#include "GKVorticity.H"
#include "GKPoissonF_F.H"
#include "FourthOrderUtil.H"
#include "SimpleDivergence.H"
#include "Directions.H"
#include "ToroidalBlockLevelExchangeCenter.H"

//#define VERIFY_MATRIX

#include "NamespaceHeader.H"
 
const char* GKVorticity::pp_name = {"GKVorticity"};

GKVorticity::GKVorticity(const ParmParse&   a_pp,
                         const ParmParse&   a_pp_base,
                         const MagGeom&     a_geom,
                         const Real         a_larmor_number,
                         const Real         a_debye_number,
                         const bool         a_second_order,
                         const bool         a_low_pollution,
                         const bool         a_include_pol_den_correction,
                         const bool         a_include_diffusion,
                         const std::string& a_model,
                         const bool         a_allocate_precond)
   : GKPoisson(a_pp_base, a_geom, a_larmor_number, a_debye_number, false, a_second_order, a_low_pollution),
     m_model(a_model),
     m_include_pol_den_correction(a_include_pol_den_correction),
     m_include_diffusion(a_include_diffusion),
     m_use_vorticity_bcs(false),
     m_precond_scale_fac(1.),
     m_dt_implicit(1.),
     m_beta(1.),
     m_flux_surface(a_geom),
     m_charge_exchange_func(NULL),
     m_parallel_conductivity_func(NULL)

{
   init(a_pp, a_pp_base, a_allocate_precond);
}
      
GKVorticity::~GKVorticity()
{
}

void
GKVorticity::init( const ParmParse&  a_pp,
                   const ParmParse&  a_pp_base,
                   const bool        a_allocate_precond)
{
   
   m_verbosity = false;
   
   parseParameters( a_pp );

   if (m_verbosity>0) {
      printParameters();
   }

   // Get distjoint box layout
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // Assign parallel conductivity profile
   m_parallel_cond_face.define( grids, 1, IntVect::Zero );
   m_parallel_cond_cell.define( grids, 1, IntVect::Zero );
   if (m_parallel_conductivity_func != NULL) {
      LevelData<FArrayBox> parallel_cond_cell( grids, 1, 2*IntVect::Unit );
      m_parallel_conductivity_func->assign( parallel_cond_cell, m_geometry, 0.0);
      m_geometry.fillInternalGhosts( parallel_cond_cell );
      fourthOrderCellToFaceCenters(m_parallel_cond_face, parallel_cond_cell);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_parallel_cond_cell[dit].copy(parallel_cond_cell[dit]);
      }
   }

   // Assign charge exchange profile
   m_charge_exchange_face.define( grids, 1, IntVect::Zero );
   if (m_charge_exchange_func != NULL ) {
      LevelData<FArrayBox> charge_exchange_cell( grids, 1, 2*IntVect::Unit );
      m_charge_exchange_func->assign( charge_exchange_cell, m_geometry, 0.0);
      m_geometry.fillInternalGhosts( charge_exchange_cell );
      fourthOrderCellToFaceCenters(m_charge_exchange_face, charge_exchange_cell);
   }
   
   // Create objects to hold plasma parameters
   m_mass_density_face.define( grids, 1, IntVect::Zero );
   m_charge_density_face.define( grids, 1, IntVect::Zero );
   m_electron_temperature_face.define( grids, 1, IntVect::Zero );
   
   // Check if we need to include high-order corrections
   // N.B. vorticity diffusion model should only solve
   // standard (2nd order) elliptic diffusion equaiton for vorticity
   m_include_high_order_corr = ((m_include_pol_den_correction || m_include_diffusion)
                                && m_model == "Vorticity") ? true : false;
   
   // Create objects to handle high-order (i.e., polarization density corrections and
   // vorticity diffusion) terms in the vorticity equation. See notes for notation.
   m_M_mapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_M_unmapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_N2_mapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_N2_unmapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   
   
   // Preconditioner initializations
   m_symmetrized_preconditioner = false;
   if ( m_include_high_order_corr ) {
      ParmParse pp( ((string)a_pp_base.prefix() + ".linear_solver.precond").c_str());
      if ( pp.contains("method") ) {
         string method;
         pp.get("method", method);
         if ( method != "MGR" ) {
            MayDay::Warning("GKVorticity::init(): include_pol_den_correction option needs MGR preconditioner");
         }
      }
      if (pp.contains("symmetrized")) {
         pp.get( "symmetrized", m_symmetrized_preconditioner );
      }
      else {
         m_symmetrized_preconditioner = false;
      }

      if ( m_symmetrized_preconditioner ) {
         // This option is highly experimental
         cout << "GKVorticity::init(): m_symmetrized_preconditioner is true!" << endl;
         exit(1);
      }
      
      if ( m_symmetrized_preconditioner ) {
         m_ne_over_Te.define(grids, 1, IntVect::Zero);
      }

      if ( pp.contains("precond_scale_fac") ) {
         pp.get("precond_scale_fac", m_precond_scale_fac);
      }
   }

   if ( m_preconditioner == NULL && a_allocate_precond) {

      int nvar = m_include_high_order_corr? 2: 1;
      m_preconditioner = allocatePreconditioner(m_geometry, m_discretization_order, nvar, m_mblx_ptr);

      ParmParse pp_precond( ((string)a_pp_base.prefix() + ".linear_solver.precond").c_str());
      m_preconditioner->setMethodParams(pp_precond);
      m_preconditioner->setConvergenceParams(pp_precond);
   }

#ifdef VERIFY_MATRIX
   m_subtract_fs_par_div = false;
#else
   m_subtract_fs_par_div = m_geometry.shearedMBGeom();
#endif

}


MBHypreSolver* GKVorticity::allocatePreconditioner( const MagGeom&                  a_geom,
                                                    const int                       a_discretization_order,
                                                    const int                       a_nvar,
                                                    MultiBlockLevelExchangeCenter*  a_mblx_ptr )
{
   MBHypreSolver* solver(NULL);
   solver = new MBHypreSolver(a_geom, a_nvar, a_discretization_order, a_mblx_ptr);
   CH_assert(solver != NULL);
   return solver;
}


void
GKVorticity::setVorticityOperatorCoefficients(const LevelData<FArrayBox>&  a_ion_mass_density,
                                              const LevelData<FArrayBox>&  a_ion_charge_density,
                                              const LevelData<FArrayBox>&  a_electron_temperature,
                                              EllipticOpBC&                a_bc,
                                              const bool                   a_update_preconditioner )
{
   /*
    This function sets base class (EllipticOp) BC data (m_bc, boundary stencils,
    and Neumann/Natural bc factors), which is needed for computeFluxDivergence()
    calls from the ImEx Vorticity model; also updates a preconditioner.
    */
   
   
   // Compute tensor coefficients
   computeVorticityCoefficients(a_ion_mass_density,
                                a_ion_charge_density,
                                a_electron_temperature,
                                a_bc,
                                m_mapped_coefficients,
                                m_unmapped_coefficients);
   
   // Set base class m_bc, boundary stencils and Neumann/Natural bc factors
   // needed for an inhomogeneous BC application of the computeFluxDivergence().
   // N.B. because this GKVorticity class only handles
   // a homogeneous/preconditioner solve we don't need to update m_bc_divergence
   // member that stores inhomogeneous RHS contributions; thus we pass
   // homogeneous_bcs = true boolean
   updateBoundaryData(m_unmapped_coefficients, a_bc, true);
   
   // Update preconditioner
   // NB: MBSolver calss that handles preconditioner
   // construct boundary stencils as part of hypre matrix construction
   if ( a_update_preconditioner ) {
      updatePreconditioner(m_preconditioner, m_volume_reciprocal, m_mapped_coefficients, a_bc);
   }
}


void
GKVorticity::computeVorticityCoefficients(const LevelData<FArrayBox>& a_ion_mass_density,
                                          const LevelData<FArrayBox>& a_ion_charge_density,
                                          const LevelData<FArrayBox>& a_electron_temperature,
                                          const EllipticOpBC&         a_bc,
                                          LevelData<FluxBox>&         a_mapped_coefficients,
                                          LevelData<FluxBox>&         a_unmapped_coefficients )
{
   CH_TIME("GKVorticity::computeCoefficients");
   CH_assert(a_mapped_coefficients.ghostVect() >= IntVect::Unit);
   CH_assert(a_unmapped_coefficients.ghostVect() >= IntVect::Unit);

   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   convertCellToFace(m_mass_density_face, a_ion_mass_density);
   convertCellToFace(m_charge_density_face, a_ion_charge_density);
   convertCellToFace(m_electron_temperature_face, a_electron_temperature);
   
   //Update parallel conductivity
   if (m_parallel_conductivity_func != NULL) {
      updateParallelConductivity(m_electron_temperature_face, m_parallel_cond_face);
   }
   
   const LevelData<FluxBox>& BFieldMag = m_geometry.getFCBFieldMag();
   
   const LevelData<FluxBox>& perp_coeff_unmapped = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& par_coeff_unmapped = m_geometry.getEllipticOpParCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();
      
   for (DataIterator dit(grids); dit.ok(); ++dit) {

      const Box& box = grids[dit];
      
      FluxBox tmp_perp_unmapped(box, SpaceDim * SpaceDim);
      tmp_perp_unmapped.copy(perp_coeff_unmapped[dit]);
      
      FluxBox tmp_par_unmapped(box, SpaceDim * SpaceDim);
      tmp_par_unmapped.copy(par_coeff_unmapped[dit]);
      
      FluxBox tmp_perp_mapped(box, SpaceDim * SpaceDim);
      tmp_perp_mapped.copy(perp_coeff_mapped[dit]);

      FluxBox tmp_par_mapped(box, SpaceDim * SpaceDim);
      tmp_par_mapped.copy(par_coeff_mapped[dit]);

      if (m_model == "Vorticity") {
      
         FluxBox polarization_fac(box, 1);
         polarization_fac.setVal(1.0);
         polarization_fac *= m_larmor_number2;
         polarization_fac.mult(m_mass_density_face[dit],box, 0, 0);
         polarization_fac.divide(BFieldMag[dit], box, 0, 0);
         polarization_fac.divide(BFieldMag[dit], box, 0, 0);
         
         FluxBox perp_fac(box, 1);
         perp_fac.copy(polarization_fac);
         perp_fac *= m_alpha;
         
         if (m_charge_exchange_func != NULL) {
            FluxBox perp_fac_chx(box, 1);
            perp_fac_chx.copy(m_charge_exchange_face[dit]);
            perp_fac_chx *= m_dt_implicit;
            perp_fac += perp_fac_chx;
         }
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_perp_unmapped.mult(perp_fac, box, 0, n);
            tmp_par_unmapped.mult(m_parallel_cond_face[dit], box, 0, n);
            
            tmp_perp_mapped.mult(perp_fac, box, 0, n);
            tmp_par_mapped.mult(m_parallel_cond_face[dit], box, 0, n);
         }

         tmp_par_unmapped *= m_dt_implicit;
         tmp_par_mapped *= m_dt_implicit;

         a_unmapped_coefficients[dit].copy(tmp_perp_unmapped);
         a_unmapped_coefficients[dit] += tmp_par_unmapped;

         a_mapped_coefficients[dit].copy(tmp_perp_mapped);
         a_mapped_coefficients[dit] += tmp_par_mapped;

         // Handle high-order terms
         m_M_mapped[dit].setVal(0.);
         m_M_unmapped[dit].setVal(0.);
         m_N2_mapped[dit].setVal(0.);
         m_N2_unmapped[dit].setVal(0.);
         
         if (m_include_high_order_corr) {
            m_M_mapped[dit] += perp_coeff_mapped[dit];
            m_M_unmapped[dit] += perp_coeff_unmapped[dit];
            for (int n=0; n<SpaceDim*SpaceDim; ++n) {
               m_M_mapped[dit].mult(polarization_fac, box, 0, n);
               m_M_unmapped[dit].mult(polarization_fac, box, 0, n);
            }
         
            // Include polarization corrections
            if (m_include_pol_den_correction) {
               m_N2_mapped[dit] += tmp_par_mapped;
               m_N2_unmapped[dit] += tmp_par_unmapped;
           
               for (int n=0; n<SpaceDim*SpaceDim; ++n) {
                  m_N2_mapped[dit].mult(m_electron_temperature_face[dit],box, 0, n);
                  m_N2_unmapped[dit].mult(m_electron_temperature_face[dit],box, 0, n);
               
                  m_N2_mapped[dit].divide(m_charge_density_face[dit], box, 0, n);
                  m_N2_unmapped[dit].divide(m_charge_density_face[dit], box, 0, n);
               }
            }
         
            // Include anomalous viscosity term
            if (m_include_diffusion) {

               FluxBox diff_coeff_unmapped(box, SpaceDim * SpaceDim);
               FluxBox diff_coeff_mapped(box, SpaceDim * SpaceDim);

               diff_coeff_unmapped.copy(m_diffusion_coeff_unmapped[dit]);
               diff_coeff_mapped.copy(m_diffusion_coeff_mapped[dit]);

               diff_coeff_unmapped *= m_dt_implicit;
               diff_coeff_mapped *= m_dt_implicit;
            
               m_N2_unmapped[dit] += diff_coeff_unmapped;
               m_N2_mapped[dit] += diff_coeff_mapped;
            }
         }
      }

      else if (m_model == "ParallelCurrent") {
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_par_unmapped.mult(m_parallel_cond_face[dit], box, 0, n);
            tmp_par_unmapped.divide(m_charge_density_face[dit],box, 0, n);
            tmp_par_unmapped.mult(m_electron_temperature_face[dit],box, 0, n);
            
            tmp_par_mapped.mult(m_parallel_cond_face[dit], box, 0, n);
            tmp_par_mapped.divide(m_charge_density_face[dit],box, 0, n);
            tmp_par_mapped.mult(m_electron_temperature_face[dit],box, 0, n);
         }

         tmp_par_unmapped *= m_dt_implicit;
         tmp_par_mapped *= m_dt_implicit;
         
         a_unmapped_coefficients[dit].copy(tmp_par_unmapped);
         a_mapped_coefficients[dit].copy(tmp_par_mapped);
      }

      else if (m_model == "VorticityDiffusion") {
         
         a_unmapped_coefficients[dit].setVal(0.);
         a_mapped_coefficients[dit].setVal(0.);
         
         if (m_include_pol_den_correction) {
            for (int n=0; n<SpaceDim*SpaceDim; ++n) {
               tmp_par_unmapped.mult(m_parallel_cond_face[dit], box, 0, n);
               tmp_par_unmapped.divide(m_charge_density_face[dit],box, 0, n);
               tmp_par_unmapped.mult(m_electron_temperature_face[dit],box, 0, n);
               
               tmp_par_mapped.mult(m_parallel_cond_face[dit], box, 0, n);
               tmp_par_mapped.divide(m_charge_density_face[dit],box, 0, n);
               tmp_par_mapped.mult(m_electron_temperature_face[dit],box, 0, n);
            }
            a_unmapped_coefficients[dit] += tmp_par_unmapped;
            a_mapped_coefficients[dit] += tmp_par_mapped;
         }
         
         if (m_include_diffusion) {
            a_unmapped_coefficients[dit] += m_diffusion_coeff_unmapped[dit];
            a_mapped_coefficients[dit] += m_diffusion_coeff_mapped[dit];
         }
      }
      
      else {
         MayDay::Error("GKVorticity:: unknown model is specified");
      }
   }

   modifyForNeumannAndNaturalBCs(a_bc, a_unmapped_coefficients, a_mapped_coefficients);
   
   if (m_model == "Vorticity" && m_include_high_order_corr) {
      modifyForNeumannAndNaturalBCs(a_bc, m_M_unmapped, m_M_mapped);
      
      if (m_use_vorticity_bcs) {
         modifyForNeumannAndNaturalBCs(*m_vorticity_bc, m_N2_unmapped, m_N2_mapped);
      }
   }

   if (!m_second_order) {

      // Fill the ghost cell data needed for conversion to face averages
      //   (a) Make sure the interior ghost cell values are consistent
      a_unmapped_coefficients.exchange();
      a_mapped_coefficients.exchange();
      //   (b) Extrapolate to transverse physical boundaries
      m_geometry.fillTransverseGhosts(a_unmapped_coefficients, false);
      m_geometry.fillTransverseGhosts(a_mapped_coefficients, false);
   
      // Convert from face-centered to face-averaged
      fourthOrderAverage(a_unmapped_coefficients);
      fourthOrderAverage(a_mapped_coefficients);

      // Update the ghosts again
      a_unmapped_coefficients.exchange();
      a_mapped_coefficients.exchange();
      m_geometry.fillTransverseGhosts(a_unmapped_coefficients, false);
      m_geometry.fillTransverseGhosts(a_mapped_coefficients, false);
   }
   
   // We don't seem to use symmetrized version, yet keep the code just in case
   if ( m_symmetrized_preconditioner ) {
      const DisjointBoxLayout& grids( m_geometry.grids() );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_ne_over_Te[dit].copy(a_ion_charge_density[dit]);
         m_ne_over_Te[dit] /= a_electron_temperature[dit];
      }
   }
   
}

void GKVorticity::setDiffusionOperatorCoeff(const LevelData<FluxBox>& a_pol_diffusion_face,
                                            const LevelData<FluxBox>& a_perp_diffusion_face)
{
   
   // Computes tensor diffusion coefficents
  
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   if (!m_diffusion_coeff_unmapped.isDefined()) {
      m_diffusion_coeff_unmapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   }
   if (!m_diffusion_coeff_mapped.isDefined()) {
     m_diffusion_coeff_mapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   }
   
   LevelData<FluxBox> poloidal_coeff(grids, SpaceDim*SpaceDim, IntVect::Unit);
   LevelData<FluxBox> poloidal_coeff_mapped(grids, SpaceDim*SpaceDim, IntVect::Unit);

   m_geometry.getCustomEllipticOpCoeff(poloidal_coeff, "poloidal");
   m_geometry.getCustomEllipticOpCoeffMapped(poloidal_coeff_mapped, "poloidal");
   
   const LevelData<FluxBox>& perp_coeff_unmapped = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
            
      FluxBox& this_coeff = m_diffusion_coeff_unmapped[dit];
      FluxBox& this_coeff_mapped = m_diffusion_coeff_mapped[dit];
      
      this_coeff.setVal(0.);
      this_coeff_mapped.setVal(0.);
      
      const Box& box = grids[dit];
      
      FluxBox diffusion_coeff_tmp(box, SpaceDim * SpaceDim);
      FluxBox diffusion_coeff_tmp_mapped(box, SpaceDim * SpaceDim);

      // Add poloidal diffusion
      diffusion_coeff_tmp.copy(poloidal_coeff[dit]);
      diffusion_coeff_tmp_mapped.copy(poloidal_coeff_mapped[dit]);

      for (int n=0; n<SpaceDim*SpaceDim; ++n) {
         diffusion_coeff_tmp.mult(a_pol_diffusion_face[dit], box, 0, n);
         diffusion_coeff_tmp_mapped.mult(a_pol_diffusion_face[dit], box, 0, n);
      }
         
      this_coeff += diffusion_coeff_tmp;
      this_coeff_mapped += diffusion_coeff_tmp_mapped;

      // Add perpendicular diffusion
      diffusion_coeff_tmp.copy(perp_coeff_unmapped[dit]);
      diffusion_coeff_tmp_mapped.copy(perp_coeff_mapped[dit]);

      for (int n=0; n<SpaceDim*SpaceDim; ++n) {
         diffusion_coeff_tmp.mult(a_perp_diffusion_face[dit], box, 0, n);
         diffusion_coeff_tmp_mapped.mult(a_perp_diffusion_face[dit], box, 0, n);
      }

      this_coeff += diffusion_coeff_tmp;
      this_coeff_mapped += diffusion_coeff_tmp_mapped;
   }
}

void
GKVorticity::setVorticityBc(EllipticOpBC& a_bc)
{
   m_vorticity_bc = a_bc.clone();
   m_use_vorticity_bcs = true;
}

void
GKVorticity::updatePreconditioner( MBSolver*              a_preconditioner,
                                   LevelData<FArrayBox>&  a_volume_reciprocal,
                                   LevelData<FluxBox>&    a_mapped_coefficients,
                                   const EllipticOpBC&    a_potential_bc )
{
   CH_TIME("GKVorticity::updatePreconditioner");

   const DisjointBoxLayout& grids = m_geometry.grids();

   a_preconditioner->initializeMatrix();

   if (m_include_high_order_corr) {

      LevelData<FluxBox> MplusN1(grids, SpaceDim*SpaceDim, IntVect::Unit);
      LevelData<FluxBox> negativeN2(grids, SpaceDim*SpaceDim, IntVect::Unit);
      LevelData<FluxBox> M(grids, SpaceDim*SpaceDim, IntVect::Unit);
      LevelData<FluxBox> zero(grids, SpaceDim*SpaceDim, IntVect::Unit);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         MplusN1[dit].copy(a_mapped_coefficients[dit]);
         M[dit].copy(m_M_mapped[dit]);
         if ( m_precond_scale_fac != 1. ) {
            MplusN1[dit] *= m_precond_scale_fac;
            M[dit] *= m_precond_scale_fac;
         }
         negativeN2[dit].copy(m_N2_mapped[dit]);
         negativeN2[dit].negate();

         zero[dit].setVal(0.);
      }

      // Upper-left block
      a_preconditioner->constructMatrixBlock(0, 0, a_volume_reciprocal, MplusN1, a_potential_bc);

      if ( m_symmetrized_preconditioner ) {

         // Upper-right block
         a_preconditioner->constructMatrixBlock(0, 1, a_volume_reciprocal, M, a_potential_bc);

         // Lower-left block
         a_preconditioner->constructMatrixBlock(1, 0, a_volume_reciprocal, M, a_potential_bc);

         // Lower-right block
         a_preconditioner->constructMatrixBlock(1, 1, a_volume_reciprocal, M, m_ne_over_Te, a_potential_bc);
      }
      else {

         LevelData<FArrayBox> Te_reciprocal(grids, 1, IntVect::Zero);
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            Te_reciprocal[dit].setVal(1.);
            // Don't need the 1/T_e factor (yet), since Te is included in the coefficient
            // Te_reciprocal[dit] /= m_electron_temperature_cell[dit];
         }

         
         RefCountedPtr<EllipticOpBC> vorticity_bc;
         
         if (m_use_vorticity_bcs) {
            vorticity_bc = m_vorticity_bc->clone();
            // Even though preconditioner solves a homogeneous problem
            // wich does not require these factors, the MBSolver requires
            // some data to be initialized
            setNeumannNaturalFactor(negativeN2, *vorticity_bc);
         }
         else {
            bool use_extrapolated_bcs = true;
            vorticity_bc = a_potential_bc.clone(use_extrapolated_bcs);
         }
         
         // Upper-right block
         a_preconditioner->constructMatrixBlock(0, 1, a_volume_reciprocal, negativeN2, *vorticity_bc);

         // Lower-left block
         a_preconditioner->constructMatrixBlock(1, 0, a_volume_reciprocal, M, a_potential_bc);

         // Lower-right block
         a_preconditioner->constructMatrixBlock(1, 1, zero, Te_reciprocal, *vorticity_bc);
      }
   }
   
   else if (m_model == "VorticityDiffusion") {
      LevelData<FArrayBox> beta(grids, 1, IntVect::Zero);
      for (DataIterator dit(beta.dataIterator()); dit.ok(); ++dit ) {
         beta[dit].setVal(m_beta);
      }
      a_preconditioner->constructMatrixBlock(0, 0, a_volume_reciprocal, a_mapped_coefficients, beta, a_potential_bc);
   }
   
   else {
      a_preconditioner->constructMatrixBlock(0, 0, a_volume_reciprocal, a_mapped_coefficients, a_potential_bc);
   }

   a_preconditioner->finalizeMatrix();

}

void
GKVorticity::solvePreconditioner( const LevelData<FArrayBox>& a_r,
                                  LevelData<FArrayBox>&       a_z )
{
   CH_TIME("GKVorticity::solvePreconditioner");

   if ( m_include_high_order_corr ) {
      const DisjointBoxLayout& grids = m_geometry.grids();

      LevelData<FArrayBox> block_in(grids, 2, IntVect::Zero);
      LevelData<FArrayBox> block_out(grids, 2, IntVect::Zero);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         block_in[dit].copy(a_r[dit],0,0,1);
         block_in[dit].setVal(0.,1);
         block_out[dit].copy(a_z[dit],0,0,1);
         block_out[dit].setVal(0.,1);

         if ( m_precond_scale_fac != 1. ) {
            block_out[dit].mult(1./m_precond_scale_fac,0);
         }
      }

      m_preconditioner->solve(block_in, block_out, true);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_z[dit].copy(block_out[dit],0,0,1);

         if ( m_precond_scale_fac != 1. ) {
            a_z[dit] *= m_precond_scale_fac;
         }
      }

      if ( m_symmetrized_preconditioner ) {
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            a_z[dit].plus(block_out[dit],1,0,1);
         }
      }
   }
   else {
      m_preconditioner->solve(a_r, a_z, true);
   }
}

void
GKVorticity::computeFluxDivergence(const LevelData<FArrayBox>& a_in,
                                   LevelData<FArrayBox>&       a_out,
                                   const bool                  a_homogeneous_bcs,
                                   const bool                  a_extrap_to_ghosts)
{
  /*
   Same as computeFluxDivergence() function of the base class,
   but also contains the option to subtruct fs-average of parallel flux
   
   Can't we use the base class function for most of the code here?
   No we can't becasue we need intermidiate calculation of flux, which
   the base class version does not provide.
   */
  CH_TIME("GKVorticity::computeFluxDivergence");

  const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

  IntVect ghosts = a_in.ghostVect();
  LevelData<FArrayBox> phi(grids, 1, ghosts);

  for (DataIterator dit(grids); dit.ok(); ++dit) {
     phi[dit].copy(a_in[dit]);
  }

  LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);

  if (!m_low_pollution) {

    if (SpaceDim == 3) {
      if (a_extrap_to_ghosts) computeField(phi, flux);
      else compute3DFieldWithBCs(phi, flux, a_homogeneous_bcs);

    }
    else {
      if (a_extrap_to_ghosts) computePoloidalField(phi, flux);
      else computePoloidalFieldWithBCs(phi, flux, a_homogeneous_bcs);
    }

    // Multiply the field by the unmapped, face-centered GKP coefficients
    multiplyCoefficients(flux, false);

    m_geometry.fillTransversePhysicalGhosts(flux);

    m_geometry.applyAxisymmetricCorrection(flux);

    // Convert to face-averaged
    if (!m_second_order) fourthOrderAverage(flux);

    m_geometry.computeMappedGridDivergence(flux, a_out, !m_second_order);
  }

  else {

    if (SpaceDim == 3) {
      if (a_extrap_to_ghosts) computeMapped3DField(phi, flux);
      else computeMapped3DFieldWithBCs(phi, flux, a_homogeneous_bcs);
   }

    else  {
      if (a_extrap_to_ghosts) computeMappedPoloidalField(phi, flux);
      else computeMappedPoloidalFieldWithBCs(phi, flux, a_homogeneous_bcs);
    }

    // Multiply the field by the mapped, face-centered GKP coefficients
    multiplyCoefficients(flux, true);

    m_geometry.fillTransversePhysicalGhosts(flux);

    // Convert to face-averaged
    if (!m_second_order) fourthOrderAverage(flux);

    m_geometry.averageAtBlockBoundaries(flux);

    LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
    for (DataIterator dit(flux.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect faceArea = block_coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; ++dir) {
         NTF_normal[dit][dir].copy(flux[dit][dir],dir,0,1);
         NTF_normal[dit][dir].mult(faceArea[dir]);
      }
    }

    RealVect fakeDx = RealVect::Unit;
    for (DataIterator dit(grids); dit.ok(); ++dit) {
      simpleDivergence(a_out[dit], NTF_normal[dit], grids[dit], fakeDx);
    }

  }

  for (DataIterator dit(grids); dit.ok(); ++dit) {
     a_out[dit] /= m_volume[dit];
  }

  if (m_subtract_fs_par_div) {
     subtractFSAverParDiv(a_out, flux);
  }
}

void
GKVorticity::subtractFSAverParDiv(LevelData<FArrayBox>&      a_div,
                                  const  LevelData<FluxBox>& a_flux) const
{
   
   const DisjointBoxLayout& grids = a_div.disjointBoxLayout();
   
   LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
   for (DataIterator dit(a_flux.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect faceArea = block_coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; ++dir) {
         if (a_flux.nComp() == 1) {
            NTF_normal[dit][dir].copy(a_flux[dit][dir]);
         }
         else {
            NTF_normal[dit][dir].copy(a_flux[dit][dir],dir,0,1);
         }
         NTF_normal[dit][dir].mult(faceArea[dir]);
         
         //zero out flux on radial faces
         if (dir == RADIAL_DIR) {
            NTF_normal[dit][dir].setVal(0.);
         }
      }
   }
   
   RealVect fakeDx = RealVect::Unit;
   LevelData<FArrayBox> parallel_div(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      simpleDivergence(parallel_div[dit], NTF_normal[dit], grids[dit], fakeDx);
   }
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      parallel_div[dit] /= m_volume[dit];
   }

   //subtract parallel divergence
   LevelData<FArrayBox> par_div_fs(grids, 1, IntVect::Zero);
   m_flux_surface.averageAndSpread(parallel_div, par_div_fs);

   //Subtract only on closed-flux-surfaces
   const MagCoordSys& coords = *m_geometry.getCoordSys();
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      int block_number = coords.whichBlock(grids[dit]);
      
      if ((typeid(coords) != typeid(SingleNullCoordSys)) ||
      ((const SingleNullCoordSys&)coords).isCORE(block_number))  {

         a_div[dit] -= par_div_fs[dit];
      }
   }
}

void
GKVorticity::updateParallelConductivity(const LevelData<FluxBox>& a_Te,
                                        LevelData<FluxBox>&       a_parallel_conductivity ) const
{
   if (m_parallel_conductivity_func != NULL) {
      for (DataIterator dit(a_parallel_conductivity.dataIterator()); dit.ok(); ++dit) {
         a_parallel_conductivity[dit].copy(m_parallel_cond_face[dit]);
      }
   }
   else {
      computeParallelConductivity(a_Te, a_parallel_conductivity);
   }
}

void
GKVorticity::updateParallelConductivity(const LevelData<FArrayBox>&  a_Te,
                                        LevelData<FArrayBox>&        a_parallel_conductivity ) const
{
   if (m_parallel_conductivity_func != NULL) {
      for (DataIterator dit(a_parallel_conductivity.dataIterator()); dit.ok(); ++dit) {
         a_parallel_conductivity[dit].copy(m_parallel_cond_cell[dit]);
      }
   }
   else {
      computeParallelConductivity(a_Te, a_parallel_conductivity);
   }
}

void
GKVorticity::computeParallelConductivity(const LevelData<FluxBox>& a_Te,
                                         LevelData<FluxBox>&       a_parallel_conductivity ) const

{
   //Get normalization parameters (units)
   double N, T, L;
   ParmParse ppunits( "units" );
   ppunits.get("number_density",N);  //[m^{-3}]
   ppunits.get("temperature",T);     //[eV]
   ppunits.get("length",L);          //[m]
   
   double Ncgs = 1.0e-6 * N;    //[cm^{-3}]
   double Tcgs = 1.602e-12 * T; //[erg]
   double Lcgs  = 1.0e2 * L;   //[cm]
   
   //Universal constants (in CGS)
   double mp = 1.6726e-24;
   double me = 9.1094e-28;
   double ech = 4.8032e-10;
   
   //In computing Coulumb logarithm we neglect variations of n and T
   //relative to thier corresponding unit values
   double Coulomb_Lg_ee = 23.5 - log( sqrt(Ncgs)/ pow(T, 5.0/4.0) )
                          - sqrt(1.0e-5 + 1.0/16.0 * pow(log(T) - 2.0, 2));

   //Braginskii tau_e [s] and sigma_parallel [1/s] (computed for the "units" values)
   double tau_e = 3.44 * 1.0e5 * pow(T, 3.0/2.0) / (Ncgs * Coulomb_Lg_ee);
   double sigma_e = ech * ech * Ncgs * tau_e / (0.51 * me);

   //Normalization coefficient
   double Vunit = sqrt(Tcgs/mp);
   double norm_coeff = Tcgs / (ech*ech * Vunit * Ncgs * Lcgs);
   double sigma_norm = sigma_e * norm_coeff;
   
   //In case we need to limit the value of sigma for numerical purposes
   double sigma_max = m_parallel_conductivity_limit;
   
   //Compute parallel conductivity
   const DisjointBoxLayout& grids = a_parallel_conductivity.disjointBoxLayout();
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         
         FArrayBox& this_sigma = a_parallel_conductivity[dit][dir];
         const FArrayBox& this_T = a_Te[dit][dir];
      
         FORT_COMPUTE_CONDUCTIVITY(CHF_BOX(this_sigma.box()),
                                  CHF_FRA1(this_sigma,0),
                                  CHF_CONST_FRA1(this_T,0),
                                  CHF_CONST_REAL(sigma_norm),
                                  CHF_CONST_REAL(sigma_max));
      
      }
   }
}

void
GKVorticity::computeParallelConductivity(const LevelData<FArrayBox>& a_Te,
                                         LevelData<FArrayBox>&       a_parallel_conductivity ) const

{
   //Get normalization parameters (units)
   double N, T, L;
   ParmParse ppunits( "units" );
   ppunits.get("number_density",N);  //[m^{-3}]
   ppunits.get("temperature",T);     //[eV]
   ppunits.get("length",L);          //[m]
   
   double Ncgs = 1.0e-6 * N;    //[cm^{-3}]
   double Tcgs = 1.602e-12 * T; //[erg]
   double Lcgs  = 1.0e2 * L;   //[cm]
   
   //Universal constants (in CGS)
   double mp = 1.6726e-24;
   double me = 9.1094e-28;
   double ech = 4.8032e-10;
   
   //In computing Coulumb logarithm we neglect variations of n and T
   //relative to thier corresponding unit values
   double Coulomb_Lg_ee = 23.5 - log( sqrt(Ncgs)/ pow(T, 5.0/4.0) )
                          - sqrt(1.0e-5 + 1.0/16.0 * pow(log(T) - 2.0, 2));

   //Braginskii tau_e [s] and sigma_parallel [1/s] (computed for the "units" values)
   double tau_e = 3.44 * 1.0e5 * pow(T, 3.0/2.0) / (Ncgs * Coulomb_Lg_ee);
   double sigma_e = ech * ech * Ncgs * tau_e / (0.51 * me);

   //Normalization coefficient
   double Vunit = sqrt(Tcgs/mp);
   double norm_coeff = Tcgs / (ech*ech * Vunit * Ncgs * Lcgs);
   double sigma_norm = sigma_e * norm_coeff;
   
   //In case we need to limit the value of sigma for numerical purposes
   double sigma_max = m_parallel_conductivity_limit;
   
   //Compute parallel conductivity
   const DisjointBoxLayout& grids = a_parallel_conductivity.disjointBoxLayout();
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
         
      FArrayBox& this_sigma = a_parallel_conductivity[dit];
      const FArrayBox& this_T = a_Te[dit];
      
      FORT_COMPUTE_CONDUCTIVITY(CHF_BOX(this_sigma.box()),
                                CHF_FRA1(this_sigma,0),
                                CHF_CONST_FRA1(this_T,0),
                                CHF_CONST_REAL(sigma_norm),
                                CHF_CONST_REAL(sigma_max));
      
   }
}

void
GKVorticity::convertCellToFace(LevelData<FluxBox>&          a_face_data,
                               const LevelData<FArrayBox>&  a_cell_data ) const
{
   CH_assert(m_num_potential_ghosts >= 2);
   
   const DisjointBoxLayout& grids = a_cell_data.disjointBoxLayout();
   
   LevelData<FArrayBox> cell_data_grown( grids, 1, m_num_potential_ghosts*IntVect::Unit );
      
   // Extrapolate to two codim 1 ghost cell layers at all block boundaries, including
   // physical boundaries
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = block_coord_sys.domain();
      
      cell_data_grown[dit].copy(a_cell_data[dit]);
      fourthOrderCellExtrapAtDomainBdry(cell_data_grown[dit], domain, grids[dit]);
   }
   
   // Fill extrablock ghosts and perform an exchange.  We then have m_num_potential_ghosts
   // layers of ghost cells filled, except at the physical boundary were we only have 2
   // codim 1 cells.
   m_geometry.fillInternalGhosts( cell_data_grown );
      
   // Average from cell averages to face centers.  We use fourth-order averaging, even
   // if m_second_order is true.  This operation only involves data from the cells along
   // the same direction (i.e., does not involve codim 2 values, and therefore it is safe
   // to use even for the second-order (since fourthOrderCellExtrapAtDomainBdry() call
   // provided all necessary information).

   fourthOrderCellToFaceCenters(a_face_data, cell_data_grown);
}

void
GKVorticity::parseParameters( const ParmParse&   a_pp )
{

   if (a_pp.contains("charge_exchange_coefficient")) {
      GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
      std::string grid_function_name;
      a_pp.get("charge_exchange_coefficient", grid_function_name );
      m_charge_exchange_func = grid_library->find( grid_function_name );
   }
   
   if (a_pp.contains("parallel_conductivity")) {
      GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
      std::string grid_function_name;
      a_pp.get("parallel_conductivity", grid_function_name );
      m_parallel_conductivity_func = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("parallel_conductivity_limit")) {
     a_pp.get("parallel_conductivity_limit",m_parallel_conductivity_limit );
   }
   else {
     m_parallel_conductivity_limit = DBL_MAX;
   }

}


void
GKVorticity::printParameters()
{
   if (procID()==0) {
      std::cout << "GKVorticity model: " << m_model << std::endl;
   }
}

#include "NamespaceFooter.H"
