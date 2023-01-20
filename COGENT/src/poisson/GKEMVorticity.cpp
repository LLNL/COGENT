#include "GKEMVorticity.H"
#include "Directions.H"
#include "SimpleDivergence.H"
#include "CONSTANTS.H"
#include <unistd.h>
#include "CellToEdge.H"
#include "EllipticOpBCFactory.H"

#undef ALFVEN_WAVE_TEST

#include "NamespaceHeader.H"
 
const char* GKEMVorticity::pp_name = {"GKEMVorticity"};


GKEMVorticity::GKEMVorticity(const ParmParse&   a_pp,
                             const ParmParse&   a_pp_base,
                             const MagGeom&     a_geom,
                             const Real         a_larmor_number,
                             const Real         a_ion_skin_depth,
                             const Real         a_electron_skin_depth,
                             const bool         a_second_order,
                             const bool         a_low_pollution,
                             const bool         a_include_pol_den_correction,
                             const bool         a_include_diffusion,
                             const std::string& a_pc_model)
   : GKVorticity(a_pp,
                 a_pp_base,
                 a_geom,
                 a_larmor_number,
                 0.,
                 a_second_order,
                 a_low_pollution,
                 a_include_pol_den_correction,
                 a_include_diffusion,
                 "Vorticity",
                 false),
     m_pc_Apar(NULL),
     m_neg_Laplacian_perp(NULL),
     m_divb_solver(NULL),
     m_pc_model(a_pc_model),
     m_pc_updated(false)
{
   m_verbosity = false;
    
   parseParameters( a_pp );

   if (m_verbosity>0) {
      printParameters();
   }
   
   m_ion_skin_depth_sq = pow(a_ion_skin_depth, 2);
   m_electron_skin_depth_sq = pow(a_electron_skin_depth, 2);

   const DisjointBoxLayout& grids( m_geometry.grids() );
   m_charge_density_cell.define(grids, 1, IntVect::Zero);
   m_electron_temperature_cell.define(grids, 1, IntVect::Zero);
   
   m_charge_density_cell.define(grids, 1, IntVect::Zero);
   m_elliptic_phi_mapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_elliptic_phi_unmapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_elliptic_Apar_mapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_elliptic_Apar_unmapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);

   // Allocated preconditioners
   if ( m_preconditioner == NULL) {
      
      ParmParse pp_precond( ((string)a_pp_base.prefix() + ".linear_solver.precond").c_str());

      int nvar;
      if (m_pc_model == "AparSolnOp") {
         nvar = 1;
      }
      else if (m_pc_model == "EMVorticity") {
         nvar = 4;
      }
      else if (m_pc_model == "AlfvenWave") {
         nvar = 2;

         // Get perpendicular Laplacian solve parameters
         ParmParse pp_perp_lapl( ((string)pp_precond.prefix() + ".perp_lapl").c_str());
         
         string method;
         if ( pp_perp_lapl.query("method", method) == 0 ) method = "AMG";
         double tol;
         if ( pp_perp_lapl.query("tol", tol) == 0 ) tol = 1.e-12;
         int max_iter;
         if ( pp_perp_lapl.query("max_iter", max_iter) == 0 ) max_iter = 100;
         bool verbose;
         if ( pp_perp_lapl.query("verbose", verbose) == 0 ) verbose = false;

         // Construct the negative perpendicular Laplacian
         m_neg_Laplacian_perp = allocatePreconditioner(m_geometry, m_discretization_order, 1, m_mblx_ptr);
         if ( method == "AMG" ) {
            m_neg_Laplacian_perp->setBasicAMGMethod();
         }
         else {
            m_neg_Laplacian_perp->setMethod(method);
         }
         m_neg_Laplacian_perp->setConvergenceParams( tol, max_iter, verbose, 0., 0, false);

         // Create potential BCs
         const std::string name("potential");
         const std::string prefix( "BC." + name );
         ParmParse ppsp( prefix.c_str() );
         EllipticOpBCFactory elliptic_op_bc_factory;

         EllipticOpBC* potential_bcs = elliptic_op_bc_factory.create(name,
                                                         ppsp,
                                                         *(m_geometry.getCoordSys()),
                                                         false);
    
         m_neg_Laplacian_perp->initializeMatrix();

         const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
         LevelData<FluxBox> tensor_coefficient(grids, SpaceDim*SpaceDim, IntVect::Unit);
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            tensor_coefficient[dit].copy(perp_coeff_mapped[dit]);
         }
         
         m_neg_Laplacian_perp->constructMatrixBlock(0, 0, m_volume_reciprocal, tensor_coefficient, *potential_bcs);

         m_neg_Laplacian_perp->finalizeMatrix();

         m_divb_solver = allocatePreconditioner(m_geometry, m_discretization_order, 1, m_mblx_ptr);
         m_divb_solver->setBasicAMGMethod();
         m_divb_solver->setConvergenceParams( tol, max_iter, verbose, 0., 0, false);

         m_divb_solver->initializeMatrix();

         const LevelData<FluxBox>& bunit = m_geometry.getFCBFieldDir();

         LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            flux[dit].copy(bunit[dit]);
         }
         LevelData<FluxBox> NTV_normal(grids, 1, IntVect::Unit);
         m_geometry.computeMetricTermProductAverage(NTV_normal, flux, false);

         RefCountedPtr<EllipticOpBC>  tmp_bcs = potential_bcs->clone(true);
         m_divb_solver->addAdvectionMatrixBlock(0, 0, m_volume_reciprocal, NTV_normal, *tmp_bcs);

         m_divb_solver->finalizeMatrix();

         // Store div(b)
         m_divb.define(grids, 1, IntVect::Zero);
         m_geometry.computeMappedGridDivergence(bunit, m_divb, false);

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            m_divb[dit] /= m_volume[dit];
         }
      }
      else if (m_pc_model == "ESVorticity") {
         nvar = m_include_high_order_corr? 2: 1;
         m_pc_Apar = allocatePreconditioner(m_geometry, m_discretization_order, 1, m_mblx_ptr);
         m_pc_Apar->setMethodParams(pp_precond);
         m_pc_Apar->setConvergenceParams(pp_precond);
      }

      if (m_pc_model == "AlfvenWave") {
         m_preconditioner = allocatePreconditioner(m_geometry, m_discretization_order, 2, m_mblx_ptr);
         m_preconditioner->setMethodParams(pp_precond);
         m_preconditioner->setConvergenceParams(pp_precond);
      }
      else {
         m_preconditioner = allocatePreconditioner(m_geometry, m_discretization_order, nvar, m_mblx_ptr);
         m_preconditioner->setMethodParams(pp_precond);
         m_preconditioner->setConvergenceParams(pp_precond);
      }
   }
   else {
      MayDay::Error("EMGKVorticity::init(): m_preconditioner is already allocated");
   }
}
      
GKEMVorticity::~GKEMVorticity()
{
   if (m_divb_solver) delete m_divb_solver;
   if (m_neg_Laplacian_perp) delete m_neg_Laplacian_perp;
   if (m_pc_Apar) delete m_pc_Apar;
}

void
GKEMVorticity::setJparallelBc(EllipticOpBC& a_bc)
{
   m_jpar_bc = a_bc.clone();
   m_use_jpar_bcs = true;
}

void
GKEMVorticity::setVorticityOperatorCoefficients(const LevelData<FArrayBox>&   a_ion_mass_density,
                                                const LevelData<FArrayBox>&   a_ion_charge_density,
                                                const LevelData<FArrayBox>&   a_electron_temperature,
                                                EllipticOpBC&                 a_bc,
                                                const bool                    a_update_preconditioner )
{
   /*
    This function sets base class (EllipticOp) BC data (m_bc, boundary stencils,
    and Neumann/Natural bc factors), which is needed for computeFluxDivergence()
    calls from the ImEx Vorticity model; also updates a preconditioner.
    */
    
   // Compute tensor coefficients
   computeEllipticOpCoefficients(a_ion_mass_density,
                                 a_ion_charge_density,
                                 a_electron_temperature,
                                 a_bc,
                                 "default",
                                 m_mapped_coefficients,
                                 m_unmapped_coefficients);
   
   // Set base class m_bc, boundary stencils and Neumann/Natural bc factors
   // needed for an inhomogeneous BC application of the computeFluxDivergence().
   // N.B. because this GKVorticity class only handles
   // a homogeneous/preconditioner solve we don't need to update m_bc_divergence
   // member that stores inhomogeneous RHS contributions; thus we pass
   // homogeneous_bcs = true boolean
   updateBoundaryData(m_unmapped_coefficients, a_bc, true);
   
   
   // Update a local copy of ion_charge_density
   // (needed to set-up an operator and a preconditioner)
   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_charge_density_cell[dit].copy(a_ion_charge_density[dit]);
   }
   
   // Update preconditioner
   if ( a_update_preconditioner) {
      updatePreconditioner(m_preconditioner, m_volume_reciprocal, m_mapped_coefficients, a_bc);
   }
   
   
   // TEMPORARILY: Now compute inhomogeneous contribution
   // needed to solve the physical-operator problem
   // move this to GKEMVorticityBE class when created
   computeBcDivergence(m_bc_divergence);
}


void
GKEMVorticity::setVorticityOperatorCoefficients(const LevelData<FArrayBox>&   a_ion_mass_density,
                                                const LevelData<FArrayBox>&   a_ion_charge_density,
                                                const LevelData<FArrayBox>&   a_electron_temperature,
                                                EllipticOpBC&                 a_bc_phi,
                                                EllipticOpBC&                 a_bc_Apar,
                                                const bool                    a_update_preconditioner )
{
   /*
    This function sets base class (EllipticOp) BC data (m_bc, boundary stencils,
    and Neumann/Natural bc factors), which is needed for computeFluxDivergence()
    calls from the ImEx Vorticity model; also updates a preconditioner.
    */
 
      
   // Compute coefficients for elliptic Apar operator
   computeEllipticOpCoefficients(a_ion_mass_density,
                                 a_ion_charge_density,
                                 a_electron_temperature,
                                 a_bc_Apar,
                                 "default",
                                 m_elliptic_Apar_mapped,
                                 m_elliptic_Apar_unmapped);
   
   setNeumannNaturalFactor(m_elliptic_Apar_unmapped, a_bc_Apar);

   // Compute coefficients for elliptic vorticity operator
   computeEllipticOpCoefficients(a_ion_mass_density,
                                 a_ion_charge_density,
                                 a_electron_temperature,
                                 a_bc_phi,
                                 "potential",
                                 m_elliptic_phi_mapped,
                                 m_elliptic_phi_unmapped);

   setNeumannNaturalFactor(m_elliptic_phi_unmapped, a_bc_phi);
   
   
   if (m_pc_model == "ESVorticity") {
      GKVorticity::setVorticityOperatorCoefficients(a_ion_mass_density,
                                                    a_ion_charge_density,
                                                    a_electron_temperature,
                                                    a_bc_phi,
                                                    a_update_preconditioner);
   }
   
   //Deal with diffusion coefficeints;
   //they are already computed in the base class but need to mofidy them
   //for the case of neuman or natural BCs here
   if (m_include_diffusion && m_use_vorticity_bcs) {
      modifyForNeumannAndNaturalBCs(*m_vorticity_bc, m_diffusion_coeff_unmapped, m_diffusion_coeff_mapped);
   }
   
   // Update plasma parameters
   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_charge_density_cell[dit].copy(a_ion_charge_density[dit]);
      m_electron_temperature_cell[dit].copy(a_electron_temperature[dit]);
      convertCellToFace(m_electron_temperature_face, m_electron_temperature_cell);
   }
   
   // Update conductivity
   if (m_parallel_conductivity_func != NULL) {
      updateParallelConductivity(a_electron_temperature, m_parallel_cond_cell);
   }
   
   // Update preconditioner
   if ( a_update_preconditioner) {
      if (m_pc_model == "EMVorticity") {
         updateEMPreconditioner(m_preconditioner, a_bc_phi, a_bc_Apar);
      }
      else if (m_pc_model == "ESVorticity") {
         updateESPreconditioner(m_preconditioner, m_pc_Apar, a_bc_phi, a_bc_Apar);
      }
      else if (m_pc_model == "AlfvenWave") {
         updateAlfvenWavePreconditioner(m_preconditioner, a_ion_mass_density, a_bc_phi, a_bc_Apar);
      }
   }
}

void
GKEMVorticity::computeEllipticOpCoefficients(const LevelData<FArrayBox>& a_ion_mass_density,
                                             const LevelData<FArrayBox>& a_ion_charge_density,
                                             const LevelData<FArrayBox>& a_electron_temperature,
                                             const EllipticOpBC&         a_bc,
                                             const std::string&          a_model,
                                             LevelData<FluxBox>&         a_mapped_coefficients,
                                             LevelData<FluxBox>&         a_unmapped_coefficients )
{
   CH_TIME("GKEMVorticity::computeCoefficients");

   CH_assert(a_mapped_coefficients.ghostVect() >= IntVect::Unit);
   CH_assert(a_unmapped_coefficients.ghostVect() >= IntVect::Unit);

   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   const LevelData<FluxBox>& perp_coeff_unmapped = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_unmapped_coefficients[dit].copy(perp_coeff_unmapped[dit]);
      a_mapped_coefficients[dit].copy(perp_coeff_mapped[dit]);
   }
   
   if (a_model == "potential") {

      convertCellToFace(m_mass_density_face, a_ion_mass_density);

      const LevelData<FluxBox>& BFieldMag = m_geometry.getFCBFieldMag();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& box = grids[dit];
          
         FluxBox polarization_fac(box, 1);
         polarization_fac.setVal(1.0);
         polarization_fac *= m_larmor_number2;
         polarization_fac.mult(m_mass_density_face[dit],box, 0, 0);
         polarization_fac.divide(BFieldMag[dit], box, 0, 0);
         polarization_fac.divide(BFieldMag[dit], box, 0, 0);
                      
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            a_mapped_coefficients[dit].mult(polarization_fac, box, 0, n);
            a_unmapped_coefficients[dit].mult(polarization_fac, box, 0, n);
         }
      }
   }
   
   modifyForNeumannAndNaturalBCs(a_bc, a_unmapped_coefficients, a_mapped_coefficients);
}


void
GKEMVorticity::updatePreconditioner(MBSolver*              a_preconditioner,
                                    LevelData<FArrayBox>&  a_volume_reciprocal,
                                    LevelData<FluxBox>&    a_mapped_coefficients,
                                    const EllipticOpBC&    a_potential_bc )
{
   CH_TIME("GKEMVorticity::updatePreconditioner");

   a_preconditioner->initializeMatrix();
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   LevelData<FArrayBox> alpha_coef(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      alpha_coef[dit].copy(a_volume_reciprocal[dit]);
      alpha_coef[dit].mult(m_electron_skin_depth_sq);
      alpha_coef[dit].divide(m_charge_density_cell[dit]);
   }
   
   LevelData<FArrayBox> beta(grids, 1, IntVect::Zero);
   for (DataIterator dit(beta.dataIterator()); dit.ok(); ++dit ) {
      beta[dit].setVal(1.0);
   }
   a_preconditioner->constructMatrixBlock(0, 0, alpha_coef, a_mapped_coefficients, beta, a_potential_bc);

   a_preconditioner->finalizeMatrix();

}


void
GKEMVorticity::updateEMPreconditioner(MBSolver*              a_preconditioner,
                                      const EllipticOpBC&    a_phi_bc,
                                      const EllipticOpBC&    a_Apar_bc)
{
   CH_TIME("GKEMVorticity::updateEMPreconditioner");

   const DisjointBoxLayout& grids = m_geometry.grids();
   LevelData<FluxBox> tensor_coefficient(grids, SpaceDim*SpaceDim, IntVect::Unit);
   LevelData<FArrayBox> alpha_coefficient(grids, 1, IntVect::Unit);
   LevelData<FArrayBox> beta_coefficient(grids, 1, IntVect::Unit);
   
   RefCountedPtr<EllipticOpBC> vorticity_bc;
   
   if (m_use_vorticity_bcs) {
      vorticity_bc = m_vorticity_bc->clone();
      // Even though preconditioner solves a homogeneous problem
      // wich does not require these factors, the MBSolver requires
      // some data to be initialized
      setNeumannNaturalFactor(m_diffusion_coeff_unmapped, *vorticity_bc);
   }
   else {
      bool use_extrapolated_bcs = true;
      vorticity_bc = a_phi_bc.clone(use_extrapolated_bcs);
   }

   a_preconditioner->initializeMatrix();
   
   // Construct zero blocks; BCs don't matter so use a_phi_bc
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      tensor_coefficient[dit].setVal(0.);
   }
   
   // Construct zero blocks; BCs don't matter so use a_phi_bc
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      beta_coefficient[dit].setVal(1.);
   }
   a_preconditioner->constructMatrixBlock(1, 1,  m_volume_reciprocal, tensor_coefficient, beta_coefficient, *vorticity_bc);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      beta_coefficient[dit].setVal(-1.0);
   }
   a_preconditioner->constructMatrixBlock(3, 3,  m_volume_reciprocal, tensor_coefficient, beta_coefficient, *m_jpar_bc);

   
   // Add diagonal contributions to advection blocks
   const LevelData<FArrayBox>& B_mag = m_geometry.getCCBFieldMag();
   const LevelData<FArrayBox>& grad_B = m_geometry.getCCGradBFieldMag();
   const LevelData<FArrayBox>& b_unit = m_geometry.getCCBFieldDir();
   
   LevelData<FArrayBox> advection_shift(grids, 1, IntVect::Unit);
   SpaceUtils::dotProduct(advection_shift, b_unit, grad_B);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      advection_shift[dit].mult(B_mag[dit]);
   }
   a_preconditioner->constructMatrixBlock(2, 0, tensor_coefficient, advection_shift, a_phi_bc);

   if (m_include_pol_den_correction) {

      LevelData<FArrayBox> neg_grad_Te(grids, SpaceDim, IntVect::Unit);
      computeField(m_electron_temperature_cell, neg_grad_Te);

      LevelData<FArrayBox> tmp(grids, 1, IntVect::Unit);
      SpaceUtils::dotProduct(tmp, neg_grad_Te, b_unit);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         advection_shift[dit].mult(m_electron_temperature_cell[dit]);
         advection_shift[dit].mult(m_charge_density_cell[dit]);
      
         tmp[dit].mult(-1.0);
         tmp[dit].divide(m_charge_density_cell[dit]);
      
         advection_shift[dit].plus(tmp[dit]);
         advection_shift[dit].mult(-1.0);
      }
      a_preconditioner->constructMatrixBlock(2, 1, tensor_coefficient, advection_shift, *vorticity_bc);
   }
      
   // Construct elliptic blocks for phi and vorticity
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      beta_coefficient[dit].setVal(-m_alpha);
      if (m_include_diffusion) {
         tensor_coefficient[dit].copy(m_diffusion_coeff_mapped[dit]);
      }
      else {
         tensor_coefficient[dit].setVal(0.);
      }
   }
   a_preconditioner->constructMatrixBlock(0, 1,  m_volume_reciprocal, tensor_coefficient, beta_coefficient, *vorticity_bc);
   a_preconditioner->constructMatrixBlock(1, 0,  m_volume_reciprocal, m_elliptic_phi_mapped, a_phi_bc);
   
   
   // Construct elliptic blocks for Apar
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      beta_coefficient[dit].setVal(m_alpha);
      
      const Box& box = grids[dit];
      FArrayBox alpha_Le2_over_ne(box, 1);

      alpha_Le2_over_ne.setVal(m_alpha);
      alpha_Le2_over_ne.divide(m_charge_density_cell[dit]);
      alpha_Le2_over_ne *= m_electron_skin_depth_sq;
      
#ifndef ALFVEN_WAVE_TEST
      alpha_coefficient[dit].setVal(m_ion_skin_depth_sq);
      alpha_coefficient[dit].divide(m_parallel_cond_cell[dit]);
      alpha_coefficient[dit].plus(alpha_Le2_over_ne);
#else
      alpha_coefficient[dit].copy(alpha_Le2_over_ne);
#endif
      
      alpha_coefficient[dit].mult(m_volume_reciprocal[dit]);
   }
   a_preconditioner->constructMatrixBlock(2, 2, alpha_coefficient, m_elliptic_Apar_mapped, beta_coefficient, a_Apar_bc);
   

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      alpha_coefficient[dit].setVal(m_ion_skin_depth_sq);
      alpha_coefficient[dit].mult(m_volume_reciprocal[dit]);
   }
   a_preconditioner->constructMatrixBlock(3, 2, alpha_coefficient, m_elliptic_Apar_mapped, a_Apar_bc);

   
   // Construct advection blocks
   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);
   const LevelData<FluxBox>& bunit = m_geometry.getFCBFieldDir();
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      flux[dit].copy(bunit[dit]);
   }
   LevelData<FluxBox> NTV_normal(grids, 1, IntVect::Unit);
   m_geometry.computeMetricTermProductAverage(NTV_normal, flux, false);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      alpha_coefficient[dit].copy(m_volume_reciprocal[dit]);
   }
   
   a_preconditioner->addAdvectionMatrixBlock(0, 3, alpha_coefficient, NTV_normal, *m_jpar_bc);
   a_preconditioner->addAdvectionMatrixBlock(2, 0, alpha_coefficient, NTV_normal, a_phi_bc);

   if (m_include_pol_den_correction) {
   
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         NTV_normal[dit] *= m_electron_temperature_face[dit];
      }

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         alpha_coefficient[dit].divide(m_charge_density_cell[dit]);
         alpha_coefficient[dit].mult(-1.0);
      }
  
      a_preconditioner->addAdvectionMatrixBlock(2, 1, alpha_coefficient, NTV_normal, *vorticity_bc);
   }
   
   a_preconditioner->finalizeMatrix();

}


void
GKEMVorticity::updateESPreconditioner(MBSolver*              a_preconditioner,
                                      MBSolver*              a_pc_Apar,
                                      const EllipticOpBC&    a_phi_bc,
                                      const EllipticOpBC&    a_Apar_bc)
{
   CH_TIME("GKEMVorticity::updateESPreconditioner");

   // Update ES vorticity preconditioner
   GKVorticity::updatePreconditioner(a_preconditioner,
                                     m_volume_reciprocal,
                                     m_mapped_coefficients,
                                     a_phi_bc);

   // Update pc for Apar solve
   const DisjointBoxLayout& grids = m_geometry.grids();
   LevelData<FArrayBox> alpha_coefficient(grids, 1, IntVect::Unit);
   LevelData<FArrayBox> beta_coefficient(grids, 1, IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      beta_coefficient[dit].setVal(m_alpha);
      
      const Box& box = grids[dit];
      FArrayBox alpha_Le2_over_ne(box, 1);

      alpha_Le2_over_ne.setVal(m_alpha);
      alpha_Le2_over_ne.divide(m_charge_density_cell[dit]);
      alpha_Le2_over_ne *= m_electron_skin_depth_sq;
      
#ifndef ALFVEN_WAVE_TEST
      alpha_coefficient[dit].setVal(m_ion_skin_depth_sq);
      alpha_coefficient[dit].divide(m_parallel_cond_cell[dit]);
      alpha_coefficient[dit].plus(alpha_Le2_over_ne);
#else
      alpha_coefficient[dit].copy(alpha_Le2_over_ne);
#endif
      
      alpha_coefficient[dit].mult(m_volume_reciprocal[dit]);
   }
   a_pc_Apar->initializeMatrix();
   a_pc_Apar->constructMatrixBlock(0, 0, alpha_coefficient, m_elliptic_Apar_mapped, beta_coefficient, a_Apar_bc);
   a_pc_Apar->finalizeMatrix();
}


void
GKEMVorticity::updateAlfvenWavePreconditioner(MBSolver*                    a_preconditioner,
                                              const LevelData<FArrayBox>&  a_ion_mass_density,
                                              const EllipticOpBC&          a_phi_bc,
                                              const EllipticOpBC&          a_Apar_bc)
{
   CH_TIME("GKEMVorticity::updateEMPreconditionerExperiment");
   
   const DisjointBoxLayout& grids = m_geometry.grids();
   
   a_preconditioner->initializeMatrix();
   
   LevelData<FluxBox> tensor_coefficient(grids, SpaceDim*SpaceDim, IntVect::Unit);
   LevelData<FArrayBox> alpha_coefficient(grids, 1, IntVect::Unit);
   LevelData<FArrayBox> beta_coefficient(grids, 1, IntVect::Unit);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      tensor_coefficient[dit].setVal(0.);
   }

   if ( !m_eta.isDefined() ) {
      m_eta.define(grids, 1, IntVect::Zero);
   }

   const LevelData<FArrayBox>& B_mag = m_geometry.getCCBFieldMag();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box& box = grids[dit];
       
      FArrayBox polarization_fac(box, 1);
      polarization_fac.setVal(m_larmor_number2);
      polarization_fac.mult(a_ion_mass_density[dit],box, 0, 0);
      polarization_fac.divide(B_mag[dit], box, 0, 0);
      polarization_fac.divide(B_mag[dit], box, 0, 0);
      polarization_fac.divide(m_ion_skin_depth_sq);
      
      m_eta[dit].copy(polarization_fac);
   }

   // Upper left block
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      alpha_coefficient[dit].setVal(0.);
      beta_coefficient[dit].setVal(1.);
      tensor_coefficient[dit].setVal(0.);
   }

   a_preconditioner->constructMatrixBlock(0, 0, alpha_coefficient, tensor_coefficient, beta_coefficient, a_phi_bc);

   // Upper right block
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      alpha_coefficient[dit].copy(m_volume_reciprocal[dit]);
      alpha_coefficient[dit] /= m_eta[dit];
      alpha_coefficient[dit] /= m_alpha;
   }

   const LevelData<FluxBox>& bunit = m_geometry.getFCBFieldDir();

   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      flux[dit].copy(bunit[dit]);
   }
   LevelData<FluxBox> NTV_normal(grids, 1, IntVect::Unit);
   m_geometry.computeMetricTermProductAverage(NTV_normal, flux, false);
   
   a_preconditioner->addAdvectionMatrixBlock(0, 1, alpha_coefficient, NTV_normal, a_Apar_bc);

   // Lower left block

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      alpha_coefficient[dit].setVal(0.);
      beta_coefficient[dit].copy(m_divb[dit]);
      beta_coefficient[dit].negate();
      tensor_coefficient[dit].setVal(0.);
   }

   a_preconditioner->constructMatrixBlock(1, 0, alpha_coefficient, tensor_coefficient, beta_coefficient, a_phi_bc);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      alpha_coefficient[dit].copy(m_volume_reciprocal[dit]);
   }

   a_preconditioner->addAdvectionMatrixBlock(1, 0, alpha_coefficient, NTV_normal, a_phi_bc);

   // Lower right block
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox alpha_Le2_over_ne(grids[dit], 1);
      alpha_Le2_over_ne.setVal(m_alpha);
      alpha_Le2_over_ne.divide(m_charge_density_cell[dit]);
      alpha_Le2_over_ne *= m_electron_skin_depth_sq;

#ifndef ALFVEN_WAVE_TEST
      alpha_coefficient[dit].setVal(m_ion_skin_depth_sq);
      alpha_coefficient[dit].divide(m_parallel_cond_cell[dit]);
      alpha_coefficient[dit].plus(alpha_Le2_over_ne);
#else
      alpha_coefficient[dit].copy(alpha_Le2_over_ne);
#endif
      alpha_coefficient[dit].mult(m_volume_reciprocal[dit]);

      beta_coefficient[dit].setVal(m_alpha);
   }

   a_preconditioner->constructMatrixBlock(1, 1, alpha_coefficient, m_elliptic_Apar_mapped, beta_coefficient, a_Apar_bc);

   a_preconditioner->finalizeMatrix();

   m_pc_updated = true;
}


void
GKEMVorticity::computeDivbFlux(const LevelData<FArrayBox>&  a_in,
                               LevelData<FArrayBox>&        a_out,
                               const bool                   a_fourth_order ) const
{
   const DisjointBoxLayout& grids = m_geometry.grids();

   LevelData<FArrayBox> in_grown( grids, 1, 2*IntVect::Unit );
      
   // This fills two codim 1 ghost cell layers at all block boundaries
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = block_coord_sys.domain();
      
      in_grown[dit].copy(a_in[dit]);

      fourthOrderCellExtrapAtDomainBdry(in_grown[dit], domain, grids[dit]);
   }
   
   // This fills extrablock ghosts and performs an exchange.  We now have 4 layers
   // of ghost cells filled, except at the physical boundary were we only have 2.
   m_geometry.fillInternalGhosts(in_grown);

   // Convert from cell averages to face centers.  We lose 2 ghost cells in this
   // process, so we now only have data on 2 ghost cell layers, except at the
   // physical boundary where we have none.
   // NB: this operation only involves data from the cells along the same direction
   // (i.e., does not involve CODIM 2 isses, and therefore it is safe to use even
   // even for the second-order (since fourthOrderCellExtrapAtDomainBdry() call provided all
   // nesesary information).

   LevelData<FluxBox> in_face( grids, 1, IntVect::Unit );
   if ( a_fourth_order ) {
      fourthOrderCellToFaceCenters(in_face, in_grown);
   }
   else {
      CellToEdge( in_grown, in_face );
   }

   const LevelData<FluxBox>& bunit = m_geometry.getFCBFieldDir();
   LevelData<FluxBox> flux( grids, SpaceDim, IntVect::Unit );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      flux[dit].copy(bunit[dit]);
      for (int dir=0; dir<SpaceDim; ++dir) {
         for (int comp=0; comp<SpaceDim; ++comp) {
            flux[dit][dir].mult(in_face[dit][dir],0,comp,1);
         }
      }
   }

   m_geometry.computeMappedGridDivergence(flux, a_out, a_fourth_order);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_out[dit] *= m_volume_reciprocal[dit];
   }
}


void GKEMVorticity::computeAltParGrad( const LevelData<FArrayBox>&  a_in,
                                       LevelData<FArrayBox>&        a_out ) const
{
   const DisjointBoxLayout& grids = m_geometry.grids();
   LevelData<FArrayBox> tmp(grids,1,IntVect::Zero);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      tmp[dit].copy(a_in[dit]);
      tmp[dit] *= m_divb[dit];
   }

   m_divb_solver->multiplyMatrix(a_in, a_out);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_out[dit] -= tmp[dit];
   }
}


void
GKEMVorticity::solvePreconditioner(const LevelData<FArrayBox>& a_r,
                                   LevelData<FArrayBox>&       a_z )
{
   CH_TIME("GKEMVorticity::solvePreconditioner");

   const DisjointBoxLayout& grids = m_geometry.grids();
   
   if (m_pc_model == "AparSolnOp" || m_pc_model == "AlfvenWaveOriginal") {
      m_preconditioner->solve(a_r, a_z, true);
   }
   else if (m_pc_model == "AlfvenWave") {

      LevelData<FArrayBox> r_phi(grids, 1, IntVect::Zero);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         r_phi[dit].copy(a_r[dit], 0, 0, 1);
      }

      // Compute -L_i^{-2} \Delta_\perp^{-1} r_phi
      LevelData<FArrayBox> neg_lperp_inv_term(grids, 1, IntVect::Zero);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         neg_lperp_inv_term[dit].setVal(0.);  // initial condition
      }
      m_neg_Laplacian_perp->solve(r_phi, neg_lperp_inv_term, true);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         neg_lperp_inv_term[dit] /= m_ion_skin_depth_sq;
      }

      LevelData<FArrayBox> rhs(grids, 2, IntVect::Zero);
      LevelData<FArrayBox> sol(grids, 2, IntVect::Zero);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs[dit].copy(neg_lperp_inv_term[dit],0,0,1);
         rhs[dit].divide(m_eta[dit],0,0,1);
         rhs[dit].divide(m_alpha,0,1);
         rhs[dit].copy(a_r[dit],1,1,1);
         sol[dit].setVal(0.);
      }         

      m_preconditioner->solve(rhs, sol, true);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_z[dit].copy(sol[dit], 0, 0, 2);
      }
   }
   else if (m_pc_model == "EMVorticity") {
   
      // Get 4 component rhs and solution vector
      // contains (phi, vorticity, Apar, jpar)
      LevelData<FArrayBox> r(grids, 4, IntVect::Zero);
      LevelData<FArrayBox> z(grids, 4, IntVect::Zero);
      for (DataIterator dit( r.dataIterator() ); dit.ok(); ++dit) {
         r[dit].setVal(0.);
         z[dit].setVal(0.);
         r[dit].copy(a_r[dit], 0, 0, 1);
         r[dit].copy(a_r[dit], 1, 2, 1);
      }
      // Solve 4x4 block system
      m_preconditioner->solve(r, z, true);
   
      // Return phi and Apar solutions
      for (DataIterator dit( z.dataIterator() ); dit.ok(); ++dit) {
         a_z[dit].copy(z[dit], 0, 0, 1);
         a_z[dit].copy(z[dit], 2, 1, 1);
      }
   }
   
   else if (m_pc_model == "ESVorticity") {

      // Set zero initial guess
      for (DataIterator dit( a_z.dataIterator() ); dit.ok(); ++dit) {
         a_z[dit].setVal(0.);
      }

      // Get phi RHS
      LevelData<FArrayBox> r(grids, 1, IntVect::Zero);
      for (DataIterator dit( a_z.dataIterator() ); dit.ok(); ++dit) {
         r[dit].copy(a_r[dit],0,0,1);
      }

      // Create aliases to the separate phi components of z
      LevelData<FArrayBox> z_phi;
      aliasLevelData(z_phi, &a_z, Interval(0,0));

      // Solve ES preconditioner
      GKVorticity::solvePreconditioner(r, z_phi);

      // Get Apar rhs
      for (DataIterator dit( a_z.dataIterator() ); dit.ok(); ++dit) {
         r[dit].copy(a_r[dit],1,0,1);
      }
      updateAparRHS(r, z_phi);

      // Create aliases to the separate Apar components of z
      LevelData<FArrayBox> z_Apar;
      aliasLevelData(z_Apar, &a_z, Interval(1,1));

      // Find Apar solution
      m_pc_Apar->solve(r, z_Apar, true);
   }
}

void
GKEMVorticity::updateAparRHS(LevelData<FArrayBox>&       a_r,
                             const LevelData<FArrayBox>& a_z_phi)
{
   const DisjointBoxLayout& grids = m_geometry.grids();
   
   LevelData<FArrayBox> neg_grad_phi(grids, SpaceDim, IntVect::Unit);
   computeField(a_z_phi, neg_grad_phi);

   LevelData<FArrayBox> par_neg_grad_phi(grids, 1, IntVect::Zero);
   m_geometry.computeParallelProjection(par_neg_grad_phi, neg_grad_phi);
   
   for (DataIterator dit( a_r.dataIterator() ); dit.ok(); ++dit) {
      a_r[dit].plus(par_neg_grad_phi[dit]);
   }
}

void
GKEMVorticity::applyOp(LevelData<FArrayBox>&       a_out,
                       const LevelData<FArrayBox>& a_in,
                       bool                        a_homogeneous )
{
   // Use base class functionality, because we don't need
   // to subtract FS-averaged parallel divergence
   EllipticOp::computeFluxDivergence(a_in, a_out, a_homogeneous);

   for (DataIterator dit( a_out.dataIterator() ); dit.ok(); ++dit) {
      a_out[dit].mult(m_electron_skin_depth_sq);
      a_out[dit].divide(m_charge_density_cell[dit]);
      a_out[dit].plus(a_in[dit]);
   }
}

void
GKEMVorticity::computeBcDivergence( LevelData<FArrayBox>& a_out )
{

   // Can't we use base class function for this -- no we can't!
   // First, base class calls computeFluxDivergence, not applyOp
   // Second, the case class computeBcDivergence() function seems
   // to call the base class computeFluxDivergence(), and not
   // its overwritten derived class version
   const DisjointBoxLayout& grids = a_out.disjointBoxLayout();
   LevelData<FArrayBox> phi(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi[dit].setVal(0.);
   }

   applyOp(a_out, phi, false);
}


void
GKEMVorticity::parseParameters( const ParmParse&   a_pp )
{

}


void
GKEMVorticity::printParameters()
{
   if (procID()==0) {
      std::cout << "GKEMVorticity pc_model: " << m_pc_model << std::endl;
   }
}

#include "NamespaceFooter.H"
