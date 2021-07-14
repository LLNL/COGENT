#include "GKPoisson.H"
#include "GKPoissonF_F.H"
#include "FourthOrderUtil.H"
#include "FluxSurface.H"
#include "SingleNullCoordSys.H"
#include "Directions.H"
#include "ToroidalBlockLevelExchangeCenter.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "Kernels.H"
#include "MomentOp.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H"
 
const char* GKPoisson::pp_name = {"gkpoisson"};


GKPoisson::GKPoisson( const ParmParse&   a_pp,
                      const MagGeom&     a_geom,
                      const Real         a_larmor_number,
                      const Real         a_debye_number,
                      const bool         a_allocate_precond )
   : EllipticOp(a_pp, a_geom),
     m_moment_op(PS::MomentOp::instance())
{
   init(a_pp, a_larmor_number, a_debye_number, a_allocate_precond);
}
      
GKPoisson::GKPoisson( const ParmParse&                  a_pp,
                      const MagGeom&                    a_geom,
                      const PS::KineticSpeciesPtrVect&  a_kin_species,
                      const Real                        a_larmor_number,
                      const Real                        a_debye_number,
                      const bool                        a_allocate_precond )
   : EllipticOp(a_pp, a_geom),
     m_moment_op(PS::MomentOp::instance())
{
   init(a_pp, a_larmor_number, a_debye_number, a_allocate_precond);
   setFLRObjects(a_kin_species);
}
      


GKPoisson::GKPoisson( const ParmParse&   a_pp,
                      const MagGeom&     a_geom,
                      const Real         a_larmor_number,
                      const Real         a_debye_number,
                      const bool         a_allocate_precond,
                      const bool         a_second_order,
                      const bool         a_low_pollution )
   : EllipticOp(a_pp, a_geom, a_second_order, a_low_pollution),
     m_moment_op(PS::MomentOp::instance())
{
   init(a_pp, a_larmor_number, a_debye_number, a_allocate_precond);
}
      

GKPoisson::~GKPoisson()
{
   if (m_preconditioner) delete m_preconditioner;
   if (m_mblx_ptr) delete m_mblx_ptr;
}


void
GKPoisson::init( const ParmParse&  a_pp,
                 const Real        a_larmor_number,
                 const Real        a_debye_number,
                 const bool        a_allocate_precond )
{
   m_debye_number2 = a_debye_number*a_debye_number;
   m_alpha = 1.;
   m_preconditioner = NULL;
   m_larmor_number2 = a_larmor_number*a_larmor_number;
   m_discretization_order = 2;
   m_mblx_ptr = NULL;
   m_include_FLR_effects = false;
   m_phase_geom = NULL;
   m_model = "GyroPoisson";
   m_verbosity = false;

   m_FLR_bc_buffer = 5;
   m_plot_FLR_switch = false;
   unsetLHSTermType();
   parseParameters( a_pp );

   if (m_verbosity>0) {
      printParameters();
   }

   // We give the coefficients one ghost cell layer so that the
   // second-order centered difference formula can be used to compute
   // the transverse gradients needed for the fourth-order formulas even
   // at box boundaries.
   const DisjointBoxLayout& grids = m_geometry.grids();
   m_mapped_coefficients.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_unmapped_coefficients.define(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   if (m_include_FLR_effects) {
     m_mapped_coeffs_FLR_1.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
     m_mapped_coeffs_FLR_2.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
     m_unmapped_coeffs_FLR_1.define(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
     m_unmapped_coeffs_FLR_2.define(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   }

   if ( m_geometry.shearedMBGeom() ) {
      m_mblx_ptr = new ToroidalBlockLevelExchangeCenter(m_geometry, m_discretization_order, m_discretization_order);
   }
   else {
      if ( m_geometry.extrablockExchange() ) {
         m_mblx_ptr = new MagFluxAlignedMBLEXCenter;
         m_mblx_ptr->define(&m_geometry, m_discretization_order, m_discretization_order);
      }
      else {
         m_mblx_ptr = NULL;
      }
   }

   if ( a_allocate_precond ) {
#ifdef with_petsc
      m_preconditioner = allocatePreconditioner(m_geometry, m_discretization_order, m_mblx_ptr);
#else
      m_preconditioner = allocatePreconditioner(m_geometry, m_discretization_order, m_mblx_ptr);
#endif
      ParmParse pp_precond( ((string)a_pp.prefix() + ".linear_solver.precond").c_str());
      m_preconditioner->setMethodParams(pp_precond);
      m_preconditioner->setConvergenceParams(pp_precond);
   }
}

void 
GKPoisson::setFLRObjects( const PS::KineticSpeciesPtrVect& a_species)
{
  if (m_include_FLR_effects) {

    int n_species(a_species.size());
    CH_assert(n_species > 0);

    m_species_vec = a_species;
    m_phase_geom = &(a_species[0]->phaseSpaceGeometry());
    const PS::LevelData<PS::FArrayBox>&  dfn(a_species[0]->distributionFunction());
    m_phase_grids = &(dfn.disjointBoxLayout());

    m_FLR_integrand_factor.resize(n_species);
    for (int s = 0; s < n_species; s++) {
      m_FLR_integrand_factor[s] 
        = new PS::LevelData<PS::FArrayBox>( *m_phase_grids, 
                                            1, 
                                            PS::IntVect::Zero);
    }

    m_FLR_bc_factor.define( m_geometry.grids(), 1, IntVect::Unit );
    if (m_FLR_switch_func != NULL) {

      m_FLR_switch_func->assign( m_FLR_bc_factor, m_geometry, 0.0 );

    } else {

      const DisjointBoxLayout& dbl(m_FLR_bc_factor.disjointBoxLayout());
      for (DataIterator dit(m_FLR_bc_factor.dataIterator()); dit.ok(); ++dit) {
        m_FLR_bc_factor[dit].setVal(-1.0);
        m_FLR_bc_factor[dit].setVal(1.0, dbl[dit], 0);
      }

      m_FLR_bc_factor.exchange();

      for (DataIterator dit(dbl.dataIterator()); dit.ok(); ++dit) {
  
        const Box& bx = dbl[dit];
        FArrayBox& fab = m_FLR_bc_factor[dit];
        const Box& bx_wg = fab.box();
  
        FArrayBox fab_copy(bx_wg, 1);
        fab_copy.copy(fab, bx_wg);
  
        for (int i = 1; i <= m_FLR_bc_buffer; i++) {
          Real xi = ((Real) i-1.0) / ((Real) m_FLR_bc_buffer);
          for (BoxIterator bit(bx); bit.ok(); ++bit) {
            IntVect iv(bit());
            for (int d = 0; d < SpaceDim; d++) {
              for (int sign = -1; sign <= 1; sign+=2) {
                IntVect iv2 = iv+sign*i*BASISV(d);
                if ( bx_wg.contains(iv2) ) {
                  if (fab_copy(iv2,0) < 0) {
                    //fab(iv,0) *= (0.5 + 0.5*tanh(5*(xi-0.5)));
                    fab(iv,0) *= xi;
                  }
                }
              }
            }
          }
        }
  
      }

    }

    if (m_plot_FLR_switch) {
      m_geometry.plotCellData("FLR_switch_var", m_FLR_bc_factor, 0.0);
    }
  }

  return;
}

   
MBHypreSolver*
GKPoisson::allocatePreconditioner( const MagGeom&                  a_geom,
                                   const int                       a_discretization_order,
                                   MultiBlockLevelExchangeCenter*  a_mblx_ptr )
{
   MBHypreSolver* solver(NULL);
   solver = new MBHypreSolver(a_geom, 1, a_discretization_order, a_mblx_ptr);
   CH_assert(solver != NULL);
   return solver;
}
                                  

void
GKPoisson::updateBoundaries( const EllipticOpBC&  a_bc )
{
   // N.B.: The operator coefficient must be computed prior to calling this function
   CH_TIME("GKPoisson::updateBoundaries");
   setBc(a_bc);
   
   computeBcDivergence( m_bc_divergence );
}


void
GKPoisson::setOperatorCoefficients( const LevelData<FArrayBox>&  a_ion_mass_density,
                                    const EllipticOpBC&          a_bc,
                                    const bool                   a_update_preconditioner )
{
   CH_TIME("GKPoisson::setOperatorCoefficients");

   computeCoefficients(a_ion_mass_density, 
                       m_mapped_coefficients, 
                       m_unmapped_coefficients );

   updateBoundaries(a_bc);

   if ( a_update_preconditioner ) {
      updatePreconditioner( m_preconditioner, 
                            m_volume_reciprocal, 
                            m_mapped_coefficients, 
                            a_bc );
   }
}


void
GKPoisson::updatePreconditioner( MBSolver*              a_preconditioner,
                                 LevelData<FArrayBox>&  a_volume_reciprocal,
                                 LevelData<FluxBox>&    a_mapped_coefficients,
                                 const EllipticOpBC&    a_bc )
{
   a_preconditioner->constructMatrix(a_volume_reciprocal, a_mapped_coefficients, a_bc);
}
                                                                      
                                 
void
GKPoisson::setOperatorCoefficients( const PS::KineticSpeciesPtrVect&  a_kin_species,
                                    const LevelData<FArrayBox>&       a_ion_mass_density,
                                    const EllipticOpBC&               a_bc,
                                    const bool                        a_update_preconditioner )
{
   CH_TIME("GKPoisson::setOperatorCoefficients");

   computeCoefficients(a_ion_mass_density, 
                       m_mapped_coefficients, 
                       m_unmapped_coefficients );

   if (m_include_FLR_effects) {
     m_species_vec = a_kin_species;
     computeFLRIntegrandFactor();
     computeCoefficientsFLR_1(a_ion_mass_density,
                              m_mapped_coeffs_FLR_1, 
                              m_unmapped_coeffs_FLR_1 );
     computeCoefficientsFLR_2(a_ion_mass_density,
                              m_mapped_coeffs_FLR_2, 
                              m_unmapped_coeffs_FLR_2 );
   }

   updateBoundaries(a_bc);

   if ( a_update_preconditioner ) {
      updatePreconditioner( m_preconditioner, 
                            m_volume_reciprocal, 
                            m_mapped_coefficients, 
                            a_bc );
   }
}


void
GKPoisson::setOperatorCoefficients( const LevelData<FArrayBox>& a_ion_mass_density,
                                    const EllipticOpBC&         a_bc,
                                    const bool                  a_update_preconditioner,
                                    double&                     a_lo_value,
                                    double&                     a_hi_value,
                                    LevelData<FArrayBox>&       a_radial_gkp_divergence_average )
{
   CH_TIME("GKPoisson::setOperatorCoefficients2");
   setOperatorCoefficients(a_ion_mass_density, a_bc, a_update_preconditioner);
   
   const DisjointBoxLayout& grids = m_geometry.grids();
   const MagCoordSys& coords = *m_geometry.getCoordSys();
   
   LevelData<FluxBox> ones( m_geometry.grids(), 1, IntVect::Zero);
   for (DataIterator dit( ones.dataIterator() ); dit.ok(); ++dit) {
      int block_number = coords.whichBlock(grids[dit]);
      
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_fab = ones[dit][dir];
         this_fab.setVal(0.);
         
         if ( dir == RADIAL_DIR ) {

	   if ((typeid(coords) != typeid(SingleNullCoordSys)) ||
               ((const SingleNullCoordSys&)coords).isCORE(block_number))  {
	     
	     ones[dit].setVal(1.,dir);
	   }
         }
      }
   }
   
   //Compute <c^2*n*m/(e*B^2)*|grad(Psi)|^2>
   computeRadialFSAverage(ones, a_lo_value, a_hi_value, a_radial_gkp_divergence_average);
}


void
GKPoisson::computeCoefficients(const LevelData<FArrayBox>& a_ion_mass_density,
                               LevelData<FluxBox>&         a_mapped_coefficients,
                               LevelData<FluxBox>&         a_unmapped_coefficients )
{
   CH_TIME("GKPoisson::computeCoefficients");

   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   LevelData<FArrayBox> density_sum_cell( grids, 1, 2*IntVect::Unit );
      
   // This fills two codim 1 ghost cell layers at all block boundaries
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = block_coord_sys.domain();
      
      // First component of a_ion_mass_density contains this quantity summed 
      // over *all* species.
      density_sum_cell[dit].copy(a_ion_mass_density[dit], 0, 0);

      fourthOrderCellExtrapAtDomainBdry(density_sum_cell[dit], domain, grids[dit]);
   }
   
   // This fills extrablock ghosts and performs an exchange.  We now have 4 layers
   // of ghost cells filled, except at the physical boundary were we only have 2.
   m_geometry.fillInternalGhosts( density_sum_cell );
   
   // Convert from cell averages to face centers.  We lose 2 ghost cells in this
   // process, so we now only have data on 2 ghost cell layers, except at the
   // physical boundary where we have none.
   // NB: this operation only involves data from the cells along the same direction
   // (i.e., does not involve CODIM 2 isses, and therefore it is safe to use even
   // even for the second-order (since fourthOrderCellExtrapAtDomainBdry() call provided all
   // nesesary information).

   LevelData<FluxBox> density_sum_face( grids, 1, IntVect::Zero );
   fourthOrderCellToFaceCenters(density_sum_face, density_sum_cell);

   const IntVect grown_ghosts( m_mapped_coefficients.ghostVect() + IntVect::Unit );
   CH_assert(grown_ghosts == 2*IntVect::Unit);
   const LevelData<FluxBox>& BFieldMag = m_geometry.getFCBFieldMag();
   CH_assert(BFieldMag.ghostVect()>=grown_ghosts);
   const LevelData<FluxBox>& BFieldDir = m_geometry.getFCBFieldDir();
   CH_assert(BFieldDir.ghostVect()>=grown_ghosts);
   
   LevelData<FluxBox> grown_mapped_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);
   LevelData<FluxBox> grown_unmapped_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);

   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& par_coeff = m_geometry.getEllipticOpParCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();

   LevelData<FluxBox> radial_coeff( grids, SpaceDim*SpaceDim, 2*IntVect::Unit );
   LevelData<FluxBox> radial_coeff_mapped( grids, SpaceDim*SpaceDim, 2*IntVect::Unit );
   if (m_model == "RadialGyroPoisson") {
      m_geometry.getEllipticOpRadCoeff(radial_coeff);
      m_geometry.getEllipticOpRadCoeffMapped(radial_coeff_mapped);
   }
  
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {

      const Box& box = grids[dit];

      FluxBox isotropic_coeff(box, SpaceDim * SpaceDim);
      isotropic_coeff.copy(perp_coeff[dit]);
      isotropic_coeff += par_coeff[dit];
      isotropic_coeff *= m_debye_number2;

      FluxBox isotropic_coeff_mapped(box, SpaceDim * SpaceDim);
      isotropic_coeff_mapped.copy(perp_coeff_mapped[dit]);
      isotropic_coeff_mapped += par_coeff_mapped[dit];
      isotropic_coeff_mapped *= m_debye_number2;
      
      FluxBox polarization_fac(box, 1);
      polarization_fac.setVal(1.0);
      polarization_fac *= m_larmor_number2;
      polarization_fac.mult(density_sum_face[dit],box, 0, 0);
      polarization_fac.divide(BFieldMag[dit], box, 0, 0);
      polarization_fac.divide(BFieldMag[dit], box, 0, 0);
      
      FluxBox tmp_perp(box, SpaceDim * SpaceDim);
      tmp_perp.copy(perp_coeff[dit]);
      
      FluxBox tmp_par(box, SpaceDim * SpaceDim);
      tmp_par.copy(par_coeff[dit]);
      
      FluxBox tmp_perp_mapped(box, SpaceDim * SpaceDim);
      tmp_perp_mapped.copy(perp_coeff_mapped[dit]);
      
      FluxBox tmp_par_mapped(box, SpaceDim * SpaceDim);
      tmp_par_mapped.copy(par_coeff_mapped[dit]);
      
      if (m_model == "GyroPoisson") {
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_perp.mult(polarization_fac, box, 0, n);
            tmp_perp_mapped.mult(polarization_fac, box, 0, n);
         }
         grown_unmapped_coefficients[dit].copy(isotropic_coeff);
         grown_unmapped_coefficients[dit] += tmp_perp;

         grown_mapped_coefficients[dit].copy(isotropic_coeff_mapped);
         grown_mapped_coefficients[dit] += tmp_perp_mapped;

      }

      else if (m_model == "PerpGyroPoisson") {
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_perp.mult(polarization_fac, box, 0, n);
            tmp_perp_mapped.mult(polarization_fac, box, 0, n);
         }
         grown_unmapped_coefficients[dit].copy(tmp_perp);
         grown_mapped_coefficients[dit].copy(tmp_perp_mapped);

      }

      else if (m_model == "RadialGyroPoisson") {

         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            radial_coeff[dit].mult(polarization_fac, box, 0, n);
            radial_coeff_mapped[dit].mult(polarization_fac, box, 0, n);
         }
         grown_unmapped_coefficients[dit].copy(radial_coeff[dit]);
         grown_mapped_coefficients[dit].copy(radial_coeff_mapped[dit]);
      }
      

      else {
         MayDay::Error("GKPoisson:: unknown model is specified");
      }
   }
   
   if (!m_second_order) {

      grown_unmapped_coefficients.exchange();
      grown_mapped_coefficients.exchange();
      
      // The mapped coefficients must now be converted to face averages, which
      // requires a layer of transverse ghost faces that we don't have at the
      // physical boundary, so we need to extrapolate them.
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);
   
      // Convert the mapped coefficients from face-centered to face-averaged
      fourthOrderAverage(grown_mapped_coefficients);
   
      m_geometry.fillTransverseGhosts(grown_unmapped_coefficients, false);
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);

   }
   
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_mapped_coefficients[dit].copy(grown_mapped_coefficients[dit]);
      a_unmapped_coefficients[dit].copy(grown_unmapped_coefficients[dit]);
   }
}

void
GKPoisson::computeCoefficientsFLR_1(const LevelData<FArrayBox>&  a_ion_mass_density,
                                    LevelData<FluxBox>&          a_mapped_coefficients,
                                    LevelData<FluxBox>&          a_unmapped_coefficients )
{
   CH_TIME("GKPoisson::computeCoefficientsFLR_1");

   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   LevelData<FArrayBox> density_sum_cell( grids, 1, 2*IntVect::Unit );
      
   // This fills two codim 1 ghost cell layers at all block boundaries
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = block_coord_sys.domain();
    
      // Second component of a_ion_mass_density contains this quantity
      // summed over all non-FLR species.
      density_sum_cell[dit].copy(a_ion_mass_density[dit], 1, 0);
      fourthOrderCellExtrapAtDomainBdry(density_sum_cell[dit], domain, grids[dit]);
   }
   
   // This fills extrablock ghosts and performs an exchange.  We now have 4 layers
   // of ghost cells filled, except at the physical boundary were we only have 2.
   m_geometry.fillInternalGhosts( density_sum_cell );
   
   // Convert from cell averages to face centers.  We lose 2 ghost cells in this
   // process, so we now only have data on 2 ghost cell layers, except at the
   // physical boundary where we have none.
   // NB: this operation only involves data from the cells along the same direction
   // (i.e., does not involve CODIM 2 isses, and therefore it is safe to use even
   // even for the second-order (since fourthOrderCellExtrapAtDomainBdry() call provided all
   // nesesary information).

   LevelData<FluxBox> density_sum_face( grids, 1, IntVect::Zero );
   fourthOrderCellToFaceCenters(density_sum_face, density_sum_cell);

   const IntVect grown_ghosts( m_mapped_coefficients.ghostVect() + IntVect::Unit );
   CH_assert(grown_ghosts == 2*IntVect::Unit);
   const LevelData<FluxBox>& BFieldMag = m_geometry.getFCBFieldMag();
   CH_assert(BFieldMag.ghostVect()>=grown_ghosts);
   const LevelData<FluxBox>& BFieldDir = m_geometry.getFCBFieldDir();
   CH_assert(BFieldDir.ghostVect()>=grown_ghosts);
   
   LevelData<FluxBox> grown_mapped_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);
   LevelData<FluxBox> grown_unmapped_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);

   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& par_coeff = m_geometry.getEllipticOpParCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();

   LevelData<FluxBox> radial_coeff( grids, SpaceDim*SpaceDim, 2*IntVect::Unit );
   LevelData<FluxBox> radial_coeff_mapped( grids, SpaceDim*SpaceDim, 2*IntVect::Unit );
   if (m_model == "RadialGyroPoisson") {
      m_geometry.getEllipticOpRadCoeff(radial_coeff);
      m_geometry.getEllipticOpRadCoeffMapped(radial_coeff_mapped);
   }
  
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {

      const Box& box = grids[dit];

      FluxBox isotropic_coeff(box, SpaceDim * SpaceDim);
      isotropic_coeff.copy(perp_coeff[dit]);
      isotropic_coeff += par_coeff[dit];
      isotropic_coeff *= m_debye_number2;

      FluxBox isotropic_coeff_mapped(box, SpaceDim * SpaceDim);
      isotropic_coeff_mapped.copy(perp_coeff_mapped[dit]);
      isotropic_coeff_mapped += par_coeff_mapped[dit];
      isotropic_coeff_mapped *= m_debye_number2;
      
      FluxBox polarization_fac(box, 1);
      polarization_fac.setVal(1.0);
      polarization_fac *= m_larmor_number2;
      polarization_fac.mult(density_sum_face[dit],box, 0, 0);
      polarization_fac.divide(BFieldMag[dit], box, 0, 0);
      polarization_fac.divide(BFieldMag[dit], box, 0, 0);
      
      FluxBox tmp_perp(box, SpaceDim * SpaceDim);
      tmp_perp.copy(perp_coeff[dit]);
      
      FluxBox tmp_par(box, SpaceDim * SpaceDim);
      tmp_par.copy(par_coeff[dit]);
      
      FluxBox tmp_perp_mapped(box, SpaceDim * SpaceDim);
      tmp_perp_mapped.copy(perp_coeff_mapped[dit]);
      
      FluxBox tmp_par_mapped(box, SpaceDim * SpaceDim);
      tmp_par_mapped.copy(par_coeff_mapped[dit]);
      
      if (m_model == "GyroPoisson") {
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_perp.mult(polarization_fac, box, 0, n);
            tmp_perp_mapped.mult(polarization_fac, box, 0, n);
         }
         grown_unmapped_coefficients[dit].copy(isotropic_coeff);
         grown_unmapped_coefficients[dit] += tmp_perp;

         grown_mapped_coefficients[dit].copy(isotropic_coeff_mapped);
         grown_mapped_coefficients[dit] += tmp_perp_mapped;

      }

      else if (m_model == "PerpGyroPoisson") {
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_perp.mult(polarization_fac, box, 0, n);
            tmp_perp_mapped.mult(polarization_fac, box, 0, n);
         }
         grown_unmapped_coefficients[dit].copy(tmp_perp);
         grown_mapped_coefficients[dit].copy(tmp_perp_mapped);

      }

      else if (m_model == "RadialGyroPoisson") {

         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            radial_coeff[dit].mult(polarization_fac, box, 0, n);
            radial_coeff_mapped[dit].mult(polarization_fac, box, 0, n);
         }
         grown_unmapped_coefficients[dit].copy(radial_coeff[dit]);
         grown_mapped_coefficients[dit].copy(radial_coeff_mapped[dit]);
      }
      

      else {
         MayDay::Error("GKPoisson:: unknown model is specified");
      }
   }
   
   if (!m_second_order) {

      grown_unmapped_coefficients.exchange();
      grown_mapped_coefficients.exchange();
      
      // The mapped coefficients must now be converted to face averages, which
      // requires a layer of transverse ghost faces that we don't have at the
      // physical boundary, so we need to extrapolate them.
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);
   
      // Convert the mapped coefficients from face-centered to face-averaged
      fourthOrderAverage(grown_mapped_coefficients);
   
      m_geometry.fillTransverseGhosts(grown_unmapped_coefficients, false);
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);

   }
   
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_mapped_coefficients[dit].copy(grown_mapped_coefficients[dit]);
      a_unmapped_coefficients[dit].copy(grown_unmapped_coefficients[dit]);
   }
}

void
GKPoisson::computeCoefficientsFLR_2(const LevelData<FArrayBox>&  a_ion_mass_density,
                                    LevelData<FluxBox>&          a_mapped_coefficients,
                                    LevelData<FluxBox>&          a_unmapped_coefficients )
{
   CH_TIME("GKPoisson::computeCoefficientsFLR_2");

   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   LevelData<FArrayBox> density_sum_cell( grids, 1, 2*IntVect::Unit );
      
   // This fills two codim 1 ghost cell layers at all block boundaries
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = block_coord_sys.domain();
    
      // First component of a_ion_mass_density contains this quantity
      // summed over *all* species.
      // Second component of a_ion_mass_density contains this quantity
      // summed over all non-FLR species.
      density_sum_cell[dit].copy(a_ion_mass_density[dit], 0, 0);
      density_sum_cell[dit].minus(a_ion_mass_density[dit], 1, 0);
      // Now we have the a_ion_mass_density summed over FLR species
      fourthOrderCellExtrapAtDomainBdry(density_sum_cell[dit], domain, grids[dit]);
   }
   
   // This fills extrablock ghosts and performs an exchange.  We now have 4 layers
   // of ghost cells filled, except at the physical boundary were we only have 2.
   m_geometry.fillInternalGhosts( density_sum_cell );
   
   // Convert from cell averages to face centers.  We lose 2 ghost cells in this
   // process, so we now only have data on 2 ghost cell layers, except at the
   // physical boundary where we have none.
   // NB: this operation only involves data from the cells along the same direction
   // (i.e., does not involve CODIM 2 isses, and therefore it is safe to use even
   // even for the second-order (since fourthOrderCellExtrapAtDomainBdry() call provided all
   // nesesary information).

   LevelData<FluxBox> density_sum_face( grids, 1, IntVect::Zero );
   fourthOrderCellToFaceCenters(density_sum_face, density_sum_cell);

   const IntVect grown_ghosts( m_mapped_coefficients.ghostVect() + IntVect::Unit );
   CH_assert(grown_ghosts == 2*IntVect::Unit);
   const LevelData<FluxBox>& BFieldMag = m_geometry.getFCBFieldMag();
   CH_assert(BFieldMag.ghostVect()>=grown_ghosts);
   const LevelData<FluxBox>& BFieldDir = m_geometry.getFCBFieldDir();
   CH_assert(BFieldDir.ghostVect()>=grown_ghosts);
   
   LevelData<FluxBox> grown_mapped_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);
   LevelData<FluxBox> grown_unmapped_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);

   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& par_coeff = m_geometry.getEllipticOpParCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();

   LevelData<FluxBox> radial_coeff( grids, SpaceDim*SpaceDim, 2*IntVect::Unit );
   LevelData<FluxBox> radial_coeff_mapped( grids, SpaceDim*SpaceDim, 2*IntVect::Unit );
   if (m_model == "RadialGyroPoisson") {
      m_geometry.getEllipticOpRadCoeff(radial_coeff);
      m_geometry.getEllipticOpRadCoeffMapped(radial_coeff_mapped);
   }
  
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {

      const Box& box = grids[dit];

      FluxBox polarization_fac(box, 1);
      polarization_fac.setVal(1.0);
      polarization_fac *= m_larmor_number2;
      polarization_fac.mult(density_sum_face[dit],box, 0, 0);
      polarization_fac.divide(BFieldMag[dit], box, 0, 0);
      polarization_fac.divide(BFieldMag[dit], box, 0, 0);
      
      FluxBox tmp_perp(box, SpaceDim * SpaceDim);
      tmp_perp.copy(perp_coeff[dit]);
      
      FluxBox tmp_par(box, SpaceDim * SpaceDim);
      tmp_par.copy(par_coeff[dit]);
      
      FluxBox tmp_perp_mapped(box, SpaceDim * SpaceDim);
      tmp_perp_mapped.copy(perp_coeff_mapped[dit]);
      
      FluxBox tmp_par_mapped(box, SpaceDim * SpaceDim);
      tmp_par_mapped.copy(par_coeff_mapped[dit]);
      
      if ( (m_model == "GyroPoisson") || (m_model == "PerpGyroPoisson") ){
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_perp.mult(polarization_fac, box, 0, n);
            tmp_perp_mapped.mult(polarization_fac, box, 0, n);
         }
         grown_unmapped_coefficients[dit].copy(tmp_perp);
         grown_mapped_coefficients[dit].copy(tmp_perp_mapped);

      }

      else if (m_model == "RadialGyroPoisson") {

         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            radial_coeff[dit].mult(polarization_fac, box, 0, n);
            radial_coeff_mapped[dit].mult(polarization_fac, box, 0, n);
         }
         grown_unmapped_coefficients[dit].copy(radial_coeff[dit]);
         grown_mapped_coefficients[dit].copy(radial_coeff_mapped[dit]);
      }

      else {
         MayDay::Error("GKPoisson:: unknown model is specified");
      }
   }
   
   if (!m_second_order) {

      grown_unmapped_coefficients.exchange();
      grown_mapped_coefficients.exchange();
      
      // The mapped coefficients must now be converted to face averages, which
      // requires a layer of transverse ghost faces that we don't have at the
      // physical boundary, so we need to extrapolate them.
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);
   
      // Convert the mapped coefficients from face-centered to face-averaged
      fourthOrderAverage(grown_mapped_coefficients);
   
      m_geometry.fillTransverseGhosts(grown_unmapped_coefficients, false);
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);

   }
   
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_mapped_coefficients[dit].copy(grown_mapped_coefficients[dit]);
      a_unmapped_coefficients[dit].copy(grown_unmapped_coefficients[dit]);
   }
}

void
GKPoisson::getMinMax(LevelData<FArrayBox>& a_data, double& a_min, double& a_max) const
{
   const DisjointBoxLayout & grids = a_data.disjointBoxLayout();

   a_min =  DBL_MAX;
   a_max = -DBL_MAX;
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      FArrayBox& this_data = a_data[dit];
      Box box(grids[dit]);
      if (this_data.min(box) < a_min) a_min = this_data.min(box);
      if (this_data.max(box) > a_max) a_max = this_data.max(box);
   }   

#ifdef CH_MPI
   double local_min = a_min;
   MPI_Allreduce(&local_min, &a_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   double local_max = a_max;
   MPI_Allreduce(&local_max, &a_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}



void
GKPoisson::getMinMax(LevelData<FluxBox>& a_data, double& a_min, double& a_max) const
{
   const DisjointBoxLayout & grids = a_data.disjointBoxLayout();

   a_min =  DBL_MAX;
   a_max = -DBL_MAX;
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         Box box(surroundingNodes(grids[dit],dir));
         FArrayBox& this_data = a_data[dit][dir];
         if (this_data.min(box) < a_min) a_min = this_data.min(box);
         if (this_data.max(box) > a_max) a_max = this_data.max(box);
      }
   }   

#ifdef CH_MPI
   double local_min = a_min;
   MPI_Allreduce(&local_min, &a_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   double local_max = a_max;
   MPI_Allreduce(&local_max, &a_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}

//void
//GKPoisson::setModel(const std::string& a_model)
//{
//   m_model = a_model;
//}



void
GKPoisson::solvePreconditioner( const LevelData<FArrayBox>& a_r,
                                LevelData<FArrayBox>&       a_z )
{
   m_preconditioner->solve(a_r, a_z, true);
}



void
GKPoisson::multiplyCoefficients( LevelData<FluxBox>& a_data,
                                 const bool a_mapped_coeff ) const
{
   CH_assert(a_data.ghostVect() <= m_unmapped_coefficients.ghostVect());

   if (m_include_FLR_effects) {

     if (lhsTermType() == gkp_ct_flr1) {

       for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
          FluxBox& this_data = a_data[dit];
          FluxBox saved_data(this_data.box(),SpaceDim);
          saved_data.copy(this_data);
          this_data.setVal(0.);
          for (int dir=0; dir<SpaceDim; ++dir) { 
             FArrayBox& this_data_dir = this_data[dir];
             FArrayBox& this_saved_data_dir = saved_data[dir];
             FArrayBox tmp(this_data_dir.box(),1);
             for (int i=0; i<SpaceDim; ++i) {
                for (int j=0; j<SpaceDim; ++j) {
                   if (a_mapped_coeff) {
                      tmp.copy( m_mapped_coeffs_FLR_1[dit][dir], 
                                SpaceDim*i+j, 
                                0, 
                                1 );
                   }
                   else {
                      tmp.copy( m_unmapped_coeffs_FLR_1[dit][dir], 
                                SpaceDim*i+j, 
                                0, 
                                1 );
                   }
                   tmp.mult(this_saved_data_dir, j, 0, 1);
                   this_data_dir.plus(tmp, 0, i, 1);
                }
             }
          }
       }

     } else if (lhsTermType() == gkp_ct_flr2) {

       for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
          FluxBox& this_data = a_data[dit];
          FluxBox saved_data(this_data.box(),SpaceDim);
          saved_data.copy(this_data);
          this_data.setVal(0.);
          for (int dir=0; dir<SpaceDim; ++dir) { 
             FArrayBox& this_data_dir = this_data[dir];
             FArrayBox& this_saved_data_dir = saved_data[dir];
             FArrayBox tmp(this_data_dir.box(),1);
             for (int i=0; i<SpaceDim; ++i) {
                for (int j=0; j<SpaceDim; ++j) {
                   if (a_mapped_coeff) {
                      tmp.copy( m_mapped_coeffs_FLR_2[dit][dir], 
                                SpaceDim*i+j,
                                0, 
                                1 );
                   }
                   else {
                      tmp.copy( m_unmapped_coeffs_FLR_2[dit][dir], 
                                SpaceDim*i+j, 
                                0, 
                                1 );
                   }
                   tmp.mult(this_saved_data_dir, j, 0, 1);
                   this_data_dir.plus(tmp, 0, i, 1);
                }
             }
          }
       }

     } else if (lhsTermType() == gkp_ct_no_flr) {

       for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
          FluxBox& this_data = a_data[dit];
          FluxBox saved_data(this_data.box(),SpaceDim);
          saved_data.copy(this_data);
          this_data.setVal(0.);
          for (int dir=0; dir<SpaceDim; ++dir) { 
             FArrayBox& this_data_dir = this_data[dir];
             FArrayBox& this_saved_data_dir = saved_data[dir];
             FArrayBox tmp(this_data_dir.box(),1);
             for (int i=0; i<SpaceDim; ++i) {
                for (int j=0; j<SpaceDim; ++j) {
                   if (a_mapped_coeff) {
                      tmp.copy( m_mapped_coefficients[dit][dir], 
                                SpaceDim*i+j,
                                0, 
                                1 );
                   }
                   else {
                      tmp.copy( m_unmapped_coefficients[dit][dir], 
                                SpaceDim*i+j, 
                                0, 
                                1 );
                   }
                   tmp.mult(this_saved_data_dir, j, 0, 1);
                   this_data_dir.plus(tmp, 0, i, 1);
                }
             }
          }
       }

     }

   } else {

      for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
         FluxBox& this_data = a_data[dit];
         FluxBox saved_data(this_data.box(),SpaceDim);
         saved_data.copy(this_data);
         this_data.setVal(0.);
         for (int dir=0; dir<SpaceDim; ++dir) { 
            FArrayBox& this_data_dir = this_data[dir];
            FArrayBox& this_saved_data_dir = saved_data[dir];
            FArrayBox tmp(this_data_dir.box(),1);
            for (int i=0; i<SpaceDim; ++i) {
               for (int j=0; j<SpaceDim; ++j) {
                  if (a_mapped_coeff) {
                     tmp.copy( m_mapped_coefficients[dit][dir], 
                               SpaceDim*i+j,
                               0, 
                               1 );
                  }
                  else {
                     tmp.copy( m_unmapped_coefficients[dit][dir], 
                               SpaceDim*i+j, 
                               0, 
                               1 );
                  }
                  tmp.mult(this_saved_data_dir, j, 0, 1);
                  this_data_dir.plus(tmp, 0, i, 1);
               }
            }
         }
      }

   }

   a_data.exchange();
}

void
GKPoisson::parseParameters( const ParmParse&   a_pp )
{
   if (a_pp.contains("model")) {
      a_pp.get( "model", m_model );
   }

   a_pp.query( "include_FLR_effects", m_include_FLR_effects );
   if (m_include_FLR_effects) {
     a_pp.query("FLR_plot_switch_func", m_plot_FLR_switch);
     if (a_pp.contains("FLR_switch_func")) {
       GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
       std::string func_name;
       a_pp.get("FLR_switch_func", func_name);
       m_FLR_switch_func = grid_library->find(func_name);
     } else {
       a_pp.query( "FLR_boundary_buffer_size", m_FLR_bc_buffer);
     }
   }

   ParmParse pp_precond( ((string)a_pp.prefix() + ".precond").c_str());
   if ( pp_precond.query("precond_method", m_precond_method) == 0 ) m_precond_method = "AMG";
   if ( pp_precond.query("precond_tol", m_precond_tol) == 0 ) m_precond_tol = 0.;
   if ( pp_precond.query("precond_max_iter", m_precond_max_iter) == 0 ) m_precond_max_iter = 1;
   if ( pp_precond.query("precond_verbose", m_precond_verbose) == 0 ) m_precond_verbose = false;

   ParmParse pp_precond_precond( ((string)pp_precond.prefix() + ".precond").c_str());
   if ( pp_precond_precond.query("precond_precond_method", m_precond_precond_method) == 0 ) m_precond_precond_method = "";
   if ( pp_precond_precond.query("precond_precond_tol", m_precond_precond_tol) == 0 ) m_precond_precond_tol = 0.;
   if ( pp_precond_precond.query("precond_precond_max_iter", m_precond_precond_max_iter) == 0 ) m_precond_precond_max_iter = 0;
   if ( pp_precond_precond.query("precond_precond_verbose", m_precond_precond_verbose) == 0 ) m_precond_precond_verbose = false;
}


void
GKPoisson::printParameters()
{
   if (procID()==0) {
      std::cout << "GKPoisson model: " << m_model << std::endl;
      std::cout << "GKPoisson include FLR effects: " 
                << (m_include_FLR_effects ? "true" : "false") << std::endl;
      if (m_include_FLR_effects) {
        if (m_FLR_switch_func != NULL) {
        } else {
          std::cout << "    FLR boundary buffer size: " << m_FLR_bc_buffer 
                    << "\n";
        }
      }
   }
}

void
GKPoisson::applyOp( LevelData<FArrayBox>&       a_out,
                    const LevelData<FArrayBox>& a_in,
                    bool                        a_homogeneous )
{
  if (!m_include_FLR_effects) {
   
    setLHSTermType(gkp_ct_no_flr);
    EllipticOp::applyOp(a_out, a_in, a_homogeneous);
    unsetLHSTermType();

  } else {

    setLHSTermType(gkp_ct_flr1);
    EllipticOp::applyOp(a_out, a_in, a_homogeneous);
    unsetLHSTermType();

    LevelData<FArrayBox> term2( a_out.disjointBoxLayout(), 
                                a_out.nComp(), 
                                a_out.ghostVect() );
    for (DataIterator dit(a_out.dataIterator()); dit.ok(); ++dit) {
      term2[dit].setVal(0.0);
    }
    setLHSTermType(gkp_ct_flr2);
    EllipticOp::applyOp(term2, a_in, a_homogeneous);
    unsetLHSTermType();

    for (DataIterator dit(a_out.dataIterator()); dit.ok(); ++dit) {
      FArrayBox factor(term2[dit].box(), 1);
      factor.setVal(1.0);
      factor.minus(m_FLR_bc_factor[dit]);

      term2[dit].mult(factor);
      a_out[dit] += term2[dit];
    }

    /* inject phi to phase space */
    PS::LevelData<PS::FArrayBox> phi;
    m_phase_geom->injectConfigurationToPhase(a_in, phi);

    const LevelData<FArrayBox>& BFieldMag = m_geometry.getCCBFieldMag();

    for (int s = 0; s < m_species_vec.size(); s++) {

      const PS::KineticSpecies& this_species = *(m_species_vec[s]);
      if (!this_species.isGyrokinetic()) continue;
      const PS::GyroaverageOperator* gyroavg_op = this_species.gyroaverageOp();
      const PS::LevelData<PS::FArrayBox>& factor = *(m_FLR_integrand_factor[s]);
   
      /* gyroaverage this phi twice */
      PS::LevelData<PS::FArrayBox> phi_tilde;
      gyroavg_op->applyOpTwice(phi_tilde, phi);

      PS::LevelData<PS::FArrayBox> integrand(*m_phase_grids, 1, PS::IntVect::Zero);
      for (PS::DataIterator dit(integrand.dataIterator()); dit.ok(); ++dit) {
        PS::FArrayBox& integrand_fab = integrand[dit];
        const PS::FArrayBox& phi_fab = phi[dit];
        const PS::FArrayBox& phi_tilde_fab = phi_tilde[dit];
        const PS::FArrayBox& factor_fab = factor[dit];
        const PS::Box& box = (*m_phase_grids)[dit];

        const PS::Box& phi_bx(phi_fab.box());
        const PS::Box& phi_tilde_bx(phi_tilde_fab.box());

        const PS::IntVect smallend_phi = phi_bx.smallEnd();
        const PS::IntVect smallend_phi_tilde = phi_tilde_bx.smallEnd();

        for (PS::BoxIterator bit(box); bit.ok(); ++bit) {
          PS::IntVect iv_phase(bit());

          PS::IntVect iv_phi(bit());
          iv_phi[VPARALLEL_DIR] = smallend_phi[VPARALLEL_DIR];
          iv_phi[MU_DIR] = smallend_phi[MU_DIR];

          PS::IntVect iv_phi_tilde(bit());
          iv_phi_tilde[VPARALLEL_DIR] = smallend_phi_tilde[VPARALLEL_DIR];

          /* multiplying by the mass because the moment op will
           * divide by the mass (see MomentOp::compute() */
          integrand_fab(iv_phase,0) =   this_species.mass() 
                                      * factor_fab(iv_phase,0)
                                      * (   phi_fab(iv_phi,0)
                                          - phi_tilde_fab(iv_phi_tilde,0) );
        }

      }

      /* compute polarization term and add it to a_out */
      LevelData<FArrayBox> n_pol_species( a_out.disjointBoxLayout(),
                                          1,
                                          IntVect::Zero );
      m_moment_op.compute(  n_pol_species, 
                            this_species, 
                            integrand, 
                            PS::DensityKernel<PS::FArrayBox>() );

      for (DataIterator dit(a_out.dataIterator()); dit.ok(); ++dit) {

        const Box& box = n_pol_species[dit].box();

        FArrayBox polarization_fac(box, 1);
        polarization_fac.setVal(-1.0);
        polarization_fac *= (2.0*this_species.charge()*this_species.charge()/this_species.mass());
        polarization_fac.divide(BFieldMag[dit], box, 0, 0);

        // multiply the near-boundary switching factor
        polarization_fac.mult(m_FLR_bc_factor[dit], box, 0, 0);

        n_pol_species[dit].mult(polarization_fac);
        a_out[dit].plus(n_pol_species[dit]);

      }

    }

  }

  return;
}

void
GKPoisson::computeFLRIntegrandFactor()
{
  CH_assert(m_include_FLR_effects);

  for (int s = 0; s < m_species_vec.size(); s++) {

    const PS::KineticSpecies& species = *(m_species_vec[s]);
    const PS::LevelData<PS::FArrayBox>& dfn = species.distributionFunction();
    PS::LevelData<PS::FArrayBox>& factor = *(m_FLR_integrand_factor[s]);
    PS::LevelData<PS::FArrayBox> grad_dfn( dfn.disjointBoxLayout(), 
                                           PDIM, PS::IntVect::Zero);
    m_phase_geom->computeGradient( dfn, grad_dfn, (m_second_order ? 2 : 4));

    for (PS::DataIterator dit(factor.dataIterator()); dit.ok(); ++dit) {

      PS::FArrayBox& factor_fab = factor[dit];
      const PS::FArrayBox& grad_dfn_fab = grad_dfn[dit];

      const PS::Box& box = (*m_phase_grids)[dit];
      for (PS::BoxIterator bit(box); bit.ok(); ++bit) {
        factor_fab(bit(),0) = grad_dfn_fab(bit(), MU_DIR);
      }
    }
  }

  return;
}

#include "NamespaceFooter.H"
