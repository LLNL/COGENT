#include "GKPoisson.H"
#include "GKPoissonF_F.H"
#include "FourthOrderUtil.H"
#include "FluxSurface.H"
#include "Directions.H"
#include "ToroidalBlockLevelExchangeCenter.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "Kernels.H"
#include "MomentOp.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "inspect.H"

#include "NamespaceHeader.H"
 
const char* GKPoisson::pp_name = {"gkpoisson"};


GKPoisson::GKPoisson( const ParmParse&   a_pp,
                      const MagGeom&     a_geom,
                      const Real         a_larmor_number,
                      const Real         a_debye_number )
   : EllipticOp(a_pp, a_geom),
     m_debye_number2(a_debye_number*a_debye_number),
     m_alpha(1.),
     m_model("GyroPoisson"),
     m_larmor_number2(a_larmor_number*a_larmor_number),
     m_mblx_ptr(NULL),
     m_include_FLR_effects(false),
     m_phase_geom(NULL),
     m_moment_op( PS::MomentOp::instance() )
{

   parseParameters( a_pp );

   if (m_verbosity>0) {
      printParameters();
   }

   // We give the coefficients one ghost cell layer so that the
   // second-order centered difference formula can be used to compute
   // the transverse gradients needed for the fourth-order formulas even
   // at box boundaries.
   const DisjointBoxLayout& grids = a_geom.grids();
   m_mapped_coefficients.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_unmapped_coefficients.define(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   if (m_include_FLR_effects) {
     m_mapped_coefficients_gk.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
     m_unmapped_coefficients_gk.define(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   }

   int discretization_order = 2;

   if ( a_geom.shearedMBGeom() ) {
      m_mblx_ptr = new ToroidalBlockLevelExchangeCenter(a_geom, discretization_order, discretization_order);
   }
   else {
      if ( a_geom.extrablockExchange() ) {
         m_mblx_ptr = new MagFluxAlignedMBLEXCenter;
         m_mblx_ptr->define(&a_geom, discretization_order, discretization_order);
      }
      else {
         m_mblx_ptr = NULL;
      }
   }

#ifdef with_petsc
   m_preconditioner = new MBPETScSolver(a_geom, discretization_order, m_mblx_ptr);
#else
   m_preconditioner = new MBHypreSolver(a_geom, discretization_order, m_mblx_ptr);
#endif
   
   m_preconditioner->setMethodParams(m_precond_method, m_precond_precond_method);
}
      


GKPoisson::~GKPoisson()
{
   if (m_preconditioner) delete m_preconditioner;
   if (m_mblx_ptr) delete m_mblx_ptr;
}


void
GKPoisson::updateBoundaries( const EllipticOpBC&  a_bc )
{
   // N.B.: The operator coefficient must computed prior to calling this function
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

   computeCoefficients( a_ion_mass_density, 
                        m_mapped_coefficients, 
                        m_unmapped_coefficients );

   updateBoundaries(a_bc);

   if ( a_update_preconditioner ) {
      m_preconditioner->constructMatrix(  m_volume_reciprocal, 
                                          m_mapped_coefficients, 
                                          a_bc);
   }
}

void
GKPoisson::setOperatorCoefficients( const PS::KineticSpeciesPtrVect&  a_kin_species,
                                    const LevelData<FArrayBox>&       a_ion_mass_density,
                                    const EllipticOpBC&               a_bc,
                                    const bool                        a_update_preconditioner )
{
   CH_TIME("GKPoisson::setOperatorCoefficients");
   static bool first_call = true;

   computeCoefficients( a_ion_mass_density, 
                        m_mapped_coefficients, 
                        m_unmapped_coefficients );

   if (m_include_FLR_effects) {

     m_species_vec = a_kin_species;

     if (first_call) {
       m_FLR_integrand_factor.resize(a_kin_species.size());
       for (int s = 0; s < a_kin_species.size(); s++) {
         m_FLR_integrand_factor[s]->define(*m_phase_grids, 1, PS::IntVect::Zero);
       }
       first_call = false;
     }

     computeFLRIntegrandFactor();

     computeCoefficientsGK( m_mapped_coefficients_gk, 
                            m_unmapped_coefficients_gk );

   }

   updateBoundaries(a_bc);

   if ( a_update_preconditioner ) {
      m_preconditioner->constructMatrix(  m_volume_reciprocal, 
                                          m_mapped_coefficients, 
                                          a_bc);
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
         
         if ( dir == RADIAL_DIR && block_number < 2 ) {
            ones[dit].setVal(1.,dir);
         }
      }
   }
   
   //Compute <c^2*n*m/(e*B^2)*|grad(Psi)|^2>
   computeRadialFSAverage(ones, a_lo_value, a_hi_value, a_radial_gkp_divergence_average);
}


void
GKPoisson::computeCoefficients( const LevelData<FArrayBox>& a_ion_mass_density,
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
      
      density_sum_cell[dit].copy(a_ion_mass_density[dit]);
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
GKPoisson::computeCoefficientsGK( LevelData<FluxBox>&         a_mapped_coefficients,
                                  LevelData<FluxBox>&         a_unmapped_coefficients )
{
   CH_TIME("GKPoisson::computeCoefficientsGK");
   CH_assert(m_include_FLR_effects);

   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   const IntVect grown_ghosts( m_mapped_coefficients.ghostVect() + IntVect::Unit );
   CH_assert(grown_ghosts == 2*IntVect::Unit);
   
   LevelData<FluxBox> grown_mapped_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);
   LevelData<FluxBox> grown_unmapped_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);
   
   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& par_coeff = m_geometry.getEllipticOpParCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {

      const Box& box = a_unmapped_coefficients[dit].box();
      
      FluxBox isotropic_coeff(box, SpaceDim * SpaceDim);
      isotropic_coeff.copy(perp_coeff[dit]);
      isotropic_coeff += par_coeff[dit];
      isotropic_coeff *= m_debye_number2;

      FluxBox isotropic_coeff_mapped(box, SpaceDim * SpaceDim);
      isotropic_coeff_mapped.copy(perp_coeff_mapped[dit]);
      isotropic_coeff_mapped += par_coeff_mapped[dit];
      isotropic_coeff_mapped *= m_debye_number2;
      
      if (m_model == "GyroPoisson") {
          grown_mapped_coefficients[dit].copy(isotropic_coeff_mapped);
          grown_unmapped_coefficients[dit].copy(isotropic_coeff);
      }

      else {
         MayDay::Error("GKPoisson:: unknown model is specified");
      }
   }

   grown_unmapped_coefficients.exchange();
   grown_mapped_coefficients.exchange();
   
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_unmapped_coefficients[dit].copy(grown_unmapped_coefficients[dit]);
   }
   
   if (!m_second_order) {
      // The mapped coefficients must now be converted to face averages, which
      // requires a layer of transverse ghost faces that we don't have at the
      // physical boundary, so we need to extrapolate them.
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);
   
      // Convert the mapped coefficients from face-centered to face-averaged
      fourthOrderAverage(grown_mapped_coefficients);
   
      m_geometry.fillTransverseGhosts(a_unmapped_coefficients, false);
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);
   }
   
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_mapped_coefficients[dit].copy(grown_mapped_coefficients[dit]);
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
GKPoisson::setPreconditionerConvergenceParams( const double a_tol,
                                               const int    a_max_iter,
                                               const double a_precond_tol,
                                               const int    a_precond_max_iter )
{
   m_preconditioner->setConvergenceParams(a_tol, a_max_iter, m_precond_verbose,
                                          a_precond_tol, a_precond_max_iter,
                                          m_precond_precond_verbose);
}



void
GKPoisson::solvePreconditioner( const LevelData<FArrayBox>& a_in,
                                LevelData<FArrayBox>&       a_out )
{
   m_preconditioner->solve(a_in, a_out, true);
}



void
GKPoisson::multiplyCoefficients( LevelData<FluxBox>& a_data,
                                 const bool a_mapped_coeff,
                                 const bool a_apply_op) const
{
   CH_assert(a_data.ghostVect() <= m_unmapped_coefficients.ghostVect());

   if (a_apply_op && m_include_FLR_effects) {

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
                    tmp.copy(m_mapped_coefficients_gk[dit][dir], SpaceDim*i+j, 0, 1);
                 }
                 else {
                    tmp.copy(m_unmapped_coefficients_gk[dit][dir], SpaceDim*i+j, 0, 1);
                 }
                 tmp.mult(this_saved_data_dir, j, 0, 1);
                 this_data_dir.plus(tmp, 0, i, 1);
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
                    tmp.copy(m_mapped_coefficients[dit][dir], SpaceDim*i+j, 0, 1);
                 }
                 else {
                    tmp.copy(m_unmapped_coefficients[dit][dir], SpaceDim*i+j, 0, 1);
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
GKPoisson::parseParameters( const ParmParse&   a_ppntr )
{
   if (a_ppntr.contains("model")) {
      a_ppntr.get( "model", m_model );
   }
   a_ppntr.query( "include_FLR_effects", m_include_FLR_effects );
}


void
GKPoisson::printParameters()
{
   if (procID()==0) {
      std::cout << "GKPoisson model: " << m_model << std::endl;
      std::cout << "GKPoisson include FLR effects: " 
                << (m_include_FLR_effects ? "true" : "false") << std::endl;
   }
}

void
GKPoisson::applyOp( LevelData<FArrayBox>&       a_out,
                    const LevelData<FArrayBox>& a_in,
                    bool                        a_homogeneous )
{
  EllipticOp::applyOp(a_out, a_in, a_homogeneous);

  if (m_include_FLR_effects) {

    /* inject phi to phase space */
    PS::LevelData<PS::FArrayBox> phi_flat;
    m_phase_geom->injectConfigurationToPhase(a_in, phi_flat);

    /* create a "full-sized" phi in phase space */
    PS::LevelData<PS::FArrayBox> phi( *m_phase_grids, 
                                      phi_flat.nComp(), 
                                      PS::IntVect::Zero );
    for (PS::DataIterator dit(m_phase_grids->dataIterator()); dit.ok(); ++dit) {
      const PS::FArrayBox& phi_flat_fab = phi_flat[dit];
      const PS::Box& flat_box = phi_flat_fab.box();
      PS::FArrayBox& phi_fab = phi[dit];
      const PS::Box& big_box = phi_fab.box();

      VEL::IntVect vel_idx;
      for (int d = VPARALLEL_DIR; d <= MU_DIR; d++) {
        int imin = flat_box.smallEnd(d);
        int imax = flat_box.bigEnd(d);
        CH_assert(imin == imax);
        vel_idx[d-VPARALLEL_DIR] = imin;
      }

      for (int n = 0; n < phi_fab.nComp(); n++) {
        for (PS::BoxIterator bit(big_box); bit.ok(); ++bit) {
          IntVect cfg_idx = m_phase_geom->config_restrict(bit());
          PS::IntVect idx = m_phase_geom->tensorProduct( cfg_idx, vel_idx );
          phi_fab(bit(), n) = phi_flat_fab(idx, n);
        }
      }
    }

    const LevelData<FArrayBox>& BFieldMag = m_geometry.getCCBFieldMag();

    for (int s = 0; s < m_species_vec.size(); s++) {

      const PS::KineticSpecies& this_species = *(m_species_vec[s]);
      if (!this_species.isGyrokinetic()) continue;
      const PS::GyroaverageOperator* gyroavg_op = this_species.gyroaverageOp();
      const PS::LevelData<PS::FArrayBox>& factor = *(m_FLR_integrand_factor[s]);
   
      /* gyroaverage this phi twice */
      PS::LevelData<PS::FArrayBox> phi_bar, phi_tilde;
      gyroavg_op->applyOp(phi_bar, phi);
      gyroavg_op->applyOp(phi_tilde, phi_bar);

      PS::LevelData<PS::FArrayBox> integrand(*m_phase_grids, 1, PS::IntVect::Zero);
      for (PS::DataIterator dit(integrand.dataIterator()); dit.ok(); ++dit) {
        PS::FArrayBox& integrand_fab = integrand[dit];
        const PS::FArrayBox& phi_fab = phi[dit];
        const PS::FArrayBox& phi_tilde_fab = phi_tilde[dit];
        const PS::FArrayBox& factor_fab = factor[dit];
        const PS::Box& box = (*m_phase_grids)[dit];

        for (PS::BoxIterator bit(box); bit.ok(); ++bit) {
          integrand_fab(bit(),0) = factor_fab(bit(),0)
                                    * ( phi_fab(bit(),0) - phi_tilde_fab(bit(),0));
        }

      }

      /* compute polarization term and add it to a_out */
      LevelData<FArrayBox> n_pol_species( a_out.disjointBoxLayout(),
                                          1,
                                          IntVect::Zero );
      m_moment_op.compute(  n_pol_species, 
                            this_species, 
                            integrand, 
                            PS::DensityKernel() );

      for (DataIterator dit(a_out.dataIterator()); dit.ok(); ++dit) {

        const Box& box = n_pol_species[dit].box();
        FArrayBox polarization_fac(box, 1);
        polarization_fac.setVal(1.0);
        //polarization_fac *= m_larmor_number2;
        polarization_fac *= ( this_species.charge() * this_species.charge() );
        polarization_fac.divide(BFieldMag[dit], box, 0, 0);

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

    const VEL::VelCoordSys& vel_coords = m_phase_geom->velSpaceCoordSys();
    const VEL::ProblemDomain& vel_domain = vel_coords.domain();
    const VEL::RealVect& vel_dx = vel_coords.dx();
    Real dmu = vel_dx[1];

    PS::LevelData<PS::FArrayBox>& factor = *(m_FLR_integrand_factor[s]);

    PS::LevelData<PS::FArrayBox> inv_bstar_par(  *m_phase_grids, 
                                                 1, 
                                                 PS::IntVect::Zero );
    for (PS::DataIterator dit(factor.dataIterator()); dit.ok(); ++dit) {
      inv_bstar_par[dit].setVal(1.0);
    }
    m_phase_geom->divideBStarParallel(inv_bstar_par);

    for (PS::DataIterator dit(factor.dataIterator()); dit.ok(); ++dit) {

      const PS::FArrayBox& dfn_fab = dfn[dit];
      const PS::FArrayBox& inv_bstar_par_fab = inv_bstar_par[dit];
      PS::FArrayBox& factor_fab = factor[dit];

      const PS::Box& box = (*m_phase_grids)[dit];

      for (PS::BoxIterator bit(box); bit.ok(); ++bit) {

        Real bstar_par = 1.0/inv_bstar_par_fab(bit(),0);
        
        Real df_dmu = (1.0 / (2*dmu) )
                      * (   dfn_fab(bit()+PS::BASISV(MU_DIR),0)
                          - dfn_fab(bit()-PS::BASISV(MU_DIR),0) );

        factor_fab(bit(),0) = bstar_par * df_dmu;

      }

    }

  }

  return;
}

#include "NamespaceFooter.H"
