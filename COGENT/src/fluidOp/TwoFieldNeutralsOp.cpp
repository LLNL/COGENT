#include "TwoFieldNeutralsOp.H"
#include "EllipticOpBCFactory.H"
#include "FluidVarBCFactory.H"
#include "FourthOrderUtil.H"
#include "SimpleDivergence.H"
#include "SpaceUtils.H.multidim"
#include "Directions.H"
#include "mappedAdvectionFlux.H"
#include "FourthOrderUtil.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "FluidOpPreconditioner.H"
//#include "altFaceAverages.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NeutralsOpF_F.H"

#include "NamespaceHeader.H" 

const char* TwoFieldNeutralsOp::pp_name = {"two_field_neutrals_op"};


TwoFieldNeutralsOp::TwoFieldNeutralsOp(const string&   a_pp_str,
                                       const string&   a_species_name,
                                       const MagGeom&  a_geometry,
                                       const Real      a_larmor,
                                       const int       a_verbosity )
   : OneFieldNeutralsOp(a_pp_str, a_species_name, a_geometry, a_larmor, a_verbosity),
     m_geometry(a_geometry),
     m_verbosity(a_verbosity),
     m_my_pc_idx(-1),
     m_bckgr_density(1.0e-7),
     m_is_time_implicit(true),
     m_include_advection_bc(true),
     m_fixed_recycling(false),
     m_first_call(true),
     m_extrapolate_density(false),
     m_advScheme("uw3"),
     m_opt_string(a_pp_str),
     m_analytic_pc(false),
     m_twofieldneutrals_pc(nullptr)
{

   // Overwrite the base class value with
   // true by default
   m_consistent_diffusion = true;
   
   ParmParse pp(a_pp_str.c_str());
   parseParameters( pp );
   if (m_verbosity>0) {
      printParameters();
   }

   std::string name_op;
   std::string prefixBC;
   EllipticOpBCFactory elliptic_op_bc_factory;

   // Construct the preconditioner
   name_op = "neutrals_pc";
   ParmParse pppc( name_op.c_str() );
   m_twofieldneutrals_pc = new TwoFieldNeutralsPC(pppc, a_geometry);

   // Construct diffusion elliptic operator
   name_op = "neutrals_diffusion";

   ParmParse ppdiff( name_op.c_str() );
   m_diffusion_op = new Diffusion(ppdiff, a_geometry);

   prefixBC = "BC." + name_op;
   ParmParse ppdiff_bc( prefixBC.c_str() );
   m_diffusion_bc = elliptic_op_bc_factory.create(name_op,
                                                  ppdiff_bc,
                                                  *(a_geometry.getCoordSys()),
                                                  false );
   
   // Construct viscosity elliptic operator
   name_op = "neutrals_viscosity";

   ParmParse ppvisc( name_op.c_str() );
   m_viscosity_op = new Diffusion(ppvisc, a_geometry);

   prefixBC = "BC." + name_op;
   ParmParse ppvisc_bc( prefixBC.c_str() );
   m_viscosity_bc = elliptic_op_bc_factory.create(name_op,
                                                  ppvisc_bc,
                                                  *(a_geometry.getCoordSys()),
                                                  false );
   

   std:: string variable_name;

   // Get initial and boundary conditions for mass density (m_g * n) variable
   variable_name = "density";
   parseInitialConditions(a_species_name, variable_name);
   
   FluidVarBCFactory fluid_var_bc_factory;
   
   FluidVarBC* dens_var_bc = fluid_var_bc_factory.create(a_species_name, variable_name, a_geometry.getCoordSys()->type(), false);
   m_fluid_bc.push_back( dens_var_bc );
      
   // Check for recyling BCs for density
   m_dens_recycling_bc = dens_var_bc->isRecyclingBC();
   
   if (m_dens_recycling_bc) {
      // Set flags to ignore standard rho_g BCs (e.g., Dirichlet)
      // set by input parameters, and use extrapolation to fill ghosts
      m_extrapolate_density = true;
      m_include_advection_bc = false;
   }
   
   // Get initial and boundary conditions for parallel momentum (m_g * n_g * V_g_par) variable
   variable_name = "momentumDensity_virtual";
   parseInitialConditions(a_species_name, variable_name);
   
   FluidVarBC* mom_var_bc = fluid_var_bc_factory.create(a_species_name, variable_name, a_geometry.getCoordSys()->type(), false);
   m_fluid_bc.push_back( mom_var_bc );
      
   //Check for recyling BCs for momentum
   m_mom_recycling_bc = mom_var_bc->isRecyclingBC();
   
   if (m_mom_recycling_bc) {
      m_include_advection_bc = false;
   }
   
   //Get disjoint box layout
   const DisjointBoxLayout& grids( a_geometry.grids() );
   
   if (m_consistent_velocity) {
      // This ghost vect is already set by the base class
      // now try to optimize expensive moment calculation
      // by computing it in valid cells only and then
      // by using extrapolation routins at phys. boundaries
      m_ghostVect = IntVect::Zero;
   }
   
   // Create atomic data objects
   m_chx_rate.define(grids, 1, m_ghostVect);
   m_ionization_rate.define(grids, 1, IntVect::Zero);
   m_recombination_rate.define(grids, 1, IntVect::Zero);
   
   
   // Consistency check
   if (!m_diffusion_op->secondOrder()  && m_dens_recycling_bc) {
     MayDay::Error("TwoFieldNeutralsOp:: recycling BCs will introduce discontinuity in 4th order diffusion ");
   }
   if (!m_geometry.secondOrder() && !m_include_advection_bc) {
      MayDay::Warning("TwoFieldNeutralsOp:: using codim-1 extrapolation with fourth-order");
   }
}


TwoFieldNeutralsOp::~TwoFieldNeutralsOp()
{

   for (std::vector<FluidVarBC*>::iterator it = m_fluid_bc.begin(); it != m_fluid_bc.end(); ++it)
   {
      delete (*it);
   }
   m_fluid_bc.clear();

   if ( m_twofieldneutrals_pc ) delete m_twofieldneutrals_pc;
   
   delete m_diffusion_op;
   delete m_diffusion_bc;
   
   delete m_viscosity_op;
   delete m_viscosity_bc;

}


void TwoFieldNeutralsOp::accumulateRHS(FluidSpeciesPtrVect&               a_rhs,
                                       const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                       const FluidSpeciesPtrVect&         a_fluid_species,
                                       const PS::ScalarPtrVect&           a_scalars,
                                       const EField&                      a_E_field,
                                       const int                          a_fluidVecComp,
                                       const Real                         a_time)
{
   // Get prescribed plasma parameters
   getPlasmaParameters(a_time);
   
   // Get atomic rates
   computeAtomicRates(m_ionization_rate, m_recombination_rate, m_chx_rate, a_time);
   
   // Get fluid rhs
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );

   // Get fluid soln
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );

   //Accumulate diffusive terms
   accumulateDiffusiveTerms(rhs_fluid, soln_fluid, a_time);

   //Accumulate parallel advection terms
   accumulateParAdvectionTerms(rhs_fluid, soln_fluid, a_time);
   
   //Accumulate perpendicular advection terms
   accumulatePerpAdvectionTerms(rhs_fluid, soln_fluid, a_time);

   //Accumulate source terms
   accumulateSourceTerms(rhs_fluid, soln_fluid, a_time);
   
   m_first_call = false;
}


void TwoFieldNeutralsOp::accumulateExplicitRHS(FluidSpeciesPtrVect&               a_rhs,
                                               const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                               const FluidSpeciesPtrVect&         a_fluid_species,
                                               const PS::ScalarPtrVect&           a_scalars,
                                               const EField&                      a_E_field,
                                               const int                          a_fluidVecComp,
                                               const Real                         a_time)
{
   if (!m_is_time_implicit) {
     accumulateRHS(a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, a_fluidVecComp, a_time);
   }
}


void TwoFieldNeutralsOp::accumulateImplicitRHS(FluidSpeciesPtrVect&               a_rhs,
                                               const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                               const FluidSpeciesPtrVect&         a_fluid_species,
                                               const PS::ScalarPtrVect&           a_scalars,
                                               const EField&                      a_E_field,
                                               const int                          a_fluidVecComp,
                                               const Real                         a_time)
{
  if (m_is_time_implicit) {
    accumulateRHS(a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, a_fluidVecComp, a_time);
  }
}


void TwoFieldNeutralsOp::accumulateDiffusiveTerms(FluidSpecies&        a_fluid_rhs,
                                                  const FluidSpecies&  a_fluid_soln,
                                                  const Real           a_time) const
{
   /* Computes the diffusive terms RHS for density and momentum eqs,
    density_rhs = \nabla_perp D \nabla_perp (rho_g * T_g), D = 1/(m_g * ne * <sigma_chx*Vti> )
    par_mom_rhs = \nabla \eta \nabla (v_par_g), \eta =rho_g * T_g/(m_g * ne * <sigma_chx*Vti> )
    */
   
   // Get rhs data
   LevelData<FArrayBox>& rhs_density( a_fluid_rhs.cell_var("density") );
   LevelData<FArrayBox>& rhs_par_mom( a_fluid_rhs.cell_var("momentumDensity_virtual") );
   
   // Get soln data
   // NB: density component is mass density rho_g = m_g * n_g
   // NB: momentumDensity_virtual is m_g * n_g * V_g
   const LevelData<FArrayBox>& soln_density( a_fluid_soln.cell_var("density") );

   //Get species mass
   const Real species_mass = a_fluid_soln.mass();
   
   // Get boxes
   const DisjointBoxLayout& grids( rhs_density.getBoxes() );

   // Create temporary object to store flux divergence
   LevelData<FArrayBox> flux_div(grids, 1, IntVect::Zero);

   // Create elliptic coeffiecient
   LevelData<FluxBox> ellip_coeff(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   LevelData<FluxBox> ellip_coeff_mapped(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);

   // Compute diffusion term
   // NB: here pressure = m_g * n_g * T_g
   LevelData<FArrayBox> neutr_pressure(grids, 1, IntVect::Zero);
   for (DataIterator dit( neutr_pressure.dataIterator() ); dit.ok(); ++dit) {
      neutr_pressure[dit].copy(soln_density[dit]);
      neutr_pressure[dit].mult(m_Tg[dit]);
   }
   computeDiffusionCoefficients(ellip_coeff, ellip_coeff_mapped, species_mass, a_time);

   m_diffusion_op->setOperatorCoefficients(ellip_coeff, ellip_coeff_mapped, *m_diffusion_bc);
   if (m_analytic_pc) {
      m_twofieldneutrals_pc->updateTopLeftCoefficients(ellip_coeff, ellip_coeff_mapped);
   }

   if (!m_dens_recycling_bc) {
      m_diffusion_op->computeFluxDivergence( neutr_pressure, flux_div, false, false);
   }
   else {
      computeDiffusionWithRecyclingBC(flux_div, *m_fluid_bc.at(0), neutr_pressure,
                                      *m_diffusion_op, m_ion_normal_flux, a_time);
   }

   // The preceding call computes a negated physical flux divergence cell average.
   // Since this function is being used to define the right-hand side of an evolution
   // equation for a mapped unknown (i.e., containing a J factor), we therefore need
   // to multiply by J.  We also flip the sign.
   m_geometry.multJonValid(flux_div);
   
   for (DataIterator dit( flux_div.dataIterator() ); dit.ok(); ++dit) {
      rhs_density[dit] -= flux_div[dit];
   }
   
   // Compute viscosity term
   LevelData<FArrayBox> parallel_velocity(grids, 1, IntVect::Zero);
   a_fluid_soln.velocity_virtual(parallel_velocity);

   computeViscosityCoefficients(ellip_coeff, ellip_coeff_mapped, soln_density, species_mass, a_time);
   // Set recycling BC for parallel velocity
   if (m_mom_recycling_bc) {
      m_fluid_bc.at(1)->applyRecyclingBC(*m_viscosity_bc, m_ion_parallel_vel_face, a_time);
   }

   m_viscosity_op->setOperatorCoefficients(ellip_coeff, ellip_coeff_mapped, *m_viscosity_bc);
   if (m_analytic_pc) {
      m_twofieldneutrals_pc->updateBottomRightCoefficients(ellip_coeff, ellip_coeff_mapped);
   }
   
   m_viscosity_op->computeFluxDivergence( parallel_velocity, flux_div, false, false);
   m_geometry.multJonValid(flux_div);
   
   for (DataIterator dit( flux_div.dataIterator() ); dit.ok(); ++dit) {
     rhs_par_mom[dit] -= flux_div[dit];
   }
}


void TwoFieldNeutralsOp::accumulateParAdvectionTerms(FluidSpecies&        a_fluid_rhs,
                                                     const FluidSpecies&  a_fluid_soln,
                                                     const Real           a_time) const
{
   
   /* Computes the parallel advection terms RHS for the density and momentum equations,
    density_rhs = div (v_par_g * bvec * rho_g)
    par_mom_rhs = div (v_par_g * bvec * rho_g * v_par_g )
   */
   
   // Get rhs data
   LevelData<FArrayBox>& rhs_density( a_fluid_rhs.cell_var("density") );
   LevelData<FArrayBox>& rhs_par_mom( a_fluid_rhs.cell_var("momentumDensity_virtual") );
   
   // Get soln data
   const LevelData<FArrayBox>& soln_density( a_fluid_soln.cell_var("density") );
   const LevelData<FArrayBox>& soln_par_mom( a_fluid_soln.cell_var("momentumDensity_virtual") );
   
   // Get boxes
   const DisjointBoxLayout& grids( rhs_density.getBoxes() );
   
   // Get parallel velocity
   LevelData<FArrayBox> parallel_vel(grids, 1, soln_par_mom.ghostVect());
   a_fluid_soln.velocity_virtual(parallel_vel);
   
   if (!m_include_advection_bc) {
      bool fourth_order = (m_geometry.secondOrder()) ? false : true;
      m_geometry.extrapolateToPhysicalGhosts(parallel_vel, fourth_order);
   }
   
   // Get parallel velocity on faces
   LevelData<FluxBox> velocity_on_faces( grids, 1, IntVect::Unit );
   convertCellToFace(velocity_on_faces, parallel_vel);
       
   // Force zero flux on radial faces for the flux-aligned grids
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      if (block_coord_sys.isFieldAligned()) {
         velocity_on_faces[dit][RADIAL_DIR].setVal(0.0);
      }
   }
   // Get unit b vector
   const LevelData<FluxBox>& b_dir = m_geometry.getFCBFieldDir();
   
   // Compute parallel velocity (physical vector)
   LevelData<FluxBox> velocity( grids, SpaceDim, IntVect::Unit );
#if CFG_DIM == 2
   m_geometry.projectPoloidalVector(b_dir, velocity);
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
         for (int n = 0; n < SpaceDim; n++) {
            velocity[dit][dir].mult(velocity_on_faces[dit][dir],0,n,1);
         }
      }
   }
#else
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      velocity[dit].copy(b_dir[dit]);
      for (int dir = 0; dir < SpaceDim; dir++) {
         for (int n = 0; n < SpaceDim; n++) {
           velocity[dit][dir].mult(velocity_on_faces[dit][dir],0,n,1);
         }
      }
   }
#endif
   
   // Apply axisymmetric corrections (for pointwsie metrics does nothing)
   m_geometry.applyAxisymmetricCorrection(velocity);
   
   // Compute parallel advection RHS in eqn for density
   LevelData<FArrayBox> rho_rhs( grids, 1, IntVect::Zero );
   computeAdvectionRHS(rho_rhs, soln_density, velocity, m_dens_recycling_bc, true, a_time);

   // Compute parallel advection RHS in eqn for parallel momentum
   LevelData<FArrayBox> par_mom_rhs( grids, 1, IntVect::Zero );
   computeAdvectionRHS(par_mom_rhs, soln_par_mom, velocity, m_dens_recycling_bc, true, a_time);
   
   // Add parallel advection to RHS of the entire system
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      rhs_density[dit].plus(rho_rhs[dit]);
      rhs_par_mom[dit].plus(par_mom_rhs[dit]);
   }
}

void TwoFieldNeutralsOp::accumulatePerpAdvectionTerms(FluidSpecies&        a_fluid_rhs,
                                                      const FluidSpecies&  a_fluid_soln,
                                                      const Real           a_time) const
{
   
   /* Computes the perpendicular advection terms RHS for the momentum equations,
    par_mom_rhs = div (v_perp_g * rho_g * v_par_g )
    v_perp_g = - 1/(m_g * rho_g * ne * <sigma_chx*Vti>) * grad_perp (rho_g * T_g)
    */
   
   // Get rhs data
   LevelData<FArrayBox>& rhs_par_mom( a_fluid_rhs.cell_var("momentumDensity_virtual") );
   
   // Get soln data
   const LevelData<FArrayBox>& soln_density( a_fluid_soln.cell_var("density") );
   const LevelData<FArrayBox>& soln_par_mom( a_fluid_soln.cell_var("momentumDensity_virtual") );
   const Real species_mass = a_fluid_soln.mass();
   
   // Get boxes
   const DisjointBoxLayout& grids( soln_density.getBoxes() );
   
   // Create a temporary for soln_density
   const IntVect& ghosts = soln_density.ghostVect();
   LevelData<FArrayBox> soln_density_tmp(grids, 1, ghosts);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
     soln_density_tmp[dit].copy(soln_density[dit]);
   }
   
   if (!m_include_advection_bc) {
      bool fourth_order = (m_geometry.secondOrder()) ? false : true;
      m_geometry.extrapolateToPhysicalGhosts(soln_density_tmp, fourth_order);
   }
   
   // Compute perpendicular velocity (physical vector)
   LevelData<FluxBox> velocity( grids, SpaceDim, IntVect::Unit );
   computePerpPressureGradient(velocity, soln_density_tmp, species_mass);

   // Get scalar factors
   LevelData<FArrayBox> fac_cell(grids, 1, 2*IntVect::Unit);
   LevelData<FluxBox> fac_face(grids, 1, IntVect::Unit);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      fac_cell[dit].copy(m_chx_rate[dit]);
      fac_cell[dit].mult(soln_density_tmp[dit]);
      fac_cell[dit].invert(-1.0);
   }
   
   if (m_ghostVect == IntVect::Zero || m_extrapolate_density) {
      bool fourth_order = (m_geometry.secondOrder()) ? false : true;
      m_geometry.extrapolateToPhysicalGhosts(fac_cell, fourth_order);
   }
   convertCellToFace(fac_face, fac_cell);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         for (int n = 0; n < SpaceDim; n++) {
            velocity[dit][dir].mult(fac_face[dit][dir],0,n,1);
         }
      }
   }
   
   // Apply axisymmetric corrections (for pointwsie metrics does nothing)
   m_geometry.applyAxisymmetricCorrection(velocity);
   
   // Compute perpendicular advection RHS in eqn for parallel momentum
   LevelData<FArrayBox> par_mom_rhs( grids, 1, IntVect::Zero );
   computeAdvectionRHS(par_mom_rhs, soln_par_mom, velocity, m_dens_recycling_bc, false, a_time);
   
   // Add perpendicular advection to RHS of the entire system
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      rhs_par_mom[dit].plus(par_mom_rhs[dit]);
   }
}

void TwoFieldNeutralsOp::accumulateSourceTerms(FluidSpecies&        a_fluid_rhs,
                                               const FluidSpecies&  a_fluid_soln,
                                               const Real           a_time) const
{
   /* Computes the source terms RHS for the density and momentum equations,
    density_rhs = rho_g * ne * (<sigma_rec*Vte>-<sigma_iz*Vte>)
    par_mom_rhs = rho_g * ne * <sigma_chx*Vti> * (v_par_i - v_par_g) - grad_par(rho_g * T_g / m_g)
    */
   
   // Get rhs data
   LevelData<FArrayBox>& rhs_density( a_fluid_rhs.cell_var("density") );
   LevelData<FArrayBox>& rhs_par_mom( a_fluid_rhs.cell_var("momentumDensity_virtual") );
   
   // Get soln data
   const LevelData<FArrayBox>& soln_density( a_fluid_soln.cell_var("density") );
   const LevelData<FArrayBox>& soln_par_mom( a_fluid_soln.cell_var("momentumDensity_virtual") );
   
   // Get boxes
   const DisjointBoxLayout& grids( rhs_density.getBoxes() );
   
   // Get species mass
   const Real species_mass = a_fluid_soln.mass();
   
   // Compute ionization and recombination terms
   LevelData<FArrayBox> dens_src(grids, 1, IntVect::Zero);
   for (DataIterator dit( dens_src.dataIterator() ); dit.ok(); ++dit) {
      dens_src[dit].copy(m_recombination_rate[dit]);
      dens_src[dit].minus(m_ionization_rate[dit]);
      dens_src[dit].mult(soln_density[dit]);

      FORT_ADD_BCKGR_SRC(CHF_BOX(dens_src[dit].box()),
                         CHF_CONST_FRA1(m_ionization_rate[dit],0),
                         CHF_CONST_FRA1(soln_density[dit],0),
                         CHF_CONST_REAL(m_bckgr_density),
                         CHF_CONST_REAL(species_mass),
                         CHF_FRA1(dens_src[dit],0));
   }
   
   m_geometry.multJonValid(dens_src);
   
   // Compute charge-exchange momentum src
   LevelData<FArrayBox> chx_friction(grids, 1, IntVect::Zero);
   for (DataIterator dit( chx_friction.dataIterator() ); dit.ok(); ++dit) {
      chx_friction[dit].copy(soln_density[dit]);
      chx_friction[dit].mult(m_ion_parallel_vel[dit]);
      chx_friction[dit].minus(soln_par_mom[dit]);
      chx_friction[dit].mult(m_chx_rate[dit]);
   }
   m_geometry.multJonValid(chx_friction);

   // Compute ionization momentum sink
   LevelData<FArrayBox> ioniz_sink(grids, 1, IntVect::Zero);
   for (DataIterator dit( ioniz_sink.dataIterator() ); dit.ok(); ++dit) {
      ioniz_sink[dit].copy(soln_par_mom[dit]);
      ioniz_sink[dit].mult(m_ionization_rate[dit]);
   }
   m_geometry.multJonValid(ioniz_sink);

   // Compute recombination momentum source
   LevelData<FArrayBox> recomb_src(grids, 1, IntVect::Zero);
   for (DataIterator dit( recomb_src.dataIterator() ); dit.ok(); ++dit) {
      recomb_src[dit].copy(soln_density[dit]);
      recomb_src[dit].mult(m_ion_parallel_vel[dit]);
      recomb_src[dit].mult(m_recombination_rate[dit]);
   }
   m_geometry.multJonValid(recomb_src);

   // Compute poloidal component of the pressure parallel force
   LevelData<FArrayBox> pressure_grad(grids, 1, IntVect::Zero);
   computeParallelPressureGradient(pressure_grad, soln_density, species_mass);

   m_geometry.multJonValid(pressure_grad);
   
   // Add sources to fluid equations rhs
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      rhs_density[dit].plus(dens_src[dit]);
      rhs_par_mom[dit].plus(chx_friction[dit]);
      rhs_par_mom[dit].plus(recomb_src[dit]);
      rhs_par_mom[dit].minus(ioniz_sink[dit]);
      rhs_par_mom[dit].minus(pressure_grad[dit]);
   }
}

void TwoFieldNeutralsOp::computeAdvectionRHS(LevelData<FArrayBox>&         a_div,
                                             const LevelData<FArrayBox>&   a_dist,
                                             const LevelData<FluxBox>&     a_velocity,
                                             const bool                    a_recycling_bc,
                                             const bool                    a_homogeneous_flux_bc,
                                             const Real                    a_time) const
{
   
   // Get boxes
   const DisjointBoxLayout& grids( a_div.getBoxes() );
   
   bool fourth_order = (m_geometry.secondOrder()) ? false : true;
   
   // Get mapped components of parallel velocity
   LevelData<FluxBox> normal_vel( grids, 1, IntVect::Unit );
   m_geometry.computeMetricTermProductAverage(normal_vel, a_velocity, false);
   
   // Get distribution on faces
   LevelData<FluxBox> faceDist( grids, 1, IntVect::Unit );
   if (m_include_advection_bc) {
      SpaceUtils::upWindToFaces(faceDist, a_dist, normal_vel, m_advScheme, fourth_order);
   }
   else {
      // Overwrite advection quantity BCs with extrapolation
      const IntVect& ghosts = a_dist.ghostVect();
      LevelData<FArrayBox> tmp(grids, 1, ghosts );
      for (DataIterator dit( a_dist.dataIterator() ); dit.ok(); ++dit) {
         tmp[dit].copy(a_dist[dit]);
      }
      m_geometry.extrapolateToPhysicalGhosts(tmp, fourth_order);
      SpaceUtils::upWindToFaces(faceDist, tmp, normal_vel, m_advScheme, fourth_order);
   }
      
   // Compute flux
   LevelData<FluxBox> flux( grids, CFG_DIM, IntVect::Unit );
   computeFlux(flux, faceDist, a_velocity, fourth_order);
   
   // Compute integrated normal flux
   LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
   m_geometry.computeMetricTermProductAverage(NTF_normal, flux, fourth_order);
   
   // Average the normal component of NTF consistent across block interfaces
   m_geometry.averageAtBlockBoundaries(NTF_normal);
   
   // Apply recycling BC
   if (a_recycling_bc) {
      
      // Apply zero flux BC
      if (a_homogeneous_flux_bc) {
         LevelData<FluxBox> zero_flux(grids, 1, IntVect::Zero);
         setZero(zero_flux);
         m_fluid_bc.at(0)->applyRecyclingBC(NTF_normal, zero_flux, a_time);
      }
      
      // Apply ion flux BC. We should only do it for the perpendicular ion momentum advection.
      else {
	LevelData<FluxBox> tmp(grids, 1, IntVect::Zero);
	for (DataIterator dit(grids); dit.ok(); ++dit) {
	  for (int dir=0; dir<SpaceDim; ++dir) {
	    tmp[dit][dir].copy(m_ion_normal_flux[dit][dir]);
	    tmp[dit][dir].mult(faceDist[dit][dir]);
	  }
	}
	m_fluid_bc.at(0)->applyRecyclingBC(NTF_normal, tmp, a_time);
      }
   }
   
   RealVect fakeDx = RealVect::Unit;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      simpleDivergence(a_div[dit], NTF_normal[dit], grids[dit], fakeDx);
   }
   
   // Divide by cell volume and negate
   for (DataIterator dit( a_div.dataIterator() ); dit.ok(); ++dit) {
      const MagBlockCoordSys& coords( m_geometry.getBlockCoordSys( grids[dit] ) );
      double fac( -1.0 / coords.getMappedCellVolume() );
      a_div[dit].mult( fac );
   }
}


void TwoFieldNeutralsOp::computeFlux(LevelData<FluxBox>&       a_flux,
                                     const LevelData<FluxBox>& a_faceDist,
                                     const LevelData<FluxBox>& a_velocity,
                                     const bool                a_fourth_order) const
{
   if ( !a_fourth_order ) {
      
      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxes( a_flux, a_faceDist, a_velocity, false );
   }
   else {
      
      const DisjointBoxLayout& grids = a_flux.disjointBoxLayout();
      LevelData<FluxBox> fourth_order_flux(grids, SpaceDim, IntVect::Zero);
      
      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxes( fourth_order_flux, a_faceDist, a_velocity, true );
      
      // Compute the second-order flux in valid plus ghost cell faces,
      // then overwrite with the fourth-order flux on the valid faces.
      CH_assert(a_flux.ghostVect() == IntVect::Unit);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_flux[dit].copy(a_velocity[dit]);
         
         Box box = grow(grids[dit],1);
         for (int dir=0; dir<SpaceDim; ++dir) {
            a_flux[dit].mult(a_faceDist[dit], box, 0, dir, 1);
         }
         
         a_flux[dit].copy(fourth_order_flux[dit], grids[dit]);
      }
   }
   
   a_flux.exchange();
}

void TwoFieldNeutralsOp::computeParallelPressureGradient(LevelData<FArrayBox>&       a_par_grad,
                                                         const LevelData<FArrayBox>& a_soln_density,
                                                         const Real                  a_mass  ) const
{
   /* Computes the parallel projection
    of the parallel pressure gradient,
    Fpar = grad_par(P) = Bpol/B * grad_pol(P)
    P= n_g * T_g = (rho_g/m_g) * T_g
    
    NB: need to generalize to 3D
    */

   // Get boxes
   const DisjointBoxLayout& grids( a_par_grad.getBoxes() );
   
   int order = (m_geometry.secondOrder()) ? 2 : 4;
   
   //Get neutrals pressure, P
   LevelData<FArrayBox> pressure(grids, 1, a_soln_density.ghostVect());
   for (DataIterator dit( pressure.dataIterator() ); dit.ok(); ++dit) {
      pressure[dit].copy(a_soln_density[dit]);
      pressure[dit].mult(m_Tg[dit]);
      pressure[dit].divide(a_mass);
   }

   // Fill ghosts if we have no data for m_Tg or do not want to use
   // density BCs (e.g., when we have recycling BCs ON)
   if (m_ghostVect == IntVect::Zero || m_extrapolate_density) {
      bool fourth_order = (m_geometry.secondOrder()) ? false : true;
      m_geometry.extrapolateToPhysicalGhosts(pressure, fourth_order);
   }
   
   // Compute 2D (R,Z) mapped components of grad(P)
   LevelData<FArrayBox> pressure_gradient_mapped(grids, CFG_DIM, IntVect::Unit);
   m_geometry.computeMappedPoloidalGradientWithGhosts(pressure, pressure_gradient_mapped, order);
   
   // Since the ghost cells of pressure_gradient_mapped are only 2nd-order
   // accurate, need to do this for 4th order implementation
   if (order == 4) m_geometry.fillInternalGhosts(pressure_gradient_mapped);
   
   // Compute 2D (R,Z) physical components of grad(P)
   LevelData<FArrayBox> pressure_gradient(grids, CFG_DIM, IntVect::Unit);
   m_geometry.unmapPoloidalGradient(pressure_gradient_mapped, pressure_gradient);

   
   
   // Get unit b vector
   const LevelData<FArrayBox>& b_dir = m_geometry.getCCBFieldDir();
   
   // Compute grad_parallel(P)
   LevelData<FArrayBox> parallel_grad(grids, 1, IntVect::Unit);
   for (DataIterator dit( parallel_grad.dataIterator() ); dit.ok(); ++dit) {
      FORT_COMPUTE_B_PROJECTION(CHF_BOX(parallel_grad[dit].box()),
                                CHF_CONST_FRA(pressure_gradient[dit]),
                                CHF_CONST_FRA(b_dir[dit]),
                                CHF_FRA1(parallel_grad[dit],0));

      if (order == 4) {
         const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         const ProblemDomain& domain = coord_sys.domain();
         fourthOrderAverageCell(parallel_grad[dit], domain, grids[dit]);
      }
   }
   
   if (order == 4) m_geometry.fillInternalGhosts(parallel_grad);
   
   //Compute parallel pressure force, Fpar = -Bpol/B * grad_pol(P)
   for (DataIterator dit( a_par_grad.dataIterator() ); dit.ok(); ++dit) {
      a_par_grad[dit].copy(parallel_grad[dit]);
   }
}

void TwoFieldNeutralsOp::computePerpPressureGradient(LevelData<FluxBox>&           a_perp_grad,
                                                     const LevelData<FArrayBox>&   a_soln_density,
                                                     const Real                    a_mass) const
{
   /* Computes the perpendicular pressure gradient,
    Fperp = grad_perp(P) = (grad(P) - grad_par(P))
    P= n_g * T_g = (rho_g/m_g) * T_g

    NB: need to generalize to 3D
    */
   
   // Get boxes
   const DisjointBoxLayout& grids( a_perp_grad.getBoxes() );

   int order = (m_geometry.secondOrder()) ? 2 : 4;
   
   //Get neutrals pressure, P
   LevelData<FArrayBox> pressure(grids, 1, a_soln_density.ghostVect());
   for (DataIterator dit( pressure.dataIterator() ); dit.ok(); ++dit) {
      pressure[dit].copy(a_soln_density[dit]);
      pressure[dit].mult(m_Tg[dit]);
      pressure[dit].divide(a_mass);
   }
   
   // Fill ghosts if we have no data for m_Tg or do not want to use
   // density BCs (e.g., when we have recycling BCs ON)
   if (m_ghostVect == IntVect::Zero || m_extrapolate_density) {
      bool fourth_order = (m_geometry.secondOrder()) ? false : true;
      m_geometry.extrapolateToPhysicalGhosts(pressure, fourth_order);
   }
   
   //Compute 2D (R,Z) mapped components of grad(P)
   LevelData<FluxBox> pressure_gradient_mapped(grids, CFG_DIM, IntVect::Unit);
   m_geometry.computeMappedPoloidalGradientWithGhosts(pressure, pressure_gradient_mapped, order);
   
   //Compute 2D (R,Z) physical components of grad(P)
   LevelData<FluxBox> pressure_gradient(grids, CFG_DIM, IntVect::Unit);
   m_geometry.unmapPoloidalGradient(pressure_gradient_mapped, pressure_gradient);
   
   //Compute physical components of grad_parallel(P)
   LevelData<FluxBox> pressure_grad_par(grids, CFG_DIM, IntVect::Unit);
   for (DataIterator dit( pressure.dataIterator() ); dit.ok(); ++dit) {
      pressure_grad_par[dit].copy(pressure_gradient[dit]);
   }
   m_geometry.projectOntoParallel(pressure_grad_par);
   
   //Subtract parallel grad from the full gradient
   for (DataIterator dit( pressure.dataIterator() ); dit.ok(); ++dit) {
      a_perp_grad[dit].copy(pressure_gradient[dit]);
      a_perp_grad[dit] -= pressure_grad_par[dit];
   }
}

void TwoFieldNeutralsOp::computeDiffusionCoefficients(LevelData<FluxBox>& a_ellip_coeff,
                                                      LevelData<FluxBox>& a_ellip_coeff_mapped,
                                                      const Real          a_mass,
                                                      const Real          a_time) const
{
   /*
    Computes the perpendicular diffusivity coefficients
    D_perp = 1/(mg * ne * <sigma_cx*Vti>)
   */
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // Get perpendicular tensor coefficeints
   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   
   // Get scalar factors
   LevelData<FArrayBox> D_cell(grids, 1, 2*IntVect::Unit);
   LevelData<FluxBox> D_face(grids, 1, 2*IntVect::Unit);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      D_cell[dit].copy(m_chx_rate[dit]);
      D_cell[dit].mult(a_mass);
      D_cell[dit].invert(1.0);
   }
   
   if (m_ghostVect == IntVect::Zero) {
      bool fourth_order = (m_geometry.secondOrder()) ? false : true;
      m_geometry.extrapolateToPhysicalGhosts(D_cell, fourth_order);
   }
   convertCellToFace(D_face, D_cell);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_ellip_coeff[dit].copy(perp_coeff[dit]);
      a_ellip_coeff_mapped[dit].copy(perp_coeff_mapped[dit]);
      for (int dir = 0; dir < SpaceDim; dir++) {
         for (int n = 0; n < SpaceDim*SpaceDim; n++) {
            a_ellip_coeff[dit][dir].mult(D_face[dit][dir],0,n,1);
            a_ellip_coeff_mapped[dit][dir].mult(D_face[dit][dir],0,n,1);
         }
      }
   }
}

void TwoFieldNeutralsOp::computeViscosityCoefficients(LevelData<FluxBox>&         a_ellip_coeff,
                                                      LevelData<FluxBox>&         a_ellip_coeff_mapped,
                                                      const LevelData<FArrayBox>& a_density,
                                                      const Real                  a_mass,
                                                      const Real                  a_time) const
{
   /*
    Computes the isotropic viscosity coefficients
    D_isotropic = T_g/(mg * ne * <sigma_cx*Vti>)
    */
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   //Parallel and perpendicular tensor coefficeints
   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& par_coeff = m_geometry.getEllipticOpParCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();
   
   // Get scalar factors
   LevelData<FArrayBox> D_cell(grids, 1, 2*IntVect::Unit);
   LevelData<FluxBox> D_face(grids, 1, 2*IntVect::Unit);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      D_cell[dit].copy(m_Tg[dit]);
      D_cell[dit].mult(a_density[dit]);
      D_cell[dit].divide(m_chx_rate[dit]);
      D_cell[dit].divide(a_mass);
   }
   
   // Get the viscosity coefficient on faces
   if (m_ghostVect == IntVect::Zero || m_extrapolate_density) {
      bool fourth_order = (m_geometry.secondOrder()) ? false : true;
      m_geometry.extrapolateToPhysicalGhosts(D_cell, fourth_order);
   }
   convertCellToFace(D_face, D_cell);
      
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_ellip_coeff[dit].copy(perp_coeff[dit]);
      a_ellip_coeff[dit] += par_coeff[dit];
      a_ellip_coeff_mapped[dit].copy(perp_coeff_mapped[dit]);
      a_ellip_coeff_mapped[dit] += par_coeff_mapped[dit];
      for (int dir = 0; dir < SpaceDim; dir++) {
         for (int n = 0; n < SpaceDim*SpaceDim; n++) {
            a_ellip_coeff[dit][dir].mult(D_face[dit][dir],0,n,1);
            a_ellip_coeff_mapped[dit][dir].mult(D_face[dit][dir],0,n,1);
         }
      }
   }
}

void TwoFieldNeutralsOp::computeAtomicRates(LevelData<FArrayBox>&   a_ionization_rate,
                                            LevelData<FArrayBox>&   a_recombination_rate,
                                            LevelData<FArrayBox>&   a_chx_rate,
                                            const Real              a_time) const
{
   /*
      This function computes ionization and charge-exchange
      frequencies (in COGENT units)
    */
   
   computeIonizationRate(a_ionization_rate, a_time);
   computeChargeExchangeRate(a_chx_rate, a_time);
   
   //Presently set the recombination rate to zero
   setZero(a_recombination_rate);
}

void TwoFieldNeutralsOp::preOpEval(const PS::KineticSpeciesPtrVect&   a_kinetic_species,
                                   const FluidSpeciesPtrVect&         a_fluid_species,
                                   const PS::ScalarPtrVect&           a_scalars,
                                   const EField&                      a_E_field,
                                   const double                       a_time )
{
   CH_TIME("TwoFieldNeutralsOp::preOpEval");
   
   const DisjointBoxLayout& grids = m_geometry.grids();
   
   // The moment calculations used below compute moment including the ghost-cell regions
   // of the kinent distritbuion function, and some of them (e.g., temperature) expect
   // (through the use of intermidiate variables) that cfg objects have the same number of ghosts.
   // This can be fixed by using direct functions from momentOp class if needed for performance.

   IntVect ghostVectMom(IntVect::Zero);
   if (m_consistent_density || m_consistent_temperature || m_consistent_velocity) {
     const PS::KineticSpecies& master_species( *(a_kinetic_species[0]) );
     const PS::PhaseGeom& phase_geom = master_species.phaseSpaceGeometry();
     const PS::LevelData<PS::FArrayBox>& master_dfn( master_species.distributionFunction() );
     PS::IntVect phase_ghosts = master_dfn.ghostVect();
     ghostVectMom = phase_geom.config_restrict(phase_ghosts);
   }
    
   // Compute electron density (by summing charge density from all speceis)
   if (m_consistent_density) {
      if (!m_ne.isDefined()) m_ne.define(grids, 1, ghostVectMom);
      computeIonChargeDensity( m_ne, a_kinetic_species );
   }
   
   //Compute temperature of the zero species (assume it to be the main species)
   if (m_consistent_temperature) {
      if (!m_Tg.isDefined()) m_Tg.define(grids, 1, ghostVectMom);
      const PS::KineticSpecies& this_species( *(a_kinetic_species[0]) );
      this_species.temperature( m_Tg );
   }
   
   // Compute parallel velocity of the zero species (assume it to be the main species)
   if (m_consistent_velocity) {
      if (!m_ion_parallel_vel.isDefined()) {
         m_ion_parallel_vel.define(grids, 1, ghostVectMom);
      }
      const PS::KineticSpecies& this_species( *(a_kinetic_species[0]) );
      this_species.ParallelVelocity( m_ion_parallel_vel );
   }
   
   // Compute normal particle flux thourgh a cell face (NT * Gamma_phys * dA_mapped )
   if (m_dens_recycling_bc && !m_fixed_recycling) {
      if (!m_ion_normal_flux.isDefined()) m_ion_normal_flux.define(grids, 1, IntVect::Zero);
      computeIntegratedNormalIonParticleFlux( m_ion_normal_flux, a_kinetic_species, a_E_field, a_time );
   }
   
   if (m_dens_recycling_bc && m_fixed_recycling) {
      if (!m_ion_normal_flux.isDefined()) m_ion_normal_flux.define(grids, 1, IntVect::Zero);
      computeIntegratedNormalPrescribedFlux( m_ion_normal_flux, m_ne, m_ion_parallel_vel);
   }
   
   // Compute parallel velocity on faces (to be used in recycling BC V_g_par = R * V_i_par)
   if (m_mom_recycling_bc) {
      bool fourth_order = (m_geometry.secondOrder()) ? false : true;
      if (!m_ion_parallel_vel_face.isDefined()) {
	m_ion_parallel_vel_face.define(grids, 1, IntVect::Zero);
      }
      LevelData<FArrayBox> tmp(grids, 1, 2*IntVect::Unit);
      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         tmp[dit].copy(m_ion_parallel_vel[dit]);
      }
      m_geometry.extrapolateToPhysicalGhosts(tmp, fourth_order);
      convertCellToFace(m_ion_parallel_vel_face, tmp);
   }
}

void
TwoFieldNeutralsOp::computeIntegratedNormalPrescribedFlux(LevelData<FluxBox>&         flux,
                                                          const LevelData<FArrayBox>& ne,
                                                          const LevelData<FArrayBox>& Vi) const
{
   MayDay::Error( "TwoFieldNeutralsOp::computeIntegratedNormalPrescribedFlux: not implemented!!! " );
}

void TwoFieldNeutralsOp::defineBlockPC(std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>& a_pc,
                                       std::vector<PS::DOFList>&                                    a_dof_list,
                                       const PS::ODEVector&                                         a_soln_vec,
                                       void*                                                        a_gkops,
                                       const std::string&                                           a_out_string,
                                       const std::string&                                           a_opt_string,
                                       bool                                                         a_im,
                                       const FluidSpecies&                                          a_fluid_species,
                                       const PS::GlobalDOFFluidSpecies&                             a_global_dofs,
                                       const int                                                    a_species_idx,
                                       const int                                                    a_id )
{
  if (a_im && m_is_time_implicit) {
    CH_assert(a_pc.size() == a_dof_list.size());
  
    if (!procID()) {
      std::cout << "  Fluid Species " << a_species_idx
                << " : "
                << " creating " << _FLUID_OP_PC_ << " preconditioner"
                << " (index = " << a_pc.size() << ").\n";
    }
  
    PS::Preconditioner<PS::ODEVector,PS::AppCtxt> *pc;
    pc = new PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt>;
    pc->setSysID(a_id);
    dynamic_cast<PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt>*>
      (pc)->define(a_soln_vec, a_gkops, *this, m_opt_string, m_opt_string, a_im);
    dynamic_cast<PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt>*>
      (pc)->speciesIndex(a_species_idx);
  
    PS::DOFList dof_list(0);
  
    for (int var=0; var<a_fluid_species.num_cell_vars(); ++var) {

       const LevelData<FArrayBox>& soln_data   (a_fluid_species.cell_var(var));
       const DisjointBoxLayout&    grids       (soln_data.disjointBoxLayout());
       const int                   n_comp      (soln_data.nComp());
       const LevelData<FArrayBox>& pMapping    (a_global_dofs.data(var)); 
  
       for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
          const Box& grid = grids[dit];
          const FArrayBox& pMap = pMapping[dit];
          for (BoxIterator bit(grid); bit.ok(); ++bit) {
             IntVect ic = bit();
             for (int n(0); n < n_comp; n++) {
                dof_list.push_back((int) pMap.get(ic ,n) - a_global_dofs.mpiOffset());
             }
          }
       }
    }
  
    m_my_pc_idx = a_pc.size();

    a_pc.push_back(pc);
    a_dof_list.push_back(dof_list);
  }

  return;
}


void TwoFieldNeutralsOp::updateBlockPC(std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>& a_pc,
                                       const PS::KineticSpeciesPtrVect&                             a_kin_species_phys,
                                       const FluidSpeciesPtrVect&                                   a_fluid_species,
                                       const Real                                                   a_time,
                                       const int                                                    a_step,
                                       const int                                                    a_stage,
                                       const Real                                                   a_shift,
                                       const bool                                                   a_im,
                                       const int                                                    a_species_idx )
{
  if (a_im && m_is_time_implicit) {
    CH_assert(m_my_pc_idx >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx);
  
    if (!procID()) {
      std::cout << "    ==> Updating " << _FLUID_OP_PC_ << " preconditioner " 
                << " for TwoFieldNeutralsOp of fluid species " << a_species_idx << ".\n";
    }
  
    PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt> *pc 
      = dynamic_cast<PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt>*>(a_pc[m_my_pc_idx]);
    CH_assert(pc != NULL);
    pc->update( a_kin_species_phys, 
                a_fluid_species, 
                a_time, 
                a_step,
                a_stage,
                a_shift, 
                a_im, 
                a_species_idx);
  }
  return;
}


void TwoFieldNeutralsOp::updatePCImEx(const FluidSpeciesPtrVect&       a_fluid_species,
                                      const PS::KineticSpeciesPtrVect& a_kinetic_species,
                                      const double                     a_time,
                                      const int                        a_step,
                                      const int                        a_stage,
                                      const double                     a_shift,
                                      const int                        a_component)
{
   CH_TIME("TwoFieldNeutralsOp::updatePCImEx");

   // Get boxes
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // Get computational fluid species data
   const FluidSpecies& fluid_species = static_cast<const FluidSpecies&>(*a_fluid_species[a_component]);
   const LevelData<FArrayBox>& density( fluid_species.cell_var("density") );
   
   // Get physical density
   LevelData<FArrayBox> density_phys(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      density_phys[dit].copy(density[dit]);
   }
   m_geometry.divideJonValid(density_phys);
   
   //Get species mass
   const Real species_mass = fluid_species.mass();
   
   // Compute linear coefficients (beta) and update density elliptic solver
   LevelData<FArrayBox> beta_diffusion(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      beta_diffusion[dit].setVal(a_shift);
      beta_diffusion[dit].plus(m_ionization_rate[dit]);
      beta_diffusion[dit].minus(m_recombination_rate[dit]);
      
      FORT_ADD_BCKGR_SRC_BETA(CHF_BOX(beta_diffusion[dit].box()),
                              CHF_CONST_FRA1(m_ionization_rate[dit],0),
                              CHF_CONST_FRA1(density_phys[dit],0),
                              CHF_CONST_REAL(m_bckgr_density),
                              CHF_CONST_REAL(species_mass),
                              CHF_FRA1(beta_diffusion[dit],0));
      
      beta_diffusion[dit].divide(m_Tg[dit]);

   }
   if (!m_analytic_pc) {
      m_diffusion_op->updateImExPreconditioner( beta_diffusion, *m_diffusion_bc );
   }
   
   // Compute linear coefficients (beta) and update velocity elliptic solver
   LevelData<FArrayBox> beta_viscosity(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      beta_viscosity[dit].setVal(a_shift);
      beta_viscosity[dit].plus(m_ionization_rate[dit]);
      beta_viscosity[dit].plus(m_chx_rate[dit]);
      beta_viscosity[dit].mult(density_phys[dit]);
   }
   
   if (!m_analytic_pc) {
      m_viscosity_op->updateImExPreconditioner( beta_viscosity, *m_viscosity_bc );
   }
   
   // Save density for the use in solvePCImEx function
   if ( !m_pc_density.isDefined() ) {
      m_pc_density.define(grids, 1, IntVect::Zero);
   }
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      m_pc_density[dit].copy(density_phys[dit]);
   }

   // If using the analytic preconditioner, we update here instead of the
   // two previous calls to m_diffusion_op and m_viscosity_op
   if (m_analytic_pc) {
      LevelData<FArrayBox> zero_beta(grids, 1, IntVect::Zero);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
	 zero_beta[dit].setVal(0.);
      }
      m_twofieldneutrals_pc->updateSourceTermComponents(beta_diffusion,
							zero_beta,
						        beta_viscosity);
      m_twofieldneutrals_pc->updatePreconditioner(*m_diffusion_bc, *m_viscosity_bc);
   }
}


void TwoFieldNeutralsOp::solvePCImEx(FluidSpeciesPtrVect&              a_fluid_species_solution,
                                     const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                     const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                     const int                         a_component )
{
   /*
    This function computes solution to the preconditioned problem.
    The block-diagonal preconditioner involves all difusion and linear source/sink terms.
   */

   CH_TIME("TwoFieldNeutralsOp::solvePCImEx");
   if (!m_is_time_implicit) {
     MayDay::Error("Error in TwoFieldNeutralsOp::solvePCImEx(): m_is_time_implicit is false. This function shouldn't have been called at all!");
   }
   
   // Get RHS data
   const FluidSpecies& rhs_species = static_cast<const FluidSpecies&>(*a_fluid_species_rhs[a_component]);
   const LevelData<FArrayBox>& rhs_density( rhs_species.cell_var("density") );
   const LevelData<FArrayBox>& rhs_par_mom( rhs_species.cell_var("momentumDensity_virtual") );
   
   // Get boxes
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // Convert the rhs to physical space
   LevelData<FArrayBox> rhs_density_phys(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> rhs_par_mom_phys(grids, 1, IntVect::Zero);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      rhs_density_phys[dit].copy(rhs_density[dit]);
      rhs_par_mom_phys[dit].copy(rhs_par_mom[dit]);
   }

   m_geometry.divideJonValid(rhs_density_phys);
   m_geometry.divideJonValid(rhs_par_mom_phys);
   
   // Get state vector components, z_density = J * rho_g; z_par_mom = J * rho_g * V_par_g
   FluidSpecies& sol_species = static_cast<FluidSpecies&>(*a_fluid_species_solution[a_component]);
   LevelData<FArrayBox>& z_density = sol_species.cell_var("density");
   LevelData<FArrayBox>& z_par_mom = sol_species.cell_var("momentumDensity_virtual");

   // Get convergence parameters for diffusion MBSolver
   double tol = m_diffusion_op->getPrecondTol();
   double max_iter = m_diffusion_op->getPrecondMaxIter();
   double precond_tol = 0.;
   double precond_max_iter = 1;

   if (m_analytic_pc) {
      m_twofieldneutrals_pc->setImExPreconditionerConvergenceParams(tol, max_iter, precond_tol, precond_max_iter);
      m_twofieldneutrals_pc->solvePreconditioner(rhs_density_phys, rhs_par_mom_phys, z_density, z_par_mom);
   }
   else
   {
      // Solve for z_density = rho_g * T_g
      m_diffusion_op->setImExPreconditionerConvergenceParams(tol, max_iter, precond_tol, precond_max_iter);
      m_diffusion_op->solveImExPreconditioner(rhs_density_phys, z_density);

      // Get convergence parameters for viscosity MBSolver
      tol = m_viscosity_op->getPrecondTol();
      max_iter = m_viscosity_op->getPrecondMaxIter();

      // Solve for z_par_mom = par_mom / rho_g
      m_viscosity_op->setImExPreconditionerConvergenceParams(tol, max_iter, precond_tol, precond_max_iter);
      m_viscosity_op->solveImExPreconditioner(rhs_par_mom_phys, z_par_mom);

   }


   // Convert the elliptic solutions to the state vector components
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      z_density[dit].divide(m_Tg[dit]);
      z_par_mom[dit].mult(m_pc_density[dit]);
   }
   m_geometry.multJonValid(z_density);
   m_geometry.multJonValid(z_par_mom);

/*
   // Use this for an identity preconditioner
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     z_par_mom[dit].copy(rhs_par_mom[dit]);
     z_density[dit].copy(rhs_density[dit]);
   }
*/

}

void TwoFieldNeutralsOp::fillGhostCells(FluidSpecies&  a_species_phys,
                                        const double   a_time )
{
   for (int n=0; n<a_species_phys.num_cell_vars(); ++n) {
      // Fill ghost cells except for those on physical boundaries
      LevelData<FArrayBox>& fld( a_species_phys.cell_var(n) );
      m_geometry.fillInternalGhosts( fld );
      // Fill ghost cells on physical boundaries
      m_fluid_bc.at(n)->apply( a_species_phys, a_time );
   }
}

void TwoFieldNeutralsOp::getPlasmaParameters(const Real a_time)
{
   // Get prescribed electron density (m_ne)
   // electron temperature (m_Te) and neutrals temperature (m_Tg)
   OneFieldNeutralsOp::getPlasmaParameters(a_time);
   
   const DisjointBoxLayout& grids = m_geometry.grids();
   
   // Get prescribed ion parallel velocity
   if (!m_ion_parallel_vel.isDefined()) {
      IntVect ghostVect = IntVect::Zero;
      if (m_mom_recycling_bc) {
         int ghosts = (m_geometry.secondOrder()) ? 1 : 2;
         ghostVect = ghosts*IntVect::Unit;
      }
      m_ion_parallel_vel.define(grids, 1, ghostVect);
   }
   
   if (!m_consistent_velocity && m_first_call) {
      m_ion_parallel_vel_func->assign( m_ion_parallel_vel, m_geometry, a_time);
   }

   if (m_dens_recycling_bc && m_fixed_recycling && m_first_call) {
      if (!m_ion_normal_flux.isDefined()) m_ion_normal_flux.define(grids, 1, IntVect::Zero);
      computeIntegratedNormalPrescribedFlux( m_ion_normal_flux, m_ne, m_ion_parallel_vel);
   }   
}

void TwoFieldNeutralsOp::parseParameters( ParmParse& a_pp )
{
   
   //NB: other important parameters are parsed
   //when OneFieldNeutralsOp base class is initialized
   
   a_pp.query( "time_implicit", m_is_time_implicit);

   a_pp.query( "use_analytic_pc", m_analytic_pc);
   
   if (a_pp.contains("advScheme")) {
      a_pp.get("advScheme", m_advScheme );
   }

   if (a_pp.contains("background_density")) {
      a_pp.get("background_density", m_bckgr_density );
   }

   GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
   std::string grid_function_name;
   
   a_pp.query("fixed_density_recycling", m_fixed_recycling);
   
   if (a_pp.contains("ion_velocity")) {
      a_pp.get( "ion_velocity", grid_function_name );
      m_ion_parallel_vel_func = grid_library->find( grid_function_name );
      m_consistent_velocity = false;
   }
   else {
      m_consistent_velocity = true;
   }
   
   //Check input for consistency
   if (!m_consistent_temperature && m_Tg_func == NULL) {
      MayDay::Error("Neutral temperature must be specified or set consistent");
   }
}

void TwoFieldNeutralsOp::getMemberVar(LevelData<FArrayBox>&        a_Var,
                                      const FluidSpecies&          a_fluid_species,
                                      const string&                a_name ) const
{
   const DisjointBoxLayout& grids( m_geometry.grids() );
   const IntVect ghostVect = IntVect::Zero;
   
   if(a_name=="parallel_velocity") {
      a_Var.define(grids, 1, ghostVect);
      a_fluid_species.velocity_virtual(a_Var);
   }
   else if(a_name=="ion_parallel_velocity") {
      a_Var.define(grids, 1, IntVect::Zero);
      // Skip first call as this function is being called before
      // m_ion_parallel_vel is defined
      if (!m_first_call) {
         for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
            a_Var[dit].copy(m_ion_parallel_vel[dit]);
         }
      }
   }
   else {
      cout << "Trying to getMemberVar() from TwoFieldNeutralsOp for a_name = " <<
              a_name << endl;
      MayDay::Error("Trying to getMemberVar() from TwoFieldNeutralsOp!!!  ");
   }
}

void TwoFieldNeutralsOp::printParameters()
{
   if (procID()==0) {
      std::cout << "TwoFieldNeutralsOp parameters:" << std::endl;
      std::cout << "  time_implicit  =  " << m_is_time_implicit << std::endl;
      std::cout << "  advScheme  =  " << m_advScheme << std::endl;
   }
}

Real TwoFieldNeutralsOp::computeDt( const FluidSpeciesPtrVect&  a_fluid_species )
{
   return DBL_MAX;
   std::cout << " TwoFieldNeutralsOp::computeDt is being called " << std::endl;

   // Need stable step based on Courant for advection
   // Need stable step for explicit diffusion if using explicit
   // Need stable step for advection-diffusion, which is diff for C2 adv vs U1 adv
   // See Pozrikidis pg 578
   // dtmaxCour = 0.5*dX*dY/(Umax_x*dZ + Umax_z*dX)
   // dtmaxDiff = 0.5/D*dX^2*dY^2/(dX^2 + dZ^2)
   // dtmaxAD_C2 = 2.0*D/(Umax_x^2 + Umax_z^2)
   // dtmaxAD_U1 = 2.0*D/dX^2 + 2*D/dZ^2 + Umax_x/dX + Umax_z/dZ
}


#include "NamespaceFooter.H"
