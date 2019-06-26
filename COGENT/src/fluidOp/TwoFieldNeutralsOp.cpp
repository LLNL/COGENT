#include "TwoFieldNeutralsOp.H"
#include "EllipticOpBCFactory.H"
#include "FluidVarBCFactory.H"
#include "FourthOrderUtil.H"
#include "SpaceUtils.H"
#include "Directions.H"
#include "mappedAdvectionFlux.H"
#include "FourthOrderUtil.H"
#include "inspect.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "GKOps.H"
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
                                       const int       a_verbosity )
   : m_verbosity(a_verbosity),
     m_geometry(a_geometry),
     m_is_time_implicit(true),
     m_advScheme("uw3"),
     m_bckgr_density(1.0e-7),
     m_opt_string(a_pp_str),
     m_my_pc_idx(-1)
{
   ParmParse pp(a_pp_str.c_str());
   parseParameters( pp );
   if (m_verbosity>0) {
      printParameters();
   }

   std::string name;
   std::string prefix;
   std::string prefixBC;
   EllipticOpBCFactory elliptic_op_bc_factory;

   // Construct diffusion elliptic operator
   name = "neutrals_diffusion";

   prefix = name ;
   ParmParse ppdiff( prefix.c_str() );
   m_diffusion_op = new Diffusion(ppdiff, a_geometry);

   prefixBC = "BC." + name ;
   ParmParse ppdiff_bc( prefixBC.c_str() );
   m_diffusion_bc = elliptic_op_bc_factory.create(name,
                                                  ppdiff_bc,
                                                  a_geometry.getCoordSys()->type(),
                                                  false );

   // Construct viscosity elliptic operator
   name = "neutrals_viscosity";

   prefix = name;
   ParmParse ppvisc( prefix.c_str() );
   m_viscosity_op = new Diffusion(ppvisc, a_geometry);

   prefixBC = "BC." + name ;
   ParmParse ppvisc_bc( prefixBC.c_str() );
   m_viscosity_bc = elliptic_op_bc_factory.create(name,
                                                  ppvisc_bc,
                                                  a_geometry.getCoordSys()->type(),
                                                  false );
   
   // Get initial conditions for cell variables and define boundary conditions instance
   const std::string species_name( "IC." + a_species_name);
   std:: string variable_name;
   ParmParse ppic(  species_name.c_str() );
   FluidVarBCFactory fluid_var_bc_factory;
   
   if (ppic.contains("density.function")) {
      variable_name = "density";
      parseInitialConditions(a_species_name, variable_name);
      m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, a_geometry.getCoordSys()->type(), false) );
   }
   
   if (ppic.contains("momentumDensity_virtual.function")) {
      variable_name = "momentumDensity_virtual";
      parseInitialConditions(a_species_name, variable_name);
      m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, a_geometry.getCoordSys()->type(), false) );
   }

   //Get disjoint box layout
   const DisjointBoxLayout& grids( a_geometry.grids() );
   
   // Set ghost vect for the initialization objects
   IntVect ghost_vec(4*IntVect::Unit);

   // Initialize neutral temperature
   m_neutral_temperature.define(grids, 1, ghost_vec);
   m_ntr_temp->assign( m_neutral_temperature, m_geometry, 0.);
   
   // Initialize atomic data
   m_ionization_rate.define(grids, 1, ghost_vec);
   m_recombination_rate.define(grids, 1, ghost_vec);
   m_chx_rate.define(grids, 1, ghost_vec);
   computeAtomicRates(m_ionization_rate, m_recombination_rate, m_chx_rate, 0.);

   // Fill internal ghosts: needed if the mapping in ghosts
   // is possibly inaccurate (e.g., flux-aligned SN)
   // and we want to overwrite assign results with extraBlockExchange
   m_geometry.fillInternalGhosts( m_neutral_temperature );
   m_geometry.fillInternalGhosts( m_ionization_rate );
   m_geometry.fillInternalGhosts( m_recombination_rate );
   m_geometry.fillInternalGhosts( m_chx_rate );

   // Initialize ion parallel velocity
   if (m_prescribed_plasma_bckrg) {
      m_ion_parallel_vel.define(grids, 1, ghost_vec);
      m_ion_velocity->assign( m_ion_parallel_vel, m_geometry, 0.);
      m_geometry.fillInternalGhosts( m_ion_parallel_vel );
   }
}


TwoFieldNeutralsOp::~TwoFieldNeutralsOp()
{

   for (std::vector<FluidVarBC*>::iterator it = m_fluid_bc.begin(); it != m_fluid_bc.end(); ++it)
   {
      delete (*it);
   }
   m_fluid_bc.clear();
   
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
    density_rhs = \nabla_perp D \nabla_perp (rho_g * T_g), D = 1/(rho_g * ne * (<sigma_chx*Vti> + <sigma_iz*Vte>) )
    par_mom_rhs = \nabla \eta \nabla (v_par_g), \eta =rho_g * T_g/(m_g * ne * <sigma_chx*Vti> )
    */
   
   // Get rhs data
   LevelData<FArrayBox>& rhs_density( a_fluid_rhs.cell_var("density") );
   LevelData<FArrayBox>& rhs_par_mom( a_fluid_rhs.cell_var("momentumDensity_virtual") );
   
   // Get soln data
   const LevelData<FArrayBox>& soln_density( a_fluid_soln.cell_var("density") );
   const LevelData<FArrayBox>& soln_par_mom( a_fluid_soln.cell_var("momentumDensity_virtual") );

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
      neutr_pressure[dit].mult(m_neutral_temperature[dit]);
   }
   computeDiffusiveCoefficients(ellip_coeff, ellip_coeff_mapped, species_mass, a_time);
   m_diffusion_op->setOperatorCoefficients(ellip_coeff, ellip_coeff_mapped, *m_diffusion_bc);
   m_diffusion_op->computeFluxDivergence( neutr_pressure, flux_div, false, false);
   
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
   m_viscosity_op->setOperatorCoefficients(ellip_coeff, ellip_coeff_mapped, *m_viscosity_bc);
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
   
   // Get parallel velocity on faces; loose two layers of ghost cells
   LevelData<FluxBox> velocity_on_faces( grids, 1, IntVect::Unit );
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; dir++) {
      
      int order = (m_geometry.secondOrder()) ? 2 : 4;
      
      IntVect grow_vect = 1*IntVect::Unit;
      grow_vect[dir] = 0;
      Box box = grow(grids[dit],grow_vect);
      
      SpaceUtils::faceInterpolate(dir,
                                  surroundingNodes(box,dir),
                                  order,
                                  parallel_vel[dit],
                                  velocity_on_faces[dit][dir] );
      }

      // Force zero flux on radial faces for the flux-aligned grids
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
   computeAdvectionRHS(rho_rhs, soln_density, velocity);

   // Compute parallel advection RHS in eqn for parallel momentum
   LevelData<FArrayBox> par_mom_rhs( grids, 1, IntVect::Zero );
   computeAdvectionRHS(par_mom_rhs, soln_par_mom, velocity);
   
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
    v_perp_g = - (1/(m_i * rho_g * ne * (<sigma_chx*Vti> + <sigma_iz*Vte>))) * grad_perp (rho_g * T_g)
    */
   
   // Get rhs data
   LevelData<FArrayBox>& rhs_par_mom( a_fluid_rhs.cell_var("momentumDensity_virtual") );
   
   // Get soln data
   const LevelData<FArrayBox>& soln_density( a_fluid_soln.cell_var("density") );
   const LevelData<FArrayBox>& soln_par_mom( a_fluid_soln.cell_var("momentumDensity_virtual") );
   const Real species_mass = a_fluid_soln.mass();
   
  
   // Get boxes
   const DisjointBoxLayout& grids( soln_density.getBoxes() );
   
   // Compute perpendicular velocity (physical vector)
   LevelData<FluxBox> velocity( grids, SpaceDim, IntVect::Unit );
   computePerpPressureGradient(velocity, soln_density, species_mass);

   // Get scalar factors
   LevelData<FArrayBox> fac_cell(grids, 1, 3*IntVect::Unit);
   LevelData<FluxBox> fac_face(grids, 1, IntVect::Unit);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      fac_cell[dit].copy(m_ionization_rate[dit]);
      fac_cell[dit].copy(m_chx_rate[dit]);
      fac_cell[dit].mult(soln_density[dit]);
      fac_cell[dit].invert(-1.0);

      for (int dir=0; dir<SpaceDim; dir++) {
         
         int order = (m_geometry.secondOrder()) ? 2 : 4;
         
         IntVect grow_vect = IntVect::Unit;
         grow_vect[dir] = 0;
         Box box = grow(grids[dit],grow_vect);
         
         SpaceUtils::faceInterpolate(dir,
                                     surroundingNodes(box,dir),
                                     order,
                                     fac_cell[dit],
                                     fac_face[dit][dir] );

         for (int n = 0; n < SpaceDim; n++) {
            velocity[dit][dir].mult(fac_face[dit][dir],0,n,1);
         }
      }
   }
   
   // Apply axisymmetric corrections (for pointwsie metrics does nothing)
   m_geometry.applyAxisymmetricCorrection(velocity);
   
   // Compute parallel advection RHS in eqn for parallel momentum
   LevelData<FArrayBox> par_mom_rhs( grids, 1, IntVect::Zero );
   computeAdvectionRHS(par_mom_rhs, soln_par_mom, velocity);
   
   // Add parallel advection to RHS of the entire system
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

void TwoFieldNeutralsOp::computeAdvectionRHS(LevelData<FArrayBox>& a_div,
                                             const LevelData<FArrayBox>& a_dist,
                                             const LevelData<FluxBox>& a_velocity) const
{
   
   // Get boxes
   const DisjointBoxLayout& grids( a_div.getBoxes() );
   
   // Get mapped components of parallel velocity
   LevelData<FluxBox> normal_vel( grids, 1, IntVect::Unit );
   m_geometry.computeMetricTermProductAverage(normal_vel, a_velocity, false);
   
   // Get distribution on faces
   LevelData<FluxBox> faceDist( grids, 1, IntVect::Unit );
   SpaceUtils::upWindToFaces(faceDist, a_dist, normal_vel, m_advScheme);
   
   // Compute flux
   bool fourth_order = (m_geometry.secondOrder()) ? false : true;
   LevelData<FluxBox> flux( grids, CFG_DIM, IntVect::Unit );
   computeFlux(flux, faceDist, a_velocity, fourth_order);
   
   // Compute mapped divergence
   const bool OMIT_NT(false);
   m_geometry.computeMappedGridDivergence( flux, a_div, OMIT_NT );
   
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
      pressure[dit].mult(m_neutral_temperature[dit]);
      pressure[dit].divide(a_mass);
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
      pressure[dit].mult(m_neutral_temperature[dit]);
      pressure[dit].divide(a_mass);
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

void TwoFieldNeutralsOp::computeDiffusiveCoefficients(LevelData<FluxBox>& a_ellip_coeff,
                                                      LevelData<FluxBox>& a_ellip_coeff_mapped,
                                                      const Real          a_mass,
                                                      const Real          a_time) const
{
   /*
    Computes the perpendicular diffusivity coefficients
    D_perp = 1/(mg * ne * (<sigma_cx*Vti> + <sigma_iz*Vte>))
   */
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // Get perpendicular tensor coefficeints
   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   
   // Get scalar factors
   LevelData<FArrayBox> D_cell(grids, 1, 4*IntVect::Unit);
   LevelData<FluxBox> D_face(grids, 1, 2*IntVect::Unit);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      D_cell[dit].copy(m_ionization_rate[dit]);
      D_cell[dit].copy(m_chx_rate[dit]);
      D_cell[dit].mult(a_mass);
      D_cell[dit].invert(1.0);

      for (int dir=0; dir<SpaceDim; dir++) {
   
         int order = (m_geometry.secondOrder()) ? 2 : 4;
         
         IntVect grow_vect = 2*IntVect::Unit;
         grow_vect[dir] = 0;
         Box box = grow(grids[dit],grow_vect);
         
         SpaceUtils::faceInterpolate(dir,
                                     surroundingNodes(box,dir),
                                     order,
                                     D_cell[dit],
                                     D_face[dit][dir] );
      }
   }

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
   LevelData<FArrayBox> D_cell(grids, 1, 4*IntVect::Unit);
   LevelData<FluxBox> D_face(grids, 1, 2*IntVect::Unit);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      D_cell[dit].copy(m_neutral_temperature[dit]);
      D_cell[dit].mult(a_density[dit]);
      D_cell[dit].divide(m_chx_rate[dit]);
      D_cell[dit].divide(a_mass);
      
      for (int dir=0; dir<SpaceDim; dir++) {
         
         int order = (m_geometry.secondOrder()) ? 2 : 4;
         
         IntVect grow_vect = 2*IntVect::Unit;
         grow_vect[dir] = 0;
         Box box = grow(grids[dit],grow_vect);
         
         SpaceUtils::faceInterpolate(dir,
                                     surroundingNodes(box,dir),
                                     order,
                                     D_cell[dit],
                                     D_face[dit][dir] );
      }
   }
   
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
   /* This function computes ionization and chanrge-exchange
    frequency (in COGENT units) from electron parameters provided
    in the input file (in COGENT units)
    */
   
   //Universal constants (in CGS)
   double mp = 1.6726e-24;
   
   //Get normalization parameters (units)
   double N, T, L;
   ParmParse ppunits( "units" );
   ppunits.get("number_density",N);  //[m^{-3}]
   ppunits.get("temperature",T);     //[eV]
   ppunits.get("length",L);          //[m]
   
   double Tcgs = 1.602e-12 * T; //[erg]
   double Lcgs  = 1.0e2 * L;   //[cm]
   
   double time_norm = Lcgs / sqrt(Tcgs/mp); //[s]
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   IntVect ghost(4*IntVect::Unit);
   
   // Get electron density and temeprature in COGENT units
   LevelData<FArrayBox> ne(grids, 1, ghost);
   m_electron_dens->assign( ne, m_geometry, a_time);
   
   LevelData<FArrayBox> Te(grids, 1, ghost);
   m_electron_temp->assign( Te, m_geometry, a_time);
   
   //Compute ionization frequency
   LevelData<FArrayBox> a2(grids, 1, ghost);
   double fac = pow(T/10.0,2);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a2[dit].copy(Te[dit]);
      a2[dit].mult(Te[dit]);
      a2[dit].mult(fac);
   }
   
   //Get <sigmaV> in m^3/s
   LevelData<FArrayBox> sigmaV(grids, 1, ghost);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      sigmaV[dit].copy(a2[dit]);
      sigmaV[dit].plus(3.0);
      sigmaV[dit].divide(a2[dit]);
      sigmaV[dit].invert(1.0);
      sigmaV[dit].mult(3.0e-14);
   }
   
   //Compute ionization frequency in COGENT time units
   double norm_iz = time_norm * N;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_ionization_rate[dit].copy(sigmaV[dit]);
      a_ionization_rate[dit].mult(ne[dit]);
      a_ionization_rate[dit].mult(norm_iz);
   }
   
   
   //Compute charge-exchange normalization factor
   double norm_chx = time_norm * 1.98e-14 * sqrt(T) * N;
   
   //Get sqrt(Ti); assume T_i=T_g
   LevelData<FArrayBox> sqrtTi(grids, 1, ghost);
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      FORT_COMPUTE_SQRT_FUNC(CHF_BOX(sqrtTi[dit].box()),
                             CHF_CONST_FRA1(m_neutral_temperature[dit],0),
                             CHF_FRA1(sqrtTi[dit],0));
   }
   
   //Compute charge-exchange frequency in COGENT time units
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_chx_rate[dit].copy(ne[dit]);
      a_chx_rate[dit].mult(sqrtTi[dit]);
      a_chx_rate[dit].mult(norm_chx);
   }
   
   //Temporary set recombination rate to 0
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_recombination_rate[dit].setVal(0.);
   }
}

void TwoFieldNeutralsOp::defineBlockPC(std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
                                       std::vector<PS::DOFList>&                                 a_dof_list,
                                       const PS::GKVector&                                       a_soln_vec,
                                       PS::GKOps&                                                a_gkops,
                                       const std::string&                                        a_out_string,
                                       const std::string&                                        a_opt_string,
                                       bool                                                      a_im,
                                       const FluidSpecies&                                       a_fluid_species,
                                       const PS::GlobalDOFFluidSpecies&                          a_global_dofs,
                                       int                                                       a_species_idx )
{
  if (a_im && m_is_time_implicit) {
    CH_assert(a_pc.size() == a_dof_list.size());
  
    if (!procID()) {
      std::cout << "  Fluid Species " << a_species_idx
                << " : "
                << " creating " << _FLUID_OP_PC_ << " preconditioner"
                << " (index = " << a_pc.size() << ").\n";
    }
  
    PS::Preconditioner<PS::GKVector, PS::GKOps> *pc;
    pc = new PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>;
    dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>
      (pc)->define(a_soln_vec, a_gkops, *this, m_opt_string, m_opt_string, a_im);
    dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>
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


void TwoFieldNeutralsOp::updateBlockPC(std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
                                       const PS::KineticSpeciesPtrVect&                          a_kin_species_phys,
                                       const FluidSpeciesPtrVect&                                a_fluid_species,
                                       const Real                                                a_time,
                                       const Real                                                a_shift,
                                       const bool                                                a_im,
                                       const int                                                 a_species_idx )
{
  if (a_im && m_is_time_implicit) {
    CH_assert(m_my_pc_idx >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx);
  
    if (!procID()) {
      std::cout << "    ==> Updating " << _FLUID_OP_PC_ << " preconditioner " 
                << " for TwoFieldNeutralsOp of fluid species " << a_species_idx << ".\n";
    }
  
    PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps> *pc 
      = dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>(a_pc[m_my_pc_idx]);
    CH_assert(pc != NULL);
    pc->update(a_kin_species_phys, a_fluid_species, a_time, a_shift, a_im, a_species_idx);
  }
  return;
}


void TwoFieldNeutralsOp::updatePCImEx(const FluidSpeciesPtrVect&       a_fluid_species,
                                      const PS::KineticSpeciesPtrVect& a_kinetic_species,
                                      const double                     a_time,
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
      
      beta_diffusion[dit].divide(m_neutral_temperature[dit]);

   }
   m_diffusion_op->updateImExPreconditioner( beta_diffusion, *m_diffusion_bc );
   
   // Compute linear coefficients (beta) and update velocity elliptic solver
   LevelData<FArrayBox> beta_viscosity(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      beta_viscosity[dit].setVal(a_shift);
      beta_viscosity[dit].plus(m_ionization_rate[dit]);
      beta_viscosity[dit].plus(m_chx_rate[dit]);
      beta_viscosity[dit].mult(density_phys[dit]);
   }
   m_viscosity_op->updateImExPreconditioner( beta_viscosity, *m_viscosity_bc );
   
   // Save density for the use in solvePCImEx function
   if ( !m_pc_density.isDefined() ) {
      m_pc_density.define(grids, 1, IntVect::Zero);
   }
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      m_pc_density[dit].copy(density_phys[dit]);
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

   // Solve for z_density = rho_g * T_g
   m_diffusion_op->setImExPreconditionerConvergenceParams(tol, max_iter, precond_tol, precond_max_iter);
   m_diffusion_op->solveImExPreconditioner(rhs_density_phys, z_density);

   // Get convergence parameters for viscosity MBSolver
   tol = m_viscosity_op->getPrecondTol();
   max_iter = m_viscosity_op->getPrecondMaxIter();

   // Solve for z_par_mom = par_mom / rho_g
   m_viscosity_op->setImExPreconditionerConvergenceParams(tol, max_iter, precond_tol, precond_max_iter);
   m_viscosity_op->solveImExPreconditioner(rhs_par_mom_phys, z_par_mom);

   // Convert the elliptic solutions to the state vector components
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      z_density[dit].divide(m_neutral_temperature[dit]);
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


void TwoFieldNeutralsOp::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "time_implicit", m_is_time_implicit);

   GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
   std::string grid_function_name;
   
   if (a_pp.contains("prescribed_plasma_bckrg")) {
      a_pp.get("prescribed_plasma_bckrg", m_prescribed_plasma_bckrg );
   }
   else {
      MayDay::Error("Neutral temperature must be specified for the two-field model");
   }

   if (a_pp.contains("neutral_temperature")) {
      a_pp.get("neutral_temperature", grid_function_name );
      m_ntr_temp = grid_library->find( grid_function_name );
   }
   else {
      m_prescribed_plasma_bckrg = false;
   }
   
   if (m_prescribed_plasma_bckrg) {
      
      if (a_pp.contains("electron_density")) {
         a_pp.get( "electron_density", grid_function_name );
         m_electron_dens = grid_library->find( grid_function_name );
      }
      else{
         MayDay::Error("Electron density must be specified for prescribed ionization");
      }
      
      if (a_pp.contains("electron_temperature")) {
         a_pp.get( "electron_temperature", grid_function_name );
         m_electron_temp = grid_library->find( grid_function_name );
      }
      else{
         MayDay::Error("Electron temperature must be specified for prescribed ionization");
      }

      if (a_pp.contains("ion_velocity")) {
         a_pp.get( "ion_velocity", grid_function_name );
         m_ion_velocity = grid_library->find( grid_function_name );
      }
      else{
         MayDay::Error("Electron temperature must be specified for prescribed ionization");
      }
   }


   if (a_pp.contains("advScheme")) {
      a_pp.get("advScheme", m_advScheme );
   }

   if (a_pp.contains("background_density")) {
      a_pp.get("background_density", m_bckgr_density );
   }
 
}

void TwoFieldNeutralsOp::printParameters()
{
   if (procID()==0) {
      std::cout << "TwoFieldNeutralsOp parameters:" << std::endl;
      std::cout << "  time_implicit  =  " << m_is_time_implicit << std::endl;
      std::cout << "  prescribed_plasma_bckgrd  =  " << m_prescribed_plasma_bckrg << std::endl;
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
