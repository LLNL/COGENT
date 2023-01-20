#include <fstream>

#include "GKOps.H"
#include "MomentOp.H"
#include "FourthOrderUtil.H"
#include "Directions.H"
#include "GKVlasovAmpere.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FourthOrderUtil.H"
#include "GKPoisson.H"
#include "PhiSelfConsistentBCOps.H"
#include "FluxSurface.H"
#include "EMVorticityOp.H"
#include "CFGVars.H"
#include "inspect.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"


void GKOps::define( const GKState&            a_state, 
                    const GKSystemParameters& a_gkparams )
{
   CH_assert( isDefined()==false );
   m_phase_geometry = a_state.geometry();
   CH_assert( m_phase_geometry != NULL );

   /* default IMEX partitioning:
    * Explicit: Vlasov, Neutrals
    * Implicit: Collisions, Transport
   */
   m_transport_imex_implicit  = false;
   m_neutrals_imex_implicit   = false;

   setParameters( a_gkparams );

   m_units = new GKUnits();

   ParmParse pp;
   m_boundary_conditions = new GKSystemBC( pp, a_state );
   m_initial_conditions = new GKSystemIC( pp, a_state );

   const KineticSpeciesPtrVect& kin_species(a_state.dataKinetic());
   int n_species = kin_species.size();
   std::vector<std::string> name_list(n_species);
   for (int n=0; n<n_species; n++) {
     const KineticSpecies& species = *(kin_species[n]);
     name_list[n] = species.name();
   }
 
   const double larmor( m_units->larmorNumber() );
   if ( a_state.dataKinetic().size() > 0 ) {
      ParmParse pp_vlasov(GKVlasov::pp_name);
      if ( m_consistent_potential_bcs && m_use_vlasov_ampere_bcs) {
	m_vlasov = new GKVlasovAmpere( pp_vlasov, larmor, name_list, !m_ampere_law );
      }
      else {
         m_vlasov = new GKVlasov( pp_vlasov, larmor, name_list );
      }
   }
   else {
      m_vlasov = NULL;
   }
   m_collisions  = new GKCollisions( m_verbosity );
   m_fluidOp     = new CFG::GKFluidOp( m_phase_geometry->magGeom(), *m_units, m_verbosity );
   m_scalarOp    = new GKScalarOp( *m_phase_geometry, *m_units, m_verbosity );
   if (m_transport_model_on) {
      m_transport = new GKTransport( m_verbosity );
   }
   if (m_neutrals_model_on) {
      m_neutrals = new GKNeutrals( m_verbosity );
   }
 
   m_Y.define(a_state);

   m_state_includes_potential = stateIncludesPotential(a_state);

   const CFG::FluidSpeciesPtrVect fluid_species( a_state.dataFluid() );
   m_trivial_solution_op = m_fluidOp->trivialSolutionOp(fluid_species);

   if ( m_old_vorticity_model ) {
      ParmParse pp( CFG::VorticityOp::pp_name );
      const double larmor( m_units->larmorNumber() );
      if (m_EM_effects) {
         const double ion_skin_depth( m_units->ionSkinNumber() );
         const double electron_skin_depth( m_units->electronSkinNumber() );
         m_VorticityOp = new CFG::EMVorticityBEOp(pp, m_phase_geometry->magGeom(), larmor, ion_skin_depth, electron_skin_depth,  false);
      }
      else {
         m_VorticityOp = new CFG::VorticityBEOp(pp, m_phase_geometry->magGeom(), larmor, false);
      }
   }
   else {
      m_VorticityOp = NULL;
   }
   
   CFG::IntVect phi_ghost_vect = m_fluid_ghosts*CFG::IntVect::Unit;
   m_initial_phi.define( m_phase_geometry->magGeom().gridsFull(), 1, phi_ghost_vect );

   GyroaverageOperatorFactory::createOps(m_gyroavg_ops, a_state.dataKinetic());
   {
     KineticSpeciesPtrVect& kin_species(m_Y.dataKinetic());
     for (int s=0; s<kin_species.size(); s++) {
       KineticSpecies& species = *(kin_species[s]);
       if (species.isGyrokinetic()) species.gyroaverageOp(m_gyroavg_ops[species.name()]);
     }
   }

   m_fluidOp->setStepConstKinCoeff( m_step_const_kin_coeff_fluidop, 
                                    fluid_species );

   m_is_defined = true;
}


void GKOps::defineShell(  const GKState&            a_state, 
                          const GKSystemParameters& a_gkparams )
{
   CH_assert( isDefined()==false );
   m_phase_geometry = a_state.geometry();
   CH_assert( m_phase_geometry != NULL );

   setParameters( a_gkparams );
   m_units = new GKUnits();
   m_fluidOp = new CFG::GKFluidOp(  m_phase_geometry->magGeom(), 
                                    *m_units, 
                                    m_verbosity );
 
   CFG::IntVect phi_ghost_vect = m_fluid_ghosts*CFG::IntVect::Unit;
   m_initial_phi.define( m_phase_geometry->magGeom().gridsFull(), 1, phi_ghost_vect );

   m_is_defined = true;
}


void GKOps::initializeElectricField( const GKState& a_state_phys,
                                     const int      a_cur_step,
                                     const double   a_cur_time )
{
   GKState initial_state( a_state_phys );
   
   const KineticSpeciesPtrVect kinetic_species( a_state_phys.dataKinetic() );
   KineticSpeciesPtrVect& initial_kinetic_species( initial_state.dataKinetic() );
   for (int s(0); s<kinetic_species.size(); s++) {
      initial_kinetic_species[s] = kinetic_species[s]->clone( IntVect::Unit );
   }

   const CFG::FluidSpeciesPtrVect fluid_species( a_state_phys.dataFluid() );
   CFG::FluidSpeciesPtrVect& initial_fluid_species( initial_state.dataFluid() );
   for (int s(0); s<fluid_species.size(); s++) {
      initial_fluid_species[s] = fluid_species[s]->clone( CFG::IntVect::Unit );
   }

   const ScalarPtrVect scalars( a_state_phys.dataScalar() );
   ScalarPtrVect& initial_scalars( initial_state.dataScalar() );
   for (int s(0); s<scalars.size(); s++) {
      initial_scalars[s] = scalars[s]->clone();
   }

   if ( a_cur_step > 0 ) {
      // If we are restarting, initial_kinetic_species contains the current
      // solution, so overwrite it with the initial condition. In particular,
      // we need this to pass the initial ion density into some electron models
      initializeState( initial_state, a_cur_step, a_cur_time );

      //NB: Although the input parametr a_state_phys is physical satate, initialize state gives
      //mapped_space (a.k.a. comp_state). Thus we need to get rid of J to bring it to physical space.
      for (int s(0); s<initial_kinetic_species.size(); s++) {
         LevelData<FArrayBox>& dfn( initial_kinetic_species[s]->distributionFunction() );
         const PhaseGeom& geometry( initial_kinetic_species[s]->phaseSpaceGeometry() );
         geometry.divideJonValid( dfn );
      }

      // Initializee the potential if the state includes it
      if ( m_state_includes_potential ) {
         const CFG::FluidSpeciesPtrVect& fluids( a_state_phys.dataFluid() );
         const CFG::LevelData<CFG::FArrayBox>& state_phi = fluids[a_state_phys.getFluidComponent("Potentials")]->cell_var("potential");
         
         for (CFG::DataIterator dit(state_phi.dataIterator()); dit.ok(); ++dit) {
            m_initial_phi[dit].copy(state_phi[dit]);
         }
      }
   }

   // Get a pointer to the electric field operation object and define it
   CFG::IntVect ghostVect = (m_phase_geometry->secondOrder()) ? CFG::IntVect::Zero : CFG::IntVect::Unit;
   if ( m_ampere_law ) {
      m_phi_ops = new CFG::PhiAmpereOps("", m_phase_geometry->magGeom());
   }
   else if ( m_consistent_potential_bcs ) {
      m_phi_ops = new CFG::PhiSelfConsistentBCOps("", m_phase_geometry->magGeom());
   }
   else {
      m_phi_ops = new CFG::PhiOps("", m_phase_geometry->magGeom());
   }
   m_phi_ops->define( m_units->larmorNumber(),
                      m_units->debyeNumber(),
                      initial_kinetic_species,
                      m_boltzmann_electron,
                      m_fixed_efield,
                      m_phase_geometry->divFreeVelocity(),
                      m_fluid_ghosts,
                      a_cur_step );
   /*
     Initialize the m_EM_fields and m_EM_fields_ImOpImEx objects by:
       (a) defining (allocating) the potential (and nodal potential if the divergence-free GK velocity
           option is being used),
       (b) defining (allocating) the E-field (cell-centered, face-centered, injected, etc.)
       (c) setting the potential in m_EM_fields to m_initial_phi (in case it's needed as an initial
           condition for an update of the potential as described in the following bullet
       (d) computing the E-field,  A new potential (and resulting E-field) is computed if 
           (1) the potential is not a state variable, and if
           (2) the old vorticity model is not being used, and if
           (3) the m_ampere_law option (which doesn't compute a potential when using Boltzmann electrons) is
               not being used, and if
           (4) the E-field is not fixed
           Otherwise, the E-field is computed from m_initial_phi.  In all cases, m_EM_fields will
           contain both the E-field and the potential used to compute it.
   */

   bool include_phi_node = m_phase_geometry->divFreeVelocity();
   m_phi_ops->definePotential(m_initial_phi, include_phi_node, m_EM_fields);
   m_phi_ops->definePotential(m_initial_phi, include_phi_node, m_EM_fields_ImOpImEx);
   m_phi_ops->defineEField(m_EM_fields);
   m_phi_ops->defineEField(m_EM_fields_ImOpImEx);

   m_EM_fields.setPhi(m_initial_phi);

   bool recompute_potential = !m_state_includes_potential && !m_old_vorticity_model && !m_ampere_law && !m_fixed_efield;

   m_phi_ops->updateEField( m_EM_fields,
                            recompute_potential,
                            kinetic_species,
                            fluid_species,
                            scalars,
                            m_boundary_conditions->getEllipticOpBC(),
                            true,
                            true );
      
   //Improve Er field calculation to take into accout the dealigment between the grid and magnetic surfaces
   //Should not be used with poloidal variations: FIX LATER!!!

   CH_assert( m_phase_geometry != NULL );
   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   const CFG::MagCoordSys& coords = *mag_geom.getCoordSys();
   if ( (typeid(coords) == typeid(CFG::SingleNullCoordSys)) && (m_dealignment_corrections)) {
      mag_geom.interpolateErFromMagFS(m_EM_fields.getEFieldFace(), m_EM_fields.getEFieldCell());
   }

   if ( m_old_vorticity_model ) {
      Vector<Real> scalar_data;
      if (m_consistent_potential_bcs) {
         scalar_data = a_state_phys.dataScalar()[m_Y.getScalarComponent("Er_boundary")]->data();
      }
      if (m_EM_effects) {
         ((CFG::EMVorticityBEOp*)m_VorticityOp)->initializeBEModel(a_state_phys.dataKinetic(), a_state_phys.dataFluid(), scalar_data);
      }
      else{
         ((CFG::VorticityBEOp*)m_VorticityOp)->initializeBEModel(a_state_phys.dataKinetic(), a_state_phys.dataFluid(), scalar_data);
      }
   }

   m_Apar_ops = new CFG::AparOps("", m_phase_geometry->magGeom());
   m_Apar_ops->defineAparData(m_EM_fields);
}


void GKOps::initializeElectricFieldShell( const GKState& a_state_phys,
                                          const int      a_cur_step,
                                          const double   a_cur_time )
{
   // Get a pointer to the electric field operation object and define it
   CFG::IntVect ghostVect = (m_phase_geometry->secondOrder()) ? CFG::IntVect::Zero : CFG::IntVect::Unit;
   if ( m_ampere_law ) {
      m_phi_ops = new CFG::PhiAmpereOps("", m_phase_geometry->magGeom());
   }
   else if ( m_consistent_potential_bcs ) {
      m_phi_ops = new CFG::PhiSelfConsistentBCOps("", m_phase_geometry->magGeom());
   }
   else {
      m_phi_ops = new CFG::PhiOps("", m_phase_geometry->magGeom());
   }
   m_phi_ops->define( m_units->larmorNumber(),
                      m_units->debyeNumber(),
                      a_state_phys.dataKinetic(),
                      m_boltzmann_electron,
                      m_fixed_efield,
                      m_phase_geometry->divFreeVelocity(),
                      m_fluid_ghosts,
                      a_cur_step );
   /*
     Initialize the m_EM_fields and m_EM_fields_ImOpImEx objects by:
       (a) defining (allocating) the potential (and nodel potential if the divergence-free GK velocity
           option is being used),
       (b) defining (allocating) the E-field (cell-centered, face-centered, injected, etc.)
   */

   bool include_phi_node = m_phase_geometry->divFreeVelocity();
   m_phi_ops->definePotential(m_initial_phi, include_phi_node, m_EM_fields);
   m_phi_ops->definePotential(m_initial_phi, include_phi_node, m_EM_fields_ImOpImEx);
   m_phi_ops->defineEField(m_EM_fields);
   m_phi_ops->defineEField(m_EM_fields_ImOpImEx);
   
   CH_assert( m_phase_geometry != NULL );
}

void GKOps::initializeKineticSpeciesPhysical( KineticSpeciesPtrVect& a_species_comp,
                                              const double a_time )
{
   CH_TIME("GKOps::initializeKineticSpeciesPhysical()");
   m_kinetic_species_phys.resize( a_species_comp.size() );
   for (int s(0); s<a_species_comp.size(); s++) {
      m_kinetic_species_phys[s] = a_species_comp[s]->clone( m_kinetic_ghosts*IntVect::Unit );
      LevelData<FArrayBox>& dfn_phys( m_kinetic_species_phys[s]->distributionFunction() );
      const PhaseGeom& geometry( a_species_comp[s]->phaseSpaceGeometry() );
      geometry.divideJonValid( dfn_phys );
   }
}

void GKOps::initializeFluidSpeciesPhysical( CFG::FluidSpeciesPtrVect&  a_species_comp,
                                            const double a_time )
{
   CH_TIME("GKOps::initializeFluidSpeciesPhysical()");
   int num_species = a_species_comp.size();
   CFG::IntVect ghost_vect_cfg = m_fluid_ghosts*CFG::IntVect::Unit;
   m_fluid_species_phys.resize(num_species);
   for (int species(0); species<num_species; species++) {
      m_fluid_species_phys[species] = a_species_comp[species]->clone(ghost_vect_cfg, true);
      m_fluid_species_phys[species]->convertToPhysical();
   }   
   m_fluidOp->fillGhostCells( m_fluid_species_phys, a_time );
   if ( a_time == 0.0 ) {
      m_fluidOp->initializeWithBC( a_species_comp, m_fluid_species_phys, a_time ); 
   }

}

GKOps::~GKOps()
{
   if (m_phi_ops) delete m_phi_ops;
   if (m_Apar_ops) delete m_Apar_ops;
   if (m_boltzmann_electron) delete m_boltzmann_electron;
   if (m_VorticityOp) delete m_VorticityOp;
   if (m_collisions) delete m_collisions;
   if (m_fluidOp) delete m_fluidOp;
   if (m_scalarOp) delete m_scalarOp;
   if (m_transport) delete m_transport;
   if (m_neutrals)  delete m_neutrals;
   if (m_vlasov) delete m_vlasov;
   for (int ihist(0); ihist<m_hist_count; ihist++) {
     FieldHist *save_hist = &m_fieldHistLists[ihist];
     delete save_hist->fieldvals;
     delete save_hist->timevals;
     delete save_hist->timestep;
   }
   if (m_initial_conditions) delete m_initial_conditions;
   if (m_boundary_conditions) delete m_boundary_conditions;
   delete m_units;

   m_kinetic_species_phys.clear();
   m_fluid_species_phys.clear();

   GyroaverageOperatorFactory::deleteOps(m_gyroavg_ops);
}

Real GKOps::stableDtExpl( const GKState& a_state, const int a_step_number )
{
   std::string who = "nobody";
   Real dt_stable( DBL_MAX );
   for (int i = 0; i < m_dt_vlasov_expl.size(); i++) {
     Real dt = m_dt_vlasov_expl[i].second;
     if (dt < dt_stable) {
       dt_stable = dt;
       who = "Vlasov (" + m_dt_vlasov_expl[i].first + ")";
     }
   }
   if (m_dt_collisions_expl > 0.) {
     if (m_dt_collisions_expl < dt_stable) {
       dt_stable = m_dt_collisions_expl;
       who = "collisions";
     }
   }

   if (m_dt_fluids_expl < dt_stable) {
     dt_stable = m_dt_fluids_expl;
     who = "fluids";
   }

   if (m_dt_scalars_expl < dt_stable) {
     dt_stable = m_dt_scalars_expl;
     who = "scalars";
   }

   if (m_transport_model_on) {
      if (m_dt_transport_expl < dt_stable) {
        dt_stable = m_dt_transport_expl;
        who = "transport";
      }
   }

   if (m_neutrals_model_on) {
      if (m_dt_neutrals_expl < dt_stable) {
        dt_stable = m_dt_neutrals_expl;
        who = "neutrals";
      }
   }

   if (!procID()) {
     if (who != "nobody") {
      printf("Stable dt for explicit time integration: %f", dt_stable);
      printf(", fastest physics: %s\n", who.c_str());
     }
   }

   return dt_stable;
}

Real GKOps::stableDtImEx( const GKState& a_state, const int a_step_number )
{
   std::string who = "nobody";
   Real dt_stable( DBL_MAX );
   for (int i = 0; i < m_dt_vlasov_imex.size(); i++) {
     Real dt = m_dt_vlasov_imex[i].second;
     if (dt < dt_stable) {
       dt_stable = dt;
       who = "Vlasov (" + m_dt_vlasov_imex[i].first + ")";
     }
   }
   if (m_dt_collisions_imex > 0.) {
     if (m_dt_collisions_imex < dt_stable) {
       dt_stable = m_dt_collisions_imex;
       who = "collisions";
     }
   }

   if (m_dt_fluids_imex < dt_stable) {
     dt_stable = m_dt_fluids_imex;
     who = "fluids";
   }

   if (m_dt_scalars_imex < dt_stable) {
     dt_stable = m_dt_scalars_imex;
     who = "scalars";
   }

   if (m_transport_model_on) {
      if (m_dt_transport_imex < dt_stable) {
        dt_stable = m_dt_transport_imex;
        who = "transport";
      }
   }

   if (m_neutrals_model_on) {
      if (m_dt_neutrals_imex < dt_stable) {
        dt_stable = m_dt_neutrals_imex;
        who = "neutrals";
      }
   }

   if (!procID()) {
     if (who != "nobody") {
       printf("Stable dt for IMEX time integration: %f", dt_stable);
       printf(", fastest physics on explicit side: %s\n", who.c_str());
     }
   }

   return dt_stable;
}

void GKOps::initializeTI  ( const int       a_step, 
                            const Real      a_time, 
                            GKRHSData&      a_rhs,
                            const GKState&  a_state_comp, 
                            const GKState&  a_state_phys )
{
   /* Many of the physics implementations allocate data structures
    * the first time they are needed. This function runs once through
    * the RHS calculation steps necessary to get all such data allocated
    * and sensibly initialized so that the time integrators can create
    * preconditioners before beginning the time steps. */

   CH_TIMERS("GKOps::initializeTI()");

   const KineticSpeciesPtrVect& soln_comp( a_state_comp.dataKinetic() );
   const KineticSpeciesPtrVect& soln_phys( a_state_phys.dataKinetic() );
   m_collisions->preTimeStep( soln_comp, a_time, soln_phys );

   preOpEval(a_state_comp, a_time, chkpt_ti_init);
   explicitOp(a_rhs, a_time, a_state_comp);

   if (!trivialSolutionOp()) {
     preSolutionOpEval(a_state_comp, a_time, chkpt_ti_init);
     solutionOp(a_rhs, a_time, a_state_comp);
   }
   
   return;
}


void GKOps::preTimeStep (const int       a_step, 
                         const Real      a_time, 
                         const GKState&  a_state_comp, 
                         const GKState&  a_state_phys )
{
   CH_TIMERS("GKOps::preTimeStep()");
   CH_TIMER("compute_dt",t_compute_dt);

   if(!m_no_efield) {
      updatePhysicalSpeciesVector( a_state_comp.dataFluid(), a_time);
      updatePhysicalSpeciesVector( a_state_comp.dataKinetic(), a_time, false);
      setEMFields( a_state_comp, a_step, m_EM_fields );
   }

   const KineticSpeciesPtrVect& soln_comp( a_state_comp.dataKinetic() );
   const KineticSpeciesPtrVect& soln_phys( a_state_phys.dataKinetic() );
   m_collisions->preTimeStep( soln_comp, a_time, soln_phys );
   
   /* compute dt for each physics */

   CH_START(t_compute_dt);
   if ( m_vlasov && m_vlasov->computeCfl()) {
      m_dt_vlasov_expl = m_vlasov->computeDtExplicitTI(m_EM_fields, soln_comp, a_time);
      m_dt_vlasov_imex = m_vlasov->computeDtImExTI(m_EM_fields, soln_comp, a_time);
      m_time_scale_vlasov = m_vlasov->computeTimeScales(m_EM_fields, soln_comp, a_time);
   }
   else {
      m_dt_vlasov_expl.clear();
      m_dt_vlasov_imex.clear();
      m_time_scale_vlasov.clear();
   }

   m_dt_collisions_expl = m_collisions->computeDtExplicitTI( soln_comp );
   m_dt_collisions_imex = m_collisions->computeDtImExTI( soln_comp );
   m_time_scale_collisions = m_collisions->computeTimeScale( soln_comp );

   const CFG::FluidSpeciesPtrVect& fluids_comp( a_state_comp.dataFluid() );
   m_dt_fluids_expl = m_fluidOp->computeDtExplicitTI( fluids_comp );
   m_dt_fluids_imex = m_fluidOp->computeDtImExTI( fluids_comp );
   m_time_scale_fluids = m_fluidOp->computeTimeScale( fluids_comp );

   const ScalarPtrVect& scalars_comp( a_state_comp.dataScalar() );
   m_dt_scalars_expl = m_scalarOp->computeDtExplicitTI( scalars_comp );
   m_dt_scalars_imex = m_scalarOp->computeDtImExTI( scalars_comp );
   m_time_scale_scalars = m_scalarOp->computeTimeScale( scalars_comp );

   if (m_transport_model_on) {
      m_transport->preTimeStep( soln_comp, a_time, soln_phys );
      m_dt_transport_expl = m_transport->computeDtExplicitTI( soln_comp );
      m_dt_transport_imex = m_transport->computeDtImExTI( soln_comp );
      m_time_scale_transport = m_transport->computeTimeScale( soln_comp );
   }
   if (m_neutrals_model_on) {
      m_neutrals->preTimeStep( soln_comp, a_time, soln_phys );
      m_dt_neutrals_expl = m_neutrals->computeDtExplicitTI( soln_phys );
      m_dt_neutrals_imex = m_neutrals->computeDtImExTI( soln_phys );
      m_time_scale_neutrals = m_neutrals->computeTimeScale( soln_phys );
   }
   CH_STOP(t_compute_dt);
   
   preOpEval(a_state_comp, a_time, chkpt_pre_time_step);

   return;
}


void GKOps::postTimeStep( const int   a_step,
                          const Real  a_dt,
                          const Real  a_time,
                          GKState&    a_state )
{
   CH_TIME("GKOps::postTimeStep");
   if ( m_state_includes_potential ) {
      const CFG::FluidSpeciesPtrVect& fluids( a_state.dataFluid() );
      const CFG::LevelData<CFG::FArrayBox>& state_phi = fluids[a_state.getFluidComponent("Potentials")]->cell_var("potential");
      m_EM_fields.setPhi(state_phi);
   }

   if ( m_ampere_post_step_update && m_ampere_law ) {

      const ScalarPtrVect& scalars( a_state.dataScalar() );
      const Vector<Real>& scalar_data = scalars[a_state.getScalarComponent("Er_boundary")]->data();

      const CFG::FluidSpeciesPtrVect& fluids( a_state.dataFluid() );
      const CFG::LevelData<CFG::FArrayBox>& Er 
         = fluids[a_state.getFluidComponent("Er_flux_surfaces")]->cell_var("Er_flux_surfaces");

      (static_cast<CFG::PhiAmpereOps*>(m_phi_ops))->updateAverageAndPerturbation(scalar_data, Er, m_kinetic_species_phys, a_time);
   }

   if ( m_old_vorticity_model ) {

      updatePhysicalSpeciesVector(a_state.dataKinetic(), a_time);
      Vector<Real> scalar_data;
      if (m_consistent_potential_bcs) {
        scalar_data = a_state.dataScalar()[m_Y.getScalarComponent("Er_boundary")]->data();
      }

      bool cbc_or_ampere = m_consistent_potential_bcs || m_ampere_law;
      if (m_EM_effects) {
         ((CFG::EMVorticityBEOp*)m_VorticityOp)->updatePotentialBEModel(m_EM_fields,
                                                                        *m_phi_ops,
                                                                        cbc_or_ampere,
                                                                        m_kinetic_species_phys,
                                                                        a_state.dataFluid(),
                                                                        a_state.dataScalar(),
                                                                        m_boundary_conditions->getEllipticOpBC(),
                                                                        a_dt, a_time);
      }
      else {
         ((CFG::VorticityBEOp*)m_VorticityOp)->updatePotentialBEModel(m_EM_fields,
                                                                      *m_phi_ops,
                                                                      cbc_or_ampere,
                                                                      m_kinetic_species_phys,
                                                                      a_state.dataFluid(),
                                                                      a_state.dataScalar(),
                                                                      m_boundary_conditions->getEllipticOpBC(),
                                                                      a_dt, a_time);
      }
   }
  
   m_fluidOp->enforcePositivity( a_state.dataFluid() ); 
   updatePhysicalSpeciesVector( a_state.dataFluid(), a_time ); 
   m_fluidOp->postTimeEval(a_state.dataFluid(), m_fluid_species_phys, a_dt, a_time);
}

void GKOps::postTimeStage( const int   a_step,
                           const Real  a_time,
                           const Real  a_dt,
                           GKState&    a_state_comp,
                           const int   a_stage )
{
   CH_TIME("GKOps::postTimeStage()");
   // Stage number follows C convention: it goes from 0, ..., n_stages-1
   // Calculation of the field is skipped for stage 0, since we assume it
   // has already been calculated to compute the velocity used in the
   // time step estimate.

   updatePhysicalSpeciesVector( a_state_comp.dataFluid(), a_time ); 
   updatePhysicalSpeciesVector( a_state_comp.dataKinetic(), a_time, false ); 
   if (a_stage && !m_no_efield && !m_skip_efield_stage_update ) {
      setEMFields( a_state_comp, a_step, m_EM_fields );
   }
   fillPhysicalKineticSpeciesGhostCells(a_time);

   const KineticSpeciesPtrVect& soln( a_state_comp.dataKinetic() );
   m_collisions->postTimeStage( soln, a_time, a_stage );

   if ( !m_ampere_post_step_update && m_ampere_law ) {

      const ScalarPtrVect& scalars( a_state_comp.dataScalar() );
      const Vector<Real>& scalar_data = scalars[a_state_comp.getScalarComponent("Er_boundary")]->data();

      const CFG::FluidSpeciesPtrVect& fluids( a_state_comp.dataFluid() );
      const CFG::LevelData<CFG::FArrayBox>& Er = fluids[a_state_comp.getFluidComponent("Er_flux_surfaces")]->cell_var("Er_flux_surfaces");

      (static_cast<CFG::PhiAmpereOps*>(m_phi_ops))->updateAverageAndPerturbation(scalar_data, Er, m_kinetic_species_phys, a_time);
   }
   
   // post stage eval used for semi-implicit advance in relaxation scheme
   m_fluidOp->enforcePositivity( a_state_comp.dataFluid() );
   updatePhysicalSpeciesVector( a_state_comp.dataFluid(), a_time );
   m_fluidOp->postStageEval(  a_state_comp.dataFluid(), 
                              m_fluid_species_phys, 
                              a_stage, 
                              a_dt, 
                              a_time  );

}


void GKOps::preSolutionOpEval( const GKState&     a_state,
                               const Real         a_time,
                               const Checkpoint&  a_chkpt )
{
  CH_assert(!trivialSolutionOp());
  CH_TIMERS("GKOps::preSolutionOpEval");
  CH_TIMER("updatePhysicalSpeciesVector",t_update_phys_vector);
  CH_TIMER("fluidpreOpEval",t_fluid_preop_eval);

  // The purpose of this function is to enable the calculation of data prior to
  // calling solutionOp().  Presently, it is used to precompute the physical
  // kinetic species from the computational kinetic species, which is somewhat
  // expensive due to the need to fill ghost cells.  By caching the result of this
  // calculation, it can be reused in function evaluations in which it is known that
  // the distribution function hasn't changed between function calls, or we don't care
  // if it did, e.g., for an approximate preconditioner operator evaluation.
  // Dependence on class state like this is dangerous, but necessary for performance.
 
  if ((!m_enable_ti_optimizations) && (!m_no_efield)) {
    setEMFields( a_state, 
                 -1, 
                 m_EM_fields_ImOpImEx, 
                 !m_step_const_kin_coeff_fluidop );
  }

  if ((a_chkpt == chkpt_ti_init) || (!m_enable_ti_optimizations)) {

    updatePhysicalSpeciesVector( a_state.dataFluid(), a_time );

    CH_START(t_update_phys_vector);
    updatePhysicalSpeciesVector( a_state.dataKinetic(), a_time );
    CH_STOP(t_update_phys_vector);

    CH_START(t_fluid_preop_eval);
    m_fluidOp->preSolutionOpEval( m_kinetic_species_phys, 
                                  a_state.dataFluid(), 
                                  a_state.dataScalar(), 
                                  m_EM_fields, 
                                  a_time );
    CH_STOP(t_fluid_preop_eval);

  } else if (a_chkpt == chkpt_ti_advance_1) {

    /* in time step, before stage calculation */

    updatePhysicalSpeciesVector( a_state.dataFluid(), a_time );

    if (!m_step_const_kin_coeff_fluidop) {

      CH_START(t_update_phys_vector);
      updatePhysicalSpeciesVector( a_state.dataKinetic(), a_time );
      CH_STOP(t_update_phys_vector);
  
      CH_START(t_fluid_preop_eval);
      m_fluidOp->preSolutionOpEval( m_kinetic_species_phys, 
                                    a_state.dataFluid(), 
                                    a_state.dataScalar(), 
                                    m_EM_fields, 
                                    a_time );
      CH_STOP(t_fluid_preop_eval);

    }

  } else if (a_chkpt == chkpt_ti_advance_3) {

    /* in time step, before step completion */

    CH_START(t_update_phys_vector);
    updatePhysicalSpeciesVector( a_state.dataKinetic(), a_time );
    CH_STOP(t_update_phys_vector);
  
    updatePhysicalSpeciesVector( a_state.dataFluid(), a_time );

    if (!m_step_const_kin_coeff_fluidop) {

      CH_START(t_fluid_preop_eval);
      m_fluidOp->preSolutionOpEval( m_kinetic_species_phys, 
                                    a_state.dataFluid(), 
                                    a_state.dataScalar(), 
                                    m_EM_fields, 
                                    a_time );
      CH_STOP(t_fluid_preop_eval);

    }

  } else {

    /* this function should not be called from any other location in the code */
    MayDay::Error("GKOps::preSolutionOpEval() - invalid value for a_chkpt");

  }

  return;
}


void GKOps::preOpEval( const GKState&     a_state,
                       const Real         a_time,
                       const Checkpoint&  a_chkpt )
{
  CH_TIMERS("GKOps::preOpEval");
  CH_TIMER("updatePhysicalSpeciesVector",t_update_phys_vector);
  CH_TIMER("fluidpreOpEval",t_fluid_preop_eval);
  
  // The purpose of this function is to enable the calculation of data prior to
  // calling the operator evaluation class members (e.g., explicitOp, explicitOpImEx,
  // implicitOpImEx, etc.).  Presently, it is used to precompute the physical
  // kinetic species from the computational kinetic species, which is somewhat
  // expensive due to the need to divide by J and fill ghost cells.  By caching
  // the result of this calculation, it can be reused in function evaluations
  // in which it is known that the distribution function hasn't changed between
  // function calls, or we don't care if it did, e.g., for an approximate preconditioner
  // operator evaluation.  Dependence on class state like this is dangerous,
  // but necessary for performance.
  
  bool  update_kinetic_phys(true),
        update_fluid_phys(true),
        compute_fluid_pre_op(true);

  if (m_enable_ti_optimizations) {

    std::vector<std::string> implicit_vlasov_species(0);
    if (m_vlasov) m_vlasov->implicitSpecies(implicit_vlasov_species);
    bool is_vlasov_implicit = (implicit_vlasov_species.size() > 0);
  
    if (a_chkpt == chkpt_naive) {
  
      update_kinetic_phys   = true;
      update_fluid_phys     = true;
      compute_fluid_pre_op  = true;
  
    } else if (a_chkpt == chkpt_ti_init) {
  
      update_kinetic_phys   = true;
      update_fluid_phys     = true;
      compute_fluid_pre_op  = true;
  
    } else if (a_chkpt == chkpt_pre_time_step) {
  
      // pre-timestep //
  
      if (m_step_const_kin_coeff_fluidop) {
        update_kinetic_phys   = true;
        update_fluid_phys     = true;
        compute_fluid_pre_op  = true;
      } else {
        update_kinetic_phys   = false;
        update_fluid_phys     = false;
        compute_fluid_pre_op  = false;
      }
  
    } else if (a_chkpt == chkpt_ti_advance_2) {
  
      // in time step, after stage calculation //
  
      update_kinetic_phys   = false; // done in postTimeStage
      update_fluid_phys     = false; // done in postTimeStage
      compute_fluid_pre_op  = !m_step_const_kin_coeff_fluidop;
  
    } else if (a_chkpt == chkpt_stage_func_0) {
  
      // stage function evaluation for first Newton iteration //
  
      update_kinetic_phys   = true;
      update_fluid_phys     = true;
      compute_fluid_pre_op  = !m_step_const_kin_coeff_fluidop;
  
    } else if (a_chkpt == chkpt_stage_func_n) {
  
      // stage function evaluation for Newton iteration > 1 //
  
      update_kinetic_phys   = is_vlasov_implicit;
      update_fluid_phys     = true;
      compute_fluid_pre_op  = false;
  
    } else if (a_chkpt == chkpt_stage_jac) {
  
      // stage Jacobian evaluation //
  
      update_kinetic_phys   = is_vlasov_implicit;
      update_fluid_phys     = true;
      compute_fluid_pre_op  = false;
  
    } else {
  
      // this function should not be called from any other location in the code //
      MayDay::Error("GKOps::preOpEval() - invalid value for a_chkpt");
  
    }

  } else {
  
    update_kinetic_phys   = true;
    update_fluid_phys     = true;
    compute_fluid_pre_op  = true;

    if(!m_no_efield) {
      setEMFields( a_state, 
                   -1, 
                   m_EM_fields_ImOpImEx, 
                   !m_step_const_kin_coeff_fluidop );
    }

  }

  if ( update_kinetic_phys ) {
    CH_START(t_update_phys_vector);
    updatePhysicalSpeciesVector( a_state.dataKinetic(), a_time );
    CH_STOP(t_update_phys_vector);
  }

  if (update_fluid_phys) {
    updatePhysicalSpeciesVector( a_state.dataFluid(), a_time );
  }

  if (compute_fluid_pre_op) {
    CH_START(t_fluid_preop_eval);
    m_fluidOp->preOpEval( m_kinetic_species_phys, 
                          a_state.dataFluid(), 
                          a_state.dataScalar(), 
                          m_EM_fields, 
                          a_time );
    CH_STOP(t_fluid_preop_eval);
  }
  
  return;
}


void GKOps::solutionOp( GKRHSData&      a_rhs,
                        const Real      a_time,
                        const GKState&  a_state )
{
   CH_TIME("GKOps::solutionOp");

   a_rhs.zero();
   applyKineticSpeciesSolutionOperator( a_rhs.dataKinetic(), 
                                        a_state.dataKinetic(), 
                                        a_time  );
   
   applyFluidSpeciesSolutionOperator( a_rhs.dataFluid(), 
                                      a_state.dataKinetic(),
                                      m_kinetic_species_phys, 
                                      a_state.dataFluid(), 
                                      m_fluid_species_phys, 
                                      a_state.dataScalar(), 
                                      a_time );
   
   applyScalarSolutionOperator( a_rhs.dataScalar(), a_state.dataScalar(), a_time );
}

void GKOps::explicitOp( GKRHSData&      a_rhs,
                        const Real      a_time,
                        const GKState&  a_state_comp )
{
   CH_TIME("GKOps::explicitOp");

   a_rhs.zero();

   applyVlasovOperator( a_rhs, m_kinetic_species_phys, unsplit, a_time );

   if (m_transport_model_on) {
      applyTransportOperator( a_rhs.dataKinetic(), m_kinetic_species_phys, a_time );
   }

   if (m_neutrals_model_on) {
      applyNeutralsOperator( a_rhs.dataKinetic(), m_kinetic_species_phys,
                             m_fluid_species_phys, a_time );
   }

   applyCollisionOperator( a_rhs.dataKinetic(), a_state_comp.dataKinetic(), unsplit, a_time );
   
   applyFluidOperator( a_rhs.dataFluid(), m_kinetic_species_phys, m_fluid_species_phys,
                       a_state_comp.dataScalar(), m_EM_fields, unsplit, a_time );

   applyScalarOperator( a_rhs, m_kinetic_species_phys, m_fluid_species_phys,
                        a_state_comp.dataScalar(), m_EM_fields, *m_phi_ops, unsplit, true, a_time );
}

void GKOps::explicitOpImEx( GKRHSData&      a_rhs,
                            const Real      a_time,
                            const GKState&  a_state_comp )
{
   CH_TIME("GKOps::explicitOpImEx");

   a_rhs.zero();

   applyVlasovOperator( a_rhs, m_kinetic_species_phys, imex_exp, a_time );

   if (m_transport_model_on && (!m_transport_imex_implicit) ) {
      applyTransportOperator( a_rhs.dataKinetic(), m_kinetic_species_phys, a_time );
   }

   applyCollisionOperator( a_rhs.dataKinetic(), a_state_comp.dataKinetic(), imex_exp, a_time );

   if (m_neutrals_model_on && (!m_neutrals_imex_implicit) ) {
      applyNeutralsOperator( a_rhs.dataKinetic(), m_kinetic_species_phys,
                             m_fluid_species_phys, a_time );
   }

   applyFluidOperator( a_rhs.dataFluid(), m_kinetic_species_phys, m_fluid_species_phys,
                       a_state_comp.dataScalar(), m_EM_fields, imex_exp, a_time );

   applyScalarOperator( a_rhs, m_kinetic_species_phys, m_fluid_species_phys,
                        a_state_comp.dataScalar(), m_EM_fields, *m_phi_ops, imex_exp, true, a_time );
}

void GKOps::implicitOpImEx( GKRHSData&      a_rhs,
                            const Real      a_time,
                            const GKState&  a_state_comp )
{
   CH_TIME("GKOps::implicitOpImEx");

   if ((!m_no_efield) && (m_enable_ti_optimizations)) {
      // The last argument is for whether E_field will be injected or not;
      // If m_step_const_kin_coeff_fluidop = true, it means that kinietc distribution is not updated;
      // so we do not need to inject E_field

      setEMFields( a_state_comp, 
                   -1, 
                   m_EM_fields_ImOpImEx, 
                   !m_step_const_kin_coeff_fluidop );
   }

   a_rhs.zero();

   applyVlasovOperator( a_rhs, m_kinetic_species_phys, imex_imp, a_time );

   if (m_transport_model_on && (m_transport_imex_implicit) ) {
      applyTransportOperator( a_rhs.dataKinetic(), m_kinetic_species_phys, a_time );
   }

   applyCollisionOperator( a_rhs.dataKinetic(), a_state_comp.dataKinetic(), imex_imp, a_time );

   if (m_neutrals_model_on && (m_neutrals_imex_implicit) ) {
      applyNeutralsOperator(  a_rhs.dataKinetic(), 
                              m_kinetic_species_phys,
                              m_fluid_species_phys, 
                              a_time );
   }

   CFG::EMFields* EMfields;
   if (m_enable_ti_optimizations) {
     EMfields = &m_EM_fields_ImOpImEx;
   } else {
     EMfields = &m_EM_fields;
   }

   applyFluidOperator(  a_rhs.dataFluid(), 
                        m_kinetic_species_phys, 
                        m_fluid_species_phys,
                        a_state_comp.dataScalar(), 
                        *EMfields, 
                        imex_imp, 
                        a_time );

   applyScalarOperator( a_rhs, 
                        m_kinetic_species_phys, 
                        m_fluid_species_phys,
                        a_state_comp.dataScalar(), 
                        *EMfields, 
                        *m_phi_ops, 
                        imex_imp, 
                        true, 
                        a_time );

   return;
}

void GKOps::solutionPC( GKRHSData&     a_rhs,
                        const Real     a_time,
                        const GKState& a_state )
{
   CH_TIME("GKOps::solutionPC");

   /* This function is the same as solutionOp, but it is called by
    * the preconditioner, not the actual time integrator or Newton solver
    * or linear solver. Thus, inaccuracies in the result of this function
    * does not affect the overall solution.
   */

   a_rhs.zero();

   applyKineticSpeciesSolutionOperator(a_rhs.dataKinetic(), a_state.dataKinetic(), a_time);
     
   updatePhysicalSpeciesVector( a_state.dataFluid(), a_time); // needed?
   
   applyFluidSpeciesSolutionOperator( a_rhs.dataFluid(), 
                                      a_state.dataKinetic(),
                                      m_kinetic_species_phys, 
                                      a_state.dataFluid(),
                                      m_fluid_species_phys, 
                                      a_state.dataScalar(), 
                                      a_time );
   
   applyScalarSolutionOperator( a_rhs.dataScalar(), a_state.dataScalar(), a_time );
}

void GKOps::explicitPC( GKRHSData&     a_rhs,
                        const Real     a_time,
                        const GKState& a_state )
{
   CH_TIME("GKOps::explicitPC");

   /* This function is the same as explicitOp, but it is called by
    * the preconditioner, not the actual time integrator or Newton solver
    * or linear solver. Thus, inaccuracies in the result of this function
    * does not affect the overall solution.
   */

   a_rhs.zero();

   applyVlasovOperator( a_rhs, m_kinetic_species_phys, unsplit, a_time );
     
   updatePhysicalSpeciesVector( a_state.dataFluid(), a_time); // needed?

   if (m_transport_model_on) {
      applyTransportOperator( a_rhs.dataKinetic(), m_kinetic_species_phys, a_time );
   }

   if (m_neutrals_model_on) {
      applyNeutralsOperator( a_rhs.dataKinetic(), m_kinetic_species_phys,
                             m_fluid_species_phys, a_time );
   }

   applyCollisionOperator( a_rhs.dataKinetic(), a_state.dataKinetic(), unsplit, a_time );
   
   applyFluidOperator( a_rhs.dataFluid(), m_kinetic_species_phys, m_fluid_species_phys, a_state.dataScalar(),
                       m_EM_fields, unsplit, a_time );

   applyScalarOperator( a_rhs, m_kinetic_species_phys, m_fluid_species_phys,
                        a_state.dataScalar(), m_EM_fields, *m_phi_ops, unsplit, true, a_time );
}

void GKOps::implicitPCImEx( GKRHSData&     a_rhs,
                            const Real     a_time,
                            const GKState& a_state_comp )
{
   CH_TIME("GKOps::implicitPCImEx");
   updatePhysicalSpeciesVector( a_state_comp.dataFluid(), a_time); // needed?
   implicitOpImEx(a_rhs, a_time, a_state_comp);
}

void GKOps::defineMultiPhysicsPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>&  a_pc,
                                  std::vector<DOFList>&                             a_dof_list,
                                  const ODEVector&                                  a_Y,
                                  const GKState&                                    a_state,
                                  void*                                             a_sys,
                                  const std::string&                                a_out_string,
                                  const std::string&                                a_opt_string,
                                  bool                                              a_im,
                                  const int                                         a_id )
{
  const GlobalDOF* global_dof(a_state.getGlobalDOF());

  if (m_vlasov) {
    m_vlasov->defineMultiPhysicsPC( a_pc,
                                    a_dof_list,
                                    a_state.dataKinetic(),
                                    global_dof->dataKinetic(),
                                    a_Y,
                                    a_sys,
                                    a_out_string,
                                    a_opt_string,
                                    a_im,
                                    a_id );
  }

  m_collisions->defineMultiPhysicsPC( a_pc, 
                                      a_dof_list,
                                      a_state.dataKinetic(),
                                      global_dof->dataKinetic(),
                                      a_Y, 
                                      a_sys, 
                                      a_out_string, 
                                      a_opt_string, 
                                      a_im,
                                      a_id );

  m_fluidOp->defineMultiPhysicsPC(  a_pc,
                                    a_dof_list,
                                    a_state.dataFluid(),
                                    global_dof->dataFluid(),
                                    a_Y,
                                    a_sys,
                                    a_out_string,
                                    a_opt_string,
                                    a_im,
                                    a_id );

  return;
}

void GKOps::updateMultiPhysicsPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>&  a_pc,
                                  const GKState&                                    a_state,
                                  const Real                                        a_time,
                                  const int                                         a_step,
                                  const int                                         a_stage,
                                  const Real                                        a_shift,
                                  const bool                                        a_im )
{
  const GlobalDOF* global_dof(a_state.getGlobalDOF());

  if (m_vlasov) {
    m_vlasov->updateMultiPhysicsPC( a_pc,
                                    a_state.dataKinetic(),
                                    m_EM_fields,
                                    global_dof->dataKinetic(),
                                    a_time, 
                                    a_step,
                                    a_stage,
                                    a_shift, 
                                    a_im  );
  }

  m_collisions->updateMultiPhysicsPC( a_pc, 
                                      a_state.dataKinetic(), 
                                      global_dof->dataKinetic(), 
                                      a_time, 
                                      a_step,
                                      a_stage,
                                      a_shift, 
                                      a_im  );

  m_fluidOp->updateMultiPhysicsPC(  a_pc,
                                    m_kinetic_species_phys,
                                    a_state.dataFluid(),
                                    a_time, 
                                    a_step,
                                    a_stage,
                                    a_shift, 
                                    a_im );
  return;
}

void GKOps::solveSolutionPC( GKRHSData&     a_rhs,
                             const GKState& a_state,
                             int            a_idx )
{
   CH_TIME("GKOps::solveSolutionPC");
   m_fluidOp->solveSolutionPC(  a_rhs.dataFluid(), 
                                a_state.dataKinetic(), 
                                a_state.dataFluid(), 
                                a_idx );
}

void GKOps::solveFluidOpPCImEx( GKRHSData&          a_rhs,
                                const GKState&      a_state,
                                const std::string&  a_name,
                                int                 a_idx )
{
   CH_TIME("GKOps::solvePCImEx");
   m_fluidOp->solvePCImEx(  a_rhs.dataFluid(), 
                            a_state.dataKinetic(), 
                            a_state.dataFluid(), 
                            a_name,
                            a_idx );
}

void GKOps::solveVlasovOpPCImEx(Preconditioner<ODEVector,AppCtxt>* const a_pc,
                                GKRHSData&                               a_rhs,
                                const GKState&                           a_state,
                                int                                      a_idx )
{
   CH_TIME("GKOps::solvePCImEx");
   m_vlasov->solvePCImEx( a_pc,
                          a_rhs.dataKinetic(), 
                          a_state.dataKinetic(), 
                          a_state.dataFluid(), 
                          a_idx );
}


inline
void GKOps::setEMFields( const GKState&   a_state_comp,
                         const int        a_step,
                         CFG::EMFields&   a_EM_fields,
                         const bool       a_phase_space_update ) const
{
   CH_TIME("GKOps::setEMFields");

   if ( m_state_includes_potential ) {
      const CFG::FluidSpeciesPtrVect& fluids( a_state_comp.dataFluid() );
      const CFG::CFGVars* species_ptr = fluids[a_state_comp.getFluidComponent("Potentials")];
      const CFG::Potentials* potential_vars = dynamic_cast<const CFG::Potentials*>(species_ptr);
      CH_assert(potential_vars);

      const CFG::LevelData<CFG::FArrayBox>& phi_state = potential_vars->cell_var("potential");
      a_EM_fields.setPhi(phi_state);
      m_phi_ops->updateEField( a_EM_fields,
                               false,
                               m_kinetic_species_phys,
                               m_fluid_species_phys,
                               a_state_comp.dataScalar(),
                               m_boundary_conditions->getEllipticOpBC(),
                               a_phase_space_update,
                               false );

      if ( potential_vars->evolvingApar() && a_phase_space_update ) {
      
         const CFG::LevelData<CFG::FArrayBox>& Apar_state = potential_vars->cell_var("A_parallel");

         m_Apar_ops->computeAparDerivs( a_EM_fields,
                                        m_kinetic_species_phys,
                                        Apar_state,
                                        m_fluidOp->fluidModel("Potentials"),
                                        true );
      }
   }
   else {

      // If the potential is not a state variable (i.e., there is not a separate equation
      // describing it's evolution) but needs to be updated as the result of a change in
      // state variables, we do so here.  Otherwise, we just update the field using the
      // potential currently in a_EM_fields.
      bool compute_potential = !(m_old_vorticity_model || m_fixed_efield);

      m_phi_ops->updateEField( a_EM_fields,
                               compute_potential,
                               m_kinetic_species_phys,
                               m_fluid_species_phys,
                               a_state_comp.dataScalar(),
                               m_boundary_conditions->getEllipticOpBC(),
                               a_phase_space_update,
                               false );
   }
}



inline
void GKOps::updatePhysicalSpeciesVector( const KineticSpeciesPtrVect&  a_species_comp,
                                         const double                  a_time,
                                         const bool                    a_fill_ghost_cells )
{
   CH_TIMERS("GKOps::updatePhysicalSpeciesVector (kinetic)");

   divideJ( a_species_comp, m_kinetic_species_phys );

   if (a_fill_ghost_cells) {
      m_boundary_conditions->fillGhostCells(m_kinetic_species_phys,
                                            m_EM_fields,
                                            a_time);
   }
}


inline
void GKOps::fillPhysicalKineticSpeciesGhostCells( const double a_time )
{
   CH_TIMERS("GKOps::fillPhysicalKineticSpeciesGhostCells");

   m_boundary_conditions->fillGhostCells( m_kinetic_species_phys,
                                          m_EM_fields,
                                          a_time);
}


inline
void GKOps::updatePhysicalSpeciesVector( const CFG::FluidSpeciesPtrVect&  a_species_comp,
                                         const double                     a_time,
                                         const bool                       a_fill_ghost_cells )
{
   CH_TIME("GKOps::updatePhysicalSpeciesVector (fluid)");

   m_fluidOp->convertToPhysical( m_fluid_species_phys, a_species_comp ); 

   if (a_fill_ghost_cells) {
     m_fluidOp->fillGhostCells( m_fluid_species_phys, a_time );
   }
}


inline
void GKOps::fillPhysicalFluidSpeciesGhostCells( const double a_time )
{
   CH_TIME("GKOps::fillPhysicalFluidSpeciesGhostCells");
   m_fluidOp->fillGhostCells( m_fluid_species_phys, a_time );
}


inline
void GKOps::applyVlasovOperator( GKRHSData&                    a_rhs,
                                 const KineticSpeciesPtrVect&  a_kinetic_phys,
                                 const opType                  a_op_type,
                                 const Real&                   a_time )
{
   CH_TIME("GKOps::applyVlasovOperator");

   if ( m_vlasov ) {
      m_count_vlasov++;
      if (a_op_type == unsplit) {
         m_vlasov->accumulateRHS(a_rhs, a_kinetic_phys, m_EM_fields, *m_phi_ops, false, a_time );
         m_vlasov->accumulateRHS(a_rhs, a_kinetic_phys, m_EM_fields,  *m_phi_ops, true, a_time );
      } else if (a_op_type == imex_exp) {
         m_vlasov->accumulateRHS(a_rhs, a_kinetic_phys, m_EM_fields,  *m_phi_ops, false, a_time );
      } else if (a_op_type == imex_imp) {
         m_vlasov->accumulateRHS(a_rhs, a_kinetic_phys, m_EM_fields,  *m_phi_ops, true, a_time );
      } else {
        MayDay::Error("GKOps::applyVlasovOperator: invalid value for argument a_op_type");
      }
   }
}

inline
void GKOps::applyKineticSpeciesSolutionOperator( KineticSpeciesPtrVect&        a_rhs,
                                                 const KineticSpeciesPtrVect&  a_soln,
                                                 const Real&                   a_time )
{
   for (int s(0); s<a_rhs.size(); s++) {
     LevelData<FArrayBox>& rhs_dfn(a_rhs[s]->distributionFunction());
     const LevelData<FArrayBox>& soln_dfn(a_soln[s]->distributionFunction());

     for (DataIterator dit(rhs_dfn.dataIterator()); dit.ok(); ++dit) {
        rhs_dfn[dit].copy(soln_dfn[dit]);
     }
   }
}

inline
void GKOps::applyFluidSpeciesSolutionOperator( CFG::FluidSpeciesPtrVect&        a_rhs,
                                               const KineticSpeciesPtrVect&     a_kinetic_species_comp,
                                               const KineticSpeciesPtrVect&     a_kinetic_species_phys,
                                               const CFG::FluidSpeciesPtrVect&  a_fluid_species_comp,
                                               const CFG::FluidSpeciesPtrVect&  a_fluid_species_phys,
                                               const ScalarPtrVect&             a_scalars,
                                               const Real&                      a_time )
{
   m_fluidOp->evalSolutionOp( a_rhs, 
                              a_kinetic_species_comp, 
                              a_kinetic_species_phys, 
                              a_fluid_species_comp, 
                              a_fluid_species_phys, 
                              a_scalars, 
                              a_time );
}

inline
void GKOps::applyScalarSolutionOperator( ScalarPtrVect&        a_rhs,
                                         const ScalarPtrVect&  a_soln,
                                         const Real&           a_time )
{
   for (int s(0); s<a_rhs.size(); s++) {
      Vector<Real>& rhs_data(a_rhs[s]->data());
      const Vector<Real>& soln_data(a_soln[s]->data());
      CH_assert(rhs_data.size() == soln_data.size());

      for (int i=0; i<rhs_data.size(); ++i) {
         rhs_data[i] = soln_data[i];
      }
   }
}

inline
void GKOps::applyCollisionOperator( KineticSpeciesPtrVect&       a_rhs,
                                    const KineticSpeciesPtrVect& a_soln,
                                    const opType                 a_op_type,
                                    const Real&                  a_time )
{
   m_count_collision++;
   if (a_op_type == imex_exp) {
     m_collisions->accumulateRHS( a_rhs, a_soln, false, a_time );
   } else if (a_op_type == imex_imp) {
     m_collisions->accumulateRHS( a_rhs, a_soln, true, a_time );
   } else if (a_op_type == unsplit) {
     m_collisions->accumulateRHS( a_rhs, a_soln, false, a_time );
     m_collisions->accumulateRHS( a_rhs, a_soln, true, a_time );
   } else {
     MayDay::Error("Unknown a_op_type!");
   }
}

inline
void GKOps::applyTransportOperator( KineticSpeciesPtrVect&       a_rhs,
                                    const KineticSpeciesPtrVect& a_soln,
                                    const Real&                  a_time )
{
   m_count_transport++;
   m_transport->accumulateRHS( a_rhs, a_soln, a_time );

}

inline
void GKOps::applyNeutralsOperator( KineticSpeciesPtrVect&            a_rhs,
                                   const KineticSpeciesPtrVect&      a_kinetic_species_phys,
                                   const CFG::FluidSpeciesPtrVect&   a_fluid_species_phys,
                                   const Real&                       a_time )
{
    m_count_neutrals++;
    m_neutrals->accumulateRHS( a_rhs, a_kinetic_species_phys, a_fluid_species_phys, a_time );
    
}


inline
void GKOps::applyFluidOperator( CFG::FluidSpeciesPtrVect&                  a_rhs,
                                const KineticSpeciesPtrVect&               a_kinetic_species_phys,
                                const CFG::FluidSpeciesPtrVect&            a_fluid_species,
                                const ScalarPtrVect&                       a_scalars,
                                //                                const CFG::EField&                         a_E_field,
                                const CFG::EMFields&                       a_EM_fields,
                                const opType                               a_op_type,
                                const Real&                                a_time)
{
   CH_TIME("GKOps::applyFluidOperator");

   m_count_fluids++;
   if (a_op_type == imex_exp) {
      m_fluidOp->accumulateRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_EM_fields,
                                false, a_time );
   } else if (a_op_type == imex_imp) {
      m_fluidOp->accumulateRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_EM_fields,
                                true, a_time );
   } else if (a_op_type == unsplit) { 
      m_fluidOp->accumulateRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_EM_fields,
                                false, a_time );
      m_fluidOp->accumulateRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_EM_fields,
                                true, a_time );
   } else {
     MayDay::Error("Unknown a_op_type!");
   }
}


inline
void GKOps::applyScalarOperator( GKRHSData&                       a_rhs,
                                 const KineticSpeciesPtrVect&     a_kinetic_species,
                                 const CFG::FluidSpeciesPtrVect&  a_fluid_species,
                                 const ScalarPtrVect&             a_scalars,
                                 const CFG::EMFields&             a_EM_fields,
                                 const CFG::PhiOps&               a_phi_ops,
                                 const opType                     a_op_type,
                                 const bool                       a_recompute_kinetic_terms,
                                 const Real&                      a_time)
{
   CH_TIME("GKOps::applyScalarOperator");

   m_count_scalars++;
   //m_count_fluids++;
   if (a_op_type == imex_exp) {
      m_scalarOp->accumulateRHS( a_rhs, a_kinetic_species, a_fluid_species, a_scalars, a_EM_fields, a_phi_ops,
                                 false, a_recompute_kinetic_terms, a_time );
   } else if (a_op_type == imex_imp) {
      m_scalarOp->accumulateRHS( a_rhs, a_kinetic_species, a_fluid_species, a_scalars, a_EM_fields, a_phi_ops,
                                 true, a_recompute_kinetic_terms, a_time );
   } else if (a_op_type == unsplit) {
      m_scalarOp->accumulateRHS( a_rhs, a_kinetic_species, a_fluid_species, a_scalars, a_EM_fields, a_phi_ops,
                                 false, a_recompute_kinetic_terms, a_time );
      m_scalarOp->accumulateRHS( a_rhs, a_kinetic_species, a_fluid_species, a_scalars, a_EM_fields, a_phi_ops,
                                 true, a_recompute_kinetic_terms, a_time );
   } else {
     MayDay::Error("Unknown a_op_type!");
   }
}


inline
void setZero( CFG::LevelData<CFG::FArrayBox>& a_data )
{
   for (CFG::DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit].setVal(0.);
   }
}


void GKOps::computeSignedChargeDensities( CFG::LevelData<CFG::FArrayBox>& a_pos_charge_density,
                                          CFG::LevelData<CFG::FArrayBox>& a_neg_charge_density,
                                          const KineticSpeciesPtrVect& a_species ) const
{
   // Container for individual species charge density
   CH_assert( m_phase_geometry != NULL );
   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   CFG::LevelData<CFG::FArrayBox> species_charge_density( mag_geom.gridsFull(),
                                                          1,
                                                          CFG::IntVect::Zero );

   setZero( a_pos_charge_density );
   setZero( a_neg_charge_density );

   for (int species(0); species<a_species.size(); species++) {
      const KineticSpecies& this_species( *(a_species[species]) );

      // Compute the charge density for this species
      this_species.chargeDensity( species_charge_density );

      if ( this_species.charge() > 0.0 ) {
         CFG::DataIterator dit( a_pos_charge_density.dataIterator() );
         for (dit.begin(); dit.ok(); ++dit) {
            a_pos_charge_density[dit].plus( species_charge_density[dit] );
         }
      }
      else {
         CFG::DataIterator dit( a_neg_charge_density.dataIterator() );
         for (dit.begin(); dit.ok(); ++dit) {
            a_neg_charge_density[dit].plus( species_charge_density[dit] );
         }
      }
   }
}


void GKOps::setErAverage( const LevelData<FluxBox>& Er_face_injected )
{
   CFG::LevelData<CFG::FluxBox> Er_average_face(m_phase_geometry->magGeom().gridsFull(), 3, CFG::IntVect::Unit);
   for (CFG::DataIterator dit(Er_average_face.disjointBoxLayout()); dit.ok(); ++dit) {
      Er_average_face[dit].setVal(0.);
   }
   m_phase_geometry->projectPhaseToConfiguration(Er_face_injected, Er_average_face);

   CFG::PhiAmpereOps& phi_ops = static_cast<CFG::PhiAmpereOps&>(*m_phi_ops);
   phi_ops.setErAverage(Er_average_face);
}

void GKOps::setErAverage( const LevelData<FArrayBox>& Er_cell_injected)
{
   CFG::LevelData<CFG::FArrayBox> Er_average_cell(m_phase_geometry->magGeom().gridsFull(), 3, CFG::IntVect::Unit);
   m_phase_geometry->projectPhaseToConfiguration(Er_cell_injected, Er_average_cell);

   CFG::PhiAmpereOps& phi_ops = static_cast<CFG::PhiAmpereOps&>(*m_phi_ops);
   phi_ops.setErAverage(Er_average_cell);
}

void GKOps::setETilde( const LevelData<FluxBox>& E_tilde_face_injected )
{
   CFG::LevelData<CFG::FluxBox> E_tilde_face(m_phase_geometry->magGeom().gridsFull(), 3, CFG::IntVect::Unit);
   for (CFG::DataIterator dit(E_tilde_face.disjointBoxLayout()); dit.ok(); ++dit) {
      E_tilde_face[dit].setVal(0.);
   }
   m_phase_geometry->projectPhaseToConfiguration(E_tilde_face_injected, E_tilde_face);

   CFG::PhiAmpereOps& phi_ops = static_cast<CFG::PhiAmpereOps&>(*m_phi_ops);
   phi_ops.setETilde(E_tilde_face);
}

void GKOps::setETilde( const LevelData<FArrayBox>& E_tilde_cell_injected)
{
   CFG::LevelData<CFG::FArrayBox> E_tilde_cell(m_phase_geometry->magGeom().gridsFull(), 3, CFG::IntVect::Unit);
   m_phase_geometry->projectPhaseToConfiguration(E_tilde_cell_injected, E_tilde_cell);

   CFG::PhiAmpereOps& phi_ops = static_cast<CFG::PhiAmpereOps&>(*m_phi_ops);
   phi_ops.setETilde(E_tilde_cell);
}

void GKOps::setInitialPhi( const LevelData<FArrayBox>& phi_injected)
{
   // This member function is needed to set the initial phi from the value obtained from
   // a restart file when using the old vorticity model
   CFG::LevelData<CFG::FArrayBox> phi_tmp(m_phase_geometry->magGeom().gridsFull(), 1, CFG::IntVect::Zero);
   m_phase_geometry->projectPhaseToConfiguration(phi_injected, phi_tmp);
   
   for (CFG::DataIterator dit(m_initial_phi.disjointBoxLayout()); dit.ok(); ++dit) {
      m_initial_phi[dit].copy(phi_tmp[dit]);
   }
}


void GKOps::divideJ( const KineticSpeciesPtrVect& a_soln_mapped,
                     KineticSpeciesPtrVect&       a_soln_physical ) const
{
   for (int species(0); species<a_soln_physical.size(); species++) {
      const KineticSpecies& soln_species_mapped( *(a_soln_mapped[species]) );
      KineticSpecies& soln_species_physical( *(a_soln_physical[species]) );

      const LevelData<FArrayBox> & dfn_mapped = soln_species_mapped.distributionFunction();
      LevelData<FArrayBox> & dfn_physical = soln_species_physical.distributionFunction();

      DataIterator dit = dfn_physical.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
         dfn_physical[dit].copy( dfn_mapped[dit] );
      }

      const PhaseGeom& geometry( soln_species_mapped.phaseSpaceGeometry() );
      geometry.divideJonValid( dfn_physical );
   }
}


void GKOps::convertToPhysical( const GKState& a_soln_mapped, GKState& a_soln_physical,
                               const Real a_cur_time )
{
   CH_TIME("GKOps::convertToPhysical()");
   divideJ( a_soln_mapped.dataKinetic(), a_soln_physical.dataKinetic() );

   
   const CFG::FluidSpeciesPtrVect& soln_mapped_fluid = a_soln_mapped.dataFluid();
   CFG::FluidSpeciesPtrVect& soln_physical_fluid = a_soln_physical.dataFluid();
   //for (int s=0; s<soln_mapped_fluid.size(); ++s) {
   //   soln_physical_fluid[s]->copy(*soln_mapped_fluid[s]);
   //   soln_physical_fluid[s]->convertToPhysical();
   //}
   m_fluidOp->convertToPhysical( soln_physical_fluid, soln_mapped_fluid );
   
 
   //CFG::FluidSpeciesPtrVect& soln_physical_fluid = a_soln_physical.dataFluid();
   //updatePhysicalSpeciesVector( a_soln_mapped.dataFluid(), a_cur_time );
   //const CFG::FluidSpeciesPtrVect& soln_phys_fluid_wghosts = m_fluid_species_phys;
   //for (int s=0; s<soln_physical_fluid.size(); ++s) {
   //   soln_physical_fluid[s]->copy(*soln_phys_fluid_wghosts[s]);
   //}

   const ScalarPtrVect& soln_mapped_scalar = a_soln_mapped.dataScalar();
   const ScalarPtrVect& soln_physical_scalar = a_soln_physical.dataScalar();
   for (int s=0; s<soln_mapped_scalar.size(); ++s) {
      Vector<Real>& this_physical_scalar_vec = soln_physical_scalar[s]->data();
      Vector<Real>& this_mapped_scalar_vec = soln_mapped_scalar[s]->data();
      for (int n=0; n<this_physical_scalar_vec.size(); ++n) {
         this_physical_scalar_vec[n] = this_mapped_scalar_vec[n];
      }
   }
}


void GKOps::setParameters( const GKSystemParameters& a_gkparams )
{
  m_fixed_efield = a_gkparams.fixedEField();
  m_no_efield = a_gkparams.noEField();

  m_transport_imex_implicit = a_gkparams.imexTransportImplicit();
  m_neutrals_imex_implicit = a_gkparams.imexNeutralsImplicit();

  m_old_vorticity_model = a_gkparams.oldVorticityModel();
  m_EM_effects = a_gkparams.includeEMeffects();
  m_skip_efield_stage_update = a_gkparams.skipEFieldStageUpdate();

  m_consistent_potential_bcs = a_gkparams.consistentPotentialBCs();
  m_use_vlasov_ampere_bcs = a_gkparams.useVlasovAmpereBCs();
  m_Esol_extrapolation = a_gkparams.ESolExtrapolation();

  m_dealignment_corrections = a_gkparams.dealignmentCorrections();

  m_ampere_law = a_gkparams.ampereLaw();
  if (m_ampere_law) {
    m_consistent_potential_bcs = true;
    m_ampere_cold_electrons = a_gkparams.ampereColdElectrons();
  } else {
    m_ampere_cold_electrons = false;
  }
  m_ampere_post_step_update = a_gkparams.amperePostStepUpdate();

  if ( m_fixed_efield && m_ampere_law ) {
     MayDay::Error("GKOps::setParameters(): Specify either fixed field or ampere law, but not both"); 
  }

  m_transport_model_on = a_gkparams.transportModelOn();
  m_neutrals_model_on = a_gkparams.neutralsModelOn();

  m_enforce_quasineutrality = a_gkparams.enforceQuasineutrality();
  m_step_const_kin_coeff_fluidop = a_gkparams.stepConstKineticCoeffFluidOp();
  m_enable_ti_optimizations = a_gkparams.enableTIOptimizations();

  m_kinetic_ghosts = a_gkparams.kineticGhosts();
  m_fluid_ghosts = a_gkparams.fluidGhosts();

  bool using_boltzmann_electrons(true);
  m_boltzmann_electron = NULL;
  ParmParse pp_be( "boltzmann_electron" );
  
  string name;
  if (using_boltzmann_electrons && pp_be.contains("name")) {
     pp_be.get( "name", name );
     using_boltzmann_electrons = (name=="electron");
  }
  else using_boltzmann_electrons = false;
  
  double mass;
  if (using_boltzmann_electrons && pp_be.contains("mass")) {
     pp_be.get( "mass", mass );
  }
  else using_boltzmann_electrons = false;
  
  double charge;
  if (using_boltzmann_electrons && pp_be.contains("charge")) {
     pp_be.get( "charge", charge );
  }
  else using_boltzmann_electrons = false;
  
  double temperature;
  if (using_boltzmann_electrons && pp_be.contains("temperature")) {
     pp_be.get( "temperature", temperature );
  }
  else using_boltzmann_electrons = false;
  
  if (using_boltzmann_electrons) {
     CH_assert( m_phase_geometry != NULL );
     const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
     CFG::LevelData<CFG::FArrayBox> temp( mag_geom.gridsFull(),
                                          1,
                                          CFG::IntVect::Unit );
     for (CFG::DataIterator dit( temp.dataIterator() ); dit.ok(); ++dit) {
        temp[dit].setVal( temperature );
     }
     
     m_boltzmann_electron = new CFG::BoltzmannElectron( mass,
                                                        charge,
                                                        mag_geom,
                                                        temp );
  }
}

void GKOps::getEField(  CFG::LevelData<CFG::FArrayBox>& a_E ) const
{
   const CFG::LevelData<CFG::FArrayBox>& E_field(m_EM_fields.getEFieldCell());
   a_E.define(E_field.disjointBoxLayout(), 3, CFG::IntVect::Zero);
   for (auto dit(a_E.dataIterator()); dit.ok(); ++dit) {
      a_E[dit].copy(E_field[dit]);
   }
   return;
}

void GKOps::getEField(  CFG::LevelData<CFG::FluxBox>& a_E ) const
{
   const CFG::LevelData<CFG::FluxBox>& E_field(m_EM_fields.getEFieldFace());
   a_E.define(E_field.disjointBoxLayout(), 3, CFG::IntVect::Zero);
   for (auto dit(a_E.dataIterator()); dit.ok(); ++dit) {
      a_E[dit].copy(E_field[dit]);
   }
   return;
}

void GKOps::getEField( LevelData<FluxBox>& a_E) const
{
   const LevelData<FluxBox>& E_field(m_EM_fields.getEField());
   a_E.define(E_field.disjointBoxLayout(), 3, IntVect::Zero);
   for (auto dit(a_E.dataIterator()); dit.ok(); ++dit) {
      a_E[dit].copy(E_field[dit]);
   }
   return;
}

void GKOps::setupFieldHistories( ParmParse& a_ppsim )
{
   //query whether we should write history files
   a_ppsim.query( "histories", m_history );

   if (m_history) {
#ifdef CH_USE_HDF5
      //query how frequently to save histories
      a_ppsim.query( "history_frequency", m_hist_freq );
      //query for indices to generate the history.   If out of bounds, will result in no history.
      std::vector<int> read_hist_indices( CFG_DIM );
      // Set up default to be in middle of array
      bool finding_indices(true);
      // look for index/field pairs numbered sequentially 1, 2, 3, ...
      // e.g., "simulation.1.history_indices"
      // with any break in the enumeration terminating the search
      while (finding_indices) {
         int count( m_hist_count + 1 );
         stringstream sind;
         sind << count << ".history_indices" << ends; // grid indices
         
         stringstream sfield;
         sfield << count << ".history_field" << ends; // e.g., "potential"

         if ( a_ppsim.contains( sind.str().c_str() )   &&
              a_ppsim.contains( sfield.str().c_str() ) ){
            a_ppsim.getarr( sind.str().c_str(), read_hist_indices, 0, CFG_DIM );
            for (int d(0); d<CFG_DIM; d++) {
               m_hist_indices[d] = read_hist_indices[d];
            }

            // query to see what field's history to accumulate
            a_ppsim.get( sfield.str().c_str(), m_hist_fieldname );
            if (m_hist_fieldname == "potential") {
               FieldHist *save_hist = new FieldHist; //create structure
               save_hist->hist_index = m_hist_count;
               save_hist->grid_indices = m_hist_indices; // grid indices
               save_hist->fieldvals = new Vector<Real>( m_expand_incr ); // allocate memory for vectors
               save_hist->timevals = new Vector<Real>( m_expand_incr );
               save_hist->timestep = new Vector<int>( m_expand_incr );
               save_hist->valsize = m_expand_incr; // length of vectors above
               save_hist->cur_entry = 0; // incremented for each entry
               save_hist->fieldname = m_hist_fieldname;
               m_fieldHistLists.push_back(*save_hist); //save structure in array
               m_hist_count++; // count of watchpoint histories in deck
            }
            else if (m_hist_fieldname == "Efield") {
               FieldHist *save_hist = new FieldHist; //create structure
               save_hist->hist_index = m_hist_count;
               save_hist->grid_indices = m_hist_indices; // grid indices
               save_hist->fieldvals = new Vector<Real>( m_expand_incr ); // allocate memory for vectors
               save_hist->timevals = new Vector<Real>( m_expand_incr );
               save_hist->timestep = new Vector<int>( m_expand_incr );
               save_hist->valsize = m_expand_incr; // length of vectors above
               save_hist->cur_entry = 0; // incremented for each entry
               save_hist->fieldname = m_hist_fieldname;
               m_fieldHistLists.push_back(*save_hist); //save structure in array
               m_hist_count++; // count of watchpoint histories in deck
            } else {
               MayDay::Error("Unimplemented field name");
            }
         }
         else {
            if (m_hist_count == 0) {
               MayDay::Error( "If histories are requested, history_field and history_indices must also be specified" );
            } else {
               finding_indices = false;
            }
         }
      }
#else
   MayDay::Error( "histories only defined with hdf5" );
#endif
   }
}


void GKOps::writeFieldHistory( const int a_cur_step,
                               const Real a_cur_time,
                               const bool a_startup_flag,
                               const std::string& a_prefix )
{
   if (a_cur_step % m_hist_freq != 0 && !a_startup_flag) {
      return;
   }
   
   for (int ihist(0); ihist<m_hist_count; ihist++) {
      FieldHist *field_hist_ptr = &m_fieldHistLists[ihist];
      std::string hist_field_name = field_hist_ptr->fieldname;
      CFG:: IntVect inode_hist = field_hist_ptr->grid_indices;
      
      std::string fname;
      std::string fname_suf;
      std::string fname_step;
      const CFG::LevelData<CFG::FArrayBox>* field = NULL;

      if (hist_field_name == "potential") {
         field = &(m_EM_fields.getPhi());
         fname = a_prefix + "potential_hist_";
         stringstream fname_temp;
         fname_temp << fname << ihist + 1 << ".curve" << ends;
         fname_suf = fname_temp.str();
         stringstream fname_step_temp;
         fname_step_temp << fname << ihist + 1 << ".txt" << ends;
         fname_step = fname_step_temp.str();
      }
      else if (hist_field_name == "Efield") {

         CFG::LevelData<CFG::FArrayBox>& E_field_cell = m_EM_fields.getEFieldCell();
         field = &E_field_cell;
         fname = a_prefix + "Efield_hist_";
         stringstream fname_temp;
         fname_temp << fname << ihist + 1 << ".curve" << ends;
         fname_suf = fname_temp.str();
         stringstream fname_step_temp;
         fname_step_temp << fname << ihist + 1 << ".txt" << ends;
         fname_step = fname_step_temp.str();
      }
      else{
         MayDay::Error( "Unimplemented field name" );
      }
      
      // Writes value of a spatial field at a specified point inode_hist to a text file
      // with name fname, adding to what is there
      Real field_val(0.0);
      for (CFG::DataIterator dit( field->dataIterator() ); dit.ok(); ++dit) {
         
        // Extract local fortran array for the field on this patch
        // This differs from syntax in many other parts of code in that a_field is now
        // a pointer, so we need to de-reference it.
         const CFG::FArrayBox& field_on_patch( (*field)[dit] );

         // now loop over nodes of the box and if the node's global indices match
         //   the specified global indices, print.  This will only happen on one
         //   patch, unless it happens in a ghost cell -- accept that duplication
         //   for now.
         CFG::Box patchBox( (field->disjointBoxLayout())[dit] );

         bool found(true);
         for (int i(0); i<CFG_DIM; i++) {
            if ((inode_hist[i]<patchBox.smallEnd(i)) ||
                (inode_hist[i]>patchBox.bigEnd(i))   ){
               found = false;
            }
         }
         if (found) {
            field_val += field_on_patch( inode_hist, 0 );
         }
      }
      
      Real field_val_sum(0.0);
#ifdef CH_MPI
      MPI_Allreduce( &field_val, &field_val_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
      field_val_sum = field_val;
#endif

      const int cur_index( field_hist_ptr->cur_entry );
      if (cur_index >= field_hist_ptr->valsize) {
         const int length( cur_index + m_expand_incr );
         (*field_hist_ptr->fieldvals).resize( length );
         (*field_hist_ptr->timevals).resize( length );
         (*field_hist_ptr->timestep).resize( length );
         field_hist_ptr->valsize = length;
      }
      (*field_hist_ptr->fieldvals)[cur_index] = field_val_sum;
      (*field_hist_ptr->timevals)[cur_index]  = a_cur_time;
      (*field_hist_ptr->timestep)[cur_index]  = a_cur_step;
      field_hist_ptr->cur_entry += 1;
      
#ifdef CH_MPI
      if (procID()==0) { // write out a file for each watchpoint
#endif
         //overwrite any existing file; this one is for Visit (with no column for the step)
         ofstream fieldout( fname_suf.c_str(), ios::out ); 
#if CFG_DIM==3
         fieldout << "# " << inode_hist[0] << ", " << inode_hist[1] << ", " << inode_hist[2] << endl;
#else
         fieldout << "# " << inode_hist[0] << ", " << inode_hist[1] << endl;
#endif

         //overwrite any existing file; this one is for human viewing (1st column for the step)
         ofstream fieldout_step( fname_step.c_str(), ios::out ); 
#if CFG_DIM==3
         fieldout_step << "# step  time  potential at ix,iy,iz = " << inode_hist[0] << ", " << inode_hist[1] << ", " << inode_hist[2]<<  endl;
#else
         fieldout_step << "# step  time  potential at ix,iy = " << inode_hist[0] << ", " << inode_hist[1] <<  endl;
#endif

         for (int j(0); j<field_hist_ptr->cur_entry; j++) {
            fieldout << (*field_hist_ptr->timevals)[j] << " "
                     << (*field_hist_ptr->fieldvals)[j] << endl;
            fieldout_step << (*field_hist_ptr->timestep)[j] << " "
                          << (*field_hist_ptr->timevals)[j] << " "
                          << (*field_hist_ptr->fieldvals)[j] << endl;
         }
         fieldout.close();
         fieldout_step.close();
#ifdef CH_MPI
      }
#endif
   }
}


void GKOps::writeCheckpointFile( HDF5Handle& a_handle ) const
{
   char buff[100];
   hsize_t flatdims[1], count[1], sizebuff[1];

   for (int ihist(0); ihist<m_hist_count; ihist++) {

      const FieldHist *field_hist_ptr = &m_fieldHistLists[ihist];

      sprintf(buff,"field_hist_%d", ihist+1);
      a_handle.setGroup(buff);

      flatdims[0] = 1;
      count[0] = 1;
      sizebuff[0] = field_hist_ptr->valsize;

      hid_t sizedataspace = H5Screate_simple(1, flatdims, NULL);
#ifdef H516
      hid_t sizedataset   = H5Dcreate(a_handle.groupID(), "size",
                        H5T_NATIVE_INT, sizedataspace, H5P_DEFAULT);
#else
      hid_t sizedataset   = H5Dcreate(a_handle.groupID(), "size",
                        H5T_NATIVE_INT, sizedataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
      hid_t szmemdataspace = H5Screate_simple(1, count, NULL);
      H5Dwrite(sizedataset, H5T_NATIVE_INT, szmemdataspace, sizedataspace,
             H5P_DEFAULT, sizebuff);
      H5Dclose(sizedataset);

      int indices[CFG_DIM];

      for (int i = 0; i < CFG_DIM; i++) {
          indices[i] = (field_hist_ptr->grid_indices)[i];
      }

      flatdims[0] =  CFG_DIM;
      count[0] =  CFG_DIM;

      hid_t indexdataspace = H5Screate_simple(1, flatdims, NULL);
#ifdef H516
      hid_t indexdataset   = H5Dcreate(a_handle.groupID(), "indices",
                        H5T_NATIVE_INT, indexdataspace, H5P_DEFAULT);
#else
      hid_t indexdataset   = H5Dcreate(a_handle.groupID(), "indices",
                        H5T_NATIVE_INT, indexdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
      hid_t imemdataspace = H5Screate_simple(1, count, NULL);
      H5Dwrite(indexdataset, H5T_NATIVE_INT, imemdataspace, indexdataspace,
             H5P_DEFAULT, indices);
      H5Dclose(indexdataset);
      flatdims[0] =  field_hist_ptr->valsize;
      count[0] =  field_hist_ptr->valsize;

      //cout << "valsize = " << field_hist_ptr->valsize << endl;

      Real *rbuff = new Real[field_hist_ptr->valsize];

      for (int i = 0; i < field_hist_ptr->valsize; i++) {
        rbuff[i] = (*field_hist_ptr->fieldvals)[i];
      }
      hid_t fielddataspace = H5Screate_simple(1, flatdims, NULL);
#ifdef H516
      hid_t fielddataset   = H5Dcreate(a_handle.groupID(), "field",
                        H5T_NATIVE_REAL, fielddataspace, H5P_DEFAULT);
#else
      hid_t fielddataset   = H5Dcreate(a_handle.groupID(), "field",
                        H5T_NATIVE_REAL, fielddataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
      hid_t fmemdataspace = H5Screate_simple(1, count, NULL);
      H5Dwrite(fielddataset, H5T_NATIVE_REAL, fmemdataspace, fielddataspace,
             H5P_DEFAULT, rbuff);
      H5Dclose(fielddataset);

      for (int i = 0; i < field_hist_ptr->valsize; i++) {
        rbuff[i] = (*field_hist_ptr->timevals)[i];
      }
      hid_t timedataspace = H5Screate_simple(1, flatdims, NULL);
#ifdef H516
      hid_t timedataset   = H5Dcreate(a_handle.groupID(), "times",
                                      H5T_NATIVE_REAL, timedataspace, 
                                      H5P_DEFAULT);
      
#else
      hid_t timedataset   = H5Dcreate(a_handle.groupID(), "times",
                        H5T_NATIVE_REAL, timedataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
      hid_t tmemdataspace = H5Screate_simple(1, count, NULL);
      H5Dwrite(timedataset, H5T_NATIVE_REAL, tmemdataspace, timedataspace,
             H5P_DEFAULT, rbuff);
      H5Dclose(timedataset);

      delete [] rbuff;

      int *ibuff = new int[field_hist_ptr->valsize];
      for (int i = 0; i < field_hist_ptr->valsize; i++) {
        ibuff[i] = (*field_hist_ptr->timestep)[i];
      }
      hid_t stepdataspace = H5Screate_simple(1, flatdims, NULL);
#ifdef H516
      hid_t stepdataset   = H5Dcreate(a_handle.groupID(), "steps",
                        H5T_NATIVE_INT, stepdataspace, H5P_DEFAULT);
#else
      hid_t stepdataset   = H5Dcreate(a_handle.groupID(), "steps",
                        H5T_NATIVE_INT, stepdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
      hid_t smemdataspace = H5Screate_simple(1, count, NULL);
      H5Dwrite(stepdataset, H5T_NATIVE_INT, smemdataspace, stepdataspace,
             H5P_DEFAULT, ibuff);
      H5Dclose(stepdataset);

      delete [] ibuff;

   }
}
   
void GKOps::readCheckpointFile( HDF5Handle& a_handle, const int& a_cur_step )
{
   char buff[100];
   hsize_t flatdims[1], count[1], sizebuff[1];

   for (int ihist(0); ihist<m_hist_count; ihist++) {
      
      FieldHist *field_hist_ptr = &m_fieldHistLists[ihist];
      
      sprintf(buff,"field_hist_%d", ihist+1);
      a_handle.setGroup(buff);
      
      flatdims[0] =  1;
      count[0] = 1;
      hsize_t dims[1], maxdims[1];

#ifdef H516
      hid_t sizedataset   = H5Dopen(a_handle.groupID(), "size");
#else
      hid_t sizedataset   = H5Dopen(a_handle.groupID(), "size", H5P_DEFAULT);
#endif
      hid_t sizedataspace = H5Dget_space(sizedataset);
      H5Sget_simple_extent_dims(sizedataspace, dims, maxdims);
      hid_t szmemdataspace = H5Screate_simple(1, dims, NULL);
      H5Dread(sizedataset, H5T_NATIVE_INT, szmemdataspace, sizedataspace,
             H5P_DEFAULT, sizebuff);
      H5Dclose(sizedataset);

      if (sizebuff[0] >= field_hist_ptr->valsize) {
         (*field_hist_ptr->fieldvals).resize(sizebuff[0]);
         (*field_hist_ptr->timevals).resize(sizebuff[0]);
         (*field_hist_ptr->timestep).resize(sizebuff[0]);
         field_hist_ptr->valsize = sizebuff[0];
      }

      int indices[CFG_DIM];
      int readin_indices[CFG_DIM];

      for (int i = 0; i < CFG_DIM; i++) {
          indices[i] = (field_hist_ptr->grid_indices)[i];
      }
      flatdims[0] =  CFG_DIM;
      count[0] =  CFG_DIM;

#ifdef H516
      hid_t indexdataset   = H5Dopen(a_handle.groupID(), "indices");
#else
      hid_t indexdataset   = H5Dopen(a_handle.groupID(), "indices", H5P_DEFAULT);
#endif
      hid_t indexdataspace = H5Dget_space(indexdataset);
      H5Sget_simple_extent_dims(indexdataspace, dims, maxdims);
      hid_t imemdataspace = H5Screate_simple(1, dims, NULL);
      H5Dread(indexdataset, H5T_NATIVE_INT, imemdataspace, indexdataspace,
             H5P_DEFAULT, readin_indices);
      H5Dclose(indexdataset);
      for (int i = 0; i < CFG_DIM; i++) {
        if (indices[i] != readin_indices[i]) {
          MayDay::Error("Grid indices for field history don't match previous run");
        }
      }

      flatdims[0] =  field_hist_ptr->valsize;
      count[0] =  field_hist_ptr->valsize;

      // cout << "valsize = " << field_hist_ptr->valsize << endl;

      Real *rbuff = new Real[field_hist_ptr->valsize];
#ifdef H516
      hid_t fielddataset   = H5Dopen(a_handle.groupID(), "field");
#else
      hid_t fielddataset   = H5Dopen(a_handle.groupID(), "field", H5P_DEFAULT);
#endif
      hid_t fielddataspace = H5Dget_space(fielddataset);
      H5Sget_simple_extent_dims(fielddataspace, dims, maxdims);
      hid_t fmemdataspace = H5Screate_simple(1, dims, NULL);
      H5Dread(fielddataset, H5T_NATIVE_REAL, fmemdataspace, fielddataspace,
             H5P_DEFAULT, rbuff);
      H5Dclose(fielddataset);
      for (int i = 0; i < field_hist_ptr->valsize; i++) {
        (*field_hist_ptr->fieldvals)[i] = rbuff[i];
      }

#ifdef H516
      hid_t timedataset   = H5Dopen(a_handle.groupID(), "times");
#else
      hid_t timedataset   = H5Dopen(a_handle.groupID(), "times", H5P_DEFAULT);
#endif
      hid_t timedataspace = H5Dget_space(timedataset);
      H5Sget_simple_extent_dims(timedataspace, dims, maxdims);
      hid_t tmemdataspace = H5Screate_simple(1, dims, NULL);
      H5Dread(timedataset, H5T_NATIVE_REAL, tmemdataspace, timedataspace,
             H5P_DEFAULT, rbuff);
      H5Dclose(timedataset);
      for (int i = 0; i < field_hist_ptr->valsize; i++) {
        (*field_hist_ptr->timevals)[i] = rbuff[i];
      }

      delete [] rbuff;

      int *ibuff = new int[field_hist_ptr->valsize];
#ifdef H516
      hid_t stepdataset   = H5Dopen(a_handle.groupID(), "steps");
#else
      hid_t stepdataset   = H5Dopen(a_handle.groupID(), "steps", H5P_DEFAULT);
#endif
      hid_t stepdataspace = H5Dget_space(stepdataset);
      H5Sget_simple_extent_dims(stepdataspace, dims, maxdims);
      hid_t smemdataspace = H5Screate_simple(1, dims, NULL);
      H5Dread(stepdataset, H5T_NATIVE_INT, smemdataspace, stepdataspace,
             H5P_DEFAULT, ibuff);
      H5Dclose(stepdataset);
      for (int i = 0; i < field_hist_ptr->valsize; i++) {
        (*field_hist_ptr->timestep)[i] = ibuff[i];
      }

      field_hist_ptr->cur_entry = field_hist_ptr->valsize; //where to add new entries
      for (int i = 0; i < field_hist_ptr->valsize; i++) {
        if (ibuff[i] > a_cur_step || (i>0 && ibuff[i] == 0)) {
          field_hist_ptr->cur_entry = i; // first location in table for next entry
          break;
        }
      }

      delete [] ibuff;

   }
}


inline
KineticSpecies& findElectronSpecies( KineticSpeciesPtrVect& a_species )
{
   int electron_index(-1);
   for (int species(0); species<a_species.size(); species++) {
      if (a_species[species]->isSpecies( "electron" )) {
         electron_index = species;
      }
   }
   CH_assert( electron_index>-1 );
   return ( *(a_species[electron_index]) );
}


void GKOps::enforceQuasiNeutrality(
   KineticSpeciesPtrVect&          a_species,
   CFG::LevelData<CFG::FArrayBox>& a_potential ) const
{
   if (m_enforce_quasineutrality) {
      
      CH_assert( m_phase_geometry != NULL );
      const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
      const CFG::DisjointBoxLayout& grids( mag_geom.gridsFull() );
      CFG::LevelData<CFG::FArrayBox> ion_density( grids, 1, CFG::IntVect::Zero );
      CFG::LevelData<CFG::FArrayBox> electron_density( grids, 1, CFG::IntVect::Zero );
      
      computeSignedChargeDensities( ion_density, electron_density, a_species );
 
      CFG::IntVect cfg_nghosts = 2*CFG::IntVect::Unit;
      CFG::LevelData<CFG::FArrayBox> quasineutral_density( grids, 1, cfg_nghosts );
      computeQuasiNeutralElectronDensity( quasineutral_density,
                                          a_potential,
                                          m_boundary_conditions->getEllipticOpBC(),
                                          ion_density );

      KineticSpecies& electron_species( findElectronSpecies( a_species ) );
      renormalizeElectronDistributionFunction( electron_species,
                                               quasineutral_density,
                                               electron_density );
   }
}


inline
void GKOps::computeQuasiNeutralElectronDensity(
   CFG::LevelData<CFG::FArrayBox>&       a_quasineutral_density,
   CFG::LevelData<CFG::FArrayBox>&       a_potential,
   CFG::EllipticOpBC&                    a_bc, 
   const CFG::LevelData<CFG::FArrayBox>& a_ion_density) const
{
   m_phi_ops->computeQuasiNeutralElectronDensity( a_quasineutral_density,
                                                  a_potential,
                                                  a_bc,
                                                  a_ion_density );
}


inline
void GKOps::renormalizeElectronDistributionFunction(
   KineticSpecies&                 a_electron_species,
   CFG::LevelData<CFG::FArrayBox>& a_quasineutral_density,
   CFG::LevelData<CFG::FArrayBox>& a_initial_density ) const
{
   const int CELL_AVERAGE_TO_POINT_VALUE(-1);

   CFG::fourthOrderAverage( a_quasineutral_density, CELL_AVERAGE_TO_POINT_VALUE );
   LevelData<FArrayBox> injected_quasineutral_density;
   m_phase_geometry->injectConfigurationToPhase( a_quasineutral_density,
                                                 injected_quasineutral_density );

   CFG::fourthOrderAverage( a_initial_density, CELL_AVERAGE_TO_POINT_VALUE );
   LevelData<FArrayBox> injected_initial_density;
   m_phase_geometry->injectConfigurationToPhase( a_initial_density, injected_initial_density );

   LevelData<FArrayBox>& electron_dfn( a_electron_species.distributionFunction() );
   fourthOrderAverage( electron_dfn, CELL_AVERAGE_TO_POINT_VALUE );
   for (DataIterator dit( electron_dfn.dataIterator() ); dit.ok(); ++dit) {
      electron_dfn[dit].divide( injected_initial_density[dit] );
      electron_dfn[dit].mult( injected_quasineutral_density[dit] );
   }
   fourthOrderAverage( electron_dfn );
}


double GKOps::getLoRadialField(const GKState& a_state) const
{
   double value = 0.;

   if ( m_consistent_potential_bcs ) {

      const ScalarPtrVect& scalars( a_state.dataScalar() );
      int scalar_component = a_state.getScalarComponent("Er_boundary");
      const Vector<Real>& scalar_data = scalars[scalar_component]->data();
      value = scalar_data[0];
   }

   return value;
}

void GKOps::setLoRadialField(GKState& a_state,
                             const double   a_Er_lo) const
{
   if ( m_consistent_potential_bcs ) {

      ScalarPtrVect& scalars( a_state.dataScalar() );
      int scalar_component = a_state.getScalarComponent("Er_boundary");
      Vector<Real>& scalar_data = scalars[scalar_component]->data();
      scalar_data[0] = a_Er_lo;
   }
}

double GKOps::getHiRadialField( const GKState& a_state ) const
{
   double value = 0.;

   if ( m_consistent_potential_bcs ) {
      const ScalarPtrVect& scalars( a_state.dataScalar() );
      int scalar_component = a_state.getScalarComponent("Er_boundary");
      const Vector<Real>& scalar_data = scalars[scalar_component]->data();
      value = scalar_data[1];
   }

   return value;
}

void GKOps::setHiRadialField( GKState&      a_state,
                              const double  a_Er_hi ) const
{
   if ( m_consistent_potential_bcs ) {
      ScalarPtrVect& scalars( a_state.dataScalar() );
      int scalar_component = a_state.getScalarComponent("Er_boundary");
      Vector<Real>& scalar_data = scalars[scalar_component]->data();
      scalar_data[1] = a_Er_hi;
   }
}

bool GKOps::stateIncludesPotential( const GKState& a_state ) const
{
   bool includes_potential = false;

   const CFG::FluidSpeciesPtrVect fluid_species( a_state.dataFluid() );
   for (int s(0); s<fluid_species.size(); s++) {
      if (fluid_species[s]->name() == "Potentials") {
         includes_potential = true;
         break;
      }
   }

   return includes_potential;
}

#include "NamespaceFooter.H"

