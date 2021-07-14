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
#include "EFieldAmpere.H"
#include "EFieldSelfConsistentBC.H"
#include "FluxSurface.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"


void GKOps::define( const GKState& a_state )
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

   ParmParse ppgksys( "gksystem" );
   parseParameters( ppgksys );

   m_units = new GKUnits( ppgksys );

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
      m_VorticityOp = new CFG::VorticityOp(pp, m_phase_geometry->magGeom(), larmor, false);
   }
   else {
      m_VorticityOp = NULL;
   }
   
   CFG::IntVect phi_ghost_vect = m_fluid_ghosts*CFG::IntVect::Unit;
   m_phi.define( m_phase_geometry->magGeom().gridsFull(), 1, phi_ghost_vect );

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
         const CFG::LevelData<CFG::FArrayBox>& phi = fluids[a_state_phys.getFluidComponent("potential")]->cell_var("potential");
         
         for (CFG::DataIterator dit(phi.dataIterator()); dit.ok(); ++dit) {
            m_phi[dit].copy(phi[dit]);
         }
      }
   }

   // Get pointers to the primary electric field and one used internally by implicitOpImEx
   string field_type;
   CFG::IntVect ghostVect = (m_phase_geometry->secondOrder()) ? CFG::IntVect::Zero : CFG::IntVect::Unit;
   if ( m_ampere_law ) {
      field_type = "EFieldAmpere";
      m_E_field          = new CFG::EFieldAmpere("", "EField", m_phase_geometry->magGeom(), ghostVect);
      m_E_field_ImOpImEx = new CFG::EFieldAmpere("", "EField", m_phase_geometry->magGeom(), ghostVect);

   }
   else if ( m_consistent_potential_bcs ) {
      field_type = "EFieldSelfConsistentBC";
      m_E_field          = new CFG::EFieldSelfConsistentBC("", "EField", m_phase_geometry->magGeom(), ghostVect);
      m_E_field_ImOpImEx = new CFG::EFieldSelfConsistentBC("", "EField", m_phase_geometry->magGeom(), ghostVect);
   }
   else {
      field_type = "EField";
      m_E_field          = new CFG::EField("", "EField", m_phase_geometry->magGeom(), ghostVect);
      m_E_field_ImOpImEx = new CFG::EField("", "EField", m_phase_geometry->magGeom(), ghostVect);
   }

   m_E_field->define( m_units->larmorNumber(),
                      m_units->debyeNumber(),
                      initial_kinetic_species,
                      m_boltzmann_electron,
                      m_fixed_efield,
                      m_phase_geometry->divFreeVelocity(),
                      a_cur_step );

   m_E_field_ImOpImEx->define( m_units->larmorNumber(),
                               m_units->debyeNumber(),
                               initial_kinetic_species,
                               m_boltzmann_electron,
                               m_fixed_efield,
                               m_phase_geometry->divFreeVelocity(),
                               -1 );

   if ( a_cur_step == 0 && m_ampere_law ) {
      // Initialize the Ampere model fields unless this is a restart, in which case these fields are set
      // by GKSystem::readCheckpointFile()

      m_E_field->fillInternalGhosts(m_phi);
      CFG::LevelData<CFG::FArrayBox>& Er_average_cell = (static_cast<CFG::EFieldAmpere*>(m_E_field))->getErAverageCell();
      m_E_field->computeField(m_phi, Er_average_cell );
      CFG::LevelData<CFG::FluxBox>& Er_average_face = (static_cast<CFG::EFieldAmpere*>(m_E_field))->getErAverageFace();
      m_E_field->computeField(m_phi, Er_average_face );
   }
   
   CH_assert( m_phase_geometry != NULL );

   // If the potential is a state variable or the old vorticity model is being used or the efield
   // is fixed, then something else is controlling the potential, so don't recompute it here.
   bool compute_potential = !(m_state_includes_potential || m_old_vorticity_model || m_fixed_efield);

   m_E_field->computeEField( m_Y,
                             kinetic_species,
                             fluid_species,
                             scalars,
                             m_phi,
                             m_boundary_conditions->getEllipticOpBC(),
                             compute_potential,
                             true,
                             true );

   //Improve Er field calculation to take into accout the dealigment between the grid and magnetic surfaces
   //Should not be used with poloidal variations: FIX LATER!!!

   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   const CFG::MagCoordSys& coords = *mag_geom.getCoordSys();
   if ( (typeid(coords) == typeid(CFG::SingleNullCoordSys)) && (m_dealignment_corrections)) {
      mag_geom.interpolateErFromMagFS(m_E_field->getFaceCenteredField(), m_E_field->getCellCenteredField());
   }

   if ( m_old_vorticity_model ) {
      Vector<Real> scalar_data;
      if (m_consistent_potential_bcs) {
	scalar_data = a_state_phys.dataScalar()[m_Y.getScalarComponent("Er_boundary")]->data();
      }
      m_VorticityOp->initializeOldModel(a_state_phys.dataKinetic(), scalar_data);
   }
   
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
   if (m_E_field_ImOpImEx) delete m_E_field_ImOpImEx;
   if (m_E_field) delete m_E_field;
   if (m_boltzmann_electron) delete m_boltzmann_electron;
   if (m_VorticityOp) delete m_VorticityOp;
   delete m_collisions;
   delete m_fluidOp;
   delete m_scalarOp;
   if (m_transport) delete m_transport;
   if (m_neutrals)  delete m_neutrals;
   if (m_vlasov) delete m_vlasov;
   for (int ihist(0); ihist<m_hist_count; ihist++) {
     FieldHist *save_hist = &m_fieldHistLists[ihist];
     delete save_hist->fieldvals;
     delete save_hist->timevals;
     delete save_hist->timestep;
   }
   delete m_initial_conditions;
   delete m_boundary_conditions;
   delete m_units;

   m_kinetic_species_phys.clear();
   m_fluid_species_phys.clear();

   GyroaverageOperatorFactory::deleteOps(m_gyroavg_ops);
}

Real GKOps::stableDtExpl( const GKState& a_state, const int a_step_number )
{
   Real dt_stable( DBL_MAX );
   if (m_dt_vlasov_expl > 0.) {
     dt_stable = Min( dt_stable, m_dt_vlasov_expl );
   }
   if (m_dt_collisions_expl > 0.) {
      dt_stable = Min( dt_stable, m_dt_collisions_expl );
   }
   dt_stable = Min( dt_stable, m_dt_fluids_expl );
   dt_stable = Min( dt_stable, m_dt_scalars_expl );
   if (m_transport_model_on) {
      dt_stable = Min( dt_stable, m_dt_transport_expl );
   }
   if (m_neutrals_model_on) {
      dt_stable = Min( dt_stable, m_dt_neutrals_expl );
   }
   return dt_stable;
}

Real GKOps::stableDtImEx( const GKState& a_state, const int a_step_number )
{
   Real dt_stable( DBL_MAX );
   if (m_dt_vlasov_imex > 0.) {
     dt_stable = Min( dt_stable, m_dt_vlasov_imex );
   }
   if (m_dt_collisions_imex > 0.) {
      dt_stable = Min( dt_stable, m_dt_collisions_imex );
   }
   dt_stable = Min( dt_stable, m_dt_fluids_imex );
   dt_stable = Min( dt_stable, m_dt_scalars_imex );
   if (m_transport_model_on) {
      dt_stable = Min( dt_stable, m_dt_transport_imex );
   }
   if (m_neutrals_model_on) {
      dt_stable = Min( dt_stable, m_dt_neutrals_imex );
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

   if(!m_zero_efield) {
      updatePhysicalSpeciesVector( a_state_comp.dataFluid(), a_time);
      updatePhysicalSpeciesVector( a_state_comp.dataKinetic(), a_time, false);
      setElectricField( a_state_comp, a_step, m_phi, *m_E_field );
   }

   const KineticSpeciesPtrVect& soln_comp( a_state_comp.dataKinetic() );
   const KineticSpeciesPtrVect& soln_phys( a_state_phys.dataKinetic() );
   m_collisions->preTimeStep( soln_comp, a_time, soln_phys );
   
   /* compute dt for each physics */

   CH_START(t_compute_dt);
   if ( m_vlasov && m_vlasov->computeCfl()) {
      m_dt_vlasov_expl = m_vlasov->computeDtExplicitTI(*m_E_field, soln_comp, a_time);
      m_dt_vlasov_imex = m_vlasov->computeDtImExTI(*m_E_field, soln_comp, a_time);
      m_time_scale_vlasov = m_vlasov->computeTimeScale(*m_E_field, soln_comp, a_time);
   }
   else {
      m_dt_vlasov_expl = -1;
      m_dt_vlasov_imex = -1;
      m_time_scale_vlasov = -1;
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
      const CFG::LevelData<CFG::FArrayBox>& phi = fluids[a_state.getFluidComponent("potential")]->cell_var("potential");
         
      for (CFG::DataIterator dit(phi.dataIterator()); dit.ok(); ++dit) {
         m_phi[dit].copy(phi[dit]);
      }
   }

   if ( m_ampere_post_step_update && m_ampere_law ) {

      const ScalarPtrVect& scalars( a_state.dataScalar() );
      const Vector<Real>& scalar_data = scalars[a_state.getScalarComponent("Er_boundary")]->data();

      const CFG::FluidSpeciesPtrVect& fluids( a_state.dataFluid() );
      const CFG::LevelData<CFG::FArrayBox>& Er = fluids[a_state.getFluidComponent("Er_flux_surfaces")]->cell_var("Er_flux_surfaces");

      (static_cast<CFG::EFieldAmpere*>(m_E_field))->updateAverageAndPerturbation(scalar_data, Er, m_kinetic_species_phys, a_time);
   }

   if ( m_old_vorticity_model ) {

      updatePhysicalSpeciesVector(a_state.dataKinetic(), a_time);
      Vector<Real> scalar_data;
      if (m_consistent_potential_bcs) {
        scalar_data = a_state.dataScalar()[m_Y.getScalarComponent("Er_boundary")]->data();
      }
      m_VorticityOp->updatePotentialOldModel(m_phi, *m_E_field, m_kinetic_species_phys, scalar_data, a_dt, a_time);
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
   if (a_stage && !m_zero_efield && !m_skip_efield_stage_update ) {
      setElectricField( a_state_comp, a_step, m_phi, *m_E_field );
   }
   fillPhysicalKineticSpeciesGhostCells(a_time);

   const KineticSpeciesPtrVect& soln( a_state_comp.dataKinetic() );
   m_collisions->postTimeStage( soln, a_time, a_stage );

   if ( !m_ampere_post_step_update && m_ampere_law ) {

      const ScalarPtrVect& scalars( a_state_comp.dataScalar() );
      const Vector<Real>& scalar_data = scalars[a_state_comp.getScalarComponent("Er_boundary")]->data();

      const CFG::FluidSpeciesPtrVect& fluids( a_state_comp.dataFluid() );
      const CFG::LevelData<CFG::FArrayBox>& Er = fluids[a_state_comp.getFluidComponent("Er_flux_surfaces")]->cell_var("Er_flux_surfaces");

      (static_cast<CFG::EFieldAmpere*>(m_E_field))->updateAverageAndPerturbation(scalar_data, Er, m_kinetic_species_phys, a_time);
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
 
  if (a_chkpt == chkpt_ti_init) {

    updatePhysicalSpeciesVector( a_state.dataFluid(), a_time );

    CH_START(t_update_phys_vector);
    updatePhysicalSpeciesVector( a_state.dataKinetic(), a_time );
    CH_STOP(t_update_phys_vector);

    CH_START(t_fluid_preop_eval);
    m_fluidOp->preSolutionOpEval(m_kinetic_species_phys, a_state.dataFluid(), a_state.dataScalar(), *m_E_field, a_time);
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
                                    *m_E_field, 
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
      m_fluidOp->preSolutionOpEval(m_kinetic_species_phys, a_state.dataFluid(), a_state.dataScalar(), *m_E_field, a_time);
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
  
  bool  update_kinetic_phys,
        update_fluid_phys,
        compute_fluid_pre_op;

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
    m_fluidOp->preOpEval(m_kinetic_species_phys, a_state.dataFluid(), a_state.dataScalar(), *m_E_field, a_time);
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
                       a_state_comp.dataScalar(), *m_E_field, unsplit, a_time );

   applyScalarOperator( a_rhs, m_kinetic_species_phys, m_fluid_species_phys,
                        a_state_comp.dataScalar(), *m_E_field, unsplit, true, a_time );
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
                       a_state_comp.dataScalar(), *m_E_field, imex_exp, a_time );

   applyScalarOperator( a_rhs, m_kinetic_species_phys, m_fluid_species_phys,
                        a_state_comp.dataScalar(), *m_E_field, imex_exp, true, a_time );
}

void GKOps::implicitOpImEx( GKRHSData&      a_rhs,
                            const Real      a_time,
                            const GKState&  a_state_comp )
{
   CH_TIME("GKOps::implicitOpImEx");

   CFG::LevelData<CFG::FArrayBox> phi( m_phase_geometry->magGeom().gridsFull(), 1, 4*CFG::IntVect::Unit );

   if(!m_zero_efield) {
      // The last argument is for whether E_field will be injected or not;
      // If m_step_const_kin_coeff_fluidop = true, it means that kinietc distribution is not updated;
      // so we do not need to inject E_field
      setElectricField( a_state_comp, 
                        -1, 
                        phi, 
                        *m_E_field_ImOpImEx, 
                        !m_step_const_kin_coeff_fluidop );
   }

   a_rhs.zero();

   applyVlasovOperator( a_rhs, m_kinetic_species_phys, imex_imp, a_time );

   if (m_transport_model_on && (m_transport_imex_implicit) ) {
      applyTransportOperator( a_rhs.dataKinetic(), m_kinetic_species_phys, a_time );
   }

   applyCollisionOperator( a_rhs.dataKinetic(), a_state_comp.dataKinetic(), imex_imp, a_time );

   if (m_neutrals_model_on && (m_neutrals_imex_implicit) ) {
      applyNeutralsOperator( a_rhs.dataKinetic(), m_kinetic_species_phys,
                             m_fluid_species_phys, a_time );
   }
    
   applyFluidOperator( a_rhs.dataFluid(), m_kinetic_species_phys, m_fluid_species_phys,
                       a_state_comp.dataScalar(), *m_E_field_ImOpImEx, imex_imp, a_time );

   applyScalarOperator( a_rhs, m_kinetic_species_phys, m_fluid_species_phys,
                        a_state_comp.dataScalar(), *m_E_field_ImOpImEx, imex_imp, true, a_time );
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
                       *m_E_field, unsplit, a_time );

   applyScalarOperator( a_rhs, m_kinetic_species_phys, m_fluid_species_phys,
                        a_state.dataScalar(), *m_E_field, unsplit, true, a_time );
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
                                  bool                                              a_im )
{
  a_pc.clear();
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
                                    a_im );
  }

  m_collisions->defineMultiPhysicsPC( a_pc, 
                                      a_dof_list,
                                      a_state.dataKinetic(),
                                      global_dof->dataKinetic(),
                                      a_Y, 
                                      a_sys, 
                                      a_out_string, 
                                      a_opt_string, 
                                      a_im );

  m_fluidOp->defineMultiPhysicsPC(  a_pc,
                                    a_dof_list,
                                    a_state.dataFluid(),
                                    global_dof->dataFluid(),
                                    a_Y,
                                    a_sys,
                                    a_out_string,
                                    a_opt_string,
                                    a_im );

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
                                    m_phi,
                                    *m_E_field,
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

void GKOps::solveFluidOpPCImEx( GKRHSData&     a_rhs,
                                const GKState& a_state,
                                int            a_idx )
{
   CH_TIME("GKOps::solvePCImEx");
   m_fluidOp->solvePCImEx(  a_rhs.dataFluid(), 
                            a_state.dataKinetic(), 
                            a_state.dataFluid(), 
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
void GKOps::setElectricField( const GKState&                   a_state_comp,
                              const int                        a_step,
                              CFG::LevelData<CFG::FArrayBox>&  a_phi,
                              CFG::EField&                     a_E_field,
                              const bool                       a_inject_field) const
{
   CH_TIME("GKOps::setElectricField");
 
   if ( m_state_includes_potential ) {
      const CFG::FluidSpeciesPtrVect& fluids( a_state_comp.dataFluid() );
      const CFG::LevelData<CFG::FArrayBox>& phi = fluids[a_state_comp.getFluidComponent("potential")]->cell_var("potential");
         
      for (CFG::DataIterator dit(phi.dataIterator()); dit.ok(); ++dit) {
         a_phi[dit].copy(phi[dit]);
      }
   }

   // If the potential is a state variable or the old vorticity model is being used or the
   // efield is fixed, then something else is controlling the potential, so don't recompute it here.
   bool compute_potential = !(m_state_includes_potential || m_old_vorticity_model || m_fixed_efield);

   a_E_field.computeEField( m_Y,
                            m_kinetic_species_phys,
                            m_fluid_species_phys,
                            a_state_comp.dataScalar(),
                            a_phi,
                            m_boundary_conditions->getEllipticOpBC(),
                            compute_potential,
                            a_inject_field,
                            false );

}


inline
void GKOps::updatePhysicalSpeciesVector( const KineticSpeciesPtrVect&  a_species_comp,
                                         const double                  a_time,
                                         const bool                    a_fill_ghost_cells )
{
   CH_TIMERS("GKOps::updatePhysicalSpeciesVector (kinetic)");

   divideJ( a_species_comp, m_kinetic_species_phys );

   if (a_fill_ghost_cells) {
     m_boundary_conditions->fillGhostCells( m_kinetic_species_phys,
                                            m_phi,
                                            m_E_field->getInjectedField(),
                                            a_time );
   }
}


inline
void GKOps::fillPhysicalKineticSpeciesGhostCells( const double a_time )
{
   CH_TIMERS("GKOps::fillPhysicalKineticSpeciesGhostCells");
   m_boundary_conditions->fillGhostCells( m_kinetic_species_phys,
                                          m_phi,
                                          m_E_field->getInjectedField(),
                                          a_time );
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
        m_vlasov->accumulateRHS(a_rhs, a_kinetic_phys, m_phi, *m_E_field, false, a_time );
        m_vlasov->accumulateRHS(a_rhs, a_kinetic_phys, m_phi, *m_E_field, true, a_time );
      } else if (a_op_type == imex_exp) {
        m_vlasov->accumulateRHS(a_rhs, a_kinetic_phys, m_phi, *m_E_field, false, a_time );
      } else if (a_op_type == imex_imp) {
        m_vlasov->accumulateRHS(a_rhs, a_kinetic_phys, m_phi, *m_E_field, true, a_time );
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
                                const CFG::EField&                         a_E_field,
                                const opType                               a_op_type,
                                const Real&                                a_time)
{
   CH_TIME("GKOps::applyFluidOperator");

   m_count_fluids++;
   if (a_op_type == imex_exp) {
     m_fluidOp->accumulateRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field,
                               false, a_time );
   } else if (a_op_type == imex_imp) {
     m_fluidOp->accumulateRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field,
                               true, a_time );
   } else if (a_op_type == unsplit) { 
     m_fluidOp->accumulateRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field,
                               false, a_time );
     m_fluidOp->accumulateRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field,
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
                                 const CFG::EField&               a_E_field,
                                 const opType                     a_op_type,
                                 const bool                       a_recompute_kinetic_terms,
                                 const Real&                      a_time)
{
   CH_TIME("GKOps::applyScalarOperator");

   m_count_scalars++;
   //m_count_fluids++;
   if (a_op_type == imex_exp) {
     m_scalarOp->accumulateRHS( a_rhs, a_kinetic_species, a_fluid_species, a_scalars, a_E_field, 
                                false, a_recompute_kinetic_terms, a_time );
   } else if (a_op_type == imex_imp) {
     m_scalarOp->accumulateRHS( a_rhs, a_kinetic_species, a_fluid_species, a_scalars, a_E_field, 
                                true, a_recompute_kinetic_terms, a_time );
   } else if (a_op_type == unsplit) {
     m_scalarOp->accumulateRHS( a_rhs, a_kinetic_species, a_fluid_species, a_scalars, a_E_field, 
                                false, a_recompute_kinetic_terms, a_time );
     m_scalarOp->accumulateRHS( a_rhs, a_kinetic_species, a_fluid_species, a_scalars, a_E_field, 
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

   CFG::EFieldAmpere& efield = static_cast<CFG::EFieldAmpere&>(*m_E_field);
   efield.setErAverage(Er_average_face);
}

void GKOps::setErAverage( const LevelData<FArrayBox>& Er_cell_injected)
{
   CFG::LevelData<CFG::FArrayBox> Er_average_cell(m_phase_geometry->magGeom().gridsFull(), 3, CFG::IntVect::Unit);
   m_phase_geometry->projectPhaseToConfiguration(Er_cell_injected, Er_average_cell);

   CFG::EFieldAmpere& efield = static_cast<CFG::EFieldAmpere&>(*m_E_field);
   efield.setErAverage(Er_average_cell);
}

void GKOps::setETilde( const LevelData<FluxBox>& E_tilde_face_injected )
{
   CFG::LevelData<CFG::FluxBox> E_tilde_face(m_phase_geometry->magGeom().gridsFull(), 3, CFG::IntVect::Unit);
   for (CFG::DataIterator dit(E_tilde_face.disjointBoxLayout()); dit.ok(); ++dit) {
      E_tilde_face[dit].setVal(0.);
   }
   m_phase_geometry->projectPhaseToConfiguration(E_tilde_face_injected, E_tilde_face);

   CFG::EFieldAmpere& efield = static_cast<CFG::EFieldAmpere&>(*m_E_field);
   efield.setETilde(E_tilde_face);
}

void GKOps::setETilde( const LevelData<FArrayBox>& E_tilde_cell_injected)
{
   CFG::LevelData<CFG::FArrayBox> E_tilde_cell(m_phase_geometry->magGeom().gridsFull(), 3, CFG::IntVect::Unit);
   m_phase_geometry->projectPhaseToConfiguration(E_tilde_cell_injected, E_tilde_cell);

   CFG::EFieldAmpere& efield = static_cast<CFG::EFieldAmpere&>(*m_E_field);
   efield.setETilde(E_tilde_cell);
}

void GKOps::setPhi( const LevelData<FArrayBox>& phi_injected)
{
   CFG::LevelData<CFG::FArrayBox> phi_tmp(m_phase_geometry->magGeom().gridsFull(), 1, CFG::IntVect::Zero);
   m_phase_geometry->projectPhaseToConfiguration(phi_injected, phi_tmp);
   
   for (CFG::DataIterator dit(m_phi.disjointBoxLayout()); dit.ok(); ++dit) {
      m_phi[dit].copy(phi_tmp[dit]);
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


void GKOps::parseParameters( ParmParse& a_ppgksys )
{
   if (a_ppgksys.contains("fixed_efield")) {
      a_ppgksys.query( "fixed_efield", m_fixed_efield );
   }
   else {
      m_fixed_efield = false;
   }
   
   if (a_ppgksys.contains("zero_efield")) {
      a_ppgksys.query( "zero_efield", m_zero_efield );
   }

   a_ppgksys.query( "imex_transport_implicit", m_transport_imex_implicit );
   a_ppgksys.query( "imex_neutrals_implicit", m_neutrals_imex_implicit );

   if (a_ppgksys.contains("old_vorticity_model")) {
      a_ppgksys.get("old_vorticity_model", m_old_vorticity_model);
   }
   else {
      m_old_vorticity_model = false;
   }

   if (a_ppgksys.contains("skip_efield_stage_update")) {
      a_ppgksys.get("skip_efield_stage_update", m_skip_efield_stage_update);
   }
   else {
      m_skip_efield_stage_update = false;
   }

   if (a_ppgksys.contains("consistent_potential_bcs")) {
      a_ppgksys.query( "consistent_potential_bcs", m_consistent_potential_bcs );
   }
   else {
      m_consistent_potential_bcs = false;
   }

  if (a_ppgksys.contains("use_vlasov_ampere_bcs")) {
      a_ppgksys.query( "use_vlasov_ampere_bcs", m_use_vlasov_ampere_bcs );
   }
   else {
      m_use_vlasov_ampere_bcs = true;
   }
   
   if (a_ppgksys.contains("extrapolated_sol_efield")) {
     a_ppgksys.query( "extrapolated_sol_efield", m_Esol_extrapolation );
   }
   else {
     m_Esol_extrapolation = true;
   }

   if (a_ppgksys.contains("efield_dealignment_corrections")) {
     a_ppgksys.query( "efield_dealignment_corrections", m_dealignment_corrections );
   }
   else {
     m_dealignment_corrections = false;
   }

   if (a_ppgksys.contains("ampere_law")) {
      a_ppgksys.query( "ampere_law", m_ampere_law );

      if (m_ampere_law) {
         m_consistent_potential_bcs = true;

         if (a_ppgksys.contains("ampere_cold_electrons")) {
            a_ppgksys.query( "ampere_cold_electrons", m_ampere_cold_electrons );
         }
         else {
            m_ampere_cold_electrons = false;
         }
      }
      else {
         m_ampere_cold_electrons = false;
      }
   }
   else {
      m_ampere_law = false;
      m_ampere_cold_electrons = false;
   }

   if (a_ppgksys.contains("ampere_post_step_update")) {
      a_ppgksys.query( "ampere_post_step_update", m_ampere_post_step_update );
   }
   else {
      m_ampere_post_step_update = false;
   }

   if ( m_fixed_efield && m_ampere_law ) {
      MayDay::Error("GKOps::parseParameters(): Specify either fixed field or ampere law, but not both"); 
   }

   if (a_ppgksys.contains("transport_model_on")) {
      a_ppgksys.get("transport_model_on", m_transport_model_on);
   }
   else {
      m_transport_model_on = false;
   }

   if (a_ppgksys.contains("neutrals_model_on")) {
      a_ppgksys.get("neutrals_model_on", m_neutrals_model_on);
   }
   else {
      m_neutrals_model_on = false;
   }

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

   a_ppgksys.query( "gksystem.enforce_quasineutrality", m_enforce_quasineutrality );

   a_ppgksys.query("step_const_coef", m_step_const_kin_coeff_fluidop);
   
   
   // Get kinetic ghost layer width
   a_ppgksys.query( "kinetic_ghost_width", m_kinetic_ghosts );

   // Get fluid ghost layer width
   a_ppgksys.query( "fluid_ghost_width", m_fluid_ghosts );
}



void GKOps::plotPotential( const std::string& a_filename,
                           const bool&        a_non_zonal_comp,
                           const double&      a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   if (!a_non_zonal_comp){
      m_phase_geometry->plotConfigurationData( a_filename.c_str(), m_phi, a_time );
   }
   else {
      const PhaseGeom& phase_geometry( *m_phase_geometry );
      const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
      CFG::LevelData<CFG::FArrayBox> phi( mag_geom.gridsFull(), 2, CFG::IntVect::Zero);
      
      CFG::LevelData<CFG::FArrayBox> phi_axisym;
      phi_axisym.define(m_phi);
      
      CFG::FluxSurface flux_surface(mag_geom);
      flux_surface.averageAndSpread(m_phi, phi_axisym);
      for (CFG::DataIterator dit(m_phi.dataIterator()); dit.ok(); ++dit) {
         phi[dit].copy(m_phi[dit],0,0,1);
         phi[dit].copy(m_phi[dit],0,1,1);
         phi[dit].minus(phi_axisym[dit],0,1,1);
      }
      m_phase_geometry->plotConfigurationData( a_filename.c_str(), phi, a_time );
   }
}

void GKOps::plotEField( const std::string& a_filename,
                        const double&      a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );

   // Plot the psi and theta projections of the physical field

   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   const CFG::DisjointBoxLayout& grids = mag_geom.gridsFull();
   CFG::LevelData<CFG::FluxBox> Efield( grids, CFG_DIM, CFG::IntVect::Unit);

#if 1
#if CFG_DIM==2
   m_E_field->computePoloidalField( m_phi, Efield );
   for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
      CFG::FluxBox& this_Efield = Efield[dit];
      mag_geom.getBlockCoordSys(grids[dit]).computePsiThetaProjections(this_Efield);
   }
#endif
#if CFG_DIM==3
   m_E_field->computeField( m_phi, Efield );
#endif
#else
   m_E_field->getFaceCenteredFieldOnCore( m_phi, Efield );
#endif
      
   if ( m_ampere_law ) {
      CFG::LevelData<CFG::FArrayBox>& E_field = m_E_field->getCellCenteredField();
      CFG::LevelData<CFG::FArrayBox> E_field_no_ghosts(grids, CFG_DIM, CFG::IntVect::Zero);
      for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
         E_field_no_ghosts[dit].copy(E_field[dit], grids[dit]);
      }
      phase_geometry.plotConfigurationData( a_filename.c_str(), E_field_no_ghosts, a_time );
   }
   else {
      CFG::LevelData<CFG::FluxBox> E_field_no_ghosts(grids, CFG_DIM, CFG::IntVect::Zero);
      for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
         E_field_no_ghosts[dit].copy(Efield[dit], grids[dit]);
      }
      phase_geometry.plotConfigurationData( a_filename.c_str(), E_field_no_ghosts, a_time );
   }
}

void GKOps::plotExBData(const std::string& a_filename,
                        const double&      a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );

   // Plot the psi and theta projections of the physical field

   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   const CFG::DisjointBoxLayout& grids = mag_geom.gridsFull();

   CFG::LevelData<CFG::FArrayBox> Efield( grids, 3, CFG::IntVect::Unit);
   m_E_field->computeField( m_phi, Efield );
   
   CFG::LevelData<CFG::FArrayBox> ExB_drift( grids, 3, CFG::IntVect::Zero);
   mag_geom.computeEXBDrift(Efield, ExB_drift);
   
   const double larmor_number( m_units->larmorNumber() );
   for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
      ExB_drift[dit].mult(larmor_number);
   }
   
   CFG::LevelData<CFG::FArrayBox> ExB_data( grids, 4, CFG::IntVect::Zero);
   CFG::LevelData<CFG::FArrayBox> tmp( grids, 1, CFG::IntVect::Zero);
   
   // Compute radial component of ExB velocity
   mag_geom.computeRadialProjection(tmp, ExB_drift);
   for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
      ExB_data[dit].copy(tmp[dit],0,0,1);
   }

   // Compute FS averaged poloidal component of ExB velocity
   mag_geom.computePoloidalProjection(tmp, ExB_drift);

   CFG::FluxSurface flux_surface(mag_geom);
   CFG::LevelData<CFG::FArrayBox> tmp_axisym( grids, 1, CFG::IntVect::Zero);
   flux_surface.averageAndSpread(tmp, tmp_axisym);
   
   for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
      ExB_data[dit].copy(tmp_axisym[dit],0,1,1);
   }

   // Compute radial component of the <ExB velocity> shear
   CFG::LevelData<CFG::FArrayBox> ExB_shear( grids, 3, CFG::IntVect::Unit);
   m_E_field->computeField(tmp_axisym, ExB_shear);
   
   mag_geom.computeRadialProjection(tmp, ExB_shear);
   for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
      ExB_data[dit].copy(tmp[dit],0,2,1);
      //Multiply by -1 to compensate for the minus sign
      // in computeField calculation that returns -grad(phi)
      ExB_data[dit].mult(-1.0,2,1);
   }

   // Compute flux surface average of the <ExB> shear
   //CFG::FluxSurface flux_surface(mag_geom);
   CFG::LevelData<CFG::FArrayBox> negative_ExB_shear_axisym( grids, 1, CFG::IntVect::Zero);
   flux_surface.averageAndSpread(tmp, negative_ExB_shear_axisym);
   for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
      ExB_data[dit].copy(negative_ExB_shear_axisym[dit],0,3,1);
      ExB_data[dit].mult(-1.0,3,1);
   }
   
   phase_geometry.plotConfigurationData( a_filename.c_str(), ExB_data, a_time );
}

void GKOps::plotDistributionFunction( const std::string&    a_filename,
                                      const KineticSpecies& a_soln_species,
                                      const double&         a_time ) const
{
   const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& soln_dfn = a_soln_species.distributionFunction();
   LevelData<FArrayBox> dfn( species_geometry.gridsFull(), 1, IntVect::Zero );
   for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      dfn[dit].copy(soln_dfn[dit]);
   }

   species_geometry.divideBStarParallel( dfn );
   species_geometry.plotData( a_filename.c_str(), dfn, a_time );
}

void GKOps::plotFluid( const std::string&       a_filename,
                       const CFG::FluidSpecies& a_fluid_species,
                       const string&            a_var_name,
                       const double&            a_time ) const
{
   CH_TIME("GKOps::plotFluid()");
   m_phase_geometry->plotConfigurationData( a_filename.c_str(), a_fluid_species.cell_var(a_var_name), a_time );
}

void GKOps::plotFluidAtCellFromFaceNorms( const std::string&        a_filename,
                                          const CFG::FluidSpecies&  a_fluid_species,
                                          const string&             a_var_name,
                                          const double&             a_time ) const
{
   CH_TIME("GKOps::plotFluidAtCellFromFaceNorms()");

   // convert contravarient vars on faces to physical vars at cells
   //
   const CFG::IntVect& ghostVect = a_fluid_species.face_var(a_var_name).ghostVect();
   CFG::LevelData<CFG::FArrayBox> this_cell_data(a_fluid_species.face_var(a_var_name).getBoxes(), CFG_DIM, ghostVect);
   a_fluid_species.interpFaceVarToCell(this_cell_data, a_var_name);
   
   //const PhaseGeom& phase_geometry( *m_phase_geometry );
   //const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   mag_geom.convertPhysToContravar(this_cell_data,1); // 1 is for inverse
   
   // plot data
   //
   m_phase_geometry->plotConfigurationData( a_filename.c_str(), this_cell_data, a_time );
   
}

void GKOps::plotFluidAtCellFromEdgeTans( const std::string&        a_filename,
                                         const CFG::FluidSpecies&  a_fluid_species,
                                         const string&             a_var_name,
                                         const double&             a_time ) const
{
   CH_TIME("GKOps::plotFluidAtCellFromEdgeTans()");
   
   // convert covariant vars on edges to physical vars at cells
   //
   const CFG::IntVect& ghostVect = a_fluid_species.edge_var(a_var_name).ghostVect();
   CFG::LevelData<CFG::FArrayBox> this_cell_data(a_fluid_species.edge_var(a_var_name).getBoxes(), CFG_DIM, ghostVect);
   a_fluid_species.interpEdgeVarToCell(this_cell_data, a_var_name);
   
   //const PhaseGeom& phase_geometry( *m_phase_geometry );
   //const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   mag_geom.convertPhysToCovar(this_cell_data,1); // 1 is for inverse

   m_phase_geometry->plotConfigurationData( a_filename.c_str(), this_cell_data, a_time );

}

void GKOps::plotFluidOpMember( const std::string&        a_filename,
                               const CFG::FluidSpecies&  a_fluid_species,
                               const std::string&        a_member_var_name,
                               const double              a_time ) const

{
   CH_TIME("GKOps::plotFluidOpMember()");
   //const CFG::IntVect& ghostVect = 1*CFG::IntVect::Unit;
   //const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   CFG::LevelData<CFG::FArrayBox> this_cell_data;
   m_fluidOp->getMemberVarForPlotting(this_cell_data, a_fluid_species, a_member_var_name);
   m_phase_geometry->plotConfigurationData( a_filename.c_str(), this_cell_data, a_time );
  
}

void GKOps::plotBStarParallel( const std::string&    a_filename,
                               const KineticSpecies& a_soln_species,
                               const double&         a_time ) const
{
   const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();

   LevelData<FArrayBox> BStarParallel( species_geometry.gridsFull(), 1, IntVect::Zero );
   species_geometry.getBStarParallel(BStarParallel);

   species_geometry.plotData( a_filename.c_str(), BStarParallel, a_time );
}

void GKOps::plotDeltaF( const std::string&    a_filename,
                        const KineticSpecies& a_soln_species,
                        const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> density( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.numberDensity( density );

   CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.ParallelMomentum( ParallelMom );

   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      ParallelMom[dit].divide(density[dit]);
   }
   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.pressureMoment(pressure, ParallelMom);

   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      pressure[dit].divide(density[dit]);
   }

   const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();
   LevelData<FArrayBox> deltaF( species_geometry.gridsFull(), 1, IntVect::Zero );
   DeltaFKernel<FArrayBox> DeltaF_Kernel(density, pressure, ParallelMom);
   DeltaF_Kernel.eval(deltaF, a_soln_species);

   species_geometry.plotData( a_filename.c_str(), deltaF, a_time );

}


void GKOps::plotDistributionFunctionAtMu( const std::string&    a_filename,
                                          const KineticSpecies& a_soln_species,
                                          const int             a_mu,
                                          const double&         a_time ) const
{
   const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& soln_dfn = a_soln_species.distributionFunction();
   LevelData<FArrayBox> dfn( species_geometry.gridsFull(), 1, IntVect::Zero );
   for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      dfn[dit].copy(soln_dfn[dit]);
   }

   species_geometry.divideBStarParallel( dfn );
   species_geometry.plotAtMuIndex( a_filename.c_str(), a_mu, soln_dfn, a_time );
}


void GKOps::plotVParallelTheta( const std::string&    a_filename,
                                const KineticSpecies& a_soln_species,
                                const int             a_radial_index,
                                const int             a_toroidal_index,
                                const int             a_mu_index,
                                const double&         a_time ) const
{
   const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& soln_dfn = a_soln_species.distributionFunction();
   LevelData<FArrayBox> dfn( species_geometry.gridsFull(), 1, IntVect::Zero );
   for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      dfn[dit].copy(soln_dfn[dit]);
   }

   species_geometry.divideBStarParallel( dfn );
   species_geometry.plotVParPoloidalData( a_filename.c_str(), a_radial_index,
                                          a_toroidal_index, a_mu_index, dfn, a_time );
}

void GKOps::plotbfVParallelTheta( const std::string&    a_filename,
                                const KineticSpecies& a_soln_species,
                                const int             a_radial_index,
                                const int             a_toroidal_index,
                                const int             a_mu_index,
                                const double&         a_time ) const
{
   const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& soln_dfn = a_soln_species.distributionFunction();
   LevelData<FArrayBox> dfn( species_geometry.gridsFull(), 1, IntVect::Zero );
   for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      dfn[dit].copy(soln_dfn[dit]);
   }

   // Same as plotVParallelTheta but not dividing by BStarParallel this time
   species_geometry.plotVParPoloidalData( a_filename.c_str(), a_radial_index,
                                          a_toroidal_index, a_mu_index, dfn, a_time );
}


void GKOps::plotRTheta( const std::string&    a_filename,
                        const KineticSpecies& a_soln_species,
                        const int             a_vpar_index,
                        const int             a_mu_index,
                        const double&         a_time ) const
{
   const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& soln_dfn = a_soln_species.distributionFunction();
   LevelData<FArrayBox> dfn( species_geometry.gridsFull(), 1, IntVect::Zero );
   for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      dfn[dit].copy(soln_dfn[dit]);
   }

   VEL::IntVect vspace_index(a_vpar_index, a_mu_index);

   species_geometry.divideBStarParallel( dfn );
   species_geometry.plotAtVelocityIndex( a_filename.c_str(), vspace_index, dfn, a_time );
}



void GKOps::plotVParallelMu( const std::string&    a_filename,
                             const KineticSpecies& a_soln_species,
                             const int             a_radial_index,
                             const int             a_poloidal_index,
                             const int             a_toroidal_index,
                             const double&         a_time ) const
{
   const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& soln_dfn = a_soln_species.distributionFunction();
   LevelData<FArrayBox> dfn( species_geometry.gridsFull(), 1, IntVect::Zero );
   for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      dfn[dit].copy(soln_dfn[dit]);
   }

#if CFG_DIM==3
   CFG::IntVect cspace_index(a_radial_index, a_poloidal_index, a_toroidal_index);
#else
   CFG::IntVect cspace_index(a_radial_index, a_poloidal_index);
#endif
   species_geometry.divideBStarParallel( dfn );
   species_geometry.plotAtConfigurationIndex( a_filename.c_str(),
                                              cspace_index, dfn, a_time );
}


void GKOps::plotChargeDensity( const std::string&    a_filename,
                              const KineticSpecies& a_soln_species,
                              const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   CFG::LevelData<CFG::FArrayBox> charge_density( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.chargeDensity( charge_density );
   
   phase_geometry.plotConfigurationData( a_filename.c_str(), charge_density, a_time );
}


void GKOps::plotChargeDensity( const std::string&     a_filename,
                               const KineticSpeciesPtrVect& a_species,
                               const double&          a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   CFG::LevelData<CFG::FArrayBox> charge_density( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   CFG::LevelData<CFG::FArrayBox> species_charge_density( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );

   for (CFG::DataIterator dit(charge_density.dataIterator()); dit.ok(); ++dit) {
      charge_density[dit].setVal(0.);
   }

   for (int species(0); species<a_species.size(); species++) {
      KineticSpecies& soln_species( *(a_species[species]) );
      soln_species.chargeDensity(species_charge_density);

      for (CFG::DataIterator dit(charge_density.dataIterator()); dit.ok(); ++dit) {
         charge_density[dit] += species_charge_density[dit];
      }
   }

   if (m_boltzmann_electron) {
      const CFG::LevelData<CFG::FArrayBox>& ne = m_boltzmann_electron->numberDensity();

      for (CFG::DataIterator dit(charge_density.dataIterator()); dit.ok(); ++dit) {
         charge_density[dit] -= ne[dit];
      }
   }

   phase_geometry.plotConfigurationData( a_filename.c_str(), charge_density, a_time );
}


void GKOps::plotMomentum( const std::string&    a_filename,
                          const KineticSpecies& a_soln_species,
                          const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const double larmor_number( m_units->larmorNumber() );
   LevelData<FluxBox> E_field_tmp;

   phase_geometry.injectConfigurationToPhase( m_E_field->getFaceCenteredField(),
   					      m_E_field->getCellCenteredField(),
					      E_field_tmp );

   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> momentum( mag_geom.gridsFull(), CFG_DIM, CFG::IntVect::Zero );
   a_soln_species.momentumDensity( momentum, E_field_tmp, larmor_number );

   phase_geometry.plotConfigurationData( a_filename.c_str(), momentum, a_time );
}


void GKOps::plotParallelMomentum( const std::string&    a_filename,
                                  const KineticSpecies& a_soln_species,
                                  const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   CFG::LevelData<CFG::FArrayBox> parallel_mom( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.ParallelMomentum(parallel_mom);
   phase_geometry.plotConfigurationData( a_filename.c_str(), parallel_mom, a_time );
}


void GKOps::plotPoloidalMomentum( const std::string&    a_filename,
                                  const KineticSpecies& a_soln_species,
                                  const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const double larmor_number( m_units->larmorNumber() );
   LevelData<FluxBox> E_field_tmp;

   phase_geometry.injectConfigurationToPhase( m_E_field->getFaceCenteredField(),
                                              m_E_field->getCellCenteredField(),
                                              E_field_tmp );

   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> poloidal_mom( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.PoloidalMomentum( poloidal_mom, E_field_tmp, larmor_number );

   phase_geometry.plotConfigurationData( a_filename.c_str(), poloidal_mom, a_time );
}

void GKOps::plotParallelVelocity( const std::string&    a_filename,
                                  const KineticSpecies& a_soln_species,
                                  const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> parallel_velocity( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.ParallelVelocity( parallel_velocity );

   phase_geometry.plotConfigurationData( a_filename.c_str(), parallel_velocity, a_time );
}


void GKOps::plotEnergyDensity(const std::string&    a_filename,
                              const KineticSpecies& a_soln_species,
                              const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> energy_density( mag_geom.gridsFull(), 3, CFG::IntVect::Zero );
   
   CFG::LevelData<CFG::FArrayBox> energy_density_loc( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.energyMoment( energy_density_loc, m_phi );

   for (CFG::DataIterator dit(energy_density.dataIterator()); dit.ok(); ++dit) {
      energy_density[dit].copy(energy_density_loc[dit], 0, 0, 1);
   }

   // Get flux-surface shell-averaged and shell-integrated components
   CFG::FluxSurface flux_surface(mag_geom);
   CFG::LevelData<CFG::FArrayBox> energy_density_tmp( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   
   flux_surface.averageAndSpread(energy_density_loc, energy_density_tmp);
   for (CFG::DataIterator dit(energy_density.dataIterator()); dit.ok(); ++dit) {
      energy_density[dit].copy(energy_density_tmp[dit], 0, 1, 1);
   }
   
   flux_surface.integrateAndSpread(energy_density_loc, energy_density_tmp);
   for (CFG::DataIterator dit(energy_density.dataIterator()); dit.ok(); ++dit) {
      energy_density[dit].copy(energy_density_tmp[dit], 0, 2, 1);
   }
   
   phase_geometry.plotConfigurationData( a_filename.c_str(), energy_density, a_time );
}

void GKOps::plotKineticEnergyDensity(const std::string&    a_filename,
                                     const KineticSpecies& a_soln_species,
                                     const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> energy_density( mag_geom.gridsFull(), 3, CFG::IntVect::Zero );
   
   CFG::LevelData<CFG::FArrayBox> energy_density_loc( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.kineticEnergyMoment( energy_density_loc );

   for (CFG::DataIterator dit(energy_density.dataIterator()); dit.ok(); ++dit) {
      energy_density[dit].copy(energy_density_loc[dit], 0, 0, 1);
   }

   // Get flux-surface shell-averaged and shell-integrated components
   CFG::FluxSurface flux_surface(mag_geom);
   CFG::LevelData<CFG::FArrayBox> energy_density_tmp( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   
   flux_surface.averageAndSpread(energy_density_loc, energy_density_tmp);
   for (CFG::DataIterator dit(energy_density.dataIterator()); dit.ok(); ++dit) {
      energy_density[dit].copy(energy_density_tmp[dit], 0, 1, 1);
   }
   
   flux_surface.integrateAndSpread(energy_density_loc, energy_density_tmp);
   for (CFG::DataIterator dit(energy_density.dataIterator()); dit.ok(); ++dit) {
      energy_density[dit].copy(energy_density_tmp[dit], 0, 2, 1);
   }
   
   phase_geometry.plotConfigurationData( a_filename.c_str(), energy_density, a_time );
}

void GKOps::plotParallelEnergyDensity(const std::string&    a_filename,
                                      const KineticSpecies& a_soln_species,
                                      const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> energy_density( mag_geom.gridsFull(), 3, CFG::IntVect::Zero );
   
   CFG::LevelData<CFG::FArrayBox> energy_density_loc( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.parallelEnergyDensity( energy_density_loc );

   for (CFG::DataIterator dit(energy_density.dataIterator()); dit.ok(); ++dit) {
      energy_density[dit].copy(energy_density_loc[dit], 0, 0, 1);
   }

   // Get flux-surface shell-averaged and shell-integrated components
   CFG::FluxSurface flux_surface(mag_geom);
   CFG::LevelData<CFG::FArrayBox> energy_density_tmp( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   
   flux_surface.averageAndSpread(energy_density_loc, energy_density_tmp);
   for (CFG::DataIterator dit(energy_density.dataIterator()); dit.ok(); ++dit) {
      energy_density[dit].copy(energy_density_tmp[dit], 0, 1, 1);
   }
   
   flux_surface.integrateAndSpread(energy_density_loc, energy_density_tmp);
   for (CFG::DataIterator dit(energy_density.dataIterator()); dit.ok(); ++dit) {
      energy_density[dit].copy(energy_density_tmp[dit], 0, 2, 1);
   }
   
   phase_geometry.plotConfigurationData( a_filename.c_str(), energy_density, a_time );
}


void GKOps::plotPerpEnergyDensity(const std::string&    a_filename,
                                  const KineticSpecies& a_soln_species,
                                  const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> energy_density( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.perpEnergyDensity( energy_density );
   
   phase_geometry.plotConfigurationData( a_filename.c_str(), energy_density, a_time );
}


void GKOps::plotPressure( const std::string&    a_filename,
                          const KineticSpecies& a_soln_species,
                          const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.pressureMoment( pressure );

   phase_geometry.plotConfigurationData( a_filename.c_str(), pressure, a_time );
}


void GKOps::plotParallelPressure(const std::string&    a_filename,
                                 const KineticSpecies& a_soln_species,
                                 const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.parallelPressure( pressure );

   phase_geometry.plotConfigurationData( a_filename.c_str(), pressure, a_time );
}


void GKOps::plotPerpPressure(const std::string&    a_filename,
                             const KineticSpecies& a_soln_species,
                             const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.perpPressure( pressure );

   phase_geometry.plotConfigurationData( a_filename.c_str(), pressure, a_time );
}

void GKOps::plotGradPoverN(const std::string&    a_filename,
                           const KineticSpecies& a_soln_species,
                           const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> density( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.numberDensity( density );
   a_soln_species.ParallelMomentum( ParallelMom );

   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      ParallelMom[dit].divide(density[dit]);
   }
   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.pressureMoment(pressure, ParallelMom);

   const CFG::DisjointBoxLayout& cfg_grids = mag_geom.gridsFull();
   CFG::LevelData<CFG::FArrayBox> pressure_grown( cfg_grids, 1, 2*CFG::IntVect::Unit );
   for (CFG::DataIterator dit(pressure.dataIterator()); dit.ok(); ++dit) {
      pressure_grown[dit].copy(pressure[dit]);
      const CFG::MagBlockCoordSys& coord_sys = mag_geom.getBlockCoordSys(cfg_grids[dit]);
      const CFG::ProblemDomain& domain = coord_sys.domain();
      secondOrderCellExtrapAtDomainBdry(pressure_grown[dit], cfg_grids[dit], domain);
   }
   mag_geom.fillInternalGhosts(pressure_grown);

#if CFG_DIM==2

   CFG::LevelData<CFG::FArrayBox> pressure_gradient_mapped( cfg_grids, CFG_DIM, CFG::IntVect::Unit );
   mag_geom.computeMappedPoloidalGradientWithGhosts(pressure_grown, pressure_gradient_mapped, 2);
   
   CFG::LevelData<CFG::FArrayBox> pressure_gradient(cfg_grids, CFG_DIM, CFG::IntVect::Unit);
   mag_geom.unmapPoloidalGradient(pressure_gradient_mapped, pressure_gradient);

   CFG::LevelData<CFG::FArrayBox> gradPoverN( cfg_grids, CFG_DIM, CFG::IntVect::Zero );
   for (CFG::DataIterator dit(cfg_grids); dit.ok(); ++dit) {
      mag_geom.getBlockCoordSys(cfg_grids[dit]).computePsiThetaProjections(pressure_gradient[dit]);
      gradPoverN[dit].copy(pressure_gradient[dit]);
      for (int n=0; n<gradPoverN.nComp(); ++n) {
         gradPoverN[dit].divide(density[dit], 0, n);
      }
   }
   
   phase_geometry.plotConfigurationData( a_filename.c_str(), gradPoverN, a_time );

#else
   MayDay::Error("GKOps::gradPover this diagnostcs is not implemented in 3D");
#endif
}

void GKOps::plotTemperature( const std::string&    a_filename,
                             const KineticSpecies& a_soln_species,
                             const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> temperature( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.temperature( temperature );

   phase_geometry.plotConfigurationData( a_filename.c_str(), temperature, a_time );
}



void GKOps::plotParallelTemperature(const std::string&    a_filename,
                                    const KineticSpecies& a_soln_species,
                                    const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> temperature( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.parallelTemperature( temperature );

   phase_geometry.plotConfigurationData( a_filename.c_str(), temperature, a_time );
}


void GKOps::plotPerpTemperature( const std::string&    a_filename,
                                 const KineticSpecies& a_soln_species,
                                 const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> temperature( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.perpTemperature( temperature );

   phase_geometry.plotConfigurationData( a_filename.c_str(), temperature, a_time );
}


void GKOps::plotParallelHeatFlux( const std::string&    a_filename,
                                  const KineticSpecies& a_soln_species,
                                  const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   
   CFG::LevelData<CFG::FArrayBox> density( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.numberDensity( density );
   a_soln_species.ParallelMomentum( ParallelMom );
   
   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      ParallelMom[dit].divide(density[dit]);
   }
   CFG::LevelData<CFG::FArrayBox> parallelHeatFlux( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.parallelHeatFluxMoment( parallelHeatFlux, ParallelMom );
   
   phase_geometry.plotConfigurationData( a_filename.c_str(), parallelHeatFlux, a_time );
}



void GKOps::plotTotalParallelHeatFlux(const std::string&    a_filename,
                                      const KineticSpecies& a_soln_species,
                                      const double&         a_time ) const
{
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   
   CFG::LevelData<CFG::FArrayBox> zero_velocity( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   for (CFG::DataIterator dit(zero_velocity.dataIterator()); dit.ok(); ++dit) {
      zero_velocity[dit].setVal(0.0);
   }

   CFG::LevelData<CFG::FArrayBox> parallelHeatFlux( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.parallelHeatFluxMoment( parallelHeatFlux, zero_velocity );
   
   phase_geometry.plotConfigurationData( a_filename.c_str(), parallelHeatFlux, a_time );
}

void GKOps::plotFourthMoment( const std::string&    a_filename,
                              const KineticSpecies& a_soln_species,
                              const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   CFG::LevelData<CFG::FArrayBox> density( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.numberDensity( density );

   CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.ParallelMomentum( ParallelMom );
   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      ParallelMom[dit].divide(density[dit]);
   }

   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );
   a_soln_species.pressureMoment(pressure, ParallelMom);

   CFG::LevelData<CFG::FArrayBox> fourthMoment( mag_geom.gridsFull(), 1, CFG::IntVect::Zero);
   a_soln_species.fourthMoment( fourthMoment );

   CFG::LevelData<CFG::FArrayBox> temp;
   temp.define(pressure);

   for (CFG::DataIterator dit(fourthMoment.dataIterator()); dit.ok(); ++dit) {
      fourthMoment[dit].divide(temp[dit]);       // fourthMom/Pressure
      temp[dit].divide(density[dit]);            // Pressure/Density
      fourthMoment[dit].divide(temp[dit]);       // fourthMom/(N*T^2)
      fourthMoment[dit].mult(4.0/15.0);          // should be unity for Maxwellian!!!
   }

   phase_geometry.plotConfigurationData(a_filename.c_str(), fourthMoment, a_time);
}


void GKOps::plotParticleFlux( const std::string&    a_filename,
                              const KineticSpecies& a_soln_species,
                              const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   LevelData<FluxBox> E_field_tmp;

   phase_geometry.injectConfigurationToPhase( m_E_field->getFaceCenteredField(),
                                              m_E_field->getCellCenteredField(),
                                              E_field_tmp );

   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   
   CFG::LevelData<CFG::FArrayBox> particle_flux( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );
   a_soln_species.ParticleFlux( particle_flux, E_field_tmp );

   phase_geometry.plotConfigurationData( a_filename.c_str(), particle_flux, a_time );
}


void GKOps::plotHeatFlux( const std::string&    a_filename,
                          const KineticSpecies& a_soln_species,
                          const double&         a_time ) const
{
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   LevelData<FluxBox> E_field_tmp;
   LevelData<FArrayBox> phi_injected_tmp;

   phase_geometry.injectConfigurationToPhase( m_E_field->getFaceCenteredField(),
                                              m_E_field->getCellCenteredField(),
                                              E_field_tmp );

   phase_geometry.injectConfigurationToPhase( m_phi, phi_injected_tmp );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> heat_flux( mag_geom.gridsFull(), 2, CFG::IntVect::Zero );
   a_soln_species.HeatFlux( heat_flux, E_field_tmp, phi_injected_tmp );

   phase_geometry.plotConfigurationData( a_filename.c_str(), heat_flux, a_time );
}


void GKOps::plotAmpereErIncrement( const std::string&                a_filename,
				   const CFG::FluidSpeciesPtrVect&   a_soln_fluid,
				   const ScalarPtrVect&              a_soln_scalar,
				   const double&                     a_time ) const
{
  CH_assert( m_phase_geometry != NULL );
  const PhaseGeom& phase_geometry( *m_phase_geometry );
  const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   
  CFG::LevelData<CFG::FArrayBox> AmpereErIncr( mag_geom.gridsFull(), 1, CFG::IntVect::Zero );

  int scalar_component = m_Y.getScalarComponent("Er_boundary");
  const Vector<Real>& Er_boundary = a_soln_scalar[scalar_component]->data();

  int fluid_component = m_Y.getFluidComponent("Er_flux_surfaces");
  const CFG::LevelData<CFG::FArrayBox>& Er = a_soln_fluid[fluid_component]->cell_var("Er_flux_surfaces");

  const CFG::MagCoordSys& coords = *mag_geom.getCoordSys();
  const CFG::MagBlockCoordSys& block0_coord_sys = (const CFG::MagBlockCoordSys&)(*(coords.getCoordSys(0)));
  int hi_radial_index = block0_coord_sys.domain().domainBox().bigEnd(RADIAL_DIR);

  for (CFG::DataIterator dit(AmpereErIncr.dataIterator()); dit.ok(); ++dit) {
    int block_number( coords.whichBlock( mag_geom.gridsFull()[dit] ) );
    CFG::FArrayBox& this_AmpereErIncr = AmpereErIncr[dit];
        
    if ( block_number < 2 ) {
            
      const CFG::FArrayBox& this_Er = Er[dit];
      for (CFG::BoxIterator bit(this_AmpereErIncr.box()); bit.ok(); ++bit) {
         CFG::IntVect iv = bit();
         CFG::IntVect iv_shift = bit();
         iv_shift[0] = iv[0]+1;
         this_AmpereErIncr(iv,0) = this_Er(iv_shift,0);
         if (iv[0]==hi_radial_index) {
            this_AmpereErIncr(iv,0) = Er_boundary[1];
         }
      }
    }
    else {
      this_AmpereErIncr.setVal(0.);
    }
  }

  phase_geometry.plotConfigurationData( a_filename.c_str(), AmpereErIncr, a_time );
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
                               const bool a_startup_flag )
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
      CFG::LevelData<CFG::FArrayBox>* field = NULL;

      if (hist_field_name == "potential") {
         field = &m_phi;
         fname = "potential_hist";
         stringstream fname_temp;
         fname_temp << "potential_hist_" << ihist + 1 << ".curve" << ends;
         fname_suf = fname_temp.str();
         stringstream fname_step_temp;
         fname_step_temp << "potential_hist_" << ihist + 1 << ".txt" << ends;
         fname_step = fname_step_temp.str();
      }
      else if (hist_field_name == "Efield") {

         CFG::LevelData<CFG::FArrayBox>& E_field_cell = m_E_field->getCellCenteredField();
         field = &E_field_cell;
         fname = "Efield_hist";
         stringstream fname_temp;
         fname_temp << "Efield_hist_" << ihist + 1 << ".curve" << ends;
         fname_suf = fname_temp.str();
         stringstream fname_step_temp;
         fname_step_temp << "Efield_hist_" << ihist + 1 << ".txt" << ends;
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
   const CFG::EllipticOpBC&              a_bc, 
   const CFG::LevelData<CFG::FArrayBox>& a_ion_density) const
{
   m_E_field->computeQuasiNeutralElectronDensity( a_quasineutral_density,
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
      if (fluid_species[s]->name() == "potential") {
         includes_potential = true;
         break;
      }
   }

   return includes_potential;
}

#include "NamespaceFooter.H"

