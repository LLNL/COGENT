#include "GKOps.H"
#include "CONSTANTS.H"
#include "FORT_PROTO.H"
#include "inspect.H"

#include <fstream>
#include <sstream>

#include "MillerPhaseCoordSys.H"
#include "SingleNullPhaseCoordSys.H"
#include "SNCorePhaseCoordSys.H"
#include "SlabPhaseCoordSys.H"
#include "Kernels.H"


#include "GKCollisions.H"

#include "LocalizedF_F.H"

#if 1  // warning, OS dependencies, will not work on all platforms
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "CH_HDF5.H"
#include "HDF5Portable.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "MagBlockCoordSys.H"
#include "MillerBlockCoordSys.H"
#include "MillerCoordSys.H"
#include "SlabBlockCoordSys.H"
#include "SlabCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "SNCoreBlockCoordSys.H"
#include "newMappedGridIO.H"
#include "FourthOrderUtil.H"
#include "DataArray.H"
#include "inspect.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM 3
#include "AMRIO.H"
#include "Slicing.H.transdim"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM


#include "FourthOrderUtil.H"
#include "Directions.H"

#include "NamespaceHeader.H"

#define EXTRAPOLATE_FIELD

using namespace CH_MultiDim;

#if 0
#ifdef USE_ARK4
void GKOps::define( const GKState& a_state,
                    const Real a_dt,
                    const Real a_dt_scale )
{
   m_dt_scale = a_dt_scale;
   define( a_state, a_dt );
}
#endif
#endif


void GKOps::define( const GKState& a_state,
                    const Real a_dt )
{
   CH_assert( isDefined()==false );
   CH_assert( a_dt>0.0 );
   m_dt = a_dt;
   m_phase_geometry = a_state.geometry();
   CH_assert( m_phase_geometry != NULL );

   ParmParse ppgksys( "gksystem" );
   parseParameters( ppgksys );

   m_units = new GKUnits( ppgksys );

   ParmParse pp;
   m_initial_conditions = new GKSystemIC( pp, *m_phase_geometry, a_state.data() );
   m_boundary_conditions = new GKSystemBC( pp, *m_phase_geometry, a_state.data() );
   
   const double larmor( m_units->larmorNumber() );
   ParmParse pp_vlasov(GKVlasov::pp_name);
   m_vlasov      = new GKVlasov( pp_vlasov, larmor );
   m_collisions  = new GKCollisions( m_verbosity );
   if (m_transport_model_on) {
      m_transport = new GKTransport( m_verbosity );
   }
   if (m_neutrals_model_on) {
      m_neutrals = new GKNeutrals( m_verbosity );
   }
 

   CFG::IntVect phi_ghost_vect( 4*CFG::IntVect::Unit );
   m_phi.define( m_phase_geometry->magGeom().grids(), 1, phi_ghost_vect );

   m_is_defined = true;
}


void
GKOps::initialize( const KineticSpeciesPtrVect& a_kinetic_species,
                   const int                    a_cur_step )
{
   KineticSpeciesPtrVect initial_kinetic_species;

   // Would use createTemporarySpeciesVector() here, but it clones with ghost cells
   // and divides by J....
   initial_kinetic_species.resize( a_kinetic_species.size() );
   for (int s(0); s<a_kinetic_species.size(); s++) {
      initial_kinetic_species[s] = a_kinetic_species[s]->clone( IntVect::Zero );
   }

   if ( a_cur_step > 0 ) {
      // If we are restarting, initial_kinetic_species contains the current
      // solution, so overwrite it with the initial condition.
      applyInitialConditions( initial_kinetic_species, 0. );

      // And get rid of J
      for (int s(0); s<initial_kinetic_species.size(); s++) {
         LevelData<FArrayBox>& dfn( initial_kinetic_species[s]->distributionFunction() );
         const PhaseGeom& geometry( initial_kinetic_species[s]->phaseSpaceGeometry() );
         geometry.divideJonValid( dfn );
      }
   }
   else if (m_consistent_potential_bcs) {
      m_Er_lo = 0.;
      m_Er_hi = 0.;
   }

   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   const CFG::DisjointBoxLayout& mag_grids = mag_geom.grids();

   if (m_consistent_potential_bcs) {

      m_radial_gkp_divergence_average.define( mag_grids, 1, CFG::IntVect::Zero );
      m_radial_flux_divergence_average.define( mag_grids, 1, CFG::IntVect::Zero );
    
      if (a_cur_step==0) {
         m_Er_average_face.define( mag_grids, 3, 2*CFG::IntVect::Unit );
         m_Er_average_cell.define( mag_grids, 3, 2*CFG::IntVect::Unit );
         m_E_tilde_face.define( mag_grids, 3, 2*CFG::IntVect::Unit );
         m_E_tilde_cell.define( mag_grids, 3, 2*CFG::IntVect::Unit );

         for (CFG::DataIterator dit(mag_grids); dit.ok(); ++dit) {
            m_Er_average_face[dit].setVal(0.);
            m_Er_average_cell[dit].setVal(0.);
            m_E_tilde_face[dit].setVal(0.);
            m_E_tilde_cell[dit].setVal(0.);
         }
      }
   }

   CFG::LevelData<CFG::FArrayBox> initial_ion_charge_density( mag_grids, 1, CFG::IntVect::Zero);
   computeIonChargeDensity( initial_ion_charge_density, initial_kinetic_species );

   createGKPoisson( initial_ion_charge_density );
    
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   m_E_field_face.define( mag_grids, 3, 2*CFG::IntVect::Unit );
   m_E_field_cell.define( mag_grids, 3, 2*CFG::IntVect::Unit );
   computeEField( m_E_field_face,
                  m_E_field_cell,
                  a_kinetic_species,
                  a_cur_step );
    
   if ( m_ampere_law && a_cur_step == 0 ) {
      //Initialize E-field
      CFG::PotentialBC& bc = m_boundary_conditions->getPotentialBC();
      m_poisson->fillGhosts( m_phi );
#ifdef EXTRAPOLATE_FIELD
      m_poisson->computeElectricField( m_phi, m_Er_average_cell );
      m_poisson->computeElectricField( m_phi, m_Er_average_face );
#else
      m_poisson->computeElectricField( m_phi, bc, m_Er_average_cell );
      m_poisson->computeElectricField( m_phi, bc, m_Er_average_face );
#endif
   }
}

GKOps::~GKOps()
{
   if (m_boltzmann_electron) delete m_boltzmann_electron;
   delete m_collisions;
   if (m_transport) delete m_transport;
   if (m_neutrals)  delete m_neutrals;
   delete m_vlasov;
   if (m_poisson) delete m_poisson;
   for (int ihist(0); ihist<m_hist_count; ihist++) {
     FieldHist *save_hist = &m_fieldHistLists[ihist];
     delete save_hist->fieldvals;
     delete save_hist->timevals;
     delete save_hist->timestep;
   }
}

Real GKOps::dtScale_Vlasov( const KineticSpeciesPtrVect& a_soln,
                             const int a_step_number )
{
   return m_TimeScale_Vlasov;
}

Real GKOps::dtScale_Collisions( const KineticSpeciesPtrVect& a_soln,
                                 const int a_step_number )
{
  return m_TimeScale_Collisions;
}

Real GKOps::stableDt( const KineticSpeciesPtrVect& a_soln,
                      const int a_step_number )
{
   CH_assert( isDefined() );
   LevelData<FluxBox> E_field;
   computeElectricField( E_field, a_soln, a_step_number );

   m_dt_Vlasov      = m_vlasov->computeDt( E_field, a_soln );
   m_dt_Collisions  = m_collisions->computeDt( a_soln );
   Real dt_stable   = Min( m_dt_Vlasov, m_dt_Collisions );

   /* compute and store the time scales */
   m_TimeScale_Vlasov     = m_vlasov->computeTimeScale(E_field, a_soln);
   m_TimeScale_Collisions = m_collisions->computeTimeScale( a_soln );

   if (m_transport_model_on) {
      Real dt_transport( m_transport->computeDt( a_soln ) );
      dt_stable = Min( dt_stable, dt_transport );
   }

   if (m_neutrals_model_on) {
      Real dt_neutrals( m_neutrals->computeDt( a_soln ) );
      dt_stable = Min( dt_stable, dt_neutrals );
   }
 
   return dt_stable;
}


void GKOps::explicitOp( GKRHSData& a_rhs,
                        const Real a_time,
                        const GKState& a_state,
                        const int a_stage )
{
   CH_assert( isDefined() );

   a_rhs.zero();

   const KineticSpeciesPtrVect& species_comp( a_state.data() );
   
   LevelData<FluxBox> E_field;
   if (m_consistent_potential_bcs) {
      if (a_stage == 0) {
         m_stage0_time = a_time;
      }

      // We're not fourth-order accurate with this option anyway,
      // so only compute the field at the beginning of the step.
      computeElectricField( E_field, species_comp, -1, -1 );
   }
   else {
      computeElectricField( E_field, species_comp, -1, a_stage );
   }

   KineticSpeciesPtrVect species_phys;
   createTemporarySpeciesVector( species_phys, species_comp );
   fillGhostCells( species_phys, E_field, a_time );
   applyVlasovOperator( a_rhs.data(), species_phys, E_field, a_time );

   if (m_transport_model_on) {
      applyTransportOperator( a_rhs.data(), species_phys, a_time );
   }
    
   if (m_neutrals_model_on) {
      applyNeutralsOperator( a_rhs.data(), species_comp, a_time );
   }

   applyCollisionOperator( a_rhs.data(), species_comp, a_time );

   if (m_consistent_potential_bcs && a_stage == 3) {

      double dt = a_time - m_stage0_time;

      m_Er_lo += dt * (-m_lo_radial_flux_divergence_average / m_lo_radial_gkp_divergence_average);
      m_Er_hi += dt * (-m_hi_radial_flux_divergence_average / m_hi_radial_gkp_divergence_average);

      setCoreBC( m_Er_lo, -m_Er_hi, m_boundary_conditions->getPotentialBC() );

      if (procID()==0) cout << "  Er_lo = " << m_Er_lo << ", Er_hi = " << m_Er_hi << endl;
       
      if (m_ampere_law) {

         CFG::LevelData<CFG::FluxBox> E_tilde_mapped_face;
         CFG::LevelData<CFG::FArrayBox> E_tilde_mapped_cell;
         CFG::LevelData<CFG::FArrayBox> phi_tilde_fs_average;
         double phi_tilde_fs_average_lo;
         double phi_tilde_fs_average_hi;

         if (!m_ampere_cold_electrons) {

            const CFG::MagGeom& mag_geom = m_phase_geometry->magGeom();
            const CFG::DisjointBoxLayout& mag_grids = mag_geom.grids();

            // Compute the poloidal variation phi_tilde.  The poloidal variation of the
            // current potential is used as the initial guess for the interative procedure.
            int num_phi_ghosts_filled = m_poisson->numPotentialGhostsFilled();
            CH_assert(num_phi_ghosts_filled >=2 );
            CFG::LevelData<CFG::FArrayBox> phi_tilde(mag_grids, 1, num_phi_ghosts_filled*CFG::IntVect::Unit);

            computePhiTilde( species_phys, *m_boltzmann_electron, phi_tilde );

            // Fill the two ghost cell layers at block boundaries as required by the
            // subsequent field calculations
            m_poisson->fillGhosts( phi_tilde );

            E_tilde_mapped_face.define(mag_grids, 3, CFG::IntVect::Zero);
            m_poisson->computeMappedElectricField( phi_tilde, E_tilde_mapped_face );

            E_tilde_mapped_cell.define(mag_grids, 3, CFG::IntVect::Zero);
            m_poisson->computeMappedElectricField( phi_tilde, E_tilde_mapped_cell );

            CFG::LevelData<CFG::FluxBox> Er_tilde_mapped(mag_grids, 1, CFG::IntVect::Zero);
            m_poisson->extractNormalComponent(E_tilde_mapped_face, RADIAL_DIR, Er_tilde_mapped);

            phi_tilde_fs_average.define(mag_grids, 1, CFG::IntVect::Zero);
            m_poisson->computeRadialFSAverage( Er_tilde_mapped, phi_tilde_fs_average_lo, phi_tilde_fs_average_hi, phi_tilde_fs_average );
         }

         updateAveragedEfield( m_Er_average_face, m_Er_average_cell,
                               m_radial_flux_divergence_average, m_lo_radial_flux_divergence_average, m_hi_radial_flux_divergence_average, dt );
          
         updateEfieldPoloidalVariation( E_tilde_mapped_face, E_tilde_mapped_cell, m_E_tilde_face, m_E_tilde_cell,
                                        phi_tilde_fs_average, phi_tilde_fs_average_lo, phi_tilde_fs_average_hi );
      }
   }
}

void GKOps::explicitOpImEx( GKRHSData& a_rhs,
                            const Real a_time,
                            const GKState& a_state,
                            const int a_stage )
{
   CH_assert( isDefined() );

   a_rhs.zero();

   const KineticSpeciesPtrVect& species_comp( a_state.data() );
   
   LevelData<FluxBox> E_field;
   if (m_consistent_potential_bcs) {
      if (a_stage == 0) {
         m_stage0_time = a_time;
      }

      // We're not fourth-order accurate with this option anyway,
      // so only compute the field at the beginning of the step.
      computeElectricField( E_field, species_comp, -1, -1 );
   }
   else {
      computeElectricField( E_field, species_comp, -1, a_stage );
   }

   KineticSpeciesPtrVect species_phys;
   createTemporarySpeciesVector( species_phys, species_comp );
   fillGhostCells( species_phys, E_field, a_time );
   applyVlasovOperator( a_rhs.data(), species_phys, E_field, a_time );

   if (m_transport_model_on) {
      applyTransportOperator( a_rhs.data(), species_phys, a_time );
   }
    
   if (m_neutrals_model_on) {
      applyNeutralsOperator( a_rhs.data(), species_comp, a_time );
   }

   /* The following is hard-coded for a 4-stage time integrator!! */
   if (m_consistent_potential_bcs && a_stage == 3) {

      double dt = a_time - m_stage0_time;

      m_Er_lo += dt * (-m_lo_radial_flux_divergence_average / m_lo_radial_gkp_divergence_average);
      m_Er_hi += dt * (-m_hi_radial_flux_divergence_average / m_hi_radial_gkp_divergence_average);

      setCoreBC( m_Er_lo, -m_Er_hi, m_boundary_conditions->getPotentialBC() );

      if (procID()==0) cout << "  Er_lo = " << m_Er_lo << ", Er_hi = " << m_Er_hi << endl;
       
      if (m_ampere_law) {

         CFG::LevelData<CFG::FluxBox> E_tilde_mapped_face;
         CFG::LevelData<CFG::FArrayBox> E_tilde_mapped_cell;
         CFG::LevelData<CFG::FArrayBox> phi_tilde_fs_average;
         double phi_tilde_fs_average_lo;
         double phi_tilde_fs_average_hi;

         if (!m_ampere_cold_electrons) {

            const CFG::MagGeom& mag_geom = m_phase_geometry->magGeom();
            const CFG::DisjointBoxLayout& mag_grids = mag_geom.grids();

            // Compute the poloidal variation phi_tilde.  The poloidal variation of the
            // current potential is used as the initial guess for the interative procedure.
            int num_phi_ghosts_filled = m_poisson->numPotentialGhostsFilled();
            CH_assert(num_phi_ghosts_filled >=2 );
            CFG::LevelData<CFG::FArrayBox> phi_tilde(mag_grids, 1, num_phi_ghosts_filled*CFG::IntVect::Unit);

            computePhiTilde( species_phys, *m_boltzmann_electron, phi_tilde );

            // Fill the two ghost cell layers at block boundaries as required by the
            // subsequent field calculations
            m_poisson->fillGhosts( phi_tilde );

            E_tilde_mapped_face.define(mag_grids, 3, CFG::IntVect::Zero);
            m_poisson->computeMappedElectricField( phi_tilde, E_tilde_mapped_face );

            E_tilde_mapped_cell.define(mag_grids, 3, CFG::IntVect::Zero);
            m_poisson->computeMappedElectricField( phi_tilde, E_tilde_mapped_cell );

            CFG::LevelData<CFG::FluxBox> Er_tilde_mapped(mag_grids, 1, CFG::IntVect::Zero);
            m_poisson->extractNormalComponent(E_tilde_mapped_face, RADIAL_DIR, Er_tilde_mapped);

            phi_tilde_fs_average.define(mag_grids, 1, CFG::IntVect::Zero);
            m_poisson->computeRadialFSAverage( Er_tilde_mapped, phi_tilde_fs_average_lo, phi_tilde_fs_average_hi, phi_tilde_fs_average );
         }

         updateAveragedEfield( m_Er_average_face, m_Er_average_cell,
                               m_radial_flux_divergence_average, m_lo_radial_flux_divergence_average, m_hi_radial_flux_divergence_average, dt );
          
         updateEfieldPoloidalVariation( E_tilde_mapped_face, E_tilde_mapped_cell, m_E_tilde_face, m_E_tilde_cell,
                                        phi_tilde_fs_average, phi_tilde_fs_average_lo, phi_tilde_fs_average_hi );
      }
   }
}

void GKOps::implicitOpImEx( GKRHSData& a_rhs,
                            const Real a_time,
                            const GKState& a_state,
                            const int a_stage )
{
   CH_assert( isDefined() );

   const KineticSpeciesPtrVect& species_comp( a_state.data() );
   a_rhs.zero();
   applyCollisionOperator( a_rhs.data(), species_comp, a_time );
}

#if 0
void GKOps::solve( GKState& a_state,
                   const GKRHSData& a_rhs,
                   const Real a_time,
                   const int a_stage )
// Solve for a_state in
// (I - s_aIdiag * m_dt * FI) (a_state) = a_rhs.
{

// So we could think of this iteratively as
//
//    [G(phi_s+eps)-G(phi_s)]/eps * dphi_s = G(phi_s),
//
// with phi_{s+1} = phi_s + dphi_s, where
//
//    G(a_state) = a_state - s_aIdiag * dt * f_I(a_state) - a_rhs = 0.
//
// The only thing that is necessary from an operator treated implicitly is then
// the implicit rhs function f_I(phi).

   CH_assert( isDefined() );

   KineticSpeciesPtrVect& species_comp( a_state.data() );
   const KineticSpeciesPtrVect& species_rhs( a_rhs.data() );

   const int num_kinetic_species( species_comp.size() );
   for (int s(0); s<num_kinetic_species; s++) {
      
      KineticSpecies& soln_species( *(species_comp[s]) );
      CLSInterface& cls_ref( m_collisions->collisionModel( soln_species.name() ) );
      MyKrook* krook = NULL;
      krook = dynamic_cast<MyKrook*>( &cls_ref );

      if (krook) {
         Real nu( krook->collisionFrequency() );
         CH_assert( nu>=0.0 );
         Real coeff( m_dt_scale * m_dt * nu );
         Real factor( 1.0 / (1.0 + coeff) );
         
         LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );
         
         const KineticSpecies& rhs_species( *(species_rhs[s]) );
         const LevelData<FArrayBox>& rhs_dfn( rhs_species.distributionFunction() );
         
         rhs_dfn.copyTo( soln_dfn );
         krook->addReferenceDfn( soln_species, a_time, coeff );
         for (DataIterator sdit(soln_dfn.dataIterator()); sdit.ok(); ++sdit) {
            soln_dfn[sdit].mult( factor );
         }
         
      }
      else {
         MayDay::Error( "Implicit solution is currently limited to linear, non-conservative Krook operator!" );
      }
   }
}
#endif

/* Solve for a_state in  (I - shift * FI) (a_state) = a_rhs. 
   So we could think of this iteratively as

      [G(phi_s+eps)-G(phi_s)]/eps * dphi_s = G(phi_s),

   with phi_{s+1} = phi_s + dphi_s, where

      G(a_state) = a_state - shift * f_I(a_state) - a_rhs = 0.

   The only thing that is necessary from an operator treated implicitly is then
   the implicit rhs function f_I(phi).
 */
void GKOps::solveImEx( GKState& a_state,
                       const GKRHSData& a_rhs,
                       const Real shift, 
                       const Real a_time,
                       const int a_stage )
{
   CH_assert( isDefined() );

   KineticSpeciesPtrVect& species_comp( a_state.data() );
   const KineticSpeciesPtrVect& species_rhs( a_rhs.data() );

   const int num_kinetic_species( species_comp.size() );
   for (int s(0); s<num_kinetic_species; s++) {
      
      KineticSpecies& soln_species( *(species_comp[s]) );
      std::string CLSModelName (m_collisions->collisionModelName( soln_species.name() ));
      CLSInterface& cls_ref( m_collisions->collisionModel( soln_species.name() ) );

      if (CLSModelName == _CLS_KROOK_) {
        Krook* clsObject = dynamic_cast<Krook*>( &cls_ref );
        Real nu( clsObject->collisionFrequency() );
        CH_assert( nu>=0.0 );
        Real coeff( shift * nu );
        Real factor( 1.0 / (1.0 + coeff) );
         
        LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );
         
        const KineticSpecies& rhs_species( *(species_rhs[s]) );
        const LevelData<FArrayBox>& rhs_dfn( rhs_species.distributionFunction() );
         
        rhs_dfn.copyTo( soln_dfn );
        clsObject->addReferenceDfn( soln_species, a_time, coeff );
        for (DataIterator sdit(soln_dfn.dataIterator()); sdit.ok(); ++sdit) {
           soln_dfn[sdit].mult( factor );
        }
      } else if (CLSModelName == _CLS_MYKROOK_) {
        MyKrook* clsObject = dynamic_cast<MyKrook*>( &cls_ref );
        Real nu( clsObject->collisionFrequency() );
        CH_assert( nu>=0.0 );
        Real coeff( shift * nu );
        Real factor( 1.0 / (1.0 + coeff) );
         
        LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );
         
        const KineticSpecies& rhs_species( *(species_rhs[s]) );
        const LevelData<FArrayBox>& rhs_dfn( rhs_species.distributionFunction() );
         
        rhs_dfn.copyTo( soln_dfn );
        clsObject->addReferenceDfn( soln_species, a_time, coeff );
        for (DataIterator sdit(soln_dfn.dataIterator()); sdit.ok(); ++sdit) {
           soln_dfn[sdit].mult( factor );
        }
      } else if (CLSModelName == _CLS_NONE_) {
        a_state.zero();
        a_state.increment(a_rhs);
      } else {
        MayDay::Error( "Implicit solution solver is not yet implemented for the chosen collision model!" );
      }
   }
}


inline
void GKOps::computeElectricField( LevelData<FluxBox>&          a_E_field,
                                  const KineticSpeciesPtrVect& a_kinetic_species,
                                  const int                    a_step_number,
                                  const int                    a_stage )
{
   // If this is the first evaluation this step, then a field solve
   // has already been performed to compute the potential needed to
   // evaluate the velocity for determining the stable time step.
   // We therefore check for this so we don't do another (redundant)
   // field solve.
   CH_assert( isDefined() );
   if ( a_stage > 0 ) {
      const int num_kinetic_species( a_kinetic_species.size() );
      KineticSpeciesPtrVect result;
      result.resize(num_kinetic_species);
      for (int species(0); species<num_kinetic_species; species++) {
         result[species] = a_kinetic_species[species]->clone( m_ghost_vect );
      }
      divideJ( a_kinetic_species, result ); /// FIXME
      computeEField( m_E_field_face,
                     m_E_field_cell,
                     result,
                     a_step_number );
   }

   CH_assert( m_phase_geometry != NULL );
   m_phase_geometry->injectConfigurationToPhase( m_E_field_face,
                                                 m_E_field_cell,
                                                 a_E_field );
}


inline
void GKOps::createTemporarySpeciesVector( KineticSpeciesPtrVect& a_out,
                                          const KineticSpeciesPtrVect& a_in )
{
   a_out.resize( a_in.size() );
   for (int s(0); s<a_in.size(); s++) {
      a_out[s] = a_in[s]->clone( m_ghost_vect );
      LevelData<FArrayBox>& out_dfn( a_out[s]->distributionFunction() );
      const PhaseGeom& geometry( a_in[s]->phaseSpaceGeometry() );
      geometry.divideJonValid( out_dfn );
   }
}


inline
void GKOps::fillGhostCells(
   KineticSpeciesPtrVect&       a_species_phys,
   const LevelData<FluxBox>&    a_E_field,
   const Real&                  a_time )
{
   CH_assert( isDefined() );
   m_boundary_conditions->fillGhostCells( a_species_phys,
                                          m_phi,
                                          a_E_field,
                                          a_time );
}


inline
void GKOps::applyVlasovOperator(
   KineticSpeciesPtrVect&       a_rhs,
   const KineticSpeciesPtrVect& a_soln,
   const LevelData<FluxBox>&    a_E_field,
   const Real&                  a_time )
{
   CH_assert( isDefined() );

   if (m_consistent_potential_bcs) {
      m_lo_radial_flux_divergence_average = 0.;
      m_hi_radial_flux_divergence_average = 0.;

      for (int s(0); s<a_rhs.size(); s++) {
         const KineticSpecies& soln_species( *(a_soln[s]) );
         KineticSpecies& rhs_species( *(a_rhs[s]) );
         double lo_value, hi_value;
         m_vlasov->evalRHS( rhs_species, lo_value, hi_value, m_radial_flux_divergence_average, soln_species, a_E_field, a_time );

         m_lo_radial_flux_divergence_average += lo_value;
         m_hi_radial_flux_divergence_average += hi_value;
      }
   }
   else {
      for (int s(0); s<a_rhs.size(); s++) {
         const KineticSpecies& soln_species( *(a_soln[s]) );
         KineticSpecies& rhs_species( *(a_rhs[s]) );
         m_vlasov->evalRHS( rhs_species, soln_species, a_E_field, a_time );
      }
   }
}


inline
void GKOps::applyCollisionOperator( KineticSpeciesPtrVect&       a_rhs,
                                    const KineticSpeciesPtrVect& a_soln,
                                    const Real&                  a_time )
{
   CH_assert( isDefined() );
   m_collisions->accumulateRHS( a_rhs, a_soln, a_time );
}


inline
void GKOps::applyTransportOperator( KineticSpeciesPtrVect&       a_rhs,
                                    const KineticSpeciesPtrVect& a_soln,
                                    const Real&                  a_time )
{
   CH_assert( isDefined() );
   m_transport->accumulateRHS( a_rhs, a_soln, a_time );

}

inline
void GKOps::applyNeutralsOperator( KineticSpeciesPtrVect&       a_rhs,
                                   const KineticSpeciesPtrVect& a_soln,
                                   const Real&                  a_time )
{
    CH_assert( isDefined() );
    m_neutrals->accumulateRHS( a_rhs, a_soln, a_time );
    
}


inline
void setZero( CFG::LevelData<CFG::FArrayBox>& a_data )
{
   for (CFG::DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit].setVal(0.);
   }
}


void
GKOps::computeTotalChargeDensity( CFG::LevelData<CFG::FArrayBox>& a_charge_density,
                                  const KineticSpeciesPtrVect&    a_species ) const
{
   // Container for individual species charge density
   CH_assert( m_phase_geometry != NULL );
   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   CFG::LevelData<CFG::FArrayBox> species_charge_density( mag_geom.grids(),
                                                          1,
                                                          CFG::IntVect::Zero );

   setZero( a_charge_density );

   for (int species(0); species<a_species.size(); species++) {
      const KineticSpecies& this_species( *(a_species[species]) );

      // Compute the charge density for this species
      this_species.chargeDensity( species_charge_density );
      
      CFG::DataIterator dit( a_charge_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         a_charge_density[dit].plus( species_charge_density[dit] );
      }
   }
}


void
GKOps::computeIonChargeDensity( CFG::LevelData<CFG::FArrayBox>& a_ion_charge_density,
                                const KineticSpeciesPtrVect&    a_species ) const
{
   // Container for individual species charge density
   CH_assert( m_phase_geometry != NULL );
   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   CFG::LevelData<CFG::FArrayBox> species_charge_density( mag_geom.grids(),
                                                          1,
                                                          CFG::IntVect::Zero );

   setZero( a_ion_charge_density );

   for (int species(0); species<a_species.size(); species++) {

      const KineticSpecies& this_species( *(a_species[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      // Compute the charge density for this species
      this_species.chargeDensity( species_charge_density );
      
      CFG::DataIterator dit( a_ion_charge_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         a_ion_charge_density[dit].plus( species_charge_density[dit] );
      }
   }
}

void
GKOps::computeIonParallelCurrentDensity( CFG::LevelData<CFG::FArrayBox>& a_ion_current_density,
                                         const KineticSpeciesPtrVect&    a_species ) const
{
   // Container for individual species charge density
   CH_assert( m_phase_geometry != NULL );
   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   CFG::LevelData<CFG::FArrayBox> species_current_density( mag_geom.grids(),
                                                          1,
                                                          CFG::IntVect::Zero );

   setZero( a_ion_current_density );

   for (int species(0); species<a_species.size(); species++) {

      const KineticSpecies& this_species( *(a_species[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      // Compute the charge density for this species
      this_species.ParallelMomentum( species_current_density );
      
      CFG::DataIterator dit( a_ion_current_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         species_current_density[dit].mult( this_species.charge() );
         a_ion_current_density[dit].plus( species_current_density[dit] );
      }
   }
}


void
GKOps::computeSignedChargeDensities( CFG::LevelData<CFG::FArrayBox>& a_pos_charge_density,
                                     CFG::LevelData<CFG::FArrayBox>& a_neg_charge_density,
                                     const KineticSpeciesPtrVect& a_species ) const
{
   // Container for individual species charge density
   CH_assert( m_phase_geometry != NULL );
   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   CFG::LevelData<CFG::FArrayBox> species_charge_density( mag_geom.grids(),
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


void
GKOps::computeIonMassDensity( CFG::LevelData<CFG::FArrayBox>& a_mass_density,
                              const KineticSpeciesPtrVect&    a_species ) const
{
   // Container for individual species charge density
   CH_assert( m_phase_geometry != NULL );
   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   CFG::LevelData<CFG::FArrayBox> species_mass_density( mag_geom.grids(),
                                                        1,
                                                        CFG::IntVect::Zero );

   setZero( a_mass_density );

   for (int species(0); species<a_species.size(); species++) {
         
      const KineticSpecies& this_species( *(a_species[species]) );
      if ( this_species.charge() < 0.0 ) continue;

      // Compute the charge density for this species
      this_species.massDensity( species_mass_density );
      
      CFG::DataIterator dit( a_mass_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         a_mass_density[dit].plus( species_mass_density[dit] );
      }
   }
}


void GKOps::computeEField( CFG::LevelData<CFG::FluxBox>&   a_E_field_face,
                           CFG::LevelData<CFG::FArrayBox>& a_E_field_cell,
                           const KineticSpeciesPtrVect&    a_soln,
                           const int                       a_step_number )
{
   CH_assert( isDefined() );

   if (m_poisson) {
      CFG::PotentialBC& bc = m_boundary_conditions->getPotentialBC();

      const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
      const CFG::DisjointBoxLayout& grids( mag_geom.grids() );

      if (!m_fixed_efield){
         CH_assert( m_phase_geometry != NULL );
         // Update the potential and field, if not fixed_efield

         CFG::LevelData<CFG::FArrayBox> ion_mass_density( grids, 1, CFG::IntVect::Zero );
         computeIonMassDensity( ion_mass_density, a_soln );

         if (m_boltzmann_electron == NULL) {
            m_poisson->setOperatorCoefficients( ion_mass_density, bc );

            CFG::LevelData<CFG::FArrayBox> total_charge_density( grids, 1, CFG::IntVect::Zero );
            computeTotalChargeDensity( total_charge_density, a_soln );
            m_poisson->computePotential( m_phi, bc, total_charge_density );

         }
         else {

            // Boltzmann electron model

            if (m_consistent_potential_bcs) {
               setCoreBC( m_Er_lo, -m_Er_hi, bc );
            }

            CFG::LevelData<CFG::FArrayBox> ion_charge_density( grids, 1, CFG::IntVect::Zero );
            computeIonChargeDensity( ion_charge_density, a_soln );
            
            bool single_null = typeid(*(mag_geom.getCoordSys())) == typeid(CFG::SingleNullCoordSys);

            if (single_null) {

               CFG::LevelData<CFG::FArrayBox> ion_parallel_current_density( grids, 1, CFG::IntVect::Zero );
               computeIonParallelCurrentDensity( ion_parallel_current_density, a_soln );

              ((CFG::NewGKPoissonBoltzmann*)m_poisson)
                 ->setDivertorBVs( ion_charge_density, ion_parallel_current_density, bc );
            }

            if (m_consistent_potential_bcs) {

               double lo_value, hi_value;

               if (single_null) {

                  double core_outer_bv = -m_Er_hi;

                  ((CFG::NewGKPoissonBoltzmann*)m_poisson)
                     ->setOperatorCoefficients( ion_mass_density, bc, core_outer_bv, lo_value, hi_value, m_radial_gkp_divergence_average );
               }
               else {
                  ((CFG::GKPoissonBoltzmann*)m_poisson)
                     ->setOperatorCoefficients( ion_mass_density, bc, lo_value, hi_value, m_radial_gkp_divergence_average );
               }

               m_lo_radial_gkp_divergence_average = lo_value;
               m_hi_radial_gkp_divergence_average = hi_value;
            }
            else {
               m_poisson->setOperatorCoefficients( ion_mass_density, bc );
            }

            if ( !m_ampere_law ) {
             
               if (single_null) {
                  ((CFG::NewGKPoissonBoltzmann*)m_poisson)
                     ->computePotentialAndElectronDensity( m_phi,
                                                           *m_boltzmann_electron,
                                                           ion_charge_density,
                                                           bc,
                                                           a_step_number==0 );
               }
               else {
                  ((CFG::GKPoissonBoltzmann*)m_poisson)
                     ->computePotentialAndElectronDensity( m_phi,
                                                           *m_boltzmann_electron,
                                                           ion_charge_density,
                                                           bc,
                                                           a_step_number==0 );
               }
            }
         }
      }

      if (!m_fixed_efield || (m_initializedE==false)){
         // only calculate E if not fixed efield or if this is the first step in the simulation.
         //  Need to take care that this is done right for a restart
         m_poisson->fillGhosts( m_phi );

         if (m_ampere_law) {
            for (CFG::DataIterator dit(a_E_field_cell.dataIterator()); dit.ok(); ++dit) {
               a_E_field_face[dit].copy(m_Er_average_face[dit]);
               a_E_field_face[dit] += m_E_tilde_face[dit];
               a_E_field_cell[dit].copy(m_Er_average_cell[dit]);
               a_E_field_cell[dit] += m_E_tilde_cell[dit];
            }
         }
         else {

#ifdef EXTRAPOLATE_FIELD
            m_poisson->computeElectricField( m_phi, a_E_field_cell );
            m_poisson->computeElectricField( m_phi, a_E_field_face );
#else
            m_poisson->computeElectricField( m_phi, bc, a_E_field_cell );
            m_poisson->computeElectricField( m_phi, bc, a_E_field_face );
#endif
         }          
         m_initializedE = true;
      }
   }
   else if (m_initializedE==false) {
      // Only calculate E if not fixed efield or if this is the first step in the simulation.
      // Need to take care that this is done right for a restart
         
      for (CFG::DataIterator dit( a_E_field_face.dataIterator() ); dit.ok(); ++dit) {
         a_E_field_cell[dit].setVal(0.);
         for (int dir(0); dir<CFG_DIM; ++dir) {
            a_E_field_face[dit].setVal(0.,dir);
         }
      }
      m_initializedE = true;
   }
}


void
GKOps::setCoreBC( const double      a_core_inner_bv,
                  const double      a_core_outer_bv,
                  CFG::PotentialBC& a_bc ) const 
{
   const CFG::MagGeom& mag_geom = m_phase_geometry->magGeom();

   if ( typeid(*(mag_geom.getCoordSys())) == typeid(CFG::MillerCoordSys) ) {
      a_bc.setBCValue(0,RADIAL_DIR,0,a_core_inner_bv);
      a_bc.setBCValue(0,RADIAL_DIR,1,a_core_outer_bv);
   }
   else if ( typeid(*(mag_geom.getCoordSys())) == typeid(CFG::SlabCoordSys) ) {
      a_bc.setBCValue(0,RADIAL_DIR,0,a_core_inner_bv);
      a_bc.setBCValue(0,RADIAL_DIR,1,a_core_outer_bv);
   }
   else if ( typeid(*(mag_geom.getCoordSys())) == typeid(CFG::SNCoreCoordSys) ) {
      a_bc.setBCValue(L_CORE,RADIAL_DIR,0,a_core_inner_bv);
      a_bc.setBCValue(L_CORE,RADIAL_DIR,1,a_core_outer_bv);
   }
   else if ( typeid(*(mag_geom.getCoordSys())) == typeid(CFG::SingleNullCoordSys) ) {
      a_bc.setBCValue(LCORE,RADIAL_DIR,0,a_core_inner_bv);
   }
   else {
      MayDay::Error("GKOps::computeEfield(): unknown geometry with consistent bcs");
   }
}


void
GKOps::computePhiTilde( const KineticSpeciesPtrVect&          a_kinetic_species,
                        const CFG::BoltzmannElectron&         a_ne,
                        CFG::LevelData<CFG::FArrayBox>&       a_phi_tilde ) const
{
   CFG::LevelData<CFG::FArrayBox> Zni(m_phase_geometry->magGeom().grids(), 1, CFG::IntVect::Zero);
   computeIonChargeDensity( Zni, a_kinetic_species );

   if ( CFG::GKPoissonBoltzmann* gkpb = dynamic_cast<CFG::GKPoissonBoltzmann*>(m_poisson) ) {
      gkpb->getPhiTilde( Zni, a_ne, a_phi_tilde );
   }
   else if ( CFG::NewGKPoissonBoltzmann* new_gkpb = dynamic_cast<CFG::NewGKPoissonBoltzmann*>(m_poisson) ) {
      new_gkpb->getPhiTilde( Zni, a_ne, a_phi_tilde );
   }
   else {
      MayDay::Error("GKOps::computeGKPPhiTildeDivergence() can only be called with Boltzmann electrons");
   }
}


void
GKOps::updateAveragedEfield( CFG::LevelData<CFG::FluxBox>&   a_Er_average_face,
                             CFG::LevelData<CFG::FArrayBox>& a_Er_average_cell,
                             CFG::LevelData<CFG::FArrayBox>& a_flux_divergence,
                             double&                         a_flux_divergence_lo,
                             double&                         a_flux_divergence_hi,
                             const double                    dt ) const
{
    //Geometry parameters
    const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
    const CFG::MagCoordSys& coords = *mag_geom.getCoordSys();
    const CFG::MagBlockCoordSys& block0_coord_sys = (const CFG::MagBlockCoordSys&)(*(coords.getCoordSys(0)));
    int lo_radial_index = block0_coord_sys.domain().domainBox().smallEnd(RADIAL_DIR);
    int hi_radial_index = block0_coord_sys.domain().domainBox().bigEnd(RADIAL_DIR);
    
    //Create face-centered increment for the mapped E-field (- grad Phi )
    CFG::LevelData<CFG::FluxBox> Er_mapped_face(mag_geom.grids(), 3, 2*CFG::IntVect::Unit );
    
    //Creating tmp arrays with extra layers of ghost cells
    CFG::LevelData<CFG::FArrayBox> gkp_divergence_tmp( mag_geom.grids(), 1, 3*CFG::IntVect::Unit );
    CFG::LevelData<CFG::FArrayBox> flux_divergence_tmp( mag_geom.grids(), 1, 3*CFG::IntVect::Unit );
    
    setZero(gkp_divergence_tmp);
    setZero(flux_divergence_tmp);
    
    for (CFG::DataIterator dit(gkp_divergence_tmp.dataIterator()); dit.ok(); ++dit) {
        gkp_divergence_tmp[dit].copy(m_radial_gkp_divergence_average[dit]);        
        flux_divergence_tmp[dit].copy(a_flux_divergence[dit]);
    }

    gkp_divergence_tmp.exchange();
    flux_divergence_tmp.exchange();
    
    for (CFG::DataIterator dit(a_Er_average_face.dataIterator()); dit.ok(); ++dit) {
        
        Er_mapped_face[dit].setVal(0.0);
        
        CFG::FArrayBox& this_gkp_divergence = gkp_divergence_tmp[dit];
        CFG::FArrayBox& this_flux_divergence = flux_divergence_tmp[dit];
        
        int block_number( coords.whichBlock( mag_geom.grids()[dit] ) );
        
        if ( block_number < 2 ) {
            
            for (int dir = 0; dir < 2; dir++) {
                
                CFG::FArrayBox& this_Er_mapped_dir = Er_mapped_face[dit][dir];
                
                CFG::Box box( this_Er_mapped_dir.box() );
                CFG::BoxIterator bit(box);
                for (bit.begin(); bit.ok(); ++bit) {
                    
                    CFG::IntVect iv = bit();
                    
                    if (iv[0]<lo_radial_index) {
                        this_Er_mapped_dir(iv,0) =   dt * (-a_flux_divergence_lo / m_lo_radial_gkp_divergence_average);
                    }
                    else if (iv[0]>hi_radial_index) {
                        this_Er_mapped_dir(iv,0) =  dt * (-a_flux_divergence_hi / m_hi_radial_gkp_divergence_average);
                    }
                    else {
                        
                        CFG::IntVect iv_tmp = bit();
                        iv_tmp[1] = mag_geom.grids()[dit].smallEnd(POLOIDAL_DIR);
                        
                        if ( dir == RADIAL_DIR ) {
                            this_Er_mapped_dir(iv,0) =   dt * ( - this_flux_divergence(iv_tmp,0) / this_gkp_divergence(iv_tmp,0) );
                        }
                        
                        if ( dir == POLOIDAL_DIR ) {
                            CFG::IntVect iv_incr = iv_tmp;
                            iv_incr[0] = iv_tmp[0] + 1;
                            if (iv[0]<hi_radial_index) {
                                this_Er_mapped_dir(iv,0) =   0.5 * dt * ( - this_flux_divergence(iv_tmp,0) / this_gkp_divergence(iv_tmp,0)
                                                                         - this_flux_divergence(iv_incr,0) / this_gkp_divergence(iv_incr,0) );
                            }
                            else  {
                                this_Er_mapped_dir(iv,0) =   0.5 * dt * ( - this_flux_divergence(iv_tmp,0) / this_gkp_divergence(iv_tmp,0)
                                                                         -a_flux_divergence_hi / m_hi_radial_gkp_divergence_average );
                            }
                            
                        }
                        
                    }
                }
            }
        }
        
        else {
            //temporarily setting E-field to zero on open field lines
            double E_open = 0.0;
            Er_mapped_face[dit].setVal(E_open,RADIAL_DIR,0,1);
            Er_mapped_face[dit].setVal(E_open,POLOIDAL_DIR,0,1);
        }
    }
    Er_mapped_face.exchange();
    
    //Compute cell-centered increment for the mapped E-field
    CFG::LevelData<CFG::FArrayBox> Er_mapped_cell(mag_geom.grids(), 3, 2*CFG::IntVect::Unit );
    CFG::DataIterator dit( Er_mapped_cell.dataIterator() );
    for (dit.begin(); dit.ok(); ++dit) {
        CFG::Box box( Er_mapped_cell[dit].box() );
        CFG::BoxIterator bit(box);
        for (bit.begin(); bit.ok(); ++bit) {
            CFG::IntVect iv = bit();
            CFG::IntVect iv_r = bit();
            CFG::IntVect iv_th = bit();
            iv_r[0] = iv[0] + 1;
            iv_th[1] = iv[1] + 1;
            Er_mapped_cell[dit](iv,0) = 0.25*(Er_mapped_face[dit][1](iv,0)+Er_mapped_face[dit][1](iv_th,0)+Er_mapped_face[dit][0](iv,0)+Er_mapped_face[dit][0](iv_r,0));
            Er_mapped_cell[dit](iv,1) = 0.25*(Er_mapped_face[dit][1](iv,1)+Er_mapped_face[dit][1](iv_th,1)+Er_mapped_face[dit][0](iv,1)+Er_mapped_face[dit][0](iv_r,1));
            Er_mapped_cell[dit](iv,2) = 0.25*(Er_mapped_face[dit][1](iv,2)+Er_mapped_face[dit][1](iv_th,2)+Er_mapped_face[dit][0](iv,2)+Er_mapped_face[dit][0](iv_r,2));
        }
    }
    
    //Multiply by NJtranspose (E_fs_average = - NJInverse * grad Phi)
    CFG::LevelData<CFG::FluxBox> Er_face_incr( mag_geom.grids(), 3, 2*CFG::IntVect::Unit );
    m_poisson->getUnmappedField(Er_mapped_face, Er_face_incr);
    
    CFG::LevelData<CFG::FArrayBox> Er_cell_incr( mag_geom.grids(), 3, 2*CFG::IntVect::Unit );
    m_poisson->getUnmappedField(Er_mapped_cell, Er_cell_incr);

    //Update E-field
    for (CFG::DataIterator dit(a_Er_average_face.dataIterator()); dit.ok(); ++dit) {
        if (dt>0) {
            a_Er_average_face[dit] += Er_face_incr[dit];
            a_Er_average_cell[dit] += Er_cell_incr[dit];
        }
        
        else {
          Real inv_dt = 1.0 / dt;
           a_Er_average_face[dit].copy(Er_face_incr[dit]);
           a_Er_average_face[dit] *= inv_dt;
           a_Er_average_cell[dit].copy(Er_cell_incr[dit]);
           a_Er_average_cell[dit].mult(inv_dt);
        }
    }
    
    //Extrapolate E_field into the open-field-line region (presently, constant extrapolation)
    if ( typeid(coords) == typeid(CFG::SingleNullCoordSys) ) {
      // ((CFG::SingleNullCoordSys&)coords).extrapolateEfield(a_Er_average_face, a_Er_average_cell);
    }
}


void
GKOps::updateEfieldPoloidalVariation( const CFG::LevelData<CFG::FluxBox>&   a_E_tilde_mapped_face,
                                      const CFG::LevelData<CFG::FArrayBox>& a_E_tilde_mapped_cell,
                                      CFG::LevelData<CFG::FluxBox>&         a_E_tilde_phys_face,
                                      CFG::LevelData<CFG::FArrayBox>&       a_E_tilde_phys_cell,
                                      CFG::LevelData<CFG::FArrayBox>&       a_phi_tilde_fs_average,
                                      double&                               a_phi_tilde_fs_average_lo,
                                      double&                               a_phi_tilde_fs_average_hi ) const
{
   //values in the ghost cells shold not matter since it's overwritten in the GK-velocity calculation
   for (CFG::DataIterator dit(a_E_tilde_phys_face.dataIterator()); dit.ok(); ++dit) {
      a_E_tilde_phys_face[dit].setVal(0.0);
   }

   for (CFG::DataIterator dit(a_E_tilde_phys_cell.dataIterator()); dit.ok(); ++dit) {
      a_E_tilde_phys_cell[dit].setVal(0.0);
   }

   if ( !m_ampere_cold_electrons ) {
      m_poisson->getUnmappedField(a_E_tilde_mapped_face, a_E_tilde_phys_face);
      m_poisson->getUnmappedField(a_E_tilde_mapped_cell, a_E_tilde_phys_cell);
      
      CFG::LevelData<CFG::FluxBox> Er_average_face_pol_contrib(m_phase_geometry->magGeom().grids(), 3, 2*CFG::IntVect::Unit);
      CFG::LevelData<CFG::FArrayBox> Er_average_cell_pol_contrib(m_phase_geometry->magGeom().grids(), 3, 2*CFG::IntVect::Unit);
 
      updateAveragedEfield( Er_average_face_pol_contrib, Er_average_cell_pol_contrib,
                            a_phi_tilde_fs_average, a_phi_tilde_fs_average_lo, a_phi_tilde_fs_average_hi, -1.0 );
       
      for (CFG::DataIterator dit(a_E_tilde_phys_face.dataIterator()); dit.ok(); ++dit) {
         a_E_tilde_phys_face[dit] -= Er_average_face_pol_contrib[dit];
         a_E_tilde_phys_cell[dit] -= Er_average_cell_pol_contrib[dit];
      }
   }
}

void
GKOps::setErAverage( const LevelData<FluxBox>& Er_face_injected )
{
   m_Er_average_face.define(m_phase_geometry->magGeom().grids(), 3, 2*CFG::IntVect::Unit);
   m_phase_geometry->projectPhaseToConfiguration(Er_face_injected, m_Er_average_face);
}

void
GKOps::setErAverage( const LevelData<FArrayBox>& Er_cell_injected)
{
   m_Er_average_cell.define(m_phase_geometry->magGeom().grids(), 3, 2*CFG::IntVect::Unit);
   m_phase_geometry->projectPhaseToConfiguration(Er_cell_injected, m_Er_average_cell);
}

void
GKOps::setETilde( const LevelData<FluxBox>& E_tilde_face_injected )
{
   m_E_tilde_face.define(m_phase_geometry->magGeom().grids(), 3, 2*CFG::IntVect::Unit);
   m_phase_geometry->projectPhaseToConfiguration(E_tilde_face_injected, m_E_tilde_face);
}

void
GKOps::setETilde( const LevelData<FArrayBox>& E_tilde_cell_injected)
{
   m_E_tilde_cell.define(m_phase_geometry->magGeom().grids(), 3, 2*CFG::IntVect::Unit);
   m_phase_geometry->projectPhaseToConfiguration(E_tilde_cell_injected, m_E_tilde_cell);
}

void GKOps::divideJ( const KineticSpeciesPtrVect& a_soln_mapped,
                     KineticSpeciesPtrVect&       a_soln_physical )
{
   CH_assert( isDefined() );
   for (int species(0); species<a_soln_physical.size(); species++) {
      const KineticSpecies& soln_species_mapped( *(a_soln_mapped[species]) );
      KineticSpecies& soln_species_physical( *(a_soln_physical[species]) );

      const LevelData<FArrayBox> & dfn_mapped = soln_species_mapped.distributionFunction();
      LevelData<FArrayBox> & dfn_physical = soln_species_physical.distributionFunction();

      DataIterator dit = dfn_physical.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
         dfn_physical[dit].copy( dfn_mapped[dit] );
      }

      CH_assert( m_phase_geometry != NULL );
      m_phase_geometry->divideJonValid( dfn_physical );
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

   if (a_ppgksys.contains("consistent_potential_bcs")) {
      a_ppgksys.query( "consistent_potential_bcs", m_consistent_potential_bcs );
   }
   else {
      m_consistent_potential_bcs = false;
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
      CFG::LevelData<CFG::FArrayBox> temp( mag_geom.grids(),
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
}


void GKOps::createGKPoisson( const CFG::LevelData<CFG::FArrayBox>& a_initial_ion_charge_density )
{
   if (m_poisson) {
      MayDay::Error( "GKPoisson has already been created" );
   }

   const double larmor( m_units->larmorNumber() );
   const double debye( m_units->debyeNumber() );
   CH_assert( m_phase_geometry != NULL );
   const CFG::MagGeom& mag_geom( m_phase_geometry->magGeom() );
   if ( m_boltzmann_electron ) {
      ParmParse pp( CFG::GKPoissonBoltzmann::pp_name );
      if ( typeid(*(mag_geom.getCoordSys())) == typeid(CFG::SingleNullCoordSys) ) {
         m_poisson = new CFG::NewGKPoissonBoltzmann( pp, mag_geom, larmor, debye, a_initial_ion_charge_density );
      }
      else {
         m_poisson = new CFG::GKPoissonBoltzmann( pp, mag_geom, larmor, debye, a_initial_ion_charge_density );
      }
   }
   else {
      ParmParse pp( CFG::GKPoisson::pp_name );
      m_poisson = new CFG::GKPoisson( pp, mag_geom, larmor, debye );
   }
}


void GKOps::plotPotential( const std::string& a_filename,
                           const double&      a_time ) const
{
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   m_phase_geometry->plotConfigurationData( a_filename.c_str(), m_phi, a_time );
}

void GKOps::plotEField( const std::string& a_filename,
                        const double&      a_time ) const
{
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );

   if (m_poisson) { // Plot the psi and theta projections of the physical field
      const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
      CFG::LevelData<CFG::FluxBox> Efield( mag_geom.grids(), 2, CFG::IntVect::Zero);

      const CFG::PotentialBC& bc = m_boundary_conditions->getPotentialBC();

#if 1
      m_poisson->computeFaceCenteredField( m_phi, bc, Efield );
#else
      ((CFG::NewGKPoissonBoltzmann*)m_poisson)->getFaceCenteredFieldOnCore( m_phi, Efield );
#endif
      
      const CFG::DisjointBoxLayout& grids = mag_geom.grids();
      for (CFG::DataIterator dit(Efield.dataIterator()); dit.ok(); ++dit) {
        CFG::FluxBox& this_Efield = Efield[dit];
        mag_geom.getBlockCoordSys(grids[dit]).computePsiThetaProjections(this_Efield);
      }

      phase_geometry.plotConfigurationData( a_filename.c_str(), Efield, a_time );
      //      phase_geometry.plotConfigurationData( a_filename.c_str(), m_E_field_cell, a_time ); 

       
   }
   else { // Plot the unmapped field
      phase_geometry.plotConfigurationData( a_filename.c_str(), m_E_field_cell, a_time );
   }
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

   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> density( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.numberDensity( density );

   CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.ParallelMomentum( ParallelMom );

   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      ParallelMom[dit].divide(density[dit]);
   }
   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.pressureMoment(pressure, ParallelMom);

   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      pressure[dit].divide(density[dit]);
   }

   const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();
   LevelData<FArrayBox> deltaF( species_geometry.gridsFull(), 1, IntVect::Zero );
   DeltaFKernel DeltaF_Kernel(density, pressure, ParallelMom);
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
                             const double&         a_time ) const
{
   const PhaseGeom& species_geometry = a_soln_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& soln_dfn = a_soln_species.distributionFunction();
   LevelData<FArrayBox> dfn( species_geometry.gridsFull(), 1, IntVect::Zero );
   for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      dfn[dit].copy(soln_dfn[dit]);
   }

   CFG::IntVect cspace_index(a_radial_index, a_poloidal_index);
   species_geometry.divideBStarParallel( dfn );
   species_geometry.plotAtConfigurationIndex( a_filename.c_str(),
                                              cspace_index, dfn, a_time );
}



void GKOps::plotChargeDensity( const std::string&    a_filename,
                               const KineticSpecies& a_soln_species,
                               const double&         a_time ) const
{
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );

   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   CFG::LevelData<CFG::FArrayBox> charge_density( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.chargeDensity( charge_density );

   phase_geometry.plotConfigurationData( a_filename.c_str(), charge_density, a_time );
}


void GKOps::plotChargeDensity( const std::string&     a_filename,
                               KineticSpeciesPtrVect& a_species,
                               const double&          a_time ) const
{
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   CFG::LevelData<CFG::FArrayBox> charge_density( mag_geom.grids(), 1, CFG::IntVect::Zero );
   CFG::LevelData<CFG::FArrayBox> species_charge_density( mag_geom.grids(), 1, CFG::IntVect::Zero );

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


void GKOps::plotParallelMomentum( const std::string&    a_filename,
                                  const KineticSpecies& a_soln_species,
                                  const double&         a_time ) const
{
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   CFG::LevelData<CFG::FArrayBox> parallel_mom( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.ParallelMomentum(parallel_mom);
   phase_geometry.plotConfigurationData( a_filename.c_str(), parallel_mom, a_time );
}


void GKOps::plotPoloidalMomentum( const std::string&    a_filename,
                                  const KineticSpecies& a_soln_species,
                                  const double&         a_time ) const
{
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const double larmor_number( m_units->larmorNumber() );
   LevelData<FluxBox> E_field_tmp;
   phase_geometry.injectConfigurationToPhase( m_E_field_face,
                                              m_E_field_cell,
                                              E_field_tmp );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> poloidal_mom( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.PoloidalMomentum( poloidal_mom, E_field_tmp, larmor_number );

   phase_geometry.plotConfigurationData( a_filename.c_str(), poloidal_mom, a_time );
}


void GKOps::plotPressure( const std::string&    a_filename,
                          const KineticSpecies& a_soln_species,
                          const double&         a_time ) const
{
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> density( mag_geom.grids(), 1, CFG::IntVect::Zero );
   CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.numberDensity( density );
   a_soln_species.ParallelMomentum( ParallelMom );

   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      ParallelMom[dit].divide(density[dit]);
   }
   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.pressureMoment(pressure, ParallelMom);

   phase_geometry.plotConfigurationData( a_filename.c_str(), pressure, a_time );
}


void GKOps::plotTemperature( const std::string&    a_filename,
                             const KineticSpecies& a_soln_species,
                             const double&         a_time ) const
{
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   CFG::LevelData<CFG::FArrayBox> density( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.numberDensity( density );

   CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.ParallelMomentum( ParallelMom );
   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      ParallelMom[dit].divide(density[dit]);
   }

   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.pressureMoment(pressure, ParallelMom);

   CFG::LevelData<CFG::FArrayBox> temperature;
   temperature.define( pressure );

   for (CFG::DataIterator dit(temperature.dataIterator()); dit.ok(); ++dit) {
      temperature[dit].divide(density[dit]);
   }

   phase_geometry.plotConfigurationData( a_filename.c_str(), temperature, a_time );
}


void GKOps::plotFourthMoment( const std::string&    a_filename,
                              const KineticSpecies& a_soln_species,
                              const double&         a_time ) const
{
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   CFG::LevelData<CFG::FArrayBox> density( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.numberDensity( density );

   CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.ParallelMomentum( ParallelMom );
   for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      ParallelMom[dit].divide(density[dit]);
   }

   CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.pressureMoment(pressure, ParallelMom);

   CFG::LevelData<CFG::FArrayBox> fourthMoment( mag_geom.grids(), 1, CFG::IntVect::Zero);
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
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   LevelData<FluxBox> E_field_tmp;
   phase_geometry.injectConfigurationToPhase( m_E_field_face,
                                              m_E_field_cell,
                                              E_field_tmp );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );
   
   CFG::LevelData<CFG::FArrayBox> particle_flux( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.ParticleFlux( particle_flux, E_field_tmp );

   phase_geometry.plotConfigurationData( a_filename.c_str(), particle_flux, a_time );
}


void GKOps::plotHeatFlux( const std::string&    a_filename,
                          const KineticSpecies& a_soln_species,
                          const double&         a_time ) const
{
   CH_assert( isDefined() );
   CH_assert( m_phase_geometry != NULL );
   const PhaseGeom& phase_geometry( *m_phase_geometry );
   LevelData<FluxBox> E_field_tmp;
   LevelData<FArrayBox> phi_injected_tmp;
   phase_geometry.injectConfigurationToPhase( m_E_field_face,
                                              m_E_field_cell,
                                              E_field_tmp );
   phase_geometry.injectConfigurationToPhase( m_phi, phi_injected_tmp );
   const CFG::MagGeom& mag_geom( phase_geometry.magGeom() );

   CFG::LevelData<CFG::FArrayBox> heat_flux( mag_geom.grids(), 1, CFG::IntVect::Zero );
   a_soln_species.HeatFlux( heat_flux, E_field_tmp, phi_injected_tmp );

   phase_geometry.plotConfigurationData( a_filename.c_str(), heat_flux, a_time );
}


void GKOps::setupFieldHistories( ParmParse& a_ppsim )
{
   CH_assert( isDefined() );
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
   CH_assert( isDefined() );
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
         field = &m_E_field_cell;
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
         fieldout << "# " << inode_hist[0] << ", " << inode_hist[1] << endl;

         //overwrite any existing file; this one is for human viewing (1st column for the step)
         ofstream fieldout_step( fname_step.c_str(), ios::out ); 
         fieldout_step << "# step  time  potential at ix,iy = " << inode_hist[0] << ", " << inode_hist[1] <<  endl;

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


void GKOps::plotVlasovDivergence( const std::string&    a_filename,
                                  const KineticSpecies& a_soln_species,
                                  const double&         a_time)
{
   CH_assert( isDefined() );
#if 0
   // This is quite broken
   
   // Make a temporary species vector for the output
   KineticSpeciesPtrVect tmp_species_vect(1);
   tmp_species_vect[0] = a_soln_species.clone( IntVect::Zero, false );

   LevelData<FluxBox> E_field;
   computeElectricField( E_field, tmp_species_vect, 0, 4 );

  // Make another temporary species with ghost cells, divide by J,
  // and apply Vlasov operator
  KineticSpeciesPtrVect species_phys;
  createTemporarySpeciesVector( species_phys, tmp_species_vect );
  fillGhostCells( species_phys, E_field, 0. );
  applyVlasovOperator( tmp_species_vect, species_phys, E_field, 0. );

    KineticSpecies& this_species( *(tmp_species_vect[0]) );

    // Get the velocity index at (v_parallel_max, 0)
    const PhaseGeom& geom = this_species.phaseSpaceGeometry();
    const Box& domainBox = geom.domain().domainBox();
    VEL::IntVect vpt(domainBox.bigEnd(VPARALLEL_DIR),0);

#if 1
    // Construct filename
    string filename("vlasov_divergence." + this_species.name());

    geom.plotAtVelocityIndex( filename.c_str(), vpt, this_species.distributionFunction() a_time );

#else

    // Slice in the mu direction at the low mu coordinate
    SliceSpec slice_mu(MU_DIR,domainBox.smallEnd(MU_DIR));

    CP1::Box sliced_box;
    sliceBox(sliced_box, domainBox, slice_mu);

    CP1::LevelData<D3::FArrayBox> sliced_divergence;
    sliceLevelData(sliced_divergence, this_species.distributionFunction(), slice_mu);

    Vector<string> names;
    names.push_back("velocity_divergence");
    Vector<CP1::DisjointBoxLayout> layouts;
    layouts.push_back(sliced_divergence.disjointBoxLayout());

    Vector<CP1::LevelData<CP1::FArrayBox>* > data;
    data.push_back(&sliced_divergence);

    double dx = 1.;
    double dt = 1.;
    Vector<int> ratio;
    for (int dir=0; dir<CFG_DIM+1; ++dir) {
      ratio.push_back(1);
    }
    int num_levels = 1;

    // Construct filename
    string filename("vlasov_divergence." + this_species.name() + ".hdf5");

    WriteAMRHierarchyHDF5(filename.c_str(), layouts, data, names, sliced_box,
                          dx, dt, a_time, ratio, num_levels);
#endif
#endif
}

void GKOps::writeCheckpointFile( HDF5Handle& a_handle ) const
{
   CH_assert( isDefined() );
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
   CH_assert( isDefined() );
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
      const CFG::DisjointBoxLayout& grids( mag_geom.grids() );
      CFG::LevelData<CFG::FArrayBox> ion_density( grids, 1, CFG::IntVect::Zero );
      CFG::LevelData<CFG::FArrayBox> electron_density( grids, 1, CFG::IntVect::Zero );
      
      computeSignedChargeDensities( ion_density, electron_density, a_species );
 
      CFG::LevelData<CFG::FArrayBox> quasineutral_density( grids, 1, m_cfg_nghosts );
      computeQuasiNeutralElectronDensity( quasineutral_density,
                                          a_potential,
                                          m_boundary_conditions->getPotentialBC(),
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
   const CFG::PotentialBC&               a_bc, 
   const CFG::LevelData<CFG::FArrayBox>& a_ion_density) const
{
   a_ion_density.copyTo( a_quasineutral_density );

   m_poisson->setOperatorCoefficients( a_ion_density, a_bc );
   const CFG::DisjointBoxLayout& grids( a_ion_density.disjointBoxLayout() );
   CFG::LevelData<CFG::FArrayBox> polarization_density( grids, 1, m_cfg_nghosts );

   const CFG::PotentialBC& bc = m_boundary_conditions->getPotentialBC();
   m_poisson->applyOperator( a_potential, bc, polarization_density );
   for (CFG::DataIterator cdit( grids.dataIterator() ); cdit.ok(); ++cdit) {
      a_quasineutral_density[cdit] -= polarization_density[cdit];
   }
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



#include "NamespaceFooter.H"

