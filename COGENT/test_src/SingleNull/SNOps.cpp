#include "SNOps.H"
#include "CONSTANTS.H"
#include "FORT_PROTO.H"
#include "inspect.H"

#include <fstream>
#include <sstream>

#include "LocalizedF_F.H"

#if 1  // warning, OS dependencies, will not work on all platforms
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "CH_HDF5.H"
#include "HDF5Portable.H"

#include "BoxIterator.H"
#include "LoadBalance.H"
#include "MagBlockCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "newMappedGridIO.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "inspect.H"

#include "AMRIO.H"

#include "EdgeToCell.H"
#include "FourthOrderUtil.H"
//#include "Directions.H"

#include "NamespaceHeader.H"



void SNOps::define( const SNState& a_state, const Real a_dt )
{
   CH_assert( isDefined()==false );
   CH_assert( a_dt>0.0 );
   m_dt = a_dt;

   ParmParse ppsnsys( "snsystem" );
   parseParameters( ppsnsys );

   ParmParse pp;
   m_boundary_conditions = new SNSystemBC( pp, a_state );
   m_initial_conditions = new SNSystemIC( pp, a_state );
   
   m_Y.define(a_state);
   m_rhs.define(a_state);

   ParmParse pp_advect("Advect");
   m_operator = new Advect(pp_advect);

   m_is_defined = true;
}



SNOps::~SNOps()
{
   delete m_operator;
}



Real SNOps::stableDtExpl( const SNState& a_state, const int a_step_number )
{
   CH_assert( isDefined() );
   Real dt_stable( DBL_MAX );
   dt_stable = Min( dt_stable, m_dt );
   //   dt_stable = Min( dt_stable, m_dt_vlasov );
   //   dt_stable = Min( dt_stable, m_dt_collisions );
   //   dt_stable = Min( dt_stable, m_dt_fields );
   //   dt_stable = Min( dt_stable, m_dt_fluids );
   //   if (m_transport_model_on) {
   //      dt_stable = Min( dt_stable, m_dt_transport );
   //   }
   //   if (m_neutrals_model_on) {
   //      dt_stable = Min( dt_stable, m_dt_neutrals );
   //   }
   return dt_stable;
}



Real SNOps::stableDtImEx( const SNState& a_state, const int a_step_number )
{
   CH_assert( isDefined() );
   Real dt_stable( DBL_MAX );
   dt_stable = Min( dt_stable, m_dt );
   //   dt_stable = Min( dt_stable, m_dt_vlasov );
   //   dt_stable = Min( dt_stable, m_dt_fields );
   //   dt_stable = Min( dt_stable, m_dt_fluids );
   //   if (m_neutrals_model_on) {
   //      dt_stable = Min( dt_stable, m_dt_neutrals );
   //   }
   return dt_stable;
}



Real SNOps::stableDt( const SNState& a_state, 
                      const int a_step_number )
{
   CH_assert( isDefined() );
   const FluidSpeciesPtrVect& species_comp( a_state.data() );
   return m_operator->computeDt( species_comp );
}



void SNOps::preTimeStep (const int a_step, 
                         const Real a_time, 
                         const SNState& a_state_comp, 
                         const SNState& a_state_phys)
{
   CH_assert( isDefined() );
   const FluidSpeciesPtrVect& soln_comp( a_state_comp.data() );
   m_dt = m_operator->computeDt( soln_comp );
   m_time_scale = m_operator->computeTimeScale( soln_comp );
}



void SNOps::postTimeStep(const int a_step, const Real a_time, const SNState& a_state)
{
}



void SNOps::postTimeStage(const int a_step, const Real a_time, const SNState& a_state, const int a_stage )
{
}



void SNOps::postTimeStage(const int a_step, const Real a_time, const SNVector& a_vector, const int a_stage )
{
}



void SNOps::setVelocities( FluidSpeciesPtrVect& a_species_vec,
                           const bool           a_mapped )
{
   for (int s=0; s<a_species_vec.size(); ++s) {
      FluidSpecies& species = *a_species_vec[s];

      m_operator->updateVelocity(species.configurationSpaceGeometry(),
                                 species.velocity(),
                                 a_mapped);
   }
}



void SNOps::explicitOp( SNRHSData& a_rhs,
                        const Real a_time,
                        const SNState& a_state,
                        const int a_stage )
{
   CH_assert( isDefined() );
   a_rhs.zero();
   const FluidSpeciesPtrVect& species_comp( a_state.data() );
      
   FluidSpeciesPtrVect species_phys;
   createTemporarySpeciesVector( species_phys, species_comp );

   // Set mapped velocities for filling the BCs
   setVelocities( species_phys, true );
   fillGhostCells( species_phys, a_time );

   // Set physical velocities for the operator evaluation
   setVelocities( species_phys, false );
   applyOperator( a_rhs.data(), species_phys, a_time );
}



void SNOps::explicitOp( SNRHSData& a_rhs,
                        const Real a_time,
                        const SNRHSData& a_state,
                        const int a_stage )
{
   SNState temp;
   temp.copy( a_state );
   explicitOp( a_rhs, a_time, temp, a_stage );
}



void SNOps::explicitOp( SNVector& a_rhs,
                        const Real a_time,
                        const SNVector& a_state,
                        const int a_stage )
{
  m_Y.copyFrom(a_state.data());
  explicitOp(m_rhs,a_time,m_Y,a_stage);
  m_rhs.copyTo(a_rhs.data());
}



inline
void SNOps::createTemporarySpeciesVector( FluidSpeciesPtrVect& a_out,
                                          const FluidSpeciesPtrVect& a_in )
{
   IntVect ghost_vect;
   for (int d(0); d<CFG_DIM; d++) {
      ghost_vect[d] = m_ghost_vect[d];
   }
   a_out.resize( a_in.size() );
   for (int s(0); s<a_in.size(); s++) {
      a_out[s] = a_in[s]->clone( ghost_vect );
      LevelData<FArrayBox>& out_data( a_out[s]->data() );
      const MagGeom& geometry( a_in[s]->configurationSpaceGeometry() );
      geometry.divideJonValid( out_data );
   }
}



inline
void SNOps::createTemporaryState( SNState& a_out, const SNState& a_in )
{
   createTemporarySpeciesVector( a_out.data(), a_in.data() );
}



inline
void SNOps::fillGhostCells( FluidSpeciesPtrVect& a_field_phys, const Real& a_time )
{
   CH_assert( isDefined() );
   m_boundary_conditions->fillFluidSpeciesGhostCells( a_field_phys, a_time );
}



inline
void SNOps::applyOperator( FluidSpeciesPtrVect&       a_rhs,
                           const FluidSpeciesPtrVect& a_soln,
                           const Real&                a_time )
{
   CH_assert( isDefined() );
   for (int s(0); s<a_rhs.size(); s++) {
      const FluidSpecies& soln_species( *(a_soln[s]) );
      FluidSpecies& rhs_species( *(a_rhs[s]) );
      m_operator->evalRHS( rhs_species, soln_species, a_time );
   }
}



void SNOps::divideJ( const FluidSpeciesPtrVect& a_soln_mapped,
                     FluidSpeciesPtrVect&       a_soln_physical )
{
   CH_assert( isDefined() );
   for (int species(0); species<a_soln_physical.size(); species++) {
      const FluidSpecies& soln_species_mapped( *(a_soln_mapped[species]) );
      FluidSpecies& soln_species_physical( *(a_soln_physical[species]) );
      
      const LevelData<FArrayBox>& dfn_mapped( soln_species_mapped.data() );
      LevelData<FArrayBox>& dfn_physical( soln_species_physical.data() );

      for (DataIterator dit( dfn_physical.dataIterator() ); dit.ok(); ++dit) {
         dfn_physical[dit].copy( dfn_mapped[dit] );
      }

      const MagGeom& geometry( soln_species_mapped.configurationSpaceGeometry() );
      geometry.divideJonValid( dfn_physical );
   }
}



void SNOps::divideJ( const SNState& a_soln_mapped, SNState& a_soln_physical )
{
   CH_assert( isDefined() );
   divideJ( a_soln_mapped.data(), a_soln_physical.data() );
}



void SNOps::parseParameters( ParmParse& a_ppgksys )
{
}



#include "NamespaceFooter.H"

