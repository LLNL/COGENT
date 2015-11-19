#include "GKTimeIntegration.H"
#include "inspect.H"

#include "NamespaceHeader.H"

void GKState::increment( const GKRHSData& a_rhs,
                         Real a_factor,
                         bool a_update_flux_register )
{
   CH_assert( isDefined() );
   const KineticSpeciesPtrVect& rhs_species( a_rhs.data() );
   CH_assert( m_species_mapped.size()==rhs_species.size() );
   for (int s(0); s<m_species_mapped.size(); s++) {
//      inspect( rhs_species[s]->distributionFunction() );
      m_species_mapped[s]->addData( *(rhs_species[s]), a_factor );
//      inspect( m_species_mapped[s]->distributionFunction() );
   }
}


void GKRHSData::increment( const GKRHSData& a_increment,
                           Real a_factor,
                           bool a_update_flux_register )
{
   CH_assert( isDefined() );
   const KineticSpeciesPtrVect& increment_species( a_increment.data() );
   CH_assert( m_species_mapped.size()==increment_species.size() );
   for (int s(0); s<m_species_mapped.size(); s++) {
//      inspect( increment_species[s]->distributionFunction() );
      m_species_mapped[s]->addData( *(increment_species[s]), a_factor );
//      inspect( m_species_mapped[s]->distributionFunction() );
   }
}

#include "NamespaceFooter.H"
