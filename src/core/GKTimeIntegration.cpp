#include "GKTimeIntegration.H"
#include "inspect.H"

#include "NamespaceHeader.H"

int GKState::getVectorSize()
{
  CH_assert( isDefined() );
  int TotalSize = 0;
  for (int s=0; s < m_species_mapped.size(); s++) {
    const LevelData<FArrayBox>& x = m_species_mapped[s]->distributionFunction();
    const DisjointBoxLayout& dbl = x.disjointBoxLayout();
    DataIterator dit = dbl.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) TotalSize += dbl[dit].numPts();
  }
  return TotalSize;
}

void GKState::copyTo(Real *Y)
{
   CH_assert( isDefined() );
   int offset = 0;
   for (int s=0; s < m_species_mapped.size(); s++) {
     const LevelData<FArrayBox>& x = m_species_mapped[s]->distributionFunction();
     const DisjointBoxLayout& dbl = x.disjointBoxLayout();
     DataIterator dit = x.dataIterator();
     for (dit.begin(); dit.ok(); ++dit) {
       FArrayBox tmp(dbl[dit],1,(Y+offset));
       tmp.copy(x[dit]);
       offset += dbl[dit].numPts();
     }
   }
}

void GKState::copyFrom(Real *Y)
{
   CH_assert( isDefined() );
   int offset = 0;
   for (int s=0; s < m_species_mapped.size(); s++) {
     LevelData<FArrayBox>& x = m_species_mapped[s]->distributionFunction();
     const DisjointBoxLayout& dbl = x.disjointBoxLayout();
     DataIterator dit = x.dataIterator();
     for (dit.begin(); dit.ok(); ++dit) {
       const FArrayBox tmp(dbl[dit],1,(Y+offset));
       x[dit].copy(tmp);
       offset += dbl[dit].numPts();
     }
   }
}

void GKState::addFrom(Real *Y, Real a_a)
{
   CH_assert( isDefined() );
   int offset = 0;
   for (int s=0; s < m_species_mapped.size(); s++) {
     LevelData<FArrayBox>& x = m_species_mapped[s]->distributionFunction();
     const DisjointBoxLayout& dbl = x.disjointBoxLayout();
     DataIterator dit = x.dataIterator();
     for (dit.begin(); dit.ok(); ++dit) {
       const FArrayBox tmp(dbl[dit],1,(Y+offset));
       x[dit].plus(tmp,a_a);
       offset += dbl[dit].numPts();
     }
   }
}

void GKRHSData::copyTo(Real *Y)
{
   CH_assert( isDefined() );
   int offset = 0;
   for (int s=0; s < m_species_mapped.size(); s++) {
     const LevelData<FArrayBox>& x = m_species_mapped[s]->distributionFunction();
     const DisjointBoxLayout& dbl = x.disjointBoxLayout();
     DataIterator dit = x.dataIterator();
     for (dit.begin(); dit.ok(); ++dit) {
       FArrayBox tmp(dbl[dit],1,(Y+offset));
       tmp.copy(x[dit]);
       offset += dbl[dit].numPts();
     }
   }
}

void GKRHSData::copyFrom(Real *Y)
{
   CH_assert( isDefined() );
   int offset = 0;
   for (int s=0; s < m_species_mapped.size(); s++) {
     LevelData<FArrayBox>& x = m_species_mapped[s]->distributionFunction();
     const DisjointBoxLayout& dbl = x.disjointBoxLayout();
     DataIterator dit = x.dataIterator();
     for (dit.begin(); dit.ok(); ++dit) {
       const FArrayBox tmp(dbl[dit],1,(Y+offset));
       x[dit].copy(tmp);
       offset += dbl[dit].numPts();
     }
   }
}

void GKState::increment( const GKState& a_state,
                         Real a_factor,
                         bool a_update_flux_register )
{
   CH_assert( isDefined() );
   const KineticSpeciesPtrVect& rhs_species( a_state.data() );
   CH_assert( m_species_mapped.size()==rhs_species.size() );
   for (int s(0); s<m_species_mapped.size(); s++) {
      m_species_mapped[s]->addData( *(rhs_species[s]), a_factor );
   }
}

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
