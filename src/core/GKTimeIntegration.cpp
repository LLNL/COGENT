#include "GKTimeIntegration.H"
#include "inspect.H"

#include "NamespaceHeader.H"

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

void GKState::computeNorms( Real& a_norm_1, Real& a_norm_2, Real& a_norm_inf)
{
  CH_assert(isDefined());
  Real local_max = 0, local_sum1 = 0, local_sum2 = 0;
  for (int s(0); s<m_species_mapped.size(); s++) {
    const LevelData<FArrayBox>& solution = m_species_mapped[s]->distributionFunction();
    const PhaseGeom& geom = m_species_mapped[s]->phaseSpaceGeometry();
    const DisjointBoxLayout& grids = solution.disjointBoxLayout();
    LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
    geom.getCellVolumes(volume);
    DataIterator dit = solution.dataIterator();
    /* L2 */
    for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox tmp(grids[dit],1);
      tmp.copy(solution[dit]); 
      tmp *= solution[dit];
      tmp *= volume[dit];
      local_sum2 += tmp.sum(grids[dit],0,1);
    }
    /* L1 and inf */
    for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox tmp(grids[dit],1);
      tmp.copy(solution[dit]);
      tmp.abs();
      double box_max = tmp.max(grids[dit]);
      if (box_max > local_max) local_max = box_max;
      tmp *= volume[dit];
      local_sum1 += tmp.sum(grids[dit],0,1);
    }
  }
  MPI_Allreduce(&local_max , &a_norm_inf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&local_sum1, &a_norm_1  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_sum2, &a_norm_2  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
  a_norm_2 = sqrt(a_norm_2);
}

void GKRHSData::computeNorms( Real& a_norm_1, Real& a_norm_2, Real& a_norm_inf)
{
  CH_assert(isDefined());
  Real local_max = 0, local_sum1 = 0, local_sum2 = 0;
  for (int s(0); s<m_species_mapped.size(); s++) {
    const LevelData<FArrayBox>& solution = m_species_mapped[s]->distributionFunction();
    const PhaseGeom& geom = m_species_mapped[s]->phaseSpaceGeometry();
    const DisjointBoxLayout& grids = solution.disjointBoxLayout();
    LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
    geom.getCellVolumes(volume);
    DataIterator dit = solution.dataIterator();
    /* L2 */
    for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox tmp(grids[dit],1);
      tmp.copy(solution[dit]); 
      tmp *= solution[dit];
      tmp *= volume[dit];
      local_sum2 += tmp.sum(grids[dit],0,1);
    }
    /* L1 and inf */
    for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox tmp(grids[dit],1);
      tmp.copy(solution[dit]);
      tmp.abs();
      double box_max = tmp.max(grids[dit]);
      if (box_max > local_max) local_max = box_max;
      tmp *= volume[dit];
      local_sum1 += tmp.sum(grids[dit],0,1);
    }
  }
  MPI_Allreduce(&local_max , &a_norm_inf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&local_sum1, &a_norm_1  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_sum2, &a_norm_2  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
  a_norm_2 = sqrt(a_norm_2);
}

#include "NamespaceFooter.H"
