#include "GKTimeIntegration.H"
#include "inspect.H"

#include "NamespaceHeader.H"

void GKState::copy( const GKState& a_state )
{
  m_species_mapped.resize( a_state.m_species_mapped.size() );
  for (int s(0); s<m_species_mapped.size(); s++) {
    m_species_mapped[s]->copy( *(a_state.m_species_mapped[s]) );
  }
}

void GKState::copy( const GKRHSData& a_rhs )
{
  const KineticSpeciesPtrVect& rhs_species_mapped( a_rhs.data() ); 
  m_species_mapped.resize( rhs_species_mapped.size() );
  for (int s(0); s<m_species_mapped.size(); s++) {
    m_species_mapped[s]->copy( *(rhs_species_mapped[s]) );
  }
}

void GKRHSData::copy( const GKRHSData& a_rhs )
{
  m_species_mapped.resize( a_rhs.m_species_mapped.size() );
  for (int s(0); s<m_species_mapped.size(); s++) {
    m_species_mapped[s]->copy( *(a_rhs.m_species_mapped[s]) );
  }
}

void GKRHSData::copy( const GKState& a_rhs )
{
  const KineticSpeciesPtrVect& rhs_species_mapped( a_rhs.data() ); 
  m_species_mapped.resize( rhs_species_mapped.size() );
  for (int s(0); s<m_species_mapped.size(); s++) {
    m_species_mapped[s]->copy( *(rhs_species_mapped[s]) );
  }
}

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


void GKRHSData::scale(const Real& a_factor)
{
   CH_assert( isDefined() );
   for (int s(0); s<m_species_mapped.size(); s++) {
      LevelData<FArrayBox>& solution = m_species_mapped[s]->distributionFunction();
      DataIterator dit = solution.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) solution[dit] *= a_factor;
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


Real GKRHSData::dotProduct(const GKRHSData& a_Y)
{
   CH_assert( isDefined() );
   Real dotProduct_local = 0, dotProduct = 0;
   const KineticSpeciesPtrVect& Y_species( a_Y.data() );
   CH_assert( m_species_mapped.size()==Y_species.size() );
   for (int s(0); s<m_species_mapped.size(); s++) {
      const LevelData<FArrayBox>& vec_a = m_species_mapped[s]->distributionFunction();
      const LevelData<FArrayBox>& vec_b = Y_species[s]->distributionFunction();
      const DisjointBoxLayout& grids = vec_a.disjointBoxLayout();
      DataIterator dit = vec_a.dataIterator();
      for (dit.begin();dit.ok();++dit) dotProduct_local += vec_a[dit].dotProduct(vec_b[dit],grids[dit]);
   }
   MPI_Allreduce(&dotProduct_local,&dotProduct,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   return dotProduct;
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

Real GKState::computeNorm( int a_ord)
{
  CH_assert(isDefined());
  CH_assert((a_ord <= 2) && (a_ord >= 0));
  Real local_norm = 0, norm = 0;

  if (a_ord == 0) {
    /* max norm */
    for (int s(0); s<m_species_mapped.size(); s++) {
      const LevelData<FArrayBox>& solution = m_species_mapped[s]->distributionFunction();
      const PhaseGeom& geom = m_species_mapped[s]->phaseSpaceGeometry();
      const DisjointBoxLayout& grids = solution.disjointBoxLayout();
      LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
      geom.getCellVolumes(volume);
      DataIterator dit = solution.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox tmp(grids[dit],1);
        tmp.copy(solution[dit]);
        tmp.abs();
        double box_max = tmp.max(grids[dit]);
        if (box_max > local_norm) local_norm = box_max;
      }
    }
    MPI_Allreduce(&local_norm,&norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  } else if (a_ord = 1) {
    /* L1 */
    for (int s(0); s<m_species_mapped.size(); s++) {
      const LevelData<FArrayBox>& solution = m_species_mapped[s]->distributionFunction();
      const PhaseGeom& geom = m_species_mapped[s]->phaseSpaceGeometry();
      const DisjointBoxLayout& grids = solution.disjointBoxLayout();
      LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
      geom.getCellVolumes(volume);
      DataIterator dit = solution.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox tmp(grids[dit],1);
        tmp.copy(solution[dit]);
        tmp.abs();
        tmp *= volume[dit];
        local_norm += tmp.sum(grids[dit],0,1);
      }
    }
    MPI_Allreduce(&local_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  } else {
    /* L2 */
    for (int s(0); s<m_species_mapped.size(); s++) {
      const LevelData<FArrayBox>& solution = m_species_mapped[s]->distributionFunction();
      const PhaseGeom& geom = m_species_mapped[s]->phaseSpaceGeometry();
      const DisjointBoxLayout& grids = solution.disjointBoxLayout();
      LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
      geom.getCellVolumes(volume);
      DataIterator dit = solution.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox tmp(grids[dit],1);
        tmp.copy(solution[dit]); 
        tmp *= solution[dit];
        tmp *= volume[dit];
        local_norm += tmp.sum(grids[dit],0,1);
      }
    }
    MPI_Allreduce(&local_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  return(norm);
}

Real GKRHSData::computeNorm( int a_ord)
{
  CH_assert(isDefined());
  CH_assert((a_ord <= 2) && (a_ord >= 0));
  Real local_norm = 0, norm = 0;

  if (a_ord == 0) {
    /* max norm */
    for (int s(0); s<m_species_mapped.size(); s++) {
      const LevelData<FArrayBox>& solution = m_species_mapped[s]->distributionFunction();
      const PhaseGeom& geom = m_species_mapped[s]->phaseSpaceGeometry();
      const DisjointBoxLayout& grids = solution.disjointBoxLayout();
      LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
      geom.getCellVolumes(volume);
      DataIterator dit = solution.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox tmp(grids[dit],1);
        tmp.copy(solution[dit]);
        tmp.abs();
        double box_max = tmp.max(grids[dit]);
        if (box_max > local_norm) local_norm = box_max;
      }
    }
    MPI_Allreduce(&local_norm,&norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  } else if (a_ord = 1) {
    /* L1 */
    for (int s(0); s<m_species_mapped.size(); s++) {
      const LevelData<FArrayBox>& solution = m_species_mapped[s]->distributionFunction();
      const PhaseGeom& geom = m_species_mapped[s]->phaseSpaceGeometry();
      const DisjointBoxLayout& grids = solution.disjointBoxLayout();
      LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
      geom.getCellVolumes(volume);
      DataIterator dit = solution.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox tmp(grids[dit],1);
        tmp.copy(solution[dit]);
        tmp.abs();
        tmp *= volume[dit];
        local_norm += tmp.sum(grids[dit],0,1);
      }
    }
    MPI_Allreduce(&local_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  } else {
    /* L2 */
    for (int s(0); s<m_species_mapped.size(); s++) {
      const LevelData<FArrayBox>& solution = m_species_mapped[s]->distributionFunction();
      const PhaseGeom& geom = m_species_mapped[s]->phaseSpaceGeometry();
      const DisjointBoxLayout& grids = solution.disjointBoxLayout();
      LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
      geom.getCellVolumes(volume);
      DataIterator dit = solution.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox tmp(grids[dit],1);
        tmp.copy(solution[dit]); 
        tmp *= solution[dit];
        tmp *= volume[dit];
        local_norm += tmp.sum(grids[dit],0,1);
      }
    }
    MPI_Allreduce(&local_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  return(norm);
}

#include "NamespaceFooter.H"
