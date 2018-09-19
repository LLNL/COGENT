#include "SNTimeIntegration.H"
#include "inspect.H"

#include "SNUtils.H"

#include "NamespaceHeader.H"

////////////////////////////////////////////////////////////////////

inline
void copyData( FluidSpeciesPtrVect& a_species_dst,
               const FluidSpeciesPtrVect& a_species_src )
{
   a_species_dst.resize( a_species_src.size() );
   for (int s(0); s<a_species_dst.size(); s++) {
      a_species_dst[s]->copy( *(a_species_src[s]) );
   } 
}
   
void SNState::copy( const SNState& a_state )
{
   copyData( m_species, a_state.data() );
}

void SNState::copy( const SNRHSData& a_rhs )
{
   copyData( m_species, a_rhs.data() );
}

void SNRHSData::copy( const SNRHSData& a_rhs )
{
   copyData( m_species, a_rhs.data() );
}

void SNRHSData::copy( const SNState& a_state )
{
   copyData( m_species, a_state.data() );
}

////////////////////////////////////////////////////////////////////

inline
int computeVectorSize( const FluidSpeciesPtrVect& a_species )
{
   int total_size(0);
   for (int s(0); s<a_species.size(); s++) {
      total_size += a_species[s]->size();
   }
   return total_size;
}

int SNState::getVectorSize() const
{
   CH_assert( isDefined() );
   return computeVectorSize( m_species );
}

int SNRHSData::getVectorSize()
{
   CH_assert( isDefined() );
   return computeVectorSize( m_species );
}

////////////////////////////////////////////////////////////////////

inline
void copyToArray( Real* a_dst, const FluidSpeciesPtrVect& a_species )
{
   int offset(0);
   for (int s(0); s<a_species.size(); s++) {
      offset += SNUtils::copyFromLevelData( a_dst + offset,
                                                 a_species[s]->data() );
   }
}

void SNState::copyTo( Real* a_vector ) const
{
   CH_assert( isDefined() );
   copyToArray( a_vector,  m_species );
}

void SNRHSData::copyTo( Real* a_vector ) const 
{
   CH_assert( isDefined() );
   copyToArray( a_vector,  m_species );
}

////////////////////////////////////////////////////////////////////

inline
void copyFromArray( FluidSpeciesPtrVect& a_species,
                    const Real* a_src )
{
   int offset(0);
   for (int s(0); s<a_species.size(); s++) {
      offset += SNUtils::copyToLevelData( a_species[s]->data(),
                                               a_src + offset );
   }
}

void SNState::copyFrom( const Real* a_vector )
{
   CH_assert( isDefined() );
   copyFromArray( m_species, a_vector );
}

void SNRHSData::copyFrom( const Real* a_vector )
{
   CH_assert( isDefined() );
   copyFromArray( m_species, a_vector );
}

////////////////////////////////////////////////////////////////////

inline
void addFromArray( FluidSpeciesPtrVect& a_species,
                   Real* a_src,
                   const Real& a_factor )
{
   int offset(0);
   for (int s(0); s<a_species.size(); s++) {
      offset += SNUtils::addToLevelData( a_species[s]->data(),
                                              a_src + offset,
                                              a_factor );
   }
}

void SNState::addFrom( Real* a_vector, const Real& a_factor )
{
   CH_assert( isDefined() );
   addFromArray( m_species, a_vector, a_factor );
}

////////////////////////////////////////////////////////////////////

inline
void scaleData( FluidSpeciesPtrVect& a_species, const Real& a_factor )
{
   for (int s(0); s<a_species.size(); s++) {
      SNUtils::scaleLevelData( a_species[s]->data(), a_factor );
   }
}

void SNState::scale( const Real& a_factor )
{
   CH_assert( isDefined() );
   scaleData( m_species, a_factor );
}

void SNRHSData::scale( const Real& a_factor )
{
   CH_assert( isDefined() );
   scaleData( m_species, a_factor );
}

////////////////////////////////////////////////////////////////////

inline
void incrementData( FluidSpeciesPtrVect& a_species,
                    const FluidSpeciesPtrVect& a_dspecies,
                    const Real& a_factor )
{
   CH_assert( a_species.size()==a_dspecies.size() );
   for (int s(0); s<a_species.size(); s++) {
      a_species[s]->addData( *(a_dspecies[s]), a_factor );
   }
}

void SNState::increment( const SNState& a_delta,
                         const Real& a_factor,
                         const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_species, a_delta.data(), a_factor );
}

void SNState::increment( const SNRHSData& a_delta,
                         const Real& a_factor,
                         const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_species, a_delta.data(), a_factor );
}

void SNRHSData::increment( const SNRHSData& a_delta,
                           const Real& a_factor,
                           const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_species, a_delta.data(), a_factor );
}

void SNRHSData::increment( const SNState& a_delta,
                           const Real& a_factor,
                           const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_species, a_delta.data(), a_factor );
}

////////////////////////////////////////////////////////////////////

Real SNRHSData::dotProduct( const SNRHSData& a_vector )
{
   CH_assert( isDefined() );
   double sum_local(0);

   const  FluidSpeciesPtrVect& b_species( a_vector.data() );
   CH_assert( m_species.size()==b_species.size() );
   for (int s(0); s<m_species.size(); s++) {
      sum_local += SNUtils::innerProductLevelData(
         m_species[s]->data(),
         b_species[s]->data() );
   }
   
   double sum(0.0);

#ifdef CH_MPI
   MPI_Allreduce( &sum_local, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
   sum = sum_local;
#endif

   return sum;
}
   
////////////////////////////////////////////////////////////////////
   
inline
Real normLp( const FluidSpeciesPtrVect& a_species, const int a_p )
{
   double accum(0);
   for (int s(0); s<a_species.size(); s++) {
      accum += SNUtils::accumulateLevelData( a_species[s]->data(), a_p );
   }
   return SNUtils::reduce( accum, a_p );
}

   
Real SNState::computeNorm( const int a_p ) const
{
   CH_assert(isDefined()); 
   CH_assert((a_p <= 2) && (a_p >= 0));
   return normLp( m_species, a_p );
}

Real SNRHSData::computeNorm( const int a_p ) const
{
   CH_assert(isDefined());
   CH_assert((a_p <= 2) && (a_p >= 0));
   return normLp( m_species, a_p );
}

#include "NamespaceFooter.H"
