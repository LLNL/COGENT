#include "BoundaryLookupTable.H"

#include "NamespaceHeader.H"

BoundaryLookupTable*
BoundaryLookupTable::s_lookup_table_instance = (BoundaryLookupTable *)NULL;

/// PUBLIC METHODS /////////////////////////////////////////////////////////////////
const BoundaryLookupTable& BoundaryLookupTable::getLookupTable()
{
   if (!s_lookup_table_instance) {
      s_lookup_table_instance = new BoundaryLookupTable();
   }
   return *(s_lookup_table_instance);
}


const Vector<int>& BoundaryLookupTable::getDirections( const int a_iloc,
                                                       const int a_codim ) const
{
   CH_assert((a_codim > 0) && (a_codim <= CH_SPACEDIM));
   CH_assert((a_iloc >= 0) && (a_iloc < m_ncase[a_codim-1]));
   int icomb( iComb( a_iloc, a_codim ) );
   return m_dirs_table[a_codim-1][icomb];
}


const Vector<Side::LoHiSide>& BoundaryLookupTable::getSides( const int a_iloc,
                                                             const int a_codim ) const
{
   CH_assert((a_codim > 0) && (a_codim <= CH_SPACEDIM));
   CH_assert((a_iloc >= 0) && (a_iloc < m_ncase[a_codim-1]));
   return m_side_table[a_codim-1][a_iloc];
}

void BoundaryLookupTable::printTable( std::ostream& a_out ) const
{
   a_out << "codim\tdir\tloc index\tcomb index\tperm index" << std::endl;
   for (int icodim(0); icodim< CH_SPACEDIM; icodim++) {
      for (int iloc(0); iloc<m_ncase[icodim]; iloc++) {
         a_out << icodim+1 << "\t";
         const Vector<int>& tuple( getDirections( iloc, icodim+1 ) );
         for (int d(0); d<tuple.size(); d++) {
            a_out << tuple[d];
         }
         a_out << "\t"
               << iloc << "\t\t"
               << iComb(iloc, icodim+1) << "\t\t"
               << iPerm(iloc, icodim+1) << "\t";
         const Vector<Side::LoHiSide>& sides( getSides( iloc, icodim+1 ) );
         for (int d(0); d<sides.size(); d++) {
            a_out << sides[d];
         }
         a_out << std::endl;
      }
   }
}
/// PUBLIC METHODS /////////////////////////////////////////////////////////////////


/// PRIVATE METHODS //////////////////////////////////////////////////////////////////
BoundaryLookupTable::BoundaryLookupTable()
{
   initializeCounts();
   buildDirsVectors();
   buildSideVectors();
}


BoundaryLookupTable::~BoundaryLookupTable()
{
   for (int i=0; i<m_side_table.size(); ++i) {
      for (int j=0; j<m_side_table[i].size(); ++j) {
         m_side_table[i][j].clear();
      }
      m_side_table[i].clear();
   }
   m_side_table.clear();
}


void BoundaryLookupTable::incrementDirTuple( Vector<int>& a_tuple,
                                             const int    a_index )
{
   if (a_index>=0) {
      const int codim(a_tuple.size());
      a_tuple[a_index]++;
      if (a_tuple[a_index]==CH_SPACEDIM - codim + 1 + a_index) {
         incrementDirTuple( a_tuple, a_index-1 );
         a_tuple[a_index] = a_tuple[a_index-1] + 1;
      }
   }
}
/// PRIVATE METHODS //////////////////////////////////////////////////////////////////


/// PRIVATE INLINE METHODS //////////////////////////////////////////////////////////
inline
int factorial( int a_n )
{
   CH_assert(a_n>=0);
   int result(1);
   for (int i(1); i<=a_n; i++) result *= i;
   return result;
}

inline
int nchoosek( int a_n, int a_k )
{
   CH_assert(a_n>=0);
   CH_assert(a_k>=0);
   return factorial(a_n) / ( factorial(a_k) * factorial(a_n-a_k) );
}

inline
bool BoundaryLookupTable::isUpper( const int a_iloc,
                                   const int a_codim,
                                   const int a_index ) const
{
   CH_assert((a_codim > 0) && (a_codim <= CH_SPACEDIM));
   CH_assert((a_iloc >= 0) && (a_iloc < m_ncase[a_codim-1]));
   CH_assert((a_index >= 0) && (a_index < a_codim));
   const int iperm( iPerm(a_iloc, a_codim) );
   return (iperm & (1 << (a_codim-1-a_index) ));
}

inline
int BoundaryLookupTable::iComb( const int a_iloc,
                                const int a_codim ) const
{
   CH_assert((a_codim > 0) && (a_codim <= CH_SPACEDIM));
   CH_assert((a_iloc >= 0) && (a_iloc < m_ncase[a_codim-1]));
   return (a_iloc / m_nperm[a_codim-1]);
}

inline
int BoundaryLookupTable::iPerm( const int a_iloc,
                                const int a_codim ) const
{
   CH_assert((a_codim > 0) && (a_codim <= CH_SPACEDIM));
   CH_assert((a_iloc >= 0) && (a_iloc < m_ncase[a_codim-1]));
   return (a_iloc % m_nperm[a_codim-1]);
}

inline
void BoundaryLookupTable::initializeCounts()
{
   m_ncomb.resize( CH_SPACEDIM );
   m_nperm.resize( CH_SPACEDIM );
   m_ncase.resize( CH_SPACEDIM );

   for (int codim(1); codim<=CH_SPACEDIM; codim++) {
      int cdm1(codim-1);
      m_ncomb[cdm1] = nchoosek( CH_SPACEDIM, codim );
      m_nperm[cdm1] = (1 << codim);
      m_ncase[cdm1] = m_ncomb[cdm1] * m_nperm[cdm1];
   }
}


inline
void BoundaryLookupTable::buildDirsVectors()
{
   m_dirs_table.resize( CH_SPACEDIM );
   for (int icodim(0); icodim<CH_SPACEDIM; icodim++) {
      const int codim(icodim+1);
      m_dirs_table[icodim].resize( m_ncomb[icodim] );

      Vector<int> work(codim);
      for (int i(0); i<codim; i++) {
         work[i] = i;
      }
      m_dirs_table[icodim][0].resize( codim );
      m_dirs_table[icodim][0] = work;

      for (int icomb(1); icomb<m_ncomb[icodim]; icomb++) {
         incrementDirTuple( work, icodim );
         m_dirs_table[icodim][icomb].resize( codim );
         m_dirs_table[icodim][icomb] = work;
      }
   }
}


inline
void BoundaryLookupTable::buildSideVectors()
{
   m_side_table.resize( CH_SPACEDIM );
   for (int icodim(0); icodim<CH_SPACEDIM; icodim++) {
      m_side_table[icodim].resize( m_ncase[icodim] );
      for (int iloc(0); iloc<m_ncase[icodim]; iloc++) {
         const int codim( icodim + 1 );
         m_side_table[icodim][iloc].resize( codim );
         for (int j(0); j<codim; j++) {
            if (isUpper(iloc, codim, j)) {
               m_side_table[icodim][iloc][j] = Side::Hi;
            }
            else {
               m_side_table[icodim][iloc][j] = Side::Lo;
            }
         }
      }
   }
}
/// PRIVATE INLINE METHODS //////////////////////////////////////////////////////////


#include "NamespaceFooter.H"
