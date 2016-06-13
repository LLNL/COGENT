#include "REAL.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "DataIterator.H"

#include "NamespaceHeader.H"

namespace GKUtils {

   inline
   int copyFromLevelData( Real* a_dst, const LevelData<FArrayBox>& a_src )
   {
      int offset(0);
      const DisjointBoxLayout& dbl( a_src.disjointBoxLayout() );
      for (DataIterator dit( a_src.dataIterator() ); dit.ok(); ++dit) {
         FArrayBox tmp( dbl[dit], a_src.nComp(), ( a_dst + offset ) );
         tmp.copy( a_src[dit] );
         offset += dbl[dit].numPts() * a_src.nComp();
      }
      return offset;
   }
   
   inline
   int copyToLevelData( LevelData<FArrayBox>& a_dst, Real* a_src )
   {
      int offset(0);
      const DisjointBoxLayout& dbl( a_dst.disjointBoxLayout() );
      for (DataIterator dit( a_dst.dataIterator() ); dit.ok(); ++dit) {
         FArrayBox tmp( dbl[dit], a_dst.nComp(), ( a_src + offset ) );
         a_dst[dit].copy( tmp );
         offset += dbl[dit].numPts() * a_dst.nComp();
      }
      return offset;
   }
   
   inline
   int addToLevelData( LevelData<FArrayBox>& a_dst,
                       Real* a_src,
                       const Real& a_factor )
   {
      int offset(0);
      const DisjointBoxLayout& dbl( a_dst.disjointBoxLayout() );
      for (DataIterator dit( a_dst.dataIterator() ); dit.ok(); ++dit) {
         FArrayBox tmp( dbl[dit], a_dst.nComp(), ( a_src + offset ) );
         a_dst[dit].plus( tmp, a_factor );
         offset += dbl[dit].numPts() * a_dst.nComp();
      }
      return offset;
   }
   
   inline
   void scaleLevelData( LevelData<FArrayBox>& a_dst, const Real& a_factor )
   {
      for (DataIterator dit( a_dst.dataIterator() ); dit.ok(); ++dit) {
         a_dst[dit] *= a_factor;
      }
   }
   
   inline
   void incrementLevelData( LevelData<FArrayBox>& a_u,
                            const LevelData<FArrayBox>& a_du,
                            const Real& a_factor )
   {
      const DisjointBoxLayout& dbl( a_u.disjointBoxLayout() );
      for (DataIterator dit( a_u.dataIterator() ); dit.ok(); ++dit) {
         FArrayBox tmp( dbl[dit], a_du.nComp() );
         tmp *= a_factor;
         a_u[dit] += tmp;
      }
   }
   
   inline
   double innerProductLevelData( const LevelData<FArrayBox>& a_vec_a,
                                 const LevelData<FArrayBox>& a_vec_b )
   {
      double sum(0.0);
      const DisjointBoxLayout& dbl( a_vec_a.disjointBoxLayout() );
      for (DataIterator dit( a_vec_a.dataIterator() ); dit.ok(); ++dit) {
         sum += a_vec_a[dit].dotProduct( a_vec_b[dit], dbl[dit] );
      }
      return sum;
   }
   
   inline
   double accumulate( const FArrayBox& a_data,
                      const Box& a_box,
                      const int a_p,
                      const int a_comp,
                      const int a_n_comps )
   {
      if (a_p==0) {
         Real max_val(0);
         for (int n(0); n<a_n_comps; n++) {
            max_val = std::max( max_val, a_data.max( a_box, n ) );
         }
         return max_val;
      }
      else if (a_p==1) {
         FArrayBox tmp( a_box, int(a_n_comps) );
         tmp.copy( a_data, a_box );
         tmp.abs();
         return tmp.sum( int(a_comp), int(a_n_comps) );
      }
      return a_data.sumPow( a_box, int(a_p), int(a_comp), int(a_n_comps) );
   }
   
   inline
   Real reduce( double& a_accum, const int a_p )
   {
      double result(0.0);

#ifdef CH_MPI
      MPI_Op type( MPI_SUM );
      if (a_p==0) {
         type = MPI_MAX;
      }
      MPI_Allreduce( &a_accum, &result, 1, MPI_DOUBLE, type, MPI_COMM_WORLD );
#else
      result = a_accum;
#endif
      
      if (a_p>1) { // Lp, p>1
         double exponent( 1.0 / a_p );
         return std::pow( result, exponent );
      }
      return result; 
   }

   inline
   double accumulateLevelData( const LevelData<FArrayBox>& a_data,
                               const int a_p )
   {
      double accum(0.0);
      const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
      const int n( a_data.nComp() );
      for (DataIterator dit( a_data.dataIterator() ); dit.ok(); ++dit) {
         accum += GKUtils::accumulate( a_data[dit], grids[dit], a_p, 0, n );
      }
      return accum;
   }
   
}
   
 
#include "NamespaceFooter.H"


