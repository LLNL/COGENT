//#ifndef _SGINTERP_H_
//#define _SGINTERP_H_

#include <vector>
#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "ProblemDomain.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "SpaceUtils.H.multidim"
#include "FourthOrderUtilsClass.H.multidim"

#include "NamespaceHeader.H"

namespace SGInterp{

  const int num_gpt = 4;

  void coarsen1DCC( FArrayBox&,
                    const Box&,
                    const FArrayBox&,
                    const Box&,
                    const int );

  void refine1DCC(  FArrayBox&,
                    const Box&,
                    const FArrayBox&,
                    const Box&,
                    const int );

  void coarsen1DFC( FArrayBox&,
                    const Box&,
                    const FArrayBox&,
                    const Box&,
                    const int );

  void refine1DFC(  FArrayBox&,
                    const Box&,
                    const FArrayBox&,
                    const Box&,
                    const int );

  void interpolate( FArrayBox&,
                    const Box&,
                    const FArrayBox&,
                    const Box& );

  void interpolate( FluxBox&,
                    const Box&,
                    const FluxBox&,
                    const Box& );

  void interpolate( EdgeDataBox&,
                    const Box&,
                    const EdgeDataBox&,
                    const Box& );

  template<typename T>
  void interpolate( LevelData<T>&       a_dst,
                    const LevelData<T>& a_src )
  {
    CH_assert(a_dst.isDefined());
    CH_assert(a_src.isDefined());
    CH_assert( a_src.nComp() == a_dst.nComp() );
  
    const DisjointBoxLayout& dst_dbl( a_dst.disjointBoxLayout() );
    const DisjointBoxLayout& src_dbl( a_src.disjointBoxLayout() );
  
    const ProblemDomain& dst_domain( dst_dbl.physDomain() );
    const ProblemDomain& src_domain( src_dbl.physDomain() );
  
    const IntVect dst_size =  dst_domain.domainBox().bigEnd() 
                              - dst_domain.domainBox().smallEnd() 
                              + IntVect::Unit;
    const IntVect src_size =  src_domain.domainBox().bigEnd() 
                              - src_domain.domainBox().smallEnd() 
                              + IntVect::Unit;
  
    IntVect refine_vec(IntVect::Unit), coarse_vec(IntVect::Unit);
    for (int dir = 0; dir < SpaceDim; dir++) {
      if (dst_size[dir] > src_size[dir]) {
        refine_vec[dir] = dst_size[dir] / src_size[dir];
      } else if (src_size[dir] > dst_size[dir]) {
        coarse_vec[dir] = src_size[dir] / dst_size[dir];
      }
    }
  
    DisjointBoxLayout tmp_dbl, tmp_tmp_dbl;
    refine( tmp_tmp_dbl, src_dbl, refine_vec );
    coarsen( tmp_dbl, tmp_tmp_dbl, coarse_vec );
    tmp_dbl.setDomain( dst_domain );
  
    LevelData<T> src_wg(src_dbl, a_src.nComp(), num_gpt*IntVect::Unit);
    SpaceUtils::copyAndFillGhostCellsSimple( src_wg, a_src );
  
    FourthOrderUtil FourthOrderOperators;
    FourthOrderOperators.setSG(true); 
    FourthOrderOperators.fourthOrderAverageGen( src_wg, -1 ); 

    LevelData<T> dst_tmp( tmp_dbl, a_dst.nComp(), IntVect::Zero );
    for (auto dit(src_dbl.dataIterator()); dit.ok(); ++dit) {
      interpolate(  dst_tmp[dit],
                    tmp_dbl[dit],
                    src_wg[dit],
                    src_dbl[dit] );
    }

    dst_tmp.copyTo(a_dst);

    LevelData<T> dst_wg( dst_dbl, a_dst.nComp(), IntVect::Unit );
    SpaceUtils::copyAndFillGhostCellsSimple( dst_wg, a_dst );
    FourthOrderOperators.fourthOrderAverageGen( dst_wg, 1 ); 

    for (auto dit(dst_dbl.dataIterator()); dit.ok(); ++dit) {
      a_dst[dit].copy( dst_wg[dit] );
    }
  
    return;
  }

  template<typename T>
  void combine( LevelData<T>* const                 a_dst,
                const std::vector< LevelData<T>* >  a_src_vec,
                const std::vector<Real>&            a_coeffs )
  {
    CH_assert( a_coeffs.size() == a_src_vec.size() );
    for (auto dit(a_dst->dataIterator()); dit.ok(); ++dit) {
      (*a_dst)[dit].setVal(0.0);
    }
    for (int s = 0; s < a_coeffs.size(); s++) {
      LevelData<T> tmp (  a_dst->disjointBoxLayout(), 
                          a_dst->nComp(),
                          IntVect::Zero );
      interpolate( tmp, *(a_src_vec[s]) );
      SpaceUtils::incrementLevelData( *a_dst, tmp, a_coeffs[s] );
    }
    return;
  }

}

#include "NamespaceFooter.H"
//#endif
