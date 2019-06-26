#ifndef _UNSTRUCTUREDINTERPUTILS_H_
#define _UNSTRUCTUREDINTERPUTILS_H_

#include <vector>
#include <utility>
#include <map>
#include "Dimensions.H"
#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "PhaseCoordSys.H"

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;

typedef pair<CFG::IntVect,Real> StencilPoint;
typedef std::vector<StencilPoint> Stencil;

namespace UnstructuredInterpUtils {

  void bilinearInterpStencil( Stencil&,
                              const CFG::FArrayBox&,
                              const std::vector<CFG::IntVect>&,
                              const CFG::RealVect& );

  void trilinearInterpStencil(  Stencil&,
                                const CFG::FArrayBox&,
                                const std::vector<CFG::IntVect>&,
                                const CFG::RealVect& );

}


/*! \def _MatrixInvert_
 * B =A^{-1} (Inverse of \a A), where \a A, \a B are square matrices of
 * size \a N stored as 1D arrays in row-major format. Uses LU decomposition.
 * If \a A is singular, \a B will contain NaN or Inf.
*/
#define _MatrixInvert_(A,B,N,status)  \
  { \
    status = 0; \
    int matopsi,matopsj,matopsk; \
    double matopsfactor, matopssum, matopsAc[(N)*(N)]; \
    \
    /* make a copy of A */ \
    for (matopsi=0; matopsi<(N)*(N); matopsi++) *(matopsAc+matopsi) = *((A)+matopsi); \
    \
    /* set B as the identity matrix */  \
    for (matopsi=0; matopsi<(N)*(N); matopsi++) *((B)+matopsi)             = 0.0; \
    for (matopsi=0; matopsi<(N)    ; matopsi++) *((B)+matopsi*(N)+matopsi) = 1.0; \
    \
    /* LU Decomposition - Forward Sweep */ \
    for (matopsi=0; matopsi<(N)-1; matopsi++) { \
      if (*(matopsAc+matopsi*(N)+matopsi) == 0) status = 1; \
      for (matopsj=matopsi+1; matopsj<(N); matopsj++) { \
        matopsfactor = *(matopsAc+matopsj*(N)+matopsi) / *(matopsAc+matopsi*(N)+matopsi); \
        for (matopsk=matopsi+1; matopsk<(N); matopsk++) *(matopsAc+matopsj*(N)+matopsk) -= (matopsfactor * *(matopsAc +matopsi*(N)+matopsk)); \
        for (matopsk=0; matopsk<matopsj; matopsk++) *((B)+matopsj*(N)+matopsk) -= (matopsfactor * *((B)+matopsi*(N)+matopsk)); \
      } \
    } \
    \
    /* LU Decomposition - Backward Sweep */ \
    for (matopsi=(N)-1; matopsi>=0; matopsi--) { \
      for (matopsk=0; matopsk<(N); matopsk++) { \
        matopssum = 0.0; \
        for (matopsj=matopsi+1; matopsj<(N); matopsj++) matopssum += (*(matopsAc+matopsi*(N)+matopsj) * *((B)+matopsj*(N)+matopsk)); \
        *((B)+matopsi*(N)+matopsk) = (*((B)+matopsi*(N)+matopsk) - matopssum) / *(matopsAc+matopsi*(N)+matopsi); \
      } \
    } \
    \
    /* Done - B contains A^{-1} now */ \
  }
#include "NamespaceFooter.H"
#endif

