#ifndef _SPACEUTILS_H_
#define _SPACEUTILS_H_

#include "Dimensions.H"
#include "BoundaryBoxLayout.H.multidim"

#include "FourthOrderUtil.H"
#include "Box.H"
#include "FArrayBox.H"

#include "NamespaceHeader.H"

namespace SpaceUtils {

  void faceInterpolate( const int,
                        const Box&,
                        const int,
                        const FArrayBox&,
                        FArrayBox& );

  void faceInterpolate( const int,
                        const Box&,
                        const Box&,
                        const int,
                        const FArrayBox&,
                        FArrayBox& );

  void cellCenteredGradientComponent( const Box&,
                                      const int,
                                      const FArrayBox&,
                                      const RealVect&,
                                      const int,
                                      FArrayBox& );

  void faceCenteredGradientComponent( const Box&,
                                      const int,
                                      const FArrayBox&,
                                      const RealVect&,
                                      const int,
                                      FArrayBox& );

  void extrapBoundaryGhostsForCC(  FArrayBox&, 
                                   const Box&, 
                                   const Box&, 
                                   const int, 
                                   const int, 
                                   const Tuple<BlockBoundary, 2*SpaceDim>&);
  
  void extrapBoundaryGhostsForFC(  FArrayBox&, 
                                   const Box&, 
                                   const Box&, 
                                   const int, 
                                   const int, 
                                   const Tuple<BlockBoundary, 2*SpaceDim>&);

  void secondOrderTransExtrapAtDomainBdry(  FArrayBox&,
                                            const int,
                                            const Box&,
                                            const ProblemDomain&,
                                            const int );

}

#include "NamespaceFooter.H"
#endif

