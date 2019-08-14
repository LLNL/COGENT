#ifndef _SPACEUTILS_H_
#define _SPACEUTILS_H_

#include "Dimensions.H"
#include "BoundaryBoxLayout.H.multidim"

#include "FourthOrderUtil.H"
#include "Box.H"
#include "FArrayBox.H"
#include "EdgeDataBox.H"

#include "NamespaceHeader.H"

namespace SpaceUtils {

  void computeLaxSplitting( LevelData<FArrayBox>&,
                            LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>& );
  
  void computeLaxSplitting( LevelData<FArrayBox>&,
                            LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>&,
                      const LevelData<FArrayBox>&,
                      const int );
  
  void computeLaxSplitting( FArrayBox&,
                            FArrayBox&,
                      const FArrayBox&,
                      const FArrayBox&,
                      const FArrayBox& );
  
  void interpToFaces( LevelData<FluxBox>&,
                const LevelData<FArrayBox>&,
                const LevelData<FArrayBox>&,
                const LevelData<FluxBox>&,
                const std::string& );
  
  void interpToFaces( FluxBox&,
                const FArrayBox&,
                const FArrayBox&,
                const FluxBox&,
                const Box&,
                const std::string& );

  void upWindToFaces( LevelData<FluxBox>&,
                const LevelData<FArrayBox>&,
                const LevelData<FluxBox>&,
                const std::string& );
  
  void upWindToFaces( FluxBox&,
                const FArrayBox&,
                const FluxBox&,
                const Box&,
                const std::string& );
  
  
  void interpCellToEdges( LevelData<EdgeDataBox>&,
                    const LevelData<FArrayBox>&,
                    const LevelData<EdgeDataBox>&,
                    const std::string& );
  
  void interpEdgesToCell( LevelData<FArrayBox>&,
                    const LevelData<EdgeDataBox>&,
                    const std::string& );
  
  void interpFacesToCell( LevelData<FArrayBox>&,
                    const LevelData<FluxBox>&,
                    const std::string& );
  
  void interpToFacesWENO( LevelData<FluxBox>&,
                    const LevelData<FArrayBox>&,
                    const LevelData<FluxBox>&,
                    const LevelData<FArrayBox>&,
                    const std::string& );
  
  void interpToFacesWENO( FluxBox&,
                    const FArrayBox&,
                    const FluxBox&,
                    const FArrayBox&,
                    const Box&,
                    const std::string& );

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

  void fillGhostCellsSimple(  FArrayBox&, const Box&, const int);
  
  void copyAndFillGhostCellsSimple( LevelData<FArrayBox>&, 
                                    const LevelData<FArrayBox>& );

  void copy( FArrayBox&,
	     const FArrayBox& );

}

#include "NamespaceFooter.H"
#endif

