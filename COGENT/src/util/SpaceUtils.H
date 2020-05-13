#ifndef _SPACEUTILS_H_
#define _SPACEUTILS_H_

#include "Dimensions.H"
#include "BoundaryBoxLayout.H.multidim"

#include "FourthOrderUtil.H"
#include "Box.H"
#include "FArrayBox.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"

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
  
  void interpEdgesToEdges( LevelData<EdgeDataBox>&,
                     const LevelData<EdgeDataBox>&,
                     const std::string& );

  void interpNodesToEdges( LevelData<EdgeDataBox>&,
                     const LevelData<NodeFArrayBox>&,
                     const std::string& );

  void interpNodesToCells( LevelData<FArrayBox>&,
                     const LevelData<NodeFArrayBox>&,
                     const std::string& );

  void interpCellsToNodes( LevelData<NodeFArrayBox>&,
                     const LevelData<FArrayBox>&,
                     const std::string& );

  void ParaEdgeGradientAtNodes( NodeFArrayBox&,
                          const EdgeDataBox&,
                          const RealVect&,
                          const Box&,
                          const std::string& );

  void PerpEdgeGradientAtCells( FArrayBox&,
                          const EdgeDataBox&,
                          const RealVect&,
                          const Box&,
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
                                   const int);
  
  void extrapBoundaryGhostsForFC(  FArrayBox&, 
                                   const Box&, 
                                   const Box&, 
                                   const int, 
                                   const int, 
                                   const int);

  void secondOrderTransExtrapAtDomainBdry(  FArrayBox&,
                                            const int,
                                            const Box&,
                                            const ProblemDomain&,
                                            const int );

  void fillGhostCellsSimple(  FArrayBox&, const Box&, const int);
  
  void copyAndFillGhostCellsSimple( LevelData<FArrayBox>&, 
                                    const LevelData<FArrayBox>& );
  
  void inspectFArrayBox( const LevelData<FArrayBox>&,
                         const int,
                         const int                   );
  
  void inspectFluxBox( const LevelData<FluxBox>&,
                       const int                 );
  
  void inspectEdgeDataBox( const LevelData<EdgeDataBox>&,
                           const int                     );

  void exchangeEdgeDataBox( LevelData<EdgeDataBox>& );

  void copyEdgeDataBox( EdgeDataBox&,
                  const EdgeDataBox& );

  void copy( FArrayBox&,
	     const FArrayBox& );

  void copy( FArrayBox&,
       const FArrayBox&,
       const Box&        );

  void copyNodeToCell( FArrayBox&,
                 const FArrayBox&,
                 const Box&,
                 const int );

  void localVectorNorm( FArrayBox&,
                  const FArrayBox& );

}

#include "NamespaceFooter.H"
#endif

