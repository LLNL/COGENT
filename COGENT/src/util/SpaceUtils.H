//#ifndef _SPACEUTILS_H_
//#define _SPACEUTILS_H_

#include "BoundaryBoxLayout.H.multidim"

#include "FourthOrderUtil.H"
#include "Box.H"
#include "FArrayBox.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"
#include "CH_HDF5.H"

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

  void interpToFaces(LevelData<FluxBox>&,
                     const LevelData<FArrayBox>&,
                     const int );
  
  void interpToFaces( FluxBox&,
                const FArrayBox&,
                const FArrayBox&,
                const FluxBox&,
                const Box&,
                const std::string& );

  void upWindToFaces( LevelData<FluxBox>&,
                const LevelData<FArrayBox>&,
                const LevelData<FluxBox>&,
                const std::string&,
                const bool& fourth_order = false);
  
  void upWindToFaces( FluxBox&,
                const FArrayBox&,
                const FluxBox&,
                const Box&,
                const std::string&,
                const bool& fourth_order = false);

  void upWindMappedGradient(FArrayBox&,
                            const FArrayBox&,
                            const FArrayBox&,
                            const RealVect&,
                            const std::string&);
  
  void interpToFacesZIP( LevelData<FluxBox>&    a_face_fV,
                   const LevelData<FArrayBox>&  a_cell_f,
                   const LevelData<FArrayBox>&  a_cell_V );
  
  void interpToFacesZIP( LevelData<FluxBox>&    a_face_fuV,
                   const LevelData<FArrayBox>&  a_cell_u,
                   const LevelData<FArrayBox>&  a_cell_fu,
                   const LevelData<FArrayBox>&  a_cell_V,
                   const LevelData<FArrayBox>&  a_cell_fV, 
                   const int                    a_udir );

  void constrainedCrossProductOnEdges( LevelData<EdgeDataBox>&  a_edge_E,
                                 const LevelData<FluxBox>&      a_face_B,
                                 const LevelData<FArrayBox>&    a_cell_V,
                                 const LevelData<FArrayBox>&    a_cell_Cspeed,
                                 const std::string&             a_method );

  void constrainedCrossProductOnNodes( LevelData<NodeFArrayBox>&  a_node_E,
                                 const LevelData<FluxBox>&        a_face_B,
                                 const LevelData<FArrayBox>&      a_cell_V,
                                 const LevelData<FArrayBox>&      a_cell_Cspeed,
                                 const std::string&               a_method );

  void interpToEdges( FArrayBox&    a_Vphi,
                const FArrayBox&    a_phi,
                const FArrayBox&    a_V,
                const Box&          a_box,
                const int           a_dir,
                const std::string&  a_method );
  
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
  
  void interpFacesToFaces( LevelData<FluxBox>&,
                     const LevelData<FluxBox>& );

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

  void interpFaceVectorToCell( LevelData<FArrayBox>&,
                         const LevelData<FluxBox>&,
                         const std::string& );
  
  void interpFaceScalarToCell( LevelData<FArrayBox>&,
                         const LevelData<FluxBox>&,
                         const std::string& );
  
  void interpFaceScalarToCell( LevelData<FArrayBox>&,
                         const LevelData<FluxBox>&,
                         const int a_dir );
  
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
                                      FArrayBox&,
                                      const int a_var_comp = 0 );

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
  
  void copyAndFillGhostCellsSimple( LevelData<FluxBox>&, 
                                    const LevelData<FluxBox>& );
  
  void copyAndFillGhostCellsSimple( LevelData<EdgeDataBox>&, 
                                    const LevelData<EdgeDataBox>& );
  
  void copyAndFillGhostCellsSimple( LevelData<NodeFArrayBox>&, 
                                    const LevelData<NodeFArrayBox>& );
  
  void inspectFArrayBox( const LevelData<FArrayBox>&,
                         const int,
                         const int                   );
  
  void inspectFluxBox( const LevelData<FluxBox>&,
                       const int,
                       const int                 );
  
  void inspectEdgeDataBox( const LevelData<EdgeDataBox>&,
                           const int,
                           const int                 );

  void exchangeFluxBox( LevelData<FluxBox>& );
  
  void exchangeEdgeDataBox( LevelData<EdgeDataBox>& );
  
  void exchangeNodeFArrayBox( LevelData<NodeFArrayBox>& );

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

  void incrementLevelData(  LevelData<FArrayBox>&,
                            const LevelData<FArrayBox>&,
                            const Real a_a = 1);

  void incrementLevelData(  LevelData<FluxBox>&,
                            const LevelData<FluxBox>&,
                            const Real a_a = 1);

  void incrementLevelData(  LevelData<EdgeDataBox>&,
                            const LevelData<EdgeDataBox>&,
                            const Real a_a = 1);

   void findMax(const LevelData<FArrayBox>& a,
                double&                     max_val,
                IntVect&                    max_iv );

   void dotProduct(LevelData<FArrayBox>& dotProduct,
                   const LevelData<FArrayBox>& data_1,
                   const LevelData<FArrayBox>& data_2);

   Real dotProduct( const LevelData<FArrayBox>&  data_1,
                    const LevelData<FArrayBox>&  data_2 );

   Real  MaxNorm( const LevelData<FArrayBox>& a );
   
   Real  MaxNorm( const LevelData<FluxBox>& a );
   
   void applyHarmonicFiltering(LevelData<FArrayBox>& phi,
                               const int& dir);

   void enforcePositivity(LevelData<FArrayBox>& phi);

   // Returns the number of bytes occupied by the data on the calling processor
   unsigned long long int getLevelDataLocalSize( const LevelData<FArrayBox>& data );

   // Returns the number of bytes occupied by the data on the calling processor
   unsigned long long int getLevelDataLocalSize( const LevelData<FluxBox>& data );

   // Returns the number of bytes occupied by the data on the calling processor
   unsigned long long int getLevelDataLocalSize( const LevelData<EdgeDataBox>& data );

   // Returns the number of bytes occupied by the data on the calling processor
   unsigned long long int getLevelDataLocalSize( const LevelData<NodeFArrayBox>& data );

   // Returns the number of bytes occupied by the data on the calling processor
   unsigned long long int getBoxLayoutDataLocalSize( const BoxLayoutData<FArrayBox>& data );

   // Returns the number of bytes occupied by the data on the calling processor
   unsigned long long int getBoxLayoutLocalSize( const BoxLayout& bl );
   
   // Returns the maximum number of bytes occupied by the data per processor
   void getBoxLayoutDataMaxSize(const std::string&               name,
                                const BoxLayoutData<FArrayBox>&  data );

   // Returns the maximum number of bytes occupied by the data per processor
   void getLevelDataMaxSize(const std::string&           name,
                            const LevelData<FArrayBox>&  data );

   // Returns the maximum number of bytes occupied by the data per processor
   void getFArrayBoxMaxSize(const std::string& name,
                            const FArrayBox&   data );

   template<typename T>
   void writeHDF5GroupData(const LevelData<T>&  data,
                           const string&        group,
                           HDF5Handle&          handle)
   {
      handle.setGroup( group );

      HDF5HeaderData header;
      
      header.m_int["is_defined"] = 0;

      if (data.isDefined()) {
         header.m_int["is_defined"] = 1;
         write(handle,data,"data",data.ghostVect());
      }

      header.writeToFile( handle );
   }

   template<typename T>
   void readHDF5GroupData(LevelData<T>&            data,
                          const string&            group,
                          const DisjointBoxLayout& dbl,
                          HDF5Handle&              handle)
   {
   
      handle.setGroup( group );
      
      HDF5HeaderData header;
      header.readFromFile( handle );
   
      int is_data_defined = header.m_int["is_defined"];

      if (is_data_defined == 1) {
         read( handle, data, "data", dbl );
      }
   }

}

#include "NamespaceFooter.H"
//#endif

