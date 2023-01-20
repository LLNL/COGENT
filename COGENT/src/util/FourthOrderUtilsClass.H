/**
 * Class for holding fourth order stencil operations, both original Chombo
 * versions and new sparse grid-friendly versions
*/

#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

/**
 * A class for holding lots of operators that move between point values, face averages, and cell averages
 * at fourth-order accuracy.  The only data member is m_useSG, which decides whether to use the (more expensive)
 * tensor product stencils that have an error expansion that is compatible with the sparse grid combination technique,
 * or to use the cheaper stencils implemented in Chombo.  
 *
 * Naming conventions: This class is FourthOrder... with a capital F to distinguish it from the numerous Chombo
 * utility functions that start with fourthOrder.  
 *
 * fourthOrder{foo}Gen: A wrapper that calls the Chombo function fourthOrder{foo} if !m_useSG, and calls fourthOrder{foo}SG 
 *    otherwise
 *
 * fourthOrder{foo}SG: Performs the analogous operation to the Chombo function fourthOrder{foo}, but using the sparse grid-friendly 
 *    stencils.  
 *
 * The hope is that with this naming convention, anywhere in COGENT a call to fourthOrder{foo} is made, it can just be replaced with
 * a call to fourthOrder{foo}Gen.
 *
 * The plan is for a class that uses FourthOrderUtil to have its own version of m_useSG that gets set by a parm parse of the input
 * parameter
 *
 * sparse_grids.useSGstencils
 *
 */

class FourthOrderUtil
{
public:

   /// Constructor with initialization
   /**
    * @param[in] pp the input database, should include the useSG boolean
    */
   //FourthOrderUtil( ParmParse& pp);

   /// Destructor.
   /**
    */
   virtual ~FourthOrderUtil() {};

   void setSG( bool SGval ) {useSG = SGval;}

   /****************************************************************************************/

   // These next three map (cell/face)-centered data
   // to (cell/face)-averaged data or vice versa depending 
   // on whether a_sgn = +1 or -1
   void fourthOrderAverageGen( LevelData<FluxBox>& a_phi, 
                               int a_sgn = 1);

   void fourthOrderAverageGen( LevelData<FArrayBox>& a_phi, 
                               int a_sgn = 1);

   void fourthOrderAverageGen( LevelData<EdgeDataBox>& a_phi, 
                               int a_sgn = 1);

   //void fourthOrderAverageGen( LevelData<FArrayBox>& a_phi, 
   //                         const ProblemDomain& a_domain, 
   //                         int a_sgn = 1);
   
   void fourthOrderAverageSG(LevelData<FArrayBox>& a_phi,
							 int a_sgn = 1);

   void fourthOrderAverageSG(LevelData<FluxBox>& a_phi,
							 int a_sgn = 1);
   
   void fourthOrderAverageSG(LevelData<EdgeDataBox>& a_phi,
							 int a_sgn = 1);
   
	// Cell-(averaged/centered) to face-(averaged/centered)
   void fourthOrderCellToFace( LevelData<FluxBox>& a_faceData, 
                                const LevelData<FArrayBox>& a_cellData);

   void fourthOrderCellToFace( FluxBox& a_faceData,
							   const FArrayBox& a_cellData );

   // This maps cell-centered values to face-centered values
   void fourthOrderCellToFaceCenters( LevelData<FluxBox>& a_faceData,
                                      const LevelData<FArrayBox>& a_cellData);

   void fourthOrderCellToFaceCenters( FluxBox& a_faceData, 
									  const FArrayBox& a_cellData );

   // Take face averages of u and v, compute face averages of u*v
   void fourthOrderMult( LevelData<FluxBox>& a_uTimesV,
                         const LevelData<FluxBox>& a_u,
                         const LevelData<FluxBox>& a_v,
                         int a_sgn = 1);

   // Take cell averages of u and v, compute cell average of u*v
   void fourthOrderCellProd( LevelData<FArrayBox>& a_uTimesV,
                             const LevelData<FArrayBox>& a_u,
                             const LevelData<FArrayBox>& a_v);

   /****************************************************************************************/

   void fourthOrderAverageCellGen( FArrayBox& a_phi, 
								   const ProblemDomain& a_domain, 
								   const Box& 			a_bx, 
								   int a_sgn = 1);
     

   // Used by fourthOrderAverage when applied to face-(centered/averaged) LevelData<FluxBox>
   void fourthOrderAverageFaceSG(FluxBox& a_phi, 
                               int a_sgn = 1);

   void fourthOrderAverageFaceSG(FArrayBox& a_phi, 
								 int a_dir, 
								 int a_sgn = 1);

   void fourthOrderAverageEdgeSG(EdgeDataBox& a_phi, 
                               int a_sgn = 1);

// Used by fourthOrderAverage when applied to cell-(centered/averaged) LevelData<FArrayBox>
   void fourthOrderAverageCellSG(FArrayBox& a_phi,
                               int a_sgn = 1);

   // Same as above, but does one-sided cell-(averaging/deaveraging) at domain boundaries
   void fourthOrderAverageCellSG(FArrayBox& a_phi,
								 const ProblemDomain& a_domain,
								 const Box& 		  a_bx,
								 int a_sgn = 1);

   // Used by fourthOrderAverageCellToFace
   void fourthOrderCellToFaceSG(FluxBox& a_faceData, 
                              const FArrayBox& a_cellData);

   // Compute <uv> given <u> and <v> to fourth order.  
   void fourthOrderMultSG( LevelData<FluxBox>& a_uTimesV, 
                           const LevelData<FluxBox>& a_u, 
                           const LevelData<FluxBox>& a_v,
                           int a_sgn = 1);

   // Called by the above
   void fourthOrderMultFaceSG( FluxBox& a_uTimesV, 
                               const FluxBox& a_u, 
                               const FluxBox& a_v, 
                               int a_sgn = 1);
   
   // Take cell averages of u and v, compute cell average of u*v
   void fourthOrderCellProdSG( LevelData<FArrayBox>& a_uTimesV,
                             const LevelData<FArrayBox>& a_u,
                             const LevelData<FArrayBox>& a_v);
   
   // Take cell averages of u and v, compute cell average of u*v
   void fourthOrderCellProdSG( FArrayBox& a_uTimesV,
                             const FArrayBox& a_u,
                             const FArrayBox& a_v);

protected:
   // Boolean that decides whether to use the SG-friendly versions of the various stencils
   bool useSG;

};

#include "NamespaceFooter.H"
