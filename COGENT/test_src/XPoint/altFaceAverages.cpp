#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "altFaceAverages.H"
#include "altFaceAveragesF_F.H"

#include <algorithm>

#include "NamespaceHeader.H"
using std::max;

void
altFaceAverages( LevelData<FluxBox>&           a_facePhi,
                 const LevelData<FArrayBox>&   a_cellPhi,
                 const LevelData<FluxBox>&     a_faceVel,
                 const MultiBlockCoordSys&     a_coords)
{
#if 0
   const int nComp( a_facePhi.nComp() );

   // this specifies the number of ghost faces we will need in
   // each transverse and normal direction
   const int transverseGrow(2);
   const int normalGrow(0);

   // may need to do exchange on facePhi
   a_facePhi.exchange();

   const DisjointBoxLayout& grids( a_facePhi.getBoxes() );
   DataIterator dit( grids.dataIterator() );

   // in order to do upwinding, need normal velocities in
   // computational space. computing Fourth-order face
   // averages just in order to upwinding is over-working the
   // issue, but it's the most convenient thing to do here.
   LevelData<FluxBox> normalVel( grids, 1, a_faceVel.ghostVect() );

   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = a_coords.whichBlock(grids[dit]);
      const PhaseBlockCoordSys* block_coord_sys
         = static_cast<const PhaseBlockCoordSys*>(a_coords.getCoordSys(block_number));
      RealVect face_area = block_coord_sys->getMappedFaceArea();

      const FluxBox& thisFaceVel( a_faceVel[dit] );
      FluxBox& thisNormalVel( normalVel[dit] );
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& thisNormalVel_dir( thisNormalVel[dir] );
         thisNormalVel_dir.copy(thisFaceVel[dir],dir,0,1);
         thisNormalVel_dir *= face_area[dir];
      }
   }

   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = a_coords.whichBlock(grids[dit]);
      const PhaseBlockCoordSys* block_coord_sys
         = static_cast<const PhaseBlockCoordSys*>(a_coords.getCoordSys(block_number));
      const RealVect& dx = block_coord_sys->dx();

      const Box& gridBox( grids[dit] );
      Box lapBox( gridBox );

      // need two laplacian ghost cells for this limiter
      int lapBoxGrow( max(normalGrow+2, transverseGrow) );
      lapBox.grow( lapBoxGrow );
      FArrayBox ccLaplacians( lapBox, nComp );
      FluxBox& thisFacePhi( a_facePhi[dit] );
      const FArrayBox& thisCellPhi( a_cellPhi[dit] );
      const FluxBox& thisNormalVel( normalVel[dit] );

      // check to be sure that we have enough phi to compute
      // all of these Laplacians
      {
         Box LapPhiBox( lapBox );
         LapPhiBox.grow(1);
         CH_assert(thisCellPhi.box().contains(LapPhiBox));
      }

      for (int dir=0; dir<SpaceDim; dir++) {

         // box of valid edges for this grid
         Box faceBox( gridBox );

         // for 4th order, we need extra faces in the transverse
         // directions. (for mapped-grid case, we actually need
         // _2_ transverse faces

         // need extra faces in the tangential directions in order
         // to handle 4th-order multiplies
         faceBox.grow(transverseGrow);
         // however, need different number of extra faces in normal dir
         faceBox.grow(dir,normalGrow - transverseGrow);

         faceBox.surroundingNodes(dir);
         FArrayBox& thisFacePhiDir( thisFacePhi[dir] );
         const FArrayBox& thisNormalVelDir( thisNormalVel[dir] );
         {
            // context for computing limited face values
            Box grownFaceBox( faceBox );
            // need an extra face's worth of the FC laplacians
            grownFaceBox.grow( dir,1 );

            FArrayBox centeredLaplacian( grownFaceBox, nComp );

            // now compute limited face value
            FORT_ALTFACEVALUES(CHF_FRA(thisFacePhiDir),
                               CHF_CONST_FRA(thisCellPhi),
                               CHF_CONST_FRA1(thisNormalVelDir,0),
                               CHF_BOX(grownFaceBox),
                               CHF_CONST_INT(dir),
                               CHF_CONST_REAL(dx[dir]) );
            // end context for computing limited face values
         }

      } // end loop over directions

   } // end loop over grids
#endif
}

#include "NamespaceFooter.H"
