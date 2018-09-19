#include "XPointCoordSys.H"
#include "XPointBlockCoordSys.H"
#include "FourthOrderUtil.H"
#include "BlockBoundary.H"
#include "SpreadingCopier.H"
#include "ReductionOps.H"
#include "Directions.H"
#include "CONSTANTS.H"


#include "NamespaceHeader.H"


XPointCoordSys::XPointCoordSys()
{
  ParmParse pp("xpsystem");

  Vector<double> decomp(SpaceDim);
  pp.getarr("block_decomp", decomp, 0, SpaceDim);
  for (int i=0; i<SpaceDim; ++i) {
    m_decomp[i] = decomp[i];
  }
}



XPointCoordSys::~XPointCoordSys()
{
}



void
XPointCoordSys::define( Vector<MagBlockCoordSys *>& a_geom_vec )
{
   // Just define the base class data for now.  We'll probably need to put more here later.
   MagCoordSys::defineCoordSystemsAndBoundaries(a_geom_vec);
}



void
XPointCoordSys::defineBoundaries()
{
   CH_assert(SpaceDim==2);
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   IntVect permDiff = IntVect(D_DECL(1, 0, 2));
   IntVect clockwise = IntVect(D_DECL(-1, 1, 1));
   IntVect anticlockwise = IntVect(D_DECL(1, -1, 1));

   for (int iblock = 0; iblock < numBlocks(); iblock++) {

      IndicesTransformation it;
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[iblock];
      const XPointBlockCoordSys* coord_sys_src = (const XPointBlockCoordSys*)m_coordSysVect[iblock];
      int block_type = coord_sys_src->blockType();
      int neighbor_block[SpaceDim];

      neighbor_block[0] = (iblock+numBlocks()-1) % numBlocks();
      neighbor_block[1] = (iblock+1) % numBlocks();

      for (int dir=0; dir<SpaceDim; ++dir) {
         const XPointBlockCoordSys* coord_sys_dst = (const XPointBlockCoordSys*)m_coordSysVect[neighbor_block[dir]];
         IntVect rotation;
         if (dir==0) {
            rotation = clockwise;
         }
         else {
            rotation = anticlockwise;
         }
         IntVect pold = coord_sys_src->domainSmallEnd();
         IntVect pnew = coord_sys_dst->domainSmallEnd();
         it.defineFromPivot(pold, pnew, permDiff, rotation);
         blockBoundaries[dir].define(it,neighbor_block[dir]);
         blockBoundaries[dir + SpaceDim].define(0);
      }
   }
}



void
XPointCoordSys::blockRemapping(RealVect&       a_xi_valid,
                               int&            a_n_valid,
                               const RealVect& a_xiSrc,
                               int             a_nSrc) const
{
   const XPointBlockCoordSys* coord_sys_Src = (const XPointBlockCoordSys*)m_coordSysVect[a_nSrc];
   RealVect X = coord_sys_Src->realCoord(a_xiSrc);

   double R = X[0] - coord_sys_Src->centralMajorRadius();
   double Z = X[1];

   if (R >= 0. && Z > 0.) {     // upper right quadrant
      if (Z > R) {
         a_n_valid = 0;
      }
      else {
         a_n_valid = 1;
      }
   }
   else if (R > 0. && Z <= 0.) { // lower right quadrant
      if (Z > -R) {
         a_n_valid = 2;
      }
      else {
         a_n_valid = 3;
      }
   }
   else if (R <= 0. && Z < 0.) { // lower left quadrant
      if (Z < R) {
         a_n_valid = 4;
      }
      else {
         a_n_valid = 5;
      }
   }
   else if (R < 0. && Z >= 0.) { // upper left quadrant
      if (Z < -R) {
         a_n_valid = 6;
      }
      else {
         a_n_valid = 7;
      }
   }

   const XPointBlockCoordSys* coord_sys_valid = (const XPointBlockCoordSys*)m_coordSysVect[a_n_valid];
   a_xi_valid = coord_sys_valid->mappedCoord(X);

   // Perform a final check that a_xi_valid is contained in the valid region of the block whose
   // index is being returned
   if ( coord_sys_valid->lowerMappedCoordinate(0) > a_xi_valid[0] ||
        coord_sys_valid->upperMappedCoordinate(0) < a_xi_valid[0] ||
        coord_sys_valid->lowerMappedCoordinate(1) > a_xi_valid[1] ||
        coord_sys_valid->upperMappedCoordinate(1) < a_xi_valid[1] ) {
      a_n_valid = -1;
      a_xi_valid = a_xiSrc;
   }
}



#include "NamespaceFooter.H"
