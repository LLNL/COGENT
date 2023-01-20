#include "XPointCoordSys.H"
#include "XPointBlockCoordSys.H"

#include "NamespaceHeader.H"


XPointCoordSys::XPointCoordSys()
{
   ParmParse pp("xpsystem");

   Vector<Box> domain_boxes;
   createBlocks(domain_boxes, m_decomp);

   int num_blocks = domain_boxes.size();

   Vector<MagBlockCoordSys *> coord_vec;
   for ( int block_number = 0; block_number < num_blocks; ++block_number ) {
      MagBlockCoordSys* geom
         = new XPointBlockCoordSys(pp, ProblemDomain(domain_boxes[block_number]), block_number );
      coord_vec.push_back(geom);
   }

   defineCoordSystemsAndBoundaries(coord_vec);
}



void
XPointCoordSys::createBlocks( Vector<Box>& a_blocks,
                              IntVect&     a_block_decomp ) const
{
   ParmParse xppp( "xpsystem" );

   Vector<int> block_size(SpaceDim);
   if (xppp.contains("block_size")) {
      xppp.getarr( "block_size", block_size, 0, SpaceDim );
   }
   else {
      MayDay::Error("XPointCoordSys::createBlocks(): xpsystem.block_size not found in input file");
   }

   for (int block=0; block<8; ++block) {

      IntVect lo_mapped_index, hi_mapped_index;
      lo_mapped_index[0] = 2 * block * block_size[0];
      lo_mapped_index[1] = 0;
      hi_mapped_index[0] = lo_mapped_index[0] + block_size[0] - 1;
      hi_mapped_index[1] = block_size[1] - 1;

      a_blocks.push_back(Box(lo_mapped_index, hi_mapped_index));
   }
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
XPointCoordSys::blockRemapping( RealVect&       a_xi_valid,
                                int&            a_n_valid,
                                const RealVect& a_xiSrc,
                                int             a_nSrc ) const
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
