#include "MagCoordSys.H"

#include "NamespaceHeader.H"


void
MagCoordSys::defineCoordSystemsAndBoundaries( Vector<MagBlockCoordSys *>& a_geom_vec )
{
   m_mappingBlocks.resize(a_geom_vec.size());
   m_coordSysVect.resize(numBlocks(), NULL);

   for (int n=0; n<numBlocks(); ++n) {
      m_coordSysVect[n] = (NewCoordSys*)a_geom_vec[n];
      m_mappingBlocks[n] = a_geom_vec[n]->domain().domainBox();
   }
   m_gotCoordSysVect = m_gotMappingBlocks = true;

   defineBoundaries();
   m_gotBoundaries = true;

   initializeBlockTransformations();
}


RealVect
MagCoordSys::blockRemappingGeneral( int              a_nDst,
                                    const RealVect&  a_xiSrc,
                                    int              a_nSrc) const
{
   return MultiBlockCoordSys::blockRemappingGeneral(a_nDst, a_xiSrc, a_nSrc);
}


bool
MagCoordSys::containsPhysicalBoundary( int                    a_block_number,
                                       int                    a_dir,
                                       const Side::LoHiSide&  a_side ) const
{
   const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = m_boundaries[a_block_number];

   bool contains_boundary = false;

   if ( a_side == Side::LoHiSide::Lo ) {
      contains_boundary = this_block_boundaries[a_dir].isDomainBoundary();
   }
   else if ( a_side == Side::LoHiSide::Hi ) {
      contains_boundary = this_block_boundaries[a_dir+SpaceDim].isDomainBoundary();
   }

   return contains_boundary;
}
   

#include "NamespaceFooter.H"
