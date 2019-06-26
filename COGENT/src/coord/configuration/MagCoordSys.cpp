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


#include "NamespaceFooter.H"
