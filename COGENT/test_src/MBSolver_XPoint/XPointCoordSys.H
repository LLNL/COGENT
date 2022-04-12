#ifndef XPOINTCOORDSYS_H
#define XPOINTCOORDSYS_H

#include "MagCoordSys.H"

#include "NamespaceHeader.H"

class XPointCoordSys
   : public MagCoordSys
{
public:

   XPointCoordSys();

   ~XPointCoordSys() {};

   virtual void blockRemapping( RealVect&       xi_valid,
                                int&            n_valid,
                                const RealVect& xiSrc,
                                int             nSrc ) const;

   virtual void defineBoundaries();

   virtual IntVect getDecomposition(int block) const {return m_decomp;}

   virtual std::string type() const {return "XPoint";}

protected:

   void createBlocks( Vector<Box>& blocks,
                      IntVect&     block_decomp ) const;

   IntVect m_decomp;
};


#include "NamespaceFooter.H"

#endif











