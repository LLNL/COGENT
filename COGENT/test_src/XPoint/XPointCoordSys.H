#ifndef XPOINTCOORDSYS_H
#define XPOINTCOORDSYS_H

#define XPOINT_BLOCK_SEPARATION 12

#include "MagBlockCoordSys.H"
#include "MagCoordSys.H"
#include "NamespaceHeader.H"

class XPointCoordSys : public MagCoordSys
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.
  */
   //  XPointCoordSys();
  XPointCoordSys();

  ///
  /** Destructor.
   */
  virtual ~XPointCoordSys();

  /*@}*/

   virtual void define( Vector<MagBlockCoordSys *>& geom_vec);

  /// block mapping conversion function
  /** given xi_k in block k, returns which block (k_valid)
      contains that location in its valid region and its
      location in that block's computational coordinate space.
  */
   virtual void blockRemapping(RealVect& xi_valid,
                               int& n_valid,
                               const RealVect& xiSrc,
                               int nSrc) const;

   virtual void defineBoundaries();

   virtual IntVect getDecomposition(int block) const {return m_decomp;}

protected:

   IntVect m_decomp;

   string m_geom_type;
};



#include "NamespaceFooter.H"

#endif











