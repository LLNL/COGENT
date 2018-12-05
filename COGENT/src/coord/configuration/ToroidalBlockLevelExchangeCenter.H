#ifndef _TOROIDAL_BLOCK_LEVEL_EXCHANGE_CENTER_
#define _TOROIDAL_BLOCK_LEVEL_EXCHANGE_CENTER_

#include "MultiBlockLevelExchangeCenter.H"
#include "MagGeom.H"

#include "NamespaceHeader.H"

class ToroidalBlockLevelExchangeCenter
  : public MultiBlockLevelExchangeCenter
{
public:

   ToroidalBlockLevelExchangeCenter( const MagGeom&  geometry,
                                     int             ghosts,
                                     int             order );

   virtual ~ToroidalBlockLevelExchangeCenter() {m_isDefined = false;}

   void interpGhosts( LevelData<FArrayBox>& data ) const;

protected:

   IntVectSet extraBlockGhosts( const Box&  box,
                                int         ghosts ) const;

   const MagGeom& m_geometry;
};

#include "NamespaceFooter.H"

#endif
