#ifndef _PHASE_TBLEX_
#define _PHASE_TBLEX_

#include "PhaseGeom.H"

#include "NamespaceHeader.H"

class PhaseTBLEX
{
public:

   PhaseTBLEX( const PhaseGeom&  geometry,
               int               ghosts,
               int               order );

   virtual ~PhaseTBLEX() {;}

   void interpGhosts( LevelData<FArrayBox>& data ) const;

   const LayoutData< IntVectSet >& ghostCells() const {return m_ghostCells;}

   const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencils() const {return m_stencils;}

protected:

   IntVectSet extraBlockGhosts( const Box&  box,
                                int         ghosts ) const;

   const PhaseGeom& m_geometry;

   LayoutData< IntVectSet > m_ghostCells;
   LayoutData< RefCountedPtr< IVSFAB<MBStencil> > > m_stencils;
};

#include "NamespaceFooter.H"

#endif
