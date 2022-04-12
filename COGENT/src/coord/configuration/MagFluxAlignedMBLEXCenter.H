#ifndef _MAGFLUXALIGNEDMBLEXCENTER_H_
#define _MAGFLUXALIGNEDMBLEXCENTER_H_

#include "MultiBlockLevelExchangeCenter.H"
#include "SingleNullCoordSys.H"

#include "NamespaceHeader.H"

/*
  This is a specialization of MultiBlockLevelExchangeCenter for use
  on flux-aligned single null geometries.  Rather that performing
  interpolation to fill extrablock ghost cells, the interpGhosts() 
  member will instead perform extrablock ghost exchanges.  The m_ghostCells
  and m_stencils data members are obtained from SingleNullCoordSys::GhostCellsUe()
  and SingleNullCoordSys::getStencilsUe() by overriding the
  MultiBlockLevelExchangeCenter::define() member.

  N.B.: Unless the geometry coordinate system is SingleNullCoordSys with
  extrablock ghost cell exchanges enabled, this class will revert to the
  standard base class behavior.
*/

class MagFluxAlignedMBLEXCenter : public MultiBlockLevelExchangeCenter
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.
  */
   MagFluxAlignedMBLEXCenter()
      : m_do_extrablock_exchanges(false) {};

  ///
  /** Destructor.
   */
   virtual ~MagFluxAlignedMBLEXCenter()
   {
      if ( m_do_extrablock_exchanges ) {
         m_isDefined = false;  // Keeps the base destructor from trying to delete
                               // stuff that was never defined by the over-riding define()
      }
   }

  /*@}*/

  virtual void define(const MultiBlockLevelGeom*  a_geomPtr,
                      int                         a_ghosts,
                      int                         a_order)
   {
      const SingleNullCoordSys* coord_sys_ptr = dynamic_cast<SingleNullCoordSys *>(a_geomPtr->coordSysPtr());

      if (coord_sys_ptr && ((MagGeom*)a_geomPtr)->extrablockExchange() ) {

         m_geomPtr = (MultiBlockLevelGeom*)a_geomPtr;

         const LayoutData< IntVectSet >& ghostCells_ue = coord_sys_ptr->getGhostCellsUe();

         m_ghostCells.define(ghostCells_ue.boxLayout());
         for (DataIterator dit(m_ghostCells.dataIterator()); dit.ok(); ++dit) {
            m_ghostCells[dit].define(ghostCells_ue[dit]);
         }

         const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencil_ue = coord_sys_ptr->getStencilsUe();

         m_stencils.define(stencil_ue.boxLayout());
         for (DataIterator dit(m_stencils.dataIterator()); dit.ok(); ++dit) {
            m_stencils[dit] = stencil_ue[dit];
         }

         m_isDefined = true;
         m_do_extrablock_exchanges = true;
      }
      else {
         MultiBlockLevelExchangeCenter::define(a_geomPtr, a_ghosts, a_order);
      }
   }

  virtual void interpGhosts(LevelData<FArrayBox>& a_data) const
   {
      if ( m_do_extrablock_exchanges ) {
         ((MagGeom*)m_geomPtr)->exchangeExtraBlockGhosts(a_data);
      }
      else {
         MultiBlockLevelExchangeCenter::interpGhosts(a_data);
      }
   }

private:

   bool m_do_extrablock_exchanges;

};


#include "NamespaceFooter.H"

#endif
