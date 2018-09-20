#include "CodimBC.H"
#include "CodimBCF_F.H"

#include "BoundaryLookupTable.H.multidim"

#include "NamespaceHeader.H"

inline
Box CodimBC::getGhostBox( const int a_codim,
                          const int& a_current_index,
                          const Box& a_box,
                          const Box& a_domain_box,
                          const Vector<int>& a_dir,
                          const Vector<Side::LoHiSide>& a_side,
                          const IntVect& a_ghosts )
{
   const int dir( a_dir[a_current_index] );
   const Side::LoHiSide side( a_side[a_current_index] );
   const int num_ghosts( a_ghosts[dir] );
   const Box gBox( adjCellBox( a_box, dir, side, num_ghosts ) );
   if (a_current_index+1<a_codim) {
      return getGhostBox( a_codim,
                          a_current_index+1,
                          gBox,
                          a_domain_box,
                          a_dir,
                          a_side,
                          a_ghosts );
   }
   Box gDomainBox( a_domain_box );
   gDomainBox.grow( a_ghosts );
   return (gBox & gDomainBox);
}


inline
void CodimBC::fillCodimGhostCells(
   FArrayBox&                    a_this_soln,
   const Box&                    a_boundary_box,
   const Vector<int>&            a_dirs,
   const Vector<Side::LoHiSide>& a_sides,
   const int&                    a_codim )
{
   Vector<int> sides( a_sides.size() );
   for (int i(0); i<sides.size(); i++)
      sides[i] = static_cast<int>(a_sides[i]);
   FORT_FILL_CODIM_GHOST_CELLS( CHF_FRA(a_this_soln),
                                CHF_BOX(a_boundary_box),
                                CHF_CONST_I1D(&(a_dirs[0]),a_dirs.size()),
                                CHF_CONST_I1D(&(sides[0]),sides.size()),
                                CHF_CONST_INT(a_codim) );
}


#if 0
void CodimBC::setCodimBoundaryValues(
   LevelData<FArrayBox>&     a_uJ,
   const MultiBlockCoordSys& a_coord_sys )
{
   const BoundaryLookupTable& boundaries( BoundaryLookupTable::getLookupTable() );
   const DisjointBoxLayout& grids( a_uJ.disjointBoxLayout() );
   for (DataIterator dit( a_uJ.dataIterator() ); dit.ok(); ++dit) {
      for (int codim(2); codim<=SpaceDim; codim++) {
         for (int m(0); m<boundaries.numberOfBoundaryCases(codim); m++) {

            const Vector<int>& dirs( boundaries.getDirections( m, codim ) );
            const Vector<Side::LoHiSide>& sides( boundaries.getSides( m, codim ) );
            const Box& box( grids[dit] );
            if (isPhysicalBoundary( a_coord_sys, box, dirs, sides )) {
               int CURRENT_INDEX(0);
               const Box& domain_box( BCUtils::blockDomain( box, a_coord_sys ) );

               Box boundary_box( CodimBC::getGhostBox( codim,
                                                       CURRENT_INDEX,
                                                       box,
                                                       domain_box,
                                                       dirs,
                                                       sides,
                                                       a_uJ.ghostVect() ) );

               if (!boundary_box.isEmpty()) {
                  fillCodimGhostCells( a_uJ[dit], boundary_box, dirs, sides, codim );
               }
            }
         }
      }
   }
}
#else
void CodimBC::setCodimBoundaryValues(
   LevelData<FArrayBox>&     a_u,
   const MultiBlockCoordSys& a_coord_sys )
{
   const BoundaryLookupTable& boundaries( BoundaryLookupTable::getLookupTable() );
   const DisjointBoxLayout& grids( a_u.disjointBoxLayout() );
   for (DataIterator dit( a_u.dataIterator() ); dit.ok(); ++dit) {
      for (int codim(2); codim<=SpaceDim; codim++) {
         for (int m(0); m<boundaries.numberOfBoundaryCases(codim); m++) {

            const Vector<int>& dirs( boundaries.getDirections( m, codim ) );
            const Vector<Side::LoHiSide>& sides( boundaries.getSides( m, codim ) );
            const Box& box( grids[dit] );
            if (isPhysicalBoundary( a_coord_sys, box, dirs, sides )) {
               int CURRENT_INDEX(0);
               const Box& domain_box( BCUtils::blockDomain( box, a_coord_sys ) );

               Box boundary_box( CodimBC::getGhostBox( codim,
                                                       CURRENT_INDEX,
                                                       box,
                                                       domain_box,
                                                       dirs,
                                                       sides,
                                                       a_u.ghostVect() ) );
               if (!boundary_box.isEmpty()) {
                  fillCodimGhostCells( a_u[dit], boundary_box, dirs, sides, codim );
               }
            }
         }
      }
   }
}
#endif

#include "NamespaceFooter.H"


