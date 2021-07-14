//#ifndef  _CODIM_BC_H_
//#define  _CODIM_BC_H_

#include "Box.H"
#include "IntVect.H"
#include "LoHiSide.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "Vector.H"
#include "MultiBlockCoordSys.H"

#include "BCUtils.H.multidim"

#include "NamespaceHeader.H"

namespace CodimBC {
   void setCodimBoundaryValues(
      LevelData<FArrayBox>&     a_soln,
      const MultiBlockCoordSys& a_coord_sys );
   
   void setCodimCornerValues(
      LevelData<FArrayBox>&     a_soln,
      const MultiBlockCoordSys& a_coord_sys );

   void fillCodimGhostCells( FArrayBox&                    a_this_soln,
                             const Box&                    a_boundary_box,
                             const Vector<int>&            a_dirs,
                             const Vector<Side::LoHiSide>& a_sides,
                             const int&                    a_codim );

   Box getGhostBox( const int a_codim,
                    const int& a_current_dir,
                    const Box& a_box,
                    const Box& a_domain_box,
                    const Vector<int>& a_dir,
                    const Vector<Side::LoHiSide>& a_side,
                    const IntVect& a_ghosts );


   inline
   bool isPhysicalBoundary( const MultiBlockCoordSys&     a_coord_sys,
                            const Box&                    a_box,
                            const Vector<int>&            a_dirs,
                            const Vector<Side::LoHiSide>& a_sides )
   {
      bool is_physical_boundary(false);
      for (int n(0); n<a_dirs.size(); n++) {
         is_physical_boundary |= BCUtils::touchesPhysicalBoundary( a_coord_sys,
                                                                   a_box,
                                                                   a_dirs[n],
                                                                   a_sides[n] );
      }
      return is_physical_boundary;
   }
   
   inline
   bool isPhysicalCorner( const MultiBlockCoordSys&     a_coord_sys,
                            const Box&                    a_box,
                            const Vector<int>&            a_dirs,
                            const Vector<Side::LoHiSide>& a_sides,
                            const int                     a_codim )
   {
      bool is_physical_boundary(false);
      int num_physical_boundaries = 0;
      for (int n(0); n<a_dirs.size(); n++) {
         is_physical_boundary = BCUtils::touchesPhysicalBoundary( a_coord_sys,
                                                                  a_box,
                                                                  a_dirs[n],
                                                                  a_sides[n] );
         if(is_physical_boundary) { // need a_codim physical boundaries for corner
            num_physical_boundaries = num_physical_boundaries + 1;
         }
      }
      
      bool is_physical_corner = false;
      if(num_physical_boundaries>=a_codim) is_physical_corner = true;
      return is_physical_corner;
   }

}

#include "NamespaceFooter.H"

//#endif
