#ifndef _PHASE_BC_UTILS_H_
#define _PHASE_BC_UTILS_H_

#include "Box.H"
#include "LoHiSide.H"
#include "BlockBoundary.H"
#include "Vector.H"
#include "MultiBlockCoordSys.H"
#include "Directions.H"
#include "DisjointBoxLayout.H"
#include "IntVect.H"
#include "KineticSpecies.H"

#include "BCUtils.H.multidim"
#include "BoundaryBoxLayout.H.multidim"
#include "FourthOrderBC.H.multidim"

#include "NamespaceHeader.H"

namespace PhaseBCUtils {

   inline
   void defineBoundaryBoxLayouts( BoundaryBoxLayoutPtrVect& a_bdry_layouts,
                                  const DisjointBoxLayout& a_grids,
                                  const MultiBlockCoordSys& a_coord_sys,
                                  const IntVect& a_ghost_vect )
   {
      const Vector<Box>& blocks( a_coord_sys.mappingBlocks() );
      for (int b(0); b<blocks.size(); b++) {
         for (int dir(RADIAL_DIR); dir<=MU_DIR; dir++) {
            for (SideIterator si; si.ok(); ++si) {
               Side::LoHiSide side( si() );
               if (BCUtils::isPhysicalBoundary( a_coord_sys, blocks[b], dir, side )) {
                  a_bdry_layouts.push_back(
                     BoundaryBoxLayoutPtr( new BoundaryBoxLayout( a_grids,
                                                                  a_coord_sys,
                                                                  blocks[b],
                                                                  dir,
                                                                  side,
                                                                  a_ghost_vect )));
               }
            }
         }
      } 
   }
 
   inline
   void defineInflowDataStorage( KineticSpeciesPtrVect& a_bdry_data,
                                 const BoundaryBoxLayoutPtrVect& a_bdry_layout,
                                 const KineticSpecies& a_prototype_species )

   {
      const LevelData<FArrayBox>& prototype_dfn( a_prototype_species.distributionFunction() );
      a_bdry_data.resize( a_bdry_layout.size() );
      for (int i(0); i<a_bdry_layout.size(); i++) {
         const DisjointBoxLayout& dbl( a_bdry_layout[i]->disjointBoxLayout() );
         a_bdry_data[i] = a_prototype_species.clone( IntVect::Zero, false );
         LevelData<FArrayBox>& dfn( a_bdry_data[i]->distributionFunction() );
         dfn.define( dbl, prototype_dfn.nComp(), IntVect::Zero );
      }
   }


   inline
   void setInflowOutflowBC( LevelData<FArrayBox>& a_BfJ,
                            const BoundaryBoxLayoutPtrVect& a_all_bdry_layouts,
                            const KineticSpeciesPtrVect& a_all_bdry_data,
                            const MultiBlockCoordSys& a_coord_sys,
                            const LevelData<FluxBox>& a_velocity )
   {
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            const Box fill_box( bdry_grids[dit] );
            
            FArrayBox& this_BfJ( a_BfJ[bdry_layout.dataIndex(dit)] );
            const KineticSpecies& bdry_species( *(a_all_bdry_data[b]) );
            const LevelData<FArrayBox>& bdry_data( bdry_species.distributionFunction() );
            const FArrayBox& this_inflow_data( bdry_data[dit] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            const FluxBox& this_face_vel( a_velocity[interior_dit] );
            FourthOrderBC::setInflowOutflowBC( this_BfJ,
                                               fill_box,
                                               this_inflow_data,
                                               this_face_vel,
                                               bdry_layout.dir(),
                                               bdry_layout.side() );
         }
      }
   }
   
}
#include "NamespaceFooter.H"

#endif
