#ifndef _FLUID_BC_UTILS_H_
#define _FLUID_BC_UTILS_H_

#include "Box.H"
#include "LoHiSide.H"
#include "BlockBoundary.H"
#include "Vector.H"
#include "MultiBlockCoordSys.H"
#include "Directions.H"
#include "DisjointBoxLayout.H"
#include "IntVect.H"
#include "FluidSpecies.H"
#include "EdgeDataBox.H"

#include "BCUtils.H.multidim"
#include "BoundaryBoxLayout.H.multidim"
#include "FourthOrderBC.H.multidim"

#include "NamespaceHeader.H"

namespace FluidBCUtils {

   inline
   void defineBoundaryBoxLayouts( BoundaryBoxLayoutPtrVect& a_bdry_layouts,
                                  const DisjointBoxLayout& a_grids,
                                  const MultiBlockCoordSys& a_coord_sys,
                                  const IntVect& a_ghost_vect )
   {
      const Vector<Box>& blocks( a_coord_sys.mappingBlocks() );
      for (int b(0); b<blocks.size(); b++) {
         for (int dir(RADIAL_DIR); dir<SpaceDim; dir++) {
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
   void defineInflowDataStorage( Vector<RefCountedPtr<LevelData<FArrayBox>>>&  a_bdry_data,
                           const BoundaryBoxLayoutPtrVect&                     a_bdry_layout,
                           const std::string&                                  a_variable_name,
                           const FluidSpecies&                                 a_prototype_species )

   {
      const LevelData<FArrayBox>& prototype_dfn( a_prototype_species.cell_var(a_variable_name) );
      for (int i(0); i<a_bdry_layout.size(); i++) {
         const DisjointBoxLayout& dbl( a_bdry_layout[i]->disjointBoxLayout() );
         a_bdry_data.push_back(RefCountedPtr<LevelData<FArrayBox> >
                               (new LevelData<FArrayBox>(dbl, prototype_dfn.nComp(), IntVect::Zero)));
      }
   }

   inline
   void setInflowOutflowBC( LevelData<FArrayBox>&                         a_BfJ,
                      const BoundaryBoxLayoutPtrVect&                     a_all_bdry_layouts,
                      const Vector<RefCountedPtr<LevelData<FArrayBox>>>&  a_all_bdry_data,
                      const Vector<std::string>&                          a_all_bc_type,
                      const MagGeom&                                      a_geometry,
                      const LevelData<FluxBox>&                           a_velocity )
   {
      CH_TIME("FluidBCUtils::setInflowOutflowBC");
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         const DisjointBoxLayout& grids = a_BfJ.getBoxes();
            
         const LevelData<FArrayBox>& bdry_data( *(a_all_bdry_data[b]) );
         
         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
            const Box fill_box( bdry_grids[dit] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            FArrayBox& this_BfJ( a_BfJ[interior_dit] );
            
            const FArrayBox& this_inflow_data( bdry_data[dit] );
            
            const FluxBox& this_face_vel( a_velocity[interior_dit] );

            // JRA stuff for rotating momentum vector and applying symmetry BC
            // This is needed for non-planar geometries, such as anodes in dpfs
            //
            if(this_bc_type=="symmetry") {
               const MagBlockCoordSys& coord_sys = a_geometry.getBlockCoordSys(grids[interior_dit]);
               Box this_box = this_BfJ.box();
               FArrayBox this_BfJ_norm(this_box, a_BfJ.nComp());
               this_BfJ_norm.copy(this_BfJ);
               //coord_sys.projectOntoFluxSurfaceDir( this_BfJ_norm, bdry_layout.dir(), 0 );
               a_geometry.projectOntoFluxSurfaceDir( this_BfJ_norm, interior_dit, bdry_layout.dir(), 0 );
               FourthOrderBC::setInflowOutflowBC( this_BfJ_norm,
                                                  fill_box,
                                                  this_inflow_data,
                                                  this_bc_type,
                                                  this_face_vel,
                                                  bdry_layout.dir(),
                                                  bdry_layout.side() );
               //coord_sys.projectOntoFluxSurfaceDir( this_BfJ_norm, bdry_layout.dir(), 1 );
               a_geometry.projectOntoFluxSurfaceDir( this_BfJ_norm, interior_dit, bdry_layout.dir(), 1 );
               //this_BfJ_norm.setVal(0.0);
               this_BfJ.copy(this_BfJ_norm, fill_box);
            }
            else if(this_bc_type=="extrapolateJ" || this_bc_type=="neumannJ") {
               const MagBlockCoordSys& coord_sys = a_geometry.getBlockCoordSys(grids[interior_dit]);
               Box this_box = this_BfJ.box();
               FArrayBox J(this_box, a_BfJ.nComp());
               FArrayBox Xi(this_box, SpaceDim);
               coord_sys.getCellCenteredMappedCoords(Xi);
               coord_sys.pointwiseJ(J, Xi, this_box);
               this_BfJ.mult(J);
               std::string this_bc_type_sub;
               this_bc_type_sub = this_bc_type.substr(0,this_bc_type.size()-1);
               FourthOrderBC::setInflowOutflowBC( this_BfJ,
                                                  fill_box,
                                                  this_inflow_data,
                                                  this_bc_type_sub,
                                                  this_face_vel,
                                                  bdry_layout.dir(),
                                                  bdry_layout.side() );
               this_BfJ.divide(J);
            } 
            else {
               FourthOrderBC::setInflowOutflowBC( this_BfJ,
                                                  fill_box,
                                                  this_inflow_data,
                                                  this_bc_type,
                                                  this_face_vel,
                                                  bdry_layout.dir(),
                                                  bdry_layout.side() );
            } 
         }
      }
   }

   inline
   void setFluxBC( LevelData<FluxBox>& a_dst,
                   const BoundaryBoxLayoutPtrVect& a_all_bdry_layouts,
                   const Vector<std::string>& a_all_bc_type,
                   const MultiBlockCoordSys& a_coord_sys,
                   const LevelData<FluxBox>& a_flux )
   {
      CH_TIME("FluidBCUtils::setFluxBC");
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         //if(this_bc_type=="symmetry") {
         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            const Box fill_box( bdry_grids[dit] );
            FluxBox& this_a_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            const FluxBox& this_flux( a_flux[interior_dit] );
            this_a_dst.copy(this_flux, fill_box);
         }
         //}
      }
   }
   
   inline
   void setEdgeBC( LevelData<EdgeDataBox>& a_dst,
                   const BoundaryBoxLayoutPtrVect& a_all_bdry_layouts,
                   const Vector<std::string>& a_all_bc_type,
                   const MultiBlockCoordSys& a_coord_sys,
                   const LevelData<EdgeDataBox>& a_edge )
   {
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         //if(this_bc_type=="symmetry") {
         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            const int dir=1;
            const Box fill_box( bdry_grids[dit] );
            EdgeDataBox& this_a_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            const EdgeDataBox& this_edge( a_edge[interior_dit] );
            Interval Cdst(0,0);
            Interval Csrc(0,0);
            this_a_dst.copy(fill_box, Cdst, fill_box, this_edge, Csrc);
         }
         //}
      }
   }

   
}
#include "NamespaceFooter.H"

#endif
