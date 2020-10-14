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
#include "BoxIterator.H"
#include "FluidSpecies.H"
#include "EdgeDataBox.H"
#include "SpaceUtils.H"

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
   
         // JRA: need to set some metric terms outside dit loop for
         // call to projectOntoFluxSurfaceDir()
         // (poor design...think how to fix)
         //
         if(this_bc_type=="symmetry") { 
            const IntVect ghosts = a_BfJ.ghostVect();   
            a_geometry.setPointwiseNJInverseOnCells(grids, ghosts); 
            a_geometry.setdXdxiOnCells(grids, ghosts);
         }
         //
         /////////////////////////////////////////////////

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
            const Box fill_box( bdry_grids[dit] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            FArrayBox& this_BfJ( a_BfJ[interior_dit] );
            
            const FArrayBox& this_inflow_data( bdry_data[dit] );
            const FluxBox& this_face_vel( a_velocity[interior_dit] );

            // JRA: Rotating momentum vector and applying symmetry BC.
            // This is needed for non-planar geometries, such as anodes in dpfs
            //
            if(this_bc_type=="symmetry") { 
               Box this_box = this_BfJ.box();
               FArrayBox this_BfJ_copy(this_box, a_BfJ.nComp());
               SpaceUtils::copy(this_BfJ_copy,this_BfJ,this_box);
               a_geometry.projectOntoFluxSurfaceDir( this_BfJ_copy, interior_dit, 
                                                     bdry_layout.dir(), 0 );
               FourthOrderBC::setInflowOutflowBC( this_BfJ_copy,
                                                  fill_box,
                                                  this_inflow_data,
                                                  this_bc_type,
                                                  this_face_vel,
                                                  bdry_layout.dir(),
                                                  bdry_layout.side() );
               a_geometry.projectOntoFluxSurfaceDir( this_BfJ_copy, interior_dit, 
                                                     bdry_layout.dir(), 1 );
               SpaceUtils::copy(this_BfJ,this_BfJ_copy,fill_box);
            }
            else if(this_bc_type=="extrapolateJ" || this_bc_type=="neumannJ") {
               Box this_box = this_BfJ.box();
               FArrayBox this_BfJ_copy(this_box, a_BfJ.nComp());
               SpaceUtils::copy(this_BfJ_copy,this_BfJ,this_box);
               a_geometry.multJ(this_BfJ_copy,interior_dit);
               std::string this_bc_type_sub;
               this_bc_type_sub = this_bc_type.substr(0,this_bc_type.size()-1);
               FourthOrderBC::setInflowOutflowBC( this_BfJ_copy,
                                                  fill_box,
                                                  this_inflow_data,
                                                  this_bc_type_sub,
                                                  this_face_vel,
                                                  bdry_layout.dir(),
                                                  bdry_layout.side() );
               a_geometry.divideJ(this_BfJ_copy,interior_dit);
               SpaceUtils::copy(this_BfJ,this_BfJ_copy,fill_box);
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
   void setFluxBC( LevelData<FluxBox>&        a_dst,
             const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
             const Vector<std::string>&       a_all_bc_type, 
             const MagGeom&                   a_geometry )
   {
      CH_TIME("FluidBCUtils::setFluxBC() from applyFluxBC()");
      
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         if( this_bc_type=="symmetry" || this_bc_type=="wall" ) {
            const DisjointBoxLayout& grids = a_dst.getBoxes();
            const IntVect ghosts = a_dst.ghostVect();   
            a_geometry.setdXdxiOnFaces(grids, ghosts);
         }
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
            const Box fill_box( bdry_grids[dit] );
            FluxBox& this_a_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            
            if(this_bc_type=="zero_on_face") {
               for (int dir=0; dir<SpaceDim; dir++) {
                  if (dir == bdry_layout.dir()) {
                     this_a_dst.setVal(0.0, fill_box, dir, 0, this_a_dst.nComp());
                  }
               }
            } 
            else if( this_bc_type=="symmetry" || this_bc_type=="wall" ) {
               
               // Passed FluxBox is contravariant on faces
               // Boundary condition is applied to covariant vector for non-dir directions
               //
               Box this_box = this_a_dst.box();
               FluxBox this_a_dst_copy(this_box, a_dst.nComp());
               this_a_dst_copy.copy(this_a_dst,this_box);
               a_geometry.convertContravarToCovar( this_a_dst_copy, interior_dit, 0 );
               const int this_dir = bdry_layout.dir(); // use contra for dir direction
               Box this_box_dir = this_a_dst[this_dir].box();
               this_a_dst_copy[this_dir].copy(this_a_dst[this_dir],this_box_dir);
               FourthOrderBC::setFluxBC( this_a_dst_copy,
               //FourthOrderBC::setFluxBC( this_a_dst,
                                         fill_box,
                                         this_bc_type,
                                         bdry_layout.dir(),
                                         bdry_layout.side() );
               const int ISIDE(bdry_layout.side());
               if(ISIDE==0) this_box_dir.setBig(this_dir,fill_box.bigEnd(this_dir));
               if(ISIDE==1) this_box_dir.setSmall(this_dir,fill_box.smallEnd(this_dir));
               this_a_dst[this_dir].copy(this_a_dst_copy[this_dir],this_box_dir);
               a_geometry.convertContravarToCovar( this_a_dst_copy, interior_dit, 1 );
               for(int dir=0; dir<SpaceDim; ++dir) {
                  if(dir!=this_dir) {
                     this_box_dir = this_a_dst[dir].box();
                     if(ISIDE==0) this_box_dir.setBig(this_dir,fill_box.bigEnd(this_dir));
                     if(ISIDE==1) this_box_dir.setSmall(this_dir,fill_box.smallEnd(this_dir));
                     this_a_dst[dir].copy(this_a_dst_copy[dir],this_box_dir);
                  }
               }
            
            } 
            else {
               cout << "Warning: No valid flux bc type specified for side "<< bdry_layout.side() << endl; 
            }
            
         }
      }
   }

   inline
   void setFluxBC( LevelData<FluxBox>& a_dst,
                   const BoundaryBoxLayoutPtrVect& a_all_bdry_layouts,
                   const LevelData<FluxBox>& a_flux )
   {
      CH_TIME("FluidBCUtils::setFluxBC() from setFluxBC()");
      CH_assert(a_dst.nComp() == a_flux.nComp());
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            const Box fill_box( bdry_grids[dit] );
            FluxBox& this_a_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            const FluxBox& this_flux( a_flux[interior_dit] );
            this_a_dst.copy(this_flux, fill_box);
         }
      }
   }
   
   inline
   void setBC( LevelData<FArrayBox>&      a_dst,
         const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
         const Vector<std::string>&       a_all_bc_type,
         const MultiBlockCoordSys&        a_coord_sys )
   {
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            //const int dir=1;
            const Box fill_box( bdry_grids[dit] );
            FArrayBox& this_a_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            //if(this_bc_type=="natural") {
                
               FourthOrderBC::setBC( this_a_dst,
                                     fill_box,
                                     this_bc_type,
                                     bdry_layout.dir(),
                                     bdry_layout.side() );
               
            //} 
         }
      }
   }
   
   inline
   void setEdgeBC( LevelData<EdgeDataBox>&    a_dst,
             const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
             const Vector<std::string>&       a_all_bc_type )
   {
      CH_TIME("FluidBCUtils::setEdgeBC() from applyEdgeBC()");
      
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
            const Box fill_box( bdry_grids[dit] );
            EdgeDataBox& this_a_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            
            if(this_bc_type=="zero_on_edge") {
               for (int dir=0; dir<SpaceDim; dir++) {
                  if (dir != bdry_layout.dir()) {
                     this_a_dst.setVal(0.0, fill_box, dir, 0, this_a_dst.nComp());
                  }
               }
            } 
            else if( this_bc_type=="neumann_on_edge" || 
                     this_bc_type=="natural"         ||
                     this_bc_type=="symmetry" ) {
               
               FourthOrderBC::setEdgeBC( this_a_dst,
                                         fill_box,
                                         this_bc_type,
                                         bdry_layout.dir(),
                                         bdry_layout.side() );
               //   int thisdir = 0;
               //   this_a_dst.setVal(0.0,fill_box,thisdir,0,this_a_dst.nComp());
            } 
            else {
               cout << "Warning: No valid edge bc type specified for side "<< bdry_layout.side() << endl; 
            }
         }
      }
   }
   
   inline
   void setEdgeBC( LevelData<EdgeDataBox>&    a_dst,
             const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
             const LevelData<EdgeDataBox>&    a_edge )
   {
      CH_TIME("FluidBCUtils::setEdgeBC() from setEdgeBC()");
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            const Box fill_box( bdry_grids[dit] );
            EdgeDataBox& this_a_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            const EdgeDataBox& this_edge( a_edge[interior_dit] );
            Interval Cdst(0,0);
            Interval Csrc(0,0);
            this_a_dst.copy(fill_box, Cdst, fill_box, this_edge, Csrc);
            /*
            for (int dir=0; dir<SpaceDim; dir++) {
               if (dir != bdry_layout.dir()) {
                  const FArrayBox& edge_on_dir = this_edge[dir]; 
                  FArrayBox& dst_on_dir = this_a_dst[dir]; 
                  const Box dst_box( dst_on_dir.box() );
                  const Box src_box( edge_on_dir.box() );
                  //dst_on_dir.copy(edge_on_dir, fill_box, 0, fill_box, 0, 1);
                  dst_on_dir.copy(src_box, Cdst, fill_box, edge_on_dir, Csrc);
               }
            }
            */
         }
      }
   }
   
}
#include "NamespaceFooter.H"

#endif
