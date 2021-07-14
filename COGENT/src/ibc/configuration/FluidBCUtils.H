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
#include "CONSTANTS.H"

#include "MagGeomF_F.H"

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
         IntVect trans_ghost_vect = prototype_dfn.ghostVect();
         trans_ghost_vect[a_bdry_layout[i]->dir()] = 0;
         const DisjointBoxLayout& dbl( a_bdry_layout[i]->disjointBoxLayout() );
         a_bdry_data.push_back(RefCountedPtr<LevelData<FArrayBox> >
                               (new LevelData<FArrayBox>(dbl, prototype_dfn.nComp(), trans_ghost_vect)));
                               //(new LevelData<FArrayBox>(dbl, prototype_dfn.nComp(), IntVect::Zero)));
      }
   }
   
   inline
   void setInflowBC( LevelData<FArrayBox>&  a_BfJ,
               const LevelData<FArrayBox>&  a_bdry_data,
               const LevelData<FluxBox>&    a_velocity,
               const BoundaryBoxLayout&     a_bdry_layout,
               const std::string&           a_bc_type )
   {
      
      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         FArrayBox& this_BfJ( a_BfJ[interior_dit] );
            
         // grow fill box to include tranverse ghosts 
         const Box fill_box( bdry_grids[dit] );
         IntVect grow_vect = a_BfJ.ghostVect();
         grow_vect[a_bdry_layout.dir()] = 0;
         Box fill_box_grown = grow(fill_box,grow_vect);
            
         const FArrayBox& this_inflow_data( a_bdry_data[dit] );
         const FluxBox& this_face_vel( a_velocity[interior_dit] );
            
         FourthOrderBC::setInflowBC( this_BfJ,
                                     fill_box_grown,
                                     this_inflow_data,
                                     a_bc_type,
                                     this_face_vel,
                                     a_bdry_layout.dir(),
                                     a_bdry_layout.side() );
      }

   }

   inline
   void projectOntoFluxSurfaceDir( FArrayBox&  a_data,
                             const FArrayBox&  a_NJinverse,
                             const FArrayBox&  a_dXdxi,
                             const int         a_dir,
                             const int         a_inverse )
   {
      CH_TIME("FluidBCUtils::projectOntoFluxSurfaceDir()");
      CH_assert(a_data.nComp() == SpaceDim);
   
      // Project the vector from phys directions to unit normal coord sys
      // with 0-component normal to the a_dir direction (or inverse)

      const FArrayBox& dxidX_on_patch( a_NJinverse );
      const FArrayBox& dXdxi_on_patch = a_dXdxi;

      const Box& box = a_data.box();
      FORT_PROJECT_PHYS_TO_MAPPED_DIR( CHF_BOX(box),
                                       CHF_CONST_FRA(dxidX_on_patch),
                                       CHF_CONST_FRA(dXdxi_on_patch),
                                       CHF_CONST_INT(a_dir),
                                       CHF_CONST_INT(a_inverse),
                                       CHF_FRA(a_data) );

   }
   
   inline
   void setCellBC( LevelData<FArrayBox>&      a_BfJ,
             const LevelData<FluxBox>&        a_boundary_values,
             const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
             const Vector<RefCountedPtr<LevelData<FArrayBox>>>&  a_all_bdry_data,
             const Vector<std::string>&       a_all_bc_type,
             const MagGeom&                   a_geometry )
   {
      CH_TIME("FluidBCUtils::setCellBC");

      // define references to geometry data used for some geometry-dependent BCs 
      const DisjointBoxLayout& grids = a_BfJ.getBoxes();
      a_geometry.setJ( grids, a_BfJ.ghostVect());
      a_geometry.setdXdxiOnCells( grids, a_BfJ.ghostVect());
      a_geometry.setPointwiseNJInverseOnCells( grids, a_BfJ.ghostVect() );
      const LevelData<FArrayBox>& Ja = a_geometry.getCellCenteredJacobian();
      //const LevelData<FArrayBox>& Xphys = a_geometry.getCellCenteredRealCoords();
      const LevelData<FArrayBox>& dxidX = a_geometry.getCellCentereddxidX();
      const LevelData<FArrayBox>& dXdxi = a_geometry.getCellCentereddXdxi();

      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
 
         const std::string this_bc_type (a_all_bc_type[b]);
         if(this_bc_type=="dirichlet" || this_bc_type=="fixed_ghosts") { // setBCs that use bdry_data
            CH_assert(!a_all_bdry_data[b].isNull())
            const LevelData<FArrayBox>& bdry_data( *(a_all_bdry_data[b]) );
            setInflowBC( a_BfJ, bdry_data, a_boundary_values, bdry_layout, this_bc_type );

         }
         else { // set BCs that don't use bdry_data

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            const FluxBox& this_face_val( a_boundary_values[interior_dit] );
            FArrayBox& this_BfJ( a_BfJ[interior_dit] );
        
            // grow fill box to include tranverse ghosts 
            const Box fill_box( bdry_grids[dit] );
            IntVect grow_vect = a_BfJ.ghostVect();
            grow_vect[bdry_layout.dir()] = 0;
            Box fill_box_grown = grow(fill_box,grow_vect);
            
            // set the BCs according to the BC type
            if(this_bc_type=="symmetry") { 
               Box this_box = this_BfJ.box();
               FArrayBox this_BfJ_copy(this_box, a_BfJ.nComp());
               SpaceUtils::copy(this_BfJ_copy,this_BfJ,this_box);
               projectOntoFluxSurfaceDir( this_BfJ_copy, dxidX[interior_dit], 
                                          dXdxi[interior_dit], bdry_layout.dir(), 0 );
               FourthOrderBC::setCellBC( this_BfJ_copy,
                                         this_face_val,
                                         fill_box_grown,
                                         this_bc_type,
                                         bdry_layout.dir(),
                                         bdry_layout.side() );
               projectOntoFluxSurfaceDir( this_BfJ_copy, dxidX[interior_dit], 
                                          dXdxi[interior_dit], bdry_layout.dir(), 1 );
               SpaceUtils::copy(this_BfJ,this_BfJ_copy,fill_box_grown);
            }
            else if(this_bc_type=="extrapolateJ" || this_bc_type=="neumannJ") {
               Box this_box = this_BfJ.box();
               FArrayBox this_BfJ_copy(this_box, a_BfJ.nComp());
               SpaceUtils::copy(this_BfJ_copy,this_BfJ,this_box);
               const FArrayBox& this_Ja( Ja[interior_dit] );
               //const FArrayBox& this_Xphys( Xphys[interior_dit] );
               const MagBlockCoordSys& coord_sys = a_geometry.getBlockCoordSys(grids[interior_dit]);
               bool twoDaxisymm = coord_sys.isAxisymmetric();
               if(twoDaxisymm) {
                  for (int n=0; n<this_BfJ_copy.nComp(); n++) {
                     this_BfJ_copy.mult(this_Ja,this_BfJ_copy.box(),0,n,1);
                     //this_BfJ_copy.mult(this_Xphys,this_BfJ_copy.box(),0,n,1);
                  }
               }
               
               std::string this_bc_type_sub;
               this_bc_type_sub = this_bc_type.substr(0,this_bc_type.size()-1);
               FourthOrderBC::setCellBC( this_BfJ_copy,
                                         this_face_val,
                                         fill_box_grown,
                                         this_bc_type_sub,
                                         bdry_layout.dir(),
                                         bdry_layout.side() );
               if(twoDaxisymm) {
                  for (int n=0; n<this_BfJ_copy.nComp(); n++) {
                     this_BfJ_copy.divide(this_Ja,this_BfJ_copy.box(),0,n,1);
                     //this_BfJ_copy.divide(this_Xphys,this_BfJ_copy.box(),0,n,1);
                  }
               }
               SpaceUtils::copy(this_BfJ,this_BfJ_copy,fill_box_grown);
            }
            else if(this_bc_type=="insulatorConductor_Bv") {
               // do nothing for this, as it will be set elsewhere
            } 
            else {
               FourthOrderBC::setCellBC( this_BfJ,
                                         this_face_val,
                                         fill_box_grown,
                                         this_bc_type,
                                         bdry_layout.dir(),
                                         bdry_layout.side() );
            } 
         }
         }
      }
   }
   
   inline
   void convertContravarToCovar( FluxBox&    a_F,
                           const FluxBox&    a_dXdxi,
                           const int         a_inverse )
   {
      CH_TIME("FluidBCUtils::convertContravarToCovar() on FluxBox");
   
      // convert a_dir comp of contravariant vector a_F
      // to covariant or vice versa if a_inverse not zero
      CH_assert(a_F.nComp() == 1);

      for (int dir=0; dir<SpaceDim; ++dir) {
         const FArrayBox& dXdxi_on_dir( a_dXdxi[dir] );
         FArrayBox& this_F_on_dir = a_F[dir];
         const Box& thisbox = a_F.box();
        
         FORT_CONTRA_TO_COVAR( CHF_BOX(thisbox),
                               CHF_CONST_FRA(dXdxi_on_dir),
                               CHF_CONST_INT(dir),
                               CHF_CONST_INT(a_inverse),
                               CHF_FRA1(this_F_on_dir,0) );
      }
   }
   
   inline
   void setFluxBC( LevelData<FluxBox>&        a_dst,
             const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
             const Vector<std::string>&       a_all_bc_type, 
             const MagGeom&                   a_geometry )
   {
      CH_TIME("FluidBCUtils::setFluxBC() from applyFluxBC()");
      
      const DisjointBoxLayout& grids = a_dst.getBoxes();
      a_geometry.setdXdxiOnFaces(grids, a_dst.ghostVect());
      const LevelData<FluxBox>& dXdxi = a_geometry.getFaceCentereddXdxi();
      
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
            const Box fill_box( bdry_grids[dit] );
            FluxBox& this_a_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            
            if(this_bc_type=="zero_on_face") {
               
               // grow fill box to include tranverse ghosts 
               const int bdry_dir = bdry_layout.dir();
               IntVect grow_vect = a_dst.ghostVect();
               grow_vect[bdry_dir] = 0;
               Box fill_box_grown = grow(fill_box,grow_vect);

               for (int dir=0; dir<SpaceDim; dir++) {
                  if (dir == bdry_layout.dir()) {
                     Box faceBox(surroundingNodes(fill_box_grown, dir));
                     this_a_dst.setVal(0.0, faceBox, dir, 0, this_a_dst.nComp());
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
               convertContravarToCovar( this_a_dst_copy, dXdxi[interior_dit], 0 );
               const int this_dir = bdry_layout.dir();
               Box this_box_dir = this_a_dst[this_dir].box();
               this_a_dst_copy[this_dir].copy(this_a_dst[this_dir],this_box_dir);
              
               FourthOrderBC::setFluxBC( this_a_dst_copy,
                                         fill_box,
                                         this_bc_type,
                                         bdry_layout.dir(),
                                         bdry_layout.side() );
              
               const int ISIDE(bdry_layout.side());
               if(ISIDE==0) this_box_dir.setBig(this_dir,fill_box.bigEnd(this_dir));
               if(ISIDE==1) this_box_dir.setSmall(this_dir,fill_box.smallEnd(this_dir));
               this_a_dst[this_dir].copy(this_a_dst_copy[this_dir],this_box_dir);
               convertContravarToCovar( this_a_dst_copy, dXdxi[interior_dit], 1 );
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
             const LevelData<FluxBox>& a_src )
   {
      CH_TIME("FluidBCUtils::setFluxBC() from setFluxBC()");
      CH_assert(a_dst.nComp() == a_src.nComp());

      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         const int bdry_dir = bdry_layout.dir();
         const int bdry_side = bdry_layout.side();

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
    
            FluxBox& this_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const FluxBox& this_src( a_src[bdry_layout.dataIndex(dit)] );
            
            if(a_src.ghostVect()==IntVect::Zero) {
               for (int dir=0; dir<CFG_DIM; dir++) {
                  FArrayBox& dst_dir(this_dst[dir]);
                  const FArrayBox& src_dir(this_src[dir]);
                  if(dir==bdry_dir) {
                     Box fill_box_mod = src_dir.box();
                     if(bdry_side==0) fill_box_mod.setBig(bdry_dir,fill_box_mod.smallEnd(bdry_dir));
                     if(bdry_side==1) fill_box_mod.setSmall(bdry_dir,fill_box_mod.bigEnd(bdry_dir));
                     dst_dir.copy(src_dir, fill_box_mod);
                  }
               }
            } 
            else {
               // grow fill box to include tranverse ghosts 
               const Box& fill_box( bdry_grids[dit] );
               IntVect grow_vect = a_dst.ghostVect();
               grow_vect[bdry_dir] = 0;
               Box fill_box_grown = grow(fill_box,grow_vect);

               const Box& src_box = this_src.box();
               if(bdry_side==1) fill_box_grown.setBig(bdry_dir,src_box.bigEnd(bdry_dir));
               if(bdry_side==0) fill_box_grown.setSmall(bdry_dir,src_box.smallEnd(bdry_dir));
               this_dst.copy(this_src, fill_box_grown);
            }
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
            
            if(this_bc_type=="zero_on_edge") {

               // grow fill box to include tranverse ghosts 
               const int bdry_dir = bdry_layout.dir();
               IntVect grow_vect = a_dst.ghostVect();
               grow_vect[bdry_dir] = 0;
               Box fill_box_grown = grow(fill_box,grow_vect);

               //this_a_dst.setVal(0.0, fill_box_grown); // set zero everywhere (does edge properly?)
               for (int dir=0; dir<SpaceDim; dir++) {
                  if (dir != bdry_layout.dir()) {
                     Box edgeBox(surroundingNodes(fill_box_grown));
                     edgeBox.enclosedCells(dir);
                     this_a_dst.setVal(0.0, edgeBox, dir, 0, this_a_dst.nComp());
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
            } 
            else if(this_bc_type=="insulatorConductor_currentDensity") {
             // set elsewhere
            }
            else {
               cout << "Warning this_bc_type = " << this_bc_type << endl; 
               cout << "Warning: No valid edge bc type specified for side "<< bdry_layout.side() << endl; 
            }
         }
      }
   }
   
   inline
   void setEdgeBC( LevelData<EdgeDataBox>& a_dst,
             const BoundaryBoxLayoutPtrVect& a_all_bdry_layouts,
             const LevelData<EdgeDataBox>& a_src )
   {
      CH_TIME("FluidBCUtils::setEdgeBC() from setEdgeBC()");
      CH_assert(a_dst.nComp() == a_src.nComp());

      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         const int bdry_dir = bdry_layout.dir();
         const int bdry_side = bdry_layout.side();

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
    
            EdgeDataBox& this_dst( a_dst[bdry_layout.dataIndex(dit)] );
            const EdgeDataBox& this_src( a_src[bdry_layout.dataIndex(dit)] );
            
            if(a_src.ghostVect()==IntVect::Zero) {
               for (int dir=0; dir<CFG_DIM; dir++) {
                  FArrayBox& dst_dir(this_dst[dir]);
                  const FArrayBox& src_dir(this_src[dir]);
                  if(dir!=bdry_dir) {
                     Box fill_box_mod = src_dir.box();
                     if(bdry_side==0) fill_box_mod.setBig(bdry_dir,fill_box_mod.smallEnd(bdry_dir));
                     if(bdry_side==1) fill_box_mod.setSmall(bdry_dir,fill_box_mod.bigEnd(bdry_dir));
                     dst_dir.copy(src_dir, fill_box_mod);
                  }
               }
            } 
            else {
               // grow fill box to include tranverse ghosts 
               const Box& fill_box( bdry_grids[dit] );
               IntVect grow_vect = a_dst.ghostVect();
               grow_vect[bdry_dir] = 0;
               Box fill_box_grown = grow(fill_box,grow_vect);

               const Box& src_box = this_src.box();
               if(bdry_side==1) fill_box_grown.setBig(bdry_dir,src_box.bigEnd(bdry_dir));
               if(bdry_side==0) fill_box_grown.setSmall(bdry_dir,src_box.smallEnd(bdry_dir));
               this_dst.copy(fill_box_grown, a_dst.interval(), fill_box_grown, this_src, a_src.interval());
            }
         }

      }
   }
   
   inline
   void setOnAxisCurlBC( LevelData<EdgeDataBox>&    a_curl,
                   const LevelData<FArrayBox>&      a_By_phys,
                   const MagGeom&                   a_geometry,
                   const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
                   const Vector<std::string>&       a_all_bc_type )
   {
      CH_TIME("FluidBCUtils::setOnAxisCurlBC()");
 
      // Compute the appropriate BC on axis for the curl of a function
      // that is odd wrt the axis for axisymmetric systems where the Jacobian is zero.
      //
      // Assumed a_curl is covariant curl(a_By_phys) on cell edges   
   
      const DisjointBoxLayout& grids = a_curl.getBoxes();

      a_geometry.setdXdxiOnEdges( grids, a_By_phys.ghostVect());
      const LevelData<EdgeDataBox>& Xphys = a_geometry.getEdgeCenteredRealCoords();
      const LevelData<EdgeDataBox>& dXdxi = a_geometry.getEdgeCentereddXdxi();
      
      for (int b(0); b<a_all_bdry_layouts.size(); b++) {
      
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
         RealVect dX_mapped;
         int bdry_dir = bdry_layout.dir();
         int curl_dir = 0;
         if(bdry_dir==0) curl_dir = 1;

         for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         
            const Box fill_box( bdry_grids[dit] );
            const DataIndex& interior_dit( bdry_layout.dataIndex(dit) );
            const FArrayBox& this_By( a_By_phys[bdry_layout.dataIndex(dit)] );
            
            const MagBlockCoordSys& coord_sys = a_geometry.getBlockCoordSys(grids[interior_dit]);
            bool twoDaxisymm = coord_sys.isAxisymmetric();
            if(!twoDaxisymm) break;

            dX_mapped = coord_sys.getMappedCellSize();
                     
            // grow fill box to include tranverse ghosts 
            IntVect grow_vect = a_curl.ghostVect();
            grow_vect[bdry_dir] = 0;
            Box fill_box_grown = grow(fill_box,grow_vect);
            
            // convert fill_box_grown from cell box to edge box
            fill_box_grown.surroundingNodes( );       // grow hi end by one in all dirs
            fill_box_grown.enclosedCells( curl_dir ); // shrink hi end by 1 in curl_dir
            
            // collapse in bdry_dir direction to 1 cell that lives on the physical boundary
            const int ISIDE(bdry_layout.side());
            if(ISIDE==1) fill_box_grown.setBig(bdry_dir,fill_box_grown.smallEnd(bdry_dir));
            if(ISIDE==0) fill_box_grown.setSmall(bdry_dir,fill_box_grown.bigEnd(bdry_dir));
            
            if( this_bc_type=="odd" || this_bc_type=="axis" ) {
               const FArrayBox& this_Xphys( Xphys[interior_dit][curl_dir] );
               const FArrayBox& this_dXdxi( dXdxi[interior_dit][curl_dir] );
                     FArrayBox& this_curl( a_curl[interior_dit][curl_dir] );

               BoxIterator bit(fill_box_grown); 
               IntVect ib, shift;
               shift = IntVect::Zero;
               shift[bdry_dir] = -ISIDE;
               Real bdry_val, By0, Xphys0, J2D, g11;
               for(bit.begin(); bit.ok(); ++bit) {
                  ib = bit();
                  g11 = this_dXdxi.get(ib,0)*this_dXdxi.get(ib,0) + this_dXdxi.get(ib,1)*this_dXdxi.get(ib,1);
                  J2D = this_dXdxi.get(ib,0)*this_dXdxi.get(ib,3) - this_dXdxi.get(ib,1)*this_dXdxi.get(ib,2);
                  Xphys0 = this_Xphys.get(ib,0);
                  By0 = this_By.get(ib+shift,0);
                  bdry_val = g11/J2D*4.0*By0/dX_mapped[bdry_dir];
                  this_curl.set(ib,0,bdry_val);
               } 
            }
         
         }

      }
      
   }
   
}
#include "NamespaceFooter.H"

#endif
