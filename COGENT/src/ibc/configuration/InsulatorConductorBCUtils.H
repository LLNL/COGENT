#ifndef _INSU_COND_BC_UTILS_H_
#define _INSU_COND_BC_UTILS_H_

#include "Box.H"
#include "LoHiSide.H"
//#include "BlockBoundary.H"
//#include "Vector.H"
//#include "MultiBlockCoordSys.H"
#include "Directions.H"
#include "DisjointBoxLayout.H"
#include "IntVect.H"
#include "BoxIterator.H"
#include "EdgeDataBox.H"
#include "SpaceUtils.H.multidim"
#include "CONSTANTS.H"
#include "LevelData.H"

#include "MagGeomF_F.H"

//#include "BCUtils.H.multidim"
#include "BoundaryBoxLayout.H.multidim"
#include "FourthOrderBC.H.multidim"

#include "NamespaceHeader.H"

namespace InsulatorConductorBCUtils {
   
   inline
   void convertContravarToCovar( FArrayBox&  a_Fdir,
                           const FArrayBox&  a_dXdxi,
                           const int         a_dir,
                           const int         a_inverse )
   {
      CH_TIME("InsulatorConductorBCUtils::convertContravarToCovar()");
   
      // orthogonal coords assumed near boundaries

      // convert a_dir comp of contravariant vector a_F
      // to covariant or vice versa if a_inverse not zero
      CH_assert(a_Fdir.nComp() == 1);

      const Box& thisbox = a_Fdir.box();
      FORT_CONTRA_TO_COVAR_ORTHO( CHF_BOX(thisbox),
                                  CHF_CONST_FRA(a_dXdxi),
                                  CHF_CONST_INT(a_dir),
                                  CHF_CONST_INT(a_inverse),
                                  CHF_FRA1(a_Fdir,0) );
   
   }

   inline
   void setInsulatorBC( LevelData<FArrayBox>&  a_dst,
                  const LevelData<FArrayBox>&  a_src,
                  const LevelData<FArrayBox>&  a_ICbinary,
                  const BoundaryBoxLayout&     a_bdry_layout )
   {
      CH_TIME("InsulatorConductorBCUtils::setInsulatorBC() overloaded");
 
      // JRA: total hack for setting a_dst at insulator to a_src
      //
      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
            
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const Box fill_box( bdry_grids[dit] );
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         FArrayBox& this_dst( a_dst[interior_dit] );
         const FArrayBox& this_src( a_src[interior_dit] );
         const FArrayBox& this_ICbinary( a_ICbinary[interior_dit] );
            
         // grow fill box to include tranverse ghosts 
         IntVect grow_vect = a_dst.ghostVect();
         grow_vect[a_bdry_layout.dir()] = 0;
         Box fill_box_grown = grow(fill_box,grow_vect);
         
         // set boundary values in insulator region
         BoxIterator bit(fill_box_grown); // grid indices for boundary box    
         IntVect ib;         // cell index
         Real bdry_ICbinary; // 1 for insulator 0 for conductor
         Real bdry_val;
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            bdry_ICbinary = this_ICbinary.get(ib,0);
            if(bdry_ICbinary==1.0) {
               bdry_val = this_src.get(ib,0);
               this_dst.set(ib,0,bdry_val);
            }
         }
 
      }
      
   }
   
   inline
   void setInsulatorConductorBC( LevelData<FArrayBox>&  a_Bv,
                           const LevelData<FArrayBox>&  a_ICbinary,
                           const BoundaryBoxLayout&     a_bdry_layout,
                           const std::string&           a_bc_type,
                           const Real                   a_ft,
                           const RealVect&              a_Xmin_insulator,
                           const MagGeom&               a_geometry )
   {
      CH_TIME("InsulatorConductorBCUtils::setInsulatorConductorBC()");
      CH_assert(a_bc_type=="insulatorConductor_Bv");
 
      // Set insulator/conductor boundary condition for
      // virtual magnetic field. a_ICbinary is binary box containing ones
      // where an insulator is and zeros where a conductor is
      //
      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const DisjointBoxLayout& grids = a_Bv.getBoxes();
      const LevelData<FArrayBox>& Xphys = a_geometry.getCellCenteredRealCoords();
      bool twoDaxisymm;
      
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
   
         // check for axisymmetric and correct the boundary value 
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         const MagBlockCoordSys& coord_sys = a_geometry.getBlockCoordSys(grids[interior_dit]);
         twoDaxisymm = coord_sys.isAxisymmetric();
         Real bdry_val = a_ft;
         if(twoDaxisymm) bdry_val = a_ft*a_Xmin_insulator[0];
            
         // grow fill box to include tranverse ghosts 
         IntVect grow_vect = a_Bv.ghostVect();
         grow_vect[a_bdry_layout.dir()] = 0;
         const Box fill_box( bdry_grids[dit] );
         Box fill_box_grown = grow(fill_box,grow_vect);
         
         // set boundary values
         FArrayBox& this_Bv( a_Bv[interior_dit] );
         const FArrayBox& this_ICbinary( a_ICbinary[interior_dit] );
         if(twoDaxisymm) {
            Box Bv_box = this_Bv.box();
            FArrayBox this_Bv_copy(Bv_box, a_Bv.nComp());
            SpaceUtils::copy(this_Bv_copy,this_Bv,Bv_box);
            const FArrayBox& this_Xphys( Xphys[interior_dit] );
            for (int n=0; n<this_Bv_copy.nComp(); n++) this_Bv_copy.mult(this_Xphys,Bv_box,0,n,1);
            FourthOrderBC::setInsulatorConductorBC( this_Bv_copy,
                                                    this_ICbinary,
                                                    bdry_val,
                                                    fill_box_grown,
                                                    a_bdry_layout.dir(),
                                                    a_bdry_layout.side() );
            for (int n=0; n<this_Bv_copy.nComp(); n++) this_Bv_copy.divide(this_Xphys,Bv_box,0,n,1);
            SpaceUtils::copy(this_Bv,this_Bv_copy,fill_box_grown);
         }
         else {
            FourthOrderBC::setInsulatorConductorBC( this_Bv,
                                                    this_ICbinary,
                                                    bdry_val,
                                                    fill_box_grown,
                                                    a_bdry_layout.dir(),
                                                    a_bdry_layout.side() );
         } 

      }
      
   }
   
   inline
   void setInsulatorConductorBC( LevelData<FluxBox>&  a_B,
                           const LevelData<FluxBox>&  a_ICbinary_fc,
                           const BoundaryBoxLayout&   a_bdry_layout,
                           const std::string&         a_bc_type,
                           const Real                 a_ft,
                           const RealVect&            a_Xmin_insulator,
                           const MagGeom&             a_geometry )
   {
      CH_TIME("InsulatorConductorBCUtils::setInsulatorConductorBC() on faces");
      CH_assert(a_bc_type=="insulatorConductor_B");
 
      // Set insulator/conductor boundary condition for
      // magnetic field. a_ICbinary is binary box containing ones
      // where an insulator is and zeros where a conductor is

      const int bdry_dir = a_bdry_layout.dir();   
      const int bdry_side(a_bdry_layout.side());
      const Real bdry_val = a_ft;
      
      const DisjointBoxLayout& grids = a_B.getBoxes();
      a_geometry.setdXdxiOnFaces(grids, a_B.ghostVect());
      const LevelData<FluxBox>& dXdxi = a_geometry.getFaceCentereddXdxi();
   
      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
   
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         const Box& fill_box( bdry_grids[dit] );
               
         for (int dir=0; dir<SpaceDim; dir++) {
         
            FArrayBox& this_Bdir( a_B[interior_dit][dir] );
            const FArrayBox& this_ICbinary( a_ICbinary_fc[interior_dit][dir] );
         
            // copy this_Bdir to tmp and convert to covar
            FArrayBox Bdir_tmp(this_Bdir.box(),a_B.nComp());
            Bdir_tmp.copy(this_Bdir,this_Bdir.box());
            convertContravarToCovar( Bdir_tmp, dXdxi[interior_dit][dir], dir, 0 );
               
            // convert fill box to face box and grow to include tranverse ghosts 
            Box face_fill_box = surroundingNodes(fill_box,dir); 
            IntVect grow_vect = a_B.ghostVect();
            grow_vect[bdry_dir] = 0;
            face_fill_box.grow(grow_vect);

            // don't change values that live right on the boundary
            if(dir==bdry_dir && bdry_side==0) face_fill_box.growHi(dir,-1);        
            if(dir==bdry_dir && bdry_side==1) face_fill_box.growLo(dir,-1);        
 
            // set boundary values
            FourthOrderBC::setInsulatorConductorBC( Bdir_tmp,
                                                    this_ICbinary,
                                                    bdry_val,
                                                    face_fill_box,
                                                    bdry_dir,
                                                    a_bdry_layout.side() );
            
            // convert tmp back to contravar and copy back to this_Bdir
            convertContravarToCovar( Bdir_tmp, dXdxi[interior_dit][dir], dir, 1 );
            this_Bdir.copy(Bdir_tmp,face_fill_box);
         
         }
      
      }
      
   }
   
   inline
   void setInsulatorConductorBC( LevelData<EdgeDataBox>&  a_Jcovar,
                           const LevelData<EdgeDataBox>&  a_ICbinary_ec,
                           const BoundaryBoxLayout&       a_bdry_layout,
                           const std::string&             a_bc_type )
   {
      CH_TIME("InsulatorConductorBCUtils::setInsulatorConductorBC() on current density");
      CH_assert(a_bc_type=="insulatorConductor_currentDensity");
      
      // Set insulator/conductor boundary condition for
      // covariant current density on edges. a_ICbinary is 
      // binary box containing ones where an insulator is and 
      // zeros where a conductor is
      //
      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
   
         // check for axisymmetric and correct boundary value 
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         const Box fill_box( bdry_grids[dit] );
         
         // set boundary values
         EdgeDataBox& this_Jcovar( a_Jcovar[interior_dit] );
         const EdgeDataBox& this_ICbinary( a_ICbinary_ec[interior_dit] );
         FourthOrderBC::setInsulatorConductorEdgeBC( this_Jcovar,
                                                     this_ICbinary,
                                                     fill_box,
                                                     a_bc_type,
                                                     a_bdry_layout.dir(),
                                                     a_bdry_layout.side() );
      }
      
   }
   
   inline
   void setInsulatorBC( LevelData<FArrayBox>&  a_Bv,
                  const LevelData<FArrayBox>&  a_ICbinary,
                  const BoundaryBoxLayout&     a_bdry_layout,
                  const std::string&           a_bc_type,
                  const Real&                  a_ft,
                  const RealVect&              a_Xmin_insulator,
                  const MagGeom&               a_geometry )
   {
      CH_TIME("InsulatorConductorBCUtils::setInsulatorBC()");
 
      // Set insulator boundary condition for virtual magnetic field
      //
      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<FArrayBox>& Xphys_cc = a_geometry.getCellCenteredRealCoords();
 
      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const Box fill_box( bdry_grids[dit] );
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         FArrayBox& this_Bv( a_Bv[interior_dit] );
         const FArrayBox& this_ICbinary( a_ICbinary[interior_dit] );
            
         // grow fill box to include tranverse ghosts 
         IntVect grow_vect = a_Bv.ghostVect();
         grow_vect[a_bdry_layout.dir()] = 0;
         Box fill_box_grown = grow(fill_box,grow_vect);
         
         // set boundary values in insulator region
         BoxIterator bit(fill_box_grown); // grid indices for boundary box    
         IntVect ib;         // cell index
         RealVect local_Xcc; // cell center coords at cell index
         Real bdry_ICbinary; // 1 for insulator 0 for conductor
         Real bdry_val;
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            bdry_ICbinary = this_ICbinary.get(ib,0);
            if(bdry_ICbinary==1.0) {
               for(int dir=0; dir<SpaceDim; dir++) {
                  local_Xcc[dir] = Xphys_cc[interior_dit].get(ib,dir);
               } 
               bdry_val = a_ft;
               if(a_bc_type=="neumannJ") bdry_val = bdry_val*a_Xmin_insulator[0]/local_Xcc[0];
               this_Bv.set(ib,0,bdry_val);
            }
         } 

      }
      
   }
   
   inline
   void defineInsulatorConductorBinary( LevelData<FArrayBox>&  a_ICbinary,
                                  const BoundaryBoxLayout&     a_bdry_layout,
                                  const RealVect&              a_Xmin_insulator,
                                  const RealVect&              a_Xmax_insulator,
                                  const MagGeom&               a_geometry )
   {
      CH_TIME("InsulatorConductorBCUtils::defineInsulatorConductorBinary()");
 
      // Define a LevelData<FArrayBox> on bdry_grids that contains zeros 
      // at coords that touch conductors and ones at physical locations that 
      // touch insulator. (What about axis?)

      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<FArrayBox>& Xphys_cc = a_geometry.getCellCenteredRealCoords();

      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
            
         const Box fill_box( bdry_grids[dit] );
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         FArrayBox& this_ICbinary( a_ICbinary[interior_dit] );

         // grow fill box to include tranverse ghosts 
         IntVect grow_vect = a_ICbinary.ghostVect();
         grow_vect[a_bdry_layout.dir()] = 0;
         Box fill_box_grown = grow(fill_box,grow_vect);
         
         BoxIterator bit(fill_box_grown); // grid indices for boundary box    
         IntVect ib;         // cell index
         RealVect local_Xcc; // cell center coords at cell index
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            for(int dir=0; dir<SpaceDim; dir++) {
               local_Xcc[dir] = Xphys_cc[interior_dit].get(ib,dir);
            } 
            if(local_Xcc[0]<=a_Xmax_insulator[0] && local_Xcc[0]>=a_Xmin_insulator[0] 
            && local_Xcc[1]<=a_Xmax_insulator[1] && local_Xcc[1]>=a_Xmin_insulator[1]) { 
               this_ICbinary.set(ib,0,1.0); // insulator
            }
         } 

      }
      
   }
   
   inline
   void defineInsulatorConductorBinary( LevelData<FluxBox>&  a_ICbinary,
                                  const BoundaryBoxLayout&   a_bdry_layout,
                                  const RealVect&            a_Xmin_insulator,
                                  const RealVect&            a_Xmax_insulator,
                                  const MagGeom&             a_geometry )
   {
      CH_TIME("InsulatorConductorBCUtils::defineInsulatorConductorBinary()");

      // JRA, need to think about to use mapped coords here. This doesnt work for r-th!
 
      // Define LevelData<FluxBox> on bdry_grids that contains zeros 
      // at coords that touch conductors and ones at physical locations that 
      // touch insulator. (What about axis?)

      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<FluxBox>& Xphys_fc = a_geometry.getFaceCenteredRealCoords();

      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         const Box fill_box( bdry_grids[dit] );
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         
         for (int dir=0; dir<SpaceDim; dir++) {         
            FArrayBox& this_ICbinary( a_ICbinary[interior_dit][dir] );

            Box face_box = surroundingNodes(fill_box,dir);

            // grow fill box to include tranverse ghosts 
            IntVect grow_vect = a_ICbinary.ghostVect();
            grow_vect[a_bdry_layout.dir()] = 0;
            Box fill_box_grown = grow(face_box,grow_vect);
                  
            BoxIterator bit(fill_box_grown); // grid indices for boundary box    
            IntVect ib;         // cell index
            RealVect local_Xfc; // face center coords at cell index
            for(bit.begin(); bit.ok(); ++bit) {
               ib = bit();
               for(int dir0=0; dir0<SpaceDim; dir0++) {
                  local_Xfc[dir0] = Xphys_fc[interior_dit][dir].get(ib,dir0);
               } 
               if(local_Xfc[0]<=a_Xmax_insulator[0] && local_Xfc[0]>=a_Xmin_insulator[0]
               && local_Xfc[1]<=a_Xmax_insulator[1] && local_Xfc[1]>=a_Xmin_insulator[1]) { 
                  this_ICbinary.set(ib,0,1.0); // insulator
               }
            } 

         }
      }
   }
   
   inline
   void defineInsulatorConductorBinary( LevelData<EdgeDataBox>&  a_ICbinary,
                                  const BoundaryBoxLayout&       a_bdry_layout,
                                  const RealVect&                a_Xmin_insulator,
                                  const RealVect&                a_Xmax_insulator,
                                  const MagGeom&                 a_geometry )
   {
      CH_TIME("InsulatorConductorBCUtils::defineInsulatorConductorBinary()");
 
      // Define LevelData<EdgeDataBox> on bdry_grids that contains zeros 
      // at coords that touch conductors and ones at physical locations that 
      // touch insulator. (What about axis?)

      const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
      const LevelData<EdgeDataBox>& Xphys_ec = a_geometry.getEdgeCenteredRealCoords();

      for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
         const Box fill_box( bdry_grids[dit] );
         const DataIndex& interior_dit( a_bdry_layout.dataIndex(dit) );
         
         for (int dir=0; dir<SpaceDim; dir++) {         
            FArrayBox& this_ICbinary( a_ICbinary[interior_dit][dir] );
            
            Box edge_box = surroundingNodes(fill_box);
            edge_box.enclosedCells(dir);

            // grow fill box to include tranverse ghosts 
            IntVect grow_vect = a_ICbinary.ghostVect();
            grow_vect[a_bdry_layout.dir()] = 0;
            Box fill_box_grown = grow(edge_box,grow_vect);
         
            BoxIterator bit(fill_box_grown); // grid indices for boundary box    
            IntVect ib;         // cell index
            RealVect local_Xce; // cell center coords at cell index
            for(bit.begin(); bit.ok(); ++bit) {
               ib = bit();
               for(int dir0=0; dir0<SpaceDim; dir0++) {
                  local_Xce[dir0] = Xphys_ec[interior_dit][dir].get(ib,dir0);
               } 
               if(local_Xce[0]<=a_Xmax_insulator[0] && local_Xce[0]>=a_Xmin_insulator[0] 
               && local_Xce[1]<=a_Xmax_insulator[1] && local_Xce[1]>=a_Xmin_insulator[1]) { 
                  this_ICbinary.set(ib,0,1.0); // insulator
               }
            } 

         }
      }     
   }
   
}
#include "NamespaceFooter.H"

#endif
