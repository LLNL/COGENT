#include "SpaceUtils.H"
#include "SpaceUtilsF_F.H"
//#include "altFaceAveragesF_F.H"
//#include "upwindSchemesF_F.H"

#include "NamespaceHeader.H"
      
void
SpaceUtils::upWindToFaces( LevelData<FluxBox>&    a_face_phi,
                     const LevelData<FArrayBox>&  a_cell_phi,
                     const LevelData<FluxBox>&    a_norm_vel,
                     const std::string&           a_method )
{
   CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );

   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FluxBox& this_face_phi( a_face_phi[dit] );
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
      const FluxBox& this_norm_vel( a_norm_vel[dit] );
      const Box& this_dbl_box( grids[dit] ); // this box has no ghost cells
      upWindToFaces( this_face_phi, this_cell_phi, this_norm_vel, this_dbl_box, a_method );
       
   } 
   //a_face_phi.exchange();

}

void
SpaceUtils::upWindToFaces( FluxBox&      this_face_phi,
                     const FArrayBox&    this_cell_phi,
                     const FluxBox&      a_norm_vel,
                     const Box&          a_dbl_box,
                     const std::string&  a_method )
{
   CH_TIME("SpaceUtils::upWindToFaces()");
   CH_assert( a_method=="c2" ||
              a_method=="TVD1st" || a_method=="TVDvanleer" ||
              a_method=="TVDminmod"  || a_method=="TVDsuperbee" ||
              a_method=="uw1" || a_method=="uw3" || a_method=="uw5" ||
              a_method=="quick" || a_method=="weno5" || a_method=="bweno");

   int this_limiter;
   if(a_method=="TVD1st")      this_limiter=0;
   if(a_method=="TVDvanleer")  this_limiter=1;
   if(a_method=="TVDminmod")   this_limiter=2;
   if(a_method=="TVDsuperbee") this_limiter=3;


   int normVelcomp = 0;
   for (int dir(0); dir<SpaceDim; dir++) {
         
      if(a_norm_vel.nComp()==SpaceDim) normVelcomp = dir; 

      Box face_box( a_dbl_box ); 
    
      // for 4th order, need extra face in mapped-grid,
      // transverse diection to handle 4th-order products  
      for (int tdir(0); tdir<SpaceDim; tdir++) {
         if (tdir!=dir) {
            const int TRANSVERSE_GROW(1);
            face_box.grow( tdir, TRANSVERSE_GROW );
         }
      }
      
      face_box.surroundingNodes( dir );
      face_box.grow( dir, 1 );

      //  now compute limited face values
         
      FArrayBox& this_face_phi_dir( this_face_phi[dir] );
      const FArrayBox& this_norm_vel_dir( a_norm_vel[dir] );
      if(a_method=="c2") {
         FORT_C2FACE( CHF_FRA( this_face_phi_dir ),
                      CHF_CONST_FRA( this_cell_phi ),
                      CHF_BOX( face_box ),
                      CHF_CONST_INT( dir ) );
      } else
      if(a_method=="uw1") {
         FORT_UW1FACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( dir ) );
      } else
      if(a_method=="uw3") {
         FORT_UW3FACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( dir ) );
      } else
      if(a_method=="uw5") {
         FORT_UW5FACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( dir ) );
      } else
      if(a_method=="quick") {
         FORT_QUICKFACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( dir ) );
      }
      if(a_method=="weno5") {
         FArrayBox this_unit_smooth(this_cell_phi.box(),SpaceDim);
         this_unit_smooth.setVal(1.0);  
         FORT_WENO5FACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_CONST_FRA( this_unit_smooth ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( dir ) );
      }
      if(a_method=="bweno") {
         FORT_BWENOFACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( dir ) );
      }
      if(a_method=="TVD1st"    || a_method=="TVDvanleer" || 
         a_method=="TVDminmod" || a_method=="TVDsuperbee") {
         FORT_TVDFACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( this_limiter ),
                       CHF_CONST_INT( dir ) );
      }
   } // end loop over directions

}

void
SpaceUtils::interpToFaces( LevelData<FluxBox>&    a_face_phi,
                     const LevelData<FArrayBox>&  a_cell_phi,
                     const LevelData<FArrayBox>&  a_cell_fun,
                     const LevelData<FluxBox>&    a_norm_vel,
                     const std::string&           a_method )
{
   CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );

   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FluxBox& this_face_phi( a_face_phi[dit] );
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
      const FArrayBox& this_cell_fun( a_cell_fun[dit] );
      const FluxBox& this_norm_vel( a_norm_vel[dit] );
      const Box& this_dbl_box( grids[dit] ); // this box has no ghost cells
      interpToFaces( this_face_phi, this_cell_phi, this_cell_fun, this_norm_vel, this_dbl_box, a_method );
       
   } 
   //a_face_phi.exchange();
}

void
SpaceUtils::interpToFaces( FluxBox&      this_face_phi,
                     const FArrayBox&    this_cell_phi,
                     const FArrayBox&    this_cell_fun,
                     const FluxBox&      a_norm_vel,
                     const Box&          a_dbl_box,
                     const std::string&  a_method )
{
   CH_TIME("SpaceUtils::interpToFaces()");
   CH_assert( a_method=="c2" || a_method == "uw1" || 
              a_method=="uw1c21st" || a_method=="uw1c2vanleer" || 
              a_method== "uw1c2minmod" || a_method == "uw1c2superbee" ||  
              a_method=="uw1c2vanAlbada1" || a_method=="uw1c2vanAlbada2" );
   int this_limiter;
   if(a_method=="uw1c21st")      this_limiter=0;
   if(a_method=="uw1c2vanleer")  this_limiter=1;
   if(a_method=="uw1c2minmod")   this_limiter=2;
   if(a_method=="uw1c2superbee") this_limiter=3;
   if(a_method=="uw1c2vanAlbada1") this_limiter=4;
   if(a_method=="uw1c2vanAlbada2") this_limiter=5;

   int normVelcomp = 0;
   for (int dir(0); dir<SpaceDim; dir++) {
         
      if(a_norm_vel.nComp()==SpaceDim) normVelcomp = dir; 

      Box face_box( a_dbl_box ); 
    
      // for 4th order, need extra face in mapped-grid,
      // transverse diection to handle 4th-order products  
      for (int tdir(0); tdir<SpaceDim; tdir++) {
         if (tdir!=dir) {
            const int TRANSVERSE_GROW(1);
            face_box.grow( tdir, TRANSVERSE_GROW );
         }
      }
      
      face_box.surroundingNodes( dir );
      face_box.grow( dir, 1 );

      //  now compute limited face values
         
      FArrayBox& this_face_phi_dir( this_face_phi[dir] );
      const FArrayBox& this_norm_vel_dir( a_norm_vel[dir] );
      if(a_method=="c2") {
         FORT_C2FACE( CHF_FRA( this_face_phi_dir ),
                      CHF_CONST_FRA( this_cell_phi ),
                      CHF_BOX( face_box ),
                      CHF_CONST_INT( dir ) );
      } else
      if(a_method=="uw1") {
         FORT_UW1FACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( dir ) );
      } 
      else {
         FORT_UW1C2FACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA( this_cell_fun ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( this_limiter ),
                         CHF_CONST_INT( dir ) );
      }
   }
}

void
SpaceUtils::interpToFacesWENO( LevelData<FluxBox>&   a_face_phi,
                         const LevelData<FArrayBox>& a_cell_phi,
                         const LevelData<FluxBox>&   a_norm_vel,
                         const LevelData<FArrayBox>& a_smooth,
                         const std::string&          a_method )
{
   CH_TIME("SpaceUtils::interpToFacesWENO()");
   CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );

   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FluxBox& this_face_phi( a_face_phi[dit] );
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
      const FluxBox& this_norm_vel( a_norm_vel[dit] );
      const FArrayBox& this_smooth( a_smooth[dit] );
      const Box& this_dbl_box( grids[dit] ); // this box has no ghost cells
      interpToFacesWENO( this_face_phi, this_cell_phi, this_norm_vel, this_smooth, this_dbl_box, a_method );
      
   } // end loop over grids
   //a_face_phi.exchange();
}

void
SpaceUtils::interpToFacesWENO( FluxBox&      this_face_phi,
                         const FArrayBox&    this_cell_phi,
                         const FluxBox&      a_norm_vel,
                         const FArrayBox&    this_smooth,
                         const Box&          a_dbl_box,
                         const std::string&  a_method )
{
   CH_TIME("SpaceUtils::interpToFacesWENO()");
   CH_assert( a_method=="weno5" || a_method=="bweno");

   
   //FluxBox normal_vel( a_norm_vel.box(), 1 );
   int normVelcomp = 0;
   for (int dir(0); dir<SpaceDim; dir++) {
         
      if(a_norm_vel.nComp()==SpaceDim) normVelcomp = dir; 
      //normal_vel.copy(a_norm_vel,normVelcomp,0,1);

      Box face_box( a_dbl_box ); 
      /*   
      for (int tdir(0); tdir<SpaceDim; tdir++) {
         if (tdir!=dir) {
            const int TRANSVERSE_GROW(1);
            face_box.grow( tdir, TRANSVERSE_GROW );
         }
      }
      */
      face_box.surroundingNodes( dir );
      face_box.grow( dir, 1 );

      // Interp fluxes from cell centers to cell faces
      //
      FArrayBox& this_face_phi_dir( this_face_phi[dir] );
      //const FArrayBox& this_norm_vel_dir( normal_vel[dir] );
      const FArrayBox& this_norm_vel_dir( a_norm_vel[dir] );
      if(a_method=="weno5") {
         FORT_WENO5FACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_CONST_FRA( this_smooth ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( dir ) );
      }
      if(a_method=="bweno") {
         FORT_BWENOFACE( CHF_FRA( this_face_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_BOX( face_box ),
                         CHF_CONST_INT( dir ) );
      }
   } // end loop over directions

}

void
SpaceUtils::interpCellToEdges( LevelData<EdgeDataBox>&   a_edge_phi,
                         const LevelData<FArrayBox>&     a_cell_phi,
                         const LevelData<EdgeDataBox>&   a_norm_vel,
                         const std::string&              a_method )
{
   CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );
   CH_assert( a_method=="c2" || 
              a_method=="TVD1st" || a_method=="TVDvanleer" ||
              a_method=="TVDminmod"  || a_method=="TVDsuperbee" ||
              a_method=="uw1" || a_method=="uw3" || a_method=="uw5" ||
              a_method=="quick" || a_method=="weno5" || a_method=="bweno");
   
   int this_limiter;
   if(a_method=="TVD1st")      this_limiter=0;
   if(a_method=="TVDvanleer")  this_limiter=1;
   if(a_method=="TVDminmod")   this_limiter=2;
   if(a_method=="TVDsuperbee") this_limiter=3;

   // Note that EdgeDataBox has extra value in each of the non-dir directions,
   // opposed to FluxBox which has extra value in dir direction. 

   const DisjointBoxLayout& grids( a_edge_phi.getBoxes() );
   int normVelcomp = 0;

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );

      EdgeDataBox& this_edge_phi( a_edge_phi[dit] );
      for (int dir(0); dir<SpaceDim; dir++) {
         
         if(a_norm_vel.nComp()==SpaceDim) normVelcomp = dir; 

         //int dir_edge = (dir + 1) % SpaceDim;
         //cout << "dir = " << dir << endl;
         //cout << "dir_edge = " << dir_edge << endl;
         
         FArrayBox& this_edge_phi_dir( this_edge_phi[dir] );
         const Box& edge_box = this_edge_phi_dir.box();
         const FArrayBox& this_norm_vel_dir( a_norm_vel[dit][dir] );
         if(a_method=="c2") {
         FORT_C2EDGE( CHF_FRA( this_edge_phi_dir ),
                      CHF_CONST_FRA( this_cell_phi ),
                      CHF_BOX( edge_box ),
                      CHF_CONST_INT( dir ) );
         } else
         if(a_method=="uw1") {
         FORT_UW1EDGE( CHF_FRA( this_edge_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( edge_box ),
                       CHF_CONST_INT( dir ) );
         } else
         if(a_method=="uw3") {
         FORT_UW3FACE( CHF_FRA( this_edge_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( edge_box ),
                       CHF_CONST_INT( dir ) );
         } else
         if(a_method=="uw5") {
         FORT_UW5FACE( CHF_FRA( this_edge_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( edge_box ),
                       CHF_CONST_INT( dir ) );
         } else
         if(a_method=="quick") {
         FORT_QUICKFACE( CHF_FRA( this_edge_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_BOX( edge_box ),
                         CHF_CONST_INT( dir ) );
         }
         if(a_method=="weno5") {
         FArrayBox this_unit_smooth(this_cell_phi.box(),SpaceDim);
         this_unit_smooth.setVal(1.0);  
         FORT_WENO5FACE( CHF_FRA( this_edge_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_CONST_FRA( this_unit_smooth ),
                         CHF_BOX( edge_box ),
                         CHF_CONST_INT( dir ) );
         }
         if(a_method=="bweno") {
         FORT_BWENOFACE( CHF_FRA( this_edge_phi_dir ),
                         CHF_CONST_FRA( this_cell_phi ),
                         CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                         CHF_BOX( edge_box ),
                         CHF_CONST_INT( dir ) );
         }
         if(a_method=="TVD1st"    || a_method=="TVDvanleer" || 
            a_method=="TVDminmod" || a_method=="TVDsuperbee") {
         FORT_TVDFACE( CHF_FRA( this_edge_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( edge_box ),
                       CHF_CONST_INT( this_limiter ),
                       CHF_CONST_INT( dir ) );
         }
      } // end loop over directions
   } // end loop over grids
   //a_edge_phi.exchange();

}

void
SpaceUtils::interpEdgesToCell( LevelData<FArrayBox>&    a_cell_phi,
                         const LevelData<EdgeDataBox>&  a_edge_phi,
                         const std::string&             a_method )
{
   CH_TIME("SpaceUtils::interpEdgesToCell()");
   CH_assert( a_cell_phi.nComp()==SpaceDim );
   //CH_assert( a_edge_phi.ghostVect()>=a_cell_phi.ghostVect() );
   CH_assert( a_method=="c2" );
   
   const DisjointBoxLayout& grids( a_edge_phi.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FArrayBox& this_cell_phi( a_cell_phi[dit] );
      //Box cell_box( grids[dit] ); // has no ghosts
      Box cell_box( this_cell_phi.box() ); // has ghosts
      for (int dir(0); dir<SpaceDim; dir++) {
         const FArrayBox& this_edge_phi_dir( a_edge_phi[dit][dir] );
         if(a_method=="c2") {
            FORT_C2CELL( CHF_BOX( cell_box ),
                         CHF_CONST_INT( dir ),
                         CHF_CONST_FRA1( this_edge_phi_dir,0 ),
                         CHF_FRA1( this_cell_phi,dir ) );
         }
      }

   }
   //a_cell_phi.exchange();

}

void
SpaceUtils::interpFacesToCell( LevelData<FArrayBox>&  a_cell_phi,
                         const LevelData<FluxBox>&    a_face_phi,
                         const std::string&           a_method )
{
   CH_TIME("SpaceUtils::interpFacesToCell()");
   CH_assert( a_cell_phi.nComp()==SpaceDim );
   //CH_assert( a_edge_phi.ghostVect()>=a_cell_phi.ghostVect() );
   CH_assert( a_method=="c2" );
   
   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FArrayBox& this_cell_phi( a_cell_phi[dit] );
      Box cell_box( grids[dit] );
      for (int dir(0); dir<SpaceDim; dir++) {
         const FArrayBox& this_face_phi_dir( a_face_phi[dit][dir] );
         if(a_method=="c2") {
            FORT_C2FACETOCELL( CHF_BOX( cell_box ),
                               CHF_CONST_INT( dir ),
                               CHF_CONST_FRA1( this_face_phi_dir,0 ),
                               CHF_FRA1( this_cell_phi,dir ) );
         }
      }

   }
   //a_cell_phi.exchange();

}

void
SpaceUtils::computeLaxSplitting( LevelData<FArrayBox>&  a_fluxR,
                                 LevelData<FArrayBox>&  a_fluxL,
                           const LevelData<FArrayBox>&  a_flux,
                           const LevelData<FArrayBox>&  a_Cspeed,
                           const LevelData<FArrayBox>&  a_fun )
{
   CH_assert( a_flux.nComp()==a_Cspeed.nComp() );

   const DisjointBoxLayout& grids( a_fun.getBoxes() );

   //   get left and right going flux at cell-center
   //   fluxR = 0.5*(flux + Cspeed*fun),
   //   fluxL = 0.5*(flux - Cspeed*fun),
   //   Cspeed = abs(max(eigenValue of Flux Jacobian))
   //
   a_fluxR.define(grids, a_flux.nComp(), a_flux.ghostVect());
   a_fluxL.define(grids, a_flux.nComp(), a_flux.ghostVect());
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_fluxR[dit].copy(a_Cspeed[dit]);
      a_fluxL[dit].copy(a_Cspeed[dit]);
      a_fluxL[dit].mult(-1.0);
      for (int n(0); n < a_flux[dit].nComp(); n++) {
         a_fluxR[dit].mult(a_fun[dit],0,n,1);
         a_fluxL[dit].mult(a_fun[dit],0,n,1);
      }
      a_fluxR[dit].plus(a_flux[dit]);
      a_fluxL[dit].plus(a_flux[dit]);
      a_fluxR[dit].mult(0.5);
      a_fluxL[dit].mult(0.5);
   }

}

void
SpaceUtils::computeLaxSplitting( FArrayBox&  a_fluxR,
                                 FArrayBox&  a_fluxL,
                           const FArrayBox&  a_flux,
                           const FArrayBox&  a_Cspeed,
                           const FArrayBox&  a_fun )
{
   CH_TIME("SpaceUtils::computeLaxSplitting()");
   CH_assert( a_flux.nComp()==a_Cspeed.nComp() );

   const Box& thisbox = a_fun.box();

   //   get left and right going flux at cell-center
   //   fluxR = 0.5*(flux + Cspeed*fun),
   //   fluxL = 0.5*(flux - Cspeed*fun),
   //   Cspeed = abs(max(eigenValue of Flux Jacobian))
   //
   a_fluxR.define(thisbox, a_flux.nComp());
   a_fluxL.define(thisbox, a_flux.nComp());

   a_fluxR.copy(a_Cspeed, a_fluxR.box());
   a_fluxL.copy(a_Cspeed, a_fluxL.box());
   //a_fluxL.mult(-1.0);
   a_fluxL.negate();
   for (int n(0); n < a_flux.nComp(); n++) {
      a_fluxR.mult(a_fun,0,n,1);
      a_fluxL.mult(a_fun,0,n,1);
   }
   a_fluxR.plus(a_flux);
   a_fluxL.plus(a_flux);
   a_fluxR.mult(0.5);
   a_fluxL.mult(0.5);

}

void
SpaceUtils::computeLaxSplitting( LevelData<FArrayBox>&  a_fluxR,
                                 LevelData<FArrayBox>&  a_fluxL,
                           const LevelData<FArrayBox>&  a_flux,
                           const LevelData<FArrayBox>&  a_Cspeed,
                           const LevelData<FArrayBox>&  a_fun,
                           const int  a_fun_comp )
{

   const DisjointBoxLayout& grids( a_fun.getBoxes() );
   LevelData<FArrayBox> this_a_fun(grids, 1, a_fun.ghostVect());
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      this_a_fun[dit].copy(a_fun[dit],a_fun_comp,0,1);
   }
   computeLaxSplitting(a_fluxR,a_fluxL,a_flux,a_Cspeed,this_a_fun);

}

void
SpaceUtils::faceInterpolate(  const int         a_dir,
                              const Box&        a_fc_box,
                              const int         a_order,
                              const FArrayBox&  a_var,
                              FArrayBox&        a_face_var )
{
  int ncomp = a_var.nComp();
  CH_assert(a_face_var.nComp() == ncomp);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert (a_dir >= 0 && a_dir < SpaceDim);

  for (int n=0; n<ncomp; n++) {
    FORT_FACE_INTERPOLATE(  CHF_CONST_INT(a_dir),
                            CHF_BOX(a_fc_box),
                            CHF_CONST_INT(a_order),
                            CHF_CONST_FRA1(a_var, n),
                            CHF_FRA1(a_face_var, n) );
  }

  return;
}


void
SpaceUtils::faceInterpolate(  const int         a_dir,
                              const Box&        a_fc_box,
                              const Box&        a_cc_box,
                              const int         a_order,
                              const FArrayBox&  a_var,
                              FArrayBox&        a_face_var )
{
  int ncomp = a_var.nComp();
  CH_assert(a_face_var.nComp() == ncomp);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert (a_dir >= 0 && a_dir < SpaceDim);

  int imin = a_cc_box.smallEnd(a_dir),
      imax = a_cc_box.bigEnd(a_dir);

  if (imin == imax) {

    /* box is flat in along this dimension,
     * so copy instead of interpolate */
    BoxIterator bit(a_cc_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      IntVect iw = iv; iw[a_dir]++;

      for (int n=0; n<ncomp; n++) {
        a_face_var(iv,n) = a_var(iv,n);
        a_face_var(iw,n) = a_var(iv,n);
      }
    }

  } else {

    faceInterpolate(a_dir, a_fc_box, a_order, a_var, a_face_var);

  }

  return;
}


void
SpaceUtils::cellCenteredGradientComponent(  const Box&        a_box,
                                            const int         a_dir,
                                            const FArrayBox&  a_var,
                                            const RealVect&   a_dx,
                                            const int         a_order,
                                            FArrayBox&        a_grad_var )
{
  CH_assert(a_var.nComp() == 1);
  CH_assert(a_grad_var.nComp() > a_dir);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  FORT_CELL_CENTERED_GRAD_COMPONENT(  CHF_BOX(a_box),
                                      CHF_CONST_INT(a_dir),
                                      CHF_CONST_FRA1(a_var,0),
                                      CHF_CONST_REALVECT(a_dx),
                                      CHF_CONST_INT(a_order),
                                      CHF_CONST_FRA1(a_grad_var, a_dir) );

  return;
}

void
SpaceUtils::faceCenteredGradientComponent(  const Box&        a_box,
                                            const int         a_dir,
                                            const FArrayBox&  a_var,
                                            const RealVect&   a_dx,
                                            const int         a_order,
                                            FArrayBox&        a_grad_var )
{
  CH_assert(a_var.nComp() == 1);
  CH_assert(a_grad_var.nComp() > a_dir);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  FORT_FACE_CENTERED_GRAD_COMPONENT(  CHF_BOX(a_box),
                                      CHF_CONST_INT(a_dir),
                                      CHF_CONST_FRA1(a_var, 0),
                                      CHF_CONST_REALVECT(a_dx),
                                      CHF_CONST_INT(a_order),
                                      CHF_FRA1(a_grad_var, a_dir) );

  return;
}

void
SpaceUtils::extrapBoundaryGhostsForCC( FArrayBox&                              a_data,
                                       const Box&                              a_interiorbox,
                                       const Box&                              a_domain_box,
                                       const int                               a_dir,
                                       const int                               a_order,
                                       const int                               a_side)
{
   // If a_order = 2, this function second-order extrapolates to fill two layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   // If a_order = 4, this function fourth-order extrapolates to fill three layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));

   int imin = a_interiorbox.smallEnd(a_dir);
   int imax = a_interiorbox.bigEnd(a_dir);

   int depth = (a_order==2)? 2: 3;
   
   Box dstbox;
   bool isbc = false;
      
   switch(a_side) 
     {
     case -1:
       // Handle the low side
       isbc = (a_interiorbox.smallEnd(a_dir) == a_domain_box.smallEnd(a_dir));
       if (isbc) {
	 dstbox = adjCellLo(a_interiorbox, a_dir, depth);
       }
       break;
     case 1:
       // Handle high side
       isbc = (a_interiorbox.bigEnd(a_dir) == a_domain_box.bigEnd(a_dir));
       if (isbc) {
	 dstbox = adjCellHi(a_interiorbox, a_dir, depth);
       }
       break;
     }

   if (isbc) {

     CH_assert(bx.contains(dstbox));
     
     if (imin == imax) {

       BoxIterator bit(dstbox);
       for (bit.begin(); bit.ok(); ++bit) {
	 IntVect i1 = bit();
	 IntVect i2 = i1; i2[a_dir] = imin;
	 for (int n=0; n<a_data.nComp(); n++) {
	   a_data(i1,n) = a_data(i2,n);
	 }
       }
       
     } else {

       CH_assert( (a_order==4 && a_interiorbox.size(a_dir)>=5)
		  || (a_order==2 && a_interiorbox.size(a_dir)>=3));
       FORT_EXTRAP_FOR_CC_OPS(CHF_CONST_INT(a_dir),
			      CHF_CONST_INT(a_side),
			      CHF_CONST_INT(a_order),
			      CHF_BOX(dstbox),
			      CHF_BOX(a_interiorbox),
			      CHF_FRA(a_data));

     }
   }

   return;
} 

void
SpaceUtils::extrapBoundaryGhostsForFC( FArrayBox&                              a_data,
                                       const Box&                              a_interiorbox,
                                       const Box&                              a_domain_box,
                                       const int                               a_dir,
                                       const int                               a_order,
				       const int                               a_side)
{
   // This function fourth-order extrapolates to fill two layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));

   int imin = a_interiorbox.smallEnd(a_dir);
   int imax = a_interiorbox.bigEnd(a_dir);

   Box dstbox;
   bool isbc = false;
      
   switch(a_side) 
     {
     case -1:
       // Handle the low side
       isbc = (a_interiorbox.smallEnd(a_dir) == a_domain_box.smallEnd(a_dir));
       if (isbc) {
	 dstbox = adjCellLo(a_interiorbox, a_dir, 2);
       }
       break;
     case 1:
       // Handle high side
       isbc = (a_interiorbox.bigEnd(a_dir) == a_domain_box.bigEnd(a_dir));
       if (isbc) {
	 dstbox = adjCellHi(a_interiorbox, a_dir, 2);
       }
       break;
     }
   
   if (isbc) {
     
     CH_assert(bx.contains(dstbox));

     if (imin == imax) {

       BoxIterator bit(dstbox);
       for (bit.begin(); bit.ok(); ++bit) {
	 IntVect i1 = bit();
	 IntVect i2 = i1; i2[a_dir] = imin;
	 for (int n=0; n<a_data.nComp(); n++) {
	   a_data(i1,n) = a_data(i2,n);
	 }
       }
       
     } else {
       
       CH_assert( (a_order==4 && a_interiorbox.size(a_dir)>=5)
		  || (a_order==2 && a_interiorbox.size(a_dir)>=3));
       FORT_EXTRAP_FOR_FC_OPS(CHF_CONST_INT(a_dir),
			      CHF_CONST_INT(a_side),
			      CHF_CONST_INT(a_order),
			      CHF_BOX(dstbox),
			      CHF_BOX(a_interiorbox),
			      CHF_FRA(a_data));
       
       
     }
   }

   return;
} 


void
SpaceUtils::secondOrderTransExtrapAtDomainBdry(FArrayBox&           a_data,
                                               const int            a_dir,
                                               const Box&           a_interiorbox,
                                               const ProblemDomain& a_domain,
                                               const int            a_maxdim )
{
   CH_TIME("secondOrderTransExtrapAtDomainBdry");

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));
   const Box& domBox = a_domain.domainBox();

   for (int tdir=0; tdir<a_maxdim; ++tdir)
   {
      if (tdir != a_dir && !a_domain.isPeriodic(tdir))
      {
         CH_assert(a_interiorbox.size(tdir)>=3);

         for (int side=-1; side<2; side+=2)
         {
            Box dstbox;
            bool isbc = false;

            switch(side)
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(tdir) == domBox.smallEnd(tdir));
               if (isbc)
               {
                  dstbox = adjCellLo(a_interiorbox,tdir,1);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(tdir) == domBox.bigEnd(tdir));
               if (isbc)
               {
                  dstbox = adjCellHi(a_interiorbox,tdir,1);
               }
               break;
            }

            if (isbc)
            {
               CH_assert(bx.contains(dstbox));
               FORT_SECOND_ORDER_EXTRAPOLATION( CHF_CONST_INT(tdir),
                                                CHF_CONST_INT(side),
                                                CHF_BOX(a_interiorbox),
                                                CHF_BOX(dstbox),
                                                CHF_FRA(a_data) );
            }
         }
      }
   }

   return;
}



void 
SpaceUtils::fillGhostCellsSimple( FArrayBox& a_phi,
                                  const Box& a_bx_int,
                                  const int  a_dir )
{
  CH_assert((a_dir >=0) && (a_dir < SpaceDim));
  const Box& phi_bx = a_phi.box();

  for (int side=-1; side<2; side+=2) {

    int ii;
    Box bdrybox;

    int n_gpt;
    if (side == -1) {
      n_gpt = a_bx_int.smallEnd(a_dir) - phi_bx.smallEnd(a_dir);
      if (n_gpt > 0) {
        ii = a_bx_int.smallEnd(a_dir);
        bdrybox = adjCellLo(a_bx_int, a_dir, n_gpt);
      }
    } else {
      n_gpt = phi_bx.bigEnd(a_dir) - a_bx_int.bigEnd(a_dir);
      if (n_gpt > 0) {
        ii = a_bx_int.bigEnd(a_dir);
        bdrybox = adjCellHi(a_bx_int, a_dir, n_gpt);
      }
    }

    if (n_gpt > 0) {
      BoxIterator bit(bdrybox);
      for (bit.begin(); bit.ok(); ++bit) {
  
        IntVect iv(bit());
        int j = (Real) (side < 0 ? ii-iv[a_dir] : iv[a_dir]-ii );
  
        IntVect i_int_0(iv); i_int_0[a_dir] = ii;
        IntVect i_int_1(iv); i_int_1[a_dir] = ii - side;
        IntVect i_int_2(iv); i_int_2[a_dir] = ii - 2*side;
        IntVect i_int_3(iv); i_int_3[a_dir] = ii - 3*side;
  
        Real c0 = ((1.0+j)*(2.0+j)*(3.0+j))/6.0;
        Real c1 = -(j*(2.0+j)*(3.0+j))/2.0;
        Real c2 = (j*(1.0+j)*(3.0+j))/2.0;
        Real c3 = -(j*(1.0+j)*(2.0+j))/6.0;
  
        for (int n=0; n<a_phi.nComp(); n++) {
          a_phi(iv,n) =   c0 * a_phi(i_int_0,n)
                        + c1 * a_phi(i_int_1,n)
                        + c2 * a_phi(i_int_2,n)
                        + c3 * a_phi(i_int_3,n);
        }
  
      }
    }
  }

  return;
}

void
SpaceUtils::copyAndFillGhostCellsSimple(  LevelData<FArrayBox>&       a_var_wg,
                                          const LevelData<FArrayBox>& a_var )
{
  const DisjointBoxLayout& var_dbl = a_var.disjointBoxLayout();

  /* interior */
  for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
    a_var_wg[dit].setVal(0.0);
    a_var_wg[dit].copy(a_var[dit], var_dbl[dit]);
  }
  a_var_wg.exchange();

  /* codim-1 boundaries */
  for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
    for (int dir=0; dir<SpaceDim; dir++) {
      fillGhostCellsSimple(a_var_wg[dit], var_dbl[dit], dir);
    }
  }
  a_var_wg.exchange();

  /* higher codim boundaries */
  for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
    Box bx_int(var_dbl[dit]);
    const Box& bx_var(a_var_wg[dit].box());
    for (int dir=0; dir<SpaceDim-1; dir++) {
      bx_int.growLo(dir, (bx_int.smallEnd(dir)-bx_var.smallEnd(dir)) );
      bx_int.growHi(dir, (bx_var.bigEnd(dir)-bx_int.bigEnd(dir)) );
      fillGhostCellsSimple(a_var_wg[dit], bx_int, dir+1);
    }
  }
  a_var_wg.exchange();

  /* done - hopefully! */
  return;
}

void 
SpaceUtils::inspectFArrayBox(const LevelData<FArrayBox>&  a_F0,
                             const int                    a_comp)
{
   CH_assert(a_comp<a_F0.nComp());
   const DisjointBoxLayout& grids( a_F0.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& F0_on_patch = a_F0[dit];   
      const Box& thisbox = F0_on_patch.box();

      FORT_INSPECT_FARRAYBOX( CHF_BOX(thisbox), 
                              CHF_CONST_FRA1(F0_on_patch,a_comp) );
   }
}

void 
SpaceUtils::inspectFluxBox(const LevelData<FluxBox>&  a_Flux,
                           const int                  a_dir)
{
   const DisjointBoxLayout& grids( a_Flux.getBoxes() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FluxBox& Flux_on_patch = a_Flux[dit]; 
      //const Box& thisbox = Flux_on_patch.box();

      //Box face_box( Flux_on_patch.box());
      //face_box.surroundingNodes( a_dir );
      //face_box.grow( a_dir, 1 );

      const FArrayBox& Flux_on_dir = Flux_on_patch[a_dir]; 
      const Box& thisbox = Flux_on_dir.box();

      FORT_INSPECT_FLUXBOX( CHF_BOX(thisbox), 
                            CHF_CONST_FRA(Flux_on_dir),
                            CHF_CONST_INT(a_dir) );
   }
}

void 
SpaceUtils::inspectEdgeDataBox(const LevelData<EdgeDataBox>&  a_Edge,
                               const int                      a_dir)
{
   const DisjointBoxLayout& grids( a_Edge.getBoxes() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const EdgeDataBox& Edge_on_patch = a_Edge[dit]; 
      const Box& thisbox = Edge_on_patch.box();

      const FArrayBox& Edge_on_dir = Edge_on_patch[a_dir]; 

      FORT_INSPECT_FLUXBOX( CHF_BOX(thisbox), 
                            CHF_CONST_FRA(Edge_on_dir),
                            CHF_CONST_INT(a_dir) );
   }
}

#if 1
void
SpaceUtils::copy(  FArrayBox& a_dst,
		   const FArrayBox& a_src )
{
  CH_assert(a_dst.nComp() == a_src.nComp());
  CH_assert(a_src.box().contains(a_dst.box()));
  FORT_COPY(CHF_BOX(a_dst.box()), 
	    CHF_FRA( a_dst ),
	    CHF_CONST_FRA( a_src ));

}
#endif

#include "NamespaceFooter.H"
