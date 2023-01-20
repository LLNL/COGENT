#include "SpaceUtils.H.multidim"
#include "SpaceUtilsF_F.H"
#include "CONSTANTS.H"
//#include "altFaceAveragesF_F.H"
//#include "upwindSchemesF_F.H"

#include "NamespaceHeader.H"
      
void
SpaceUtils::upWindToFaces( LevelData<FluxBox>&    a_face_phi,
                     const LevelData<FArrayBox>&  a_cell_phi,
                     const LevelData<FluxBox>&    a_norm_vel,
                     const std::string&           a_method,
                     const bool&                  a_fourth_order)
{
   CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );
   
   if (a_fourth_order) {
      CH_assert( a_cell_phi.ghostVect() >= 2*IntVect::Unit );
      CH_assert( a_face_phi.ghostVect() >= IntVect::Unit );
   }

   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FluxBox& this_face_phi( a_face_phi[dit] );
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
      const FluxBox& this_norm_vel( a_norm_vel[dit] );
      const Box& this_dbl_box( grids[dit] ); // this box has no ghost cells
      upWindToFaces( this_face_phi, this_cell_phi, this_norm_vel, this_dbl_box, a_method, a_fourth_order );
       
   } 

}

void
SpaceUtils::upWindToFaces( FluxBox&      this_face_phi,
                     const FArrayBox&    this_cell_phi,
                     const FluxBox&      a_norm_vel,
                     const Box&          a_dbl_box,
                     const std::string&  a_method,
                     const bool&         a_fourth_order)
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
      // transverse direction to handle 4th-order products
      if (a_fourth_order) {
         for (int tdir(0); tdir<SpaceDim; tdir++) {
            if (tdir!=dir) {
               const int TRANSVERSE_GROW(1);
               face_box.grow( tdir, TRANSVERSE_GROW );
            }
         }
      }
      
      face_box.surroundingNodes( dir );

      //  now compute limited face values
         
      FArrayBox& this_face_phi_dir( this_face_phi[dir] );
      const FArrayBox& this_norm_vel_dir( a_norm_vel[dir] );
      if(a_method=="c2") {
         FORT_C2FACE( CHF_FRA( this_face_phi_dir ),
                      CHF_CONST_FRA( this_cell_phi ),
                      CHF_BOX( face_box ),
                      CHF_CONST_INT( dir ) );
      }
      if(a_method=="uw1") {
         FORT_UW1FACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( dir ) );
      }
      if(a_method=="uw3") {
         FORT_UW3FACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( dir ) );
      }
      if(a_method=="uw5") {
         FORT_UW5FACE( CHF_FRA( this_face_phi_dir ),
                       CHF_CONST_FRA( this_cell_phi ),
                       CHF_CONST_FRA1( this_norm_vel_dir, normVelcomp ),
                       CHF_BOX( face_box ),
                       CHF_CONST_INT( dir ) );
      }
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
SpaceUtils::constrainedCrossProductOnEdges( LevelData<EdgeDataBox>&  a_edge_E,
                                      const LevelData<FluxBox>&      a_face_B,
                                      const LevelData<FArrayBox>&    a_cell_V,
                                      const LevelData<FArrayBox>&    a_cell_Cspeed,
                                      const std::string&             a_method )
{
   CH_TIME("SpaceUtils::constrainedCrossProductOnEdges()");
 
   // Compute the cross product, E=BxV, on cell edges in such a way to ensure div of curl(E)
   // is identically zero while properly upwinding. 
   //
   // 1) Interpolate non-dir components of vector V at cell center to dir face
   // 2) Multiply by dir component of vector B on face. 
   // 3) Upwind both VB terms on each face to the appropriate surrounding edges.
   //

   //CH_assert( a_cell_V.ghostVect()>=2*IntVect::Unit );
   //CH_assert( a_face_B.ghostVect()>=2*IntVect::Unit );
   CH_assert( a_face_B.nComp()==1 );
   CH_assert( a_cell_V.nComp()==SpaceDim );
   CH_assert( a_edge_E.nComp()==1 );

   int numGhosts = 2;
   if(a_method == "weno5" || a_method == "bweno" || a_method == "uw5") numGhosts = 3;

   const DisjointBoxLayout& grids( a_cell_V.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
            EdgeDataBox&  this_edge_E( a_edge_E[dit] );
            this_edge_E.setVal(0.0);
      const FluxBox&   this_face_B( a_face_B[dit] );
      const FArrayBox& this_cell_V( a_cell_V[dit] );
      const FArrayBox& this_cell_Cspeed( a_cell_Cspeed[dit] );
            FluxBox    VB_faces(this_face_B.box(),1); // container for face center VB comps
      Box face_box; 
      Box edge_box; 
      int signB=-1, signV=-1, sign;
      for(int Bdir=0; Bdir<SpaceDim; Bdir++) {
         
         signB = -1*signB;
         for(int Vdir=0; Vdir<SpaceDim; Vdir++) {
            if(Vdir!=Bdir) {
               signV = -1*signV;
               // Interpolate Vdir comp of V at cell to Bdir face and 
               // multiply by Bdir comp of B
               //
               face_box = grids[dit];             // no ghosts
               face_box.surroundingNodes(Bdir);   // grow hi end by one in Bdir direction
               face_box.grow(Vdir,numGhosts);     // grow each end by numGhosts in Vdir direction
               FORT_CELLFACEPRODUCT( CHF_BOX( face_box ),
                                     CHF_CONST_INT( Bdir ),
                                     CHF_CONST_FRA1( this_face_B[Bdir],0 ),
                                     CHF_CONST_FRA1( this_cell_V,Vdir ),
                                     CHF_FRA1( VB_faces[Bdir],0 ) );
               
               // perform Lax flux splitting
               //
               FArrayBox fluxR_face(face_box,1);  
               FArrayBox fluxL_face(face_box,1);  
               FArrayBox Cspeed_Vdir_face(face_box,1); 
               fluxR_face.setVal(1.0); 
               //Cspeed_Vdir_face.setVal(1.0); 
               FORT_CELLFACEPRODUCT( CHF_BOX( face_box ),
                                     CHF_CONST_INT( Bdir ), 
                                     CHF_CONST_FRA1( fluxR_face,0 ), 
                                     CHF_CONST_FRA1( this_cell_Cspeed,Vdir ), 
                                     CHF_FRA1( Cspeed_Vdir_face,0 ) );
               computeLaxSplitting( fluxR_face, fluxL_face, VB_faces[Bdir],
                                    Cspeed_Vdir_face, this_face_B[Bdir] );

               // upwind Left/Right VB on Bdir face in Vdir direction to Edir edges
               //
               for(int Edir=0; Edir<SpaceDim; Edir++) {
                  if(Edir!=Vdir && Edir!=Bdir) {
                     edge_box = grids[dit]; // no ghosts
                     edge_box.surroundingNodes( );   // grow hi end by one in all dirs
                     edge_box.enclosedCells( Edir ); // shrink hi end by 1 in Edir
                     FArrayBox fluxR_edges(edge_box,1); 
                     FArrayBox fluxL_edges(edge_box,1); 
                     FArrayBox Cspeed_Vdir_edges(edge_box,1); 
                     interpToEdges(Cspeed_Vdir_edges, Cspeed_Vdir_face, Cspeed_Vdir_edges, edge_box, Vdir, "c2");
                     //Cspeed_Vdir_edges.setVal(1.0);
                     interpToEdges(fluxR_edges, fluxR_face, Cspeed_Vdir_edges, edge_box, Vdir, a_method);
                     Cspeed_Vdir_edges.negate();
                     interpToEdges(fluxL_edges, fluxL_face, Cspeed_Vdir_edges, edge_box, Vdir, a_method);
                     sign = signV*signB;              
                     if(sign==1) {
                        this_edge_E[Edir].plus(fluxR_edges,0,0,1);
                        this_edge_E[Edir].plus(fluxL_edges,0,0,1);
                     }
                     if(sign==-1) {
                        this_edge_E[Edir].minus(fluxR_edges,0,0,1);
                        this_edge_E[Edir].minus(fluxL_edges,0,0,1);
                     }
                  }
               }
            }
         }
      }
   } 

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
}


void
SpaceUtils::interpToFaces( LevelData<FluxBox>&    a_face_phi,
                     const LevelData<FArrayBox>&  a_cell_phi,
                     const int                    a_order )
{
   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );

   if (a_order == 4) {
      CH_assert(a_cell_phi.ghostVect() >= 2*IntVect::Unit );
   }
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {

       for (int dir=0; dir<SpaceDim; dir++) {
       
       Box box(grids[dit]);
          
       if (a_order == 4 && (a_face_phi.ghostVect() >= IntVect::Unit)) {
          for (int tdir(0); tdir<SpaceDim; tdir++) {
             if (tdir!=dir) {
                const int TRANSVERSE_GROW(1);
                box.grow( tdir, TRANSVERSE_GROW );
             }
          }
       }
       
       SpaceUtils::faceInterpolate(dir,
                                   surroundingNodes(box,dir),
                                   a_order,
                                   a_cell_phi[dit],
                                   a_face_phi[dit][dir] );
       }
    }
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

      Box face_box( a_dbl_box );         // no ghosts 
      face_box.surroundingNodes( dir );  // grow hi end by 1 in dir direction
    
      // compute limited face values
      //
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
      
   }

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

   
   int normVelcomp = 0;
   for (int dir(0); dir<SpaceDim; dir++) {
         
      if(a_norm_vel.nComp()==SpaceDim) normVelcomp = dir; 

      Box face_box( a_dbl_box ); 
      face_box.surroundingNodes( dir );

      // Interp fluxes from cell centers to cell faces
      //
      FArrayBox& this_face_phi_dir( this_face_phi[dir] );
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
   }

}

void
SpaceUtils::interpToEdges( FArrayBox&    a_Vphi,
                     const FArrayBox&    a_phi,
                     const FArrayBox&    a_V,
                     const Box&          a_box,
                     const int           a_dir,
                     const std::string&  a_method )
{
   CH_TIME("SpaceUtils::interpToEdges()");
   CH_assert( a_method=="c2" ||
              a_method=="TVD1st" || a_method=="TVDvanleer" ||
              a_method=="TVDminmod"  || a_method=="TVDsuperbee" ||
              a_method=="uw1" || a_method=="uw3" || a_method=="uw5" ||
              a_method=="quick" || a_method=="weno5" );

   int this_limiter;
   if(a_method=="TVD1st")      this_limiter=0;
   if(a_method=="TVDvanleer")  this_limiter=1;
   if(a_method=="TVDminmod")   this_limiter=2;
   if(a_method=="TVDsuperbee") this_limiter=3;
     
   int normVelcomp = 0; 
   if(a_method=="c2") {
      FORT_C2FACE( CHF_FRA( a_Vphi ),
                   CHF_CONST_FRA( a_phi ),
                   CHF_BOX( a_box ),
                   CHF_CONST_INT( a_dir ) );
   }
   if(a_method=="uw1") {
      FORT_UW1FACE( CHF_FRA( a_Vphi ),
                    CHF_CONST_FRA( a_phi ),
                    CHF_CONST_FRA1( a_V, normVelcomp ),
                    CHF_BOX( a_box ),
                    CHF_CONST_INT( a_dir ) );
   }
   if(a_method=="uw3") {
      FORT_UW3FACE( CHF_FRA( a_Vphi ),
                    CHF_CONST_FRA( a_phi ),
                    CHF_CONST_FRA1( a_V, normVelcomp ),
                    CHF_BOX( a_box ),
                    CHF_CONST_INT( a_dir ) );
   }
   if(a_method=="uw5") {
      FORT_UW5FACE( CHF_FRA( a_Vphi ),
                    CHF_CONST_FRA( a_phi ),
                    CHF_CONST_FRA1( a_V, normVelcomp ),
                    CHF_BOX( a_box ),
                    CHF_CONST_INT( a_dir ) );
   }
   if(a_method=="quick") {
      FORT_QUICKFACE( CHF_FRA( a_Vphi ),
                      CHF_CONST_FRA( a_phi ),
                      CHF_CONST_FRA1( a_V, normVelcomp ),
                      CHF_BOX( a_box ),
                      CHF_CONST_INT( a_dir ) );
   }
   if(a_method=="weno5") {
         FArrayBox this_unit_smooth(a_phi.box(),SpaceDim);
         this_unit_smooth.setVal(1.0);  
         FORT_WENO5FACE( CHF_FRA( a_Vphi ),
                         CHF_CONST_FRA( a_phi ),
                         CHF_CONST_FRA1( a_V, normVelcomp ),
                         CHF_CONST_FRA( this_unit_smooth ),
                         CHF_BOX( a_box ),
                         CHF_CONST_INT( a_dir ) );
   }
   if(a_method=="TVD1st"    || a_method=="TVDvanleer" || 
      a_method=="TVDminmod" || a_method=="TVDsuperbee") {
      FORT_TVDFACE( CHF_FRA( a_Vphi ),
                    CHF_CONST_FRA( a_phi ),
                    CHF_CONST_FRA1( a_V, normVelcomp ),
                    CHF_BOX( a_box ),
                    CHF_CONST_INT( this_limiter ),
                    CHF_CONST_INT( a_dir ) );
   }

}

void
SpaceUtils::interpCellToEdges( LevelData<EdgeDataBox>&   a_edge_phi,
                         const LevelData<FArrayBox>&     a_cell_phi,
                         const LevelData<EdgeDataBox>&   a_norm_vel,
                         const std::string&              a_method )
{
   const IntVect cellgv = a_cell_phi.ghostVect();
   const IntVect edgegv = a_edge_phi.ghostVect();
   CH_assert( cellgv >= edgegv || cellgv >= IntVect::Unit );
   //CH_assert( a_cell_phi.ghostVect() >= a_edge_phi.ghostVect() );
   //CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );  
   CH_assert( a_method=="c2" ); 
           
   IntVect grow_vect = IntVect::Zero; 
   if(edgegv == cellgv) grow_vect = edgegv-IntVect::Unit;
  
   // Note that EdgeDataBox has extra value in each of the non-dir directions,
   // opposed to FluxBox which has extra value in dir direction. 

   const DisjointBoxLayout& grids( a_edge_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
      EdgeDataBox& this_edge_phi( a_edge_phi[dit] );
      
      for (int dir(0); dir<SpaceDim; dir++) {
         Box edge_box;
         if(cellgv > edgegv )
            edge_box = this_edge_phi[dir].box();    // has ghost 
         else {
            edge_box = grow(grids[dit],grow_vect);  // grown cell box 
            edge_box.surroundingNodes( );  // grow hi end by one in all dirs
            edge_box.enclosedCells( dir ); // shrink hi end by 1 in dir
         }
         
         
         FArrayBox& this_edge_phi_dir( this_edge_phi[dir] );
         if(SpaceDim==2) {
            FORT_C2_CELL_TO_EDGE_2D( CHF_FRA( this_edge_phi_dir ),
                                     CHF_CONST_FRA( this_cell_phi ),
                                     CHF_BOX( edge_box ),
                                     CHF_CONST_INT( dir ) );
         }
         else {
            int dirj = dir + 1;
                dirj = dirj % SpaceDim;
            int dirk = dir + 2;
                dirk = dirk % SpaceDim;
            FORT_C2_CELL_TO_EDGE_3D( CHF_BOX( edge_box ),
                                     CHF_CONST_INT( dirj ),
                                     CHF_CONST_INT( dirk ),
                                     CHF_CONST_FRA1( this_cell_phi,dir ),
                                     CHF_FRA1( this_edge_phi_dir,0 ) );
         }
      } // end loop over directions
   } // end loop over grids

}

void
SpaceUtils::interpEdgesToCell( LevelData<FArrayBox>&    a_cell_phi,
                         const LevelData<EdgeDataBox>&  a_edge_phi,
                         const std::string&             a_method )
{
   CH_TIME("SpaceUtils::interpEdgesToCell()");
   const int cell_ncomp = a_cell_phi.nComp();
   const int edge_ncomp = a_edge_phi.nComp();
   const IntVect cellgv = a_cell_phi.ghostVect();
   const IntVect edgegv = a_edge_phi.ghostVect();
   CH_assert( (cell_ncomp==SpaceDim && edge_ncomp==1 )
            || cell_ncomp==edge_ncomp );
   CH_assert( a_method=="c2" );

   const DisjointBoxLayout& grids( a_edge_phi.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& this_cell_phi( a_cell_phi[dit] );
      Box cell_box( this_cell_phi.box() ); // has ghosts
      if(cellgv > edgegv ) cell_box = grids[dit];    // no ghost 

      if(cell_ncomp==SpaceDim && edge_ncomp==1) { // edge_phi is vector
         for (int dir(0); dir<SpaceDim; dir++) {
            const FArrayBox& this_edge_phi_dir( a_edge_phi[dit][dir] );
            if(a_method=="c2" && SpaceDim<3) {
               FORT_C2_EDGE_TO_CELL_2D( CHF_BOX( cell_box ),
                                        CHF_CONST_INT( dir ),
                                        CHF_CONST_FRA1( this_edge_phi_dir,0 ),
                                        CHF_FRA1( this_cell_phi,dir ) );
            } 
            else {
               int dirj = dir + 1;
                   dirj = dirj % SpaceDim;
               int dirk = dir + 2;
                   dirk = dirk % SpaceDim;
               FORT_C2_EDGE_TO_CELL_3D( CHF_BOX( cell_box ),
                                        CHF_CONST_INT( dirj ),
                                        CHF_CONST_INT( dirk ),
                                        CHF_CONST_FRA1( this_edge_phi_dir,0 ),
                                        CHF_FRA1( this_cell_phi,dir ) );
            }
         }
      }
      else { // edge_phi is scalar
         this_cell_phi.setVal(0.0);
         for (int dir(0); dir<SpaceDim; dir++) {
            const FArrayBox& this_edge_phi_dir( a_edge_phi[dit][dir] );
            if(a_method=="c2") {
               FORT_EDGE_SCALAR_TO_CELL( CHF_BOX( cell_box ),
                                         CHF_CONST_INT( dir ),
                                         CHF_CONST_FRA( this_edge_phi_dir ),
                                         CHF_FRA( this_cell_phi ) );
            }
         }
      }

   }

}

void
SpaceUtils::interpEdgesToEdges( LevelData<EdgeDataBox>&  a_edge_out,
                          const LevelData<EdgeDataBox>&  a_edge_in,
                          const std::string&             a_method )
{
   CH_TIME("SpaceUtils::interpEdgesToEdges()");
   CH_assert( a_edge_in.nComp()==1 );
   CH_assert( a_edge_in.ghostVect() >= 1*IntVect::Unit);
   CH_assert( a_method=="c2" );
   
   // interpolate values that live on one edge to the other edges
   // example: Ein = [Ex_onz, Ez_onx], Eout = [Ex_onx, Ez_onz]
   //
   // Hack for 2D right now... need to think about 3D
   //
   
   int dir0=1;
   const DisjointBoxLayout& grids( a_edge_out.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      for (int dir(0); dir<SpaceDim; dir++) {
         
         Box edge_box( grids[dit] ); // no ghosts
         edge_box.surroundingNodes( );  // grow hi end by one in all dirs
         edge_box.enclosedCells( dir ); // shrink hi end by 1 in dir
        
         if(dir==1) dir0=0;
         const FArrayBox& this_edge_in(  a_edge_in[dit][dir0] );
               FArrayBox& this_edge_out( a_edge_out[dit][dir] );
        
         FORT_C2_EDGES_TO_EDGES( CHF_BOX( edge_box ),
                                 CHF_CONST_FRA1( this_edge_in,0 ),
                                 CHF_FRA1( this_edge_out,0 ),
                                 CHF_CONST_INT( dir ) );
      }

   }

}

void
SpaceUtils::interpNodesToEdges( LevelData<EdgeDataBox>&    a_edgePhi_out,
                          const LevelData<NodeFArrayBox>&  a_nodePhi_in,
                          const std::string&               a_method )
{
   CH_TIME("SpaceUtils::interpNodesToEdges()");
   CH_assert( a_edgePhi_out.nComp()==a_nodePhi_in.nComp() );
   CH_assert( a_method=="c2" );
   
   // interpolate values that live cell nodes to cell edges
   //
   const DisjointBoxLayout& grids( a_edgePhi_out.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      for (int dir(0); dir<SpaceDim; dir++) {
         
         Box edgebox( grids[dit] ); // no ghosts
         edgebox.surroundingNodes( );  // grow hi end by one in all dirs
         edgebox.enclosedCells( dir ); // shrink hi end by 1 in dir
        
         const FArrayBox& this_node_in(  a_nodePhi_in[dit].getFab() );
               FArrayBox& this_edge_out( a_edgePhi_out[dit][dir] );
        
         FORT_C2_NODES_TO_EDGES( CHF_BOX( edgebox ),
                                 CHF_CONST_INT( dir ),
                                 CHF_CONST_FRA( this_node_in ),
                                 CHF_FRA( this_edge_out ) );
      }

   }
   
}

void
SpaceUtils::interpCellsToNodes( LevelData<NodeFArrayBox>&  a_node_phi,
                          const LevelData<FArrayBox>&      a_cell_phi,
                          const std::string&               a_method )
{
   const IntVect cellgv = a_cell_phi.ghostVect();
   const IntVect nodegv = a_node_phi.ghostVect();
   CH_assert( cellgv > nodegv || cellgv >= IntVect::Unit );
   CH_assert( a_cell_phi.nComp() == a_node_phi.nComp() );
   CH_assert( a_method=="c2" );
   
   
   const DisjointBoxLayout& grids( a_node_phi.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
            FArrayBox& this_node_phi( a_node_phi[dit].getFab() );
      
      Box node_box;
      if(cellgv > nodegv )
         node_box = this_cell_phi.box();    // has ghost 
      else {  
         node_box = grids[dit];    // no ghost 
      }
      node_box.surroundingNodes( );  // grow hi end by one in all dirs

      FORT_C2_CELLS_TO_NODES( CHF_BOX( node_box ),
                             CHF_CONST_FRA( this_cell_phi ),
                             CHF_FRA( this_node_phi ) );
   
   }

}

void
SpaceUtils::interpNodesToCells( LevelData<FArrayBox>&      a_cellPhi_out,
                          const LevelData<NodeFArrayBox>&  a_nodePhi_in,
                          const std::string&               a_method )
{
   CH_TIME("SpaceUtils::interpNodesToCells()");
   CH_assert( a_cellPhi_out.nComp()==a_nodePhi_in.nComp() );
   CH_assert( a_nodePhi_in.ghostVect() >= a_cellPhi_out.ghostVect() );
   CH_assert( a_method=="c2" );
   
   // interpolate values that live cell nodes to cell centers
   //
   const DisjointBoxLayout& grids( a_cellPhi_out.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& this_node_in(  a_nodePhi_in[dit].getFab() );
            FArrayBox& this_cell_out( a_cellPhi_out[dit] );
      
      //Box cellbox( grids[dit] ); // no ghosts
      //cellbox.grow( 1 );         // add one layer of ghost cells
      const Box& cellbox( this_cell_out.box() );        

      FORT_C2_NODES_TO_CELLS( CHF_BOX( cellbox ),
                              CHF_CONST_FRA( this_node_in ),
                              CHF_FRA( this_cell_out ) );
   }

}

void
SpaceUtils::ParaEdgeGradientAtNodes( NodeFArrayBox&  a_gradPhi,
                               const EdgeDataBox&    a_edgePhi,
                               const RealVect&       a_dx,
                               const Box&            a_box,
                               const std::string&    a_method )
{
   CH_TIME("SpaceUtils::ParaEdgeGradientAtNodes()");
   CH_assert( a_method=="c2" );
   
   // calculate normal gradient at nodes from edge data
   // i.e.
   // gradPhi(comp=0) = d(edgePhi0)/dX0
   // gradPhi(comp=1) = d(edgePhi1)/dX1
   //
      
   FArrayBox& this_gradPhi( a_gradPhi.getFab() );
   for (int dir(0); dir<SpaceDim; dir++) {
         
      double thisdx = a_dx[dir];
      const FArrayBox& this_edgePhi( a_edgePhi[dir] );
        
      FORT_EDGE_GRAD_AT_NODES( CHF_BOX( a_box ),
                               CHF_CONST_INT( dir ),
                               CHF_CONST_REAL( thisdx ),
                               CHF_CONST_FRA1( this_edgePhi,0 ),
                               CHF_FRA1( this_gradPhi,dir ) );
   }

}

void
SpaceUtils::PerpEdgeGradientAtCells( FArrayBox&    a_gradPhi,
                               const EdgeDataBox&  a_edgePhi,
                               const RealVect&     a_dx,
                               const Box&          a_box,
                               const std::string&  a_method )
{
   CH_TIME("SpaceUtils::PerpEdgeGradientAtCells()");
   CH_assert( a_gradPhi.nComp()==SpaceDim );
   CH_assert( a_method=="c2" );
   
   // calculate cross gradient comps at cells from edge data
   // i.e.
   // gradPhi(comp=0) = d(edgePhi0)/dX1
   // gradPhi(comp=1) = d(edgePhi1)/dX0
   //
          
   int dir0 = 1;
   for (int dir(0); dir<SpaceDim; dir++) {
         
      const FArrayBox& this_edgePhi( a_edgePhi[dir] );

      if(dir==1) dir0 = 0;
      double thisdx = a_dx[dir0];

      FORT_EDGE_GRAD_AT_CELLS( CHF_BOX( a_box ),
                               CHF_CONST_INT( dir0 ),
                               CHF_CONST_REAL( thisdx ),
                               CHF_CONST_FRA1( this_edgePhi,0 ),
                               CHF_FRA1( a_gradPhi,dir ) );
   }
      
}

void
SpaceUtils::interpFaceVectorToCell( LevelData<FArrayBox>&  a_cell_phi,
                              const LevelData<FluxBox>&    a_face_phi,
                              const std::string&           a_method )
{
   CH_TIME("SpaceUtils::interpFaceVectorToCell()");
   CH_assert( a_cell_phi.nComp() == SpaceDim && a_face_phi.nComp() == 1 );
   CH_assert( a_method=="c2" );

   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& this_cell_phi( a_cell_phi[dit] );
      Box cell_box( grids[dit] );
      if(a_face_phi.ghostVect()>=a_cell_phi.ghostVect()) {
         cell_box = this_cell_phi.box();
      }
      for (int dir(0); dir<SpaceDim; dir++) {
         const FArrayBox& this_face_phi_dir( a_face_phi[dit][dir] );
         FORT_C2_FACE_TO_CELL( CHF_BOX( cell_box ),
                               CHF_CONST_INT( dir ),
                               CHF_CONST_FRA1( this_face_phi_dir,0 ),
                               CHF_FRA1( this_cell_phi,dir ) );
      }
      
   } 

}

void
SpaceUtils::interpFaceScalarToCell( LevelData<FArrayBox>&  a_cell_phi,
                              const LevelData<FluxBox>&    a_face_phi,
                              const std::string&           a_method )
{
   CH_TIME("SpaceUtils::interpFaceScalarToCell()");
   CH_assert( a_cell_phi.nComp() == a_face_phi.nComp() );
   CH_assert( a_method=="c2" );

   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {

      FArrayBox& this_cell_phi( a_cell_phi[dit] );
      Box cell_box( grids[dit] );
      if(a_face_phi.ghostVect()>=a_cell_phi.ghostVect()) {
         cell_box = this_cell_phi.box();
      }
      for (int dir(0); dir<SpaceDim; dir++) {
         const FArrayBox& this_face_phi_dir( a_face_phi[dit][dir] );
         FORT_FACE_SCALAR_TO_CELL( CHF_BOX( cell_box ),
                                   CHF_CONST_INT( dir ),
                                   CHF_CONST_FRA( this_face_phi_dir ),
                                   CHF_FRA( this_cell_phi ) );
      }  
      
   } 

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
      
   if (a_side == -1) { 
      // Handle the low side
      isbc = (a_interiorbox.smallEnd(a_dir) == a_domain_box.smallEnd(a_dir));
      if (isbc) {
         dstbox = adjCellLo(a_interiorbox, a_dir, 2);
      }
   } else if (a_side == 1) {
      // Handle high side
      isbc = (a_interiorbox.bigEnd(a_dir) == a_domain_box.bigEnd(a_dir));
      if (isbc) {
         dstbox = adjCellHi(a_interiorbox, a_dir, 2);
      }
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
       CH_assert(    (a_order==4 && a_interiorbox.size(a_dir)>=5)
		              || (a_order==2 && a_interiorbox.size(a_dir)>=3));
       FORT_EXTRAP_FOR_FC_OPS( CHF_CONST_INT(a_dir),
                               CHF_CONST_INT(a_side),
                               CHF_CONST_INT(a_order),
                               CHF_BOX(dstbox),
                               CHF_BOX(a_interiorbox),
                               CHF_FRA(a_data) );
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

   for (int tdir=0; tdir<a_maxdim; ++tdir) {
      if (tdir != a_dir && !a_domain.isPeriodic(tdir)) {

         int imin = a_interiorbox.smallEnd(tdir);
         int imax = a_interiorbox.bigEnd(tdir);

         for (int side=-1; side<2; side+=2) {
            Box dstbox;
            bool isbc = false;

            if (side == -1) {
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(tdir) == domBox.smallEnd(tdir));
               if (isbc) {
                  dstbox = adjCellLo(a_interiorbox,tdir,1);
               }
            } else if (side == 1) {
               // Handle high side
               isbc = (a_interiorbox.bigEnd(tdir) == domBox.bigEnd(tdir));
               if (isbc) {
                  dstbox = adjCellHi(a_interiorbox,tdir,1);
               }
            } 
            
            if (isbc) {
               CH_assert(bx.contains(dstbox));
               if (imin == imax) {
                  BoxIterator bit(dstbox);
                  for (bit.begin(); bit.ok(); ++bit) {
                     IntVect i1 = bit();
                     IntVect i2 = i1; i2[tdir] = imin;
                     for (int n=0; n<a_data.nComp(); n++) {
                        a_data(i1,n) = a_data(i2,n);
                     }
                  }
               } else {
                  CH_assert(a_interiorbox.size(tdir)>=3);
                  FORT_SECOND_ORDER_EXTRAPOLATION( CHF_CONST_INT(tdir),
                                                   CHF_CONST_INT(side),
                                                   CHF_BOX(a_interiorbox),
                                                   CHF_BOX(dstbox),
                                                   CHF_FRA(a_data) );
               }
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
SpaceUtils::copyAndFillGhostCellsSimple(  LevelData<FluxBox>&       a_var_wg,
                                          const LevelData<FluxBox>& a_var )
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
      for (int dim = 0; dim < SpaceDim; dim++) {
        fillGhostCellsSimple(a_var_wg[dit][dim], var_dbl[dit], dir);
      }
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
      for (int dim = 0; dim < SpaceDim; dim++) {
        fillGhostCellsSimple(a_var_wg[dit][dim], bx_int, dir+1);
      }
    }
  }
  a_var_wg.exchange();

  /* done - hopefully! */
  return;
}

void
SpaceUtils::copyAndFillGhostCellsSimple(  LevelData<EdgeDataBox>&       a_var_wg,
                                          const LevelData<EdgeDataBox>& a_var )
{
  const DisjointBoxLayout& var_dbl = a_var.disjointBoxLayout();

  /* interior */
  for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
    for (int dim = 0; dim < SpaceDim; dim++) {
      a_var_wg[dit][dim].setVal(0.0);
      a_var_wg[dit][dim].copy(a_var[dit][dim], var_dbl[dit]);
    }
  }
  a_var_wg.exchange();

  /* codim-1 boundaries */
  for (DataIterator dit(a_var.dataIterator()); dit.ok(); ++dit) {
    for (int dir=0; dir<SpaceDim; dir++) {
      for (int dim = 0; dim < SpaceDim; dim++) {
        fillGhostCellsSimple(a_var_wg[dit][dim], var_dbl[dit], dir);
      }
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
      for (int dim = 0; dim < SpaceDim; dim++) {
        fillGhostCellsSimple(a_var_wg[dit][dim], bx_int, dir+1);
      }
    }
  }
  a_var_wg.exchange();

  /* done - hopefully! */
  return;
}

void 
SpaceUtils::inspectFArrayBox(const LevelData<FArrayBox>&  a_F0,
                             const int                    a_comp, 
                             const int                    a_ghosts )
{
   CH_assert(a_comp<a_F0.nComp());
   CH_assert(a_F0.ghostVect() >= a_ghosts*IntVect::Unit);
   const DisjointBoxLayout& grids( a_F0.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& F0_on_patch = a_F0[dit];   
      Box cellbox( grids[dit] ); // no ghosts
      if(a_ghosts) cellbox.grow( a_ghosts );
      //cellbox.setSmall(0,cellbox.bigEnd(0)-3);

      FORT_INSPECT_FARRAYBOX( CHF_BOX(cellbox), 
                              CHF_CONST_FRA1(F0_on_patch,a_comp) );
   }
}

void 
SpaceUtils::inspectFluxBox(const LevelData<FluxBox>&  a_Flux,
                           const int                  a_dir,
                           const int                  a_ghosts )
{
   const DisjointBoxLayout& grids( a_Flux.getBoxes() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FluxBox& Flux_on_patch = a_Flux[dit]; 

      Box facebox( grids[dit] );
      facebox.surroundingNodes( a_dir );
      //facebox.setSmall(0,facebox.bigEnd(0)-3);

      const FArrayBox& Flux_on_dir = Flux_on_patch[a_dir]; 
      if(a_ghosts) facebox = Flux_on_dir.box();

      FORT_INSPECT_FLUXBOX( CHF_BOX(facebox), 
                            CHF_CONST_FRA(Flux_on_dir),
                            CHF_CONST_INT(a_dir) );
   }
}

void 
SpaceUtils::inspectEdgeDataBox( const LevelData<EdgeDataBox>&  a_Edge,
                                const int                      a_dir,
                                const int                      a_ghosts )
{
   const DisjointBoxLayout& grids( a_Edge.getBoxes() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const EdgeDataBox& Edge_on_patch = a_Edge[dit]; 
      const FArrayBox& Edge_on_dir = Edge_on_patch[a_dir]; 

      Box edgebox( grids[dit] );
      edgebox.surroundingNodes( );  // grow hi end by one in all dirs
      edgebox.enclosedCells( a_dir ); // shrink hi end by 1 in dir
      if(a_ghosts) edgebox.grow( a_ghosts );
      //edgebox.setBig(1,0);

      //const Box& edgebox = Edge_on_dir.box();
      FORT_INSPECT_FLUXBOX( CHF_BOX(edgebox), 
                            CHF_CONST_FRA(Edge_on_dir),
                            CHF_CONST_INT(a_dir) );
   }
}

void
SpaceUtils::copyEdgeDataBox( EdgeDataBox& a_dst,
		       const EdgeDataBox& a_src )
{
  
  CH_assert(a_dst.nComp() == a_src.nComp());
  for (int dir=0; dir<SpaceDim; dir++) {
     FArrayBox& dst_dir = a_dst[dir];
     const FArrayBox& src_dir = a_src[dir];
   
     Box thisbox(dst_dir.box()); // use the smaller box
     if(!src_dir.box().contains(dst_dir.box())) thisbox = src_dir.box();

     copy( dst_dir, src_dir, thisbox );
  }

}

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

void
SpaceUtils::copy(  FArrayBox& a_dst,
             const FArrayBox& a_src,
             const Box&       a_box )
{
  CH_assert(a_dst.nComp() == a_src.nComp());
  CH_assert(a_dst.box().contains(a_box));
  CH_assert(a_src.box().contains(a_box));
  FORT_COPY(CHF_BOX( a_box ),
            CHF_FRA( a_dst ),
            CHF_CONST_FRA( a_src ));

}

void
SpaceUtils::copyNodeToCell(  FArrayBox&  a_dst,
                       const FArrayBox&  a_src,
                       const Box&        a_box,
                       const int         a_dir )
{ 
  CH_assert(a_dst.nComp() == a_src.nComp());
  CH_assert(a_dst.box().type(a_dir) == 0);
  CH_assert(a_src.box().type(a_dir) == 1);
  //cout << "copyNodeToCell: src.box().type() = " << a_src.box().type() << endl;
  FORT_COPY(CHF_BOX( a_box ),
            CHF_FRA( a_dst ),
            CHF_CONST_FRA( a_src ));

}

void
SpaceUtils::localVectorNorm(  FArrayBox& a_dst,
                        const FArrayBox& a_src )
{ 
  CH_assert(a_dst.nComp() == 1);
  CH_assert(a_src.box().contains(a_dst.box()));
  FORT_VECTOR_NORM( CHF_BOX(a_dst.box()),
                    CHF_FRA1( a_dst,0 ), 
                    CHF_CONST_FRA( a_src ) );

}

void
SpaceUtils::exchangeFluxBox( LevelData<FluxBox>& a_Face )
{
   CH_TIME("SpaceUtils::exchangeFluxBox()");
   
   const DisjointBoxLayout& grids( a_Face.getBoxes() );
   DataIterator dit = grids.dataIterator();

   // copy to a temporary
   //
   LevelData<FluxBox> tmp_Face(grids, a_Face.nComp());
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) { 
         const FArrayBox& Face_on_dir     = a_Face[dit][dir]; 
               FArrayBox& tmp_Face_on_dir = tmp_Face[dit][dir]; 
         copy(tmp_Face_on_dir,Face_on_dir);
      }
   }

   // call exchange on passed FluxBox
   //
   a_Face.exchange();

   // now use the original data stored in temp to modify
   // the exchanged Face data. In this example, I'm copying data
   // on the low-side grid edge from temp->Face, which has the
   // effect of declaring the high-side value on the shared edge
   // to be the "correct" value

   for (dit.begin(); dit.ok(); ++dit) {

      FluxBox& thisTemp = tmp_Face[dit];
      FluxBox& thisFace = a_Face[dit];

      for (int dir=0; dir<SpaceDim; dir++) {
        
         // define a face-centered box on the low-side of this
         // grid box
         Box loFaceBox=thisTemp[dir].box();
         loFaceBox.setBig(dir,loFaceBox.smallEnd(dir));

         copy(thisFace[dir],thisTemp[dir],loFaceBox);

      } // end loop over directions

   } // end loop over grids on this level

}

void
SpaceUtils::exchangeEdgeDataBox( LevelData<EdgeDataBox>& a_Edge )
{
   CH_TIME("SpaceUtils::exchangeEdgeDataBox()");
   
   const DisjointBoxLayout& grids( a_Edge.getBoxes() );
   DataIterator dit = grids.dataIterator();

   // copy to a temporary
   //
   LevelData<EdgeDataBox> tmp_Edge(grids, a_Edge.nComp());
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) { 
         const FArrayBox& Edge_on_dir     = a_Edge[dit][dir]; 
               FArrayBox& tmp_Edge_on_dir = tmp_Edge[dit][dir]; 
         copy(tmp_Edge_on_dir,Edge_on_dir);
      }
   }

   // call exchange on passed EdgeDataBox
   //
   a_Edge.exchange();

   // now use the original data stored in temp to modify
   // the exchanged Edge data. In this example, I'm copying data
   // on the low-side grid edge from temp->Edge, which has the
   // effect of declaring the high-side value on the shared edge
   // to be the "correct" value

   for (dit.begin(); dit.ok(); ++dit) {
      EdgeDataBox& thisTemp = tmp_Edge[dit];
      EdgeDataBox& thisEdge = a_Edge[dit];

      for (int dir=0; dir<SpaceDim; dir++) {
         // define an edge-centered box on the low-side of this
         // grid box
         for (int dir0=0; dir0<SpaceDim; dir0++) {
             if(dir0!=dir) {
                Box loEdgeBox=thisTemp[dir].box();
                loEdgeBox.setBig(dir0,loEdgeBox.smallEnd(dir0));

                // now copy from temp->Edge on the low edge.
                // we can be a bit cavalier about this, because
                // on the edges which are not shared between grids,
                // exchange shouldn't have changed anything anyway
                //thisEdge[dir].copy(thisTemp[dir],loEdgeBox);
                copy(thisEdge[dir],thisTemp[dir],loEdgeBox);
             }
         }

      } // end loop over directions

   } // end loop over grids on this level

}

void
SpaceUtils::incrementLevelData( LevelData<FArrayBox>&       a_dst,
                                const LevelData<FArrayBox>& a_src,
                                const Real                  a_a )
{
  for (auto dit(a_dst.dataIterator()); dit.ok(); ++dit) {
    a_dst[dit].plus(a_src[dit], a_a);
  }
}

void
SpaceUtils::incrementLevelData( LevelData<FluxBox>&       a_dst,
                                const LevelData<FluxBox>& a_src,
                                const Real                a_a )
{
  for (auto dit(a_dst.dataIterator()); dit.ok(); ++dit) {
    for (int dim(0); dim<SpaceDim; dim++) {
      a_dst[dit][dim].plus(a_src[dit][dim], a_a);
    }
  }
}

void
SpaceUtils::incrementLevelData( LevelData<EdgeDataBox>&       a_dst,
                                const LevelData<EdgeDataBox>& a_src,
                                const Real                    a_a )
{
  for (auto dit(a_dst.dataIterator()); dit.ok(); ++dit) {
    for (int dim(0); dim<SpaceDim; dim++) {
      a_dst[dit][dim].plus(a_src[dit][dim], a_a);
    }
  }
}

void
SpaceUtils::findMax( const LevelData<FArrayBox>& a,
                    double&                     a_max_val,
                    IntVect&                    a_max_iv )
{
   CH_assert(a.nComp()==1);
   const DisjointBoxLayout& grids = a.disjointBoxLayout();
   
   double local_max = -DBL_MAX;
   IntVect local_max_loc = IntVect::Zero;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      IntVect max_loc = a[dit].maxIndex(0);
      double max_val = a[dit](max_loc,0);
      if (max_val > local_max) {
         local_max = max_val;
         local_max_loc = max_loc;
      }
   }

   struct {
      double val;
      int   rank;
   } in, out;
   int myrank;
 
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   in.val = local_max;
   in.rank = myrank;
   MPI_Allreduce( &in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD );

   a_max_val = out.val;

   int iv[CH_SPACEDIM];

   if (myrank == out.rank) {
      for (int i=0; i<SpaceDim; ++i) iv[i] = local_max_loc[i];
   }

   MPI_Bcast(iv, SpaceDim, MPI_INT, out.rank, MPI_COMM_WORLD);

   for (int i=0; i<SpaceDim; ++i) a_max_iv[i] = iv[i];
}


Real
SpaceUtils::MaxNorm( const LevelData<FArrayBox>& a )
{
   const DisjointBoxLayout& grids = a.disjointBoxLayout();

   double local_max = -DBL_MAX;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      double this_max = a[dit].max();
      if (this_max > local_max) local_max = this_max;
   }

   double global_max;
#ifdef CH_MPI
   MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
   global_max = local_max;
#endif

   return global_max;
}


Real
SpaceUtils::MaxNorm( const LevelData<FluxBox>& a )
{
   const DisjointBoxLayout& grids = a.disjointBoxLayout();

   double local_max = -DBL_MAX;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         double this_max = a[dit][dir].max();
         if (this_max > local_max) local_max = this_max;
      }
   }

   double global_max;
#ifdef CH_MPI
   MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
   global_max = local_max;
#endif

   return global_max;
}


void
SpaceUtils::dotProduct(LevelData<FArrayBox>& a_dotProduct,
                       const LevelData<FArrayBox>& a_data_1,
                       const LevelData<FArrayBox>& a_data_2)
{
   CH_assert( a_data_1.nComp()==a_data_2.nComp() );

   const DisjointBoxLayout& grids = a_dotProduct.getBoxes();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     
      a_dotProduct[dit].setVal(0.);
      FArrayBox tmp(a_dotProduct[dit].box(), 1);
      for (int n=0; n<a_data_1.nComp(); ++n) {
         tmp.copy(a_data_1[dit],n, 0, 1);
         tmp.mult(a_data_2[dit],n, 0, 1);
         a_dotProduct[dit].plus(tmp);
      }
   }
}

Real
SpaceUtils::dotProduct( const LevelData<FArrayBox>&  a_1,
                        const LevelData<FArrayBox>&  a_2 )
{
   const DisjointBoxLayout& grids = a_1.disjointBoxLayout();

   double local_sum = 0.;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      local_sum += a_1[dit].dotProduct(a_2[dit],grids[dit]);
   }

   double global_sum;
#ifdef CH_MPI
   MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
   global_sum = local_sum;
#endif

   return global_sum;
}

void
SpaceUtils::applyHarmonicFiltering(LevelData<FArrayBox>& a_phi,
                                   const int& a_dir)
{
   CH_TIME("SpaceUtils::::applyHarmonicFiltering");
  
   const DisjointBoxLayout& grids = a_phi.getBoxes();
   
   LevelData<FArrayBox> tmp;
   tmp.define(a_phi);
   
   Box domain_box = grids.physDomain().domainBox();
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      int Npts = grids[dit].size(a_dir);
      
      if (Npts != domain_box.size(a_dir)) {
         MayDay::Error("SpaceUtils::applyHarmonicFiltering(): the direction of filtering must be single-block and has no domain decomposition");
      }
      
      for (BoxIterator bit( grids[dit] ); bit.ok(); ++bit) {
         IntVect iv( bit() );
         
         double A = 0.0;
         double B = 0.0;
         for (int n=0; n<Npts; ++n) {
            IntVect ivSum(iv);
            ivSum[a_dir] = n;
            double phase = 2.0 * Pi * (n+0.5)/Npts;
            A += (2.0/Npts) * tmp[dit](ivSum,0) * cos(phase);
            B += (2.0/Npts) * tmp[dit](ivSum,0) * sin(phase);
         }
         double phase = 2.0 * Pi * (iv[a_dir]+0.5)/Npts;
         a_phi[dit](iv,0) = A*cos(phase)+B*sin(phase);
      }
   }
}

void
SpaceUtils::enforcePositivity(LevelData<FArrayBox>& a_phi)
{
   // Replaces all negative values with zeroes
  
   CH_TIME("SpaceUtils::::enforcePositivity");
  
   const DisjointBoxLayout& grids = a_phi.getBoxes();
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      FORT_ENFORCE_POSITIVITY(CHF_FRA( a_phi[dit] ),
                              CHF_BOX( a_phi[dit].box() ));
   }
}


unsigned long long int
SpaceUtils::getLevelDataLocalSize( const LevelData<FArrayBox>& a_data )
{
   // Returns the number of bytes occupied by the data on the calling processor
   unsigned long long int local = 0;

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      local += a_data[dit].box().numPts();
   }

   return local * a_data.nComp() * sizeof(Real);
}

unsigned long long int
SpaceUtils::getLevelDataLocalSize( const LevelData<FluxBox>& a_data )
{
   // Returns the number of bytes occupied by the data on the calling processor
   unsigned long long int local = 0;

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         local += a_data[dit][dir].box().numPts();
      }
   }

   return local * a_data.nComp() * sizeof(Real);
}

unsigned long long int
SpaceUtils::getLevelDataLocalSize( const LevelData<EdgeDataBox>& a_data )
{
   // Returns the number of bytes occupied by the data on the calling processor
   unsigned long long int local = 0;

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         local += a_data[dit][dir].box().numPts();
      }
   }

   return local * a_data.nComp() * sizeof(Real);
}

unsigned long long int
SpaceUtils::getLevelDataLocalSize( const LevelData<NodeFArrayBox>& a_data )
{
   // Returns the number of bytes occupied by the data on the calling processor
   unsigned long long int local = 0;

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      local += a_data[dit].getFab().box().numPts();
   }

   return local * a_data.nComp() * sizeof(Real);
}

unsigned long long int
SpaceUtils::getBoxLayoutDataLocalSize( const BoxLayoutData<FArrayBox>& a_data )
{
   // Returns the number of bytes occupied by the data on the calling processor
   unsigned long long int local = 0;

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      local += a_data[dit].box().numPts();
   }

   return local * a_data.nComp() * sizeof(Real);
}

unsigned long long int
SpaceUtils::getBoxLayoutLocalSize( const BoxLayout& a_bl )
{
   // Returns the number of cells occupied by the data on the calling processor
   unsigned long long int local = 0;

   if ( a_bl.size() > 0 ) {
      for (DataIterator dit(a_bl); dit.ok(); ++dit) {
         local += a_bl[dit].numPts();
      }
   }

   return local;
}

void
SpaceUtils::getBoxLayoutDataMaxSize(const std::string&               a_name,
                                    const BoxLayoutData<FArrayBox>&  a_data )
{
   // Returns the maximum number of bytes occupied by the data per processor
   unsigned long long int local_bytes;
   local_bytes = SpaceUtils::getBoxLayoutDataLocalSize(a_data);
   
   unsigned long long int max_bytes;
#ifdef CH_MPI
   MPI_Allreduce(&local_bytes, &max_bytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
#else
   max_bytes = local_bytes;
#endif
   if ( procID() == 0 ) {
      cout << a_name << " maximum per process memory = " << max_bytes << " bytes " << endl;
   }
}

void
SpaceUtils::getLevelDataMaxSize(const std::string&          a_name,
                                const LevelData<FArrayBox>& a_data )
{
   // Returns the maximum number of bytes occupied by the data per processor
   unsigned long long int local_bytes;
   local_bytes = SpaceUtils::getLevelDataLocalSize(a_data);
   
   unsigned long long int max_bytes;
#ifdef CH_MPI
   MPI_Allreduce(&local_bytes, &max_bytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
#else
   max_bytes = local_bytes;
#endif
   if ( procID() == 0 ) {
      cout << a_name << " maximum per process memory = " << max_bytes << " bytes " << endl;
   }
}

void
SpaceUtils::getFArrayBoxMaxSize(const std::string& a_name,
                                const FArrayBox&   a_data )
{
   // Returns the maximum number of bytes occupied by the data per processor
   
   unsigned long long int local_bytes;
   local_bytes = a_data.box().numPts() * a_data.nComp() * sizeof(Real);
   
   unsigned long long int max_bytes;
#ifdef CH_MPI
   MPI_Allreduce(&local_bytes, &max_bytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
#else
   max_bytes = local_bytes;
#endif
   if ( procID() == 0 ) {
      cout << a_name << " maximum per process memory = " << max_bytes << " bytes " << endl;
   }
}

#include "NamespaceFooter.H"
