#include "altFaceAverages.H"
#include "altFaceAveragesF_F.H"

#include "NamespaceHeader.H"


AltFaceAverages::AltFaceAverages(const PhaseGeom& a_phase_geometry)
{
   IntVect ghostVect = (a_phase_geometry.secondOrder()) ? IntVect::Zero : IntVect::Unit;
   m_normal_vel.define( a_phase_geometry.gridsFull(), 1, ghostVect );
}

AltFaceAverages::~AltFaceAverages()
{
}

void
AltFaceAverages::computeFaceAverages(LevelData<FluxBox>&         a_face_phi,
                                     const LevelData<FArrayBox>& a_cell_phi,
                                     const LevelData<FluxBox>&   a_face_vel,
                                     const PhaseGeom&            a_geom,
                                     const std::string&          a_method,
                                     const bool                  a_mapped_velocity )
{
   CH_TIME("AltFaceAverages::computeFaceAverages");
   
   CH_assert( a_cell_phi.ghostVect()>=IntVect::Unit );

   // In order to do upwinding, need normal velocities in
   // computational space.  We just need a sign, so we don't need fourth
   // order for this.
   const DisjointBoxLayout& grids( a_face_phi.getBoxes() );
   if ( a_mapped_velocity ) {
      // The passed in velocity is in the mapped frame
      int nvcomp = a_face_vel.nComp();
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const PhaseBlockCoordSys& block_coord_sys = a_geom.getBlockCoordSys(grids[dit]);
         RealVect face_area = block_coord_sys.getMappedFaceArea();

         const FluxBox& thisFaceVel( a_face_vel[dit] );
         FluxBox& thisNormalVel( m_normal_vel[dit] );

         for (int dir=0; dir<SpaceDim; ++dir) {
            int src_comp = (nvcomp==1)? 0: dir;
            FArrayBox& thisNormalVel_dir( thisNormalVel[dir] );
            thisNormalVel_dir.copy(thisFaceVel[dir],src_comp,0,1);
            thisNormalVel_dir *= face_area[dir];
         }
      }
   }
   else {
      // The passed in velocity is in the physical frame (default)
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_normal_vel[dit].setVal( 0.0 );
      }
      a_geom.computeMetricTermProductAverage( m_normal_vel, a_face_vel, false );
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FArrayBox& this_cell_phi( a_cell_phi[dit] );
      const FluxBox& this_normal_vel( m_normal_vel[dit] );
      FluxBox& this_face_phi( a_face_phi[dit] );

      for (int dir(0); dir<SpaceDim; dir++) {
         Box face_box( grids[dit] );
         // for 4th order, we need an extra face in the mapped-grid,
         // transverse directions to handle 4th-order products
         if (!a_geom.secondOrder()) {
            for (int tdir(0); tdir<SpaceDim; tdir++) {
               if (tdir!=dir) {
                  const int TRANSVERSE_GROW(1);
                  face_box.grow( tdir, TRANSVERSE_GROW );
               }
            }
         }
	 
	 // This will turn cell-centered box into
	 //face-centered box in the direction dir
         face_box.surroundingNodes( dir );

         // now compute limited face value
         FArrayBox& this_face_phi_dir( this_face_phi[dir] );
         const FArrayBox& this_normal_vel_dir( this_normal_vel[dir] );
         
         if (a_method == "bweno") {
            FORT_BWENOFACEVALUES(CHF_FRA( this_face_phi_dir ),
                                 CHF_CONST_FRA( this_cell_phi ),
                                 CHF_CONST_FRA1( this_normal_vel_dir, 0 ),
                                 CHF_BOX( face_box ),
                                 CHF_CONST_INT( dir ) );
         }
         else if (a_method == "uw5") {
            FORT_UW5FACEVALUES(CHF_FRA( this_face_phi_dir ),
                               CHF_CONST_FRA( this_cell_phi ),
                               CHF_CONST_FRA1( this_normal_vel_dir, 0 ),
                               CHF_BOX( face_box ),
                               CHF_CONST_INT( dir ) );
         
         }
         else if (a_method == "uw3") {
            FORT_UW3FACEVALUES(CHF_FRA( this_face_phi_dir ),
                               CHF_CONST_FRA( this_cell_phi ),
                               CHF_CONST_FRA1( this_normal_vel_dir, 0 ),
                               CHF_BOX( face_box ),
                               CHF_CONST_INT( dir ) );
         }
         else if (a_method == "uw1") {
            FORT_UW1FACEVALUES(CHF_FRA( this_face_phi_dir ),
                               CHF_CONST_FRA( this_cell_phi ),
                               CHF_CONST_FRA1( this_normal_vel_dir, 0 ),
                               CHF_BOX( face_box ),
                               CHF_CONST_INT( dir ) );
         }
         else if (a_method == "c2") {
            FORT_C2FACEVALUES(CHF_FRA( this_face_phi_dir ),
                               CHF_CONST_FRA( this_cell_phi ),
                               CHF_BOX( face_box ),
                               CHF_CONST_INT( dir ) );
         }
         else if (a_method == "weno5") {
            FORT_WENO5FACEVALUES(CHF_FRA( this_face_phi_dir ),
                               CHF_CONST_FRA( this_cell_phi ),
                               CHF_CONST_FRA1( this_normal_vel_dir, 0 ),
                               CHF_BOX( face_box ),
                               CHF_CONST_INT( dir ) );
         }
         else {
            const std::string msg( "altFaceAverages: unknown method ");
            MayDay::Error( msg.c_str() );
         }
      } // end loop over directions
   } // end loop over grids

   if ( !a_geom.secondOrder() ) {
      a_face_phi.exchange();
   }
}

#include "NamespaceFooter.H"
