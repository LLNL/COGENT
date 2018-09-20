#include "InOutFlowBC.H"
#include "InOutFlowBCF_F.H"
#include "MultiBlockCoordSys.H"
#include "MagBlockCoordSys.H"

#include "NamespaceHeader.H"

enum directions {RADIAL_DIR, POLOIDAL_DIR, TOROIDAL_DIR};

InOutFlowBC::InOutFlowBC()
   : m_face_vel(NULL)
{
}


void InOutFlowBC::ghostCellBC( LevelData<FArrayBox>&      a_soln,
                               const MultiBlockLevelGeom& a_geometry,
                               const Real                 a_dx,
                               const Real                 time )
{
   if ( m_face_vel==NULL ) {
      MayDay::Error( "InOutFlowBC(): Attempt to set ghost cell values without valide velocity!" );
   }

   const MultiBlockCoordSys* coord_sys( a_geometry.coordSysPtr() );

   for (int block(0); block<coord_sys->numBlocks(); ++block) {
      const MagBlockCoordSys* block_coord_sys(
         static_cast<const MagBlockCoordSys*>(
            coord_sys->getCoordSys( block ) ) );
      const ProblemDomain& domain( block_coord_sys->domain() );
      RealVect face_area( block_coord_sys->getMappedFaceArea() );
      setRadialBoundaryValues( a_soln, domain, *m_face_vel, face_area );
   }

   clearFaceVelocity();
}


Box InOutFlowBC::getGhostBox( const Box& a_domain_box,
                              const Box& a_box,
                              const int& a_dir,
                              const Side::LoHiSide& a_side,
                              const int& a_num_ghosts ) const
{
   const Box gBox( adjCellBox( a_box, a_dir, a_side, a_num_ghosts ) );
   const Box gDomainBox( adjCellBox( a_domain_box, a_dir, a_side, a_num_ghosts ) );
   return (gBox & gDomainBox);
}


void InOutFlowBC::computeNormalVelocity( FArrayBox&      a_normalFaceVel,
                                         const FluxBox&  a_faceVel,
                                         const RealVect& a_faceArea,
                                         const int&      a_dir)
{
   const int SRC_COMP(a_dir);
   const int DST_COMP(0);
   const int NUM_COMP(1);
   a_normalFaceVel.copy( a_faceVel[a_dir], SRC_COMP, DST_COMP, NUM_COMP );
   a_normalFaceVel *= a_faceArea[a_dir];
}


void InOutFlowBC::setRadialBoundaryValues( LevelData<FArrayBox>&     a_soln,
                                           const ProblemDomain&      a_domain,
                                           const LevelData<FluxBox>& a_face_vel,
                                           const RealVect&           a_face_area )
{
   // If direction dir is periodic, don't do anything
   const int DIR(RADIAL_DIR);
   if (!a_domain.isPeriodic(DIR)) {

      const DisjointBoxLayout& dbl( a_soln.getBoxes() );
      const Box& domain_box( a_domain.domainBox() );

      const IntVect& ghostVect( a_soln.ghostVect() );
      int num_ghosts_dir( ghostVect[DIR] );
      IntVect transverseGhosts( ghostVect );
      transverseGhosts[DIR] = 0;

      DataIterator dit( a_soln.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         FArrayBox& this_soln( a_soln[dit] );
         const Box& this_box( dbl[dit] );
         const FluxBox& this_face_vel( a_face_vel[dit] );

         for (SideIterator si; si.ok(); ++si) {
            Side::LoHiSide side( si() );

            Box boundary_box( getGhostBox( domain_box,
                                           this_box,
                                           DIR,
                                           side,
                                           num_ghosts_dir) );

            if (!boundary_box.isEmpty()) {
//             boundary_box.grow( transverseGhosts ); // FIX: fills corners?

               // inflow
               this_soln.setVal( -1.0, boundary_box, 0 );

               //outflow
               Box vel_box( boundary_box );
               vel_box.grow( DIR, 1 );
               vel_box &= this_box;
               vel_box.shiftHalf( DIR, sign( side ) );
               FArrayBox this_normal_face_vel( vel_box, 1 );
               computeNormalVelocity( this_normal_face_vel,
                                      this_face_vel,
                                      a_face_area,
                                      DIR );
               const int ISIDE(side);
               FORT_FOURTH_ORDER_OUTFLOW_BC( CHF_FRA(this_soln),
                                             CHF_BOX(boundary_box),
                                             CHF_CONST_INTVECT(ghostVect),
                                             CHF_CONST_FRA(this_normal_face_vel),
                                             CHF_CONST_INT(DIR),
                                             CHF_CONST_INT(ISIDE) );
            }
         }
      }
   }
}

#include "NamespaceFooter.H"
