#include "FourthOrderBC.H"
#include "RealVect.H"
#include "Directions.H"
#include "inspect.H"

#include "FourthOrderBCF_F.H"

#include "NamespaceHeader.H"

void computeNormalVelocity( FArrayBox&      a_normalFaceVel,
                            const FluxBox&  a_faceVel,
                            const int&      a_dir)
{
   CH_TIMERS("computeNormalVelocity");
   const int SRC_COMP(a_dir);
   const int DST_COMP(0);
   const int NUM_COMP(1);
   a_normalFaceVel.copy( a_faceVel[a_dir], SRC_COMP, DST_COMP, NUM_COMP );
}


void computeNormalVelocity( FArrayBox&      a_normalFaceVel,
                            const FluxBox&  a_faceVel,
                            const RealVect& a_faceArea,
                            const int&      a_dir)
{
   computeNormalVelocity( a_normalFaceVel, a_faceVel, a_dir);
   a_normalFaceVel *= a_faceArea[a_dir];
}


void FourthOrderBC::setOutflowBC( FArrayBox&               a_this_soln,
                                  const Box&               a_boundary_box,
                                  const FluxBox&           a_this_face_vel,
                                  const int&               a_dir,
                                  const Side::LoHiSide&    a_side,
				  const int&               a_extrapolated_bc)
{
   CH_TIMERS("FourthOrderBC::setOutflowBC");

   // assumes the inflow data is already set
   // define the face box for this boundary
   Box face_box( a_boundary_box );
   face_box.grow( a_dir, 1 );
   const Box& this_box( a_this_soln.box() );
   face_box &= this_box;
   face_box.shiftHalf( a_dir, sign( a_side ) );

   //outflow
   FArrayBox face_velocity( face_box, 1 );
   computeNormalVelocity( face_velocity,
                          a_this_face_vel,
                          a_dir );

   const int ISIDE(a_side);
   FORT_FOURTH_ORDER_OUTFLOW_BC( CHF_FRA(a_this_soln),
                                 CHF_BOX(a_boundary_box),
                                 CHF_CONST_FRA(face_velocity),
                                 CHF_CONST_INT(a_dir),
                                 CHF_CONST_INT(ISIDE),
                                 CHF_CONST_INT(a_extrapolated_bc));
}


void FourthOrderBC::setInflowOutflowBC( FArrayBox&               a_this_soln,
                                        const Box&               a_boundary_box,
                                        const FArrayBox&         a_inflow_data,
                                        const std::string&       a_bc_type,
                                        const FluxBox&           a_this_face_vel,
                                        const int&               a_dir,
                                        const Side::LoHiSide&    a_side )
{

   CH_TIMERS("FourthOrderBC::setInflowOutflowBC");
   CH_TIMER("copy_inflow_func", t_copy_inflow_func);

   CH_assert( a_this_soln.nComp()==a_inflow_data.nComp() );
   CH_assert( a_inflow_data.box().contains( a_boundary_box ) );

   // inflow
   CH_START(t_copy_inflow_func);
   a_this_soln.copy( a_inflow_data, a_boundary_box );
   CH_STOP(t_copy_inflow_func);

   // outflow
   if ((a_bc_type != "dirichlet") &&
       (a_bc_type != "neumann")   &&
       (a_bc_type != "odd")   &&
       (a_bc_type != "axis")   &&
       (a_bc_type != "symmetry")   &&
       (a_bc_type != "fixed_ghosts") ) {

      int extrapolated_bc(0);
      if (a_bc_type == "extrapolate") extrapolated_bc = 1;
      
      FourthOrderBC::setOutflowBC( a_this_soln,
                                   a_boundary_box,
                                   a_this_face_vel,
                                   a_dir,
                                   a_side,
                                   extrapolated_bc);
      
   }

   else if (a_bc_type == "dirichlet") {

     const int ISIDE(a_side);
     FORT_SECOND_ORDER_DIRICHLET_BC(CHF_FRA(a_this_soln),
                                    CHF_BOX(a_boundary_box),
                                    CHF_CONST_FRA(a_inflow_data),
                                    CHF_CONST_INT(a_dir),
                                    CHF_CONST_INT(ISIDE) );
   } 
   
   else if (a_bc_type == "odd") {

     const int ISIDE(a_side);
     IntVect ghosts = a_boundary_box.size();
     
     int evenodd = -1; // -1 for odd
     FORT_FOURTH_ORDER_EVENODD_BC( CHF_FRA(a_this_soln),
                               CHF_BOX(a_boundary_box),
                               CHF_CONST_INT(evenodd),
                               CHF_CONST_INT(a_dir),
                               CHF_CONST_INT(ISIDE) );
   }
   
   else if (a_bc_type == "axis") { // special BC for Btheta on axis

     const int ISIDE(a_side);
     IntVect ghosts = a_boundary_box.size();
     
     FORT_FOURTH_ORDER_AXIS_BC( CHF_FRA(a_this_soln),
                                CHF_BOX(a_boundary_box),
         	                CHF_CONST_INTVECT(ghosts),
                                CHF_CONST_INT(a_dir),
                                CHF_CONST_INT(ISIDE) );
   }

   else if (a_bc_type == "neumann") {

     const int ISIDE(a_side);
     IntVect ghosts = a_boundary_box.size();

     FORT_FOURTH_ORDER_EVEN_BC( CHF_FRA(a_this_soln),
                                CHF_BOX(a_boundary_box),
				CHF_CONST_INTVECT(ghosts),
                                CHF_CONST_INT(a_dir),
                                CHF_CONST_INT(ISIDE) );
   }
   
   else if (a_bc_type == "symmetry") {

     CH_assert( a_this_soln.nComp()==SpaceDim );
     const int ISIDE(a_side);
     IntVect ghosts = a_boundary_box.size();
     
     FORT_FOURTH_ORDER_SYMMETRY_BC( CHF_FRA(a_this_soln),
                                    CHF_BOX(a_boundary_box),
	              	            CHF_CONST_INTVECT(ghosts),
                                    CHF_CONST_INT(a_dir),
                                    CHF_CONST_INT(ISIDE) );

   } 
}

void FourthOrderBC::setBC( FArrayBox&       a_this_soln,
                     const Box&             a_boundary_box,
                     const std::string&     a_bc_type,
                     const int&             a_dir,
                     const Side::LoHiSide&  a_side )
{
   CH_TIMERS("FourthOrderBC::setBC()");
   
   if (a_bc_type == "symmetry" || a_bc_type == "natural" ) {

     const int ISIDE(a_side);
     IntVect ghosts = a_boundary_box.size();
     
     for (int n=0; n<a_this_soln.nComp(); ++n) {
     
        FORT_FOURTH_ORDER_EVEN_BC( CHF_FRA(a_this_soln),
                                   CHF_BOX(a_boundary_box),
	                	   CHF_CONST_INTVECT(ghosts),
                                   CHF_CONST_INT(a_dir),
                                   CHF_CONST_INT(ISIDE) );
        
     }
  
   } 

}

void FourthOrderBC::setFluxBC( FluxBox&           a_dst,
                         const Box&               a_fill_box,
                         const std::string&       a_bc_type,
                         const int                a_dir,
                         const Side::LoHiSide&    a_side )
{
   CH_TIMERS("FourthOrderBC::setFluxBC()");

   const int ISIDE(a_side);
   for (int dir=0; dir<SpaceDim; dir++) {
      
      // adjust box appropriately depending on dir
      // Note that dir data is stag in dir direction for FluxBox
      //
      FArrayBox& this_dst_dir(a_dst[dir]);
      const Box& this_dst_dir_box( this_dst_dir.box() );
      Box a_fill_box_grown = this_dst_dir_box;
      if(ISIDE==0) a_fill_box_grown.setBig(a_dir,a_fill_box.bigEnd(a_dir));
      if(ISIDE==1) a_fill_box_grown.setSmall(a_dir,a_fill_box.smallEnd(a_dir));
      if(dir==a_dir && ISIDE==1) a_fill_box_grown.growLo(a_dir,-1); // dont change value on boundary

      int evenodd = 1;
      if (a_bc_type == "symmetry") {
         if(dir == a_dir) {
            evenodd = -1; // -1 for odd
            FORT_FOURTH_ORDER_EVENODD_BC_STAG( CHF_FRA(this_dst_dir),
                                               CHF_BOX(a_fill_box_grown),
                                               CHF_CONST_INT(evenodd),
                                               CHF_CONST_INT(a_dir),
                                               CHF_CONST_INT(ISIDE) );
         } else {
            evenodd = 1; // 1 for even
            FORT_FOURTH_ORDER_EVENODD_BC( CHF_FRA(this_dst_dir),
                                          CHF_BOX(a_fill_box_grown),
                                          CHF_CONST_INT(evenodd),
                                          CHF_CONST_INT(a_dir),
                                          CHF_CONST_INT(ISIDE) );
         }
      }
      if (a_bc_type == "wall") {
         if (dir == a_dir) {
            int order = 4;
            FORT_FOURTH_ORDER_EXTRAP_BC( CHF_FRA(this_dst_dir),
                                         CHF_BOX(a_fill_box_grown),
                                         CHF_CONST_INT(a_dir),
                                         CHF_CONST_INT(ISIDE),
                                         CHF_CONST_INT(order) );
         } 
         else {
            FORT_FOURTH_ORDER_NEUMANN_BC( CHF_FRA(this_dst_dir),
                                          CHF_BOX(a_fill_box_grown),
                                          CHF_CONST_INT(a_dir),
                                          CHF_CONST_INT(ISIDE) );
         }
      }
      
   }

}
  
void FourthOrderBC::setEdgeBC( EdgeDataBox&     a_dst,
                         const Box&             a_fill_box,
                         const std::string&     a_bc_type,
                         const int              a_dir,
                         const Side::LoHiSide&  a_side )
{
   CH_TIMERS("FourthOrderBC::setEdgeBC()");

   const int ISIDE(a_side);
   for (int dir=0; dir<SpaceDim; dir++) {
      
      // adjust box appropriately depending on dir
      // Note that dir data is stag in all other dirs for EdgeDataBox
      // (compliment of FluxBox)
      //
      FArrayBox& this_dst_dir(a_dst[dir]);
      const Box& this_dst_dir_box( this_dst_dir.box() );
      Box a_fill_box_grown = this_dst_dir_box;
      if(ISIDE==0) a_fill_box_grown.setBig(a_dir,a_fill_box.bigEnd(a_dir));
      if(ISIDE==1) a_fill_box_grown.setSmall(a_dir,a_fill_box.smallEnd(a_dir));
      if(dir!=a_dir && ISIDE==1) a_fill_box_grown.growLo(a_dir,-1); // dont change value on boundary

      if (a_bc_type == "natural") {
         int order = 2;
         FORT_FOURTH_ORDER_EXTRAP_BC( CHF_FRA(this_dst_dir),
                                      CHF_BOX(a_fill_box_grown),
                                      CHF_CONST_INT(a_dir),
                                      CHF_CONST_INT(ISIDE),
                                      CHF_CONST_INT(order) );
      }
      if (a_bc_type == "symmetry") {
         int evenodd;
         if(dir == a_dir) {
            evenodd = -1; // -1 for odd
            FORT_FOURTH_ORDER_EVENODD_BC( CHF_FRA(this_dst_dir),
                                          CHF_BOX(a_fill_box_grown),
                                          CHF_CONST_INT(evenodd),
                                          CHF_CONST_INT(a_dir),
                                          CHF_CONST_INT(ISIDE) );
         } 
         else {
            evenodd = 1; // 1 for even
            FORT_FOURTH_ORDER_EVENODD_BC_STAG( CHF_FRA(this_dst_dir),
                                               CHF_BOX(a_fill_box_grown),
                                               CHF_CONST_INT(evenodd),
                                               CHF_CONST_INT(a_dir),
                                               CHF_CONST_INT(ISIDE) );
         }
      }
      if (a_bc_type == "neumann_on_edge") {
         if (dir != a_dir) {
            FORT_FOURTH_ORDER_NEUMANN_BC_ON_EDGE( CHF_FRA(this_dst_dir),
                                                  CHF_BOX(a_fill_box_grown),
                                                  CHF_CONST_INT(a_dir),
                                                  CHF_CONST_INT(ISIDE) );
         }
      }
      
   }

}


void FourthOrderBC::setNeumannBC( FArrayBox&               a_this_soln,
                                  const IntVect&           a_ghosts,
                                  const Box&               a_boundary_box,
                                  const int&               a_dir,
                                  const Side::LoHiSide&    a_side )
{
   const int ISIDE(a_side);
   FORT_FOURTH_ORDER_EVEN_BC( CHF_FRA(a_this_soln),
                              CHF_BOX(a_boundary_box),
                              CHF_CONST_INTVECT(a_ghosts),
                              CHF_CONST_INT(a_dir),
                              CHF_CONST_INT(ISIDE) );
}

void FourthOrderBC::setDirichletBC( FArrayBox&               a_this_soln,
                                    const IntVect&           a_ghosts,
                                    const Box&               a_boundary_box,
                                    const FluxBox&           a_face_values,
                                    const int&               a_dir,
                                    const Side::LoHiSide&    a_side )

{
   const int ISIDE(a_side);
   FORT_FOURTH_ORDER_DIRICHLET_BC( CHF_FRA(a_this_soln),
                                   CHF_BOX(a_boundary_box),
                                   CHF_CONST_INTVECT(a_ghosts),
                                   CHF_CONST_FRA(a_face_values[a_dir]),
                                   CHF_CONST_INT(a_dir),
                                   CHF_CONST_INT(ISIDE) );
}

#include "NamespaceFooter.H"
