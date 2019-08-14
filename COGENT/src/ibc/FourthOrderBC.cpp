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
     
     FORT_FOURTH_ORDER_ODD_BC( CHF_FRA(a_this_soln),
                               CHF_BOX(a_boundary_box),
         	               CHF_CONST_INTVECT(ghosts),
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

     FORT_FOURTH_ORDER_NEUMANN_BC( CHF_FRA(a_this_soln),
                                   CHF_BOX(a_boundary_box),
				   CHF_CONST_INTVECT(ghosts),
                                   CHF_CONST_INT(a_dir),
                                   CHF_CONST_INT(ISIDE) );
   }
   
   else if (a_bc_type == "symmetry") {

     CH_assert( a_this_soln.nComp()==SpaceDim );
     const int ISIDE(a_side);
     IntVect ghosts = a_boundary_box.size();
     
     FluxBox this_face_val;
     this_face_val.define(a_this_face_vel);
     this_face_val.copy(a_this_face_vel);
     this_face_val.setVal(0.);
      
     Box box = a_this_soln.box();
     FArrayBox this_soln_copy(box,1);
     
     for (int n=0; n<a_this_soln.nComp(); ++n) {
     
     if(n==a_dir) {
     
        // apply zero value bc to normal components
        //
        this_soln_copy.copy(a_this_soln,n,0,1);
        FORT_FOURTH_ORDER_ODD_BC( CHF_FRA(this_soln_copy),
                                  CHF_BOX(a_boundary_box),
	              	          CHF_CONST_INTVECT(ghosts),
                                  CHF_CONST_INT(a_dir),
                                  CHF_CONST_INT(ISIDE) );
        a_this_soln.copy(this_soln_copy,0,n,1);
     
     }
     else {
     
        // apply zero gradient bc to tangential components
        //
        this_soln_copy.copy(a_this_soln,n,0,1);
        FORT_FOURTH_ORDER_NEUMANN_BC( CHF_FRA(this_soln_copy),
                                      CHF_BOX(a_boundary_box),
	        	      	      CHF_CONST_INTVECT(ghosts),
                                      CHF_CONST_INT(a_dir),
                                      CHF_CONST_INT(ISIDE) );
        a_this_soln.copy(this_soln_copy,0,n,1);
        
     }
     }
  
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
     
        FORT_FOURTH_ORDER_NEUMANN_BC( CHF_FRA(a_this_soln),
                                      CHF_BOX(a_boundary_box),
	        	      	      CHF_CONST_INTVECT(ghosts),
                                      CHF_CONST_INT(a_dir),
                                      CHF_CONST_INT(ISIDE) );
        
     }
  
   } 

}

void FourthOrderBC::setFluxBc( FluxBox&           a_dst,
                         const Box&               a_fill_box,
                         const std::string&       a_bc_type,
                         const FluxBox&           a_src,
                         const int                a_dir,
                         const Side::LoHiSide&    a_side )
{
   CH_assert( a_dst.nComp()==a_src.nComp() );
   CH_assert( a_src.box().contains( a_fill_box ) );

   // copy at boundary
   a_dst.copy( a_src, a_fill_box );

}
  
void FourthOrderBC::setEdgeBC( EdgeDataBox&     a_dst,
                         const Box&             a_fill_box,
                         const std::string&     a_bc_type,
                         const int              a_dir,
                         const Side::LoHiSide&  a_side )
{
   // apply neumann bc to non-dir comps (tangential to edge)
   //
   IntVect ghosts = a_fill_box.size();
   const int ISIDE(a_side);
   const int GROW(1);
   Box a_fill_box_grown = a_fill_box;
   //cout << "JRA: a_fill_box = " << a_fill_box << endl;
   if (a_bc_type == "neumann_on_edge") {
      a_fill_box_grown.grow(a_dir,GROW);
      
      for (int dir=0; dir<SpaceDim; dir++) {
         if (dir != a_dir) {
            FArrayBox& this_dst_dir(a_dst[dir]);
            FORT_FOURTH_ORDER_NEUMANN_BC_ON_EDGE( CHF_FRA(this_dst_dir),
                                                  CHF_BOX(a_fill_box_grown),
             	              	                  CHF_CONST_INTVECT(ghosts),
                                                  CHF_CONST_INT(a_dir),
                                                  CHF_CONST_INT(ISIDE) );
            //this_dst_dir.setVal(0.0,a_fill_box,0);
         }
      }

   }
   for (int dir=0; dir<SpaceDim; dir++) {

      // adjust box appropriately depending on dir
      // Note that dir data is stag in all other dirs
      //
      if (dir == a_dir) {
         a_fill_box_grown.growDir(a_dir,a_side,1);
         for (int dir0=0; dir0<SpaceDim; dir0++) {
            if (dir0 != a_dir) a_fill_box_grown.growHi(dir0,1);
         }
      }
      else {
         if(ISIDE==1) a_fill_box_grown.shift(a_dir,1);
      }
      FArrayBox& this_dst_dir(a_dst[dir]);
   
      if (a_bc_type == "natural") {
         //cout << "natural for iside = " << ISIDE << endl;
         FORT_FOURTH_ORDER_EXTRAP_BC_ON_EDGE( CHF_FRA(this_dst_dir),
                                              CHF_BOX(a_fill_box_grown),
                                              CHF_CONST_INT(a_dir),
                                              CHF_CONST_INT(ISIDE) );
      }
      if (a_bc_type == "symmetry") {
         int evenodd = 1; // 1 for even
         if(dir == a_dir) evenodd = -1; // -1 for odd
         FORT_FOURTH_ORDER_EVENODD_BC_ON_EDGE( CHF_FRA(this_dst_dir),
                                               CHF_BOX(a_fill_box_grown),
                                               CHF_CONST_INT(evenodd),
                                               CHF_CONST_INT(a_dir),
                                               CHF_CONST_INT(ISIDE) );
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
   FORT_FOURTH_ORDER_NEUMANN_BC( CHF_FRA(a_this_soln),
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
