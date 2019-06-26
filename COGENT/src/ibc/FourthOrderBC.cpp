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
                                  const Side::LoHiSide&    a_side )
{
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
                                 CHF_CONST_INT(ISIDE) );
}


void FourthOrderBC::setInflowOutflowBC( FArrayBox&               a_this_soln,
                                        const Box&               a_boundary_box,
                                        const FArrayBox&         a_inflow_data,
                                        const std::string&       a_bc_type,
                                        const FluxBox&           a_this_face_vel,
                                        const int&               a_dir,
                                        const Side::LoHiSide&    a_side )
{
   CH_assert( a_this_soln.nComp()==a_inflow_data.nComp() );
   CH_assert( a_inflow_data.box().contains( a_boundary_box ) );

   // inflow
   a_this_soln.copy( a_inflow_data, a_boundary_box );
 
   // outflow
   if ((a_bc_type != "dirichlet") &&
       (a_bc_type != "neumann")   &&
       (a_bc_type != "odd")   &&
       (a_bc_type != "symmetry")   &&
       (a_bc_type != "fixed_ghosts") ) {

      FluxBox this_face_vel;
      this_face_vel.define(a_this_face_vel);
      this_face_vel.copy(a_this_face_vel);
      
      if (a_bc_type == "extrapolate" || a_bc_type == "symmetry") {
         this_face_vel.setVal(0.);
      }
      
      FourthOrderBC::setOutflowBC( a_this_soln,
                                   a_boundary_box,
                                   this_face_vel,
                                   a_dir,
                                   a_side );
      

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
     
        // overwrite 0 component (normal to surface) with dirichlet bc
        //
        //JRA, origally I just set all to NEUMANN, and then applied
        //Dirichlet to n=a_dir. However, this led to the first ghost 
        //cell being set equal to the negative of the second cell
        //inside the domain rather than the first for some reason.
        //Does Neumann_BC change size of FArrayBox?
        //SAME ISSUE when using 4th order dirichlet!!!
        this_soln_copy.copy(a_this_soln,n,0,1);
        //FORT_SECOND_ORDER_DIRICHLET_BC( CHF_FRA(this_soln_copy),
	//            	                CHF_BOX(a_boundary_box),
	//   			        CHF_CONST_INT(a_dir),
	//     			        CHF_CONST_INT(ISIDE) );
        FORT_FOURTH_ORDER_ODD_BC( CHF_FRA(this_soln_copy),
                                  CHF_BOX(a_boundary_box),
	              	          CHF_CONST_INTVECT(ghosts),
                                  CHF_CONST_INT(a_dir),
                                  CHF_CONST_INT(ISIDE) );
        //FORT_FOURTH_ORDER_DIRICHLET_BC( CHF_FRA(this_soln_copy),
        //                                CHF_BOX(a_boundary_box),
        //                                CHF_CONST_INTVECT(ghosts),
        //                                CHF_CONST_FRA(this_face_val[a_dir]),
        //                                CHF_CONST_INT(a_dir),
        //                                CHF_CONST_INT(ISIDE) );
        a_this_soln.copy(this_soln_copy,0,n,1);
     
     }
     else {
     
        // apply neumann bc to tangential components
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

   /*
   if ( a_bc_type=="dirichlet" ) {
      //FluxBox this_face_vel;
      //this_face_vel.define(a_this_face_vel);
      //this_face_vel.copy(a_this_face_vel);
      FArrayBox& this_dst_dir = a_dist[a_dir];
      this_dst_dir.setVal(0.0, a_fill_box, 0);
   }
   */

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
