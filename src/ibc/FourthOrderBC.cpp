#include "FourthOrderBC.H"
#include "RealVect.H"
#include "Directions.H"

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
   //   if ( a_dir == RADIAL_DIR ) {

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
//   }
}


void FourthOrderBC::setInflowOutflowBC( FArrayBox&               a_this_soln,
                                        const Box&               a_boundary_box,
                                        const FArrayBox&         a_inflow_data,
                                        const FluxBox&           a_this_face_vel,
                                        const int&               a_dir,
                                        const Side::LoHiSide&    a_side )
{
   CH_assert( a_this_soln.nComp()==a_inflow_data.nComp() );
   CH_assert( a_inflow_data.box().contains( a_boundary_box ) );

   // inflow
   a_this_soln.copy( a_inflow_data, a_boundary_box );

   // outflow
   FourthOrderBC::setOutflowBC( a_this_soln,
                                a_boundary_box,
                                a_this_face_vel,
                                a_dir,
                                a_side );
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
