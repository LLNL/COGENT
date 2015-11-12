//#ifndef  _FOURTH_ORDER_BC_H_
//#define  _FOURTH_ORDER_BC_H_

#include "IntVect.H"
#include "Box.H"
#include "LoHiSide.H"
#include "FArrayBox.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"

namespace FourthOrderBC {

   void setOutflowBC( FArrayBox&               a_this_soln,
                      const Box&               a_boundary_box,
                      const FluxBox&           a_this_face_vel,
                      const int&               a_dir,
                      const Side::LoHiSide&    a_side );

   void setInflowOutflowBC( FArrayBox&               a_this_soln,
                            const Box&               a_boundary_box,
                            const FArrayBox&         a_inflow_data,
                            const FluxBox&           a_this_face_vel,
                            const int&               a_dir,
                            const Side::LoHiSide&    a_side );

   void setNeumannBC( FArrayBox&               a_this_soln,
                      const IntVect&           a_ghosts,
                      const Box&               a_boundary_box,
                      const int&               a_dir,
                      const Side::LoHiSide&    a_side );

   void setDirichletBC( FArrayBox&               a_this_soln,
                        const IntVect&           a_ghosts,
                        const Box&               a_boundary_box,
                        const FluxBox&           a_face_values,
                        const int&               a_dir,
                        const Side::LoHiSide&    a_side );

}

#include "NamespaceFooter.H"

//#endif
