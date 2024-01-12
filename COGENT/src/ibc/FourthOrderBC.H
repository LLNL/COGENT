//#ifndef  _FOURTH_ORDER_BC_H_
//#define  _FOURTH_ORDER_BC_H_

#include "IntVect.H"
#include "Box.H"
#include "LoHiSide.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"

#include "NamespaceHeader.H"

namespace FourthOrderBC {

   void setInflowBC( FArrayBox&       a_this_soln,
               const Box&             a_boundary_box,
               const FArrayBox&       a_inflow_data,
               const std::string&     a_this_bc_type,
               const FluxBox&         a_this_face_vel,
               const int&             a_dir,
               const Side::LoHiSide&  a_side );
   
   void setInflowOutflowBC( FArrayBox&       a_this_soln,
                      const Box&             a_boundary_box,
                      const FArrayBox&       a_inflow_data,
                      const std::string&     a_this_bc_type,
                      const FluxBox&         a_this_face_vel,
                      const int&             a_dir,
                      const Side::LoHiSide&  a_side );
   
   void setOutflowBC( FArrayBox&       a_this_soln,
                const Box&             a_boundary_box,
                const FluxBox&         a_this_face_vel,
                const int&             a_dir,
                const Side::LoHiSide&  a_side,
                const int&             a_extrapolated_bc );

   void setCellBC( FArrayBox&       a_this_soln,
             const FluxBox&         a_boundary_vals,
             const Box&             a_boundary_box,
             const std::string&     a_bc_type,
             const int&             a_dir,
             const Side::LoHiSide&  a_side );
   
   void setInsulatorConductorBC( FArrayBox&  a_Bv,
                           const FArrayBox&  a_ICBinary, 
                           const Real        a_ft,
                           const Box&        a_boundary_box,
                           const int         a_dir,
                           const Side::LoHiSide&  a_side );

   void setInsulatorConductorBC( FArrayBox&  a_Bdir,
                           const FArrayBox&  a_ICbinary, 
                           const int         a_dir,
                           const Real        a_ft,
                           const Box&        a_bdry_box,
                           const int         a_bdry_dir,
                           const Side::LoHiSide&  a_bdry_side );

   void setInsulatorConductorEdgeBC( EdgeDataBox&     a_dst,
                               const EdgeDataBox&     a_ICbinary_ce,
                               const Box&             a_fill_box,
                               const std::string&     a_bc_type,
                               const int              a_dir,
                               const Side::LoHiSide&  a_side );
   
   void setFluxBC( FluxBox&           a_dst,
             const Box&               a_fill_box,
             const std::string&       a_bc_type,
             const int                a_dir,
             const Side::LoHiSide&    a_side );

   void setEdgeBC( EdgeDataBox&       a_dst,
             const Box&               a_fill_box,
             const std::string&       a_bc_type,
             const int                a_dir,
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
