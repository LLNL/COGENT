#ifndef _ROSENBLUTHPOTENTIALS_H_
#define _ROSENBLUTHPOTENTIALS_H_

#include "ParmParse.H"
#include "PhaseGeom.H"
#include <sstream>

#include "FArrayBox.H"
#include "LevelData.H"

#include "HYPRE_struct_ls.h"
#include "HYPRE_struct_mv.h"

#include "NamespaceHeader.H"


/**
 * RosenbluthPotentials collisions operator class.
*/
class RosenbluthPotentials
{
public:

   /// Constructor with initialization.

   RosenbluthPotentials( LevelData<FArrayBox>&      phi_one,
                         LevelData<FArrayBox>&       phi_two,
                         const LevelData<FArrayBox>& rhs,
                         const PhaseGeom&            phase_geom,
                         const Real                  mass,
                         const Real                  pcg_tol,
                         const Real                  pcg_maxiter,
                         const int                   mult_num,   
                         const int                   verbocity );

   /// Destructor.
   /**
    */
   virtual ~RosenbluthPotentials();

   /// Solves the Laplace equation 
   void solve( LevelData<FArrayBox>&       solution,
               const LevelData<FArrayBox>& rhs ) const;

   ///Computes multipole coefficients
   void computeMultipoleCoeff( CFG::LevelData<CFG::FArrayBox>&  multipole_coeff,
                               const LevelData<FArrayBox>&      rho ) const;

   ///Imposes multipole BC for the Poisson eq-n for phi_one
   void imposeBC1( LevelData<FArrayBox>&       rhsBC,
                   const LevelData<FArrayBox>& rho ) const;

   ///Imposes multipole BC for the Poisson eq-n for phi_two
   /* Here, the situation is more comples, since solution for phi_one
      does not rapidly decay inside the computational domain. Therefore, 
      we fisrt separate the inner part of phi_one, and use it for multipole 
      solution as in imposeBC1. For the outer part of phi_one and its 
      multipole contributions we use the analytic expression in the 
      entire velocity space 
   */
   void imposeBC2( LevelData<FArrayBox>&       rhsBC,
                   const LevelData<FArrayBox>& phi_one,
                   const LevelData<FArrayBox>& rho ) const;

   ///Hypre routines

   void createHypreData( const MPI_Comm&       comm,
                         const List<VEL::Box>& boxes,
                         HYPRE_StructGrid&     grid,
                         HYPRE_StructStencil&  stencil,
                         HYPRE_StructMatrix&   matrix,
                         HYPRE_StructVector&   x,
                         HYPRE_StructVector&   b,
                         HYPRE_StructSolver&   solver ) const;

   void destroyHypreData( HYPRE_StructSolver&  solver,
                          HYPRE_StructVector&  b,
                          HYPRE_StructVector&  x,
                          HYPRE_StructMatrix&  matrix,
                          HYPRE_StructStencil& stencil,
                          HYPRE_StructGrid&    grid ) const;

   void constructMatrix( const VEL::ProblemDomain& domain,
                         const VEL::RealVect&      dx,
                         const List<VEL::Box>&     boxes,
                         const double&             Bfield,
                         HYPRE_StructMatrix&       matrix ) const;

   /// Utility function (used to make proper block arrangement)
   void getSliceBoxes( const CFG::IntVect     cfg_iv,
                       const List<VEL::Box>&  vel_boxes,
                       List<Box>&             phase_boxes ) const;

private:
 
   int m_verbosity;
   ParmParse m_ppcls;

   const PhaseGeom& m_phase_geom;
   Real m_mass;
   Real m_pcg_tol;
   Real m_pcg_maxiter;
   int m_mult_num;   

   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "CLS.species" section
    * of ParmParse database, where "species" is replaced by the species name
    *
    * \param[in] ppcls CLS.species input database.
    */

};


#include "NamespaceFooter.H"

#endif
