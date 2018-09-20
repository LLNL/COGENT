#ifndef _MAPPED_ADVECTION_OP_H_
#define _MAPPED_ADVECTION_OP_H_

#include "REAL.H"
#include "ParmParse.H"
#include "FourthOrderCoordSys.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "Vector.H"

#include "NamespaceHeader.H"

/// Generic fourth-order mapped grid advection operator.
/**
 * This operator class is a generic implementation of limited,
 * positivity-preserving, fourth-order, finite volume advection in a mapped
 * coordinate system.
*/
class MappedAdvectionOp
{
   public:

      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       *  @param[in] boxes the box layout.
       *  @param[in] coords the coordinate system mapping.
       *  @param[in] dx the computational mesh size.
       */
      void MappedAdvectionOp( ParmParse                  pp,
                              const DisjointboxLayouts&  boxes,
                              const FourthOrderCoordSys& coords,
                              const Real                 dx );

      /// Destructor.
      /**
       */
      ~MappedAdvectionOp() {;}

      /// Evaluates the RHS of the semi-discrete advection operator.
      /**
       * Evaluates the RHS of the semi-discrete advection operator given
       * the face-averaged physical components of the velocity vector.
       *
       * @param[out] rhs data holder for rhs (without ghost cells).
       * @param[in] soln current solution with which to evaluate RHS.
       * @param[in] velocity the face-averaged velocity vectors.
       * @param[in] time current time.
       */
      void evalRHS( LevelData<FArrayBox>&       rhs,
                    const LevelData<FArrayBox>& soln,
                    const LevelData<FluxBox>&   velocity,
                    const Real                  time );

      /// Computes the maximum velocity.
      /**
       * Given the face-averaged physical components of the velocity vector,
       * computes the maximum cell-centered computational space components
       * of the velocity vector.
       *
       * @param[out] max_vel vector of the maximum velocity components.
       * @param[in] velocity the face-averaged velocity vectors.
       */
      void maximumVelocity( Vector<Real>&             max_vel,
                            const LevelData<FluxBox>& velocity );

};

#include "NamespaceFooter.H"

#endif
