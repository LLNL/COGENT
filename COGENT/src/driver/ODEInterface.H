#ifndef _ODE_INTEGRATOR_H_
#define _ODE_INTEGRATOR_H_

#include "REAL.H"

#include "NamespaceHeader.H"

/// ODE interface for RK4Integrator.
/**
 * Abstract base class that defines the interface needed for the Runge-Kutta
 * integrator class RK4Integrator.
 */
template <class TSoln>
class ODEInterface
{
   public:

      /// Virtual destructor.
      /**
       */
      virtual ~ODEInterface() {;}

      /// Define TSoln data holder for data.
      /**
       * Defines the TSoln data holder for data without ghost cells.
       *
       * @param[out] rhs data holder for rhs (without ghost
       * cells).
       * @param[in] prototype prototype data holder.
       */
      virtual void defineRHSData( TSoln& rhs, const TSoln& prototype ) = 0;

      /// Defines the TSoln data holder for data with ghost cells.
      /**
       * Defines the TSoln data holder for data with ghost cells.
       *
       * @param[out] soln data holder for solution with ghost cells.
       * @param[in] prototype prototype data holder.
       */
      virtual void defineSolnData( TSoln& soln, const TSoln& prototype ) = 0;

      /// Copies all data from srcSoln to dstSoln.
      /**
       * Copies all data from srcSoln to dstSoln.
       *
       * @param[out] dstSoln destination data holder.
       * @param[in] srcSoln source data holder.
       */
      virtual void copySolnData( TSoln& dstSoln, const TSoln& srcSoln ) = 0;

      /// Evaluates the RHS of the ODE.
      /**
       * Given the current time and solution, evaluates the RHS of the ODE.
       *
       * @param[out] rhs data holder for rhs (without ghost cells).
       * @param[in] soln current solution with which to evaluate RHS.
       * @param[in] step_number current step number
       * @param[in] time current time.
       * @param[in] stage stage number
       */
      virtual void evalRHS( TSoln&       rhs,
                            const TSoln& soln,
                            const int    step_number,
                            const Real   time,
                            const int    stage) = 0;

      /// Checks for conformity between soln and protoSoln.
      /**
       * Checks for conformity between soln and protoSoln.
       *
       * @param[out] soln data holder.
       * @param[in] protoSoln prototype soln data holder.
       */
      virtual bool validSolnData( const TSoln& soln, const TSoln& protoSoln ) = 0;

      /// Checks for conformity between rhs and protoSoln.
      /**
       * Checks for conformity between soln and protoSoln ignoring ghost cells.
       *
       * @param[out] rhs data holder.
       * @param[in] protoSoln prototype soln data holder.
       */
      virtual bool validRHSData( const TSoln& rhs, const TSoln& protoSoln ) = 0;

      /// Sets the TSoln data to zero.
      /**
       * Sets the TSoln data to zero.
       *
       * @param[out] soln data holder.
       */
      virtual void zeroSolnData( TSoln& soln ) = 0;

      /// Adds (scale * increment) to the TSoln data holder.
      /**
       * Adds (scale * increment) to the TSoln data holder.
       *
       * @param[out] soln solution data to which is added.
       * @param[in] increment amount to add.
       * @param[in] scale multiplier by which to scale increment.
       */
      virtual void addSolnData( TSoln& soln,
                                const TSoln& increment,
                                const Real  scale = 1.0 ) = 0;

      /// Callback to allow for additional operations after stage advance.
      /**
       * Callback to allow for additional operations after stage advance.
       *
       * @param[in,out] soln data holder.
       */
      virtual void postStageAdvance( TSoln& soln ) = 0;

};

#include "NamespaceFooter.H"

#endif
