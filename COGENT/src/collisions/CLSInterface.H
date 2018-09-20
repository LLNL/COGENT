#ifndef  _CLSINTERFACE_H_
#define  _CLSINTERFACE_H_

#include "REAL.H"
#include "ParmParse.H"
#include "KineticSpecies.H"
#include "GlobalDOF.H"

#include "NamespaceHeader.H"

/**
 * Collision interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set collision on species KineticSpecies.
*/
class CLSInterface
{
   public:

      /// Destructor.
      /**
       */
      virtual ~CLSInterface() {;}

      /// Evaluate collison RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  collision RHS on the gven species.
       *
       *  @param[in,out] rhs - collision contribution to be added to rhs.
       *  @param[in] soln - kinetic species.
       *  @param[in] species - species to which collision model is applied.
       *  @param[in] time the time at which the collison RHS is to be evaluated
       */
      virtual void evalClsRHS( KineticSpeciesPtrVect&       rhs,
                               const KineticSpeciesPtrVect& soln,
                               const int                    species,
                               const Real                   time ) = 0;

      virtual void evalClsApproxRHS(  KineticSpeciesPtrVect&       rhs,
                                      const KineticSpeciesPtrVect& soln,
                                      const int                    species,
                                      const Real                   time )
      {
        /* default implementation */
        evalClsRHS(rhs, soln, species, time);
      }

      virtual inline bool isLinear() { return false; }
      
      virtual Real computeDt(const KineticSpeciesPtrVect& soln, const int i) { return DBL_MAX; }
      virtual Real TimeScale(const KineticSpeciesPtrVect& soln, const int i) { return -1; }

      virtual int precondMatrixBands() { return(0); }

      virtual void assemblePrecondMatrix( void*,
                                          const KineticSpecies&,
                                          const GlobalDOFKineticSpecies& ) { }

      virtual void preTimeStep  (const KineticSpeciesPtrVect&,
                                 const int, const Real,
                                 const KineticSpeciesPtrVect&) { }

      virtual void postTimeStage(const KineticSpeciesPtrVect& soln,
                                 const int, const Real, const int ) { }
};

#include "NamespaceFooter.H"
#endif
