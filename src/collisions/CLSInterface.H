#ifndef  _CLSINTERFACE_H_
#define  _CLSINTERFACE_H_

#include "REAL.H"
#include "ParmParse.H"
#include "KineticSpecies.H"
#include "Mapping.H"

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
                               const Real                   time,
                               const int                    flag=1) = 0;

      virtual inline bool isLinear() { return false; }
      
      virtual Real computeDt(const KineticSpeciesPtrVect& soln) { return DBL_MAX; }
      virtual Real TimeScale(const KineticSpeciesPtrVect& soln) { return DBL_MAX; }

//      virtual bool setupPrecondMatrix(void *a_P,int a_N,int a_bs) { return false; }
      virtual bool setupPrecondMatrix( void *a_P, int a_N )
   {
      return false;
   }

   virtual void assemblePrecondMatrix( void *a_P,
                                          const KineticSpeciesPtrVect& a_soln,
                                          int a_i,
                                          const Mapping& a_mapping) { }

      virtual void preTimeStep  (const KineticSpeciesPtrVect&,
                                 const int, const Real,
                                 const KineticSpeciesPtrVect&) { }

      virtual void postTimeStage(const KineticSpeciesPtrVect& soln,
                                 const int, const Real, const int ) { }
};

#include "NamespaceFooter.H"
#endif
