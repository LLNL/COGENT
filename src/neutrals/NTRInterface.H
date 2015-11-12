#ifndef  _NTRINTERFACE_H_
#define  _NTRINTERFACE_H_

#include "REAL.H"
#include "ParmParse.H"
#include "KineticSpecies.H"

#include "NamespaceHeader.H"

/**
 * Neutral interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set neutral on species KineticSpecies.
*/
class NTRInterface
{
   public:

      /// Destructor.
      /**
       */
      virtual ~NTRInterface() {;}

      /// Evaluate collison RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  neutral RHS on the gven species.
       *
       *  @param[in,out] rhs - neutral contribution to be added to rhs.
       *  @param[in] soln - kinetic species.
       *  @param[in] species - species to which neutral model is applied.
       *  @param[in] time the time at which the collison RHS is to be evaluated
       */
      virtual void evalNtrRHS( KineticSpeciesPtrVect&       rhs,
                               const KineticSpeciesPtrVect& soln,
                               const int                    species,
                               const Real                   time ) = 0;


};

#include "NamespaceFooter.H"
#endif
