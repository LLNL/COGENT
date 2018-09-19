#ifndef  _TPMINTERFACE_H_
#define  _TPMINTERFACE_H_

#include "REAL.H"
#include "ParmParse.H"
#include "KineticSpecies.H"

#include "NamespaceHeader.H"

/**
 * Transport Model Interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set transport on species KineticSpecies.
*/
class TPMInterface
{
   public:

      /// Destructor.
      /**
       */
      virtual ~TPMInterface() {;}

      /// Evaluate transport RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  transport RHS on the gven species.
       *
       *  @param[in,out] rhs - transport contribution to be added to rhs.
       *  @param[in] soln - kinetic species.
       *  @param[in] species - species to which transport model is applied.
       *  @param[in] time the time at which the transport RHS is to be evaluated
       */
      virtual void evalTpmRHS( KineticSpecies&              rhs_species,
                               const KineticSpeciesPtrVect& soln,
                               const int                    species,
                               const Real                   time ) = 0;

      virtual Real computeStableDt(const KineticSpecies& soln) = 0;
      //virtual Real computeDt(const KineticSpecies& soln, const int i) = 0;

};

#include "NamespaceFooter.H"
#endif
