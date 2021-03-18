#ifndef  _NTRINTERFACE_H_
#define  _NTRINTERFACE_H_

#include "REAL.H"
#include "ParmParse.H"
#include "KineticSpecies.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluidSpecies.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

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
      virtual void evalNtrRHS( KineticSpecies&                 rhs_species,
                               const KineticSpeciesPtrVect&    kinetic_species_phys,
                               const CFG::FluidSpeciesPtrVect& fluid_species_phys,
                               const int                       species,
                               const Real                      time ) = 0;

      virtual Real computeDtExplicitTI(const KineticSpeciesPtrVect& soln, const int i) { return DBL_MAX; }
      virtual Real computeDtImExTI(const KineticSpeciesPtrVect& soln, const int i) { return DBL_MAX; }
      virtual Real TimeScale(const KineticSpeciesPtrVect& soln, const int i) { return -1; }


};

#include "NamespaceFooter.H"
#endif
