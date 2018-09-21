#ifndef  _FLUIDOPINTERFACE_H_
#define  _FLUIDOPINTERFACE_H_

#include "FluidSpecies.H"
#include "EField.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "KineticSpecies.H"
#include "Scalar.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H"
namespace PS = PS_NAMESPACE;

/**
 * FluidOp interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set operators on fluid species.
*/
class FluidOpInterface
{
   public:

      /// Destructor.
      /**
       */
      virtual ~FluidOpInterface() {;}

      /// Evaluates the RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  fluid RHS on the given species.
       *
       *  @param[in,out] rhs           -  contribution to be added to rhs.
       *  @param[in] kinetic_specties  -  current solution for kinetic species.
       *  @param[in] fluid_species     -  current solution for fluid species.
       *  @param[in] component         -  fluid component to be computed
       *  @param[in] time              -  the time at which the fluid RHS is to be evaluated.
       */
      virtual void accumulateExplicitRHS( FluidSpeciesPtrVect&               rhs,
                                          const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                          const FluidSpeciesPtrVect&         fluid_species,
                                          const PS::ScalarPtrVect&           scalars,
                                          const EField&                      E_field,
                                          const int                          component,
                                          const Real                         time) = 0;

      virtual void accumulateImplicitRHS( FluidSpeciesPtrVect&               rhs,
                                          const PS::KineticSpeciesPtrVect&   kinetic_species,
                                          const FluidSpeciesPtrVect&         fluid_species,
                                          const PS::ScalarPtrVect&           scalars,
                                          const EField&                      E_field,
                                          const int                          component,
                                          const bool                         recompute_kinetic_terms,
                                          const Real                         time) = 0;

      virtual void evalSolutionOp( FluidSpeciesPtrVect&               rhs,
                                   const PS::KineticSpeciesPtrVect&   kinetic_species,
                                   const FluidSpeciesPtrVect&         fluid_species,
                                   const PS::ScalarPtrVect&           scalars,
                                   const int                          component,
                                   const Real                         time )
      {
         rhs[component]->copy(*(fluid_species[component]));
         MayDay::Error("FluidOpInterface::evalSolutionOp() Error: this default implementation seems to be incorrect since GKOps will pass the *physical* fluid species while rhs needs to be set to the *computational* fluid species. This line will be removed once this issue has been sorted.\n");
      }

      virtual void updatePCImEx( const PS::KineticSpeciesPtrVect& kinetic_species,
                                 const double                     mshift ) {};

      virtual void solveSolutionPC( FluidSpeciesPtrVect&              fluid_species_solution,
                                    const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                                    const FluidSpeciesPtrVect&        fluid_species_rhs,
                                    const int                         component )
      {
         // The default solution operator is just the identity
         fluid_species_solution[component]->copy(*(fluid_species_rhs[component]));
      }

      virtual void updatePC( const PS::KineticSpeciesPtrVect&  kinetic_species,
                             const EllipticOpBC&               bc ) {};

      virtual void solvePCImEx( FluidSpeciesPtrVect&              fluid_species_solution,
                                const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                                const FluidSpeciesPtrVect&        fluid_species_rhs,
                                const int                         component )
      {
         // The default solution operator is just the identity
         fluid_species_solution[component]->copy(*(fluid_species_rhs[component]));
      }

      virtual Real computeDt(const FluidSpeciesPtrVect& fluids)
      {
         return DBL_MAX;
      }
   
      virtual Real TimeScale(const FluidSpeciesPtrVect& fluids)
      {
         return DBL_MAX;
      }

      virtual bool trivialSolutionOp() const = 0;
};

#include "NamespaceFooter.H"
#endif
