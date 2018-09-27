#ifndef  _FLUIDOPINTERFACE_H_
#define  _FLUIDOPINTERFACE_H_

#include "FluidSpecies.H"
#include "EField.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "KineticSpecies.H"
#include "Scalar.H"
#include "GlobalDOF.H"
#include "MultiPhysicsPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "NamespaceHeader.H"
class GKVector;
class GKOps;
#include "NamespaceFooter.H"
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

      virtual void defineBlockPC( std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>&,
                                  std::vector<PS::DOFList>&,
                                  const PS::GKVector&,
                                  PS::GKOps&,
                                  const std::string&,
                                  const std::string&,
                                  bool,
                                  const FluidSpecies&,
                                  const PS::GlobalDOFFluidSpecies&,
                                  int ) { }

      virtual void updateBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>&,
                                   const PS::KineticSpeciesPtrVect&,
                                   const FluidSpeciesPtrVect&,
                                   const Real,
                                   const bool,
                                   const int ) { }

      virtual void evalSolutionOp( FluidSpeciesPtrVect&               rhs,
                                   const PS::KineticSpeciesPtrVect&   kinetic_species_comp,
                                   const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                   const FluidSpeciesPtrVect&         fluid_species_comp,
                                   const FluidSpeciesPtrVect&         fluid_species_phys,
                                   const PS::ScalarPtrVect&           scalars,
                                   const int                          component,
                                   const Real                         time )
      {
         // The default solution operator is just the identity
         rhs[component]->copy(*(fluid_species_comp[component]));
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
