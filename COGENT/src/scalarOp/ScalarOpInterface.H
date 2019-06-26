#ifndef  _SCALAROPINTERFACE_H_
#define  _SCALAROPINTERFACE_H_

#include "KineticSpecies.H"
#include "Scalar.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluidSpecies.H"
#include "EField.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

/**
 * ScalarOp interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set operators on scalar species.
*/
class ScalarOpInterface
{
   public:

      /// Destructor.
      /**
       */
      virtual ~ScalarOpInterface() {;}

      /// Accumulates the scalar RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  scalar RHS on the given species.
       *
       *  @param[in,out] rhs           -  contribution to be added to rhs.
       *  @param[in] kinetic_specties  -  current solution for kinetic species.
       *  @param[in] fluid_species     -  current solution for fluid species.
       *  @param[in] scalars           -  current solution for scalars.
       *  @param[in] E_field           -  current electric field
       *  @param[in] component         -  scalar component to which operator is applied.
       *  @param[in] time              -  the time at which the scalar RHS is to be evaluated.
       */
      virtual void accumulateExplicitRHS( GKRHSData&                         rhs,
                                          const KineticSpeciesPtrVect&       kinetic_species,
                                          const CFG::FluidSpeciesPtrVect&    fluid_species,
                                          const ScalarPtrVect&               scalars,
                                          const CFG::EField&                 E_field,
                                          const int                          component,
                                          const Real                         time ) = 0;

      virtual void accumulateImplicitRHS( GKRHSData&                         rhs,
                                          const KineticSpeciesPtrVect&       kinetic_species,
                                          const CFG::FluidSpeciesPtrVect&    fluid_species,
                                          const ScalarPtrVect&               scalars,
                                          const CFG::EField&                 E_field,
                                          const int                          component,
                                          const bool                         recompute_kinetic_terms,
                                          const Real                         time ) = 0;

      virtual Real computeDtExplicitTI(const ScalarPtrVect& scalars)
      {
         return DBL_MAX;
      }
   
      virtual Real computeDtImExTI(const ScalarPtrVect& scalars)
      {
         return DBL_MAX;
      }
   
   virtual Real TimeScale(const ScalarPtrVect& scalars)
      {
         return DBL_MAX;
      }
};

#include "NamespaceFooter.H"
#endif
