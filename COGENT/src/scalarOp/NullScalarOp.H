#ifndef  _NULLSCALAROP_H_
#define  _NULLSCALAROP_H_

#include "ScalarOpInterface.H"

#include "NamespaceHeader.H"

/**
 * NullScalarOp scalar operator class.
*/
class NullScalarOp
   : public ScalarOpInterface
{
   public:
      
      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      NullScalarOp() {;}
      
      /// Destructor.
      /**
       */
      virtual ~NullScalarOp() {;}
  
      /// Accumulates the scalar RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  scalar RHS on the given species.
       *
       *  @param[in,out] rhs           -  contribution to be added to rhs.
       *  @param[in] kinetic_specties  -  current solution for kinetic species.
       *  @param[in] fluid_species     -  current solution for fluid species.
       *  @param[in] scalars           -  current solution for scalars
       *  @param[in] component         -  scalar component to which operator is applied.
       *  @param[in] time              -  the time at which the scalar RHS is to be evaluated
       */
      virtual void accumulateExplicitRHS( GKRHSData&                         rhs,
                                          const KineticSpeciesPtrVect&       kinetic_species,
                                          const CFG::FluidSpeciesPtrVect&    fluid_species,
                                          const ScalarPtrVect&               scalars,
                                          const CFG::EMFields&               EM_fields,
                                          const CFG::PhiOps&                 phi_ops,
                                          const int                          component,
                                          const Real                         time) {;}

      virtual void accumulateImplicitRHS( GKRHSData&                         rhs,
                                          const KineticSpeciesPtrVect&       kinetic_species,
                                          const CFG::FluidSpeciesPtrVect&    fluid_species,
                                          const ScalarPtrVect&               scalars,
                                          const CFG::EMFields&               EM_fields,
                                          const int                          component,
                                          const bool                         recompute_kinetic_terms,
                                          const Real                         time) {;}

};


#include "NamespaceFooter.H"

#endif
