#ifndef  _NULLFLUIDOP_H_
#define  _NULLFLUIDOP_H_

#include "FluidOpInterface.H"

#include "NamespaceHeader.H"

/**
 * NullFluidOp fluid operator class.
*/
class NullFluidOp
   : public FluidOpInterface
{
   public:
      
      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      NullFluidOp() {;}
      
      /// Destructor.
      /**
       */
      virtual ~NullFluidOp() {;}
  
      /// Evaluates the RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  fluid RHS on the given species.
       *
       *  @param[in,out] rhs           -  contribution to be added to rhs.
       *  @param[in] kinetic_specties  -  current solution for kinetic species.
       *  @param[in] fluid_species     -  current solution for fluid species.
       *  @param[in] scalars           -  current solution for scalars.
       *  @param[in] component         -  fluid component to be computed
       *  @param[in] time              -  the time at which the fluid RHS is to be evaluated
       */
      virtual void accumulateExplicitRHS( FluidSpeciesPtrVect&               rhs,
                                          const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                          const FluidSpeciesPtrVect&         fluid_species,
                                          const PS::ScalarPtrVect&           scalars,
                                          const EMFields&                    EM_fields,
                                          const int                          component,
                                          const Real                         time) {;}

      virtual void accumulateImplicitRHS( FluidSpeciesPtrVect&               rhs,
                                          const PS::KineticSpeciesPtrVect&   kinetic_species,
                                          const FluidSpeciesPtrVect&         fluid_species,
                                          const PS::ScalarPtrVect&           scalars,
                                          const EMFields&                    EM_fields,
                                          const int                          component,
                                          const Real                         time) {;}

      virtual bool trivialSolutionOp() const {return true;};

      virtual CFGVars* createStateVariable( const string&       pp_prefix,
                                            const std::string&  name,
                                            const MagGeom&      geometry,
                                            const IntVect&      ghost_vect ) const {return NULL;}

      virtual void fillGhostCells( FluidSpecies&  species_phys,
                                   const double   time ) {;}

      const char* pp_name = {"null_op"};
};


#include "NamespaceFooter.H"

#endif
