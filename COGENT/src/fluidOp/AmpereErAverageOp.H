#ifndef  _AMPEREERAVERAGEOP_H_
#define  _AMPEREERAVERAGEOP_H_

#include "FluidOpInterface.H"

#include "NamespaceHeader.H"

/**
 * AmpereErAverageOp fluid operator class.  This is essentially a dummy clase that is only needed
 * for its parseInitialConditions() and createStateVariable() members.  Everything else normally
 * provided by a FluidOp is actually performed by the VlasovAmpereOp, for which this special class
 * is a companion
*/
class AmpereErAverageOp
   : public FluidOpInterface
{
   public:
      
      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      AmpereErAverageOp()
      {parseInitialConditions("", "Er_flux_surfaces");}
      
      /// Destructor.
      /**
       */
      virtual ~AmpereErAverageOp() {;}
  
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
                                          const EField&                      E_field,
                                          const int                          component,
                                          const Real                         time) {;}

      virtual void accumulateImplicitRHS( FluidSpeciesPtrVect&               rhs,
                                          const PS::KineticSpeciesPtrVect&   kinetic_species,
                                          const FluidSpeciesPtrVect&         fluid_species,
                                          const PS::ScalarPtrVect&           scalars,
                                          const EField&                      E_field,
                                          const int                          component,
                                          const bool                         recompute_kinetic_terms,
                                          const Real                         time) {;}

      virtual bool trivialSolutionOp() const {return true;};

      virtual CFGVars* createStateVariable( const string&       a_pp_prefix,
                                            const std::string&  a_name,
                                            const MagGeom&      a_geometry,
                                            const IntVect&      a_ghost_vect ) const
      {
         return new CellVar(a_pp_prefix, a_name, a_geometry, 1, a_ghost_vect);
      }

      virtual void fillGhostCells( FluidSpecies&  species_phys,
                                   const double   time ) {;}
};


#include "NamespaceFooter.H"

#endif
