#ifndef  _NULLFLUIDOP_H_
#define  _NULLFLUIDOP_H_

#include "Field.H"
#include "FluidSpecies.H"
#include "REAL.H"
#include "FluidOpInterface.H"
#include "ParmParse.H"
#include <sstream>

#include "NamespaceHeader.H"

/**
 * NullFieldOp fluid operator class.
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
  
      /// Evaluates fluid RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  fluid RHS on the given species.
       *
       *  @param[in,out] rhs           -  contribution to be added to rhs.
       *  @param[in] fields            -  current solution for fields.
       *  @param[in] fluids            -  current solution for fluids.
       *  @param[in] kinetic_specties  -  current solution for kinetic species.
       *  @param[in] species           -  fluid species to which operator is applied.
       *  @param[in] time              -  the time at which the fluid RHS is to be evaluated
       */
      virtual void evalFluidRHS( FluidSpeciesPtrVect&               rhs,
                                 const FieldPtrVect&                fields,
                                 const FluidSpeciesPtrVect&         fluids,
                                 const PS::KineticSpeciesPtrVect&   kinetic_species,
                                 const LevelData<FluxBox>&          E_field,
                                 const int                          species,
                                 const Real                         time) {;}

};


#include "NamespaceFooter.H"

#endif
