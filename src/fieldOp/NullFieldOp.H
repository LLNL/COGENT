#ifndef  _NULLFIELDOP_H_
#define  _NULLFIELDOP_H_

#include "Field.H"
#include "FluidSpecies.H"
#include "REAL.H"
#include "FieldOpInterface.H"
#include "ParmParse.H"
#include <sstream>

#include "NamespaceHeader.H"

/**
 * NullFieldOp field operator class.
*/
class NullFieldOp
   : public FieldOpInterface
{
   public:
      
      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      NullFieldOp() {;}
      
      /// Destructor.
      /**
       */
      virtual ~NullFieldOp() {;}
  
      /// Evaluates field RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  field RHS on the given field.
       *
       *  @param[in,out] rhs           -  contribution to be added to rhs.
       *  @param[in] fields            -  current solution for fields.
       *  @param[in] fluids            -  current solution for fluids.
       *  @param[in] kinetic_specties  -  current solution for kinetic species.
       *  @param[in] fieldVecComp      -  component of the field vector to which operator is applied.
       *  @param[in] time              -  the time at which the field RHS is to be evaluated.
       */
      virtual void evalFieldRHS( FieldPtrVect&                      rhs,
                                 const FieldPtrVect&                fields,
                                 const FluidSpeciesPtrVect&         fluids,
                                 const PS::KineticSpeciesPtrVect&   kinetic_species,
                                 const LevelData<FluxBox>&          E_field,
                                 const int                          fieldVecComp,
                                 const Real                         time) {;}

};


#include "NamespaceFooter.H"

#endif
