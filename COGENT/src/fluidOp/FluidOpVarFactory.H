#ifndef _FLUIDOPVARFACTORY_H_
#define _FLUIDOPVARFACTORY_H_

#include "CellVar.H"
#include "OneFieldOp.H"
#include "BurgersOp.H"
#include "IdealMhdOp.H"
#include "FullFluidOp.H"
#include "TwoFieldNeutralsOp.H"
#include "AmpereErAverageOp.H"

#include <string>

#include "NamespaceHeader.H"

class FluidOpVarFactory
{
   public:

      /// Constructor.
      /**
       */
      FluidOpVarFactory() {;}

      /// Destructor.
      /**
       */
      ~FluidOpVarFactory() {;}

      /// Creates a new CFGVars object
      /**
       * Returns a dynamically allocated new Field object of the
       * type indicated in the passed in string a_type.  The calling
       * routine is responsible for de-allocating the memory.
       *
       * @param[in] name String identifier of this Field variable
       */
     inline CFGVars* create( const string&      pp_prefix,
                             const std::string& name,
                             const std::string& type,
                             const MagGeom&     geometry,
                             const IntVect&     ghost_vect );
};

inline
CFGVars* FluidOpVarFactory::create( const string&      a_pp_prefix,
                                    const std::string& a_species_name,
                                    const std::string& a_op_type,
                                    const MagGeom&     a_geometry,
                                    const IntVect&     a_ghost_vect )
{
   CFGVars* result(NULL);

   const std::string prefix( "fluid_model." + a_species_name );

   if ( a_op_type == "VorticityOp" ) {
      VorticityOp op(prefix, a_geometry, 0., false);
      result = op.createStateVariable(a_pp_prefix, a_species_name, a_geometry, a_ghost_vect);
   }
   else if (a_op_type == "OneFieldOp" ) {
      OneFieldOp op(prefix, a_species_name, a_geometry, false);
      result = op.createStateVariable(a_pp_prefix, a_species_name, a_geometry, a_ghost_vect);
   }
   else if (a_op_type == "BurgersOp" ) {
      BurgersOp op(prefix, a_species_name, a_geometry, true);
      result = op.createStateVariable(a_pp_prefix, a_species_name, a_geometry, a_ghost_vect);
   }
   else if (a_op_type == "IdealMhdOp" ) {
      IdealMhdOp op(prefix, a_species_name, a_geometry, true);
      result = op.createStateVariable(a_pp_prefix, a_species_name, a_geometry, a_ghost_vect);
   }
   else if (a_op_type == "FullFluidOp" ) {
      FullFluidOp op(prefix, a_species_name, a_geometry, true);
      result = op.createStateVariable(a_pp_prefix, a_species_name, a_geometry, a_ghost_vect);
   }
   else if (a_op_type == "TwoFieldNeutralsOp" ) {
      TwoFieldNeutralsOp op(prefix, a_species_name, a_geometry, true);
      result = op.createStateVariable(a_pp_prefix, a_species_name, a_geometry, a_ghost_vect);
   }
   else if ( a_op_type == "AmpereErAverageOp" ) {
      AmpereErAverageOp op;
      result = op.createStateVariable(a_pp_prefix, a_species_name, a_geometry, a_ghost_vect);
   }
   else {
      MayDay::Error("FluidOpVarFactory: Unrecognized op_type");
   }

   return result;
}

#include "NamespaceFooter.H"

#endif
