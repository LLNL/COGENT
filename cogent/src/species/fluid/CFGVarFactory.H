#ifndef _CFGVARFACTORY_H_
#define _CFGVARFACTORY_H_

//#include "ParmParse.H"
#include "EField.H"
#include "EFieldAmpere.H"
#include "EFieldSelfConsistentBC.H"
#include "CellVar.H"
#include "FluidSpecies.H"

#include <string>

#include "NamespaceHeader.H"

class CFGVarFactory
{
   public:

      /// Constructor.
      /**
       */
      CFGVarFactory() {;}

      /// Destructor.
      /**
       */
      ~CFGVarFactory() {;}

      /// Creates a new CFGVar object
      /**
       * Returns a dynamically allocated new Field object of the
       * type indicated in the passed in string a_type.  The calling
       * routine is responsible for de-allocating the memory.
       *
       * @param[in] name String identifier of this Field variable
       */
     inline CFGVar* create( const string&      pp_prefix,
                            const std::string& name,
                            const std::string& type,
                            const MagGeom&     geometry,
                            const IntVect&     ghost_vect );
};

inline
CFGVar* CFGVarFactory::create( const string&      a_pp_prefix,
                               const std::string& a_name,
                               const std::string& a_type,
                               const MagGeom&     a_geometry,
                               const IntVect&     a_ghost_vect )
{
   CFGVar* result(NULL);

   if (a_type == "VorticityOp" || a_type == "AmpereErAverageOp" ) {
      result = new CellVar(a_pp_prefix, a_name, a_geometry, 1, a_ghost_vect);
   }
   else if (a_type == "FluidSpecies" ) {
      result = new FluidSpecies(a_pp_prefix, a_name, a_geometry, a_ghost_vect);
   }
   else if (a_type == "OneFieldOp" ) {
      result = new FluidSpecies(a_pp_prefix, a_name, a_geometry, a_ghost_vect);
   }
   else if (a_type == "EField" ) {
      result = new EField(a_pp_prefix, a_name, a_geometry, a_ghost_vect);
   }
   else if (a_type == "EFieldAmpere" ) {
      result = new EFieldAmpere(a_pp_prefix, a_name, a_geometry, a_ghost_vect);
   }
   else if (a_type == "EFieldSelfConsistentBC" ) {
      result = new EFieldSelfConsistentBC(a_pp_prefix, a_name, a_geometry, a_ghost_vect);
   }
   else {
      MayDay::Error("CFGVarFactory: Unrecognized type");
   }

   return result;
}

#include "NamespaceFooter.H"

#endif
