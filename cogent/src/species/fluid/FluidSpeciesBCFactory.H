#ifndef _FLUID_SPECIES_BC_FACTORY_H_
#define _FLUID_SPECIES_BC_FACTORY_H_

#include "ParmParse.H"
#include "FluidSpecies.H"
#include "FluidSpeciesBC.H"
#include "SingleNullConfigurationBC.H"
#include "SNCoreConfigurationBC.H"
#include "LogRectConfigurationBC.H"

#include <string>

#include "NamespaceHeader.H"

class FluidSpeciesBCFactory
{
   public:

      /// Constructor.
      /**
       */
      FluidSpeciesBCFactory() {;}

      /// Destructor.
      /**
       */
      ~FluidSpeciesBCFactory() {;}

      /// Create new InitialCondition object
      /**
       * Returns a dynamically allocated new InitialCondition object of the
       * type indicated in the passed ParmParse database.  The calling
       * routine is responsible for de-allocating the memory.
       *
       * @param[in] name String identifier of this initial condition.
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new InitialCondition object.
       * @param[in] verbosity Flag to turn on additional levels of output
       *
       * The following represents a sample input entry:
       *
       * \verbatim
       * type = "localized"
       */
      inline FluidSpeciesBC* create( const std::string& name,
                                     ParmParse& pp,
                                     const std::string& coord_sys_type,
                                     const int& verbosity = 0 );
};

inline
FluidSpeciesBC* FluidSpeciesBCFactory::create(
   const std::string& a_name,
   ParmParse& a_pp,
   const std::string& a_coord_sys_type,
   const int& a_verbosity )
{
   FluidSpeciesBC* result(NULL);
   if (a_coord_sys_type == "SingleNull" ) {
      result = new SingleNullConfigurationBC( a_name, a_pp, a_verbosity );
   }
   else if (a_coord_sys_type == "SNCore" ) {
      result = new SNCoreConfigurationBC( a_name, a_pp, a_verbosity );
   }
   else if (a_coord_sys_type == "LogicallyRectangular" ) {
      result = new LogRectConfigurationBC( a_name, a_pp, a_verbosity );
   }
   //   else if (a_coord_sys_type == ... ) {
   //      result = new ...
   //   }
   else {
      MayDay::Error( "FluidSpeciesBCFactory: Unimplemented boundary condition for given coordinate system type" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
