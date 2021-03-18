#ifndef _FLUID_VAR_BC_FACTORY_H_
#define _FLUID_VAR_BC_FACTORY_H_

#include "FluidVarBC.H"
#include "SingleNullConfigurationBC.H"
#include "SNCoreConfigurationBC.H"
#include "LogRectConfigurationBC.H"
#include "SingleNullBlockCoordSys.H"

#include <string>

#include "NamespaceHeader.H"

class FluidVarBCFactory
{
   public:

      /// Constructor.
      /**
       */
      FluidVarBCFactory() {;}

      /// Destructor.
      /**
       */
      ~FluidVarBCFactory() {;}

      /// Create new fluid variable boundary condition object
      /**
       * Returns a dynamically allocated new FluidVarBC object of the
       * type indicated in the coord_sys_type argument.  The calling
       * routine is responsible for de-allocating the memory.
       *
       * @param[in] species_name String identifier of the species name
       * @param[in] variable_name String identifier of the variable name
       * @param[in] coord_sys_type String identified of the coordinate system type
       * @param[in] verbosity Flag to turn on additional levels of output
       *
       */

      inline FluidVarBC* create( const std::string& species_name,
                                 const std::string& variable_name,
                                 const std::string& coord_sys_type,
                                 const int& verbosity = 0 );
};

inline
FluidVarBC* FluidVarBCFactory::create( const std::string&  a_species_name,
                                       const std::string&  a_variable_name,
                                       const std::string&  a_coord_sys_type,
                                       const int&          a_verbosity )
{
   FluidVarBC* result(NULL);

   if (a_coord_sys_type == "SingleNull" ) {
      int poloidal_blocks_per_sector = SingleNullBlockCoordSys::NUM_POLOIDAL_BLOCKS;
      result = new SingleNullConfigurationBC( a_species_name, a_variable_name, poloidal_blocks_per_sector, a_verbosity );
   }
   else if (a_coord_sys_type == "SNCore" ) {
      result = new SNCoreConfigurationBC( a_species_name, a_variable_name, a_verbosity );
   }
   else if (a_coord_sys_type == "LogicallyRectangular" ) {
      result = new LogRectConfigurationBC( a_species_name, a_variable_name, a_verbosity );
   }
   else {
      MayDay::Error( "FluidVarBCFactory: Unimplemented boundary condition for given coordinate system type" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
