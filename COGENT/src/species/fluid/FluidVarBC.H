#ifndef _FLUID_VAR_BC_H_
#define _FLUID_VAR_BC_H_

#include "FluidSpecies.H"

#include "NamespaceHeader.H"

class FluidVarBC
{
   public:

      virtual ~FluidVarBC() {;}

      virtual void apply( FluidSpecies&  a_species_phys,
                          const Real&    a_time ) = 0;

      virtual bool isForVariable( const std::string& variable_name ) const = 0;
};

#include "NamespaceFooter.H"

#endif
