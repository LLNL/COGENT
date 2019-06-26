#ifndef _FLUID_VAR_BC_H_
#define _FLUID_VAR_BC_H_

#include "FluidSpecies.H"
#include "EdgeDataBox.H"

#include "NamespaceHeader.H"

class FluidVarBC
{
   public:

      virtual ~FluidVarBC() {;}
      
      virtual void apply( FluidSpecies&  a_species_phys,
                          const Real&    a_time ) = 0;
      
      virtual void applyFluxBC( const FluidSpecies&       a_species_phys,
                                LevelData<FluxBox>& a_dst,
                          const LevelData<FluxBox>& a_src,
                          const Real&               a_time ) = 0;
      
      virtual void applyEdgeBC( const FluidSpecies&       a_species_phys,
                                LevelData<EdgeDataBox>& a_dst,
                          const LevelData<EdgeDataBox>& a_src,
                          const Real&               a_time ) = 0;

      virtual bool isForVariable( const std::string& variable_name ) const = 0;
      
};


#include "NamespaceFooter.H"

#endif
