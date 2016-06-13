#ifndef _GKFLUIDOP_H_
#define _GKFLUIDOP_H_

#include "REAL.H"
#include "FluidOpInterface.H"
#include "Field.H"
#include "FluidSpecies.H"
#include "ParmParse.H"

#include "NullFluidOp.H"

#include "NamespaceHeader.H"

#define _FLUID_NONE_   "None"

/**
 * Gyrokinetic fluid operator class.
*/
class GKFluidOp
{
   public:

      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      GKFluidOp( const int verbose );

      /// Destructor.
      /**
       */
      virtual ~GKFluidOp();

      /// Accumulates the RHS of the collision operator.
      /**
       * @param[out] rhs               - data holder for rhs.
       * @param[in]  fields            - current solution for fields.
       * @param[in]  fluids            - current solution for fluids.
       * @param[in]  kinetic_specties  - current solution for kinetic species.
       * @param[in]  time current time.
       */
      virtual void accumulateRHS( FluidSpeciesPtrVect&               rhs,
                                  const FieldPtrVect&                fields,
                                  const FluidSpeciesPtrVect&         fluids,
                                  const PS::KineticSpeciesPtrVect&   kinetic_species,
                                  const LevelData<FluxBox>&          E_field,
                                  const Real                         time);

      /// Compute a stable time step.
      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDt( const FieldPtrVect&         fields,
                      const FluidSpeciesPtrVect&  fluids);
   

      /// Compute time scale.
      /**
       * Computes and returns an estimate of the time scale of the collision term.
       */
      Real computeTimeScale( const FieldPtrVect&         fields,
                             const FluidSpeciesPtrVect&  fluids);

      /// returns the fluid model associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      FluidOpInterface& fluidModel( const std::string& name );

      /// returns the fluid model name associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      std::string fluidModelName( const std::string& name );

   private:

      // prevent copying
      GKFluidOp( const GKFluidOp& );
      const GKFluidOp& operator=( const GKFluidOp& );

      inline
      int getFluidModelIndex( std::map<std::string,int>::iterator& it )
      {
         return it->second;
      }
   
      bool m_verbose;
      std::map<std::string,int> m_species_map;
      std::vector<FluidOpInterface*> m_fluid_model;
      std::map<std::string,int> m_fluid_model_name;
};

#include "NamespaceFooter.H"

#endif
