#ifndef _GKFIELDOP_H_
#define _GKFIELDOP_H_

#include "REAL.H"
#include "FluidOpInterface.H"
#include "Field.H"
#include "FluidSpecies.H"
#include "ParmParse.H"

#include "NullFieldOp.H"

#include "NamespaceHeader.H"

#define _FIELD_NONE_   "None"
#define _VORTICITY_   "Vorticity"

/**
 * Gyrokinetic field operator class.
*/
class GKFieldOp
{
   public:

      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      GKFieldOp( const int verbose );

      /// Destructor.
      /**
       */
      virtual ~GKFieldOp();

      /// Accumulates the RHS of the field operator.
      /**
       * @param[out] rhs               - data holder for rhs.
       * @param[in]  fields            - current solution for fields.
       * @param[in]  fluids            - current solution for fluids.
       * @param[in]  kinetic_specties  - current solution for kinetic species.
       * @param[in]  time              - current time.
       */
      virtual void accumulateRHS( FieldPtrVect&                      rhs,
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

      /// returns the field model associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      FieldOpInterface& fieldModel( const std::string& name );

      /// returns the field model name associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      std::string fieldModelName( const std::string& name );

   private:

      // prevent copying
      GKFieldOp( const GKFieldOp& );
      const GKFieldOp& operator=( const GKFieldOp& );

      inline
      int getfieldModelIndex( std::map<std::string,int>::iterator& it )
      {
         return it->second;
      }
   
      bool m_verbose;
      std::map<std::string,int> m_fieldVecComp_map;
      std::vector<FieldOpInterface*> m_field_model;
      std::map<std::string,int> m_field_model_name;
};

#include "NamespaceFooter.H"

#endif
