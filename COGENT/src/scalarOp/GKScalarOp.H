#ifndef _GKSCALAROP_H_
#define _GKSCALAROP_H_

#include "ScalarOpInterface.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluidSpecies.H"
#include "EField.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

/**
 * Scalar operator class.
*/
class GKScalarOp
{
   public:

      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
     GKScalarOp( const int verbose );

      /// Destructor.
      /**
       */
      virtual ~GKScalarOp();

      /// Accumulates the RHS of the scalar operator.
      /**
       * @param[out] rhs               - data holder for rhs.
       * @param[in]  fluids            - current solution for fluids.
       * @param[in]  kinetic_specties  - current solution for kinetic species.
       * @param[in]  time              - current time.
       */

   virtual void accumulateRHS( GKRHSData&                     rhs,
                                  const KineticSpeciesPtrVect&       kinetic_species,
                                  const CFG::FluidSpeciesPtrVect&    fluid_species,
                                  const ScalarPtrVect&               scalars,
                                  const CFG::EField&                 E_field,
                                  const bool                         implicit,
                                  const bool                         recompute_kinetic_terms,
                                  const Real                         time);

      /// Compute a stable time step.
      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDt( const ScalarPtrVect& scalars);

      /// Compute time scale.
      /**
       * Computes and returns an estimate of the time scale of the collision term.
       */
      Real computeTimeScale( const ScalarPtrVect&    scalars);

      /// returns the operator associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      ScalarOpInterface& scalarOp( const std::string& name );

      /// returns the scalar operator name associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      std::string scalarOpName( const std::string& name );

   private:

      // prevent copying
      GKScalarOp( const GKScalarOp& );
      const GKScalarOp& operator=( const GKScalarOp& );

      inline
      int getscalarOpIndex( std::map<std::string,int>::iterator& it )
      {
         return it->second;
      }
   
      bool m_verbose;
      std::map<std::string,int> m_scalarVecComp_map;
      std::vector<ScalarOpInterface*> m_scalar_op;
      std::map<std::string,int> m_scalar_op_name;
};

#include "NamespaceFooter.H"

#endif
