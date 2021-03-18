#ifndef _GKNEUTRALS_H_
#define _GKNEUTRALS_H_

#include "REAL.H"
#include "NTRInterface.H"
#include "KineticSpecies.H"
#include "ParmParse.H"

#include <map>

#include "NamespaceHeader.H"

/**
 * Gyrokinetic neutrals operator class.
*/
class GKNeutrals
{
   public:

      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      GKNeutrals( const int verbose );

      /// Destructor.
      /**
       */
      virtual ~GKNeutrals();

      /// Accumulates the RHS of the neutral operator.
      /**
       * @param[out] rhs data holder for rhs (without ghost cells).
       * @param[in] soln current solution with which to evaluate RHS.
       * @param[in] time current time.
       */
      virtual void accumulateRHS( KineticSpeciesPtrVect&          rhs,
                                  const KineticSpeciesPtrVect&    kinetic_species_phys,
                                  const CFG::FluidSpeciesPtrVect& fluid_species_phys,
                                  const Real                      time );

      /// Compute a stable time step.
      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDtExplicitTI( const KineticSpeciesPtrVect& soln );
      Real computeDtImExTI( const KineticSpeciesPtrVect& soln );

      /// Compute the time scale.
      /**
       * Computes and returns an estimate of the time scale.
       */
      Real computeTimeScale( const KineticSpeciesPtrVect& soln );

      /// returns the neutral model associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      NTRInterface& neutralModel( const std::string& name );

   private:

      // prevent copying
      GKNeutrals( const GKNeutrals& );
      const GKNeutrals& operator=( const GKNeutrals& );

      bool m_verbose;
      std::map<std::string,int> m_species_map;
      std::vector<NTRInterface*> m_neutral_model;
      std::map<std::string,int> m_neutral_model_name;
};

#include "NamespaceFooter.H"

#endif
