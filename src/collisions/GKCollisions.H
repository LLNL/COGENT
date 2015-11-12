#ifndef _GKCOLLISIONS_H_
#define _GKCOLLISIONS_H_

#include "REAL.H"
#include "CLSInterface.H"
#include "KineticSpecies.H"
#include "ParmParse.H"

#include <map>

#include "NamespaceHeader.H"

/**
 * Gyrokinetic collisions operator class.
*/
class GKCollisions
{
   public:

      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      GKCollisions( const int verbose );

      /// Destructor.
      /**
       */
      virtual ~GKCollisions();

      /// Accumulates the RHS of the collision operator.
      /**
       * @param[out] rhs data holder for rhs (without ghost cells).
       * @param[in] soln current solution with which to evaluate RHS.
       * @param[in] time current time.
       */
      virtual void accumulateRHS( KineticSpeciesPtrVect&       rhs,
                                  const KineticSpeciesPtrVect& soln,
                                  const Real                   time );

      /// Compute a stable time step.
      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDt( const KineticSpeciesPtrVect& soln );

      /// returns the collision model associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      CLSInterface& collisionModel( const std::string& name );

   private:

      // prevent copying
      GKCollisions( const GKCollisions& );
      const GKCollisions& operator=( const GKCollisions& );

      bool m_verbose;
      std::map<std::string,int> m_species_map;
      std::vector<CLSInterface*> m_collision_model;
};

#include "NamespaceFooter.H"

#endif
