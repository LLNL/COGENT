#ifndef _GKCOLLISIONS_H_
#define _GKCOLLISIONS_H_

#include "REAL.H"
#include "CLSInterface.H"
#include "KineticSpecies.H"
#include "ParmParse.H"
#include "GlobalDOF.H"
#include "AppCtxt.H"

#include "Krook.H"
#include "Lorentz.H"
#include "Linearized.H"
#include "LinearizedUnlike.H"
#include "FokkerPlanck.H"
#include "ConsDragDiff.H"
#include "NullCLS.H"

#include <string>
#include <map>

#include "NamespaceHeader.H"

#define _CLS_KROOK_              "Krook"
#define _CLS_LORENTZ_            "Lorentz"
#define _CLS_LINEARIZED_         "Linearized"
#define _CLS_LINEARIZEDUNLIKE_   "LinearizedUnlike"
#define _CLS_FOKKERPLANCK_       "FokkerPlanck"
#define _CLS_CONSDRAGDIFF_       "ConsDragDiff"
#define _CLS_NONE_               "None"

class ODEVector;

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
                                  const bool                   implicit,
                                  const Real                   time );

      /// Compute a stable time step.
      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDtExplicitTI( const KineticSpeciesPtrVect& soln );
      Real computeDtImExTI( const KineticSpeciesPtrVect& soln );

      /// Compute time scale.
      /**
       * Computes and returns an estimate of the time scale of the collision term.
       */
      Real computeTimeScale( const KineticSpeciesPtrVect& soln );

      int precondMatrixBands();

      void defineMultiPhysicsPC(std::vector<Preconditioner<ODEVector,AppCtxt>*>&,
                                std::vector<DOFList>&,
                                const KineticSpeciesPtrVect&,
                                const GlobalDOFKineticSpeciesPtrVect&,
                                const ODEVector&,
                                void*,
                                const std::string&,
                                const std::string&,
                                bool,
                                const int );

      void updateMultiPhysicsPC(std::vector<Preconditioner<ODEVector,AppCtxt>*>&,
                                const KineticSpeciesPtrVect&,
                                const GlobalDOFKineticSpeciesPtrVect&,
                                const Real, 
                                const int,
                                const int,
                                const Real, 
                                const bool);


      void assemblePrecondMatrix( void*,
                                  const KineticSpeciesPtrVect&,
                                  const GlobalDOFKineticSpeciesPtrVect&,
                                  const int,
                                  const int,
                                  const Real );

      /// returns the collision model associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      CLSInterface& collisionModel( const std::string& name,
                                    const std::string& name_bkgr);

      /// returns the collision model name associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      std::string collisionModelName( const std::string& name,
                                      const std::string& name_bkgr);

      bool isLinear();

      void preTimeStep    (const KineticSpeciesPtrVect&, const Real, const KineticSpeciesPtrVect&);
      void postTimeStage  (const KineticSpeciesPtrVect&, const Real, const int);

   private:

      // prevent copying
      GKCollisions( const GKCollisions& );
      const GKCollisions& operator=( const GKCollisions& );

      inline
      int getCollisionModelIndex( std::map<std::string,int>::iterator& it )
      {
         return it->second;
      }
   
      bool m_verbose;
      int m_num_species;
      std::map<std::string,int> m_species_map;
      std::vector<CLSInterface*> m_collision_model;
      std::vector<std::string>  m_collision_model_name;
};

#include "NamespaceFooter.H"

#endif
