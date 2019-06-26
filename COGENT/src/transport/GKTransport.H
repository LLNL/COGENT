#ifndef _GKTRANSPORT_H_
#define _GKTRANSPORT_H_

#include "REAL.H"
#include "TPMInterface.H"
#include "KineticSpecies.H"
#include "ParmParse.H"

#include <map>
#include "NamespaceHeader.H"

/**
 * Gyrokinetic transport operator class.
*/
class GKTransport
{
   public:

      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      GKTransport( const int verbose );

      /// Destructor.
      /**
       */
      virtual ~GKTransport();

      /// Accumulates the RHS of the transport operator.
      /**
       * @param[out] rhs data holder for rhs (without ghost cells).
       * @param[in] soln current solution with which to evaluate RHS.
       * @param[in] time current time.
       */
      virtual void accumulateRHS( KineticSpeciesPtrVect&       rhs,
                                  const KineticSpeciesPtrVect& soln,
                                  const Real                   time );

      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDtExplicitTI( const KineticSpeciesPtrVect& soln );
      Real computeDtImExTI( const KineticSpeciesPtrVect& soln );
      Real computeTimeScale( const KineticSpeciesPtrVect& soln );

   private:

      // prevent copying
      GKTransport( const GKTransport& );
      const GKTransport& operator=( const GKTransport& );

      /// returns the transport model associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      TPMInterface& transportModel( const std::string& name );

      bool m_verbose;
      std::map<std::string,int> m_species_map;
      std::vector<TPMInterface*> m_transport_model;
      std::map<std::string,int> m_transport_model_name;


};

#include "NamespaceFooter.H"

#endif
