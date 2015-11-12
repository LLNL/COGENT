#ifndef  _VIRTUALTPM_H_
#define  _VIRTUALTPM_H_

#include "KineticSpecies.H"
#include "REAL.H"
#include "TPMInterface.H"
#include "ParmParse.H"
#include <sstream>

#include "NamespaceHeader.H"

/**
 * NullTPM transport operator class.
*/
class NullTPM
   : public TPMInterface
{
   public:
      
      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      NullTPM() {;}
      
      /// Destructor.
      /**
       */
      virtual ~NullTPM() {;}
  
      /// Accumulates the RHS of the transport operator.
      /**
       * @param[out] rhs data holder for rhs (without ghost cells).
       * @param[in] soln current solution with which to evaluate RHS.
       * @param[in] species index of species to which transport is applied.
       * @param[in] time current time.
       */
      virtual void evalTpmRHS( KineticSpeciesPtrVect&       rhs,
                               const KineticSpeciesPtrVect& soln,
                               const int                    species,
                               const Real                   time ) {;}
  
      /// Compute a stable time step.
      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDt( const KineticSpeciesPtrVect& soln ) { return DBL_MAX; }
};


#include "NamespaceFooter.H"

#endif
