#ifndef  _VIRTUALNTR_H_
#define  _VIRTUALNTR_H_

#include "KineticSpecies.H"
#include "REAL.H"
#include "NTRInterface.H"
#include "ParmParse.H"
#include <sstream>

#include "NamespaceHeader.H"

/**
 * NullNTR neutrals operator class.
*/
class NullNTR
   : public NTRInterface
{
   public:
      
      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      NullNTR() {;}
      
      /// Destructor.
      /**
       */
      virtual ~NullNTR() {;}
  
      /// Accumulates the RHS of the neutral operator.
      /**
       * @param[out] rhs data holder for rhs (without ghost cells).
       * @param[in] soln current solution with which to evaluate RHS.
       * @param[in] species index of species to which neutral is applied.
       * @param[in] time current time.
       */
      virtual void evalNtrRHS( KineticSpecies&              rhs_species,
                               const KineticSpeciesPtrVect& soln,
                               const int                    species,
                               const Real                   time ) {;}
  
      /// Compute a stable time step.
      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDtExplicitTI( const KineticSpeciesPtrVect& soln ) { return DBL_MAX; }
      Real computeDtImExTI( const KineticSpeciesPtrVect& soln ) { return DBL_MAX; }
};


#include "NamespaceFooter.H"

#endif
