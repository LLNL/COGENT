#ifndef  _VIRTUALCLS_H_
#define  _VIRTUALCLS_H_

#include "KineticSpecies.H"
#include "REAL.H"
#include "CLSInterface.H"
#include "ParmParse.H"
#include <sstream>

#include "NamespaceHeader.H"

/**
 * NullCLS collisions operator class.
*/
class NullCLS
   : public CLSInterface
{
   public:
      
      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      NullCLS() {;}
      
      /// Destructor.
      /**
       */
      virtual ~NullCLS() {;}
  
      /// Accumulates the RHS of the collision operator.
      /**
       * @param[out] rhs data holder for rhs (without ghost cells).
       * @param[in] soln current solution with which to evaluate RHS.
       * @param[in] species index of species to which collision is applied.
       * @param[in] time current time.
       */
      virtual void evalClsRHS( KineticSpeciesPtrVect&       rhs,
                               const KineticSpeciesPtrVect& soln,
                               const int                    species,
                               const Real                   time ) {}

      inline
      virtual void evalClsRHSExplicit( KineticSpeciesPtrVect&        a_rhs,
                                       const KineticSpeciesPtrVect&  a_soln,
                                       const int                     a_species,
                                       const Real                    a_time ) {}
    
      inline
      virtual void evalClsRHSImplicit( KineticSpeciesPtrVect&        a_rhs,
                                       const KineticSpeciesPtrVect&  a_soln,
                                       const int                     a_species,
                                       const Real                    a_time ) {}

};


#include "NamespaceFooter.H"

#endif
