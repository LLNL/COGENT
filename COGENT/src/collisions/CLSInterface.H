#ifndef  _CLSINTERFACE_H_
#define  _CLSINTERFACE_H_

#include <string>

#include "REAL.H"
#include "ParmParse.H"
#include "KineticSpecies.H"
#include "GlobalDOF.H"
#include "MultiPhysicsPreconditioner.H"

#include "NamespaceHeader.H"

class GKVector;
class GKOps;

/**
 * Collision interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set collision on species KineticSpecies.
*/
class CLSInterface
{
   public:

      /// Destructor.
      /**
       */
      virtual ~CLSInterface() {;}

      /// Evaluate collison RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  collision RHS on the gven species.
       *
       *  @param[in,out] rhs - collision contribution to be added to rhs.
       *  @param[in] soln - kinetic species.
       *  @param[in] species - species to which collision model is applied.
       *  @param[in] species_bkgr - background species, which the species collides with 
       *  @param[in] time the time at which the collison RHS is to be evaluated
       */
      virtual void evalClsRHS( KineticSpeciesPtrVect&       rhs,
                               const KineticSpeciesPtrVect& soln,
                               const int                    species,
                               const int                    species_bkgr,
                               const Real                   time ) = 0;

      virtual void evalClsRHSImplicit( KineticSpeciesPtrVect&       rhs,
                                       const KineticSpeciesPtrVect& soln,
                                       const int                    species,
                                       const int                    species_bkgr,
                                       const Real                   time ) = 0;

      virtual void evalClsRHSExplicit( KineticSpeciesPtrVect&       rhs,
                                       const KineticSpeciesPtrVect& soln,
                                       const int                    species,
                                       const int                    species_bkgr,
                                       const Real                   time ) = 0;

      virtual inline bool isLinear() { return false; }
      
      virtual Real computeDtExplicitTI(const KineticSpeciesPtrVect& soln,
                                       const int                    species) = 0;

      virtual Real computeDtImExTI(const KineticSpeciesPtrVect& soln,
                                   const int                    species) = 0;

      virtual Real computeTimeScale(const KineticSpeciesPtrVect& soln,
                                    const int                    species) = 0;

      virtual int precondMatrixBands() { return(0); }

      virtual void assemblePrecondMatrix( void*,
                                          const KineticSpecies&,
                                          const GlobalDOFKineticSpecies&,
                                          const Real ) { }

      virtual void defineBlockPC( std::vector<Preconditioner<GKVector,GKOps>*>&,
                                  std::vector<DOFList>&,
                                  const GKVector&,
                                  GKOps&,
                                  const std::string&,
                                  const std::string&,
                                  bool,
                                  const KineticSpecies&,
                                  const GlobalDOFKineticSpecies&,
                                  int ) { }

      virtual void updateBlockPC( std::vector<Preconditioner<GKVector,GKOps>*>&,
                                  const KineticSpecies&,
                                  const GlobalDOFKineticSpecies&,
                                  const Real,
                                  const Real,
                                  const bool,
                                  const int ) { }

      virtual void preTimeStep  (const KineticSpeciesPtrVect&,
                                 const int, const Real,
                                 const KineticSpeciesPtrVect&) { }

      virtual void postTimeStage(const KineticSpeciesPtrVect& soln,
                                 const int, const Real, const int ) { }
};

#include "NamespaceFooter.H"
#endif
