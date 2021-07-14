#ifndef  _SELFCONSISTENTBCOP_H_
#define  _SELFCONSISTENTBCOP_H_

#include "ScalarOpInterface.H"
#include "GKVlasovAmpere.H"

#include "NamespaceHeader.H"

/**
 * SelfConsistentBCOp  class.
 *
*/
class SelfConsistentBCOp
   : public ScalarOpInterface
{
  public:

    SelfConsistentBCOp(const ParmParse&    pp,
                       const PhaseGeom&    geometry,
                       const double        larmor,
                       const int           verbosity);

    virtual ~SelfConsistentBCOp();

/// Accumulates the scalar RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  scalar RHS on the given species.
       *
       *  @param[in,out] rhs           -  contribution to be added to rhs.
       *  @param[in] kinetic_specties  -  current solution for kinetic species.
       *  @param[in] fluid_species     -  current solution for fluid species.
       *  @param[in] scalars           -  current solution for scalars.
       *  @param[in] component         -  scalar component to which operator is applied.
       *  @param[in] time              -  the time at which the scalar RHS is to be evaluated.
       */
      virtual void accumulateExplicitRHS( GKRHSData&                         rhs,
                                          const KineticSpeciesPtrVect&       kinetic_species,
                                          const CFG::FluidSpeciesPtrVect&    fluid_species,
                                          const ScalarPtrVect&               scalars,
                                          const CFG::EField&                 E_field,
                                          const int                          component,
                                          const Real                         time);

      virtual void accumulateImplicitRHS( GKRHSData&                         rhs,
                                          const KineticSpeciesPtrVect&       kinetic_species,
                                          const CFG::FluidSpeciesPtrVect&    fluid_species,
                                          const ScalarPtrVect&               scalars,
                                          const CFG::EField&                 E_field,
                                          const int                          component,
                                          const bool                         recompute_kinetic_terms,
                                          const Real                         time);
  static const char* pp_name;
  
  private:

    GKVlasovAmpere* m_vlasov;
    int m_verbosity;
    const PhaseGeom& m_geometry;
    GKRHSData        m_rhs;
    Real m_larmor;
    bool m_is_vlasov_implicit;
   
};

#include "NamespaceFooter.H"
#endif
