#ifndef  _PRESCRIBEDSOURCES_H_
#define  _PRESCRIBEDSOURCES_H_

#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "KineticFunctionLibrary.H"
#include "NTRInterface.H"
#include "ParmParse.H"
#include <sstream>

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunction.H"
#include "GridFunctionLibrary.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

/**
 * PrescribedSources neutrals operator class.
*/
class PrescribedSources
   : public NTRInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  PrescribedSources(ParmParse&  a_ppntr,
                    const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~PrescribedSources();

  /// Accumulates the RHS of the collision operator.
  /**
   * @param[out] rhs data holder for rhs (without ghost cells).
   * @param[in] soln current solution with which to evaluate RHS.
   * @param[in] species index of species to which collision is applied.
   * @param[in] time current time.
   */
  virtual void evalNtrRHS( KineticSpecies&                  rhs_species,
                           const KineticSpeciesPtrVect&     soln,
                           const CFG::FluidSpeciesPtrVect&  fluid_species_phys,
                           const int                        species,
                           const Real                       time );

  /// Compute a stable time step.
  /**
   * Computes and returns an estimate of the maximum stable time step.
   */
  Real computeDtExplicitTI(const KineticSpeciesPtrVect& a_soln, const int a_species);
  Real computeDtImExTI(const KineticSpeciesPtrVect& a_soln, const int a_species);
  Real TimeScale(const KineticSpeciesPtrVect& a_soln, const int a_species);
   
private:
  
   void computeHeatSrc(LevelData<FArrayBox>&  heat_src,
                       const KineticSpecies&  kinetic_species,
                       const Real&            mass,
                       const Real&            time);
  
  void computeThermalBathSrc(LevelData<FArrayBox>&  heat_src,
                             const KineticSpecies&  kinetic_species,
                             const Real&            time);
  
   void getNormalizationFactors();
  
   int m_verbosity;
   bool m_first_call;
   bool m_thermal_bath;
   bool m_diagnostics;

   double m_power;
   double m_T_src;
   double m_heat_src_norm;
   
   RefCountedPtr<KineticFunction> m_ref_func;
   RefCountedPtr<CFG::GridFunction> m_heat_src_spatial_func;
   
   LevelData<FArrayBox> m_heat_src;
   LevelData<FArrayBox> m_ref_dfn;
   LevelData<FArrayBox> m_spatial_profile_inj;

  
   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "NTR.species" section
    * of ParmParse database, where "species" is replaced by the species name
    *
    * \param[in] ppntr NTR.species input database.
    */
   void parseParameters( ParmParse& pp );
   
   /// Print parameters.
   /**
    * Private method to print parameters.
    */
   void printParameters();

   //Plots particle, parallel momentum and energy sources
   void diagnostics(const LevelData<FArrayBox>& a_rhs,
                    const KineticSpecies&       a_rhs_species,
                    const double                a_time) const;

};


#include "NamespaceFooter.H"

#endif
