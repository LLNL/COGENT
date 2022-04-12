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
  
  inline void preTimeStep  ( const KineticSpeciesPtrVect& a_soln_mapped,
                             const int a_species,
                             const Real a_time,
                             const KineticSpeciesPtrVect& a_soln_physical);
   
private:
  
   void computeHeatSrc(LevelData<FArrayBox>&  src,
                       const KineticSpecies&  kinetic_species,
                       const Real&            mass,
                       const Real&            time);
  
   void computeRelaxationSrc(LevelData<FArrayBox>&  src,
                            const KineticSpecies&  kinetic_species,
                            const Real&            time);

   void computeThermalizationSrc(LevelData<FArrayBox>&  src,
                                 const KineticSpecies&  kinetic_species,
                                 const Real&            time);
  
   void getNormalizationFactors();
  
   int m_verbosity;
   int m_it_counter;
   int m_update_freq;
   bool m_first_call;
   bool m_include_relaxation;
   bool m_include_thermalization;
   bool m_include_heat_src;
   bool m_diagnostics;

   double m_power;
   double m_T_src;
   double m_heat_src_norm;
   
   RefCountedPtr<KineticFunction>   m_ref_func;
   RefCountedPtr<CFG::GridFunction> m_heat_src_profile_func;
   RefCountedPtr<CFG::GridFunction> m_relaxation_profile_func;
   RefCountedPtr<CFG::GridFunction> m_thermalization_profile_func;
   RefCountedPtr<CFG::GridFunction> m_maxw_density_func;
   RefCountedPtr<CFG::GridFunction> m_maxw_vparallel_func;
  
   LevelData<FArrayBox> m_heat_src;
   LevelData<FArrayBox> m_relaxation_src;
   LevelData<FArrayBox> m_thermalization_src;
   LevelData<FArrayBox> m_ref_dfn;
   LevelData<FArrayBox> m_maxw_dfn;
   LevelData<FArrayBox> m_relaxation_profile_inj;
   LevelData<FArrayBox> m_thermalization_profile_inj;
   LevelData<FArrayBox> m_heat_src_profile_inj;
   CFG::LevelData<CFG::FArrayBox> m_maxw_density;
   CFG::LevelData<CFG::FArrayBox> m_maxw_vparallel;

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
