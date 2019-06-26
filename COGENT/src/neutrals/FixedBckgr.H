#ifndef  _FIXEDBCKGR_H_
#define  _FIXEDBCKGR_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "REAL.H"
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
 * FixedBckgr neutrals operator class.
*/
class FixedBckgr
   : public NTRInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  FixedBckgr( ParmParse&  a_ppntr,
        const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~FixedBckgr();

  /// Accumulates the RHS of the collision operator.
  /**
   * @param[out] rhs data holder for rhs (without ghost cells).
   * @param[in] soln current solution with which to evaluate RHS.
   * @param[in] species index of species to which collision is applied.
   * @param[in] time current time.
   */
  virtual void evalNtrRHS( KineticSpecies&               rhs_species,
                           const KineticSpeciesPtrVect&  soln,
                           const int                     species,
                           const Real                    time );

  /// Compute a stable time step.
  /**
   * Computes and returns an estimate of the maximum stable time step.
   */
  Real computeDtExplicitTI(const KineticSpeciesPtrVect& a_soln, const int a_species);
  Real computeDtImExTI(const KineticSpeciesPtrVect& a_soln, const int a_species);
  Real TimeScale(const KineticSpeciesPtrVect& a_soln, const int a_species);
   
private:

   int m_verbosity;

   bool m_include_chx, m_chx_model_friction, m_first_call;
   bool m_SI_input, m_fixed_source_dfn;

   double m_ionization_norm;
   double m_chx_norm;

   RefCountedPtr<KineticFunction>   m_neutr_func;   
   RefCountedPtr<CFG::GridFunction> m_neutr_dens;
   RefCountedPtr<CFG::GridFunction> m_neutr_vel;
   RefCountedPtr<CFG::GridFunction> m_neutr_temp;
   RefCountedPtr<CFG::GridFunction> m_ioniz_rate;

   CFG::LevelData<CFG::FArrayBox> m_neutral_density_cfg;
   CFG::LevelData<CFG::FArrayBox> m_neutral_velocity_cfg;
   
   
   LevelData<FArrayBox> m_ionization_rate;
   LevelData<FArrayBox> m_neutral_density;
   LevelData<FArrayBox> m_neutral_dfn;

   
   /// Computes the ionization contribution to RHS.
   /*
      rhs = <sigmaV_ionization> * neutr_dens * fB
    */
   void computeIonizationRHS( LevelData<FArrayBox>& rhs,
                              const LevelData<FArrayBox>& soln_dfn,
                              const LevelData<FArrayBox>& neutral_density) const;
   

   /// Computes the charge-exchange contribution to RHS.
   /*
    rhs = <sigmaV_chx> (n_ion * fB_neutral - n_neutral * fB)
    */

   void computeChargeExchangeRHS(LevelData<FArrayBox>& rhs,
                                 const LevelData<FArrayBox>& soln_dfn,
                                 const LevelData<FArrayBox>& ion_density,
                                 const LevelData<FArrayBox>& ion_temperature,
                                 const LevelData<FArrayBox>& neutral_density,
                                 const LevelData<FArrayBox>& neutral_dfn) const;

   /// Computes the model friction charge-exchange contribution to RHS.
   /*
    rhs = <sigmaV_chx> * n_neutral * vpar * upar * mass_ion * fB_unshifted_Maxw / Ti)
    where 
      upar = neutral_par_vel - ion_par_vel
      fB_unshifted_Maxw is an unshifted Maxwellian with ion density and temperature
    */
   void computeModelChargeExchangeRHS(LevelData<FArrayBox>&       rhs,
                                      const LevelData<FArrayBox>& neutral_density,
                                      const LevelData<FArrayBox>& ion_temperature,
                                      const LevelData<FArrayBox>& par_vel_shift,
                                      const LevelData<FArrayBox>& ion_unshifted_maxw,
                                      const double                mass,
                                      const PhaseGeom&            phase_geom) const;

   
   
   /// Computes normalizations to handle the case of SI_input data
   void computeChxNormalization(double&       ioniz_norm,
                                double&       chx_norm,
                                double&       dens_norm,
                                double&       vel_norm,
                                double&       temp_norm,
                                const bool   SI_input) const;

  
   
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
   
   //Plots particle and parallel momentum sources
   void diagnostics(const LevelData<FArrayBox>& a_rhs,
		    const KineticSpecies&       a_rhs_species,
		    const double                a_time) const;


};


#include "NamespaceFooter.H"

#endif
