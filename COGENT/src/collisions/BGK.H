#ifndef  _BGK_H_
#define  _BGK_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "REAL.H"
#include "CLSInterface.H"
#include "ParmParse.H"
#include <sstream>

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunction.H"
#include "GridFunctionLibrary.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "SpaceUtils.H.multidim"

#include "NamespaceHeader.H"

/**
 * BGK collisions operator class.
*/
class BGK
   : public CLSInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  BGK(const std::string&  a_ppcls,
      const int           a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~BGK();

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
                           const int                    species_bkgr,
                           const Real                   time );

  inline
  virtual void evalClsRHSExplicit( KineticSpeciesPtrVect&        a_rhs,
                                   const KineticSpeciesPtrVect&  a_soln,
                                   const int                     a_species,
                                   const int                     a_species_bkgr,
                                   const Real                    a_time )
  {
    if (!m_time_implicit) evalClsRHS(a_rhs, a_soln, a_species, a_species_bkgr, a_time);
  }

  inline
  virtual void evalClsRHSImplicit( KineticSpeciesPtrVect&        a_rhs,
                                   const KineticSpeciesPtrVect&  a_soln,
                                   const int                     a_species,
                                   const int                     a_species_bkgr,
                                   const Real                    a_time )
  {
    if (m_time_implicit) evalClsRHS(a_rhs, a_soln, a_species, a_species_bkgr, a_time);
  }

  /// Compute a stable time step.
  /**
   * Computes and returns an estimate of the maximum stable time step.
   */
  Real computeDtExplicitTI(const KineticSpeciesPtrVect& soln, const int);
  Real computeDtImExTI(const KineticSpeciesPtrVect& soln, const int);

  /// Compute the time scale.
  /**
   * Computes and returns an estimate of the collision term time scale.
   */
  Real computeTimeScale(const KineticSpeciesPtrVect& soln, const int );

  Real collisionFrequency ();

  /// Report that collision operator is linear
  /**
   * Reports that this collision operator is linear
   */
  inline bool isLinear() {return true; }

private:

   void computeReferenceSolution(KineticSpecies& ref_species,
                                 const KineticSpecies& soln_species,
                                 const Real time);
   
   void computeTwoConservationCoeff(CFG::LevelData<CFG::FArrayBox>& coeff);
   
   void computeThreeConservationCoeff(CFG::LevelData<CFG::FArrayBox>& coeff,
                                      const KineticSpecies& soln_species);
   

   inline void preTimeStep  ( const KineticSpeciesPtrVect& soln_mapped,
                              const int species,
                              const Real time,
                              const KineticSpeciesPtrVect& soln_physical);
   
  ///Computes the self-consistent cls_freq from the density and temperature profiles
  /**
   * Computes the dimensionless normalization cls_norm * n/T^(3/2), wher n and T are the 
   * Maxwellian density and temeprature profiles specified by grid functions, and
   * 
   * cls_norm[s] = (4*Pi*Z^4*Coulomb_Lg/2^(3/2)) * N[1/m^3]*10^(-6) * a_charge^4 *ech[CGS]^4 
   *                * L[m]*10^2 / (a_mass^(1/2) * T[ev]^2 * (1.6*10^(-12))^2)
   *
   * or, equvivalently in SI units
   *
   * cls_norm[s] = (1/(4*Pi*eps0^2)*Z^4*Coulomb_Lg/2^(3/2)) * N[1/m^3] * a_charge^4 * ech[SI]^2 
   *                * L[m] / (a_mass^(1/2) * T[eV]^2) 
   * 
   */

  void computeSelfConsistFreq(LevelData<FArrayBox>&       cls_freq,
                              const LevelData<FArrayBox>& density,
                              const LevelData<FArrayBox>& temperature,
                              const double                mass,
                              const double                charge      ) const;


   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "CLS.species" section
    * of ParmParse database, where "species" is replaced by the species name
    *
    * \param[in] ppcls CLS.species input database.
    */
   void parseParameters( ParmParse& pp );
   
   /// Print parameters.
   /**
    * Private method to print parameters.
    */
   void printParameters();
   
   
   int m_verbosity;
   int m_update_freq;
   int m_it_counter;
   
   bool m_time_implicit;

   Real m_cls_freq;
   bool m_conserve_two_moments;
   bool m_conserve_three_moments;
   bool m_fixed_cls_freq;
   bool m_first_call;
   bool m_first_stage;

   RefCountedPtr<KineticSpecies> m_ref_species;
   
   LevelData<FArrayBox> m_density_inj;
   LevelData<FArrayBox> m_temperature_inj;
   LevelData<FArrayBox> m_sc_cls_freq;
   LevelData<FArrayBox> m_coeff_inj;
   LevelData<FArrayBox> m_kernels;
   CFG::LevelData<CFG::FArrayBox> m_density;
   CFG::LevelData<CFG::FArrayBox> m_temperature;
   CFG::LevelData<CFG::FArrayBox> m_parallel_part_flux;
   
};


#include "NamespaceFooter.H"

#endif
