#ifndef  _LinearizedUnlike_H_
#define  _LinearizedUnlike_H_

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

#include "NamespaceHeader.H"

/**
 * Linearized collisions operator class for unlike species collisions.
*/
class LinearizedUnlike
   : public CLSInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  LinearizedUnlike( const std::string&  a_ppcls_str,
                    const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~LinearizedUnlike();

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
  Real computeTimeScale(const KineticSpeciesPtrVect& soln, const int);

  /// Report that collision operator is linear
  /**
   * Reports that this collision operator is linear
   */
  inline bool isLinear() {return true; }

private:

  /// Compute the test-particle collision RHS.
  /**
   * Computes the exact test-particle collision RHS.
   */
  void testPartCollRHS(LevelData<FArrayBox>&                   rhs_coll,
                       const LevelData<FArrayBox>&             delta_dfn_a,
                       const LevelData<FArrayBox>&             Nb_inj,
                       const LevelData<FArrayBox>&             Tb_inj,
                       const KineticSpecies&                   species,
                       const double                            ma,
                       const double                            mb) const;

  /// Compute the field-particle collision RHS.
  /**
   * Computes Kolesnikov et al (2010) model for field-particle collision RHS.
   */
  void fieldPartCollRHS(LevelData<FArrayBox>&       rhs_coll,
                        const LevelData<FArrayBox>& delta_dfn_b,
                        const KineticSpecies&       species_b) const;

   /// Compute the Rab * vpar and Qab terms in Eqs. (30) and (31) of Kolesnikov (2010).
   /**
    * Computes the Rab * vpar and Qab terms in Eqs. (30) and (31) of Kolesnikov (2010).
    */
   void computeRandQTerms(LevelData<FArrayBox>& R,
                          LevelData<FArrayBox>& Q,
                          const KineticSpecies& species) const;
   
   /// Compute a Maxwellian fit
   void computeMaxwellianFit(LevelData<FArrayBox>&             F0,
                             CFG::LevelData<CFG::FArrayBox>&   temperature,
                             CFG::LevelData<CFG::FArrayBox>&   density,
			     LevelData<FArrayBox>&             temperature_inj,
                             LevelData<FArrayBox>&             density_inj,
			     const KineticSpecies&             species,
			     const LevelData<FArrayBox>&       dfn,
                             const Real                        time) const;

   /// Assign Maxwellian background profile data, when it is specified from input
   void getReferenceData(LevelData<FArrayBox>&           F0,
			 CFG::LevelData<CFG::FArrayBox>& T,
			 CFG::LevelData<CFG::FArrayBox>& N,
			 LevelData<FArrayBox>&  T_inj,
			 LevelData<FArrayBox>&  N_inj,
			 const KineticSpecies&  species,
			 const RefCountedPtr<KineticFunction>& ref_func,
			 const RefCountedPtr<CFG::GridFunction>& ref_T,
			 const RefCountedPtr<CFG::GridFunction>& ref_N,
			 const Real time) const;
  
  ///Computes the self-consistent cls_freq from the density and temperature profiles
  /**
   * Computes the dimensionless normalization (to be multipled  by nb/sqrt(m_a), wher nb is the
   * Maxwellian background density (in COGENT units) and m_a is the speeces a mass in COGENT mass units
   *
   * cls_norm[s] = (4*Pi*Z^4*Coulomb_Lg/2^(3/2)) * N[1/m^3]*10^(-6) * charge_a^2 * charge_b^2 *ech[CGS]^4
   *                * L[m]*10^2 /  (T[ev]^2 * (1.6*10^(-12))^2)
   *
   * or, equvivalently in SI units
   *
   * cls_norm[s] = (1/(4*Pi*eps0^2)*Z^4*Coulomb_Lg/2^(3/2)) * N[1/m^3] * charge_a^2 * charge_b^2 * ech[SI]^2
   *                * L[m] /  T[eV]^2
   *
   */

  void getClsFreqNorm(LevelData<FArrayBox>&             cls_norm,
                      CFG::LevelData<CFG::FArrayBox>&   Ta,
                      CFG::LevelData<CFG::FArrayBox>&   Tb,
                      CFG::LevelData<CFG::FArrayBox>&   Na,
                      CFG::LevelData<CFG::FArrayBox>&   Nb,
		      const PhaseGeom&                  phase_geom) const;



  void convertToCellAverage(const PhaseGeom&  coord_sys,
                            LevelData<FArrayBox>&   dfn) const;

  void convertToCellCenters(const PhaseGeom&  coord_sys,
                            LevelData<FArrayBox>&  dfn) const;
  
  void convertToCellFaces(LevelData<FluxBox>& faceData,
                          const LevelData<FArrayBox>& cellData) const;

   
   inline void preTimeStep  ( const KineticSpeciesPtrVect& soln_mapped,
                              const int species,
                              const Real time,
                              const KineticSpeciesPtrVect& soln_physical);

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
   
   /// Diagnostics capability (left for debugging purposes)
   void diagnostics(const KineticSpecies&  a_rhs_species,
		    const int              a_species) const;
   
   
   int m_verbosity;
   bool m_time_implicit;

   Real m_mass_a;
   Real m_mass_b;
   Real m_charge_a;
   Real m_charge_b;
   int m_update_freq;
   int m_it_counter;
   bool m_second_order;
   bool m_conserve_momentum;
   bool m_conserve_energy;
   bool m_fixed_bkgr;
   bool m_fixed_cls_freq;
   bool m_first_call;
   bool m_first_stage;
   bool m_include_bkgr_cls;

   RefCountedPtr<KineticFunction> m_ref_func_a;
   RefCountedPtr<KineticFunction> m_ref_func_b;
   RefCountedPtr<CFG::GridFunction> m_ref_Ta;
   RefCountedPtr<CFG::GridFunction> m_ref_Tb;
   RefCountedPtr<CFG::GridFunction> m_ref_Na;
   RefCountedPtr<CFG::GridFunction> m_ref_Nb;
   
   LevelData<FArrayBox> m_F0a;
   LevelData<FArrayBox> m_F0b;
   LevelData<FArrayBox> m_delta_dfn_a;
   LevelData<FArrayBox> m_delta_dfn_b;
   LevelData<FArrayBox> m_Rab;
   LevelData<FArrayBox> m_Qab;
   
   LevelData<FArrayBox>  m_rhs_tmp;
   mutable LevelData<FArrayBox>  m_tmp_dfn_withGhosts;
   mutable LevelData<FArrayBox>  m_tmp_dfn_noGhosts;
   mutable LevelData<FArrayBox>  m_flux_vmu_cc;
   mutable LevelData<FluxBox>    m_flux;
   mutable LevelData<FluxBox>    m_flux_vmu;
   
   CFG::LevelData<CFG::FArrayBox> m_Ta;
   CFG::LevelData<CFG::FArrayBox> m_Tb;
   CFG::LevelData<CFG::FArrayBox> m_Na;
   CFG::LevelData<CFG::FArrayBox> m_Nb;
   LevelData<FArrayBox> m_Ta_inj;
   LevelData<FArrayBox> m_Tb_inj;
   LevelData<FArrayBox> m_Na_inj;
   LevelData<FArrayBox> m_Nb_inj;
   LevelData<FArrayBox> m_cls_norm;
   
   MomentOp& m_moment_op;
   
};


#include "NamespaceFooter.H"

#endif
