#ifndef  _Lorentz_H_
#define  _Lorentz_H_

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
 * Lorentz collisions operator class.
*/
class Lorentz
   : public CLSInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  Lorentz( const std::string&  a_ppcls_str,
        const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~Lorentz();

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
                           const Real                   time );

  inline
  virtual void evalClsRHSExplicit( KineticSpeciesPtrVect&        a_rhs,
                                   const KineticSpeciesPtrVect&  a_soln,
                                   const int                     a_species,
                                   const Real                    a_time )
  {
    if (!m_time_implicit) evalClsRHS(a_rhs, a_soln, a_species, a_time);
  }

  inline
  virtual void evalClsRHSImplicit( KineticSpeciesPtrVect&        a_rhs,
                                   const KineticSpeciesPtrVect&  a_soln,
                                   const int                     a_species,
                                   const Real                    a_time )
  {
    if (m_time_implicit) evalClsRHS(a_rhs, a_soln, a_species, a_time);
  }

  /// Compute a stable time step.
  /**
   * Computes and returns an estimate of the maximum stable time step.
   */
  Real computeDtExplicitTI(const KineticSpeciesPtrVect& soln, const int);
  Real computeDtImExTI(const KineticSpeciesPtrVect& soln, const int);
  Real computeTimeScale(const KineticSpeciesPtrVect& soln, const int);


private:

  int m_verbosity;
  bool m_time_implicit;

  Real m_cls_freq;
  bool m_constant_freq;
  bool m_conserve_momentum;
  bool m_fixed_cls_freq;
  bool m_first_step;

  RefCountedPtr<KineticFunction> m_ref_func;
  RefCountedPtr<CFG::GridFunction> m_ref_temp;
  RefCountedPtr<CFG::GridFunction> m_ref_dens;

  LevelData<FArrayBox> m_temperature;
  LevelData<FArrayBox> m_sc_cls_freq;
  CFG::LevelData<CFG::FArrayBox> m_norm_momentum;

  inline bool enforceConservation()
  {
     return (m_conserve_momentum);
  }

  void addConservativeTerms( KineticSpecies& rhs_species,
                             const KineticSpecies& soln_species,
                             const LevelData<FArrayBox>& tp_rhs_coll,
                             const LevelData<FArrayBox>& delta_F,
                             const Real time );


  /// Computes the Lorentz collisions.
  /**
   * Computes the Lorentz collisions.
   */

  void testPartCollRHS(LevelData<FArrayBox>& rhs_coll,
                       const LevelData<FArrayBox>& delta_dfn,
                       const PhaseGeom& phase_geom,
                       const double mass ) const;


  /// Computes the momentum restoring term.
  /**
   * Computes the momentum restoring term.
   */

  void fieldPartCollRHS(LevelData<FArrayBox>& rhs_coll,
                        const LevelData<FArrayBox>& delta_dfn,
                        const PhaseGeom& phase_geom,
                        const CFG::LevelData<CFG::FArrayBox>& moment_rest,
                        const double mass ) const;


  /// Computes the kernels for the restoring terms
  void collKernels(LevelData<FArrayBox>& kern_moment,
                     const LevelData<FArrayBox>& test_part_coll_RHS,
                       const PhaseGeom& phase_geom,
                       const double mass ) const;


  ///Computes the kernels for the normalization factors
  void evaluateNormKern( LevelData<FArrayBox>& kern_moment_norm,
                         const PhaseGeom& phase_geom,
                         const double mass ) const;


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

  void computeSelfConsistFreq(LevelData<FArrayBox>&       a_cls_freq,
                              const LevelData<FArrayBox>& a_density,
                              const double                a_mass,
                              const double                a_charge      ) const;

  void convertToCellAverage(const PhaseGeom&  coord_sys,
                            LevelData<FArrayBox>&   dfn) const;

  void convertToCellCenters(const PhaseGeom&  coord_sys,
                            LevelData<FArrayBox>&  dfn) const;
  
  void convertToCellFaces(LevelData<FluxBox>& faceData,
                          const LevelData<FArrayBox>& cellData) const;


  /// Compute linear collision frequencies (notaion from Abel et al, 2008 (3)-(7)).
  /**
   * Computes and returns coefficients in collision operator fro maxwellian background.
   */
   double nu_D(double x) const;

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

};


#include "NamespaceFooter.H"

#endif
