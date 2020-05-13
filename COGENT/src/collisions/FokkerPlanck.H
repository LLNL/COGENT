#ifndef  _FOKKERPLANCK_H_
#define  _FOKKERPLANCK_H_

#include "BandedMatrix.H"
#include "FArrayBox.H"
#include "FaceBox.H"
#include "CornerBox.H"
#include "LevelData.H"
#include "RefCountedPtr.H"
#include "Kernels.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "REAL.H"
#include "GlobalDOF.H"
#include "CLSInterface.H"
#include "ParmParse.H"
#include "BasicGKPreconditioner.H"

#include <sstream>
#include <iostream>
#include <string>

#include "NamespaceHeader.H"

#if CH_SPACEDIM==5
#define _VPAR_DIM_ 3
#define _MU_DIM_ 4
#else
#define _VPAR_DIM_ 2
#define _MU_DIM_ 3
#endif

typedef struct fokkerplankcoeffs
{

   LevelData<FaceBox>   D_v, D_mu, D_vv, D_mumu;
   LevelData<CornerBox> D_vmu;

} FokkerPlanckCoeffs;

/**
 * FokkerPlanck collisions operator class.
*/
class FokkerPlanck
   : public CLSInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  FokkerPlanck( const std::string&  a_ppcls_str,
                const int           a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~FokkerPlanck();

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
  Real computeDtExplicitTI(const KineticSpeciesPtrVect&, const int);
  Real computeDtImExTI(const KineticSpeciesPtrVect&, const int);
  Real computeTimeScale(const KineticSpeciesPtrVect&, const int);

  /// Report that collision operator is nonlinear
  /**
   * Reports that this collision operator is nonlinear
   */
  inline bool isLinear() { return m_rosenbluth_skip_stage; }

  inline int precondMatrixBands() 
  { 
    if (m_time_implicit)  return(m_nbands); 
    else                  return(0);
  }

  void defineBlockPC( std::vector<Preconditioner<GKVector,GKOps>*>&,
                      std::vector<DOFList>&,
                      const GKVector&,
                      GKOps&,
                      const std::string&,
                      const std::string&,
                      bool,
                      const KineticSpecies&,
                      const GlobalDOFKineticSpecies&,
                      int );

  virtual void updateBlockPC( std::vector<Preconditioner<GKVector,GKOps>*>&,
                              const KineticSpecies&,
                              const GlobalDOFKineticSpecies&,
                              const Real,
                              const Real,
                              const bool,
                              const int );

  void assemblePrecondMatrix( void*,
                              const KineticSpecies&,
                              const GlobalDOFKineticSpecies&,
                              const Real );

  inline void preTimeStep  ( const KineticSpeciesPtrVect& a_soln_mapped, 
                             const int a_species, 
                             const Real a_time,
                             const KineticSpeciesPtrVect& a_soln_physical);

  inline void postTimeStage( const KineticSpeciesPtrVect& a_soln, 
                             const int a_species, 
                             const Real a_time,
                             const int a_stage)
    { 
      return;
    }

private:

   bool m_verbosity;
   bool m_time_implicit;
   bool m_limiters;
   ParmParse m_ppcls;

   bool m_fixed_cls_freq;
   Real m_cls_freq;
   Real m_cls_norm;
   Real m_pcg_tol;
   Real m_pcg_maxiter;
   int  m_mult_num;
   int  m_nD;

   int  m_update_freq;
   int  m_it_counter;

   LevelData<FArrayBox> m_phi;
   FokkerPlanckCoeffs m_coeffs;

   int m_nbands;
   int m_my_pc_idx;

   bool m_rosenbluth_skip_stage;

   bool m_subtract_background;
   bool m_ref_updated;
   bool m_compute_maxwellian;
   bool m_conserve_energy;
   Real m_fp_energy_cons_epsilon;
   RefCountedPtr<KineticFunction> m_ref_func;

   LevelData<FArrayBox> m_F0;
   LevelData<FArrayBox> m_phi_F0;
   FokkerPlanckCoeffs   m_coeffs_F0;

   CFG::LevelData<CFG::FArrayBox> m_fp_kinetic_energy;
   CFG::LevelData<CFG::FArrayBox> m_energy_cons;

   std::string m_opt_string;

   /// Computes the Rosenbluth potentials
   void evalRosenbluthPotentials( LevelData<FArrayBox>& phi,
                                  const RefCountedPtr<PhaseGeom>& phase_geom,
                                  const LevelData<FArrayBox>& dfn,
                                  const double mass ) const;

  ///Computes the normalization factor 
  /**
   * Computes the normalization factor 
   * 
   * cls_norm = dfn_norm[s^3/m^6] * tau_norm[s] * 
   *            Coulomb_Lg * a_charge^4 * ech[SI]^4 / (a_mass * M[kg] * eps0 )^2
   *
   * dfn_norm = (1/Pi) * N[1/m^3] / (T[J]/M[kg])^(3/2)
   *
   * tau_norm = L[m] / Sqrt(T[J]/M[kg])
   * 
   */
  void computeClsNorm(Real&                       a_cls_norm,
                      const double                a_mass,
                      const double                a_charge      ) const;


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

   void convertToCellCenters(const PhaseGeom&, LevelData<FArrayBox>&) const;

   void convertToCellAverages(const PhaseGeom&, LevelData<FArrayBox>&) const;

   void fillGhostCellsPhi(const PhaseGeom&, LevelData<FArrayBox>&) const;

   void fillGhostCellsDfn     (const PhaseGeom&, LevelData<FArrayBox>&) const;
   void fillGhostCellsDfnVPar (const PhaseGeom&, LevelData<FArrayBox>&) const;
   void fillGhostCellsDfnMu   (const PhaseGeom&, LevelData<FArrayBox>&) const;

   void computeReferenceSolution(const KineticSpeciesPtrVect&, const int, const Real);

   void evalCoefficients( FokkerPlanckCoeffs&, 
                          const LevelData<FArrayBox>&,
                          const PhaseGeom&,
                          const Real, const Real) const;

   void computePotentialsAndCoeffs( const KineticSpecies&,
                                    const LevelData<FArrayBox>&, 
                                    const DisjointBoxLayout&,
                                    const RefCountedPtr<PhaseGeom>& );

   void computeFokkerPlanckFlux( LevelData<FaceBox>&,
                                 LevelData<FaceBox>&,
                                 const FokkerPlanckCoeffs&,
                                 const LevelData<FArrayBox>&, 
                                 const PhaseGeom& );

   void computeFokkerPlanckRHS( LevelData<FArrayBox>&,
                                const LevelData<FaceBox>&,
                                const LevelData<FaceBox>&,
                                const PhaseGeom& );

   void computeEnergyConservationFactor(const KineticSpecies&,
                                        const LevelData<FaceBox>&,
                                        const LevelData<FaceBox>&,
                                        const DisjointBoxLayout&,
                                        const PhaseGeom& );

};

#include "NamespaceFooter.H"

#endif
