#ifndef  _CONSDRAGDIFF_H_
#define  _CONSDRAGDIFF_H_

#include "KineticSpecies.H"
#include "REAL.H"
#include "CLSInterface.H"
#include "ParmParse.H"
#include "BasicGKPreconditioner.H"
#include "GlobalDOF.H"
#include <sstream>

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunction.H"
#include "GridFunctionLibrary.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

class MomentOp;

/**
 * Fully conservative drag-diffusion collisions operator class.
*/
class ConsDragDiff
   : public CLSInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  ConsDragDiff( const string&  a_species_name,
                const string&  a_ppcls_str,
                const int      a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~ConsDragDiff();

  /// Accumulates the RHS of the collision operator.
  /**
   * @param[out] rhs data holder for rhs (without ghost cells).
   * @param[in] soln current solution with which to evaluate RHS.
   * @param[in] species index of species to which collision is applied.
   * @param[in] time current time.
   */
  /// Compute a stable time step.
  /**
   * Computes and returns an estimate of the maximum stable time step.
   */
  Real computeDtExplicitTI( const KineticSpeciesPtrVect&  soln, const int species);
  Real computeDtImExTI( const KineticSpeciesPtrVect&  soln, const int species);
  Real computeTimeScale( const KineticSpeciesPtrVect&  soln, const int species);

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

  inline int precondMatrixBands()
  {
    if (m_time_implicit)  return(m_nbands);
    else                  return(0);
  }

  void defineBlockPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>&,
                      std::vector<DOFList>&,
                      const ODEVector&,
                      void*,
                      const std::string&,
                      const std::string&,
                      bool,
                      const KineticSpecies&,
                      const GlobalDOFKineticSpecies&,
                      const int,
                      const int );

  virtual void updateBlockPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>&,
                              const KineticSpecies&,
                              const GlobalDOFKineticSpecies&,
                              const Real,
                              const int,
                              const int,
                              const Real,
                              const bool,
                              const int );

  void assemblePrecondMatrix( BandedMatrix&,
                              const KineticSpecies&,
                              const GlobalDOFKineticSpecies&,
                              const int,
                              const int,
                              const Real );
  
  
private:

  
  void computeVelFluxesDiv(LevelData<FArrayBox>&        Jpsi,
                           const LevelData<FArrayBox>&  fBJ,
                           const PhaseGeom&             phase_geom,
                           const double                 mass);
  
  void computeConsUparTemp(LevelData<FArrayBox>& Jpsi,
                           const KineticSpecies& soln_species,
                           const PhaseGeom&      phase_geom );
  
  void computeClsFreq(const KineticSpecies&     soln_species,
                      const PhaseGeom&          phase_geom,
                      const DisjointBoxLayout&  grids,
                      const Real                time);
  
  void computeSelfConsistFreq(LevelData<FArrayBox>&       cls_freq,
                              const LevelData<FArrayBox>& density,
                              const LevelData<FArrayBox>& temperature,
                              const double                mass,
                              const double                charge      ) const;

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
  void ParseParameters(const ParmParse& a_pp);

  /// Print parameters.
  /**
   * Private method to print parameters.
   */
  void printParameters(const KineticSpeciesPtrVect& soln,
                       const int                    species);
  
  
  void diagnostics(const LevelData<FArrayBox>& rhs,
                   const KineticSpecies&       rhs_species,
                   const double                time) const;
  
  RefCountedPtr<CFG::GridFunction> m_cls_freq_func;

  LevelData<FArrayBox> m_cls_freq;
  LevelData<FArrayBox> m_fBJ_vel_ghost;
  LevelData<FArrayBox> m_Jpsi;
  LevelData<FArrayBox> m_rhs_cls;
  
  LevelData<FluxBox> m_fluxes;
  
  LevelData<FArrayBox> m_inj_vpar_moms;
  LevelData<FArrayBox> m_inj_pres_moms;
  
  LevelData<FArrayBox> m_dens_inj;
  LevelData<FArrayBox> m_temp_inj;
  
  LevelData<FArrayBox> m_Upar;
  LevelData<FArrayBox> m_Temp;

  string m_species_name;
  MomentOp& moment_op;

  bool m_cls_only;
  bool m_first_call;
  bool m_skip_stage_update;
  int m_it_counter;
  int m_update_cls_freq;
  int m_update_operator_coeff;
  bool m_time_implicit;
  bool m_diagnostics;

  int m_verbosity;
  ParmParse m_ppcls;
  std::string m_opt_string;
  
  int m_nbands;
  int m_my_pc_idx;

};


#include "NamespaceFooter.H"

#endif
