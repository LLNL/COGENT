#ifndef  _ANOMALOUS_H_
#define  _ANOMALOUS_H_

#include "KineticSpecies.H"
#include "REAL.H"
#include "TPMInterface.H"
#include "ParmParse.H"
#include <sstream>

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunction.H"
#include "GridFunctionLibrary.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

class Kernel;
class MomentOp;

/**
 * Anomalous transport operator class.
 * See Justin Angus' PhD thesis Chapter 5 for detailed explaination
 * and verification tests
*/
class Anomalous
   : public TPMInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  Anomalous( const string& a_species_name,
         ParmParse& a_pptpm,
         const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~Anomalous();

  /// Accumulates the RHS of the tranpsort operator.
  /**
   * @param[out] rhs data holder for rhs (without ghost cells).
   * @param[in] soln current solution with which to evaluate RHS.
   * @param[in] species index of species to which tranpsport is applied.
   * @param[in] time current time.
   */
  void evalTpmRHS( KineticSpecies&              rhs_species,
                           const KineticSpeciesPtrVect& soln,
                           const int                    species,
                           const Real                   time );

  /// Compute a stable time step.
  /**
  * Computes and returns an estimate of the maximum stable time step.
  */
  Real computeStableDt( const KineticSpecies& soln_species );
   
   
  const CFG::LevelData<CFG::FluxBox>& getPrecondCoeff() const {return m_precond_D;}
   
   
private:
  MomentOp& m_moment_op;
  Vector<Real> m_D_fluid;
  Vector<Real> m_D_kinet;
  
  double m_DN0;
   
  RefCountedPtr<CFG::GridFunction> m_shape_function;
  RefCountedPtr<CFG::GridFunction> m_shape_function_D0, m_shape_function_D1;
  RefCountedPtr<CFG::GridFunction> m_shape_function_D2, m_shape_function_D3;
  
  bool m_model_only;
  bool m_const_coeff;
  bool m_first_call;
  bool m_first_stage;
  bool m_arbitrary_grid;
  bool m_simple_diffusion;
  int m_update_freq;
  int m_it_counter;
  int m_verbosity;

   
  IntVect m_ghostVect;
  ParmParse pptpm;
  string species_name;

  LevelData<FArrayBox> m_fB, m_fB_beta, m_dfBdmu_cc, m_rhs_transport;
  LevelData<FArrayBox> m_hr, m_beta;
  LevelData<FluxBox> m_lame_faces, m_inj_bunit, m_inj_pointwiseNJinv;
  LevelData<FluxBox> m_D_kinet_faces;
  LevelData<FluxBox> m_flux, m_fluxNorm;

  LevelData<FArrayBox> m_density, m_density_beta;
  LevelData<FArrayBox> m_Upar, m_Upar_beta;
  LevelData<FArrayBox> m_temperature, m_temperature_beta;
  LevelData<FArrayBox> m_fourth_coef, m_fourth_coef_beta;
  LevelData<FArrayBox> m_perp_coef, m_perp_coef_beta;

  CFG::LevelData<CFG::FluxBox> m_precond_D;


  /// Parse parameters.
  /**
   * Private method to obtain control parameters from "TPM.species" section
   * of ParmParse database, where "species" is replaced by the species name
   *
   * \param[in] pptpm TPM.species input database.
   */
  void ParseParameters();

  /// Print parameters.
  /**
   * Private method to print parameters.
   */
  void printParameters(const KineticSpecies& soln_species);

  /**
   * Private method to get lame coefficients on cell faces.
   */
  void getFaceCenteredLameCoefficients( LevelData<FluxBox>&      a_lame_faces,
                                        const PhaseGeom&         a_phase_geom,
                                        const DisjointBoxLayout& a_dbl);

  /**
   * Private method to get lame coefficients on cell centers.
   */
  void getCellCenteredLameCoefficients( FArrayBox&             a_lame_cells,
                                        const PhaseGeom&       a_phase_geom,
                                        const Box&             a_dbl);


  /**
   * Private method for setting the spatial dependence of the transport coefficients
   */
  void setTransportCoeffSpatialDependence( CFG::LevelData<CFG::FArrayBox>&   a_D_kinet_cfg,
                                           const CFG::MagGeom&             a_mag_geom,
                                           const Real                      a_time );

  
  /// Compute the preconditioner coeffficient for the implicit solver
  void computePrecondCoefficient( CFG::LevelData<CFG::FluxBox>&            a_D,
                                  const CFG::MagGeom&                      a_phase_geom,
                                  const CFG::LevelData<CFG::FluxBox>&      a_D_faces );
 
   /// Compute cell-centered beta = U_psi^2/D_psi
   void computeBeta(LevelData<FArrayBox>& beta,
                    const KineticSpecies& soln_species );
   
   /// Pre time step calculation
   inline void preTimeStep  ( const KineticSpeciesPtrVect& soln_mapped,
                              const int species,
                              const Real time,
                              const KineticSpeciesPtrVect& soln_physical);

   
};


#include "NamespaceFooter.H"

#endif
