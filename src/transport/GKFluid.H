#ifndef  _GKFLUID_H_
#define  _GKFLUID_H_

#include "KineticSpecies.H"
#include "REAL.H"
#include "TPMInterface.H"
#include "ParmParse.H"
#include <sstream>

#include "NamespaceHeader.H"

class Kernel;
class MomentOp;

/**
 * Anomalous transport operator class.
*/
class GKFluid
   : public TPMInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  GKFluid( const string& a_species_name,
         ParmParse& a_pptpm,
         const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~GKFluid();

  /// Accumulates the RHS of the tranpsort operator.
  /**
   * @param[out] rhs data holder for rhs (without ghost cells).
   * @param[in] soln current solution with which to evaluate RHS.
   * @param[in] species index of species to which tranpsport is applied.
   * @param[in] time current time.
   */
  virtual void evalTpmRHS( KineticSpeciesPtrVect&       rhs,
                           const KineticSpeciesPtrVect& soln,
                           const int                    species,
                           const Real                   time );

private:

  Vector<Real> D_fluid;
  Vector<Real> D_kinet;
  bool model_only;
  MomentOp& moment_op;
  bool m_first_step;

  int verbosity;
  ParmParse pptpm;
  string species_name;

  LevelData<FluxBox> metrics_faces;

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
  void printParameters(const KineticSpeciesPtrVect& soln);

  /// Get Metrics hr, htheta, and hphi
  /**
   * Private method to get metric components on cell faces.
   */
  void getFaceCenteredMetrics(LevelData<FluxBox>&                 metrics_faces,
                              const PhaseGeom&         a_phase_geom,
                              const DisjointBoxLayout& a_dbl);

  /// Get Metrics hr, htheta, and hphi
  /**
   * Private method to get metric components on cell centers.
   */
  void getCellCenteredMetrics(FArrayBox&               metrics_cells,
                              const PhaseGeom&         a_phase_geom,
                              const Box& a_dbl);

};


#include "NamespaceFooter.H"

#endif
