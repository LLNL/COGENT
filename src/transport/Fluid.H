#ifndef  _FLUID_H_
#define  _FLUID_H_

#include "KineticSpecies.H"
#include "REAL.H"
#include "TPMInterface.H"
#include "ParmParse.H"
#include <sstream>

#include "NamespaceHeader.H"

/**
 * Fluid transport operator class.
*/
class Fluid
   : public TPMInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  Fluid( const string& a_species_name,
         ParmParse& a_pptpm,
         const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~Fluid();

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

  Vector<Real> D_m;
  bool m_model_only;

  int m_verbosity;
  ParmParse m_pptpm;
  string m_species_name;
  bool m_first_step;

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
  void getFaceCenteredMetrics(FluxBox&                 metrics_faces,
                              const PhaseGeom&         a_phase_geom,
                              const Box& a_dbl);

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
