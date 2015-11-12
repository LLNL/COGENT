#ifndef  _CONSDRAGDIFF_H_
#define  _CONSDRAGDIFF_H_

#include "KineticSpecies.H"
#include "REAL.H"
#include "CLSInterface.H"
#include "ParmParse.H"
#include <sstream>

#include "NamespaceHeader.H"

class Kernel;
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
                ParmParse&     a_ppcls,
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
  virtual void evalClsRHS( KineticSpeciesPtrVect&        rhs,
                           const KineticSpeciesPtrVect&  soln,
                           const int                     species,
                           const Real                    time );

  /// Compute a stable time step.
  /**
   * Computes and returns an estimate of the maximum stable time step.
   */
  Real computeDt( const KineticSpeciesPtrVect&  soln);


private:

  Real cls_freq;
  bool cls_only;
  string species_name;
  MomentOp& moment_op;

  int verbosity;
  ParmParse ppcls;

  /// Parse parameters.
  /**
   * Private method to obtain control parameters from "CLS.species" section
   * of ParmParse database, where "species" is replaced by the species name
   *
   * \param[in] ppcls CLS.species input database.
   */
  void ParseParameters();

  /// Print parameters.
  /**
   * Private method to print parameters.
   */
  void printParameters( const KineticSpeciesPtrVect&  soln );

};


#include "NamespaceFooter.H"

#endif
