#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "REAL.H"

#include "NamespaceHeader.H"

/**
 * Interface which represents time integrators
 *
 * @tparam AppVec   A vector storing the ODE system state
 *
 * @tparam System   The type of the ODE system to solve.
 *                  See class AppCtxt for the interface definition.
 */
template <class AppVec, class System>
class Simulation
{
  public: 
    Simulation<AppVec, System>() {}
    virtual ~Simulation<AppVec, System>() {}

    /// Set up the time integrator.
    /// This must be called once before `solve() is called.
    virtual void initialize(ParmParse&, System* const) = 0;

    /// Run the simulation
    virtual void solve() = 0;

    /// Tidy up before deleting
    virtual void finalize() = 0;
};

#include "NamespaceFooter.H"
#endif
