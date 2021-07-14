#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "REAL.H"

#include "NamespaceHeader.H"

template <class AppVec, class System>
class Simulation
{
  public: 
    Simulation<AppVec, System>() {}
    virtual ~Simulation<AppVec, System>() {}

    virtual void initialize(ParmParse&, System* const) = 0;
    virtual void solve() = 0;
    virtual void finalize() = 0;
};

#include "NamespaceFooter.H"
#endif
