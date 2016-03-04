#ifndef _TIME_INTEGRATOR_H_
#define _TIME_INTEGRATOR_H_

#include <string>

#include "REAL.H"
#include "ParmParse.H"
#include "parstream.H"

#include "NamespaceHeader.H"

#define _TI_EXTERNAL_ -1
#define _TI_EXPLICIT_ 0
#define _TI_IMEX_     1

template <class Solution, class RHS, class Ops>
class TimeIntegrator
{

  public:
    virtual Ops&  getOperators() = 0;
    virtual void  define(ParmParse& a_pp,std::string a_name, Solution& a_state, Real a_dt) = 0;
    virtual void  advance(Real a_time, Solution& Y) = 0;

    virtual void setTimeStep(int a_n)         = 0;
    virtual void getTimeStep(int *a_n)        = 0;
    virtual void setTimeStepSize(Real a_dt)   = 0;
    virtual void getTimeStepSize(Real *a_dt)  = 0;
    virtual void setCurrentTime(Real a_time)  = 0;
    virtual void getCurrentTime(Real *a_time) = 0;

    virtual ~TimeIntegrator<Solution,RHS,Ops>() {}

    virtual short int getType() { return _TI_EXPLICIT_; }
    virtual void getCounts(int *a,int *b,int *c) = 0;
};

#include "NamespaceFooter.H"

#endif
