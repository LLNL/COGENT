#ifndef _TiARK_H_
#define _TiARK_H_

#include <iostream>
#include <string>

#include "REAL.H"
#include "TimeIntegrator.H"

#include "NamespaceHeader.H"

template <class Solution, class RHS, class Ops>
class TiARK : public TimeIntegrator<Solution,RHS,Ops>
{

  public:

    /// Constructor
    /**
     * Constructor: set m_is_Defined to false.
     */ 
    TiARK<Solution,RHS,Ops>() { m_is_Defined = false; }

    /// Destructor
    /*
     * Clean up allocations
     */
    ~TiARK();

    /// Define the specific ARK method
    /**
     * define the specific ARK method (eg, "1bee", "2a", "3", "4", etc)
     *
     * @param[in] a_name string containing the method name
     */
    void define(std::string a_name, Solution& a_state, Real a_dt);

    /// Advance one time step
    /**
     * Advance one time step.
     *
     * @param[in] a_dt time step
     * @param[in] a_time current simulation time
     * @param[out] a_Y solution
     */
    void advance(Real a_time, Solution& a_Y);

    /// Check if method is defined
    /**
     * Returns the value of m_is_Defined
     */
    bool isDefined() const { return m_is_Defined; }
    
    /// Get the operators for the time integrator
    /**
     * get the operators for the time integrator
     */
    Ops& getOperators() { return m_Operators; }

    /// Set the time step size
    /**
     * set the time step size for the time integrator
     *
     * @param[in] a_dt the specified time step size
     */
    void setTimeStep(Real a_dt) { m_dt = a_dt; }

    /// Get the time step size
    /**
     * get the time step size of the time integrator
     *
     * @param[out] a_dt the time step size
     */ 
    void getTimeStep(Real *a_dt) { *a_dt = m_dt; }

    /// Set the current simulation time 
    /**
     * set the current simulation time
     *
     * @param[in] a_time the specified simulation time
     */ 
    void setCurrentTime(Real a_time) { m_time = a_time; }

    /// Get the current simulation time
    /*
     * get the current simulation time
     *
     * @param[out] a_time the current simulation time
     */ 
    void getCurrentTime(Real *a_time) { *a_time = m_time; }

  protected:

  private:
    bool        m_is_Defined;
    std::string m_name;
    int         m_nstages;
    Real        *m_Ae, *m_be, *m_ce;
    Real        *m_Ai, *m_bi, *m_ci;
    Solution    m_YStage;
    RHS         *m_rhsStage_exp, *m_rhsStage_imp, m_work;
    Ops         m_Operators;
    Real        m_time;
    Real        m_dt;
};


template <class Solution, class RHS, class Ops>
void TiARK<Solution,RHS,Ops>::define(std::string a_name, Solution& a_state, Real a_dt)
{
  m_dt = a_dt;

  if (a_name == "1bee") {
    m_name = a_name;
    m_nstages = 3;
  } else if (a_name == "rk4") {
    /* explicit RK4, disguised as an ARK method -- for debugging */
    m_name = a_name;
    m_nstages = 4;
  } else {
    /* default: forward Euler */
    m_name = "1bee";
    m_nstages = 3;
    //MayDay.Warning("Warning: Unrecognized ARK method. Defaulting to 1bee.");
  }

  /* allocate Butcher tableaux coefficients */
  m_Ae = new Real[m_nstages*m_nstages];
  m_be = new Real[m_nstages];
  m_ce = new Real[m_nstages];
  m_Ai = new Real[m_nstages*m_nstages];
  m_bi = new Real[m_nstages];
  m_ci = new Real[m_nstages];

  /* allocate RHS */
  m_rhsStage_exp  = new RHS[m_nstages];
  m_rhsStage_imp  = new RHS[m_nstages];

  /* set Butcher coefficients */
  for (int i=0; i<m_nstages*m_nstages; i++) m_Ae[i] = m_Ai[i] = 0; 
  for (int i=0; i<m_nstages; i++) m_be[i] = m_bi[i] = m_ce[i] = m_ci[i] = 0; 
  if (m_name == "1bee") {

    m_Ae[7] = 0.5;

    m_Ai[0] = 1.0;
    m_Ai[4] = 0.5;
    m_Ai[7] = 0.5;
    m_Ai[8] = 0.5;

    m_ce[0] = 0.0;
    m_ce[1] = 0.0;
    m_ce[2] = 0.5;

    m_ci[0] = 1.0;
    m_ci[1] = 0.5;
    m_ci[2] = 1.0;

    m_be[1] = m_be[2] = m_bi[1] = m_bi[2]  = 0.5;

  } else if (m_name == "rk4") {

    m_Ae[4]  = 0.5;
    m_Ae[9]  = 0.5;
    m_Ae[14] = 1.0;

    m_Ai[4]  = 0.5;
    m_Ai[9]  = 0.5;
    m_Ai[14] = 1.0;

    m_ce[1] = m_ce[2] = 0.5;
    m_ce[3] = 1.0;

    m_ci[1] = m_ci[2] = 0.5;
    m_ci[3] = 1.0;

    m_be[0] = 1.0/6.0; 
    m_be[1] = m_be[2] = 1.0/3.0; 
    m_be[3] = 1.0/6.0;

    m_bi[0] = 1.0/6.0; 
    m_bi[1] = m_bi[2] = 1.0/3.0; 
    m_bi[3] = 1.0/6.0;

  }

  m_YStage.define(a_state);
  m_work.define(a_state);
  for (int i=0; i<m_nstages; i++) {
    m_rhsStage_exp[i].define(a_state);
    m_rhsStage_imp[i].define(a_state);
  }
  m_Operators.define(a_state, m_dt);
  m_is_Defined = true;

  if (!procID()) cout << "Time integration method: ark (" << m_name << ")\n" ;
}

template <class Solution, class RHS, class Ops>
TiARK<Solution, RHS, Ops>::~TiARK()
{
  delete[] m_Ae;
  delete[] m_be;
  delete[] m_ce;
  delete[] m_Ai;
  delete[] m_bi;
  delete[] m_ci;
  delete[] m_rhsStage_exp;
  delete[] m_rhsStage_imp;
}

template <class Solution, class RHS, class Ops>
void TiARK<Solution, RHS, Ops>::advance(Real a_time, Solution& a_Y)
{
  //CH_TIMERS("TiARK::advance");
  CH_assert(isDefined());
  CH_assert(m_time == a_time);

  /* Stage calculations */
  int i, j;
  for (i = 0; i < m_nstages; i++) {
    m_work.copy(a_Y);
    for (j=0; j<i; j++) m_work.increment(m_rhsStage_exp[j],(m_dt*m_Ae[i*m_nstages+j]));
    for (j=0; j<i; j++) m_work.increment(m_rhsStage_imp[j],(m_dt*m_Ai[i*m_nstages+j]));

    /* implicit stage */
    Real shift = m_Ai[i*m_nstages+j] * m_dt;
    if (shift != 0.0) m_Operators.solveImEx(m_YStage,m_work,shift,(m_time+m_ci[i]*m_dt),i);

    m_Operators.explicitOpImEx(m_rhsStage_exp[i],(m_time+m_ce[i]*m_dt),m_YStage,i);
    m_Operators.implicitOpImEx(m_rhsStage_imp[i],(m_time+m_ci[i]*m_dt),m_YStage,i);
  }
  /* Step completion */
  for (i = 0; i < m_nstages; i++) {
    a_Y.increment(m_rhsStage_exp[i],(m_dt*m_be[i]));
    a_Y.increment(m_rhsStage_imp[i],(m_dt*m_bi[i]));
  }

  /* update current time */
  m_time += m_dt; 
}

#include "NamespaceFooter.H"

#endif
