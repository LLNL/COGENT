#ifndef _TiRK_H_
#define _TiRK_H_

#include <string>

#include "REAL.H"
#include "TimeIntegrator.H"

#include "NamespaceHeader.H"

template <class Solution, class RHS, class Ops>
class TiRK : public TimeIntegrator<Solution,RHS,Ops>
{

  public:

    /// Constructor
    /**
     * Constructor: set m_is_Defined to false.
     */ 
    TiRK<Solution,RHS,Ops>() { m_is_Defined = false; }

    /// Destructor
    /*
     * Clean up allocations
     */
    ~TiRK();

    /// Define the specific RK method
    /**
     * define the specific RK method (eg, "1fe", "2a", "3", "4", etc)
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
     *
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
    Real        *m_A, *m_b, *m_c;
    Solution    m_YStage;
    RHS         *m_rhsStage;
    Ops         m_Operators;
    Real        m_time;
    Real        m_dt;
};


template <class Solution, class RHS, class Ops>
void TiRK<Solution,RHS,Ops>::define(std::string a_name, Solution& a_state, Real a_dt)
{
  m_dt = a_dt;

  if (a_name == "1fe") {
    /* forward Euler */
    m_name = a_name;
    m_nstages = 1;
  } else if (a_name == "2a") {
    /* 2nd order, 2-stage Runge-Kutta */
    m_name = a_name;
    m_nstages = 2;
  } else if (a_name == "3") {
    /* 3rd order, 3-stage Runge-Kutta */
    m_name = a_name;
    m_nstages = 3;
  } else if (a_name == "4") {
    /* 4th order, 4-stage Runge-Kutta */
    m_name = a_name;
    m_nstages = 4;
  } else {
    /* default: forward Euler */
    m_name = "1fe";
    m_nstages = 1;
    //MayDay.Warning("Warning: Unrecognized RK method. Defaulting to forward Euler.");
  }

  /* allocate Butcher tableaux coefficients */
  m_A = new Real[m_nstages*m_nstages];
  m_b = new Real[m_nstages];
  m_c = new Real[m_nstages];

  /* allocate RHS */
  m_rhsStage  = new RHS[m_nstages];

  /* set Butcher coefficients */
  for (int i=0; i<m_nstages*m_nstages; i++) m_A[i] = 0; 
  for (int i=0; i<m_nstages; i++) m_b[i] = m_c[i] = 0; 
  if (m_name == "1fe") {
    m_b[0] = 1.0;
  } else if (m_name == "2a") {
    m_A[2] = 1.0;
    m_c[1] = 1.0;
    m_b[0] = m_b[1] = 0.5;
  } else if (m_name == "3") {
    m_A[3] = 2.0/3.0;
    m_A[6] = 2.0/3.0 - 1.0/4.0;
    m_A[7] = 1.0/4.0;
    m_c[1] = 2.0/3.0;
    m_c[2] = 2.0/3.0;
    m_b[0] = 1.0/4.0;
    m_b[1] = -1.0/4.0;
    m_b[2] = 1.0;
  } else if (m_name == "4") {
    m_A[4] = 0.5;
    m_A[9] = 0.5;
    m_A[14] = 1.0;
    m_c[1] = m_c[2] = 0.5;
    m_c[3] = 1.0;
    m_b[0] = 1.0/6.0; 
    m_b[1] = m_b[2] = 1.0/3.0; 
    m_b[3] = 1.0/6.0;
  }

  m_YStage.define(a_state);
  for (int i=0; i<m_nstages; i++) {
    m_rhsStage[i].define(a_state);
  }
  m_Operators.define(a_state, m_dt);
  m_is_Defined = true;
}

template <class Solution, class RHS, class Ops>
TiRK<Solution, RHS, Ops>::~TiRK()
{
  delete[] m_A;
  delete[] m_b;
  delete[] m_c;
  delete[] m_rhsStage;
}

template <class Solution, class RHS, class Ops>
void TiRK<Solution, RHS, Ops>::advance(Real a_time, Solution& a_Y)
{
  //CH_TIMERS("TiRK::advance");
  CH_assert(isDefined());
  CH_assert(m_time == a_time);

  /* Stage calculations */
  int i, j;
  for (i = 0; i < m_nstages; i++) {
    m_YStage.copy(a_Y);
    for (j=0; j<i; j++) m_YStage.increment(m_rhsStage[j],(m_dt*m_A[i*m_nstages+j]));
    m_Operators.explicitOp(m_rhsStage[i],(m_time+m_c[i]*m_dt),m_YStage,i);
  }
  /* Step completion */
  for (i = 0; i < m_nstages; i++) a_Y.increment(m_rhsStage[i],(m_dt*m_b[i]));

  /* update current time */
  m_time += m_dt; 
}

#include "NamespaceFooter.H"

#endif
