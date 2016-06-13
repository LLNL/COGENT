#ifndef _TiRK_H_
#define _TiRK_H_

#include <iostream>
#include <string>

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
   TiRK<Solution,RHS,Ops>() : m_is_Defined(false) {}

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
    virtual void define( ParmParse& a_pp, std::string a_name, Solution& a_state, Real a_dt );

    /// Advance one time step
    /**
     * Advance one time step.
     *
     * @param[in] a_dt time step
     * @param[in] a_time current simulation time
     * @param[out] a_Y solution
     */
    virtual void advance( const Real& a_time, Solution& a_Y );

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
    virtual Ops& getOperators() { return m_Operators; }

    /// Set the time step size
    /**
     * set the time step size for the time integrator
     *
     * @param[in] a_dt the specified time step size
     */
    virtual void setTimeStepSize( const Real& a_dt ) { m_dt = a_dt; }

    /// Get the time step size
    /**
     * get the time step size of the time integrator
     *
     * @param[out] a_dt the time step size
     */ 
    virtual void getTimeStepSize( Real& a_dt ) const { a_dt = m_dt; }

    /// Set the time step
    /**
     * set the time step for the time integrator
     *
     * @param[in] a_n the specified time step
     */
    virtual void setTimeStep( const int& a_n ) { m_cur_step = a_n; }

    /// Get the time step
    /**
     * get the time step of the time integrator
     *
     * @param[out] a_n the time step
     */ 
    virtual void getTimeStep( int& a_n ) const { a_n = m_cur_step; }

    /// Set the current simulation time 
    /**
     * set the current simulation time
     *
     * @param[in] a_time the specified simulation time
     */ 
    virtual void setCurrentTime( const Real& a_time ) { m_time = a_time; }

    /// Get the current simulation time
    /*
     * get the current simulation time
     *
     * @param[out] a_time the current simulation time
     */ 
    virtual void getCurrentTime( Real& a_time ) const { a_time = m_time; }

    /// Set Butcher tableaux coefficients
    /*
     * set the Butcher tableaux coefficients for the method
     *
     * @param[in] a_nstages Number of stages
     * @param[in] a_A Stage calculation coefficients
     * @param[in] a_b Step completion coefficients
     */ 
    virtual void setCoefficients(int a_nstages, const Real* a_A, const Real* a_b);

    virtual bool isExplicit() const { return true; }
   
    virtual bool isImEx() const { return false; }

    virtual bool isExternal() const { return false; }

    virtual void printCounts() const
      {
        if (!procID()) {
          cout << "  Time integrator counts:-\n";
          cout << "    Time steps: " << m_count << "\n";
        }
      }

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
    int         m_cur_step, m_count;
};

template <class Solution, class RHS, class Ops>
void TiRK<Solution,RHS,Ops>::define( ParmParse& a_pp,
                                     std::string a_name, Solution& a_state,
                                     Real a_dt )
{
  m_dt = a_dt;

  /* set Butcher coefficients */
  if (a_name == "1fe") {

    /* forward Euler */
    m_name = a_name;
    m_nstages = 1;

    const Real
      A[1][1] = {{0.0}},
      b[1]    = {1.0};
    setCoefficients(m_nstages,&A[0][0],&b[0]);

  } else if (a_name == "2a") {

    /* 2nd order, 2-stage Runge-Kutta */
    m_name = a_name;
    m_nstages = 2;
    
    const Real
      A[2][2]     = {{0.0,0.0},
                    {1.0,0.0}},
      b[2]        = {0.5,0.5};
    setCoefficients(m_nstages,&A[0][0],&b[0]); 

  } else if (a_name == "3") {

    /* 3rd order, 3-stage Runge-Kutta */
    m_name = a_name;
    m_nstages = 3;
    
    const Real
      A[3][3] = {{0,0,0},
                 {2.0/3.0,0,0},
                 {-1.0/3.0,1.0,0}},
      b[3]    = {0.25,0.5,0.25};
    setCoefficients(m_nstages,&A[0][0],&b[0]); 

  } else if (a_name == "4") {
    
    /* 4th order, 4-stage Runge-Kutta */
    m_name = a_name;
    m_nstages = 4;
    
    const Real
      A[4][4] = {{0,0,0,0},
                 {0.5,0,0,0},
                 {0,0.5,0,0},
                 {0,0,1.0,0}},
      b[4]    = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
    setCoefficients(m_nstages,&A[0][0],&b[0]); 

  } else {
    
    /* default: RK4 */
    if (!procID()) cout << "Warning: unknown RK method specified " << a_name << ". Using default.\n";
    m_name = "4";
    m_nstages = 4;
    
    const Real
      A[4][4] = {{0,0,0,0},
                 {0.5,0,0,0},
                 {0,0.5,0,0},
                 {0,0,1.0,0}},
      b[4]    = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
    setCoefficients(m_nstages,&A[0][0],&b[0]); 

  }

  /* allocate RHS */
  m_rhsStage  = new RHS[m_nstages];

  m_YStage.define(a_state);
  for (int i=0; i<m_nstages; i++) {
    m_rhsStage[i].define(a_state);
  }
  m_Operators.define(a_state, m_dt);
  m_count = 0;
  m_is_Defined = true;

  if (!procID()) cout << "Time integration method: rk (" << m_name << ")\n" ;
}

template <class Solution, class RHS, class Ops>
void TiRK<Solution, RHS, Ops>::setCoefficients( int a_nstages,
                                                const Real* a_A, 
                                                const Real* a_b
                                              )
{
  CH_assert(!isDefined());
  CH_assert(a_nstages == m_nstages);

  /* allocate Butcher tableaux coefficients 
   * deallocated in destructor */
  m_A = new Real[m_nstages*m_nstages];
  m_b = new Real[m_nstages];
  m_c = new Real[m_nstages];

  int i, j;
  for (i=0; i<m_nstages*m_nstages; i++) m_A[i] = a_A[i];
  for (i=0; i<m_nstages;           i++) m_b[i] = a_b[i];
  for (i=0; i<m_nstages; i++) {
    m_c[i] = 0.0; for(j=0; j<m_nstages; j++) m_c[i] += a_A[i*m_nstages+j]; 
  }
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
void TiRK<Solution, RHS, Ops>::advance( const Real& a_time, Solution& a_Y )
{
  //CH_TIMERS("TiRK::advance");
  CH_assert(isDefined());
  CH_assert(m_time == a_time);

  /* Stage calculations */
  int i, j;
  for (i = 0; i < m_nstages; i++) {
    m_YStage.copy(a_Y);
    for (j=0; j<i; j++) m_YStage.increment(m_rhsStage[j],(m_dt*m_A[i*m_nstages+j]));
    Real stage_time = m_time+m_c[i]*m_dt;
    m_Operators.postTimeStage(m_cur_step,stage_time,m_YStage,i);
    m_Operators.explicitOp(m_rhsStage[i],stage_time,m_YStage,i);
  }
  /* Step completion */
  for (i = 0; i < m_nstages; i++) a_Y.increment(m_rhsStage[i],(m_dt*m_b[i]));

  /* update current time and step number */
  m_cur_step++;
  m_time += m_dt; 
  m_count++;
}

#include "NamespaceFooter.H"

#endif
