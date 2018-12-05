#ifndef _TiRK_H_
#define _TiRK_H_

#include <iostream>
#include <string>

#include "TimeIntegrator.H"
#include "ExplicitStageFunction.H"
#include "ExplicitStageJacobian.H"
#include "NewtonSolver.H"

#define _TI_RK_ "rk"

#include "NamespaceHeader.H"

template <class Vector, class Ops>
class TiRK : public TimeIntegrator<Vector,Ops>
{

  public:

    /// Constructor
    /**
     * Constructor: set m_is_Defined to false.
     */ 
   TiRK<Vector,Ops>() : m_is_Defined(false) {}

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
    virtual void define( ParmParse&, std::string, Vector&, Real );

    /// Advance one time step
    /**
     * Advance one time step.
     *
     * @param[in] a_dt time step
     * @param[in] a_time current simulation time
     * @param[out] a_Y solution
     */
    virtual void advance( const Real& a_time, Vector& a_Y );

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

    virtual void printCounts() const
      {
        if (!procID()) {
          cout << "  Time integrator counts:-\n";
          cout << "    Time steps          : " << m_count << "\n";
          if (!m_Operators.trivialSolutionOp()) {
            cout << "    Nonlinear iterations: " 
                 << m_ExplicitStageSolver.getCount()
                 << "\n";
            cout << "    Linear iterations   : " 
                 << m_ExplicitStageSolver.getLinearCount()
                 << "\n";
          }
        }
      }


private:
    bool        m_is_Defined;
    bool        m_isLinear;
    std::string m_name;
    int         m_nstages;
    Real        *m_A, *m_b, *m_c;
    Vector      m_YStage, *m_rhsStage, m_R;
    Ops         m_Operators;
    Real        m_time;
    Real        m_dt;
    Real        m_epsJFNK;
    int         m_cur_step, m_count;
    bool        m_step_const_coef;

    ExplicitStageFunction<Vector,Ops>  *m_Efunction;
    ExplicitStageJacobian<Vector,Ops>  *m_EJacobian;
    NewtonSolver<Vector,Ops>           m_ExplicitStageSolver;

};

template <class Vector, class Ops>
void TiRK<Vector,Ops>::define( ParmParse& a_pp,
                               std::string a_name, Vector& a_Y,
                               Real a_dt )
{
  m_dt = a_dt;
  m_epsJFNK = 1e-6;
  m_step_const_coef = false;
  
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

  } else if (a_name == "5dp") {

    /* 5th order, 7-stage Dormand-Prince */
    m_name = a_name;
    m_nstages = 7;

    const Real
      A[7][7]   = {{0,0,0,0,0,0,0},
                   {1.0/5.0,0,0,0,0,0,0},
                   {3.0/40.0,9.0/40.0,0,0,0,0,0},
                   {44.0/45.0,-56.0/15.0,32.0/9.0,0,0,0,0},
                   {19372.0/6561.0,-25360.0/2187.0,64448.0/6561.0,-212.0/729.0,0,0,0},
                   {9017.0/3168.0,-355.0/33.0,46732.0/5247.0,49.0/176.0,-5103.0/18656.0,0,0},
                   {35.0/384.0,0,500.0/1113.0,125.0/192.0,-2187.0/6784.0,11.0/84.0,0}},
      b[7]      = {35.0/384.0,0,500.0/1113.0,125.0/192.0,-2187.0/6784.0,11.0/84.0,0};
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
  m_rhsStage  = new Vector[m_nstages];

  m_YStage.define(a_Y);
  m_R.define(a_Y);
  for (int i=0; i<m_nstages; i++) {
    m_rhsStage[i].define(a_Y);
  }
  m_Operators.define((*a_Y.getState()), m_dt);
  m_isLinear = false; 

  /* define string prefixes for options parsing */
  std::string esolverOutPrefix = "TiRK.e";
  std::string solverOptPrefix = "rk";

  if (!m_Operators.trivialSolutionOp()) {
    Function<Vector,Ops>       
      **ExplicitStageSolverFunction = m_ExplicitStageSolver.getFunction();
    LinearFunction<Vector,Ops> 
      **ExplicitStageSolverJacobian = m_ExplicitStageSolver.getJacobian();
  
    *ExplicitStageSolverFunction = new ExplicitStageFunction<Vector,Ops>;
    *ExplicitStageSolverJacobian = new ExplicitStageJacobian<Vector,Ops>;
    
    m_Efunction = dynamic_cast<ExplicitStageFunction<Vector,Ops>*>
                    (*ExplicitStageSolverFunction);
    m_EJacobian = dynamic_cast<ExplicitStageJacobian<Vector,Ops>*>
                    (*ExplicitStageSolverJacobian);
  
    m_Efunction->define(m_R,m_Operators);
    m_EJacobian->define(m_R,m_Operators,esolverOutPrefix,solverOptPrefix);
    m_Efunction->setJac(m_EJacobian);
    m_EJacobian->setIsLinear(m_isLinear);
    m_EJacobian->setJFNKEps(m_epsJFNK);
  
    m_ExplicitStageSolver.define(a_pp,m_R,esolverOutPrefix,solverOptPrefix);
    m_ExplicitStageSolver.setIsLinear(m_isLinear);
  }

  m_count = 0;

  a_pp.query("step_const_coef", m_step_const_coef); 

  /* done */
  m_is_Defined = true;
  if (!procID()) cout << "Time integration method: rk (" << m_name << ")\n" ;
}

template <class Vector, class Ops>
void TiRK<Vector, Ops>::setCoefficients( int a_nstages,
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

template <class Vector, class Ops>
TiRK<Vector, Ops>::~TiRK()
{
  delete[] m_A;
  delete[] m_b;
  delete[] m_c;
  delete[] m_rhsStage;
  if (!m_Operators.trivialSolutionOp()) {
    delete m_Efunction;
    delete m_EJacobian;
  }
}

template <class Vector, class Ops>
void TiRK<Vector, Ops>::advance( const Real& a_time, Vector& a_Y )
{
  CH_TIMERS("TiRK::advance");
  CH_assert(isDefined());
  CH_assert(m_time == a_time);
  m_YStage.copy(a_Y);

  /* Stage calculations */
  int i, j;
  for (i = 0; i < m_nstages; i++) {

    Real stage_time = m_time+m_c[i]*m_dt;

    if (m_Operators.trivialSolutionOp()) {
      if ( m_step_const_coef && i==0 ) {
         // Initialize the coefficients for the right-hand side operators,
         // which are assumed to be constant over the step
         m_Operators.preOpEval(a_Y,m_time);
      }
      m_YStage.copy(a_Y);
      for (j=0; j<i; j++) m_YStage.increment(m_rhsStage[j],(m_dt*m_A[i*m_nstages+j]));

    } else {

      m_Efunction->setStageTime(stage_time);
      m_EJacobian->setStageTime(stage_time);
  
      /* first stage - a_Y
       * subsequent stages - m_YStage from previous stage */
      m_EJacobian->updatePreCondMat(m_YStage);
  
      if ( m_step_const_coef ) {
         if ( i==0 ) {
            // Initialize the coefficients for both the solution and right-hand side
            // operators, which are assumed to be constant over the step
            m_Operators.preOpEval(a_Y,stage_time);
         }
      }
      else {
         m_Operators.preSolutionOpEval(a_Y,stage_time);
      }
      m_Operators.solutionOp(m_R,stage_time,a_Y);
      for (j=0; j<i; j++) m_R.increment(m_rhsStage[j],(m_dt*m_A[i*m_nstages+j]));
  
      m_ExplicitStageSolver.solve(m_YStage,m_R,!m_step_const_coef);

    }

    m_Operators.postTimeStage(m_cur_step,stage_time,m_YStage,i);
    if ( !m_step_const_coef ) m_Operators.preOpEval(m_YStage,stage_time);
    m_Operators.explicitOp(m_rhsStage[i],stage_time,m_YStage);
  }

  /* Step completion */

  if (m_Operators.trivialSolutionOp()) {

    for (i = 0; i < m_nstages; i++) {
      a_Y.increment(m_rhsStage[i],(m_dt*m_b[i]));
    }

  } else {

    m_Efunction->setStageTime(m_time+m_dt);
    m_EJacobian->setStageTime(m_time+m_dt);
    m_EJacobian->updatePreCondMat(a_Y);

    if ( !m_step_const_coef ) m_Operators.preSolutionOpEval(a_Y,(m_time+m_dt));
    m_Operators.solutionOp(m_R,(m_time+m_dt),a_Y);
    for (i = 0; i < m_nstages; i++) {
      m_R.increment(m_rhsStage[i],(m_dt*m_b[i]));
    }

    m_ExplicitStageSolver.solve(a_Y,m_R,!m_step_const_coef);

  }

  /* update current time and step number */
  m_cur_step++;
  m_time += m_dt; 
  m_count++;
}

#include "NamespaceFooter.H"

#endif
