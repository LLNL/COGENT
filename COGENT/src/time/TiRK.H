#ifndef _TiRK_H_
#define _TiRK_H_

#include <iostream>
#include <string>

#include "TimeIntegrator.H"
#include "ExplicitStageFunction.H"
#include "ExplicitStageJacobian.H"
#include "NewtonSolver.H"
#include "TiDefinitions.H"
#include "TiButcherTables.H"

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
    TiRK<Vector,Ops>() : m_is_Defined(false)
    {
      m_A = m_b = m_c = m_b_emb = NULL;
      m_rhsStage = NULL;
      m_Efunction = NULL;
      m_EJacobian = NULL;
    }

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
    virtual void define(  ParmParse&,
                          const std::string&,
                          const Vector&,
                          Ops* const );

    /// Initialize
    virtual void initialize( const Real, const Vector& );

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
    virtual Ops& getOperators() { return *m_Operators; }

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

    virtual TimeIntegratorType type() const { return ti_explicit; }

    virtual void printCounts() const
      {
        if (!procID()) {
          cout << "  Time integrator counts:-\n";
          cout << "    Time steps          : " << m_count << "\n";
          if (!m_Operators->trivialSolutionOp()) {
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
    int         m_order, m_order_emb;
    Real        *m_A, *m_b, *m_c, *m_b_emb;
    Vector      m_YStage, *m_rhsStage, m_R;
    Ops*        m_Operators;
    Real        m_time;
    Real        m_dt;
    Real        m_epsJFNK;
    int         m_cur_step, m_count;

    ExplicitStageFunction<Vector,Ops>  *m_Efunction;
    ExplicitStageJacobian<Vector,Ops>  *m_EJacobian;
    NewtonSolver<Vector,Ops>           m_ExplicitStageSolver;

};

template <class Vector, class Ops>
void TiRK<Vector,Ops>::define( ParmParse&         a_pp,
                               const std::string& a_name,
                               const Vector&      a_Y,
                               Ops* const         a_ops )
{
  m_dt = 1.0;
  m_epsJFNK = 1e-6;

  m_name = a_name;
  TiButcherTables::getButcherTableRK( a_name,
                                      m_nstages,
                                      m_order,
                                      m_order_emb,
                                      &m_A,
                                      &m_b,
                                      &m_c,
                                      &m_b_emb );

  /* allocate RHS */
  m_rhsStage  = new Vector[m_nstages];

  m_YStage.define(a_Y);
  m_R.define(a_Y);
  for (int i=0; i<m_nstages; i++) {
    m_rhsStage[i].define(a_Y);
  }
  m_Operators = a_ops;
  m_isLinear = false;

  /* define string prefixes for options parsing */
  std::string esolverOutPrefix = "TiRK.e";
  std::string solverOptPrefix = "rk";

  if (!m_Operators->trivialSolutionOp()) {
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

    m_Efunction->define(m_R,*m_Operators);
    m_EJacobian->define(m_R,*m_Operators,esolverOutPrefix,solverOptPrefix);
    m_Efunction->setJac(m_EJacobian);
    m_EJacobian->setIsLinear(m_isLinear);
    m_EJacobian->setJFNKEps(m_epsJFNK);

    m_ExplicitStageSolver.define(a_pp,m_R,esolverOutPrefix,solverOptPrefix);
    m_ExplicitStageSolver.setIsLinear(m_isLinear);
  }

  m_count = 0;

  /* done */
  m_is_Defined = true;
  if (!procID()) cout << "Time integration method: rk (" << m_name << ")\n" ;
}

template <class Vector, class Ops>
void TiRK<Vector,Ops>::initialize(const Real a_time, const Vector& a_Y)
{
  if (!m_Operators->trivialSolutionOp()) {
    m_Efunction->setStageTime(a_time);
    m_EJacobian->setStageTime(a_time);
    m_EJacobian->updatePreCondMat(a_Y, -1, -1);
  }
  return;
}

template <class Vector, class Ops>
TiRK<Vector, Ops>::~TiRK()
{
  if (m_A) delete[] m_A;
  if (m_b) delete[] m_b;
  if (m_c) delete[] m_c;
  if (m_b_emb) delete[] m_b_emb;
  if (m_rhsStage) delete[] m_rhsStage;
  if (!m_Operators->trivialSolutionOp()) {
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

  if (!m_Operators->trivialSolutionOp()) {
    m_Efunction->setStageTime(a_time);
    m_EJacobian->setStageTime(a_time);
    m_EJacobian->updatePreCondMat(m_YStage, m_cur_step, -1);
  }

  /* Stage calculations */

  int i, j;
  for (i = 0; i < m_nstages; i++) {

    Real stage_time = m_time+m_c[i]*m_dt;

    if (m_Operators->trivialSolutionOp()) {

      m_YStage.copy(a_Y);
      for (j=0; j<i; j++) m_YStage.increment(m_rhsStage[j],(m_dt*m_A[i*m_nstages+j]));

    } else {

      m_Efunction->setStageTime(stage_time);
      m_EJacobian->setStageTime(stage_time);

      /* first stage - a_Y
       * subsequent stages - m_YStage from previous stage */
      m_EJacobian->updatePreCondMat(m_YStage, m_cur_step, i);

      m_Operators->preSolutionOpEval(a_Y,m_time);
      m_Operators->solutionOp(m_R,stage_time,a_Y);
      for (j=0; j<i; j++) m_R.increment(  m_rhsStage[j],
                                          (m_dt*m_A[i*m_nstages+j]) );
      m_ExplicitStageSolver.solve(m_YStage,m_R);

    }

    m_Operators->postTimeStage(m_cur_step,stage_time,m_c[i]*m_dt,m_YStage,i);
    m_Operators->preOpEval(m_YStage,stage_time,chkpt_ti_advance_2);
    m_Operators->explicitOp(m_rhsStage[i],stage_time,m_YStage);
  }

  /* Step completion */
  if (m_Operators->trivialSolutionOp()) {

    for (i = 0; i < m_nstages; i++) {
      a_Y.increment(m_rhsStage[i],(m_dt*m_b[i]));
    }

  } else {

    m_Efunction->setStageTime(m_time+m_dt);
    m_EJacobian->setStageTime(m_time+m_dt);
    m_EJacobian->updatePreCondMat(a_Y, m_cur_step, m_nstages);

    m_Operators->preSolutionOpEval(a_Y,(m_time+m_dt));
    m_Operators->solutionOp(m_R,(m_time+m_dt),a_Y);
    for (i = 0; i < m_nstages; i++) {
      m_R.increment(m_rhsStage[i],(m_dt*m_b[i]));
    }

    m_ExplicitStageSolver.solve(a_Y,m_R);

  }

  /* update current time and step number */
  m_cur_step++;
  m_time += m_dt;
  m_count++;
}

#include "NamespaceFooter.H"

#endif
