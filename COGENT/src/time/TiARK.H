#ifndef _TiARK_H_
#define _TiARK_H_

#include <iostream>
#include <string>

#include "ParmParse.H"
#include "parstream.H"
#include "TimeIntegrator.H"
#include "ImplicitStageFunction.H"
#include "ImplicitStageJacobian.H"
#include "ExplicitStageFunction.H"
#include "ExplicitStageJacobian.H"
#include "NewtonSolver.H"
#include "TiDefinitions.H"
#include "TiButcherTables.H"

#define _TI_ARK_ "ark"

#include "NamespaceHeader.H"

template <class Vector, class Ops>
class TiARK : public TimeIntegrator<Vector,Ops>
{

  public:

    /// Constructor
    /**
     * Constructor: set m_is_Defined to false.
     */ 
    TiARK<Vector,Ops>() : m_is_Defined(false) 
    {
      m_Ae = m_be = m_c = m_Ai = m_bi = NULL;
      m_be_emb = m_bi_emb = NULL;
      m_binterpe = m_binterpi = NULL;
      m_rhsStage_exp = m_rhsStage_imp = NULL;
      m_rhsStage_exp_prev = m_rhsStage_imp_prev = NULL;
      m_Efunction = NULL;
      m_Ifunction = NULL;
      m_EJacobian = NULL;
      m_IJacobian = NULL;
    }

    /// Destructor
    /*
     * Clean up allocations
     */
    virtual ~TiARK();

    /// Define the specific ARK method
    /**
     * define the specific ARK method (eg, "1bee", "2a", "3", "4", etc)
     *
     * @param[in] a_name string containing the method name
    */
    virtual void define(  ParmParse&, 
                          const std::string& a_name,
                          const Vector&, 
                          Ops* const );

    /// Initialize
    /**
     * @param[in] a_time   Simulation time
     * @param[in] a_Y      Vector passed to ImplicitStageJacobian::updatePreCondMat
     */
    virtual void initialize( const Real a_time, const Vector& a_Y);

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

    virtual TimeIntegratorType type() const { return ti_imex; }

    virtual void printCounts() const
      {
        if (!procID()) {
          cout << "  Time integrator counts:-\n";
          cout << "    Time steps          : " << m_count << "\n";
          if (m_Operators->trivialSolutionOp()) {
            cout << "    Nonlinear iterations: " 
                 << m_ImplicitStageSolver.getCount()
                 << "\n";
            cout << "    Linear iterations   : " 
                 << m_ImplicitStageSolver.getLinearCount()
                 << "\n";
          } else {
            cout << "    Nonlinear iterations: " 
                 << m_ImplicitStageSolver.getCount() + m_ExplicitStageSolver.getCount()
                 << "\n";
            cout << "    Linear iterations   : " 
                 << m_ImplicitStageSolver.getLinearCount() + m_ExplicitStageSolver.getLinearCount()
                 << "\n";
          }
        }
      }

  protected:

  private:
    bool          m_is_Defined;
    bool          m_isLinear;
    bool          m_stagePredictor;
    std::string   m_name;
    int           m_nstages;
    int           m_order, m_order_emb;
    int           m_pinterp;
    Real          *m_Ae, *m_be, *m_c;
    Real          *m_Ai, *m_bi;
    Real          *m_be_emb, *m_bi_emb;
    Real          *m_binterpe, *m_binterpi;
    Vector        m_YStage, m_YPrev,
                  *m_rhsStage_exp, *m_rhsStage_imp, 
                  *m_rhsStage_exp_prev, *m_rhsStage_imp_prev,
                  m_R, m_Z, m_W;
    Ops*          m_Operators;
    Real          m_time;
    Real          m_dt;
    Real          m_epsJFNK;
    int           m_cur_step,
                  m_count;
    bool          m_use_scaled_norms;

    ExplicitStageFunction<Vector,Ops>*  m_Efunction;
    ExplicitStageJacobian<Vector,Ops>*  m_EJacobian;
    ImplicitStageFunction<Vector,Ops>*  m_Ifunction;
    ImplicitStageJacobian<Vector,Ops>*  m_IJacobian;

    NewtonSolver<Vector,Ops>            m_ImplicitStageSolver,
                                        m_ExplicitStageSolver;

    void extrapolate( Real, Vector& );
    void parseParameters( ParmParse& );
};

template <class Vector, class Ops>
void TiARK<Vector,Ops>::define( ParmParse&          a_pp, 
                                const std::string&  a_name, 
                                const Vector&       a_Y, 
                                Ops* const          a_ops )
{
  m_dt = 1.0;
  m_stagePredictor = false;
  m_epsJFNK = 1e-6;
  m_use_scaled_norms = false;

  m_name = a_name;
  TiButcherTables::getButcherTableARK(  m_name,
                                        m_nstages,
                                        m_order,
                                        m_order_emb,
                                        &m_Ai, 
                                        &m_Ae,
                                        &m_bi,
                                        &m_be,
                                        &m_c,
                                        &m_bi_emb,
                                        &m_be_emb,
                                        &m_binterpi,
                                        &m_binterpe );

  ParmParse ppARK("ark");
  parseParameters( ppARK );
  if ((m_pinterp == 0) && (m_stagePredictor)) {
    if (!procID()) {
      cout << "ARK method \"" << m_name << "\" does not have dense output. ";
      cout << "Cannot use stage predictors. Set ark.stage_predictor to false ";
      cout << "or remove it from input file (default: false).\n";
    }
    return;
  }

  /* allocate RHS */
  m_rhsStage_exp  = new Vector[m_nstages];
  m_rhsStage_imp  = new Vector[m_nstages];

  /* define the work vectors and operators */
  m_YStage.define(a_Y);
  m_R.define(a_Y);
  m_Z.define(a_Y);
  m_W.define(a_Y);
  for (int i=0; i<m_nstages; i++) {
    m_rhsStage_exp[i].define(a_Y);
    m_rhsStage_imp[i].define(a_Y);
  }
  m_Operators = a_ops;
  m_Operators->useScales(m_use_scaled_norms);

  m_isLinear = m_Operators->isLinear();

  /* define vectors for stage predictors */
  if (m_stagePredictor) {
    m_YPrev.define(a_Y);
    m_rhsStage_exp_prev  = new Vector[m_nstages];
    m_rhsStage_imp_prev  = new Vector[m_nstages];
    for (int i=0; i<m_nstages; i++) {
      m_rhsStage_exp_prev[i].define(a_Y);
      m_rhsStage_imp_prev[i].define(a_Y);
    }
  } else m_rhsStage_exp_prev = m_rhsStage_imp_prev = NULL;

  /* define string prefixes for options parsing */
  std::string isolverOutPrefix = "TiARK.i";
  std::string esolverOutPrefix = "TiARK.e";
  std::string solverOptPrefix = "ark";

  /* define the Newton solvers */
  Function<Vector,Ops>       
    **ImplicitStageSolverFunction = m_ImplicitStageSolver.getFunction();
  LinearFunction<Vector,Ops> 
    **ImplicitStageSolverJacobian = m_ImplicitStageSolver.getJacobian();

  *ImplicitStageSolverFunction = new ImplicitStageFunction<Vector,Ops>;
  *ImplicitStageSolverJacobian = new ImplicitStageJacobian<Vector,Ops>;
  
  m_Ifunction = dynamic_cast<ImplicitStageFunction<Vector,Ops>*>
                  (*ImplicitStageSolverFunction);
  m_IJacobian = dynamic_cast<ImplicitStageJacobian<Vector,Ops>*>
                  (*ImplicitStageSolverJacobian);

  m_Ifunction->define(m_R,*m_Operators);
  m_IJacobian->define(m_R,*m_Operators,isolverOutPrefix,solverOptPrefix);
  m_Ifunction->setJac(m_IJacobian);
  m_IJacobian->setIsLinear(m_isLinear);
  m_IJacobian->setJFNKEps(m_epsJFNK);

  m_ImplicitStageSolver.define(a_pp,m_R,isolverOutPrefix,solverOptPrefix);
  m_ImplicitStageSolver.setIsLinear(m_isLinear);
  
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

  m_Operators->setScaleTolerance(m_ImplicitStageSolver.getAbsTol());

  m_count = 0;

  /* done */
  m_is_Defined = true;
  if (!procID()) cout << "Time integration method: ark (" << m_name << ")\n" ;
}

template <class Vector, class Ops>
void TiARK<Vector,Ops>::initialize(const Real a_time, const Vector& a_Y)
{
  Real shift = 0;
  for (int i = 0; i < m_nstages; i++) {
    if (m_Ai[i*m_nstages+i] != 0.0) {
      shift = 1.0 / (m_Ai[i*m_nstages+i] * m_dt);
      break;
    }
  }
  CH_assert(shift > 0.0);
  m_Ifunction->setShift(shift);
  m_Ifunction->setStageTime(a_time);
  m_IJacobian->setShift(shift);
  m_IJacobian->setStageTime(a_time);
  m_IJacobian->updatePreCondMat(a_Y, -1, -1);

  if (!m_Operators->trivialSolutionOp()) {
    m_Efunction->setStageTime(a_time);
    m_EJacobian->setStageTime(a_time);
    m_EJacobian->updatePreCondMat(a_Y, -1, -1);
  }

  return;
}

template <class Vector, class Ops>
TiARK<Vector, Ops>::~TiARK()
{
  if (m_Ae) delete[] m_Ae;
  if (m_be) delete[] m_be;
  if (m_c) delete[] m_c;
  if (m_Ai) delete[] m_Ai;
  if (m_bi) delete[] m_bi;
  if (m_bi_emb) delete[] m_bi_emb;
  if (m_be_emb) delete[] m_be_emb;
  if (m_pinterp != 0) {
    delete[] m_binterpe;
    delete[] m_binterpi;
  }
  delete[] m_rhsStage_exp;
  delete[] m_rhsStage_imp;
  if (m_stagePredictor) {
    delete[] m_rhsStage_exp_prev;
    delete[] m_rhsStage_imp_prev;
  }
  delete m_Ifunction;
  delete m_IJacobian;
  if (!m_Operators->trivialSolutionOp()) {
    delete m_Efunction;
    delete m_EJacobian;
  }
}

template <class Vector, class Ops>
void TiARK<Vector, Ops>::advance( const Real& a_time, Vector& a_Y )
{
  CH_TIMERS("TiARK::advance");
  CH_TIMER("implicitStageSolver",t_implicit_stage_solver);
  CH_TIMER("explicitStageSolver",t_explicit_stage_solver);
  CH_assert(isDefined());
  CH_assert(m_time == a_time);

  m_YStage = a_Y;

  Real shift = 0;
  for (int i = 0; i < m_nstages; i++) {
    if (m_Ai[i*m_nstages+i] != 0.0) {
      shift = 1.0 / (m_Ai[i*m_nstages+i] * m_dt);
      break;
    }
  }
  CH_assert(shift > 0.0);
  m_Ifunction->setShift(shift);
  m_Ifunction->setStageTime(a_time);
  m_IJacobian->setShift(shift);
  m_IJacobian->setStageTime(a_time);
  m_IJacobian->updatePreCondMat(m_YStage, m_cur_step, -1);
  if (!m_Operators->trivialSolutionOp()) {
    m_Efunction->setStageTime(a_time);
    m_EJacobian->setStageTime(a_time);
    m_EJacobian->updatePreCondMat(m_YStage, m_cur_step, -1);
  }

  /* Assert explicit-first-stage - post-time-stage function implementation 
   * in GKOps.cpp assumes that YStage(1) = y_n! */
  CH_assert(m_Ai[0] == 0.0);

  /* Stage calculations */
  int i, j;
  for (i = 0; i < m_nstages; i++) {

    Real stage_time = m_time+m_c[i]*m_dt;

    /* explicit stage solution */
    m_Z = a_Y;
    for (j=0; j<i; j++) m_Z.increment(  m_rhsStage_exp[j],
                                        (m_dt*m_Ae[i*m_nstages+j])  );
    for (j=0; j<i; j++) m_Z.increment(  m_rhsStage_imp[j],
                                        (m_dt*m_Ai[i*m_nstages+j])  );

    if (!m_Operators->trivialSolutionOp()) {
      m_Operators->preSolutionOpEval(a_Y,m_time,chkpt_ti_advance_1);
      m_Operators->solutionOp(m_W,m_time,a_Y);
      for (j=0; j<i; j++) m_W.increment(  m_rhsStage_exp[j],
                                          (m_dt*m_Ae[i*m_nstages+j])  );
      for (j=0; j<i; j++) m_W.increment(  m_rhsStage_imp[j],
                                          (m_dt*m_Ai[i*m_nstages+j])  );
    }

    if (m_Ai[i*m_nstages+i] != 0.0) {
      
      /* implicit stage */
      Real shift = 1.0 / (m_Ai[i*m_nstages+i] * m_dt);

      m_Ifunction->setShift(shift);
      m_Ifunction->setStageTime(stage_time);

      m_IJacobian->setShift(shift);
      m_IJacobian->setStageTime(stage_time);

      /* Build precond matrix from m_YStage, which is
       * first stage - a_Y
       * subsequent stages - m_YStage from previous stage */
      m_IJacobian->updatePreCondMat(m_YStage, m_cur_step, i);

      /* RHS for stage solve */
      if (m_Operators->trivialSolutionOp()) {
        m_R = m_Z;
      } else {
        m_R = m_W;
      }
      m_R *= shift;

      /* Before the Newton solve begins, m_Z contains the
       * initial guess. If using stage predictor, initial 
       * guess is computed by extrapolating from previous step */
      if (m_stagePredictor && (m_cur_step > 0)) extrapolate(m_c[i],m_Z);

      CH_START(t_implicit_stage_solver);
      m_ImplicitStageSolver.solve(m_Z,m_R);
      CH_STOP(t_implicit_stage_solver);
      m_YStage = m_Z;

    } else if (!m_Operators->trivialSolutionOp()) {

      m_Efunction->setStageTime(stage_time);
      m_EJacobian->setStageTime(stage_time);

      /* Build precond matrix from m_YStage, which is
       * first stage - a_Y
       * subsequent stages - m_YStage from previous stage */
      m_EJacobian->updatePreCondMat(m_YStage, m_cur_step, i);

      /* RHS for stage solve */
      m_R = m_W;

      /* initial guess */
      m_Z = m_YStage;

      CH_START(t_explicit_stage_solver);
      m_ExplicitStageSolver.solve(m_Z, m_R);
      CH_STOP(t_explicit_stage_solver);
      m_YStage = m_Z;

    } else {

      m_YStage = m_Z;

    }

    m_Operators->postTimeStage(m_cur_step,stage_time,m_c[i]*m_dt,m_YStage,i);
    m_Operators->preOpEval(m_YStage,stage_time,chkpt_ti_advance_2);
    m_Operators->explicitOpImEx(m_rhsStage_exp[i],stage_time,m_YStage);
    m_Operators->implicitOpImEx(m_rhsStage_imp[i],stage_time,m_YStage);
  }

  /* save stuff for stage predictor */
  if (m_stagePredictor) {
    m_YPrev.copy(a_Y);
    for (i = 0; i< m_nstages; i++) {
      m_rhsStage_exp_prev[i].copy(m_rhsStage_exp[i]);
      m_rhsStage_imp_prev[i].copy(m_rhsStage_imp[i]);
    }
  }

  /* Step completion */
  if (m_Operators->trivialSolutionOp()) {

    for (i = 0; i < m_nstages; i++) {
      a_Y.increment(m_rhsStage_exp[i],(m_dt*m_be[i]));
      a_Y.increment(m_rhsStage_imp[i],(m_dt*m_bi[i]));
    }

  } else {

    m_Efunction->setStageTime(m_time+m_dt);
    m_EJacobian->setStageTime(m_time+m_dt);
    m_EJacobian->updatePreCondMat(a_Y, m_cur_step, m_nstages);

    m_Operators->preSolutionOpEval(a_Y,(m_time+m_dt),chkpt_ti_advance_3);
    m_Operators->solutionOp(m_R,(m_time+m_dt),a_Y);
    for (i = 0; i < m_nstages; i++) {
      m_R.increment(m_rhsStage_exp[i],(m_dt*m_be[i]));
      m_R.increment(m_rhsStage_imp[i],(m_dt*m_bi[i]));
    }

    CH_START(t_explicit_stage_solver);
    m_ExplicitStageSolver.solve(a_Y,m_R);
    CH_STOP(t_explicit_stage_solver);
  }

  /* update current time and step number */
  m_cur_step++;
  m_time += m_dt; 
  m_count++;
}

/* Use dense output to extrapolate from previous stage */
template <class Vector, class Ops>
void TiARK<Vector, Ops>::extrapolate(Real c, Vector& a_Y)
{
  Real h, tt, t, *bt,*b;
  int i,j;
  if (!m_binterpi || !m_binterpe) {
    if (!procID()) cout << "Error: ARK method \"" << m_name << "\" does not have an interpolation formula.\n";
    return;
  }

  t = 1.0 + c;
  h = m_dt;
  bt = (Real*) calloc (m_nstages,sizeof(Real));
  b  = (Real*) calloc (m_nstages,sizeof(Real));

  for (i=0; i<m_nstages; i++) bt[i] = b[i] = 0;
  for (j=0,tt=t; j<m_pinterp; j++,tt*=t) {
    for (i=0; i<m_nstages; i++) {
      bt[i] += h * m_binterpi[i*m_pinterp+j] * tt;
      b[i]  += h * m_binterpe[i*m_pinterp+j] * tt;
    }
  }

  a_Y.copy(m_YPrev);
  for (i=0; i<m_nstages; i++) {
    a_Y.increment(m_rhsStage_imp_prev[i],bt[i]);
    a_Y.increment(m_rhsStage_exp_prev[i],b[i] );
  }

  free(bt);
  free(b);
}

template <class Vector, class Ops>
void TiARK<Vector,Ops>::parseParameters( ParmParse& a_pp)
{
  /* use stage predictor based on dense output */
  a_pp.query("stage_predictor", m_stagePredictor); 
  a_pp.query("jfnk_epsilon", m_epsJFNK); 
  a_pp.query("scaled_norms", m_use_scaled_norms); 
}

#include "NamespaceFooter.H"

#endif
