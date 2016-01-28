#ifndef _TiARK_H_
#define _TiARK_H_

#include <iostream>
#include <string>

#include "TimeIntegrator.H"
#include "ImplicitStageFunction.H"
#include "ImplicitStageJacobian.H"
#include "NewtonSolver.H"

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
    void define(ParmParse& a_pp, std::string a_name, Solution& a_state, Real a_dt);

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
    void setTimeStepSize(Real a_dt) { m_dt = a_dt; }

    /// Get the time step size
    /**
     * get the time step size of the time integrator
     *
     * @param[out] a_dt the time step size
     */ 
    void getTimeStepSize(Real *a_dt) { *a_dt = m_dt; }

    /// Set the time step 
    /**
     * set the time step for the time integrator
     *
     * @param[in] a_n the specified time step
     */
    void setTimeStep(int a_n) { m_cur_step = a_n; }

    /// Get the time step
    /**
     * get the time step of the time integrator
     *
     * @param[out] a_n the time step
     */ 
    void getTimeStep(int *a_n) { *a_n = m_cur_step; }

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

    /// Set Butcher tableaux coefficients
    /*
     * set the Butcher tableaux coefficients for the method
     *
     * @param[in] a_nstages Number of stages
     * @param[in] a_Ae Explicit stage calculation coefficients
     * @param[in] a_Ai Implicit stage calculation coefficients
     * @param[in] a_be Explicit step completion coefficients (if NULL, equal to a_be)
     * @param[in] a_bi Implicit step completion coefficients (if NULL, last row of a_Ai) 
     */ 
    void setCoefficients(int a_nstages, const Real* a_Ae, const Real* a_Ai, const Real* a_be, const Real* a_bi);

  protected:

  private:
    bool          m_is_Defined;
    std::string   m_name;
    int           m_nstages;
    Real          *m_Ae, *m_be, *m_ce;
    Real          *m_Ai, *m_bi, *m_ci;
    Solution      m_YStage;
    RHS           *m_rhsStage_exp, *m_rhsStage_imp, m_R, m_Z;
    Ops           m_Operators;
    Real          m_time;
    Real          m_dt;
    int           m_cur_step;

    ImplicitStageFunction<RHS,Ops>  *m_Ifunction;
    ImplicitStageJacobian<RHS,Ops>  *m_IJacobian;
    NewtonSolver<RHS,Ops>           m_NewtonSolver;
};


template <class Solution, class RHS, class Ops>
void TiARK<Solution,RHS,Ops>::define(ParmParse& a_pp, std::string a_name, Solution& a_state, Real a_dt)
{
  m_dt = a_dt;

  /* set Butcher coefficients */
  if (a_name == "1bee") {

    m_name = a_name;
    m_nstages = 3;

    const Real
      Ae[3][3] = {{0.0,0.0,0.0},
                  {0.0,0.0,0.0},
                  {0.0,0.5,0.0}},
      Ai[3][3] = {{1.0,0.0,0.0},
                  {0.0,0.5,0.0},
                  {0.0,0.5,0.5}},
      be[3]    =  {0.0,0.5,0.5};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],&be[0],NULL);

  } else if (a_name == "4") {

    /* Fourth order ARK scheme with L-stable implicit part.
       This method has one explicit stage and four implicit stages.
       References: Kennedy and Carpenter 2003.
    */
    m_name = a_name;
    m_nstages = 6;

    const Real 
      Ae[6][6] = {{0,0,0,0,0,0},
                 {1./2,0,0,0,0,0},
                 {13861./62500.,6889./62500.,0,0,0,0},
                 {-116923316275./2393684061468.,-2731218467317./15368042101831.,9408046702089./11113171139209.,0,0,0},
                 {-451086348788./2902428689909.,-2682348792572./7519795681897.,12662868775082./11960479115383.,3355817975965./11060851509271.,0,0},
                 {647845179188./3216320057751.,73281519250./8382639484533.,552539513391./3454668386233.,3354512671639./8306763924573.,4040./17871.,0}},
      Ai[6][6] = {{0,0,0,0,0,0},
                  {1./4,1./4,0,0,0,0},
                  {8611./62500.,-1743./31250.,1./4,0,0,0},
                  {5012029./34652500.,-654441./2922500.,174375./388108.,1./4,0,0},
                  {15267082809./155376265600.,-71443401./120774400.,730878875./902184768.,2285395./8070912.,1./4,0},
                  {82889./524892.,0,15625./83664.,69875./102672.,-2260./8211,1./4}};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],NULL,NULL);

  } else if (a_name == "rk4") {

    /* explicit RK4, disguised as an ARK method -- for debugging */
    m_name = a_name;
    m_nstages = 4;

    const Real
      A[4][4] = {{0,0,0,0},
                 {0.5,0,0,0},
                 {0,0.5,0,0},
                 {0,0,1.0,0}},
      b[4]    = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
    setCoefficients(m_nstages,&A[0][0],&A[0][0],&b[0],&b[0]);

  } else {

    /* default: ARK4 */
    if (!procID()) cout << "Error: unknown ARK method specified " << a_name << ". Using default.\n";
    m_name = "4";
    m_nstages = 6;

    const Real 
      Ae[6][6] = {{0,0,0,0,0,0},
                 {1./2,0,0,0,0,0},
                 {13861./62500.,6889./62500.,0,0,0,0},
                 {-116923316275./2393684061468.,-2731218467317./15368042101831.,9408046702089./11113171139209.,0,0,0},
                 {-451086348788./2902428689909.,-2682348792572./7519795681897.,12662868775082./11960479115383.,3355817975965./11060851509271.,0,0},
                 {647845179188./3216320057751.,73281519250./8382639484533.,552539513391./3454668386233.,3354512671639./8306763924573.,4040./17871.,0}},
      Ai[6][6] = {{0,0,0,0,0,0},
                  {1./4,1./4,0,0,0,0},
                  {8611./62500.,-1743./31250.,1./4,0,0,0},
                  {5012029./34652500.,-654441./2922500.,174375./388108.,1./4,0,0},
                  {15267082809./155376265600.,-71443401./120774400.,730878875./902184768.,2285395./8070912.,1./4,0},
                  {82889./524892.,0,15625./83664.,69875./102672.,-2260./8211,1./4}};
    setCoefficients(m_nstages,&Ae[0][0],&Ai[0][0],NULL,NULL);

  }

  /* allocate RHS */
  m_rhsStage_exp  = new RHS[m_nstages];
  m_rhsStage_imp  = new RHS[m_nstages];

  /* define the work vectors and operators */
  m_YStage.define(a_state);
  m_R.define(a_state);
  m_Z.define(a_state);
  for (int i=0; i<m_nstages; i++) {
    m_rhsStage_exp[i].define(a_state);
    m_rhsStage_imp[i].define(a_state);
  }
  m_Operators.define(a_state, m_dt);

  /* define the Newton solver */
  Function<RHS,Ops>       **NewtonSolverFunction = m_NewtonSolver.getFunction();
  LinearFunction<RHS,Ops> **NewtonSolverJacobian = m_NewtonSolver.getJacobian();
  *NewtonSolverFunction = new ImplicitStageFunction<RHS,Ops>;
  *NewtonSolverJacobian = new ImplicitStageJacobian<RHS,Ops>;
  
  m_Ifunction = dynamic_cast<ImplicitStageFunction<RHS,Ops>*>(*NewtonSolverFunction);
  m_IJacobian = dynamic_cast<ImplicitStageJacobian<RHS,Ops>*>(*NewtonSolverJacobian);
  m_Ifunction->define(m_R,m_Operators);
  m_IJacobian->define(m_R,m_Operators);
  m_Ifunction->setJac(m_IJacobian);

  std::string solverOutPrefix = "TiARK";
  std::string solverOptPrefix = "ark";
  m_NewtonSolver.define(a_pp,m_R,solverOutPrefix,solverOptPrefix);

  /* done */
  m_is_Defined = true;
  if (!procID()) cout << "Time integration method: ark (" << m_name << ")\n" ;
}

template <class Solution, class RHS, class Ops>
void TiARK<Solution, RHS, Ops>::setCoefficients( int a_nstages,
                                                 const Real* a_Ae, const Real* a_Ai, 
                                                 const Real* a_be, const Real* a_bi
                                               )
{
  CH_assert(!isDefined());
  CH_assert(a_nstages == m_nstages);

  /* allocate Butcher tableaux coefficients 
   * deallocated in destructor */
  m_Ae = new Real[m_nstages*m_nstages];
  m_be = new Real[m_nstages];
  m_ce = new Real[m_nstages];
  m_Ai = new Real[m_nstages*m_nstages];
  m_bi = new Real[m_nstages];
  m_ci = new Real[m_nstages];

  int i, j;
  for (i=0; i<m_nstages*m_nstages; i++) {
    m_Ae[i] = a_Ae[i];
    m_Ai[i] = a_Ai[i];
  }
  if (a_bi) for (i=0; i<m_nstages; i++) m_bi[i] = a_bi[i];
  else      for (i=0; i<m_nstages; i++) m_bi[i] = a_Ai[(m_nstages-1)*m_nstages+i];
  if (a_be) for (i=0; i<m_nstages; i++) m_be[i] = a_be[i];
  else      for (i=0; i<m_nstages; i++) m_be[i] = m_bi[i];

  for (i=0; i<m_nstages; i++) {
    m_ce[i] = m_ci[i] = 0.0;
    for(j=0; j<m_nstages; j++) { 
      m_ce[i] += a_Ae[i*m_nstages+j]; 
      m_ci[i] += a_Ai[i*m_nstages+j]; 
    }
  }
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
  delete m_Ifunction;
  delete m_IJacobian;
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
    m_YStage.copy(a_Y);
    for (j=0; j<i; j++) m_YStage.increment(m_rhsStage_exp[j],(m_dt*m_Ae[i*m_nstages+j]));
    for (j=0; j<i; j++) m_YStage.increment(m_rhsStage_imp[j],(m_dt*m_Ai[i*m_nstages+j]));

    /* implicit stage */
    Real shift = m_Ai[i*m_nstages+j] * m_dt;
    if (shift != 0.0) {
      m_Ifunction->setShift(shift);
      m_Ifunction->setStageNumber(i);
      m_Ifunction->setStageTime(m_time+m_ce[i]*m_dt);
      m_IJacobian->setShift(shift);
      m_IJacobian->setStageNumber(i);
      m_IJacobian->setStageTime(m_time+m_ce[i]*m_dt);

      m_R.copy(m_YStage);  /* right-hand side */
      m_Z.copy(a_Y);       /* initial guess */
      m_NewtonSolver.solve(m_Z,m_R);
      m_YStage.copy(m_Z);
    }

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
