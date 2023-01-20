#ifndef _NewtonSolver_H_
#define _NewtonSolver_H_

#include <iostream>
#include <string>

#include "REAL.H"
#include "ParmParse.H"
#include "parstream.H"
#include "Function.H"
#include "LinearFunction.H"
#include "BiCGStab.H"
#include "GMRES.H"

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class NewtonSolver
{
  public:

    /// Constructor
    /**
     * Constructor: set m_is_Defined to false.
     */ 
    NewtonSolver<T,Ops>()
       : m_is_Defined(false),
         m_norm(DBL_MAX),
         m_norm0(DBL_MAX),
         m_norm_prev(DBL_MAX),
         m_Function(NULL),
         m_Jacobian(NULL),
         m_LinearSolver(NULL) {};

    /// Destructor
    /*
     * Clean up allocations
     */
    ~NewtonSolver();

    /// Check if object is defined
    /**
     * Returns the value of m_is_Defined
     */
    bool isDefined() const { return m_is_Defined; }

    /// Define the Newton solver
    /**
     * define and set up the Newton solver
     *
     * @param[in] a_state solution object
     */
    void define(ParmParse&, T&, std::string, std::string);

    /// Get the number of iterations
    /**
     * get the number of Newton iterations
     */
    int getConvergedIterations() const { return m_its; }

    /// Solve the equation using Newton's method
    /**
     * solve F(x) = b using Newton's method
     *
     * @param[in] a_Y initial guess; contains the solution at exit
     * @param[in] a_b right-hand side
     */
    void solve(T& a_Y, T& a_b);

    /// Get converged reason
    /**
     * Get the reason for Newton solver exit
     */
    std::string getConvergedReason() const {return m_convergedReason; }

    /// Get exit status
    /*
     * Get the exit status for the Newton solve
     */
    int getExitStatus() const { return m_exitStatus; }

    /// Get the function object
    /**
     * Get the pointer to the function object
     */
    Function<T,Ops>** getFunction() { return &m_Function; }

    /// Get the Jacobian object
    /**
     * Get the pointer to the Jacobian object
     */
    LinearFunction<T,Ops>** getJacobian() { return &m_Jacobian; }

    inline void setIsLinear(bool a_flag) { m_isLinear = a_flag; }

    inline int getCount() const { return m_count; }
    inline int getLinearCount() const { return m_LinearSolver->getCount(); }

    inline Real getAbsTol() const { return m_atol; }
    inline Real getRelTol() const { return m_rtol; }

  protected:

  private:

    void parseParameters( ParmParse& );

    bool      m_is_Defined;
    int       m_its, m_maxits, m_exitStatus, m_count;
    Real      m_norm, m_norm0, m_norm_prev,
              m_rtol, m_atol, m_stol;
    T         m_dY, m_F;
    bool      m_verbose, m_isLinear;

    std::string           m_convergedReason,
                          m_outPrefix,
                          m_optPrefix,
                          m_linear_solver;

    Function<T,Ops>             *m_Function;
    LinearFunction<T,Ops>       *m_Jacobian;
    LinearFunctionSolver<T,Ops> *m_LinearSolver;
};

template <class T, class Ops>
NewtonSolver<T,Ops>::~NewtonSolver()
{
  /*
   * The calling module using this solver is responsible for
   * deallocating the objects pointed to by m_Function and
   * m_Jacobian (since it will also create these objects).
   * m_LinearSolver is, however, owned by this class
   */

   if (m_LinearSolver) delete m_LinearSolver;
   
  return;
}

template <class T, class Ops>
void NewtonSolver<T,Ops>::define(ParmParse& a_pp, T& a_state, std::string a_outPrefix, std::string a_optPrefix)
{
  CH_assert(!isDefined());

  /* default values */
  m_maxits  = 100;
  m_rtol    = 1e-6;
  m_atol    = 1e-10;
  m_stol    = 1e-14;
  m_verbose = true;
  m_count   = 0;
  m_isLinear = false;
  m_linear_solver = "GMRES";

  m_outPrefix = a_outPrefix;
  m_optPrefix = a_optPrefix;
  std::string optString = "_newton";
  optString = m_optPrefix + optString;

  ParmParse ppNewton( optString.c_str() );
  parseParameters( ppNewton );

  if ( m_linear_solver == "GMRES" || 
       m_linear_solver == "gmres" ) {
     m_LinearSolver = new GMRESSolver<T,Ops>;
  }
  else if ( m_linear_solver == "BICGSTAB" || 
            m_linear_solver == "bicgstab" ||
            m_linear_solver == "BiCGStab" ) {
     m_LinearSolver = new BiCGStabSolver<T,Ops>;
  }
  else {
     MayDay::Error("NewtonSolver<T,Ops>::define(): Unrecognized linear solver");
  }
  
  m_F.define(a_state); 
  m_dY.define(a_state);
  m_LinearSolver->define(m_Jacobian,a_outPrefix,a_optPrefix,a_state);

  m_is_Defined = true;
}

template <class T, class Ops>
void NewtonSolver<T,Ops>::solve(T& a_Y, T& a_b)
{
  CH_assert(isDefined());
  CH_assert(m_Function);

  m_its = 0;
  if (m_isLinear) m_maxits = 0; /* Take only 1 Newton iteration */

  while(1) {

    m_Function->evalFunction(m_F, a_Y, a_b, m_its); /* F = F(Y) - b */

    /* compute norm of m_F */
    m_norm_prev = m_norm;
    m_norm = m_F.computeNorm();
    if (m_its == 0) m_norm0 = m_norm;
    if (m_norm0 == 0) m_norm0 = 1e-16;

    if (m_verbose && (!procID())) {
      printf("    --> (%s.Newton) iteration %3d, residual = %1.6e (abs.), %1.6e (rel.)\n",
              m_outPrefix.c_str(), m_its, m_norm, m_norm/m_norm0 );
    }

    /* check for exit conditions */
    if (m_norm < m_atol) {
      m_convergedReason = "norm less than absolute tolerance";
      m_exitStatus = 0;
      if (m_verbose && (!procID())) {
        printf("    (%s.Newton) exiting: satisfied absolute tolerance (%1.3e).\n", 
               m_outPrefix.c_str(), m_atol);
      }
      break;
    }
    if (m_norm/m_norm0 < m_rtol) {
      m_convergedReason = "norm less than relative tolerance";
      m_exitStatus = 0;
      if (m_verbose && (!procID())) {
        printf("    (%s.Newton) exiting: satisfied relative tolerance (%1.3e).\n", 
               m_outPrefix.c_str(), m_rtol);
      }
      break;
    }
    if (m_its >= m_maxits) {
      m_convergedReason = "number of iterations exceeded maximum";
      m_exitStatus = 1;
      if (m_verbose && (!procID())) {
        printf("    (%s.Newton) exiting: iterations exceed max iterations (%d).\n", 
               m_outPrefix.c_str(), m_maxits);
      }
      break;
    }
    if (m_norm > 100*m_norm0) {
      m_convergedReason = "Newton solve diverged";
      m_exitStatus = -1;
      if (m_verbose && (!procID())) {
        printf( "    (%s.Newton) exiting: newton solver DIVERGED!\n",
                m_outPrefix.c_str() );
      }
      break;
    }

    /* Compute step size [Jac]dY = F */
    m_dY.zero();
    m_LinearSolver->solve(m_dY,m_F);
    Real step_norm = m_Function->computeStepNorm(m_dY);

    /* Update solution Y = Y - dY */
    a_Y.increment(m_dY,-1);
    m_its++;

    if (step_norm < m_stol) {
      m_convergedReason = "step size norm less than tolerance";
      if (m_verbose && (!procID())) {
        printf( "    (%s.Newton) exiting: step norm (%1.3e) satisfied step size tolerance (%1.3e).\n",
                m_outPrefix.c_str(), step_norm, m_stol );
      }
      break;
    }

  }
  m_count += m_its;
}

template <class T, class Ops>
void NewtonSolver<T,Ops>::parseParameters( ParmParse& a_pp)
{
  a_pp.query("verbose",       m_verbose);        /* verbosity */
  a_pp.query("atol",          m_atol);           /* absolute tolerance */
  a_pp.query("rtol",          m_rtol);           /* relative tolerance */
  a_pp.query("stol",          m_stol);           /* step size tolerance */
  a_pp.query("maxits",        m_maxits);         /* maximum iterations */
  a_pp.query("linear_solver", m_linear_solver);  /* linear solver */
}

#include "NamespaceFooter.H"

#endif
