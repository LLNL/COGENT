#ifndef _NewtonSolver_H_
#define _NewtonSolver_H_

#include <iostream>
#include <string>

#include "REAL.H"
#include "ParmParse.H"
#include "parstream.H"
#include "Function.H"
#include "LinearFunction.H"
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
    NewtonSolver<T,Ops>() { m_is_Defined = false; m_Function = NULL; m_Jacobian = NULL; }

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

  protected:

  private:

    void parseParameters( ParmParse& );

    bool      m_is_Defined;
    int       m_its, m_maxits, m_exitStatus;
    Real      m_norm, m_norm0,
              m_rtol, m_atol, m_stol;
    T         m_dY, m_F;
    bool      m_verbose;

    std::string           m_convergedReason,
                          m_outPrefix,
                          m_optPrefix;

    Function<T,Ops>       *m_Function;
    LinearFunction<T,Ops> *m_Jacobian;
    GMRESSolver<T,Ops>    m_LinearSolver;
};

template <class T, class Ops>
NewtonSolver<T,Ops>::~NewtonSolver()
{
  /*
   * The calling module using this solver is responsible for
   * deallocating the objects pointed to by m_Function and
   * m_Jacobian (since it will also create these objects).
   */
  return;
}

template <class T, class Ops>
void NewtonSolver<T,Ops>::define(ParmParse& a_pp, T& a_state, std::string a_outPrefix, std::string a_optPrefix)
{
  CH_assert(!isDefined());

  /* default values */
  m_maxits  = 100;
  m_rtol    = 1e-6;
  m_atol    = 1e-6;
  m_stol    = 1e-14;
  m_verbose = true;

  m_outPrefix = a_outPrefix;
  m_optPrefix = a_optPrefix;
  std::string optString = "_newton";
  optString = m_optPrefix + optString;

  ParmParse ppNewton( optString.c_str() );
  parseParameters( ppNewton );

  m_F.define(a_state); 
  m_dY.define(a_state);
  m_LinearSolver.define(a_pp,m_Jacobian,a_outPrefix,a_optPrefix);

  m_is_Defined = true;
}

template <class T, class Ops>
void NewtonSolver<T,Ops>::solve(T& a_Y, T& a_b)
{
  CH_assert(isDefined());
  CH_assert(m_Function);

  m_its = 0;
  while(1) {
    m_Function->evalFunction(m_F,a_Y,a_b); /* F = F(Y) - b */
    m_its++;

    /* compute norm of m_F */
    m_norm = m_F.computeNorm(2);
    if (m_its == 1) m_norm0 = m_norm;

    if (m_verbose && (!procID())) {
      cout << "    --> (" << m_outPrefix << ") Newton iteration " << m_its << ", absolute residual norm " << m_norm;
      cout << ", relative residual norm " << m_norm/m_norm0 << "\n"; 
    }

    /* check for exit conditions */
    if (m_norm < m_atol) {
      m_convergedReason = "norm less than absolute tolerance";
      m_exitStatus = 0;
      break;
    }
    if (m_norm/m_norm0 < m_rtol) {
      m_convergedReason = "norm less than relative tolerance";
      m_exitStatus = 0;
      break;
    }
    if (m_its > m_maxits) {
      m_convergedReason = "number of iterations exceeded maximum";
      m_exitStatus = 1;
      break;
    }
    if (m_norm > m_norm0) {
      m_convergedReason = "Newton solve diverged";
      m_exitStatus = -1;
      break;
    }

    /* Compute step size [Jac]dY = F */
    m_dY.zero();
    m_LinearSolver.solve(m_dY,m_F);
    Real step_norm = m_dY.computeNorm(2);

    if (step_norm < m_stol) {
      m_convergedReason = "step size norm less than tolerance";
      break;
    }

    /* Update solution Y = Y - dY */
    a_Y.increment(m_dY,-1);
  }
}

template <class T, class Ops>
void NewtonSolver<T,Ops>::parseParameters( ParmParse& a_pp)
{
  a_pp.query("verbose", m_verbose); /* verbosity */
  a_pp.query("atol",    m_atol);    /* absolute tolerance */
  a_pp.query("rtol",    m_rtol);    /* relative tolerance */
  a_pp.query("stol",    m_stol);    /* step size tolerance */
  a_pp.query("maxits",  m_maxits);  /* maximum iterations */
}

#include "NamespaceFooter.H"

#endif
