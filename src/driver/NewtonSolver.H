#ifndef _NewtonSolver_H_
#define _NewtonSolver_H_

#include <iostream>
#include <string>

#include "REAL.H"
#include "Function.H"

#include "NamespaceHeader.H"

template <class Solution, class RHS, class Ops>
class NewtonSolver
{
  public:

    /// Constructor
    /**
     * Constructor: set m_is_Defined to false.
     */ 
    NewtonSolver<Solution,RHS,Ops>() { m_is_Defined = false; }

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
    void define(Solution& a_state);

    /// Get the number of iterations
    /**
     * get the number of Newton iterations
     */
    int getNewtonIterations() const { return m_its; }

    /// Solve the equation using Newton's method
    /**
     * solve F(x) = b using Newton's method
     *
     * @param[in] a_Y initial guess; contains the solution at exit
     * @param[in] a_b right-hand side
     */
    void solve(Solution& a_Y, RHS& a_b, Function<Solution,RHS,Ops>& a_Func);

    /// Get converged reason
    /**
     * Get the reason for Newton solver exit
     */
    std::string getConvergedReason() const {return m_convergedReason; }

  protected:

  private:

    bool      m_is_Defined;
    int       m_its, m_maxits;
    Real      m_norm, m_norm0,
              m_rtol, m_atol;
    Solution  m_dY;
    RHS       m_F;
    bool      m_verbose;

    std::string m_convergedReason;
};

template <class Solution, class RHS, class Ops>
NewtonSolver<Solution,RHS,Ops>::~NewtonSolver()
{
  return;
}

template <class Solution, class RHS, class Ops>
void NewtonSolver<Solution,RHS,Ops>::define(Solution& a_state)
{
  CH_assert(!isDefined());

  m_F.define(a_state); 
  m_dY.define(a_state);

  /* default values */
  m_maxits  = 10000;
  m_rtol    = 1e-6;
  m_atol    = 1e-6;
  m_verbose = true;

  m_is_Defined = true;
}

template <class Solution, class RHS, class Ops>
void NewtonSolver<Solution,RHS,Ops>::solve(Solution& a_Y, RHS& a_b, Function<Solution,RHS,Ops>& a_Func)
{
  CH_assert(isDefined());

  m_its = 0;
  while(1) {
    a_Func.evalFunction(m_F,a_Y,a_b); /* F = F(Y) - b */
    m_its++;

    /* compute norm of m_F */
    Real norm_1,norm_2,norm_inf;
    m_F.computeNorms(norm_1,norm_2,norm_inf);
    m_norm = norm_2; /* using the 2-norm */
    if (m_its == 0) m_norm0 = m_norm;

    if (m_verbose && (!procID())) {
      cout << "    Newton iteration " << m_its << ", absolute residual norm " << m_norm;
      cout << ", relative residual norm " << m_norm/m_norm0 << "\n"; 
    }

    /* check for exit conditions */
    if (m_norm < m_atol) {
      m_convergedReason = "norm less than absolute tolerance";
      break;
    }
    if (m_norm/m_norm0 < m_rtol) {
      m_convergedReason = "norm less than relative tolerance";
      break;
    }
    if (m_its > m_maxits) {
      m_convergedReason = "number of iterations exceeded maximum";
      break;
    }

    /* Compute step size */

    /* Update solution */
    a_Y.increment(m_dY);
  }

}

#include "NamespaceFooter.H"

#endif
