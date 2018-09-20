#ifndef _PointJacobiSolver_H_
#define _PointJacobiSolver_H_

#include "LinearFunction.H"
#include "LinearFunctionSolver.H"
#include "ParmParse.H"
#include "parstream.H"

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class PointJacobiSolver : public LinearFunctionSolver<T, Ops>
{
  public:

    PointJacobiSolver<T,Ops>()   { m_is_defined = false; m_op = NULL; }
    ~PointJacobiSolver<T,Ops>()  {}

    void solve  (T& a_x, const T& a_b, int a_its = -1);
    void define (LinearFunction<T,Ops>* const, const std::string&, const std::string&, const T&);

    inline bool isDefined() const { return m_is_defined; }

    inline void setTolerances(Real a_atol, Real a_rtol)
    {
      m_atol = a_atol;
      m_rtol = m_rtol;
    }

  private:
    
    bool m_is_defined;
    bool m_verbosity;
    std::string m_outPrefix, m_optPrefix;
    int m_its; /*!< Number of Jacobi iterations */
    Real m_omega;
    T  m_y;     /*!< Work vector */
    Real m_atol, m_rtol;

    LinearFunction<T,Ops>* m_op;

    void parseParameters(ParmParse&);
};

template <class T, class Ops>
void PointJacobiSolver<T,Ops>::define(LinearFunction<T,Ops> * const a_op, 
                                      const std::string&            a_out_prefix, 
                                      const std::string&            a_opt_prefix, 
                                      const T&                      a_x )
{
  /* default values */
  m_verbosity = true;
  m_omega = 0.9;
  m_its = 200;
  m_atol = 1e-10;
  m_rtol = 1e-10;

  m_outPrefix = a_out_prefix + ".Jacobi";
  m_optPrefix = a_opt_prefix + "_jacobi";

  m_op = a_op;
  CH_assert(m_op != NULL);

  ParmParse pp_Jac( m_optPrefix.c_str() );
  parseParameters( pp_Jac );

  m_y.define(a_x);
  m_is_defined = true;
}

template <class T, class Ops>
void PointJacobiSolver<T,Ops>::parseParameters( ParmParse& a_pp )
{
  a_pp.query("its", m_its);
  a_pp.query("omega", m_omega);
  a_pp.query("verbose", m_verbosity);
  a_pp.query("atol", m_atol);
  a_pp.query("rtol", m_rtol);
  return;
}

template <class T, class Ops>
void PointJacobiSolver<T,Ops>::solve(T& a_x, const T& a_b, int a_its)
{
  CH_assert(m_is_defined);
  if (a_its < 0) a_its = m_its;
  Real norm, norm0;

  std::vector<Real> diag;
  m_op->getDiagonal(diag);

  int n;
  for (n = 0; n < a_its; n++) {

    m_op->residual(m_y, a_x, a_b);
    norm = m_y.computeNorm(2);
    if (n == 0) norm0 = norm;
    if ((norm < m_atol) || (norm/norm0 < m_rtol)) break;

    if (m_verbosity && (!procID())) {
      cout  << "      --> (" << m_outPrefix << ") iteration " << n
            << ", residual = " << norm << endl;
    }

    //m_op->applyOp(m_y, a_x);
    m_y += a_b;
    for (int i(0); i < diag.size(); i++) {
      a_x[i] += m_omega * (a_b[i] - m_y[i]) / diag[i];
    }

  }

  if (m_verbosity && (!procID()) && (n > 0) && (n < a_its) ) {
    cout  << "      --> (" << m_outPrefix << ") iteration " << n
          << ", residual = " << norm << endl;
  }

  return;
}

#include "NamespaceFooter.H"
#endif
