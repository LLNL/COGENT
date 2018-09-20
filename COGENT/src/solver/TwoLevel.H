#ifndef _TwoLevelSolver_H_
#define _TwoLevelSolver_H_

#include "LinearFunction.H"
#include "LinearFunctionSolver.H"
#include "ParmParse.H"
#include "parstream.H"
#include "BandedMatrixSolver.H"
#include "BandedMatrixGaussSeidelSolver.H"
#include "BandedMatrixJacobiSolver.H"

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class TwoLevelSolver : public LinearFunctionSolver<T, Ops>
{
  public:

    TwoLevelSolver<T,Ops>() 
    {
      m_is_defined = false; 
      m_op = NULL; 
      m_smoother = NULL;
      m_coarse_solver = NULL;
    }
    ~TwoLevelSolver<T,Ops>()  
    {
      delete m_smoother;
      delete m_coarse_solver;
    }

    void solve  (T& a_x, const T& a_b, int a_its = -1);
    void define (LinearFunction<T,Ops>* const, const std::string&, const std::string&, const T&);

    inline bool isDefined() const { return m_is_defined; }

    inline void setTolerances(Real a_atol, Real a_rtol)
    {
      m_atol = a_atol;
      m_rtol = m_rtol;
    }

  private:
    
    bool  m_is_defined;
    bool  m_verbosity;
    bool  m_norm_based_exit;

    int   m_its; /*!< Number of two level iterations */
    int   m_its_presmoothing;
    int   m_its_coarse;
    int   m_its_postsmoothing;

    Real  m_atol, m_rtol;

    std::string m_outPrefix, m_optPrefix;
    std::string m_smoother_type;
    std::string m_coarse_solver_type;

    T m_y;   /*!< Work vector */
    T m_r;   /*!< Work vector */

    LinearFunction<T,Ops>* m_op;
    LinearFunctionSolver<T,Ops>* m_smoother;
    BandedMatrixSolver<T>* m_coarse_solver;

    void parseParameters(ParmParse&);
};

template <class T, class Ops>
void TwoLevelSolver<T,Ops>::define(LinearFunction<T,Ops> * const  a_op, 
                                   const std::string&             a_out_prefix, 
                                   const std::string&             a_opt_prefix, 
                                   const T&                       a_x )
{
  /* default values */
  m_its = 10;
  m_its_presmoothing = 10;
  m_its_coarse = 100;
  m_its_postsmoothing = 10;
  m_verbosity = true;
  m_smoother_type = "jacobi";
  m_coarse_solver_type = "gauss_seidel";
  m_atol = 1e-10;
  m_rtol = 1e-10;
  m_norm_based_exit = false;

  m_outPrefix = a_out_prefix + ".TwoLevel";
  m_optPrefix = a_opt_prefix + "_twolevel";

  m_op = a_op;
  CH_assert(m_op != NULL);

  ParmParse pp_TL( m_optPrefix.c_str() );
  parseParameters( pp_TL );

  m_y.define(a_x);
  m_r.define(a_x);

  if (m_smoother_type == "jacobi") {
    m_smoother = new PointJacobiSolver<T,Ops>;
  } else {
    if (!procID()) {
      cout << "Error in TwoLevelSolver::define() - " << m_smoother_type
           << " is not a valid type.\n";
    }
    m_smoother = NULL;
  }
  CH_assert(m_smoother != NULL);
  std::string opt_string_smoother = m_optPrefix + "_smoother";
  std::string out_string_smoother = m_outPrefix + ".Smoother";
  m_smoother->define(m_op, out_string_smoother, opt_string_smoother, a_x);

  if (m_coarse_solver_type == "gauss_seidel") {
    m_coarse_solver = new BandedMatrixGaussSeidelSolver<T>;
  } else if (m_coarse_solver_type == "jacobi") {
    m_coarse_solver = new BandedMatrixJacobiSolver<T>;
  } else {
    if (!procID()) {
      cout << "Error in TwoLevelSolver::define() - " << m_coarse_solver_type
           << " is not a valid type.\n";
    }
    m_coarse_solver = NULL;
  }
  CH_assert(m_coarse_solver != NULL);
  std::string opt_string_coarse = m_optPrefix + "_coarse_solver";
  std::string out_string_coarse = m_outPrefix + ".CoarseSolver";
  m_coarse_solver->define(opt_string_coarse, out_string_coarse, a_x);

  m_is_defined = true;
}

template <class T, class Ops>
void TwoLevelSolver<T,Ops>::parseParameters( ParmParse& a_pp )
{
  a_pp.query("its", m_its);
  a_pp.query("presmoothing_its", m_its_presmoothing);
  a_pp.query("coarse_solve_its", m_its_coarse);
  a_pp.query("postsmoothing_its",m_its_postsmoothing);
  a_pp.query("coarse_solver_type",m_coarse_solver_type);
  a_pp.query("smoother_type",m_smoother_type);
  a_pp.query("verbose", m_verbosity);
  a_pp.query("atol", m_atol);
  a_pp.query("rtol", m_rtol);
  a_pp.query("norm_based_exit", m_norm_based_exit);
  return;
}

template <class T, class Ops>
void TwoLevelSolver<T,Ops>::solve(T& a_x, const T& a_b, int a_its)
{
  CH_assert(m_is_defined);
  if (a_its < 0) a_its = m_its;
  Real norm, norm0;

  int n;
  for (n = 0; n < a_its; n++) {

    if (m_norm_based_exit) {
      m_op->residual(m_r, a_x, a_b);
      norm = m_r.computeNorm(2);
      if (n == 0) norm0 = norm;
      if ((norm < m_atol) || (norm/norm0 < m_rtol)) break;
    }

    if (m_verbosity && (!procID())) {
      cout  << "      --> (" << m_outPrefix << ") iteration " << n
            << ", residual = " << norm << endl;
    }

    /* pre-smoothing */
    m_smoother->solve(a_x, a_b, m_its_presmoothing);

    /* coarse solve */
    m_op->residual(m_r, a_x, a_b);
    m_r *= -1.0;
    m_y.zero();
    m_coarse_solver->solve(m_y,m_op->getBandedMatrixForm(),m_r,m_its_coarse);
    a_x += m_y;

    /* post-smoothing */
    m_smoother->solve(a_x, a_b, m_its_postsmoothing);

  }

  if ((m_verbosity) && (!procID()) && (n > 0) && (n < a_its)) {
    cout  << "      --> (" << m_outPrefix << ") iteration " << n
          << ", residual = " << norm << endl;
  }

  return;
}

#include "NamespaceFooter.H"
#endif
