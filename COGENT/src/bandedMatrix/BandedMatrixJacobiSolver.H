#ifndef _BandedMatrixJacobiSolver_H_
#define _BandedMatrixJacobiSolver_H_

#include <iostream>
#include "BandedMatrix.H"
#include "BandedMatrixSolver.H"

#include "NamespaceHeader.H"

template <class Vec>
class BandedMatrixJacobiSolver : public BandedMatrixSolver<Vec>
{
  public:

    BandedMatrixJacobiSolver<Vec>()   { m_is_defined = false; }
    ~BandedMatrixJacobiSolver<Vec>()  {}

    void solve  (Vec&, const BandedMatrix&, const Vec&, int a_its = -1);
    void define (const std::string&, const std::string&, const Vec&);

    bool isDefined() const { return m_is_defined; }

  private:
    
    bool m_is_defined;
    bool m_verbosity;
    int  m_its; /*!< Number of Jacobi iterations */
    Vec  m_y;   /*!< Work vector */
    std::string m_opt_string, m_out_string;

    std::unordered_map<int,Real> m_nonlocal_data;

    void parseParameters(ParmParse&);
};

template <class Vec>
void BandedMatrixJacobiSolver<Vec>::define( const std::string&  a_opt_prefix, 
                                            const std::string&  a_out_prefix,
                                            const Vec&          a_x)
{
  /* default values */
  m_its = 10;
  m_verbosity = true;

  m_opt_string = a_opt_prefix + "_jacobi";
  m_out_string = a_out_prefix + ".Jacobi";
  ParmParse pp_Jac( m_opt_string.c_str() );
  parseParameters( pp_Jac );

  m_y.define(a_x);
  m_is_defined = true;
}

template <class Vec>
void BandedMatrixJacobiSolver<Vec>::parseParameters( ParmParse& a_pp )
{
  a_pp.query("its", m_its);
  a_pp.query("verbose", m_verbosity);
  return;
}

template <class Vec>
void BandedMatrixJacobiSolver<Vec>::solve(Vec& a_x, const BandedMatrix& a_A, const Vec& a_b, int a_its)
{
  CH_assert(m_is_defined);
  CH_assert(a_A.isAssembled());

  if (a_its < 0) a_its = m_its;

  std::vector<Real> diag;
  a_A.getDiagonal(diag);
  CH_assert(diag.size() == a_A.nRows());

  for (int n(0); n < a_its; n++) {

    BandedMatrixSolver<Vec>::applyOp(m_y, a_x, a_A);
    for (int i(0); i < a_A.nRows(); i++) {
      m_y[i] -= diag[i] * a_x[i];
      a_x[i] = (a_b[i] - m_y[i]) / diag[i];
    }

    if (m_verbosity) {
      BandedMatrixSolver<Vec>::applyOp(m_y, a_x, a_A);
      m_y *= -1.0;
      m_y += a_b;
      Real norm = m_y.computeNorm(2);
      if (!procID()) {
        std::cout << "      --> (" << m_out_string << ") iteration " << n+1
                  << ", residual = " << norm << std::endl;
      }
    }

  }

  return;
}

#include "NamespaceFooter.H"
#endif
