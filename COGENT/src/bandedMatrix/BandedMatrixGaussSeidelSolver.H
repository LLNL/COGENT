#ifndef _BandedMatrixGaussSeidelSolver_H_
#define _BandedMatrixGaussSeidelSolver_H_

#include "BandedMatrix.H"
#include "BandedMatrixSolver.H"

#include "NamespaceHeader.H"

template <class Vec>
class BandedMatrixGaussSeidelSolver : public BandedMatrixSolver<Vec>
{
  public:

    BandedMatrixGaussSeidelSolver<Vec>()   { m_is_defined = false; }
    ~BandedMatrixGaussSeidelSolver<Vec>()  {}

    void solve  (Vec&, const BandedMatrix&, const Vec&, int a_its = -1);
    void define (const std::string&, const std::string&, const Vec&);

    bool isDefined() const { return m_is_defined; }

  private:
    
    bool m_is_defined;
    bool m_verbosity;
    int  m_its; /*!< Number of GaussSeidel iterations */
    std::string m_opt_string, m_out_string;
    std::unordered_map<int,Real> m_nonlocal_data;

    void parseParameters(ParmParse&);
};

template <class Vec>
void BandedMatrixGaussSeidelSolver<Vec>::define(const std::string& a_opt_prefix, 
                                                const std::string& a_out_prefix,
                                                const Vec&)
{
  /* default values */
  m_its = 100;
  m_verbosity = false;

  m_opt_string = a_opt_prefix + "_gauss_seidel";
  m_out_string = a_out_prefix + ".GaussSeidel";
  ParmParse pp_GS( m_opt_string.c_str() );
  parseParameters( pp_GS );

  m_is_defined = true;
}

template <class Vec>
void BandedMatrixGaussSeidelSolver<Vec>::parseParameters( ParmParse& a_pp )
{
  a_pp.query("its", m_its);
  a_pp.query("verbose", m_verbosity);
  return;
}

template <class Vec>
void BandedMatrixGaussSeidelSolver<Vec>::solve( Vec&                a_x, 
                                                const BandedMatrix& a_A, 
                                                const Vec&          a_b, 
                                                int                 a_its)
{
  CH_assert(m_is_defined);
  CH_assert(a_A.isAssembled());

  if (a_its < 0) a_its = m_its;

  for (int n(0); n < a_its; n++) {

    BandedMatrixSolver<Vec>::exchangeData(a_A, a_x, m_nonlocal_data);

    for (int i(0); i < a_A.nRows(); i++) {

      Real *data_ptr      = a_A.data()  + i*a_A.getNBands();
      int  *col_ptr       = a_A.iCols() + i*a_A.getNBands();

      Real diag;
      Real rhs = a_b[i];

      for (int k(0); k < a_A.nCols(i); k++) {

        int icol = col_ptr[k];

        if (a_A.isLocal(icol)) {
          
          int j = a_A.globalToLocal(icol);
          if (j == i) {
            diag = data_ptr[k];
          } else {
            rhs -= ( data_ptr[k] * a_x[j] );
          }

        } else {

          rhs -= (data_ptr[k] * m_nonlocal_data.at(icol) );

        }
      }

      a_x[i] = rhs / diag;
    }

    if (m_verbosity) {
      Vec tmp;
      tmp.define(a_x);
      BandedMatrixSolver<Vec>::applyOp(tmp, a_x, a_A);
      tmp *= -1.0;
      tmp += a_b;
      Real norm = tmp.computeNorm();
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
