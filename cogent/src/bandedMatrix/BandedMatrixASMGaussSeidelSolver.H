#ifndef _BandedMatrixASMGaussSeidelSolver_H_
#define _BandedMatrixASMGaussSeidelSolver_H_

#include "BandedMatrix.H"
#include "BandedMatrixSolver.H"

#include "NamespaceHeader.H"

template <class Vec>
class BandedMatrixASMGaussSeidelSolver : public BandedMatrixSolver<Vec>
{
  public:

    BandedMatrixASMGaussSeidelSolver<Vec>()   { m_is_defined = false; }
    ~BandedMatrixASMGaussSeidelSolver<Vec>()  {}

    void solve  (Vec&, const BandedMatrix&, const Vec&, int a_its = -1);
    void define (const std::string&, const std::string&, const Vec&);

    bool isDefined() const { return m_is_defined; }

  private:
    
    bool m_is_defined;
    bool m_verbosity;
    int  m_its; /*!< Number of GaussSeidel iterations */

    std::string m_opt_string, m_out_string;

    void parseParameters(ParmParse&);
};

template <class Vec>
void BandedMatrixASMGaussSeidelSolver<Vec>::define( const std::string& a_opt_prefix, 
                                                    const std::string& a_out_prefix,
                                                    const Vec&)
{
  /* default values */
  m_its = 5;

  m_opt_string = a_opt_prefix + "_asm_gauss_seidel";
  m_out_string = a_out_prefix + ".ASMGaussSeidel";
  ParmParse pp_ASMGS( m_opt_string.c_str() );
  parseParameters( pp_ASMGS );

  m_is_defined = true;
}

template <class Vec>
void BandedMatrixASMGaussSeidelSolver<Vec>::parseParameters( ParmParse& a_pp )
{
  a_pp.query("its", m_its);
  a_pp.query("verbose", m_verbosity);
  return;
}

template <class Vec>
void BandedMatrixASMGaussSeidelSolver<Vec>::solve(Vec& a_x, const BandedMatrix& a_A, const Vec& a_b, int a_its)
{
  CH_assert(m_is_defined);
  CH_assert(a_A.isAssembled());

  if (a_its < 0) a_its = m_its;

  /* GaussSeidel iterations with the local part of the matrix */
  for (int n(0); n < a_its; n++) {
    for (int i(0); i < a_A.nRows(); i++) {
      Real *data_ptr      = a_A.data()        + i*a_A.getNBands();
      int  *col_ptr       = a_A.iCols()       + i*a_A.getNBands();
      int  *local_col_ptr = a_A.iColsLocal()  + i*a_A.getNBands();

      Real diag;
      Real rhs = a_b[i];

      for (int k(0); k < a_A.nColsLocal(i); k++) {

        int icol = local_col_ptr[k];
        int j    = a_A.globalToLocal(col_ptr[icol]);

        CH_assert(a_A.isLocal(col_ptr[icol]));

        if (j == i) {
          diag = data_ptr[icol];
        } else {
          rhs -= (data_ptr[icol]*a_x[j]);
        }
      }

      a_x[i] = rhs / diag;
    }
  }
  return;
}

#include "NamespaceFooter.H"
#endif
