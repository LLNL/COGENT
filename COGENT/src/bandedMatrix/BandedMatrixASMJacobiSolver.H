#ifndef _BandedMatrixASMJacobiSolver_H_
#define _BandedMatrixASMJacobiSolver_H_

#include "BandedMatrix.H"
#include "BandedMatrixSolver.H"

#include "NamespaceHeader.H"

template <class Vec>
class BandedMatrixASMJacobiSolver : public BandedMatrixSolver<Vec>
{
  public:

    BandedMatrixASMJacobiSolver<Vec>()   { m_is_defined = false; }
    ~BandedMatrixASMJacobiSolver<Vec>()  {}

    void solve  (Vec&, const BandedMatrix&, const Vec&, int a_its = -1);
    void define (const std::string&, const std::string&,const Vec&);

    bool isDefined() const { return m_is_defined; }

  private:
    
    bool m_is_defined;
    bool m_verbosity;
    int  m_its; /*!< Number of Jacobi iterations */
    Vec  m_y;   /*!< Work vector */
    std::string m_opt_string, m_out_string;

    void parseParameters(ParmParse&);
};

template <class Vec>
void BandedMatrixASMJacobiSolver<Vec>::define(const std::string&  a_opt_prefix,
                                              const std::string&  a_out_prefix,
                                              const Vec&          a_x)
{
  /* default values */
  m_its = 5;

  m_opt_string = a_opt_prefix + "_asm_jacobi";
  m_out_string = a_out_prefix + ".ASMJacobi";
  ParmParse pp_ASMJac( m_opt_string.c_str() );
  parseParameters( pp_ASMJac );

  m_y.define(a_x);
  m_is_defined = true;
}

template <class Vec>
void BandedMatrixASMJacobiSolver<Vec>::parseParameters( ParmParse& a_pp )
{
  a_pp.query("its", m_its);
  a_pp.query("verbose", m_verbosity);
  return;
}

template <class Vec>
void BandedMatrixASMJacobiSolver<Vec>::solve(Vec& a_x, const BandedMatrix& a_A, const Vec& a_b, int a_its)
{
  CH_assert(m_is_defined);
  CH_assert(a_A.isAssembled());

  if (a_its < 0) a_its = m_its;

  /* Jacobi iterations with the local part of the matrix */
  for (int n(0); n < a_its; n++) {

    m_y.zero();

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

      m_y[i] = rhs / diag;
    }

    a_x.copy(m_y);
  }
  return;
}

#include "NamespaceFooter.H"
#endif
