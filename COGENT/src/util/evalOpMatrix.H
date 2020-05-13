#ifndef _EvalOpMatrix_H_
#define _EvalOpMatrix_H_

#include <stdio.h>
#include <iostream>
#include <string>
#include <time.h>
#include "BandedMatrix.H"

#include "NamespaceHeader.H"

template<class Vector, class Ops>
class EvalFunction
{
  public: 
    virtual ~EvalFunction<Vector,Ops>() {}
    virtual void define(Ops*) = 0;
    virtual void eval(Vector&, const Vector&, Real) const = 0;
};

template<class Vector, class Ops>
class EvalExplicitFunction : public EvalFunction<Vector,Ops>
{
  public:

    EvalExplicitFunction<Vector,Ops>()
    {
      m_ops = NULL;
      m_is_defined = false;
    }

    ~EvalExplicitFunction<Vector,Ops>() {}

    inline
    void define(Ops *a_ops)
    {
      m_ops = a_ops;
      m_is_defined = true;
      return;
    }

    inline
    void eval(Vector& a_f, const Vector& a_u, Real a_time) const
    {
      CH_assert(m_is_defined);
      m_ops->explicitOp(a_f,a_time,a_u);
      return;
    }

  private:

    bool  m_is_defined;
    Ops*  m_ops;
};

template<class Vector, class Ops>
class EvalExplicitFunctionImEx : public EvalFunction<Vector,Ops>
{
  public:

    EvalExplicitFunctionImEx<Vector,Ops>()
    {
      m_ops = NULL;
      m_is_defined = false;
    }

    ~EvalExplicitFunctionImEx<Vector,Ops>() {}

    inline
    void define(Ops *a_ops)
    {
      m_ops = a_ops;
      m_is_defined = true;
      return;
    }

    inline
    void eval(Vector& a_f, const Vector& a_u, Real a_time) const
    {
      CH_assert(m_is_defined);
      m_ops->explicitOpImEx(a_f,a_time,a_u);
      return;
    }

  private:

    bool  m_is_defined;
    Ops*  m_ops;
};

template<class Vector, class Ops>
class EvalImplicitFunctionImEx : public EvalFunction<Vector,Ops>
{
  public:

    EvalImplicitFunctionImEx<Vector,Ops>()
    {
      m_ops = NULL;
      m_is_defined = false;
    }

    ~EvalImplicitFunctionImEx<Vector,Ops>() {}

    inline
    void define(Ops *a_ops)
    {
      m_ops = a_ops;
      m_is_defined = true;
      return;
    }

    inline
    void eval(Vector& a_f, const Vector& a_u, Real a_time) const
    {
      CH_assert(m_is_defined);
      m_ops->implicitOpImEx(a_f,a_time,a_u);
      return;
    }

  private:

    bool  m_is_defined;
    Ops*  m_ops;
};

template<class Vector, class Ops>
class EvalOpMatrix
{
  public:
  
    EvalOpMatrix<Vector,Ops>()
    {
      m_is_defined = false;
      m_ops = NULL;
    }

    ~EvalOpMatrix<Vector,Ops>() {}

    inline
    void define(Ops *a_ops, const Vector& a_vec, Real a_tolerance = 1e-6)
    {
      m_ops = a_ops;

      m_Y.define(a_vec);
      m_F0.define(a_vec);
      m_F.define(a_vec);

      m_explicit.define(a_ops);
      m_explicit_imex.define(a_ops);
      m_implicit_imex.define(a_ops);

      m_fname_explicit      = "matOPExplicit.dat";
      m_fname_explicit_imex = "matOPExplicitImEx.dat";
      m_fname_implicit_imex = "matOPImplicitImEx.dat";
      m_fname_pc_imex       = "matPCImplicitImEx.dat";

      m_tolerance = a_tolerance;
      m_epsJFNK = 1e-6;

      m_is_defined = m_ops->preallocatePCMatrixImEx((void*)&m_precon_matrix_imex, a_vec); 
      return;
    }

    void computeMatrices(const Vector& a_X, Real a_t, int a_i)
    {
      int op_flag = 1;
      int ndof = a_X.getVectorSize();
      if (!procID()) {
        std::cout << "  Computing operator matrix for RHS (NDOF = " << ndof << ").\n";
      }
      computeOpMatrix(a_X, m_explicit, a_t, a_i, m_fname_explicit, op_flag);
      if (!procID()) {
        std::cout << "  Computing operator matrix for explicit part of RHS (IMEX) (NDOF = " 
                  << ndof << ").\n";
      }
      computeOpMatrix(a_X, m_explicit_imex, a_t, a_i, m_fname_explicit_imex, op_flag);
      if (!procID()) {
        std::cout << "  Computing operator matrix for implicit part of RHS (IMEX) (NDOF = " 
                  << ndof << ").\n";
      }
      computeOpMatrix(a_X, m_implicit_imex, a_t, a_i, m_fname_implicit_imex, op_flag);
      if (!procID()) {
        std::cout << "  Assembling preconditioning matrix for implicit part of RHS (IMEX) (NDOF = " 
                  << m_precon_matrix_imex.nRows() << ").\n";
      }
      m_ops->assemblePCImEx(true, (void*)&m_precon_matrix_imex, a_X, a_t, 0.);
      m_precon_matrix_imex.writeToFile(m_fname_pc_imex, op_flag);
    }

  private:

    bool    m_is_defined;
    Ops*    m_ops;
    Vector  m_Y, m_F0, m_F;
    Real    m_tolerance, m_epsJFNK;

    BandedMatrix  m_precon_matrix_imex;

    EvalExplicitFunction<Vector,Ops>      m_explicit;
    EvalExplicitFunctionImEx<Vector,Ops>  m_explicit_imex;
    EvalImplicitFunctionImEx<Vector,Ops>  m_implicit_imex;

    std::string m_fname_explicit, m_fname_explicit_imex, m_fname_implicit_imex,
                m_fname_pc_imex;

    inline Real abs(Real x) { return(x < 0 ? -x : x); }

    void computeOpMatrix(const Vector&, const EvalFunction<Vector,Ops>&,
                         Real, int, std::string&, int);

};

template <class Vector, class Ops>
void EvalOpMatrix<Vector,Ops>::computeOpMatrix( const Vector& a_Y0,
                                                const EvalFunction<Vector,Ops>& a_function,
                                                Real a_t, 
                                                int a_i, 
                                                std::string& a_fname,
                                                int a_op_flag )
{
  /* This operation is implemented only for 1 MPI rank */
  CH_assert(procID() == 0);
  /* It's too cumbersome to implement this for multiple MPI ranks.
   * Besides, do not **ever** use this O(N^2) function for large
   * problems */

  clock_t timer_start = clock();

  int ndof = a_Y0.getVectorSize();
  Real eps, norm = a_Y0.computeNorm(2);
  if (norm < m_tolerance) {
    eps = m_epsJFNK;
  } else {
    eps = m_epsJFNK * a_Y0.computeNorm(2);
  }

  if (!procID()) {
    std::cout << "  --> EvalOpMatrix<Vector,Ops>::computeOpMatrix: Writing to sparse matrix file "
              << a_fname << ".\n";
  }

  FILE *fout;
  if (a_op_flag == 0) {
    fout = fopen(a_fname.c_str(),"wb");
    fwrite(&ndof,sizeof(int),1,fout);
  } else {
    fout = fopen(a_fname.c_str(),"w");
    fprintf(fout,"%d\n",ndof);
  }

  m_F0.zero();
  a_function.eval(m_F0, a_Y0, a_t);

  for (int i=0; i<ndof; i++) {
    m_Y = a_Y0;
    m_Y[i] += eps;

    m_F.zero();
    a_function.eval(m_F, m_Y, a_t);

    m_F -= m_F0;

    for (int j=0; j<ndof; j++) {
      Real mat_elem = m_F[j] / eps;
      if (abs(mat_elem) > m_tolerance) {
        int ii = i+1, jj = j+1;
        if (a_op_flag == 0) {
          fwrite(&jj,sizeof(int),1,fout);
          fwrite(&ii,sizeof(int),1,fout);
          fwrite(&mat_elem,sizeof(double),1,fout);
        } else {
          fprintf(fout,"%10d  %10d  %+1.16e\n",j+1,i+1,mat_elem);
        }
      }
    }
  }

  clock_t timer_end = clock();

  double walltime = ((double)(timer_end - timer_start)) / CLOCKS_PER_SEC;
  if (!procID()) {
    std::cout << "  -->" << walltime << " seconds.\n";
  }

  fclose(fout);
  return;
}

#include "NamespaceFooter.H"

#endif
