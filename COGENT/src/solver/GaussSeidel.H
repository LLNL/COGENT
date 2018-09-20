#ifndef _GaussSeidelSolver_H_
#define _GaussSeidelSolver_H_

#include "LinearFunction.H"
#include "LinearFunctionSolver.H"
#include "ParmParse.H"
#include "parstream.H"
#include "BandedMatrixSolver.H"
#include "BandedMatrixGaussSeidelSolver.H"
#include "BandedMatrixJacobiSolver.H"

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class GaussSeidelSolver : public LinearFunctionSolver<T, Ops>
{
  public:

    GaussSeidelSolver<T,Ops>() 
    {
      m_is_defined = false; 
      m_op = NULL; 
      m_solver = NULL;
    }
    ~GaussSeidelSolver<T,Ops>()  
    {
      delete m_solver;
    }

    void solve  (T& a_x, const T& a_b, int a_its = -1);
    void define (LinearFunction<T,Ops>* const, const std::string&, const std::string&, const T&);

    inline bool isDefined() const { return m_is_defined; }

    inline void setTolerances(Real a_atol, Real a_rtol)
    {
      //m_solver->setTolerances(a_atol, a_rtol);
    }

  private:
    
    bool  m_is_defined;

    LinearFunction<T,Ops>* m_op;
    BandedMatrixSolver<T>* m_solver;
};

template <class T, class Ops>
void GaussSeidelSolver<T,Ops>::define(LinearFunction<T,Ops> * const a_op, 
                                      const std::string&            a_out_prefix, 
                                      const std::string&            a_opt_prefix, 
                                      const T&                      a_x )
{
  /* default values */
  m_op = a_op;
  CH_assert(m_op != NULL);

  m_solver = new BandedMatrixGaussSeidelSolver<T>;
  m_solver->define(a_opt_prefix, a_out_prefix, a_x);

  m_is_defined = true;
}

template <class T, class Ops>
void GaussSeidelSolver<T,Ops>::solve(T& a_x, const T& a_b, int a_its)
{
  CH_assert(m_is_defined);
  m_solver->solve(a_x,m_op->getBandedMatrixForm(),a_b,a_its);
  return;
}

#include "NamespaceFooter.H"
#endif
