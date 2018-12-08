#ifndef _ImplicitStageFunction_H_
#define _ImplicitStageFunction_H_

#include "Function.H"
#include "ImplicitStageJacobian.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class ImplicitStageFunction: public Function<T,Ops>
{
  public:
    ImplicitStageFunction<T,Ops>() {m_is_Defined = false; m_shift = 0; }
    ~ImplicitStageFunction<T,Ops>() {}

   void evalFunction(T& a_F, const T& a_Y, const T& a_b, int a_its, bool a_update_fluid_ops);
    void define(const T& a_state, Ops& a_ops);

    inline void setShift        (Real a_shift )       { m_shift = a_shift;    }
    inline void setStageTime    (Real a_time  )       { m_StageTime = a_time; }
    inline bool isDefined       ()              const { return m_is_Defined;  }

    inline void setJac(ImplicitStageJacobian<T,Ops>* a_jac) { m_jac = a_jac; }

  private:
    bool m_is_Defined;
    Real m_shift, m_StageTime;
    T    m_Z, m_W;
    Ops  *m_ops;
    ImplicitStageJacobian<T,Ops> *m_jac;
};

template <class T, class Ops>
void ImplicitStageFunction<T,Ops>::define(const T& a_state, Ops& a_ops)
{
  m_Z.define(a_state);
  m_W.define(a_state);
  m_ops = &a_ops;
  m_is_Defined = true;
}

/*
 * To evaluate: alpha*M[Y(i)] - F[Y(i)] - R
 * where
 *  alpha = 1 / (dt*a_ii) (m_shift)
 *  Y(i) is the stage solution to solve for
 *  M(Y) is the mass matrix function
 *  F(Y) is the implicitly-treated function
 *  R    is the right-hand-side
*/
template <class T, class Ops>
void ImplicitStageFunction<T,Ops>::evalFunction(T& a_F, const T& a_Y, const T& a_b, int a_its, bool a_update_fluid_ops)
{
  CH_assert(isDefined());

  if ( a_its==0 ) m_ops->preOpEval(a_Y,m_StageTime,a_update_fluid_ops);
  if (m_ops->trivialSolutionOp()){
    m_W.copy(a_Y);
  } else{
    m_ops->solutionOp(m_W,m_StageTime,a_Y);
  }
  m_ops->implicitOpImEx(m_Z,m_StageTime,a_Y);

  /* set base Y and G(Y) for Jacobian */
  m_jac->setBaseSolution(a_Y);
  m_jac->setBaseRHS     (m_Z, m_W);

  /* Compute ([(1/shift)*M(Y) - G(Y)] - b) */
  a_F.copy(m_W);
  a_F.scale(m_shift);
  a_F.increment(m_Z,-1);
  a_F.increment(a_b,-1);
}

#include "NamespaceFooter.H"
#endif
