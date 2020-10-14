#ifndef _ExplicitStageFunction_H_
#define _ExplicitStageFunction_H_

#include "Function.H"
#include "ExplicitStageJacobian.H"
#include "Checkpoints.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class ExplicitStageFunction: public Function<T,Ops>
{
  public:
    ExplicitStageFunction<T,Ops>() {m_is_Defined = false; m_shift = 0; }
    ~ExplicitStageFunction<T,Ops>() {}

    void evalFunction(T& a_F, const T& a_Y, const T& a_b, int a_its);
    void define(const T& a_state, Ops& a_ops);

    inline void setShift        (Real a_shift )       { m_shift = a_shift;    }
    inline void setStageTime    (Real a_time  )       { m_StageTime = a_time; }
    inline bool isDefined       ()              const { return m_is_Defined;  }

    inline void setJac(ExplicitStageJacobian<T,Ops>* a_jac) { m_jac = a_jac; }

    inline Real computeStepNorm(const GKVector& a_vec) const
    {
      return m_ops->scaledVecNormU(a_vec);
    }

  private:
    bool m_is_Defined;
    Real m_shift, m_StageTime;
    T    m_W;
    Ops  *m_ops;
    ExplicitStageJacobian<T,Ops> *m_jac;
};

template <class T, class Ops>
void ExplicitStageFunction<T,Ops>::define(const T& a_state, Ops& a_ops)
{
  m_W.define(a_state);
  m_ops = &a_ops;
  m_is_Defined = true;
}

/*
 * To evaluate: M[Y(i)] - R
 * where
 *  alpha = 1 / (dt*a_ii) (m_shift)
 *  Y(i) is the stage solution to solve for
 *  M(Y) is the mass matrix function
 *  R    is the right-hand-side
*/
template <class T, class Ops>
void ExplicitStageFunction<T,Ops>::evalFunction(  T&        a_F, 
                                                  const T&  a_Y, 
                                                  const T&  a_b, 
                                                  int       a_its )
{
  CH_TIME("ExplicitStageFunction::evalFunction()");
  CH_assert(isDefined());

  m_ops->preOpEval( a_Y, 
                    m_StageTime, 
                    (a_its == 0 ? chkpt_stage_func_0 : chkpt_stage_func_n) );

  m_ops->solutionOp(m_W,m_StageTime,a_Y);

  /* set base Y and G(Y) for Jacobian */
  m_jac->setBaseSolution(a_Y);
  m_jac->setBaseRHS     (m_W);

  /* Compute (M(Y) - b) */
  a_F.copy(m_W);
  a_F.increment(a_b,-1);

  if (a_its == 0) {
    m_ops->computeUScale(a_Y);
    m_ops->computeLHSOpScale(a_F);
  }

  m_ops->scaleLHSOp(a_F);

  return;
}

#include "NamespaceFooter.H"
#endif
