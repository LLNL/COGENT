#ifndef _ImplicitStage_H_
#define _ImplicitStage_H_

#include "Function.H"

#include "NamespaceHeader.H"

template <class Solution, class RHS, class Ops>
class ImplicitStage: public Function<Solution, RHS, Ops>
{
  public:
    ImplicitStage<Solution,RHS,Ops>() {m_is_Defined = false; m_shift = 0; }
    ~ImplicitStage<Solution,RHS,Ops>() {}

    void evalFunction(RHS& a_F, Solution& a_Y, RHS& a_b);

    void setShift(Real a_shift) { m_shift = a_shift; }
    void setStageTime(Real a_time) { m_StageTime = a_time; }
    void setStageNumber(int a_Stage) { m_Stage = a_Stage; }

    bool isDefined() const{ return m_is_Defined; }
    void define(Solution& a_state, Ops& a_ops);

  private:
    bool m_is_Defined;
    int  m_Stage;
    Real m_shift, m_StageTime;
    RHS  m_Z;
    Ops  *m_ops;
};

template <class Solution, class RHS, class Ops>
void ImplicitStage<Solution,RHS,Ops>::define(Solution& a_state, Ops& a_ops)
{
  m_Z.define(a_state);
  m_ops = &a_ops;
  m_is_Defined = true;
}

template <class Solution, class RHS, class Ops>
void ImplicitStage<Solution,RHS,Ops>::evalFunction(RHS& a_F, Solution& a_Y, RHS& a_b)
{
  CH_assert(isDefined());
  /* Compute ([Y - shift*G(Y)] - b) */
  a_F.copy(a_Y);
  m_ops->implicitOpImEx(m_Z,m_StageTime,a_Y,m_Stage);
  a_F.increment(m_Z,-m_shift);
  a_F.increment(a_b,-1);
}

#include "NamespaceFooter.H"
#endif
