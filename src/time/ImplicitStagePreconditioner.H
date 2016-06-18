#ifndef _ImplicitStagePreconditioner_H_
#define _ImplicitStagePreconditioner_H_

#include "BandedMatrix.H"

#include "NamespaceHeader.H"

template <class T,class Ops>
class ImplicitStagePreconditioner
{
  public:
    ImplicitStagePreconditioner<T,Ops>() { m_is_Defined = false; m_shift = 0; }
    ~ImplicitStagePreconditioner<T,Ops>() {}

    void define   (T&,Ops&);
    void applyPinv(T&, const T&);
    void updateP  (const T&);

    inline void setShift  (Real a_a)  { m_shift = a_a; }

  private:
    bool  m_is_Defined;
    Real  m_shift;
    Ops   *m_ops;

    BandedMatrix m_P;
};

template <class T,class Ops>
void ImplicitStagePreconditioner<T,Ops>::define(T& a_state,Ops& a_ops)
{
  m_ops = &a_ops;
  m_is_Defined = m_ops->setupPCImEx((void*)&m_P,a_state);
}

template <class T,class Ops>
void ImplicitStagePreconditioner<T,Ops>::updateP(const T& a_x)
{
  if (!m_is_Defined) {
    m_P.zeroEntries();
    m_ops->assemblePCImEx((void*)&m_P,a_x);
    m_P.scaleEntries(-1.0);
    m_P.shift(m_shift);
  }
  return;
}

template <class T,class Ops>
void ImplicitStagePreconditioner<T,Ops>::applyPinv(T& a_y,const T& a_x)
{
  //if (!m_is_Defined) a_y.copy(a_x);
  //else {
  //}
  a_y.copy(a_x);
  return;
}

#include "NamespaceFooter.H"
#endif
