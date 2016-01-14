#ifndef _ImplicitStageJacobian_H_
#define _ImplicitStageJacobian_H_

#include "LinearFunction.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class ImplicitStageJacobian: public LinearFunction<T,Ops>
{
  public:

    ImplicitStageJacobian<T,Ops>() {m_is_Defined = false; m_shift = 0; }
    ~ImplicitStageJacobian<T,Ops>() {}

    void residual   (T& a_residual, const T& a_y, const T& a_rhs, bool a_homogeneous);
    void preCond    (T& a_cor, const T& a_residual);
    void applyOp    (T& a_F, const T& a_y, bool a_homogeneous);
    void axby       (T& a_Z, const T& a_X, const T& a_Y, Real a_a, Real a_b);

    inline void create    (T& a_Z, const T& a_Y)              { a_Z.define(a_Y); }
    inline void assign    (T& a_Z, const T& a_Y)              { a_Z.copy(a_Y); }
    inline void incr      (T& a_Z, const T& a_Y, Real a_scale){ a_Z.increment(a_Y,a_scale); } 
    inline void scale     (T& a_Y, const Real& a_scale)       { a_Y.scale(a_scale); }
    inline void setToZero (T& a_Y)                            { a_Y.zero(); }
    inline Real dotProduct(T& a_X, T& a_Y)                    { return(a_X.dotProduct(a_Y)); }
    inline Real norm      (T& a_Y, int a_ord)                 { return(a_Y.computeNorm(a_ord)); }

    inline void setShift        (Real a_shift )       { m_shift = a_shift;    }
    inline void setStageTime    (Real a_time  )       { m_StageTime = a_time; }
    inline void setStageNumber  (int a_Stage  )       { m_Stage = a_Stage;    }
    inline bool isDefined       ()              const { return m_is_Defined;  }
    inline void setBaseSolution (T& a_Y)              { m_Y.copy(a_Y); }

    void define(T& a_state, Ops& a_ops);

  private:
  
    bool m_is_Defined;
    int  m_Stage;
    Real m_shift, m_StageTime;
    T    m_Z, m_Y, m_F1, m_F2;
    Ops  *m_ops;
};

template <class T, class Ops>
void ImplicitStageJacobian<T,Ops>::define(T& a_state, Ops& a_ops)
{
  m_Z.define (a_state);
  m_Y.define (a_state);
  m_F1.define(a_state);
  m_F2.define(a_state);
  m_ops = &a_ops;
  m_is_Defined = true;
}

/* a_residual = [L]a_Y - a_rhs */
template <class T, class Ops>
void ImplicitStageJacobian<T,Ops>::residual(
                                              T& a_residual,
                                              const T& a_Y,
                                              const T& a_rhs,
                                              bool a_homogeneous = false
                                           )
{
  applyOp(a_residual,a_Y);
  a_residual.increment(a_rhs,-1);
}

template <class T, class Ops>
void ImplicitStageJacobian<T,Ops>::preCond(T& a_cor, const T& a_residual)
{
}

/* To compute: a_F = [L]a_y where [L] is the Jacobian of an implicit ARK stage,
 * i.e., [L] = [I - alpha*F'(Y)] where F(Y) is the part of the right-hand-side
 * being treated implicitly, and alpha is m_shift.
 * F'(Y) is approximated using a direction derivative, and its action on a vector
 * y is computed as: [F'(Y)]y ~ (1/eps)[F(Y+eps*y) - F(Y)].
 * Thus, the action of the Jacobian operator [L] on a vector y is:
 * [L]y = y - alpha * (1/eps) * [F(Y+eps*y) - F(Y)]
 *
 * For linear F(Y), eps = 1.
*/
template <class T, class Ops>
void ImplicitStageJacobian<T,Ops>::applyOp(T& a_F, const T& a_y, bool a_homogeneous=false)
{
  Real eps = 1.0;

  m_Z.copy(m_Y); m_Z.increment(a_y,eps);
  m_ops->implicitOpImEx(m_F1,m_StageTime,m_Y,m_Stage);
  m_ops->implicitOpImEx(m_F2,m_StageTime,m_Z,m_Stage);

  m_Z.copy(m_F2); m_Z.increment(m_F1,-1);
  a_F.copy(a_y); a_F.increment(m_Z,(-m_shift/eps));

}

template <class T, class Ops>
void ImplicitStageJacobian<T,Ops>::axby(T& a_Z, const T& a_X, const T& a_Y, Real a_a, Real a_b)
{
  a_Z.zero();
  a_Z.increment(a_X,a_a);
  a_Z.increment(a_Y,a_b);
}

#include "NamespaceFooter.H"
#endif

