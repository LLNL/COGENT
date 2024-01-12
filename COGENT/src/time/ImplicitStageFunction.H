#ifndef _ImplicitStageFunction_H_
#define _ImplicitStageFunction_H_

#include "Function.H"
#include "ImplicitStageJacobian.H"
#include "TiDefinitions.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class ImplicitStageFunction: public Function<T,Ops>
{
  public:
    ImplicitStageFunction<T,Ops>() {m_is_Defined = false; m_shift = 0; }
    ~ImplicitStageFunction<T,Ops>() {}

    void evalFunction(T& a_F, const T& a_Y, const T& a_b, int a_its );
    void define(const T& a_state, Ops& a_ops);

    bool checkLinearity( const T& );

    inline void setShift        (Real a_shift )       { m_shift = a_shift;    }
    inline void setStageTime    (Real a_time  )       { m_StageTime = a_time; }
    inline bool isDefined       ()              const { return m_is_Defined;  }

    inline void setJac(ImplicitStageJacobian<T,Ops>* a_jac) { m_jac = a_jac; }

    inline Real computeStepNorm(const ODEVector& a_vec) const
    {
      return m_ops->scaledVecNormU(a_vec);
    }

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
void ImplicitStageFunction<T,Ops>::evalFunction(  T&        a_F, 
                                                  const T&  a_Y, 
                                                  const T&  a_b, 
                                                  int       a_its )
{
  CH_TIME("ImplicitStageFunction::evalFunction()");
  CH_assert(isDefined());

  m_ops->preOpEval( a_Y, 
                    m_StageTime, 
                    (a_its == 0 ? chkpt_stage_func_0 : chkpt_stage_func_n) );

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

  if (a_its == 0) {
    m_ops->computeUScale(a_Y);
    m_ops->computeRHSOpScale(a_F);
  }

  m_ops->scaleRHSOp(a_F);

  return;
}

template <class T, class Ops>
bool ImplicitStageFunction<T,Ops>::checkLinearity(  const T& a_Y)
{
  const double pi = 4.0*atan(1.0);

  T rhsvec(a_Y);
  rhsvec.zero();

  double Ynorm = a_Y.computeNorm();
  if (!procID()) {
    printf("    |ImplicitStageFunction::checkLinearity():\n");
  }

  for (double eps = 1e-2; eps >= 1e-14; eps *= 1e-2) {
    double reps = eps;
    if (Ynorm > 1e-14) reps *= Ynorm;

    T Z(a_Y), Fy(a_Y), Fz(a_Y), Fyz(a_Y);

    Z = a_Y;
    for (int i=0; i<Z.getVectorSize(); i++) {
      Z[i] += reps * cos(2*pi*(((double)i)/((double)Z.getVectorSize())));
    }

    evalFunction(Fy, a_Y, rhsvec, 0);
    evalFunction(Fz,   Z, rhsvec, 0);

    Z += a_Y;
    evalFunction(Fyz, Z, rhsvec, 0);

    double Fnorm = Fy.computeNorm();
    Fyz -= Fy;
    Fyz -= Fz;
    double deltaNorm = Fyz.computeNorm();

    if (!procID()) {
      printf("    |  eps=%1.1e, norm of F(Y): %1.4e, norm of F(Y+Z)-[F(Y)+F(Z)]: %1.4e\n",
             reps, Fnorm, deltaNorm );
    }
  }

  return false;
}

#include "NamespaceFooter.H"
#endif
