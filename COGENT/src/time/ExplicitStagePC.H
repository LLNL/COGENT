#ifndef _ExplicitStagePC_H_
#define _ExplicitStagePC_H_

#include <string>

#include "IMEXStagePC.H"
#include "BandedMatrix.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class ExplicitStagePC: public IMEXStagePCFunction<T,Ops>
{
  public:

    ExplicitStagePC<T,Ops>() {m_is_Defined = false; m_ops = NULL; }
    ~ExplicitStagePC<T,Ops>() {}

    void residual   (T& a_residual, const T& a_y, const T& a_rhs);
    void applyOp    (T& a_F, const T& a_y);
    void axby       (T& a_Z, const T& a_X, const T& a_Y, Real a_a, Real a_b);

    inline void preCond(T& a_Y, const T& a_X) { a_Y.copy(a_X); }

    inline void create    (T& a_Z, const T& a_Y)              { a_Z.define(a_Y); }
    inline void assign    (T& a_Z, const T& a_Y)              { a_Z.copy(a_Y); }
    inline void incr      (T& a_Z, const T& a_Y, Real a_scale){ a_Z.increment(a_Y,a_scale); } 
    inline void scale     (T& a_Y, const Real& a_scale)       { a_Y.scale(a_scale); }
    inline void setToZero (T& a_Y)                            { a_Y.zero(); }
    inline Real dotProduct(T& a_X, const T& a_Y)              { return(a_X.dotProduct(a_Y)); }
    inline Real norm      (T& a_Y, int a_ord)                 { return(a_Y.computeNorm(a_ord)); }

    inline void setStageTime    (Real a_time  )       { m_StageTime = a_time; }
    inline bool isDefined       ()              const { return m_is_Defined;  }

    inline void setBaseSolution (const T& a_Y)
    { 
      m_Y0.copy(a_Y); 
      m_normY0 = m_Y0.computeNorm(2); 
    }

    inline void setBaseRHS (const T& a_M)              
    { 
      m_M0.copy(a_M);
    }

    inline void setJFNKEps      (Real a_eps)          { m_epsJFNK = a_eps; }

    void define(const T&, Ops&, const std::string&, const std::string&);
    inline void setIsLinear(bool a_isLinear) { m_isLinear = a_isLinear; }

    inline void update(const T& a_X, const Real a_a = 0.0) 
    {
      m_P.zeroEntries();
      m_ops->assembleSolPC(true, (void*)&m_P, a_X);
    }

    inline void getDiagonal(std::vector<Real>& a_D) const
    {
      return m_P.getDiagonal(a_D);
    }

    inline const BandedMatrix& getBandedMatrixForm() const { return m_P; }
    inline BandedMatrix& getBandedMatrixForm() { return m_P; }

  private:
  
    bool m_is_Defined, m_isLinear;
    Real m_StageTime, m_epsJFNK, m_normY0;
    T    m_Z, m_Y0, m_M0, m_M;

    Ops *m_ops;
    BandedMatrix m_P;
    std::string m_pc_type;

    void parseParameters( ParmParse& );
};

template <class T, class Ops>
void ExplicitStagePC<T,Ops>::define(  const T&            a_state, 
                                      Ops&                a_ops, 
                                      const std::string&  a_outPrefix, 
                                      const std::string&  a_optPrefix )
{
  m_Z.define  (a_state);
  m_Y0.define (a_state);
  m_M0.define (a_state);
  m_M.define  (a_state);
  m_ops = &a_ops;
  CH_assert(m_ops != NULL);

  /* default values */
  m_epsJFNK = 1e-6;

  ParmParse pp_IPreCond( a_optPrefix.c_str() );
  parseParameters( pp_IPreCond );

  m_is_Defined = m_ops->preallocateSolPCMatrix((void*)&m_P, a_state);
}

template <class T,class Ops>
void ExplicitStagePC<T,Ops>::parseParameters( ParmParse& a_pp )
{
}

/* a_residual = [L]a_Y - a_rhs */
template <class T, class Ops>
void ExplicitStagePC<T,Ops>::residual(
                                       T& a_residual,
                                       const T& a_Y,
                                       const T& a_rhs )
{
  applyOp(a_residual,a_Y);
  a_residual.increment(a_rhs,-1);
}

template <class T, class Ops>
void ExplicitStagePC<T,Ops>::applyOp(T& a_F, const T& a_y)
{
  Real normY = a_y.computeNorm(2);
  if (normY < 1e-15) a_F.zero();
  else {
    Real eps;
    if (m_isLinear) {
      eps = 1.0;
    } else {
      eps = m_epsJFNK * sqrt(1 + m_normY0) / normY;
    }

    m_Z.copy(m_Y0); m_Z.increment(a_y,eps);
    m_ops->solutionPC(m_M,m_StageTime,m_Z);

    a_F.copy(m_M); 
    a_F.increment(m_M0,-1);
    a_F.scale(1.0/eps);
  }
}

template <class T, class Ops>
void ExplicitStagePC<T,Ops>::axby(T& a_Z, const T& a_X, const T& a_Y, Real a_a, Real a_b)
{
  a_Z.zero();
  a_Z.increment(a_X,a_a);
  a_Z.increment(a_Y,a_b);
}

#include "NamespaceFooter.H"
#endif

