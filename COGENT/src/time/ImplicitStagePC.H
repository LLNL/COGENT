#ifndef _ImplicitStagePC_H_
#define _ImplicitStagePC_H_

#include <string>

#include "IMEXStagePC.H"
#include "BandedMatrix.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class ImplicitStagePC: public IMEXStagePCFunction<T,Ops>
{
  public:

    ImplicitStagePC<T,Ops>() 
    {
      m_is_Defined = false; 
      m_shift = 0; 
      m_ops = NULL; 
    }
    ~ImplicitStagePC<T,Ops>() {}

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

    inline bool isDefined() const { return m_is_Defined;  }

    inline void setJFNKEps(Real a_eps) { m_epsJFNK = a_eps; }

    inline void setBaseSolution (const T& a_Y)              
    { 
      m_Y0.copy(a_Y); 
      m_normY0 = m_Y0.computeNorm(2); 
    }

    inline void setBaseRHS (const T& a_F, const T& a_M)      
    { 
      m_F0.copy(a_F);
      m_M0.copy(a_M); 
    }

    void define(const T&, Ops&, const std::string&, const std::string&);
    inline void setIsLinear(bool a_isLinear) { m_isLinear = a_isLinear; }

    inline void update(const T& a_X, const Real a_time, const Real a_shift) 
    {
       m_shift = a_shift;
       m_stage_time = a_time;
       m_ops->assemblePCImEx(true, (void*)&m_Pmat, a_X, a_time, a_shift);
    }

    inline void getDiagonal(std::vector<Real>& a_D) const
    {
      return m_Pmat.getDiagonal(a_D);
    }

    inline const BandedMatrix& getBandedMatrixForm() const { return m_Pmat; }
    inline BandedMatrix& getBandedMatrixForm() { return m_Pmat; }

  private:
  
    bool m_is_Defined, m_isLinear;
    Real m_shift, m_stage_time, m_epsJFNK, m_normY0;
    T    m_Z, m_Y0, m_F0, m_F, m_M0, m_M;

    Ops *m_ops;
    BandedMatrix m_Pmat;
   
    void parseParameters( ParmParse& );
};

template <class T, class Ops>
void ImplicitStagePC<T,Ops>::define(  const T&            a_state, 
                                      Ops&                a_ops, 
                                      const std::string&  a_outPrefix, 
                                      const std::string&  a_optPrefix )
{
  m_Z.define  (a_state);
  m_Y0.define (a_state);
  m_F0.define (a_state);
  m_F.define  (a_state);
  m_M0.define (a_state);
  m_M.define  (a_state);
  m_ops = &a_ops;
  CH_assert(m_ops != NULL);

  /* default values */
  m_epsJFNK = 1e-6;

  ParmParse pp_IPreCond( a_optPrefix.c_str() );
  parseParameters( pp_IPreCond );

  m_ops->preallocatePCMatrixImEx((void*)&m_Pmat, a_state);
  m_is_Defined = true;
  return;
}

template <class T,class Ops>
void ImplicitStagePC<T,Ops>::parseParameters( ParmParse& a_pp )
{
}

/* a_residual = [L]a_Y - a_rhs */
template <class T, class Ops>
void ImplicitStagePC<T,Ops>::residual(
                                        T& a_residual,
                                        const T& a_Y,
                                        const T& a_rhs
                                     )
{
  applyOp(a_residual,a_Y);
  a_residual.increment(a_rhs,-1);
}

template <class T, class Ops>
void ImplicitStagePC<T,Ops>::applyOp(T& a_F, const T& a_y)
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
    m_ops->implicitPCImEx(m_F,m_stage_time,m_Z);
    if (m_ops->trivialSolutionOp()) {
      m_M.copy(m_Z);
    } else {
      m_ops->solutionPC(m_M,m_stage_time,m_Z);
    }

    m_Z.copy(m_F); 
    m_Z.increment(m_F0,-1);

    a_F.copy(m_M);
    a_F.increment(m_M0,-1);
    a_F.scale(m_shift/eps);
    a_F.increment(m_Z,(-1.0/eps));
  }
}

template <class T, class Ops>
void ImplicitStagePC<T,Ops>::axby(T& a_Z, const T& a_X, const T& a_Y, Real a_a, Real a_b)
{
  a_Z.zero();
  a_Z.increment(a_X,a_a);
  a_Z.increment(a_Y,a_b);
}

#include "NamespaceFooter.H"
#endif

