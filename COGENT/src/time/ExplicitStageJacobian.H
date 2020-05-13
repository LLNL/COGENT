#ifndef _ExplicitStageJacobian_H_
#define _ExplicitStageJacobian_H_

#include "PreconditionerLibrary.H"
#include "Checkpoints.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class ExplicitStageJacobian: public LinearFunction<T,Ops>
{
  public:

    ExplicitStageJacobian<T,Ops>() {m_is_Defined = false; }
    ~ExplicitStageJacobian<T,Ops>() 
    {
      delete m_preCond;
      return;
    }

    void residual   (T& a_residual, const T& a_y, const T& a_rhs);
    void applyOp    (T& a_F, const T& a_y);
    void axby       (T& a_Z, const T& a_X, const T& a_Y, Real a_a, Real a_b);

    inline void preCond(T& a_Y, const T& a_X)     
    { 
      if (m_usePreCond) {
        a_Y.zero();
        m_preCond->apply(a_Y, a_X);
      } else {
        a_Y.copy(a_X);
      }
      return;
    }

    inline void updatePreCondMat(const T& a_X)
    {
      if (m_usePreCond) m_preCond->update(a_X, m_stage_time);
      return;
    }

    inline void create    (T& a_Z, const T& a_Y)              { a_Z.define(a_Y); }
    inline void assign    (T& a_Z, const T& a_Y)              { a_Z.copy(a_Y); }
    inline void incr      (T& a_Z, const T& a_Y, Real a_scale){ a_Z.increment(a_Y,a_scale); } 
    inline void scale     (T& a_Y, const Real& a_scale)       { a_Y.scale(a_scale); }
    inline void setToZero (T& a_Y)                            { a_Y.zero(); }
    inline Real dotProduct(T& a_X, const T& a_Y)              { return(a_X.dotProduct(a_Y)); }
    inline Real norm      (T& a_Y, int a_ord)                 { return(a_Y.computeNorm(a_ord)); }

    inline void setStageTime (Real a_time  ) { m_stage_time = a_time; }

    inline bool isDefined() const { return m_is_Defined;  }

    inline void setBaseSolution (const T& a_Y)
    { 
      m_Y0.copy(a_Y); 
      m_normY0 = m_Y0.computeNorm(2); 
      if (m_usePreCond) m_preCond->setBaseSolution(a_Y);
    }

    inline void setBaseRHS (const T& a_M)
    { 
      m_M0.copy(a_M); 
      if (m_usePreCond) m_preCond->setBaseRHS(a_M);
    }

    inline void setJFNKEps (Real a_eps)
    { 
      m_epsJFNK = a_eps; 
      if (m_usePreCond) m_preCond->setJFNKEps(a_eps);
    }

    void setIsLinear (bool a_isLinear) 
    { 
      m_isLinear = a_isLinear; 
      if (m_usePreCond) m_preCond->setIsLinear(true);
    }

    void define(const T&, Ops&, const std::string&, const std::string&);

  private:
  
    bool m_is_Defined, m_isLinear, m_usePreCond;
    Real m_stage_time, m_epsJFNK, m_normY0;
    T    m_Z, m_Y0, m_M0, m_M;

    Preconditioner<T,Ops> *m_preCond;

    Ops *m_ops;
    std::string m_optPrefix, m_outPrefix, m_pc_type;
    
    void parseParameters( ParmParse& );
};

template <class T, class Ops>
void ExplicitStageJacobian<T,Ops>::define(  const T&            a_state, 
                                            Ops&                a_ops,
                                            const std::string&  a_outPrefix,
                                            const std::string&  a_optPrefix )
{
  m_Z.define  (a_state);
  m_Y0.define (a_state);
  m_M0.define (a_state);
  m_M.define  (a_state);
  m_ops = &a_ops;

  /* default values */
  m_epsJFNK = 1e-6;
  m_usePreCond = true;
  m_pc_type = _MULTIPHYSICS_PC_;

  m_outPrefix = a_outPrefix;
  m_optPrefix = a_optPrefix;
  std::string optString = a_optPrefix + "_jacobian";
  ParmParse ppIJac( optString.c_str() );
  parseParameters( ppIJac );

  m_preCond = NULL;
  if (m_usePreCond) {

     std::string pc_outPrefix = a_outPrefix + ".PC";
     std::string pc_optPrefix = a_optPrefix + "_pc";

     ParmParse pp_pc(pc_optPrefix.c_str());
     pp_pc.query("type", m_pc_type);

    if (m_pc_type == _BASIC_GK_PC_) {
      m_preCond = new BasicGKPreconditioner<T,Ops>;
    } else if (m_pc_type == _GK_VLASOV_PC_) {
      m_preCond = new GKVlasovPreconditioner<T,Ops>;
    } else if (m_pc_type == _MULTIPHYSICS_PC_) {
      m_preCond = new MultiPhysicsPreconditioner<T,Ops>;
    } else {
      std::cout << "Error in ExplicitStageJacobian::define() - "
                << m_pc_type << " "
                << "is not a valid type of preconditioner.\n";
      m_preCond = NULL;
    }
    CH_assert(m_preCond != NULL);
    m_preCond->define(a_state, a_ops, pc_outPrefix, pc_optPrefix, false);
  }

  m_is_Defined = true;
  return;
}

template <class T, class Ops>
void ExplicitStageJacobian<T,Ops>::parseParameters( ParmParse& a_pp)
{
  a_pp.query("jfnk_eps", m_epsJFNK);
  a_pp.query("with_pc",  m_usePreCond);
}

/* a_residual = [L]a_Y - a_rhs */
template <class T, class Ops>
void ExplicitStageJacobian<T,Ops>::residual(
                                              T& a_residual,
                                              const T& a_Y,
                                              const T& a_rhs
                                           )
{
  applyOp(a_residual,a_Y);
  a_residual.increment(a_rhs,-1);
}

/* To compute: a_F = [L]a_y where 
 * [L] = M'(Y),
 * M(Y) is the mass matrix function.
 * M'(Y) is approximated using a direction derivative, 
 * and its action on a vector y is computed as: 
 * [M'(Y)]y ~ (1/eps)[M(Y+eps*y) - M(Y)].
 * Thus, the action of the Jacobian operator [L] on a vector y is:
 * [L]y = (1/eps) * [M(Y+eps*y)-M(Y)]
 *
 * For linear problems, eps = 1.*/
template <class T, class Ops>
void ExplicitStageJacobian<T,Ops>::applyOp(T& a_F, const T& a_y)
{
  CH_TIME("ExplicitStageJacobian::apllyOp()");
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
    m_ops->preOpEval(m_Z,m_stage_time,chkpt_stage_jac);
    m_ops->solutionOp(m_M,m_stage_time,m_Z);

    a_F.copy(m_M); 
    a_F.increment(m_M0,-1);
    a_F.scale(1.0/eps);
  }
}

template <class T, class Ops>
void ExplicitStageJacobian<T,Ops>::axby(T& a_Z, const T& a_X, const T& a_Y, Real a_a, Real a_b)
{
  a_Z.zero();
  a_Z.increment(a_X,a_a);
  a_Z.increment(a_Y,a_b);
}

#include "NamespaceFooter.H"
#endif

