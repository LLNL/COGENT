#ifndef _GK_VLASOV_PC_H_
#define _GK_VLASOV_PC_H_

#define _GK_VLASOV_PC_ "GKVlasovPC"

#include "NamespaceHeader.H"


template <class T, class Ops = void>
class GKVlasovPreconditioner : public Preconditioner<T,Ops>
{
  public:

    GKVlasovPreconditioner<T,Ops>()
       : m_ops(NULL),
         m_shift(0.0),
         m_species_idx(-1),
         m_is_defined(false) {};

    ~GKVlasovPreconditioner<T,Ops>() {};

    inline bool isDefined() const 
    { 
      return (m_is_defined && (m_ops != NULL)); 
    }

    void define(const T&, Ops&, const std::string&, const std::string&, bool);

    inline void update(const T& a_X, const Real a_time, const Real a_shift = 0.0)
    {
      CH_assert(isDefined());
      if (m_im) m_ops->assembleVlasovOpPCImEx(NULL, a_X, a_time, a_shift);
      else      m_ops->assembleSolPC(false, NULL, a_X, a_time);
      return;
    }

    inline void apply(T& a_Y, const T& a_X)
    {
      CH_assert(isDefined());
      if (m_im) m_ops->solveVlasovOpPCImEx(a_Y, a_X, m_species_idx);
      else      m_ops->solveSolutionPC(a_Y, a_X, m_species_idx);
      return;
    }

    inline const std::string typeName() { return _GK_VLASOV_PC_; }

    inline void speciesIndex(int a_idx) { m_species_idx = a_idx; }

  protected:

  private:

    bool  m_is_defined;
    int   m_species_idx;
    Ops*  m_ops;
    Real  m_shift;
    bool  m_im;
};


template <class T, class Ops>
void GKVlasovPreconditioner<T,Ops>::define(const T&            a_x,
                                           Ops&                a_ops,
                                           const std::string&  a_out_string,
                                           const std::string&  a_opt_string,
                                           bool                a_im )
{
  CH_assert(!isDefined());

  m_ops = &a_ops;
  m_im = a_im;

  m_is_defined = true;
}


#include "NamespaceFooter.H"
#endif
