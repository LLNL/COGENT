#ifndef _GK_VLASOV_PC_H_
#define _GK_VLASOV_PC_H_

#include "Preconditioner.H"
#include "GKVlasovAMG.H"

#define _GK_VLASOV_PC_ "GKVlasovPC"

#include "NamespaceHeader.H"


template <class T, class Ops = void>
class GKVlasovPreconditioner : public Preconditioner<T,Ops>
{
  public:

    GKVlasovPreconditioner<T,Ops>()
       : m_is_defined(false),
         m_species_idx(-1),
         m_ops(NULL),
         m_shift(0.0),
         m_amg_pc(NULL) {}

    ~GKVlasovPreconditioner<T,Ops>() 
    {
      delete m_amg_pc;
    }

    inline bool isDefined() const 
    { 
      return (      m_is_defined 
                &&  (m_ops != NULL) 
                &&  (m_species_idx >= 0) ); 
    }

    void define(  const T&, 
                  void*, 
                  const std::string&, 
                  const std::string&, 
                  bool  );

    void defineAMGPC(const ParmParse&, const PhaseGeom& );

    inline void apply(T& a_Y, const T& a_X)
    {
      CH_assert(isDefined());
      CH_assert(m_im);
      m_ops->solveVlasovOpPCImEx( this,
                                  a_Y, 
                                  a_X, 
                                  Preconditioner<T,Ops>::getSysID(),
                                  m_species_idx );
      return;
    }

    inline const std::string typeName() { return _GK_VLASOV_PC_; }

    inline void speciesIndex(int a_idx) { m_species_idx = a_idx; }

    inline GKVlasovAMG* getAMGPC() { return m_amg_pc; }

  protected:

  private:

    bool  m_is_defined;
    int   m_species_idx;
    Ops*  m_ops;
    Real  m_shift;
    bool  m_im;

    GKVlasovAMG *m_amg_pc;
};


template <class T, class Ops>
void GKVlasovPreconditioner<T,Ops>::define(const T&            a_x,
                                           void*               a_ops,
                                           const std::string&  a_out_string,
                                           const std::string&  a_opt_string,
                                           bool                a_im )
{
  CH_assert(!isDefined());

  m_ops = (Ops*) a_ops;
  m_im = a_im;

  m_is_defined = false; //not yet
}

template <class T, class Ops>
void GKVlasovPreconditioner<T,Ops>::defineAMGPC(const ParmParse&  a_pp,
                                                const PhaseGeom&  a_geom )
{
  CH_assert(!isDefined());
  m_amg_pc = new GKVlasovAMG(a_pp, a_geom);
  m_is_defined = true; //now it's defined
}


#include "NamespaceFooter.H"
#endif
