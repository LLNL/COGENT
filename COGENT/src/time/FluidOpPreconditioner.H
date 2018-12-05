#ifndef _FLUID_OP_PC_H_
#define _FLUID_OP_PC_H_

#include <iostream>
#include <string>

#include "Dimensions.H"
#include "Preconditioner.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluidOpInterface.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#define _FLUID_OP_PC_ "FluidOpPC"

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;

template <class T, class Ops = void>
class FluidOpPreconditioner : public Preconditioner<T,Ops>
{
  public:

    FluidOpPreconditioner<T,Ops>()
    {
      m_ops = NULL;
      m_fluid_ops = NULL;
      m_shift = 0.0;
      m_species_idx = -1;
      m_is_defined = false;
      return;
    }

    ~FluidOpPreconditioner<T,Ops>()
    {
      return;
    }

    inline bool isDefined() const 
    { 
      return (m_is_defined && (m_ops != NULL)); 
    }

    void define(const T&, Ops&, const std::string&, const std::string&, bool)
    {
      MayDay::Error("Error in FluidOpPreconditioner::define() - Using disabled interface.\n");
      return;
    }

    void define(const T&, Ops&, CFG::FluidOpInterface&, const std::string&, const std::string&, bool);

    inline void update(const T& a_X, const Real a_time, const Real a_shift = 0.0)
    {
      CH_assert(isDefined());
      if (m_im) m_ops->assemblePCImEx(false, NULL, a_X, a_time, a_shift);
      else      m_ops->assembleSolPC(false, NULL, a_X, a_time);
      return;
    }

    inline void update( const KineticSpeciesPtrVect&    a_kin_species_phys,
                        const CFG::FluidSpeciesPtrVect& a_fluid_species,
                        const Real                      a_time,
                        const Real                      a_shift,
                        const bool                      a_im,
                        const int                       a_idx )
    {
      CH_assert(a_idx == m_species_idx);
      CH_assert(a_im == m_im);

      if (m_im) {
        m_fluid_ops->updatePCImEx(a_kin_species_phys, a_time, a_shift);
      }
      return;
    }

    inline void apply(T& a_Y, const T& a_X)
    {
      CH_assert(isDefined() && (m_species_idx >= 0));
      if (m_im) m_ops->solvePCImEx(a_Y, a_X, m_species_idx);
      else      m_ops->solveSolutionPC(a_Y, a_X, m_species_idx);
      return;
    }

    inline const std::string typeName() { return _FLUID_OP_PC_; }

    inline void speciesIndex(int a_idx) { m_species_idx = a_idx; }

  protected:

  private:

    bool                    m_is_defined;
    int                     m_species_idx;
    Ops*                    m_ops;
    CFG::FluidOpInterface*  m_fluid_ops;
    Real                    m_shift;
    bool                    m_im;

};

template <class T, class Ops>
void FluidOpPreconditioner<T,Ops>::define(const T&                a_x,
                                          Ops&                    a_ops,
                                          CFG::FluidOpInterface&  a_fluid_ops,
                                          const std::string&      a_out_string,
                                          const std::string&      a_opt_string,
                                          bool                    a_im )
{
  CH_assert(!isDefined());

  m_ops = &a_ops;
  m_fluid_ops = &a_fluid_ops;
  m_im = a_im;

  m_is_defined = true;
  return;
}

#include "NamespaceFooter.H"
#endif
