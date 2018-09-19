#ifndef _AMG_GK_PC_H_
#define _AMG_GK_PC_H_

#include <iostream>
#include <string>

#include "Preconditioner.H"

#define _AMG_GK_PC_ "amgGKPC"

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class AMGGKPreconditioner : public Preconditioner<T,Ops>
{
  public:

    AMGGKPreconditioner<T,Ops>()
    {
      m_ops = NULL;
      m_shift = 0.0;
      m_is_defined = false;
      return;
    }

    ~AMGGKPreconditioner<T,Ops>()
    {
      return;
    }

    inline bool isDefined() const 
    { 
      return (m_is_defined && (m_ops != NULL)); 
    }

    void define(const T&, Ops&, const std::string&, const std::string&, bool);

    inline void update(const T& a_X, const Real a_shift = 0.0)
    {
      CH_assert(isDefined());
      if (m_im) m_ops->assemblePCImEx(false,NULL,a_X,a_shift);
      else      m_ops->assembleSolPC(false,NULL,a_X);
      return;
    }

    inline void apply(T& a_Y, const T& a_X)
    {
      CH_assert(isDefined());
      if (m_im) m_ops->solvePCImEx(a_Y, a_X);
      else      m_ops->solveSolutionPC(a_Y, a_X);
      return;
    }

  protected:

  private:

    bool  m_is_defined;
    Ops*  m_ops;
    Real  m_shift;
    bool  m_im;

};

template <class T, class Ops>
void AMGGKPreconditioner<T,Ops>::define(const T&            a_x,
                                        Ops&                a_ops,
                                        const std::string&  a_out_string,
                                        const std::string&  a_opt_string,
                                        bool                a_im )
{
  CH_assert(!isDefined());

  m_ops = &a_ops;
  m_im = a_im;

  m_is_defined = true;
  return;
}

#include "NamespaceFooter.H"
#endif
