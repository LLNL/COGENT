#ifndef _MultiPhysicsPreconditioner_H_
#define _MultiPhysicsPreconditioner_H_

#include <vector>
#include <utility>
#include <string>
#include "Preconditioner.H"

#define _MULTIPHYSICS_PC_ "multiphysicsPC"

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class MultiPhysicsPreconditioner : public Preconditioner<T,Ops>
{
  public:

    MultiPhysicsPreconditioner<T,Ops>() 
    { 
      m_is_defined = false; 
      m_im = false;
      m_ops = NULL;
      m_shift = -DBL_MAX;
      m_dof_list.clear();
      m_block_pc.clear(); 
    }
    ~MultiPhysicsPreconditioner<T,Ops>();

    virtual void define(const T&,
                        Ops&,
                        const std::string&,
                        const std::string&,
                        bool );

    virtual void update(const T&, const Real, const Real a_a = 0.0);
    virtual void apply(T&, const T&);
    virtual bool isDefined() const { return m_is_defined; }

    inline virtual void setBaseSolution(const T& a_x)
    {
      for (int i=0; i<m_block_pc.size(); i++) {
        m_block_pc[i]->setBaseSolution(a_x);
      }
    }

    inline virtual void setBaseRHS(const T& a_F, const T& a_M)
    {
      for (int i=0; i<m_block_pc.size(); i++) {
        m_block_pc[i]->setBaseRHS(a_F, a_M);
      }
    }

    inline virtual void setBaseRHS(const T& a_M)
    {
      for (int i=0; i<m_block_pc.size(); i++) {
        m_block_pc[i]->setBaseRHS(a_M);
      }
    }

    inline virtual void setJFNKEps(Real a_e)
    {
      for (int i=0; i<m_block_pc.size(); i++) {
        m_block_pc[i]->setJFNKEps(a_e);
      }
    }

    inline virtual void setIsLinear(bool a_b)
    {
      for (int i=0; i<m_block_pc.size(); i++) {
        m_block_pc[i]->setIsLinear(a_b);
      }
    }

    inline const std::string typeName() { return _MULTIPHYSICS_PC_; }

  protected:

    bool m_is_defined;
    bool m_im;
    Real m_shift;
    Real m_time;

    std::vector<DOFList> m_dof_list;
    std::vector<Preconditioner<T,Ops>*> m_block_pc;

    std::string m_out_string, m_opt_string;

    Ops* m_ops;

  private:
};

template <class T, class Ops>
MultiPhysicsPreconditioner<T,Ops>::~MultiPhysicsPreconditioner()
{
  for (int i=0; i<m_block_pc.size(); i++) {
    delete(m_block_pc[i]);
  }
  m_block_pc.clear();
  m_dof_list.clear();
  return;
}

template <class T, class Ops>
void MultiPhysicsPreconditioner<T,Ops>::define( const T&            a_x,
                                                Ops&                a_ops,
                                                const std::string&  a_out_string,
                                                const std::string&  a_opt_string,
                                                bool                a_im )
{
  CH_assert(!isDefined());

  m_out_string = a_out_string + "." + _MULTIPHYSICS_PC_;
  m_opt_string = a_opt_string + "_" + _MULTIPHYSICS_PC_;

  m_ops = &a_ops;
  m_im = a_im;

  m_ops->defineMultiPhysicsPC(m_block_pc, m_dof_list, a_x, m_out_string, m_opt_string, a_im);
  CH_assert(m_block_pc.size() == m_dof_list.size());

  m_is_defined = true;
  for (int i=0; i<m_block_pc.size(); i++) {
    if (m_block_pc[i] == NULL) {
      MayDay::Error("Error in MultiPhysicsPreconditioner::define(): NULL pointer to preconditioner!");
    }
    m_is_defined = (m_is_defined && m_block_pc[i]->isDefined());
  }

  CH_assert(isDefined());

  if (!procID()) {
    std::cout << "\n";
    std::cout << _MULTIPHYSICS_PC_ << " defined "
              << "for " << (a_im ? "RHS " : "LHS ")
              << "with "
              << m_block_pc.size() << " block preconditioners: "
              << "\n";
    for (int i=0; i<m_block_pc.size(); i++) {
      std::cout << "  " << i << ": " << m_block_pc[i]->typeName()
                << "\n";
    }
    std::cout << "\n";
  }
}

template <class T, class Ops>
void MultiPhysicsPreconditioner<T,Ops>::update( const T&    a_x,
                                                const Real  a_t,
                                                const Real  a_a )
{
  CH_assert(isDefined());
  if (m_im) CH_assert(a_a > 0);

  m_time = a_t;
  m_shift = a_a;

  m_ops->updateMultiPhysicsPC(m_block_pc, a_x, a_t, a_a, m_im);

  return;
}

template <class T, class Ops>
void MultiPhysicsPreconditioner<T,Ops>::apply(T&        a_x,
                                              const T&  a_b )
{
  CH_assert(isDefined());

  /* default: identity matrix as preconditioner */
  a_x = a_b; 

  if (m_block_pc.size() > 1) {

    if (m_im) {
  
      CH_assert(m_shift > 0);
  
      /* On the RHS, apply the block preconditioners
       * in a general way (no assumptions of overlap)*/
  
      T z(a_x); /* work vector */
      for (int i=0; i<m_block_pc.size(); i++) {
    
        Preconditioner<T,Ops> *pc(m_block_pc[i]);
    
        if (i == 0) {
  
          a_x = a_b;
          a_x.scale(1.0/m_shift);
  
          z = a_b;
          z.scale(1.0/m_shift);
          m_ops->solveSolutionPC(a_x, z);
  
          pc->apply(a_x, a_b);
  
        } else {
  
          z = a_x;
          m_ops->solutionPC(z, m_time, a_x);
          z.scale(m_shift);
  
          pc->apply(a_x, z);
  
        }
    
      }
  
    } else {
  
      /* On the LHS, assume nontrivial solution 
       * operators are nonoverlapping */
  
      T z(a_x);
      for (int i=0; i<m_block_pc.size(); i++) {
        
        Preconditioner<T,Ops> *pc(m_block_pc[i]);
        const DOFList& dof_list(m_dof_list[i]);
        
        z = a_b;
        pc->apply(z, a_b);
        
        a_x.copy(z, dof_list);
      }
  
    }

  } else if (m_block_pc.size() == 1) {

    Preconditioner<T,Ops> *pc(m_block_pc[0]);
    pc->apply(a_x, a_b);

  }

  return;
}

#include "NamespaceFooter.H"

#endif
