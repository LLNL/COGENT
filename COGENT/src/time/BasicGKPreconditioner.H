#ifndef _BasicGKPC_H_
#define _BasicGKPC_H_

#include <iostream>
#include <string>

#include "Preconditioner.H"
#include "ExplicitStagePC.H"
#include "ImplicitStagePC.H"
#include "LinearFunctionSolverLibrary.H"

#define _BASIC_GK_PC_ "basicGKPC"

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class BasicGKPreconditioner : public Preconditioner<T,Ops>
{
  public:

    BasicGKPreconditioner<T,Ops>()
    {
      m_dof_list = DOFList(0);
      m_is_defined = false;
      return;
    }
    ~BasicGKPreconditioner<T,Ops>()
    {
      delete m_pc_object;
      delete m_pc_solver;
      return;
    }

    inline bool isDefined() const { return m_is_defined; }

    inline void define( const T&            a_X, 
                        Ops&                a_ops, 
                        const std::string&  a_out_str, 
                        const std::string&  a_opt_str, 
                        bool                a_im )
    {
      DOFList tmp = DOFList(0);
      define(a_X,a_ops,a_out_str,a_opt_str,a_im,tmp);
    }
    void define(const T&, Ops&, const std::string&, const std::string&, bool, const DOFList&);

    inline void update(const T& a_X, const Real a_time, const Real a_shift = 0.0)
    {
      CH_assert(isDefined());
      if (!procID()) {
        std::cout << "    ==> (" 
                  << m_out_string << ") Updating preconditioning matrix.\n";
      }
      m_pc_object->update(a_X, a_time, a_shift);
      return;
    }

    inline void apply(T& a_Y, const T& a_X)
    {
      CH_assert(isDefined());
      if (m_dof_list.size() == 0) {
        m_pc_solver->solve(a_Y, a_X);
      } else {
        m_Z = a_Y;
        m_pc_solver->solve(m_Z,a_X);
        a_Y.copy(m_Z, m_dof_list);
      }
      return;
    }

    inline void setShift (Real a_a) { m_pc_object->setShift(a_a); }
    inline void setTime(Real a_t) { m_pc_object->setStageTime(a_t); }

    inline void setBaseSolution(const T& a_x) 
    { 
      m_pc_object->setBaseSolution(a_x); 
    }

    inline void setBaseRHS(const T& a_F, const T& a_M) 
    { 
      m_pc_object->setBaseRHS(a_F, a_M); 
    }

    inline void setBaseRHS(const T& a_M) 
    { 
      m_pc_object->setBaseRHS(a_M); 
    }

    inline void setJFNKEps(Real a_eps) { m_pc_object->setJFNKEps(a_eps); }
    inline void setIsLinear(bool a_b) { m_pc_object->setIsLinear(a_b); }

    inline const std::string typeName() { return _BASIC_GK_PC_; }

    inline void* getBandedMatrix() { return((void*)(&(m_pc_object->getBandedMatrixForm()) ) ); }

  protected:

  private:
    
    bool m_is_defined;
    std::string m_pc_type;
    std::string m_out_string, m_opt_string;
    DOFList m_dof_list;
    T m_Z;

    IMEXStagePCFunction<T,Ops>  *m_pc_object;
    LinearFunctionSolver<T,Ops> *m_pc_solver;
};

template <class T, class Ops>
void BasicGKPreconditioner<T,Ops>::define(const T&            a_x,
                                          Ops&                a_ops,
                                          const std::string&  a_out_string,
                                          const std::string&  a_opt_string,
                                          bool                a_im,
                                          const DOFList&      a_dof_list)
{
  CH_assert(!isDefined());
  /* default values */
  m_pc_type = "gauss_seidel";

  m_out_string = a_out_string + "." + _BASIC_GK_PC_;
  m_opt_string = a_opt_string + "_" + _BASIC_GK_PC_;

  /* read inputs */
  ParmParse pp_pc(m_opt_string.c_str());
  pp_pc.query("solver_type", m_pc_type);

  if (a_im) m_pc_object = new ImplicitStagePC<T,Ops>;
  else      m_pc_object = new ExplicitStagePC<T,Ops>;
  m_pc_object->define(a_x, a_ops, m_out_string, m_opt_string);

  if      (m_pc_type == "gmres" )       m_pc_solver = new GMRESSolver<T,Ops>;
  else if (m_pc_type == "jacobi")       m_pc_solver = new PointJacobiSolver<T,Ops>;
  else if (m_pc_type == "twolevel")     m_pc_solver = new TwoLevelSolver<T,Ops>;
  else if (m_pc_type == "gauss_seidel") m_pc_solver = new GaussSeidelSolver<T,Ops>;
  else {
     std::cout << "Error in BasicGKPreconditioner::define() - " 
               << m_pc_type << " "
               << "is not a valid type for preconditioning solver.\n";
    m_pc_solver = NULL;
  }
  CH_assert(m_pc_solver != NULL);
  m_pc_solver->define(m_pc_object, m_out_string, m_opt_string, a_x);
  
  m_Z.define(a_x);
  m_dof_list = a_dof_list;
  m_is_defined = true;
  return;
}

#include "NamespaceFooter.H"
#endif

