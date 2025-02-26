#ifndef _BasicGKPC_H_
#define _BasicGKPC_H_

#include <iostream>
#include <string>

#include "BandedMatrix.H"
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
                        void*               a_ops, 
                        const std::string&  a_out_str, 
                        const std::string&  a_opt_str, 
                        bool                a_im )
    {
      DOFList tmp = DOFList(0);
      define( a_X,
              a_ops,
              a_out_str,
              a_opt_str,
              a_im,
              tmp  );
    }

    void define(  const T&, 
                  void*, 
                  const std::string&, 
                  const std::string&, 
                  bool, 
                  const DOFList& );

    inline void apply(T& a_Y, const T& a_X)
    {
      CH_assert(isDefined());
      T X, Y;
      m_ops->getVec( a_X, Preconditioner<T,Ops>::getSysID(), X);
      m_ops->getVec( a_Y, Preconditioner<T,Ops>::getSysID(), Y);

      if (m_dof_list.size() == 0) {
        m_pc_solver->solve(Y, X);
      } else {
        m_Z = Y;
        m_pc_solver->solve(m_Z, X);
        Y.copy(m_Z, m_dof_list);
      }
      return;
    }

    inline void setShift (Real a_a) { m_pc_object->setShift(a_a); }
    inline void setTime(Real a_t) { m_pc_object->setStageTime(a_t); }

    inline void setBaseSolution(const T& a_x) 
    { 
      T x;
      m_ops->getVec( a_x, Preconditioner<T,Ops>::getSysID(), x );
      m_pc_object->setBaseSolution(x); 
    }

    inline void setBaseRHS(const T& a_F, const T& a_M) 
    { 
      T F, M;
      m_ops->getVec( a_F, Preconditioner<T,Ops>::getSysID(), F);
      m_ops->getVec( a_M, Preconditioner<T,Ops>::getSysID(), M);
      m_pc_object->setBaseRHS(F, M); 
    }

    inline void setBaseRHS(const T& a_M) 
    { 
      T M;
      m_ops->getVec( a_M, Preconditioner<T,Ops>::getSysID(), M);
      m_pc_object->setBaseRHS(M); 
    }

    inline void setJFNKEps(Real a_eps) { m_pc_object->setJFNKEps(a_eps); }
    inline void setIsLinear(bool a_b) { m_pc_object->setIsLinear(a_b); }

    inline const std::string typeName() { return _BASIC_GK_PC_; }

    inline BandedMatrix& getBandedMatrix() 
    { 
      return(m_pc_object->getBandedMatrixForm()); 
    }

  protected:

  private:
    
    bool m_is_defined;
    std::string m_pc_type;
    std::string m_out_string, m_opt_string;
    DOFList m_dof_list;
    T m_Z;
    Ops* m_ops;

    IMEXStagePCFunction<T,Ops>  *m_pc_object;
    LinearFunctionSolver<T,Ops> *m_pc_solver;
};

template <class T, class Ops>
void BasicGKPreconditioner<T,Ops>::define(const T&            a_x,
                                          void*               a_ops,
                                          const std::string&  a_out_string,
                                          const std::string&  a_opt_string,
                                          bool                a_im,
                                          const DOFList&      a_dof_list)
{
  CH_assert(!isDefined());

  m_ops = (Ops*) a_ops;

  if (a_im) m_pc_object = new ImplicitStagePC<T,Ops>;
  else      m_pc_object = new ExplicitStagePC<T,Ops>;
  m_pc_object->define(a_x, *m_ops, m_out_string, m_opt_string);

  m_Z.define(a_x);
  m_dof_list = a_dof_list;

  /* read inputs */
  m_pc_type = "gauss_seidel";
  
  m_opt_string = a_opt_string + "_basicGKPC";
  m_out_string = a_out_string + ".basicGKPC";
  
  ParmParse pp_pc(m_opt_string.c_str());
  pp_pc.query("solver_type", m_pc_type);

  if      (m_pc_type == "gmres" ) {
    m_pc_solver = new GMRESSolver<T,Ops>;
  } else if (m_pc_type == "jacobi") {
    m_pc_solver = new PointJacobiSolver<T,Ops>;
  } else if (m_pc_type == "twolevel") {
    m_pc_solver = new TwoLevelSolver<T,Ops>;
  } else if (m_pc_type == "gauss_seidel") {
    m_pc_solver = new GaussSeidelSolver<T,Ops>;
  } else {
     std::cout << "Error in BasicGKPreconditioner::define() - " 
               << m_pc_type << " "
               << "is not a valid type for preconditioning solver.\n";
    m_pc_solver = NULL;
  }
  CH_assert(m_pc_solver != NULL);

  m_pc_solver->define(m_pc_object, m_out_string, m_opt_string, m_Z);
  m_is_defined = true;
  return;
}

#include "NamespaceFooter.H"
#endif

