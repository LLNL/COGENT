#ifndef _Preconditioner_H_
#define _Preconditioner_H_

#include<string>

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class Preconditioner
{
  public:

    Preconditioner() { m_sys_id = -1; }
    virtual ~Preconditioner() {;}

    virtual void define(const T&, 
                        void*, 
                        const std::string&, 
                        const std::string&,
                        bool ) = 0;

    virtual void update(  const T&, 
                          const Real, 
                          const int, 
                          const int, 
                          const Real a_a = 0.0)
    {
      MayDay::Error("Preconditioner::update() - not implemented!");
    }

    virtual void apply(T&, const T&) = 0;
    virtual bool isDefined() const = 0;

    /* additional stuff for matrix-free types */
    virtual void setBaseSolution(const T& a_x) {}
    virtual void setBaseRHS(const T& a_F, const T& a_M) {}
    virtual void setBaseRHS(const T& a_M) {}
    virtual void setJFNKEps(Real a_e) {}
    virtual void setIsLinear(bool a_b) {}

    virtual const std::string typeName() { return "BaseClass"; }

    inline void setSysID(const int a_id) { m_sys_id = a_id; }
    inline int getSysID() const { return m_sys_id; }

  protected:

    int m_sys_id;

  private:

};

typedef std::vector<int> DOFList;

#include "NamespaceFooter.H"
#endif

