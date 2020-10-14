#ifndef _LINEARFUNCTIONSOLVER_H_
#define _LINEARFUNCTIONSOLVER_H_

#include "LinearFunction.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class LinearFunctionSolver
{
  public:

    virtual ~LinearFunctionSolver() {}

    virtual void define(LinearFunction<T,Ops>* const, const std::string&, const std::string& ,const T&) = 0;
    virtual void solve(T& a_x, const T& a_r, int a_its = -1) = 0;
    virtual void setTolerances(Real, Real) = 0;
    virtual int getCount() const {return 0;}
};

#include "NamespaceFooter.H"
#endif
