#ifndef _Function_H_
#define _Function_H_

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class Function
{
  public:
    virtual ~Function() {;}
    virtual void evalFunction(T&, const T&, const T&, int) = 0;
    virtual void define(const T&, Ops&) = 0;

    virtual bool checkLinearity(const T&)
    {
      return false;
    }

    virtual Real computeStepNorm(const T&) const = 0;
};

#include "NamespaceFooter.H"
#endif
