#ifndef _Function_H_
#define _Function_H_

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class Function
{
  public:
    virtual ~Function() {;}
    virtual void evalFunction(T&, const T&, const T&, int, bool) = 0;
    virtual void define(const T&, Ops&) = 0;
};

#include "NamespaceFooter.H"
#endif
