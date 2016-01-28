#ifndef _Function_H_
#define _Function_H_

#include "NamespaceHeader.H"

template <class T, class Ops = void>
class Function
{
  public:
    virtual void evalFunction(T&, T&, T&) = 0;
    virtual void define(T&,Ops&) = 0;
};

#include "NamespaceFooter.H"
#endif
