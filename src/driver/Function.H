#ifndef _Function_H_
#define _Function_H_

#include "NamespaceHeader.H"

template <class Solution, class RHS, class Ops>
class Function
{
  public:
    virtual void evalFunction(RHS&, Solution&, RHS&) = 0;
};

#include "NamespaceFooter.H"
#endif
