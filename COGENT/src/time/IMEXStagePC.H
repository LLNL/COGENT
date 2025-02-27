#ifndef _IMEX_STAGE_PC_H_
#define _IMEX_STAGE_PC_H_

#include "LinearFunction.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class IMEXStagePCFunction : public LinearFunction<T,Ops>
{
  public:

    virtual void setBaseSolution(const T&) = 0;

    virtual void setBaseRHS(const T&, const T&)
    {
      MayDay::Error("IMEXStagePCFunction::setBaseRHS() called on base class!\n");
    }
    virtual void setBaseRHS(const T&)
    {
      MayDay::Error("IMEXStagePCFunction::setBaseRHS() called on base class!\n");
    }

    virtual void setJFNKEps(Real) = 0;
    virtual void setIsLinear(bool) = 0;

    virtual const BandedMatrix& getBandedMatrixForm() const = 0;
    virtual BandedMatrix& getBandedMatrixForm() = 0;

    virtual bool isDefined() const = 0;

  protected:

  private:
    
};

#include "NamespaceFooter.H"
#endif
