/*
 * Copied from Chombo's LinearSolver.H, but added
 * an additional template parameter for operators.
*/

#ifndef _LINEARFUNCTION_H_
#define _LINEARFUNCTION_H_

#include <cmath>
#include "REAL.H"
#include "BandedMatrix.H"
#include "NamespaceHeader.H"

template <class T, class Ops = void>
class LinearFunction
{
  public:
  
    virtual ~LinearFunction() {}
    virtual void clear(T& a_lhs) {}
  
    virtual void define(const T&, Ops&, const std::string&, const std::string&)
    {
      MayDay::Error("Calling define on base class!");
    }
  
    virtual void residual(T& a_lhs, const T& a_phi, const T& a_rhs) = 0;
    virtual void preCond(T& a_cor, const T& a_residual) = 0;
    virtual void applyOp(T& a_lhs, const T& a_phi) = 0;
    virtual void create(T& a_lhs, const T& a_rhs) = 0;
    virtual void assign(T& a_lhs, const T& a_rhs) = 0;
    virtual Real dotProduct(T& a_1, const T& a_2) = 0;
    virtual void incr(T& a_lhs, const T& a_x, Real a_scale) = 0;
    virtual void axby(T& a_lhs, const T& a_x, const T& a_y, Real a_a, Real a_b) = 0;
    virtual void scale(T& a_lhs, const Real& a_scale)  = 0;
    virtual Real norm(T& a_rhs, int a_ord) = 0;
    virtual void setToZero(T& a_lhs) = 0;
    
    virtual void assignLocal(T& a_lhs, const T& a_rhs)
    {
      this->assign(a_lhs, a_rhs);
    }
    virtual void mDotProduct(T& a_1, const int a_sz, const T a_2[], Real a_mdots[])
    {
      for (int j=0; j<a_sz; j++)
        {
          a_mdots[j] = dotProduct(a_1, a_2[j]);
        }
    }
    virtual Real dx() const
    {
      MayDay::Warning(" calling dx on base class\n");
      return 0.;
    }
    virtual void write(const T* a, const char* filename)
    {
      MayDay::Warning("LinearFunction::write not implemented");
    }
  
    virtual void getDiagonal(std::vector<Real>& a_x) const
    {
      MayDay::Warning("LinearFunction::getDiagonal() not implemented");
      a_x.resize(0);
    }
  
    virtual const BandedMatrix& getBandedMatrixForm() const
    {
      MayDay::Warning("LinearFunction::getMatrix() not implemented. The code should crash now!");
      return(m_bandedmatrix_dummy);
    }

  private:
    
    BandedMatrix m_bandedmatrix_dummy;

};

#include "NamespaceFooter.H"

#endif /*_LINEARFUNCTION_H_*/
