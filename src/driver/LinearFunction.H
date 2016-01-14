/*
 * Copied from Chombo's LinearSolver.H, but added
 * an additional template parameter for operators.
*/

#ifndef _LINEARFUNCTION_H_
#define _LINEARFUNCTION_H_

#include "REAL.H"
#include "Box.H"
#include <cmath>
#include "NamespaceHeader.H"

/**
   Operator class for Linear solver for solving L(phi) = rhs.
 */
template <class T, class Ops>
class LinearFunction
{
public:
  virtual ~LinearFunction()
  {
  }

  ///
  /**
     Say you are  solving L(phi) = rhs.   Make a_lhs = L(a_phi) - a_rhs.   If a_homogeneous is true,
     evaluate the operator using homogeneous boundary conditions.
   */
  virtual void residual(  T& a_lhs, const T& a_phi, const T& a_rhs, bool a_homogeneous = false) = 0;

  ///
  /**
     Given the current state of the residual the correction, apply your preconditioner to a_cor.
   */
  virtual void preCond(   T& a_cor, const T& a_residual)       = 0;

  ///
  /**
     In the context of solving L(phi) = rhs, set a_lhs = L(a_phi).  If a_homogeneous is true,
     evaluate the operator using homogeneous boundary conditions.
   */
  virtual void applyOp(   T& a_lhs, const T& a_phi, bool a_homogeneous = false) = 0;

  ///
  /**
     Creat data holder a_lhs that mirrors a_rhs.   You do not need to copy the data of a_rhs,
     just  make a holder the same size.
   */
  virtual void create(    T& a_lhs, const T& a_rhs) = 0;

  ///
  /**
     Opposite of create -- perform any operations required before lhs goes
     out of scope. In general, the only time this needs to be defined in
     a derived class is if the create() function called new. Otherwise, the
     default no-op function is sufficient.
  */
  virtual void clear(T& a_lhs)
  {
  }

  ///
  /**
     Set a_lhs  equal to a_rhs.
   */
  virtual void assign(    T& a_lhs, const T& a_rhs)       = 0;

  virtual void assignLocal(T& a_lhs, const T& a_rhs)
  {
    this->assign(a_lhs, a_rhs);
  }
  ///
  /**
     Compute and return the dot product of a_1 and a_2.   In most contexts, this
     means return the sum over all data points of a_1*a_2.
   */
  virtual Real dotProduct(T& a_1, T& a_2)     = 0;
  /* multiple dot products (for GMRES) */
  virtual void mDotProduct(T& a_1, const int a_sz, T a_2[], Real a_mdots[])
  {
    for (int j=0; j<a_sz; j++)
      {
        a_mdots[j] = dotProduct(a_1, a_2[j]);
      }
  }

  ///
  /**
     Increment by scaled amount (a_lhs += a_scale*a_x).
   */
  virtual void incr  (    T& a_lhs, const T& a_x, Real a_scale) = 0;

  ///
  /**
     Set input to a scaled sum (a_lhs = a_a*a_x + a_b*a_y).
   */
  virtual void axby(      T& a_lhs, const T& a_x, const T& a_y, Real a_a, Real a_b) = 0;

  ///
  /**
     Multiply the input by a given scale (a_lhs *= a_scale).
   */
  virtual void scale(     T& a_lhs, const Real& a_scale)  = 0;

  ///
  /**
     Return the norm of  a_rhs.
     a_ord == 0  max norm, a_ord == 1 sum(abs(a_rhs)), else, L(a_ord) norm.
   */
  virtual Real norm(T& a_rhs, int a_ord) = 0;

  ///
  /**
     Return dx at this level of refinement
  */
  virtual Real dx() const
  {
    MayDay::Warning(" calling dx on base class\n");
    return 0.;
  }
  ///
  /**
     Set a_lhs to zero.
   */
  virtual void setToZero(T& a_lhs) = 0;

  ///
  /**
     Debugging aid for solvers.  Print out a "T" to a file named "filename"
     default implementation is to print out a message saying "LinearFunction::write not implemented"
  */
  virtual void write(const T* a, const char* filename)
  {
    MayDay::Warning("LinearFunction::write not implemented");
  }

};

#include "NamespaceFooter.H"

#endif /*_LINEARFUNCTION_H_*/
