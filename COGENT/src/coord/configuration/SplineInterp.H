#ifndef _SPLINEINTERP_H_
#define _SPLINEINTERP_H_

#include "ParmParse.H"
#include "Interp.H"

#include "NamespaceHeader.H"

class SplineInterp
  : public Interp
{
  public:

    SplineInterp( const ParmParse& parm_parse,
                  const FArrayBox& coords,
                  const FArrayBox& data );

    ~SplineInterp();

    virtual double interpolate(const int       component,
                               const int       mode,
                               const RealVect& xi) const;

    virtual void interpolate(const int        src_component,
                             const int        dst_component,
                             const int        mode,
                             const FArrayBox& xi,
                             const Box&       box,
                             FArrayBox&       dXdXi) const;

  private:

    double **m_wk;
};

#include "NamespaceFooter.H"

#endif
