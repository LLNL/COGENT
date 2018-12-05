#ifndef _HERMITEINTERP_H_
#define _HERMITEINTERP_H_

#include "ParmParse.H"
#include "Interp.H"

#include "NamespaceHeader.H"

class HermiteInterp
  : public Interp
{
  public:

    HermiteInterp( const ParmParse& parm_parse,
                   const FArrayBox& coords,
                   const FArrayBox& data );

    virtual ~HermiteInterp();

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
