#ifndef _BILINEARINTERP_H_
#define _BILINEARINTERP_H_

#include "ParmParse.H"
#include "Interp.H"

#include "NamespaceHeader.H"

class BilinearInterp
  : public Interp
{
  public:

    BilinearInterp( const ParmParse& parm_parse,
                  const FArrayBox& coords,
                  const FArrayBox& data );

    ~BilinearInterp();

    virtual double interpolate(const int       component,
                               const int       mode,
                               const RealVect& xi) const;

    virtual void interpolate(const int        a_src_component,
                             const int        a_dst_component,
                             const int        a_mode,
                             const FArrayBox& a_xi,
                             const Box&       a_box,
                             FArrayBox&       a_dXdXi) const;

  private:

    FArrayBox m_coords;

    bool m_print_diagnostics;

};

#include "NamespaceFooter.H"

#endif
