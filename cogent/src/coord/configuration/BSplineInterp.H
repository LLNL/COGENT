#ifndef _BSPLINEINTERP_H_
#define _BSPLINEINTERP_H_

#include "ParmParse.H"
#include "Interp.H"

#include "NamespaceHeader.H"

class BSplineInterp
  : public Interp
{
  public:

    BSplineInterp( const ParmParse& parm_parse,
                   const FArrayBox& coords,
                   const FArrayBox& data );

    ~BSplineInterp();

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

    int m_bspline_order;
    bool m_print_diagnostics;

    double *m_radial_knots;
    double *m_poloidal_knots;
    double **m_wk;
    mutable double** m_db2val_work;
    mutable int *m_icont;
    mutable int **m_iwork;
};

#include "NamespaceFooter.H"

#endif
