#ifndef _BSPLINEINTERP3D_H_
#define _BSPLINEINTERP3D_H_

#include "ParmParse.H"
#include "Interp3D.H"


#include "NamespaceHeader.H"

class BSplineInterp3D
   : public Interp3D
{
  public:

    BSplineInterp3D( const ParmParse&  pp,
                     const FArrayBox&  nodes,
                     const FArrayBox&  data );

    ~BSplineInterp3D();

    virtual double interpolate( const int        component,
                                const int        mode,
                                const RealVect&  xi ) const;

    virtual void interpolate( const int         src_component,
                              const int         dst_component,
                              const int         mode,
                              const FArrayBox&  xi,
                              const Box&        box,
                              FArrayBox&        dXdXi ) const;

  private:

    void setIDs(const int mode, int& idx, int& idy, int& idz) const;
   
    int m_bspline_order;
    bool m_print_diagnostics;

    double *m_x_knots;
    double *m_y_knots;
    double *m_z_knots;
    double **m_wk;
    mutable double** m_db3val_work;
    mutable int *m_icont;
    mutable int **m_iwork;
};


#include "NamespaceFooter.H"

#endif
