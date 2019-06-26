#ifndef _INTERP3D_H_
#define _INTERP3D_H_

#include "RealVect.H"
#include "FArrayBox.H"

#include "NamespaceHeader.H"


class Interp3D
{
  public:

    Interp3D( const FArrayBox&  nodes,
              const FArrayBox&  data );

    virtual ~Interp3D();

    virtual double interpolate( const int        component,
                                const int        mode,
                                const RealVect&  xi ) const = 0;

    virtual void interpolate( const int         src_component,
                              const int         dst_component,
                              const int         mode,
                              const FArrayBox&  xi,
                              const Box&        box,
                              FArrayBox&        dXdXi ) const = 0;

    int nComp() const {return m_data.nComp();}

    int dim(int dir) const {return m_data.box().size(dir);}

    enum interpolation_mode {FUNCTION, X_DERIVATIVE, Y_DERIVATIVE, Z_DERIVATIVE, XX_SECOND_DERIVATIVE,
                             XY_SECOND_DERIVATIVE, XZ_SECOND_DERIVATIVE, YY_SECOND_DERIVATIVE,
                             YZ_SECOND_DERIVATIVE, ZZ_SECOND_DERIVATIVE};

  protected:

    double *m_x;
    double *m_y;
    double *m_z;

    FArrayBox m_data;
};


#include "NamespaceFooter.H"

#endif
