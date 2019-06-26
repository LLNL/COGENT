#ifndef _INTERP_H_
#define _INTERP_H_

#include "RealVect.H"
#include "FArrayBox.H"

#include "NamespaceHeader.H"

class Interp
{
  public:

    Interp( const FArrayBox& coords,
            const FArrayBox& data );

    virtual ~Interp();

    virtual double interpolate(const int       component,
                               const int       mode,
                               const RealVect& xi) const = 0;

    virtual void interpolate(const int        src_component,
                             const int        dst_component,
                             const int        mode,
                             const FArrayBox& xi,
                             const Box&       box,
                             FArrayBox&       dXdXi) const = 0;



    int nComp() const {return m_data.nComp();}

    int dim(int dir) const {return m_data.box().size(dir);}

    enum interpolation_mode {FUNCTION, RADIAL_DERIVATIVE, POLOIDAL_DERIVATIVE, SECOND_RADIAL_DERIVATIVE,
                             RADIAL_POLOIDAL_DERIVATIVE, SECOND_POLOIDAL_DERIVATIVE};

    protected:

    double *m_x;
    double *m_y;

    FArrayBox m_data;
};

#include "NamespaceFooter.H"

#endif
