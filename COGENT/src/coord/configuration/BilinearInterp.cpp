#include "BilinearInterp.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

using namespace std;

#define max(a,b) ((a) >= (b) ? (a) : (b))

extern "C" {
  void db2int_(const double* x, const int& nx, const double* y, const int& ny,
               const int& kx, const int& ky, double* tx, double* ty, double* fcn,
               const int& ldf, double* work, int& iflag);

  double db2val_(const double& x, const double& y, const int& idx, const int& idy,
                 const double* tx, const double *ty, const int& nx, const int& ny,
                 const int& kx, const int& ky, const double *fcn, const int& ldf,
                 int *icont, int *iwork, double* work, int& iflag);
}


BilinearInterp::BilinearInterp( const ParmParse& a_parm_parse,
                                const FArrayBox& a_coords,
                                const FArrayBox& a_data )
   : Interp(a_coords, a_data)
{
  if (SpaceDim != 2) {
    MayDay::Error("BilinearInterp is only two-dimensional");
  }

  m_coords.define(a_coords.box(), a_coords.nComp());
  m_coords.copy(a_coords);
}



BilinearInterp::~BilinearInterp()
{
}


#if CFG_DIM ==3
double BilinearInterp::interpolate(const int       a_component,
                                   const int       a_mode,
                                   const RealVect& a_xi) const
{
    return 0.0;
}
#else
double BilinearInterp::interpolate(const int       a_component,
                                   const int       a_mode,
                                   const RealVect& a_xi) const
{
  CH_assert(a_component >= 0 && a_component < nComp());

  const Box& box = m_coords.box();

  IntVect lo = box.smallEnd();
  IntVect hi = box.bigEnd();

  double h0 = (m_coords(hi,0) - m_coords(lo,0)) / (box.size(0) - 1);
  double h1 = (m_coords(hi,1) - m_coords(lo,1)) / (box.size(1) - 1);

  double x = a_xi[0];
  double y = a_xi[1];

  int i;
  bool found_cell0 = false;
  for (int n=0; n<box.size(0)-1; ++n) {
     if ( m_x[n] <= x && x <m_x[n+1] ) {
        i = n;
        found_cell0 = true;
        break;
     }
  }
  i += lo[0];

  int j;
  bool found_cell1 = false;
  for (int n=0; n<box.size(1)-1; ++n) {
     if ( m_y[n] <= y && y <m_y[n+1] ) {
        j = n;
        found_cell1 = true;
        break;
     }
  }
  j += lo[1];

  CH_assert(found_cell0 && found_cell1);

  IntVect ll = IntVect(i,j);
  IntVect ul = ll; ul[1]++;
  IntVect lr = ll; lr[0]++;
  IntVect ur = ll + IntVect::Unit;

  double xlo = m_coords(ll,0);
  double ylo = m_coords(ll,1);
  double xhi = m_coords(ur,0);
  double yhi = m_coords(ur,1);

  CH_assert(x >= xlo && x <= xhi && y >= ylo && y <= yhi);

  double val;

  switch( a_mode ) {
  case 0:
     val = (x - xlo)*(y - ylo)*m_data(ur,a_component)
         + (x - xlo)*(yhi - y)*m_data(lr,a_component)
         + (xhi - x)*(y - ylo)*m_data(ul,a_component)
         + (xhi - x)*(yhi - y)*m_data(ll,a_component);
     break;
  case 1:
     val = (y - ylo)*m_data(ur,a_component)
         + (yhi - y)*m_data(lr,a_component)
         - (y - ylo)*m_data(ul,a_component)
         - (yhi - y)*m_data(ll,a_component);
     break;
  case 2:
     val = (x - xlo)*m_data(ur,a_component)
         - (x - xlo)*m_data(lr,a_component)
         + (xhi - x)*m_data(ul,a_component)
         - (xhi - x)*m_data(ll,a_component);
     break;
  default:
     MayDay::Error("BilinearInterp::interpolate(): Invalid a_mode argument");
  }

  return val / (h0*h1);
}
#endif


void BilinearInterp::interpolate(const int        a_src_component,
                                 const int        a_dst_component,
                                 const int        a_mode,
                                 const FArrayBox& a_xi,
                                 const Box&       a_box,
                                 FArrayBox&       a_dXdXi) const
{
  CH_assert(a_src_component >= 0 && a_src_component < nComp());
  CH_assert(a_dst_component >= 0 && a_dst_component < a_dXdXi.nComp());

  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    IntVect iv = bit();

    RealVect xi;
    for (int n=0; n<SpaceDim; ++n) {
       xi[n] = a_xi(iv,n);
    }

    a_dXdXi(iv,a_dst_component) = interpolate(a_src_component, a_mode, xi);
  }
}



#include "NamespaceFooter.H"


