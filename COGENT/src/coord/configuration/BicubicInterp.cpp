#include "BicubicInterp.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"


static double cubic( double f[4],
                     bool   derivative,
                     double t )
{
   /*
     Evaluates the cubic determined by the four input "f" values at unit-spaced
     coordinates x0 through x3 a point x located between x1 and x2.  t is the
     distance between the point and x1 and must be contained in the unit interval:

         f[0]        f[1]        f[2]        f[3]


          x0          x1     x    x2          x3

                    -->|  t  |<--

     If derivative is true, the value of the first derivative is returned instead.

     (https://en.wikipedia.org/wiki/Bicubic_interpolation)
   */
   
   CH_assert(t >= 0. && t <= 1.);

   double a =           2.*f[1];
   double b =   -f[0]           +    f[2];
   double c = 2.*f[0] - 5.*f[1] + 4.*f[2] - f[3];
   double d =   -f[0] + 3.*f[1] - 3.*f[2] + f[3];

   double value;

   if (derivative) {
      value = 0.5 * (b + t*(2.*c + t*3.*d)); 
   }
   else {
      value = 0.5 * (a + t*(b + t*(c + t*d))); 
   }

   return value;
}


BicubicInterp::BicubicInterp( const ParmParse& a_parm_parse,
                                const FArrayBox& a_coords,
                                const FArrayBox& a_data )
   : Interp(a_coords, a_data)
{
  if (SpaceDim != 2) {
    MayDay::Error("BicubicInterp is only two-dimensional");
  }

  m_coords.define(a_coords.box(), a_coords.nComp());
  m_coords.copy(a_coords);
}



BicubicInterp::~BicubicInterp()
{
}


#if CFG_DIM ==3
double BicubicInterp::interpolate(const int       a_component,
                                  const int       a_mode,
                                  const RealVect& a_xi) const
{
    return 0.0;
}
#else
double BicubicInterp::interpolate(const int       a_component,
                                  const int       a_mode,
                                  const RealVect& a_xi) const
{
  CH_assert(a_component >= 0 && a_component < nComp());
  CH_assert(a_mode >= 0 && a_mode <= 2);

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

  double tx = (x - m_coords(IntVect(i,j),0)) / h0;

  double b[4];
  for (int m=0; m<4; ++m) {
     double a[4];

     IntVect iv1 = IntVect(i,j+m-1);
     IntVect iv0 = iv1; iv0[0]--;
     IntVect iv2 = iv1; iv2[0]++;
     IntVect iv3 = iv2; iv3[0]++;

     a[0] = m_data(iv0,a_component);
     a[1] = m_data(iv1,a_component);
     a[2] = m_data(iv2,a_component);
     a[3] = m_data(iv3,a_component);

     b[m] = cubic(a, a_mode == 1, tx);
  }

  double ty = (y - m_coords(IntVect(i,j),1)) / h1;

  return cubic(b, a_mode == 2, ty);
}
#endif


void BicubicInterp::interpolate(const int        a_src_component,
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


