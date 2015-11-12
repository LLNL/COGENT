#include "SplineInterp.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

using namespace std;

extern "C" {
  void compton_(const int& md, const int& nxd, const int& nyd,
                const double *xd, const double *yd, const double *zd,
                const int& nip, const double* xi, const double* yi,
                double* zi, int& ier, const int& eval_mode, double* wk);
}

SplineInterp::SplineInterp( const ParmParse& a_parm_parse,
                            const FArrayBox& a_coords,
                            const FArrayBox& a_data )
   : Interp(a_coords, a_data)
{
  if (SpaceDim != 2) {
    MayDay::Error("SplineInterp is only two-dimensional");
  }

  int num_components = nComp();
  m_wk = new double*[num_components];
  int work_length = 3*dim(0)*dim(1);

  for (int component=0; component<num_components; ++component) {
    m_wk[component] = new double[work_length];

    double dummy;
    int ier;

    compton_(1, dim(0), dim(1), m_x, m_y, m_data.dataPtr(component), 1,
             &m_x[0], &m_y[0], &dummy, ier, 0, m_wk[component]);
    if ( ier != 0 ) {
      cout << "ier = " << ier << endl;
      MayDay::Error("SplineInterp(): RGBI3P returned a nonzero value");
    }
  }
}



SplineInterp::~SplineInterp()
{
  for (int component=0; component<nComp(); ++component) {
    if (m_wk[component]) delete [] m_wk[component];
  }
  delete [] m_wk;
}



double SplineInterp::interpolate(const int       a_component,
                                 const int       a_mode,
                                 const RealVect& a_xi) const
{
  CH_assert(a_component >= 0 && a_component < nComp());
  double val;
  int ier;

  compton_(2, dim(0), dim(1), m_x, m_y, m_data.dataPtr(a_component), 1,
           &a_xi[0], &a_xi[1], &val, ier, a_mode, m_wk[a_component]);
  if ( ier != 0 ) {
    cout << "ier = " << ier << endl;
    MayDay::Error("SplineInterp::interpolate(): compton returned a nonzero value");
  }

  return val;
}



void SplineInterp::interpolate(const int        a_src_component,
                               const int        a_dst_component,
                               const int        a_mode,
                               const FArrayBox& a_xi,
                               const Box&       a_box,
                               FArrayBox&       a_dXdXi) const
{
  CH_assert(a_src_component >= 0 && a_src_component < nComp());
  CH_assert(a_dst_component >= 0 && a_dst_component < a_dXdXi.nComp());
  int num_points = a_box.numPts();

  double *buffer = new double[3*num_points];
  double *x = buffer;
  double *y = &buffer[num_points];
  double *values = &buffer[2*num_points];

  BoxIterator bit(a_box);
  int k = 0;
  for (bit.begin(); bit.ok(); ++bit) {
    IntVect iv = bit();

    x[k] = a_xi(iv,0);
    y[k] = a_xi(iv,1);

    k++;
  }

  int ier;
  compton_(2, dim(0), dim(1), m_x, m_y, m_data.dataPtr(a_src_component), num_points,
           x, y, values, ier, a_mode, m_wk[a_src_component]);

  if ( ier != 0 ) {
    cout << "ier = " << ier << endl;
    MayDay::Error("SplineInterp::interpolate(): compton_ returned a nonzero value");
  }

  k = 0;
  for (bit.begin(); bit.ok(); ++bit) {
    a_dXdXi(bit(),a_dst_component) = values[k++];
  }

  delete [] buffer;
}



#include "NamespaceFooter.H"


