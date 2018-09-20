#include "BSplineInterp.H"
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


BSplineInterp::BSplineInterp( const ParmParse& a_parm_parse,
                              const int        a_order,
                              const FArrayBox& a_coords,
                              const FArrayBox& a_data )
   : Interp(a_coords, a_data)
{
  if (SpaceDim != 2) {
    MayDay::Error("BSplineInterp is only two-dimensional");
  }

#if 0
  if (a_parm_parse.contains("bspline_order")) {
    a_parm_parse.get("bspline_order", m_bspline_order);
  }
  else {
     m_bspline_order = 6;
  }
#else
  m_bspline_order = a_order;
#endif

  if (a_parm_parse.contains("bspline_diagnostics")) {
    a_parm_parse.get("bspline_diagnostics", m_print_diagnostics);
  }
  else {
    m_print_diagnostics = false;
  }

  m_radial_knots   = new double[dim(0) + m_bspline_order];
  m_poloidal_knots = new double[dim(1) + m_bspline_order];

  int work_length = dim(0)*dim(1) + 2*m_bspline_order*(max(dim(0),dim(1))+1);

  int num_components = nComp();
  m_wk          = new double*[num_components];
  m_db2val_work = new double*[num_components];
  m_icont       = new int[num_components];
  m_iwork       = new int*[num_components];

  for (int component=0; component<num_components; ++component) {

    m_wk[component] = new double[work_length];
    m_db2val_work[component] = new double[4*m_bspline_order+1];
    m_iwork[component] = new int[6];

    int iflag = 1;
    db2int_(m_x, dim(0), m_y, dim(1), m_bspline_order, m_bspline_order,
            m_radial_knots, m_poloidal_knots, m_data.dataPtr(component),
            dim(0), m_wk[component], iflag);
    if (iflag != 0) {
      cout << "db2int returned " << iflag << endl;
      exit(1);
    }

    m_icont[component] = 0;
    double dummy = db2val_(m_x[0], m_y[0], 0, 0, m_radial_knots, m_poloidal_knots,
                           dim(0), dim(1), m_bspline_order, m_bspline_order,
                           m_data.dataPtr(component), dim(0), &m_icont[component],
                           m_iwork[component], m_db2val_work[component], iflag );
    if (iflag != 0) {
      cout << "db2int returned " << iflag << endl;
      exit(1);
    }
  }
}



BSplineInterp::~BSplineInterp()
{
  for (int component=0; component<nComp(); ++component) {
    if (m_wk[component]) delete [] m_wk[component];
    if (m_db2val_work[component]) delete [] m_db2val_work[component];
    if (m_iwork[component]) delete [] m_iwork[component];
  }
  delete [] m_wk;
  delete [] m_db2val_work;
  delete [] m_icont;
  delete [] m_iwork;

  if (m_radial_knots) delete [] m_radial_knots;
  if (m_poloidal_knots) delete [] m_poloidal_knots;
}



double BSplineInterp::interpolate(const int       a_component,
                                  const int       a_mode,
                                  const RealVect& a_xi) const
{
  CH_assert(a_component >= 0 && a_component < nComp());
  double val;

  if ( m_print_diagnostics && procID() == 0 ) {
    if (a_xi[0] < m_x[0] || a_xi[0] > m_x[dim(0)-1] ) {
      cout << "BSplineInterp::interpolate(): input point out of range in direction 0" << endl;
      cout << "min, max = " << m_x[0] << "  " << m_x[dim(0)-1] << ", coord = " << a_xi[0] << endl;
    }

    if (a_xi[1] < m_y[0] || a_xi[1] > m_y[dim(1)-1] ) {
      cout << "BSplineInterp::interpolate(): input point out of range in direction 1" << endl;
      cout << "min, max = " << m_y[0] << "  " << m_y[dim(1)-1] << ", coord = " << a_xi[1] << endl;
    }
  }

  int idx = 0;
  int idy = 0;
  if (a_mode == 1) {
    idx = 1;
  }
  else if (a_mode == 2) {
    idy = 1;
  }

  int iflag;
  val = db2val_(a_xi[0], a_xi[1], idx, idy, m_radial_knots, m_poloidal_knots,
                dim(0), dim(1), m_bspline_order, m_bspline_order,
                m_data.dataPtr(a_component), dim(0),
                &m_icont[a_component], m_iwork[a_component], m_db2val_work[a_component], iflag );
  if (iflag != 0) {
    cout << "db2val returned " << iflag << endl;
    exit(1);
  }

  return val;
}



void BSplineInterp::interpolate(const int        a_src_component,
                                const int        a_dst_component,
                                const int        a_mode,
                                const FArrayBox& a_xi,
                                const Box&       a_box,
                                FArrayBox&       a_dXdXi) const
{
  CH_assert(a_src_component >= 0 && a_src_component < nComp());
  CH_assert(a_dst_component >= 0 && a_dst_component < a_dXdXi.nComp());

  int idx = 0;
  int idy = 0;
  if (a_mode == 1) {
    idx = 1;
  }
  else if (a_mode == 2) {
    idy = 1;
  }

  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit) {
    IntVect iv = bit();

    int iflag;
    a_dXdXi(iv,a_dst_component)
      = db2val_(a_xi(iv,0), a_xi(iv,1), idx, idy, m_radial_knots, m_poloidal_knots,
                dim(0), dim(1), m_bspline_order, m_bspline_order,
                m_data.dataPtr(a_src_component), dim(0), &m_icont[a_src_component],
                m_iwork[a_src_component], m_db2val_work[a_src_component], iflag );
    if (iflag != 0) {
      cout << "db2val returned " << iflag << endl;
      exit(1);
    }
  }
}



#include "NamespaceFooter.H"


