#include "BSplineInterp3D.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

#define max(a,b) ((a) >= (b) ? (a) : (b))

extern "C" {
   void db3int_(const double* x, const int& nx, const double* y, const int& ny, const double* z, const int& nz,
                const int& kx, const int& ky, const int& kz, double* tx, double* ty, double* tz,
                double* fcn, const int& ldf1, const int& ldf2, double* work, int& iflag);

   double db3val_(const double& x, const double& y, const double& z, const int& idx, const int& idy, const int& idz,
                  const double* tx, const double* ty, const double* tz, const int& nx, const int& ny, const int& nz,
                  const int& kx, const int& ky, const int& kz, const double *fcn, const int& ldf1, const int& ldf2,
                  int *icont, int *iwork, double* work, int& iflag);
}



BSplineInterp3D::BSplineInterp3D( const ParmParse&  a_pp,
                                  const FArrayBox&  a_nodes,
                                  const FArrayBox&  a_data )
   : Interp3D(a_nodes, a_data)
{
   if (SpaceDim != 3) {
      MayDay::Error("BSplineInterp3D is three-dimensional");
   }

   if (a_pp.contains("bspline_order")) {
      a_pp.get("bspline_order", m_bspline_order);
   }
   else {
      m_bspline_order = 6;
   }

   if (a_pp.contains("bspline_diagnostics")) {
      a_pp.get("bspline_diagnostics", m_print_diagnostics);
   }
   else {
      m_print_diagnostics = false;
   }

   m_x_knots = new double[dim(0) + m_bspline_order];
   m_y_knots = new double[dim(1) + m_bspline_order];
   m_z_knots = new double[dim(2) + m_bspline_order];

   int n_max = max( max(dim(0),dim(1)), dim(2) );
   int work_length = dim(0)*dim(1)*dim(2) + 2*m_bspline_order*(n_max+1);

   int num_components = nComp();
   m_wk          = new double*[num_components];
   m_db3val_work = new double*[num_components];
   m_icont       = new int[num_components];
   m_iwork       = new int*[num_components];

   for (int component=0; component<num_components; ++component) {

      m_wk[component] = new double[work_length];
      m_db3val_work[component] = new double[m_bspline_order*(m_bspline_order + 4) + 2];
      m_iwork[component] = new int[10];

      int iflag = 1; // Knot sequence chosen by db3int
      db3int_(m_x, dim(0), m_y, dim(1), m_z, dim(2),
              m_bspline_order, 2, m_bspline_order,
              m_x_knots, m_y_knots, m_z_knots,
              m_data.dataPtr(component), dim(0), dim(1),
              m_wk[component], iflag);
      if (iflag != 0) {
         cout << "db3int returned " << iflag << endl;
         exit(1);
      }

      m_icont[component] = 0; // First call
      db3val_(m_x[0], m_y[0], m_z[0], 0, 0, 0, m_x_knots, m_y_knots, m_z_knots,
              dim(0), dim(1), dim(2), m_bspline_order, 2, m_bspline_order,
              m_data.dataPtr(component), dim(0), dim(1), &m_icont[component],
              m_iwork[component], m_db3val_work[component], iflag );
      if (iflag != 0) {
         cout << "db3int returned " << iflag << endl;
         exit(1);
      }
   }
}


BSplineInterp3D::~BSplineInterp3D()
{
   for (int component=0; component<nComp(); ++component) {
      if (m_wk[component]) delete [] m_wk[component];
      if (m_db3val_work[component]) delete [] m_db3val_work[component];
      if (m_iwork[component]) delete [] m_iwork[component];
   }
   delete [] m_wk;
   delete [] m_db3val_work;
   delete [] m_icont;
   delete [] m_iwork;

   if (m_x_knots) delete [] m_x_knots;
   if (m_y_knots) delete [] m_y_knots;
   if (m_z_knots) delete [] m_z_knots;
}


double BSplineInterp3D::interpolate( const int        a_component,
                                     const int        a_mode,
                                     const RealVect&  a_xi ) const
{
   CH_assert(a_component >= 0 && a_component < nComp());
   double val;

   if ( m_print_diagnostics && procID() == 0 ) {
      if (a_xi[0] < m_x[0] || a_xi[0] > m_x[dim(0)-1] ) {
         cout << "BSplineInterp3D::interpolate(): input point out of range in direction 0" << endl;
         cout << "min, max = " << m_x[0] << "  " << m_x[dim(0)-1] << ", coord = " << a_xi[0] << endl;
      }

      if (a_xi[1] < m_y[0] || a_xi[1] > m_y[dim(1)-1] ) {
         cout << "BSplineInterp3D::interpolate(): input point out of range in direction 1" << endl;
         cout << "min, max = " << m_y[0] << "  " << m_y[dim(1)-1] << ", coord = " << a_xi[1] << endl;
      }

      if (a_xi[2] < m_z[0] || a_xi[2] > m_z[dim(2)-1] ) {
         cout << "BSplineInterp3D::interpolate(): input point out of range in direction 2" << endl;
         cout << "min, max = " << m_z[0] << "  " << m_z[dim(2)-1] << ", coord = " << a_xi[2] << endl;
      }
   }

   int idx, idy, idz;
   setIDs(a_mode, idx, idy, idz);

   int iflag;
   val = db3val_(a_xi[0], a_xi[1], a_xi[2], idx, idy, idz, m_x_knots, m_y_knots, m_z_knots,
                 dim(0), dim(1), dim(2), m_bspline_order, 2, m_bspline_order,
                 m_data.dataPtr(a_component), dim(0), dim(1),
                 &m_icont[a_component], m_iwork[a_component], m_db3val_work[a_component], iflag );
   if (iflag != 0) {
      cout << "db3val returned " << iflag << endl;
      exit(1);
   }

   return val;
}


void BSplineInterp3D::interpolate( const int         a_src_component,
                                   const int         a_dst_component,
                                   const int         a_mode,
                                   const FArrayBox&  a_xi,
                                   const Box&        a_box,
                                   FArrayBox&        a_dXdXi ) const
{
   CH_assert(a_src_component >= 0 && a_src_component < nComp());
   CH_assert(a_dst_component >= 0 && a_dst_component < a_dXdXi.nComp());

   int idx, idy, idz;
   setIDs(a_mode, idx, idy, idz);

   for (BoxIterator bit(a_box); bit.ok(); ++bit) {
      IntVect iv = bit();

      int iflag;
      a_dXdXi(iv,a_dst_component)
         = db3val_(a_xi(iv,0), a_xi(iv,1), a_xi(iv,2), idx, idy, idz, m_x_knots, m_y_knots, m_z_knots,
                   dim(0), dim(1), dim(2), m_bspline_order, m_bspline_order, m_bspline_order,
                   m_data.dataPtr(a_src_component), dim(0), dim(1), &m_icont[a_src_component],
                   m_iwork[a_src_component], m_db3val_work[a_src_component], iflag );
      if (iflag != 0) {
         cout << "db3val returned " << iflag << endl;
         exit(1);
      }
   }
}


void BSplineInterp3D::setIDs( const int a_mode,
                              int&      a_idx,
                              int&      a_idy,
                              int&      a_idz ) const
{
   a_idx = 0;
   a_idy = 0;
   a_idz = 0;
   switch ( a_mode )
      {
      case X_DERIVATIVE:
         a_idx = 1;
         break;
      case Y_DERIVATIVE:
         a_idy = 1;
         break;
      case Z_DERIVATIVE:
         a_idz = 1;
         break;
      case XX_SECOND_DERIVATIVE:
         a_idx = 2;
         break;
      case XY_SECOND_DERIVATIVE:
         a_idx = 1;
         a_idy = 1;
         break;
      case XZ_SECOND_DERIVATIVE:
         a_idx = 1;
         a_idz = 1;
         break;
      case YY_SECOND_DERIVATIVE:
         a_idy = 2;
         break;
      case YZ_SECOND_DERIVATIVE:
         a_idy = 1;
         a_idz = 1;
         break;
      case ZZ_SECOND_DERIVATIVE:
         a_idz = 2;
         break;
      default:
         break;
      }
}



#include "NamespaceFooter.H"


