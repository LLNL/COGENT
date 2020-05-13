#ifndef _RADIALUTILS_H_
#define _RADIALUTILS_H_

#include <tuple>
#include <string>
#include <list>
#include <vector>
#include <iostream>
#include "NamespaceHeader.H"


class HermiteSpline
{
  /** The class contains three arrays of the coefficients for cubic interpolation by Hermite splines.
   *  There are a_N points given: x[i] and y[i], including end-points. Therefore there are a_N-1 domains.
   *  In each domain the function is approximated by a third order polynomial
   *  y = y[i] + a[i]*(x-x[i]) + b[i]*(x-x[i])**2 + c[i]*(x-x[i])**3.
   *  The approximation, its derivative, and second derivative are continuous everywhere.
   *  Hermite spline uses first derivatives as boundary conditions: m_a_0 and m_a_N,
   *  instead zero second derivatives, as it is realized in a conventional cubic spline.
   *  If not provided, the derivatives are obtained from cubic interpolation at the boundaries,
   *  therefore there should be at least 4 points a_N>=4. Otherwise, for the boundary conditions provied,  a_N>=2.
   *  The constructor is aborted for a_N=1, and the spline is not created.
   *  The arrays of the coefficients a[i], b[i], and c[i] are created and computed at by the constructor of the class.
   */
public:
  /// data_type are currently either double*, std::vector<double>&,
  /// std::list<std::pair<double,double>>&, std::vector<std::pair<double,double>>&  -- NO
  HermiteSpline();
  ///          (data_x, data_y, a_N, derivatives, deriv_0, deriv_N,  tolerance)
  HermiteSpline(double*, double*, int, bool=false, double=0, double=0, double = 1.0e-10);
  bool buildSpline(double*, double*, int, bool=false, double=0, double=0, double = 1.0e-10);
  /// Constructors for std::list<std::pair<double,double>>& and std::vector<std::pair<double,double>>&
  
  ~HermiteSpline();
  
  /// Return functions
  double getValAt(int a_ind, double a_x) const;
  void getValAt(double* a_res, double a_x) const;

  int m_N;                  /// Number of grid points
  bool m_initialized;       /// Flag if the spline was initialized
  double m_ERR;             /// Accuracy (tolerance)
  double m_a_0, m_a_N;      /// Derivatives at x[0] and x[N]

  /// Arrays for coefficients [N]
  double* m_a;
  double* m_b;
  double* m_c;

  double* m_h;              /// Array for grid spacing [m_N-1]
  double* m_arr_x;          /// Array for the grid [m_N]
  double* m_arr_y;          /// Array for the values on the grid [m_N]
};

double cubicFit(double* a_p_abc, double a_x_f, double* a_x, double* a_y, double a_ERR = 1.0e-10);

bool threeDiagSolver(double* a_x, double* a_a, double* a_b, double* a_c, double* a_f, int a_N, double a_ERR = 1.0e-10);

bool fillDataBRadial(HermiteSpline& a_spline, std::string& a_str);

bool readDataBRadial(std::list< std::pair<double, double> >& a_lst,  std::string& a_str);

#include "NamespaceFooter.H"

#endif
