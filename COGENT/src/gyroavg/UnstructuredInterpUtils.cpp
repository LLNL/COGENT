#include "UnstructuredInterpUtils.H"
#include "LUP.H"

#include "NamespaceHeader.H"

static void createMatrix( const int a_npts,
                          Real***   a_Amat )
{
  (*a_Amat) = (Real**) calloc (a_npts, sizeof(Real*));
  for (int i=0; i<a_npts; i++) {
    (*a_Amat)[i] = (Real*) calloc (a_npts, sizeof(Real));
  }
}

static void setMatrixBicubic( const std::vector<Real>&  a_x,
                              const std::vector<Real>&  a_y,
                              const int                 a_npts,
                              Real**                    a_Amat )
{
  CH_assert(a_npts == 9);
  for (int i=0; i<a_npts; i++) {
    a_Amat[i][ 0] = 1.0;
    a_Amat[i][ 1] = a_x[i];
    a_Amat[i][ 2] = a_y[i];
    a_Amat[i][ 3] = a_x[i]*a_y[i];
    a_Amat[i][ 4] = a_x[i]*a_x[i];
    a_Amat[i][ 5] = a_y[i]*a_y[i];
    a_Amat[i][ 6] = a_x[i]*a_x[i]*a_y[i];
    a_Amat[i][ 7] = a_x[i]*a_y[i]*a_y[i];
    a_Amat[i][ 8] = a_x[i]*a_x[i]*a_y[i]*a_y[i];
  }
}

static void setMatrixBiquartic( const std::vector<Real>&  a_x,
                                const std::vector<Real>&  a_y,
                                const int                 a_npts,
                                Real**                    a_Amat )
{
  CH_assert(a_npts == 16);
  for (int i=0; i<a_npts; i++) {
    a_Amat[i][ 0] = 1.0;
    a_Amat[i][ 1] = a_x[i];
    a_Amat[i][ 2] = a_y[i];
    a_Amat[i][ 3] = a_x[i]*a_y[i];
    a_Amat[i][ 4] = a_x[i]*a_x[i];
    a_Amat[i][ 5] = a_y[i]*a_y[i];
    a_Amat[i][ 6] = a_x[i]*a_x[i]*a_y[i];
    a_Amat[i][ 7] = a_x[i]*a_y[i]*a_y[i];
    a_Amat[i][ 8] = a_x[i]*a_x[i]*a_y[i]*a_y[i];
    a_Amat[i][ 9] = a_x[i]*a_x[i]*a_x[i];
    a_Amat[i][10] = a_y[i]*a_y[i]*a_y[i];
    a_Amat[i][11] = a_x[i]*a_x[i]*a_x[i]*a_y[i];
    a_Amat[i][12] = a_x[i]*a_y[i]*a_y[i]*a_y[i];
    a_Amat[i][13] = a_x[i]*a_x[i]*a_y[i]*a_y[i]*a_y[i];
    a_Amat[i][14] = a_x[i]*a_x[i]*a_x[i]*a_y[i]*a_y[i];
    a_Amat[i][15] = a_x[i]*a_x[i]*a_x[i]*a_y[i]*a_y[i]*a_y[i];
  }
}

static void setMatrixTrilinear( const std::vector<Real>&  a_x,
                                const std::vector<Real>&  a_y,
                                const std::vector<Real>&  a_z,
                                const int                 a_npts,
                                Real**                    a_Amat )
{
  CH_assert(a_npts == 8);
  for (int i=0; i<a_npts; i++) {
    a_Amat[i][0] = 1.0;
    a_Amat[i][1] = a_x[i];
    a_Amat[i][2] = a_y[i];
    a_Amat[i][3] = a_z[i];
    a_Amat[i][4] = a_x[i]*a_y[i];
    a_Amat[i][5] = a_y[i]*a_z[i];
    a_Amat[i][6] = a_x[i]*a_z[i];
    a_Amat[i][7] = a_x[i]*a_y[i]*a_z[i];
  }
}

static void setMatrixTricubic(  const std::vector<Real>&  a_x,
                                const std::vector<Real>&  a_y,
                                const std::vector<Real>&  a_z,
                                const int                 a_npts,
                                Real**                    a_Amat )
{
  CH_assert(a_npts == 27);
  for (int i=0; i<a_npts; i++) {
    double xc = a_x[i];
    double yc = a_y[i];
    double zc = a_z[i];
    a_Amat[i][ 0] = 1.0;
    a_Amat[i][ 1] = xc;
    a_Amat[i][ 2] = yc;
    a_Amat[i][ 3] = zc;
    a_Amat[i][ 4] = xc*yc;
    a_Amat[i][ 5] = yc*zc;
    a_Amat[i][ 6] = xc*zc;
    a_Amat[i][ 7] = xc*yc*zc;
    a_Amat[i][ 8] = xc*xc;
    a_Amat[i][ 9] = yc*yc;
    a_Amat[i][10] = zc*zc;
    a_Amat[i][11] = xc*xc*yc;
    a_Amat[i][12] = xc*xc*zc;
    a_Amat[i][13] = yc*yc*xc;
    a_Amat[i][14] = yc*yc*zc;
    a_Amat[i][15] = zc*zc*xc;
    a_Amat[i][16] = zc*zc*yc;
    a_Amat[i][17] = xc*xc*yc*yc;
    a_Amat[i][18] = xc*xc*zc*zc;
    a_Amat[i][19] = yc*yc*zc*zc;
    a_Amat[i][20] = xc*xc*yc*zc;
    a_Amat[i][21] = xc*yc*yc*zc;
    a_Amat[i][22] = xc*yc*zc*zc;
    a_Amat[i][23] = xc*xc*yc*yc*zc;
    a_Amat[i][24] = xc*xc*yc*zc*zc;
    a_Amat[i][25] = xc*yc*yc*zc*zc;
    a_Amat[i][26] = xc*xc*yc*yc*zc*zc;
  }
}

static void destroyMatrix(  const int a_npts,
                            Real**    a_Amat )
{
  for (int i=0; i<a_npts; i++) free(a_Amat[i]);
  free(a_Amat);
}

static void rotatePoints2D( std::vector<Real>& a_x,
                            std::vector<Real>& a_y,
                            Real a_xc,
                            Real a_yc )
{
  /* Set the angles */
  Real sin_t = 0.1;
  Real cos_t = sqrt(1.0-sin_t*sin_t);

  int npts = a_x.size();
  for (int n=0; n<npts; n++) {
    Real x, y;
    x = a_x[n]*cos_t - a_y[n]*sin_t;
    y = a_x[n]*sin_t + a_y[n]*cos_t;
    a_x[n] = x;
    a_y[n] = y;
  }
  {
    Real x, y;
    x = a_xc*cos_t - a_yc*sin_t;
    y = a_xc*sin_t + a_yc*cos_t;
    a_xc = x;
    a_yc = y;
  }
}

static void rotatePoints3D( std::vector<Real>& a_x,
                            std::vector<Real>& a_y,
                            std::vector<Real>& a_z,
                            Real a_xc,
                            Real a_yc,
                            Real a_zc )
{
  /* Set the angles */
  Real sin_a = 0.1;
  Real sin_b = 0.1;
  Real sin_c = 0.1;
  Real cos_a = sqrt(1.0-sin_a*sin_a);
  Real cos_b = sqrt(1.0-sin_b*sin_b);
  Real cos_c = sqrt(1.0-sin_c*sin_c);

  int npts = a_x.size();
  for (int i=0; i<npts; i++) {
    Real x, y, z;
    x =   cos_a*cos_b*a_x[i] 
        + (cos_a*sin_b*sin_c - sin_a*cos_c)*a_y[i]
        + (cos_a*sin_b*cos_c + sin_a*sin_c)*a_z[i];
    y =   sin_a*cos_b*a_x[i] 
        + (sin_a*sin_b*sin_c + cos_a*cos_c)*a_y[i]
        + (sin_a*sin_b*cos_c - cos_a*sin_c)*a_z[i];
    z = - sin_b*a_x[i] 
        + cos_b*sin_c*a_y[i] 
        + cos_b*cos_c*a_z[i];
    a_x[i] = x;
    a_y[i] = y;
    a_z[i] = z;
  }
  {
    Real x, y, z;
    x =   cos_a*cos_b*a_xc 
        + (cos_a*sin_b*sin_c - sin_a*cos_c)*a_yc
        + (cos_a*sin_b*cos_c + sin_a*sin_c)*a_zc;
    y =   sin_a*cos_b*a_xc 
        + (sin_a*sin_b*sin_c + cos_a*cos_c)*a_yc
        + (sin_a*sin_b*cos_c - cos_a*sin_c)*a_zc;
    z = - sin_b*a_xc 
        + cos_b*sin_c*a_yc 
        + cos_b*cos_c*a_zc;
    a_xc = x;
    a_yc = y;
    a_zc = z;
  }
}

Real UnstructuredInterpUtils::determinant(  const std::vector<Real>& a_x,
                                            const std::vector<Real>& a_y )
{
  CH_assert(a_x.size() == a_y.size());

  Real retval = 0.0;
  if (a_x.size() == 4) {

    /* 2D */

    Real x0 = a_x[0];
    Real x1 = a_x[1];
    Real x2 = a_x[2];
    Real x3 = a_x[3];
    Real y0 = a_y[0];
    Real y1 = a_y[1];
    Real y2 = a_y[2];
    Real y3 = a_y[3];
  
    retval = (  x0*(x3*(y1 - y2)*(y0 - y3) - x2*(y0 - y2)*(y1 - y3) + x1*(y0 - y1)*(y2 - y3)) 
           + x1*x2*(y1 - y2)*(y0 - y3) - x1*x3*(y0 - y2)*(y1 - y3) + x2*x3*(y0 - y1)*(y2 - y3) );

  } else {

    MayDay::Error("UnstructuredInterpUtils::determinant() not yet implemented for the size of input vectors.");

  }

  return retval;
}


void UnstructuredInterpUtils::bilinearInterpStencil( Stencil&                          a_stencil,
                                                     const CFG::FArrayBox&             a_coords,
                                                     const std::vector<CFG::IntVect>&  a_idx_set,
                                                     const CFG::RealVect&              a_x,
                                                     const Real                        a_tol )
{
  a_stencil.clear();

  /* 
   * Bilinear interpolation using 4 arbitrary points in 2D space 
  */

  const int n_pts = 4;
  CH_assert(a_idx_set.size() == n_pts);

  /* extract the (x,y) coords of the interpolation point */
  Real  xc = a_x[0],
        yc = a_x[1];

  /* extract the (x,y) coords of the four points */
  std::vector<Real> x(n_pts), y(n_pts);
  for (int n=0; n<n_pts; n++) {
    x[n] = a_coords(a_idx_set[n],0);
    y[n] = a_coords(a_idx_set[n],1);
  }

  /* let's shift everything such that (x[0], y[0]) is the origin */
  {
    Real  x0 = x[0],
          y0 = y[0];
    xc -= x0;
    yc -= y0;
    for (int n=0; n<n_pts; n++) {
      x[n] -= x0;
      y[n] -= y0;
    }
  }

  /* let's find the max distance and scale by it so that everything is O(1) */
  {
    Real dmax = 0;
    for (int n=1; n<n_pts; n++) {
      Real d2 =   (x[n]-x[0])*(x[n]-x[0]) 
                + (y[n]-y[0])*(y[n]-y[0]);
      Real d = sqrt(d2);
      if (d > dmax) dmax = d;
    }
    CH_assert(dmax > 0);
    xc /= dmax;
    yc /= dmax;
    for (int n=0; n<n_pts; n++) {
      x[n] /= dmax;
      y[n] /= dmax;
    }
  }

  /* if the quadrilateral is "degenerate", let's rotate everything 
   * such that one side is aligned with the x-axis */
  Real den = determinant(x, y);
  while ((den*den) < a_tol) {
    rotatePoints2D(x, y, xc, yc);
    den = determinant(x, y);
  }

  /* bilinear interp coefficients - computed using Mathematica */
  Real c0, c1, c2, c3;
  {
    Real x0 = x[0];
    Real x1 = x[1];
    Real x2 = x[2];
    Real x3 = x[3];
    Real y0 = y[0];
    Real y1 = y[1];
    Real y2 = y[2];
    Real y3 = y[3];
  
    c0 = (x1*(-(xc*(y2 - y3)*(y1 - yc)) + x3*(y1 - y3)*(y2 - yc) - x2*(y1 - y2)*(y3 - yc)) - x2*x3*(y2 - y3)*(y1 - yc) + x2*xc*(y1 - y3)*(y2 - yc) - x3*xc*(y1 - y2)*(y3 - yc));
    c1 = (x0*(xc*(y2 - y3)*(y0 - yc) - x3*(y0 - y3)*(y2 - yc) + x2*(y0 - y2)*(y3 - yc)) + x2*x3*(y2 - y3)*(y0 - yc) - x2*xc*(y0 - y3)*(y2 - yc) + x3*xc*(y0 - y2)*(y3 - yc));
    c2 = (x0*(-(xc*(y1 - y3)*(y0 - yc)) + x3*(y0 - y3)*(y1 - yc) - x1*(y0 - y1)*(y3 - yc)) - x1*x3*(y1 - y3)*(y0 - yc) + x1*xc*(y0 - y3)*(y1 - yc) - x3*xc*(y0 - y1)*(y3 - yc));
    c3 = (x0*(xc*(y1 - y2)*(y0 - yc) - x2*(y0 - y2)*(y1 - yc) + x1*(y0 - y1)*(y2 - yc)) + x1*x2*(y1 - y2)*(y0 - yc) - x1*xc*(y0 - y2)*(y1 - yc) +  x2*xc*(y0 - y1)*(y2 - yc));
  
    c0 /= den;
    c1 /= den;
    c2 /= den;
    c3 /= den;
  }

  /* push these coefficients and coordinates to the stencil */
  a_stencil.push_back(StencilPoint(a_idx_set[0],c0));
  a_stencil.push_back(StencilPoint(a_idx_set[1],c1));
  a_stencil.push_back(StencilPoint(a_idx_set[2],c2));
  a_stencil.push_back(StencilPoint(a_idx_set[3],c3));

  /* done */
  return;
}

void UnstructuredInterpUtils::bicubicInterpStencil(Stencil&                          a_stencil,
                                                   const CFG::FArrayBox&             a_coords,
                                                   const std::vector<CFG::IntVect>&  a_idx_set,
                                                   const CFG::RealVect&              a_x,
                                                   const Real                        a_tol )
{
  a_stencil.clear();

  /* 
   * Bicubic interpolation using 9 arbitrary points in 2D space 
  */

  const int n_pts = 9;
  CH_assert(a_idx_set.size() == n_pts);

  /* extract the (x,y) coords of the interpolation point */
  Real  xc = a_x[0],
        yc = a_x[1];

  /* extract the (x,y) coords of the four points */
  std::vector<Real> x(n_pts), y(n_pts);
  for (int n=0; n<n_pts; n++) {
    x[n] = a_coords(a_idx_set[n],0);
    y[n] = a_coords(a_idx_set[n],1);
  }

  /* let's shift everything such that (x[0], y[0]) is the origin */
  {
    Real  x0 = x[0],
          y0 = y[0];
    xc -= x0;
    yc -= y0;
    for (int n=0; n<n_pts; n++) {
      x[n] -= x0;
      y[n] -= y0;
    }
  }

  /* let's find the max distance and scale by it so that everything is O(1) */
  {
    Real dmax = 0;
    for (int n=1; n<n_pts; n++) {
      Real d2 =   (x[n]-x[0])*(x[n]-x[0]) 
                + (y[n]-y[0])*(y[n]-y[0]);
      Real d = sqrt(d2);
      if (d > dmax) dmax = d;
    }
    CH_assert(dmax > 0);
    xc /= dmax;
    yc /= dmax;
    for (int n=0; n<n_pts; n++) {
      x[n] /= dmax;
      y[n] /= dmax;
    }
  }

  /* create LHS matrix */
  Real **Amat = NULL;
  createMatrix(n_pts, &Amat);
  setMatrixBicubic(x, y, n_pts, Amat);

  /* allocate memory for pivot array */
  int *pivot = (int*) calloc (n_pts+1, sizeof(int));

  /* decompose A */
  int attempt = 0;
  int status = 0;
  status = LUP::decompose(Amat, n_pts, 1e-14, pivot);

  while ((!status) && (attempt < 60)) {
    rotatePoints2D(x, y, xc, yc);
    setMatrixBicubic(x, y, n_pts, Amat);
    status = LUP::decompose(Amat, n_pts, 1e-14, pivot);
    attempt++;
  }

  if (!status) {
    /* matrix was singular */
    pout()    << "Error in UnstructuredInterpUtils::bicubicInterpStencil():\n"
              << "  The nine points provided are degenerate.\n";
    std::cout << "Rank " << procID() << ": "
              << "Error in UnstructuredInterpUtils::bicubicInterpStencil() - see pout.<n> for details.\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  /* allocate memory for Ainv */
  Real **Ainv = NULL;
  createMatrix(n_pts, &Ainv);

  /* compute inverse of A */
  LUP::invert(Amat, pivot, n_pts, Ainv);

  /* compute the coefficients */
  std::vector<Real> interp_coeffs(n_pts);
  for (int n=0; n<n_pts; n++) {
    /* evaluate the polynomial at xc, yc, zc */
    interp_coeffs[n] =      Ainv[ 0][n] * 1.0
                        +   Ainv[ 1][n] * xc
                        +   Ainv[ 2][n] * yc
                        +   Ainv[ 3][n] * xc*yc
                        +   Ainv[ 4][n] * xc*xc
                        +   Ainv[ 5][n] * yc*yc
                        +   Ainv[ 6][n] * xc*xc*yc
                        +   Ainv[ 7][n] * xc*yc*yc
                        +   Ainv[ 8][n] * xc*xc*yc*yc;
  }

  /* push these coefficients and coordinates to the stencil */
  for (int n=0; n<n_pts; n++) {
    a_stencil.push_back(StencilPoint(a_idx_set[n],interp_coeffs[n]));
  }

  /* done */
  destroyMatrix(n_pts, Amat);
  destroyMatrix(n_pts, Ainv);
  free(pivot);
  return;
}

void UnstructuredInterpUtils::biquarticInterpStencil(Stencil&                          a_stencil,
                                                     const CFG::FArrayBox&             a_coords,
                                                     const std::vector<CFG::IntVect>&  a_idx_set,
                                                     const CFG::RealVect&              a_x,
                                                     const Real                        a_tol )
{
  a_stencil.clear();

  /* 
   * Biquartic interpolation using 16 arbitrary points in 2D space 
  */

  const int n_pts = 16;
  CH_assert(a_idx_set.size() == n_pts);

  /* extract the (x,y) coords of the interpolation point */
  Real  xc = a_x[0],
        yc = a_x[1];

  /* extract the (x,y) coords of the four points */
  std::vector<Real> x(n_pts), y(n_pts);
  for (int n=0; n<n_pts; n++) {
    x[n] = a_coords(a_idx_set[n],0);
    y[n] = a_coords(a_idx_set[n],1);
  }

  /* let's shift everything such that (x[0], y[0]) is the origin */
  {
    Real  x0 = x[0],
          y0 = y[0];
    xc -= x0;
    yc -= y0;
    for (int n=0; n<n_pts; n++) {
      x[n] -= x0;
      y[n] -= y0;
    }
  }

  /* let's find the max distance and scale by it so that everything is O(1) */
  {
    Real dmax = 0;
    for (int n=1; n<n_pts; n++) {
      Real d2 =   (x[n]-x[0])*(x[n]-x[0]) 
                + (y[n]-y[0])*(y[n]-y[0]);
      Real d = sqrt(d2);
      if (d > dmax) dmax = d;
    }
    CH_assert(dmax > 0);
    xc /= dmax;
    yc /= dmax;
    for (int n=0; n<n_pts; n++) {
      x[n] /= dmax;
      y[n] /= dmax;
    }
  }

  /* create LHS matrix */
  Real **Amat = NULL;
  createMatrix(n_pts, &Amat);
  setMatrixBiquartic(x, y, n_pts, Amat);

  /* allocate memory for pivot array */
  int *pivot = (int*) calloc (n_pts+1, sizeof(int));

  /* decompose A */
  int attempt = 0;
  int status = 0;
  status = LUP::decompose(Amat, n_pts, 1e-14, pivot);

  while ((!status) && (attempt < 60)) {
    rotatePoints2D(x, y, xc, yc);
    setMatrixBiquartic(x, y, n_pts, Amat);
    status = LUP::decompose(Amat, n_pts, 1e-14, pivot);
    attempt++;
  }

  if (!status) {
    /* matrix was singular */
    pout()    << "Error in UnstructuredInterpUtils::biquarticInterpStencil():\n"
              << "  The four points provided are degenerate.\n";
    std::cout << "Rank " << procID() << ": "
              << "Error in UnstructuredInterpUtils::bilinearInterpStencil() - see pout.<n> for details.\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  /* allocate memory for Ainv */
  Real **Ainv = NULL;
  createMatrix(n_pts, &Ainv);

  /* compute inverse of A */
  LUP::invert(Amat, pivot, n_pts, Ainv);

  /* compute the coefficients */
  std::vector<Real> interp_coeffs(n_pts);
  for (int n=0; n<n_pts; n++) {
    /* evaluate the polynomial at xc, yc, zc */
    interp_coeffs[n] =      Ainv[ 0][n] * 1.0
                        +   Ainv[ 1][n] * xc
                        +   Ainv[ 2][n] * yc
                        +   Ainv[ 3][n] * xc*yc
                        +   Ainv[ 4][n] * xc*xc
                        +   Ainv[ 5][n] * yc*yc
                        +   Ainv[ 6][n] * xc*xc*yc
                        +   Ainv[ 7][n] * xc*yc*yc
                        +   Ainv[ 8][n] * xc*xc*yc*yc
                        +   Ainv[ 9][n] * xc*xc*xc
                        +   Ainv[10][n] * yc*yc*yc
                        +   Ainv[11][n] * xc*xc*xc*yc
                        +   Ainv[12][n] * xc*yc*yc*yc
                        +   Ainv[13][n] * xc*xc*yc*yc*yc
                        +   Ainv[14][n] * xc*xc*xc*yc*yc
                        +   Ainv[15][n] * xc*xc*xc*yc*yc*yc;
  }

  /* push these coefficients and coordinates to the stencil */
  for (int n=0; n<n_pts; n++) {
    a_stencil.push_back(StencilPoint(a_idx_set[n],interp_coeffs[n]));
  }

  /* done */
  destroyMatrix(n_pts, Amat);
  destroyMatrix(n_pts, Ainv);
  free(pivot);
  return;
}

void UnstructuredInterpUtils::trilinearInterpStencil( Stencil&                          a_stencil,
                                                      const CFG::FArrayBox&             a_coords,
                                                      const std::vector<CFG::IntVect>&  a_idx_set,
                                                      const CFG::RealVect&              a_x,
                                                      const Real                        a_tol )
{
  a_stencil.clear();

  /* 
   * Trilinear interpolation using 8 arbitrary points in 3D space 
  */

  const int n_pts = 8;
  CH_assert(a_idx_set.size() == n_pts);

  /* extract the (x,y) coords of the interpolation point */
  Real  xc = a_x[0],
        yc = a_x[1],
        zc = a_x[2];

  /* extract the (x,y) coords of the four points */
  std::vector<Real> x(n_pts), y(n_pts), z(n_pts);
  for (int n=0; n<n_pts; n++) {
    x[n] = a_coords(a_idx_set[n],0);
    y[n] = a_coords(a_idx_set[n],1);
    z[n] = a_coords(a_idx_set[n],2);
  }

  /* let's shift everything such that (x[0], y[0], z[0]) is the origin */
  {
    Real  x0 = x[0],
          y0 = y[0],
          z0 = z[0];
    xc -= x0;
    yc -= y0;
    zc -= z0;
    for (int n=0; n<n_pts; n++) {
      x[n] -= x0;
      y[n] -= y0;
      z[n] -= z0;
    }
  }
  /* let's find the max distance and scale by it so that everything is O(1) */
  {
    Real dmax = 0;
    for (int n=1; n<n_pts; n++) {
      Real d2 =   (x[n]-x[0])*(x[n]-x[0]) 
                + (y[n]-y[0])*(y[n]-y[0]) 
                + (z[n]-z[0])*(z[n]-z[0]);
      Real d = sqrt(d2);
      if (d > dmax) dmax = d;
    }
    CH_assert(dmax > 0);
    xc /= dmax;
    yc /= dmax;
    zc /= dmax;
    for (int n=0; n<n_pts; n++) {
      x[n] /= dmax;
      y[n] /= dmax;
      z[n] /= dmax;
    }
  }

  /* create LHS matrix and a copy */
  Real **Amat(NULL);
  createMatrix(n_pts, &Amat ); 
  setMatrixTrilinear(x, y, z, n_pts, Amat ); 

  /* allocate memory for pivot array */
  int *pivot = (int*) calloc (n_pts+1, sizeof(int));

  /* decompose A */
  int attempt = 0;
  int status = 0;
  status = LUP::decompose(Amat, n_pts, 1e-14, pivot);

  while ((!status) && (attempt < 60)) {
    rotatePoints3D(x, y, z, xc, yc, zc);
    setMatrixTrilinear(x, y, z, n_pts, Amat);
    status = LUP::decompose(Amat, n_pts, 1e-14, pivot);
    attempt++;
  }

  if (!status) {
    /* matrix was singular */
    pout()    << "Error in UnstructuredInterpUtils::trilinearInterpStencil():\n"
              << "  The eight points provided are degenerate. More info:-\n"
              << "  Coordinates of the eight points are:\n"
              << setiosflags(ios::showpoint)
              << setiosflags(ios::scientific)
              << setprecision(16)
              << "    " << x[0] << ", " << y[0] << ", " << z[0] << "\n"
              << "    " << x[1] << ", " << y[1] << ", " << z[1] << "\n"
              << "    " << x[2] << ", " << y[2] << ", " << z[2] << "\n"
              << "    " << x[3] << ", " << y[3] << ", " << z[3] << "\n"
              << "    " << x[4] << ", " << y[4] << ", " << z[4] << "\n"
              << "    " << x[5] << ", " << y[5] << ", " << z[5] << "\n"
              << "    " << x[6] << ", " << y[6] << ", " << z[6] << "\n"
              << "    " << x[7] << ", " << y[7] << ", " << z[7] << "\n"
              << "  Coordinates of interp point: " << xc << ", " << yc << ", " << zc
              << resetiosflags(ios::scientific)
              << "\n";
    std::cout << "Rank " << procID() << ": "
              << "Error in UnstructuredInterpUtils::trilinearInterpStencil() - see pout.<n> for details.\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  /* allocate memory for Ainv */
  Real **Ainv(NULL);
  createMatrix(n_pts, &Ainv);

  /* compute inverse of A */
  LUP::invert(Amat, pivot, n_pts, Ainv);

  /* compute the coefficients */
  std::vector<Real> interp_coeffs(n_pts);
  for (int n=0; n<n_pts; n++) {
    /* polynomial coeffs */
    Real  a = Ainv[0][n],
          b = Ainv[1][n],
          c = Ainv[2][n],
          d = Ainv[3][n],
          e = Ainv[4][n],
          f = Ainv[5][n],
          g = Ainv[6][n],
          h = Ainv[7][n];

    /* evaluate the polynomial at xc, yc, zc */
    interp_coeffs[n] =      a * 1.0
                        +   b * xc
                        +   c * yc
                        +   d * zc
                        +   e * xc*yc
                        +   f * yc*zc
                        +   g * xc*zc
                        +   h * xc*yc*zc;
  }

  /* push these coefficients and coordinates to the stencil */
  for (int n=0; n<n_pts; n++) {
    a_stencil.push_back(StencilPoint(a_idx_set[n],interp_coeffs[n]));
  }

  /* done */
  destroyMatrix(n_pts, Amat);
  destroyMatrix(n_pts, Ainv);
  free(pivot);
  return;
}

void UnstructuredInterpUtils::tricubicInterpStencil(Stencil&                          a_stencil,
                                                    const CFG::FArrayBox&             a_coords,
                                                    const std::vector<CFG::IntVect>&  a_idx_set,
                                                    const CFG::RealVect&              a_x,
                                                    const Real                        a_tol )
{
  a_stencil.clear();

  /* 
   * Tricubic interpolation using 8 arbitrary points in 3D space 
  */

  const int n_pts = 27;
  CH_assert(a_idx_set.size() == n_pts);

  /* extract the (x,y,z) coords of the interpolation point */
  Real  xc = a_x[0],
        yc = a_x[1],
        zc = a_x[2];

  /* extract the (x,y,z) coords of the stencil points */
  std::vector<Real> x(n_pts), y(n_pts), z(n_pts);
  for (int n=0; n<n_pts; n++) {
    x[n] = a_coords(a_idx_set[n],0);
    y[n] = a_coords(a_idx_set[n],1);
    z[n] = a_coords(a_idx_set[n],2);
  }

  /* let's shift everything such that (x[0], y[0], z[0]) is the origin */
  {
    Real  x0 = x[0],
          y0 = y[0],
          z0 = z[0];
    xc -= x0;
    yc -= y0;
    zc -= z0;
    for (int n=0; n<n_pts; n++) {
      x[n] -= x0;
      y[n] -= y0;
      z[n] -= z0;
    }
  }
  /* let's find the max distance and scale by it so that everything is O(1) */
  {
    Real dmax = 0;
    for (int n=1; n<n_pts; n++) {
      Real d2 =   (x[n]-x[0])*(x[n]-x[0]) 
                + (y[n]-y[0])*(y[n]-y[0]) 
                + (z[n]-z[0])*(z[n]-z[0]);
      Real d = sqrt(d2);
      if (d > dmax) dmax = d;
    }
    CH_assert(dmax > 0);
    xc /= dmax;
    yc /= dmax;
    zc /= dmax;
    for (int n=0; n<n_pts; n++) {
      x[n] /= dmax;
      y[n] /= dmax;
      z[n] /= dmax;
    }
  }

  /* create LHS matrix and a copy */
  Real **Amat(NULL);
  createMatrix(n_pts, &Amat ); 
  setMatrixTricubic(x, y, z, n_pts, Amat ); 

  /* allocate memory for pivot array */
  int *pivot = (int*) calloc (n_pts+1, sizeof(int));

  /* decompose A */
  int attempt = 0;
  int status = 0;
  status = LUP::decompose(Amat, n_pts, 1e-14, pivot);

  while ((!status) && (attempt < 60)) {
    rotatePoints3D(x, y, z, xc, yc, zc);
    setMatrixTricubic(x, y, z, n_pts, Amat);
    status = LUP::decompose(Amat, n_pts, 1e-14, pivot);
    attempt++;
  }

  if (!status) {
    /* matrix was singular */
    pout()    << "Error in UnstructuredInterpUtils::tricubicInterpStencil():\n"
              << "  The points provided are degenerate.\n";
    std::cout << "Rank " << procID() << ": "
              << "Error in UnstructuredInterpUtils::tricubicInterpStencil() - see pout.<n> for details.\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  /* allocate memory for Ainv */
  Real **Ainv(NULL);
  createMatrix(n_pts, &Ainv);

  /* compute inverse of A */
  LUP::invert(Amat, pivot, n_pts, Ainv);

  /* compute the coefficients */
  std::vector<Real> interp_coeffs(n_pts);
  for (int n=0; n<n_pts; n++) {
    /* evaluate the polynomial at xc, yc, zc */
    interp_coeffs[n] =      Ainv[ 0][n] * 1.0
                        +   Ainv[ 1][n] * xc
                        +   Ainv[ 2][n] * yc
                        +   Ainv[ 3][n] * zc
                        +   Ainv[ 4][n] * xc*yc
                        +   Ainv[ 5][n] * yc*zc
                        +   Ainv[ 6][n] * xc*zc
                        +   Ainv[ 7][n] * xc*yc*zc
                        +   Ainv[ 8][n] * xc*xc
                        +   Ainv[ 9][n] * yc*yc
                        +   Ainv[10][n] * zc*zc
                        +   Ainv[11][n] * xc*xc*yc
                        +   Ainv[12][n] * xc*xc*zc
                        +   Ainv[13][n] * yc*yc*xc
                        +   Ainv[14][n] * yc*yc*zc
                        +   Ainv[15][n] * zc*zc*xc
                        +   Ainv[16][n] * zc*zc*yc
                        +   Ainv[17][n] * xc*xc*yc*yc
                        +   Ainv[18][n] * xc*xc*zc*zc
                        +   Ainv[19][n] * yc*yc*zc*zc
                        +   Ainv[20][n] * xc*xc*yc*zc
                        +   Ainv[21][n] * xc*yc*yc*zc
                        +   Ainv[22][n] * xc*yc*zc*zc
                        +   Ainv[23][n] * xc*xc*yc*yc*zc
                        +   Ainv[24][n] * xc*xc*yc*zc*zc
                        +   Ainv[25][n] * xc*yc*yc*zc*zc
                        +   Ainv[26][n] * xc*xc*yc*yc*zc*zc;
  }

  /* push these coefficients and coordinates to the stencil */
  for (int n=0; n<n_pts; n++) {
    a_stencil.push_back(StencilPoint(a_idx_set[n],interp_coeffs[n]));
  }

  /* done */
  destroyMatrix(n_pts, Amat);
  destroyMatrix(n_pts, Ainv);
  free(pivot);
  return;
}

void StencilUtils::compress(Stencil& a_stencil) {

  std::map<CFG::IntVect, Real> tmp;

  for (int i = 0; i < a_stencil.size(); i++) {
    auto it = tmp.find(a_stencil[i].first);
    if (it == tmp.end()) {
      tmp[a_stencil[i].first] = a_stencil[i].second;
    } else {
      tmp[a_stencil[i].first] += a_stencil[i].second;
    }
  }

  a_stencil.clear();
  a_stencil.resize(tmp.size());
  int i = 0;
  for (auto it = tmp.cbegin(); it != tmp.cend(); ++it) {
    a_stencil[i].first = it->first;
    a_stencil[i].second = it->second;
    i++;
  }
  CH_assert(i == a_stencil.size());
}

#include "NamespaceFooter.H"
