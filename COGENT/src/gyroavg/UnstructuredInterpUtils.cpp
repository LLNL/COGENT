#include "UnstructuredInterpUtils.H"
#include "LUP.H"

#include "NamespaceHeader.H"


void UnstructuredInterpUtils::bilinearInterpStencil( Stencil&                          a_stencil,
                                                     const CFG::FArrayBox&             a_coords,
                                                     const std::vector<CFG::IntVect>&  a_idx_set,
                                                     const CFG::RealVect&              a_x )
{
  a_stencil.clear();

  /* 
   * Bilinear interpolation using 4 arbitrary points in 2D space 
  */

  CH_assert(a_idx_set.size() == 4);

  /* extract the (x,y) coords of the interpolation point */
  Real  xc = a_x[0],
        yc = a_x[1];

  /* extract the (x,y) coords of the four points */
  Real  x0 = a_coords(a_idx_set[0],0),
        y0 = a_coords(a_idx_set[0],1);
  Real  x1 = a_coords(a_idx_set[1],0),
        y1 = a_coords(a_idx_set[1],1);
  Real  x2 = a_coords(a_idx_set[2],0),
        y2 = a_coords(a_idx_set[2],1);
  Real  x3 = a_coords(a_idx_set[3],0),
        y3 = a_coords(a_idx_set[3],1);

  /* bilinear interp coefficients - computed using Mathematica */
  Real c0, c1, c2, c3;
  c0 = (x1*(-(xc*(y2 - y3)*(y1 - yc)) + x3*(y1 - y3)*(y2 - yc) - x2*(y1 - y2)*(y3 - yc)) - 
       x2*x3*(y2 - y3)*(y1 - yc) + x2*xc*(y1 - y3)*(y2 - yc) - x3*xc*(y1 - y2)*(y3 - yc))/
          (x0*(x3*(y1 - y2)*(y0 - y3) - x2*(y0 - y2)*(y1 - y3) + x1*(y0 - y1)*(y2 - y3)) + 
               x1*x2*(y1 - y2)*(y0 - y3) - x1*x3*(y0 - y2)*(y1 - y3) + x2*x3*(y0 - y1)*(y2 - y3));
  c1 = (x0*(xc*(y2 - y3)*(y0 - yc) - x3*(y0 - y3)*(y2 - yc) + x2*(y0 - y2)*(y3 - yc)) + 
       x2*x3*(y2 - y3)*(y0 - yc) - x2*xc*(y0 - y3)*(y2 - yc) + x3*xc*(y0 - y2)*(y3 - yc))/
          (x0*(x3*(y1 - y2)*(y0 - y3) - x2*(y0 - y2)*(y1 - y3) + x1*(y0 - y1)*(y2 - y3)) + 
               x1*x2*(y1 - y2)*(y0 - y3) - x1*x3*(y0 - y2)*(y1 - y3) + x2*x3*(y0 - y1)*(y2 - y3));
  c2 = (x0*(-(xc*(y1 - y3)*(y0 - yc)) + x3*(y0 - y3)*(y1 - yc) - x1*(y0 - y1)*(y3 - yc)) - 
       x1*x3*(y1 - y3)*(y0 - yc) + x1*xc*(y0 - y3)*(y1 - yc) - x3*xc*(y0 - y1)*(y3 - yc))/
          (x0*(x3*(y1 - y2)*(y0 - y3) - x2*(y0 - y2)*(y1 - y3) + x1*(y0 - y1)*(y2 - y3)) + 
               x1*x2*(y1 - y2)*(y0 - y3) - x1*x3*(y0 - y2)*(y1 - y3) + x2*x3*(y0 - y1)*(y2 - y3));
  c3 = (x0*(xc*(y1 - y2)*(y0 - yc) - x2*(y0 - y2)*(y1 - yc) + x1*(y0 - y1)*(y2 - yc)) + 
       x1*x2*(y1 - y2)*(y0 - yc) - x1*xc*(y0 - y2)*(y1 - yc) + x2*xc*(y0 - y1)*(y2 - yc))/
          (x0*(x3*(y1 - y2)*(y0 - y3) - x2*(y0 - y2)*(y1 - y3) + x1*(y0 - y1)*(y2 - y3)) + 
               x1*x2*(y1 - y2)*(y0 - y3) - x1*x3*(y0 - y2)*(y1 - y3) + x2*x3*(y0 - y1)*(y2 - y3));

  /* push these coefficients and coordinates to the stencil */
  a_stencil.push_back(StencilPoint(a_idx_set[0],c0));
  a_stencil.push_back(StencilPoint(a_idx_set[1],c1));
  a_stencil.push_back(StencilPoint(a_idx_set[2],c2));
  a_stencil.push_back(StencilPoint(a_idx_set[3],c3));

  /* done */
  return;
}

void UnstructuredInterpUtils::trilinearInterpStencil( Stencil&                          a_stencil,
                                                      const CFG::FArrayBox&             a_coords,
                                                      const std::vector<CFG::IntVect>&  a_idx_set,
                                                      const CFG::RealVect&              a_x )
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

  /* create LHS matrix */
  Real **Amat = (Real**) calloc (n_pts, sizeof(Real*));
  for (int i=0; i<n_pts; i++) {
    Amat[i] = (Real*) calloc (n_pts, sizeof(Real));
    Amat[i][0] = 1.0;
    Amat[i][1] = x[i];
    Amat[i][2] = y[i];
    Amat[i][3] = z[i];
    Amat[i][4] = x[i]*y[i];
    Amat[i][5] = y[i]*z[i];
    Amat[i][6] = x[i]*z[i];
    Amat[i][7] = x[i]*y[i]*z[i];
  }

  /* allocate memory for pivot array */
  int *pivot = (int*) calloc (n_pts+1, sizeof(int));

  /* decompose A */
  int status = LUP::decompose(Amat, n_pts, 1e-14, pivot);
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
  Real **Ainv = (Real**) calloc (n_pts, sizeof(Real*));
  for (int i=0; i<n_pts; i++) {
    Ainv[i] = (Real*) calloc (n_pts, sizeof(Real));
  }

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
  for (int n=0; n<n_pts; n++) a_stencil.push_back(StencilPoint(a_idx_set[n],interp_coeffs[n]));

  /* done */
  for (int i=0; i<n_pts; i++) {
    free(Amat[i]);
    free(Ainv[i]);
  }
  free(Amat);
  free(Ainv);
  free(pivot);
  return;
}

#include "NamespaceFooter.H"
