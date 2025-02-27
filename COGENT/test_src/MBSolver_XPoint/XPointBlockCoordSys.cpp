#include "XPointBlockCoordSys.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

const char* XPointBlockCoordSys::pp_name = {"XPoint"};

static double affine[8][4] = { {1,0,1,1}, {1,1,0,1}, {1,1,-1,0}, {0,1,-1,-1}, {-1,0,-1,-1}, {-1,-1,0,-1}, {-1,-1,1,0}, {0,-1,1,1} };

XPointBlockCoordSys::XPointBlockCoordSys( ParmParse&               a_parm_parse,
                                          const ProblemDomain&     a_domain,
                                          const int                a_block_type)
   : MagBlockCoordSys(a_parm_parse),
     m_block_type(a_block_type),
     m_alpha(3.-sqrt(5.))
     //     m_alpha(0.)
{
   if (SpaceDim != 2) {
      MayDay::Error("XPointBlockCoordSys is only two-dimensional");
   }

    a_parm_parse.getarr("origin", m_origin, 0, SpaceDim);

    // Compute the mesh size
    IntVect dimensions = a_domain.size();
    RealVect cellSpacing;
    for (int dir=0; dir<SpaceDim; ++dir) {
       cellSpacing[dir] = 1. / (double)dimensions[dir];
    }

    // Finish defining the object now that we also have the mesh spacing
    define( a_domain, cellSpacing );

#ifdef PLOT_INVERSE_POINTS
   char file_name[80];
   sprintf(file_name, "Block%d_%d", m_block_type, procID());
   m_ipt_file[m_block_type] = fopen(file_name, "w");

   if (procID()==0) {

     char file_name2[80];
     sprintf(file_name2, "Pts%d", m_block_type);
     FILE *fptr = fopen(file_name2, "w");

     Box box = a_domain.domainBox();
     box.surroundingNodes();

     BoxIterator bit(box);
     for (bit.begin(); bit.ok(); ++bit) {
        IntVect iv = bit();
        RealVect xi;
        for (int dir=0; dir<SpaceDim; ++dir) {
          xi[dir] = iv[dir]*cellSpacing[dir];
        }
        RealVect X = realCoord(xi);
        fprintf(fptr, "%e %e\n", X[0], X[1]);
     }

     fclose(fptr);
   }


#endif
}



XPointBlockCoordSys::~XPointBlockCoordSys()
{
#ifdef PLOT_INVERSE_POINTS
   for (int block=0; block<8; ++block) {
      fclose(m_ipt_file[block]);
   }
#endif
}



RealVect XPointBlockCoordSys::realCoord( const RealVect& a_xi ) const
{
   RealVect xi_loc;
   xi_loc[0] = a_xi[0] - lowerMappedCoordinate(0);
   xi_loc[1] = a_xi[1] - lowerMappedCoordinate(1);

   double R_loc = 0.25 * xi_loc[0] * (3. + pow(1. - m_alpha*xi_loc[1],2));
   double Z_loc = 0.25 * xi_loc[1] * (3. + pow(1. - m_alpha*xi_loc[0],2));

   double R = affine[m_block_type][0] * R_loc + affine[m_block_type][1] * Z_loc;
   double Z = affine[m_block_type][2] * R_loc + affine[m_block_type][3] * Z_loc;

   RealVect x;
   D_TERM(x[0] = R + m_origin[0];,
          x[1] = Z + m_origin[1];,
          x[2] = 0.;)

   return x;
}



RealVect XPointBlockCoordSys::mappedCoord( const RealVect& a_X ) const
{
#ifdef PLOT_INVERSE_POINTS
   fprintf(m_ipt_file[m_block_type], "%20.12e %20.12e\n", a_X[0], a_X[1]);
#endif
   /*
     Use Newton iteration to evaluate the mapping of physical to computational
     coordinates by inverting the mapping of computational to physical coordinates.
   */
   double tol = 1.e-12;         // Read from input if we find some reason to later
   double max_iter = 100;

   // Use the mapped block centroid as the initial guess
   RealVect xi;
   for (int dir=0; dir<SpaceDim; ++dir) {
      xi[dir] = 0.5 * ( lowerMappedCoordinate(dir) + upperMappedCoordinate(dir) );
   }

   double conv_tol = (xi.vectorLength() > 0.)? tol * xi.vectorLength(): tol;
   bool converged = false;
   int num_iters = 0;

   while ( !converged ) {

      if (num_iters > max_iter) {
         MayDay::Error("XPointBlockCoordSys::mappedCoord(): Newton iteration failed to converge");
      }

      RealVect residual = realCoord(xi) - a_X;

      // Invert the Jacobian to get the update
      double a = dXdXi(xi,0,0);
      double b = dXdXi(xi,0,1);
      double c = dXdXi(xi,1,0);
      double d = dXdXi(xi,1,1);
      double J = a*d - b*c;

      RealVect delta;
      delta[0] = ( d * residual[0] - b * residual[1] ) / J;
      delta[1] = (-c * residual[0] + a * residual[1] ) / J;

      // Update the current solution
      xi -= delta;

      converged = delta.vectorLength() <= conv_tol;
      num_iters++;
   }

   return xi;
}



Real XPointBlockCoordSys::dXdXi( const RealVect& a_Xi,
                                 int             a_dirX,
                                 int             a_dirXi ) const
{
   RealVect xi_loc;
   xi_loc[0] = a_Xi[0] - lowerMappedCoordinate(0);
   xi_loc[1] = a_Xi[1] - lowerMappedCoordinate(1);

   double d00 = 0.75 + 0.25 * pow(1. - m_alpha * xi_loc[1], 2);
   double d01 = -0.5 * m_alpha * xi_loc[0] * (1. - m_alpha * xi_loc[1]);
   double d10 = -0.5 * m_alpha * xi_loc[1] * (1. - m_alpha * xi_loc[0]);
   double d11 = 0.75 + 0.25 * pow(1. - m_alpha * xi_loc[0], 2);

   double value;

   if (a_dirX == 0) {
      if (a_dirXi == 0) {
         value = affine[m_block_type][0] * d00 + affine[m_block_type][1] * d10;
      }
      else if (a_dirXi == 1) {
         value = affine[m_block_type][0] * d01 + affine[m_block_type][1] * d11;
      }
   }
   else if (a_dirX == 1) {
      if (a_dirXi == 0) {
         value = affine[m_block_type][2] * d00 + affine[m_block_type][3] * d10;
      }
      else if (a_dirXi == 1) {
         value = affine[m_block_type][2] * d01 + affine[m_block_type][3] * d11;
      }
   }
   else {
      MayDay::Error("Bad dirX in XPointBlockCoordSys::dXdXi");
   }

   return value;
}



void XPointBlockCoordSys::dXdXi(FArrayBox&       a_dXdXi,
                                const FArrayBox& a_Xi,
                                int              a_destComp,
                                int              a_dirX,
                                int              a_dirXi,
                                const Box&       a_box) const
{
   BoxIterator bit(a_box);
   for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect xi;
      for (int dir=0; dir<SpaceDim; ++dir) {
         xi[dir] = a_Xi(iv,dir);
      }
      a_dXdXi(iv,a_destComp) = dXdXi(xi, a_dirX, a_dirXi);
   }
}



Real XPointBlockCoordSys::dX2dXi2( const RealVect& a_Xi,
                                   int             a_dirX,
                                   int             a_dirXi0,
                                   int             a_dirXi1 ) const
{
   return 0.;
}



void
XPointBlockCoordSys::getRBpoloidal( const int  a_dir,
                                    FArrayBox& a_RBpol,
                                    FArrayBox& a_dRBpoldt ) const
{
   a_RBpol.setVal(1.);
   a_dRBpoldt.setVal(0.);
}


#include "NamespaceFooter.H"
