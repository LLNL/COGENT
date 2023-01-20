#include <array>
#include <cmath>
#include <float.h>
#include <fstream>
#include "OneBlockCoordSys.H"
#include "HermiteInterp.H"
#include "SplineInterp.H"
#include "BSplineInterp.H"
#include "Directions.H"
#include "BoxIterator.H"
#include "ConstFact.H"
#include "RadialUtils.H"
#include "CONSTANTS.H"

#undef REPORT_NEWTON_FAILURE

enum interpolation_package {HERMITE, SPLINE, BSPLINE};
enum interpolation_mode {FUNCTION, RADIAL_DERIVATIVE, POLOIDAL_DERIVATIVE, SECOND_RADIAL_DERIVATIVE, RADIAL_POLOIDAL_DERIVATIVE, SECOND_POLOIDAL_DERIVATIVE};
//enum field_vars {RBR_VAR, RBZ_VAR, RB_VAR, UNIT_BR_VAR, UNIT_BPHI_VAR, UNIT_BZ_VAR, NUM_FIELD_VARS};
enum rz_vars {R_VAR, Z_VAR};

#include "NamespaceHeader.H"

const std::string OneBlockCoordSys::pp_name = "oneblock";

OneBlockCoordSys::OneBlockCoordSys( ParmParse&      a_parm_parse,
                                            const ProblemDomain&  a_domain )
   : MagBlockCoordSys(a_parm_parse),
     m_RZ_interp(NULL)
{

   /// Read the input data specific to this geometry

   /// Get grid information values
   a_parm_parse.get("r_min", m_rmin);
   a_parm_parse.get("r_max", m_rmax);
   
   a_parm_parse.get("theta_min", m_thetamin);
   a_parm_parse.get("theta_max", m_thetamax);
   //   CH_assert(m_thetamin>=0.0 && m_thetamax<=3.141592653589793); 

#if CFG_DIM==3
   a_parm_parse.get("phi_min", m_phimin);
   a_parm_parse.get("phi_max", m_phimax);
#endif
   
   std::string conf_output; /// Output message after magnetic geometry is created
   /// Get magnetic field type
   /// Check for "B_type" entry
   if (a_parm_parse.contains("B_type")){
     std::string B_type;
     a_parm_parse.get("B_type", B_type);

     /// Constant magnetic type: B = B0
     if (B_type=="constant") {
       this->calcB = &OneBlockCoordSys::calcBConst;
       if (a_parm_parse.contains("Btor_0")) {
	 a_parm_parse.get("Btor_0", m_Btor_0);
	 conf_output = "Constant B field: Btor_0 = ";
	 conf_output.append(std::to_string(m_Btor_0));
       }
       else { MayDay::Error("ERROR: OneBlockCoordSys::OneBlockCoordSys() failed! (\"Btor_0\" is missing)");
       }
     }
     else { // Unknown type of magnetic geometry is specified
       conf_output = "ERROR: OneBlockCoordSys::OneBlockCoordSys() failed! Unknown magnetic geometry \"";
       conf_output.append(B_type);
       conf_output.append("\"");
       const char* cstr = conf_output.c_str();
       MayDay::Error(cstr);
     }
   }
   else { MayDay::Error("ERROR: OneBlockCoordSys::OneBlockCoordSys() failed! (\"B_type\" is missing)");
   }

   // Print geometry data to stdout
   if (m_verbose && !procID()) {
      cout << "Constructing oneblock geometry..." << endl;
      cout << "rmin = " << m_rmin << ", rmax = " << m_rmax<<endl;
      cout << "thetamin = " << m_thetamin << ", thetamax = " << m_thetamax<<endl;
#if CFG_DIM==3
      cout << "phi_min = " << m_phimin << ", phi_max = " << m_phimax<<endl;
#endif
      cout<< conf_output << endl;
   }

   IntVect dimensions = a_domain.size();
   double dr = (m_rmax - m_rmin)/(double)dimensions[0];
   double dtheta = (m_thetamax - m_thetamin)/(double)dimensions[1];
#if CFG_DIM==3
   double dphi = (m_phimax-m_phimin)/(double)dimensions[2];
   RealVect cellSpacing(dr,dtheta,dphi);
#endif

#if CFG_DIM==2
   RealVect cellSpacing(dr,dtheta);
#endif

   // Finish defining the object now that we also have the mesh spacing
   define( a_domain, cellSpacing );

   init( a_parm_parse );

   if (m_verbose && procID()==0) {
     cout << "Done constructing oneblock geometry" << endl;
   }
}

OneBlockCoordSys::~OneBlockCoordSys()
{
  if(m_RZ_interp) delete m_RZ_interp;
}

void
OneBlockCoordSys::init( ParmParse& a_pp )
{
   int n_radial_points, n_poloidal_points;
   int n_radial_extend, n_poloidal_extend;
   int num_block_elements;
   ifstream inFile;

#ifdef CH_MPI
   if (procID() == 0) {
#endif

      string geometry_file;
      if (a_pp.contains("geometry_file")) {
         a_pp.get("geometry_file", geometry_file);
      }
      else {
         MayDay::Error("OneBlockCoordSys::init(): No geometry mapping file specified");
      }

      inFile.open( geometry_file.c_str() );

      if (!inFile) {
         cout << "Unable to open geometry mapping file" << endl;
         exit(1);
      }

      cout << "Reading geometry mapping file" << endl;
 
      inFile >> n_radial_points;
      inFile >> n_radial_extend;
      inFile >> n_poloidal_points;
      inFile >> n_poloidal_extend;

      num_block_elements = (n_radial_points + 2*n_radial_extend) * (n_poloidal_points + 2*n_poloidal_extend);

#ifdef CH_MPI
   }

   MPI_Bcast(&n_radial_points, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n_radial_extend, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n_poloidal_points, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n_poloidal_extend, 1, MPI_INT, 0, MPI_COMM_WORLD);

#endif // end if MPI

   m_mapping_block_size[RADIAL_DIR]   = n_radial_points   + 2*n_radial_extend;
   m_mapping_block_size[POLOIDAL_DIR] = n_poloidal_points + 2*n_poloidal_extend;

   num_block_elements = m_mapping_block_size[RADIAL_DIR] * m_mapping_block_size[POLOIDAL_DIR];

   double *input_values = new double[2*num_block_elements];

#ifdef CH_MPI
   if (procID() == 0) {
#endif

      for (int i=0; i<2*num_block_elements; ++i) {
         inFile >> input_values[i];
      }
      inFile.close();
      cout << "Done reading geometry mapping file" << endl;

#ifdef CH_MPI
   }

   MPI_Bcast(input_values, 2*num_block_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif // end if MPI

   int num_radial_cells = n_radial_points - 1;

   double dr = ( upperMappedCoordinate(RADIAL_DIR)
             -   lowerMappedCoordinate(RADIAL_DIR) ) / num_radial_cells;

   double *r_pts = new double[m_mapping_block_size[RADIAL_DIR]];
            
   r_pts[0] = lowerMappedCoordinate(RADIAL_DIR) - n_radial_extend*dr;
   for (int i=1; i<m_mapping_block_size[RADIAL_DIR]; ++i) {
      r_pts[i] = r_pts[i-1] + dr;
   }
   
   int num_poloidal_cells = n_poloidal_points - 1;

   double dtheta = ( upperMappedCoordinate(POLOIDAL_DIR)
                 -   lowerMappedCoordinate(POLOIDAL_DIR) ) / num_poloidal_cells;
   
   double *theta_pts = new double[m_mapping_block_size[POLOIDAL_DIR]];
   
   theta_pts[0] = lowerMappedCoordinate(POLOIDAL_DIR) - n_poloidal_extend*dtheta;
   for (int i=1; i<m_mapping_block_size[POLOIDAL_DIR]; ++i) {
      theta_pts[i] = theta_pts[i-1] + dtheta;
   }

   // Determine interpolation method
   string interp_method = "bspline";
   if (a_pp.contains("interpolation_method")) {
      a_pp.get("interpolation_method", interp_method);
   }

   Box interp_box(IntVect::Zero,m_mapping_block_size-IntVect::Unit);
   
   FArrayBox interp_node_coords(interp_box, SpaceDim);

   // Create an interpolator for the cylindrical coordinates

   FArrayBox RZ_data(interp_box, 2);

   BoxIterator bit(interp_box);
   for (bit.begin();bit.ok();++bit) {
      IntVect iv = bit();
      int k = m_mapping_block_size[RADIAL_DIR]*iv[1] + iv[0];

      RZ_data(iv,0) = input_values[2*k];
      RZ_data(iv,1) = input_values[2*k+1];

      interp_node_coords(iv,0) = r_pts[iv[0] - interp_box.smallEnd(0)];
      interp_node_coords(iv,1) = theta_pts[iv[1] - interp_box.smallEnd(1)];
   }

   if(SpaceDim<3) {
      if (interp_method == "bspline") {
         m_RZ_interp = new BSplineInterp(a_pp, interp_node_coords, RZ_data);
      }
      else if (interp_method == "hermite") {
         m_RZ_interp = new HermiteInterp(a_pp, interp_node_coords, RZ_data);
      }
      else if (interp_method == "spline") {
         m_RZ_interp = new SplineInterp(a_pp, interp_node_coords, RZ_data);
      }
      else {
         MayDay::Error("OneBlockCoordSys::init(): Interp methods not setup for 3D configuration ");
      }
   }

   delete [] r_pts;
   delete [] theta_pts;
   delete [] input_values;

   getCellCenterRealCoords(); 

}

RealVect OneBlockCoordSys::realCoord( const RealVect& a_xi ) const
{
   RealVect x;
   D_TERM(x[0] = m_RZ_interp->interpolate(R_VAR, FUNCTION, a_xi);,
          x[1] = m_RZ_interp->interpolate(Z_VAR, FUNCTION, a_xi);,
          x[2] = 0.;)  
 
   return x;
}

RealVect OneBlockCoordSys::mappedCoord( const RealVect& a_X ) const
{
   /*
     Search for the cell whose centroid real space coordinate is closest to
     the passed argument and use the cooresponding mapped coordinate as the
     initial guess for the Newton solve.  This should be replaced by something
     more efficient.
   */

   RealVect xi_initial = RealVect::Zero;

   RealVect rc_dx = m_dx * m_rc_coarsen_ratio;
   RealVect offset = 0.5*RealVect::Unit;
   offset *= rc_dx;

   double max_dist = DBL_MAX;

   IntVect iv_closest;

   BoxIterator bit(m_realCoords.box());
   for (bit.begin();bit.ok();++bit) {
      IntVect iv = bit();
      RealVect center_xi = rc_dx*iv + offset;
      RealVect center_X;
      for (int dir=0; dir<SpaceDim; ++dir) {
         center_X[dir] = m_realCoords(iv,dir);
      }
      RealVect diff = center_X;
      diff -= a_X;
      double dist = diff.vectorLength();
      if (dist < max_dist) {
         xi_initial = center_xi;
         max_dist = dist;
         iv_closest = iv;
      }
   }
 
   return mappedCoordNewton( a_X, xi_initial, iv_closest );

}

RealVect
OneBlockCoordSys::mappedCoordNewton( const RealVect& a_X,
                                            const RealVect& a_xi_initial,
                                            const IntVect& a_iv_initial ) const
{
   /*
     Use Newton iteration to evaluate the mapping of physical to computational
     coordinates by inverting the mapping of computational to physical coordinates.
   */
   double tol = 1.e-10;         // Read from input if we find some reason to later
   int max_iter = 40;

   RealVect xi = a_xi_initial;

   RealVect residual;

   double R = m_RZ_interp->interpolate(R_VAR, FUNCTION, xi);
   double Z = m_RZ_interp->interpolate(Z_VAR, FUNCTION, xi);
   double dRdr = m_RZ_interp->interpolate(R_VAR, RADIAL_DERIVATIVE, xi);
   double dRdtheta = m_RZ_interp->interpolate(R_VAR, POLOIDAL_DERIVATIVE, xi);
   double dZdr = m_RZ_interp->interpolate(Z_VAR, RADIAL_DERIVATIVE, xi);
   double dZdtheta = m_RZ_interp->interpolate(Z_VAR, POLOIDAL_DERIVATIVE, xi);

   residual[0] = R - a_X[0];
   residual[1] = Z - a_X[1];

   double Fnorm = residual.vectorLength();

   bool converged = Fnorm <= tol;
   int num_iters = 0;
   int bt_steps;

   RealVect xi_saved[100];
   double Fnorm_saved[100];

   xi_saved[0] = xi;
   Fnorm_saved[0] = Fnorm;

   while ( !converged && num_iters < max_iter) {

      // Invert the Jacobian to get the update
      double a = dRdr;
      double b = dRdtheta;
      double c = dZdr;
      double d = dZdtheta;
      double J = a*d - b*c;

      RealVect delta;
      delta[0] = - ( d * residual[0] - b * residual[1] ) / J;
      delta[1] = - (-c * residual[0] + a * residual[1] ) / J;

      double s = 1.;

      RealVect xi_proposed = xi + delta;

      R = m_RZ_interp->interpolate(R_VAR, FUNCTION, xi_proposed);
      Z = m_RZ_interp->interpolate(Z_VAR, FUNCTION, xi_proposed);
      dRdr = m_RZ_interp->interpolate(R_VAR, RADIAL_DERIVATIVE, xi_proposed);
      dRdtheta = m_RZ_interp->interpolate(R_VAR, POLOIDAL_DERIVATIVE, xi_proposed);
      dZdr = m_RZ_interp->interpolate(Z_VAR, RADIAL_DERIVATIVE, xi_proposed);
      dZdtheta = m_RZ_interp->interpolate(Z_VAR, POLOIDAL_DERIVATIVE, xi_proposed);

      residual[0] = R - a_X[0];
      residual[1] = Z - a_X[1];

      double Fnorm_proposed = residual.vectorLength();

#if 1
      /*
        Backtracking algorithm from R. P. Pawlowski, J. N. Shadid, J. P. Simonis and H. F. Walker,
        "Globalization Techniques for Newton-Krylov Methods and Applications to the Fully Coupled
        Solution of the Navier-Stokes Equations", SIAM Review, 48 (4), 2006, pp. 700-721.
      */

      double t = 1.e-4;
      double eta = 0.;
      bt_steps = 0;
      int bt_steps_max = 15;

      bool backtracking = Fnorm_proposed > (1. - t*(1. - eta)) * Fnorm && s > 0.;

      while ( backtracking ) {

         double theta = 0.1;
         delta *= theta;
         eta = 1. - theta*(1. - eta);

         xi_proposed = xi + delta;

         R = m_RZ_interp->interpolate(R_VAR, FUNCTION, xi_proposed);
         Z = m_RZ_interp->interpolate(Z_VAR, FUNCTION, xi_proposed);
         dRdr = m_RZ_interp->interpolate(R_VAR, RADIAL_DERIVATIVE, xi_proposed);
         dRdtheta = m_RZ_interp->interpolate(R_VAR, POLOIDAL_DERIVATIVE, xi_proposed);
         dZdr = m_RZ_interp->interpolate(Z_VAR, RADIAL_DERIVATIVE, xi_proposed);
         dZdtheta = m_RZ_interp->interpolate(Z_VAR, POLOIDAL_DERIVATIVE, xi_proposed);

         residual[0] = R - a_X[0];
         residual[1] = Z - a_X[1];

         Fnorm_proposed = residual.vectorLength();

         backtracking = (Fnorm_proposed > (1. - t*(1. - eta)) * Fnorm) && (bt_steps < bt_steps_max);

         bt_steps++;
      }
#endif

      // Update the current solution and resisual norm
      xi = xi_proposed;
      Fnorm = Fnorm_proposed;

      // Test convergence
      converged = Fnorm <= tol;

      num_iters++;

      xi_saved[num_iters] = xi;
      Fnorm_saved[num_iters] = Fnorm;
   }

#ifdef REPORT_NEWTON_FAILURE
   if ( !converged && num_iters >= max_iter ) {
      cout << "Newton solve did not converge at " << a_X << ", a_xi_initial = " << a_xi_initial << ", iv_init = " << a_iv_initial << ", Fnorm = " << Fnorm << endl;
#if 0
      for (int i=0; i<num_iters; ++i) {
         cout << i << " xi = " << xi_saved[i] << ", Fnorm = " << Fnorm_saved[i] << endl;
      }
#endif
   }
#endif

   return xi;

}


void
OneBlockCoordSys::getCellCenterRealCoords()
{
   Box box = m_domain.domainBox();
   
   m_rc_coarsen_ratio = 8;
   for (int dir=0; dir<SpaceDim; ++dir) {
     if (box.size(dir) < 64) {
       m_rc_coarsen_ratio = 1;
       break;
     }
   }

   box.coarsen(m_rc_coarsen_ratio);

   m_realCoords.define(box, 2);

   RealVect rc_dx = m_dx * m_rc_coarsen_ratio;

   RealVect offset = 0.5*RealVect::Unit;
   offset *= rc_dx;

   BoxIterator bit(box);
   for (bit.begin();bit.ok();++bit) {
      IntVect iv = bit();
      RealVect xi = rc_dx*iv + offset;
      RealVect X = realCoord(xi);
      for (int dir=0; dir<SpaceDim; ++dir) {
         m_realCoords(iv,dir) = X[dir];
      }
   }
   
}

Real 
OneBlockCoordSys::dXdXi( const RealVect&  a_Xi,
                                  int            a_dirX,
                                  int            a_dirXi ) const
{
  if (!m_subgrid_geometry) return m_RZ_interp->interpolate(a_dirX, a_dirXi+1, a_Xi);
  else {
    RealVect a_Xi_tmp(a_Xi);
    RealVect x = realCoord(a_Xi_tmp);

    RealVect Xi_hi(a_Xi_tmp);
    RealVect Xi_lo(a_Xi_tmp);

    Xi_hi[a_dirXi] = a_Xi[a_dirXi] +  m_dx[a_dirXi]/m_N_mr[a_dirXi];
    Xi_lo[a_dirXi] = a_Xi[a_dirXi] -  m_dx[a_dirXi]/m_N_mr[a_dirXi];

    RealVect x_hi = realCoord(Xi_hi);
    RealVect x_lo = realCoord(Xi_lo);

    RealVect metrics;
    for (int dir(0); dir<SpaceDim; ++dir) {
      metrics[dir] = (x_hi[dir] - x_lo[dir])/m_dx[a_dirXi] * m_N_mr[a_dirXi]/2.0;
    }

    return metrics[a_dirX];
  }


}

void 
OneBlockCoordSys::dXdXi( FArrayBox&  a_dXdXi,
                          const FArrayBox&  a_Xi,
                                int         a_destComp,
                                int         a_dirX,
                                int         a_dirXi,
                          const Box&        a_box ) const
{
  if (!m_subgrid_geometry) m_RZ_interp->interpolate(a_dirX, a_destComp, a_dirXi+1, a_Xi, a_box, a_dXdXi);

  else {
     BoxIterator bit(a_box);
     for (bit.begin(); bit.ok(); ++bit) {
        const IntVect& iv = bit();
        RealVect this_Xi;
        for (int dir=0; dir<SpaceDim; ++dir) {
           this_Xi[dir] = a_Xi(iv,dir);
        }
        a_dXdXi(iv,a_destComp) = dXdXi(this_Xi,a_dirX,a_dirXi);
     }
  }
}

bool
OneBlockCoordSys::isValid(const RealVect& xi) const
{
  return xi[RADIAL_DIR]   >= lowerMappedCoordinate(RADIAL_DIR)   &&
         xi[RADIAL_DIR]   <= upperMappedCoordinate(RADIAL_DIR)   &&
         xi[POLOIDAL_DIR] >= lowerMappedCoordinate(POLOIDAL_DIR) &&
         xi[POLOIDAL_DIR] <= upperMappedCoordinate(POLOIDAL_DIR);
}

void
OneBlockCoordSys::computeFieldData( const int   a_dir,
                                        FArrayBox&  a_BField,
                                        FArrayBox&  a_BFieldMag,
                                        FArrayBox&  a_BFieldDir,
                                        FArrayBox&  a_gradBFieldMag,
                                        FArrayBox&  a_curlBFieldDir,
                                        FArrayBox&  a_BFieldDirdotcurlBFieldDir,
                                  const bool        a_derived_data_only ) const
// NEED TO FOLLOW AND SEE WHAT IS GOTTEN IN FORTRAN
{
   //MayDay::Error("OneBlockCoordSys::computeFieldData() not implemented");
  
   // Get box intersection
   Box box = a_BField.box();
   box &= a_BFieldMag.box();
   box &= a_BFieldDir.box();
   box &= a_gradBFieldMag.box();
   box &= a_curlBFieldDir.box();
   box &= a_BFieldDirdotcurlBFieldDir.box();

   FArrayBox RZ(box,SpaceDim);
   if (a_dir>=0 && a_dir<SpaceDim) {
      getFaceCenteredRealCoords(a_dir, RZ);
   }
   else {
      getCellCenteredRealCoords(RZ);
   }

   for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();
      
      double R = RZ(iv,0);
      double Z = RZ(iv,1);
      double p_B[17];	/// Array for magnetic filed configuration
	  
      // Calculate magnetic geometry components
      (this->*calcB)(R, Z, p_B);
	
      a_BField(iv,0) = p_B[0]; /// B[0]
      a_BField(iv,1) = p_B[1]; /// B[1]
      a_BField(iv,2) = p_B[2]; /// B[2]
      
      a_BFieldMag(iv,0) = p_B[3]; /// abs(B)
      
      a_BFieldDir(iv,0) = p_B[4]; /// b[0], b=b/abs(B)
      a_BFieldDir(iv,1) = p_B[5]; /// b[1]
      a_BFieldDir(iv,2) = p_B[6]; /// b[2]

      a_gradBFieldMag(iv,0) = p_B[7]; /// gB[0], gB = \nabla abs(B)
      a_gradBFieldMag(iv,1) = p_B[8]; /// gB[1]
      a_gradBFieldMag(iv,2) = p_B[9]; /// gB[2]
      
      a_curlBFieldDir(iv,0) = p_B[10]; /// cb[1], cb = \nabla \times b
      a_curlBFieldDir(iv,1) = p_B[11]; /// cb[2]
      a_curlBFieldDir(iv,2) = p_B[12]; /// cb[3]
      
      a_BFieldDirdotcurlBFieldDir(iv,0) = p_B[13]; /// b \cdot \nabla \times b
   }
}

void
OneBlockCoordSys::calcBConst( double a_R, double a_Z, double* a_p_B ) const
{
  double sign_b; if (m_Btor_0>0) {sign_b = 1.0;} else {sign_b = -1.0;}
  
  a_p_B[0] = 0.0;
  a_p_B[1] = m_Btor_0;
  a_p_B[2] = 0.0;
      
  a_p_B[3] = fabs( a_p_B[1] );
      
  a_p_B[4] = 0.0;
  a_p_B[5] = sign_b;
  a_p_B[6] = 0.0;

  a_p_B[7] = 0.0;
  a_p_B[8] = 0.0;
  a_p_B[9] = 0.0;
     
  a_p_B[10] = 0.0;
  a_p_B[11] = 0.0;
  a_p_B[12] = sign_b / a_R;
      
  a_p_B[13] = 0.0;
  
  a_p_B[14] = a_Z * a_p_B[1];
  a_p_B[15] = 0.0;
  a_p_B[16] = 0.0;
}

array<double,3>
OneBlockCoordSys::computeBField(const RealVect& a_X) const
{
   MayDay::Error("OneBlockCoordSys::computeBField() not implemented");
   
   array<double,3> result;
   
   double p_B[17];
   (this->*calcB)(a_X[0], a_X[1], p_B);
   
   result[0] =  p_B[0];
   result[1] =  p_B[1];
   result[2] =  p_B[2];
   
   return result;
}

void OneBlockCoordSys::getMagneticFlux( const FArrayBox& physical_coordinates,
                                                FArrayBox&       magnetic_flux ) const
{
   MayDay::Error("OneBlockCoordSys::getMagneticFlux() not implemented");
}



double OneBlockCoordSys::getMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   MayDay::Error("OneBlockCoordSys::getMagneticFlux() not implemented");
   return(0.0);
}


void OneBlockCoordSys::getNodalFieldData(FArrayBox& points,
                                                 FArrayBox& A,
                                                 FArrayBox& b,
                                                 FArrayBox& Bmag) const
{
   MayDay::Error("OneBlockCoordSys::getNodalFieldData() not implemented");
  
  for (BoxIterator bit(points.box()); bit.ok(); ++bit) {
      IntVect iv = bit();

      RealVect xi;
      for (int n=0; n<SpaceDim; ++n) {
         xi[n] = points(iv,n);
      }

      RealVect X = realCoord(xi);

      double p_B[17];
      (this->*calcB)(X[0], X[1], p_B);

      A(iv,0) = p_B[14];    // A_r
      A(iv,1) = p_B[15];    // A_\theta
      A(iv,2) = p_B[16];    // A_z

      Bmag(iv,0) = p_B[3];  // abs(B)
      
      b(iv,0) = p_B[4];     // b_r
      b(iv,1) = p_B[5];     // b_\theta
      b(iv,2) = p_B[6];     // b_z

   }

}

double OneBlockCoordSys:: getRBtoroidal() const
{
    MayDay::Error("ERROR: Magnetic geometry is not toroidal; divfree_velocity is defined for toroidal geometry only");
    return 0.0;
}

#include "NamespaceFooter.H"
