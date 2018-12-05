#include "MillerBlockCoordSys.H"
#include "MagBlockCoordSysF_F.H"
#include "MillerBlockCoordSysF_F.H"

#include "Directions.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

//char* MillerBlockCoordSys::pp_name = {"miller"};

const std::string MillerBlockCoordSys::pp_name = "miller";

MillerBlockCoordSys::MillerBlockCoordSys( ParmParse&               a_parm_parse,
                                          const ProblemDomain&     a_domain )
   : MagBlockCoordSys(a_parm_parse)
{
  /*     if (SpaceDim != 2) {
      MayDay::Error("MillerBlockCoordSys is only two-dimensional");
      }
  */
  // ABOVE WAS OLD TEST BEFORE CONVERT TO 3D

   // Read the input data specific to this geometry

   if (a_parm_parse.contains("inner_radial_bdry")) {
      a_parm_parse.get("inner_radial_bdry", m_inner_radial_boundary);
   }
   
   if (a_parm_parse.contains("outer_radial_bdry")) {
      a_parm_parse.get("outer_radial_bdry", m_outer_radial_boundary);
   }

   m_average_radius = 0.5*(m_inner_radial_boundary + m_outer_radial_boundary);

   // Get the flag to use a modified geometry where the minor radius is replaced by
   // its average, to enable radial periodicity (e.g. for GAM tests)
   if (a_parm_parse.contains("l_const_minorrad")) {
     a_parm_parse.get("l_const_minorrad",m_const_minorrad);
   }
   else {
      m_const_minorrad = 0;
   }

   // If the constant minor radius option is selected, make sure that the radial
   // boundaries are periodic.  This restriction is due to the need to multiply
   // and divide by J in the application of radial boundary conditions.  Currently,
   // the code that does that does not know about the constant minor radius option.

   if ( m_const_minorrad != 0 && !a_domain.isPeriodic(RADIAL_DIR ) ) {
       MayDay::Error("Radial direction must be periodic with constant minor radius option");
   }
   
   // Get the coordinate system offset
   a_parm_parse.getarr("origin", m_origin, 0, SpaceDim);

   // Get the triangularity parameter
   double delta;
   a_parm_parse.get("delta", delta);
   m_beta = asin(delta);

   // Get the elongation parameter
   a_parm_parse.get("kappa", m_kappa);

   // Get the other parameters
   a_parm_parse.get("dpsidr", m_dpsidr);
   a_parm_parse.get("drR0", m_drR0);
   a_parm_parse.get("s_kappa", m_s_kappa);
   a_parm_parse.get("s_delta", m_s_delta);

   // Get the toroidal field component scale factor
   if (a_parm_parse.contains("Btor_scale")) {
      a_parm_parse.get("Btor_scale", m_Btor_scale);
   }
   else {
      m_Btor_scale = 0.;  // default
   }

   // Print geometry data to stdout
   if (m_verbose && !procID()) {
      cout << "Constructing Miller equilibrium geometry..." << endl;
      cout << "   delta = " << delta << endl;
      cout << "   kappa = " << m_kappa << endl;
      cout << "   drR0 = " << m_drR0 << ", s_kappa = " << m_s_kappa << ", m_dpsidr = " << m_dpsidr
           << ", s_delta = " << m_s_delta << endl;
      cout << "   origin =";
      for (int n=0; n<SpaceDim; n++) {
         cout << " " << m_origin[n];
      }
      cout << endl;
      cout << "   inner_radial_bdry = " << m_inner_radial_boundary << endl;
      cout << "   outer_radial_bdry = " << m_outer_radial_boundary << endl;
   }

   // Compute the mesh size
   IntVect dimensions = a_domain.size();
   double dr = (m_outer_radial_boundary - m_inner_radial_boundary) / (double)dimensions[0];
   double dtheta = 2.*Pi/(double)dimensions[1];
#if CFG_DIM==3
   double m_phimax = 2.*Pi; //for quick debug
   double m_phimin = 0.0;
   double dphi = (m_phimax - m_phimin)/(double)dimensions[2];
   RealVect cellSpacing(dr,dtheta,dphi);
#endif

#if CFG_DIM==2
   RealVect cellSpacing(dr,dtheta);
#endif

   // Finish defining the object now that we also have the mesh spacing
   define( a_domain, cellSpacing );

   if (m_verbose && procID()==0) {
      cout << "Done constructing Miller equilibrium geometry" << endl;
   }
}



RealVect MillerBlockCoordSys::realCoord( const RealVect& a_xi ) const
{
   RealVect x;
   double r = a_xi[0] + m_inner_radial_boundary;

   double theta = a_xi[1];
   D_TERM(x[0] = m_origin[0] + r*cos(theta+m_beta*sin(theta));,
          x[1] = m_origin[1] + r*m_kappa*sin(theta);,
          x[2] = a_xi[2];)

   return x;
}



RealVect MillerBlockCoordSys::mappedCoord( const RealVect& a_x ) const
{
   double x_translated = a_x[0] - m_origin[0];
   double y_translated = a_x[1] - m_origin[1];
   double r;
   double theta;
   double two_pi = 2.*Pi;

   /*
      We need to solve for theta, since it cannot be expressed
      analytically in terms of x and y for nonzero m_beta.

      We use a Newton iteration.  The initial guess for theta is exact
      when delta = m_beta = 0, and it should also be pretty good for
      nonzero values of the latter.  The iteration will be performed
      in the domain -pi < theta <= pi.  Before returning, theta is
      converted to the interval 0 <= theta < 2pi.
   */
   theta = atan2(y_translated, x_translated);

   double st, ct, t1, st1, ct1;

   while (true) {

      st  = sin(theta);
      ct  = cos(theta);
      t1  = theta + m_beta*st;
      st1 = sin(t1);
      ct1 = cos(t1);

      double g = m_kappa*x_translated*st - y_translated*ct1;

      if (fabs(g) > 1.e-11) {
         double gprime = m_kappa*x_translated*ct + y_translated*st1*(1. + m_beta*ct);
         theta -= g/gprime;
      }
      else break;
   }

   if( fabs(st) > 1.e-11) {
      r = y_translated / (m_kappa * st);
   }
   else {
      r = x_translated / ct1;
   }

   // Now shift theta to the range 0 <= theta < 2pi.
   if (y_translated < 0) theta += two_pi;

   RealVect xi;
   D_TERM(xi[0] = r - m_inner_radial_boundary;,
          xi[1] = theta;,
          xi[2] = a_x[2];)

   return xi;
}



Real MillerBlockCoordSys::dXdXi( const RealVect& a_Xi,
                                 int             a_dirX,
                                 int             a_dirXi ) const
{
  Real value = 0.0;
  // this is for test purposes

  double r = a_Xi[0] + m_inner_radial_boundary;
  double theta = a_Xi[1];

  if (a_dirX == 0)
    {
      if (a_dirXi == 0)
        {
          // dx/dr = cos(theta + beta*sin(theta))
          value = cos(theta + m_beta*sin(theta));
        }
      else if (a_dirXi == 1)
        {
          // dx/dtheta = -rsin(theta + beta*sin(theta))*(1 + beta*cos(theta))
          value = -r*sin(theta + m_beta*sin(theta))*(1. + m_beta*cos(theta));
        }
      // else value stays 0
    }
  else if (a_dirX == 1)
    {
      if (a_dirXi == 0)
        {
          // dy/dr = kappa*sin(theta)
          value = m_kappa * sin(theta);
        }
      else if (a_dirXi == 1)
        {
          // dy/dtheta = r kappa cos(theta)
          value = r * m_kappa * cos(theta);
        }
      // else value stays 0
    }
  else if (a_dirX == 2)
    {
      if (a_dirXi == 2)
        {
          value = 1.0;
        }
    }
  else
    {
       MayDay::Error("Bad dirX in MillerBlockCoordSys::dXdXi");
    }

  return value;
}



void MillerBlockCoordSys::dXdXi(FArrayBox&       a_dXdXi,
                                const FArrayBox& a_Xi,
                                int              a_destComp,
                                int              a_dirX,
                                int              a_dirXi,
                                const Box&       a_box) const
{
   FArrayBox Xi_shifted(a_Xi.box(), a_Xi.nComp());
   Xi_shifted.copy(a_Xi);
   Xi_shifted.plus(m_inner_radial_boundary,0,1);

   FORT_MILLER_DXDXI(CHF_BOX(a_box),
                     CHF_CONST_REAL(m_beta),
                     CHF_CONST_REAL(m_kappa),
                     CHF_CONST_FRA(Xi_shifted),
                     CHF_CONST_INT(a_destComp),
                     CHF_CONST_INT(a_dirX),
                     CHF_CONST_INT(a_dirXi),
                     CHF_FRA(a_dXdXi));
}



Real MillerBlockCoordSys::dX2dXi2( const RealVect& a_Xi,
                                   int             a_dirX,
                                   int             a_dirXi0,
                                   int             a_dirXi1) const
{
  Real value = 0.0;
  // this is for test purposes

  double r = a_Xi[0] + m_inner_radial_boundary;
  double theta = a_Xi[1];

  if (a_dirX == 0) {
     if (a_dirXi0 == 1 && a_dirXi1 == 1) {
        //    dR2/dtheta2
        double stpbs = sin(theta + m_beta*sin(theta));
        double ctpbs = cos(theta + m_beta*sin(theta));
        double opct = 1. + m_beta * cos(theta);

        value = r * (- ctpbs * opct*opct + stpbs * m_beta * sin(theta));
     }
     else {
        MayDay::Error("MillerBlockCoordSys::dXdXi2(): branch not implemented");
     }
  }
  else if (a_dirX == 1) {
     if (a_dirXi0 == 1 && a_dirXi1 == 1) {
        //    dZ2/dtheta2
        value = - r * m_kappa * sin(theta);
     }
     else {
        MayDay::Error("MillerBlockCoordSys::dXdXi2(): branch not implemented");
     }
    }
  else {
     MayDay::Error("MillerBlockCoordSys::dXdXi2(): branch not implemented");
  }

  return value;
}



Real
MillerBlockCoordSys::majorRadiusConstMinorRadAdj(const RealVect a_xi) const
{
  RealVect xi(a_xi);

  if (m_const_minorrad == 1) {
    xi[RADIAL_DIR] = m_average_radius - m_inner_radial_boundary;
  }

  return majorRadius(xi);
}


void
MillerBlockCoordSys::getRBpoloidal( const int  a_dir,
                                    FArrayBox& a_RBpol,
                                    FArrayBox& a_dRBpoldt ) const
{
   const Box& box = a_RBpol.box();
   CH_assert(box == a_dRBpoldt.box());

   FORT_GET_RBPOL_MILLER(CHF_CONST_INT(a_dir),
                         CHF_BOX(box),
                         CHF_CONST_REALVECT(m_dx),
                         CHF_CONST_REAL(m_beta),
                         CHF_CONST_REAL(m_kappa),
                         CHF_CONST_REAL(m_dpsidr),
                         CHF_CONST_REAL(m_origin[0]),
                         CHF_CONST_REAL(m_drR0),
                         CHF_CONST_REAL(m_s_kappa),
                         CHF_CONST_REAL(m_s_delta),
                         CHF_FRA1(a_RBpol,0),
                         CHF_FRA1(a_dRBpoldt,0));
}


#if 0
void
MillerBlockCoordSys::computeFieldData( const int  a_dir,
                                       FArrayBox& a_BField,
                                       FArrayBox& a_BFieldMag,
                                       FArrayBox& a_BFieldDir,
                                       FArrayBox& a_gradBFieldMag,
                                       FArrayBox& a_curlBFieldDir,
                                       FArrayBox& a_BFieldDirdotcurlBFieldDir,
                                       const bool a_derived_data_only ) const
{
   // Get box intersection
   Box box = a_BField.box();
   box &= a_BFieldMag.box();
   box &= a_BFieldDir.box();
   box &= a_gradBFieldMag.box();
   box &= a_curlBFieldDir.box();
   box &= a_BFieldDirdotcurlBFieldDir.box();

   FArrayBox RBpol(box, 1);
   FArrayBox dRBpoldt(box, 1);

   getRBpoloidal(a_dir, RBpol, dRBpoldt);

   FArrayBox R(box, 1);
   FArrayBox dRdr(box, 1);
   FArrayBox dRdt(box, 1);
   FArrayBox dZdr(box, 1);
   FArrayBox dZdt(box, 1);
   FArrayBox dR2dt2(box, 1);
   FArrayBox dZ2dt2(box, 1);

   // Compute centering offset
   RealVect offset = 0.5*RealVect::Unit;
   offset *= m_dx;
   for (int dir=0; dir<SpaceDim; ++dir) {
      if (dir == a_dir) offset[dir] = 0.;
   }

   BoxIterator bit(box);
   for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect xi = m_dx*iv + offset;

      if (m_const_minorrad == 1) {
         xi[RADIAL_DIR] = m_average_radius - m_inner_radial_boundary;
      }

      R(iv)    = realCoord(xi)[RADIAL_DIR];
      dRdr(iv) = dXdXi(xi, 0, 0);
      dRdt(iv) = dXdXi(xi, 0, 1);
      dZdr(iv) = dXdXi(xi, 1, 0);
      dZdt(iv) = dXdXi(xi, 1, 1);
      dR2dt2(iv) = dX2dXi2(xi, 0, 1, 1);
      dZ2dt2(iv) = dX2dXi2(xi, 1, 1, 1);
   }

   double RBtor = getRBtoroidal();

   FORT_GET_MILLER_FIELD_DATA(CHF_BOX(box),
                              CHF_CONST_FRA1(RBpol,0),
                              CHF_CONST_FRA1(dRBpoldt,0),
                              CHF_CONST_REAL(RBtor),
                              CHF_CONST_FRA1(R,0),
                              CHF_CONST_FRA1(dRdr,0),
                              CHF_CONST_FRA1(dRdt,0),
                              CHF_CONST_FRA1(dZdr,0),
                              CHF_CONST_FRA1(dZdt,0),
                              CHF_CONST_FRA1(dR2dt2,0),
                              CHF_CONST_FRA1(dZ2dt2,0),
                              CHF_FRA(a_BField),
                              CHF_FRA1(a_BFieldMag,0),
                              CHF_FRA(a_BFieldDir),
                              CHF_FRA(a_gradBFieldMag),
                              CHF_FRA(a_curlBFieldDir),
                              CHF_FRA1(a_BFieldDirdotcurlBFieldDir,0));
}
#else

void
MillerBlockCoordSys::computeFieldData( const int  a_dir,
                                       FArrayBox& a_BField,
                                       FArrayBox& a_BFieldMag,
                                       FArrayBox& a_BFieldDir,
                                       FArrayBox& a_gradBFieldMag,
                                       FArrayBox& a_curlBFieldDir,
                                       FArrayBox& a_BFieldDirdotcurlBFieldDir,
                                       const bool a_derived_data_only ) const
{
   // This option computes the Miller field FOR AN ANNULUS ONLY from the the flux derivative dpsi_dr

   if ( m_beta != 0. || m_kappa != 1. || m_drR0 != 0. || m_s_kappa != 0. || m_s_delta != 0.) {
      MayDay::Error("MillerBlockCoordSys::computeFieldData(): Not implemented for specified geometry parameters");
   }

   // Get box intersection
   Box box = a_BField.box();
   box    &= a_BFieldMag.box();
   box    &= a_BFieldDir.box();
   box    &= a_gradBFieldMag.box();
   box    &= a_curlBFieldDir.box();
   box    &= a_BFieldDirdotcurlBFieldDir.box();

   FArrayBox RZ(box,SpaceDim);
   if (a_dir>=0 && a_dir<SpaceDim) {
     getFaceCenteredRealCoords(a_dir, RZ);
   }
   else {
     getCellCenteredRealCoords(RZ);
   }

   FArrayBox RB(box,2);
   FArrayBox dRBdR(box,2);
   FArrayBox dRBdZ(box,2);

   for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();
      
      double R_shift = RZ(iv,0) - m_origin[0];
      double Z_shift = RZ(iv,1) - m_origin[1];

      double r = sqrt( pow(R_shift,2) + pow(Z_shift,2) );
      double theta = atan2(Z_shift, R_shift);
      if (theta < 0.) theta += 2. * Pi;
      double st = sin(theta);
      double ct = cos(theta);

      RB(iv,0) = -m_dpsidr * st;  // R*B_R = -dpsi_dZ
      RB(iv,1) =  m_dpsidr * ct;  // R*B_Z =  dpsi_dR

      dRBdZ(iv,0) = -m_dpsidr * ct *ct / r;  // d(R*B_R)/dZ = -d2psi_dZ2
      dRBdR(iv,1) =  m_dpsidr * st *st / r;  // d(R*B_Z)/dR =  d2psi_dR2

      dRBdR(iv,0) = m_dpsidr * st * ct / r;  // d(R*B_R)/dR = -d2psi_dRdZ
      dRBdZ(iv,1) = -dRBdR(iv,0);            // d(R*B_Z)/dZ =  d2psi_dRdZ
   }

   double RBtor = getRBtoroidal();

   FORT_GET_FIELD_DATA( CHF_BOX(box),
                        CHF_CONST_FRA(RZ),
                        CHF_CONST_FRA(RB),
                        CHF_CONST_FRA(dRBdR),
                        CHF_CONST_FRA(dRBdZ),
                        CHF_CONST_REAL(RBtor),
                        CHF_FRA(a_BField),
                        CHF_FRA1(a_BFieldMag,0),
                        CHF_FRA(a_BFieldDir),
                        CHF_FRA(a_gradBFieldMag),
                        CHF_FRA(a_curlBFieldDir),
                        CHF_FRA1(a_BFieldDirdotcurlBFieldDir,0) );
}

#endif

Vector<Real>
MillerBlockCoordSys::computeBField(const RealVect& a_X) const
{
   
   if ( m_beta != 0. || m_kappa != 1. || m_drR0 != 0. || m_s_kappa != 0. || m_s_delta != 0.) {
      MayDay::Error("MillerBlockCoordSys::computeBField(): Not implemented for specified geometry parameters");
   }
   
   Vector<Real> result(3,0);
   
   double R_shift = a_X[0] - m_origin[0];
   double Z_shift = a_X[1] - m_origin[1];
   
   //   double r = sqrt( pow(R_shift,2) + pow(Z_shift,2) );
   double theta = atan2(Z_shift, R_shift);
   if (theta < 0.) theta += 2. * Pi;
   double st = sin(theta);
   double ct = cos(theta);
   
   result[0] = -m_dpsidr * st / a_X[0];
   result[1] =  m_Btor_scale/a_X[0];
   result[2] =  m_dpsidr * ct / a_X[0];
   
   return result;
   
}


void
MillerBlockCoordSys::getMagneticFlux( const FArrayBox& a_physical_coordinates,
                                      FArrayBox&       a_magnetic_flux ) const
{
   for (BoxIterator bit(a_physical_coordinates.box()); bit.ok(); ++bit) {
      IntVect iv = bit();
      
      RealVect coord;
      for (int n=0; n<SpaceDim; ++n) {
         coord[n] = a_physical_coordinates(iv,n);
      }

      a_magnetic_flux(iv,0) = getMagneticFlux(coord);
   }
}



double
MillerBlockCoordSys::getMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   double psi = 0.;

   // The general non-annulus case is not yet implemented.
   if ( m_beta != 0. || m_kappa != 1. || m_drR0 != 0. || m_s_kappa != 0. || m_s_delta != 0.) {
      MayDay::Error("MillerBlockCoordSys::getMagneticFlux(): Not implemented for specified geometry parameters");
   }
   else {
      double R_shift = a_physical_coordinate[0] - m_origin[0];
      double Z_shift = a_physical_coordinate[1] - m_origin[1];

      psi = m_dpsidr * sqrt( pow(R_shift,2) + pow(Z_shift,2) );
   }

   return psi;
}


void
MillerBlockCoordSys::getNodalFieldData( FArrayBox& a_points,
                                        FArrayBox& a_A,
                                        FArrayBox& a_b,
                                        FArrayBox& a_Bmag ) const
{
   if ( m_beta != 0. || m_kappa != 1. || m_drR0 != 0. || m_s_kappa != 0. || m_s_delta != 0.) {
      MayDay::Error("MillerBlockCoordSys::getNodalFieldData(): Not implemented for specified geometry parameters");
   }

   Box box = a_points.box();

   FArrayBox RZ(box,SpaceDim);
   for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();

      RealVect xi;
      for (int n=0; n<SpaceDim; ++n) {
         xi[n] = a_points(iv,n);
      }

      RealVect X = realCoord(xi);
      for (int n=0; n<SpaceDim; ++n) {
         RZ(iv,n) = X[n];
      }
   }

   FArrayBox psi(box,1);
   FArrayBox RB(box,2);

   for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();
      
      double R_shift = RZ(iv,0) - m_origin[0];
      double Z_shift = RZ(iv,1) - m_origin[1];

      double r = sqrt( pow(R_shift,2) + pow(Z_shift,2) );
      double theta = atan2(Z_shift, R_shift);
      if (theta < 0.) theta += 2. * Pi;

      psi(iv,0) = m_dpsidr * r;

      RB(iv,0) = -m_dpsidr * sin(theta);  // R*B_R = -dpsi_dZ
      RB(iv,1) =  m_dpsidr * cos(theta);  // R*B_Z =  dpsi_dR
   }

   FORT_GET_NODAL_FIELD_DATA( CHF_BOX(box),
                              CHF_CONST_FRA(RZ),
                              CHF_CONST_FRA1(psi,0),
                              CHF_CONST_FRA(RB),
                              CHF_CONST_REAL(m_Btor_scale),
                              CHF_FRA(a_A),
                              CHF_FRA(a_b),
                              CHF_FRA(a_Bmag));
}

void MillerBlockCoordSys::getNodeCenteredMappedCoords(FArrayBox& a_xi) const
{
   double rmid = getAvMinorRad() - m_inner_radial_boundary;
   int const_minor_radius = getConstMinorrad();

   FORT_GET_MILLER_NC_MAPPED_COORDS(CHF_BOX(a_xi.box()),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_CONST_REAL(rmid),
                                    CHF_CONST_INT(const_minor_radius),
                                    CHF_FRA(a_xi));
}



void MillerBlockCoordSys::getCellCenteredMappedCoords(FArrayBox& a_xi) const
{
   double rmid = getAvMinorRad() - m_inner_radial_boundary;
   int const_minor_radius = getConstMinorrad();

   FORT_GET_MILLER_CC_MAPPED_COORDS(CHF_BOX(a_xi.box()),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_CONST_REAL(rmid),
                                    CHF_CONST_INT(const_minor_radius),
                                    CHF_FRA(a_xi));
}



void MillerBlockCoordSys::getFaceCenteredMappedCoords(const int a_dir,
                                                      FArrayBox&  a_xi) const
{
   double rmid = getAvMinorRad() - m_inner_radial_boundary;
   int const_minor_radius = getConstMinorrad();

   FORT_GET_MILLER_FC_MAPPED_COORDS(CHF_BOX(a_xi.box()),
                                    CHF_CONST_INT(a_dir),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_CONST_REAL(rmid),
                                    CHF_CONST_INT(const_minor_radius),
                                    CHF_FRA(a_xi));
}



#include "NamespaceFooter.H"
