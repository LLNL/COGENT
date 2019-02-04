#if CFG_DIM==3
#include "ToroidalBlockCoordSys.H"
#include "MagBlockCoordSysF_F.H"

#include "Directions.H"
#include "BoxIterator.H"
//#include "ConstFact.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

const std::string ToroidalBlockCoordSys::pp_name = "toroidal";

ToroidalBlockCoordSys::ToroidalBlockCoordSys( ParmParse&               a_parm_parse,
                                              const ProblemDomain&     a_domain,
                                              const int                a_block,
                                              const int                a_num_blocks,
                                              const int                a_block_separation)
   : MagBlockCoordSys(a_parm_parse),
     m_block(a_block),
     m_num_blocks(a_num_blocks),
     m_is_flux_defined(false),
     m_eps(1.0e-10)
{
   if (SpaceDim != 3){
      MayDay::Error("ToroidalBlockCoordSys:: 3D space is expected; for axisymmetric 2D case use Miller");
   }
   
   // Read the input data specific to this geometry
   parseParameters (a_parm_parse);
   
   // Compute the mesh size in computational coordinates
   // radial (minor radius) coordinate runs from rmin to rmax
   // NB: the mapped radial coordinate runs from 0 to rmax-rmin
   // poloidal coordinate runs from 0 to  2*pi
   // toroidal coordinate runs for 0 to phi_max
   // for num_blocks>1 the block-normal mapped coordinate goes
   // from a_block*a_block_seprattion + phi_max/num_blocks * m_block
   // to a_block*a_block_sepration + phi_max/num_blocks * (m_block + 1)
   

   double twopi = 2. * Pi;
   m_mapped_block_size[RADIAL_DIR] = m_rmax - m_rmin;
   m_mapped_block_size[POLOIDAL_DIR] = twopi;
   m_mapped_block_size[TOROIDAL_DIR] = twopi * m_phi_frac/(double)m_num_blocks;
   
   IntVect dimensions = a_domain.size();
   RealVect dx;
   for (int dir=0; dir<SpaceDim; ++dir) {
      dx[dir] = m_mapped_block_size[dir]/(double)dimensions[dir];
   }

   //This defines the block_coord_sys object (e.g., sets m_dx)
   define( a_domain, dx );

   //Physical coordinate of a block center in the multiblock dir
   //here we assume that the block size in the toroidal direction
   //is the same for the mapped and physical spaces
   m_mb_center = m_mapped_block_size[TOROIDAL_DIR] * ((double)m_block + 0.5);
   
   //Block separation measured in mapped space
   m_separation = dx[TOROIDAL_DIR] * a_block_separation;
   
   // Print geometry data
   if (m_verbose) {
      printParameters();
   }
}


RealVect ToroidalBlockCoordSys::realCoord( const RealVect& a_xi ) const
{
   RealVect x;
  
   double r_norm = (m_rmin + a_xi[RADIAL_DIR]) / m_R0;
   double phi = a_xi[TOROIDAL_DIR] - m_separation * m_block;
   double theta = a_xi[POLOIDAL_DIR];
   
   if (m_field_aligned) {
     theta = computeFieldAlignedTheta(r_norm, phi, a_xi[POLOIDAL_DIR]);
   }
   
   x[RADIAL_DIR] = m_R0*(1.0 + r_norm*cos(theta))*cos(phi);
   x[TOROIDAL_DIR] = m_R0*(1.0 + r_norm*cos(theta))*sin(phi);
   x[POLOIDAL_DIR] = m_R0*r_norm*sin(theta);
   
   return x;
}

RealVect ToroidalBlockCoordSys::mappedCoord( const RealVect& a_x ) const
{

   RealVect xi;
   
   double x = a_x[RADIAL_DIR];
   double y = a_x[TOROIDAL_DIR];
   double z = a_x[POLOIDAL_DIR];
   
   double r_norm = sqrt(pow(sqrt(x*x+y*y)-m_R0,2)+z*z)/m_R0;
   double phi = computePhi(x,y);
   double r_xy_norm = sqrt(x*x + y*y)/m_R0;
   double z_norm = z/m_R0;
   double theta = computeTheta(z_norm, r_norm, r_xy_norm, phi);
   
   xi[RADIAL_DIR] = r_norm * m_R0 - m_rmin;
   xi[TOROIDAL_DIR] = phi + m_separation * m_block;
   xi[POLOIDAL_DIR] = theta;
   
   return xi;
}

Real ToroidalBlockCoordSys::dXdXi( const RealVect& a_Xi,
                                   int             a_dirX,
                                   int             a_dirXi ) const
{
  Real value = 0.0;
  
  double r = a_Xi[RADIAL_DIR] + m_rmin;
  double phi = a_Xi[TOROIDAL_DIR] - m_separation * m_block;
  double theta;
  if (!m_field_aligned) {
    theta = a_Xi[POLOIDAL_DIR];
  }
  else {
    double rnorm = r / m_R0;
    theta = computeFieldAlignedTheta(rnorm, phi, a_Xi[POLOIDAL_DIR]);
  }

  if (a_dirX == RADIAL_DIR) {
    if (a_dirXi == RADIAL_DIR) {
      value = cos(theta)*cos(phi);
    }
    else if (a_dirXi == TOROIDAL_DIR) {
      value = -(m_R0 + r*cos(theta))*sin(phi);
    }
    else if (a_dirXi == POLOIDAL_DIR) {
      value = -r*sin(theta)*cos(phi);
    }
  }
  
  if (a_dirX == TOROIDAL_DIR) {
    if (a_dirXi == RADIAL_DIR) {
      value = cos(theta)*sin(phi);
    }
    else if (a_dirXi == TOROIDAL_DIR) {
      value = (m_R0 + r*cos(theta))*cos(phi);
    }
    else if (a_dirXi == POLOIDAL_DIR) {
      value = -r*sin(theta)*sin(phi);
    }
  }

  if (a_dirX == POLOIDAL_DIR) {
    if (a_dirXi == RADIAL_DIR) {
      value = sin(theta);
    }
    else if (a_dirXi == TOROIDAL_DIR) {
      value = 0.0;
    }
    else if (a_dirXi == POLOIDAL_DIR) {
      value = r*cos(theta);
    }
  }

  if (m_field_aligned) {
    
    double rnorm = r / m_R0;
    RealVect dtheta_dXi = getThetaDerivatives(rnorm, phi, a_Xi[POLOIDAL_DIR]);
    
    if (a_dirX == RADIAL_DIR) {
      if (a_dirXi == RADIAL_DIR) {
        value -= r*sin(theta)*dtheta_dXi[RADIAL_DIR]*cos(phi);
      }
      else if (a_dirXi == TOROIDAL_DIR) {
        value -= r*sin(theta)*dtheta_dXi[TOROIDAL_DIR]*cos(phi);
      }
      else if (a_dirXi == POLOIDAL_DIR) {
        value *= dtheta_dXi[POLOIDAL_DIR];
      }
    }
    
    if (a_dirX == TOROIDAL_DIR) {
      if (a_dirXi == RADIAL_DIR) {
        value -= r*sin(theta)*dtheta_dXi[RADIAL_DIR]*sin(phi);
      }
      else if (a_dirXi == TOROIDAL_DIR) {
        value -= r*sin(theta)*dtheta_dXi[TOROIDAL_DIR]*sin(phi);
      }
      else if (a_dirXi == POLOIDAL_DIR) {
        value *= dtheta_dXi[POLOIDAL_DIR];
      }
    }
    
    if (a_dirX == POLOIDAL_DIR) {
      if (a_dirXi == RADIAL_DIR) {
        value += r*cos(theta)*dtheta_dXi[RADIAL_DIR];
      }
      else if (a_dirXi == TOROIDAL_DIR) {
        value = r*cos(theta)*dtheta_dXi[TOROIDAL_DIR];
      }
      else if (a_dirXi == POLOIDAL_DIR) {
        value *= dtheta_dXi[POLOIDAL_DIR];
      }
    }
  }
  
  return value;
}

void ToroidalBlockCoordSys::dXdXi(FArrayBox&       a_dXdXi,
                              const FArrayBox& a_Xi,
                              int              a_destComp,
                              int              a_dirX,
                              int              a_dirXi,
                              const Box&       a_box) const
{
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

void
ToroidalBlockCoordSys::computeFieldData(const int  a_dir,
                                        FArrayBox& a_BField,
                                        FArrayBox& a_BFieldMag,
                                        FArrayBox& a_BFieldDir,
                                        FArrayBox& a_gradBFieldMag,
                                        FArrayBox& a_curlBFieldDir,
                                        FArrayBox& a_BFieldDirdotcurlBFieldDir,
                                        const bool a_derived_data_only ) const
{
  // This option computes the magnetic field from the flux derivative dpsi_dr
  
  // Get box intersection
  Box box = a_BField.box();
  box    &= a_BFieldMag.box();
  box    &= a_BFieldDir.box();
  box    &= a_gradBFieldMag.box();
  box    &= a_curlBFieldDir.box();
  box    &= a_BFieldDirdotcurlBFieldDir.box();
  
  FArrayBox X(box,SpaceDim);
  if (a_dir>=0 && a_dir<SpaceDim) {
    getFaceCenteredRealCoords(a_dir, X);
  }
  else {
    getCellCenteredRealCoords(X);
  }
  
  FArrayBox RB(box,2);
  FArrayBox dRBdR(box,2);
  FArrayBox dRBdZ(box,2);
  
  FArrayBox RZ(box,SpaceDim-1);

  //compute Cylindrical vector components 
  //as done for axisymmetric (2D) geometries 
  for (BoxIterator bit(box); bit.ok(); ++bit) {
    IntVect iv = bit();
    
    RZ(iv,0) = sqrt(pow(X(iv,0),2)+pow(X(iv,1),2));
    RZ(iv,1) = X(iv,2);
    
    double R_shift = RZ(iv,0) - m_R0;
    double Z_shift = RZ(iv,1);
    
    double r = sqrt( pow(R_shift,2) + pow(Z_shift,2) );
    double theta = atan2(Z_shift, R_shift);
    if (theta < 0.) theta += 2. * Pi;
    double st = sin(theta);
    double ct = cos(theta);
    
    RB(iv,0) = -dpsidr(r) * st;  // R*B_R = -dpsi_dZ
    RB(iv,1) =  dpsidr(r) * ct;  // R*B_Z =  dpsi_dR
    
    dRBdZ(iv,0) = -dpsidr(r) * ct *ct / r - d2psidr2(r) * st * st;   // d(R*B_R)/dZ = -d2psi_dZ2
    dRBdR(iv,1) =  dpsidr(r) * st *st / r + d2psidr2(r) * ct * ct;   // d(R*B_Z)/dR =  d2psi_dR2
    
    dRBdR(iv,0) = dpsidr(r) * st * ct / r - d2psidr2(r) * st * ct;   // d(R*B_R)/dR = -d2psi_dRdZ
    dRBdZ(iv,1) = -dRBdR(iv,0);                                      // d(R*B_Z)/dZ =  d2psi_dRdZ
  }
  
  double RBtor = getRBtoroidal();
  
  FORT_GET_FIELD_DATA(CHF_BOX(box),
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

  //get physical (Cartesian) vector components
  convertCylindricalToCartesian(a_BField,X);
  convertCylindricalToCartesian(a_BFieldDir,X);
  convertCylindricalToCartesian(a_gradBFieldMag,X);
  convertCylindricalToCartesian(a_curlBFieldDir,X);
}

Vector<Real>
ToroidalBlockCoordSys::computeBField(const RealVect& a_X) const
{
  
  Vector<Real> result(3,0);
  
  double R_shift = sqrt(pow(a_X[0],2)+pow(a_X[1],2)) - m_R0;
  double Z_shift = a_X[2];
  
  double r = sqrt( pow(R_shift,2) + pow(Z_shift,2) );
  double theta = atan2(Z_shift, R_shift);
  if (theta < 0.) theta += 2. * Pi;
  double st = sin(theta);
  double ct = cos(theta);
  
  
  double R = R_shift + m_R0;
  
  //compute Cylindrical components
  result[0] = -dpsidr(r) * st / R;
  result[1] =  m_Btor_scale / R;
  result[2] =  dpsidr(r) * ct / R;

  //get physical (Cartesian) components
  RealVect tmp;
  tmp[0] = result[0];
  tmp[1] = result[1];
  tmp[2] = result[2];

  convertCylindricalToCartesian(tmp, a_X);

  result[0] = tmp[0];
  result[1] = tmp[1];
  result[2] = tmp[2];
  
  return result;
  
}


void
ToroidalBlockCoordSys::getMagneticFlux(const FArrayBox& a_physical_coordinates,
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
ToroidalBlockCoordSys::getMagneticFlux( const RealVect& a_physical_coordinate ) const
{
  double psi = 0.;
  
  RealVect X(a_physical_coordinate);
  
  double R_shift = sqrt(pow(X[0],2)+pow(X[1],2)) - m_R0;
  double Z_shift = X[2];
  
  double r = sqrt( pow(R_shift,2) + pow(Z_shift,2) );
  
  psi = getMagneticFlux(r);
  
  return psi;
}

Real
ToroidalBlockCoordSys::getMagneticFlux( const Real a_r) const
{
   
   //Set the region where magnetic flux will be defined
   Real r_in = m_rmin - 4.0 * m_dx[RADIAL_DIR];
   Real r_out = m_rmax + 4.0 * m_dx[RADIAL_DIR];

   if (!m_is_flux_defined) {
      initializeMagneticFluxFunction(r_in, r_out);
      m_is_flux_defined = true;
   }
   
   Real psi;
   
   // use analytic expression if q2=0
   if (m_q[2] < 1.0 + m_eps && m_q[2] > 1.0 - m_eps) {
      psi = (m_q[1]*a_r/m_a + m_q[0]*log(m_q[0]/(m_q[0] + m_q[1]*a_r/m_a)))/pow(m_q[1],2);
      psi *= pow(m_a,2) * m_Btor_scale / m_R0;

   }
   
   else {

      //m_magnetic_flux_inter vector contains data at
      //cell centeres and face centeres, so in priciple,
      //we don't need to add any interpolation function, since
      //we don't expect other locations in the code. However,
      //just in case we use here first-order interpolation
      //(to deal with possible round-off errors)
      
      Real dr = m_dx[RADIAL_DIR] / 2.0;
      Real r_lo_bnd = r_in + dr;
      int n = (a_r > r_lo_bnd) ?  floor((a_r - r_lo_bnd)/dr) : 0;
      
      Real r_lo = r_lo_bnd + n * dr;
      Real coeff_lo = 1.0 - (a_r - r_lo)/dr;
      Real coeff_hi = (a_r - r_lo)/dr;
      
      psi = coeff_lo * m_magnetic_flux_interp[n] + coeff_hi * m_magnetic_flux_interp[n+1];
   }
   
 
   return psi;
}

void
ToroidalBlockCoordSys::getNodalFieldData(FArrayBox& a_points,
                                         FArrayBox& a_A,
                                         FArrayBox& a_b,
                                         FArrayBox& a_Bmag ) const
{
  Box box = a_points.box();
  
  FArrayBox coord(box,SpaceDim);
  for (BoxIterator bit(box); bit.ok(); ++bit) {
    IntVect iv = bit();
    
    RealVect xi;
    for (int n=0; n<SpaceDim; ++n) {
      xi[n] = a_points(iv,n);
    }
    
    RealVect X = realCoord(xi);
    for (int n=0; n<SpaceDim; ++n) {
      coord(iv,n) = X[n];
    }
  }
  
  //Compute cylindrical field component
  //as done for axisymmetric geometries
  FArrayBox psi(box,1);
  FArrayBox RB(box,2);
  
  FArrayBox RZ(box,SpaceDim-1);
  for (BoxIterator bit(box); bit.ok(); ++bit) {
    IntVect iv = bit();
    
    RZ(iv,0) = sqrt(pow(coord(iv,0),2)+pow(coord(iv,1),2));
    RZ(iv,1) = coord(iv,2);
    
    double R_shift = RZ(iv,0) - m_R0;
    double Z_shift = RZ(iv,1);
    
    double r = sqrt( pow(R_shift,2) + pow(Z_shift,2) );
    double theta = atan2(Z_shift, R_shift);
    if (theta < 0.) theta += 2. * Pi;
    
    psi(iv,0) = getMagneticFlux(r);
    
    RB(iv,0) = -dpsidr(r) * sin(theta);  // R*B_R = -dpsi_dZ
    RB(iv,1) =  dpsidr(r) * cos(theta);  // R*B_Z =  dpsi_dR
  }
  
  FORT_GET_NODAL_FIELD_DATA( CHF_BOX(box),
                             CHF_CONST_FRA(RZ),
                             CHF_CONST_FRA1(psi,0),
                             CHF_CONST_FRA(RB),
                             CHF_CONST_REAL(m_Btor_scale),
                             CHF_FRA(a_A),
                             CHF_FRA(a_b),
                             CHF_FRA(a_Bmag));

  //Convert cylindrical vector components to Cartisian (physical) 
  convertCylindricalToCartesian(a_A,coord);
  convertCylindricalToCartesian(a_b,coord);
}


Real ToroidalBlockCoordSys::getSafetyFactor( const Real a_rnorm ) const
{
  Real q = m_q[0] + m_q[1]*exp(m_q[2]*log(a_rnorm * m_R0/m_a));
  return q;
}


Real ToroidalBlockCoordSys::getSafetyFactorDerivative( const Real a_rnorm ) const
{
  Real dqdr = m_q[1] * exp(m_q[2]*log(a_rnorm * m_R0/m_a)) * m_q[2] / a_rnorm;
  return dqdr;
}

Real ToroidalBlockCoordSys::dpsidr( const Real a_r ) const
{
  Real rnorm = a_r / m_R0;
  return m_Btor_scale * rnorm / getSafetyFactor(rnorm);
}

Real ToroidalBlockCoordSys::d2psidr2( const Real a_r ) const
{
   Real rnorm = a_r / m_R0;
   Real q = getSafetyFactor(rnorm);
   
   Real result = 1.0/q;
   
   result -= rnorm*getSafetyFactorDerivative(rnorm)/pow(q,2);
   
   result *= m_Btor_scale / m_R0;
   
   return result;
}

double ToroidalBlockCoordSys::computePhi( const Real a_x,
                                          const Real a_y ) const
{
   double r0 = sqrt(a_x*a_x + a_y*a_y);
   double phi = acos(a_x/r0);
   if (a_y<0) phi = 2. * Pi - phi;
   return phi;
}

double ToroidalBlockCoordSys::computeTheta( const Real a_z,
                                            const Real a_r,
                                            const Real a_rxy,
                                            const Real a_phi) const
{
   double theta = asin(a_z/a_r);
   if (a_rxy < 1.0) theta = Pi - theta;
   if (theta < 0) theta += 2. * Pi;

   //slight shift to avoid singularities
   if (theta == Pi) theta -= m_eps;
  
   if (m_field_aligned) {
      
      double q = getSafetyFactor(a_r);
      
      //return theta limits to {-Pi,Pi} so that
      //atan(theta/2) is uniquely defined
      if (theta > Pi) theta -= 2. * Pi;
      double angle = atan(sqrt(1.0-a_r*a_r)/(1.0+a_r)*tan(theta/2.0));
      theta = 2.0 * atan((1.0+a_r)/sqrt(1.0-a_r*a_r)*tan(angle - sqrt(1.0-a_r*a_r)*(a_phi-m_mb_center)/(2.0*q)));
      
      //restore {0,2Pi} theta limits
      if (theta < 0) theta += 2. * Pi;
      
   }
   
   return theta;
}

double ToroidalBlockCoordSys::computeFieldAlignedTheta(double a_r,
                                                       double a_phi,
                                                       double a_xi2) const
{
  
  double theta = a_xi2;

  //slight shift to avoid singularities
  if (theta == Pi) theta -= m_eps;
   
  //return xi2 limits to {-Pi,Pi} so that
  //atan(theta/2) is uniquely defined
  if (a_xi2 > Pi) theta -= 2. * Pi;

  double q = getSafetyFactor(a_r);
  
  double angle = 0.5*sqrt(1.0-a_r*a_r)*(a_phi-m_mb_center)/q;
  angle += atan(sqrt(1.0-a_r*a_r)/(1.0 + a_r)*tan(theta/2.0));
  
  theta = 2.0 * atan((1.0 + a_r)/sqrt(1.0-a_r*a_r)*tan(angle));

  return theta;
}


RealVect ToroidalBlockCoordSys::getThetaDerivatives(double a_r,
                                                    double a_phi,
                                                    double a_xi2) const
{
  RealVect result;
  
  double theta = a_xi2;
  double phi = a_phi - m_mb_center;
  double r = a_r;

  //slight shift to avoid singularities
  if (theta == Pi) theta -= m_eps;
  
  //return xi2 limits to {-Pi,Pi} so that
  //atan(theta/2) is uniquely defined
  if (a_xi2 > Pi) theta -= 2. * Pi;
  
  double q = getSafetyFactor(r);
  double dqdr = getSafetyFactorDerivative(r);
  
  double angle = 0.5*sqrt(1.0-r*r)*phi/q;
  angle += atan(sqrt(1.0-r*r)/(1.0 + r)*tan(theta/2.0));
  
  //get dtheta/dr
  double numerator = -q*phi*r*(1.0+r*cos(a_xi2));
  numerator -= q*q*sin(a_xi2);
  numerator += dqdr*phi*(r*r-1.0)*(1.0+r*cos(a_xi2));

  numerator *= -1.0*(r-1.0)*pow(r+1,2)/pow(cos(angle),2);
  numerator /= sqrt(1.0-r*r)*q*q*(1.0+r*cos(a_xi2));
  numerator += 2.0*(1.0+r)*tan(angle);
  
  double denominator = pow(1.0-r*r,(3.0/2.0))*(1.0+(1.0+r)/(1.0-r)*pow(tan(angle),2));
  
  result[RADIAL_DIR] = (1/m_R0)*numerator/denominator;
  
  
  //get dtheta/dphi
  result[TOROIDAL_DIR] = (1.0-r*r)/(1.0-r*cos(2.0*angle))/q;

  //get dtheta/dxi2
  result[POLOIDAL_DIR] = (1.0-r*r)/(1.0-r*cos(2.0*angle))/(1.0+r*cos(a_xi2));
  
  return result;
  
}

void ToroidalBlockCoordSys::initializeMagneticFluxFunction(const Real& a_rmin,
                                                           const Real& a_rmax) const
{
   //normalization factor
   Real fac = m_Btor_scale / m_R0;
   
   // refine factor to define a step
   // used in numerical integration
   // must be one over a power of two
   int refine_factor = 4;
   
   //compute psi at the inner boundary
   Real dr = a_rmin / round(a_rmin / m_dx[RADIAL_DIR]);
   dr /= double(refine_factor);
   Real r = 0.5*dr;
   Real psi_shift = 0.0;
   while ( r < a_rmin ) {
      Real r_norm = r / m_R0;
      psi_shift += fac * r / getSafetyFactor(r_norm) * dr;
      r += dr;
   }
   
   //create magnetic flux object to hold data
   int Npts = int(round((a_rmax - a_rmin) / m_dx[RADIAL_DIR]));
   Npts *= 2;
   m_magnetic_flux_interp.resize(Npts);

   //define step for numerical integration
   dr = m_dx[RADIAL_DIR]/double(refine_factor);
   r = a_rmin + 0.5*dr;
   Real psi = psi_shift;
   for (int n=0; n<Npts; ++n) {
      for (int k=0; k<refine_factor/2; ++k) {
         Real r_norm = r / m_R0;
         psi += fac * r / getSafetyFactor(r_norm) * dr;
         r += dr;
      }
      m_magnetic_flux_interp[n] = psi;
   }
}

void ToroidalBlockCoordSys::convertCartesianToToroidal(RealVect& a_vect,
                                                       const RealVect& a_coord) const
{
  //get Cartisian coordinates
  Real x = a_coord[0];
  Real y = a_coord[1];
  Real z = a_coord[2];

  //get XYZ components of the physical vector
  Real F_x = a_vect[0];
  Real F_y = a_vect[1];
  Real F_z = a_vect[2];

  //compute the corresponding toroidal components 
  Real r_xy = sqrt(x*x + y*y);
  Real F_R = (F_x*x + F_y * y)/r_xy;
  Real F_phi = (-F_x * y + F_y * x)/r_xy;  

  Real r_RZ = sqrt(pow(r_xy-m_R0,2) + z*z);
  Real F_rad = (F_R * (r_xy-m_R0) + F_z * z)/r_RZ;
  Real F_theta = (-F_R * z + F_z * (r_xy-m_R0))/r_RZ;

  a_vect[RADIAL_DIR] = F_rad;
  a_vect[TOROIDAL_DIR] = F_phi;
  a_vect[POLOIDAL_DIR] = F_theta;
}

void ToroidalBlockCoordSys::convertCartesianToToroidal(RealVect& a_coord) const
{
  //get Cartisian coordinates
  Real x = a_coord[0];
  Real y = a_coord[1];
  Real z = a_coord[2];

  //compute toroidal coordinates
  Real phi = atan2(y,x);
  if (phi < 0.) phi += 2.*Pi;

  Real r = sqrt(pow(sqrt(x*x+y*y)-m_R0,2)+z*z);
  Real r_xy = sqrt(x*x + y*y);
  Real theta = asin(z/r);
  if (r_xy < m_R0) theta = Pi - theta;
  if (theta < 0) theta += 2. * Pi;

  a_coord[RADIAL_DIR] = r;
  a_coord[TOROIDAL_DIR] = phi;
  a_coord[POLOIDAL_DIR] = theta;

}

void ToroidalBlockCoordSys::convertCylindricalToCartesian(RealVect& a_vect,
                                                          const RealVect& a_coord) const
{
  //get toroidal coordinates
  Real x = a_coord[0];
  Real y = a_coord[1];
  Real r_xy = sqrt(x*x + y*y);

  //get R,Phi components of the physical vector
  Real F_R = a_vect[0];
  Real F_Phi = a_vect[1];

  //compute the corresponding X,Y components
  Real F_X = (x/r_xy) * F_R - (y/r_xy) * F_Phi;
  Real F_Y = (y/r_xy) * F_R + (x/r_xy) * F_Phi;
  
  a_vect[0] = F_X;
  a_vect[1] = F_Y;
}

void ToroidalBlockCoordSys::convertCylindricalToCartesian(FArrayBox& a_vect,
                                                          const FArrayBox& a_coord) const
{
  Box box = a_coord.box();

  RealVect vect;
  RealVect coord;

  for (BoxIterator bit(box); bit.ok(); ++bit) {
    IntVect iv = bit();

    for (int dir=0; dir<SpaceDim; ++dir) {       
      vect[dir] = a_vect(iv,dir);
      coord[dir] = a_coord(iv,dir);
    }

    convertCylindricalToCartesian(vect,coord);

    for (int dir=0; dir<SpaceDim; ++dir) {
      a_vect(iv,dir) = vect[dir];
    }
  }
}

void ToroidalBlockCoordSys::getToroidalCoords(FArrayBox& a_coords) const
{
  Box box = a_coords.box();

  FArrayBox X(box,SpaceDim);
  getCellCenteredRealCoords(X);

  RealVect coord;
  for (BoxIterator bit(box); bit.ok(); ++bit) {
    IntVect iv = bit();

    for (int dir=0; dir<SpaceDim; ++dir) {       
      coord[dir] = X(iv,dir);
    }

    convertCartesianToToroidal(coord);

    for (int dir=0; dir<SpaceDim; ++dir) {
      a_coords(iv,dir) = coord[dir];
    }
  }
}

void ToroidalBlockCoordSys::parseParameters(ParmParse& a_parm_parse)
{
   
  // Get the wedge dimensions
  a_parm_parse.get("inner_radial_bdry", m_rmin);
  a_parm_parse.get("outer_radial_bdry", m_rmax);
  a_parm_parse.get("toroidal_fraction", m_phi_frac);
  
  a_parm_parse.get("major_radius", m_R0);
  a_parm_parse.get("minor_radius", m_a);
  
  // Get the toroidal field component scale factor
  if (a_parm_parse.contains("Btor_scale")) {
    a_parm_parse.get("Btor_scale", m_Btor_scale);
  }
  else {
    m_Btor_scale = 0.;  // default
  }
  
  // Get the magnetic safety factor parameters
  m_q.resize(3);
  a_parm_parse.getarr("safety_factor_parm", m_q, 0, 3);
  
}
void ToroidalBlockCoordSys::printParameters()
{
   // Print geometry data to stdout
   if (procID()==0 && m_block == 0) {
      cout << "Toroidal geometry parameters" << endl;
      cout << "rmin = " << m_rmin << ", rmax = " << m_rmax << ", phimax = " << m_phi_frac <<endl;
      cout << "majorRadius = " << m_R0 << ", minorRadius = " << m_a  << endl;
     
      cout << "Btor_scale = " << m_Btor_scale << endl;
      cout << "safety_fact_parm = " << m_q << endl;

      cout << "Number of blocks =" << m_num_blocks << endl;
   }
}

#include "NamespaceFooter.H"

#endif
