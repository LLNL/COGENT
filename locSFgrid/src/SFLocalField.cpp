#include <fstream>
#include <float.h>
#include "SFLocalField.H"
#include "BoxIterator.H"
#include "parstream.H"

#include "NamespaceHeader.H"

enum SingleNullBlockType {RCORE,LCORE,LSOL,LDIVLEG,LPFLEG1,LPFLEG2,LPFLEG3,RPFLEG3,RPFLEG2,RPFLEG1,RDIVLEG,RSOL,NUM_SNOWFLAKE_BLOCKS};
enum field_vars {BR_VAR, BZ_VAR, NUM_FIELD_VARS};
enum rz_vars {R_VAR, Z_VAR};
enum directions {RADIAL_DIR, POLOIDAL_DIR};

SFLocalField::SFLocalField(const double& a_Rmaj, const double& a_zmax, const double& a_zm0,
                     const string& a_block_name)
{
   if ( a_block_name == "lcore" ) {
      m_block_number = LCORE;
   }
   else if ( a_block_name == "rcore" ) {
      m_block_number = RCORE;
   }
   else if ( a_block_name == "lsol" ) {
      m_block_number = LSOL;
   }
   else if ( a_block_name == "ldivleg" ) {
      m_block_number = LDIVLEG;
   }
   else if ( a_block_name == "lpfleg1" ) {
      m_block_number = LPFLEG1;
   }
   else if ( a_block_name == "lpfleg2" ) {
      m_block_number = LPFLEG2;
   }
   else if ( a_block_name == "lpfleg3" ) {
      m_block_number = LPFLEG3;
   }
   else if ( a_block_name == "rpfleg3" ) {
      m_block_number = RPFLEG3;
   }
   else if ( a_block_name == "rpfleg2" ) {
      m_block_number = RPFLEG2;
   }
   else if ( a_block_name == "rpfleg1" ) {
      m_block_number = RPFLEG1;
   }
   else if ( a_block_name == "rdivleg" ) {
      m_block_number = RDIVLEG;
   }
   else if ( a_block_name == "rsol" ) {
      m_block_number = RSOL;
   }
   else {
      MayDay::Error("SFLocalField: Unrecognized block name");
   }

   m_zmax = a_zmax;
   m_zm0 = a_zm0;
   m_Rmaj = a_Rmaj;

   setIdealSFCoefs();  // Sets the coefficients of the cubic for an ideal snowflake

   /* Get maximum value of Phi */
   double zm0sq = m_zm0*m_zm0;
   m_phimax = m_l2*m_zm0 + m_q3*zm0sq + m_c4*m_zm0*zm0sq;
   m_beta = m_phimax/(m_zmax*m_zmax*m_zmax);

   third = 1.0/3.0;
   
   // print out inputs and computed variables
   pout() << "zmax = " << m_zmax << endl;
   pout() << "zm0 = " << m_zm0 << endl;
   pout() << "Rmaj = " << m_Rmaj << endl;
   pout() << "phimax = " << m_phimax << endl;
   pout() << "beta = " << m_beta << endl;
   pout() << "l1 = " << m_l1 << endl;
   pout() << "l2 = " << m_l2 << endl;
   pout() << "q1 = " << m_q1 << endl;
   pout() << "q2 = " << m_q2 << endl;
   pout() << "q3 = " << m_q3 << endl;
   pout() << "c1 = " << m_c1 << endl;
   pout() << "c2 = " << m_c2 << endl;
   pout() << "c3 = " << m_c3 << endl;
   pout() << "c4 = " << m_c4 << endl;
   pout() << "block_number = " << m_block_number << endl;
   pout() << "sin, cos rot angle = " << m_sinth[m_block_number] << ", " << m_costh[m_block_number] << endl;

}

SFLocalField::~SFLocalField()
{
}

void
SFLocalField::setIdealSFCoefs()   // Sets coefficients for ideal snowflake: coefs of cubic, and cos/sin of rotation angles
{
  /* Set sin and cos of angles of rotation for transformation from block orientation to upright. These 
     angles are, respectively, 0, 0, -pi/3, -pi/3, -2pi/3,-2pi/3, pi, pi, -2pi/3, -2pi/3   */
  double rut3o2 = sqrt(3.)/2.;
  /*  const double cosarr [] = {1.0, 1.0, 0.5, 0.5, -0.5, -0.5, -1.0, -1.0, -0.5, -0.5, 0.5, 0.5};
  const double sinarr [] = {0.0, 0.0, -rut3o2, -rut3o2,-rut3o2, -rut3o2, 0.0, 0.0, rut3o2, rut3o2, rut3o2, rut3o2};
  The above are the actual sin, cos's for the rotatin angles, but there is also reflection for the 
  even-numbered block numbers, so absorb reflection coefficients in cosarr, sinarr */
  const double cosarr [] = {1.0, -1.0, 0.5, -0.5, -0.5, 0.5, -1.0, 1.0, -0.5, 0.5, 0.5, -0.5};
  const double sinarr [] = {0.0, 0.0, -rut3o2, rut3o2,-rut3o2, rut3o2, 0.0, 0.0, rut3o2, -rut3o2, rut3o2, -rut3o2};
  memcpy (m_costh, cosarr, 12);
  memcpy (m_sinth, sinarr, 12);

  /* Now assign values for the coefficients for the cubic expansion for ideal snowflake */  
  /* Exact snowflake has 2nd order null, so m_l's and m_q's are 0 */
  m_l1 = 0.0;
  m_l2 = 0.0;
  m_q1 = 0.0;
  m_q2 = 0.0;
  m_q3 = 0.0;

  /* Remaining coefficients: per Ryutov PPCF 2010, in vacc approx Phi satisfies Laplace's eq, which implies
     that m_c3 = -3 m_c1, and m_c2 = -3 m_c4.   The condition that the ideal snowflake is oriented so as to have
     x-reflection symmetry (z axis bisects the "core" region) is that terms odd in x vanish, c1=c3=0.)
     SO we only set m_c4, and Phi = c4*(-3x^2 z + z^3).   Without loss of generality we can take m_c4=1. */
  m_c1 = 0.0;
  m_c3 = 0.0;
  m_c4 = 1.0;
  m_c2 = -3.0;
}

RealVect
SFLocalField::rotateToVert( const RealVect& a_X) const
{
  /* Rotate the vector a_X from its position in its block to where it would be in block 0 */
  RealVect X_rot;
  X_rot[0] = a_X[0]*m_costh[m_block_number] - a_X[1]*m_sinth[m_block_number];
  X_rot[1] = a_X[1]*m_costh[m_block_number] - a_X[0]*m_sinth[m_block_number];
  return X_rot;
}

RealVect
SFLocalField::rotateFromVert( const RealVect& a_X) const
{
  /* Rotate the vector a_X to its position in its block from where it would be in block 0 */
  RealVect X_rot;
  X_rot[0] = a_X[0]*m_costh[m_block_number] + a_X[1]*m_sinth[m_block_number];
  X_rot[1] = a_X[1]*m_costh[m_block_number] + a_X[0]*m_sinth[m_block_number];
  return X_rot;
}



RealVect
SFLocalField::getMappedCoordinate( const RealVect& a_X ) const
{
   /*
     Construct mapped coordinates (Ryutov's phi = cubic, and x' which is abs value of x after
     reflection and/or rotation to position of outer core)
   */
   RealVect mapped_x;
   RealVect Xprime;
   /* For now the coding for Xprime only correct for properly oriented ideal snowflake; generalize later */
   Xprime = rotateToVert( a_X );
   double xx = Xprime[R_VAR];
   double zz = Xprime[Z_VAR];
   double xsq = xx*xx;
   double zsq = zz*zz;
   /* The following will be general once we code general map to vertical orientation */
   double phi = m_l1*xx + m_l2*zz + m_q1*xsq + 2.0*m_q2*xx*zz + m_q3*zsq + m_c1*xx*xsq + m_c2*xsq*zz + m_c3*xx*zsq * m_c4*zz*zsq;
   mapped_x[RADIAL_DIR] = pow(phi/m_phimax,third);
   mapped_x[POLOIDAL_DIR] = xx/(xx+m_zmax-zz);
   return mapped_x;
}


void
SFLocalField::getPhysicalCoordinates(const FArrayBox& a_xi,
                                  FArrayBox&       a_physical_coordinates) const
{
   // Fills the components of a_physical_coordinates, evaluated at the (initial guess) mapped coordinates a_xi.
   const Box& box(a_xi.box());
   CH_assert(a_physical_coordinates.box().contains(box));
   CH_assert(a_physical_coordinates.nComp() == 2);
   const double zm0cubed = m_zm0*m_zm0*m_zm0;
   const double zmaxcubed = m_zmax*m_zmax*m_zmax;
   const double rut3m1 = sqrt(3.)-1.;
   const double rut3 = sqrt(3.);

   // loop over box cells, one by one, first finding x, z as if in block 0, then transforming to the proper box
   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();
      RealVect XC;    // Computational coordinate for a cell
      for (int dir=0; dir<SpaceDim; ++dir) {
         XC[dir] = a_xi(iv,dir);
      }
      const double eta = XC[RADIAL_DIR];
      const double xi = XC[POLOIDAL_DIR];
      double zz;
      double xx;

      /* For now specialize coding to ideal snowflake, normalized flux eta = (-3x^2 z + z^3)^(1/3)/zm0 
	and 2nd computational variable xi = x/(x+zmax-z)
	Substitute x = alpha*(zmax-z) in expression for eta, gives cubic for y = z/zmax:
	y^3(1 -  3 alpha^2) + 6 alpha^2 y^2 - 3 alpha^2 y - eta^3 beta = 0, where alpha = xi/(1-xi) and beta = zm0^2/zmax^2. */

      const double etacubed = eta*eta*eta;
      /* First do specialzed cases for 3 of the boundaries: */
      if (abs(1.0 - xi) < 1e-6) {
        zz = m_zmax;
        xx = sqrt((zmaxcubed - etacubed*zm0cubed)/(3.*zz));
      }
      else if (abs(xi) < 1e-6){
        xx = 0.0;
        zz = m_zm0*eta;
      }
      else if (abs(eta) < 1e-6){
        xx = xi*m_zmax/(1+rut3m1*xi);
        zz = rut3*xx;
      }
      else{          // set up the cubic
        const double alpha = xi/(1-xi);
        const double alphasq = alpha*alpha;
        const double a = 1.-3.*alphasq;
        const double b = 6.*alphasq;
        const double c = -3.*alphasq;
        const double d = -etacubed*m_beta;
	double guess;
	// for initial guess, note that d is typically very small while c is negative.  So a good initial guess is
	if (abs(a) > 1.e-6)
	  {guess = (-b + sqrt(b*b - 4.*a*c))/(2.*a);
	  }
	else 
	  {
	    guess = -c/b;
	  }
	
        zz  = m_zmax*newtonSolveCubic(a,b,c,d,guess);
        xx = alpha*(m_zmax-zz);
	// pout() << "  eta,xi,a,b,c,d,guess = " << eta << ", " << xi << ", " << a << ", " << b << ", " << c << ", " << d << ", " << guess << endl;
        // pout() << "eta,xi,a,b,c,d,guess,xx,zz = " << eta << ", " << xi << ", " << a << ", " << b << ", " << c << ", " << d << ", " << guess << ", " << xx << ", " << zz << endl;
      }

      RealVect XPV;
      XPV[R_VAR] = xx;
      XPV[Z_VAR] = zz;
      const RealVect XV = rotateFromVert(XPV);

      for (int n=0; n<2; ++n) {
	a_physical_coordinates(iv,n) = XV[n];
      }
   }
}


double
SFLocalField::newtonSolveCubic(const double a_a, const double a_b, const double a_c, 
			       const double a_d, double a_guess)  const {

  /* Root of cubic a_a x^3 + a_b*x^2 + a_c*x + a_d = 0 ; find root by Newton starting from a_guess */

   double tol = 1.e-12;         // Read from input if we find some reason to later
   int max_iter = 40;
   int num_iters = 0;
   double xold = a_guess;
   double xsq = xold*xold;
   double xnew;
   double Fprime;
   double resid;
   double F = a_a*xold*xsq + a_b*xsq + a_c*xold + a_d;
   resid = abs(F);
   while ( resid > tol && num_iters < max_iter) {
     Fprime = 3.0*a_a*xsq + 2.*a_b*xold + a_c;
     xnew = xold - F/Fprime;
     xsq = xnew*xnew;
     F = a_a*xnew*xsq + a_b*xsq + a_c*xnew + a_d;
     resid = abs(F);
     xold = xnew;
     num_iters++;
   }
   // pout() << " resid, tol " << resid << ", " << tol << endl;
   // pout() << "a,b,c,d,xnew,F,num_iters " << a_a << ", " << a_b << ", " << a_c << ", " << a_d << ", " << xnew << ", " << F << ", " << num_iters << endl;
   // pout() << " resid, tol " << resid << ", " << tol << endl;
   return xnew;
}

void
SFLocalField::getFieldUnitVector(const FArrayBox& a_physical_coordinates,
                              FArrayBox&       a_field_unit_vector) const
{
   /*
     Fill the components of a_field_unit_vector as follows:

         Component 0: First component of b
         Component 1: Second component of b
         Component 2: first derivative wrt xi[0] of the first component of b
         Component 3: first derivative wrt xi[1] of the first component of b
         Component 4: first derivative wrt xi[0] of the second component of b
         Component 5: first derivative wrt xi[1] of the second component of b
       Don't need the remaining 6 for now.
         Component 6: second derivative wrt xi[0] of the first component of b
         Component 7: mixed second derivative wrt xi[0] and xi[1] of the first component of b
         Component 8: second derivative wrt xi[1] of the first component of b
         Component 9: second derivative wrt xi[0] of the second component of b
         Component 10: mixed second derivative wrt xi[0] and xi[1] of the second component of b
         Component 11: second derivative wrt xi[1] of the second component of b

     evaluated at the physical coordinates a_xi.

     For the snowflake model, we have explicit expressions for Bx(x,z) and Bz(x,z).  The field unit vector b is
     (Bx e_x + Bz e_z)/B, where B=sqrt(Bx^2+Bz^2).  
     So we need to find Bx, Bz, dBxdx, dBxdz, dBzdx, dBzdz, d^2Bx/dx^2, d^2Bx/dxdz, d^2Bx/dz^2,
            d^2Bz/dx^2, d^2Bz/dx/dz, d^2Bz/dz^2

   */

   const Box& box(a_physical_coordinates.box());
   CH_assert(a_field_unit_vector.box().contains(box));
   CH_assert(a_field_unit_vector.nComp() == 6);

   // Loop through boxes obtaining physical coords x and z for each, and evaluate needed derivs from analytic expressions
   FArrayBox xi(box, SpaceDim);
   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();
      RealVect X;
      for (int dir=0; dir<SpaceDim; ++dir) {
         X[dir] = a_physical_coordinates(iv,dir);
      }
      double xx = X[R_VAR];  // physical x coordinate, using xx in case of case-insenstive compile
      double zz = X[Z_VAR];
      double invRpx = 1./(m_Rmaj+xx);
      // Get field components, directly using Ryutov expressions
      double Bx = -(m_l2 + 2.0*m_q2*xx + 2.0*m_q3*zz + m_c2*xx*xx + 2.0*m_c3*xx*zz + 3.*m_c4*zz*zz)*invRpx;
      double Bz = (m_l1 + 2.0*m_q1*xx + 2.0*m_q2*zz + 3.0*m_c1*xx*xx + 2.0*m_c2*xx*zz + m_c3*zz*zz)*invRpx;
      double Bmagsq = Bx*Bx + Bz*Bz;
      double Bmag = sqrt(Bmagsq);
      if (Bmag > 1.e-6) {
	// Get derivatives of the above
        double dBxdx = -(2.0*m_q2 + 2.0*m_c2*xx + 2.0*m_c3*zz)*invRpx - Bx*invRpx;
        double dBxdz = -(2.0*m_q3 + 2.0*m_c3*xx + 6.0*m_c4*zz)*invRpx;
        double dBzdx = (2.0*m_q1 + 6.0*m_c1*xx + 2.0*m_c2*zz)*invRpx - Bz*invRpx;
        double dBzdz = (2.0*m_q2 + 2.0*m_c2*xx + 2.0*m_c3*zz)*invRpx;

	/* The following if we ever need 2nd derivs of unit vector, so far only for Dmitri's secondary expansion
	   double d2Bxdxdx = -6.0*m_c4*xx*invRpx - 2.0*dBxdx*invRpx;
	   double d2Bxdxdz = -6.0*m_c1*invRpx - dBxdz*invRpx;
	   double d2Bxdzdz = -6.0*m_c4*invRpx;
	   double d2Bzdxdx = 6.0*m_c1*invRpx - 2.0*dBzdx*invRpx;
	   double d2Bzdxdz = -6.0*m_c4*invRpx - dBzdx*invRpx;
	   double d2Bzdzdz = -6.0*m_c1*invRpx;
	*/
        double dBmagidx = -(1.0/(Bmag*Bmagsq))*(Bx*dBxdx + Bz*dBzdx);
        double dBmagidz = -(1.0/(Bmag*Bmagsq))*(Bx*dBxdz + Bz*dBzdz);

	// Combine these to get unit vector and derivatives
	// unit vector:
        a_field_unit_vector(iv,0) = Bx/Bmag;
        a_field_unit_vector(iv,1) = Bz/Bmag;
        a_field_unit_vector(iv,2) = dBxdx/Bmag - Bx*dBmagidx;
        a_field_unit_vector(iv,3) = dBxdz/Bmag - Bx*dBmagidz;
        a_field_unit_vector(iv,4) = dBzdx/Bmag - Bz*dBmagidx;
        a_field_unit_vector(iv,5) = dBzdz/Bmag - Bz*dBmagidz;
      
      }
      else {
         for (int comp=0; comp<6; ++comp) {
            a_field_unit_vector(iv,comp) = 0.0;
         }
      }
   }
}

void
SFLocalField::writePhysicalCoordinates(const FArrayBox& a_physical_coordinates) const
{
   char file_name[80];
   sprintf(file_name, "coords%d", m_block_number);
   FILE* fd = fopen(file_name, "w");

   const Box& box = a_physical_coordinates.box();

   fprintf(fd, "%d %d\n", box.size(0), box.size(1));

   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();
      fprintf(fd, "%20.12e %20.12e \n", a_physical_coordinates(iv,0), a_physical_coordinates(iv,1));
   }

   fclose(fd);
}


void
SFLocalField::writeVectors(const FArrayBox& a_physical_coordinates,
                        const FArrayBox& a_vector_data) const
{
   char file_name[80];
   sprintf(file_name, "vectors%d", m_block_number);
   FILE* fd = fopen(file_name, "w");

   const Box& box = a_physical_coordinates.box();
   CH_assert(box==a_vector_data.box());

   fprintf(fd, "%d %d 0 0\n", box.size(0), box.size(1));

   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();
      fprintf(fd, "%20.12e %20.12e %20.12e %20.12e\n", a_physical_coordinates(iv,0), a_physical_coordinates(iv,1),
              a_vector_data(iv,0), a_vector_data(iv,1));
   }

   fclose(fd);
}



#include "NamespaceFooter.H"


