#include <array>
#include <cmath>
#include "CylindricalBlockCoordSys.H"
#include "Directions.H"
#include "CylindricalBlockCoordSys.H"

#include "Directions.H"
#include "BoxIterator.H"
#include "ConstFact.H"
//#include "CONSTANTS.H"

#include "NamespaceHeader.H"

//char* CylindricalBlockCoordSys::pp_name = {"cylindrical"};

const std::string CylindricalBlockCoordSys::pp_name = "cylindrical";

CylindricalBlockCoordSys::CylindricalBlockCoordSys( ParmParse&               a_parm_parse,
                                                    const ProblemDomain&     a_domain )
   : MagBlockCoordSys(a_parm_parse)
{

   // Read the input data specific to this geometry

   // Get max values of R, Z, (and if 5D, Phi)
   a_parm_parse.get("r_min", m_rmin);
   a_parm_parse.get("r_max", m_rmax);

   a_parm_parse.get("z_min", m_zmin);
   a_parm_parse.get("z_max", m_zmax);

#if CFG_DIM==3
   a_parm_parse.get("phi_min", m_phimin);
   a_parm_parse.get("phi_max", m_phimax);
#endif

   std::string conf_output; // output message after magnetic geometry is created
   // Get magnetic field type
   // Check for "B_type" entry
   if (a_parm_parse.contains("B_type")){
     std::string B_type;
     a_parm_parse.get("B_type", B_type);

     // Constant magnetic type: B = B0
     if (B_type=="constant") {
       this->calcB = &CylindricalBlockCoordSys::calcBConst;
       if (a_parm_parse.contains("Btor_0")) {
	 a_parm_parse.get("Btor_0", m_Btor_0);
	 conf_output = "Constant B field: Btor_0 = ";
	 conf_output.append(std::to_string(m_Btor_0));
       }
       else {MayDay::Error("ERROR: No magnetic field amplitude is specidied (\"Btor_0\" is missing)");}
     }
     // Toroidal magnetic type: B = m_Btor_scale / R
     else if (B_type=="toroidal") {
       this->calcB = &CylindricalBlockCoordSys::calcBTor;
       if (a_parm_parse.contains("Btor_scale")) {
	 a_parm_parse.get("Btor_scale", m_Btor_scale);
	 conf_output = "Toroidal B field: Btor_scale = ";
	 conf_output.append(std::to_string(m_Btor_scale));
       }
       else {MayDay::Error("ERROR: No magnetic field amplitude is specidied (\"Btor_scale\" is missing)");}
     }
     // Bennett magnetic type: B = 2 B0 *r/m_r0 / (1 + r*r/m_r0/m_r0)
     else if (B_type=="bennett") {
       this->calcB = &CylindricalBlockCoordSys::calcBBennett;
       if (a_parm_parse.contains("r_0")) {a_parm_parse.get("r_0", m_r0);}
       else {MayDay::Error("ERROR: No Bennett profile radius is specidied (\"r_0\" is missing)");}
       if (a_parm_parse.contains("Btor_0")) {
	 a_parm_parse.get("Btor_0", m_Btor_0);
	 conf_output = "Bennett B field: Btor_0 = ";
	 conf_output.append(std::to_string(m_Btor_0));
	 conf_output.append("; r0 = ");
	 conf_output.append(std::to_string(m_r0));
       }
       else {MayDay::Error("ERROR: No magnetic field amplitude is specidied (\"Btor_0\" is missing)");}
     }
     // Unknown type of magnetic geometry is specified
     else {
       conf_output = "ERROR: Unknown magnetic geometry \"";
       conf_output.append(B_type);
       conf_output.append("\"");
       const char* cstr = conf_output.c_str();
       MayDay::Error(cstr);
     }
   }
   else {MayDay::Error("ERROR: No magnetic geometry type is specified (\"B_type\" is missing)");}

   // Print geometry data to stdout
   if (m_verbose && !procID()) {
      cout << "Constructing cylindrical geometry..." << endl;
      cout << "rmin = " << m_rmin << ", rmax = " << m_rmax<<endl;
      cout << "zmin = " << m_zmin << ", zmax = " << m_zmax<<endl;
#if CFG_DIM==3
      cout << "phi_min = " << m_phimin << ", phi_max = " << m_phimax<<endl;
#endif
      cout<< conf_output << endl;
   }

   IntVect dimensions = a_domain.size();
   double dr = (m_rmax - m_rmin)/(double)dimensions[0];
   double dz = (m_zmax - m_zmin)/(double)dimensions[1];
#if CFG_DIM==3
   double dphi = (m_phimax-m_phimin)/(double)dimensions[2];
   RealVect cellSpacing(dr,dz,dphi);
#endif

#if CFG_DIM==2
   RealVect cellSpacing(dr,dz);
#endif

   // Finish defining the object now that we also have the mesh spacing
   define( a_domain, cellSpacing );

   if (m_verbose && procID()==0) {
      cout << "Done constructing cylindrical geometry" << endl;
   }
}


RealVect CylindricalBlockCoordSys::realCoord( const RealVect& a_xi ) const
{
   RealVect x;
   D_TERM(x[0] = m_rmin + a_xi[0];,
          x[1] = m_zmin + a_xi[1];,
          x[2] = m_phimin + a_xi[2];)
   return x;
}

RealVect CylindricalBlockCoordSys::mappedCoord( const RealVect& a_x ) const
{

   RealVect xi;
   D_TERM(xi[0] = a_x[0] - m_rmin;,
          xi[1] = a_x[1] - m_zmin;,
          xi[2] = a_x[2] - m_phimin;)
   return xi;
}

Real CylindricalBlockCoordSys::dXdXi( const RealVect& a_Xi,
                                     int             a_dirX,
                                     int             a_dirXi ) const
{

  //Presently, consider 2D case only
  Real value = 0.0;
  if (a_dirX != a_dirXi) {
     value = 0.0;
  }
  else{
     value = 1.0;
  }

  return value;
}

void CylindricalBlockCoordSys::dXdXi(FArrayBox&       a_dXdXi,
                                     const FArrayBox& a_Xi,
                                     int              a_destComp,
                                     int              a_dirX,
                                     int              a_dirXi,
                                     const Box&       a_box) const
{
  
  //Presently, consider 2D case only
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit) {
      const IntVect& iv = bit();
      if (a_dirX != a_dirXi) {
        a_dXdXi(iv,a_destComp) = 0;
      }
      else{
        a_dXdXi(iv,a_destComp) = 1.0;
      }
   }
}

void
CylindricalBlockCoordSys::computeFieldData(const int  a_dir,
                                           FArrayBox& a_BField,
                                           FArrayBox& a_BFieldMag,
                                           FArrayBox& a_BFieldDir,
                                           FArrayBox& a_gradBFieldMag,
                                           FArrayBox& a_curlBFieldDir,
                                           FArrayBox& a_BFieldDirdotcurlBFieldDir,
                                           const bool a_derived_data_only ) const
// NEED TO FOLLOW AND SEE WHAT IS GOTTEN IN FORTRAN
{
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
      double p_B[17];	//array for magnetic filed configuration
	  
      // Calculate magnetic geometry components
      (this->*calcB)(R, Z, p_B);
	
      a_BField(iv,0) = p_B[0]; // B[0]
      a_BField(iv,1) = p_B[1]; // B[1]
      a_BField(iv,2) = p_B[2]; // B[2]
      
      a_BFieldMag(iv,0) = p_B[3]; // abs(B)
      
      a_BFieldDir(iv,0) = p_B[4]; // b[0], b=b/abs(B)
      a_BFieldDir(iv,1) = p_B[5]; // b[1]
      a_BFieldDir(iv,2) = p_B[6]; // b[2]

      a_gradBFieldMag(iv,0) = p_B[7]; // gB[0], gB = \nabla abs(B)
      a_gradBFieldMag(iv,1) = p_B[8]; // gB[1]
      a_gradBFieldMag(iv,2) = p_B[9]; // gB[2]
      
      a_curlBFieldDir(iv,0) = p_B[10]; // cb[1], cb = \nabla \times b
      a_curlBFieldDir(iv,1) = p_B[11]; // cb[2]
      a_curlBFieldDir(iv,2) = p_B[12]; // cb[3]
      
      a_BFieldDirdotcurlBFieldDir(iv,0) = p_B[13]; // b \cdot \nabla \times b
   }
}

void
CylindricalBlockCoordSys::calcBConst( double a_R, double a_Z, double* a_p_B ) const
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

void
CylindricalBlockCoordSys::calcBTor( double a_R, double a_Z, double* a_p_B ) const
{
  double sign_b; if (m_Btor_scale>0) {sign_b = 1.0;} else {sign_b = -1.0;}
  a_p_B[0] = 0.0;
  a_p_B[1] = m_Btor_scale / a_R;
  a_p_B[2] = 0.0;
      
  a_p_B[3] = fabs( a_p_B[1] );
      
  a_p_B[4] = 0.0;
  a_p_B[5] = sign_b;
  a_p_B[6] = 0.0;

  a_p_B[7] = -a_p_B[1] / a_R;
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

void
CylindricalBlockCoordSys::calcBBennett( double a_R, double a_Z, double* a_p_B ) const
{
	
  double sign_b; if (m_Btor_0>0) {sign_b = 1.0;} else {sign_b = -1.0;}
  double xi = a_R / m_r0;
  double B_tmp = 2.0 * xi / (1.0 + xi * xi);
	
  a_p_B[0] = 0.0;
  a_p_B[1] = m_Btor_0 * B_tmp;
  a_p_B[2] = 0.0;
      
  a_p_B[3] = fabs( a_p_B[1] );
      
  a_p_B[4] = 0.0;
  a_p_B[5] = sign_b;
  a_p_B[6] = 0.0;

  a_p_B[7] = a_p_B[1] * (1.0 - xi * B_tmp) / a_R;
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
CylindricalBlockCoordSys::computeBField(const RealVect& a_X) const
{
   array<double,3> result;
   
   double p_B[17];
   (this->*calcB)(a_X[0], a_X[1], p_B);
   
   result[0] =  p_B[0];
   result[1] =  p_B[1];
   result[2] =  p_B[2];
   
   
   return result;
}


void CylindricalBlockCoordSys::getMagneticFlux( const FArrayBox& physical_coordinates,
                                                FArrayBox&       magnetic_flux ) const
{
   MayDay::Error("CylindricalBlockCoordSys::getMagneticFlux() not implemented");
}



double CylindricalBlockCoordSys::getMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   MayDay::Error("CylindricalBlockCoordSys::getMagneticFlux() not implemented");
   return(0.0);
}


void CylindricalBlockCoordSys::getNodalFieldData(FArrayBox& points,
                                                 FArrayBox& A,
                                                 FArrayBox& b,
                                                 FArrayBox& Bmag) const
{
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

double CylindricalBlockCoordSys:: getRBtoroidal() const
{
  if( this->calcB==&CylindricalBlockCoordSys::calcBTor ) {return m_Btor_scale;}
  else {
    MayDay::Error("ERROR: Magnetic geometry is not toroidal; divfree_velocity is defined for toroidal geometry only");
    return 0.0;
  }
}

#include "NamespaceFooter.H"
