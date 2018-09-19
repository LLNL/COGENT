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

   // Get RBphi=const
   if (a_parm_parse.contains("Btor_scale")) {
      a_parm_parse.get("Btor_scale", m_Btor_scale);
   }
   else {
      m_Btor_scale = 0.;  // default
   }

   // Print geometry data to stdout
   if (m_verbose && !procID()) {
      cout << "Constructing cylindrical geometry..." << endl;
      cout << "rmin = " << m_rmin << ", rmax = " << m_rmax<<endl;
      cout << "zmin = " << m_zmin << ", zmax = " << m_zmax<<endl;
#if CFG_DIM==3
      cout << "phi_min = " << m_phimin << ", phi_max = " << m_phimax<<endl;
#endif
      cout << "Btor_scale = " << m_Btor_scale<< endl;
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
      
      a_BField(iv,0) = 0.0;
      a_BField(iv,1) = m_Btor_scale / R;
      a_BField(iv,2) = 0.0;
      
      a_BFieldMag(iv,0) = abs(m_Btor_scale / R);
      
      a_BFieldDir(iv,0) = 0.0;
      a_BFieldDir(iv,1) = m_Btor_scale/abs(m_Btor_scale);
      a_BFieldDir(iv,2) = 0.0;

      a_gradBFieldMag(iv,0) = -m_Btor_scale/(R * R);
      a_gradBFieldMag(iv,1) = 0.0;
      a_gradBFieldMag(iv,2) = 0.0;
      
      a_curlBFieldDir(iv,0) = 0.0;
      a_curlBFieldDir(iv,1) = 0.0;
      a_curlBFieldDir(iv,2) = (m_Btor_scale/abs(m_Btor_scale)) / R;
      
      a_BFieldDirdotcurlBFieldDir(iv,0) = 0.0;

   }
   
}

Vector<Real>
CylindricalBlockCoordSys::computeBField(const RealVect& a_X) const
{
   
   Vector<Real> result(3,0);
   
   result[0] =  0.0;
   result[1] =  m_Btor_scale / a_X[0];
   result[2] =  0.0;
   
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
   const Box& points_box = points.box();

   for (BoxIterator bit(points.box()); bit.ok(); ++bit) {
      IntVect iv = bit();

      RealVect xi;
      for (int n=0; n<SpaceDim; ++n) {
         xi[n] = points(iv,n);
      }

      RealVect X = realCoord(xi);

      A(iv,0) = X[1] * m_Btor_scale / X[0];
      A(iv,1) = 0.;
      A(iv,2) = 0.;

      Bmag(iv,0) = abs (m_Btor_scale / X[0]);
      
      b(iv,0) = 0.0;
      b(iv,1) = m_Btor_scale/abs(m_Btor_scale);
      b(iv,2) = 0.0;

   }
}

#include "NamespaceFooter.H"
