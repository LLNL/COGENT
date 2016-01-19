#include "RectangularTorusBlockCoordSys.H"
#include "RectangularTorusBlockCoordSysF_F.H"

#include "RectangularTorusDirections.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

//char* RectangularTorusBlockCoordSys::pp_name = {"rectangulartorus"};

const std::string RectangularTorusBlockCoordSys::pp_name = "rectangulartorus";

RectangularTorusBlockCoordSys::RectangularTorusBlockCoordSys( ParmParse&               a_parm_parse,
                                          const ProblemDomain&     a_domain )
   : MagBlockCoordSys(a_parm_parse)
{
   // Read the input data specific to this geometry

   m_average_radius = 0.5*(m_inner_radial_boundary + m_outer_radial_boundary);

   // Get the flag to use a modified geometry where the major radius is infinite
   if (a_parm_parse.contains("inf_majorrad")) {
     a_parm_parse.get("inf_majorrad",m_inf_majorrad);
   }
   else {
      m_inf_majorrad = 0;
   }


   // Get the coordinate system offset
   a_parm_parse.getarr("origin", m_origin, 0, SpaceDim);

   // Get the Z boundaries
   a_parm_parse.get("lower_z_bdry", m_lower_z_boundary);
   a_parm_parse.get("upper_z_bdry", m_upper_z_boundary);

#if CFG_DIM==3
   // Get the phi boundaries
   a_parm_parse.get("min_phi_bdry", m_min_phi_boundary);
   a_parm_parse.get("max_phi_bdry", m_max_phi_boundary);
#endif

   // Get the toroidal field component scale factor
   if (a_parm_parse.contains("Btor_scale")) {
      a_parm_parse.get("Btor_scale", m_Btor_scale);
   }
   else {
      m_Btor_scale = 0.;  // default
   }

   // Get the vertical field
   a_parm_parse.get("Bz_inner", m_bz_inner);
   a_parm_parse.get("Bz_outer", m_bz_outer);

   // Print geometry data to stdout
   if (m_verbose && !procID()) {
      cout << "Constructing RectangularTorus equilibrium geometry..." << endl;
      cout << "   origin =";
      for (int n=0; n<SpaceDim; n++) {
         cout << " " << m_origin[n];
      }
      cout << endl;
      cout << "   inner_radial_bdry = " << m_inner_radial_boundary << endl;
      cout << "   outer_radial_bdry = " << m_outer_radial_boundary << endl;
      cout << "   lower_z_bdry = " << m_lower_z_boundary << endl;
      cout << "   upper_z_bdry = " << m_upper_z_boundary << endl;
#if CFG_DIM==3
      cout << "   min_phi_bdry = " << m_min_phi_boundary << endl;
      cout << "   max_phi_bdry = " << m_max_phi_boundary << endl;
#endif
   }

   // Compute the mesh size
   IntVect dimensions = a_domain.size();
   double dr = (m_outer_radial_boundary - m_inner_radial_boundary) / (double)dimensions[R_DIR];
#if CFG_DIM==3
   double dphi = (m_max_phi_boundary - m_min_phi_boundary) / (double)dimensions[PHI_DIR];
   double dz = (m_upper_z_boundary - m_lower_z_boundary) / (double) dimensions[Z_DIR];
   RealVect cellSpacing(dr,dphi,dz);
#elif CFG_DIM==2
   double dz = (m_upper_z_boundary - m_lower_z_boundary) / (double) dimensions[Z_DIR];
   RealVect cellSpacing(dr,dz);
#endif

   // Finish defining the object now that we also have the mesh spacing
   define( a_domain, cellSpacing );

   if (m_verbose && procID()==0) {
      cout << "Done constructing RectangularTorus equilibrium geometry" << endl;
   }
}



RealVect RectangularTorusBlockCoordSys::realCoord( const RealVect& a_xi ) const
{
   RealVect x;
   double r = a_xi[R_DIR] + m_inner_radial_boundary;

#if CFG_DIM==3
   double phi = a_xi[PHI_DIR];
   x[0] = r*cos(phi);
   x[1] = r*sin(phi);
   x[2] = a_xi[Z_DIR];
#elif CFG_DIM==2
   x[0] = r;
   x[1] = a_xi[Z_DIR];
#endif

   return x;
}



RealVect RectangularTorusBlockCoordSys::mappedCoord( const RealVect& a_x ) const
{
   RealVect xi;

   double x_translated = a_x[0];
#if CFG_DIM==3
   double y_translated = a_x[1];
   double z_translated = a_x[2];
#elif CFG_DIM==2
   double z_translated = a_x[1];
#endif

   double r;
   double two_pi = 2.*Pi;

#if CFG_DIM==3
   double phi = atan2(y_translated, x_translated);

   double st, ct;

   while (true) {

      st  = sin(phi);
      ct  = cos(phi);

      double g = x_translated*st - y_translated*ct;

      if (fabs(g) > 1.e-11) {
         double gprime = x_translated*ct + y_translated*st*(1. + ct);
         phi -= g/gprime;
      }
      else break;
   }

   if( fabs(st) > 1.e-11) {
      r = y_translated / (st);
   }
   else {
      r = x_translated / ct;
   }

   // Now shift phi to the range 0 <= phi < 2pi.
   if (y_translated < 0) phi += two_pi;

   xi[0] = r - m_inner_radial_boundary;
   xi[1] = phi;
   xi[2] = z_translated;
#elif CFG_DIM==2
   xi[0] = x_translated - m_inner_radial_boundary;
   xi[1] = z_translated;
#endif
   return xi;
}



Real RectangularTorusBlockCoordSys::dXdXi( const RealVect& a_Xi,
                                 int             a_dirX,
                                 int             a_dirXi ) const
{
  //calculate d(x,y,z)/d(r,phi,z) where a_Xi is for (translated r from zero, phi, z)
  Real value = 0.0;
  // this is for test purposes

  double r = a_Xi[0] + m_inner_radial_boundary;
#if CFG_DIM==3
  double phi = a_Xi[1];
  double z = a_Xi[2];

  if (a_dirX == 0)
    {
      if (a_dirXi == 0)
        {
          // dx/dr = cos(phi)
          value = cos(phi);
        }
      else if (a_dirXi == 1)
        {
          // dx/dphi = -rsin(phi)
          value = -r*sin(phi);
        }
      // else value stays 0
    }
  else if (a_dirX == 1)
    {
      if (a_dirXi == 0)
        {
          // dy/dr = sin(phi)
          value =  sin(phi);
        }
      else if (a_dirXi == 1)
        {
          // dy/dphi = r cos(phi)
          value = r * cos(phi);
        }
      // else value stays 0
    }
  else if (a_dirX == 2)
    {
      if (a_dirXi == 2)
        {
          // dz/dz = 1.0
          value = 1.0;
        }
      // else value stays 0
    }
  else
    {
       MayDay::Error("Bad dirX in RectangularTorusBlockCoordSys::dXdXi");
    }

#elif CFG_DIM==2
  double z = a_Xi[1];
  if (a_dirX == 0)
    {
      if (a_dirXi == 0)
        {
          // dx/dr = 1.0
          value = 1.0;
        }
      // else value stays 0 eg. dx/dz=0
    }
  else if (a_dirX == 1)
    {
      if (a_dirXi == 1)
        {
          // dz/dz = 1.0
          value = 1.0;
        }
      // else value stays 0 eg. dz/dx=0
    }
  else
    {
       MayDay::Error("Bad dirX in RectangularTorusBlockCoordSys::dXdXi");
    }

#endif
  return value;
}



void RectangularTorusBlockCoordSys::dXdXi(FArrayBox&       a_dXdXi,
                                const FArrayBox& a_Xi,
                                int              a_destComp,
                                int              a_dirX,
                                int              a_dirXi,
                                const Box&       a_box) const
{
   FArrayBox Xi_shifted(a_Xi.box(), a_Xi.nComp());
   Xi_shifted.copy(a_Xi);
   Xi_shifted.plus(m_inner_radial_boundary,0,1);

   FORT_RECTANGULARTORUS_DXDXI(CHF_BOX(a_box),
                     CHF_CONST_FRA(Xi_shifted),
                     CHF_CONST_INT(a_destComp),
                     CHF_CONST_INT(a_dirX),
                     CHF_CONST_INT(a_dirXi),
                     CHF_FRA(a_dXdXi));
}


void
RectangularTorusBlockCoordSys::getRBvertical( const int  a_dir,
                                    FArrayBox& a_RBvert) const
{
   const Box& box = a_RBvert.box();
   int inf_major_radius = getInfMajorrad();

   FORT_GET_RBVERT_RECTANGULARTORUS(CHF_CONST_INT(a_dir),
                         CHF_BOX(box),
                         CHF_CONST_REALVECT(m_dx),
                         CHF_CONST_REAL(m_bz_inner),
                         CHF_CONST_REAL(m_bz_outer),
                         CHF_CONST_INT(inf_major_radius),
                         CHF_CONST_REAL(m_inner_radial_boundary),
                         CHF_CONST_REAL(m_outer_radial_boundary),
                         CHF_FRA1(a_RBvert,0));
}



void
RectangularTorusBlockCoordSys::computeFieldData( const int  a_dir,
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

   FArrayBox RBvert(box, 1);
//   FArrayBox dRBpoldt(box, 1);

   getRBvertical(a_dir, RBvert);

   FArrayBox R(box, 1);
#if CFG_DIM==3
   FArrayBox dXdr(box, 1);
   FArrayBox dXdp(box, 1);
   FArrayBox dXdz(box, 1);
   FArrayBox dYdr(box, 1);
   FArrayBox dYdp(box, 1);
   FArrayBox dYdz(box, 1);
   FArrayBox dZdr(box, 1);
   FArrayBox dZdp(box, 1);
   FArrayBox dZdz(box, 1);
#elif CFG_DIM==2
   FArrayBox dXdr(box, 1);
   FArrayBox dXdz(box, 1);
   FArrayBox dZdr(box, 1);
   FArrayBox dZdz(box, 1);
#endif
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

      R(iv)    = realCoord(xi)[R_DIR];
#if CFG_DIM==3
      dXdr(iv) = dXdXi(xi, 0, R_DIR);
      dXdp(iv) = dXdXi(xi, 0, PHI_DIR);
      dXdz(iv) = dXdXi(xi, 0, Z_DIR);
      dYdr(iv) = dXdXi(xi, 1, R_DIR);
      dYdp(iv) = dXdXi(xi, 1, PHI_DIR);
      dYdz(iv) = dXdXi(xi, 1, Z_DIR);
      dZdr(iv) = dXdXi(xi, 2, R_DIR);
      dZdp(iv) = dXdXi(xi, 2, PHI_DIR);
      dZdz(iv) = dXdXi(xi, 2, Z_DIR);
#elif CFG_DIM==2
      dXdr(iv) = dXdXi(xi, 0, R_DIR);
      dXdz(iv) = dXdXi(xi, 0, Z_DIR);
      dZdr(iv) = dXdXi(xi, 1, R_DIR);
      dZdz(iv) = dXdXi(xi, 1, Z_DIR);
#endif
   }

   double RBtor = getRBtoroidal();
   int inf_major_radius = getInfMajorrad();
   double rmidreal = getAvMinorRad();
   FORT_GET_RECTANGULARTORUS_FIELD_DATA(
                              CHF_BOX(box),
                              CHF_CONST_INT(inf_major_radius),
                              CHF_CONST_REAL(rmidreal),
                              CHF_CONST_FRA1(RBvert,0),
                              CHF_CONST_REAL(RBtor),
                              CHF_CONST_FRA1(R,0),
                              CHF_CONST_FRA1(dXdr,0),
                              CHF_CONST_FRA1(dXdz,0),
                              CHF_CONST_FRA1(dZdr,0),
                              CHF_CONST_FRA1(dZdz,0),
                              CHF_FRA(a_BField),
                              CHF_FRA1(a_BFieldMag,0),
                              CHF_FRA(a_BFieldDir),
                              CHF_FRA(a_gradBFieldMag),
                              CHF_FRA(a_curlBFieldDir),
                              CHF_FRA1(a_BFieldDirdotcurlBFieldDir,0));
}



void RectangularTorusBlockCoordSys::getCellCenteredMappedCoords(FArrayBox& a_xi) const
{
   double rmid = getAvMinorRad() - m_inner_radial_boundary;
   int inf_major_radius = getInfMajorrad();

   FORT_GET_RECTANGULARTORUS_CC_MAPPED_COORDS(CHF_BOX(a_xi.box()),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_CONST_REAL(rmid),
                                    CHF_CONST_INT(inf_major_radius),
                                    CHF_FRA(a_xi));
}



void RectangularTorusBlockCoordSys::getFaceCenteredMappedCoords(const int a_dir,
                                                      FArrayBox&  a_xi) const
{
   double rmid = getAvMinorRad() - m_inner_radial_boundary;
   int inf_major_radius = getInfMajorrad();

   FORT_GET_RECTANGULARTORUS_FC_MAPPED_COORDS(CHF_BOX(a_xi.box()),
                                    CHF_CONST_INT(a_dir),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_CONST_REAL(rmid),
                                    CHF_CONST_INT(inf_major_radius),
                                    CHF_FRA(a_xi));
}



#include "NamespaceFooter.H"
