#include "SlabBlockCoordSys.H"
#include "SlabBlockCoordSysF_F.H"

#include "Directions.H"
#include "BoxIterator.H"
#include "ConstFact.H"
//#include "CONSTANTS.H"

#include "NamespaceHeader.H"

//char* SlabBlockCoordSys::pp_name = {"slab"};

const std::string SlabBlockCoordSys::pp_name = "slab";

SlabBlockCoordSys::SlabBlockCoordSys( ParmParse&               a_parm_parse,
                                      const ProblemDomain&     a_domain )
   : MagBlockCoordSys(a_parm_parse)
{

   // Read the input data specific to this geometry

   // Get max values of x, y, (and if 5D, z)
   a_parm_parse.get("x_max", m_xmax);
   // convention: all computational coordinates run from 0 to 2 pi (different from others where x coord runs from
   // flux_min to flux_max.  
   m_twopi = 2.*Constants::PI;
   m_outer_radial_boundary = m_twopi;
   m_inner_radial_boundary = 0.;
   a_parm_parse.get("y_max", m_ymax);
#if CFG_DIM==3
   a_parm_parse.get("z_max", m_zmax);
#endif

   // Get Bz at inner (x=0) and outer boundaries
   a_parm_parse.get("Bz_inner", m_BzInner);
   a_parm_parse.get("Bz_outer", m_BzOuter);

   // Get By at inner boundary; set at outer boundary by zero sheer. Default is zero
   if (a_parm_parse.contains("By_inner")) {
     a_parm_parse.get("By_inner", m_ByInner);
   }
   else {
     m_ByInner = 0.; //default
   }

   // Print geometry data to stdout
   if (m_verbose && !procID()) {
      cout << "Constructing slab geometry..." << endl;
      cout << "xmax = " << m_xmax << ", ymax = " << m_ymax 
#if CFG_DIM==3
           << ", zmax = " << m_zmax 
#endif
           << endl;
      cout << "Bz_inner = " << m_BzInner << ", Bz_outer = " << m_BzOuter << endl;
      cout << "By_inner  = " << m_ByInner << endl;
   }

   // Compute the mesh size in computational coordinates
   // 5/4/15: make y and z computational extents run from 0 to 2 pi so as to be consistent with the way periodic i.c.'s 
   // are set for Miller.  1/12/16, introduce m_ximax into normalization in case we change this choice in the future
   m_ximax = m_twopi;
   IntVect dimensions = a_domain.size();
   double dx = m_ximax/(double)dimensions[0];
   double dy = m_ximax/(double)dimensions[1];
#if CFG_DIM==3
   double dz = m_ximax/(double)dimensions[2];
   RealVect cellSpacing(dx,dy,dz);
#endif

#if CFG_DIM==2
   RealVect cellSpacing(dx,dy);
#endif

   cout << "**** cellSp[0] = " << cellSpacing[0] << ",  cellSp[1] = " << cellSpacing[1] 
#if CFG_DIM==3
        << ", cellSp[2] = " << cellSpacing[2]
#endif  
        << endl;

   cout << "**** dimensions = " << dimensions[0] <<", " << dimensions[1] 
#if CFG_DIM==3
        << ", " << dimensions[2]
#endif  
        << endl;

   D_TERM(m_xmaxvec[0] = m_xmax;,
          m_xmaxvec[1] = m_ymax;,
          m_xmaxvec[2] = m_zmax;)

   // Finish defining the object now that we also have the mesh spacing
   define( a_domain, cellSpacing );

   if (m_verbose && procID()==0) {
      cout << "Done constructing slab geometry" << endl;
   }
}


RealVect SlabBlockCoordSys::realCoord( const RealVect& a_xi ) const
{
   RealVect x;
   D_TERM(x[0] = m_xmax*a_xi[0]/m_ximax;,
          x[1] = m_ymax*a_xi[1]/m_ximax;,
          x[2] = m_zmax*a_xi[2]/m_ximax;)
   return x;
}

RealVect SlabBlockCoordSys::mappedCoord( const RealVect& a_x ) const
{

   RealVect xi;
   D_TERM(xi[0] = m_ximax*a_x[0]/m_xmax;,
          xi[1] = m_ximax*a_x[1]/m_ymax;,
          xi[2] = m_ximax*a_x[2]/m_zmax;)

   return xi;
}

Real SlabBlockCoordSys::dXdXi( const RealVect& a_Xi,
                                 int             a_dirX,
                                 int             a_dirXi ) const
{
  Real value = 0.0;
#if CFG_DIM==3
  double z = a_Xi[2];
#endif
      if (a_dirX != a_dirXi) {
        value = 0.0;
      }
      else{
        value = m_xmaxvec[a_dirX]/m_ximax;
      }

  return value;
}

void SlabBlockCoordSys::dXdXi(FArrayBox&       a_dXdXi,
                                const FArrayBox& a_Xi,
                                int              a_destComp,
                                int              a_dirX,
                                int              a_dirXi,
                                const Box&       a_box) const
{
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit) {
      const IntVect& iv = bit();
      if (a_dirX != a_dirXi) {
        a_dXdXi(iv,a_destComp) = 0;
      }
      else{
        a_dXdXi(iv,a_destComp) = m_xmaxvec[a_dirX]/m_ximax;
      }
   }
}

void
SlabBlockCoordSys::computeFieldData( const int  a_dir,
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

   //   FArrayBox RBpol(box, 1);
   //   FArrayBox dRBpoldt(box, 1);

   //   getRBpoloidal(a_dir, RBpol, dRBpoldt);

   FArrayBox xix(box, 1);

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

      xix(iv)    = xi[RADIAL_DIR];
   } 

   FORT_GET_SLAB_FIELD_DATA(CHF_BOX(box),
                              CHF_CONST_FRA1(xix,0),
                              CHF_CONST_REAL(m_BzInner),
                              CHF_CONST_REAL(m_BzOuter),
                              CHF_CONST_REAL(m_ByInner),
                              CHF_CONST_REAL(m_xmax),
                              CHF_CONST_REAL(m_ximax),
                              CHF_FRA(a_BField),
                              CHF_FRA1(a_BFieldMag,0),
                              CHF_FRA(a_BFieldDir),
                              CHF_FRA(a_gradBFieldMag),
                              CHF_FRA(a_curlBFieldDir),
                              CHF_FRA1(a_BFieldDirdotcurlBFieldDir,0));
}



void SlabBlockCoordSys::getMagneticFlux( const FArrayBox& physical_coordinates,
                                         FArrayBox&       magnetic_flux ) const
{
   MayDay::Error("SlabBlockCoordSys::getMagneticFlux() not implemented");
}



double SlabBlockCoordSys::getMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   MayDay::Error("SlabBlockCoordSys::getMagneticFlux() not implemented");
   return(0.0);
}



void SlabBlockCoordSys::getCellCenteredMappedCoords(FArrayBox& a_xi) const
{

   FORT_GET_SLAB_CC_MAPPED_COORDS(CHF_BOX(a_xi.box()),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_FRA(a_xi));
}



void SlabBlockCoordSys::getFaceCenteredMappedCoords(const int a_dir,
                                                      FArrayBox&  a_xi) const
{

   FORT_GET_SLAB_FC_MAPPED_COORDS(CHF_BOX(a_xi.box()),
                                    CHF_CONST_INT(a_dir),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_FRA(a_xi));
}


#include "NamespaceFooter.H"
