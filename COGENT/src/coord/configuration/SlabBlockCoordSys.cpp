#include "SlabBlockCoordSys.H"
#include "SlabBlockCoordSysF_F.H"

#include "Directions.H"
#include "BoxIterator.H"
#include "ConstFact.H"
//#include "CONSTANTS.H"

#include "NamespaceHeader.H"

const std::string SlabBlockCoordSys::pp_name = "slab";

SlabBlockCoordSys::SlabBlockCoordSys( ParmParse&               a_parm_parse,
                                      const ProblemDomain&     a_domain,
                                      const int                a_block,
                                      const int                a_num_blocks,
                                      const int                a_block_separation)
   : MagBlockCoordSys(a_parm_parse),
     m_is_field_aligned(false),
     m_block(a_block),
     m_num_blocks(a_num_blocks)
{

   // Read the input data specific to this geometry
   parseParameters (a_parm_parse);
   
   // Compute the mesh size in computational coordinates
   // for slab geometry all computational coordinates run from 0 to 2 pi
   // however, for num_blocks>1 the block-normal mapped coordinate goes
   // from 2pi/num_blocks * m_block to 2pi/num_blocks * (m_block + 1)
   
#if CFG_DIM == 3
   m_mb_dir = TOROIDAL_DIR;
#else
   //this options has not been fully developed yet
   m_mb_dir = POLOIDAL_DIR;
#endif

   double twopi = 2.*Constants::PI;
   m_mapped_block_size = twopi * RealVect::Unit;
   m_mapped_block_size[m_mb_dir] = twopi/(double)m_num_blocks;
   
   IntVect dimensions = a_domain.size();
   RealVect dx;
   for (int dir=0; dir<SpaceDim; ++dir) {
      dx[dir] = m_mapped_block_size[dir]/(double)dimensions[dir];
   }

   //This defines the block_coord_sys object (e.g., sets m_dx)
   define( a_domain, dx );

   //Define the block size in the physical coordinates
   m_phys_block_size[0] = m_xmax;
#if CFG_DIM==2
   m_phys_block_size[1] = m_zmax;
#else
   m_phys_block_size[1] = m_ymax;
   m_phys_block_size[2] = m_zmax;
#endif
   
   //Physical block size in the multiblock dir
   m_phys_block_size[m_mb_dir] /= (double)m_num_blocks;

   //Physical coordinate of a block center in the multiblock dir
   m_mb_center = m_phys_block_size[m_mb_dir] * ((double)m_block + 0.5);
   
   //Block separation measured in mapped space
   m_separation = dx[m_mb_dir] * a_block_separation;
   
   // Print geometry data
   if (m_verbose) {
      printParameters();
   }
}


RealVect SlabBlockCoordSys::realCoord( const RealVect& a_xi ) const
{
   RealVect x = a_xi;
   
   //Remove block separations
   x[m_mb_dir] -= m_separation * m_block;
   
   //Do the scaling
   x *= m_phys_block_size;
   x /= m_mapped_block_size;
 
#if CFG_DIM ==3
   if (m_is_field_aligned) {
      double pitch = getPitch(x);
      x[POLOIDAL_DIR] += (x[m_mb_dir] - m_mb_center) * pitch;
   }
#endif

   return x;
}

RealVect SlabBlockCoordSys::mappedCoord( const RealVect& a_x ) const
{

   RealVect xi = a_x;
   
   //Do the scaling
   xi *= m_mapped_block_size;
   xi /= m_phys_block_size;
   
   //Add block separations
   xi[m_mb_dir] += m_separation * m_block;

#if CFG_DIM ==3
   if (m_is_field_aligned) {
      //Because pitch is constant along the field-line
      double pitch = getPitch(a_x);
      double shift_phys = (a_x[m_mb_dir] - m_mb_center) * pitch;
      double shift_mapped = shift_phys * m_mapped_block_size[POLOIDAL_DIR]
                                       / m_phys_block_size[POLOIDAL_DIR];
      xi[POLOIDAL_DIR] -=  shift_mapped;
   }
#endif
   
   return xi;
}

Real SlabBlockCoordSys::dXdXi( const RealVect& a_Xi,
                               int             a_dirX,
                               int             a_dirXi ) const
{
  Real value = 0.0;

  if (a_dirX != a_dirXi) {
     value = 0.0;
  }
  else{
     value = m_phys_block_size[a_dirX]/m_mapped_block_size[a_dirX];
  }

#if CFG_DIM ==3
   if (m_is_field_aligned) {

      RealVect X = realCoord(a_Xi);

      if (a_dirX == POLOIDAL_DIR && a_dirXi == TOROIDAL_DIR) {
         double pitch = getPitch(X);
         value += pitch * m_phys_block_size[TOROIDAL_DIR]
                        / m_mapped_block_size[TOROIDAL_DIR];
      }

      else if (a_dirX == POLOIDAL_DIR && a_dirXi == RADIAL_DIR) {
         double pitch_shear = getPitchShear(X);
         value += pitch_shear * ( X[TOROIDAL_DIR] - m_mb_center);
      }

   }
#endif
  
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
     RealVect this_Xi;
     for (int dir=0; dir<SpaceDim; ++dir) {
        this_Xi[dir] = a_Xi(iv,dir);
     }
     a_dXdXi(iv,a_destComp) = dXdXi(this_Xi,a_dirX,a_dirXi);
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
                            CHF_CONST_REAL(m_ByInner),
                            CHF_CONST_REAL(m_ByOuter),
                            CHF_CONST_REAL(m_BzInner),
                            CHF_CONST_REAL(m_BzOuter),
                            CHF_CONST_REAL(m_phys_block_size[0]),
                            CHF_CONST_REAL(m_mapped_block_size[0]),
                            CHF_FRA(a_BField),
                            CHF_FRA1(a_BFieldMag,0),
                            CHF_FRA(a_BFieldDir),
                            CHF_FRA(a_gradBFieldMag),
                            CHF_FRA(a_curlBFieldDir),
                            CHF_FRA1(a_BFieldDirdotcurlBFieldDir,0));
}

Vector<Real>
SlabBlockCoordSys::computeBField(const RealVect& a_X) const
{
   
   Vector<Real> result(3,0);
   
   result[0] = 0;
   result[1] = m_ByInner + (m_ByOuter - m_ByInner) * a_X[0]/m_xmax;
   result[2] = m_BzInner + (m_BzOuter - m_BzInner) * a_X[0]/m_xmax;
   
   return result;
   
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

double SlabBlockCoordSys::getPitch( const RealVect& a_physical_coordinate ) const
{
   Vector<Real> B = computeBField(a_physical_coordinate);
   double pitch = B[2] / B[1];
   return pitch;
}


double SlabBlockCoordSys::getPitchShear( const RealVect& a_physical_coordinate ) const
{
   
   double dBydX = (m_ByOuter - m_ByInner) /m_xmax;
   double dBzdX = (m_BzOuter - m_BzInner) /m_xmax;
   
   Vector<Real> B = computeBField(a_physical_coordinate);
   double Bz = B[2];
   double By = B[1];
   
   double pitch_shear = (dBzdX * By - dBydX * Bz) / (By * By);
   double mapped_pitch_shear = pitch_shear;
   mapped_pitch_shear *= m_phys_block_size[RADIAL_DIR]/m_mapped_block_size[RADIAL_DIR];
   
   return mapped_pitch_shear;
}


void SlabBlockCoordSys::getNodalFieldData(FArrayBox& points,
                                          FArrayBox& A,
                                          FArrayBox& b,
                                          FArrayBox& Bmag) const
{
   const Box& points_box = points.box();

   FArrayBox BField(points_box,3);
   FArrayBox BFieldMag(points_box,1);
   FArrayBox BFieldDir(points_box,3);
   FArrayBox gradBFieldMag(points_box,3);
   FArrayBox curlBFieldDir(points_box,3);
   FArrayBox BFieldDirdotcurlBFieldDir(points_box,1);

   computeFieldData( -1, BField, BFieldMag, BFieldDir, gradBFieldMag,
                     curlBFieldDir, BFieldDirdotcurlBFieldDir, false );

   for (BoxIterator bit(points.box()); bit.ok(); ++bit) {
      IntVect iv = bit();

      RealVect xi;
      for (int n=0; n<SpaceDim; ++n) {
         xi[n] = points(iv,n);
      }

      RealVect X = realCoord(xi);

      A(iv,0) = 0.;
      A(iv,1) =  X[0] * BField(iv,2);
      A(iv,2) = -X[0] * BField(iv,1);

      for (int n=0; n<3; ++n) {
         b(iv,n) = BFieldDir(iv,n);
      }

      Bmag(iv,0) = BFieldMag(iv,0);
   }
}

void SlabBlockCoordSys::parseParameters(ParmParse& a_parm_parse)
{
   
   // Get max values of x, z, (and if 5D, y)
   a_parm_parse.get("x_max", m_xmax);
   a_parm_parse.get("z_max", m_zmax);
#if CFG_DIM==3
   a_parm_parse.get("y_max", m_ymax);
#endif
   
   // Get By at inner (x=0) and outer boundaries
   a_parm_parse.get("By_inner", m_ByInner);
   a_parm_parse.get("By_outer", m_ByOuter);
   
   // Get Bz at inner and outer boundaries
   m_BzInner = 0.; //default
   m_BzOuter = 0.; //default
   a_parm_parse.query("Bz_inner", m_BzInner);
   a_parm_parse.query("Bz_outer", m_BzOuter);
   
   //Is the coordinate system field aligned
   a_parm_parse.query("field_aligned", m_is_field_aligned);
   
}
void SlabBlockCoordSys::printParameters()
{
   // Print geometry data to stdout
   if (procID()==0 && m_block == 0) {
      cout << "Slab geometry parameters" << endl;
      cout << "xmax = " << m_xmax << ", zmax = " << m_zmax
#if CFG_DIM==3
      << ", ymax = " << m_ymax
#endif
      << endl;
      cout << "Number of blocks =" << m_num_blocks << endl;
      cout<<"Physical block size:" << m_phys_block_size << endl;
      cout<<"Mapped block size:" << m_mapped_block_size << endl;
      
      cout << "By_inner = " << m_ByInner << ", By_outer = " << m_ByOuter << endl;
      cout << "Bz_inner = " << m_BzInner << ", Bz_outer = " << m_BzOuter << endl;

      cout << "Field aligned = "<< m_is_field_aligned << endl;
   }
}

#include "NamespaceFooter.H"
