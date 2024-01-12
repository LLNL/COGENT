#include <array>
#include <cmath>
#include <float.h>
#include <fstream>
#include "OneBlockCoordSys.H"
#include "Directions.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"

#include "inspect.H"

#undef REPORT_NEWTON_FAILURE

#include "NamespaceHeader.H"

const std::string OneBlockCoordSys::pp_name = "oneblock";

OneBlockCoordSys::OneBlockCoordSys(ParmParse&            a_parm_parse,
                                   const ProblemDomain&  a_domain,
                                   const int             a_numcells_multiblock,
                                   const int             a_toroidal_sector,
                                   const int             a_toroidal_block_separation)
   : MagBlockCoordSys(a_parm_parse),
     m_poloidal_block(0),
     m_toroidal_sector(a_toroidal_sector),
     m_poloidal_util(NULL),
     m_field_aligned_mapping(false),
     m_mapping_includes_RBdata(true),
     m_phys_coord_type(CARTESIAN),
#if CFG_DIM==3
     m_interp(NULL),
#endif
     m_toroidal_ghosts(4),
     m_toroidal_mapping_refinement(2),
     m_field_trace_step_num(10)
{

   /// Read the input data specific to this geometry

   /// Get grid information values
   a_parm_parse.get("r_min", m_rmin);
   a_parm_parse.get("r_max", m_rmax);
   
   a_parm_parse.get("theta_min", m_thetamin);
   a_parm_parse.get("theta_max", m_thetamax);

#if CFG_DIM==3
   a_parm_parse.get( "toroidal_width_over_2pi", m_toroidal_width);
   CH_assert(m_toroidal_width > 0. && m_toroidal_width <= 1.);
   m_toroidal_width *= 2.*Pi;
#endif
   
   // Print geometry data to stdout
   if (m_verbose && !procID()) {
      cout << "Constructing oneblock geometry..." << endl;
      cout << "rmin = " << m_rmin << ", rmax = " << m_rmax<<endl;
      cout << "thetamin = " << m_thetamin << ", thetamax = " << m_thetamax<<endl;
#if CFG_DIM==3
      cout << "toroidal_width = " << m_toroidal_width <<endl;
#endif
   }

   // Domain object that is passed here is a block domain,
   // which contains all cells in poloidal directions, but only
   // a fraction of cells in the toroidal direction.
   // a_numcells_multiblock corresponds to a total number of valid
   // toroidal cells in the entire problem domain
   IntVect dimensions = a_domain.size();
   double dr = (m_rmax - m_rmin)/(double)dimensions[RADIAL_DIR];
   double dtheta = (m_thetamax - m_thetamin)/(double)dimensions[POLOIDAL_DIR];
#if CFG_DIM==3
   double dphi = m_toroidal_width/(double)a_numcells_multiblock;
   RealVect cellSpacing(dr,dphi,dtheta);
#endif

#if CFG_DIM==2
   RealVect cellSpacing(dr,dtheta);
#endif

   // Finish defining the object now that we also have the mesh spacing
   define( a_domain, cellSpacing );
   
   
   // Get the toroidal field component scale factor
   if (a_parm_parse.contains("Btor_scale")) {
      a_parm_parse.get("Btor_scale", m_RB_toroidal);
   }
   else {
      m_RB_toroidal = 0.;  // default
   }

   if (a_parm_parse.contains("field_aligned_mapping") && SpaceDim == 3) {
      a_parm_parse.get("field_aligned_mapping", m_field_aligned_mapping);
   }
   
   // Some applications compute field self-consistently, in which case
   // we set dummy Bfield to zero to avoid code breaking
   a_parm_parse.query("mapping_includes_RBdata", m_mapping_includes_RBdata);
   
   // If a field coefficients file was specified in the input, use it
   // to define the field.  Otherwise, we read it from the mapping file.
   m_spectral_field = a_parm_parse.contains("field_coefficients_file");

   m_poloidal_util = new POL::SingleNullPoloidalUtil(a_parm_parse, m_spectral_field, m_RB_toroidal);

   m_provides_flux = m_spectral_field;
   
#if CFG_DIM==3
    m_toroidal_mapped_offset = a_toroidal_sector * a_toroidal_block_separation * dphi;
#endif
   
   readFiles( a_parm_parse );

   if (m_verbose && procID()==0) {
     cout << "Done constructing poloidal oneblock geometry" << endl;
   }
}

OneBlockCoordSys::~OneBlockCoordSys()
{
   if (m_poloidal_util) delete m_poloidal_util;
#if CFG_DIM==3
   if ( m_interp ) delete m_interp;
#endif
}

void
OneBlockCoordSys::readFiles( ParmParse& a_pp )
{

   // Number of elements in each string of mapping file
   // is either 4 (R,Z,RBr,RBz) or 2 (R,Z)
   int num_string_elements = (m_mapping_includes_RBdata) ? 4 : 2;
   
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
         MayDay::Error("OneBlockCoordSys::readFiles(): No geometry mapping file specified");
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

   POL::IntVect mapping_block_size;
   mapping_block_size[0] = n_radial_points   + 2*n_radial_extend;
   mapping_block_size[1] = n_poloidal_points + 2*n_poloidal_extend;
   num_block_elements = mapping_block_size[0] * mapping_block_size[1];
   double *input_values = new double[num_string_elements*num_block_elements];

#ifdef CH_MPI
   if (procID() == 0) {
#endif

      for (int i=0; i<num_string_elements*num_block_elements; ++i) {
         inFile >> input_values[i];
      }
      inFile.close();
      cout << "Done reading geometry mapping file" << endl;

#ifdef CH_MPI
   }

   MPI_Bcast(input_values, num_string_elements*num_block_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif // end if MPI

   int num_radial_cells = n_radial_points - 1;

   double dr = ( upperMappedCoordinate(RADIAL_DIR)
             -   lowerMappedCoordinate(RADIAL_DIR) ) / num_radial_cells;

   double *r_pts = new double[mapping_block_size[0]];
            
   r_pts[0] = lowerMappedCoordinate(RADIAL_DIR) - n_radial_extend*dr;
   for (int i=1; i<mapping_block_size[0]; ++i) {
      r_pts[i] = r_pts[i-1] + dr;
   }
   
   int num_poloidal_cells = n_poloidal_points - 1;

   double dtheta = ( upperMappedCoordinate(POLOIDAL_DIR)
                 -   lowerMappedCoordinate(POLOIDAL_DIR) ) / num_poloidal_cells;
   
   double *theta_pts = new double[mapping_block_size[1]];
   
   theta_pts[0] = lowerMappedCoordinate(POLOIDAL_DIR) - n_poloidal_extend*dtheta;
   for (int i=1; i<mapping_block_size[1]; ++i) {
      theta_pts[i] = theta_pts[i-1] + dtheta;
   }

   // Determine interpolation method
   string interp_method = "bspline";
   if (a_pp.contains("interpolation_method")) {
      a_pp.get("interpolation_method", interp_method);
   }

#if CFG_DIM==3
   // Since the mapping is only defined using bspline interpolation for now (and maybe forever),
   // we need to restrict the poloidal interpolation to bspline also for consistency
   if ( interp_method != "bspline" ) {
      MayDay::Error("OnelBlockCoordSys::readFiles(): interp_method is limited to bspline in 3D for now");
   }
#endif

   POL::Box interp_box(POL::IntVect::Zero,mapping_block_size-POL::IntVect::Unit);
   POL::FArrayBox interp_node_coords(interp_box, POL_DIM);

   // Create an interpolator for the cylindrical coordinates
   POL::FArrayBox RZ_data(interp_box, POL_DIM);

   POL::BoxIterator bit(interp_box);
   for (bit.begin();bit.ok();++bit) {
      POL::IntVect iv = bit();
      int k = mapping_block_size[0]*iv[1] + iv[0];

      RZ_data(iv,0) = input_values[num_string_elements*k];
      RZ_data(iv,1) = input_values[num_string_elements*k+1];

      interp_node_coords(iv,0) = r_pts[iv[0] - interp_box.smallEnd(0)];
      interp_node_coords(iv,1) = theta_pts[iv[1] - interp_box.smallEnd(1)];
   }
   
   m_poloidal_util->setRZInterp(a_pp, interp_method, interp_node_coords, RZ_data);
#if CFG_DIM==3
   setInterp(a_pp, interp_node_coords, RZ_data);
#endif
   
   
   if ( !m_spectral_field ) {

      if (procID() == 0) {
         cout << "Reading magnetic field data from mapping file" << endl;
      }

      // Derive R times the magnetic field magnitude and the cylindrical components
      // of the field unit vector.  We need these quantities and their derivatives
      // at various points, so we create an interpolator for them.

      POL::FArrayBox field_data(interp_box, POL::SingleNullPoloidalUtil::NUM_FIELD_VARS);

      for (bit.begin();bit.ok();++bit) {
         POL::IntVect iv = bit();
         int k = mapping_block_size[0]*iv[1] + iv[0];

         double RBr, RBz, RB;

         if (m_mapping_includes_RBdata) {
            field_data(iv,POL::SingleNullPoloidalUtil::RBR_VAR) = RBr = input_values[4*k+2];
            field_data(iv,POL::SingleNullPoloidalUtil::RBZ_VAR) = RBz = input_values[4*k+3];
         }
         else {
            field_data(iv,POL::SingleNullPoloidalUtil::RBR_VAR) = RBr = 0.;
            field_data(iv,POL::SingleNullPoloidalUtil::RBZ_VAR) = RBz = 0.;
         }
         field_data(iv,POL::SingleNullPoloidalUtil::RB_VAR) = RB = sqrt( RBr*RBr + RBz*RBz + m_RB_toroidal*m_RB_toroidal );
         field_data(iv,POL::SingleNullPoloidalUtil::UNIT_BR_VAR) = RBr / RB;
         field_data(iv,POL::SingleNullPoloidalUtil::UNIT_BPHI_VAR) = m_RB_toroidal / RB;
         field_data(iv,POL::SingleNullPoloidalUtil::UNIT_BZ_VAR) = RBz / RB;
      }

      m_poloidal_util->setFieldInterp(a_pp, interp_method, interp_node_coords, field_data);
   }
   
   delete [] r_pts;
   delete [] theta_pts;
   delete [] input_values;

   getCellCenterRealCoords(); 

}

void
OneBlockCoordSys::getCellCenterRealCoords()
{
#if CFG_DIM==2
   POL::Box box_pol = restrictToPoloidal(m_domain.domainBox());
   m_rpc_coarsen_ratio = 8*POL::IntVect::Unit;
   for (int n=0; n<POL_DIM; ++n) {
      if (box_pol.size(n) < 64) m_rpc_coarsen_ratio[n] = 1;
   }
   box_pol.coarsen(m_rpc_coarsen_ratio);
   m_realPoloidalCoords.define(box_pol, POL_DIM);
   POL::RealVect rpc_dx = restrictToPoloidal(m_dx) * m_rpc_coarsen_ratio;
   POL::RealVect offset_pol = 0.5*POL::RealVect::Unit;
   offset_pol *= rpc_dx;
   for (POL::BoxIterator bit(box_pol); bit.ok(); ++bit) {
      POL::IntVect iv = bit();
      POL::RealVect xi = rpc_dx*iv + offset_pol;
      POL::RealVect X_pol = realCoord(xi);

      for (int n=0; n<POL_DIM; ++n) {
         m_realPoloidalCoords(iv,n) = X_pol[n];
      }
   }
#endif
#if CFG_DIM==3
   Box box = m_domain.domainBox();
   m_rc_coarsen_ratio = 8*IntVect::Unit;
   for (int n=0; n<SpaceDim; ++n) {
      if (box.size(n) < 64) m_rc_coarsen_ratio[n] = 1;
   }
   box.coarsen(m_rc_coarsen_ratio);
   m_realCoords.define(box, SpaceDim);
   RealVect rc_dx = m_dx * m_rc_coarsen_ratio;
   RealVect offset = 0.5*RealVect::Unit;
   offset *= rc_dx;
   for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect xi = rc_dx*iv + offset;
      RealVect X = realCoord(xi);

      for (int n=0; n<SpaceDim; ++n) {
         m_realCoords(iv,n) = X[n];
      }
   }
#endif
}

RealVect
OneBlockCoordSys::realCoord( const RealVect& a_Xi ) const
{
   POL::RealVect Xi_pol = restrictToPoloidal(a_Xi);

   RealVect X;

#if CFG_DIM==2
   X[0] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::R_VAR, POL::Interp::FUNCTION, Xi_pol);
   X[1] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::Z_VAR, POL::Interp::FUNCTION, Xi_pol);
#endif
#if CFG_DIM==3
    double R   = m_interp->interpolate(0, Interp3D::FUNCTION, a_Xi);
    double Phi = m_interp->interpolate(1, Interp3D::FUNCTION, a_Xi) - m_toroidal_mapped_offset;
    double Z   = m_interp->interpolate(2, Interp3D::FUNCTION, a_Xi);

    if ( m_phys_coord_type == CYLINDRICAL ) {
       D_TERM(X[0] = R;,
              X[1] = Phi;,
              X[2] = Z;)
          }
    else if ( m_phys_coord_type == CARTESIAN ) {
       D_TERM(X[0] = R*cos(Phi);,
              X[1] = R*sin(Phi);,
              X[2] = Z;)
          }
    else {
       MayDay::Error("OneBlockCoordSys::realCoord(): unrecognized m_phys_coord_type");
    }
#endif
   
   return X;
}



RealVect
OneBlockCoordSys::mappedCoord( const RealVect& a_X ) const
{
   /*
     Search for the cell whose centroid real space coordinate is closest to
     the passed argument and use the corresponding mapped coordinate as the
     initial guess for the Newton solve.  This should be replaced by something
     more efficient.
   */

   double max_dist = DBL_MAX;

#if CFG_DIM==2
   POL::RealVect X_pol = restrictPhysCoordToPoloidal(a_X);
   POL::RealVect xi_initial = POL::RealVect::Zero;

   POL::RealVect rpc_dx = restrictToPoloidal(m_dx) * m_rpc_coarsen_ratio;
   POL::RealVect offset = 0.5*POL::RealVect::Unit;
   offset *= rpc_dx;

   POL::IntVect iv_closest;

   for (POL::BoxIterator bit(m_realPoloidalCoords.box()); bit.ok(); ++bit) {
      POL::IntVect iv = bit();
      POL::RealVect center_xi = rpc_dx*iv + offset;
      POL::RealVect center_X;
      for (int n=0; n<POL_DIM; ++n) {
         center_X[n] = m_realPoloidalCoords(iv,n);
      }
      POL::RealVect diff = center_X;
      diff -= X_pol;
      double dist = diff.vectorLength();
      if (dist < max_dist) {
         xi_initial = center_xi;
         max_dist = dist;
         iv_closest = iv;
      }
   }

   return m_poloidal_util->mappedCoordNewton(a_X, xi_initial, iv_closest);
   
#endif

#if CFG_DIM==3
   RealVect xi_initial = RealVect::Zero;

   RealVect rc_dx = m_dx * m_rc_coarsen_ratio;
   RealVect offset = 0.5*RealVect::Unit;
   offset *= rc_dx;

   IntVect iv_closest;

   for (BoxIterator bit(m_realCoords.box()); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect center_xi = rc_dx*iv + offset;
      RealVect center_X;
      for (int n=0; n<SpaceDim; ++n) {
         center_X[n] = m_realCoords(iv,n);
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

   double R;
   double Phi;
   double Z = a_X[2];
   
   if ( m_phys_coord_type == CYLINDRICAL ) {
      R = a_X[RADIAL_DIR];
      Phi = a_X[TOROIDAL_DIR];
   }
   else if ( m_phys_coord_type == CARTESIAN ) {
      R = sqrt(a_X[0]*a_X[0] + a_X[1]*a_X[1]);
      Phi = atan2(a_X[1],a_X[0]);
   }
   else {
      MayDay::Error("OneBlockCoordSys::mappedCoord(): unrecognized m_phys_coord_type");
   }
   xi_initial[1] = Phi + m_toroidal_mapped_offset;

   RealVect xi = mappedCoordNewton(RealVect(R,Phi,Z), xi_initial, iv_closest);
   
   return xi;

#endif
}


Real
OneBlockCoordSys::dXdXi(const RealVect&  a_Xi,
                        int              a_dirX,
                        int              a_dirXi ) const
{
   CH_assert(a_dirX>=0 && a_dirX<SpaceDim);
   CH_assert(a_dirXi>=0 && a_dirXi<SpaceDim);

   if ( !m_subgrid_geometry ) {
#if CFG_DIM==2
      return m_poloidal_util->interpolate(a_dirX, a_dirXi+1, a_Xi);
#endif
#if CFG_DIM==3
      Real dXdxi;
      
      double R   = m_interp->interpolate(0, Interp3D::FUNCTION, a_Xi);
      double Phi = m_interp->interpolate(1, Interp3D::FUNCTION, a_Xi) - m_toroidal_mapped_offset;

      if ( m_phys_coord_type == CARTESIAN ) {
         if ( a_dirX == 0 ) {
            if ( a_dirXi == 0 ) {
               dXdxi = m_interp->interpolate(0, Interp3D::X_DERIVATIVE, a_Xi) * cos(Phi);
            }
            else if ( a_dirXi == 1 ) {
               dXdxi = m_interp->interpolate(0, Interp3D::Y_DERIVATIVE, a_Xi) * cos(Phi) - R * sin(Phi);
            }
            else {
               dXdxi = m_interp->interpolate(0, Interp3D::Z_DERIVATIVE, a_Xi) * cos(Phi);
            }
         }
         else if ( a_dirX == 1 ) {
            if ( a_dirXi == 0 ) {
               dXdxi = m_interp->interpolate(0, Interp3D::X_DERIVATIVE, a_Xi) * sin(Phi);
            }
            else if ( a_dirXi == 1 ) {
               dXdxi = m_interp->interpolate(0, Interp3D::Y_DERIVATIVE, a_Xi) * sin(Phi) + R * cos(Phi);
            }
            else {
               dXdxi = m_interp->interpolate(0, Interp3D::Z_DERIVATIVE, a_Xi) * sin(Phi);
            }
         }
         else if ( a_dirX == 2 ) {
            if ( a_dirXi == 0 ) {
               dXdxi = m_interp->interpolate(2, Interp3D::X_DERIVATIVE, a_Xi);
            }
            else if ( a_dirXi == 1 ) {
               dXdxi = m_interp->interpolate(2, Interp3D::Y_DERIVATIVE, a_Xi);
            }
            else {
               dXdxi = m_interp->interpolate(2, Interp3D::Z_DERIVATIVE, a_Xi);
            }
         }
      }
      else if ( m_phys_coord_type == CYLINDRICAL ) {
         if ( a_dirXi == 0 ) {
            dXdxi = m_interp->interpolate(a_dirX, Interp3D::X_DERIVATIVE, a_Xi);
         }
         else if ( a_dirXi == 1 ) {
            dXdxi = m_interp->interpolate(a_dirX, Interp3D::Y_DERIVATIVE, a_Xi);
         }
         else {
            dXdxi = m_interp->interpolate(a_dirX, Interp3D::Z_DERIVATIVE, a_Xi);
         }
      }
      else {
         MayDay::Error("OneBlockCoordSys::dXdXi(): unrecognized m_phys_coord_type");
      }

      return dXdxi;
#endif
   }
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
OneBlockCoordSys::dXdXi(FArrayBox&       a_dXdXi,
                        const FArrayBox& a_Xi,
                        int              a_destComp,
                        int              a_dirX,
                        int              a_dirXi,
                        const Box&       a_box ) const
{
   for (BoxIterator bit(a_box); bit.ok(); ++bit ) {
      IntVect iv = bit();

      RealVect xi;
      for (int n=0; n<SpaceDim; ++n) {
         xi[n] = a_Xi(iv,n);
      }

      a_dXdXi(iv,a_destComp) = dXdXi(xi, a_dirX, a_dirXi);
   }
}


bool
OneBlockCoordSys::isValid(const RealVect& a_xi,
                          const bool      a_omit_toroidal ) const
{
   bool valid = a_xi[RADIAL_DIR]   >= lowerMappedCoordinate(RADIAL_DIR)   &&
                a_xi[RADIAL_DIR]   <= upperMappedCoordinate(RADIAL_DIR)   &&
                a_xi[POLOIDAL_DIR] >= lowerMappedCoordinate(POLOIDAL_DIR) &&
                a_xi[POLOIDAL_DIR] <= upperMappedCoordinate(POLOIDAL_DIR);

#if CFG_DIM==3
   if ( !a_omit_toroidal ) {
      valid &=  a_xi[TOROIDAL_DIR] >= lowerMappedCoordinate(TOROIDAL_DIR) &&
                a_xi[TOROIDAL_DIR] <= upperMappedCoordinate(TOROIDAL_DIR);
   }
#endif

   return valid;
}


void
OneBlockCoordSys::getMagneticFlux(const FArrayBox&  a_physical_coordinates,
                                  FArrayBox&         a_magnetic_flux ) const
{
   if (m_spectral_field) {
#if CFG_DIM==2
      m_poloidal_util->getMagneticFluxFromDCT(a_physical_coordinates, a_magnetic_flux);
#endif
#if CFG_DIM==3
      POL::FArrayBox physical_coordinates_pol(restrictToPoloidal(a_physical_coordinates.box()), SpaceDim);
      restrictPhysCoordsToPoloidal(a_physical_coordinates, physical_coordinates_pol);

      const Box& box = a_magnetic_flux.box();

      POL::FArrayBox magnetic_flux_pol(restrictToPoloidal(box), 1);

      for (int i=box.smallEnd(TOROIDAL_DIR); i<=box.bigEnd(TOROIDAL_DIR); ++i) {
         restrictPhysCoordsToPoloidal(a_physical_coordinates, i, physical_coordinates_pol);

         m_poloidal_util->getMagneticFluxFromDCT(physical_coordinates_pol, magnetic_flux_pol);

         injectToroidal(magnetic_flux_pol, i, a_magnetic_flux);
      }
#endif
   }
}


double
OneBlockCoordSys::getMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   double psi;

   if (m_spectral_field) {
      psi = m_poloidal_util->getMagneticFluxFromDCT(restrictPhysCoordToPoloidal(a_physical_coordinate));
   }
   else {
      MayDay::Error("OneBlockCoordSys::getMagneticFlux() only implemented for DCT field option");
      psi = 0.;
   }

   return psi;
}


void
OneBlockCoordSys::getNormMagneticFlux(const FArrayBox& a_physical_coordinates,
                                      FArrayBox&       a_magnetic_flux ) const

{
   double physFluxOnAxis = getMagneticFlux(injectToroidal(m_poloidal_util->getMagAxis()));
   double physFluxOnSep = getMagneticFlux(injectToroidal(m_poloidal_util->getXpoint()));
   
   getMagneticFlux(a_physical_coordinates, a_magnetic_flux);
   a_magnetic_flux.plus(-physFluxOnAxis);
   a_magnetic_flux.divide(physFluxOnSep - physFluxOnAxis);
   
}


double
OneBlockCoordSys::getNormMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   double physFluxOnAxis = getMagneticFlux(injectToroidal(m_poloidal_util->getMagAxis()));
   double physFluxOnSep = getMagneticFlux(injectToroidal(m_poloidal_util->getXpoint()));
   double physFlux = getMagneticFlux(a_physical_coordinate);
   
   return (physFlux - physFluxOnAxis)/(physFluxOnSep - physFluxOnAxis);
}


double
OneBlockCoordSys::getOuterMidplaneCoord( const double a_psiNorm ) const
{
   bool residual_tolerance_satisfied = false;
   int max_iterations = 10000;
   double tol = 1.0e-5;
   
   RealVect magAxis = injectToroidal(m_poloidal_util->getMagAxis());
   double tmp_lo_R = magAxis[RADIAL_DIR];
   
   double tmp_hi_R, half_R, f_lo_R, f_hi_R, f_half_R;
   
   if (m_spectral_field) {
      tmp_hi_R = m_poloidal_util->getRmax(); // MAYBE WILL NEED TO FIX LATER IN CASE THERE ARE MULTIVALUE SOLUTIONS NEAR M_RMAX
   }
   else MayDay::Error("OneBlockCoordSys::getOuterMidplaneCoord only implemented for DCT field option");
   
   int iter = 0;
   while ( iter++ < max_iterations && !residual_tolerance_satisfied ) {
      
      half_R = 0.5 * (tmp_lo_R + tmp_hi_R);

#if CH_SPACEDIM == 2
      RealVect lo_R_vect(tmp_lo_R,magAxis[POLOIDAL_DIR]);
      RealVect hi_R_vect(tmp_hi_R,magAxis[POLOIDAL_DIR]);
      RealVect half_R_vect(half_R,magAxis[POLOIDAL_DIR]);
#endif

#if CH_SPACEDIM == 3

      //Assume that the third direction is a toroidal angle
      //and that the geoemtry is axisymmetric
      double phi = 0.0;
      
      RealVect lo_R_vect(tmp_lo_R, phi, magAxis[POLOIDAL_DIR]);
      RealVect hi_R_vect(tmp_hi_R, phi, magAxis[POLOIDAL_DIR]);
      RealVect half_R_vect(half_R, phi, magAxis[POLOIDAL_DIR]);
#endif
      
      f_lo_R = getNormMagneticFlux(lo_R_vect) - a_psiNorm;
      f_hi_R = getNormMagneticFlux(hi_R_vect) - a_psiNorm;
      f_half_R = getNormMagneticFlux(half_R_vect) - a_psiNorm;
      
      if (f_lo_R * f_half_R < 0) {
         tmp_hi_R = half_R;
      }
      
      else {
         tmp_lo_R = half_R;
      }

#if CH_SPACEDIM == 2
      RealVect R_hi_test(half_R + tol,magAxis[POLOIDAL_DIR]);
      RealVect R_lo_test(half_R - tol,magAxis[POLOIDAL_DIR]);
#endif

#if CH_SPACEDIM == 3
      RealVect R_hi_test(half_R + tol, phi, magAxis[POLOIDAL_DIR]);
      RealVect R_lo_test(half_R - tol, phi, magAxis[POLOIDAL_DIR]);
#endif

      if ( (getNormMagneticFlux(R_hi_test) - a_psiNorm) * (getNormMagneticFlux(R_lo_test) - a_psiNorm) < 0 ) {
         residual_tolerance_satisfied = true;
      }
      
   }
   
   return half_R;
}


double
OneBlockCoordSys::getOuterRsep() const
{
   double psiNorm_atSep = getNormMagneticFlux(injectToroidal(m_poloidal_util->getXpoint()));
   
   return getOuterMidplaneCoord(psiNorm_atSep);
}


void
OneBlockCoordSys::computeFieldData(const int  a_dir,
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
   box    &= a_BFieldMag.box();
   box    &= a_BFieldDir.box();
   box    &= a_gradBFieldMag.box();
   box    &= a_curlBFieldDir.box();
   box    &= a_BFieldDirdotcurlBFieldDir.box();
   
#if CFG_DIM==2

   if ( m_spectral_field ) {

      FArrayBox RZ(box,SpaceDim);
      if (a_dir>=0 && a_dir<SpaceDim) {
         getFaceCenteredRealCoords(a_dir, RZ);
      }
      else {
         getCellCenteredRealCoords(RZ);
      }

      m_poloidal_util->computeFieldDataSpectral(RZ, a_dir, a_BField, a_BFieldMag, a_BFieldDir, a_gradBFieldMag,
                                                 a_curlBFieldDir, a_BFieldDirdotcurlBFieldDir, a_derived_data_only);
   }
   else {

      FArrayBox Xi(box,SpaceDim);
      if (a_dir>=0 && a_dir<SpaceDim) {
         getFaceCenteredMappedCoords(a_dir, Xi);
      }
      else {
         getCellCenteredMappedCoords(Xi);
      }

      m_poloidal_util->computeFieldDataFromMappingFile(Xi, a_dir, a_BField, a_BFieldMag,
                                                        a_BFieldDir, a_gradBFieldMag, a_curlBFieldDir,
                                                        a_BFieldDirdotcurlBFieldDir, a_derived_data_only );
   }
#endif

#if CFG_DIM==3

   // Create poloidal arrays
   POL::FArrayBox BField_pol(restrictToPoloidal(a_BField.box()), a_BField.nComp());
   POL::FArrayBox BFieldMag_pol(restrictToPoloidal(a_BFieldMag.box()), a_BFieldMag.nComp());;
   POL::FArrayBox BFieldDir_pol(restrictToPoloidal(a_BFieldDir.box()), a_BFieldDir.nComp());;
   POL::FArrayBox gradBFieldMag_pol(restrictToPoloidal(a_gradBFieldMag.box()), a_gradBFieldMag.nComp());;
   POL::FArrayBox curlBFieldDir_pol(restrictToPoloidal(a_curlBFieldDir.box()), a_curlBFieldDir.nComp());;
   POL::FArrayBox BFieldDirdotcurlBFieldDir_pol(restrictToPoloidal(a_BFieldDirdotcurlBFieldDir.box()), a_BFieldDirdotcurlBFieldDir.nComp());;

   FArrayBox X(box,SpaceDim);
   if (a_dir>=0 && a_dir<SpaceDim) {
      getFaceCenteredRealCoords(a_dir, X);
   }
   else {
      getCellCenteredRealCoords(X);
   }

   if ( m_spectral_field ) {

      POL::FArrayBox RZ_pol(restrictToPoloidal(X.box()), POL_DIM);

      for (int toroidal_index=box.smallEnd(TOROIDAL_DIR); toroidal_index<=box.bigEnd(TOROIDAL_DIR); ++toroidal_index) {
         restrictPhysCoordsToPoloidal(X, toroidal_index, RZ_pol);

         m_poloidal_util->computeFieldDataSpectral(RZ_pol, a_dir, BField_pol, BFieldMag_pol, BFieldDir_pol, gradBFieldMag_pol,
                                                   curlBFieldDir_pol, BFieldDirdotcurlBFieldDir_pol, a_derived_data_only);

         if( !a_derived_data_only ) {
            injectToroidal(BField_pol, toroidal_index, a_BField);
         }
         injectToroidal(BFieldMag_pol, toroidal_index, a_BFieldMag);
         injectToroidal(BFieldDir_pol, toroidal_index, a_BFieldDir);
         injectToroidal(gradBFieldMag_pol, toroidal_index, a_gradBFieldMag);
         injectToroidal(curlBFieldDir_pol, toroidal_index, a_curlBFieldDir);
         injectToroidal(BFieldDirdotcurlBFieldDir_pol, toroidal_index, a_BFieldDirdotcurlBFieldDir);
      }
   }
   else {
      MayDay::Error("OneBlockCoordSys::computeFieldData(): Only implemented in 3D for spectral field option");
   }

   if ( m_phys_coord_type == CARTESIAN ) {

      if( !a_derived_data_only ) {
         convertCylindricalToCartesian(X, a_BField);
      }
      convertCylindricalToCartesian(X, a_BFieldDir);
      convertCylindricalToCartesian(X, a_gradBFieldMag);
      convertCylindricalToCartesian(X, a_curlBFieldDir);
   }

#endif

}


array<double,3>
OneBlockCoordSys::computeBField(const RealVect& a_X) const
{
   if ( !m_spectral_field ) {
      MayDay::Error("OneBlockCoordSys::computeBField(): Only implemented for DCT field option");
   }

   POL::RealVect X_pol = restrictPhysCoordToPoloidal(a_X);

   array <double,3> result = m_poloidal_util->computeBField(X_pol);
   
#if CFG_DIM==3
   if ( m_phys_coord_type == CARTESIAN ) {
      convertCylindricalToCartesian(a_X, result);
   }
#endif
  
   return result;
}


void
OneBlockCoordSys::getNodalFieldData(FArrayBox& a_points,
                                    FArrayBox& a_A,
                                    FArrayBox& a_b,
                                    FArrayBox& a_Bmag) const
{
   if ( m_spectral_field ) {

      Box box = a_points.box();
      POL::Box box_pol = restrictToPoloidal(box);
      POL::FArrayBox RZ(box_pol, POL_DIM);
      
#if CFG_DIM==2

      for (POL::BoxIterator bit(box_pol); bit.ok(); ++bit) {
         POL::IntVect iv_pol = bit();
      
         RealVect xi;
         for (int n=0; n<SpaceDim; ++n) {
            xi[n] = a_points(iv_pol,n);
         }
         
         POL::RealVect X_pol = realCoord(xi);
         for (int n=0; n<POL_DIM; ++n) {
            RZ(iv_pol,n) = X_pol[n];
         }
      }

      m_poloidal_util->getNodalFieldData(RZ, a_A, a_b, a_Bmag);

#elif CFG_DIM==3

      POL::FArrayBox A_pol(restrictToPoloidal(a_A.box()),a_A.nComp());
      POL::FArrayBox b_pol(restrictToPoloidal(a_b.box()),a_b.nComp());
      POL::FArrayBox Bmag_pol(restrictToPoloidal(a_Bmag.box()),a_Bmag.nComp());
      
      // Loop over poloidal planes
      for (int toroidal_index=box.smallEnd(TOROIDAL_DIR); toroidal_index<=box.bigEnd(TOROIDAL_DIR); ++toroidal_index) {

         for (POL::BoxIterator bit(box_pol); bit.ok(); ++bit) {
            POL::IntVect iv_pol = bit();
            IntVect iv(iv_pol[0],toroidal_index,iv_pol[1]);

            RealVect xi;
            for (int n=0; n<SpaceDim; ++n) {
               xi[n] = a_points(iv,n);
            }
         
            POL::RealVect X_pol = restrictPhysCoordToPoloidal(realCoord(xi));
            for (int n=0; n<POL_DIM; ++n) {
               RZ(iv_pol,n) = X_pol[n];
            }

            m_poloidal_util->getNodalFieldData(RZ, A_pol, b_pol, Bmag_pol);

            injectToroidal(A_pol, toroidal_index, a_A);
            injectToroidal(b_pol, toroidal_index, a_b);
            injectToroidal(Bmag_pol, toroidal_index, a_Bmag);
         }
      }
#endif

   }
   else {
      MayDay::Error("OneBlockCoordSys::getNodalFieldData(): Only implemented for DCT field option");
   }
}


#if CFG_DIM==3
void
OneBlockCoordSys::setInterp(const ParmParse&       a_pp,
                            const POL::FArrayBox&  a_interp_node_coords_pol,
                            const POL::FArrayBox&  a_RZ_data )
{
   // Make an index space for the interpolation points by expanding the poloidal index
   // space in the toroidal direction.  In the new toroidal direction, the middle n_toroidal
   // points span the valid block.  There are an additional m_toroidal_ghosts points in the
   // lower and upper toroidal directions to provide a ghost cell region.

   CH_TIME("OneBlockCoordSys::setInterp");
  
   const POL::Box& box_pol = a_interp_node_coords_pol.box();
   POL::IntVect lo = box_pol.smallEnd();
   POL::IntVect hi = box_pol.bigEnd();

   int n_toroidal = m_toroidal_mapping_refinement * m_domain.domainBox().size(TOROIDAL_DIR) + 1;
   int toroidal_ghosts = m_toroidal_mapping_refinement * m_toroidal_ghosts;
   
   Box box(IntVect(lo[0],-toroidal_ghosts,lo[1]), IntVect(hi[0],n_toroidal+toroidal_ghosts-1,hi[1]));

   FArrayBox interp_node_coords(box,3);
   FArrayBox interp_node_data(box,3);
      
   double toroidal_node_spacing
      = (upperMappedCoordinate(TOROIDAL_DIR) - lowerMappedCoordinate(TOROIDAL_DIR)) / (n_toroidal - 1);
   
   if ( m_field_aligned_mapping ) {

      // The following code generates a 3D grid from a poloidal grid by tracing field lines.
      // Since this can be expensive (depending upon the grid size and requested ray tracing
      // step size), we provide an option to store and re-read this data.

      bool read_3D_mapping_data = false;
      if ( a_pp.contains("read_3D_mapping_data")) {
         a_pp.get("read_3D_mapping_data", read_3D_mapping_data);
      }

      string file_name;
      if ( a_pp.contains("3D_mapping_data_file")) {
         a_pp.get("3D_mapping_data_file", file_name);
      }

      if ( read_3D_mapping_data && file_name.empty() ) {
         MayDay::Error("OneBlockCoordSys::setInterp(): No 3D mapping data file specified");
      }

      fstream file_stream;

      if ( !read_3D_mapping_data ) {
      
         if ( procID()==0 ) {
            cout << "Constructing 3D mapping data in toroidal sector "
                  << m_toroidal_sector << endl;
         }

         bool centered_poloidal_mapping = true;
         if ( a_pp.contains("centered_poloidal_mapping")) {
            a_pp.get("centered_poloidal_mapping", centered_poloidal_mapping);
         }

         int base_index;
         double base_coord;
         if ( centered_poloidal_mapping ) {
            CH_assert( (n_toroidal-1)%2 == 0 );  // n_toroidal is assumed to be odd for this option
            base_index = (n_toroidal-1)/2 + 1;
            base_coord = 0.5 * (lowerMappedCoordinate(TOROIDAL_DIR) + upperMappedCoordinate(TOROIDAL_DIR));
         }
         else {
            base_index = 0;
            base_coord = lowerMappedCoordinate(TOROIDAL_DIR);
         }

         // On the index space constructed above, the node coordinates are obtained by starting
         // from the physical (cylindrical) coordinates at the toroidal center of the block if
         // centered_poloidal_mapping = true; otherwise, the lower toroidal block face.

         Box start_box(IntVect(lo[0],base_index,lo[1]), IntVect(hi[0],base_index,hi[1]));
         for (BoxIterator bit(start_box); bit.ok(); ++bit) {
            IntVect iv = bit();
            POL::IntVect iv_pol = restrictToPoloidal(iv);

            interp_node_coords(iv,0) = a_interp_node_coords_pol(iv_pol,0);
            interp_node_coords(iv,1) = base_coord;
            interp_node_coords(iv,2) = a_interp_node_coords_pol(iv_pol,1);

            interp_node_data(iv,0) = a_RZ_data(iv_pol,0);
            interp_node_data(iv,1) = base_coord;
            interp_node_data(iv,2) = a_RZ_data(iv_pol,1);
         }

         // The remaining node coordinates are obtained by tracing the field starting
         // from the coordinates just set.  Since we are using a BoxIterator, which
         // iterates in a positive index direction, we need to do this in two steps,
         // corresponding to the indices less than and greater than the indices
         // defining the base of the ray tracing.

         // To speed-up ray tracing we decompose the poloidal cross-section
         DisjointBoxLayout grids;
         Vector<Box> grid_boxes;
         int decomp_num;
         getPoloidalDisjointBoxLayout(grids, grid_boxes, decomp_num, box);
         
         LevelData<FArrayBox> interp_node_data_decomp(grids, 3);
         LevelData<FArrayBox> interp_node_coords_decomp(grids, 3);
         
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            interp_node_data_decomp[dit].copy(interp_node_data, grids[dit]);
            interp_node_coords_decomp[dit].copy(interp_node_coords, grids[dit]);
         }

         Box lower_box(IntVect(lo[0],-toroidal_ghosts,lo[1]), IntVect(hi[0],base_index-1,hi[1]));
         for (DataIterator dit(grids); dit.ok(); ++dit) {

            Box this_lower_box = grids[dit];
            this_lower_box &= lower_box;
            
            FArrayBox& this_interp_node_data = interp_node_data_decomp[dit];
            FArrayBox& this_interp_node_coords = interp_node_coords_decomp[dit];
            
            for (BoxIterator bit(this_lower_box); bit.ok(); ++bit) {
               IntVect iv = bit();
               iv[1] = - (toroidal_ghosts + iv[1] + 1) + base_index;

               POL::IntVect iv_pol = restrictToPoloidal(iv);
               this_interp_node_coords(iv,0) = a_interp_node_coords_pol(iv_pol,0);
               this_interp_node_coords(iv,1) = base_coord + (iv[1] - base_index)*toroidal_node_spacing;
               this_interp_node_coords(iv,2) = a_interp_node_coords_pol(iv_pol,1);

               RealVect X_old;
               IntVect iv_old = iv; iv_old[1]++;
               for (int n=0; n<3; ++n) {
                  X_old[n] = this_interp_node_data(iv_old,n);
               }

               RealVect X_new = traceField(X_old, -toroidal_node_spacing);

               for (int n=0; n<3; ++n) {
                  this_interp_node_data(iv,n) = X_new[n];
               }
            }
         }

         Box upper_box(IntVect(lo[0],base_index+1,lo[1]), IntVect(hi[0],n_toroidal+toroidal_ghosts-1,hi[1]));

         for (DataIterator dit(grids); dit.ok(); ++dit) {

            Box this_upper_box = grids[dit];
            this_upper_box &= upper_box;
            
            FArrayBox& this_interp_node_data = interp_node_data_decomp[dit];
            FArrayBox& this_interp_node_coords = interp_node_coords_decomp[dit];

            for (BoxIterator bit(this_upper_box); bit.ok(); ++bit) {
               IntVect iv = bit();

               POL::IntVect iv_pol = restrictToPoloidal(iv);
               this_interp_node_coords(iv,0) = a_interp_node_coords_pol(iv_pol,0);
               this_interp_node_coords(iv,1) = base_coord + (iv[1] - base_index)*toroidal_node_spacing;
               this_interp_node_coords(iv,2) = a_interp_node_coords_pol(iv_pol,1);

               RealVect X_old;
               IntVect iv_old = iv; iv_old[1]--;
               for (int n=0; n<3; ++n) {
                  X_old[n] = this_interp_node_data(iv_old,n);
               }

               RealVect X_new = traceField(X_old, toroidal_node_spacing);

               for (int n=0; n<3; ++n) {
                  this_interp_node_data(iv,n) = X_new[n];
               }
            }
         }

#ifdef CH_MPI
         
         assembleDecomposedData(interp_node_data, interp_node_coords, interp_node_data_decomp, interp_node_coords_decomp, grid_boxes, decomp_num);
         
         if ( procID() == 0 ) {
#endif
            if ( !file_name.empty() ) {

               file_stream.open(file_name.c_str(), ios::app);

               cout << "Adding 3D mapping data to file " << file_name << endl << endl;
               
               file_stream.write((char*)&m_poloidal_block, sizeof(int));
               file_stream.write((char*)&m_toroidal_sector, sizeof(int));
               int data_length = 6 * box.numPts();
               file_stream.write((char*)&data_length, sizeof(int));

               for (BoxIterator bit(box); bit.ok(); ++bit) {
                  IntVect iv = bit();

                  for (int n=0; n<3; ++n) {
                     file_stream.write((char*)&interp_node_coords(iv,n), sizeof(double));
                     file_stream.write((char*)&interp_node_data(iv,n), sizeof(double));
                  }
               }

               file_stream.close();
            }
#ifdef CH_MPI
         }
#endif

      }
      else {

         file_stream.open(file_name.c_str(), ios::in);

         unsigned buffer_length = 6 * box.numPts();
         double* buffer = new double[buffer_length];

#ifdef CH_MPI
         if ( procID() == 0 ) {
#endif
            bool found_data = false;
            
            while( !file_stream.eof() ) {
               int poloidal_block;
               file_stream.read((char*)&poloidal_block, sizeof(int));

               int toroidal_sector;
               file_stream.read((char*)&toroidal_sector, sizeof(int));

               int data_length;
               file_stream.read((char*)&data_length, sizeof(int));

               if ( poloidal_block == m_poloidal_block &&
                    toroidal_sector == m_toroidal_sector ) {

                  if ( data_length != buffer_length ) {
                     MayDay::Error("OneBlockCoordSys::setInterp(): interpolation data for block is wrong size");
                  }

                  found_data = true;
                  break;
               }
               else {
                  // Skip the rest of this block
                  file_stream.seekg(data_length * sizeof(double), ios::cur);
               }
            }

            if ( !found_data ) {
               MayDay::Error("OneBlockCoordSys::setInterp(): interpolation data for block was not found");
            }
            
            cout << "Reading interpolation data from file " << file_name << " in toroidal sector " << m_toroidal_sector << endl << endl;

            file_stream.read((char*)buffer, buffer_length * sizeof(double));

            file_stream.close();

#ifdef CH_MPI
         }

         MPI_Bcast(buffer, buffer_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

         double* buffer_ptr = buffer;
         for (BoxIterator bit(box); bit.ok(); ++bit) {
            IntVect iv = bit();

            for (int n=0; n<3; ++n) {
               interp_node_coords(iv,n) = *buffer_ptr++;
               interp_node_data(iv,n)   = *buffer_ptr++;
            }
         }

         delete [] buffer;
      }
   }
   else {

      for (BoxIterator bit(box); bit.ok(); ++bit) {
         IntVect iv = bit();
         POL::IntVect iv_pol = restrictToPoloidal(iv);
      
         interp_node_coords(iv,0) = a_interp_node_coords_pol(iv_pol,0);
         interp_node_coords(iv,1) = lowerMappedCoordinate(TOROIDAL_DIR) + iv[TOROIDAL_DIR]*toroidal_node_spacing;
         interp_node_coords(iv,2) = a_interp_node_coords_pol(iv_pol,1);

         interp_node_data(iv,0) = a_RZ_data(iv_pol,0);
         interp_node_data(iv,1) = interp_node_coords(iv,1);
         interp_node_data(iv,2) = a_RZ_data(iv_pol,1);
      }
   }
   
   m_interp = new BSplineInterp3D(a_pp, interp_node_coords, interp_node_data);
}

RealVect
OneBlockCoordSys::traceField(const RealVect&  a_X,
                             const double&    a_toroidal_increment ) const
{

   CH_TIME("OneBlockCoordSys::traceField");
   
   double step = a_toroidal_increment / double(m_field_trace_step_num);
   RealVect X = a_X;
   RealVect dX;
   
   for (int n=0; n<m_field_trace_step_num; n++) {

      // Get the field unit vector
      array<double,3> b = m_poloidal_util->computeBUnit(restrictToPoloidal(X));
      
      // Set the tracing direction
      RealVect direction;
      for (int n=0; n<3; ++n) {
         direction[n] = b[n];
      }
      if ( b[TOROIDAL_DIR] < 0. ) direction *= -1.;

      // Get the step increment
      dX[TOROIDAL_DIR] = step;
      dX[RADIAL_DIR] = direction[RADIAL_DIR]/direction[TOROIDAL_DIR] * step * X[RADIAL_DIR];
      dX[POLOIDAL_DIR] = direction[POLOIDAL_DIR]/direction[TOROIDAL_DIR] * step * X[RADIAL_DIR];
      
      X += dX;
   }

   return X;
}

void
OneBlockCoordSys::getPoloidalDisjointBoxLayout(DisjointBoxLayout&   a_grids,
                                               Vector<Box>&         a_boxes,
                                               int&                 a_decomp_num,
                                               const Box&           a_box) const
{

   int n_proc = numProc();

   int n_rad = a_box.size(RADIAL_DIR);
   int n_pol = a_box.size(POLOIDAL_DIR);


   // The actual number of cells per
   // direction we will end up using
   int n_loc_rad;
   int n_loc_pol;

   // Assumption for what we would like to have
   // for number of cells per direction
   int n_opt = 4;
   CH_assert(n_rad >= n_opt && n_pol >= n_opt);
   
   // Too many processors available (hope for that situation, since
   // we often have extra procs for velocity space)
   if ((n_rad / n_opt) * (n_pol / n_opt) <= n_proc) {
      n_loc_rad = n_opt;
      n_loc_pol = n_opt;
   }
   
   // Too few processors available: start with an estimate
   else {
      int num_iter = int(floor(log2(double(n_proc))));
      int num_iter_rad;
      //Begin with dividing radial direction
      for(int n = 0; n < num_iter+1; n++) {
         int test = int(floor(n_rad/pow(2,n)));
         if (test < n_opt) {
            n_loc_rad = int(floor(n_rad/pow(2,n-1)));
            num_iter_rad = n-1;
            break;
         }
         num_iter_rad = n;
         n_loc_rad = int(floor(n_rad/pow(2,n)));
      }

      //Now divide poloidal direction
      n_loc_pol = n_pol;
      for(int n = 0; n < num_iter - num_iter_rad + 1; n++) {
         int test = int(floor(n_pol/pow(2,n)));
         if (test  < n_opt) {
            n_loc_pol = int(floor(n_pol/pow(2,n-1)));
            break;
         }
         n_loc_pol = int(floor(n_pol/pow(2,n)));
      }
   }

   Vector<int> proc_ids;

   IntVect domain_hi_end = a_box.bigEnd();
   IntVect domain_lo_end = a_box.smallEnd();

   int n_rad_dec = n_rad / n_loc_rad;
   int n_pol_dec = n_pol / n_loc_pol;
   
   a_decomp_num = n_rad_dec * n_pol_dec;

   //Now do fine tuning
   if (a_decomp_num < numProc()) {
     while (a_decomp_num < numProc()) {
       n_loc_pol = n_loc_pol - 1;
       n_pol_dec = n_pol / n_loc_pol;
       a_decomp_num = n_rad_dec * n_pol_dec;
     }
   }

   if (a_decomp_num > numProc()) {
     while (a_decomp_num > numProc()) {
       n_loc_rad = n_loc_rad + 1;
       n_rad_dec = n_rad / n_loc_rad;
       a_decomp_num = n_rad_dec * n_pol_dec;
     }
   }

   for(int j = 0; j < n_pol_dec; j++) {
      for(int i = 0; i < n_rad_dec; i++) {

         int rad_ind_lo = (domain_lo_end[RADIAL_DIR] + i * n_loc_rad);
         int rad_ind_hi = (domain_lo_end[RADIAL_DIR] + (i+1) * n_loc_rad - 1);
            
         if (i == n_rad_dec - 1) {
            rad_ind_hi = domain_hi_end[RADIAL_DIR];
         }

         int pol_ind_lo = (domain_lo_end[POLOIDAL_DIR] + j * n_loc_pol);
         int pol_ind_hi = (domain_lo_end[POLOIDAL_DIR] + (j+1) * n_loc_pol - 1);
            
         if (j == n_pol_dec - 1) {
            pol_ind_hi = domain_hi_end[POLOIDAL_DIR];
         }
            
         IntVect lo_end(rad_ind_lo, domain_lo_end[TOROIDAL_DIR], pol_ind_lo);
         IntVect hi_end(rad_ind_hi, domain_hi_end[TOROIDAL_DIR], pol_ind_hi);
            
         Box box(lo_end, hi_end);
         int procID = j * n_rad_dec + i;

         a_boxes.push_back( box );
         proc_ids.push_back( procID );
      }
   }
   
   a_grids.define( a_boxes, proc_ids );

}

void
OneBlockCoordSys::assembleDecomposedData(FArrayBox&                    a_interp_node_data,
                                         FArrayBox&                    a_interp_node_coords,
                                         const LevelData<FArrayBox>&   a_interp_node_data_decomp,
                                         const LevelData<FArrayBox>&   a_interp_node_coords_decomp,
                                         const Vector<Box>&            a_grid_boxes,
                                         const int&                    a_decomp_num ) const
{
   // We assume that each proc contains only one box (we created grids that way)
   // Also, the processor IDs assigned to grids correspond to 0 .. a_decomp_num-1

   CH_TIME("OneBlockCoordSys::assembleDecomposedData");

   CH_assert(numProc() >= a_decomp_num);

   FArrayBox this_interp_node_data;
   FArrayBox this_interp_node_coords;
   
   Vector<FArrayBox> all_interp_node_data(numProc());
   Vector<FArrayBox> all_interp_node_coords(numProc());
   
   // Fill local data
   for (int ivec = 0; ivec < numProc(); ivec++) {
      if (procID() == ivec && ivec < a_decomp_num) {
         this_interp_node_data.define(a_grid_boxes[ivec], 3);
         this_interp_node_coords.define(a_grid_boxes[ivec], 3);
         
         // Iterate over boxes of decomposed object to find the box needed
         for (DataIterator dit( a_interp_node_data_decomp.dataIterator() ); dit.ok(); ++dit) {
            if (a_interp_node_data_decomp[dit].box() == a_grid_boxes[ivec]) {
               this_interp_node_data.copy(a_interp_node_data_decomp[dit]);
               this_interp_node_coords.copy(a_interp_node_coords_decomp[dit]);
       }
         }
      }
      else {
         if (procID() == ivec) {
            // Or perhas should use the empty constructor to create invalid FArrayBox
            Box empty_box;
            this_interp_node_data.define(empty_box,1);
            this_interp_node_coords.define(empty_box, 1);
         }
      }
   }
   
   // Gather data across all processors and put it onto all processors
   for (int ivec = 0; ivec < numProc(); ivec++) {
      gather(all_interp_node_data, this_interp_node_data, ivec);
      gather(all_interp_node_coords, this_interp_node_coords, ivec);
   }
   
   // Create global interpolation FArrayBox objects on all procs
   for (int ivec = 0; ivec < a_decomp_num; ivec++) {
      a_interp_node_data.copy(all_interp_node_data[ivec], a_grid_boxes[ivec]);
      a_interp_node_coords.copy(all_interp_node_coords[ivec], a_grid_boxes[ivec]);
   }

   inspect(a_interp_node_data);
}


RealVect
OneBlockCoordSys::mappedCoordNewton(const RealVect&  a_X,
                                    const RealVect&  a_xi_initial,
                                    const IntVect&   a_iv_initial ) const
{
   /*
     Use Newton iteration to evaluate the mapping of physical to computational
     coordinates by inverting the mapping of computational to physical coordinates.
     N.B: The physical coordinate is assumed to be cylindrical.  This function is
     used to invert a coordinate mapping of the form:

     R   = R(xi_1,xi_2,xi_3)
     Phi = Phi(xi_2) = xi_2
     Z   = Z(xi_1,xi_2,xi_3)

   */
   double tol = 1.e-10;         // Read from input if we find some reason to later
   int max_iter = 40;

   RealVect xi = a_xi_initial;

   POL::RealVect residual;

   double R   = m_interp->interpolate(0, Interp3D::FUNCTION, xi);
   //   double Phi = xi[1];
   double Z   = m_interp->interpolate(2, Interp3D::FUNCTION, xi);
   double R1  = m_interp->interpolate(0, Interp3D::X_DERIVATIVE, xi);
   double R2  = m_interp->interpolate(0, Interp3D::Y_DERIVATIVE, xi);
   double R3  = m_interp->interpolate(0, Interp3D::Z_DERIVATIVE, xi);
   double Z1  = m_interp->interpolate(2, Interp3D::X_DERIVATIVE, xi);
   double Z2  = m_interp->interpolate(2, Interp3D::Y_DERIVATIVE, xi);
   double Z3  = m_interp->interpolate(2, Interp3D::Z_DERIVATIVE, xi);

   residual[0] = R - a_X[0];
   residual[1] = Z - a_X[2];

   double Fnorm = residual.vectorLength();

   bool converged = Fnorm <= tol;
   int num_iters = 0;
   int bt_steps;

   RealVect xi_saved[100];
   POL::RealVect residual_saved[100];
   double Fnorm_saved[100];

   xi_saved[0] = xi;
   Fnorm_saved[0] = Fnorm;
   residual_saved[0] = residual;

   while ( !converged && num_iters < max_iter) {

      // Invert the Jacobian to get the update
      double J = R1*Z3 - R3*Z1;

      RealVect delta;
      delta[0] = - (  Z3*residual[0] - R3*residual[1] ) / J;
      delta[1] = 0.;
      delta[2] = - ( -Z1*residual[0] + R1*residual[1] ) / J;
      
      double s = 1.;

      RealVect xi_proposed = xi + delta;

      R   = m_interp->interpolate(0, Interp3D::FUNCTION, xi_proposed);
      //      Phi = xi_proposed[1];
      Z   = m_interp->interpolate(2, Interp3D::FUNCTION, xi_proposed);
      R1  = m_interp->interpolate(0, Interp3D::X_DERIVATIVE, xi_proposed);
      R2  = m_interp->interpolate(0, Interp3D::Y_DERIVATIVE, xi_proposed);
      R3  = m_interp->interpolate(0, Interp3D::Z_DERIVATIVE, xi_proposed);
      Z1  = m_interp->interpolate(2, Interp3D::X_DERIVATIVE, xi_proposed);
      Z2  = m_interp->interpolate(2, Interp3D::Y_DERIVATIVE, xi_proposed);
      Z3  = m_interp->interpolate(2, Interp3D::Z_DERIVATIVE, xi_proposed);

      residual[0] = R - a_X[0];
      residual[1] = Z - a_X[2];

      double Fnorm_proposed = residual.vectorLength();

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

         R   = m_interp->interpolate(0, Interp3D::FUNCTION, xi_proposed);
         //         Phi = xi_proposed[1];
         Z   = m_interp->interpolate(2, Interp3D::FUNCTION, xi_proposed);
         R1  = m_interp->interpolate(0, Interp3D::X_DERIVATIVE, xi_proposed);
         R2  = m_interp->interpolate(0, Interp3D::Y_DERIVATIVE, xi_proposed);
         R3  = m_interp->interpolate(0, Interp3D::Z_DERIVATIVE, xi_proposed);
         Z1  = m_interp->interpolate(2, Interp3D::X_DERIVATIVE, xi_proposed);
         Z2  = m_interp->interpolate(2, Interp3D::Y_DERIVATIVE, xi_proposed);
         Z3  = m_interp->interpolate(2, Interp3D::Z_DERIVATIVE, xi_proposed);

         residual[0] = R - a_X[0];
         residual[1] = Z - a_X[2];

         Fnorm_proposed = residual.vectorLength();

         backtracking = (Fnorm_proposed > (1. - t*(1. - eta)) * Fnorm) && (bt_steps < bt_steps_max);

         bt_steps++;
      }

      // Update the current solution and residual norm
      xi = xi_proposed;
      Fnorm = Fnorm_proposed;

      // Test convergence
      converged = Fnorm <= tol;

      num_iters++;

      xi_saved[num_iters] = xi;
      Fnorm_saved[num_iters] = Fnorm;
      residual_saved[num_iters] = residual;
   }

#ifdef REPORT_NEWTON_FAILURE
   if ( !converged && num_iters >= max_iter && procID()==0) {
      cout << "On block (" << m_toroidal_sector << "): Newton solve did not converge at "
           << a_X << ", a_xi_initial = " << a_xi_initial << ", iv_init = " << a_iv_initial << ", Fnorm = " << Fnorm << endl;
      for (int i=0; i<num_iters; ++i) {
         cout << i << " xi = " << xi_saved[i] << ", Fnorm = " << Fnorm_saved[i] << " " << residual_saved[i] << endl;
      }
      exit(1);
   }
#endif

   return xi;
}
#endif


void
OneBlockCoordSys::convertCylindricalToCartesian(const RealVect&   a_cartesian_vect,
                                                array<double,3>&  a_vect) const
{
#if CFG_DIM==3
   // This function converts a_vect from the cylindrical to the cartesian frame
   // at the point a_cartesian_vect

   double x = a_cartesian_vect[0];
   double y = a_cartesian_vect[1];
   double r = sqrt(x*x + y*y);

   // Get the R,Phi components of the physical coordinate
   double F_R   = a_vect[0];
   double F_Phi = a_vect[1];

   // Replace them by the corresponding X,Y components
   a_vect[0] = (x/r) * F_R - (y/r) * F_Phi;
   a_vect[1] = (y/r) * F_R + (x/r) * F_Phi;
#endif
}


void
OneBlockCoordSys::convertCylindricalToCartesian(const FArrayBox&  a_cartesian_coord,
                                                FArrayBox&        a_coord ) const
{
#if CFG_DIM==3
   for (BoxIterator bit(a_cartesian_coord.box()); bit.ok(); ++bit) {
      IntVect iv = bit();

      array<double,3> vect;
      RealVect cartesian_vect;
      for (int dir=0; dir<SpaceDim; ++dir) {
         vect[dir]           = a_coord(iv,dir);
         cartesian_vect[dir] = a_cartesian_coord(iv,dir);
      }

      convertCylindricalToCartesian(cartesian_vect, vect);

      for (int dir=0; dir<SpaceDim; ++dir) {
         a_coord(iv,dir) = vect[dir];
      }
   }
#endif
}

POL::IntVect
OneBlockCoordSys::restrictToPoloidal( const IntVect& a_v ) const
{
#if CFG_DIM==2
   return a_v;
#endif
#if CFG_DIM==3
   return POL::IntVect(a_v[RADIAL_DIR], a_v[POLOIDAL_DIR]);
#endif
}


POL::RealVect
OneBlockCoordSys::restrictToPoloidal( const RealVect& a_v ) const
{
#if CFG_DIM==2
   return a_v;
#endif
#if CFG_DIM==3
   return POL::RealVect(a_v[RADIAL_DIR], a_v[POLOIDAL_DIR]);
#endif
}


POL::Box
OneBlockCoordSys::restrictToPoloidal( const Box& a_box ) const
{
#if CFG_DIM==2
   return a_box;
#endif
#if CFG_DIM==3
   return POL::Box(restrictToPoloidal(a_box.smallEnd()),restrictToPoloidal(a_box.bigEnd()));
#endif
}


POL::RealVect
OneBlockCoordSys::restrictPhysCoordToPoloidal( const RealVect& a_X ) const
{
   POL::RealVect X_pol;

#if CFG_DIM==2
   X_pol = a_X;
#endif
#if CFG_DIM==3
   if ( m_phys_coord_type == CYLINDRICAL ) {
      X_pol = restrictToPoloidal(a_X);
   }
   else if ( m_phys_coord_type == CARTESIAN ) {
      X_pol = POL::RealVect(sqrt(a_X[0]*a_X[0] + a_X[1]*a_X[1]),a_X[2]);
   }
   else {
      MayDay::Error("OneBlockCoordSys::restrictPhysCoordToPoloidal(): unrecognized m_phys_coord_type");
   }
#endif

   return X_pol;
}


void
OneBlockCoordSys::restrictPhysCoordsToPoloidal(const FArrayBox&  a_coords,
                                               POL::FArrayBox&   a_coords_pol ) const
{
#if CFG_DIM==2
   a_coords_pol.copy(a_coords);
#endif
#if CFG_DIM==3
   const Box& box = a_coords.box();
   int lo_toroidal = box.smallEnd(TOROIDAL_DIR);
   const POL::Box& box_pol = a_coords_pol.box();

   for (POL::BoxIterator bit(box_pol); bit.ok(); ++bit) {
      POL::IntVect iv_pol = bit();
      IntVect iv(iv_pol[0],lo_toroidal,iv_pol[1]);

      RealVect X;
      for (int n=0; n<CFG_DIM; ++n) {
         X[n] = a_coords(iv,n);
      }

      POL::RealVect X_pol = restrictPhysCoordToPoloidal(X);
      for (int n=0; n<POL_DIM; ++n) {
         a_coords_pol(iv_pol,n) = X_pol[n];
      }
   }
#endif
}


void
OneBlockCoordSys::restrictPhysCoordsToPoloidal(const FArrayBox&  a_coords,
                                               const int         a_toroidal_index,
                                               POL::FArrayBox&   a_coords_pol ) const
{
#if CFG_DIM==2
   a_coords_pol.copy(a_coords);
#endif
#if CFG_DIM==3
   const POL::Box& box_pol = a_coords_pol.box();

   for (POL::BoxIterator bit(box_pol); bit.ok(); ++bit) {
      POL::IntVect iv_pol = bit();
      IntVect iv(iv_pol[0],a_toroidal_index,iv_pol[1]);

      RealVect X;
      for (int n=0; n<CFG_DIM; ++n) {
         X[n] = a_coords(iv,n);
      }

      POL::RealVect X_pol = restrictPhysCoordToPoloidal(X);
      for (int n=0; n<POL_DIM; ++n) {
         a_coords_pol(iv_pol,n) = X_pol[n];
      }
   }
#endif
}


void
OneBlockCoordSys::restrictMappedCoordsToPoloidal(const FArrayBox&  a_coords,
                                                 POL::FArrayBox&   a_coords_pol ) const
{
#if CFG_DIM==2
   a_coords_pol.copy(a_coords);
#endif
#if CFG_DIM==3
   const Box& box = a_coords.box();
   int lo_toroidal = box.smallEnd(TOROIDAL_DIR);
   const POL::Box& box_pol = a_coords_pol.box();

   for (POL::BoxIterator bit(box_pol); bit.ok(); ++bit) {
      POL::IntVect iv_pol = bit();
      IntVect iv(iv_pol[0],lo_toroidal,iv_pol[1]);

      a_coords_pol(iv_pol,0) = a_coords(iv,RADIAL_DIR);
      a_coords_pol(iv_pol,1) = a_coords(iv,POLOIDAL_DIR);
   }
#endif
}


RealVect
OneBlockCoordSys::injectToroidal( const POL::RealVect  a_v ) const
{
  /*
    Converts a poloidal physical coordinate (R,Z) to
    (X,Y,Z) with Y=0
   */
#if CFG_DIM==2
   return a_v;
#endif
#if CFG_DIM==3
   //return RealVect(a_v[0], lowerMappedCoordinate(TOROIDAL_DIR), a_v[1]);
   return RealVect(a_v[0], 0.0, a_v[1]);
#endif
}


void
OneBlockCoordSys::injectToroidal(const POL::FArrayBox&  a_fab_pol,
                                 const int              a_toroidal_index,
                                 FArrayBox&             a_fab ) const
{
#if CFG_DIM==2
   a_fab.copy(a_fab_pol);
#endif
#if CFG_DIM==3
   for (POL::BoxIterator bit(a_fab_pol.box()); bit.ok(); ++bit) {
      POL::IntVect iv_pol = bit();
      IntVect iv(iv_pol[0], a_toroidal_index, iv_pol[1]);
      
      for (int n=0; n<a_fab_pol.nComp(); ++n) {
         a_fab(iv,n) = a_fab_pol(iv_pol,n);
      }
   }
#endif
}

#include "NamespaceFooter.H"
