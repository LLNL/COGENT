#include <float.h>
#include <fstream>
#include "SingleNullBlockCoordSys.H"
#include "MagBlockCoordSysF_F.H"
#include "SingleNullBlockCoordSysF_F.H"
#include "HermiteInterp.H"
#include "SplineInterp.H"
#include "BSplineInterp.H"
#include "Directions.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"

#undef REPORT_NEWTON_FAILURE

enum interpolation_package {HERMITE, SPLINE, BSPLINE};
enum interpolation_mode {FUNCTION, RADIAL_DERIVATIVE, POLOIDAL_DERIVATIVE, SECOND_RADIAL_DERIVATIVE, RADIAL_POLOIDAL_DERIVATIVE, SECOND_POLOIDAL_DERIVATIVE};
enum field_vars {RBR_VAR, RBZ_VAR, RB_VAR, UNIT_BR_VAR, UNIT_BPHI_VAR, UNIT_BZ_VAR, NUM_FIELD_VARS};
enum rz_vars {R_VAR, Z_VAR};

#include "NamespaceHeader.H"

const std::string SingleNullBlockCoordSys::pp_name = "singlenullRealisticGeom";

SingleNullBlockCoordSys::SingleNullBlockCoordSys( ParmParse&               a_parm_parse,
                                                  const ProblemDomain&     a_domain,
                                                  const RealVect&          a_dx,
                                                  const int                a_block_type )
   : MagBlockCoordSys(a_parm_parse),
     m_block_type(a_block_type),
     m_poloidally_truncated(false)
{
   if (SpaceDim != 2) {
      MayDay::Error("SingleNullBlockCoordSys is only two-dimensional");
   }

   if (m_verbose && procID()==0) {
      cout << "Constructing single null ";
      switch (a_block_type)
         {
         case LCORE:
           cout << "left core";
           break;
         case MCORE:
           cout << "middle core";
           break;
         case RCORE:
           cout << "right core";
           break;
         case LCSOL:
           cout << "left center scrape-off layer";
           break;
         case MCSOL:
           cout << "middle center scrape-off layer";
           break;
         case RCSOL:
           cout << "right center scrape-off layer";
           break;
         case LSOL:
           cout << "left scrape_off layer";
           break;
         case RSOL:
           cout << "right scrape-off layer";
           break;
         case LPF:
           cout << "left private flux";
           break;
         case RPF:
           cout << "right private flux";
           break;
         default:
           MayDay::Error("SingleNullBlockCoordSys::SingleNullBlockCoordSys(): Invalid block_type encountered");
         }

      cout << " block with global index space domain box = " << a_domain.domainBox() << endl;
   }

#ifdef PLOT_INVERSE_POINTS
   char file_name[80];
   sprintf(file_name, "Block%d_%d", m_block_type, procID());
   m_ipt_file[m_block_type] = fopen(file_name, "w");
#endif

   // Get the toroidal field component scale factor
   if (a_parm_parse.contains("Btor_scale")) {
     a_parm_parse.get("Btor_scale", m_RB_toroidal);
   }
   else {
     m_RB_toroidal = 0.;  // default
   }

   define( a_domain, a_dx );

   if (m_verbose && procID()==0) {
      switch (a_block_type)
         {
         case LCORE:
           cout << "Left core mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0) << ", "
                << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case MCORE:
           cout << "Middle core mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0) << ", "
                << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case RCORE:
            cout << "Right core mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0) << ", "
                 << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case LCSOL:
            cout << "Left center scrape-off layer mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case MCSOL:
            cout << "Middle center scrape-off layer mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case RCSOL:
            cout << "Right center scrape-off layer mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case LSOL:
            cout << "Left scrape_off layer mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case RSOL:
            cout << "Right scrape-off layer mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case LPF:
            cout << "Left private flux mapped domain: " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         case RPF:
            cout << "Right private flux mapped domain " << lowerMappedCoordinate(0) << " < xi_0 < " << upperMappedCoordinate(0)
                 << ", " << lowerMappedCoordinate(1) << " < xi_1 < " << upperMappedCoordinate(1);
            break;
         default:
           MayDay::Error("SingleNullBlockCoordSys::SingleNullBlockCoordSys(): Invalid block_type encountered");
         }

      cout << endl;
   }

   init( a_parm_parse );
}


SingleNullBlockCoordSys::~SingleNullBlockCoordSys()
{
#ifdef PLOT_INVERSE_POINTS
   for (int block=0; block<NUM_SINGLE_NULL_BLOCKS; ++block) {
      fclose(m_ipt_file[block]);
   }
#endif
}


void
SingleNullBlockCoordSys::getCellCenterRealCoords()
{
   Box box = m_domain.domainBox();;

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


void
SingleNullBlockCoordSys::init( ParmParse& a_pp )
{
  int n_radial_points, n_poloidal_points;
  int n_radial_extend, n_poloidal_extend;
  int num_block_elements;
  ifstream inFile;

#ifdef CH_MPI
  if (procID() == 0) {
#endif

    string geometry_file;
    if ( a_pp.contains("geometry_file")) {
      a_pp.get("geometry_file", geometry_file);
    }
    else {
      MayDay::Error("SingleNullBlockCoordSys::init(): No geometry mapping file specified");
    }

    inFile.open( geometry_file.c_str() );

    if (!inFile) {
      cout << "Unable to open geometry mapping file" << endl;
      exit(1);
    }

    cout << "Reading geometry mapping file" << endl;

    bool found_block = false;
    string block_name;

    while ( (inFile >> block_name) ) {

      inFile >> n_radial_points;
      inFile >> n_radial_extend;
      inFile >> n_poloidal_points;
      inFile >> n_poloidal_extend;

      num_block_elements = (n_radial_points + 2*n_radial_extend) * (n_poloidal_points + 2*n_poloidal_extend);

      if ( (block_name == "mcore" && m_block_type == MCORE) ||
           (block_name == "lcore" && m_block_type == LCORE) ||
           (block_name == "rcore" && m_block_type == RCORE) ||
           (block_name == "mcsol" && m_block_type == MCSOL) ||
           (block_name == "lcsol" && m_block_type == LCSOL) ||
           (block_name == "rcsol" && m_block_type == RCSOL) ||
           (block_name == "lpf" && m_block_type == LPF) ||
           (block_name == "rpf" && m_block_type == RPF) ||
           (block_name == "lsol" && m_block_type == LSOL) ||
           (block_name == "rsol" && m_block_type == RSOL) ) {
         found_block = true;
         break;
      }
      else {
         // Move the file pointer to the next block
         double dummy;
         for (int i=0; i<4*num_block_elements; ++i) {
            inFile >> dummy;
         }
      }
    }

    if ( !found_block ) {
       MayDay::Error("SingleNullBlockCoordSys::init(): Block not found in mapping file");
    }

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

  double *input_values = new double[4*num_block_elements];

#ifdef CH_MPI
  if (procID() == 0) {
#endif

    for (int i=0; i<4*num_block_elements; ++i) {
      inFile >> input_values[i];
    }

    inFile.close();

    cout << "Done reading geometry mapping file" << endl;

#ifdef CH_MPI
  }


  MPI_Bcast(input_values, 4*num_block_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif // end if MPI

  int num_radial_cells = n_radial_points - 1;

  double dr = ( upperMappedCoordinate(RADIAL_DIR)
              - lowerMappedCoordinate(RADIAL_DIR) ) / num_radial_cells;

  double *r_pts = new double[m_mapping_block_size[RADIAL_DIR]];

  r_pts[0] = lowerMappedCoordinate(RADIAL_DIR) - n_radial_extend*dr;
  for (int i=1; i<m_mapping_block_size[RADIAL_DIR]; ++i) {
    r_pts[i] = r_pts[i-1] + dr;
  }

  int num_poloidal_cells = n_poloidal_points - 1;

  double dtheta = ( upperMappedCoordinate(POLOIDAL_DIR)
                  - lowerMappedCoordinate(POLOIDAL_DIR) ) / num_poloidal_cells;

  double *theta_pts = new double[m_mapping_block_size[POLOIDAL_DIR]];

  int block_full_poloidal = num_poloidal_cells;
  int block_poloidal = block_full_poloidal;

  if (m_block_type == LCORE || m_block_type == LCSOL) {
    if (a_pp.contains("core_mapping_truncate")) {
      a_pp.get("core_mapping_truncate", block_poloidal);
      CH_assert(block_poloidal <= m_mapping_block_size[POLOIDAL_DIR]);
    }

    dtheta *= (double)block_full_poloidal / (double)block_poloidal;

    theta_pts[m_mapping_block_size[POLOIDAL_DIR]-1] =
      upperMappedCoordinate(POLOIDAL_DIR) + n_poloidal_extend*dtheta;
    for (int i=m_mapping_block_size[POLOIDAL_DIR]-1; i>0; --i) {
      theta_pts[i-1] = theta_pts[i] - dtheta;
    }
  }
  else if (m_block_type == RCORE || m_block_type == RCSOL) {
    if (a_pp.contains("core_mapping_truncate")) {
      a_pp.get("core_mapping_truncate", block_poloidal);
      CH_assert(block_poloidal <= m_mapping_block_size[POLOIDAL_DIR]);
    }

    dtheta *= (double)block_full_poloidal / (double)block_poloidal;

    theta_pts[0] = lowerMappedCoordinate(POLOIDAL_DIR) - n_poloidal_extend*dtheta;
    for (int i=1; i<m_mapping_block_size[POLOIDAL_DIR]; ++i) {
      theta_pts[i] = theta_pts[i-1] + dtheta;
    }
  }
  else if (m_block_type == LPF || m_block_type == LSOL) {
    if (a_pp.contains("pf_mapping_truncate")) {
      a_pp.get("pf_mapping_truncate", block_poloidal);
      CH_assert(block_poloidal <= m_mapping_block_size[POLOIDAL_DIR]);
    }

    dtheta *= (double)block_full_poloidal / (double)block_poloidal;

    theta_pts[0] = lowerMappedCoordinate(POLOIDAL_DIR) - n_poloidal_extend*dtheta;
    for (int i=1; i<m_mapping_block_size[POLOIDAL_DIR]; ++i) {
      theta_pts[i] = theta_pts[i-1] + dtheta;
    }
  }
  else if (m_block_type == RPF || m_block_type == RSOL) {
    if (a_pp.contains("pf_mapping_truncate")) {
      a_pp.get("pf_mapping_truncate", block_poloidal);
      CH_assert(block_poloidal <= m_mapping_block_size[POLOIDAL_DIR]);
    }

    dtheta *= (double)block_full_poloidal / (double)block_poloidal;

    theta_pts[m_mapping_block_size[POLOIDAL_DIR]-1] =
      upperMappedCoordinate(POLOIDAL_DIR) + n_poloidal_extend*dtheta;
    for (int i=m_mapping_block_size[POLOIDAL_DIR]-1; i>0; --i) {
      theta_pts[i-1] = theta_pts[i] - dtheta;
    }
  }
  else if (m_block_type == MCORE || m_block_type == MCSOL) {
    dtheta *= (double)block_full_poloidal / (double)block_poloidal;

    theta_pts[m_mapping_block_size[POLOIDAL_DIR]-1] =
      upperMappedCoordinate(POLOIDAL_DIR) + n_poloidal_extend*dtheta;
    for (int i=m_mapping_block_size[POLOIDAL_DIR]-1; i>0; --i) {
      theta_pts[i-1] = theta_pts[i] - dtheta;
    }
  }
  else {
    MayDay::Error("SingleNullBlockCoordSys::init(): Invalid block_type encountered");
  }

  m_poloidally_truncated = (block_poloidal < block_full_poloidal);

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

    RZ_data(iv,0) = input_values[4*k];
    RZ_data(iv,1) = input_values[4*k+1];

    interp_node_coords(iv,0) = r_pts[iv[0] - interp_box.smallEnd(0)];
    interp_node_coords(iv,1) = theta_pts[iv[1] - interp_box.smallEnd(1)];
  }

  if (interp_method == "bspline") {
     m_RZ_interp = new BSplineInterp(a_pp, interp_node_coords, RZ_data);
  }
  else if (interp_method == "hermite") {
     m_RZ_interp = new HermiteInterp(a_pp, interp_node_coords, RZ_data);
  }
  else if (interp_method == "spline") {
     m_RZ_interp = new SplineInterp(a_pp, interp_node_coords, RZ_data);
  }

  // If a field coefficients file was specified in the input, use it
  // to define the field.  Otherwise, we read it from the mapping file.
  m_spectral_field = a_pp.contains("field_coefficients_file");

  if ( !m_spectral_field ) {

     if (procID() == 0) {
        cout << "Reading magnetic field data from mapping file" << endl;
     }

     // Derive R times the magnetic field magnitude and the cylindrical components
     // of the field unit vector.  We need these quantities and their derivatives
     // at various points, so we create an interpolator for them.

     FArrayBox field_data(interp_box, NUM_FIELD_VARS);

     for (bit.begin();bit.ok();++bit) {
        IntVect iv = bit();
        int k = m_mapping_block_size[RADIAL_DIR]*iv[1] + iv[0];

        double RBr, RBz, RB;

        field_data(iv,RBR_VAR) = RBr = input_values[4*k+2];
        field_data(iv,RBZ_VAR) = RBz = input_values[4*k+3];
        field_data(iv,RB_VAR) = RB = sqrt( RBr*RBr + RBz*RBz + m_RB_toroidal*m_RB_toroidal );
        field_data(iv,UNIT_BR_VAR) = RBr / RB;
        field_data(iv,UNIT_BPHI_VAR) = m_RB_toroidal / RB;
        field_data(iv,UNIT_BZ_VAR) = RBz / RB;
     }

     if (interp_method == "bspline") {
        m_field_interp = new BSplineInterp(a_pp, interp_node_coords, field_data);
     }
     else if (interp_method == "hermite") {
        m_field_interp = new HermiteInterp(a_pp, interp_node_coords, field_data);
     }
     else if (interp_method == "spline") {
        m_field_interp = new SplineInterp(a_pp, interp_node_coords, field_data);
     }

  }
  else {

     ifstream fieldCoefficientsFile;

     int NR, NZ;

#ifdef CH_MPI
     if (procID() == 0) {
#endif
        string field_coefficients_file;
        if ( a_pp.contains("field_coefficients_file")) {
           a_pp.get("field_coefficients_file", field_coefficients_file);
        }
        else {
           MayDay::Error("SingleNullBlockCoordSys::init(): No field coefficients file specified");
        }

        fieldCoefficientsFile.open( field_coefficients_file.c_str() );

        if (!fieldCoefficientsFile) {
           cout << "Unable to open field coefficients file" << endl;
           exit(1);
        }

        cout << "Reading field coefficients from " << field_coefficients_file << endl;

        fieldCoefficientsFile >> m_Rmin;
        fieldCoefficientsFile >> m_Rmax;
        fieldCoefficientsFile >> m_Zmin;
        fieldCoefficientsFile >> m_Zmax;
        fieldCoefficientsFile >> NR;
        fieldCoefficientsFile >> NZ;
        fieldCoefficientsFile >> m_magAxis[0];
        fieldCoefficientsFile >> m_magAxis[1];
        fieldCoefficientsFile >> m_Xpoint[0];
        fieldCoefficientsFile >> m_Xpoint[1];

        
#ifdef CH_MPI
     }
#endif

     MPI_Bcast(&m_Rmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     MPI_Bcast(&m_Rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     MPI_Bcast(&m_Zmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     MPI_Bcast(&m_Zmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     MPI_Bcast(&NR, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Bcast(&NZ, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Bcast(&m_magAxis, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     MPI_Bcast(&m_Xpoint, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     
#if CFG_DIM==3
     Box kbox = Box(IntVect::Zero,IntVect(NR-1,NZ-1,1));
#else
     Box kbox = Box(IntVect::Zero,IntVect(NR-1,NZ-1));
#endif

     m_psi_coefs.define(kbox,1);

#ifdef CH_MPI
     if (procID() == 0) {
#endif
        for (BoxIterator bit(kbox); bit.ok(); ++bit) {
           fieldCoefficientsFile >> m_psi_coefs(bit());
        }
        
#ifdef CH_MPI
     }
#endif // end if MPI

     int num_values = kbox.numPts();
     MPI_Bcast(m_psi_coefs.dataPtr(), num_values, MPI_DOUBLE, 0, MPI_COMM_WORLD);

     fieldCoefficientsFile.close();

     if ( a_pp.contains("psi_expansion_order")) {
        a_pp.get("psi_expansion_order", m_psi_expansion_order);
        if (m_psi_expansion_order > NR || m_psi_expansion_order > NZ) {
           MayDay::Error("SingleNullBlockCoordSys::init(): Insufficient number of DCT coefficients for requested psi_expansion_order");
        }
     }
     else {
        m_psi_expansion_order = NR;
        if (NZ < NR) m_psi_expansion_order = NZ;
     }
  }

  delete [] r_pts;
  delete [] theta_pts;
  delete [] input_values;

  getCellCenterRealCoords();

}



RealVect
SingleNullBlockCoordSys::realCoord( const RealVect& a_xi ) const
{
  RealVect x;
  D_TERM(x[0] = m_RZ_interp->interpolate(R_VAR, FUNCTION, a_xi);,
         x[1] = m_RZ_interp->interpolate(Z_VAR, FUNCTION, a_xi);,
         x[2] = 0.;)

  return x;
}



RealVect
SingleNullBlockCoordSys::mappedCoord( const RealVect& a_X ) const
{
#ifdef PLOT_INVERSE_POINTS
   fprintf(m_ipt_file[m_block_type], "%20.12e %20.12e\n", a_X[0], a_X[1]);
#endif

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
SingleNullBlockCoordSys::mappedCoordNewton( const RealVect& a_X,
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

      // Update the current solution and residual norm
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
      cout << "On block " << m_block_type << ": Newton solve did not converge at " << a_X << ", a_xi_initial = " << a_xi_initial << ", iv_init = " << a_iv_initial << ", Fnorm = " << Fnorm << endl;
#if 0
      for (int i=0; i<num_iters; ++i) {
         cout << i << " xi = " << xi_saved[i] << ", Fnorm = " << Fnorm_saved[i] << endl;
      }
#endif
   }
#endif

   return xi;
}


Real
SingleNullBlockCoordSys::dXdXi( const RealVect& a_Xi,
                                int             a_dirX,
                                int             a_dirXi ) const
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

    double metric_factor;

    RealVect metrics;
    for (int dir(0); dir<SpaceDim; ++dir) {
      metrics[dir] = (x_hi[dir] - x_lo[dir])/m_dx[a_dirXi] * m_N_mr[a_dirXi]/2.0;
    }

    return metrics[a_dirX];
  }
}

void
SingleNullBlockCoordSys::dXdXi( FArrayBox&       a_dXdXi,
                                const FArrayBox& a_Xi,
                                int              a_destComp,
                                int              a_dirX,
                                int              a_dirXi,
                                const Box&       a_box ) const
{
  if (!m_subgrid_geometry) m_RZ_interp->interpolate(a_dirX, a_destComp, a_dirXi+1, a_Xi, a_box, a_dXdXi);

  else {
    BoxIterator bit(a_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect this_Xi;
      for (int dir=0; dir<SpaceDim; ++dir) {
         this_Xi[dir] = a_Xi(iv,dir);
      }
      a_dXdXi(iv,a_destComp) = dXdXi(this_Xi,a_dirX,a_dirXi);
    }
  }
}



bool
SingleNullBlockCoordSys::isValid(const RealVect& xi) const
{
  return xi[RADIAL_DIR]   >= lowerMappedCoordinate(RADIAL_DIR)   &&
         xi[RADIAL_DIR]   <= upperMappedCoordinate(RADIAL_DIR)   &&
         xi[POLOIDAL_DIR] >= lowerMappedCoordinate(POLOIDAL_DIR) &&
         xi[POLOIDAL_DIR] <= upperMappedCoordinate(POLOIDAL_DIR);
}


void
SingleNullBlockCoordSys::getMagneticFlux( const FArrayBox& a_physical_coordinates,
                                         FArrayBox&       a_magnetic_flux ) const
{
   if (m_spectral_field) {
      getMagneticFluxFromDCT(a_physical_coordinates, a_magnetic_flux);
   }
   else {
      MayDay::Error("SingleNullBlockCoordSys::getMagneticFlux() only implemented for DCT field option");
   }
}


double
SingleNullBlockCoordSys::getMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   double psi;
   if (m_spectral_field) {
      psi = getMagneticFluxFromDCT(a_physical_coordinate);
   }
   else {
      MayDay::Error("SingleNullBlockCoordSys::getMagneticFlux() only implemented for DCT field option");
      psi = 0.;
   }
   return psi;
}


void SingleNullBlockCoordSys::getNormMagneticFlux( const FArrayBox& a_physical_coordinates,
                                                  FArrayBox&       a_magnetic_flux ) const

{
   double physFluxOnAxis = getMagneticFlux(m_magAxis);
   double physFluxOnSep = getMagneticFlux(m_Xpoint);
   
   getMagneticFlux(a_physical_coordinates, a_magnetic_flux);
   a_magnetic_flux.plus(-physFluxOnAxis);
   a_magnetic_flux.divide(physFluxOnSep - physFluxOnAxis);
   
}

double SingleNullBlockCoordSys::getNormMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   
   double physFluxOnAxis = getMagneticFlux(m_magAxis);
   double physFluxOnSep = getMagneticFlux(m_Xpoint);
   double physFlux = getMagneticFlux(a_physical_coordinate);
   
   return (physFlux - physFluxOnAxis)/(physFluxOnSep - physFluxOnAxis);
}


double SingleNullBlockCoordSys::getOuterMidplaneCoord(const double psiNorm) const
{
   bool residual_tolerance_satisfied = false;
   int max_iterations = 10000;
   double tol = 1.0e-5;
   
   double tmp_lo_R = m_magAxis[0];
   
   double tmp_hi_R, half_R, f_lo_R, f_hi_R, f_half_R;
   
   if (m_spectral_field) tmp_hi_R = m_Rmax; // MAYBE WILL NEED TO FIX LATER IN CASE THERE ARE MULTIVALUE SOLUTIONS NEAR M_RMAX
   else MayDay::Error("SingleNullBlockCoordSys::getOuterMidplaneCoord only implemented for DCT field option");
   
   int iter = 0;
   while ( iter++ < max_iterations && !residual_tolerance_satisfied ) {
      
      half_R = 0.5 * (tmp_lo_R + tmp_hi_R);

#if CH_SPACEDIM == 2
      RealVect lo_R_vect(tmp_lo_R,m_magAxis[1]);
      RealVect hi_R_vect(tmp_hi_R,m_magAxis[1]);
      RealVect half_R_vect(half_R,m_magAxis[1]);
#endif

#if CH_SPACEDIM == 3

      //Assume that the third direction is a toroidal angle
      //and that the geoemtry is axisymmetric
      double phi = 0.0;
      
      RealVect lo_R_vect(tmp_lo_R,m_magAxis[1],phi);
      RealVect hi_R_vect(tmp_hi_R,m_magAxis[1],phi);
      RealVect half_R_vect(half_R,m_magAxis[1],phi);
#endif
      
      f_lo_R = getNormMagneticFlux(lo_R_vect) - psiNorm;
      f_hi_R = getNormMagneticFlux(hi_R_vect) - psiNorm;
      f_half_R = getNormMagneticFlux(half_R_vect) - psiNorm;
      
      if (f_lo_R * f_half_R < 0) {
         tmp_hi_R = half_R;
      }
      
      else {
         tmp_lo_R = half_R;
      }

#if CH_SPACEDIM == 2
      RealVect R_hi_test(half_R + tol,m_magAxis[1]);
      RealVect R_lo_test(half_R - tol,m_magAxis[1]);
#endif

#if CH_SPACEDIM == 3
      RealVect R_hi_test(half_R + tol,m_magAxis[1],phi);
      RealVect R_lo_test(half_R - tol,m_magAxis[1],phi);
#endif

      if ( (getNormMagneticFlux(R_hi_test) - psiNorm) * (getNormMagneticFlux(R_lo_test) - psiNorm) < 0 ) {
         residual_tolerance_satisfied = true;
      }
      
   }
   
   return half_R;
}

double
SingleNullBlockCoordSys::getOuterRsep() const
{
   double psiNorm_atSep = getNormMagneticFlux(m_Xpoint);
   
   return getOuterMidplaneCoord(psiNorm_atSep);
}


void
SingleNullBlockCoordSys::computeFieldData( const int  a_dir,
                                           FArrayBox& a_BField,
                                           FArrayBox& a_BFieldMag,
                                           FArrayBox& a_BFieldDir,
                                           FArrayBox& a_gradBFieldMag,
                                           FArrayBox& a_curlBFieldDir,
                                           FArrayBox& a_BFieldDirdotcurlBFieldDir,
                                           const bool a_derived_data_only ) const
{
   if ( m_spectral_field ) {
      computeFieldDataSpectral(a_dir, a_BField, a_BFieldMag, a_BFieldDir, a_gradBFieldMag,
                               a_curlBFieldDir, a_BFieldDirdotcurlBFieldDir, a_derived_data_only);
   }
   else {
      computeFieldDataFromMappingFile(a_dir, a_BField, a_BFieldMag, a_BFieldDir, a_gradBFieldMag,
                                      a_curlBFieldDir, a_BFieldDirdotcurlBFieldDir, a_derived_data_only);

   }
}


void
SingleNullBlockCoordSys::computeFieldDataSpectral( const int  a_dir,
                                                   FArrayBox& a_BField,
                                                   FArrayBox& a_BFieldMag,
                                                   FArrayBox& a_BFieldDir,
                                                   FArrayBox& a_gradBFieldMag,
                                                   FArrayBox& a_curlBFieldDir,
                                                   FArrayBox& a_BFieldDirdotcurlBFieldDir,
                                                   const bool a_derived_data_only ) const
{
   CH_TIME("SingleNullBlockCoordSys::computeFieldDataSpectral");

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

   int NR = m_psi_coefs.box().size(0);
   int NZ = m_psi_coefs.box().size(1);

   int llen = NR>NZ? NR: NZ;

   double * temp = new double[7*llen];
   double * lambda = temp;
   double * facR = lambda + llen;
   double * facZ = facR + llen;
   double * sinfacR = facZ + llen;
   double * cosfacR = sinfacR + llen;
   double * sinfacZ = cosfacR + llen;
   double * cosfacZ = sinfacZ + llen;

   lambda[0] = 1. / sqrt(2.);
   for (int l=1; l<llen; ++l) {
      lambda[l] = 1.;
   }

   for (int i=0; i<NR; ++i) {
      facR[i] = i * Pi / NR;
   }

   for (int j=0; j<NZ; ++j) {
      facZ[j] = j * Pi / NZ;
   }

   double Rscale = (NR-1)/(m_Rmax - m_Rmin);
   double Zscale = (NZ-1)/(m_Zmax - m_Zmin);

   for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();
      
      double sR = (RZ(iv,0) - m_Rmin) * Rscale;
      double sZ = (RZ(iv,1) - m_Zmin) * Zscale;

      for (int i=0; i<NR; ++i) {
         double t = facR[i] * (sR + 0.5);
         sinfacR[i] = sin(t);
         cosfacR[i] = cos(t);
      }

      for (int j=0; j<NZ; ++j) {
         double t = facZ[j] * (sZ + 0.5);
         sinfacZ[j] = sin(t);
         cosfacZ[j] = cos(t);
      }

      double value;

      // Compute dpsi_dR

      int derivR = 1;
      int derivZ = 0;

      FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                       CHF_CONST_INT(m_psi_expansion_order),
                       CHF_CONST_INT(derivR),
                       CHF_CONST_INT(derivZ),
                       CHF_R1D(facR,NR),
                       CHF_R1D(facZ,NZ),
                       CHF_R1D(sinfacR,NR),
                       CHF_R1D(cosfacR,NR),
                       CHF_R1D(sinfacZ,NZ),
                       CHF_R1D(cosfacZ,NZ),
                       CHF_R1D(lambda,llen),
                       CHF_REAL(value) );

      value *= Rscale;

      RB(iv,1) = value;  // R*B_Z

      // Compute dpsi_dZ

      derivR = 0;
      derivZ = 1;

      FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                       CHF_CONST_INT(m_psi_expansion_order),
                       CHF_CONST_INT(derivR),
                       CHF_CONST_INT(derivZ),
                       CHF_R1D(facR,NR),
                       CHF_R1D(facZ,NZ),
                       CHF_R1D(sinfacR,NR),
                       CHF_R1D(cosfacR,NR),
                       CHF_R1D(sinfacZ,NZ),
                       CHF_R1D(cosfacZ,NZ),
                       CHF_R1D(lambda,llen),
                       CHF_REAL(value) );

      value *= Zscale;

      RB(iv,0) = -value;  // R*B_R

      // Compute d2psi_dR2

      derivR = 2;
      derivZ = 0;

      FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                       CHF_CONST_INT(m_psi_expansion_order),
                       CHF_CONST_INT(derivR),
                       CHF_CONST_INT(derivZ),
                       CHF_R1D(facR,NR),
                       CHF_R1D(facZ,NZ),
                       CHF_R1D(sinfacR,NR),
                       CHF_R1D(cosfacR,NR),
                       CHF_R1D(sinfacZ,NZ),
                       CHF_R1D(cosfacZ,NZ),
                       CHF_R1D(lambda,llen),
                       CHF_REAL(value) );

      value *= (Rscale * Rscale);

      dRBdR(iv,1) = value;  // d(R*BZ)/dR

      // Compute d2psi_dZ2

      derivR = 0;
      derivZ = 2;

      FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                       CHF_CONST_INT(m_psi_expansion_order),
                       CHF_CONST_INT(derivR),
                       CHF_CONST_INT(derivZ),
                       CHF_R1D(facR,NR),
                       CHF_R1D(facZ,NZ),
                       CHF_R1D(sinfacR,NR),
                       CHF_R1D(cosfacR,NR),
                       CHF_R1D(sinfacZ,NZ),
                       CHF_R1D(cosfacZ,NZ),
                       CHF_R1D(lambda,llen),
                       CHF_REAL(value) );

      value *= (Zscale * Zscale);

      dRBdZ(iv,0) = -value;  // d(R*B_R)/dZ

      // Compute d2psi_dRdZ

      derivR = 1;
      derivZ = 1;

      FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                       CHF_CONST_INT(m_psi_expansion_order),
                       CHF_CONST_INT(derivR),
                       CHF_CONST_INT(derivZ),
                       CHF_R1D(facR,NR),
                       CHF_R1D(facZ,NZ),
                       CHF_R1D(sinfacR,NR),
                       CHF_R1D(cosfacR,NR),
                       CHF_R1D(sinfacZ,NZ),
                       CHF_R1D(cosfacZ,NZ),
                       CHF_R1D(lambda,llen),
                       CHF_REAL(value) );

      value *= (Rscale * Zscale);

      dRBdR(iv,0) = -value;  // d(R*B_R)/dR
      dRBdZ(iv,1) =  value;  // d(R*B_Z)/dZ
   }

   delete[] temp;

   FORT_GET_FIELD_DATA( CHF_BOX(box),
                        CHF_CONST_FRA(RZ),
                        CHF_CONST_FRA(RB),
                        CHF_CONST_FRA(dRBdR),
                        CHF_CONST_FRA(dRBdZ),
                        CHF_CONST_REAL(m_RB_toroidal),
                        CHF_FRA(a_BField),
                        CHF_FRA1(a_BFieldMag,0),
                        CHF_FRA(a_BFieldDir),
                        CHF_FRA(a_gradBFieldMag),
                        CHF_FRA(a_curlBFieldDir),
                        CHF_FRA1(a_BFieldDirdotcurlBFieldDir,0) );

}


Vector<Real>
SingleNullBlockCoordSys::computeBField(const RealVect& a_X) const
{
   CH_TIME("SingleNullBlockCoordSys::computeBField");

   if ( !m_spectral_field ) {
      MayDay::Error("SingleNullBlockCoordSys::computeBField(): Only implemented for DCT field option");
   }
   
   int NR = m_psi_coefs.box().size(0);
   int NZ = m_psi_coefs.box().size(1);
   
   int llen = NR>NZ? NR: NZ;
   
   double * temp = new double[7*llen];
   double * lambda = temp;
   double * facR = lambda + llen;
   double * facZ = facR + llen;
   double * sinfacR = facZ + llen;
   double * cosfacR = sinfacR + llen;
   double * sinfacZ = cosfacR + llen;
   double * cosfacZ = sinfacZ + llen;
   
   lambda[0] = 1. / sqrt(2.);
   for (int l=1; l<llen; ++l) {
      lambda[l] = 1.;
   }
   
   for (int i=0; i<NR; ++i) {
      facR[i] = i * Pi / NR;
   }
   
   for (int j=0; j<NZ; ++j) {
      facZ[j] = j * Pi / NZ;
   }
   
   double Rscale = (NR-1)/(m_Rmax - m_Rmin);
   double Zscale = (NZ-1)/(m_Zmax - m_Zmin);
   
   
   double sR = (a_X[0] - m_Rmin) * Rscale;
   double sZ = (a_X[1] - m_Zmin) * Zscale;
      
   for (int i=0; i<NR; ++i) {
      double t = facR[i] * (sR + 0.5);
      sinfacR[i] = sin(t);
      cosfacR[i] = cos(t);
   }
      
   for (int j=0; j<NZ; ++j) {
      double t = facZ[j] * (sZ + 0.5);
      sinfacZ[j] = sin(t);
      cosfacZ[j] = cos(t);
   }
      
   double value;
   
   RealVect RB;
   
   // Compute dpsi_dR
      
   int derivR = 1;
   int derivZ = 0;
   
   FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                    CHF_CONST_INT(m_psi_expansion_order),
                    CHF_CONST_INT(derivR),
                    CHF_CONST_INT(derivZ),
                    CHF_R1D(facR,NR),
                    CHF_R1D(facZ,NZ),
                    CHF_R1D(sinfacR,NR),
                    CHF_R1D(cosfacR,NR),
                    CHF_R1D(sinfacZ,NZ),
                    CHF_R1D(cosfacZ,NZ),
                    CHF_R1D(lambda,llen),
                    CHF_REAL(value) );
      
   value *= Rscale;
      
   RB[1] = value;  // R*B_Z
      
   // Compute dpsi_dZ
      
   derivR = 0;
   derivZ = 1;
      
   FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                    CHF_CONST_INT(m_psi_expansion_order),
                    CHF_CONST_INT(derivR),
                    CHF_CONST_INT(derivZ),
                    CHF_R1D(facR,NR),
                    CHF_R1D(facZ,NZ),
                    CHF_R1D(sinfacR,NR),
                    CHF_R1D(cosfacR,NR),
                    CHF_R1D(sinfacZ,NZ),
                    CHF_R1D(cosfacZ,NZ),
                    CHF_R1D(lambda,llen),
                    CHF_REAL(value) );
      
   value *= Zscale;
      
   RB[0] = -value;  // R*B_R
      
   Vector<Real> result(3,0);
   
   result[0] = RB[0]/a_X[0];
   result[1] = m_RB_toroidal/ a_X[0];
   result[2] = RB[1]/a_X[0];
   
   delete[] temp;
   
   return result;
}


void
SingleNullBlockCoordSys::getNodalFieldData(FArrayBox& a_points,
                                           FArrayBox& a_A,
                                           FArrayBox& a_b,
                                           FArrayBox& a_Bmag) const
{
   CH_TIME("SingleNullBlockCoordSys::getNodalFieldData");
   
   if ( m_spectral_field ) {

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
      FArrayBox dRBdR(box,2);
      FArrayBox dRBdZ(box,2);

      int NR = m_psi_coefs.box().size(0);
      int NZ = m_psi_coefs.box().size(1);

      int llen = NR>NZ? NR: NZ;

      double * temp = new double[7*llen];
      double * lambda = temp;
      double * facR = lambda + llen;
      double * facZ = facR + llen;
      double * sinfacR = facZ + llen;
      double * cosfacR = sinfacR + llen;
      double * sinfacZ = cosfacR + llen;
      double * cosfacZ = sinfacZ + llen;

      lambda[0] = 1. / sqrt(2.);
      for (int l=1; l<llen; ++l) {
         lambda[l] = 1.;
      }

      for (int i=0; i<NR; ++i) {
         facR[i] = i * Pi / NR;
      }

      for (int j=0; j<NZ; ++j) {
         facZ[j] = j * Pi / NZ;
      }

      double Rscale = (NR-1)/(m_Rmax - m_Rmin);
      double Zscale = (NZ-1)/(m_Zmax - m_Zmin);

      for (BoxIterator bit(box); bit.ok(); ++bit) {
         IntVect iv = bit();
      
         double sR = (RZ(iv,0) - m_Rmin) * Rscale;
         double sZ = (RZ(iv,1) - m_Zmin) * Zscale;

         for (int i=0; i<NR; ++i) {
            double t = facR[i] * (sR + 0.5);
            sinfacR[i] = sin(t);
            cosfacR[i] = cos(t);
         }

         for (int j=0; j<NZ; ++j) {
            double t = facZ[j] * (sZ + 0.5);
            sinfacZ[j] = sin(t);
            cosfacZ[j] = cos(t);
         }

         double value;

         // Compute psi

         int derivR = 0;
         int derivZ = 0;

         FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                          CHF_CONST_INT(m_psi_expansion_order),
                          CHF_CONST_INT(derivR),
                          CHF_CONST_INT(derivZ),
                          CHF_R1D(facR,NR),
                          CHF_R1D(facZ,NZ),
                          CHF_R1D(sinfacR,NR),
                          CHF_R1D(cosfacR,NR),
                          CHF_R1D(sinfacZ,NZ),
                          CHF_R1D(cosfacZ,NZ),
                          CHF_R1D(lambda,llen),
                          CHF_REAL(value) );

         psi(iv,0) = value;

         // Compute dpsi_dR

         derivR = 1;
         derivZ = 0;

         FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                          CHF_CONST_INT(m_psi_expansion_order),
                          CHF_CONST_INT(derivR),
                          CHF_CONST_INT(derivZ),
                          CHF_R1D(facR,NR),
                          CHF_R1D(facZ,NZ),
                          CHF_R1D(sinfacR,NR),
                          CHF_R1D(cosfacR,NR),
                          CHF_R1D(sinfacZ,NZ),
                          CHF_R1D(cosfacZ,NZ),
                          CHF_R1D(lambda,llen),
                          CHF_REAL(value) );

         value *= Rscale;

         RB(iv,1) = value;  // R*B_Z

         // Compute dpsi_dZ

         derivR = 0;
         derivZ = 1;

         FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                          CHF_CONST_INT(m_psi_expansion_order),
                          CHF_CONST_INT(derivR),
                          CHF_CONST_INT(derivZ),
                          CHF_R1D(facR,NR),
                          CHF_R1D(facZ,NZ),
                          CHF_R1D(sinfacR,NR),
                          CHF_R1D(cosfacR,NR),
                          CHF_R1D(sinfacZ,NZ),
                          CHF_R1D(cosfacZ,NZ),
                          CHF_R1D(lambda,llen),
                          CHF_REAL(value) );

         value *= Zscale;

         RB(iv,0) = -value;  // R*B_R
      }

      delete[] temp;

      FORT_GET_NODAL_FIELD_DATA( CHF_BOX(box),
                                 CHF_CONST_FRA(RZ),
                                 CHF_CONST_FRA1(psi,0),
                                 CHF_CONST_FRA(RB),
                                 CHF_CONST_REAL(m_RB_toroidal),
                                 CHF_FRA(a_A),
                                 CHF_FRA(a_b),
                                 CHF_FRA(a_Bmag));
   }
   else {
      MayDay::Error("SingleNullBlockCoordSys::getNodalFieldData(): Only implemented for DCT field option");
   }
}


void
SingleNullBlockCoordSys::computeFieldDataFromMappingFile( const int  a_dir,
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

   FArrayBox Xi(box,SpaceDim);
   if (a_dir>=0 && a_dir<SpaceDim) {
     getFaceCenteredMappedCoords(a_dir, Xi);
   }
   else {
     getCellCenteredMappedCoords(Xi);
   }

   FArrayBox RZ_data(box,5);
   m_RZ_interp->interpolate(R_VAR, 0, FUNCTION,            Xi, box, RZ_data);
   m_RZ_interp->interpolate(R_VAR, 1, RADIAL_DERIVATIVE,   Xi, box, RZ_data);
   m_RZ_interp->interpolate(R_VAR, 2, POLOIDAL_DERIVATIVE, Xi, box, RZ_data);
   m_RZ_interp->interpolate(Z_VAR, 3, RADIAL_DERIVATIVE,   Xi, box, RZ_data);
   m_RZ_interp->interpolate(Z_VAR, 4, POLOIDAL_DERIVATIVE, Xi, box, RZ_data);

   if( !a_derived_data_only ) {

      FArrayBox RBpoloidal(box, 2);
      m_field_interp->interpolate(RBR_VAR, 0, FUNCTION, Xi, box, RBpoloidal);
      m_field_interp->interpolate(RBZ_VAR, 1, FUNCTION, Xi, box, RBpoloidal);

      // Construct the field from the mapping data, RBpol and RBtor
      FORT_GET_FIELD_FROM_RBPOL_RBTOR(CHF_BOX(box),
                                      CHF_CONST_REAL(m_RB_toroidal),
                                      CHF_CONST_FRA(RBpoloidal),
                                      CHF_CONST_FRA1(RZ_data,0),
                                      CHF_FRA(a_BField));
   }

   // Compute the field magnitude and unit vector
   FORT_GET_FIELD_MAGNITUDE_AND_UNITVECTOR(CHF_BOX(box),
                                           CHF_CONST_FRA(a_BField),
                                           CHF_FRA1(a_BFieldMag,0),
                                           CHF_FRA(a_BFieldDir));

   // Compute RB = the major radius times the field magnitude
   FArrayBox RB(box, 1);
   RB.copy(RZ_data,0,0,1);
   RB *= a_BFieldMag;
   
   // Differentiate RB and the field unit vector components
   FArrayBox derivative_data(box, 8);
   if( a_derived_data_only ) {

      FArrayBox diffed_data(box, 4);
      BoxIterator bit(box);
      for (bit.begin();bit.ok();++bit) {
         IntVect iv = bit();
         diffed_data(iv,0) = RB(iv);
         for (int m=1; m<4; ++m) {
            diffed_data(iv,m) = a_BFieldDir(iv,m-1);
         }
      }

      ParmParse pp( pp_name.c_str() );
      Interp* interp = new HermiteInterp(pp, Xi, diffed_data);

      for( int i=0; i<4; ++i ) {
         interp->interpolate(i, 2*i  , RADIAL_DERIVATIVE,   Xi, box, derivative_data);
         interp->interpolate(i, 2*i+1, POLOIDAL_DERIVATIVE, Xi, box, derivative_data);
      }

      delete interp;
   }
   else {
      m_field_interp->interpolate(RB_VAR,        0, RADIAL_DERIVATIVE,   Xi, box, derivative_data);
      m_field_interp->interpolate(RB_VAR,        1, POLOIDAL_DERIVATIVE, Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BR_VAR,   2, RADIAL_DERIVATIVE,   Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BR_VAR,   3, POLOIDAL_DERIVATIVE, Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BPHI_VAR, 4, RADIAL_DERIVATIVE,   Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BPHI_VAR, 5, POLOIDAL_DERIVATIVE, Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BZ_VAR,   6, RADIAL_DERIVATIVE,   Xi, box, derivative_data);
      m_field_interp->interpolate(UNIT_BZ_VAR,   7, POLOIDAL_DERIVATIVE, Xi, box, derivative_data);
   }

   int axisymmetric = m_axisymmetric? 1: 0;

   // Compute the differentiated field quantities: grad(Bmag), curl(b) and b dot curl(b)
   FORT_GET_FIELD_DERIVATIVE_DATA(CHF_BOX(box),
                                  CHF_CONST_INT(axisymmetric),
                                  CHF_CONST_FRA1(RB,0),                      // RB
                                  CHF_CONST_FRA1(derivative_data,0),         // dRBdr
                                  CHF_CONST_FRA1(derivative_data,1),         // dRBdt
                                  CHF_CONST_FRA(a_BFieldDir),                // bunit
                                  CHF_CONST_FRA1(derivative_data,2),         // dbunitRdr
                                  CHF_CONST_FRA1(derivative_data,3),         // dbunitRdt
                                  CHF_CONST_FRA1(derivative_data,4),         // dbunitphidr
                                  CHF_CONST_FRA1(derivative_data,5),         // dbunitphdt
                                  CHF_CONST_FRA1(derivative_data,6),         // dbunitZdr
                                  CHF_CONST_FRA1(derivative_data,7),         // dbunitZdt
                                  CHF_CONST_FRA1(RZ_data,0),                 // R
                                  CHF_CONST_FRA1(RZ_data,1),                 // Rr
                                  CHF_CONST_FRA1(RZ_data,2),                 // Rt
                                  CHF_CONST_FRA1(RZ_data,3),                 // Zr
                                  CHF_CONST_FRA1(RZ_data,4),                 // Zt
                                  CHF_FRA(a_gradBFieldMag),                  // gradB
                                  CHF_FRA(a_curlBFieldDir),                  // curlbunit
                                  CHF_FRA1(a_BFieldDirdotcurlBFieldDir,0));  // bdotcurlbunit
}



void
SingleNullBlockCoordSys::getMagneticFluxFromDCT( const FArrayBox& a_physical_coordinates,
                                                 FArrayBox&       a_magnetic_flux ) const
{
   CH_TIME("SingleNullBlockCoordSys::getMagneticFluxFromDCT");

   /*
     Given the FArrayBox a_physical_coordinates of physical coordinates, fill the
     components of a_magnetic_flux with the magnetic flux by evaluating its 
     DCT expansion.
   */

   const Box& box(a_physical_coordinates.box());
   CH_assert(a_magnetic_flux.box().contains(box));

   int NR = m_psi_coefs.box().size(0);
   int NZ = m_psi_coefs.box().size(1);

   int llen = NR>NZ? NR: NZ;

   double * temp = new double[7*llen];
   double * lambda = temp;
   double * facR = lambda + llen;
   double * facZ = facR + llen;
   double * sinfacR = facZ + llen;
   double * cosfacR = sinfacR + llen;
   double * sinfacZ = cosfacR + llen;
   double * cosfacZ = sinfacZ + llen;

   lambda[0] = 1. / sqrt(2.);
   for (int l=1; l<llen; ++l) {
      lambda[l] = 1.;
   }

   for (int i=0; i<NR; ++i) {
      facR[i] = i * Pi / NR;
   }

   for (int j=0; j<NZ; ++j) {
      facZ[j] = j * Pi / NZ;
   }

   double Rscale = (NR-1)/(m_Rmax - m_Rmin);
   double Zscale = (NZ-1)/(m_Zmax - m_Zmin);

   for (BoxIterator bit(box); bit.ok(); ++bit) {
      IntVect iv = bit();
      
      double sR = (a_physical_coordinates(iv,0) - m_Rmin) * Rscale;
      double sZ = (a_physical_coordinates(iv,1) - m_Zmin) * Zscale;

      for (int i=0; i<NR; ++i) {
         double t = facR[i] * (sR + 0.5);
         sinfacR[i] = sin(t);
         cosfacR[i] = cos(t);
      }

      for (int j=0; j<NZ; ++j) {
         double t = facZ[j] * (sZ + 0.5);
         sinfacZ[j] = sin(t);
         cosfacZ[j] = cos(t);
      }

      // Compute psi

      int derivR = 0;
      int derivZ = 0;

      FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                       CHF_CONST_INT(m_psi_expansion_order),
                       CHF_CONST_INT(derivR),
                       CHF_CONST_INT(derivZ),
                       CHF_R1D(facR,NR),
                       CHF_R1D(facZ,NZ),
                       CHF_R1D(sinfacR,NR),
                       CHF_R1D(cosfacR,NR),
                       CHF_R1D(sinfacZ,NZ),
                       CHF_R1D(cosfacZ,NZ),
                       CHF_R1D(lambda,llen),
                       CHF_REAL(a_magnetic_flux(iv,0)) );
   }

   delete[] temp;
}



double
SingleNullBlockCoordSys::getMagneticFluxFromDCT( const RealVect& a_physical_coordinate ) const
{
   CH_TIME("SingleNullBlockCoordSys::getMagneticFluxFromDCT");

   /*
     Given the input physical coordinate a_physical_coordinate, return the magnetic flux
     at that point by evaluating its DCT expansion.
   */

   int NR = m_psi_coefs.box().size(0);
   int NZ = m_psi_coefs.box().size(1);

   int llen = NR>NZ? NR: NZ;

   double * temp = new double[7*llen];
   double * lambda = temp;
   double * facR = lambda + llen;
   double * facZ = facR + llen;
   double * sinfacR = facZ + llen;
   double * cosfacR = sinfacR + llen;
   double * sinfacZ = cosfacR + llen;
   double * cosfacZ = sinfacZ + llen;

   lambda[0] = 1. / sqrt(2.);
   for (int l=1; l<llen; ++l) {
      lambda[l] = 1.;
   }

   for (int i=0; i<NR; ++i) {
      facR[i] = i * Pi / NR;
   }

   for (int j=0; j<NZ; ++j) {
      facZ[j] = j * Pi / NZ;
   }

   double Rscale = (NR-1)/(m_Rmax - m_Rmin);
   double Zscale = (NZ-1)/(m_Zmax - m_Zmin);

   double psi;

   double sR = (a_physical_coordinate[0] - m_Rmin) * Rscale;
   double sZ = (a_physical_coordinate[1] - m_Zmin) * Zscale;

   for (int i=0; i<NR; ++i) {
      double t = facR[i] * (sR + 0.5);
      sinfacR[i] = sin(t);
      cosfacR[i] = cos(t);
   }

   for (int j=0; j<NZ; ++j) {
      double t = facZ[j] * (sZ + 0.5);
      sinfacZ[j] = sin(t);
      cosfacZ[j] = cos(t);
   }

   // Compute psi

   int derivR = 0;
   int derivZ = 0;

   FORT_DCT_INTERP( CHF_CONST_FRA1(m_psi_coefs,0),
                    CHF_CONST_INT(m_psi_expansion_order),
                    CHF_CONST_INT(derivR),
                    CHF_CONST_INT(derivZ),
                    CHF_R1D(facR,NR),
                    CHF_R1D(facZ,NZ),
                    CHF_R1D(sinfacR,NR),
                    CHF_R1D(cosfacR,NR),
                    CHF_R1D(sinfacZ,NZ),
                    CHF_R1D(cosfacZ,NZ),
                    CHF_R1D(lambda,llen),
                    CHF_REAL(psi) );

   delete[] temp;

   return psi;
}

#include "NamespaceFooter.H"
