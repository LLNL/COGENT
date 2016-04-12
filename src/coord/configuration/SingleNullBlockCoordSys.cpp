#include <float.h>
#include <fstream>
#include "SingleNullBlockCoordSys.H"
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

const std::string SingleNullBlockCoordSys::pp_name = "singlenull";

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


#ifndef MODEL_GEOM

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

#ifdef CH_MPI
     }
#endif

     MPI_Bcast(&m_Rmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     MPI_Bcast(&m_Rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     MPI_Bcast(&m_Zmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     MPI_Bcast(&m_Zmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     MPI_Bcast(&NR, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Bcast(&NZ, 1, MPI_INT, 0, MPI_COMM_WORLD);

     Box kbox = Box(IntVect::Zero,IntVect(NR-1,NZ-1));

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
  return m_RZ_interp->interpolate(a_dirX, a_dirXi+1, a_Xi);
}



void
SingleNullBlockCoordSys::dXdXi( FArrayBox&       a_dXdXi,
                                const FArrayBox& a_Xi,
                                int              a_destComp,
                                int              a_dirX,
                                int              a_dirXi,
                                const Box&       a_box ) const
{
  m_RZ_interp->interpolate(a_dirX, a_destComp, a_dirXi+1, a_Xi, a_box, a_dXdXi);
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

   FORT_GET_FIELD_DATA_SPECTRAL( CHF_BOX(box),
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



void
SingleNullBlockCoordSys::getMagneticFluxFromDCT( const FArrayBox& a_physical_coordinates,
                                                 FArrayBox&       a_magnetic_flux ) const
{
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


#endif


//======================================= BEGIN: THE CASE OF ANALYTIC MODEL GEOMETRY ==================================


#ifdef MODEL_GEOM

void SingleNullBlockCoordSys::getN(FluxBox& a_N, const Box& a_box) const
{
  for (int dir=0; dir<SpaceDim; ++dir) {
    FArrayBox& this_N = a_N[dir];
    const Box& box = this_N.box();
        
    FArrayBox xi(box,SpaceDim);
    getFaceCenteredMappedCoords(dir, xi);
        
    dXdXi(this_N, xi, 0, 1, 1, box);
    dXdXi(this_N, xi, 1, 0, 1, box);
    dXdXi(this_N, xi, 2, 1, 0, box);
    dXdXi(this_N, xi, 3, 0, 0, box);
    this_N.negate(1,2);
        
    if ( isAxisymmetric() ) {
       BoxIterator bit(box);
       for (bit.begin();bit.ok();++bit) {
         const IntVect& iv = bit();
                
          RealVect this_xi;
          this_xi[RADIAL_DIR] = xi(iv,RADIAL_DIR);
          this_xi[POLOIDAL_DIR] = xi(iv,POLOIDAL_DIR);
                
          double TwoPiRmaj = 2. * Pi * majorRadius(this_xi);
                
          for (int comp=0; comp<4; ++comp) {
             this_N(iv,comp) *= TwoPiRmaj;
          }
       }
    }
  }
    
  RealVect faceArea = getMappedFaceArea();
  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& this_N = a_N[dir];
      this_N.mult(faceArea[dir]);
    }

}

void
SingleNullBlockCoordSys::getFaceCenteredRealCoords( const int a_dir, FArrayBox& a_x ) const
{
   const Box& box( a_x.box() );
   FArrayBox xi_array( box, SpaceDim );
   getFaceCenteredMappedCoords( a_dir, xi_array );
   
   for (BoxIterator bit( box ); bit.ok(); ++bit) {
      IntVect iv( bit() );
      RealVect xi;
      for (int dir(0); dir<SpaceDim; ++dir) {
         xi[dir] = xi_array( iv, dir );
      }
      
      RealVect xi_hi(xi);
      RealVect xi_lo(xi);
      int perp_dir = (a_dir + 1) % 2;
      xi_hi[perp_dir] = xi[perp_dir] + 0.5 * m_dx[perp_dir];
      xi_lo[perp_dir] = xi[perp_dir] - 0.5 * m_dx[perp_dir];
      
      RealVect real_loc_hi( realCoord( xi_hi ) );
      RealVect real_loc_lo( realCoord( xi_lo ) );
      for (int dir(0); dir<SpaceDim; ++dir) {
         a_x( iv, dir ) = 0.5 * (real_loc_hi[dir] + real_loc_lo[dir]);
      }
   }
}

void
SingleNullBlockCoordSys::cellVol(FArrayBox&     a_vol,
                                 const FluxBox& a_N,
                                 const Box&     a_box ) const
{
   
   FArrayBox J(a_box,1);
   FArrayBox Xi(a_box,2);
   getCellCenteredMappedCoords(Xi);
   
   BoxIterator bit(a_box);
   for (bit.begin();bit.ok();++bit) {
      const IntVect& iv = bit();
      RealVect this_xi;
      this_xi[RADIAL_DIR] = Xi(iv,RADIAL_DIR);
      this_xi[POLOIDAL_DIR] = Xi(iv,POLOIDAL_DIR);
      
      RealVect Xi00(this_xi[0]-0.5*m_dx[0], this_xi[1]-0.5*m_dx[1]);
      RealVect X00(realCoord(Xi00));
      
      RealVect Xi10(this_xi[0]+0.5*m_dx[0], this_xi[1]-0.5*m_dx[1]);
      RealVect X10(realCoord(Xi10));
      
      RealVect Xi01(this_xi[0]-0.5*m_dx[0], this_xi[1]+0.5*m_dx[1]);
      RealVect X01(realCoord(Xi01));
      
      RealVect Xi11(this_xi[0]+0.5*m_dx[0], this_xi[1]+0.5*m_dx[1]);
      RealVect X11(realCoord(Xi11));
      
      double area_lo = 0.5 * fabs(X00[0]*(X01[1]-X10[1])+X01[0]*(X10[1]-X00[1])+X10[0]*(X00[1]-X01[1]));
      double area_hi = 0.5 * fabs(X11[0]*(X10[1]-X01[1])+X10[0]*(X01[1]-X11[1])+X01[0]*(X11[1]-X10[1]));
      
      
      a_vol(iv,0) = (area_lo + area_hi) * 2. * Pi * majorRadius(this_xi);
      
   }
   
}

void
SingleNullBlockCoordSys::init( ParmParse& a_pp )
{
   
   //Default values
   m_a =  1.2;
   m_b =  0.9;
   m_c =  0.7;
   m_R0 = 1.6;
   m_Zx = -acos(m_c/m_b);
   m_Zc  = acos(m_c/m_b);
   m_Z0  = 0.4-m_Zx;
   m_core_width = 0.066;
   m_sol_width  = 0.066;
   m_pf_width   = 0.066;
   m_pol_ref_length =  PI;
   m_pol_ref_length_mp = PI/2.0;
   m_div_leg_length = 1.0/4.0 * PI;
   m_pol_ref_z_lo = m_Zx+0.2;
   
   //Get the geometry parameters
   a_pp.query( "shape_factor_a",                   m_a );
   a_pp.query( "shape_factor_b",                   m_b );
   a_pp.query( "shape_factor_c",                   m_c );
   a_pp.query( "core_width",                       m_core_width );
   a_pp.query( "sol_width",                        m_sol_width );
   a_pp.query( "pf_width",                         m_pf_width );
   a_pp.query( "divertor_leg_poloidal_length",     m_div_leg_length );
   a_pp.query( "reference_poloidal_plane",         m_pol_ref_z_lo );
   
   //Important to define this solver parameters before doing any calculations!
   m_max_iterations = 50000;
   m_tol = 1.1e-7;
   
   // Get the poloidal field component scale factor
   if (a_pp.contains("Bpol_scale")) {
      a_pp.get("Bpol_scale", m_BpScale);
   }
   else {
      m_BpScale = m_RB_toroidal / m_R0;
   }
   
   normalizeBpScale();
   
   //Derived parameters
   double psiVal = psiAtR0(m_Zx+m_Z0);
   m_Zsep_hi = getZhi(psiVal) - m_Z0;
   m_pol_ref_psi = psiAtR0(m_pol_ref_z_lo + m_Z0);
   
   //Create an interpolator for the cylindrical coordinates
   //Here, we just use this interpolation utility as an effective storage for cell vertices and cell ceneters values (computed analytically)
   //(otherwise we need to compute coordinates (by solving algebraic equations) each time step, and it takes time)
   string interp_method = "bspline";
   if (a_pp.contains("interpolation_method")) {
      a_pp.get("interpolation_method", interp_method);
   }
   
   int m_ghost = 4;
   Box interp_box(2 * (m_domain.domainBox().smallEnd() - m_ghost * IntVect::Unit), 2 * ( m_domain.domainBox().bigEnd() + (m_ghost + 1) * IntVect::Unit));
   
   FArrayBox interp_node_coords(interp_box, SpaceDim);
   
   FArrayBox RZ_data(interp_box, 2);
   
   BoxIterator bit(interp_box);
   for (bit.begin();bit.ok();++bit) {
      IntVect iv = bit();
      RealVect mapped_coord(0.5 * m_dx[0] * iv[0], 0.5 *  m_dx[1] * iv[1]);
      
      RZ_data(iv,0) = realCoordPointwise(mapped_coord)[0];
      RZ_data(iv,1) = realCoordPointwise(mapped_coord)[1];
      
      interp_node_coords(iv,0) = mapped_coord[0];
      interp_node_coords(iv,1) = mapped_coord[1];
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
SingleNullBlockCoordSys::realCoordPointwise( const RealVect& a_xi ) const
{
   
   double epsilon = m_tol;
   RealVect x;
   
   //In computations of dXdXi we may slighly (with the machine presicion accuracy) step over the block boundaries (e.g., lower core boundary) in the mapped coordinate.
   //Becasue realCoord() function have difficulty near all block boundaries (the theta function is undefined at the block boundaries), we want to make sure that the
   //small overstepping is anihilated and special treatment of the block boundaries (see below) is used
   RealVect a_xi_tmp(a_xi);
   RealVect vicinity = epsilon * m_dx;
   for (int dir = 0; dir < SpaceDim; dir++) {
      if (a_xi[dir] > (lowerMappedCoordinate(dir) - vicinity[dir]) && a_xi[dir] < (lowerMappedCoordinate(dir) + vicinity[dir])) a_xi_tmp[dir]=lowerMappedCoordinate(dir);
      if (a_xi[dir] > (upperMappedCoordinate(dir) - vicinity[dir]) && a_xi[dir] < (upperMappedCoordinate(dir) + vicinity[dir])) a_xi_tmp[dir]=upperMappedCoordinate(dir);
   }
   //Same argumets apply for the midplane where theta function is undefined
   double midPlaneMappedPolCoord = 0.5*(lowerMappedCoordinate(POLOIDAL_DIR) + upperMappedCoordinate(POLOIDAL_DIR));
   if (a_xi[1] > (midPlaneMappedPolCoord - vicinity[1]) && a_xi[1] < (midPlaneMappedPolCoord + vicinity[1])) a_xi_tmp[1]=midPlaneMappedPolCoord;
   
   RealVect mapped_block_width;
   mapped_block_width[RADIAL_DIR] = upperMappedCoordinate(RADIAL_DIR) - lowerMappedCoordinate(RADIAL_DIR);
   mapped_block_width[POLOIDAL_DIR] = upperMappedCoordinate(POLOIDAL_DIR) - lowerMappedCoordinate(POLOIDAL_DIR);
   
   if (m_block_type == RCORE || m_block_type == LCORE) {
      
      double z_hi = (m_Zsep_hi - m_core_width + m_Z0) + m_core_width * (a_xi_tmp[RADIAL_DIR]-lowerMappedCoordinate(RADIAL_DIR))/mapped_block_width[RADIAL_DIR];
      double psiVal = psiAtR0(z_hi);
      
      double length = m_pol_ref_length * (a_xi_tmp[POLOIDAL_DIR]-lowerMappedCoordinate(POLOIDAL_DIR)) / mapped_block_width[POLOIDAL_DIR];
      if (m_block_type == LCORE) length = m_pol_ref_length * (upperMappedCoordinate(POLOIDAL_DIR) - a_xi_tmp[POLOIDAL_DIR]) / mapped_block_width[POLOIDAL_DIR];
      
      if (length < 0) length = -length;
      if (length > m_pol_ref_length) length = m_pol_ref_length - (length - m_pol_ref_length);
      
      if (length < epsilon) {
         if (a_xi_tmp[RADIAL_DIR] <= upperMappedCoordinate(RADIAL_DIR) + epsilon) {
            x[0] = m_R0 + epsilon;
            x[1] = getZlo(psiVal);
         }
         else {
            x[0] = 2.0*m_R0 - getXlo(psiVal);
            x[1] = getZlo(psiVal);
         }
      }
      
      else if ((2.0 * length > (m_pol_ref_length - epsilon)) && (2.0 * length < (m_pol_ref_length + epsilon))) {
         x[0] = m_R0 + acos(psiVal - m_b * sin(m_Zc) + m_c * m_Zc)/m_a;
         x[1] = m_Zc + m_Z0;
      }
      
      else if (length > m_pol_ref_length - epsilon) {
         x[0] = m_R0 + epsilon;
         x[1] = z_hi;
      }
      
      else {
         
         RealVect refPolPlaneInters = arcLengthInverse(length);
         double thetaVal = theta(refPolPlaneInters);
         int side = (length < m_pol_ref_length_mp) ? -1 : 1;
         x = gridLinesIntersection(psiVal, thetaVal, side);
         
         if ( (a_xi_tmp[POLOIDAL_DIR]>upperMappedCoordinate(POLOIDAL_DIR)) || (a_xi_tmp[POLOIDAL_DIR]<lowerMappedCoordinate(POLOIDAL_DIR)) ) {
            x[0] = m_R0 - (x[0]-m_R0);
         }
      }
      
      if ( m_block_type == LCORE ) x[0] = m_R0 - (x[0]-m_R0);
      
   }
   
   if (m_block_type == RCSOL || m_block_type == LCSOL) {
      
      double z_hi = (m_Zsep_hi + m_Z0) + m_sol_width * (a_xi_tmp[RADIAL_DIR]-lowerMappedCoordinate(RADIAL_DIR))/mapped_block_width[RADIAL_DIR];
      double psiVal = psiAtR0(z_hi);
      
      double length = m_pol_ref_length * (a_xi_tmp[POLOIDAL_DIR]-lowerMappedCoordinate(POLOIDAL_DIR)) / mapped_block_width[POLOIDAL_DIR];
      if (m_block_type == LCSOL) length = m_pol_ref_length * (upperMappedCoordinate(POLOIDAL_DIR) - a_xi_tmp[POLOIDAL_DIR]) / mapped_block_width[POLOIDAL_DIR];
      
      if (length < 0) length = - length;
      if (length > m_pol_ref_length) length = m_pol_ref_length - (length - m_pol_ref_length);
      
      if (length < epsilon) {
         if (a_xi_tmp[RADIAL_DIR] >= lowerMappedCoordinate(RADIAL_DIR) - epsilon) {
            x[1] = getZlo(psiVal);
            x[0] = getXlo(psiVal);
         }
         else {
            x[0] = m_R0 + epsilon;
            x[1] = getZlo(psiVal);
         }
         
         if ( m_block_type == LCSOL ) x[0] = m_R0 - (x[0]-m_R0);
      }

      else if ((2.0 * length > (m_pol_ref_length - epsilon)) && (2.0 * length < (m_pol_ref_length + epsilon))) {
         
         x[0] = m_R0 + acos(psiVal - m_b * sin(m_Zc) + m_c * m_Zc)/m_a;
         x[1] = m_Zc + m_Z0;
         
         if ( m_block_type == LCSOL ) x[0] = m_R0 - (x[0]-m_R0);
      }
      
      else if (length > m_pol_ref_length - epsilon) {
         x[0] = m_R0 + epsilon;
         x[1] = z_hi;
         
         if ( m_block_type == LCSOL ) x[0] = m_R0 - (x[0]-m_R0);
      }
      
      else {
         
         RealVect refPolPlaneInters = arcLengthInverse(length);
         double thetaVal = theta(refPolPlaneInters);
         int side = (length < m_pol_ref_length_mp) ? -1 : 1;
         x = gridLinesIntersection(psiVal, thetaVal, side);
         
         if ( ((a_xi_tmp[POLOIDAL_DIR] <= upperMappedCoordinate(POLOIDAL_DIR)) && (a_xi_tmp[POLOIDAL_DIR]>=lowerMappedCoordinate(POLOIDAL_DIR))) &&  m_block_type == LCSOL) {
            x[0] = m_R0 - (x[0]-m_R0);
         }
         
         if ( (a_xi_tmp[POLOIDAL_DIR] > upperMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == RCSOL) x[0] = m_R0 - (x[0]-m_R0);
         if ( (a_xi_tmp[POLOIDAL_DIR] < lowerMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == RCSOL) x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx);
         if ( (a_xi_tmp[POLOIDAL_DIR] > upperMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == LCSOL) {
            x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx);
            x[0] = m_R0 - (x[0]-m_R0);
         }
      }
   }
   
   
   if (m_block_type == RSOL || m_block_type == LSOL) {
      
      double z_hi = (m_Zsep_hi + m_Z0) + m_sol_width * (a_xi_tmp[RADIAL_DIR]-lowerMappedCoordinate(RADIAL_DIR))/mapped_block_width[RADIAL_DIR];
      double psiVal = psiAtR0(z_hi);
      
      double length = m_div_leg_length * (upperMappedCoordinate(POLOIDAL_DIR) - a_xi_tmp[POLOIDAL_DIR]) / mapped_block_width[POLOIDAL_DIR];
      if (m_block_type == LSOL) length = m_div_leg_length * (a_xi_tmp[POLOIDAL_DIR] - lowerMappedCoordinate(POLOIDAL_DIR)) / mapped_block_width[POLOIDAL_DIR];
      if (length < 0) length = - length;
      
      RealVect refPolPlaneInters = arcLengthInverse(length);
      double thetaVal = theta(refPolPlaneInters);
      
      if (length < epsilon) {
         if (a_xi_tmp[RADIAL_DIR] >= lowerMappedCoordinate(RADIAL_DIR) - epsilon) {
            x[1] = getZlo(psiVal);
            x[0] = getXlo(psiVal);
         }
         else {
            x[1] = getZlo(psiVal);
            x[0] = m_R0 + epsilon;
         }
         
         if ( m_block_type == LSOL ) x[0] = m_R0 - (x[0]-m_R0);
         
         x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx);
      }
      
      else if ((2.0 * length > (m_pol_ref_length - epsilon)) && (2.0 * length < (m_pol_ref_length + epsilon))) {
	    x[0] = m_R0 + acos(psiVal - m_b * sin(m_Zc) + m_c * m_Zc)/m_a;
	    x[1] = m_Zc + m_Z0;
	    if ( m_block_type == LSOL ) x[0] = m_R0 - (x[0]-m_R0);
	    x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx); 
      }

      else {
         
         RealVect refPolPlaneInters = arcLengthInverse(length);
         double thetaVal = theta(refPolPlaneInters);
         int side = (length < m_pol_ref_length_mp) ? -1 : 1;
         x = gridLinesIntersection(psiVal, thetaVal, side);
         
         x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx);
         
         if ( (a_xi_tmp[POLOIDAL_DIR]>=lowerMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == LSOL) {
            x[0] = m_R0 - (x[0]-m_R0);
         }
         
         else if ((a_xi_tmp[POLOIDAL_DIR]<lowerMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == LSOL) {
            x[0] = m_R0 - (x[0]-m_R0);
            x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx); //another reflection to bring the vertical coordinate back to upper blocks
         }
         else if ((a_xi_tmp[POLOIDAL_DIR]>upperMappedCoordinate(POLOIDAL_DIR)) &&  m_block_type == RSOL) {
            x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx); //another reflection to bring the vertical coordinate back to upper blocks
         }
         
      }
      
   }
   
   if (m_block_type == RPF || m_block_type == LPF) {
      
      double z_hi = (m_Zsep_hi - m_pf_width + m_Z0) + m_pf_width * (a_xi_tmp[RADIAL_DIR]-lowerMappedCoordinate(RADIAL_DIR))/mapped_block_width[RADIAL_DIR];
      double psiVal = psiAtR0(z_hi);
      
      double length = m_div_leg_length * (upperMappedCoordinate(POLOIDAL_DIR) - a_xi_tmp[POLOIDAL_DIR]) / mapped_block_width[POLOIDAL_DIR];
      if (m_block_type == LPF) length = m_div_leg_length * (a_xi_tmp[POLOIDAL_DIR] - lowerMappedCoordinate(POLOIDAL_DIR)) / mapped_block_width[POLOIDAL_DIR];
      if (length < 0) length = - length;
      
      if (length < epsilon) {
         if (a_xi_tmp[RADIAL_DIR] <= upperMappedCoordinate(RADIAL_DIR) + epsilon) {
            x[0] = m_R0 + epsilon;
            x[1] = getZlo(psiVal);	    
         }
         else {
            x[0] = 2.0*m_R0 - getXlo(psiVal);
            x[1] = getZlo(psiVal);
         }         
      }

      else if ((2.0 * length > (m_pol_ref_length - epsilon)) && (2.0 * length < (m_pol_ref_length + epsilon))) {
         x[0] = m_R0 + acos(psiVal - m_b * sin(m_Zc) + m_c * m_Zc)/m_a;
         x[1] = m_Zc + m_Z0;
      }
      
      else {
         
         RealVect refPolPlaneInters = arcLengthInverse(length);
         double thetaVal = theta(refPolPlaneInters);
         int side = (length < m_pol_ref_length_mp) ? -1 : 1;
         x = gridLinesIntersection(psiVal, thetaVal, side);
         
         if ( ((a_xi_tmp[POLOIDAL_DIR]>upperMappedCoordinate(POLOIDAL_DIR)) && m_block_type == RPF)
             || ((a_xi_tmp[POLOIDAL_DIR]<lowerMappedCoordinate(POLOIDAL_DIR)) && m_block_type == LPF ) ) {
            
            x[0] = m_R0 - (x[0]-m_R0);
         }
         
      }
      
      x[1] = (m_Z0 +m_Zx) - (x[1]-m_Z0-m_Zx);
      
      if ( m_block_type == LPF ) x[0] = m_R0 - (x[0]-m_R0);
      
   }
   
   return x;
}


//Keep it so we not break the code structure, but make trivial
RealVect
SingleNullBlockCoordSys::mappedCoord( const RealVect& a_X ) const
{
#ifdef PLOT_INVERSE_POINTS
   fprintf(m_ipt_file[m_block_type], "%20.12e %20.12e\n", a_X[0], a_X[1]);
#endif
   
   return a_X;
}

Real
SingleNullBlockCoordSys::dXdXi( const RealVect& a_Xi,
                               int             a_dirX,
                               int             a_dirXi ) const
{
   
   RealVect a_Xi_tmp(a_Xi);
   RealVect x = realCoord(a_Xi_tmp);
   
   RealVect Xi_hi(a_Xi_tmp);
   RealVect Xi_lo(a_Xi_tmp);
   
   Xi_hi[a_dirXi] = a_Xi[a_dirXi] + 0.5 * m_dx[a_dirXi];
   Xi_lo[a_dirXi] = a_Xi[a_dirXi] - 0.5 * m_dx[a_dirXi];
   
   RealVect x_hi = realCoord(Xi_hi);
   RealVect x_lo = realCoord(Xi_lo);
   
   double metric_factor;
   //Correct metrics to preserve the freestream property
   if ((isAxisymmetric()) && a_dirX == RADIAL_DIR) {
      
      double fac = 0.5 * (x_hi[a_dirX] + x_lo[a_dirX])/x[a_dirX];
      metric_factor = (x_hi[a_dirX] - x_lo[a_dirX])/m_dx[a_dirXi] * fac;
   }
   
   else {
      
      metric_factor = (x_hi[a_dirX] - x_lo[a_dirX])/m_dx[a_dirXi];
   }
   
   return metric_factor;
}


void
SingleNullBlockCoordSys::dXdXi( FArrayBox&       a_dXdXi,
                               const FArrayBox& a_Xi,
                               int              a_destComp,
                               int              a_dirX,
                               int              a_dirXi,
                               const Box&       a_box ) const
{
   BoxIterator bit(a_box);
   for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect this_Xi(a_Xi(iv,0),a_Xi(iv,1));
      a_dXdXi(iv,a_destComp) = dXdXi(this_Xi,a_dirX,a_dirXi);
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
SingleNullBlockCoordSys::computeFieldData( const int  a_dir,
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
   
   FArrayBox RZ(box,SpaceDim);
   if (a_dir>=0 && a_dir<SpaceDim) {
      getFaceCenteredRealCoords(a_dir, RZ);
   }
   else {
      getCellCenteredRealCoords(RZ);
   }
   
   
   for (BoxIterator bit( box ); bit.ok(); ++bit) {
      IntVect iv( bit() );
      
      //Compute unshifted coordinate
      double r = RZ(iv,0)-m_R0;
      double z = ((RZ(iv,1)-m_Z0)>m_Zx) ? (RZ(iv,1)-m_Z0) : (m_Zx + (m_Zx+m_Z0 - RZ(iv,1)));
      
      //We only have psi and thus B defined in the upper blocks (z>m_Zx), and need to use symmetry for the lower blocks.
      double sign = ((RZ(iv,1)-m_Z0)>m_Zx) ? 1.0 : -1.0;
      
      //Compute field data
      a_BField(iv,0) = sign * m_BpScale * (m_c - m_b * cos(z))/(r + m_R0);
      a_BField(iv,1) = m_RB_toroidal/(r + m_R0);
      a_BField(iv,2) = - m_a * m_BpScale * sin(m_a * r) /(r + m_R0);
      
      a_BFieldMag(iv,0) = sqrt(pow(a_BField(iv,0),2) + pow(a_BField(iv,1),2) + pow(a_BField(iv,2),2));
      
      a_BFieldDir(iv,0) = a_BField(iv,0)/a_BFieldMag(iv,0);
      a_BFieldDir(iv,1) = a_BField(iv,1)/a_BFieldMag(iv,0);
      a_BFieldDir(iv,2) = a_BField(iv,2)/a_BFieldMag(iv,0);
      
      double coeff = pow(m_BpScale*m_c,2) + pow(m_RB_toroidal,2) - 2.0 * m_b * pow(m_BpScale,2) * m_c * cos(z) + pow(m_b * m_BpScale * cos(z), 2) + pow(m_a * m_BpScale * sin(m_a * r),2);
      
      a_gradBFieldMag(iv,0) = (-coeff + 1.0/2.0 * pow(m_a,3) * pow(m_BpScale,2) * (r + m_R0) * sin(2.0 * m_a * r)) / (pow(r+m_R0,2)*sqrt(coeff));
      a_gradBFieldMag(iv,1) = 0.0;
      a_gradBFieldMag(iv,2) = sign * m_b * pow(m_BpScale,2) * (m_c - m_b * cos(z)) * sin(z)/ ((r+m_R0)*sqrt(coeff));
      
      a_curlBFieldDir(iv, 0) = sign * m_b * pow(m_BpScale,2) * m_RB_toroidal * (m_c - m_b * cos(z)) * sin(z) / pow(coeff,3.0/2.0);
      
      a_curlBFieldDir(iv, 1) = m_BpScale * (m_a * m_a * cos(m_a * r) * (pow(m_BpScale * m_c,2) + pow(m_RB_toroidal,2) - 2.0 * m_b * pow(m_BpScale,2) * m_c * cos(z) + pow(m_b * m_BpScale * cos(z),2)) + m_b * (pow(m_RB_toroidal,2) + pow(m_a * m_BpScale * sin(m_a * r),2)) * sin(z)) / pow(coeff,3.0/2.0);
      
      a_curlBFieldDir(iv, 2) = m_RB_toroidal * (coeff - 1.0/2.0 * pow(m_a,3) * pow(m_BpScale,2) * (r+m_R0) * sin(2.0 * m_a * r)) / ((r+m_R0) * pow(coeff,3.0/2.0));
      
      a_BFieldDirdotcurlBFieldDir(iv, 0) =  m_BpScale * m_RB_toroidal * (m_a * m_a * (r+m_R0) * cos(m_a *r) - m_a * sin(m_a * r) + m_b * (r+m_R0) * sin(z)) / ((r+m_R0) * coeff);
      
   }
}



void
SingleNullBlockCoordSys::getMagneticFlux( const FArrayBox& a_physical_coordinates,
                                          FArrayBox&       a_magnetic_flux ) const
{
   MayDay::Error("SingleNullBlockCoordSys::getMagneticFlux() not yet implemented for model geometry");
}



double
SingleNullBlockCoordSys::getMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   MayDay::Error("SingleNullBlockCoordSys::getMagneticFlux() not yet implemented for model geometry");

   return 0.;
}


RealVect
SingleNullBlockCoordSys::gridLinesIntersection( const double psiVal, const double thetaVal, const int side ) const
{
   
   // We do not solve here for the intersections on the "critical edges"
   //(i.e., vertical, horizontal -X , and horizonta - midPlane)
   // The critical edges are directly handled in the realCoord function
   double epsilon = m_tol;
   
   bool residual_tolerance_satisfied = false;
   
   RealVect half_Z_sol, tmp_lo_Z_sol, tmp_hi_Z_sol;
   double f_lo_Z, f_hi_Z, f_half_Z, half_Z;
   
   double tmp_lo_Z = getZlo(psiVal);
   double tmp_hi_Z = getZhi(psiVal);
   
   if (side==-1) tmp_hi_Z = m_Zc + m_Z0 - epsilon;
   if (side==1) tmp_lo_Z = m_Zc + m_Z0 + epsilon;
   
   int iter = 0;
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
      
      half_Z = 0.5 * (tmp_lo_Z + tmp_hi_Z);
      half_Z_sol[0] = m_R0 + acos(psiVal - m_b * sin(half_Z - m_Z0) + m_c * (half_Z - m_Z0)) / m_a;
      half_Z_sol[1] = half_Z;
      
      tmp_lo_Z_sol[0] = m_R0 + acos(psiVal - m_b * sin(tmp_lo_Z - m_Z0) + m_c * (tmp_lo_Z - m_Z0)) / m_a;
      tmp_lo_Z_sol[1] = tmp_lo_Z;
      
      tmp_hi_Z_sol[0] = m_R0 + acos(psiVal - m_b * sin(tmp_hi_Z - m_Z0) + m_c * (tmp_hi_Z - m_Z0)) / m_a;
      tmp_hi_Z_sol[1] = tmp_hi_Z;
      
      
      f_lo_Z = theta(tmp_lo_Z_sol) - thetaVal;
      f_hi_Z = theta(tmp_hi_Z_sol) - thetaVal;
      f_half_Z = theta(half_Z_sol) - thetaVal;
      
      if (f_lo_Z * f_half_Z < 0) {
         tmp_hi_Z = half_Z;
      }
      
      else {
         tmp_lo_Z = half_Z;
      }
      
      if ( (tmp_hi_Z - tmp_lo_Z) < m_tol ) residual_tolerance_satisfied = true;
      
   }
   
   
   RealVect x;
   
   D_TERM(x[0] = acos(psiVal - m_b * sin(half_Z - m_Z0) + m_c * (half_Z - m_Z0)) / m_a + m_R0;,
          x[1] = half_Z;,
          x[2] = 0.;)
   
   return x;
}

double SingleNullBlockCoordSys::getXlo(double psiVal ) const
{
   CH_assert(psiVal <= psiAtR0(m_Zsep_hi + m_Z0));
   double epsilon = m_tol;
   double x;
   
   if (psiVal == psiAtR0(m_Zsep_hi + m_Z0)) {
      x = m_R0 + epsilon;
   }
   
   else {
      x= m_R0 + acos(psiVal - m_b * sin(m_Zx) + m_c * m_Zx)/m_a;
   }
   
   return x;
}

double SingleNullBlockCoordSys::getZlo(double psiVal ) const
{
   
   bool residual_tolerance_satisfied = false;
   
   double epsilon = m_tol;
   double half_Z;
   
   if (psiVal < psiAtR0(m_Zsep_hi + m_Z0)) {
      half_Z = m_Zx + m_Z0;
      residual_tolerance_satisfied = true;
   }
   
   double tmp_lo_Z = m_Zx + m_Z0;
   double tmp_hi_Z = m_Zc + m_Z0;
   double f_lo_Z, f_hi_Z, f_half_Z;
   
   int iter = 0;
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
      
      half_Z = 0.5 * (tmp_lo_Z + tmp_hi_Z);
      f_lo_Z = psiAtR0(tmp_lo_Z) - psiVal;
      f_hi_Z = psiAtR0(tmp_hi_Z) - psiVal;
      f_half_Z = psiAtR0(half_Z) - psiVal;
      
      if (f_lo_Z * f_half_Z < 0) {
         tmp_hi_Z = half_Z;
      }
      
      else {
         tmp_lo_Z = half_Z;
      }
      
      if ( (psiAtR0(half_Z + m_tol) - psiVal) * (psiAtR0(half_Z - m_tol) - psiVal) < 0 ) {
         residual_tolerance_satisfied = true;
      }
      
   }
   
   // add epsilon to move away from z=Zx (where theta is undefined) and also to make sure that for x_unshifted~0, there will be a solution for x.
   return (half_Z + epsilon);
}

double SingleNullBlockCoordSys::getZhi(double psiVal ) const
{
   bool residual_tolerance_satisfied = false;
   
   double tmp_lo_Z = m_Zc + m_Z0;
   double tmp_hi_Z = (3.0*m_Zc + m_Z0); //assume the simulation domain with the ghost extensions be smaller than tmp_hi_Z
   double half_Z, f_lo_Z, f_hi_Z, f_half_Z;
   
   int iter = 0;
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
      
      half_Z = 0.5 * (tmp_lo_Z + tmp_hi_Z);
      f_lo_Z = psiAtR0(tmp_lo_Z) - psiVal;
      f_hi_Z = psiAtR0(tmp_hi_Z) - psiVal;
      f_half_Z = psiAtR0(half_Z) - psiVal;
      
      if (f_lo_Z * f_half_Z < 0) {
         tmp_hi_Z = half_Z;
      }
      
      else {
         tmp_lo_Z = half_Z;
      }
      
      if ( (psiAtR0(half_Z + m_tol) - psiVal) * (psiAtR0(half_Z - m_tol) - psiVal) < 0 ) {
         residual_tolerance_satisfied = true;
      }
      
   }
   
   //Check if our solution is outside the analytical limit (i.e., the subsequent evalution of theta fails)
   if ( (psiAtR0(half_Z) - psiVal) < 0.0 ) half_Z += -m_tol;
   
   return half_Z;
}

RealVect
SingleNullBlockCoordSys::arcLengthInverse(const double length) const
{
   //Here, length is the poloidal angle that goes counter clockwise with zero being at the bottom core boundary
   
   RealVect x;
   
   double epsilon = m_tol;
   
   double psiVal = m_pol_ref_psi;
   double z_polar_axis = m_Zc + m_Z0; //shifted from the toakamak axis where theta is not defined
   
   bool residual_tolerance_satisfied = false;
   
   double half_r, tmp_lo_r, tmp_hi_r, f_lo_r, f_hi_r, f_half_r, x_lo, x_hi, half_x, z_lo, z_hi, half_z;
   
   if ( length > (PI/2.0 - epsilon) && length < (PI/2.0 + epsilon) ) {
      x[0] = acos(m_pol_ref_psi - m_b * sin(z_polar_axis) + m_c * z_polar_axis)/m_a + m_R0;
      x[1] = z_polar_axis;
      residual_tolerance_satisfied = true;
   }
   
   tmp_lo_r = 0.0;
   if ( length > PI/2.0 ) {
      tmp_hi_r = 3.0*m_Zc; //assume that this radius encompasses the top part of the domain
   }
   
   if ( length < PI/2.0 ) {
      tmp_hi_r = fabs(m_Zx + m_Z0 - z_polar_axis); //assume that this radius encompasses the lower half of the domain
   }
   
   int iter = 0;
   while ( iter++ < m_max_iterations && !residual_tolerance_satisfied ) {
      
      x_lo = tmp_lo_r * sin(length) + m_R0;
      z_lo = z_polar_axis - tmp_lo_r * cos(length);
      f_lo_r = psi(RealVect(x_lo,z_lo)) - psiVal;
      
      x_hi = tmp_hi_r * sin(length) + m_R0;
      z_hi = z_polar_axis - tmp_hi_r * cos(length);
      f_hi_r = psi(RealVect(x_hi,z_hi)) - psiVal;
      
      half_r = 0.5 * (tmp_lo_r + tmp_hi_r);
      half_x = half_r * sin(length) + m_R0;
      half_z = z_polar_axis - half_r * cos(length);
      f_half_r = psi(RealVect(half_x,half_z)) - psiVal;
      
      if (f_lo_r * f_half_r < 0) {
         tmp_hi_r = half_r;
      }
      
      else {
         tmp_lo_r = half_r;
      }
      
      double x_test_hi = (half_r+m_tol) * sin(length) + m_R0;
      double x_test_lo = (half_r-m_tol) * sin(length) + m_R0;
      double z_test_hi = z_polar_axis - (half_r + m_tol) * cos(length);
      double z_test_lo = z_polar_axis - (half_r - m_tol) * cos(length);
      
      if ( (psi(RealVect(x_test_hi,z_test_hi)) - psiVal) * (psi(RealVect(x_test_lo,z_test_lo)) - psiVal) < 0 ) {
         residual_tolerance_satisfied = true;
      }
      
      x[0] = half_r * sin(length) + m_R0;
      x[1] = z_polar_axis - half_r * cos(length);
      
   }
   
   return x;
   
}


//NB: psi, theta and derivatives experssions are used only for z>m_Zx and symmetry is used to obtain those for z<m_Zx (we assume that psi is symmetric)
double SingleNullBlockCoordSys::psi(const RealVect& a_x) const
{
   double x = a_x[0] - m_R0;
   double z = a_x[1] - m_Z0;
   // CH_assert(z >= m_Zx);
   double psiMin = m_b * sin(m_Zc) - m_c * m_Zc;
   double psi = cos(m_a * x) + m_b * sin(z) - m_c * z;
   
   //Bacause psi is periodic (though periodicity is seen only far outside the domain)
   //we flatten it out here, baising to a constant after the "first period boundary".
   //This is the extra cautionary measure for various Newton solvers.
   double psiTruncated = (psi < psiMin) ? psiMin : psi;
   return psi; //psiTruncated;
}

double SingleNullBlockCoordSys::theta(const RealVect& a_x) const
{
   double x = a_x[0] - m_R0;
   double z = a_x[1] - m_Z0;
   
   //CH_assert(z >= m_Zx);
   
   double Fx = -log(tan((m_a * abs(x))/2.0))/pow(m_a,2);
   
   double sqr = sqrt(m_b * m_b - m_c * m_c);
   double FyMath = (2.0/sqr) * atanh((m_b + m_c) * tan(z/2.0) / sqr);
   double FyGR = (1.0/sqr) * log((sqr*tan(z/2.0) + m_b - m_c)/(sqr*tan(z/2.0) - m_b + m_c));
   double Fy = ( abs(z) < abs(m_Zx) ) ? FyMath : FyGR;
   
   return Fx - Fy;
   
}

double SingleNullBlockCoordSys::psiAtR0(const double a_Z) const
{
   double z = a_Z - m_Z0;
   //CH_assert(z >= m_Zx);
   return 1.0 + m_b*sin(z) - m_c*(z);
}

double SingleNullBlockCoordSys::psiZ(const RealVect& a_x) const
{
   double z = a_x[1] - m_Z0;
   //CH_assert(z >= m_Zx);
   return -m_c + m_b * cos(z);
}

double SingleNullBlockCoordSys::psiX(const RealVect& a_x) const
{
   double x = a_x[0] - m_R0;
   return -m_a * sin(m_a * x);
}

double SingleNullBlockCoordSys::thetaZ(const RealVect& a_x) const
{
   double z = a_x[1] - m_Z0;
   //CH_assert(z >= m_Zx);
   return 1.0/(m_c - m_b * cos(z));
}

double SingleNullBlockCoordSys::thetaX(const RealVect& a_x) const
{
   double x = a_x[0] - m_R0;
   return -1.0/(m_a*sin(m_a*x));
}

void SingleNullBlockCoordSys::normalizeBpScale() {
   
   double psiSep = psiAtR0(m_Zx + m_Z0);
   double sepTopZ = getZhi(psiSep);
   RealVect sepTop(m_R0, sepTopZ);
   m_BpScale *= m_R0/sqrt(pow(psiX(sepTop),2) + pow(psiZ(sepTop),2));
}


#endif

//======================================= END: THE CASE OF ANALYTIC MODEL GEOMETRY ==================================

#include "NamespaceFooter.H"
