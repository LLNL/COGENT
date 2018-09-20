#include <fstream>
#include <float.h>
#include "FieldData.H"
#include "FieldDataF_F.H"
#include "BoxIterator.H"
#include "IntVect.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

enum SingleNullBlockType {LCORE,RCORE,LCSOL,RCSOL,LSOL,RSOL,LPF,RPF,MCORE,MCSOL,NUM_SINGLE_NULL_BLOCKS};
enum interpolation_mode {FUNCTION, RADIAL_DERIVATIVE, POLOIDAL_DERIVATIVE};
enum field_vars {BR_VAR, BZ_VAR, NUM_FIELD_VARS};
enum rz_vars {R_VAR, Z_VAR};
enum directions {RADIAL_DIR, POLOIDAL_DIR};


FieldData::FieldData(const string& a_RZ_mapping_file_name,
                     const string& a_field_file_name,
                     const string& a_block_name)
   : m_RZ_interpolator(NULL),
     m_field_RZ_interpolator(NULL),
     m_field_interpolator(NULL)
{
   if ( a_block_name == "lcore" ) {
      m_block_number = LCORE;
   }
   else if ( a_block_name == "rcore" ) {
      m_block_number = RCORE;
   }
   else if ( a_block_name == "lcsol" ) {
      m_block_number = LCSOL;
   }
   else if ( a_block_name == "rcsol" ) {
      m_block_number = RCSOL;
   }
   else if ( a_block_name == "lsol" ) {
      m_block_number = LSOL;
   }
   else if ( a_block_name == "rsol" ) {
      m_block_number = RSOL;
   }
   else if ( a_block_name == "lpf" ) {
      m_block_number = LPF;
   }
   else if ( a_block_name == "rpf" ) {
      m_block_number = RPF;
   }
   else if ( a_block_name == "mcore" ) {
      m_block_number = MCORE;
   }
   else if ( a_block_name == "mcsol" ) {
      m_block_number = MCSOL;
   }
   else {
      MayDay::Error("FieldData: Unrecognized block name");
   }

   double RB_toroidal = 0.;

   readRZMapping(a_RZ_mapping_file_name, m_block_number);

#ifdef USE_DCT_FIELD
   readFieldFromDCTFile(a_field_file_name);
#else
   readFieldFromMappingFile(a_field_file_name, m_block_number);

   // Precompute the physical coordinates of a subset of points
   // for use in the Newton iteration used to compute the inverse
   // mapping
   precomputeCoarseMappingData(10);
#endif

}



FieldData::~FieldData()
{
   if (m_field_interpolator) delete m_field_interpolator;
   if (m_field_RZ_interpolator) delete m_field_RZ_interpolator;
   if (m_RZ_interpolator) delete m_RZ_interpolator;
}


void FieldData::readRZMapping(const string&   a_file_name,
                              const int       a_block_number)
{
   int num_block_elements;
   IntVect mapping_block_size;
   RealVect lower_mapped_coordinate = RealVect::Zero;
   RealVect upper_mapped_coordinate = RealVect::Unit;
   ifstream inFile;

   ParmParse pp;

#ifdef CH_MPI
   if (procID() == 0) {
#endif

      inFile.open( a_file_name.c_str() );

      if (!inFile) {
         cout << "Unable to open RZ mapping file" << endl;
         exit(1);
      }

      // Move the file pointer to the block of interest

      string block_name;
      for (int block=0; block<a_block_number; ++block) {

         if ( (inFile >> block_name) == NULL ) {
            cout << "EOF reached before finding block number " << a_block_number << endl;
            exit(1);
         }
         inFile >> n_radial_points;
         inFile >> n_radial_extend;
         inFile >> n_poloidal_points;
         inFile >> n_poloidal_extend;

         num_block_elements = (n_radial_points+2*n_radial_extend)*(n_poloidal_points+2*n_poloidal_extend);

         double dummy;
         for (int i=0; i<4*num_block_elements; ++i) {
            inFile >> dummy;
         }
      }

      if ( (inFile >> block_name) == NULL ) {
         cout << "EOF reached before finding block number " << a_block_number << endl;
         exit(1);
      }
      inFile >> n_radial_points;
      inFile >> n_radial_extend;
      inFile >> n_poloidal_points;
      inFile >> n_poloidal_extend;

#ifdef CH_MPI
   }

   MPI_Bcast(&n_radial_points, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n_radial_extend, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n_poloidal_points, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n_poloidal_extend, 1, MPI_INT, 0, MPI_COMM_WORLD);

#endif // end if MPI

   mapping_block_size[RADIAL_DIR]   = n_radial_points   + 2*n_radial_extend;
   mapping_block_size[POLOIDAL_DIR] = n_poloidal_points + 2*n_poloidal_extend;

   num_block_elements = mapping_block_size[RADIAL_DIR] * mapping_block_size[POLOIDAL_DIR];

   double *input_values = new double[4*num_block_elements];

#ifdef CH_MPI
   if (procID() == 0) {
#endif

      for (int i=0; i<4*num_block_elements; ++i) {
         inFile >> input_values[i];
      }

      inFile.close();

#ifdef CH_MPI
   }

   MPI_Bcast(input_values, 4*num_block_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif // end if MPI

   Box interp_box(IntVect::Zero,mapping_block_size-IntVect::Unit);

   FArrayBox interp_node_coords(interp_box, SpaceDim);

   // Create an interpolator for the cylindrical coordinates

   FArrayBox RZ_data(interp_box, 2);

   for (BoxIterator bit(interp_box);bit.ok();++bit) {
      IntVect iv = bit();
      int k = mapping_block_size[RADIAL_DIR]*iv[1] + iv[0];

      RZ_data(iv,0) = input_values[4*k];
      RZ_data(iv,1) = input_values[4*k+1];
   }

   Box interior_interp_box = interp_box;
   
   interior_interp_box.grow(RADIAL_DIR, -n_radial_extend);
   interior_interp_box.grow(POLOIDAL_DIR, -n_poloidal_extend);

   FArrayBox interior_RZ_data(interior_interp_box, 2);
   interior_RZ_data.copy(RZ_data);

   double dr = ( upper_mapped_coordinate[RADIAL_DIR]
                 - lower_mapped_coordinate[RADIAL_DIR] ) / (n_radial_points - 1);

   double dtheta = ( upper_mapped_coordinate[POLOIDAL_DIR]
                     - lower_mapped_coordinate[POLOIDAL_DIR] ) / (n_poloidal_points - 1);

   FArrayBox interior_interp_node_coords(interior_interp_box, SpaceDim);

   for (BoxIterator bit(interior_interp_box);bit.ok();++bit) {
      IntVect iv = bit();

      interior_interp_node_coords(iv,0) = (iv[0] - interior_interp_box.smallEnd(0)) * dr;
      interior_interp_node_coords(iv,1) = (iv[1] - interior_interp_box.smallEnd(1)) * dtheta;
   }

   m_RZ_interpolator = new BSplineInterp(pp, 2, interior_interp_node_coords, interior_RZ_data);

   delete [] input_values;
}


void FieldData::readFieldFromMappingFile(const string&   a_file_name,
                                         const int       a_block_number)
{
   int num_block_elements;
   IntVect mapping_block_size;
   RealVect lower_mapped_coordinate = RealVect::Zero;
   RealVect upper_mapped_coordinate = RealVect::Unit;
   ifstream inFile;

   ParmParse pp;

#ifdef CH_MPI
   if (procID() == 0) {
#endif

      inFile.open( a_file_name.c_str() );

      if (!inFile) {
         cout << "Unable to open field file" << endl;
         exit(1);
      }

      // Move the file pointer to the block of interest

      string block_name;
      for (int block=0; block<a_block_number; ++block) {

         if ( (inFile >> block_name) == NULL ) {
            cout << "EOF reached before finding block number " << a_block_number << endl;
            exit(1);
         }
         inFile >> n_radial_points;
         inFile >> n_radial_extend;
         inFile >> n_poloidal_points;
         inFile >> n_poloidal_extend;

         num_block_elements = (n_radial_points+2*n_radial_extend)*(n_poloidal_points+2*n_poloidal_extend);

         double dummy;
         for (int i=0; i<4*num_block_elements; ++i) {
            inFile >> dummy;
         }
      }

      if ( (inFile >> block_name) == NULL ) {
         cout << "EOF reached before finding block number " << a_block_number << endl;
         exit(1);
      }
      inFile >> n_radial_points;
      inFile >> n_radial_extend;
      inFile >> n_poloidal_points;
      inFile >> n_poloidal_extend;

#ifdef CH_MPI
   }

   MPI_Bcast(&n_radial_points, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n_radial_extend, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n_poloidal_points, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n_poloidal_extend, 1, MPI_INT, 0, MPI_COMM_WORLD);

#endif // end if MPI

   mapping_block_size[RADIAL_DIR]   = n_radial_points   + 2*n_radial_extend;
   mapping_block_size[POLOIDAL_DIR] = n_poloidal_points + 2*n_poloidal_extend;

   num_block_elements = mapping_block_size[RADIAL_DIR] * mapping_block_size[POLOIDAL_DIR];

   double *input_values = new double[4*num_block_elements];

#ifdef CH_MPI
   if (procID() == 0) {
#endif

      for (int i=0; i<4*num_block_elements; ++i) {
         inFile >> input_values[i];
      }

      inFile.close();

#ifdef CH_MPI
   }

   MPI_Bcast(input_values, 4*num_block_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif // end if MPI

   int num_radial_cells = n_radial_points - 1;

   double dr = ( upper_mapped_coordinate[RADIAL_DIR]
                 - lower_mapped_coordinate[RADIAL_DIR] ) / num_radial_cells;

   double *r_pts = new double[mapping_block_size[RADIAL_DIR]];

   r_pts[0] = lower_mapped_coordinate[RADIAL_DIR] - n_radial_extend*dr;
   for (int i=1; i<mapping_block_size[RADIAL_DIR]; ++i) {
      r_pts[i] = r_pts[i-1] + dr;
   }

   int num_poloidal_cells = n_poloidal_points - 1;

   double dtheta = ( upper_mapped_coordinate[POLOIDAL_DIR]
                     - lower_mapped_coordinate[POLOIDAL_DIR] ) / num_poloidal_cells;

   double *theta_pts = new double[mapping_block_size[POLOIDAL_DIR]];

   if (a_block_number == LCORE || a_block_number == LCSOL) {
      theta_pts[mapping_block_size[POLOIDAL_DIR]-1] =
         upper_mapped_coordinate[POLOIDAL_DIR] + n_poloidal_extend*dtheta;
      for (int i=mapping_block_size[POLOIDAL_DIR]-1; i>0; --i) {
         theta_pts[i-1] = theta_pts[i] - dtheta;
      }
   }
   else if (a_block_number == RCORE || a_block_number == RCSOL) {
      theta_pts[0] = lower_mapped_coordinate[POLOIDAL_DIR] - n_poloidal_extend*dtheta;
      for (int i=1; i<mapping_block_size[POLOIDAL_DIR]; ++i) {
         theta_pts[i] = theta_pts[i-1] + dtheta;
      }
   }
   else if (a_block_number == LPF || a_block_number == LSOL) {
      theta_pts[0] = lower_mapped_coordinate[POLOIDAL_DIR] - n_poloidal_extend*dtheta;
      for (int i=1; i<mapping_block_size[POLOIDAL_DIR]; ++i) {
         theta_pts[i] = theta_pts[i-1] + dtheta;
      }
   }
   else if (a_block_number == RPF || a_block_number == RSOL) {
      theta_pts[mapping_block_size[POLOIDAL_DIR]-1] =
         upper_mapped_coordinate[POLOIDAL_DIR] + n_poloidal_extend*dtheta;
      for (int i=mapping_block_size[POLOIDAL_DIR]-1; i>0; --i) {
         theta_pts[i-1] = theta_pts[i] - dtheta;
      }
   }
   else if (a_block_number == MCORE || m_block_number == MCSOL) {
      theta_pts[mapping_block_size[POLOIDAL_DIR]-1] =
         upper_mapped_coordinate[POLOIDAL_DIR] + n_poloidal_extend*dtheta;
      for (int i=mapping_block_size[POLOIDAL_DIR]-1; i>0; --i) {
         theta_pts[i-1] = theta_pts[i] - dtheta;
      }
   }
   else {
      MayDay::Error("SingleNullBlockCoordSys::init(): Invalid block_type encountered");
   }

   Box interp_box(IntVect::Zero,mapping_block_size-IntVect::Unit);

   FArrayBox interp_node_coords(interp_box, SpaceDim);

   // Create an interpolator for the cylindrical coordinates

   FArrayBox RZ_data(interp_box, 2);

   BoxIterator bit(interp_box);
   for (bit.begin();bit.ok();++bit) {
      IntVect iv = bit();
      int k = mapping_block_size[RADIAL_DIR]*iv[1] + iv[0];

      RZ_data(iv,0) = input_values[4*k];
      RZ_data(iv,1) = input_values[4*k+1];

      interp_node_coords(iv,0) = r_pts[iv[0] - interp_box.smallEnd(0)];
      interp_node_coords(iv,1) = theta_pts[iv[1] - interp_box.smallEnd(1)];
   }

   m_field_RZ_interpolator = new BSplineInterp(pp, 6, interp_node_coords, RZ_data);

   // Derive R times the magnetic field magnitude and the cylindrical components
   // of the field unit vector.  We need these quantities and their derivatives
   // at various points, so we create an interpolator for them.

   FArrayBox field_data(interp_box, NUM_FIELD_VARS);

   for (bit.begin();bit.ok();++bit) {
      IntVect iv = bit();
      int k = mapping_block_size[RADIAL_DIR]*iv[1] + iv[0];

      field_data(iv,BR_VAR) = input_values[4*k+2] / RZ_data(iv,0);
      field_data(iv,BZ_VAR) = input_values[4*k+3] / RZ_data(iv,0);
   }

   m_field_interpolator = new BSplineInterp(pp, 6, interp_node_coords, field_data);

   delete [] r_pts;
   delete [] theta_pts;
   delete [] input_values;
}



void FieldData::readFieldFromDCTFile(const string& a_file_name)
{
  ifstream fieldCoefficientsFile;

  int NR, NZ;

#ifdef CH_MPI
  if (procID() == 0) {
#endif

     fieldCoefficientsFile.open( a_file_name.c_str() );

     if (!fieldCoefficientsFile) {
        cout << "Unable to open field DCT coefficients file" << endl;
        exit(1);
     }

     //        cout << "Reading field coefficients from " << a_field_file_name.c_str() << endl;

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

#ifdef CH_MPI
  MPI_Bcast(&m_Rmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m_Rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m_Zmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m_Zmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&NR, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&NZ, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
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

#ifdef CH_MPI
  int num_values = kbox.numPts();
  MPI_Bcast(m_psi_coefs.dataPtr(), num_values, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  fieldCoefficientsFile.close();
}


void
FieldData::precomputeCoarseMappingData(const int n)
{
   // Precompute mapping data on a coarse grid for use in setting
   // the initial condition for the Newton solve used to compute
   // the inverse mapping.

   Box coarse_box(IntVect::Zero, IntVect(D_DECL6(n_radial_points - 1, n_poloidal_points - 1,0,0,0,0)));
   coarse_box.grow(RADIAL_DIR, n_radial_extend);
   coarse_box.grow(POLOIDAL_DIR, n_poloidal_extend);

   // N.B.: The input argument is intended to define a coarsening factor to be used
   // in constructing a coarser grid.  We are ignoring it for now, but can re-enable
   // it if performance becomes an issue.

   RealVect dx;
   dx[RADIAL_DIR] = 1. / (n_radial_points - 1);
   dx[POLOIDAL_DIR] = 1. / (n_poloidal_points - 1);

   m_physical_coords_coarse.define(coarse_box,SpaceDim);
   m_mapped_coords_coarse.define(coarse_box,SpaceDim);

   for (BoxIterator bit(coarse_box);bit.ok();++bit) {
      IntVect iv = bit();
      RealVect xi;
      for (int dir=0; dir<SpaceDim; ++dir) {
         xi[dir] = iv[dir]*dx[dir];
         m_mapped_coords_coarse(iv,dir) = xi[dir];
      }
      m_physical_coords_coarse(iv,0) = m_field_RZ_interpolator->interpolate(R_VAR, FUNCTION, xi);
      m_physical_coords_coarse(iv,1) = m_field_RZ_interpolator->interpolate(Z_VAR, FUNCTION, xi);
   }
}


RealVect
FieldData::getMappedCoordinate( const RealVect& a_X ) const
{
   /*
     Search for the coarse point that is closest to the passed argument
     and use the cooresponding mapped coordinate as the initial guess for the
     Newton solve.
   */

   RealVect xi_initial = RealVect::Zero;

   double max_dist = DBL_MAX;
   for (BoxIterator bit(m_physical_coords_coarse.box());bit.ok();++bit) {
      IntVect iv = bit();
      RealVect center_X;
      RealVect center_xi;
      for (int dir=0; dir<SpaceDim; ++dir) {
         center_X[dir] = m_physical_coords_coarse(iv,dir);
         center_xi[dir] = m_mapped_coords_coarse(iv,dir);
      }
      RealVect diff = center_X;
      diff -= a_X;
      double dist = diff.vectorLength();
      if (dist < max_dist) {
         xi_initial = center_xi;
         max_dist = dist;
      }
   }

   return mappedCoordNewton( a_X, xi_initial );
}



RealVect
FieldData::mappedCoordNewton( const RealVect& a_X,
                              const RealVect& a_xi_initial ) const
{
   /*
     Use Newton iteration to evaluate the mapping of physical to computational
     coordinates by inverting the mapping of computational to physical coordinates.
   */
   double tol = 1.e-6;         // Read from input if we find some reason to later
   int max_iter = 40;

   RealVect xi = a_xi_initial;

   RealVect residual;

   double R = m_field_RZ_interpolator->interpolate(R_VAR, FUNCTION, xi);
   double Z = m_field_RZ_interpolator->interpolate(Z_VAR, FUNCTION, xi);

   //   pout()<< "R = "<<R                      <<endl;
   //pout()<< "Z = "<<Z                      <<endl;
   //pout()<< "a_xi_initial = "<<a_xi_initial<<endl;

   double dRdr = m_field_RZ_interpolator->interpolate(R_VAR, RADIAL_DERIVATIVE, xi);
   double dRdtheta = m_field_RZ_interpolator->interpolate(R_VAR, POLOIDAL_DERIVATIVE, xi);
   double dZdr = m_field_RZ_interpolator->interpolate(Z_VAR, RADIAL_DERIVATIVE, xi);
   double dZdtheta = m_field_RZ_interpolator->interpolate(Z_VAR, POLOIDAL_DERIVATIVE, xi);

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

      R = m_field_RZ_interpolator->interpolate(R_VAR, FUNCTION, xi_proposed);
      Z = m_field_RZ_interpolator->interpolate(Z_VAR, FUNCTION, xi_proposed);
      dRdr = m_field_RZ_interpolator->interpolate(R_VAR, RADIAL_DERIVATIVE, xi_proposed);
      dRdtheta = m_field_RZ_interpolator->interpolate(R_VAR, POLOIDAL_DERIVATIVE, xi_proposed);
      dZdr = m_field_RZ_interpolator->interpolate(Z_VAR, RADIAL_DERIVATIVE, xi_proposed);
      dZdtheta = m_field_RZ_interpolator->interpolate(Z_VAR, POLOIDAL_DERIVATIVE, xi_proposed);

      residual[0] = R - a_X[0];
      residual[1] = Z - a_X[1];

      double Fnorm_proposed = residual.vectorLength();

#if 0
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

         R = m_field_RZ_interpolator->interpolate(R_VAR, FUNCTION, xi_proposed);
         Z = m_field_RZ_interpolator->interpolate(Z_VAR, FUNCTION, xi_proposed);
         dRdr = m_field_RZ_interpolator->interpolate(R_VAR, RADIAL_DERIVATIVE, xi_proposed);
         dRdtheta = m_field_RZ_interpolator->interpolate(R_VAR, POLOIDAL_DERIVATIVE, xi_proposed);
         dZdr = m_field_RZ_interpolator->interpolate(Z_VAR, RADIAL_DERIVATIVE, xi_proposed);
         dZdtheta = m_field_RZ_interpolator->interpolate(Z_VAR, POLOIDAL_DERIVATIVE, xi_proposed);

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

#if 1
   if ( !converged && num_iters >= max_iter ) {
      cout << "On block " << m_block_number << ": Newton solve did not converge at " << a_X << ", Fnorm = " << Fnorm << endl;
      for (int i=0; i<num_iters; ++i) {
         cout << i << " xi = " << xi_saved[i] << ", Fnorm = " << Fnorm_saved[i] << endl;
      }

      // iteration failed. Return initial guess.
      xi = a_xi_initial;
   }
#endif

   return xi;
}



void
FieldData::getPhysicalCoordinates(const FArrayBox& a_xi,
                                  FArrayBox&       a_physical_coordinates) const
{
   // Fills the components of a_physical_coordinates, evaluated at the mapped coordinates a_xi.
   const Box& box(a_xi.box());
   CH_assert(a_physical_coordinates.box().contains(box));
   CH_assert(a_physical_coordinates.nComp() == 2);

   m_RZ_interpolator->interpolate(R_VAR, 0, FUNCTION, a_xi, box, a_physical_coordinates);
   m_RZ_interpolator->interpolate(Z_VAR, 1, FUNCTION, a_xi, box, a_physical_coordinates);
}



void
FieldData::getFieldUnitVectorFromMappingFile(const FArrayBox& a_physical_coordinates,
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

     evaluated at the mapped coordinates a_xi.         
   */

   const Box& box(a_physical_coordinates.box());
   CH_assert(a_field_unit_vector.box().contains(box));
   CH_assert(a_field_unit_vector.nComp() == 6);

   // Get the corresponding mapped coordinates
   FArrayBox xi(box, SpaceDim);
   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();
      RealVect X;
      for (int dir=0; dir<SpaceDim; ++dir) {
         X[dir] = a_physical_coordinates(iv,dir);
      }
      RealVect xi_vect = getMappedCoordinate( X );
      for (int dir=0; dir<SpaceDim; ++dir) {
         xi(iv,dir) = xi_vect[dir];
      }
   }

   // Get the field components
   FArrayBox components(box,2);
   m_field_interpolator->interpolate(BR_VAR, 0, FUNCTION, xi, box, components);
   m_field_interpolator->interpolate(BZ_VAR, 1, FUNCTION, xi, box, components);

   // Compute the derivatives of the field components with respect to the mapped coordinates
   FArrayBox derivatives(box, 4);
   m_field_interpolator->interpolate(BR_VAR, 0, RADIAL_DERIVATIVE,   xi, box, derivatives);
   m_field_interpolator->interpolate(BR_VAR, 1, POLOIDAL_DERIVATIVE, xi, box, derivatives);
   m_field_interpolator->interpolate(BZ_VAR, 2, RADIAL_DERIVATIVE,   xi, box, derivatives);
   m_field_interpolator->interpolate(BZ_VAR, 3, POLOIDAL_DERIVATIVE, xi, box, derivatives);

   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();

      double Br = components(iv,0);
      double Bz = components(iv,1);
      double Bmag = sqrt( Br*Br + Bz*Bz );

      if (Bmag > 1.e-6) {

         double dBmag[2];
         for (int n=0; n<2; ++n) {
            a_field_unit_vector(iv,n) = components(iv,n) / Bmag;
            dBmag[n] = (Br*derivatives(iv,n) + Bz*derivatives(iv,n+2)) / Bmag;
         }

         a_field_unit_vector(iv,2) = (Bmag*derivatives(iv,0) - Br*dBmag[0]) / (Bmag*Bmag);
         a_field_unit_vector(iv,3) = (Bmag*derivatives(iv,1) - Br*dBmag[1]) / (Bmag*Bmag);

         a_field_unit_vector(iv,4) = (Bmag*derivatives(iv,2) - Bz*dBmag[0]) / (Bmag*Bmag);
         a_field_unit_vector(iv,5) = (Bmag*derivatives(iv,3) - Bz*dBmag[1]) / (Bmag*Bmag);

      }
      else {

         for (int comp=0; comp<6; ++comp) {
            a_field_unit_vector(iv,comp) = 0.;
         }
      }
   }
}

void
FieldData::convertToMappedDerivatives( const FArrayBox& a_mapped_coordinates,
                                       FArrayBox&       a_field_unit_vector ) const
{
   /*
     Given the FArrayBox a_mapped_coordinates of mapped coordinates, convert
     the six-component FArrayBox containing

         Component 0: First component of b
         Component 1: Second component of b
         Component 2: first derivative wrt X[0] of the first component of b
         Component 3: first derivative wrt X[1] of the first component of b
         Component 4: first derivative wrt X[0] of the second component of b
         Component 5: first derivative wrt X[1] of the second component of b,

     where X denotes the physical coordinate, to the six-component FArrayBox containing

         Component 0: First component of b
         Component 1: Second component of b
         Component 2: first derivative wrt xi[0] of the first component of b
         Component 3: first derivative wrt xi[1] of the first component of b
         Component 4: first derivative wrt xi[0] of the second component of b
         Component 5: first derivative wrt xi[1] of the second component of b,

     where xi denotes the mapped coordinates.  In other words, the first two
     components containing the field unit vectors are unchanged, and the remaining
     four are converted from physical to mapped coordinate derivatives using
     the chain rule.
   */


   CH_assert(a_mapped_coordinates.box() == a_field_unit_vector.box());

   for (BoxIterator bit(a_field_unit_vector.box()); bit.ok(); ++bit) {
      IntVect iv = bit();

      IntVect mapped_coord;
      for (int dir=0; dir<SpaceDim; ++dir) {
         mapped_coord[dir] = a_mapped_coordinates(iv,dir);
      }

      // Get the derivatives of the physical coordinates wrt the mapped coordinates
      double dRdxi0 = m_RZ_interpolator->interpolate(R_VAR, RADIAL_DERIVATIVE, mapped_coord);
      double dRdxi1 = m_RZ_interpolator->interpolate(R_VAR, POLOIDAL_DERIVATIVE, mapped_coord);
      double dZdxi0 = m_RZ_interpolator->interpolate(Z_VAR, RADIAL_DERIVATIVE, mapped_coord);
      double dZdxi1 = m_RZ_interpolator->interpolate(Z_VAR, POLOIDAL_DERIVATIVE, mapped_coord);

      // Get the derivatives of the field unit vector components wrt the physical coordinates
      double dbRdR = a_field_unit_vector(iv,2);
      double dbRdZ = a_field_unit_vector(iv,3);
      double dbZdR = a_field_unit_vector(iv,4);
      double dbZdZ = a_field_unit_vector(iv,5);

      // Replace the derivatives of the field unit vector components wrt the physical coordinate
      // by the derivatives wrt mapped coordinates
      a_field_unit_vector(iv,2) = dbRdR * dRdxi0 + dbRdZ * dZdxi0;
      a_field_unit_vector(iv,3) = dbRdR * dRdxi1 + dbRdZ * dZdxi1;
      a_field_unit_vector(iv,4) = dbZdR * dRdxi0 + dbZdZ * dZdxi0;
      a_field_unit_vector(iv,5) = dbZdR * dRdxi1 + dbZdZ * dZdxi1;
   }
}



void
FieldData::getFieldUnitVectorFromDCT( const FArrayBox& a_physical_coordinates,
                                      FArrayBox&       a_field_unit_vector ) const
{
   /*
     Given the FArrayBox a_physical_coordinates of physical coordinates, fill the
     components of a_field_unit_vector as follows:

         Component 0: First component of b
         Component 1: Second component of b
         Component 2: first derivative wrt X[0] of the first component of b
         Component 3: first derivative wrt X[1] of the first component of b
         Component 4: first derivative wrt X[0] of the second component of b
         Component 5: first derivative wrt X[1] of the second component of b

     where b is the field unit vector and X denotes the physical coordinate.
   */

   const Box& box(a_physical_coordinates.box());
   CH_assert(a_field_unit_vector.box().contains(box));
   CH_assert(a_field_unit_vector.nComp() >= 6);

   double RB[2];
   double dRBdR[2];
   double dRBdZ[2];

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

      RB[1] = value;  // R*B_Z

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

      RB[0] = -value;  // R*B_R

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

      dRBdR[1] = value;  // d(R*BZ)/dR

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

      dRBdZ[0] = -value;  // d(R*B_R)/dZ

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

      dRBdR[0] = -value;  // d(R*B_R)/dR
      dRBdZ[1] =  value;  // d(R*B_Z)/dZ

      double R = a_physical_coordinates(iv,0);
      double Br = RB[0] / R;
      double Bz = RB[1] / R;
      double Bmag = sqrt( Br*Br + Bz*Bz );
      double RBmag = Bmag * R;

      if (Bmag > 1.e-6) {

         double b[2];
         for (int n=0; n<2; ++n) {
            b[n] = RB[n] / RBmag;
         }

         double dBmag[2];
         dBmag[0] = (b[0] * dRBdR[0] + b[1] * dRBdR[1] - Bmag) / R;
         dBmag[1] = (b[0] * dRBdZ[0] + b[1] * dRBdZ[1]       ) / R;

         a_field_unit_vector(iv,0) = b[0];
         a_field_unit_vector(iv,1) = b[1];

         a_field_unit_vector(iv,2) = (dRBdR[0] - R * b[0] * dBmag[0]  - Br) / RBmag;
         a_field_unit_vector(iv,3) = (dRBdZ[0] - R * b[0] * dBmag[1]) / RBmag;
         a_field_unit_vector(iv,4) = (dRBdR[1] - R * b[1] * dBmag[0]  - Bz) / RBmag;
         a_field_unit_vector(iv,5) = (dRBdZ[1] - R * b[1] * dBmag[1]) / RBmag;

      }
      else {

         for (int comp=0; comp<6; ++comp) {
            a_field_unit_vector(iv,comp) = 0.;
         }
      }
   }

   delete[] temp;
}



RealVect
FieldData::getFieldUnitVectorFromDCT( const RealVect& a_physical_coordinate ) const
{
   /*
     Given the physical coordinate a_physical_coordinate, return field unit vector.  If the
     field is zero, the zero vector is returned;
   */

   double RB[2];
   double dRBdR[2];
   double dRBdZ[2];

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

   RB[1] = value;  // R*B_Z

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

   RB[0] = -value;  // R*B_R

   double R = a_physical_coordinate[0];
   double Br = RB[0] / R;
   double Bz = RB[1] / R;
   double Bmag = sqrt( Br*Br + Bz*Bz );
   double RBmag = Bmag * R;

   RealVect field_unit_vector;

   if (Bmag > 1.e-6) {
      for (int n=0; n<2; ++n) {
         field_unit_vector[n] = RB[n] / RBmag;
      }
   }
   else {
      field_unit_vector = RealVect::Zero;
   }

   delete[] temp;

   return field_unit_vector;
}



void
FieldData::getMagneticFluxFromDCT( const FArrayBox& a_physical_coordinates,
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
FieldData::getMagneticFluxFromDCT( const RealVect& a_physical_coordinate ) const
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



void
FieldData::writePhysicalCoordinates(const FArrayBox& a_physical_coordinates) const
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
FieldData::writeVectors(const FArrayBox& a_physical_coordinates,
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


void
FieldData::writeMagneticFlux(const FArrayBox& a_physical_coordinates,
                             const FArrayBox& a_magnetic_flux) const
{
   char file_name[80];
   sprintf(file_name, "psi%d", m_block_number);
   FILE* fd = fopen(file_name, "w");

   const Box& box = a_physical_coordinates.box();

   fprintf(fd, "%d %d 0\n", box.size(0), box.size(1));

   for (BoxIterator bit(box);bit.ok();++bit) {
      IntVect iv = bit();
      fprintf(fd, "%20.12e %20.12e %20.12e \n", a_physical_coordinates(iv,0), a_physical_coordinates(iv,1), a_magnetic_flux(iv,0));
   }

   fclose(fd);
}



#include "NamespaceFooter.H"


