#include <float.h>
#include <fstream>
#include <array>
#include "CONSTANTS.H"
#include "SingleNullBlockCoordSys.H"
#include "Directions.H"

#include "inspect.H"

#undef REPORT_NEWTON_FAILURE

#include "NamespaceHeader.H"

const std::string SingleNullBlockCoordSys::pp_name = "singlenullRealisticGeom";


SingleNullBlockCoordSys::SingleNullBlockCoordSys( ParmParse&            a_parm_parse,
                                                  const ProblemDomain&  a_domain,
                                                  const RealVect&       a_dx,
                                                  const int             a_poloidal_block,
                                                  const int             a_toroidal_sector,
                                                  const int             a_toroidal_block_separation )
    : MagBlockCoordSys(a_parm_parse),
      m_poloidal_block(a_poloidal_block),
      m_toroidal_sector(a_toroidal_sector),
      m_poloidal_util(NULL),
      m_field_aligned_mapping(false),
      m_poloidally_truncated(false),
      m_include_extrablock_mapping(true),
      m_phys_coord_type(CARTESIAN),
      m_toroidal_ghosts(4),
      m_toroidal_mapping_refinement(2),
      m_field_trace_step_num(10)
 {
 #ifdef PLOT_INVERSE_POINTS
    char file_name[80];
    sprintf(file_name, "Block%d_%d", a_poloidal_block, procID());
    m_ipt_file[a_poloidal_block] = fopen(file_name, "w");
#endif

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

    a_parm_parse.query("include_extrablock_mapping", m_include_extrablock_mapping);

    // If a field coefficients file was specified in the input, use it
    // to define the field.  Otherwise, we read it from the mapping file.
    m_spectral_field = a_parm_parse.contains("field_coefficients_file");

    m_poloidal_util = new POL::SingleNullPoloidalUtil(a_parm_parse, m_spectral_field, m_RB_toroidal);

    define( a_domain, a_dx );

#if CFG_DIM==3
    m_toroidal_mapped_offset = a_toroidal_sector * a_toroidal_block_separation * a_dx[TOROIDAL_DIR];
#endif
 }


SingleNullBlockCoordSys::~SingleNullBlockCoordSys()
{
#ifdef PLOT_INVERSE_POINTS
   for (int block=0; block<NUM_POLOIDAL_BLOCKS; ++block) {
      fclose(m_ipt_file[block]);
   }
#endif
   if (m_poloidal_util) delete m_poloidal_util;
}


string
SingleNullBlockCoordSys::blockName( int a_poloidal_block ) const
{
   CH_assert(a_poloidal_block < NUM_POLOIDAL_BLOCKS);

   string name;

   switch ( a_poloidal_block )
      {
      case LCORE:
         name = "left core";
         break;
      case MCORE:
         name = "middle core";
         break;
      case RCORE:
         name = "right core";
         break;
      case LCSOL:
         name = "left center scrape-off layer";
         break;
      case MCSOL:
         name = "middle center scrape-off layer";
         break;
      case RCSOL:
         name = "right center scrape-off layer";
         break;
      case LSOL:
         name = "left scrape_off layer";
         break;
      case RSOL:
         name = "right scrape-off layer";
         break;
      case LPF:
         name = "left private flux";
         break;
      case RPF:
         name = "right private flux";
         break;
      default:
         MayDay::Error("SingleNullBlockCoordSys::blockName(): Invalid block_type encountered");
      }

   return name;
}


void
SingleNullBlockCoordSys::printInit() const
{
   if ( m_verbose && procID()==0 ) {
      cout << "Constructing single null " << blockName(m_poloidal_block) << " block";
#if CFG_DIM==3
      cout << " in toroidal sector " << m_toroidal_sector;
#endif
      cout << " with global index space domain box = " << m_domain.domainBox() << endl;

      switch ( m_poloidal_block )
         {
         case LCORE:
            cout << "Left core mapped domain: "
                 << lowerMappedCoordinate(RADIAL_DIR) << " < xi_" << RADIAL_DIR << " < " << upperMappedCoordinate(RADIAL_DIR) << ", " 
#if CFG_DIM==3
                 << lowerMappedCoordinate(TOROIDAL_DIR) << " < xi_" << TOROIDAL_DIR << " < " << upperMappedCoordinate(TOROIDAL_DIR) << ", " 
#endif
                 << lowerMappedCoordinate(POLOIDAL_DIR) << " < xi_" << POLOIDAL_DIR << " < " << upperMappedCoordinate(POLOIDAL_DIR);
            break;
         case MCORE:
            cout << "Middle core mapped domain: " 
                 << lowerMappedCoordinate(RADIAL_DIR) << " < xi_" << RADIAL_DIR << " < " << upperMappedCoordinate(RADIAL_DIR) << ", " 
#if CFG_DIM==3
                 << lowerMappedCoordinate(TOROIDAL_DIR) << " < xi_" << TOROIDAL_DIR << " < " << upperMappedCoordinate(TOROIDAL_DIR) << ", " 
#endif
                 << lowerMappedCoordinate(POLOIDAL_DIR) << " < xi_" << POLOIDAL_DIR << " < " << upperMappedCoordinate(POLOIDAL_DIR);
            break;
         case RCORE:
            cout << "Right core mapped domain: " 
                 << lowerMappedCoordinate(RADIAL_DIR) << " < xi_" << RADIAL_DIR << " < " << upperMappedCoordinate(RADIAL_DIR) << ", " 
#if CFG_DIM==3
                 << lowerMappedCoordinate(TOROIDAL_DIR) << " < xi_" << TOROIDAL_DIR << " < " << upperMappedCoordinate(TOROIDAL_DIR) << ", " 
#endif
                 << lowerMappedCoordinate(POLOIDAL_DIR) << " < xi_" << POLOIDAL_DIR << " < " << upperMappedCoordinate(POLOIDAL_DIR);
            break;
         case LCSOL:
            cout << "Left center scrape-off layer mapped domain: " 
                 << lowerMappedCoordinate(RADIAL_DIR) << " < xi_" << RADIAL_DIR << " < " << upperMappedCoordinate(RADIAL_DIR) << ", " 
#if CFG_DIM==3
                 << lowerMappedCoordinate(TOROIDAL_DIR) << " < xi_" << TOROIDAL_DIR << " < " << upperMappedCoordinate(TOROIDAL_DIR) << ", " 
 #endif
                 << lowerMappedCoordinate(POLOIDAL_DIR) << " < xi_" << POLOIDAL_DIR << " < " << upperMappedCoordinate(POLOIDAL_DIR);
            break;
         case MCSOL:
             cout << "Middle center scrape-off layer mapped domain: " 
                  << lowerMappedCoordinate(RADIAL_DIR) << " < xi_" << RADIAL_DIR << " < " << upperMappedCoordinate(RADIAL_DIR) << ", " 
#if CFG_DIM==3
                  << lowerMappedCoordinate(TOROIDAL_DIR) << " < xi_" << TOROIDAL_DIR << " < " << upperMappedCoordinate(TOROIDAL_DIR) << ", " 
#endif
                  << lowerMappedCoordinate(POLOIDAL_DIR) << " < xi_" << POLOIDAL_DIR << " < " << upperMappedCoordinate(POLOIDAL_DIR);
             break;
         case RCSOL:
            cout << "Right center scrape-off layer mapped domain: " 
                 << lowerMappedCoordinate(RADIAL_DIR) << " < xi_" << RADIAL_DIR << " < " << upperMappedCoordinate(RADIAL_DIR) << ", " 
#if CFG_DIM==3
                 << lowerMappedCoordinate(TOROIDAL_DIR) << " < xi_" << TOROIDAL_DIR << " < " << upperMappedCoordinate(TOROIDAL_DIR) << ", " 
#endif
                 << lowerMappedCoordinate(POLOIDAL_DIR) << " < xi_" << POLOIDAL_DIR << " < " << upperMappedCoordinate(POLOIDAL_DIR);
            break;
         case LSOL:
            cout << "Left scrape_off layer mapped domain: " 
                 << lowerMappedCoordinate(RADIAL_DIR) << " < xi_" << RADIAL_DIR << " < " << upperMappedCoordinate(RADIAL_DIR) << ", " 
#if CFG_DIM==3
                 << lowerMappedCoordinate(TOROIDAL_DIR) << " < xi_" << TOROIDAL_DIR << " < " << upperMappedCoordinate(TOROIDAL_DIR) << ", " 
#endif
                 << lowerMappedCoordinate(POLOIDAL_DIR) << " < xi_" << POLOIDAL_DIR << " < " << upperMappedCoordinate(POLOIDAL_DIR);
            break;
         case RSOL:
            cout << "Right scrape-off layer mapped domain: " 
                 << lowerMappedCoordinate(RADIAL_DIR) << " < xi_" << RADIAL_DIR << " < " << upperMappedCoordinate(RADIAL_DIR) << ", " 
#if CFG_DIM==3
                 << lowerMappedCoordinate(TOROIDAL_DIR) << " < xi_" << TOROIDAL_DIR << " < " << upperMappedCoordinate(TOROIDAL_DIR) << ", " 
#endif
                 << lowerMappedCoordinate(POLOIDAL_DIR) << " < xi_" << POLOIDAL_DIR << " < " << upperMappedCoordinate(POLOIDAL_DIR);
            break;
         case LPF:
            cout << "Left private flux mapped domain: " 
                 << lowerMappedCoordinate(RADIAL_DIR) << " < xi_" << RADIAL_DIR << " < " << upperMappedCoordinate(RADIAL_DIR) << ", " 
#if CFG_DIM==3
                 << lowerMappedCoordinate(TOROIDAL_DIR) << " < xi_" << TOROIDAL_DIR << " < " << upperMappedCoordinate(TOROIDAL_DIR) << ", " 
#endif
                 << lowerMappedCoordinate(POLOIDAL_DIR) << " < xi_" << POLOIDAL_DIR << " < " << upperMappedCoordinate(POLOIDAL_DIR);
            break;
         case RPF:
            cout << "Right private flux mapped domain " 
                 << lowerMappedCoordinate(RADIAL_DIR) << " < xi_" << RADIAL_DIR << " < " << upperMappedCoordinate(RADIAL_DIR) << ", " 
#if CFG_DIM==3
                 << lowerMappedCoordinate(TOROIDAL_DIR) << " < xi_" << TOROIDAL_DIR << " < " << upperMappedCoordinate(TOROIDAL_DIR) << ", " 
#endif
                 << lowerMappedCoordinate(POLOIDAL_DIR) << " < xi_" << POLOIDAL_DIR << " < " << upperMappedCoordinate(POLOIDAL_DIR);
            break;
         default:
            MayDay::Error("SingleNullBlockCoordSys::SingleNullBlockCoordSys(): Invalid block_type encountered");
         }
      cout << endl;
   }
}

void
SingleNullBlockCoordSys::getCellCenterRealCoords()
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


bool
SingleNullBlockCoordSys::blockNameIsValid( const string& a_block_name ) const
{
   return (a_block_name == "mcore" && m_poloidal_block == MCORE) ||
      (a_block_name == "lcore" && m_poloidal_block == LCORE) ||
      (a_block_name == "rcore" && m_poloidal_block == RCORE) ||
      (a_block_name == "mcsol" && m_poloidal_block == MCSOL) ||
      (a_block_name == "lcsol" && m_poloidal_block == LCSOL) ||
      (a_block_name == "rcsol" && m_poloidal_block == RCSOL) ||
      (a_block_name == "lpf"   && m_poloidal_block == LPF) ||
      (a_block_name == "rpf"   && m_poloidal_block == RPF) ||
      (a_block_name == "lsol"  && m_poloidal_block == LSOL) ||
      (a_block_name == "rsol"  && m_poloidal_block == RSOL);
}


void
SingleNullBlockCoordSys::readFiles( ParmParse& a_pp )
{
   printInit();

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
         if ( blockNameIsValid(block_name) ) {
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

   POL::IntVect mapping_block_size;
   mapping_block_size[0] = n_radial_points   + 2*n_radial_extend;
   mapping_block_size[1] = n_poloidal_points + 2*n_poloidal_extend;
   num_block_elements = mapping_block_size[0] * mapping_block_size[1];
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

   double dr = ( upperMappedCoordinate(RADIAL_DIR)
                 - lowerMappedCoordinate(RADIAL_DIR) ) / num_radial_cells;
   
   double *r_pts = new double[mapping_block_size[0]];

   r_pts[0] = lowerMappedCoordinate(RADIAL_DIR) - n_radial_extend*dr;
   for (int i=1; i<mapping_block_size[0]; ++i) {
      r_pts[i] = r_pts[i-1] + dr;
   }

   int num_poloidal_cells = n_poloidal_points - 1;

   double dtheta = ( upperMappedCoordinate(POLOIDAL_DIR)
                     - lowerMappedCoordinate(POLOIDAL_DIR) ) / num_poloidal_cells;

   double *theta_pts = new double[mapping_block_size[1]];

   int block_full_poloidal = num_poloidal_cells;
   int block_poloidal = block_full_poloidal;

   definePoints(a_pp, block_poloidal, block_full_poloidal, mapping_block_size, n_poloidal_extend, dtheta, theta_pts); 

   m_poloidally_truncated = (block_poloidal < block_full_poloidal);

   // Determine interpolation method
   string interp_method = "bspline";
   if (a_pp.contains("interpolation_method")) {
      a_pp.get("interpolation_method", interp_method);
   }

#if CFG_DIM==3
   // Since the mapping is only defined using bspline interpolation for now (and maybe forever),
   // we need to restrict the poloidal interpolation to bspline also for consistency
   if ( interp_method != "bspline" ) {
      MayDay::Error("SingleNullBlockCoordSys::readFiles(): interp_method is limited to bspline in 3D for now");
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

      RZ_data(iv,0) = input_values[4*k];
      RZ_data(iv,1) = input_values[4*k+1];

      interp_node_coords(iv,0) = r_pts[iv[0] - interp_box.smallEnd(0)];
      interp_node_coords(iv,1) = theta_pts[iv[1] - interp_box.smallEnd(1)];
   }

   if (!m_include_extrablock_mapping) {

      POL::FArrayBox interp_node_coords_noGhost;
      POL::FArrayBox RZ_data_noGhost;
      removeExtraBlockGhosts(interp_node_coords_noGhost,
                             RZ_data_noGhost,
                             interp_node_coords,
                             RZ_data,
                             n_radial_extend,
                             n_poloidal_extend);
      
      m_poloidal_util->setRZInterp(a_pp,
                                   interp_method,
                                   interp_node_coords_noGhost,
                                   RZ_data_noGhost);
#if CFG_DIM==3
      setInterp(a_pp,
                interp_node_coords_noGhost,
                RZ_data_noGhost);
#endif
   }
   
   else {

     m_poloidal_util->setRZInterp(a_pp, interp_method, interp_node_coords, RZ_data);
#if CFG_DIM==3
     setInterp(a_pp, interp_node_coords, RZ_data);
#endif
   }
   
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

         field_data(iv,POL::SingleNullPoloidalUtil::RBR_VAR) = RBr = input_values[4*k+2];
         field_data(iv,POL::SingleNullPoloidalUtil::RBZ_VAR) = RBz = input_values[4*k+3];
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
SingleNullBlockCoordSys::definePoints( const ParmParse&     a_pp,
                                       int&                 a_block_poloidal,
                                       const int            a_block_full_poloidal,
                                       const POL::IntVect&  a_mapping_block_size,
                                       const int&           a_n_poloidal_extend,
                                       double &             a_dtheta,
                                       double *             a_theta_pts ) const
{
   if ( m_poloidal_block == LCORE || m_poloidal_block == LCSOL ) {
      if ( a_pp.contains("core_mapping_truncate") ) {
         a_pp.get("core_mapping_truncate", a_block_poloidal);
         CH_assert(a_block_poloidal <= a_mapping_block_size[1]);
      }

      a_dtheta *= (double)a_block_full_poloidal / (double)a_block_poloidal;

      a_theta_pts[a_mapping_block_size[1]-1] =
         upperMappedCoordinate(POLOIDAL_DIR) + a_n_poloidal_extend*a_dtheta;
      for (int i=a_mapping_block_size[1]-1; i>0; --i) {
         a_theta_pts[i-1] = a_theta_pts[i] - a_dtheta;
      }
   }
   else if ( m_poloidal_block == RCORE || m_poloidal_block == RCSOL ) {
      if ( a_pp.contains("core_mapping_truncate") ) {
         a_pp.get("core_mapping_truncate", a_block_poloidal);
         CH_assert(a_block_poloidal <= a_mapping_block_size[1]);
      }

      a_dtheta *= (double)a_block_full_poloidal / (double)a_block_poloidal;

      a_theta_pts[0] = lowerMappedCoordinate(POLOIDAL_DIR) - a_n_poloidal_extend*a_dtheta;
      for (int i=1; i<a_mapping_block_size[1]; ++i) {
         a_theta_pts[i] = a_theta_pts[i-1] + a_dtheta;
      }
   }
   else if ( m_poloidal_block == LPF || m_poloidal_block == LSOL ) {
      if ( a_pp.contains("pf_mapping_truncate") ) {
         a_pp.get("pf_mapping_truncate", a_block_poloidal);
         CH_assert(a_block_poloidal <= a_mapping_block_size[1]);
      }

      a_dtheta *= (double)a_block_full_poloidal / (double)a_block_poloidal;

      a_theta_pts[0] = lowerMappedCoordinate(POLOIDAL_DIR) - a_n_poloidal_extend*a_dtheta;
      for (int i=1; i<a_mapping_block_size[1]; ++i) {
         a_theta_pts[i] = a_theta_pts[i-1] + a_dtheta;
      }
   }
   else if ( m_poloidal_block == RPF || m_poloidal_block == RSOL ) {
      if ( a_pp.contains("pf_mapping_truncate") ) {
         a_pp.get("pf_mapping_truncate", a_block_poloidal);
         CH_assert(a_block_poloidal <= a_mapping_block_size[1]);
      }

      a_dtheta *= (double)a_block_full_poloidal / (double)a_block_poloidal;

      a_theta_pts[a_mapping_block_size[1]-1] =
         upperMappedCoordinate(POLOIDAL_DIR) + a_n_poloidal_extend*a_dtheta;
      for (int i=a_mapping_block_size[1]-1; i>0; --i) {
         a_theta_pts[i-1] = a_theta_pts[i] - a_dtheta;
      }
   }
   else if ( m_poloidal_block == MCORE || m_poloidal_block == MCSOL ) {
      a_dtheta *= (double)a_block_full_poloidal / (double)a_block_poloidal;

      a_theta_pts[a_mapping_block_size[1]-1] =
         upperMappedCoordinate(POLOIDAL_DIR) + a_n_poloidal_extend*a_dtheta;
      for (int i=a_mapping_block_size[1]-1; i>0; --i) {
         a_theta_pts[i-1] = a_theta_pts[i] - a_dtheta;
      }
   }
   else {
      MayDay::Error("SingleNullBlockCoordSys::init(): Invalid block_type encountered");
   }
}

void
SingleNullBlockCoordSys::removeExtraBlockGhosts( POL::FArrayBox& a_coords,
                                                 POL::FArrayBox& a_data,
                                           const POL::FArrayBox& a_coords_full,
                                           const POL::FArrayBox& a_data_full,
                                           const int             a_radial_extent,
                                           const int             a_poloidal_extent  ) const
{

   //Strip extrablock ghost cells (but not physical ghosts) from the mapping data
   
   const POL::Box& box_full = a_data_full.box();
   const POL::IntVect box_full_hi = box_full.bigEnd();
   const POL::IntVect box_full_lo = box_full.smallEnd();
   
   POL::Box inter_box(box_full);
   
   if ( m_poloidal_block == LCORE || m_poloidal_block == RCORE || m_poloidal_block == MCORE) {
      inter_box.setBig(0, box_full_hi[0] - a_radial_extent);
      inter_box.setBig(1, box_full_hi[1] - a_poloidal_extent);
      inter_box.setSmall(1, box_full_lo[1] + a_poloidal_extent);
   }
   
   if ( m_poloidal_block == LCSOL || m_poloidal_block == RCSOL || m_poloidal_block == MCSOL) {
      inter_box.setSmall(0, box_full_lo[0] + a_radial_extent);
      inter_box.setBig(1, box_full_hi[1] - a_poloidal_extent);
      inter_box.setSmall(1, box_full_lo[1] + a_poloidal_extent);
   }

   if ( m_poloidal_block == LSOL) {
      inter_box.setSmall(0, box_full_lo[0] + a_radial_extent);
      inter_box.setSmall(1, box_full_lo[1] + a_poloidal_extent);
   }

   if ( m_poloidal_block == RSOL) {
      inter_box.setSmall(0, box_full_lo[0] + a_radial_extent);
      inter_box.setBig(1, box_full_hi[1] - a_poloidal_extent);
   }

   if ( m_poloidal_block == LPF) {
      inter_box.setBig(0, box_full_hi[0] - a_radial_extent);
      inter_box.setSmall(1, box_full_lo[1] + a_poloidal_extent);
   }

   if ( m_poloidal_block == RPF) {
      inter_box.setBig(0, box_full_hi[0] - a_radial_extent);
      inter_box.setBig(1, box_full_hi[1] - a_poloidal_extent);
   }

   a_coords.define(inter_box,POL_DIM);
   a_coords.copy(a_coords_full);
   
   a_data.define(inter_box, POL_DIM);
   a_data.copy(a_data_full);
   
}

RealVect
SingleNullBlockCoordSys::realCoord( const RealVect& a_Xi ) const
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
       MayDay::Error("SingleNullBlockCoordSys::realCoord(): unrecognized m_phys_coord_type");
    }
#endif

    return X;
}



RealVect
SingleNullBlockCoordSys::mappedCoord( const RealVect& a_X ) const
{
#ifdef PLOT_INVERSE_POINTS
   fprintf(m_ipt_file[m_poloidal_block], "%20.12e %20.12e\n", a_X[0], a_X[1]);
#endif

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
      MayDay::Error("SingleNullBlockCoordSys::mappedCoord(): unrecognized m_phys_coord_type");
   }
   xi_initial[1] = Phi + m_toroidal_mapped_offset;

   RealVect xi = mappedCoordNewton(RealVect(R,Phi,Z), xi_initial, iv_closest);

   double Phi_lo = lowerMappedCoordinate(TOROIDAL_DIR) - m_toroidal_ghosts*m_dx[TOROIDAL_DIR];
   double Phi_hi = upperMappedCoordinate(TOROIDAL_DIR) + m_toroidal_ghosts*m_dx[TOROIDAL_DIR];

   if ( m_phys_coord_type == CARTESIAN ) {
#if 0
      // The toroidal angle above lies between -pi and pi; shift it into the mapped variable range
      if ( xi[1] < Phi_lo ) {
         xi[1] += 2.*Pi;
         CH_assert(xi[1] >= Phi_lo && xi[1] <= Phi_hi);
      }
      else if ( xi[1] > Phi_hi ) {
         xi[1] -= 2.*Pi;
         CH_assert(xi[1] >= Phi_lo && xi[1] <= Phi_hi);
      }
#else
      //      if ( xi[1] < Phi_lo || xi[1] > Phi_hi ) {
      //         cout << m_poloidal_block << " " << m_toroidal_sector << " " << xi[1] << " " 
      //       <<lowerMappedCoordinate(TOROIDAL_DIR) << " " << upperMappedCoordinate(TOROIDAL_DIR) << endl;
      //      }
#endif
   }
   
   return xi;

#endif
}


Real
SingleNullBlockCoordSys::dXdXi( const RealVect&  a_Xi,
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
         MayDay::Error("SingleNullBlockCoordSys::dXdXi(): unrecognized m_phys_coord_type");
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
SingleNullBlockCoordSys::dXdXi( FArrayBox&       a_dXdXi,
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
SingleNullBlockCoordSys::isValid( const RealVect& a_xi,
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
SingleNullBlockCoordSys::getMagneticFlux( const FArrayBox&  a_physical_coordinates,
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
SingleNullBlockCoordSys::getMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   double psi;

   if (m_spectral_field) {
      psi = m_poloidal_util->getMagneticFluxFromDCT(restrictPhysCoordToPoloidal(a_physical_coordinate));
   }
   else {
      MayDay::Error("SingleNullBlockCoordSys::getMagneticFlux() only implemented for DCT field option");
      psi = 0.;
   }

   return psi;
}


void
SingleNullBlockCoordSys::getNormMagneticFlux( const FArrayBox& a_physical_coordinates,
                                              FArrayBox&       a_magnetic_flux ) const

{
   double physFluxOnAxis = getMagneticFlux(injectToroidal(m_poloidal_util->getMagAxis()));
   double physFluxOnSep = getMagneticFlux(injectToroidal(m_poloidal_util->getXpoint()));
   
   getMagneticFlux(a_physical_coordinates, a_magnetic_flux);
   a_magnetic_flux.plus(-physFluxOnAxis);
   a_magnetic_flux.divide(physFluxOnSep - physFluxOnAxis);
   
}


double
SingleNullBlockCoordSys::getNormMagneticFlux( const RealVect& a_physical_coordinate ) const
{
   double physFluxOnAxis = getMagneticFlux(injectToroidal(m_poloidal_util->getMagAxis()));
   double physFluxOnSep = getMagneticFlux(injectToroidal(m_poloidal_util->getXpoint()));
   double physFlux = getMagneticFlux(a_physical_coordinate);
   
   return (physFlux - physFluxOnAxis)/(physFluxOnSep - physFluxOnAxis);
}


double
SingleNullBlockCoordSys::getOuterMidplaneCoord( const double a_psiNorm ) const
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
   else MayDay::Error("SingleNullBlockCoordSys::getOuterMidplaneCoord only implemented for DCT field option");
   
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
SingleNullBlockCoordSys::getOuterRsep() const
{
   double psiNorm_atSep = getNormMagneticFlux(injectToroidal(m_poloidal_util->getXpoint()));
   
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
      MayDay::Error("SingleNullBlockCoordSys::computeFieldData(): Only implemented in 3D for spectral field option");
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
SingleNullBlockCoordSys::computeBField(const RealVect& a_X) const
{
   if ( !m_spectral_field ) {
      MayDay::Error("SingleNullBlockCoordSys::computeBField(): Only implemented for DCT field option");
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
SingleNullBlockCoordSys::getNodalFieldData(FArrayBox& a_points,
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
      MayDay::Error("SingleNullBlockCoordSys::getNodalFieldData(): Only implemented for DCT field option");
   }
}


#if CFG_DIM==3
void
SingleNullBlockCoordSys::setInterp( const ParmParse&       a_pp,
                                    const POL::FArrayBox&  a_interp_node_coords_pol,
                                    const POL::FArrayBox&  a_RZ_data )
{
   // Make an index space for the interpolation points by expanding the poloidal index
   // space in the toroidal direction.  In the new toroidal direction, the middle n_toroidal
   // points span the valid block.  There are an additional m_toroidal_ghosts points in the
   // lower and upper toroidal directions to provide a ghost cell region.

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
         MayDay::Error("SingleNullBlockCoordSys::setInterp(): No 3D mapping data file specified");
      }

      fstream file_stream;

      if ( !read_3D_mapping_data ) {
      
         if ( procID()==0 ) {
            cout << "Constructing 3D mapping data for " << blockName(m_poloidal_block)
                 << " in toroidal sector " << m_toroidal_sector << endl;
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
                     MayDay::Error("SingleNullBlockCoordSys::setInterp(): interpolation data for block is wrong size");
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
               MayDay::Error("SingleNullBlockCoordSys::setInterp(): interpolation data for block was not found");
            }
            
            cout << "Reading interpolation data from file " << file_name << " for "
                 << blockName(m_poloidal_block) << " in toroidal sector " << m_toroidal_sector << endl << endl;

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
SingleNullBlockCoordSys::traceField( const RealVect&  a_X,
                                     const double&    a_toroidal_increment ) const
{
   
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
SingleNullBlockCoordSys::getPoloidalDisjointBoxLayout(DisjointBoxLayout&   a_grids,
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
   
   // Too few processors available (the approach below is a strating point; not
   // very effective as we may end up not using about half of available proc)
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
      }
   }

   Vector<int> proc_ids;

   IntVect domain_hi_end = a_box.bigEnd();
   IntVect domain_lo_end = a_box.smallEnd();

   int n_rad_dec = n_rad / n_loc_rad;
   int n_pol_dec = n_pol / n_loc_pol;
   
   a_decomp_num = n_rad_dec * n_pol_dec;
   
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
SingleNullBlockCoordSys::assembleDecomposedData(FArrayBox&                    a_interp_node_data,
                                                FArrayBox&                    a_interp_node_coords,
                                                const LevelData<FArrayBox>&   a_interp_node_data_decomp,
                                                const LevelData<FArrayBox>&   a_interp_node_coords_decomp,
                                                const Vector<Box>&            a_grid_boxes,
                                                const int&                    a_decomp_num ) const
{
   // We assume that each proc contains only one box (we created grids that way)
   // Also, the processor IDs assigned to grids correspond to 0 .. a_decomp_num-1
   
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
	       //inspect(this_interp_node_data);
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
SingleNullBlockCoordSys::mappedCoordNewton( const RealVect&  a_X,
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
      cout << "On block (" << m_poloidal_block << "," << m_toroidal_sector << "): Newton solve did not converge at " 
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

inline double sgn(double val) {
   return val >=0.? 1.: -1.;
}

bool
SingleNullBlockCoordSys::containsPoloidalPoint( const POL::RealVect& a_X,
                                                const bool           a_verbose ) const
{
   // Determine if the physical point a_X (in (R,Z) coordinates) is contained in the block

   double xi_0_lo = lowerMappedCoordinate(RADIAL_DIR); 
   double xi_0_hi = upperMappedCoordinate(RADIAL_DIR);
   double xi_1_lo = lowerMappedCoordinate(POLOIDAL_DIR);
   double xi_1_hi = upperMappedCoordinate(POLOIDAL_DIR);
   
   // Calculation of the winding number involves tracing the block boundary
   // taking discrete steps.  We use the mapping grid size, arbitrarily
   // multiplied by 10, to set the step size.  The block boundaries are pretty
   // smooth, but if any problems arise requiring more resolution, the
   // multiplier can be increased here.
   POL::IntVect mapping_dims = m_poloidal_util->getMappingDims();
   int multiplier = 10;

   POL::IntVect nsteps = multiplier * mapping_dims;

   double inc_0 = (xi_0_hi - xi_0_lo) / nsteps[0];
   double inc_1 = (xi_1_hi - xi_1_lo) / nsteps[1];

   double theta_total = 0.;

   POL::RealVect xi(xi_0_lo,xi_1_lo);
   POL::RealVect old_ray;
   old_ray[0] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::R_VAR, POL::Interp::FUNCTION, xi) - a_X[0];
   old_ray[1] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::Z_VAR, POL::Interp::FUNCTION, xi) - a_X[1];
   double len = old_ray.vectorLength();
   if (len == 0.) {
      return true;
   }
   else {
      old_ray /= len;
   }

   // Bottom
   for (int i=0; i<nsteps[0]; ++i) {
      xi[0] += inc_0;
      POL::RealVect new_ray;
      new_ray[0] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::R_VAR, POL::Interp::FUNCTION, xi) - a_X[0];
      new_ray[1] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::Z_VAR, POL::Interp::FUNCTION, xi) - a_X[1];
      len = new_ray.vectorLength();
      if (len == 0.) {
         return true;
      }
      else {
         new_ray /= len;
      }
      double product = new_ray.dotProduct(old_ray);
      if ( fabs(product) < 1. ) {
         theta_total += acos(product) * sgn(old_ray[0]*new_ray[1] - old_ray[1]*new_ray[0]);
      }
      if ( a_verbose ) cout << "bottom " << i << " " << theta_total << endl;
      old_ray = new_ray;
   }

   // Right
   for (int i=0; i<nsteps[1]; ++i) {
      xi[1] += inc_1;
      POL::RealVect new_ray;
      new_ray[0] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::R_VAR, POL::Interp::FUNCTION, xi) - a_X[0];
      new_ray[1] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::Z_VAR, POL::Interp::FUNCTION, xi) - a_X[1];
      len = new_ray.vectorLength();
      if (len == 0.) {
         return true;
      }
      else {
         new_ray /= len;
      }
      double product = new_ray.dotProduct(old_ray);
      if ( fabs(product) < 1. ) {
         theta_total += acos(product) * sgn(old_ray[0]*new_ray[1] - old_ray[1]*new_ray[0]);
      }
      if ( a_verbose ) cout << "right " << i << " " << theta_total << endl;
      old_ray = new_ray;
   }

   // Top
   for (int i=0; i<nsteps[0]; ++i) {
      xi[0] -= inc_0;
      POL::RealVect new_ray;
      new_ray[0] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::R_VAR, POL::Interp::FUNCTION, xi) - a_X[0];
      new_ray[1] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::Z_VAR, POL::Interp::FUNCTION, xi) - a_X[1];
      len = new_ray.vectorLength();
      if (len == 0.) {
         return true;
      }
      else {
         new_ray /= len;
      }
      double product = new_ray.dotProduct(old_ray);
      if ( fabs(product) < 1. ) {
         theta_total += acos(product) * sgn(old_ray[0]*new_ray[1] - old_ray[1]*new_ray[0]);
      }
      if ( a_verbose ) cout << "top " << i << " " << theta_total << endl;
      old_ray = new_ray;
   }

   // Left
   for (int i=0; i<nsteps[1]; ++i) {
      xi[1] -= inc_1;
      POL::RealVect new_ray;
      new_ray[0] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::R_VAR, POL::Interp::FUNCTION, xi) - a_X[0];
      new_ray[1] = m_poloidal_util->interpolate(POL::SingleNullPoloidalUtil::Z_VAR, POL::Interp::FUNCTION, xi) - a_X[1];
      len = new_ray.vectorLength();
      if (len == 0.) {
         return true;
      }
      else {
         new_ray /= len;
      }
      double product = new_ray.dotProduct(old_ray);
      if ( fabs(product) < 1. ) {
         theta_total += acos(product) * sgn(old_ray[0]*new_ray[1] - old_ray[1]*new_ray[0]);
      }
      if ( a_verbose ) cout << "left " << i << " " << theta_total << endl;
      old_ray = new_ray;
   }

   theta_total += 1.e-2;  // Accomodate a little roundoff

   return theta_total >= 2.*Pi;
}


void
SingleNullBlockCoordSys::convertCylindricalToCartesian( const RealVect&   a_cartesian_vect,
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
SingleNullBlockCoordSys::convertCylindricalToCartesian( const FArrayBox&  a_cartesian_coord,
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

void SingleNullBlockCoordSys::getToroidalCoords(FArrayBox& a_coords) const
{
  Box box = a_coords.box();

  FArrayBox X(box,SpaceDim);
  getCellCenteredRealCoords(X);

  RealVect coord;
  for (BoxIterator bit(box); bit.ok(); ++bit) {
    IntVect iv = bit();

    for (int dir=0; dir<SpaceDim; ++dir) {       
      coord[dir] = X(iv,dir);
    }

    convertCartesianToToroidal(coord);

    for (int dir=0; dir<SpaceDim; ++dir) {
      a_coords(iv,dir) = coord[dir];
    }
  }
}

void SingleNullBlockCoordSys::convertCartesianToToroidal(RealVect& a_coord) const
{
#if CFG_DIM == 3
  //get Cartesian coordinates
  Real x = a_coord[0];
  Real y = a_coord[1];
  Real z = a_coord[2];
   
  //compute toroidal coordinates

  Real flux_norm = getNormMagneticFlux(a_coord);
   
  POL::RealVect axis = m_poloidal_util->getMagAxis();
   
  Real R0 = axis[0];
  Real Z0 = axis[1];
   
  Real phi = atan2(y,x);
  if (phi < 0.) phi += 2.*Pi;
   
  Real r = sqrt(pow(sqrt(x*x+y*y)-R0,2)+pow(z-Z0,2));
  Real r_xy = sqrt(x*x + y*y);
  Real theta = asin((z-Z0)/r);
  if (r_xy < R0) theta = Pi - theta;
  if (theta < 0) theta += 2. * Pi;
   
  a_coord[RADIAL_DIR] = flux_norm;
  a_coord[TOROIDAL_DIR] = phi;
  a_coord[POLOIDAL_DIR] = theta;
#endif
}


POL::IntVect
SingleNullBlockCoordSys::restrictToPoloidal( const IntVect& a_v ) const
{
#if CFG_DIM==2
   return a_v;
#endif
#if CFG_DIM==3
   return POL::IntVect(a_v[RADIAL_DIR], a_v[POLOIDAL_DIR]);
#endif
}


POL::RealVect
SingleNullBlockCoordSys::restrictToPoloidal( const RealVect& a_v ) const
{
#if CFG_DIM==2
   return a_v;
#endif
#if CFG_DIM==3
   return POL::RealVect(a_v[RADIAL_DIR], a_v[POLOIDAL_DIR]);
#endif
}


POL::Box
SingleNullBlockCoordSys::restrictToPoloidal( const Box& a_box ) const
{
#if CFG_DIM==2
   return a_box;
#endif
#if CFG_DIM==3
   return POL::Box(restrictToPoloidal(a_box.smallEnd()),restrictToPoloidal(a_box.bigEnd()));
#endif
}


POL::RealVect
SingleNullBlockCoordSys::restrictPhysCoordToPoloidal( const RealVect& a_X ) const
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
      MayDay::Error("SingleNullBlockCoordSys::restrictPhysCoordToPoloidal(): unrecognized m_phys_coord_type");
   }
#endif

   return X_pol;
}


void
SingleNullBlockCoordSys::restrictPhysCoordsToPoloidal( const FArrayBox&  a_coords,
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
SingleNullBlockCoordSys::restrictPhysCoordsToPoloidal( const FArrayBox&  a_coords,
                                                       const int         a_toroidal_index,
                                                       POL::FArrayBox&   a_coords_pol ) const
{
#if CFG_DIM==2
   a_coords_pol.copy(a_coords);
#endif
#if CFG_DIM==3
   const Box& box = a_coords.box();
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
SingleNullBlockCoordSys::restrictMappedCoordsToPoloidal( const FArrayBox&  a_coords,
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
SingleNullBlockCoordSys::injectToroidal( const POL::RealVect  a_v ) const
{
#if CFG_DIM==2
   return a_v;
#endif
#if CFG_DIM==3
   return RealVect(a_v[0], lowerMappedCoordinate(TOROIDAL_DIR), a_v[1]);
#endif
}


void
SingleNullBlockCoordSys::injectToroidal( const POL::FArrayBox&  a_fab_pol,
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
