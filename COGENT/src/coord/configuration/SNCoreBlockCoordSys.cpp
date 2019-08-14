#include <array>
#include "SNCoreBlockCoordSys.H"
#include "Directions.H"

#include "NamespaceHeader.H"


const std::string SNCoreBlockCoordSys::pp_name = "sncoreRealisticGeom";


void
SNCoreBlockCoordSys::printInit() const
{
   if ( m_verbose && procID()==0 ) {
      cout << "Constructing single null ";
      switch (m_poloidal_block)
         {
         case MCORE:
           cout << "middle core";
           break;
         case LCORE:
           cout << "left core";
           break;
         case RCORE:
           cout << "right core";
           break;
         default:
           MayDay::Error("SNCoreBlockCoordSys::SNCoreBlockCoordSys(): Invalid block_type encountered");
         }

      cout << " block with global index space domain box = " << m_domain.domainBox() << endl;

      switch (m_poloidal_block)
         {
         case MCORE:
           cout << "Middle core mapped domain: " 
                << lowerMappedCoordinate(RADIAL_DIR) << " < xi_" << RADIAL_DIR << " < " << upperMappedCoordinate(RADIAL_DIR) << ", "
#if CFG_DIM==3
                << lowerMappedCoordinate(TOROIDAL_DIR) << " < xi_" << TOROIDAL_DIR << " < " << upperMappedCoordinate(TOROIDAL_DIR) << ", "
#endif
                << lowerMappedCoordinate(POLOIDAL_DIR) << " < xi_" << POLOIDAL_DIR << " < " << upperMappedCoordinate(POLOIDAL_DIR);
            break;
         case LCORE:
           cout << "Left core mapped domain: " 
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
         default:
           MayDay::Error("SNCoreBlockCoordSys::SNCoreBlockCoordSys(): Invalid block_type encountered");
         }

      cout << endl;
   }
}


bool
SNCoreBlockCoordSys::blockNameIsValid( const string& a_block_name ) const
{
   return (a_block_name == "mcore" && m_poloidal_block == MCORE) ||
          (a_block_name == "lcore" && m_poloidal_block == LCORE) ||
          (a_block_name == "rcore" && m_poloidal_block == RCORE);
}


void
SNCoreBlockCoordSys::definePoints( const ParmParse&  a_pp,
                                   int&              a_block_poloidal,
                                   const int         a_block_full_poloidal,
                                   const IntVect&    a_mapping_block_size,
                                   const int&        a_n_poloidal_extend,
                                   double &          a_dtheta,
                                   double *          a_theta_pts ) const
{
   a_dtheta *= (double)a_block_full_poloidal / (double)a_block_poloidal;

   if ( m_poloidal_block == MCORE ) {
      a_theta_pts[a_mapping_block_size[1]-1] =
         upperMappedCoordinate(POLOIDAL_DIR) + a_n_poloidal_extend*a_dtheta;
      for (int i=a_mapping_block_size[1]-1; i>0; --i) {
         a_theta_pts[i-1] = a_theta_pts[i] - a_dtheta;
      }
   }
   else if ( m_poloidal_block == LCORE ) {
      a_theta_pts[a_mapping_block_size[1]-1] =
         upperMappedCoordinate(POLOIDAL_DIR) + a_n_poloidal_extend*a_dtheta;
      for (int i=a_mapping_block_size[1]-1; i>0; --i) {
         a_theta_pts[i-1] = a_theta_pts[i] - a_dtheta;
      }
   }
   else if ( m_poloidal_block == RCORE ) {
      a_theta_pts[0] = lowerMappedCoordinate(POLOIDAL_DIR) - a_n_poloidal_extend*a_dtheta;
      for (int i=1; i<a_mapping_block_size[1]; ++i) {
         a_theta_pts[i] = a_theta_pts[i-1] + a_dtheta;
      }
   }
   else {
      MayDay::Error("SNCoreBlockCoordSys::definePoints(): Invalid poloidal_block encountered");
   }
}


bool
SNCoreBlockCoordSys::isValid( const RealVect&  a_xi,
                              const bool       a_omit_toroidal ) const
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


#include "NamespaceFooter.H"
