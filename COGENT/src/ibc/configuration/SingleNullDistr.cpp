#include "SingleNullDistr.H"

#include "Directions.H"
#include "SingleNullCoordSys.H"
#include "MagGeom.H"

#include "NamespaceHeader.H"

SingleNullDistr::SingleNullDistr(ParmParse& a_pp,
                                 const int& a_verbosity )
   : GridFunction(a_verbosity),
     m_core_value(1.0),
     m_sep_value(0.5),
     m_pf_value(0.1),
     m_subtype("Lorentz"),
     m_inner_radial_value(0.0),
     m_outer_radial_value(0.0),
     m_midpoint_fraction(0.5),
     m_radial_width(1.0),
     m_amplitude(1.0),
     m_floor(0.0),
     m_boltzmann_equilibrium(false),
     m_mag_surf_on(false)

{
   parseParameters( a_pp );
}


inline
void SingleNullDistr::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "subtype",            m_subtype );
   a_pp.query( "core_value",         m_core_value );
   a_pp.query( "sep_value",          m_sep_value );
   a_pp.query( "pf_value",           m_pf_value );
   a_pp.query( "inner_radial_value", m_inner_radial_value );
   a_pp.query( "outer_radial_value", m_outer_radial_value );
   a_pp.query( "midpoint_fraction",  m_midpoint_fraction );
   a_pp.query( "radial_width",       m_radial_width );
   a_pp.query( "amplitude",          m_amplitude );
   a_pp.query( "floor",              m_floor );
   a_pp.query( "boltzmann_equilibrium",    m_boltzmann_equilibrium );
   a_pp.query( "dealignment_correction",   m_mag_surf_on );



   if (m_verbosity) {
      printParameters();
   }
}


inline
void SingleNullDistr::checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const
{
   const MultiBlockCoordSys& sn_coord_sys( *(a_geometry.coordSysPtr()) );
   bool not_single_null( typeid(sn_coord_sys) != typeid(SingleNullCoordSys) );

   if ( not_single_null ) {
      const std::string msg( "SingleNullDistr: Attempt to use a not single null geometry. ");
      MayDay::Error( msg.c_str() );
   }
}

void SingleNullDistr::setPointwise(FArrayBox&                 a_data,
                                   const MultiBlockLevelGeom& a_geometry,
                                   const FArrayBox&           a_real_coords,
                                   const FArrayBox&           a_normalized_flux,
                                   const int                  a_block_number) const
{
   const MagGeom& mag_geom = (const MagGeom&) a_geometry;
   const MultiBlockCoordSys& coords( *(a_geometry.coordSysPtr()) );
   const MagBlockCoordSys& block_coord_sys = getCoordSys(a_geometry, a_block_number);
   
   Box box( a_data.box() );
   FArrayBox cell_center_coords( box, SpaceDim );
   block_coord_sys.getCellCenteredMappedCoords( cell_center_coords );
    
   FArrayBox cell_center_real_coords( box, SpaceDim );
   block_coord_sys.getCellCenteredRealCoords( cell_center_real_coords );


   //Set the distribution for (LCORE,RCORE,LCSOL,RCSOL,LSOL,RSOL) blocks 
   if (a_block_number <= SingleNullBlockCoordSys::RSOL) {
     BoxIterator bit(box);
     for (bit.begin(); bit.ok(); ++bit) {
       IntVect iv = bit();

       RealVect mapped_coord;
       mapped_coord[0] = cell_center_coords (iv,0);
       mapped_coord[1] = ((const MagBlockCoordSys*)coords.getCoordSys(SingleNullBlockCoordSys::LCORE))
                                 ->lowerMappedCoordinate(1);

       //Get a physical Z coordinate at the top of a tokamak corresponding to this flux label  
       Real this_Z;

       if ((a_block_number == SingleNullBlockCoordSys::RCORE)||(a_block_number == SingleNullBlockCoordSys::LCORE)) {
         this_Z = ((const MagBlockCoordSys*)coords.getCoordSys(SingleNullBlockCoordSys::LCORE))
                         ->realCoord(mapped_coord)[1];
       }

       else {
         this_Z = ((const MagBlockCoordSys*)coords.getCoordSys(SingleNullBlockCoordSys::LCSOL))
                         ->realCoord(mapped_coord)[1];
       }
         
       if (m_mag_surf_on && a_block_number<=SingleNullBlockCoordSys::RSOL) {
           RealVect real_coord;
           real_coord[0] = cell_center_real_coords(iv,0);
           real_coord[1] = cell_center_real_coords(iv,1);
           this_Z = mag_geom.getMagFS(real_coord);
       }

       a_data(iv,0) = radialCoreSol(this_Z, coords);

     }
   }

   //Set the distribution for LPF and RPF blocks
   else {
     BoxIterator bit(box);
     for (bit.begin(); bit.ok(); ++bit) {
       IntVect iv = bit();
       RealVect mapped_coord;
       mapped_coord[0] = cell_center_coords (iv,0);
       mapped_coord[1] = ((const MagBlockCoordSys*)coords.getCoordSys(SingleNullBlockCoordSys::LPF))
                                 ->lowerMappedCoordinate(1);

       //Get a physical Z coordinate at the top of a tokamak corresponding to this flux label  
       Real this_Z = ((const MagBlockCoordSys*)coords.getCoordSys(SingleNullBlockCoordSys::LPF))
                         ->realCoord(mapped_coord)[1];

       a_data(iv,0) = radialPF(this_Z, coords);

     }

   }

}
inline
Real SingleNullDistr::radialCoreSol(const Real                 a_Z,
                                    const MultiBlockCoordSys&  a_coords)  const

{

   //Get necessary global geometry parameters
   RealVect mapped_coord;

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->lowerMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->lowerMappedCoordinate(1);
   Real Z_core = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->upperMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->lowerMappedCoordinate(1);
   Real Z_sep = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCSOL))->upperMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCSOL))->lowerMappedCoordinate(1);
   Real Z_sol = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCSOL))->realCoord(mapped_coord)[1];

   double val;

   if (m_subtype=="Lorentz") {

        double a = (Z_sep - Z_core)/sqrt(m_core_value/m_sep_value - 1.0);
        val = m_core_value /((a_Z - Z_core) * (a_Z - Z_core)/(a*a) + 1.0 );        
   }   

   if (m_subtype=="Tanh") {

        val = formTanh(a_Z, Z_core, Z_sol);
   }

  
   if (m_subtype=="Localized") {
      
       double radial_midpoint = Z_core + (Z_sol-Z_core) * m_midpoint_fraction;
       val = m_amplitude * exp(- (a_Z - radial_midpoint) * (a_Z - radial_midpoint) / (m_radial_width * m_radial_width) ) + m_floor; 
   }

   return val;

}

inline
Real SingleNullDistr::radialPF(const Real                 a_Z,
                               const MultiBlockCoordSys&  a_coords)  const


{
   //Get necessary global geometry parameters
   RealVect mapped_coord;

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->upperMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->upperMappedCoordinate(1);
   Real Z_xpt = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LPF))->lowerMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LPF))->lowerMappedCoordinate(1);
   Real Z_pf = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LPF))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->lowerMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->lowerMappedCoordinate(1);
   Real Z_core = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->upperMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->lowerMappedCoordinate(1);
   Real Z_sep = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCORE))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCSOL))->upperMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCSOL))->lowerMappedCoordinate(1);
   Real Z_sol = ((const MagBlockCoordSys*)a_coords.getCoordSys(SingleNullBlockCoordSys::LCSOL))->realCoord(mapped_coord)[1];

   double val;

   if (m_subtype=="Lorentz") {

        val = m_sep_value + (m_sep_value - m_pf_value) * ((a_Z - Z_xpt)/(Z_xpt - Z_pf) );

   }
    
   if (m_subtype=="Tanh") {

        double tanh_sep_value = formTanh(Z_sep, Z_core, Z_sol);
        val = tanh_sep_value + (tanh_sep_value - m_pf_value) * ((a_Z - Z_xpt)/(Z_xpt - Z_pf) );

   }

   if (m_subtype=="Localized") {

        double radial_midpoint = Z_core + (Z_sol-Z_core) * m_midpoint_fraction;

        double sep_value = m_amplitude * exp(- (Z_sep - radial_midpoint) * (Z_sep - radial_midpoint)
                                               / (m_radial_width * m_radial_width) ) + m_floor;

        val = sep_value + (sep_value - m_pf_value) * ((a_Z - Z_xpt)/(Z_xpt - Z_pf) );

   }


   return val;


}


void SingleNullDistr::printParameters() const
{
   if (procID()==0) {
      std::cout << "SingleNullDistr grid function parameters:" << std::endl;
      std::cout << "  subtype: "   << m_subtype   << std::endl;
      if (m_subtype == "Lorentz") {
        std::cout << "  core_value: "   << m_core_value   << std::endl;
        std::cout << "  sep_value: "    << m_sep_value    << std::endl;
        std::cout << "  pf_value: "     << m_pf_value     << std::endl;
      }

      if (m_subtype == "Tanh") {
        std::cout << "  inner_radial_value: "   << m_inner_radial_value   << std::endl;
        std::cout << "  outer_radial_value: "   << m_outer_radial_value   << std::endl;
        std::cout << "  midpoint_fraction: "    << m_midpoint_fraction    << std::endl;
        std::cout << "  radial_width: "         << m_radial_width         << std::endl;
      }

      std::cout << std::endl;
   }
}

double SingleNullDistr::formTanh(double x, double x1, double x2) const 
{

      double radial_midpoint = x1 + (x2-x1) * m_midpoint_fraction;
      double invwidth = 1.0 / m_radial_width;
      double tanhx1 = tanh( (x1 - radial_midpoint) * invwidth );
      double tanhx2 = tanh( (x2 - radial_midpoint) * invwidth );
      double a = (m_inner_radial_value - m_outer_radial_value) / (tanhx1 - tanhx2);
      double b = m_outer_radial_value - a * tanhx2;
      double arg = (x - radial_midpoint) * invwidth;

      if (!m_boltzmann_equilibrium) {
         return a * tanh( arg ) + b;
      }

      else {
         return -log(a * tanh( arg ) + b);
      }

}


#include "NamespaceFooter.H"
