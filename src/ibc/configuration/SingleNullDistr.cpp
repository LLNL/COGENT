#include <math.h>
#include "SingleNullDistr.H"

#include <iostream>
#include <typeinfo>
#include <string>


#include "DataIterator.H"
#include "Directions.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "FourthOrderUtil.H"
#include "LevelData.H"
#include "MayDay.H"
#include "MillerCoordSys.H"
#include "MillerBlockCoordSys.H"
#include "SingleNullCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "MultiBlockCoordSys.H"
#include "Vector.H"
#include "MagGeom.H"

#include "NamespaceHeader.H"

inline
const MagBlockCoordSys& getCoordSys( const MultiBlockLevelGeom& a_geometry,
                                     const Box& a_box )
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   const int block_number( coord_sys.whichBlock( a_box ) );
   const NewCoordSys* block_coord_sys( coord_sys.getCoordSys( block_number ) );
   return static_cast<const MagBlockCoordSys&>( *block_coord_sys );
}


SingleNullDistr::SingleNullDistr( //const std::string& a_name,
                                 ParmParse& a_pp,
                                 const int& a_verbosity )
   : //m_name(a_name),
     m_verbosity(a_verbosity),
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
     m_boltzmann_equilibrium(false)


{
   parseParameters( a_pp );
}


void SingleNullDistr::assign( LevelData<FArrayBox>& a_data,
                              const MultiBlockLevelGeom& a_geometry,
                              const Real& a_time,
                              const bool& a_cell_averages ) const
{
   checkGeometryValidity( a_geometry );

   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
 
      int block_number( coord_sys.whichBlock( grids[dit] ) );

      if (a_cell_averages) {
         setCellAverages( a_data[dit], coord_sys, getCoordSys( a_geometry, grids[dit] ), block_number );
      }
      else {
         setPointwise( a_data[dit], coord_sys, getCoordSys( a_geometry, grids[dit] ), block_number );
      }
   }
   a_data.exchange();
}


void SingleNullDistr::assign( FArrayBox& a_data,
                              const MultiBlockLevelGeom& a_geometry,
                              const Box& a_box, // interior box
                              const Real& a_time,
                              const bool& a_cell_averages ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   const int block_number( coord_sys.whichBlock( a_box ) );

   if (a_cell_averages) {
      setCellAverages( a_data, coord_sys, getCoordSys( a_geometry, a_box ), block_number );
   }
   else {
      setPointwise( a_data, coord_sys, getCoordSys( a_geometry, a_box ), block_number );
   }
}


void SingleNullDistr::assign( LevelData<FArrayBox>& a_data,
                              const MultiBlockLevelGeom& a_geometry,
                              const BoundaryBoxLayout& a_bdry_layout,
                              const Real& a_time ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );

   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   // NB: This is a cheat - there's one too many cells at the (dir,side) face
   // of the boundary box, but it doesn't matter because one-sided difference
   // will be used at that face to construct the cell average.  We do want the
   // extra cell in all other directions.
   LevelData<FArrayBox> data_tmp( grids, a_data.nComp(), IntVect::Unit );
   for (DataIterator dit( grids ); dit.ok(); ++dit) {
      const Box box( a_bdry_layout.interiorBox( dit ) );
      const MagBlockCoordSys& block_coord_sys( ((MagGeom&)a_geometry).getBlockCoordSys( box ) );
      const int block_number( coord_sys.whichBlock( box ) );

      setPointwise( data_tmp[dit], coord_sys, block_coord_sys, block_number );
   }
   for (DataIterator dit( grids ); dit.ok(); ++dit) {
      Box domain_box( data_tmp[dit].box() );
      domain_box.growDir( a_bdry_layout.dir(), a_bdry_layout.side(), -1 );
      ProblemDomain domain( domain_box );
      fourthOrderAverageCell( data_tmp[dit], domain, grids[dit] );
   }
   data_tmp.copyTo( a_data );
   a_data.exchange();
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



   if (m_verbosity) {
      printParameters();
   }
}


inline
void SingleNullDistr::checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const
{
   const MultiBlockCoordSys& sn_coord_sys( *(a_geometry.coordSysPtr()) );
   bool not_single_null( typeid(sn_coord_sys) != typeid(SingleNullCoordSys) );
   not_single_null &= (typeid(sn_coord_sys) != typeid(SingleNullBlockCoordSys));

   if ( not_single_null ) {
      const std::string msg( "SingleNullDistr: Attempt to use unknown geometry. ");
      MayDay::Error( msg.c_str() );
   }
}


inline
void SingleNullDistr::setCellAverages( FArrayBox&        a_data,
                                       const MultiBlockCoordSys& a_coords,
                                       const MagBlockCoordSys&   a_block_coord,
                                       const int               block_number ) const

{
   Box box( a_data.box() );
   Box tmp_box( box );
   tmp_box.grow( IntVect::Unit );
   FArrayBox tmp( tmp_box, a_data.nComp() );

   setPointwise( tmp, a_coords, a_block_coord, block_number );

   fourthOrderAverageCell( tmp, a_block_coord.domain(), box );

   a_data.copy( tmp, box );
}


inline
void SingleNullDistr::setPointwise( FArrayBox&          a_data,
                                    const MultiBlockCoordSys& a_coords,
                                    const MagBlockCoordSys&   a_block_coords,
                                    const int                 block_number)  const

{
   Box box( a_data.box() );
   FArrayBox cell_center_coords( box, SpaceDim );
   a_block_coords.getCellCenteredMappedCoords( cell_center_coords );

   //Set the distribution for (LCORE,RCORE,LCSOL,RCSOL,LSOL,RSOL) blocks 
   if (block_number <= RSOL) { 
     BoxIterator bit(box);
     for (bit.begin(); bit.ok(); ++bit) {
       IntVect iv = bit();

       RealVect mapped_coord;
       mapped_coord[0] = cell_center_coords (iv,0);
       mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))
                                 ->lowerMappedCoordinate(1);

       //Get a physical Z coordinate at the top of a tokamak corresponding to this flux label  
       Real this_Z;

       if ((block_number == RCORE)||(block_number == LCORE)) {
         this_Z = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))
                         ->realCoord(mapped_coord)[1];
       }

       else {
         this_Z = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCSOL))
                         ->realCoord(mapped_coord)[1];
       }

       a_data(iv,0) = radialCoreSol(this_Z, a_coords);

     }
   }

   //Set the distribution for LPF and RPF blocks
   else {
     BoxIterator bit(box);
     for (bit.begin(); bit.ok(); ++bit) {
       IntVect iv = bit();
       RealVect mapped_coord;
       mapped_coord[0] = cell_center_coords (iv,0);
       mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LPF))
                                 ->lowerMappedCoordinate(1);

       //Get a physical Z coordinate at the top of a tokamak corresponding to this flux label  
       Real this_Z = ((const MagBlockCoordSys*)a_coords.getCoordSys(LPF))
                         ->realCoord(mapped_coord)[1];

       a_data(iv,0) = radialPF(this_Z, a_coords);

     }

   }

}
inline
Real SingleNullDistr::radialCoreSol(const Real                 a_Z,
                                    const MultiBlockCoordSys&  a_coords)  const

{

   //Get necessary global geometry parameters
   RealVect mapped_coord;

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->lowerMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->lowerMappedCoordinate(1);
   Real Z_core = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->upperMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->lowerMappedCoordinate(1);
   Real Z_sep = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCSOL))->upperMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCSOL))->lowerMappedCoordinate(1);
   Real Z_sol = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCSOL))->realCoord(mapped_coord)[1];

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

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->upperMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->upperMappedCoordinate(1);
   Real Z_xpt = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LPF))->lowerMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LPF))->lowerMappedCoordinate(1);
   Real Z_pf = ((const MagBlockCoordSys*)a_coords.getCoordSys(LPF))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->lowerMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->lowerMappedCoordinate(1);
   Real Z_core = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->upperMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->lowerMappedCoordinate(1);
   Real Z_sep = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCORE))->realCoord(mapped_coord)[1];

   mapped_coord[0] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCSOL))->upperMappedCoordinate(0);
   mapped_coord[1] = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCSOL))->lowerMappedCoordinate(1);
   Real Z_sol = ((const MagBlockCoordSys*)a_coords.getCoordSys(LCSOL))->realCoord(mapped_coord)[1];

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
