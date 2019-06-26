#include "Tanh.H"

#include "TanhF_F.H"

#include <iostream>
#include <typeinfo>
#include <string>

#include "ConstFact.H"
#include "Directions.H"
#include "LogRectCoordSys.H"

#include "NamespaceHeader.H"

Tanh::Tanh(ParmParse& a_pp,
           const int& a_verbosity )

   : GridFunction(a_verbosity),
     m_inner_radial_value(0.0),
     m_outer_radial_value(0.0),
     m_radial_midpoint(0.0),
     m_radial_width(1.0)
{
   parseParameters( a_pp );
}


inline
void Tanh::parseParameters( ParmParse& a_pp )
{
   a_pp.get( "inner_radial_value", m_inner_radial_value );
   a_pp.get( "outer_radial_value", m_outer_radial_value );
   a_pp.get( "radial_midpoint",    m_radial_midpoint );
   a_pp.get( "radial_width",       m_radial_width );
   CH_assert( m_radial_width>=0 );

   if (m_verbosity) {
      printParameters();
   }
}

void Tanh::checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   bool not_single_block( typeid(coord_sys) != typeid(LogRectCoordSys) );
   if ( not_single_block ) {
      const std::string msg( "Tanh: Attempt to use non-single-block geometry. ");
      MayDay::Error( msg.c_str() );
   }
}

void Tanh::setPointwise(FArrayBox&                 a_data,
                        const MultiBlockLevelGeom& a_geometry,
                        const FArrayBox&           a_real_coords,
                        const FArrayBox&           a_normalized_flux,
                        const int                  a_block_number) const
{
   const MagBlockCoordSys& coord_sys = getCoordSys(a_geometry, a_block_number);

   Box box( a_data.box() );
   FArrayBox cell_center_coords( box, SpaceDim );
   coord_sys.getCellCenteredMappedCoords( cell_center_coords );
   RealVect dx = coord_sys.getMappedCellSize();

   Real inner_radial_coord( coord_sys.lowerMappedCoordinate(RADIAL_DIR) );
   Real outer_radial_coord( coord_sys.upperMappedCoordinate(RADIAL_DIR) );
   
   Real rad_midpoint = inner_radial_coord + m_radial_midpoint * (outer_radial_coord - inner_radial_coord);
   
   FORT_SET_TANH( CHF_FRA(a_data),
                  CHF_BOX(box),
                  CHF_CONST_FRA(cell_center_coords),
                  CHF_CONST_REALVECT(dx),
                  CHF_CONST_REAL(m_inner_radial_value),
                  CHF_CONST_REAL(inner_radial_coord),
                  CHF_CONST_REAL(m_outer_radial_value),
                  CHF_CONST_REAL(outer_radial_coord),
                  CHF_CONST_REAL(rad_midpoint),
                  CHF_CONST_REAL(m_radial_width) );
}


void Tanh::printParameters() const
{
   if (procID()==0) {
      std::cout << "Tanh grid function parameters:"           << std::endl;
      std::cout << "  inner_radial_value: " << m_inner_radial_value << std::endl;
      std::cout << "  outer_radial_value: " << m_outer_radial_value << std::endl;
      std::cout << "  radial_midpoint: "    << m_radial_midpoint    << std::endl;
      std::cout << "  radial_width: "       << m_radial_width       << std::endl;
      std::cout << std::endl;
   }
}


#include "NamespaceFooter.H"
