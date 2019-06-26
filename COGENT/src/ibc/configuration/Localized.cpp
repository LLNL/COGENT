#include "Localized.H"

#include "LogRectCoordSys.H"
#include "SNCoreCoordSys.H"
#include "SingleNullCoordSys.H"

#include "LocalizedF_F.H"
#include "NamespaceHeader.H"

Localized::Localized(ParmParse& a_pp,
                     const int& a_verbosity )
   : GridFunction(a_verbosity),
     m_amplitude(1.0)
{
   parseParameters( a_pp );
}


inline
void Localized::parseParameters( ParmParse& a_pp )
{
   a_pp.get( "amplitude", m_amplitude );

   Vector<Real> temp( SpaceDim );

   temp.assign( 0.0 );
   a_pp.getarr( "location", temp, 0, SpaceDim );
   m_location = RealVect( temp );

   temp.assign( 1.0 );
   a_pp.getarr( "width", temp, 0, SpaceDim );
   m_width = RealVect( temp );

   a_pp.get( "floor", m_floor );

   if (m_verbosity) {
      printParameters();
   }
}


void Localized::checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   bool unknown_geometry( typeid(coord_sys) != typeid(LogRectCoordSys) );
   unknown_geometry &= (typeid(coord_sys) != typeid(SingleNullCoordSys));
   unknown_geometry &= (typeid(coord_sys) != typeid(SNCoreCoordSys));
   
   if ( unknown_geometry ) {
      const std::string msg( "Localized: Attempt to use unknown geometry. ");
      MayDay::Error( msg.c_str() );
   }
}

void Localized::setPointwise(FArrayBox&                 a_data,
                             const MultiBlockLevelGeom& a_geometry,
                             const FArrayBox&           a_real_coords,
                             const FArrayBox&           a_normalized_flux,
                             const int                  a_block_number) const
{
   
   const MagBlockCoordSys& coord_sys = getCoordSys(a_geometry, a_block_number);

   Box box( a_data.box() );
   FArrayBox cell_center_coords( box, SpaceDim );
   coord_sys.getCellCenteredRealCoords( cell_center_coords );
   FORT_SET_LOCALIZED( CHF_FRA(a_data),
                       CHF_BOX(box),
                       CHF_CONST_FRA(cell_center_coords),
                       CHF_CONST_REAL(m_amplitude),
                       CHF_CONST_REALVECT(m_location),
                       CHF_CONST_REALVECT(m_width),
                       CHF_CONST_REAL(m_floor));
}


void Localized::printParameters() const
{
   if (procID()==0) {
      std::cout << "Localized grid function parameters:" << std::endl;
      std::cout << "  amplitude: "   << m_amplitude   << std::endl;
      std::cout << "  location: "    << m_location    << std::endl;
      std::cout << "  width: "       << m_width       << std::endl;
      std::cout << "  floor: "       << m_floor       << std::endl;
      std::cout << std::endl;
   }
}


#include "NamespaceFooter.H"
