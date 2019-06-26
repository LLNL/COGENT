#include "Cosine.H"

#include "CosineF_F.H"
#include "Directions.H"
#include "ConstFact.H"
#include "LogRectCoordSys.H"

#include "NamespaceHeader.H"

Cosine::Cosine( ParmParse& a_pp,
                const int& a_verbosity )
   : GridFunction(a_verbosity),
     m_constant(0.0),
     m_amplitude(0.0)
{
   parseParameters( a_pp );
}


inline
void Cosine::parseParameters( ParmParse& a_pp )
{
   bool enforce_positivity(false);
   a_pp.query( "enforce_positivity", enforce_positivity );

   a_pp.get( "constant", m_constant );
   if (enforce_positivity) {
      CH_assert( m_constant>=0 );
   }

   a_pp.get( "amplitude", m_amplitude );
   if (enforce_positivity) {
      CH_assert( (m_constant - fabs( m_amplitude ))>=0 );
   }

   Vector<Real> temp( CFG_DIM );
   temp.assign( 0.0 );
   a_pp.getarr( "mode", temp, 0, CFG_DIM );
   m_mode = RealVect( temp );

   temp.assign( 0.0 );
   a_pp.getarr( "phase", temp, 0, CFG_DIM );
   m_phase = RealVect( temp );

   if (m_verbosity) {
      printParameters();
   }
}


void Cosine::checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   bool not_single_block( typeid(coord_sys) != typeid(LogRectCoordSys) );
   if ( not_single_block ) {
      const std::string msg( "Cosine: Attempt to use non-single-block geometry. ");
      MayDay::Error( msg.c_str() );
   }
}


void Cosine::setPointwise(FArrayBox&                 a_data,
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

   RealVect k( m_mode );
   k[RADIAL_DIR] *= 2.0 * Constants::PI;
   k[RADIAL_DIR] /= ( coord_sys.upperMappedCoordinate(RADIAL_DIR) - coord_sys.lowerMappedCoordinate(RADIAL_DIR) );

   Real shift = -0.5*( coord_sys.upperMappedCoordinate(RADIAL_DIR) + coord_sys.lowerMappedCoordinate(RADIAL_DIR));
   
   FORT_SET_COSINE( CHF_FRA(a_data),
                    CHF_BOX(box),
                    CHF_CONST_FRA(cell_center_coords),
                    CHF_CONST_REALVECT(dx),
                    CHF_CONST_REAL(shift),
                    CHF_CONST_REAL(m_constant),
                    CHF_CONST_REAL(m_amplitude),
                    CHF_CONST_REALVECT(k),
                    CHF_CONST_REALVECT(m_phase) );
}


void Cosine::printParameters() const
{
   if (procID()==0) {
      std::cout << "Cosine grid function parameters:" << std::endl;
      std::cout << "  constant: "  << m_constant  << std::endl;
      std::cout << "  amplitude: " << m_amplitude << std::endl;
      std::cout << "  mode: "      << m_mode      << std::endl;
      std::cout << "  phase: "     << m_phase     << std::endl;
      std::cout << std::endl;
   }
}


#include "NamespaceFooter.H"
