#include "LocalizedPinch.H"

#include "LogRectCoordSys.H"
#include "SNCoreCoordSys.H"
#include "SingleNullCoordSys.H"

#include "LocalizedPinchF_F.H"
#include "NamespaceHeader.H"

LocalizedPinch::LocalizedPinch(ParmParse& a_pp,
                     const int& a_verbosity )
   : GridFunction(a_verbosity),
     m_PreCoeffs(4,0.0),
     m_expCoeffs(4,0.0)
{
   parseParameters( a_pp );
}


inline
void LocalizedPinch::parseParameters( ParmParse& a_pp )
{
      
   a_pp.get("varType", m_varType );
   
   if(m_varType != "currentDensity" && m_varType != "magneticField") {
      cout << "m_varType = " << m_varType << " is not a valid type " << endl;
   }

   a_pp.getarr( "PreCoeffs", m_PreCoeffs, 0, 4 );

   a_pp.getarr( "expCoeffs", m_expCoeffs, 0, 4 );

   if (m_verbosity) {
      printParameters();
   }

}


void LocalizedPinch::checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   bool unknown_geometry( typeid(coord_sys) != typeid(LogRectCoordSys) );
   unknown_geometry &= (typeid(coord_sys) != typeid(SingleNullCoordSys));
   unknown_geometry &= (typeid(coord_sys) != typeid(SNCoreCoordSys));
   
   if ( unknown_geometry ) {
      const std::string msg( "LocalizedPinch: Attempt to use unknown geometry. ");
      MayDay::Error( msg.c_str() );
   }
}

void LocalizedPinch::setPointwise(FArrayBox&                 a_data,
                             const MultiBlockLevelGeom& a_geometry,
                             const FArrayBox&           a_real_coords,
                             const FArrayBox&           a_normalized_flux,
                             const int                  a_block_number) const
{
   
   const MagBlockCoordSys& coord_sys = getCoordSys(a_geometry, a_block_number);

   Box box( a_data.box() );
   FArrayBox cell_center_coords( box, SpaceDim );
   coord_sys.getCellCenteredRealCoords( cell_center_coords );
 
   if(m_varType=="magneticField") {
      FORT_SET_LOCALIZED_PINCHB( CHF_FRA(a_data),
                                 CHF_BOX(box),
                                 CHF_CONST_FRA1(cell_center_coords,0),
                                 CHF_CONST_VR(m_PreCoeffs),
                                 CHF_CONST_VR(m_expCoeffs) );
   }
   
   if(m_varType=="currentDensity") {
      FORT_SET_LOCALIZED_PINCHJ( CHF_FRA(a_data),
                                 CHF_BOX(box),
                                 CHF_CONST_FRA1(cell_center_coords,0),
                                 CHF_CONST_VR(m_PreCoeffs),
                                 CHF_CONST_VR(m_expCoeffs) );
   }


}


void LocalizedPinch::printParameters() const
{
   if (procID()==0) {
      std::cout << "LocalizedPinch grid function parameters:" << std::endl;
      std::cout << "  varType: "    << m_varType   << std::endl;
      std::cout << "  P00,P01,P10,P00: "   << m_PreCoeffs   << std::endl;
      std::cout << "  a,b,c,d: "    << m_expCoeffs    << std::endl;
      std::cout << std::endl;
   }
}


#include "NamespaceFooter.H"
