#include "Arbitrary.H"
#include "MagGeom.H"
#include "SingleNullCoordSys.H"
#include "ParsingCore.H"

#include "NamespaceHeader.H"


Arbitrary::Arbitrary(ParmParse& a_pp,
                     const int& a_verbosity )
   : GridFunction(a_verbosity),
     m_pscore2(NULL),
     m_function("UNDEFINED"),
     m_function2("UNDEFINED"),
     m_coord_type("mapped")
{
   parseParameters( a_pp );

   const char *userFormular=m_function.c_str();
   
   m_pscore = new ParsingCore(userFormular);

   if (m_function2.compare("UNDEFINED") != 0) {
     const char *userFormular2=m_function2.c_str();
     m_pscore2 = new ParsingCore(userFormular2);
   }
}


inline
void Arbitrary::parseParameters( ParmParse& a_pp )
{
   a_pp.get( "function", m_function);
   a_pp.query( "second_function", m_function2);
   a_pp.query( "coordinate_type", m_coord_type );

   if (m_verbosity) {
      printParameters();
   }
}

void Arbitrary::checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const
{
   //const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
 
   //const MagGeom& mag_geom = (const MagGeom&) a_geometry;

   // do something when needed
}

void Arbitrary::setPointwise( FArrayBox&                 a_data,
                              const MultiBlockLevelGeom& a_geometry,
                              const FArrayBox&           a_real_coords,
                              const FArrayBox&           a_normalized_flux,
                              const int                  a_block_number) const

{
   CH_TIMERS("Arbitrary::setPointwise");
   CH_TIMER("assignParseFunction", t_assignParseFunction);
   CH_TIMER("getCoords", t_getCoords);

   const MagGeom& mag_geom = (const MagGeom&) a_geometry;
   const MagCoordSys& coord_sys = *(mag_geom.getCoordSys());
   
   const MagBlockCoordSys& block_coord_sys = getCoordSys(a_geometry, a_block_number);
   
   Box box( a_data.box() );
   FArrayBox cc_coords( box, SpaceDim );
   
   CH_START(t_getCoords);
   if (m_coord_type == "physical") {
      cc_coords.copy(a_real_coords);
   }
   else if (m_coord_type == "toroidal") {
      block_coord_sys.getToroidalCoords( cc_coords, a_real_coords, false);
   }
   else if (m_coord_type == "flux") {
      block_coord_sys.getToroidalCoords( cc_coords, a_real_coords, true);
   }
   else {
      block_coord_sys.getCellCenteredMappedCoords( cc_coords );
   }
   CH_STOP(t_getCoords);
      
   CH_START(t_assignParseFunction);
   
   BoxIterator bit(a_data.box());
   for (bit.begin(); bit.ok(); ++bit)
   {
      IntVect iv = bit();
      RealVect loc(iv);
      Real val;
      
      for (int dir=0; dir<SpaceDim; dir++) {
         loc[dir] = cc_coords(iv, dir);
      }

      if (SpaceDim == 3) {
         val = m_pscore->calc3d(loc[0],loc[1],loc[2]);
      }
      else {
         val = m_pscore->calc2d(loc[0],loc[1]); // JRA, this is really slow
      }
      
      // Set data in the private flux region
      if (m_function2.compare("UNDEFINED") != 0 && ((const SingleNullCoordSys&)coord_sys).isPF(a_block_number) ) {
         if (SpaceDim == 3) {
            val = m_pscore2->calc3d(loc[0],loc[1],loc[2]);
         }
         else {
            val = m_pscore2->calc2d(loc[0],loc[1]);
         }
      }
      a_data(iv,0) = val;
   }
   
   CH_STOP(t_assignParseFunction);
}


void Arbitrary::printParameters() const
{
   if (procID()==0) {
      std::cout << "Arbitrary grid function parameters:" << std::endl;
      std::cout << "  function: "  << m_function << std::endl;
      std::cout << "  function2: "  << m_function2 << std::endl;
      std::cout << "  coordinate type: "  << m_coord_type << std::endl;
      //      std::cout << "  translated form: "  << m_pscore->getManipStr()<< std::endl;
      //      std::cout << "  postfix form: "  << m_pscore->getPostStr()<< std::endl;
      //      std::cout << "  formula form:" << m_pscore->getFormula() << std::endl;
      std::cout << std::endl;
   }
}


#include "NamespaceFooter.H"
