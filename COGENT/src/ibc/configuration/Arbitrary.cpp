#include "Arbitrary.H"

#include "ToroidalBlockCoordSys.H"
#include "SNCoreCoordSys.H"
#include "SingleNullCoordSys.H"
#include "SNCoreBlockCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "SNCoreBlockCoordSysModel.H"
#include "SingleNullBlockCoordSysModel.H"

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
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   
   if (m_coord_type == "flux" || m_coord_type == "outer_midplane") {
      bool not_sn( typeid(coord_sys) != typeid(SNCoreCoordSys) );
      not_sn &= (typeid(coord_sys) != typeid(SingleNullCoordSys));
      if ( not_sn ) {
         const std::string msg( "Arbitrary: Attempt to use not a single-null geometry with the flux or outer_midplane options. ");
         MayDay::Error( msg.c_str() );
      }
   }
}

void Arbitrary::setPointwise( FArrayBox&                 a_data,
                              const MultiBlockLevelGeom& a_geometry,
                              const FArrayBox&           a_real_coords,
                              const FArrayBox&           a_normalized_flux,
                              const int                  a_block_number) const

{
   CH_TIMERS("Arbitrary::setPointwise");
   CH_TIMER("convertCartesianToToroidal", t_convert_cartesian);
   CH_TIMER("calc3d", t_calc3d);
   CH_TIMER("getflux", t_getflux);

   const MagBlockCoordSys& coord_sys = getCoordSys(a_geometry, a_block_number);
   
   Box box( a_data.box() );
   FArrayBox cc_mapped_coords( box, SpaceDim );
   coord_sys.getCellCenteredMappedCoords( cc_mapped_coords );

   a_data.setVal(0.0);
   
   BoxIterator bit(a_data.box());
   for (bit.begin(); bit.ok(); ++bit)
   {
      IntVect iv = bit();
      RealVect loc(iv);
      RealVect phys_coordinate(iv);
        
      for (int dir=0; dir<SpaceDim; dir++) {
         phys_coordinate[dir] = a_real_coords(iv, dir);
         if (m_coord_type == "physical" || m_coord_type ==  "toroidal")  loc[dir] = a_real_coords(iv, dir);
         if (m_coord_type == "mapped" || m_coord_type == "flux" || m_coord_type == "outer_midplane")  loc[dir] = cc_mapped_coords(iv, dir);
      }
#if CFG_DIM == 3      
      CH_START(t_convert_cartesian);
      if (m_coord_type == "toroidal") {
         if (typeid(coord_sys) == typeid(ToroidalBlockCoordSys)) {
            ((const ToroidalBlockCoordSys&)coord_sys).convertCartesianToToroidal(loc);
         }
        if (typeid(coord_sys) == typeid(SingleNullBlockCoordSys)) {
          ((const SingleNullBlockCoordSys&)coord_sys).convertCartesianToToroidal(loc);
        }
      }
      CH_STOP(t_convert_cartesian);
#endif
      if (m_coord_type == "flux") {
         CH_START(t_getflux);
         loc[0] = a_normalized_flux(iv,0);
         CH_STOP(t_getflux);
      }
      
      if (m_coord_type == "outer_midplane") {
         
          if (typeid(coord_sys) == typeid(SNCoreBlockCoordSys)) {
             double fluxNorm  = ((const SNCoreBlockCoordSys&)coord_sys).getNormMagneticFlux(phys_coordinate);
             double Rsep  = ((const SNCoreBlockCoordSys&)coord_sys).getOuterRsep();
             loc[0] = ((const SNCoreBlockCoordSys&)coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;

          }

         if (typeid(coord_sys) == typeid(SingleNullBlockCoordSys)) {
            double fluxNorm  = ((const SingleNullBlockCoordSys&)coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SingleNullBlockCoordSys&)coord_sys).getOuterRsep();
            loc[0] = ((const SingleNullBlockCoordSys&)coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
         }

         if (typeid(coord_sys) == typeid(SNCoreBlockCoordSysModel)) {
            double fluxNorm  = ((const SNCoreBlockCoordSysModel&)coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SNCoreBlockCoordSysModel&)coord_sys).getOuterRsep();
            loc[0] = ((const SNCoreBlockCoordSysModel&)coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
            
         }
         
         if (typeid(coord_sys) == typeid(SingleNullBlockCoordSysModel)) {
            double fluxNorm  = ((const SingleNullBlockCoordSysModel&)coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SingleNullBlockCoordSysModel&)coord_sys).getOuterRsep();
            loc[0] = ((const SingleNullBlockCoordSysModel&)coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
         }

      }
      
#if CFG_DIM==3
       CH_START(t_calc3d);
       Real val = m_pscore->calc3d(loc[0],loc[1],loc[2]);
       CH_STOP(t_calc3d);
#else
       Real val = m_pscore->calc2d(loc[0],loc[1]); // JRA, this is really slow
#endif

       if (m_function2.compare("UNDEFINED") != 0 &&
	   (a_block_number == SingleNullBlockCoordSys::LPF || a_block_number == SingleNullBlockCoordSys::RPF) ) {

#if CFG_DIM==3
       CH_START(t_calc3d);
	 val = m_pscore2->calc3d(loc[0],loc[1],loc[2]);
       CH_STOP(t_calc3d);
#else
	 val = m_pscore2->calc2d(loc[0],loc[1]);
#endif
       }

       a_data(iv,0) += val;

   }
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
