#include "Arbitrary.H"

#include <iostream>
#include <typeinfo>
#include <string>

#include "ConstFact.H"
#include "DataIterator.H"
#include "Directions.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "FourthOrderUtil.H"
#include "LevelData.H"
#include "MayDay.H"
#include "SNCoreBlockCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "SNCoreBlockCoordSysModel.H"
#include "SingleNullBlockCoordSysModel.H"
#include "MultiBlockCoordSys.H"
#include "Vector.H"
#include "MagGeom.H"

#include "ParsingCore.H"

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


Arbitrary::Arbitrary( ParmParse& a_pp,
                    const int& a_verbosity )
   : m_verbosity(a_verbosity),
     m_function("UNDEFINED"),
     m_coord_type("mapped"),
     m_pscore2(NULL),
     m_function2("UNDEFINED")
{
   parseParameters( a_pp );

   const char *userFormular=m_function.c_str();
   
   m_pscore = new ParsingCore(userFormular);

   if (m_function2.compare("UNDEFINED") != 0) {
     const char *userFormular2=m_function2.c_str();
     m_pscore2 = new ParsingCore(userFormular2);
   }

#if DEBUG_5D==1///DEBUG_5D
      cout << " Arbitrary::Arbitrary() m_pscore->getFormula():" << m_pscore->getFormula() << std::endl;
      cout << " Arbitrary::Arbitrary() m_pscore->getManipStr():" << m_pscore->getManipStr()<< std::endl;
      cout << " Arbitrary::Arbitrary() m_pscore->getPostStr():" << m_pscore->getPostStr()<< std::endl;
#endif//DEBUG_5D


}


void Arbitrary::assign( LevelData<FArrayBox>& a_data,
                        const MultiBlockLevelGeom& a_geometry,
                        const Real& a_time,
                        const bool& a_cell_averages ) const
{

   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const int block_number( coord_sys.whichBlock( grids[dit] ) );
      if (a_cell_averages) {
 	 setCellAverages( a_data[dit], getCoordSys( a_geometry, grids[dit] ), block_number );
      }
      else {
	 setPointwise( a_data[dit], getCoordSys( a_geometry, grids[dit] ), block_number );
      }
   }
   a_data.exchange();
}


void Arbitrary::assign( LevelData<FluxBox>& a_data,
                        const MultiBlockLevelGeom& a_geometry,
                        const Real& a_time,
                        const bool& a_cell_averages ) const
{

   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const int block_number( coord_sys.whichBlock( grids[dit] ) );
      for (int dir=0; dir<SpaceDim; ++dir) {
         if (a_cell_averages) {
            setCellAverages( a_data[dit][dir], getCoordSys( a_geometry, grids[dit] ), block_number );
         }
         else {
            setPointwise( a_data[dit][dir], getCoordSys( a_geometry, grids[dit] ), block_number );
         }
      }
   }
   a_data.exchange();
}


void Arbitrary::assign( FArrayBox& a_data,
                        const MultiBlockLevelGeom& a_geometry,
                        const Box& a_box,
                        const Real& a_time,
                        const bool& a_cell_averages ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   const int block_number( coord_sys.whichBlock( a_box ) );

   if (a_cell_averages) {
      setCellAverages( a_data, getCoordSys( a_geometry, a_box ), block_number );
   }
   else {
      setPointwise( a_data, getCoordSys( a_geometry, a_box ), block_number );
   }
}


void Arbitrary::assign( FluxBox& a_data,
                        const MultiBlockLevelGeom& a_geometry,
                        const Box& a_box,
                        const Real& a_time,
                        const bool& a_cell_averages ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   const int block_number( coord_sys.whichBlock( a_box ) );

   if (a_cell_averages) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         setCellAverages( a_data[dir], getCoordSys( a_geometry, a_box ), block_number );
      }
   }
   else {
      for (int dir=0; dir<SpaceDim; ++dir) {
         setPointwise( a_data[dir], getCoordSys( a_geometry, a_box ), block_number );
      }
   }
}


void Arbitrary::assign( LevelData<FArrayBox>& a_data,
                        const MultiBlockLevelGeom& a_geometry,
                        const BoundaryBoxLayout& a_bdry_layout,
                        const Real& a_time ) const
{
   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const MultiBlockCoordSys& coords( *(a_geometry.coordSysPtr()) );

   // NB: This is a cheat - there's one too many cells at the (dir,side) face
   // of the boundary box, but it doesn't matter because one-sided difference
   // will be used at that face to construct the cell average.  We do want the
   // extra cell in all other directions.
   LevelData<FArrayBox> data_tmp( grids, a_data.nComp(), IntVect::Unit );
   for (DataIterator dit( grids ); dit.ok(); ++dit) {
      const Box box( a_bdry_layout.interiorBox( dit ) );
      const MagBlockCoordSys& coord_sys( ((MagGeom&)a_geometry).getBlockCoordSys( box ) );
      const int block_number( coords.whichBlock( grids[dit] ) );
      setPointwise( data_tmp[dit], coord_sys, block_number );
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
void Arbitrary::parseParameters( ParmParse& a_pp )
{
   a_pp.get( "function", m_function);
   a_pp.query( "second_function", m_function2);
   a_pp.query( "coordinate_type", m_coord_type );

   
   if (m_verbosity) {
      printParameters();
   }
}

inline
void Arbitrary::setCellAverages( FArrayBox&              a_data,
				 const MagBlockCoordSys& a_coord_sys, 
				 const int               a_block_number) const
{
   Box box( a_data.box() );
   Box tmp_box( box );
   tmp_box.grow( IntVect::Unit );
   FArrayBox tmp( tmp_box, a_data.nComp() );

   setPointwise( tmp, a_coord_sys, a_block_number );

   fourthOrderAverageCell( tmp, a_coord_sys.domain(), box );

   a_data.copy( tmp, box );
}


inline
void Arbitrary::setPointwise( FArrayBox&              a_data,
			      const MagBlockCoordSys& a_coord_sys, 
			      const int               a_block_number) const

{
   //Check geometry validity for "flux" or "outer_midplane" types
   if (m_coord_type == "flux" || m_coord_type == "outer_midplane") {
      bool not_sn( typeid(a_coord_sys) != typeid(SNCoreBlockCoordSys) );
      not_sn &= (typeid(a_coord_sys) != typeid(SNCoreBlockCoordSysModel));
      not_sn &= (typeid(a_coord_sys) != typeid(SingleNullBlockCoordSys));
      not_sn &= (typeid(a_coord_sys) != typeid(SingleNullBlockCoordSysModel));
      if ( not_sn ) {
         const std::string msg( "Arbitrary: Attempt to use not a single-null geometry with the flux or outer_midplane options. ");
         MayDay::Error( msg.c_str() );
      }
   }
   
   Box box( a_data.box() );
   FArrayBox cc_mapped_coords( box, SpaceDim );
   a_coord_sys.getCellCenteredMappedCoords( cc_mapped_coords );

   FArrayBox cc_phys_coords( box, SpaceDim );
   a_coord_sys.getCellCenteredRealCoords( cc_phys_coords );
   
   a_data.setVal(0.0);

   BoxIterator bit(a_data.box());
   for (bit.begin(); bit.ok(); ++bit)
   {
      IntVect iv = bit();
      RealVect loc(iv);
      RealVect phys_coordinate(iv);
      
      for (int dir=0; dir<SpaceDim; dir++) {
         phys_coordinate[dir] = cc_phys_coords(iv, dir);
         if (m_coord_type == "physical")  loc[dir] = cc_phys_coords(iv, dir);
         if (m_coord_type == "mapped" || m_coord_type == "flux" || m_coord_type == "outer_midplane")  loc[dir] = cc_mapped_coords(iv, dir);
      }
      
      if (m_coord_type == "flux") {
         if (typeid(a_coord_sys) == typeid(SNCoreBlockCoordSys))  loc[0] = ((const SNCoreBlockCoordSys&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
         if (typeid(a_coord_sys) == typeid(SingleNullBlockCoordSys))  loc[0] = ((const SingleNullBlockCoordSys&)a_coord_sys).getNormMagneticFlux(phys_coordinate);

         if (typeid(a_coord_sys) == typeid(SNCoreBlockCoordSysModel))  loc[0] = ((const SNCoreBlockCoordSysModel&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
         if (typeid(a_coord_sys) == typeid(SingleNullBlockCoordSysModel))  loc[0] = ((const SingleNullBlockCoordSysModel&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
      }
      
      if (m_coord_type == "outer_midplane") {
         
          if (typeid(a_coord_sys) == typeid(SNCoreBlockCoordSys)) {
             double fluxNorm  = ((const SNCoreBlockCoordSys&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
             double Rsep  = ((const SNCoreBlockCoordSys&)a_coord_sys).getOuterRsep();
             loc[0] = ((const SNCoreBlockCoordSys&)a_coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;

          }

         if (typeid(a_coord_sys) == typeid(SingleNullBlockCoordSys)) {
            double fluxNorm  = ((const SingleNullBlockCoordSys&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SingleNullBlockCoordSys&)a_coord_sys).getOuterRsep();
            loc[0] = ((const SingleNullBlockCoordSys&)a_coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
         }

         if (typeid(a_coord_sys) == typeid(SNCoreBlockCoordSysModel)) {
            double fluxNorm  = ((const SNCoreBlockCoordSysModel&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SNCoreBlockCoordSysModel&)a_coord_sys).getOuterRsep();
            loc[0] = ((const SNCoreBlockCoordSysModel&)a_coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
            
         }
         
         if (typeid(a_coord_sys) == typeid(SingleNullBlockCoordSysModel)) {
            double fluxNorm  = ((const SingleNullBlockCoordSysModel&)a_coord_sys).getNormMagneticFlux(phys_coordinate);
            double Rsep  = ((const SingleNullBlockCoordSysModel&)a_coord_sys).getOuterRsep();
            loc[0] = ((const SingleNullBlockCoordSysModel&)a_coord_sys).getOuterMidplaneCoord(fluxNorm) - Rsep;
         }

      }


      
#if CFG_DIM==3
       Real val = m_pscore->calc3d(loc[0],loc[1],loc[2]);
#else
       Real val = m_pscore->calc2d(loc[0],loc[1]);
#endif

       if (m_function2.compare("UNDEFINED") != 0 &&
	   (a_block_number == LPF || a_block_number == RPF) ) {

#if CFG_DIM==3
	 val = m_pscore2->calc3d(loc[0],loc[1],loc[2]);
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
//      std::cout << "  translated form: "  << m_pscore->getManipStr()<< std::endl;
//      std::cout << "  postfix form: "  << m_pscore->getPostStr()<< std::endl;
      std::cout << "  coordinate type: "  << m_coord_type << std::endl;
      std::cout << std::endl;
   }
}


#include "NamespaceFooter.H"
