#include "Localized.H"

#include <iostream>
#include <typeinfo>
#include <string>

#include "LocalizedF_F.H"

#include "DataIterator.H"
#include "Directions.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "FourthOrderUtil.H"
#include "LevelData.H"
#include "MayDay.H"
#include "MagGeom.H"
#include "LogRectCoordSys.H"
#include "SNCoreCoordSys.H"
#include "SingleNullCoordSys.H"
#include "MultiBlockCoordSys.H"
#include "Vector.H"

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


Localized::Localized( //const std::string& a_name,
                          ParmParse& a_pp,
                          const int& a_verbosity )
   : //m_name(a_name),
     m_verbosity(a_verbosity),
     m_amplitude(1.0)
{
   parseParameters( a_pp );
}


void Localized::assign( LevelData<FArrayBox>& a_data,
                        const MultiBlockLevelGeom& a_geometry,
                        const Real& a_time,
                        const bool& a_cell_averages ) const
{
   checkGeometryValidity( a_geometry );

   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      if (a_cell_averages) {
         setCellAverages( a_data[dit], getCoordSys( a_geometry, grids[dit] ) );
      }
      else {
         setPointwise( a_data[dit], getCoordSys( a_geometry, grids[dit] ) );
      }
   }
   a_data.exchange();
}


void Localized::assign( LevelData<FluxBox>& a_data,
                        const MultiBlockLevelGeom& a_geometry,
                        const Real& a_time,
                        const bool& a_cell_averages ) const
{
   checkGeometryValidity( a_geometry );

   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );

   if (a_cell_averages) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            setCellAverages( a_data[dit][dir], getCoordSys( a_geometry, grids[dit] ) );
         }
      }
   }
   else {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            setPointwise( a_data[dit][dir], getCoordSys( a_geometry, grids[dit] ) );
         }
      }
   }
   a_data.exchange();
}


void Localized::assign( FArrayBox& a_data,
                        const MultiBlockLevelGeom& a_geometry,
                        const Box& a_box, // interior box
                        const Real& a_time,
                        const bool& a_cell_averages ) const
{
   if (a_cell_averages) {
      setCellAverages( a_data, getCoordSys( a_geometry, a_box ) );
   }
   else {
      setPointwise( a_data, getCoordSys( a_geometry, a_box ) );
   }
}


void Localized::assign( FluxBox& a_data,
                        const MultiBlockLevelGeom& a_geometry,
                        const Box& a_box, // interior box
                        const Real& a_time,
                        const bool& a_cell_averages ) const
{
   if (a_cell_averages) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         setCellAverages( a_data[dir], getCoordSys( a_geometry, a_box ) );
      }
   }
   else {
      for (int dir=0; dir<SpaceDim; ++dir) {
         setPointwise( a_data[dir], getCoordSys( a_geometry, a_box ) );
      }
   }
}


void Localized::assign( LevelData<FArrayBox>& a_data,
                        const MultiBlockLevelGeom& a_geometry,
                        const BoundaryBoxLayout& a_bdry_layout,
                        const Real& a_time ) const
{
   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   // NB: This is a cheat - there's one too many cells at the (dir,side) face
   // of the boundary box, but it doesn't matter because one-sided difference
   // will be used at that face to construct the cell average.  We do want the
   // extra cell in all other directions.
   LevelData<FArrayBox> data_tmp( grids, a_data.nComp(), IntVect::Unit );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const Box box( a_bdry_layout.interiorBox( dit ) );
      const MagBlockCoordSys& coord_sys( ((MagGeom&)a_geometry).getBlockCoordSys( box ) );
      setPointwise( data_tmp[dit], coord_sys );
   }
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      Box domain_box( data_tmp[dit].box() );
      domain_box.growDir( a_bdry_layout.dir(), a_bdry_layout.side(), -1 );
      ProblemDomain domain( domain_box );
      fourthOrderAverageCell( data_tmp[dit], domain, grids[dit] );
   }
   data_tmp.copyTo( a_data );
   a_data.exchange();
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


inline
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


inline
void Localized::setCellAverages( FArrayBox&              a_data,
                                 const MagBlockCoordSys& a_coord_sys ) const

{
   Box box( a_data.box() );
   Box tmp_box( box );
   tmp_box.grow( IntVect::Unit );
   FArrayBox tmp( tmp_box, a_data.nComp() );

   setPointwise( tmp, a_coord_sys );

#if 0
   fourthOrderAverageCell( tmp, a_coord_sys.domain(), box );
#else
   Box domain_box( a_coord_sys.domain().domainBox() );
   domain_box.grow(4*IntVect::Unit);
   ProblemDomain domain(domain_box);
   fourthOrderAverageCell( tmp, domain, box );
#endif

   a_data.copy( tmp, box );
}


inline
void Localized::setPointwise( FArrayBox&              a_data,
                              const MagBlockCoordSys& a_coord_sys ) const

{
   Box box( a_data.box() );
   FArrayBox cell_center_coords( box, SpaceDim );
   a_coord_sys.getCellCenteredRealCoords( cell_center_coords );
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
