#include "Tanh.H"

#include "TanhF_F.H"

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
#include "LogRectCoordSys.H"
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


Tanh::Tanh( //const std::string& a_name,
                ParmParse& a_pp,
                const int& a_verbosity )
   : //m_name(a_name),
     m_verbosity(a_verbosity),
     m_inner_radial_value(0.0),
     m_outer_radial_value(0.0),
     m_radial_midpoint(0.0),
     m_radial_width(1.0)
{
   parseParameters( a_pp );
}


void Tanh::assign( LevelData<FArrayBox>& a_data,
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


void Tanh::assign( LevelData<FluxBox>& a_data,
                   const MultiBlockLevelGeom& a_geometry,
                   const Real& a_time,
                   const bool& a_cell_averages ) const
{
   checkGeometryValidity( a_geometry );

   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      if (a_cell_averages) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            setCellAverages( a_data[dit][dir], getCoordSys( a_geometry, grids[dit] ) );
         }
      }
      else {
         for (int dir=0; dir<SpaceDim; ++dir) {
            setPointwise( a_data[dit][dir], getCoordSys( a_geometry, grids[dit] ) );
         }
      }
   }
   a_data.exchange();
}


void Tanh::assign( FArrayBox& a_data,
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


void Tanh::assign( FluxBox& a_data,
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


void Tanh::assign( LevelData<FArrayBox>& a_data,
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


inline
void Tanh::checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   bool not_single_block( typeid(coord_sys) != typeid(LogRectCoordSys) );
   if ( not_single_block ) {
      const std::string msg( "Tanh: Attempt to use non-single-block geometry. ");
      MayDay::Error( msg.c_str() );
   }
}

inline
void Tanh::setCellAverages( FArrayBox&              a_data,
                            const MagBlockCoordSys& a_coord_sys ) const

{
   Box box( a_data.box() );
   Box tmp_box( box );
   tmp_box.grow( IntVect::Unit );
   FArrayBox tmp( tmp_box, a_data.nComp() );

   setPointwise( tmp, a_coord_sys );

   fourthOrderAverageCell( tmp, a_coord_sys.domain(), box );

   a_data.copy( tmp, box );
}


inline
void Tanh::setPointwise( FArrayBox&              a_data,
                         const MagBlockCoordSys& a_coord_sys ) const

{
   Box box( a_data.box() );
   FArrayBox cell_center_coords( box, SpaceDim );
   a_coord_sys.getCellCenteredMappedCoords( cell_center_coords );
   RealVect dx = a_coord_sys.getMappedCellSize();

   Real inner_radial_coord( a_coord_sys.lowerMappedCoordinate(RADIAL_DIR) );
   Real outer_radial_coord( a_coord_sys.upperMappedCoordinate(RADIAL_DIR) );
   
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
