#include "Cosine.H"

#include "CosineF_F.H"

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
#include "MillerCoordSys.H"
#include "MillerBlockCoordSys.H"
#include "SlabCoordSys.H"
#include "SlabBlockCoordSys.H"
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


Cosine::Cosine( //const std::string& a_name,
                    ParmParse& a_pp,
                    const int& a_verbosity )
   : //m_name(a_name),
     m_verbosity(a_verbosity),
     m_constant(0.0),
     m_amplitude(0.0)
{
   parseParameters( a_pp );
}


void Cosine::assign( LevelData<FArrayBox>& a_data,
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


void Cosine::assign( FArrayBox& a_data,
                     const MultiBlockLevelGeom& a_geometry,
                     const Box& a_box,
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


void Cosine::assign( LevelData<FArrayBox>& a_data,
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


inline
void Cosine::checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   bool not_annular( typeid(coord_sys) != typeid(MillerCoordSys) );
   not_annular &= (typeid(coord_sys) != typeid(MillerBlockCoordSys));
   not_annular &= (typeid(coord_sys) != typeid(SlabCoordSys));
   not_annular &= (typeid(coord_sys) != typeid(SlabBlockCoordSys));
   if ( not_annular ) {
      const std::string msg( "Cosine: Attempt to use non-annular geometry. ");
      MayDay::Error( msg.c_str() );
   }
}


inline
void Cosine::setCellAverages( FArrayBox&              a_data,
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
void Cosine::setPointwise( FArrayBox&              a_data,
                             const MagBlockCoordSys& a_coord_sys ) const

{
   Box box( a_data.box() );
   FArrayBox cell_center_coords( box, SpaceDim );
   a_coord_sys.getCellCenteredFluxCoords( cell_center_coords );

   RealVect k( m_mode );
   k[RADIAL_DIR] *= 2.0 * Constants::PI;
   k[RADIAL_DIR] /= ( a_coord_sys.getOuterFluxLabel() - a_coord_sys.getInnerFluxLabel() );

   cell_center_coords.plus( -0.5*( a_coord_sys.getOuterFluxLabel() + a_coord_sys.getInnerFluxLabel() ),
                            RADIAL_DIR, 1 );

   FORT_SET_COSINE( CHF_FRA(a_data),
                    CHF_BOX(box),
                    CHF_CONST_FRA(cell_center_coords),
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
