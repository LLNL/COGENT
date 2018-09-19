#include "TensorProductSineKineticFunction.H"

#include <iostream>
#include <typeinfo>

#include "DataIterator.H"
#include "Directions.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "FourthOrderUtil.H"
#include "LevelData.H"
#include "TensorProductSineF_F.H"
#include "MayDay.H"
#include "LogRectPhaseCoordSys.H"
#include "SingleNullPhaseCoordSys.H"
#include "SNCorePhaseCoordSys.H"
#include "MultiBlockCoordSys.H"
#include "PhaseBlockCoordSys.H"
#include "Vector.H"
#include "KineticFunctionUtils.H"

#include "NamespaceHeader.H"


TensorProductSineKineticFunction::TensorProductSineKineticFunction(
   ParmParse& a_pp,
   const int& a_verbosity )
   : m_verbosity(a_verbosity),
     m_amplitude(1.0)
{
   parseParameters( a_pp );
}


void TensorProductSineKineticFunction::assign( KineticSpecies& a_species,
                                               const Real& a_time ) const
{
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   checkGeometryValidity( geometry );

   LevelData<FArrayBox>& dfn( a_species.distributionFunction() );
   const DisjointBoxLayout& grids( dfn.disjointBoxLayout() );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      setPointValues( dfn[dit],
                      geometry.getBlockCoordSys( grids[dit] ),
                      a_time  );
   }
//   geometry.multBStarParallel( dfn );
//   geometry.multPointwiseJ( dfn, a_bdry_layout );
   KineticFunctionUtils::convertToCellAverage( geometry, dfn );
   geometry.multJonValid( dfn );
   dfn.exchange();
}


void TensorProductSineKineticFunction::assign( KineticSpecies& a_species,
                                               const BoundaryBoxLayout& a_bdry_layout,
                                               const Real& a_time ) const
{
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   checkGeometryValidity( geometry );

   LevelData<FArrayBox>& dfn( a_species.distributionFunction() );
   const DisjointBoxLayout& grids( dfn.disjointBoxLayout() );
   // NB: This is a cheat - there's one too many cells at the (dir,side) face
   // of the boundary box, but it doesn't matter because one-sided difference
   // will be used at that face to construct the cell average.  We do want the
   // extra cell in all other directions.
   LevelData<FArrayBox> dfn_tmp( grids, dfn.nComp(), IntVect::Unit );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const Box box( a_bdry_layout.interiorBox( dit ) );
      const PhaseBlockCoordSys& coord_sys( geometry.getBlockCoordSys( box ) );
      setPointValues( dfn_tmp[dit], coord_sys, a_time );
   }
//   geometry.multBStarParallel( dfn_tmp, a_bdry_layout );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      Box domain_box( dfn_tmp[dit].box() );
      domain_box.growDir( a_bdry_layout.dir(), a_bdry_layout.side(), -1 );
      ProblemDomain domain( domain_box );
      fourthOrderAverageCell( dfn_tmp[dit], domain, grids[dit] );
   }
   dfn_tmp.copyTo( dfn );
   dfn.exchange();
}


inline
void TensorProductSineKineticFunction::parseParameters( ParmParse& a_pp )
{
   a_pp.get( "amplitude", m_amplitude );

   Vector<Real> temp( PDIM );

   temp.assign( 0.0 );
   a_pp.queryarr( "wavevector", temp, 0, PDIM );
   m_wavevector = RealVect( temp );

   if (m_verbosity) {
      printParameters();
   }
}


inline
void TensorProductSineKineticFunction::checkGeometryValidity( const PhaseGeom& a_geometry ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   bool unknown_geom( typeid(coord_sys) != typeid(LogRectPhaseCoordSys) );
   unknown_geom &= (typeid(coord_sys) != typeid(SingleNullPhaseCoordSys));
   unknown_geom &= (typeid(coord_sys) != typeid(SNCorePhaseCoordSys));
   
   if ( unknown_geom ) {
      const std::string msg( "TensorProductSineKineticFunction: Attempt to use unknown geometry. ");
      MayDay::Error( msg.c_str() );
   }
}


inline
void TensorProductSineKineticFunction::setPointValues( FArrayBox&                a_dfn,
                                         const PhaseBlockCoordSys& a_coord_sys,
                                         const Real&               a_time) const
{
   const Box& box( a_dfn.box() );
   FArrayBox cell_center_coords( box, PDIM );
   a_coord_sys.getCellCenteredRealCoords( cell_center_coords );
   FORT_SET_TENSORPRODUCTSINE( CHF_FRA(a_dfn),
                               CHF_BOX(box),
                               CHF_CONST_FRA(cell_center_coords),
                               CHF_CONST_REAL(a_time),
                               CHF_CONST_REAL(m_amplitude),
                               CHF_CONST_REALVECT(m_wavevector) );
}


inline
void TensorProductSineKineticFunction::printParameters() const
{
   if (procID()==0) {
      std::cout << "TensorProductSine kinetic function parameters:" << std::endl;
      std::cout << "  amplitude:  " << m_amplitude  << std::endl;
      std::cout << "  wavevector: " << m_wavevector << std::endl;
      std::cout << std::endl;
   }
}

#include "NamespaceFooter.H"
