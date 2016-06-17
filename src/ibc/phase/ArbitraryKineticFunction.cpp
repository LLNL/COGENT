#include "ArbitraryKineticFunction.H"

#include <iostream>
#include <typeinfo>

#include "DataIterator.H"
#include "Directions.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "FourthOrderUtil.H"
#include "LevelData.H"
#include "LocalizedF_F.H"
#include "MayDay.H"
#include "MillerPhaseCoordSys.H"
#include "SlabPhaseCoordSys.H"
#include "MultiBlockCoordSys.H"
#include "PhaseBlockCoordSys.H"
#include "SingleNullPhaseCoordSys.H"
#include "Vector.H"
#include "KineticFunctionUtils.H"

#include "NamespaceHeader.H"


ArbitraryKineticFunction::ArbitraryKineticFunction( ParmParse& a_pp,
                                                    const int& a_verbosity )
   : m_verbosity(a_verbosity),
     m_function("UNDEFINED")
{
   parseParameters( a_pp );

   const char *userFormular=m_function.c_str();
   
   m_pscore = new ParsingCore(userFormular);

#if DEBUG_5D==1///DEBUG_5D
      cout << " ArbitraryKineticFunction::ArbitraryKineticFunction() m_pscore->getFormula():" << m_pscore->getFormula() << std::endl;
      cout << " ArbitraryKineticFunction::ArbitraryKineticFunction() m_pscore->getManipStr():" << m_pscore->getManipStr()<< std::endl;
      cout << " ArbitraryKineticFunction::ArbitraryKineticFunction() m_pscore->getPostStr():" << m_pscore->getPostStr()<< std::endl;
#endif//DEBUG_5D
   
}


void ArbitraryKineticFunction::assign( KineticSpecies& a_species,
                                       const Real& a_time ) const
{
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
//   checkGeometryValidity( geometry );

   LevelData<FArrayBox>& dfn( a_species.distributionFunction() );
   const DisjointBoxLayout& grids( dfn.disjointBoxLayout() );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      setPointValues( dfn[dit],
                      geometry.getBlockCoordSys( grids[dit] ),
                      a_time  );
   }
   geometry.multBStarParallel( dfn );
   KineticFunctionUtils::convertToCellAverage( *geometry.coordSysPtr(), dfn );
   geometry.multJonValid( dfn );
   dfn.exchange();
}


void ArbitraryKineticFunction::assign( KineticSpecies& a_species,
                                       const BoundaryBoxLayout& a_bdry_layout,
                                       const Real& a_time ) const
{
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
//   checkGeometryValidity( geometry );

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
   geometry.multBStarParallel( dfn_tmp, a_bdry_layout );
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
void ArbitraryKineticFunction::parseParameters( ParmParse& a_pp )
{
   a_pp.get( "function", m_function );

   if (m_verbosity) {
      printParameters();
   }
}


inline
void ArbitraryKineticFunction::checkGeometryValidity( const PhaseGeom& a_geometry ) const
{
   const MultiBlockCoordSys& an_coord_sys( *(a_geometry.coordSysPtr()) );
   bool not_annular( typeid(an_coord_sys) != typeid(MillerPhaseCoordSys) );
   not_annular &= (typeid(an_coord_sys) != typeid(SlabPhaseCoordSys));

   const MultiBlockCoordSys& sn_coord_sys( *(a_geometry.coordSysPtr()) );
   bool not_single_null( typeid(sn_coord_sys) != typeid(SingleNullPhaseCoordSys) );

   if ( not_annular && not_single_null ) {
      const std::string msg( "ArbitraryKineticFunction: Attempt to use unknown geometry. ");
      MayDay::Error( msg.c_str() );
   }
}


inline
void ArbitraryKineticFunction::setPointValues( FArrayBox&                a_dfn,
                                         const PhaseBlockCoordSys& a_coord_sys,
                                         const Real&               a_time) const
{
   const Box& box( a_dfn.box() );
   FArrayBox cell_center_coords( box, PDIM );
   a_coord_sys.getCellCenteredRealCoords( cell_center_coords );

   RealVect a_amrDx = a_coord_sys.dx();
// rescale a_amrDx so that x, y, z, vpar and mu span from 0 to 2*pi
   IntVect hi_index = a_coord_sys.domain().domainBox().bigEnd();
   hi_index = hi_index+1;
   a_amrDx = a_amrDx*2.0*M_PI/(a_amrDx*hi_index); 

   a_dfn.setVal(0.0);

   BoxIterator bit(a_dfn.box());
   for (bit.begin(); bit.ok(); ++bit)
   {
       IntVect iv = bit();
       RealVect loc(iv);
       loc *= a_amrDx;

#if PDIM==5
       Real val = m_pscore->calc5d(loc[0],loc[1],loc[2],loc[3],loc[4]);
#else
       Real val = m_pscore->calc4d(loc[0],loc[1],loc[2],loc[3]);
#endif
       a_dfn(iv,0) += val;

   }


}


inline
void ArbitraryKineticFunction::printParameters() const
{
   if (procID()==0) {
      std::cout << "Arbitrary kinetic function parameters:" << std::endl;
      std::cout << "  function: "  << m_function << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
   }
}

#include "NamespaceFooter.H"
