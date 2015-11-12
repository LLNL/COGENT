#include "ConstantKineticFunction.H"

#include "DataIterator.H"
#include "DisjointBoxLayout.H"
#include "FourthOrderUtil.H"
#include "PhaseGeom.H"

#include "NamespaceHeader.H"

void ConstantKineticFunction::assign( KineticSpecies& a_species,
                                      const Real& a_time ) const
{
   LevelData<FArrayBox>& data( a_species.distributionFunction() );
   const DisjointBoxLayout& grid( data.getBoxes() );
   for (DataIterator pdit( grid.dataIterator() ); pdit.ok(); ++pdit) {
      data[pdit].setVal( m_value );
   }
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   geometry.multBStarParallel( data );
   convertToCellAverage( *geometry.coordSysPtr(), data );
   geometry.multJonValid( data );
}

void ConstantKineticFunction::assign( KineticSpecies& a_species,
                                      const BoundaryBoxLayout& a_bdry_layout,
                                      const Real& a_time ) const
{
   LevelData<FArrayBox>& data( a_species.distributionFunction() );
   const DisjointBoxLayout& grids( data.disjointBoxLayout() );
   // NB: This is a cheat - there's one too many cells at the (dir,side) face
   // of the boundary box, but it doesn't matter because one-sided difference
   // will be used at that face to construct the cell average.  We do want the
   // extra cell in all other directions.
   LevelData<FArrayBox> data_tmp( grids, data.nComp(), IntVect::Unit );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      data_tmp[dit].setVal( m_value );
   }

   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   geometry.multBStarParallel( data_tmp, a_bdry_layout );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      Box domain_box( data_tmp[dit].box() );
      domain_box.growDir( a_bdry_layout.dir(), a_bdry_layout.side(), -1 );
      ProblemDomain domain( domain_box );
      fourthOrderAverageCell( data_tmp[dit], domain, grids[dit] );
   }
   data_tmp.copyTo( data );
   data.exchange();
}

void ConstantKineticFunction::convertToCellAverage(
   const MultiBlockCoordSys&  a_coord_sys,
   LevelData<FArrayBox>&      a_dfn ) const
{
   LevelData<FArrayBox> dfn_tmp(a_dfn.disjointBoxLayout(),
                                a_dfn.nComp(),
                                a_dfn.ghostVect()+IntVect::Unit);

   const DisjointBoxLayout& grids( a_dfn.disjointBoxLayout() );

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      dfn_tmp[dit].copy( a_dfn[dit] );
   }
   dfn_tmp.exchange();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const int block_number( a_coord_sys.whichBlock( grids[dit] ) );
      const PhaseBlockCoordSys* coord_sys
         = dynamic_cast<const PhaseBlockCoordSys*>( a_coord_sys.getCoordSys( block_number ) );

      fourthOrderAverageCell( dfn_tmp[dit], coord_sys->domain(), grids[dit] );
   }
   dfn_tmp.exchange();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_dfn[dit].copy( dfn_tmp[dit] );
   }
}

#include "NamespaceFooter.H"

