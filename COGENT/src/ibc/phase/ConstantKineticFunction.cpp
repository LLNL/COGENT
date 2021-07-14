#include "ConstantKineticFunction.H"

#include "DataIterator.H"
#include "DisjointBoxLayout.H"
#include "FourthOrderUtil.H"
#include "PhaseGeom.H"
#include "KineticFunctionUtils.H"

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
   if ( !(geometry.secondOrder()) )  {
     KineticFunctionUtils::convertToCellAverage( geometry, data, m_useSG );
   }
   geometry.multJonValid( data );
   data.exchange();
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
   if ( !(geometry.secondOrder()) )  {
     FourthOrderUtil FourthOrderOperators; //Object that holds various fourth-order operatiosns 
     FourthOrderOperators.setSG(m_useSG); //Whether to use the SG versions of fourth order stencils
     for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
       Box domain_box( data_tmp[dit].box() );
       domain_box.growDir( a_bdry_layout.dir(), a_bdry_layout.side(), -1 );
       ProblemDomain domain( domain_box );
       FourthOrderOperators.fourthOrderAverageCellGen(data_tmp[dit], domain, grids[dit]);
     }
   }
   data_tmp.copyTo( data );
   data.exchange();
}

#include "NamespaceFooter.H"

