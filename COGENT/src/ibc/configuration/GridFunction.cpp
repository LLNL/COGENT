#include "GridFunction.H"
#include "MagGeom.H"
#include "FourthOrderUtil.H"
#include "FourthOrderUtilsClass.H.multidim"

#include "NamespaceHeader.H"

GridFunction::GridFunction(const int& a_verbosity )

   : m_verbosity(a_verbosity)

{
}


void GridFunction::assign(LevelData<FArrayBox>&       a_data,
                          const MultiBlockLevelGeom&  a_geometry,
                          const Real&                 a_time,
                          const bool&                 a_cell_averages ) const
{
   checkGeometryValidity( a_geometry );
   
   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   
   const LevelData<FArrayBox>& real_coords =  ((MagGeom&)a_geometry).getCellCenteredRealCoords();
   const LevelData<FArrayBox>& normalized_flux =  ((MagGeom&)a_geometry).getNormalizedMagneticFluxCell();
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const int block_number( coord_sys.whichBlock( grids[dit] ) );
      
      if (a_cell_averages && !((MagGeom&)a_geometry).secondOrder()) {
         setCellAverages( a_data[dit], a_geometry, real_coords[dit], normalized_flux[dit], block_number );
      }
      else {
         setPointwise( a_data[dit], a_geometry, real_coords[dit], normalized_flux[dit], block_number );
         
      }
   }
   a_data.exchange();
}


void GridFunction::assign(LevelData<FluxBox>&         a_data,
                          const MultiBlockLevelGeom&  a_geometry,
                          const Real&                 a_time,
                          const bool&                 a_face_averages ) const
{
   checkGeometryValidity( a_geometry );

   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   
   const LevelData<FluxBox>& real_coords =  ((MagGeom&)a_geometry).getFaceCenteredRealCoords();
   const LevelData<FluxBox>& normalized_flux =  ((MagGeom&)a_geometry).getNormalizedMagneticFluxFace();
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const int block_number( coord_sys.whichBlock( grids[dit] ) );
      
      if (a_face_averages && !((MagGeom&)a_geometry).secondOrder()) {
         setFaceAverages( a_data[dit], a_geometry, real_coords[dit], normalized_flux[dit], block_number );
      }
      else {
         for (int dir=0; dir<SpaceDim; ++dir) {
            setPointwise( a_data[dit][dir], a_geometry, real_coords[dit][dir], normalized_flux[dit][dir], block_number );
         }
      }
   }
   a_data.exchange();
}


void GridFunction::assign(LevelData<EdgeDataBox>&     a_data,
                          const MultiBlockLevelGeom&  a_geometry,
                          const Real&                 a_time,
                          const bool&                 a_face_averages ) const
{
   checkGeometryValidity( a_geometry );

   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   
   const LevelData<EdgeDataBox>& real_coords =  ((MagGeom&)a_geometry).getEdgeCenteredRealCoords();
   const LevelData<FluxBox>& normalized_flux =  ((MagGeom&)a_geometry).getNormalizedMagneticFluxFace();
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const int block_number( coord_sys.whichBlock( grids[dit] ) );
      //if (a_face_averages && !((MagGeom&)a_geometry).secondOrder()) {
      //   setFaceAverages( a_data[dit], a_geometry, real_coords[dit], normalized_flux[dit], block_number );
      //}
      //else {
         for (int dir=0; dir<SpaceDim; ++dir) {
            setPointwise( a_data[dit][dir], a_geometry, real_coords[dit][dir], normalized_flux[dit][dir], block_number );
         }
      //}
   }
   a_data.exchange();
}


void GridFunction::assign(FArrayBox&                  a_data,
                          const MultiBlockLevelGeom&  a_geometry,
                          const FArrayBox&            a_real_coords,
                          const FArrayBox&            a_normalized_flux,
                          const int                   a_block_number,
                          const Real&                 a_time,
                          const bool&                 a_cell_averages ) const
{
   checkGeometryValidity( a_geometry );
   
   if (a_cell_averages && !((MagGeom&)a_geometry).secondOrder()) {
      setCellAverages( a_data, a_geometry, a_real_coords, a_normalized_flux, a_block_number );
   }
   else {
      setPointwise( a_data, a_geometry, a_real_coords, a_normalized_flux, a_block_number );
   }
}


void GridFunction::assign(FluxBox&                    a_data,
                          const MultiBlockLevelGeom&  a_geometry,
                          const FluxBox&              a_real_coords,
                          const FluxBox&              a_normalized_flux,
                          const int                   a_block_number,
                          const Real&                 a_time,
                          const bool&                 a_face_averages ) const
{
   checkGeometryValidity( a_geometry );
   
   if (a_face_averages && !((MagGeom&)a_geometry).secondOrder()) {
         setFaceAverages( a_data, a_geometry, a_real_coords, a_normalized_flux, a_block_number );
   }
   else {
      for (int dir=0; dir<SpaceDim; ++dir) {
         setPointwise( a_data[dir], a_geometry, a_real_coords[dir], a_normalized_flux[dir], a_block_number );
      }
   }
}


void GridFunction::assign(LevelData<FArrayBox>&       a_data,
                          const MultiBlockLevelGeom&  a_geometry,
                          const BoundaryBoxLayout&    a_bdry_layout,
                          const Real&                 a_time ) const
{
   checkGeometryValidity( a_geometry );
   
   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   const MultiBlockCoordSys& coords( *(a_geometry.coordSysPtr()) );
   
   const LevelData<FArrayBox>& real_coords = ((MagGeom&)a_geometry).getCellCenteredRealCoords();
   const LevelData<FArrayBox>& normalized_flux = ((MagGeom&)a_geometry).getNormalizedMagneticFluxCell();
   
   // NB: This is a cheat - there's one too many cells at the (dir,side) face
   // of the boundary box, but it doesn't matter because one-sided difference
   // will be used at that face to construct the cell average.  We do want the
   // extra cell in all other directions.
   LevelData<FArrayBox> data_tmp( grids, a_data.nComp(), IntVect::Unit );
   for (DataIterator dit( grids ); dit.ok(); ++dit) {
      const Box box( a_bdry_layout.interiorBox( dit ) );
      //const MagBlockCoordSys& coord_sys( ((MagGeom&)a_geometry).getBlockCoordSys( box ) );
      const int block_number( coords.whichBlock( box ) );
      const DataIndex& internal_dit( a_bdry_layout.dataIndex( dit ) );
      setPointwise( data_tmp[dit], a_geometry, real_coords[internal_dit], normalized_flux[internal_dit], block_number );
   }
   
   if (!((MagGeom&)a_geometry).secondOrder()) {
      FourthOrderUtil FourthOrderOperators; //Object that holds various fourth-order operatiosns 
      const bool useSG = ((MagGeom&)a_geometry).useSG(); // decide whether to use SG based on MagGeom
      FourthOrderOperators.setSG(useSG); //Whether to use the SG versions of fourth order stencils
      for (DataIterator dit( grids ); dit.ok(); ++dit) {
         Box domain_box( data_tmp[dit].box() );
         domain_box.growDir( a_bdry_layout.dir(), a_bdry_layout.side(), -1 );
         ProblemDomain domain( domain_box );
         fourthOrderAverageCell( data_tmp[dit], domain, grids[dit] );
         FourthOrderOperators.fourthOrderAverageCellGen( data_tmp[dit], domain, grids[dit] );
      }
   }
   data_tmp.copyTo( a_data );
   a_data.exchange();
}

void GridFunction::setCellAverages(FArrayBox&                  a_data,
                                   const MultiBlockLevelGeom&  a_geometry,
                                   const FArrayBox&            a_real_coords,
                                   const FArrayBox&            a_normalized_flux,
                                   const int                   a_block_number ) const
{
   Box box( a_data.box() );
   Box tmp_box( box );
   tmp_box.grow( IntVect::Unit );
   FArrayBox tmp( tmp_box, a_data.nComp() );
   
   setPointwise( tmp, a_geometry, a_real_coords, a_normalized_flux, a_block_number );
   
   const MagBlockCoordSys& coord_sys = getCoordSys(a_geometry, a_block_number);
   FourthOrderUtil FourthOrderOperators; //Object that holds various fourth-order operatiosns
   //const bool useSG = ((MagGeom&)a_geometry).useSG(); // decide whether to use SG stenciles based on MagGeom 
   FourthOrderOperators.setSG(false);
   //fourthOrderAverageCell( tmp, coord_sys.domain(), box );
   FourthOrderOperators.fourthOrderAverageCellGen( tmp, coord_sys.domain(), box );
   
   a_data.copy( tmp, box );
}

void GridFunction::setFaceAverages(FluxBox&                    a_data,
                                   const MultiBlockLevelGeom&  a_geometry,
                                   const FluxBox&              a_real_coords,
                                   const FluxBox&              a_normalized_flux,
                                   const int                   a_block_number ) const
{
   const std::string msg( "GridFunction: face averaging is not yet implemented ");
   MayDay::Error( msg.c_str() );
}

const MagBlockCoordSys&
GridFunction::getCoordSys(const MultiBlockLevelGeom& a_geometry,
                          const Box&                 a_box ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   const int block_number( coord_sys.whichBlock( a_box ) );
   const NewCoordSys* block_coord_sys( coord_sys.getCoordSys( block_number ) );
   return static_cast<const MagBlockCoordSys&>( *block_coord_sys );
}

const MagBlockCoordSys&
GridFunction::getCoordSys(const MultiBlockLevelGeom& a_geometry,
                          const int                  a_block_number ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   const NewCoordSys* block_coord_sys( coord_sys.getCoordSys( a_block_number ) );
   return static_cast<const MagBlockCoordSys&>( *block_coord_sys );
}


#include "NamespaceFooter.H"
