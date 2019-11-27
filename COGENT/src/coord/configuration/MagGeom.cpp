#include <array>
#include "MagGeom.H"
#include "MagGeomF_F.H"
#include "Poisson.H"
#include "Directions.H"
#include "CONSTANTS.H"
#include "FourthOrderUtil.H"
#include "newMappedGridIO.H"
#include "SimpleDivergence.H"
#include "EdgeToCell.H"
#include "SingleNullBlockCoordSys.H"
#include "SingleNullCoordSys.H"

#include "CoDimCopyManager.H"

#define NEW_EXTRABLOCK_EXCHANGE

#include "NamespaceHeader.H"

#undef TEST_DIVERGENCE_CLEANING

inline array<double,3> computeUnitVector( const array<double,3>& a_X )
{
   array<double,3> result;
   
   double magnitude = 0.;
   for (int n=0; n<3; ++n) {
      magnitude += a_X[n]*a_X[n];
   }
   magnitude = sqrt(magnitude);
   
   for (int n=0; n<3; ++n) {
      result[n] = a_X[n]/magnitude;
   }
   
   return result;
}

inline Vector<Real> matrixMultiplication(Vector<Real>& a, Vector<Real>& b)
{
   
   CH_assert(a.size()==SpaceDim*SpaceDim);
   CH_assert(b.size()==SpaceDim*SpaceDim);
   
   Vector<Real> result(SpaceDim * SpaceDim,0);
   
   for (int i=0; i<SpaceDim; ++i) {
      for (int j=0; j<SpaceDim; ++j) {
         for (int n=0; n<SpaceDim; ++n) {
            result[SpaceDim * i + j] += a[SpaceDim*i + n] * b[SpaceDim*n +j];
         }
      }
   }

   return result;
}

MagGeom::MagGeom(ParmParse&                         a_pp,
                 const RefCountedPtr<MagCoordSys>&  a_coord_sys,
                 const DisjointBoxLayout&           a_grids,
                 int                                a_ghosts)
   : MultiBlockLevelGeom(a_coord_sys, a_grids, a_ghosts),
     m_coord_sys(a_coord_sys),
     m_velocity_data_cached(false),
     m_correct_field(false),
     m_pp_prefix(a_pp.prefix())
{
   if (a_pp.contains("correct_field")) {
     a_pp.get("correct_field", m_correct_field);
   }

   m_extrablock_exchange = false;
   a_pp.query( "extrablock_exchange", m_extrablock_exchange);
   
   m_dealignment_corrections = false;
   a_pp.query( "dealignment_corrections", m_dealignment_corrections);
   
   m_model_flux_function = false;
   a_pp.query( "model_flux_function", m_model_flux_function);

   m_sheared_mb_geom = false;
   a_pp.query( "sheared_multiblock_geometry", m_sheared_mb_geom);

   m_second_order = false;
   a_pp.query( "second_order", m_second_order);

#if CFG_DIM ==3
   m_mb_dir = TOROIDAL_DIR;
   if (m_sheared_mb_geom) {
      initializeShearedMBGeom(a_grids);
   }
#endif

   // Check for full periodicity
   m_fully_periodic = true;
   for (int block=0; block<a_coord_sys->numBlocks(); ++block) {
     const ProblemDomain& domain = getBlockCoordSys(block).domain();
     for (int dir=0; dir<SpaceDim; dir++) {
       m_fully_periodic = m_fully_periodic && domain.isPeriodic(dir);
     }
   }
   
   if ( (typeid(*m_coord_sys) == typeid(SingleNullCoordSys)) && (m_extrablock_exchange) ) {
      int nghost = 1;
      //      ((SingleNullCoordSys*)m_coord_sys)->defineStencilsUe(a_grids, nghost);
      MagCoordSys* const coord_ptr = m_coord_sys.getRefToThePointer();
      ((SingleNullCoordSys* const)coord_ptr)->defineStencilsUe(a_grids, nghost);
   }

   // Compute and cache the physical cell centers
   setRealCoords();
   
   // Compute and cache the normalized magnetic flux
   setNormalizedMagneticFlux();

   // Precompute the face-centered field data
   IntVect field_ghosts = a_ghosts*IntVect::Unit;
   m_BField_fc.define(a_grids, 3, field_ghosts);
   m_BFieldMag_fc.define(a_grids, 1, field_ghosts);
   m_BFieldDir_fc.define(a_grids, 3, field_ghosts);
   m_gradBFieldMag_fc.define(a_grids, 3, field_ghosts);
   m_curlBFieldDir_fc.define(a_grids, 3, field_ghosts);
   m_BFieldDirdotcurlBFieldDir_fc.define(a_grids, 1, field_ghosts);

   computeFieldData( m_BField_fc,
                     m_BFieldMag_fc,
                     m_BFieldDir_fc,
                     m_gradBFieldMag_fc,
                     m_curlBFieldDir_fc,
                     m_BFieldDirdotcurlBFieldDir_fc );

   m_BField_cc.define(a_grids, 3, field_ghosts);
   m_BFieldMag_cc.define(a_grids, 1, field_ghosts);
   m_BFieldDir_cc.define(a_grids, 3, field_ghosts);
   m_gradBFieldMag_cc.define(a_grids, 3, field_ghosts);
   m_curlBFieldDir_cc.define(a_grids, 3, field_ghosts);
   m_BFieldDirdotcurlBFieldDir_cc.define(a_grids, 1, field_ghosts);

   // Precompute the cell-centered field data.  If the field has been corrected
   // to clean its divergence (in the preceding call of computeFieldData), then
   // we need to average the corrected values to cell centers (which is only
   // done to second order by EdgeToCell(), unfortunately).  Otherwise, we compute
   // the actual cell-centered data.
   if ( m_correct_field ) {
      cellCenter( m_BField_fc, m_BField_cc );
      cellCenter( m_BFieldMag_fc, m_BFieldMag_cc );
      cellCenter( m_BFieldDir_fc, m_BFieldDir_cc );
      cellCenter( m_gradBFieldMag_fc, m_gradBFieldMag_cc );
      cellCenter( m_curlBFieldDir_fc, m_curlBFieldDir_cc );
      cellCenter( m_BFieldDirdotcurlBFieldDir_fc, m_BFieldDirdotcurlBFieldDir_cc );
   }
   else {
      computeFieldData( m_BField_cc,
                        m_BFieldMag_cc,
                        m_BFieldDir_cc,
                        m_gradBFieldMag_cc,
                        m_curlBFieldDir_cc,
                        m_BFieldDirdotcurlBFieldDir_cc );
   }

   // Check the mapping consistency at interblock interfaces
   checkMultiblockMappingConsistency();

   // Run diagnostics
   plotInitializationData(a_pp, a_grids, 0.0);

   //Create multiblock exchange objects
   if ( (m_coord_sys->numBlocks() > 1) && !m_extrablock_exchange && !m_sheared_mb_geom) {
      m_mblexPtr = new MultiBlockLevelExchangeAverage();

      int spaceOrder = m_second_order? 2: 4;
      m_mblexPtr->define( this, m_ghosts, spaceOrder );

      m_exchange_transverse_block_register = new BlockRegister(m_coord_sys, a_grids, 1);
   }
   else {
      m_mblexPtr = NULL;
      m_exchange_transverse_block_register = NULL;
   }
    
    // Initializations for dealignment-correction utilities (for SN geometry only)
    if ( (typeid(*m_coord_sys) == typeid(SingleNullCoordSys)) && m_dealignment_corrections ) {
        m_magFS_mapping_cell.define(a_grids, 3, IntVect::Zero);
        m_magFS_mapping_face.define(a_grids, 3, IntVect::Zero);
        computeMagFluxMappingCell( m_magFS_mapping_cell);
        computeMagFluxMappingFace( m_magFS_mapping_face);
    }

   // Initialize geometry coeffients used in elliptic field solvers
   m_perp_coeff.define(a_grids, SpaceDim*SpaceDim, field_ghosts);
   m_par_coeff.define(a_grids, SpaceDim*SpaceDim , field_ghosts);
   computeEllipticOpCoefficients(m_perp_coeff, m_par_coeff);
   
   m_perp_coeff_mapped.define(a_grids, SpaceDim*SpaceDim, field_ghosts);
   m_par_coeff_mapped.define(a_grids, SpaceDim*SpaceDim, field_ghosts);
   computeEllipticOpCoefficientsMapped(m_perp_coeff_mapped, m_par_coeff_mapped );

   // Initialize exchangeCopier with 2 ghost layers
   //const IntVect ghostVect = 2*IntVect::Unit;
   //m_exchangeCopier2nd.define(m_gridsFull, m_gridsFull, ghostVect, true);
   //m_exchangeCopier.trimEdges(m_gridsFull, ghostVect);

   //
   ///////////////////////////////////////////////////////////


   double mapping_error = maxMappingError();
   if (procID()==0) {
      cout << "Mapping error = " << mapping_error << endl;
   }

   if (mapping_error > 0.01) {
      MayDay::Error("MagGeom::MagGeom(): Difference between forward and inverse mappings is greater than one percent");
   }
}


MagGeom::~MagGeom()
{
   if (m_exchange_transverse_block_register) delete m_exchange_transverse_block_register;
   if (m_mblexPtr) delete m_mblexPtr;
}

void
MagGeom::setMetricTerms( const DisjointBoxLayout&  a_grids,
                         const IntVect&            a_ghostVect ) const
{

   if ( !m_metrics.isDefined() || (m_metrics.ghostVect() < a_ghostVect) ) {
      
      m_metrics.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);
      
      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);
         if (!coord_sys.isPointwiseMetrics()) {
            coord_sys.getN(m_metrics[dit], m_metrics[dit].box());
         }
         else {
            coord_sys.getPointwiseMetrics(m_metrics[dit]);
         }
      }
      
      m_coord_sys->postProcessMetricData(m_metrics);
   }
   
}

void
MagGeom::getMetricTerms( LevelData<FluxBox>& a_N ) const
{

   const DisjointBoxLayout& grids = a_N.disjointBoxLayout();
   const IntVect& ghosts = a_N.ghostVect();
   setMetricTerms(grids, ghosts);

   DataIterator dit = a_N.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_N[dit].copy(m_metrics[dit]);
   }
   
}

void
MagGeom::getMetricTerms( LevelData<FluxBox>& a_N,
                         LevelData<FluxBox>& a_tanGradN ) const
{
   CH_assert(a_N.getBoxes()==a_tanGradN.getBoxes());

   const DisjointBoxLayout& grids = a_N.disjointBoxLayout();

   for (DataIterator dit(a_N.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
      
      if (!coord_sys.isPointwiseMetrics()) {
         coord_sys.getN(a_N[dit], a_N[dit].box());
      }
      else {
         coord_sys.getPointwiseMetrics(a_N[dit]);
      }
   }

   for (DataIterator dit(a_tanGradN.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
      coord_sys.computeTangentialGrad(a_tanGradN[dit], a_N[dit], a_tanGradN[dit].box());
   }

   m_coord_sys->postProcessMetricData(a_N);
   m_coord_sys->postProcessMetricData(a_tanGradN);
}

void
MagGeom::setCellVolumes( const DisjointBoxLayout&  a_grids,
                         const IntVect&            a_ghostVect ) const
{
   
   if ( !m_cell_volume.isDefined() || (m_cell_volume.ghostVect() < a_ghostVect) ) {

      m_cell_volume.define(a_grids, 1, a_ghostVect);

      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);

         Box grown_box(grow(m_cell_volume[dit].box(),1));
         FluxBox N(grown_box, coord_sys.getNumN());
	 // Chombo's getN can be expensive, especially
	 // if dxdxi is arithmetically intensive 
	 // (as in field-aliened toroidal geom)
	 // thus we only compute N if needed. 
	 // For pointwiseMetrics and sub-grid methods cellVol 
	 // is not using N, it's just passed for consistency
         if (!coord_sys.isPointwiseMetrics() && !coord_sys.isSubGridGeom()) {
	   coord_sys.getN(N, grown_box);
         }
         coord_sys.cellVol(m_cell_volume[dit], N, m_cell_volume[dit].box());
      }

      m_coord_sys->postProcessMetricData(m_cell_volume);

      m_cell_volume.exchange();

      if (m_extrablock_exchange) fillInternalGhosts(m_cell_volume);
   }

}

void
MagGeom::getCellVolumes( LevelData<FArrayBox>& a_volume ) const
{

   const DisjointBoxLayout& grids = a_volume.disjointBoxLayout();
   const IntVect& ghosts = a_volume.ghostVect();
   setCellVolumes(grids, ghosts);

   DataIterator dit = a_volume.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_volume[dit].copy(m_cell_volume[dit]);
   }

}

void
MagGeom::setJ( const DisjointBoxLayout&  a_grids,
               const IntVect&            a_ghostVect ) const
{
   if ( !m_J.isDefined() || (m_J.ghostVect() < a_ghostVect) ) {

      m_J.define(a_grids, 1, a_ghostVect);

      setCellVolumes(a_grids, a_ghostVect);

      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);
         m_J[dit].copy(m_cell_volume[dit]);
         double fac( 1.0 / coord_sys.getMappedCellVolume() );
         m_J[dit].mult(fac);
      }

      m_J.exchange();
      
      if (m_extrablock_exchange) fillInternalGhosts(m_J);
   }

}

void
MagGeom::getJ( LevelData<FArrayBox>& a_J ) const
{
   const DisjointBoxLayout& grids = a_J.disjointBoxLayout();
   const IntVect& ghosts = a_J.ghostVect();
   setJ(grids, ghosts);

   for (DataIterator dit(a_J.dataIterator()); dit.ok(); ++dit) {
      a_J[dit].copy(m_J[dit]);
   }
}

void
MagGeom::setJonFaces( const DisjointBoxLayout&  a_grids,
                      const IntVect&            a_ghostVect ) const
{
   if ( !m_JonFaces.isDefined() || (m_JonFaces.ghostVect() < a_ghostVect) ) {

      m_JonFaces.define(a_grids, 1, a_ghostVect);

      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);
         coord_sys.pointwiseJ(m_JonFaces[dit]);
      }

      m_JonFaces.exchange();
      
      //if (m_extrablock_exchange) fillInternalGhosts(m_JonFaces);
   }

}

void
MagGeom::getJ( LevelData<FluxBox>& a_J ) const
{
   const DisjointBoxLayout& grids = a_J.disjointBoxLayout();
   const IntVect& ghosts = a_J.ghostVect();
   setJonFaces(grids, ghosts);

   for (DataIterator dit(a_J.dataIterator()); dit.ok(); ++dit) {
      a_J[dit].copy(m_JonFaces[dit]);
   }
}

void
MagGeom::setJonEdges( const DisjointBoxLayout&  a_grids,
                      const IntVect&            a_ghostVect ) const
{
   if ( !m_JonEdges.isDefined() || (m_JonEdges.ghostVect() < a_ghostVect) ) {

      m_JonEdges.define(a_grids, 1, a_ghostVect);

      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);
         coord_sys.pointwiseJ(m_JonEdges[dit]);
      }

      m_JonEdges.exchange();
      
      //if (m_extrablock_exchange) fillInternalGhosts(m_JonEdges);
   }

}

void
MagGeom::getJ( LevelData<EdgeDataBox>& a_J ) const
{
   const DisjointBoxLayout& grids = a_J.disjointBoxLayout();
   const IntVect& ghosts = a_J.ghostVect();
   setJonEdges(grids, ghosts);

   for (DataIterator dit(a_J.dataIterator()); dit.ok(); ++dit) {
      a_J[dit].copy(m_JonEdges[dit],0,0,1);
   }
}

void
MagGeom::setN( const DisjointBoxLayout&  a_grids, 
               const IntVect&            a_ghostVect ) const
{
   if ( !m_N.isDefined() || (m_N.ghostVect() < a_ghostVect) ) {

      m_N.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);

      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);
         coord_sys.getPointwiseN(m_N[dit]);
      }

      m_N.exchange();
      
      if (m_extrablock_exchange) fillInternalGhosts(m_N);
   }

}

void
MagGeom::getN( LevelData<FArrayBox>& a_N ) const
{
   const DisjointBoxLayout& grids = a_N.disjointBoxLayout();
   const IntVect& ghosts = a_N.ghostVect();
   setN(grids, ghosts);

   for (DataIterator dit(a_N.dataIterator()); dit.ok(); ++dit) {
      a_N[dit].copy(m_N[dit]);
   }
}


void
MagGeom::setdXdxiOnCells( const DisjointBoxLayout&  a_grids,
                          const IntVect&            a_ghostVect ) const
{

   if ( !m_dXdxi.isDefined() || (m_dXdxi.ghostVect() < a_ghostVect) ) {

      m_dXdxi.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);

      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);
         coord_sys.getCellCentereddXdxi(m_dXdxi[dit]);
      }

      m_dXdxi.exchange();
      
      if (m_extrablock_exchange) fillInternalGhosts(m_dXdxi);
   }

}

void
MagGeom::getdXdxi( LevelData<FArrayBox>& a_dXdxi ) const
{
   CH_assert(a_dXdxi.nComp() == SpaceDim*SpaceDim);
   
   const DisjointBoxLayout& grids = a_dXdxi.disjointBoxLayout();
   const IntVect& ghosts = a_dXdxi.ghostVect();
   setdXdxiOnCells(grids, ghosts);  

   for (DataIterator dit(a_dXdxi.dataIterator()); dit.ok(); ++dit) {
      a_dXdxi[dit].copy(m_dXdxi[dit]);
   }
}

void
MagGeom::setdXdxiOnFaces( const DisjointBoxLayout&  a_grids,
                          const IntVect&            a_ghostVect ) const
{

   if ( !m_dXdxi_onFaces.isDefined() || (m_dXdxi_onFaces.ghostVect() < a_ghostVect) ) {

      m_dXdxi_onFaces.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);

      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);
         coord_sys.getFaceCentereddXdxi(m_dXdxi_onFaces[dit]);
      }

      m_dXdxi_onFaces.exchange();
      
      //if (m_extrablock_exchange) fillInternalGhosts(m_dXdxi_onFaces);
   }

}

void
MagGeom::getdXdxi( LevelData<FluxBox>& a_dXdxi ) const
{
   CH_assert(a_dXdxi.nComp() == SpaceDim*SpaceDim);

   const DisjointBoxLayout& grids = a_dXdxi.disjointBoxLayout();
   const IntVect& ghosts = a_dXdxi.ghostVect();
   setdXdxiOnFaces(grids, ghosts); 

   for (DataIterator dit(a_dXdxi.dataIterator()); dit.ok(); ++dit) {
      a_dXdxi[dit].copy(m_dXdxi_onFaces[dit]);
   }
}

void
MagGeom::setdXdxiOnEdges( const DisjointBoxLayout&  a_grids,
                          const IntVect&            a_ghostVect ) const
{
   if ( !m_dXdxi_onEdges.isDefined() || (m_dXdxi_onEdges.ghostVect() < a_ghostVect) ) {

      m_dXdxi_onEdges.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);

      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);
         coord_sys.getEdgeCentereddXdxi(m_dXdxi_onEdges[dit]);
      }

      m_dXdxi_onEdges.exchange();
      
      //if (m_extrablock_exchange) fillInternalGhosts(m_dXdxi_onEdges);
   }

}

void
MagGeom::getdXdxi( LevelData<EdgeDataBox>& a_dXdxi ) const
{
   CH_assert(a_dXdxi.nComp() == SpaceDim*SpaceDim);
   
   const DisjointBoxLayout& grids = a_dXdxi.disjointBoxLayout();
   const IntVect& ghosts = a_dXdxi.ghostVect();
   setdXdxiOnEdges(grids, ghosts);

   for (DataIterator dit(a_dXdxi.dataIterator()); dit.ok(); ++dit) {
      a_dXdxi[dit].copy(m_dXdxi_onEdges[dit]);
   }
}

void
MagGeom::setXphys( const DisjointBoxLayout&  a_grids,
                   const IntVect&            a_ghostVect ) const
{   
   
   if ( !m_Xphys.isDefined() || (m_Xphys.ghostVect() < a_ghostVect) ) {

      m_Xphys.define(a_grids, SpaceDim, a_ghostVect);

      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);
         coord_sys.getCellCenteredRealCoords(m_Xphys[dit]);
      }

      m_Xphys.exchange();
      
      if (m_extrablock_exchange) fillInternalGhosts(m_Xphys);
   }

}

void
MagGeom::getXphys( LevelData<FArrayBox>& a_Xphys ) const
{   
   CH_assert(a_Xphys.nComp() == SpaceDim);
   
   const DisjointBoxLayout& grids = a_Xphys.disjointBoxLayout();
   const IntVect& ghosts = a_Xphys.ghostVect();
   setXphys(grids, ghosts);

   for (DataIterator dit(a_Xphys.dataIterator()); dit.ok(); ++dit) {
      a_Xphys[dit].copy(m_Xphys[dit]);
   }
}


void
MagGeom::getPointwisePoloidalJ( LevelData<FArrayBox>& a_data ) const
{

   DataIterator dit = a_data.dataIterator();

   if ( !m_poloidal_J.isDefined() || !(a_data.ghostVect() <= m_poloidal_J.ghostVect()) ) {

      m_poloidal_J.define(a_data);

      const DisjointBoxLayout& grids = m_poloidal_J.disjointBoxLayout();

      for (dit.begin(); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
         coord_sys.getPointwisePoloidalJ(m_poloidal_J[dit]);
      }

      m_coord_sys->postProcessMetricData(m_poloidal_J);

      m_poloidal_J.exchange();
   }

   for (dit.begin(); dit.ok(); ++dit) {
      a_data[dit].copy(m_poloidal_J[dit]);
   }
}


void
MagGeom::getPointwiseN( LevelData<FluxBox>& a_N ) const
{

   DataIterator dit = a_N.dataIterator();

   if ( !m_N_face_centered.isDefined() || !(a_N.ghostVect() <= m_N_face_centered.ghostVect()) ) {

      m_N_face_centered.define(a_N);

      const DisjointBoxLayout& grids = m_N_face_centered.disjointBoxLayout();

      for (dit.begin(); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
         coord_sys.getPointwiseN(m_N_face_centered[dit]);
      }

      m_coord_sys->postProcessMetricData(m_N_face_centered);
   }

   for (dit.begin(); dit.ok(); ++dit) {
      a_N[dit].copy(m_N_face_centered[dit]);
   }
}

void
MagGeom::setPointwiseNJInverseOnFaces( const DisjointBoxLayout&  a_grids,
                                       const IntVect&            a_ghostVect ) const
{

   if ( !m_NJinverse_face_centered.isDefined() || 
        (m_NJinverse_face_centered.ghostVect() < a_ghostVect) ) {

      m_NJinverse_face_centered.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);

      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);
         coord_sys.getPointwiseNJInverse(m_NJinverse_face_centered[dit]);
      }

      m_coord_sys->postProcessMetricData(m_NJinverse_face_centered);

   }

}

void
MagGeom::getPointwiseNJInverse( LevelData<FluxBox>& a_NJinverse ) const
{
   const DisjointBoxLayout& grids = a_NJinverse.disjointBoxLayout();
   const IntVect& ghosts = a_NJinverse.ghostVect();
   setPointwiseNJInverseOnFaces(grids, ghosts);

   DataIterator dit = a_NJinverse.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_NJinverse[dit].copy(m_NJinverse_face_centered[dit]);
   }
}

void
MagGeom::setPointwiseNJInverseOnCells( const DisjointBoxLayout&  a_grids,
                                       const IntVect&            a_ghostVect ) const
{
   
   if ( !m_NJinverse_cell_centered.isDefined() ||
        (m_NJinverse_cell_centered.ghostVect() < a_ghostVect) ) {
      
      m_NJinverse_cell_centered.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);
      
      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(a_grids[dit]);
         coord_sys.getPointwiseNJInverse(m_NJinverse_cell_centered[dit]);
      }
      
      m_coord_sys->postProcessMetricData(m_NJinverse_cell_centered);
      
      m_NJinverse_cell_centered.exchange();

   }
   
}

void
MagGeom::getPointwiseNJInverse( LevelData<FArrayBox>& a_NJinverse ) const
{
   const DisjointBoxLayout& grids = a_NJinverse.disjointBoxLayout();
   const IntVect& ghosts = a_NJinverse.ghostVect();
   setPointwiseNJInverseOnCells(grids, ghosts);

   DataIterator dit = a_NJinverse.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_NJinverse[dit].copy(m_NJinverse_cell_centered[dit]);
   }
}


void
MagGeom::divideJonValid( LevelData<FArrayBox>& a_data ) const
{
   // Multiplies the argument by J on valid cells

   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   
   if (!m_second_order) {
     IntVect grown_vect = IntVect::Unit;
     setJ(grids, grown_vect);

     LevelData<FArrayBox> grown_data(grids, a_data.nComp(), grown_vect);
     for (DataIterator dit(grown_data.dataIterator()); dit.ok(); ++dit) {
       grown_data[dit].copy(a_data[dit]);
     }
     grown_data.exchange(); // JRA, exchange needed? copy works on ghost cells

     // Compute the fourth-order quotient
     for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
       const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
       cellFGToCellF(a_data[dit], grown_data[dit], m_J[dit], grids[dit], coord_sys.domain(), true);
     }
   }  
 
   else {
     setJ(grids, a_data.ghostVect());
     for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
       for (int n=0; n<a_data.nComp(); ++n) {
	 a_data[dit].divide(m_J[dit],0,n,1);
       }
     }
   }
}


void
MagGeom::multJonValid( LevelData<FArrayBox>& a_data ) const
{
   // Multiplies the argument by J on valid cells

   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   if (!m_second_order) {
     IntVect grown_vect = IntVect::Unit;
     setJ(grids, grown_vect);
     
     LevelData<FArrayBox> grown_data(grids, a_data.nComp(), grown_vect);
     for (DataIterator dit(grown_data.dataIterator()); dit.ok(); ++dit) {
       grown_data[dit].copy(a_data[dit]);
     }
     grown_data.exchange(); // JRA, exchange needed?
     //fourthOrderCellProd(a_data, J, grown_data); // JRA, requires all inputs to have same nComps()
     
     // Compute the fourth-order product
     for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
       const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
       fourthOrderCellProd(a_data[dit], grown_data[dit], m_J[dit], grids[dit], coord_sys.domain(), true);
     }
   }

   else {
     setJ(grids, a_data.ghostVect());
     for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
       for (int n=0; n<a_data.nComp(); ++n) {
         a_data[dit].mult(m_J[dit],0,n,1);
       }
     }
   }
}

void
MagGeom::multJonFaces( LevelData<FluxBox>& a_data ) const
{
   // Multiplies the argument by J on faces
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   const IntVect& ghosts = a_data.ghostVect();
   setJonFaces(grids, ghosts);
   
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const Box& box = a_data[dit].box();
      for (int n=0; n<a_data.nComp(); ++n) {
         a_data[dit].mult(m_JonFaces[dit],box,0,n,1);
      }
   }

}

void
MagGeom::multJonEdges( LevelData<EdgeDataBox>& a_data ) const
{
   // Multiplies the argument by J on edges
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   const IntVect& ghosts = a_data.ghostVect();
   setJonEdges(grids, ghosts);
   
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      for (int n=0; n<a_data.nComp(); ++n) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            FArrayBox& data_on_dir(a_data[dit][dir]);
            data_on_dir.mult(m_JonEdges[dit][dir],0,n,1);
         }
      }
   }

}

void
MagGeom::divideJonFaces( LevelData<FluxBox>& a_data ) const
{
   // Divides the argument by J on faces
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   const IntVect& ghosts = a_data.ghostVect();
   setJonFaces(grids, ghosts);
   
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const Box& box = a_data[dit].box();
      for (int n=0; n<a_data.nComp(); ++n) {
         a_data[dit].divide(m_JonFaces[dit],box,0,n,1); // r=0?
      }
   }

}

void
MagGeom::divideJonEdges( LevelData<EdgeDataBox>& a_data ) const
{
   // Divides the argument by J on edges
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   const IntVect& ghosts = a_data.ghostVect();
   setJonEdges(grids, ghosts);
   
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      for (int n=0; n<a_data.nComp(); ++n) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            FArrayBox& data_on_dir(a_data[dit][dir]);
            data_on_dir.divide(m_JonEdges[dit][dir],0,n,1); // r=0?
         }
      }
   }

}

void
MagGeom::applyOnAxisCurlCorrection( LevelData<EdgeDataBox>& a_Curl,
                              const LevelData<FArrayBox>&   a_val) const
{
   // set values at r=0 to be a_Curl(r=0)[dir=1] = 2*a_val(r=dr/2)/(dr/2);
   const DisjointBoxLayout& grids = a_val.disjointBoxLayout();

   RealVect dX_mapped;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
      dX_mapped = coord_sys.getMappedCellSize();
      if ( coord_sys.isAxisymmetric() ) {
         const int dir0 = 0;
         const int dir1 = 1;
         const FArrayBox& X_phys = m_edge_centered_real_coords[dit][dir1];
         const FArrayBox& a_val_on_patch = a_val[dit];
         FArrayBox& a_Curl1_on_patch = a_Curl[dit][dir1];
         const Box& thisedgebox = a_Curl1_on_patch.box();
         
         FORT_ON_AXIS_CURL_CORRECTION( CHF_BOX(thisedgebox),
                                       CHF_CONST_REAL(dX_mapped[dir0]),
                                       CHF_CONST_FRA1(X_phys,0),
                                       CHF_CONST_FRA1(a_val_on_patch,0),
                                       CHF_FRA1(a_Curl1_on_patch,0) );   
         
      } 
   }
}


void MagGeom::multPointwiseJ( LevelData<FArrayBox>& a_u,
                              const BoundaryBoxLayout& a_bdry_layout ) const
{
   const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
   for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
      const Box interior_box( a_bdry_layout.interiorBox( dit ) );
      const MagBlockCoordSys& coord_sys( getBlockCoordSys( interior_box ) );
      const Box& fill_box( a_u[dit].box() );
      FArrayBox Xi( fill_box, SpaceDim );
      coord_sys.getCellCenteredMappedCoords( Xi );
      FArrayBox J( fill_box, 1 );
      coord_sys.pointwiseJ( J, Xi, fill_box );
      a_u[dit] *= J;
   }
}


void
MagGeom::getFaceAreas( LevelData<FluxBox>& a_areas ) const
{
  const DisjointBoxLayout& grids = a_areas.disjointBoxLayout();

  DataIterator dit = a_areas.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FluxBox& this_areas = a_areas[dit];
    const Box& box = this_areas.box();

    const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);

    FluxBox N(box, SpaceDim*SpaceDim);

    // Force the base class implementation to avoid the 2piR axisymmetric modifications
    coord_sys.NewFourthOrderCoordSys::getN(N, box);

    // This gets the ratio of physical areas to computational areas
    coord_sys.magnitudeN(this_areas, N, box);
  }

  m_coord_sys->postProcessMetricData(a_areas);
}



void
MagGeom::computeMetricTermProductAverage( LevelData<FluxBox>&       a_product,
                                          const LevelData<FluxBox>& a_data,
                                          const bool                a_fourthOrder ) const
{
   CH_TIME("MagGeom::computeMetricTermProductAverage()");
   
   const DisjointBoxLayout& grids = a_product.disjointBoxLayout();
   const IntVect& ghosts = IntVect::Unit;
   setMetricTerms(grids, ghosts); 

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box& box = grids[dit];
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(box);

      block_coord_sys.computeMetricTermProductAverage( a_product[dit], a_data[dit], 
                                                       m_metrics[dit], box, a_fourthOrder );
   }

}

void
MagGeom::computedxidXProductAverage( LevelData<FluxBox>&  a_product,
                               const LevelData<FluxBox>&  a_data,
                               const bool                 a_fourthOrder ) const
{
   CH_TIME("MagGeom::computedxidXProductAverage()");
   
   // takes physical a_data and converts to contravar on faces
   // and multiplies by face area. 
   //
   // Different from computeMetricTermProductAverage
   // in that the result is not multiplied by Jacobian
   //
   
   const DisjointBoxLayout& grids = a_product.disjointBoxLayout();
   const IntVect& ghosts = IntVect::Unit;
   setPointwiseNJInverseOnFaces(grids, ghosts); 

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box& box = grids[dit];
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(box);

      block_coord_sys.computeMetricTermProductAverage( a_product[dit], a_data[dit], 
                                                       m_NJinverse_face_centered[dit], box, a_fourthOrder );
      
      const RealVect& faceArea = block_coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& this_product = a_product[dit][dir];
         this_product.mult(faceArea[dir]);
      }
   }

}



void
MagGeom::computeMappedGridDivergence( const LevelData<FluxBox>& a_F,
                                      LevelData<FArrayBox>&     a_divergence,
                                      const bool                a_fourthOrder ) const
{
   CH_assert(a_F.ghostVect() >= IntVect::Unit);
   const DisjointBoxLayout& grids = a_divergence.disjointBoxLayout();

   LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
   computeMetricTermProductAverage(NTF_normal, a_F, a_fourthOrder);

   // Average the normal component of NTF consistent across block interfaces
   averageAtBlockBoundaries(NTF_normal);

   RealVect fakeDx = RealVect::Unit;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      simpleDivergence(a_divergence[dit], NTF_normal[dit], grids[dit], fakeDx);
   }
}

void
MagGeom::mappedGridDivergenceFromFluxNorms( const LevelData<FluxBox>&   a_F,
                                                  LevelData<FArrayBox>& a_divergence) const
{
   CH_TIME("MagGeom::mappedGridDivergenceFromFluxNorms()");
   
   //CH_assert(a_F.ghostVect() >= IntVect::Unit);
   CH_assert(a_F.nComp() == a_divergence.nComp());
   const DisjointBoxLayout& grids = a_divergence.disjointBoxLayout();

   // Average the normal component of a_F consistent across block interfaces
   LevelData<FluxBox> NTF_normal(grids, a_F.nComp(), IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      NTF_normal[dit].copy(a_F[dit]);
   }
   averageAtBlockBoundaries(NTF_normal);
   //averageAtBlockBoundaries(a_F);

   RealVect fakeDx = RealVect::Unit;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      simpleDivergence(a_divergence[dit], NTF_normal[dit], grids[dit], fakeDx);
      //simpleDivergence(a_divergence[dit], a_F[dit], grids[dit], fakeDx);
   }
}

void
MagGeom::mappedGridDivergenceFromFluxNorms( const FluxBox&    a_F,
                                                  FArrayBox&  a_divergence) const
{
   CH_TIME("MagGeom::mappedGridDivergenceFromFluxNorms()");
   CH_assert(a_F.nComp() == a_divergence.nComp());
   const Box& this_box = a_divergence.box();

   // Average the normal component of a_F consistent across block interfaces
   FluxBox NTF_normal(this_box, a_F.nComp());
   NTF_normal.copy(a_F);
   //averageAtBlockBoundaries(NTF_normal);

   RealVect fakeDx = RealVect::Unit;
   simpleDivergence(a_divergence, NTF_normal, this_box, fakeDx);
}

void
MagGeom::mappedGridCurlFromEdgeTans( const LevelData<EdgeDataBox>&   a_E,
                                     const int                       a_face_dir,
                                           LevelData<FArrayBox>&     a_face_curl) const
{
   CH_TIME("MagGeom::mappedGridCurlFromEdgeTans()");

   // Assumes a_E is line-average of covarient component of vector E
   // at cell edges

   CH_assert(a_E.nComp() == a_face_curl.nComp());
   const DisjointBoxLayout& grids = a_face_curl.disjointBoxLayout();

   int dirj, dirk;
   dirj = a_face_dir + 1; // a_face_dir + 1
   dirj = dirj % SpaceDim;
   dirk = a_face_dir + 2; // a_face_dir + 2
   dirk = a_face_dir % SpaceDim;

   //cout << "dirj = " << dirj << endl;
   //cout << "dirk = " << dirk << endl;

   RealVect dX_mapped;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
   
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
      dX_mapped = coord_sys.getMappedCellSize();

      const FArrayBox& E_dirj_on_patch = a_E[dit][dirj]; // Ez (on x-edges) 
      const FArrayBox& E_dirk_on_patch = a_E[dit][dirk]; // Ex (on z-edges)
      FArrayBox& face_curl_on_patch = a_face_curl[dit]; 
      const Box& thisbox = face_curl_on_patch.box();

      FORT_SIMPLE_CURL( CHF_BOX(thisbox),
                        CHF_CONST_INT(dirj), // 1
                        CHF_CONST_INT(dirk), // 0
                        CHF_CONST_REAL(dX_mapped[dirj]), // dZ
                        CHF_CONST_REAL(dX_mapped[dirk]), // dX
                        CHF_CONST_FRA(E_dirj_on_patch),  // Ez
                        CHF_CONST_FRA(E_dirk_on_patch),  // Ex
                        CHF_FRA(face_curl_on_patch) );   
   
   }
   //a_face_curl.exchange();

}

void
MagGeom::mappedGridCurlofVirtComp( const LevelData<FArrayBox>&    a_F,
                                         LevelData<EdgeDataBox>&  a_edge_curl) const
{
   
   // a_F is covariant at cell center
   // a_edge_curl is Jacobian*curl(a_F)_contravar at cell edges 
   // ( Ja*\curl(a_F)\cdot g^l )

   // Special case for 2D with vector F in out-of-plane direction (y-direction)
   // Assumes a_F is covariant y-component of vector F
   // For 2D, dir=0 is x and dir=1 is z

   // Also need a_F.ghostVect() == a_edge_curl.ghostVect()...May need to rethink 
   // iteration boxes for this. Note that using a_F.box() wont do last edge in each
   // dir, which shouldn't matter since it's ghost cell, but yea..

   CH_assert(SpaceDim==2);
   CH_assert(a_F.nComp() == a_edge_curl.nComp());
   const DisjointBoxLayout& grids = a_edge_curl.disjointBoxLayout();

   RealVect dX_mapped;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
   
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
      dX_mapped = coord_sys.getMappedCellSize();

      //FArrayBox& Jx_on_patch = a_edge_curl[dit][0]; // lives on z-faces
      //FArrayBox& Jz_on_patch = a_edge_curl[dit][1]; // lives on x-faces
      const FArrayBox& F_on_patch = a_F[dit]; 
      int dXdir = 1;
      for (int dir(0); dir<SpaceDim; dir++) {
         if(dir==1) dXdir = 0;
         FArrayBox& edge_on_patch = a_edge_curl[dit][dir];
         const Box& thisedgebox = edge_on_patch.box();
  
         FORT_SIMPLE_CURL_VIRT( CHF_BOX(thisedgebox),
                                CHF_CONST_INT( dir ),
                                CHF_CONST_REAL(dX_mapped[dXdir]),
                                CHF_CONST_FRA(F_on_patch),
                                CHF_FRA(edge_on_patch) );   
      }
   
   }
   //a_edge_curl.exchange();

}

void
MagGeom::convertPhysToCovar( LevelData<FArrayBox>&  a_F,
                       const int                    a_inverse ) const
{
   CH_TIME("MagGeom::convertPhysToCovar()");
   // convert physical vector to covariant at cell center
   // or covariant to physical if a_inverse not zero

   CH_assert(a_F.nComp() == SpaceDim);
   const DisjointBoxLayout& grids = a_F.disjointBoxLayout();
   const IntVect ghosts = a_F.ghostVect();   
   setPointwiseNJInverseOnCells(grids, ghosts);
   setdXdxiOnCells(grids, ghosts);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
   
      const FArrayBox& dXdxi_on_patch( m_dXdxi[dit] );
      const FArrayBox& dxidX_on_patch( m_NJinverse_cell_centered[dit] );
      
      FArrayBox& F_on_patch = a_F[dit]; 
      const Box& thisbox = F_on_patch.box();

      FORT_PHYS_TO_COVAR( CHF_BOX(thisbox),
                          CHF_CONST_FRA(dXdxi_on_patch),
                          CHF_CONST_FRA(dxidX_on_patch),
                          CHF_CONST_INT(a_inverse),
                          CHF_FRA(F_on_patch) );   
         
   }

}

void
MagGeom::convertPhysToCovar( LevelData<FluxBox>&  a_F,
                       const int                  a_inverse ) const
{
   CH_TIME("MagGeom::convertPhysToCovar()");
   // convert physical vector to covariant at cell face
   // or covariant to physical if a_inverse not zero

   CH_assert(a_F.nComp() == SpaceDim);
   const DisjointBoxLayout& grids = a_F.disjointBoxLayout();
   const IntVect ghosts = a_F.ghostVect();   
   setPointwiseNJInverseOnFaces(grids, ghosts);
   setdXdxiOnFaces(grids, ghosts);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
   
      for (int dir=0; dir<SpaceDim; ++dir) {
      
         const FArrayBox& dXdxi_on_dir( m_dXdxi_onFaces[dit][dir] );
         const FArrayBox& dxidX_on_dir( m_NJinverse_face_centered[dit][dir] );
      
         FArrayBox& F_on_dir = a_F[dit][dir]; 
         const Box& thisbox = F_on_dir.box();
      
         FORT_PHYS_TO_COVAR( CHF_BOX(thisbox),
                             CHF_CONST_FRA(dXdxi_on_dir),
                             CHF_CONST_FRA(dxidX_on_dir),
                             CHF_CONST_INT(a_inverse),
                             CHF_FRA(F_on_dir) );   

      }
   
   }

}

void
MagGeom::convertPhysToContravar( LevelData<FArrayBox>&  a_F,
                           const int                    a_inverse ) const
{
   CH_TIME("MagGeom::convertPhysToContravar()");
   // convert physical vector to contravariant at cell center
   // or contravariant to physical if a_inverse not zero

   CH_assert(a_F.nComp() == SpaceDim);
   const DisjointBoxLayout& grids = a_F.disjointBoxLayout();
   const IntVect ghosts = a_F.ghostVect();   
   setPointwiseNJInverseOnCells(grids, ghosts); 
   setdXdxiOnCells(grids, ghosts);

   for (DataIterator dit(grids); dit.ok(); ++dit) {

      const FArrayBox& dXdxi_on_patch( m_dXdxi[dit] );
      const FArrayBox& dxidX_on_patch( m_NJinverse_cell_centered[dit] );

      FArrayBox& F_on_patch = a_F[dit]; 
      const Box& thisbox = F_on_patch.box();

      FORT_PHYS_TO_CONTRAVAR( CHF_BOX(thisbox),
                              CHF_CONST_FRA(dXdxi_on_patch),
                              CHF_CONST_FRA(dxidX_on_patch),
                              CHF_CONST_INT(a_inverse),
                              CHF_FRA(F_on_patch) );   
         
   }

}

void
MagGeom::convertPhysToContravar( LevelData<FluxBox>&  a_F,
                           const int                  a_inverse ) const
{
   CH_TIME("MagGeom::convertPhysToContravar()");
   // convert physical vector to contravariant at cell face
   // or contravariant to physical if a_inverse not zero

   CH_assert(a_F.nComp() == SpaceDim);
   const DisjointBoxLayout& grids = a_F.disjointBoxLayout();
   const IntVect ghosts = a_F.ghostVect();   
   setPointwiseNJInverseOnFaces(grids, ghosts); 
   setdXdxiOnFaces(grids, ghosts);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
   
      for (int dir=0; dir<SpaceDim; ++dir) {
      
         const FArrayBox& dXdxi_on_dir( m_dXdxi_onFaces[dit][dir] );
         const FArrayBox& dxidX_on_dir( m_NJinverse_face_centered[dit][dir] );
      
         FArrayBox& F_on_dir = a_F[dit][dir]; 
         const Box& thisbox = F_on_dir.box();
      
         FORT_PHYS_TO_CONTRAVAR( CHF_BOX(thisbox),
                                 CHF_CONST_FRA(dXdxi_on_dir),
                                 CHF_CONST_FRA(dxidX_on_dir),
                                 CHF_CONST_INT(a_inverse),
                                 CHF_FRA(F_on_dir) );   

      }
   
   }

}

void
MagGeom::convertContravarToCovar( LevelData<EdgeDataBox>&  a_F,
                            const int                      a_inverse ) const
{
   CH_TIME("MagGeom::convertContravarToCovar()");
   // convert contravariant vector on edges
   // to covariant or vice versa if a_inverse not zero

   CH_assert(a_F.nComp() == 1);
   const DisjointBoxLayout& grids = a_F.disjointBoxLayout();
   const IntVect& ghosts = a_F.ghostVect();
   setdXdxiOnEdges(grids, ghosts);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
         
      for (int dir=0; dir<SpaceDim; ++dir) {
      
         const FArrayBox& dXdxi_on_dir( m_dXdxi_onEdges[dit][dir] );
      
         FArrayBox& F_on_dir = a_F[dit][dir]; 
         const Box& thisbox = F_on_dir.box();

         FORT_CONTRA_TO_COVAR( CHF_BOX(thisbox),
                               CHF_CONST_FRA(dXdxi_on_dir),
                               CHF_CONST_INT(dir),
                               CHF_CONST_INT(a_inverse),
                               CHF_FRA1(F_on_dir,0) );   
      
      } 

   }

}


void MagGeom::computeJaPhysStrainFlux( LevelData<FluxBox>&    a_JaW0_Flux_phys, 
                                       LevelData<FluxBox>&    a_JaW1_Flux_phys,
                                 const LevelData<FArrayBox>&  a_velocity_phys ) const
{
   CH_TIME("MagGeom::computeJaPhysStrainFlux()");

   // calculate physical rate of strain tensor on cell faces
   // (multiplied by Jacobian)
   //
   // W = \nabla(V) + (\nabla(V))^T - 2/3*I\nabla\cdot V
   //
   // JaW0 = Ja*{W_00, W_01}
   // JaW1 = Ja*{W_10, W_11}
   //

   const DisjointBoxLayout& grids( a_velocity_phys.getBoxes() );

   LevelData<FluxBox> gradV0_cf(grids, SpaceDim, 2*IntVect::Unit);
   LevelData<FluxBox> gradV1_cf(grids, SpaceDim, 2*IntVect::Unit);
   LevelData<FArrayBox> dummyFArray0(grids, 1, 2*IntVect::Unit);
   LevelData<FArrayBox> dummyFArray(grids, SpaceDim, 2*IntVect::Unit);


   // compute divV on faces in a similar way as gradient in order to use
   // tighter stencil. Opposed to compuing divV at cc and interpolaring to faces
   //


   // convert velocity to contravariant
   // and calc gradients on cell faces
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      dummyFArray[dit].copy(a_velocity_phys[dit]);
   }
   convertPhysToContravar(dummyFArray,0);
   for (int n=0; n<SpaceDim; ++n) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         dummyFArray0[dit].copy(dummyFArray[dit],n,0,1);
         dummyFArray0[dit].mult(m_J[dit],0,0,1); 
      }
      if(n==0) {
         //computeMappedPoloidalGradientWithGhosts( dummyFArray0, gradV0_cf, 2 );
         computeMappedGradient( dummyFArray0, gradV0_cf, 2 );
      }
      if(n==1) {
         //computeMappedPoloidalGradientWithGhosts( dummyFArray0, gradV1_cf, 2 );
         computeMappedGradient( dummyFArray0, gradV1_cf, 2 );
      }
   }
   
   // compute Ja*divU on cell faces
   //
   LevelData<FluxBox> JadivU_cf(grids, 1, IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      JadivU_cf[dit].copy(gradV0_cf[dit],0,0,0,1);
      JadivU_cf[dit].copy(gradV1_cf[dit],1,1,0,1);
      const FArrayBox& gradV0_on_ydir = gradV0_cf[dit][1];
      const FArrayBox& gradV1_on_xdir = gradV1_cf[dit][0];
            FArrayBox& divU_on_ydir = JadivU_cf[dit][1];
            FArrayBox& divU_on_xdir = JadivU_cf[dit][0];
      divU_on_ydir.plus(gradV0_on_ydir,0,0,1);
      divU_on_xdir.plus(gradV1_on_xdir,1,0,1);
   }

   // compute mapped \nabla Vphys matrix on cell faces
   // and then convert to physical
   //
   // NOTE: Redefining gradV_cf containers used above
   //       as will be used in the following FORT call
   //
   for (int n=0; n<SpaceDim; ++n) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         dummyFArray0[dit].copy(a_velocity_phys[dit],n,0,1);
      }
      if(n==0) {
         //computeMappedPoloidalGradientWithGhosts( dummyFArray0, gradV0_cf, 2 );
         computeMappedGradient( dummyFArray0, gradV0_cf, 2 );
         convertPhysToCovar(gradV0_cf,1);
      }
      if(n==1) {
         //computeMappedPoloidalGradientWithGhosts( dummyFArray0, gradV1_cf, 2 );
         computeMappedGradient( dummyFArray0, gradV1_cf, 2 );
         convertPhysToCovar(gradV1_cf,1);
      }
   }
   multJonFaces(gradV0_cf);
   multJonFaces(gradV1_cf);   

   // compute physical stress tensor components on cell faces
   // JaW0 = Ja*{W_00, W_01}
   // JaW1 = Ja*{W_10, W_11}
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; ++dir) {

         const FArrayBox& divU_on_dir   = JadivU_cf[dit][dir];
         const FArrayBox& gradV0_on_dir = gradV0_cf[dit][dir];
         const FArrayBox& gradV1_on_dir = gradV1_cf[dit][dir];
         FArrayBox& Pi0_on_dir = a_JaW0_Flux_phys[dit][dir];
         FArrayBox& Pi1_on_dir = a_JaW1_Flux_phys[dit][dir];

         const Box& thisbox = Pi0_on_dir.box();
         FORT_EVAL_STRAIN_TENSOR_COMPS( CHF_BOX(thisbox),
                                        CHF_CONST_FRA1(divU_on_dir,0),
                                        CHF_CONST_FRA(gradV0_on_dir),
                                        CHF_CONST_FRA(gradV1_on_dir),
                                        CHF_FRA(Pi0_on_dir),
                                        CHF_FRA(Pi1_on_dir) );

      }

   }

}

void
MagGeom::projectOntoFluxSurfaceDir( FArrayBox&  a_data,
                              const DataIndex&  a_dit,
                              const int         a_dir,
                              const int         a_inverse ) const
{
#if 0
   CH_TIME("MagGeom::projectOntoFluxSurfaceDir()");
   CH_assert(a_data.nComp() == SpaceDim);
  
   const DisjointBoxLayout& grids = m_J.disjointBoxLayout();

   const FArrayBox& N_on_patch     = m_N[a_dit];
   const FArrayBox& dXdxi_on_patch = m_dXdxi[a_dit];
   
   const Box& box = a_data.box();
   const FArrayBox& dxidX_on_patch( m_NJinverse_cell_centered[a_dit] );

   FORT_PROJECT_PHYS_TO_MAPPED_DIR( CHF_BOX(box),
                                    CHF_CONST_FRA(dxidX_on_patch),
                                    CHF_CONST_FRA(dXdxi_on_patch),
                                    CHF_CONST_INT(a_dir),
                                    CHF_CONST_INT(a_inverse),
                                    CHF_FRA(a_data) );
#endif   
}


void
MagGeom::averageAtBlockBoundaries(LevelData<FluxBox>& a_data) const
{
   if ( m_coord_sys->numBlocks() > 1 && !m_sheared_mb_geom ) {

      const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

      BlockRegister blockRegister(m_coordSysPtr, grids, 0);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int idir = 0; idir < SpaceDim; idir++) {
            for (SideIterator sit; sit.ok(); ++sit) {
               Side::LoHiSide side = sit();
               if (blockRegister.hasInterface(dit(), idir, side)) {
                  FArrayBox flux_comp(a_data[dit][idir].box(), a_data.nComp());
                  flux_comp.copy(a_data[dit][idir]);
                  blockRegister.storeFlux(flux_comp, dit(), idir, side);
               }
            }
         }
      }
      blockRegister.close();

      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = m_coord_sys->boundaries();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& baseBox = grids[dit];
         int block_number = m_coord_sys->whichBlock(grids[dit]);
         int faceID = 0;
         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
            for (int idir = 0; idir < SpaceDim; idir++) {
               if (blockRegister.hasInterface(dit(), idir, side)) {
                  // maybe better if this is done inside BlockRegister
                  const BlockBoundary& bb = boundaries[block_number][faceID];
                  int reorientFace = bb.reorientFace(idir);
                  Box faceBox = adjCellBox(baseBox, idir, side, 1);
                  // if Lo, then shift +1; if Hi, then shift -1
                  faceBox.shiftHalf(idir, -sign(side));
                  Side::LoHiSide sideOther = flip(side);
                  // Need to define these FABs.
                  FArrayBox fluxThisFab(faceBox, a_data.nComp());
                  FArrayBox fluxOtherFab(faceBox, a_data.nComp());
                  blockRegister.getFlux(fluxThisFab, dit(),
                                        idir, side, side);
                  fluxThisFab.mult(reorientFace * 0.5);
                  blockRegister.getFlux(fluxOtherFab, dit(),
                                        idir, side, sideOther);
                  fluxOtherFab.mult(0.5);
                  fluxThisFab += fluxOtherFab;
                  a_data[dit][idir].copy(fluxThisFab);
               }
               faceID++;
            } // iterate over dimensions
         } // iterate over sides
      }
   }
}



void
MagGeom::exchangeTransverseAtBlockBoundaries( LevelData<FluxBox>& a_data ) const
{
   if ( m_coordSysPtr->numBlocks() > 1 ) {

      const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

      for ( int tdir=0; tdir<SpaceDim; ++tdir ) {   // Transverse directions

         m_exchange_transverse_block_register->setToZero( a_data.nComp() );

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            for (int idir = 0; idir < SpaceDim; idir++) {
               if ( idir != tdir ) {
                  for (SideIterator sit; sit.ok(); ++sit) {
                     Side::LoHiSide side = sit();
                     if (m_exchange_transverse_block_register->hasInterface(dit(), idir, side) ) {
                        FArrayBox flux_comp(a_data[dit][tdir].box(), a_data.nComp());
                        flux_comp.copy(a_data[dit][tdir]);
                        flux_comp.shiftHalf(tdir,1);
                        flux_comp.shiftHalf(idir, sign(side));
                        m_exchange_transverse_block_register->storeAux(flux_comp, dit(), idir, side);
                     }
                  }
               }
            }
         }
         m_exchange_transverse_block_register->close();

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            for (int idir = 0; idir < SpaceDim; idir++) {
               if ( idir != tdir ) {
                  for (SideIterator sit; sit.ok(); ++sit) {
                     Side::LoHiSide side = sit();
                     if (m_exchange_transverse_block_register->hasInterface(dit(), idir, side) ) {
                        Box faceBox = adjCellBox(grids[dit], idir, side, 1);
                        faceBox.growHi(tdir,1);
                        faceBox.shiftHalf(idir, -sign(side));
                        FArrayBox fluxOtherFab(faceBox, a_data.nComp());
                        m_exchange_transverse_block_register->getAux(fluxOtherFab, dit(),
                                                                     idir, side, side);
                        fluxOtherFab.shiftHalf(idir, sign(side));
                        fluxOtherFab.shiftHalf(tdir,-1);
                        a_data[dit][tdir].copy(fluxOtherFab);
                     }
                  } // iterate over sides
               }
            } // iterate over dimensions
         }
      }

      a_data.exchange();
   }
}


void
MagGeom::applyAxisymmetricCorrection( LevelData<FluxBox>& a_data ) const
{
   CH_TIME("MagGeom::applyAxisymmetricCorrection()");

   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);

     if ( coord_sys.isAxisymmetric() && (!coord_sys.isPointwiseMetrics()) ) {

       FluxBox& this_data = a_data[dit];
       Box box(this_data.box());
       FluxBox R(box, 1);

       coord_sys.getPointwiseMajorRadius(R);

       for (int dir=0; dir<SpaceDim; ++dir) {
          this_data[dir].mult(R[dir],0,RADIAL_DIR,1);
       }
     }
   }
}


void MagGeom::computeFieldData( LevelData<FArrayBox>& a_BField,
                                LevelData<FArrayBox>& a_BFieldMag,
                                LevelData<FArrayBox>& a_BFieldDir,
                                LevelData<FArrayBox>& a_gradBFieldMag,
                                LevelData<FArrayBox>& a_curlBFieldDir,
                                LevelData<FArrayBox>& a_BFieldDirdotcurlBFieldDir ) const
{
   const DisjointBoxLayout& grids = a_BField.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);

      coord_sys.computeFieldData(5,
                                 a_BField[dit],
                                 a_BFieldMag[dit],
                                 a_BFieldDir[dit],
                                 a_gradBFieldMag[dit],
                                 a_curlBFieldDir[dit],
                                 a_BFieldDirdotcurlBFieldDir[dit]);
   }
}



void MagGeom::computeFieldData( LevelData<FluxBox>& a_BField,
                                LevelData<FluxBox>& a_BFieldMag,
                                LevelData<FluxBox>& a_BFieldDir,
                                LevelData<FluxBox>& a_gradBFieldMag,
                                LevelData<FluxBox>& a_curlBFieldDir,
                                LevelData<FluxBox>& a_BFieldDirdotcurlBFieldDir ) const
{
   const DisjointBoxLayout& grids = a_BField.disjointBoxLayout();

   // Compute the field cell averages
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);

      FluxBox & thisB = a_BField[dit];
      FluxBox & thisBmag = a_BFieldMag[dit];
      FluxBox & thisbunit = a_BFieldDir[dit];
      FluxBox & thisgradB = a_gradBFieldMag[dit];
      FluxBox & thiscurlB = a_curlBFieldDir[dit];
      FluxBox & thisbdotcurlb = a_BFieldDirdotcurlBFieldDir[dit];

      for (int dir=0; dir<SpaceDim; ++dir) {
         coord_sys.computeFieldData(dir,
                                    thisB[dir],
                                    thisBmag[dir],
                                    thisbunit[dir],
                                    thisgradB[dir],
                                    thiscurlB[dir],
                                    thisbdotcurlb[dir]);
      }
   }

   if ( m_correct_field ) {

      if( procID() == 0 ) {
         cout << endl << "Performing a Poisson solve to correct the magnetic field" << endl << endl;
      }

#ifdef TEST_DIVERGENCE_CLEANING
      plotFaceData( string("before_correction_BField"), a_BField, 0. );
      plotFaceData( string("before_correction_BFieldMag"), a_BFieldMag, 0. );
      plotFaceData( string("before_correction_BFieldDir"), a_BFieldDir, 0. );
      plotFaceData( string("before_correction_gradB"), a_gradBFieldMag, 0. );
      plotFaceData( string("before_correction_curlBFieldDir"), a_curlBFieldDir, 0. );
      plotFaceData( string("before_correction_BFieldDirdotcurlBFieldDir"), a_BFieldDirdotcurlBFieldDir, 0. );
#endif

      if ( !m_BFieldCorrection.isDefined() ) {
         CH_assert(a_BField.ghostVect() >= IntVect::Unit);
         LevelData<FluxBox> flux( grids, SpaceDim, IntVect::Unit );
         LevelData<FArrayBox> uncorrected_divergence( grids, 1, IntVect::Zero );
#ifdef TEST_DIVERGENCE_CLEANING
         Box domain_box = grids.physDomain().domainBox();
#endif

         //LevelData<FArrayBox> cell_volume(grids, 1, IntVect::Zero);
         //getCellVolumes(cell_volume);
         setCellVolumes(grids, IntVect::Zero);

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            FluxBox& this_flux = flux[dit];
            FluxBox& this_BField = a_BField[dit];
            for (int dir=0; dir<SpaceDim; ++dir) {
               FArrayBox& this_flux_dir = this_flux[dir];
               FArrayBox& this_BField_dir = this_BField[dir];

               this_flux_dir.copy( this_BField_dir, 0, 0, 1 );
               this_flux_dir.copy( this_BField_dir, 2, 1, 1 );
            }
         }

         applyAxisymmetricCorrection( flux );

         fourthOrderAverage(flux);

         computeMappedGridDivergence( flux, uncorrected_divergence, true );

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            uncorrected_divergence[dit] /= m_cell_volume[dit];
         }
         
#ifdef TEST_DIVERGENCE_CLEANING
         WriteMappedUGHDF5("uncorrected_divergence", grids, uncorrected_divergence, *m_coord_sys, domain_box, 0. );
#endif

         ParmParse pp( "field_correction" );
         Poisson * poisson = new Poisson( pp, *this );
         CH_assert(!poisson->secondOrder());   // This function assumes fourth-order

         poisson->setOperatorCoefficients( m_coord_sys->getDivergenceCleaningBC() );

         string method;  // dummy variable

         // Defaults; over-ridden if in ParmParse object
         double tol = 1.e-12;
         int max_iter = 50;
         bool verbose = false;
         ParmParse pp_linear_solver( ((string)pp.prefix() + ".linear_solver").c_str());
         poisson->parseMethodAndParams(pp_linear_solver, method, tol, max_iter, verbose);

         // Defaults; over-ridden if in ParmParse object
         double precond_tol = 0.;
         int precond_max_iter = 1;
         bool precond_verbose = false;
         ParmParse pp_ls_precond( ((string)pp_linear_solver.prefix() + ".precond").c_str());
         poisson->parseMethodAndParams(pp_ls_precond, method, precond_tol, precond_max_iter, precond_verbose);

         // Defaults; over-ridden if in ParmParse object
         double precond_precond_tol = 0.;
         int precond_precond_max_iter = 0;
         bool precond_precond_verbose = false;
         ParmParse pp_ls_precond_precond( ((string)pp_ls_precond.prefix() + ".precond").c_str());
         poisson->parseMethodAndParams(pp_ls_precond_precond, method, precond_precond_tol,
                                       precond_precond_max_iter, precond_precond_verbose);

         poisson->setConvergenceParams(tol, max_iter, verbose,
                                       precond_tol, precond_max_iter, precond_verbose,
                                       precond_precond_tol, precond_precond_max_iter, precond_precond_verbose);

         LevelData<FArrayBox> phi( grids, 1, 4*IntVect::Unit );
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            phi[dit].setVal(0.);
         }

         poisson->solve( uncorrected_divergence, phi );

#ifdef TEST_DIVERGENCE_CLEANING
         WriteMappedUGHDF5("field_correction_potential", grids, phi, *m_coord_sys, domain_box, 0. );
#endif

         poisson->fillInternalGhosts(phi);

         // Compute cell-centered -grad phi
         poisson->computePoloidalFieldWithBCs( phi, flux, false );

         fillTransversePhysicalGhosts(flux);

#ifdef TEST_DIVERGENCE_CLEANING
         m_BFieldCorrection.define( grids, 3, IntVect::Zero);
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            FluxBox& this_BFieldCorrection = m_BFieldCorrection[dit];
            FluxBox& this_flux = flux[dit];

            this_BFieldCorrection.setVal(0.);
            for (int dir=0; dir<SpaceDim; ++dir) {
               FArrayBox& this_BFieldCorrection_dir = this_BFieldCorrection[dir];
               FArrayBox& this_flux_dir = this_flux[dir];

               Box valid_box_dir = surroundingNodes(grids[dit],dir);
               this_BFieldCorrection_dir.copy( this_flux_dir, valid_box_dir, 0, valid_box_dir, 0, 1 );
               this_BFieldCorrection_dir.copy( this_flux_dir, valid_box_dir, 1, valid_box_dir, 2, 1 );
            }
         }
#endif

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            FluxBox& this_BField = a_BField[dit];
            FluxBox& this_flux = flux[dit];
            for (int dir=0; dir<SpaceDim; ++dir) {
               FArrayBox& this_BField_dir = this_BField[dir];
               FArrayBox& this_flux_dir = this_flux[dir];

               this_flux_dir.negate();
                  
               IntVect grow_vect = IntVect::Unit;
               grow_vect[dir] = 0;
               Box box = surroundingNodes(grow(grids[dit],grow_vect),dir);

               this_flux_dir.plus(this_BField_dir,box,0,0,1);
               this_flux_dir.plus(this_BField_dir,box,2,1,1);

               this_BField_dir.copy( this_flux_dir, box, 0, box, 0, 1 );
               this_BField_dir.copy( this_flux_dir, box, 1, box, 2, 1 );
            }
         }

#ifdef TEST_DIVERGENCE_CLEANING

         {
            applyAxisymmetricCorrection( flux );

            fourthOrderAverage(flux);

            LevelData<FArrayBox> corrected_divergence( grids, 1, IntVect::Zero );
            computeMappedGridDivergence( flux, corrected_divergence, true );

            for (DataIterator dit(grids); dit.ok(); ++dit) {
               corrected_divergence[dit] /= cell_volume[dit];
            }

            WriteMappedUGHDF5("corrected_divergence", grids, corrected_divergence, *m_coord_sys, domain_box, 0.);
         }
#endif

         delete poisson;

#ifdef TEST_DIVERGENCE_CLEANING
         LevelData<FluxBox> BFieldMagChange(grids, a_BFieldMag.nComp(), a_BFieldMag.ghostVect());
         LevelData<FluxBox> BFieldDirChange(grids, a_BFieldDir.nComp(), a_BFieldDir.ghostVect());
         LevelData<FluxBox> gradBFieldMagChange(grids, a_gradBFieldMag.nComp(), a_gradBFieldMag.ghostVect());
         LevelData<FluxBox> curlBFieldDirChange(grids, a_curlBFieldDir.nComp(), a_curlBFieldDir.ghostVect());
         LevelData<FluxBox> BFieldDirdotcurlBFieldDirChange(grids, a_BFieldDirdotcurlBFieldDir.nComp(),
                                                            a_BFieldDirdotcurlBFieldDir.ghostVect());
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            BFieldMagChange[dit].copy(a_BFieldMag[dit]);
            BFieldDirChange[dit].copy(a_BFieldDir[dit]);
            gradBFieldMagChange[dit].copy(a_gradBFieldMag[dit]);
            curlBFieldDirChange[dit].copy(a_curlBFieldDir[dit]);
            BFieldDirdotcurlBFieldDirChange[dit].copy(a_BFieldDirdotcurlBFieldDir[dit]);
         }
#endif

         // This needs to be redone.  The current implementation is too approximate
         // in its calculation of the differentiated field quantities (gradB, curlb,
         // bdotcurlb) and causes problems with fourth-order convergence.  We can, however,
         // at least correct the non-differentiated quantities.
#if 0
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);

            FluxBox & thisB = a_BField[dit];
            FluxBox & thisBmag = a_BFieldMag[dit];
            FluxBox & thisbunit = a_BFieldDir[dit];
            FluxBox & thisgradB = a_gradBFieldMag[dit];
            FluxBox & thiscurlB = a_curlBFieldDir[dit];
            FluxBox & thisbdotcurlb = a_BFieldDirdotcurlBFieldDir[dit];

            // Given the just corrected field, recompute the other field quantities (setting
            // the last argument to true prevents the field itself from being recomputed).
            for (int dir=0; dir<SpaceDim; ++dir) {
               coord_sys.computeFieldData(dir,
                                          thisB[dir],
                                          thisBmag[dir],
                                          thisbunit[dir],
                                          thisgradB[dir],
                                          thiscurlB[dir],
                                          thisbdotcurlb[dir],
                                          true);
            }
         }
#else
         // We can at least correct the non-differentiated field quantities

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            FluxBox & thisB = a_BField[dit];
            FluxBox & thisBmag = a_BFieldMag[dit];
            FluxBox & thisbunit = a_BFieldDir[dit];

            for (int dir=0; dir<SpaceDim; ++dir) {
               FArrayBox& thisB_dir = thisB[dir];
               FArrayBox& thisBmag_dir = thisBmag[dir];
               FArrayBox& thisbunit_dir = thisbunit[dir];

               for (BoxIterator bit(thisB_dir.box()); bit.ok(); ++bit) {
                  IntVect iv = bit();

                  double Bmag = 0.;
                  for (int comp=0; comp<3; ++comp) {
                     Bmag += thisB_dir(iv,comp) * thisB_dir(iv,comp);
                  }
                  Bmag = sqrt(Bmag);

                  thisBmag_dir(iv,0) = Bmag;
                  for (int comp=0; comp<3; ++comp) {
                     thisbunit_dir(iv,comp) = thisB_dir(iv,comp) / Bmag;
                  }
               }
            }
         }
#endif

         // The field (and therefore its magnitude and unit vector) are already consistent
         // across block boundaries.  Do the same for the remaining quantities.
         averageAtBlockBoundaries( a_gradBFieldMag );
         averageAtBlockBoundaries( a_curlBFieldDir );
         averageAtBlockBoundaries( a_BFieldDirdotcurlBFieldDir );
         
#ifdef TEST_DIVERGENCE_CLEANING
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            FluxBox& this_change = BFieldMagChange[dit];

            FluxBox denom(this_change.box(), this_change.nComp());
            denom.copy(this_change);
            this_change -= a_BFieldMag[dit];
            for (int dir=0; dir<SpaceDim; ++dir) {
               this_change[dir] /= denom[dir];
            }
         }

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            BFieldMagChange[dit] -= a_BFieldMag[dit];
            BFieldDirChange[dit] -= a_BFieldDir[dit];
            gradBFieldMagChange[dit] -= a_gradBFieldMag[dit];
            curlBFieldDirChange[dit] -= a_curlBFieldDir[dit];
            BFieldDirdotcurlBFieldDirChange[dit] -= a_BFieldDirdotcurlBFieldDir[dit];
         }

         plotFaceData( string("BFieldChange"), m_BFieldCorrection, 0. );
         plotFaceData( string("BFieldMagChange"), BFieldMagChange, 0. );
         plotFaceData( string("BFieldDirChange"), BFieldDirChange, 0. );
         plotFaceData( string("gradBFieldMagChange"), gradBFieldMagChange, 0. );
         plotFaceData( string("curlBFieldDirChange"), curlBFieldDirChange, 0. );
         plotFaceData( string("BFieldDirdotcurlBFieldDirChange"), BFieldDirdotcurlBFieldDirChange, 0. );

         plotFaceData( string("after_correction_BField"), a_BField, 0. );
         plotFaceData( string("after_correction_BFieldMag"), a_BFieldMag, 0. );
         plotFaceData( string("after_correction_BFieldDir"), a_BFieldDir, 0. );
         plotFaceData( string("after_correction_gradB"), a_gradBFieldMag, 0. );
         plotFaceData( string("after_correction_curlBFieldDir"), a_curlBFieldDir, 0. );
         plotFaceData( string("after_correction_BFieldDirdotcurlBFieldDir"), a_BFieldDirdotcurlBFieldDir, 0. );
#endif

         if( procID() == 0 ) {
            cout << endl;
         }
      }
   }
}



void MagGeom::getMagneticFlux( LevelData<FArrayBox>& a_psi ) const
{
   const DisjointBoxLayout& grids = a_psi.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      FArrayBox& this_psi = a_psi[dit];
      
      FArrayBox phys_coords(this_psi.box(), SpaceDim);

      block_coord_sys.getCellCenteredRealCoords(phys_coords);
      block_coord_sys.getMagneticFlux(phys_coords, this_psi);
   }
}



void MagGeom::getMagneticFlux( LevelData<FluxBox>& a_psi ) const
{
   const DisjointBoxLayout& grids = a_psi.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      FluxBox& this_psi = a_psi[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_psi_dir = this_psi[dir];
         FArrayBox phys_coords(this_psi_dir.box(), SpaceDim);

         block_coord_sys.getFaceCenteredRealCoords(dir, phys_coords);
         block_coord_sys.getMagneticFlux(phys_coords, this_psi_dir);
      }
   }
}



void MagGeom::plotCellData( const string&               a_file_name,  
                            const LevelData<FArrayBox>& a_data,
                            const double&               a_time ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   Box domain_box = grids.physDomain().domainBox();
   domain_box.grow(a_data.ghostVect());
   WriteMappedUGHDF5(a_file_name.c_str(), grids, a_data, *m_coord_sys, domain_box, a_time);
}


#if 0
void MagGeom::plotFaceData( const string&             a_file_name,  
                            const LevelData<FluxBox>& a_data,
                            const double&             a_time ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   // Since the passed data is a LevelData<FluxBox>, we average the face
   // values to cell centers
   LevelData<FArrayBox> tmp(grids, SpaceDim*a_data.nComp(), a_data.ghostVect());
   EdgeToCell(a_data, tmp);

   // Next, average the directions
   LevelData<FArrayBox> data_cell(grids, a_data.nComp(), a_data.ghostVect());

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox& this_tmp = tmp[dit];
      FArrayBox& this_data_cell = data_cell[dit];

      this_data_cell.setVal(0.);
      for (int comp=0; comp<this_data_cell.nComp(); comp++) {
        for (int dir=0; dir<SpaceDim; dir++) {
          this_data_cell.plus(this_tmp,comp*SpaceDim + dir,comp,1);
        }
      }
      this_data_cell.mult(1./SpaceDim);
   }

   Box domain_box = grids.physDomain().domainBox();
   domain_box.grow(a_data.ghostVect());
   WriteMappedUGHDF5(a_file_name.c_str(), grids, data_cell, *m_coord_sys, domain_box, a_time);
}
#endif

void MagGeom::plotFaceData( const string&             a_file_name,
                            const LevelData<FluxBox>& a_data,
                            const double&             a_time ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   LevelData<FArrayBox> data_cell(grids, a_data.nComp(), a_data.ghostVect());

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const FluxBox& this_data = a_data[dit];
      FArrayBox& this_data_cell = data_cell[dit];
      this_data_cell.setVal(0.);
      for (int dir=0; dir<SpaceDim; dir++) {
         const FArrayBox& this_data_dir = this_data[dir];
         FArrayBox shift_fab(this_data_dir.box(), this_data_dir.nComp());
         shift_fab.copy(this_data_dir);

         shift_fab.shiftHalf(dir,-1);
         this_data_cell.plus(shift_fab);
         shift_fab.shiftHalf(dir,2);
         this_data_cell.plus(shift_fab);
      }
      this_data_cell.mult(0.5/SpaceDim);
   }

   Box domain_box = grids.physDomain().domainBox();
   //domain_box.grow(a_data.ghostVect());
   WriteMappedUGHDF5(a_file_name.c_str(), grids, data_cell, *m_coord_sys, domain_box, a_time);
}


void
MagGeom::writeGeometryData( const DisjointBoxLayout& a_grids,
                            const double&            a_time ) const
{
  IntVect geom_data_ghosts = 4*IntVect::Unit;
#if CFG_DIM==2
  LevelData<FArrayBox> geom_data(a_grids, 6, geom_data_ghosts);

  for (DataIterator dit(a_grids); dit.ok(); ++dit) {
    const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(a_grids[dit]);

    RealVect dx = block_coord_sys.dx();
    RealVect offset = 0.5*RealVect::Unit;
    offset *= dx;

    FArrayBox& this_geom_data = geom_data[dit];
    BoxIterator bit(geom_data[dit].box());
    for (bit.begin();bit.ok();++bit) {
      IntVect iv = bit();
      RealVect xi = dx*iv + offset;
      RealVect X = block_coord_sys.realCoord(xi);
      this_geom_data(iv,0) = X[RADIAL_DIR];
      this_geom_data(iv,1) = X[POLOIDAL_DIR];

      this_geom_data(iv,2) = block_coord_sys.dXdXi(xi,0,0);
      this_geom_data(iv,3) = block_coord_sys.dXdXi(xi,0,1);
      this_geom_data(iv,4) = block_coord_sys.dXdXi(xi,1,0);
      this_geom_data(iv,5) = block_coord_sys.dXdXi(xi,1,1);
    }
  }

#elif CFG_DIM==3

  LevelData<FArrayBox> geom_data(a_grids, 12, geom_data_ghosts);

  LevelData<FArrayBox> geom_data_cyl(a_grids, 9, geom_data_ghosts);

  for (DataIterator dit(a_grids); dit.ok(); ++dit) {
    const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(a_grids[dit]);

    RealVect dx = block_coord_sys.dx();
    RealVect offset = 0.5*RealVect::Unit;
    offset *= dx;

    FArrayBox& this_geom_data = geom_data[dit];
    BoxIterator bit(geom_data[dit].box());
    for (bit.begin();bit.ok();++bit) {
      IntVect iv = bit();
      RealVect xi = dx*iv + offset;
      RealVect X = block_coord_sys.realCoord(xi);
      this_geom_data(iv,0) = X[RADIAL_DIR];
      this_geom_data(iv,1) = X[TOROIDAL_DIR];
      this_geom_data(iv,2) = X[POLOIDAL_DIR];

      this_geom_data(iv,3) = block_coord_sys.dXdXi(xi,0,0);
      this_geom_data(iv,4) = block_coord_sys.dXdXi(xi,0,1);
      this_geom_data(iv,5) = block_coord_sys.dXdXi(xi,0,2);
      this_geom_data(iv,6) = block_coord_sys.dXdXi(xi,1,0);
      this_geom_data(iv,7) = block_coord_sys.dXdXi(xi,1,1);
      this_geom_data(iv,8) = block_coord_sys.dXdXi(xi,1,2);
      this_geom_data(iv,9) = block_coord_sys.dXdXi(xi,2,0);
      this_geom_data(iv,10) = block_coord_sys.dXdXi(xi,2,1);
      this_geom_data(iv,11) = block_coord_sys.dXdXi(xi,2,2);

      RealVect mappedLoc = dx*iv;
      mappedLoc += offset;

      double R = sqrt(X[0]*X[0] + X[1]*X[1]);
      double Phi = mappedLoc[TOROIDAL_DIR];

      for (int i=0; i<3; ++i) {
         geom_data_cyl[dit](iv,3*i) =           cos(Phi) * this_geom_data(iv,i+3) +        sin(Phi) * this_geom_data(iv,i+6);
         geom_data_cyl[dit](iv,3*i+1) = -(1./R)*sin(Phi) * this_geom_data(iv,i+3) + (1./R)*cos(Phi) * this_geom_data(iv,i+6);
         geom_data_cyl[dit](iv,3*i+2) =                    this_geom_data(iv,i+9);
      }
    }
  }

#endif

  Box domain_box = a_grids.physDomain().domainBox();
  domain_box.grow(geom_data_ghosts);

  //  WriteMappedUGHDF5("config_geom_data_cyl", a_grids, geom_data_cyl, *m_coord_sys, domain_box, a_time);

  WriteMappedUGHDF5("config_geom_data", a_grids, geom_data, *m_coord_sys, domain_box, a_time);

  setCellVolumes(a_grids, geom_data_ghosts);
  WriteMappedUGHDF5("volume", a_grids, m_cell_volume, *m_coord_sys, domain_box, a_time);

  setJ(a_grids, geom_data_ghosts);
  WriteMappedUGHDF5("J", a_grids, m_J, *m_coord_sys, domain_box, a_time);

}



void
MagGeom::checkMultiblockMappingConsistency() const
{
   if ( m_coord_sys->numBlocks() > 1 ) {

      const DisjointBoxLayout& grids = gridsFull();

      LevelData<FluxBox> coords(grids, SpaceDim, IntVect::Zero);
      LevelData<FluxBox> metrics(grids, (SpaceDim-1)*SpaceDim, IntVect::Zero);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& box = grids[dit];
         int block_number = m_coord_sys->whichBlock(box);
         const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
         const Box domain_box = block_coord_sys.domain().domainBox();
         RealVect dx = block_coord_sys.dx();

         for (int dir=0; dir<SpaceDim; ++dir) {

            RealVect offset = 0.5*RealVect::Unit;
            offset *= dx;
            offset[dir] = 0.;

            FArrayBox& this_coord = coords[dit][dir];
            FArrayBox& this_metric = metrics[dit][dir];

            if ( (box.smallEnd(dir) == domain_box.smallEnd(dir) &&
                  !m_coord_sys->containsPhysicalBoundary(block_number, dir, Side::LoHiSide::Lo)) ||
                 (box.bigEnd(dir) == domain_box.bigEnd(dir) &&
                  !m_coord_sys->containsPhysicalBoundary(block_number, dir, Side::LoHiSide::Hi)) ) {

               for (BoxIterator bit(this_coord.box()); bit.ok(); ++bit) {
                  IntVect iv = bit();

                  RealVect xi = dx*iv + offset;
                  RealVect X = block_coord_sys.realCoord(xi);
                  for (int n=0; n<SpaceDim; ++n) {
                     this_coord(iv,n) = X[n];
                  }

                  // only makes sense to check transverse derivatives
                  for (int tdir=0,m=0; tdir<SpaceDim; ++tdir) {
                     if ( tdir != dir ) {
                        for (int n=0; n<SpaceDim; ++n) {
                           this_metric(iv,n+m) = block_coord_sys.dXdXi(xi,n,tdir);
                        }
                        m += SpaceDim;
                     }
                  }
               }
            }
         }
      }

      double max_coord = maxBlockBoundaryDifference(coords);
      double max_metrics = maxBlockBoundaryDifference(metrics);

      if (procID()==0) {
         cout << "Max difference at block interfaces: " << max_coord << " (coordinates), "
              << max_metrics << " (pointwise metrics)" << endl;
      }
   }
}



double
MagGeom::maxBlockBoundaryDifference(LevelData<FluxBox>& a_data) const
{
   double max_diff = 0.;

   if ( m_coord_sys->numBlocks() > 1 ) {

      const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

      BlockRegister blockRegister(m_coordSysPtr, grids, 0);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int idir = 0; idir < SpaceDim; idir++) {
            for (SideIterator sit; sit.ok(); ++sit) {
               Side::LoHiSide side = sit();
               if (blockRegister.hasInterface(dit(), idir, side)) {
                  FArrayBox flux_comp(a_data[dit][idir].box(), a_data.nComp());
                  flux_comp.copy(a_data[dit][idir]);
                  blockRegister.storeFlux(flux_comp, dit(), idir, side);
               }
            }
         }
      }
      blockRegister.close();

      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = m_coord_sys->boundaries();

      double local_max = 0.;

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& baseBox = grids[dit];
         int block_number = m_coord_sys->whichBlock(grids[dit]);
         int faceID = 0;
         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
            for (int idir = 0; idir < SpaceDim; idir++) {
               if (blockRegister.hasInterface(dit(), idir, side)) {
                  // maybe better if this is done inside BlockRegister
                  const BlockBoundary& bb = boundaries[block_number][faceID];
                  int reorientFace = bb.reorientFace(idir);
                  Box faceBox = adjCellBox(baseBox, idir, side, 1);
                  // if Lo, then shift +1; if Hi, then shift -1
                  faceBox.shiftHalf(idir, -sign(side));
                  Side::LoHiSide sideOther = flip(side);
                  // Need to define these FABs.
                  FArrayBox fluxThisFab(faceBox, a_data.nComp());
                  FArrayBox fluxOtherFab(faceBox, a_data.nComp());
                  blockRegister.getFlux(fluxThisFab, dit(),
                                        idir, side, side);
                  fluxThisFab.mult(reorientFace);
                  blockRegister.getFlux(fluxOtherFab, dit(),
                                        idir, side, sideOther);
                  fluxThisFab -= fluxOtherFab;
                  fluxThisFab.abs();
                  double this_max = fluxThisFab.max();
                  if (this_max > local_max) local_max = this_max;
               }
               faceID++;
            } // iterate over dimensions
         } // iterate over sides
      }

#ifdef CH_MPI
      MPI_Allreduce(&local_max, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
      max_diff = local_max;
#endif
   }

   return max_diff;
}


double
MagGeom::maxMappingError() const
{
   double local_max = 0.;
   
   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      const Box& box = m_gridsFull[dit];
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(box);

      double diff = coord_sys.mappingError(box);

      if (diff > local_max) local_max = diff;
   }
      
   double max_diff;

#ifdef CH_MPI
   MPI_Allreduce(&local_max, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
   max_diff = local_max;
#endif

   return max_diff;
}


void
MagGeom::plotInitializationData(const ParmParse&  a_pp,
                                const DisjointBoxLayout& a_grids,
                                const double& a_time) const
{
   if (a_pp.contains("plot_geometry_data")) {
      bool plot_geometry_data;
      a_pp.get("plot_geometry_data", plot_geometry_data);
      
      if (plot_geometry_data) {
         writeGeometryData(a_grids, a_time);
      }
      
   }
   
   if (a_pp.contains("plot_magnetic_field_data")) {
      bool plot_magnetic_field_data;
      a_pp.get("plot_magnetic_field_data", plot_magnetic_field_data);
      
      if (plot_magnetic_field_data) {
         plotMagneticFieldData(a_time);
      }
      
   }
   
   if (a_pp.contains("plot_field_alignment")) {
      bool plot_field_alignment;
      a_pp.get("plot_field_alignment", plot_field_alignment);
      
      if (plot_field_alignment) {
         // Plot the angle (in degrees) between the poloidal field component and the
         // mapped coordinate system vector in the poloidal direction.
         plotFieldAlignment(a_time);
      }
   }
   
   if (a_pp.contains("plot_field_divergence")) {
      bool plot_field_divergence;
      a_pp.get("plot_field_divergence", plot_field_divergence);
      
      if (plot_field_divergence) {
         plotFieldDivergence(m_BField_fc, a_time);
      }
   }
}


void
MagGeom::plotMagneticFieldData(const double& a_time) const
{
   plotCellData( string("BField"), m_BField_cc, a_time );
   plotCellData( string("BFieldMag"), m_BFieldMag_cc, a_time );
   plotCellData( string("BFieldDir"), m_BFieldDir_cc, a_time );
   plotCellData( string("gradB"), m_gradBFieldMag_cc, a_time );
   plotCellData( string("curlBFieldDir"), m_curlBFieldDir_cc, a_time );
   plotCellData( string("BFieldDirdotcurlBFieldDir"), m_BFieldDirdotcurlBFieldDir_cc, a_time );
}

void
MagGeom::plotFieldAlignment(const double& a_time) const
{
   const DisjointBoxLayout& grids = gridsFull();
   LevelData<FluxBox> dotprod(grids, 1, IntVect::Zero);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
      const FluxBox& this_b = m_BFieldDir_fc[dit];
      FluxBox& this_dotprod = dotprod[dit];

      for (int dir=0; dir<SpaceDim; ++dir) {
         const FArrayBox& this_b_dir = this_b[dir];
         FArrayBox& this_dotprod_dir = this_dotprod[dir];
         const Box& box = this_dotprod_dir.box();
         FArrayBox xi(box,SpaceDim);
         coord_sys.getFaceCenteredMappedCoords(dir, xi);

         FArrayBox dXdtheta(box,2);

         coord_sys.dXdXi(dXdtheta, xi, 0, 0, POLOIDAL_DIR, box);
         coord_sys.dXdXi(dXdtheta, xi, 1, 1, POLOIDAL_DIR, box);

         BoxIterator bit(box);
         for (bit.begin();bit.ok();++bit) {
            IntVect iv = bit();
               
            double num = -(dXdtheta(iv,0)*this_b_dir(iv,0) + dXdtheta(iv,1)*this_b_dir(iv,2));
            double norm1 = sqrt( dXdtheta(iv,0) * dXdtheta(iv,0) + dXdtheta(iv,1) * dXdtheta(iv,1));
            double norm2 = sqrt(this_b_dir(iv,0) * this_b_dir(iv,0) + this_b_dir(iv,2) * this_b_dir(iv,2));

            if (norm1 * norm2 == 0.) {
               MayDay::Warning("MagGeom::plotFieldAlignment(): Zero divisor");
               return;
            }

            this_dotprod_dir(iv,0) = acos(num / (norm1 * norm2)) * 180. / Pi;
         }
      }
   }

   plotFaceData( string("field_alignment"), dotprod, a_time );
}



void
MagGeom::plotFieldDivergence( const LevelData<FluxBox>& a_field,
                              const double              a_time ) const
{
   const DisjointBoxLayout& grids = gridsFull();

   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);
   LevelData<FArrayBox> divergence(grids, 1, IntVect::Zero);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FluxBox& this_B = a_field[dit];
      FluxBox& this_flux = flux[dit];

      for (int dir=0; dir<SpaceDim; ++dir) {
         const FArrayBox& this_B_dir = this_B[dir];
         FArrayBox& this_flux_dir = this_flux[dir];

         this_flux_dir.copy( this_B_dir, 0, 0, 1 );
         this_flux_dir.copy( this_B_dir, 2, 1, 1 );
      }
   }

   applyAxisymmetricCorrection( flux );

   fourthOrderAverage(flux);

   computeMappedGridDivergence( flux, divergence, true );

   //LevelData<FArrayBox> cell_volume(grids, 1, IntVect::Zero);
   //getCellVolumes(cell_volume);
   setCellVolumes(grids, IntVect::Zero);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      divergence[dit] /= m_cell_volume[dit];
   }
         
   plotCellData( string("field_divergence"), divergence, a_time );
}


void
MagGeom::fillInternalGhosts( LevelData<FArrayBox>& a_data ) const
{
   const IntVect& nghosts = a_data.ghostVect();

   if ( nghosts.lexGT(IntVect::Zero) ) {
   
      if ( m_mblexPtr && !m_extrablock_exchange && !m_sheared_mb_geom ) {

         const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
         LevelData<FArrayBox> tmp(grids, a_data.nComp(), m_ghosts*IntVect::Unit);
      
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            tmp[dit].copy(a_data[dit]);
         }

         m_mblexPtr->interpGhosts(tmp);

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            a_data[dit].copy(tmp[dit]);
         }
      
         a_data.exchange();
      }
      else if ( m_extrablock_exchange ) {

         // Fill the internal ghosts
         a_data.exchange();

#ifdef NEW_EXTRABLOCK_EXCHANGE
         IntVect include_dir = IntVect::Unit;
#if CFG_DIM == 3
         if ( m_sheared_mb_geom ) {
            // Don't fill extrablock ghosts in the toroidal direction, since this will be
            // handled by interpolateFromShearedGhosts() below
            include_dir[TOROIDAL_DIR] = 0;
         }
#endif
         // Fill the codim ghost cells
         for (int codim=1; codim<=SpaceDim; ++codim) {
            exchangeExtraBlockGhosts(a_data, codim, include_dir);
         }

#if CFG_DIM==3
         if ( m_sheared_mb_geom && nghosts[TOROIDAL_DIR] > 0 ) {
            // Fill the codim 1 extrablock ghost cells at the toroidal block boundaries
            interpolateFromShearedGhosts(a_data);

            // Fill the codim 2 and 3 extrablock ghost cells at the toroidal block boundaries
            IntVect toroidal_dir_only = IntVect::Zero;
            toroidal_dir_only[TOROIDAL_DIR] = 1;

            for (int codim=2; codim<=SpaceDim; ++codim) {
               exchangeExtraBlockGhosts(a_data, codim, toroidal_dir_only);
            }
         }
#endif

         // Fill the corner ghost cells
         const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
         CoDimCopyManager<FArrayBox> manager(grids, grids, nghosts, true);
         manager.manageExchanges(a_data);
#else
         exchangeExtraBlockGhosts(a_data);

         //Need to call again to fill corners
         exchangeExtraBlockGhosts(a_data);
#endif      

#if CFG_DIM==3
         if ( m_sheared_mb_geom && nghosts[TOROIDAL_DIR] > 0 ) {
            interpolateFromShearedGhosts(a_data);
         }
#endif
      }
#if CFG_DIM ==3 
      else if ( m_sheared_mb_geom ) {
         if ( nghosts[TOROIDAL_DIR] > 0 ) {
            interpolateFromShearedGhosts(a_data);
         }
         a_data.exchange();

         //perform data exchange in corner ghost cells
         //might be needed for high-order calculations 
         const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

         CoDimCopyManager<FArrayBox> manager(grids, grids, nghosts, true);
         manager.manageExchanges(a_data);
      }
#endif
      else {
         a_data.exchange();
      }
   }
}


void
MagGeom::exchangeExtraBlockGhosts( LevelData<FArrayBox>& a_data ) const
{
    /*
     Fills extrablock ghost cells using a pure exchange from valid data
     rather than interpolation (for experimental purposes only).
     Presenly only works for objects with a non-zero ghost-layer, thus                                                                                                                   
     the test below.  
     */
    
    IntVect ghost_vect = a_data.ghostVect();
    int max_ghost = -1;
    for (int dir=0; dir<SpaceDim; ++dir) {
      if (ghost_vect[dir] > max_ghost) max_ghost = ghost_vect[dir];
    }

    if ( m_coord_sys->numBlocks() > 1 && max_ghost > 0 ) {
        
        const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
        
        int ncomp = a_data.nComp();
        
        BlockRegister blockRegister(m_coordSysPtr, grids, max_ghost);
        
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            for (int dir = 0; dir < SpaceDim; dir++) {
                IntVect grow_vect = ghost_vect;
                grow_vect[dir] = 0;
                int nghost = ghost_vect[dir];
                for (SideIterator sit; sit.ok(); ++sit) {
                    Side::LoHiSide side = sit();
                    if (blockRegister.hasInterface(dit(), dir, side)) {
#if 1
                       Box fill_box = adjCellBox(grids[dit], dir, side, 1);
                        fill_box.grow(grow_vect);
                        FArrayBox temp(fill_box, ncomp*nghost);
                        for (int n=0; n<nghost; ++n) {
                            temp.shift(dir, -sign(side));
                            temp.copy(a_data[dit],0,n*ncomp,ncomp);
                        }
                        temp.shiftHalf(dir, (2*nghost-1)*sign(side));
                        blockRegister.storeAux(temp, dit(), dir, side);
#else
                        Box fill_box = adjCellBox(grids[dit], dir, side, -nghost);
                        fill_box.grow(grow_vect);
                        FArrayBox temp(fill_box, ncomp);
                        temp.copy(a_data[dit],0,0,ncomp);
                        temp.shiftHalf(dir, sign(side));
                        blockRegister.storeAux(temp, dit(), dir, side);
#endif
                    }
                }
            }
        }
        blockRegister.close();
        
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            for (int dir = 0; dir < SpaceDim; dir++) {
                IntVect grow_vect = ghost_vect;
                grow_vect[dir] = 0;
                int nghost = ghost_vect[dir];
                for (SideIterator sit; sit.ok(); ++sit) {
                    Side::LoHiSide side = sit();
                    if (blockRegister.hasInterface(dit(), dir, side)) {
#if 1
                       Box fill_box = adjCellBox(grids[dit], dir, side, -1);
                        fill_box.grow(grow_vect);
                        FArrayBox temp(fill_box, ncomp*nghost);
                        temp.shiftHalf(dir, sign(side));
                        blockRegister.getAux(temp, dit(),
                                             dir, side, side);
                        temp.shiftHalf(dir, -sign(side));
                        for (int n=0; n<nghost; ++n) {
                            temp.shift(dir, sign(side));
                            a_data[dit].copy(temp,n*ncomp,0,ncomp);
                        }
#else
                  Box fill_box = adjCellBox(grids[dit], dir, side, nghost);
                  fill_box.grow(grow_vect);
                  FArrayBox temp(fill_box, ncomp);
                  temp.shiftHalf(dir, -sign(side));
                  blockRegister.getAux(temp, dit(),
                                       dir, side, side);
                  temp.shiftHalf(dir, sign(side));
                  a_data[dit].copy(temp,0,0,ncomp);
#endif
                    }
                } // iterate over dimensions
            } // iterate over sides
        }
    }
}


void
MagGeom::exchangeExtraBlockGhosts( LevelData<FArrayBox>& a_data,
                                   const int             a_codim,
                                   const IntVect&        a_include_dir ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   // Get box layout for extrablock ghost cells
   BoxLayout ghostLayout;
   LayoutData<DataIndex> index_map;
   LayoutData<IntVect> shift_map;
   getExtrablockExchangeLayout(a_codim,
                               grids,
                               a_data.ghostVect(),
                               a_include_dir,
                               ghostLayout,
                               index_map,
                               shift_map);
   
   // Fill extrablock ghost cells
   BoxLayoutData<FArrayBox> ghost_data(ghostLayout, a_data.nComp());

   if ( a_codim == 1 ) {
      a_data.copyTo(ghost_data);
   }
   else {
      Copier copier;
      const ProblemDomain domain;
      
      copier.define(ghostLayout,
                    grids,
                    domain,
                    a_data.ghostVect(),
                    false,
                    IntVect::Zero);
      copier.reverse();

      a_data.copyTo(ghost_data, copier);
   }
   
   // Fill ghosts of a_data
   for (DataIterator dit(ghostLayout.dataIterator()); dit.ok(); ++dit) {

      FArrayBox this_ghosts(Interval(0,ghost_data.nComp()-1), ghost_data[dit]);
      this_ghosts.shift(-shift_map[dit]);

      CH_assert( (a_data[index_map[dit]].box()).contains(this_ghosts.box()));
      
      a_data[index_map[dit]].copy(this_ghosts);
   }
}


void
MagGeom::getExtrablockExchangeLayout( const int                 a_codim,
                                      const DisjointBoxLayout&  a_grids,
                                      const IntVect&            a_ghost_vect,
                                      const IntVect&            a_include_dir,
                                      BoxLayout&                a_ghosts_layout,
                                      LayoutData<DataIndex>&    a_index_map,
                                      LayoutData<IntVect>&      a_shift_map ) const
{
   CH_TIME("PhaseGeom::getExtrablockExchangeLayout");
   CH_assert(a_codim>=1 && a_codim<=3);

   // Define BoxLayout and the iterator map for extrablock ghost cells

   Vector<Box> boxes; Vector<int> proc_ids;
   std::map<std::tuple<Box,IntVect>,LayoutIndex> local_data_map;
   typedef std::map<std::tuple<Box,IntVect>,LayoutIndex>::value_type dataMapType;

   std::list<std::tuple<Box,IntVect,LayoutIndex> > local_list;
   
   for (LayoutIterator lit(a_grids.layoutIterator()); lit.ok(); ++lit) {  // NB: This loops over all boxes, not just those on this processor

      int block_number = m_coord_sys->whichBlock(a_grids[lit]);
      const Vector< Tuple <BlockBoundary, 2*SpaceDim> >& boundaries = m_coord_sys->boundaries();
      const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = boundaries[block_number];

      const MagBlockCoordSys& block_coords = getBlockCoordSys(a_grids[lit]);
      const ProblemDomain& domain = block_coords.domain();
      const Box& domain_box = domain.domainBox();

      for (int dir1=0; dir1<SpaceDim; ++dir1) {
         if ( a_include_dir[dir1] != 0 ) {
            for (SideIterator sit1; sit1.ok(); ++sit1) {
               Side::LoHiSide side1 = sit1();
           
               if (((side1 == Side::LoHiSide::Lo && a_grids[lit].smallEnd(dir1) == domain_box.smallEnd(dir1))
                    || (side1 == Side::LoHiSide::Hi && a_grids[lit].bigEnd(dir1) == domain_box.bigEnd(dir1)))
                   && !(this_block_boundaries[dir1 + side1*SpaceDim].isDomainBoundary())) {
            
                  Box codim_1_ghostBox = adjCellBox(a_grids[lit], dir1, side1, a_ghost_vect[dir1]);

                  IntVect translation1 = this_block_boundaries[dir1 + side1*SpaceDim].getTransformation().getTranslation();

                  if ( a_codim == 1 ) {
                     codim_1_ghostBox.shift(translation1);

                     bool local_add_box = (a_grids.procID(lit()) == procID())
                        && local_data_map.insert( dataMapType( make_tuple(codim_1_ghostBox,translation1), lit() )).second;

                     bool add_box;
#ifdef CH_MPI
                     MPI_Allreduce(&local_add_box, &add_box, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
#else
                     add_box = local_add_box;
#endif

                     if ( add_box ) {
                        boxes.push_back( codim_1_ghostBox );
                        proc_ids.push_back( a_grids.procID( lit() ) );
                     }
                  }
                  else {
                     for (int dir2=0; dir2<SpaceDim; ++dir2) {
                        if ( a_include_dir[dir2] != 0 && dir2 != dir1 ) {
                           for (SideIterator sit2; sit2.ok(); ++sit2) {
                              Side::LoHiSide side2 = sit2();
                                 
                              if (((side2 == Side::LoHiSide::Lo && a_grids[lit].smallEnd(dir2) == domain_box.smallEnd(dir2))
                                   || (side2 == Side::LoHiSide::Hi && a_grids[lit].bigEnd(dir2) == domain_box.bigEnd(dir2)))
                                  && !(this_block_boundaries[dir2 + side2*SpaceDim].isDomainBoundary())) {
                                    
                                 Box codim_2_ghostBox = adjCellBox(codim_1_ghostBox, dir2, side2, a_ghost_vect[dir2]);
                                    
                                 IntVect translation2 = this_block_boundaries[dir2 + side2*SpaceDim].getTransformation().getTranslation() + translation1;

                                 if ( a_codim == 2 ) {
                                    codim_2_ghostBox.shift(translation2);
                  
                                    bool local_add_box = (a_grids.procID(lit()) == procID())
                                       && local_data_map.insert( dataMapType( make_tuple(codim_2_ghostBox,translation2), lit() )).second;

                                    bool add_box;
#ifdef CH_MPI
                                    MPI_Allreduce(&local_add_box, &add_box, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
#else
                                    add_box = local_add_box;
#endif
                                    if ( add_box ) {
                                       boxes.push_back( codim_2_ghostBox );
                                       proc_ids.push_back( a_grids.procID( lit() ) );
                                    }
                                 }
                                 else {

                                    for (int dir3=0; dir3<SpaceDim; ++dir3) {
                                       if ( a_include_dir[dir3] != 0 && dir3 != dir1 && dir3 != dir2 ) {
                                          for (SideIterator sit3; sit3.ok(); ++sit3) {
                                             Side::LoHiSide side3 = sit3();
            
                                             if (((side3 == Side::LoHiSide::Lo && a_grids[lit].smallEnd(dir3) == domain_box.smallEnd(dir3))
                                                  || (side3 == Side::LoHiSide::Hi && a_grids[lit].bigEnd(dir3) == domain_box.bigEnd(dir3)))
                                                 && !(this_block_boundaries[dir3 + side3*SpaceDim].isDomainBoundary())) {
                                             
                                                Box codim_3_ghostBox = adjCellBox(codim_2_ghostBox, dir3, side3, a_ghost_vect[dir3]);

                                                IntVect translation3 = this_block_boundaries[dir3 + side3*SpaceDim].getTransformation().getTranslation() + translation2;
                                                codim_3_ghostBox.shift(translation3);

                                                bool local_add_box = (a_grids.procID(lit()) == procID()) 
                                                   && local_data_map.insert( dataMapType( make_tuple(codim_3_ghostBox,translation3), lit() )).second;

                                                bool add_box;
#ifdef CH_MPI
                                                MPI_Allreduce(&local_add_box, &add_box, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
#else
                                                add_box = local_add_box;
#endif
                                                if ( add_box ) {
                                                   boxes.push_back( codim_3_ghostBox );
                                                   proc_ids.push_back( a_grids.procID( lit() ) );
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   a_ghosts_layout.define( boxes, proc_ids );

   // Create the index and shift maps
   a_index_map.define( a_ghosts_layout );
   a_shift_map.define( a_ghosts_layout );

   for (DataIterator dit(a_ghosts_layout.dataIterator()); dit.ok(); ++dit) {

      for ( std::map<std::tuple<Box,IntVect>,LayoutIndex >::iterator it = local_data_map.begin(); it != local_data_map.end(); ) {

         std::tuple<Box,IntVect> t = it->first;         

         const Box& list_box = get<0>(t);
         const IntVect& translation = get<1>(t);
         const LayoutIndex& layout_index = it->second;

         if ( list_box == a_ghosts_layout[dit] ) {
            a_index_map[dit] = DataIndex(layout_index);
            a_shift_map[dit] = translation;

            it = local_data_map.erase(it);
         }
         else {
            ++it;
         }
      }
   }
}


void
MagGeom::fillTransversePhysicalGhosts( LevelData<FArrayBox>& a_data ) const
{
   CH_assert(a_data.ghostVect() >= IntVect::Unit);

   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   a_data.exchange();

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      const ProblemDomain& block_domain = block_coord_sys.domain();
      const int block_number = m_coord_sys->whichBlock(grids[dit]);

      for (int dir=0; dir<SpaceDim; ++dir) {
         IntVect grow_vect = a_data.ghostVect();
         grow_vect[dir] = 0;
         Box interior = grow(grids[dit],grow_vect);

         // If the low or high boundaries in this transverse direction is
         // a block interface, then grow the interior box to trick the
         // extrapolation utility into ignoring those sides
         if ( !m_coord_sys->containsPhysicalBoundary(block_number, dir, Side::LoHiSide::Lo) ) {
            interior.growLo(dir,1);
         }
         if ( !m_coord_sys->containsPhysicalBoundary(block_number, dir, Side::LoHiSide::Hi) ) {
            interior.growHi(dir,1);
         }

         secondOrderTransExtrapAtDomainBdry(a_data[dit],
                                            dir,
                                            interior,
                                            block_domain);
      }
   }
}



void
MagGeom::fillTransversePhysicalGhosts( LevelData<FluxBox>& a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   a_data.exchange();

   IntVect ghostVect = a_data.ghostVect();

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      const ProblemDomain& block_domain = block_coord_sys.domain();
      const int block_number = m_coord_sys->whichBlock(grids[dit]);

      for (int dir=0; dir<SpaceDim; ++dir) {
         Box interior = surroundingNodes(grids[dit],dir);
         interior.grow(dir,ghostVect[dir]);

         for (int tdir=0; tdir<SpaceDim; ++tdir) {
            if (tdir != dir) {  // Transverse directions only

               // If the low or high boundaries in this transverse direction is
               // a block interface, then grow the interior box to trick the
               // extrapolation utility into ignoring those sides
               if ( !m_coord_sys->containsPhysicalBoundary(block_number, tdir, Side::LoHiSide::Lo) ) {
                  interior.growLo(tdir,1);
               }
               if ( !m_coord_sys->containsPhysicalBoundary(block_number, tdir, Side::LoHiSide::Hi) ) {
                  interior.growHi(tdir,1);
               }
            }
         }

         secondOrderTransExtrapAtDomainBdry(a_data[dit][dir],
                                            dir,
                                            interior,
                                            block_domain);
      }
   }
}



void
MagGeom::fillTransverseGhosts( LevelData<FluxBox>& a_data,
                               const bool          a_do_block_exchange ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      const ProblemDomain& block_domain = block_coord_sys.domain();
      const int block_number = m_coord_sys->whichBlock(grids[dit]);

      for (int dir=0; dir<SpaceDim; ++dir) {
         Box interior = surroundingNodes(grids[dit],dir);

         if ( a_do_block_exchange ) {
            for (int tdir=0; tdir<SpaceDim; ++tdir) {
               if (tdir != dir) {  // Transverse directions only

                  // If the low or high boundaries in this transverse direction is
                  // a block interface, then grow the interior box to trick the
                  // extrapolation utility into ignoring those sides
                  if ( !m_coord_sys->containsPhysicalBoundary(block_number, tdir, Side::LoHiSide::Lo) ) {
                     interior.growLo(tdir,1);
                  }
                  if ( !m_coord_sys->containsPhysicalBoundary(block_number, tdir, Side::LoHiSide::Hi) ) {
                     interior.growHi(tdir,1);
                  }
               }
            }
         }

         secondOrderTransExtrapAtDomainBdry(a_data[dit][dir],
                                            dir,
                                            interior,
                                            block_domain);
      }
   }

   a_data.exchange();

   if ( a_do_block_exchange ) {
      exchangeTransverseAtBlockBoundaries(a_data);
   }
}



const MagBlockCoordSys&
MagGeom::getBlockCoordSys(const Box& box) const
{
   int block_number = m_coord_sys->whichBlock(box);
   return getBlockCoordSys(block_number);
}



const MagBlockCoordSys&
MagGeom::getBlockCoordSys(const int block_number) const
{
   return *m_coord_sys->getCoordSys(block_number);
}



void
MagGeom::cellCenter( const LevelData<FluxBox>& a_fc_data,
                     LevelData<FArrayBox>&     a_cc_data ) const
{
   CH_assert(a_fc_data.ghostVect() >= a_cc_data.ghostVect());

   const DisjointBoxLayout& grids = a_cc_data.disjointBoxLayout();

   LevelData<FArrayBox> tmp(grids, SpaceDim*a_fc_data.nComp(), a_fc_data.ghostVect());
   EdgeToCell(a_fc_data, tmp);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FArrayBox& this_tmp = tmp[dit];
      FArrayBox& this_data = a_cc_data[dit];

      this_data.setVal(0.);
      for (int comp=0; comp<a_fc_data.nComp(); comp++) {
        for (int dir=0; dir<SpaceDim; dir++) {
          this_data.plus(this_tmp,comp*SpaceDim + dir,comp,1);
        }
      }
      this_data.mult(1./SpaceDim);
   }
}



void
MagGeom::getIntegralsForVelocity( const LevelData<FArrayBox>& a_Efield_cell,
                                  const LevelData<FArrayBox>& a_phi_node,
                                  const bool                  a_fourth_order,
                                  LevelData<FArrayBox>&       a_nodal_integrals,
                                  LevelData<FluxBox>&         a_face_integrals,
                                  LevelData<FArrayBox>&       a_volume_integrals ) const
{
   CH_assert(SpaceDim==2);
   CH_assert(a_nodal_integrals.ghostVect() == 2*IntVect::Unit);
   CH_assert(a_face_integrals.ghostVect()  == IntVect::Unit);
   CH_assert(a_volume_integrals.ghostVect() == IntVect::Unit);
   CH_assert(a_nodal_integrals.nComp()  == 2);
   CH_assert(a_face_integrals.nComp()   == 2);
   CH_assert(a_volume_integrals.nComp() == 2);

   const DisjointBoxLayout& grids = gridsFull();

   if ( !m_velocity_data_cached ) {

      m_nodal_integrals.define(grids, 2, 2*IntVect::Unit);
      LevelData<FArrayBox> volume_integral_data(grids, 1, IntVect::Unit);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);

         // Get jacobian at cell centers
         FArrayBox& this_J = volume_integral_data[dit];
         Box this_Jbox = this_J.box();

         FArrayBox xi_cell(this_Jbox, SpaceDim);
         coord_sys.getCellCenteredMappedCoords(xi_cell);
         coord_sys.pointwiseJ(this_J, xi_cell, this_Jbox);
      }

      m_coord_sys->postProcessMetricData(volume_integral_data);

      LevelData<FArrayBox> Bmag(grids, 1, 2*IntVect::Unit);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
         FArrayBox& this_nodal_integral_data = m_nodal_integrals[dit];
         const Box& nodal_box = this_nodal_integral_data.box();

         FArrayBox xi_node(nodal_box, SpaceDim);
         coord_sys.getNodeCenteredMappedCoords(xi_node);

         FArrayBox A(nodal_box,3);
         FArrayBox b(nodal_box,3);
         FArrayBox Bmagnitude(nodal_box,1);
         coord_sys.getNodalFieldData(xi_node, A, b, Bmagnitude);
         Bmag[dit].copy(Bmagnitude);

         this_nodal_integral_data.copy(A,1,0,1);
         this_nodal_integral_data.copy(b,1,1,1);
         this_nodal_integral_data.negate();

         if ( coord_sys.isAxisymmetric()) {

            // Get major radius at nodes
            FArrayBox this_majorR_node(nodal_box, 1);
            coord_sys.getNodeCenteredMappedCoords(xi_node);
            coord_sys.getPointwiseMajorRadius(xi_node, this_majorR_node);

            for (int n=0; n<2; ++n) {
               this_nodal_integral_data.mult(this_majorR_node,0,n,1);
            }
            this_nodal_integral_data *= 2.*Pi;
         }
      }

      m_nodal_integrals.exchange();
      volume_integral_data.exchange();

      m_face_bXgradB_integrals.define(grids, 1, IntVect::Unit);
      m_volume_B_integrals.define(grids, 3, IntVect::Unit);
      m_volume_BdotGradB_integrals.define(grids, 1, IntVect::Unit);

      m_face_Binverse_integrals.define(grids, 1, IntVect::Unit);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
         const RealVect& dx = coord_sys.dx();

         FluxBox& this_face_Binverse_integrals = m_face_Binverse_integrals[dit];
         FluxBox& this_face_bXgradB_integrals = m_face_bXgradB_integrals[dit];

         for (int dir=0; dir<SpaceDim; ++dir) {

            FArrayBox& this_face_Binverse_integrals_dir = this_face_Binverse_integrals[dir];

            FORT_COMPUTE_FACE_BINVERSE_INTEGRALS(CHF_BOX(this_face_Binverse_integrals_dir.box()),
                                                 CHF_CONST_INT(dir),
                                                 CHF_REALVECT(dx),
                                                 CHF_CONST_FRA1(m_BFieldMag_fc[dit][dir],0),
                                                 CHF_FRA1(this_face_Binverse_integrals_dir,0));

            FArrayBox& this_face_bXgradB_integrals_dir = this_face_bXgradB_integrals[dir];
            
            FORT_COMPUTE_FACE_BXGRADB_INTEGRALS(CHF_BOX(this_face_bXgradB_integrals_dir.box()),
                                                CHF_CONST_INT(dir),
                                                CHF_CONST_FRA1(Bmag[dit],0),
                                                CHF_FRA1(this_face_bXgradB_integrals_dir,0));

            if ( !coord_sys.isAxisymmetric() ) {
               // Multiply by By (assumed uniform)
               this_face_Binverse_integrals_dir.mult(m_BField_fc[dit][dir], 1, 0, 1);
               this_face_bXgradB_integrals_dir.mult(m_BField_fc[dit][dir], 1, 0, 1);
            }
            else {
               double two_pi_RBtor = 2. * Pi * coord_sys.getRBtoroidal();
               this_face_Binverse_integrals_dir *= two_pi_RBtor;
               this_face_bXgradB_integrals_dir  *= two_pi_RBtor;
            }
         }

         FArrayBox& this_vol_integrals = a_volume_integrals[dit];

         FORT_COMPUTE_VOLUME_B_INTEGRALS(CHF_BOX(this_vol_integrals.box()),
                                         CHF_CONST_REALVECT(dx),
                                         CHF_CONST_FRA(m_BField_cc[dit]),
                                         CHF_CONST_FRA1(volume_integral_data[dit],0),
                                         CHF_FRA(m_volume_B_integrals[dit]));

         FORT_COMPUTE_VOLUME_B_DOT_GRADB_INTEGRALS(CHF_BOX(this_vol_integrals.box()),
                                                   CHF_CONST_REALVECT(dx),
                                                   CHF_CONST_FRA(m_BField_cc[dit]),
                                                   CHF_CONST_FRA(m_gradBFieldMag_cc[dit]),
                                                   CHF_CONST_FRA1(volume_integral_data[dit],0),
                                                   CHF_FRA1(m_volume_BdotGradB_integrals[dit],0));
      }

      if ( a_fourth_order ) {
         // We need to convert to a fourth-order integral here before the multiplication
         // by the average of the the transverse electric field.  Converting the volume
         // integrals to fourth-order is done below after everything is added up.
         fourthOrderAverage(m_face_Binverse_integrals);
      }

      m_face_Binverse_integrals.exchange();

      m_face_bXgradB_integrals.exchange();
      m_volume_B_integrals.exchange();
      m_volume_BdotGradB_integrals.exchange();

      m_velocity_data_cached = true;
   }

   // Construct the b X E integrals
   LevelData<FluxBox> face_bXE_integrals(grids, 1, IntVect::Unit);
   computeBxEIntegrals(a_phi_node, a_fourth_order, face_bXE_integrals);

   LevelData<FArrayBox> exb_divergence(grids, 1, IntVect::Zero);

   const IntVect ghosts( a_volume_integrals.ghostVect());
   setPointwiseNJInverseOnCells(grids, ghosts);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {

      // Copy the cached nodal data
      a_nodal_integrals[dit].copy(m_nodal_integrals[dit]);

      // Copy the b X E integrals and the cached b X gradB integrals
      a_face_integrals[dit].copy(face_bXE_integrals[dit], 0, 0, 1);
      a_face_integrals[dit].copy(m_face_bXgradB_integrals[dit], 0, 1, 1);

      // Construct the B dot E integrals

      FArrayBox& this_vol_integrals = a_volume_integrals[dit];

      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      
      if (!block_coord_sys.isFieldAligned()) {
         FORT_COMPUTE_VOLUME_B_DOT_E_INTEGRALS(CHF_BOX(this_vol_integrals.box()),
                                               CHF_CONST_FRA(a_Efield_cell[dit]),
                                               CHF_CONST_FRA(m_volume_B_integrals[dit]),
                                               CHF_FRA1(this_vol_integrals,0));
      }

      else {
         RealVect dx = block_coord_sys.dx();
         FORT_COMPUTE_VOLUME_B_DOT_E_INTEGRALS_FIELD_ALIGNED(CHF_BOX(this_vol_integrals.box()),
                                                             CHF_CONST_FRA1(a_phi_node[dit],0),
                                                             CHF_CONST_REALVECT(dx),
                                                             CHF_CONST_FRA(m_NJinverse_cell_centered[dit]),
                                                             CHF_CONST_FRA(m_volume_B_integrals[dit]),
                                                             CHF_FRA1(this_vol_integrals,0));
      }

      // Copy the cached B dot gradB integrals
      this_vol_integrals.copy(m_volume_BdotGradB_integrals[dit], 0, 1, 1);
   }

   if ( a_fourth_order ) {
      fourthOrderAverage(a_volume_integrals);
   }

   a_face_integrals.exchange();
   a_volume_integrals.exchange();
}



void
MagGeom::computeBxEIntegrals( const LevelData<FArrayBox>& a_phi_node,
                              bool                      a_fourth_order,
                              LevelData<FluxBox>&       a_bXE_integral ) const
{
   CH_assert(a_phi_node.ghostVect() == 2*IntVect::Unit);
   CH_assert(m_face_Binverse_integrals.isDefined());
   const DisjointBoxLayout& grids = gridsFull();

   LevelData<FluxBox> E_transverse(grids, 1, IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();

      FluxBox& this_E_transverse = E_transverse[dit];
      const FArrayBox& this_phi_node = a_phi_node[dit];

      for (int dir=0; dir<SpaceDim; ++dir) {
         Box face_box = grow(surroundingNodes(grids[dit],dir), 1);

         FORT_COMPUTE_FACE_TRANSVERSE_E(CHF_BOX(face_box),
                                        CHF_CONST_INT(dir),
                                        CHF_CONST_REALVECT(dx),
                                        CHF_CONST_FRA1(this_phi_node,0),
                                        CHF_FRA1(this_E_transverse[dir],0));
      }
   }

   E_transverse.exchange();
   fillTransversePhysicalGhosts(E_transverse);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FluxBox& this_Binverse_integral = m_face_Binverse_integrals[dit];
      FluxBox& this_E_transverse = E_transverse[dit];
      FluxBox& this_bXE_integral = a_bXE_integral[dit];

      if ( a_fourth_order ) {
         fourthOrderMultFace(this_bXE_integral, this_Binverse_integral, this_E_transverse);
      }
      else {
         for (int dir=0; dir<SpaceDim; ++dir) {
            this_bXE_integral[dir].copy(this_Binverse_integral[dir]);
            this_bXE_integral[dir] *= this_E_transverse[dir];
         }
      }
   }

   a_bXE_integral.exchange();
}

void
MagGeom::computeEXBDrift( const LevelData<FluxBox>& a_E_field,
                          LevelData<FluxBox>&       a_ExB_drift ) const
{
   const DisjointBoxLayout& grids = gridsFull();
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {

      
      const FluxBox& this_E_field = a_E_field[dit];
      const FluxBox& this_B_field = m_BField_fc[dit];
      const FluxBox& this_B_mag = m_BFieldMag_fc[dit];
      
      for (int dir=0; dir<SpaceDim; ++dir) {
         
         FArrayBox& this_ExB_drift_dir = a_ExB_drift[dit][dir];

         FORT_COMPUTE_EXB_DRIFT(CHF_BOX(this_ExB_drift_dir.box()),
                                CHF_CONST_FRA(this_E_field[dir]),
                                CHF_CONST_FRA(this_B_field[dir]),
				CHF_CONST_FRA1(this_B_mag[dir],0),
                                CHF_FRA(this_ExB_drift_dir));
      }
   }
}

#if 0
void
MagGeom::computeNormalExBDrift( const LevelData<FArrayBox>&  a_Efield_cell,
                                const LevelData<FArrayBox>&  a_phi_node,
                                const bool                   a_fourth_order_Efield,
                                LevelData<FluxBox>&          a_normal_ExB_drift )
{
   const DisjointBoxLayout& grids = a_normal_ExB_drift.disjointBoxLayout();

   LevelData<FArrayBox> nodal_integrals(grids, 2, 2*IntVect::Unit);
   LevelData<FluxBox>   face_integrals(grids, 2, IntVect::Unit);
   LevelData<FArrayBox> volume_integrals(grids, 2, IntVect::Unit);
   getIntegralsForVelocity(a_Efield_cell, a_phi_node, a_fourth_order_Efield,
                           nodal_integrals, face_integrals, volume_integrals);



}
#endif


void
MagGeom::setRealCoords()
{
   m_cell_centered_real_coords.define(m_gridsFull, SpaceDim, (m_ghosts+1)*IntVect::Unit);
   m_face_centered_real_coords.define(m_gridsFull, SpaceDim, (m_ghosts+1)*IntVect::Unit);
   m_edge_centered_real_coords.define(m_gridsFull, SpaceDim, (m_ghosts+1)*IntVect::Unit);

   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(m_gridsFull[dit]);
      block_coord_sys.getCellCenteredRealCoords(m_cell_centered_real_coords[dit]);
      for ( int dir=0; dir<SpaceDim; ++dir ) {
         block_coord_sys.getFaceCenteredRealCoords(dir, m_face_centered_real_coords[dit][dir]);
         block_coord_sys.getEdgeCenteredRealCoords(dir, m_edge_centered_real_coords[dit][dir]);
      }
   }
}


void
MagGeom::setNormalizedMagneticFlux()
{
   m_normalized_magnetic_flux_cell.define(m_gridsFull, 1, (m_ghosts+1)*IntVect::Unit);
   m_normalized_magnetic_flux_face.define(m_gridsFull, 1, (m_ghosts+1)*IntVect::Unit);

   if ( m_coord_sys->providesFlux() ) {
      for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(m_gridsFull[dit]);
         block_coord_sys.getNormMagneticFlux(m_cell_centered_real_coords[dit], m_normalized_magnetic_flux_cell[dit]);
         for (int dir=0; dir<SpaceDim; ++dir) {
            block_coord_sys.getNormMagneticFlux(m_face_centered_real_coords[dit][dir], m_normalized_magnetic_flux_face[dit][dir]);
         }
      }
   }
   else {
      for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
         m_normalized_magnetic_flux_cell[dit].setVal(0.);
         for (int dir=0; dir<SpaceDim; ++dir) {
            m_normalized_magnetic_flux_face[dit].setVal(0.);
         }
      }
   }
}


void MagGeom::interpolateErFromMagFS(LevelData<FluxBox>&   a_Er_face,
                                     LevelData<FArrayBox>& a_Er_cell) const
{
    //Here, we consider the model geometry, in which Ez along the top cut is the radial component
    //Needs to be generalized for the case of a DIII-D geometry (by extracting the actual radial component)
    //Also needs to be generalized for the 10 block structure
    
    //X-point vicinity, in which we address the dealignemnt
    double Xpt_vicinity = DBL_MAX;
   
    const SingleNullBlockCoordSys& lcore_coord_sys =
       (SingleNullBlockCoordSys&)getBlockCoordSys(SingleNullBlockCoordSys::LCORE);
   
    double R0 = (lcore_coord_sys.getMagAxis())[0];
    double Zaxis = (lcore_coord_sys.getMagAxis())[1];
   
    if (m_model_flux_function) {
      Zaxis = 1.75935;
      R0 = 1.6;
    }
    
    const DisjointBoxLayout& grids = a_Er_cell.disjointBoxLayout();
    
    LevelData<FArrayBox> Ez_cell( grids, 1, IntVect::Zero );
    LevelData<FluxBox> Ez_face( grids, 1, IntVect::Zero );
    
    //Get the Er
    DataIterator dit( a_Er_cell.dataIterator() );
    for (dit.begin(); dit.ok(); ++dit) {
        Box box( Ez_cell[dit].box() );
        Ez_cell[dit].copy(a_Er_cell[dit], box, 2 , box, 0, 1);
        Ez_face[dit].copy(a_Er_face[dit], 2 , 0, 1);
    }
    
    //Interpolate the  top-Er flux function onto the grid
    LevelData<FArrayBox> Ez_cell_interp( grids, 1, IntVect::Zero );
    LevelData<FluxBox> Ez_face_interp( grids, 1, IntVect::Zero );
    interpolateFromMagFS(Ez_face, Ez_cell, Ez_cell_interp);
    interpolateFromMagFS(Ez_face, Ez_face_interp);
    
    
    //Get the poloidal magneitc field
    LevelData<FArrayBox> BpolMag_cell( grids, 1, IntVect::Zero );
    for (dit.begin(); dit.ok(); ++dit) {
        Box box( BpolMag_cell[dit].box() );
        BoxIterator bit(box);
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            BpolMag_cell[dit](iv,0) = sqrt(pow(m_BField_cc[dit](iv,0),2)+pow(m_BField_cc[dit](iv,2),2));
        }
    }
    
    LevelData<FluxBox> BpolMag_face( grids, 1, IntVect::Zero );
    for (dit.begin(); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; dir++) {
            Box box( BpolMag_face[dit][dir].box() );
            BoxIterator bit(box);
            for (bit.begin(); bit.ok(); ++bit) {
                IntVect iv = bit();
                BpolMag_face[dit][dir](iv,0) = sqrt(pow(m_BField_fc[dit][dir](iv,0),2)+pow(m_BField_fc[dit][dir](iv,2),2));
            }
        }
    }
    
    //Interpolate the top-Bp flux function onto the grid
    LevelData<FArrayBox> Bpol_cell_interp( grids, 1, IntVect::Zero );
    LevelData<FluxBox> Bpol_face_interp( grids, 1, IntVect::Zero );
    interpolateFromMagFS(BpolMag_face, BpolMag_cell, Bpol_cell_interp);
    interpolateFromMagFS(BpolMag_face, Bpol_face_interp);
    
    MagCoordSys* const coord_sys_ptr = m_coord_sys.getRefToThePointer();

      //Compute E-field interpolation on face-centers
    for (dit.begin(); dit.ok(); ++dit) {
        int block_number = m_coord_sys->whichBlock(grids[dit]);
        const MagBlockCoordSys& block_coord_sys =  getBlockCoordSys(grids[dit]);
        
        for (int dir = 0; dir < SpaceDim; dir++) {
            FArrayBox& this_E_face = a_Er_face[dit][dir];
            Box box( Ez_face_interp[dit][dir].box() );
            FArrayBox face_center_coords_dir( box, 2 );
            block_coord_sys.getFaceCenteredRealCoords( dir, face_center_coords_dir );
            
            BoxIterator bit(box);
            for (bit.begin(); bit.ok(); ++bit) {
                IntVect iv = bit();
                double delta_R = face_center_coords_dir(iv,0) - ((SingleNullCoordSys* const)coord_sys_ptr)->m_Xpoint[0];
                double delta_Z = face_center_coords_dir(iv,1) - ((SingleNullCoordSys* const)coord_sys_ptr)->m_Xpoint[1];
                double delta_Z_axis = face_center_coords_dir(iv,1) - Zaxis;
                double distance = sqrt(delta_Z*delta_Z + delta_R*delta_R);
                
                if (distance < Xpt_vicinity) {
                    
                    double Bp = BpolMag_face[dit][dir](iv, 0);
                    double BR = abs(m_BField_fc[dit][dir](iv,0));
                    double BZ = abs(m_BField_fc[dit][dir](iv,2));
                    double Bp_top = Bpol_face_interp[dit][dir](iv,0);
                    double Ez_top = Ez_face_interp[dit][dir](iv,0);
                    double R = face_center_coords_dir(iv,0);
                    double signEZ = (delta_Z_axis < 0) ? -1 : (delta_Z_axis > 0);
                    double signER = (delta_R < 0) ? -1 : (delta_R > 0);
                    
                    this_E_face(iv,0) = Ez_top * (Bp * R)/(Bp_top * R0) * BZ/Bp * signER;
                    this_E_face(iv,2) = Ez_top * (Bp * R)/(Bp_top * R0) * BR/Bp * signEZ;
                    
                    if (block_number>SingleNullBlockCoordSys::RCSOL) this_E_face(iv,2) *= -1.0;
                    
                }
            }
        }
    }
    
    //Compute E-field interpolation on cell-centers
    for (dit.begin(); dit.ok(); ++dit) {
        int block_number = m_coord_sys->whichBlock(grids[dit]);
        const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
        
        FArrayBox& this_E_cell = a_Er_cell[dit];
        Box box( Ez_cell_interp[dit].box() );
        FArrayBox cell_center_coords( box, 2 );
        block_coord_sys.getCellCenteredRealCoords( cell_center_coords );
        
        BoxIterator bit(box);
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            double delta_R = cell_center_coords(iv,0) - ((SingleNullCoordSys* const)coord_sys_ptr)->m_Xpoint[0];
            double delta_Z = cell_center_coords(iv,1) - ((SingleNullCoordSys* const)coord_sys_ptr)->m_Xpoint[1];
            double delta_Z_axis = cell_center_coords(iv,1) - Zaxis;
            double distance = sqrt(delta_Z*delta_Z + delta_R*delta_R);
            
            if (distance < Xpt_vicinity) {
                
                double Bp = BpolMag_cell[dit](iv, 0);
                double BR = abs(m_BField_cc[dit](iv,0));
                double BZ = abs(m_BField_cc[dit](iv,2));
                double Bp_top = Bpol_cell_interp[dit](iv,0);
                double Ez_top = Ez_cell_interp[dit](iv,0);
                double R = cell_center_coords(iv,0);
                double signEZ = (delta_Z_axis < 0) ? -1 : (delta_Z_axis > 0);
                double signER = (delta_R < 0) ? -1 : (delta_R > 0);
                
                this_E_cell(iv,0) = Ez_top * (Bp * R)/(Bp_top * R0) * BZ/Bp * signER;
                this_E_cell(iv,2) = Ez_top * (Bp * R)/(Bp_top * R0) * BR/Bp * signEZ;
                
                if (block_number>3) this_E_cell(iv,2) *= -1.0;
            }
        }
    }
    
    a_Er_face.exchange();
    a_Er_cell.exchange();
    
}


void MagGeom::computeMagFluxMappingCell( LevelData<FArrayBox>& a_magFS_mapping_cell) const


{
    const MagBlockCoordSys& lcore_coord_sys = getBlockCoordSys(SingleNullBlockCoordSys::LCORE);
    RealVect dx_core = lcore_coord_sys.dx();
    int nrad_core = lcore_coord_sys.domain().domainBox().size(RADIAL_DIR);
    
    const MagBlockCoordSys& lcsol_coord_sys = getBlockCoordSys(SingleNullBlockCoordSys::LCSOL);
    RealVect dx_sol = lcsol_coord_sys.dx();
    int nrad_sol = lcsol_coord_sys.domain().domainBox().size(RADIAL_DIR);
    
    //Finding cell-ceneters at the top of a tokamak
    int size = nrad_core + nrad_sol;
    double *z_top_faces = new double[size];
    RealVect mapped_coord;
    //Core region
    mapped_coord[1] = lcore_coord_sys.lowerMappedCoordinate(POLOIDAL_DIR);
    for (int i=0; i<nrad_core; ++i) {
        mapped_coord[0] = dx_core[RADIAL_DIR]*(i+0.5);
        z_top_faces[i] = lcore_coord_sys.realCoord(mapped_coord)[1];
    }
    //SOL region
    mapped_coord[1] = lcsol_coord_sys.lowerMappedCoordinate(POLOIDAL_DIR);
    for (int i=nrad_core; i<size; ++i) {
        mapped_coord[0] = lcsol_coord_sys.lowerMappedCoordinate(0) + dx_sol[RADIAL_DIR]*(i - nrad_core + 0.5);
        z_top_faces[i] = lcsol_coord_sys.realCoord(mapped_coord)[1];
    }
    
    //Computing mapping data that relates the input cell-data to the magnetic-flux surfaces
    const DisjointBoxLayout& grids = a_magFS_mapping_cell.disjointBoxLayout();
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        
        a_magFS_mapping_cell[dit].setVal(0.0);
        FArrayBox& this_mapping_cell = a_magFS_mapping_cell[dit];
        const MagBlockCoordSys& block_coord_sys =  getBlockCoordSys(grids[dit]);
        
        Box box( this_mapping_cell.box() );
        
        FArrayBox cell_center_coords( box, SpaceDim );
        block_coord_sys.getCellCenteredRealCoords( cell_center_coords );
        
        for (BoxIterator bit(box); bit.ok(); ++bit) {
            IntVect iv = bit();
            RealVect this_coord;
            this_coord[0] = cell_center_coords(iv,0);
            this_coord[1] = cell_center_coords(iv,1);
            
            //Computing the corresponding flux surface (here, given by the z-coordinate at a tokamak top)
            this_mapping_cell(iv,0) = getMagFS(this_coord);
            
            //Computing the closest upper-bound cell-center and interpolating coefficients
            if ( getMagFS(this_coord) > z_top_faces[size-1]) {
                this_mapping_cell(iv,1) = size;
                this_mapping_cell(iv,2) = (getMagFS(this_coord) - z_top_faces[size-1]) /(z_top_faces[size-1]-z_top_faces[size-2]);
            }
            
            else if ( getMagFS(this_coord) < z_top_faces[0]) {
                this_mapping_cell(iv,1) = 0;
                this_mapping_cell(iv,2) = (z_top_faces[0]-getMagFS(this_coord)) / (z_top_faces[1]-z_top_faces[0]);
            }
            
            else {
                int irad=0;
                bool index_found = false;
                while ( !index_found ) {
                    irad++;
                    if (getMagFS(this_coord) < z_top_faces[irad]) {index_found = true;}
                }
                this_mapping_cell(iv,1) = irad;
                this_mapping_cell(iv,2) = (z_top_faces[irad] - getMagFS(this_coord)) / (z_top_faces[irad]-z_top_faces[irad-1]);
            }
        }
    }
    
    delete [] z_top_faces;
    
}


void MagGeom::computeMagFluxMappingFace( LevelData<FluxBox>& a_magFS_mapping_face) const

{
    const MagBlockCoordSys& lcore_coord_sys = getBlockCoordSys(SingleNullBlockCoordSys::LCORE);
    RealVect dx_core = lcore_coord_sys.dx();
    int nrad_core = lcore_coord_sys.domain().domainBox().size(RADIAL_DIR);
    
    const MagBlockCoordSys& lcsol_coord_sys = getBlockCoordSys(SingleNullBlockCoordSys::LCSOL);
    RealVect dx_sol = lcsol_coord_sys.dx();
    int nrad_sol = lcsol_coord_sys.domain().domainBox().size(RADIAL_DIR);
    
    //Finding cell-ceneters at the top of a tokamak
    int size = nrad_core + nrad_sol;
    double *z_top_faces = new double[size];
    RealVect mapped_coord;
    //Core region
    mapped_coord[1] = lcore_coord_sys.lowerMappedCoordinate(POLOIDAL_DIR);
    for (int i=0; i<nrad_core; ++i) {
        mapped_coord[0] = dx_core[RADIAL_DIR]*(i+0.5);
        z_top_faces[i] = lcore_coord_sys.realCoord(mapped_coord)[1];
    }
    //SOL region
    mapped_coord[1] = lcsol_coord_sys.lowerMappedCoordinate(POLOIDAL_DIR);
    for (int i=nrad_core; i<size; ++i) {
        mapped_coord[0] = lcsol_coord_sys.lowerMappedCoordinate(0) + dx_sol[RADIAL_DIR]*(i - nrad_core + 0.5);
        z_top_faces[i] = lcsol_coord_sys.realCoord(mapped_coord)[1];
    }
    
    //Computing mapping data that relates the input cell-data to the magnetic-flux surfaces
    const DisjointBoxLayout& grids = a_magFS_mapping_face.disjointBoxLayout();
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        
        a_magFS_mapping_face[dit].setVal(0.0);
        for (int dir(0); dir<SpaceDim; ++dir) {
            
            FArrayBox& this_mapping_face_dir = a_magFS_mapping_face[dit][dir];
            const MagBlockCoordSys& block_coord_sys =  getBlockCoordSys(grids[dit]);
            
            Box box( this_mapping_face_dir.box() );
            
            FArrayBox face_center_coords_dir( box, SpaceDim );
            block_coord_sys.getFaceCenteredRealCoords( dir, face_center_coords_dir );
            
            for (BoxIterator bit(box); bit.ok(); ++bit) {
                IntVect iv = bit();
                RealVect this_coord;
                this_coord[0] = face_center_coords_dir(iv,0);
                this_coord[1] = face_center_coords_dir(iv,1);
                
                //Computing the corresponding flux surface (here, given by the z-coordinate at a tokamak top)
                this_mapping_face_dir(iv,0) = getMagFS(this_coord);
                
                //Computing the closest upper-bound cell-center and interpolating coefficients
                if ( getMagFS(this_coord) > z_top_faces[size-1]) {
                    this_mapping_face_dir(iv,1) = size;
                    this_mapping_face_dir(iv,2) = (getMagFS(this_coord) - z_top_faces[size-1]) /(z_top_faces[size-1]-z_top_faces[size-2]);
                }
                
                else if ( getMagFS(this_coord) < z_top_faces[0]) {
                    this_mapping_face_dir(iv,1) = 0;
                    this_mapping_face_dir(iv,2) = (z_top_faces[0]-getMagFS(this_coord)) / (z_top_faces[1]-z_top_faces[0]);
                }
                
                else {
                    int irad=0;
                    bool index_found = false;
                    while ( !index_found ) {
                        irad++;
                        if (getMagFS(this_coord) < z_top_faces[irad]) {index_found = true;}
                    }
                    this_mapping_face_dir(iv,1) = irad;
                    this_mapping_face_dir(iv, 2) = (z_top_faces[irad] - getMagFS(this_coord)) / (z_top_faces[irad]-z_top_faces[irad-1]);
                }
            }
        }
    }
    
    delete [] z_top_faces;
    
}

void MagGeom::interpolateFromMagFS( const LevelData<FluxBox>& a_data,
                                   LevelData<FluxBox>& a_interp) const

{
    CH_assert(a_data.nComp() == 1);
    CH_assert(a_data.ghostVect() == IntVect::Zero);
    
    const MagBlockCoordSys& lcore_coord_sys = getBlockCoordSys(SingleNullBlockCoordSys::LCORE);
    int nrad_core  = lcore_coord_sys.domain().domainBox().size(RADIAL_DIR);
    int lo_pol_LCORE = lcore_coord_sys.domain().domainBox().smallEnd(POLOIDAL_DIR);
    
    const MagBlockCoordSys& lcsol_coord_sys = getBlockCoordSys(SingleNullBlockCoordSys::LCSOL);
    int nrad_sol = lcsol_coord_sys.domain().domainBox().size(RADIAL_DIR);
    int lo_pol_LCSOL = lcsol_coord_sys.domain().domainBox().smallEnd(POLOIDAL_DIR);
    int lo_rad_LCSOL = lcsol_coord_sys.domain().domainBox().smallEnd(RADIAL_DIR);
    
    //Finding cell-ceneters at the top of a tokamak
    int size = nrad_core + nrad_sol;
    double *data_Z = new double[size];
    double *data_Z_loc = new double[size];
    
    for (int i=0; i<size; ++i) {
        data_Z[i]=0.0;
        data_Z_loc[i]=0.0;
    }
    
    //Extracting the data along the LCORE boundary (top),
    //and creating the corresponding magneitc flux function
    const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        int block_number = m_coord_sys->whichBlock(grids[dit]);
        if (block_number == SingleNullBlockCoordSys::LCORE) {
            const FArrayBox& this_data_pol_dir = a_data[dit][POLOIDAL_DIR];
            const Box box( this_data_pol_dir.box() );
            for (BoxIterator bit(box); bit.ok(); ++bit) {
                IntVect iv = bit();
                if (iv[POLOIDAL_DIR]==lo_pol_LCORE) {
                    int irad = iv[0];
                    data_Z_loc[irad] = this_data_pol_dir(iv,0);
                }
            }
        }
        
        if (block_number == SingleNullBlockCoordSys::LCSOL) {
            const FArrayBox& this_data_pol_dir = a_data[dit][POLOIDAL_DIR];
            const Box box( this_data_pol_dir.box() );
            
            for (BoxIterator bit(box); bit.ok(); ++bit) {
                IntVect iv = bit();
                if (iv[POLOIDAL_DIR]==lo_pol_LCSOL) {
                    int irad = iv[0] - lo_rad_LCSOL + nrad_core;
                    data_Z_loc[irad] = this_data_pol_dir(iv,0);
                }
            }
        }
    }
    
    MPI_Allreduce(data_Z_loc, data_Z, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    //Interpolating magnetic flux function onto the grid
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        a_interp[dit].copy(a_data[dit]);
        
        for (int dir(0); dir<SpaceDim; ++dir) {
            FArrayBox& this_interp_dir = a_interp[dit][dir];
            const FArrayBox& this_magFS_mapping_face = m_magFS_mapping_face[dit][dir];
            Box box( this_interp_dir.box() );
            
            for (BoxIterator bit(box); bit.ok(); ++bit) {
                IntVect iv = bit();
                
                int irad = this_magFS_mapping_face(iv, 1);
                double interp_coeff = this_magFS_mapping_face(iv, 2);
                
                if (irad == 0) {
                    this_interp_dir(iv,0) = data_Z[0]-(data_Z[1]-data_Z[0])*interp_coeff;
                }
                
                else if (irad == size) {
                    this_interp_dir(iv,0) = data_Z[size-1]+(data_Z[size-1]-data_Z[size-2])*interp_coeff;
                }
                
                else  {
                    this_interp_dir(iv,0) = data_Z[irad] - (data_Z[irad]-data_Z[irad-1]) * interp_coeff;
                }
            }
        }
    }
    
    
    delete [] data_Z;
    delete [] data_Z_loc;
}


void MagGeom::interpolateFromMagFS(const LevelData<FluxBox>& a_data_face,
                                   const LevelData<FArrayBox>& a_data_cell,
                                   LevelData<FArrayBox>& a_interp) const

{
    CH_assert(a_data_face.nComp() == 1);
    CH_assert(a_data_face.ghostVect() == IntVect::Zero);
    CH_assert(a_data_cell.nComp() == 1);
    CH_assert(a_data_cell.ghostVect() == IntVect::Zero);
    
    const MagBlockCoordSys& lcore_coord_sys = getBlockCoordSys(SingleNullBlockCoordSys::LCORE);
    int nrad_core  = lcore_coord_sys.domain().domainBox().size(RADIAL_DIR);
    int lo_pol_LCORE = lcore_coord_sys.domain().domainBox().smallEnd(POLOIDAL_DIR);
    
    const MagBlockCoordSys& lcsol_coord_sys = getBlockCoordSys(SingleNullBlockCoordSys::LCSOL);
    int nrad_sol = lcsol_coord_sys.domain().domainBox().size(RADIAL_DIR);
    int lo_pol_LCSOL = lcsol_coord_sys.domain().domainBox().smallEnd(POLOIDAL_DIR);
    int lo_rad_LCSOL = lcsol_coord_sys.domain().domainBox().smallEnd(RADIAL_DIR);
    
    //Finding cell-centers at the top of a tokamak
    int size = nrad_core + nrad_sol;
    double *data_Z = new double[size];
    double *data_Z_loc = new double[size];
    
    for (int i=0; i<size; ++i) {
        data_Z[i]=0.0;
        data_Z_loc[i]=0.0;
    }
    
    
    //Extracting the data along the LCORE and LCSOL boundaries (top cut),
    //and creating the corresponding magneitc flux function
    const DisjointBoxLayout& grids = a_data_cell.disjointBoxLayout();
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        int block_number = m_coord_sys->whichBlock(grids[dit]);
        if (block_number == SingleNullBlockCoordSys::LCORE) {
            const FArrayBox& this_data_pol_dir = a_data_face[dit][POLOIDAL_DIR];
            const Box box( this_data_pol_dir.box() );
            for (BoxIterator bit(box); bit.ok(); ++bit) {
                IntVect iv = bit();
                if (iv[POLOIDAL_DIR]==lo_pol_LCORE) {
                    int irad = iv[0];
                    data_Z_loc[irad] = this_data_pol_dir(iv,0);
                }
            }
        }
        
        if (block_number == SingleNullBlockCoordSys::LCSOL) {
            const FArrayBox& this_data_pol_dir = a_data_face[dit][POLOIDAL_DIR];
            const Box box( this_data_pol_dir.box() );
            
            for (BoxIterator bit(box); bit.ok(); ++bit) {
                IntVect iv = bit();
                if (iv[POLOIDAL_DIR]==lo_pol_LCSOL) {
                    int irad = iv[0] - lo_rad_LCSOL + nrad_core;
                    data_Z_loc[irad] = this_data_pol_dir(iv,0);
                }
            }
        }
    }
    
    MPI_Allreduce(data_Z_loc, data_Z, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    //Interpolating the magnetic flux function onto the grid
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        a_interp[dit].copy(a_data_cell[dit]);
        FArrayBox& this_interp = a_interp[dit];
        const FArrayBox& this_magFS_mapping_cell = m_magFS_mapping_cell[dit];
        Box box( this_interp.box() );
        
        for (BoxIterator bit(box); bit.ok(); ++bit) {
            IntVect iv = bit();
            
            int irad = this_magFS_mapping_cell(iv, 1);
            double interp_coeff = this_magFS_mapping_cell(iv, 2);
            
            if (irad == 0) {
                this_interp(iv,0) = data_Z[0]-(data_Z[1]-data_Z[0])*interp_coeff;
            }
            
            else if (irad == size) {
                this_interp(iv,0) = data_Z[size-1]+(data_Z[size-1]-data_Z[size-2])*interp_coeff;
            }
            
            else  {
                this_interp(iv,0) = data_Z[irad] - (data_Z[irad]-data_Z[irad-1]) * interp_coeff;
            }
        }
    }
    
    
    delete [] data_Z;
    delete [] data_Z_loc;
}


double
MagGeom::getMagFS( const RealVect& a_x ) const
{
    //Here, we assume that magnetic surfaces are labled by the z-coordinate at the tokamak top
    //(coinsides with the boundary between LCORE and RCORE) for the 8-block model geometry.
    //The algorithm should be straightforward to generalize, e.g., that is we can lable magnetic surfaces as
    //the z-coordinate along any radial block boundary. The key function to modify will be magFluxZ()
    //that needs to return the magnetic flux function magFlux() at the radial block boundary for a given z.
    
    
    //Solver parameters
    double tol = 1.0e-8;
    int max_iterations = 10000;
    
    //Getting min (lo_Z) and maximum (hi_Z) vertical coordinates that bound the solution.
    //Since we call this function to set the initial distribution (including some number ghost cells, presently 2)
    //we set lo_Z and hi_Z coresponding to the 3 ghost layer extensions (total we may 4 ghost layers).
    //(OLD COMMENT: Offset should not be too large, so that lo_Z and hi_Z can still be computed)
    const MagBlockCoordSys& lcore_coord_sys =  getBlockCoordSys(SingleNullBlockCoordSys::LCORE);
    RealVect dx_core = lcore_coord_sys.dx();
    double offset_core = 3.0*dx_core[0];
    RealVect mapped_coord_loZ;
    mapped_coord_loZ[1] = lcore_coord_sys.lowerMappedCoordinate(1);
    mapped_coord_loZ[0] = lcore_coord_sys.lowerMappedCoordinate(0) - offset_core;
    Real lo_Z = lcore_coord_sys.realCoord(mapped_coord_loZ)[1];
    
    const MagBlockCoordSys& lcsol_coord_sys =  getBlockCoordSys(SingleNullBlockCoordSys::LCSOL);
    RealVect dx_sol = lcsol_coord_sys.dx();
    double offset_sol = 3.0*dx_sol[0];
    RealVect mapped_coord_hiZ;
    mapped_coord_hiZ[1] = lcsol_coord_sys.lowerMappedCoordinate(1);
    mapped_coord_hiZ[0] = lcsol_coord_sys.upperMappedCoordinate(0) + offset_sol;
    Real hi_Z = lcsol_coord_sys.realCoord(mapped_coord_hiZ)[1];
    
    //Solving for the Z coordinate at the top of the core block, which designates the flux surface
    double tmp_lo_Z = lo_Z;
    double tmp_hi_Z = hi_Z;
    double f_lo_Z, f_hi_Z, f_half_Z, half_Z;
    
    bool residual_tolerance_satisfied = false;
    double magFluxValue;
    if (!m_model_flux_function) magFluxValue = lcore_coord_sys.getMagneticFlux(a_x);
    else magFluxValue = modelMagFlux(a_x);

    int iter = 0;
    while ( iter++ < max_iterations && !residual_tolerance_satisfied ) {
        
        half_Z = 0.5 * (tmp_lo_Z + tmp_hi_Z);
        f_lo_Z = magFluxZ(tmp_lo_Z) - magFluxValue;
        f_hi_Z = magFluxZ(tmp_hi_Z) - magFluxValue;
        f_half_Z = magFluxZ(half_Z) - magFluxValue;
        
        if (f_lo_Z * f_half_Z < 0) {
            tmp_hi_Z = half_Z;
        }
        
        else {
            tmp_lo_Z = half_Z;
        }
        
        if ( (magFluxZ(half_Z + tol) - magFluxValue) * (magFluxZ(half_Z - tol) - magFluxValue) < 0 ) {
            residual_tolerance_satisfied = true;
        }
        
    }
    return half_Z;
    
}

double MagGeom::magFluxZ(const double a_Z) const
{

   //Presently returns the magnetic flux function for the model geometry case
   //LATER: generalize by making use of the DCT representation of psi
   
   const SingleNullBlockCoordSys& lcore_coord_sys =  (SingleNullBlockCoordSys&)getBlockCoordSys(SingleNullBlockCoordSys::LCORE);
   
   RealVect coord;
   coord[0] = (lcore_coord_sys.getMagAxis())[0];
   coord[1] = a_Z;

   if (m_model_flux_function) coord[0] = 1.6;
 
   double result;

   if (!m_model_flux_function) {
     result =  lcore_coord_sys.getMagneticFlux(coord);
   }

   else {
     result =  modelMagFlux(coord);
   }

   return result;
}


double MagGeom::modelMagFlux(const RealVect& a_x) const
{
  //Presently returns the magnetic flux function for the model geometry case                                                                                                                      
  //LATER: generalize by making use of the DCT representation of psi                                                                                                                              
  double a=1.2;
  double b=0.9;
  double c=0.7;
  double R0=1.6;
  double Z0=1.07967;
  double Zxpt = 0.4;
  double val;

  if (a_x[1]>Zxpt) {
    val = cos(a*(a_x[0]-R0)) + b*sin((a_x[1]-Z0)) - c*(a_x[1]-Z0);
  }
  else {
    val = cos(a*(a_x[0]-R0)) + b*sin((2.0*Zxpt-a_x[1]-Z0)) - c*(2.0*Zxpt-a_x[1]-Z0);
  }

  return val;
}

void
MagGeom::unmapGradient( const LevelData<FArrayBox>& a_mapped_gradient,
                        LevelData<FArrayBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 3);
   CH_assert(a_gradient.nComp() == 3);

#if CFG_DIM ==3
   unmap3DGradient(a_mapped_gradient, a_gradient);
#else

   LevelData<FArrayBox> mapped_poloidal_gradient(m_gridsFull, 2, a_mapped_gradient.ghostVect());
   projectPoloidalVector(a_mapped_gradient, mapped_poloidal_gradient);

   LevelData<FArrayBox> poloidal_gradient(m_gridsFull, 2, a_gradient.ghostVect());
   unmapPoloidalGradient(mapped_poloidal_gradient, poloidal_gradient);

   injectPoloidalVector(poloidal_gradient, a_gradient);

#endif
}



void
MagGeom::unmapGradient( const LevelData<FluxBox>& a_mapped_gradient,
                        LevelData<FluxBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 3);
   CH_assert(a_gradient.nComp() == 3);
#if CFG_DIM ==3
   unmap3DGradient(a_mapped_gradient, a_gradient);
#else

   LevelData<FluxBox> mapped_poloidal_gradient(m_gridsFull, 2, a_mapped_gradient.ghostVect());
   projectPoloidalVector(a_mapped_gradient, mapped_poloidal_gradient);

   LevelData<FluxBox> poloidal_gradient(m_gridsFull, 2, a_gradient.ghostVect());
   unmapPoloidalGradient(mapped_poloidal_gradient, poloidal_gradient);

   injectPoloidalVector(poloidal_gradient, a_gradient);
#endif
}


void
MagGeom::unmap3DGradient( const LevelData<FArrayBox>& a_mapped_gradient,
                          LevelData<FArrayBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 3);
   CH_assert(a_gradient.nComp() == 3);

   // Multiply by NJInverse
   const DisjointBoxLayout& grids = a_mapped_gradient.disjointBoxLayout();
   const IntVect ghosts( a_mapped_gradient.ghostVect());
   setPointwiseNJInverseOnCells(grids, ghosts);
   
   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      FORT_MULT_NJINVERSE(CHF_BOX(a_mapped_gradient[dit].box()),
                          CHF_CONST_FRA(a_mapped_gradient[dit]),
                          CHF_CONST_FRA(m_NJinverse_cell_centered[dit]),
                          CHF_FRA(a_gradient[dit]));
   }
}


void
MagGeom::unmap3DGradient( const LevelData<FluxBox>& a_mapped_gradient,
                          LevelData<FluxBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 3);
   CH_assert(a_gradient.nComp() == 3);

   // Multiply by NJInverse
   const DisjointBoxLayout& grids = a_mapped_gradient.disjointBoxLayout();
   const IntVect ghosts( a_mapped_gradient.ghostVect());
   setPointwiseNJInverseOnFaces(grids, ghosts); 

   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         FORT_MULT_NJINVERSE(CHF_BOX(a_mapped_gradient[dit][dir].box()),
                             CHF_CONST_FRA(a_mapped_gradient[dit][dir]),
                             CHF_CONST_FRA(m_NJinverse_face_centered[dit][dir]),
                             CHF_FRA(a_gradient[dit][dir]));
      }
   }
}


void
MagGeom::unmapPoloidalGradient( const LevelData<FArrayBox>& a_mapped_gradient,
                                LevelData<FArrayBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 2);
   CH_assert(a_gradient.nComp() == 2);
   CH_assert(SpaceDim == 2);

   // Multiply by NJInverse
   const DisjointBoxLayout& grids = a_mapped_gradient.disjointBoxLayout();
   const IntVect ghosts( a_mapped_gradient.ghostVect());
   setPointwiseNJInverseOnCells(grids, ghosts);

   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      FORT_MULT_NJINVERSE(CHF_BOX(a_mapped_gradient[dit].box()),
                          CHF_CONST_FRA(a_mapped_gradient[dit]),
                          CHF_CONST_FRA(m_NJinverse_cell_centered[dit]),
                          CHF_FRA(a_gradient[dit]));
   }
}



void
MagGeom::unmapPoloidalGradient( const LevelData<FluxBox>& a_mapped_gradient,
                                LevelData<FluxBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 2);
   CH_assert(a_gradient.nComp() == 2);
   CH_assert(SpaceDim == 2);

   // Multiply by NJInverse
   const DisjointBoxLayout& grids = a_mapped_gradient.disjointBoxLayout();
   const IntVect ghosts( a_mapped_gradient.ghostVect());
   setPointwiseNJInverseOnFaces(grids, ghosts); 
   
   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         FORT_MULT_NJINVERSE(CHF_BOX(a_mapped_gradient[dit][dir].box()),
                             CHF_CONST_FRA(a_mapped_gradient[dit][dir]),
                             CHF_CONST_FRA(m_NJinverse_face_centered[dit][dir]),
                             CHF_FRA(a_gradient[dit][dir]));
      }
   }
}



void
MagGeom::injectPoloidalVector( const LevelData<FArrayBox>& a_poloidal_vector,
                               LevelData<FArrayBox>&       a_vector ) const
{
   CH_assert(a_poloidal_vector.nComp() == 2);
   CH_assert(a_vector.nComp() == 3);

   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      const FArrayBox& this_poloidal_vector = a_poloidal_vector[dit];
      FArrayBox& this_vector = a_vector[dit];

      this_vector.copy(this_poloidal_vector, 0, 0, 1);
      this_vector.setVal(0., 1);
      this_vector.copy(this_poloidal_vector, 1, 2, 1);
   }   
}



void
MagGeom::injectPoloidalVector( const LevelData<FluxBox>& a_poloidal_vector,
                               LevelData<FluxBox>&       a_vector ) const
{
   CH_assert(a_poloidal_vector.nComp() == 2);
   CH_assert(a_vector.nComp() == 3);

   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      const FluxBox& this_poloidal_vector = a_poloidal_vector[dit];
      FluxBox& this_vector = a_vector[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         const FArrayBox& this_poloidal_vector_dir = this_poloidal_vector[dir];
         FArrayBox& this_vector_dir = this_vector[dir];

         this_vector_dir.copy(this_poloidal_vector_dir, 0, 0, 1);
         this_vector_dir.setVal(0., 1);
         this_vector_dir.copy(this_poloidal_vector_dir, 1, 2, 1);
      }
   }   
}



void
MagGeom::projectPoloidalVector( const LevelData<FArrayBox>& a_vector,
                                LevelData<FArrayBox>&       a_poloidal_vector ) const
{
   CH_assert(a_vector.nComp() == 3);
   CH_assert(a_poloidal_vector.nComp() == 2);

   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      const FArrayBox& this_vector = a_vector[dit];
      FArrayBox& this_poloidal_vector = a_poloidal_vector[dit];

      this_poloidal_vector.copy(this_vector, 0, 0, 1);
      this_poloidal_vector.copy(this_vector, 2, 1, 1);
   }   
}



void
MagGeom::projectPoloidalVector( const LevelData<FluxBox>& a_vector,
                                LevelData<FluxBox>&       a_poloidal_vector ) const
{
   CH_assert(a_vector.nComp() == 3);
   CH_assert(a_poloidal_vector.nComp() == 2);

   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      const FluxBox& this_vector = a_vector[dit];
      FluxBox& this_poloidal_vector = a_poloidal_vector[dit];
      for (int dir=0; dir<2; ++dir) {
         const FArrayBox& this_vector_dir = this_vector[dir];
         FArrayBox& this_poloidal_vector_dir = this_poloidal_vector[dir];

         this_poloidal_vector_dir.copy(this_vector_dir, 0, 0, 1);
         this_poloidal_vector_dir.copy(this_vector_dir, 2, 1, 1);
      }
   }   
}

void
MagGeom::projectOntoParallel( LevelData<FArrayBox>& a_vector ) const
{
   CH_assert(a_vector.nComp() == CFG_DIM);

   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      
      FORT_PROJECT_ONTO_PARALLEL(CHF_BOX(a_vector[dit].box()),
                                 CHF_CONST_FRA(m_BFieldDir_cc[dit]),
                                 CHF_FRA(a_vector[dit]));

   }
}

void
MagGeom::projectOntoParallel( LevelData<FluxBox>& a_vector ) const
{
   CH_assert(a_vector.nComp() == CFG_DIM);
   
   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
      
         FORT_PROJECT_ONTO_PARALLEL(CHF_BOX(a_vector[dit][dir].box()),
                                    CHF_CONST_FRA(m_BFieldDir_fc[dit][dir]),
                                    CHF_FRA(a_vector[dit][dir]));
      }
   }
}

void
MagGeom::computeRadialProjection( LevelData<FArrayBox>& a_radComp,
                                  const LevelData<FArrayBox>& a_vector) const
{

   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
         
      FORT_COMPUTE_RADIAL_PROJECTION(CHF_BOX(a_radComp[dit].box()),
                                     CHF_CONST_FRA(m_BFieldDir_cc[dit]),
                                     CHF_CONST_FRA(a_vector[dit]),
                                     CHF_FRA1(a_radComp[dit],0));
      
   }
}

void
MagGeom::computeRadialProjection( LevelData<FluxBox>& a_radComp,
                                  const LevelData<FluxBox>& a_vector) const
{
   //Assumes that radComp box structure is within vector box structure
   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
            
         FORT_COMPUTE_RADIAL_PROJECTION(CHF_BOX(a_radComp[dit][dir].box()),
                                        CHF_CONST_FRA(m_BFieldDir_fc[dit][dir]),
                                        CHF_CONST_FRA(a_vector[dit][dir]),
                                        CHF_FRA1(a_radComp[dit][dir],0));
      }
   }
}

void
MagGeom::computePoloidalProjection( LevelData<FArrayBox>& a_polComp,
                                  const LevelData<FArrayBox>& a_vector) const
{

   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
         
      FORT_COMPUTE_POLOIDAL_PROJECTION(CHF_BOX(a_polComp[dit].box()),
                                       CHF_CONST_FRA(m_BFieldDir_cc[dit]),
                                       CHF_CONST_FRA(a_vector[dit]),
                                       CHF_FRA1(a_polComp[dit],0));
         
   }
}
   
void
MagGeom::computePoloidalProjection( LevelData<FluxBox>& a_polComp,
                                    const LevelData<FluxBox>& a_vector) const
{

   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
            
         FORT_COMPUTE_POLOIDAL_PROJECTION(CHF_BOX(a_polComp[dit][dir].box()),
                                          CHF_CONST_FRA(m_BFieldDir_fc[dit][dir]),
                                          CHF_CONST_FRA(a_vector[dit][dir]),
                                          CHF_FRA1(a_polComp[dit][dir],0));
      }
   }
}
   
void
MagGeom::computeMappedPoloidalGradientWithGhosts( const LevelData<FArrayBox>& a_phi,
                                                  LevelData<FArrayBox>&       a_field,
                                                  const int                   a_order ) const
{
   CH_TIME("MagGeom::computeMappedPoloidalGradientWithGhosts()");
   
   CH_assert(a_phi.ghostVect() >= 2*IntVect::Unit);
   CH_assert(a_field.nComp() == 2);
   //CH_assert(a_field.ghostVect() == IntVect::Unit);
   CH_assert(a_field.ghostVect() >= IntVect::Unit); // JRA ?
   CH_assert(a_order == 2 || a_order == 4);
      
   // The following assumes that we have potential values in at least two layers of ghost cells
   
   // Compute the field to second_order including one layer of ghost cells
      
   const DisjointBoxLayout& grids = a_phi.disjointBoxLayout();
      
   int tmp_order = 2;
      
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();
      Box box = grow(grids[dit],1);
         
      // Initializing to NaN to ensure that the field is actually computed everywhere it's used
      a_field[dit].setVal(1./0.);
         
      for (int dir=0; dir<2; dir++) {
         SpaceUtils::cellCenteredGradientComponent( box,
                                                    dir,
                                                    a_phi[dit],
                                                    dx,
                                                    tmp_order,
                                                    a_field[dit] );
      }
   }
      
   a_field.exchange();
      
   // If fourth-order, recompute the field at valid cell centers to fourth-order
      
   if ( a_order == 4 ) {
         
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
         RealVect dx = block_coord_sys.dx();
            
         for (int dir=0; dir<2; ++dir) {
            SpaceUtils::cellCenteredGradientComponent(  grids[dit],
                                                        dir,
                                                        a_phi[dit],
                                                        dx,
                                                        a_order,
                                                        a_field[dit] );
         }
      }
   }
}

void
MagGeom::computeMappedPoloidalGradientWithGhosts( const LevelData<FArrayBox>& a_phi,
                                                  LevelData<FluxBox>&         a_field,
                                                  const int                   a_order ) const
{
   CH_TIME("MagGeom::computeMappedPoloidalGradientWithGhosts()");
   
   CH_assert(a_phi.ghostVect() >= 2*IntVect::Unit);
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() >= IntVect::Unit);
   CH_assert(a_order == 2 || a_order == 4);
   
   // The following assumes that we have potential values in at least two layers of ghost cells at
   // all block boundaries.
   
   const DisjointBoxLayout& grids = a_phi.disjointBoxLayout();
   LevelData<FluxBox> phi_face(grids, 1, 2*IntVect::Unit);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<2; dir++) {
         IntVect grow_vect = 2*IntVect::Unit;
         grow_vect[dir] = 0;
         Box box = grow(grids[dit],grow_vect);
         
         SpaceUtils::faceInterpolate( dir,
                                      surroundingNodes(box,dir),
                                      a_order,
                                      a_phi[dit],
                                      phi_face[dit][dir] );
      }
   }
   phi_face.exchange();
   
   // Compute the field to second-order including one layer of transverse faces
   
   int tmp_order = 2;
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();
      Box box = grow(grids[dit],1);
      
      // Initializing to NaN to ensure that the field is actually computed everywhere it's used
      a_field[dit].setVal(1./0.);
      
      for (int dir=0; dir<2; dir++) {
         Box box_dir = surroundingNodes(box,dir);
         SpaceUtils::faceCenteredGradientComponent( box_dir,
                                                    dir,
                                                    a_phi[dit],
                                                    dx,
                                                    tmp_order,
                                                    a_field[dit][dir] );
         
         for (int tdir=0; tdir<2; ++tdir) {
            if (tdir != dir) {
               SpaceUtils::cellCenteredGradientComponent( box_dir,
                                                          tdir,
                                                          phi_face[dit][dir],
                                                          dx,
                                                          tmp_order,
                                                          a_field[dit][dir] );
            }
         }
      }
   }
   
   // Second-order extrapolate the field to the transverse physical boundary ghosts
   fillTransversePhysicalGhosts(a_field);
   
   // If fourth-order, recompute the field on valid cell faces
   
   if ( a_order == 4 ) {
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
         RealVect dx = block_coord_sys.dx();
         
         for (int dir=0; dir<2; dir++) {
            Box box_dir = surroundingNodes(grids[dit],dir);
            SpaceUtils::faceCenteredGradientComponent( box_dir,
                                                       dir,
                                                       a_phi[dit],
                                                       dx,
                                                       a_order,
                                                       a_field[dit][dir] );
            
            for (int tdir=0; tdir<2; ++tdir) {
               if (tdir != dir) {
                  SpaceUtils::cellCenteredGradientComponent(  box_dir,
                                                              tdir,
                                                              phi_face[dit][dir],
                                                              dx,
                                                              a_order,
                                                              a_field[dit][dir] );
               }
            }
         }
      }
      a_field.exchange();
   
   }
   
}

void
MagGeom::computeMappedGradient( const LevelData<FArrayBox>&  a_phi,
                                      LevelData<FluxBox>&    a_field,
                                const int                    a_order ) const
{
   CH_TIME("MagGeom::computeMappedGradient()");
   
   // Supposed to be general method for mapped gradients on cell faces for
   // 2D and 3D.
   // For now is copy of computeMappedPoloidalGradientWithGhosts, but without 
   // the ghosts
   //

   CH_assert(a_phi.ghostVect() >= 2*IntVect::Unit);
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() >= IntVect::Unit);
   CH_assert(a_order == 2);
   
   // The following assumes that we have potential values in at least two layers of ghost cells at
   // all block boundaries.
   
   const DisjointBoxLayout& grids = a_phi.disjointBoxLayout();
   LevelData<FluxBox> phi_face(grids, 1, 2*IntVect::Unit);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<2; dir++) {
         IntVect grow_vect = 2*IntVect::Unit;
         grow_vect[dir] = 0;
         Box box = grow(grids[dit],grow_vect);
         
         SpaceUtils::faceInterpolate( dir,
                                      surroundingNodes(box,dir),
                                      a_order,
                                      a_phi[dit],
                                      phi_face[dit][dir] );
      }
   }
   //phi_face.exchange(m_exchangeCopier2nd);
   
   // Compute the field to second-order including one layer of transverse faces
   
   int tmp_order = 2;
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();
      Box box = grow(grids[dit],1);
      
      // Initializing to NaN to ensure that the field is actually computed everywhere it's used
      a_field[dit].setVal(1./0.);
      
      for (int dir=0; dir<2; dir++) {
         Box box_dir = surroundingNodes(box,dir);
         SpaceUtils::faceCenteredGradientComponent( box_dir,
                                                    dir,
                                                    a_phi[dit],
                                                    dx,
                                                    tmp_order,
                                                    a_field[dit][dir] );
         
         for (int tdir=0; tdir<2; ++tdir) {
            if (tdir != dir) {
               SpaceUtils::cellCenteredGradientComponent( box_dir,
                                                          tdir,
                                                          phi_face[dit][dir],
                                                          dx,
                                                          tmp_order,
                                                          a_field[dit][dir] );
            }
         }
      }
   }

}

void
MagGeom::getEllipticOpRadCoeff(LevelData<FluxBox>& a_rad_coeff) const
{

   DataIterator dit = a_rad_coeff.dataIterator();
   
   if ( !m_rad_coeff.isDefined() || !(a_rad_coeff.ghostVect() <= m_rad_coeff.ghostVect()) ) {
      
      m_rad_coeff.define(a_rad_coeff);
      
      for (dit.begin(); dit.ok(); ++dit) {
         const FluxBox& this_bunit = m_BFieldDir_fc[dit];
         FluxBox& this_rad_coeff = m_rad_coeff[dit];
      
         for (int dir=0; dir<SpaceDim; ++dir) {
            const Box& box_dir = this_rad_coeff[dir].box();
         
            FORT_COMPUTE_RADIAL_ELLIPTIC_OP_COEFF(CHF_BOX(box_dir),
                                                  CHF_CONST_FRA(this_bunit[dir]),
                                                  CHF_FRA(this_rad_coeff[dir])
                                                  );
         
         }
      }
   }
   
   for (dit.begin(); dit.ok(); ++dit) {
      a_rad_coeff[dit].copy(m_rad_coeff[dit]);
   }
}

void
MagGeom::getEllipticOpRadCoeffMapped(LevelData<FluxBox>& a_rad_coeff) const
{
   CH_TIME("MagGeom::getEllipticOpRadCoefMapped()");
   
   DataIterator dit = a_rad_coeff.dataIterator();
   
   if ( !m_rad_coeff_mapped.isDefined() || !(a_rad_coeff.ghostVect() <= m_rad_coeff_mapped.ghostVect()) ) {
      
      m_rad_coeff_mapped.define(a_rad_coeff);
      
      const DisjointBoxLayout& grids = a_rad_coeff.disjointBoxLayout();
      
      const IntVect ghosts( a_rad_coeff.ghostVect());
      CH_assert(m_BFieldDir_fc.ghostVect()>=ghosts);
      
      LevelData<FluxBox> N(grids, SpaceDim*SpaceDim, ghosts);
      getPointwiseN(N);
      
      setPointwiseNJInverseOnFaces(grids, ghosts); 
      
      
      for (dit.begin(); dit.ok(); ++dit) {
         const FluxBox& this_bunit = m_BFieldDir_fc[dit];
         FluxBox& this_N = N[dit];
         FluxBox& this_NJinverse = m_NJinverse_face_centered[dit];
         FluxBox& this_rad_coeff = m_rad_coeff_mapped[dit];
         
         for (int dir=0; dir<SpaceDim; ++dir) {
            const Box& box_dir = this_rad_coeff[dir].box();
            
            FORT_COMPUTE_RADIAL_ELLIPTIC_OP_COEFF(CHF_BOX(box_dir),
                                                  CHF_CONST_FRA(this_bunit[dir]),
                                                  CHF_FRA(this_rad_coeff[dir])
                                                  );
            
            FORT_COMPUTE_ELLIPTIC_OP_COEFF_MAPPED(CHF_BOX(box_dir),
                                                  CHF_CONST_FRA(this_N[dir]),
                                                  CHF_CONST_FRA(this_NJinverse[dir]),
                                                  CHF_FRA(this_rad_coeff[dir])
                                                  );
            
         }
      }
   }
   
   for (dit.begin(); dit.ok(); ++dit) {
      a_rad_coeff[dit].copy(m_rad_coeff_mapped[dit]);
   }
}


void
MagGeom::computeEllipticOpCoefficients(LevelData<FluxBox>& a_perp_coeff,
                                       LevelData<FluxBox>& a_par_coeff ) const
{
   const DisjointBoxLayout& grids = a_perp_coeff.disjointBoxLayout();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const FluxBox& this_bunit = m_BFieldDir_fc[dit];
      FluxBox& this_perp_coeff = a_perp_coeff[dit];
      FluxBox& this_par_coeff = a_par_coeff[dit];
      
      for (int dir=0; dir<SpaceDim; ++dir) {
         const Box& box_dir = this_perp_coeff[dir].box();
         
         FORT_COMPUTE_ELLIPTIC_OP_COEFF(CHF_BOX(box_dir),
                                        CHF_CONST_FRA(this_bunit[dir]),
                                        CHF_FRA(this_perp_coeff[dir]),
                                        CHF_FRA(this_par_coeff[dir])
                                       );
         
      }
   }

}

void
MagGeom::computeEllipticOpCoefficientsMapped(LevelData<FluxBox>& a_perp_coeff,
                                             LevelData<FluxBox>& a_par_coeff ) const
{
   CH_TIME("MagGeom::computeEllipticOpCoefficientsMapped()");
   
   const DisjointBoxLayout& grids = a_perp_coeff.disjointBoxLayout();
  
   const IntVect ghosts( a_perp_coeff.ghostVect());
   CH_assert(m_BFieldDir_fc.ghostVect()>=ghosts);

   LevelData<FluxBox> N(grids, SpaceDim*SpaceDim, ghosts);
   getPointwiseN(N);
   
   setPointwiseNJInverseOnFaces(grids, ghosts);

   
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
      FluxBox& this_perp_coeff = a_perp_coeff[dit];
      FluxBox& this_par_coeff = a_par_coeff[dit];
      FluxBox& this_N = N[dit];
      FluxBox& this_NJinverse = m_NJinverse_face_centered[dit];

      if (!coord_sys.isSubGridGeom()) {
         
         const FluxBox& this_bunit = m_BFieldDir_fc[dit];

         for (int dir=0; dir<SpaceDim; ++dir) {
            const Box& box_dir = this_perp_coeff[dir].box();

            
            FORT_COMPUTE_ELLIPTIC_OP_COEFF(CHF_BOX(box_dir),
                                           CHF_CONST_FRA(this_bunit[dir]),
                                           CHF_FRA(this_perp_coeff[dir]),
                                           CHF_FRA(this_par_coeff[dir])
                                           );

            FORT_COMPUTE_ELLIPTIC_OP_COEFF_MAPPED(CHF_BOX(box_dir),
                                                  CHF_CONST_FRA(this_N[dir]),
                                                  CHF_CONST_FRA(this_NJinverse[dir]),
                                                  CHF_FRA(this_perp_coeff[dir])
                                                  );
         
            FORT_COMPUTE_ELLIPTIC_OP_COEFF_MAPPED(CHF_BOX(box_dir),
                                                  CHF_CONST_FRA(this_N[dir]),
                                                  CHF_CONST_FRA(this_NJinverse[dir]),
                                                  CHF_FRA(this_par_coeff[dir])
                                                  );

         }
      }

      //Computes mapped face-averaged coefficients by subgrid integrations
      else if (coord_sys.isSubGridGeom()) {
         computeSubGridEllipticOpCoefficientsMapped(this_perp_coeff, this_par_coeff, coord_sys);
      }

      //For the case of field-aligned coordinates all matrix entries except for
      //the diagonal entry corresponding to the parallel dir should be zero
      //This code forces it to be discrete zero to suppress truncation error
      //for the case of an approximate mapping
      if (coord_sys.isFieldAligned()) {

         int parallel_dir_entry;

         if (SpaceDim == 2) {
            //parallel_dir = poloidal direction
            parallel_dir_entry = 3;
         }
         else {
            //parallel_dir = toroidal direction
            parallel_dir_entry = 4;
         }
         
         for (int dir=0; dir<SpaceDim; ++dir) {
            for (int n=0; n<SpaceDim*SpaceDim; ++n) {
               if (n != parallel_dir_entry) {
                  this_par_coeff[dir].setVal(0.0, n);
               }
            }
         }
      }
   }
}
   
void
MagGeom::computeSubGridEllipticOpCoefficientsMapped(FluxBox& a_perp_coeff,
                                                    FluxBox& a_par_coeff,
                                                    const MagBlockCoordSys& a_coord_sys) const
{
   
   //get the subgrid parameters
   IntVect n_mr = a_coord_sys.getSubGridData();
   
   //check the validity of the subgrid parameters
   for (int dir=0; dir<SpaceDim; dir++) {
      if (n_mr[dir] < 2) {
         MayDay::Error("MagGeom::computeSubGridEllipticOpCoefficientsMapped refinement must be not less than 2 in all dirs");
      }
      
      if (n_mr[dir] % 2 != 0) {
         MayDay::Error("MagGeom::computeSubGridEllipticOpCoefficientsMapped refinement must be even in all dirs");
      }
   }

   
   RealVect dx = a_coord_sys.dx();
   
   for (int dir=0; dir<SpaceDim; ++dir) {
      const Box& box = a_perp_coeff[dir].box();
      
      FArrayBox Xi(box,SpaceDim);
      a_coord_sys.getFaceCenteredMappedCoords(dir, Xi);
      
      BoxIterator bit(box);
      for (bit.begin();bit.ok();++bit) {
         const IntVect& iv = bit();
         RealVect Xi_loc;
         
         Vector<Real> averaged_perp_coeff(SpaceDim*SpaceDim,0);
         Vector<Real> averaged_par_coeff(SpaceDim*SpaceDim,0);
         
         
         Vector<int> tanDir(SpaceDim-1,0);
         for (int tdir=0; tdir<SpaceDim-1; ++tdir) {
            tanDir[tdir] = dir + tdir + 1;
            if (tanDir[tdir] >= SpaceDim) tanDir[tdir] -= SpaceDim;
         }
         
         Xi_loc[dir] = Xi(iv,dir);
         
#if CFG_DIM == 2
         
         for (int i = -n_mr[tanDir[0]]/2; i < n_mr[tanDir[0]]/2; i++) {
            Xi_loc[tanDir[0]] = Xi(iv,tanDir[0]) + (i + 0.5) * dx[tanDir[0]]/n_mr[tanDir[0]];
#else
         for (int i = -n_mr[tanDir[0]]/2; i < n_mr[tanDir[0]]/2; i++) {
            for (int j = -n_mr[tanDir[1]]/2; j < n_mr[tanDir[1]]/2; j++) {
               Xi_loc[tanDir[0]] = Xi(iv,tanDir[0]) + (i + 0.5) * dx[tanDir[0]]/n_mr[tanDir[0]];
               Xi_loc[tanDir[1]] = Xi(iv,tanDir[1]) + (j + 0.5) * dx[tanDir[1]]/n_mr[tanDir[1]];
#endif
                  
               RealVect X_loc(a_coord_sys.realCoord(Xi_loc));
                  
               Vector<Real> NJinv = a_coord_sys.getPointwiseNJInverse(Xi_loc);
                  
               array<double,3> Bfield = a_coord_sys.computeBField(X_loc);
               array<double,3> b = computeUnitVector(Bfield);
                  
               Vector<Real> NT = a_coord_sys.getPointwiseNT(Xi_loc);
                  
               Vector<Real> perp_fac(SpaceDim*SpaceDim,0);
               if (SpaceDim == 2) {
                  perp_fac[0] = (1.0 - b[0] * b[0]);
                  perp_fac[1] = (    - b[0] * b[2]);
                  perp_fac[2] = (    - b[2] * b[0]);
                  perp_fac[3] = (1.0 - b[2] * b[2]);
               }
                  
               else {
                  perp_fac[0] = (1.0 - b[0] * b[0]);
                  perp_fac[1] = (    - b[0] * b[1]);
                  perp_fac[2] = (    - b[0] * b[2]);
                  perp_fac[3] = (    - b[1] * b[0]);
                  perp_fac[4] = (1.0 - b[1] * b[1]);
                  perp_fac[5] = (    - b[1] * b[2]);
                  perp_fac[6] = (    - b[2] * b[0]);
                  perp_fac[7] = (    - b[2] * b[1]);
                  perp_fac[8] = (1.0 - b[2] * b[2]);
               }
                  
               Vector<Real> par_fac(SpaceDim*SpaceDim,0);
               if (SpaceDim == 2 ) {
                  par_fac[0] = ( b[0] * b[0]);
                  par_fac[1] = ( b[0] * b[2]);
                  par_fac[2] = ( b[2] * b[0]);
                  par_fac[3] = ( b[2] * b[2]);
               }
               else {
                  par_fac[0] = ( b[0] * b[0]);
                  par_fac[1] = ( b[0] * b[1]);
                  par_fac[2] = ( b[0] * b[2]);
                  par_fac[3] = ( b[1] * b[0]);
                  par_fac[4] = ( b[1] * b[1]);
                  par_fac[5] = ( b[1] * b[2]);
                  par_fac[6] = ( b[2] * b[0]);
                  par_fac[7] = ( b[2] * b[1]);
                  par_fac[8] = ( b[2] * b[2]);
               }
                  
               //Get local mapped elliptic coefficients
               Vector<Real> loc_perp_coeff_tmp = matrixMultiplication(perp_fac, NJinv);
               Vector<Real> loc_perp_coeff = matrixMultiplication(NT, loc_perp_coeff_tmp);
                  
               Vector<Real> loc_par_coeff_tmp = matrixMultiplication(par_fac, NJinv);
               Vector<Real> loc_par_coeff = matrixMultiplication(NT, loc_par_coeff_tmp);
                  
               //Get the normalization factor
               double norm_fac = 1.0;
               for (int n=0; n<SpaceDim-1; ++n) {
                  norm_fac /= n_mr[tanDir[n]];
               }
                  
               //Average over all local values
               for (int nComp=0; nComp<SpaceDim*SpaceDim; ++nComp) {
                  averaged_perp_coeff[nComp] = averaged_perp_coeff[nComp] + norm_fac * loc_perp_coeff[nComp];
                  averaged_par_coeff[nComp]  = averaged_par_coeff[nComp]  + norm_fac * loc_par_coeff[nComp];
               }
#if CFG_DIM == 3
            }
#endif
         }
            
         for (int nComp=0; nComp<SpaceDim*SpaceDim; ++nComp) {
            a_perp_coeff[dir](iv,nComp) = averaged_perp_coeff[nComp];
            a_par_coeff[dir](iv,nComp)  = averaged_par_coeff[nComp];
         }
      }
   }
}
   

#if CFG_DIM == 3
void
MagGeom::interpolateFromShearedGhosts(LevelData<FArrayBox>& a_data) const
{

   int nComp = a_data.nComp();

   int ghost_mb_dir = (a_data.ghostVect()[TOROIDAL_DIR] < m_sheared_ghosts) ?
                       a_data.ghostVect()[TOROIDAL_DIR] : m_sheared_ghosts;

   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   
   for (SideIterator sit; sit.ok(); ++sit) {
      Side::LoHiSide side = sit();
      
      BoxLayoutData<FArrayBox> ghosts;
      
      Copier copier;
      const ProblemDomain& domain = grids.physDomain();
      
      IntVect included_ghosts = IntVect::Zero;
      if (m_extrablock_exchange) {
         included_ghosts[POLOIDAL_DIR] = a_data.ghostVect()[POLOIDAL_DIR];
      }
      
      //store ghosts
      if (side == Side::LoHiSide::Lo) {
         ghosts.define(m_shearedGhostBLLoEnd, nComp);

         copier.define(m_shearedGhostBLLoEnd,
                       grids,
                       domain,
                       included_ghosts,
                       false,
                       IntVect::Zero);
      }
      else {
         ghosts.define(m_shearedGhostBLHiEnd, nComp);

         copier.define(m_shearedGhostBLHiEnd,
                       grids,
                       domain,
                       included_ghosts,
                       false,
                       IntVect::Zero);
      }
      copier.reverse();

      // The reverse copier perform copyTo not only from valid data (as typycal copiers do)
      // but also from ghost data of a_data designated by included_ghosts GhostVect object
      a_data.copyTo(ghosts, copier);

      //iterate over boxes, fill the ghosts at the block boundaries
      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
	
         const MagBlockCoordSys& block_coords = getBlockCoordSys(grids[dit]);
         const ProblemDomain& domain = block_coords.domain();
         const Box& domain_box = domain.domainBox();
         
         const Box& base_box = grids[dit];
         
         if (    (side == Side::LoHiSide::Lo && base_box.smallEnd(m_mb_dir) == domain_box.smallEnd(m_mb_dir))
              || (side == Side::LoHiSide::Hi && base_box.bigEnd(m_mb_dir) == domain_box.bigEnd(m_mb_dir)) ) {
   
            Box bndryBox = adjCellBox(base_box, m_mb_dir, side, ghost_mb_dir);
            
            for (BoxIterator bit(bndryBox);bit.ok();++bit) {
               const IntVect& iv = bit();

               //get global index of the central ghost cell
               IntVect iv_ghost;
               for (int dir=0; dir<SpaceDim; ++dir) {
                  iv_ghost[dir] = m_sheared_remapped_index[dit](iv,dir);
               }

               //interpolate ghost data
               for (int comp=0; comp<nComp; ++comp) {
                  
                  double ghost_val = 0.0;
                  
                  // Get the coefficient that also designates whether the ghost cell is
                  //(a) fully emerged inside the saw-tooth BC (fac>1); physical BC handles that.
                  //(b) partially emerged in the saw-tooth (0<fac<1)
                  //(c) belongs to the poloidal interior (fac < 0)
                  
                  Real fac = m_sheared_interp_stencil[dit](iv,m_sheared_interp_order + 1);
                  
                  if (fac < 1.0) {
                     for (int n=0; n<m_sheared_interp_order + 1; ++n) {
                        IntVect iv_offset(iv_ghost);
                        iv_offset[POLOIDAL_DIR] += (int)(m_sheared_interp_stencil_offsets)[dit](iv,n);
                        ghost_val += ghosts[dit](iv_offset,comp) * m_sheared_interp_stencil[dit](iv,n);
                     }
                  
                     if (fac < 0.) {
                        a_data[dit](iv,comp) = ghost_val;
                     }
                     else {
                        a_data[dit](iv,comp) = fac*a_data[dit](iv,comp) + (1.0-fac) * ghost_val;
                     }
                  }
               }
            }
         }
      }
   }
}

void
MagGeom::initializeShearedMBGeom(const DisjointBoxLayout& a_grids)
{
   //Here, set number of ghost layers, which is later used by phaseGeom
   //to be consistent with advection scheme (uw5 needs 3 layers, uw3 needs 2)
   m_sheared_ghosts = 3;
   m_sheared_interp_order = 2;
   //   IntVect ghostVect(0,m_sheared_ghosts,0);
   IntVect ghostVect(1,m_sheared_ghosts,1);
   m_sheared_remapped_index.define(a_grids, SpaceDim, ghostVect);
   
   m_sheared_interp_stencil_offsets.define(a_grids, m_sheared_interp_order + 1, ghostVect);

   // Interpolation stenticil has one extra component to store info for saw-tooth boundary
   m_sheared_interp_stencil.define(a_grids, m_sheared_interp_order + 2, ghostVect);

   for (SideIterator sit; sit.ok(); ++sit) {
      Side::LoHiSide side = sit();
      getMagShearInterpCoeff(m_sheared_remapped_index,
                             m_sheared_interp_stencil,
                             m_sheared_interp_stencil_offsets,
                             a_grids,
                             side);
   }
   
   // The preceding has filled codim 1 extrablock ghosts on toroidal block boundaries.
   // We now fill the codim >=2 toroidal extrablock ghost cells which are needed by
   // ToroidalBlockLevelExchange in the MBSolvers.

   IntVect toroidal_dir_only = IntVect::Zero;
   toroidal_dir_only[TOROIDAL_DIR] = 1;

   for (int codim=2; codim<=SpaceDim; ++codim) {
      exchangeExtraBlockGhosts(m_sheared_remapped_index, codim, toroidal_dir_only);
      exchangeExtraBlockGhosts(m_sheared_interp_stencil, codim, toroidal_dir_only);
      exchangeExtraBlockGhosts(m_sheared_interp_stencil_offsets, codim, toroidal_dir_only);
   }

   getShearedGhostBoxLayout();
}

void
MagGeom::getShearedGhostBoxLayout()
{
   const DisjointBoxLayout& grids = m_sheared_remapped_index.disjointBoxLayout();
   
   m_shearedGhostBLHiEnd.deepCopy(grids);
   m_shearedGhostBLLoEnd.deepCopy(grids);
   
   for (int ivec = 0; ivec < grids.rawPtr()->size(); ++ivec)
   {
    
      const Box& base_box = *const_cast<Box*>(&((*grids.rawPtr())[ivec].box));
      
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(base_box);
      const ProblemDomain& domain = block_coord_sys.domain();
      const Box& domain_box = domain.domainBox();
      
      for (SideIterator sit; sit.ok(); ++sit) {
         Side::LoHiSide side = sit();
   
         if (    (side == Side::LoHiSide::Lo && base_box.smallEnd(m_mb_dir) == domain_box.smallEnd(m_mb_dir))
              || (side == Side::LoHiSide::Hi && base_box.bigEnd(m_mb_dir) == domain_box.bigEnd(m_mb_dir)) ) {
            
            Box bndryBox = adjCellBox(base_box, m_mb_dir, side, m_sheared_ghosts);
            
            double *iv_lo = new double[SpaceDim];
            double *iv_hi = new double[SpaceDim];
            double *iv_lo_loc = new double[SpaceDim];
            double *iv_hi_loc = new double[SpaceDim];
            
            for (int dir=0; dir<SpaceDim; ++dir) {
               iv_lo[dir] = 0;
               iv_hi[dir] = 0;
               iv_lo_loc[dir] = 0;
               iv_hi_loc[dir] = 0;
            }
            
            for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
               if (grids[dit].contains(base_box)) {
                  for (int dir=0; dir<SpaceDim; ++dir) {
                     iv_lo_loc[dir] = int(m_sheared_remapped_index[dit].min(bndryBox,dir));
                     iv_hi_loc[dir] = int(m_sheared_remapped_index[dit].max(bndryBox,dir));
                  }
               }
            }

            MPI_Allreduce(iv_lo_loc, iv_lo, SpaceDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(iv_hi_loc, iv_hi, SpaceDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            
            IntVect hiEnd_remapped;
            IntVect loEnd_remapped;
            
            for (int dir=0; dir<SpaceDim; ++dir) {
               loEnd_remapped[dir] = iv_lo[dir];
               hiEnd_remapped[dir] = iv_hi[dir];
            }
            
            // For the case of core-geometry the offset is only m_sheared_interp_order/2
            // However, the saw-tooth BCs requires larger offset near poloidal boundaries
            // If affects the performance, treat core and SN geom separately
            loEnd_remapped[POLOIDAL_DIR] -= (m_sheared_interp_order + 1);
            hiEnd_remapped[POLOIDAL_DIR] += (m_sheared_interp_order + 1);
            
            Box remapped_box(loEnd_remapped, hiEnd_remapped);
            
            if (side == Side::LoHiSide::Lo) {
               Box& new_box = *const_cast<Box*>(&((*m_shearedGhostBLLoEnd.rawPtr())[ivec].box));
               new_box = remapped_box;
            }
            else {
               Box& new_box = *const_cast<Box*>(&((*m_shearedGhostBLHiEnd.rawPtr())[ivec].box));
               new_box = remapped_box;
            }
            
            delete [] iv_lo;
            delete [] iv_hi;
            delete [] iv_lo_loc;
            delete [] iv_hi_loc;

         }
      }
   }

   m_shearedGhostBLLoEnd.closeNoSort();
   m_shearedGhostBLHiEnd.closeNoSort();
}

void
MagGeom::getMagShearInterpCoeff(LevelData<FArrayBox>& a_remapped_iv,
                                LevelData<FArrayBox>& a_interp_stencil,
                                LevelData<FArrayBox>& a_interp_stencil_offsets,
                                const DisjointBoxLayout& a_grids,
                                const Side::LoHiSide& a_side) const
{
   
   const DisjointBoxLayout& grids = a_remapped_iv.disjointBoxLayout();
  
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      
      const MagBlockCoordSys& coord_sys_src = getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = coord_sys_src.domain();
      const Box& domain_box = domain.domainBox();
      
      int nSrc = m_coord_sys->whichBlock(grids[dit]);

      const Box& base_box = grids[dit];
      
      if (    (a_side == Side::LoHiSide::Lo && base_box.smallEnd(m_mb_dir) == domain_box.smallEnd(m_mb_dir))
           || (a_side == Side::LoHiSide::Hi && base_box.bigEnd(m_mb_dir) == domain_box.bigEnd(m_mb_dir)) ) {
         
         Box bndryBox = adjCellBox(base_box, m_mb_dir, a_side, m_sheared_ghosts);

         bndryBox.grow(POLOIDAL_DIR,1);
         bndryBox.grow(RADIAL_DIR,1);
         
         FArrayBox mappedBndryCoord(bndryBox,SpaceDim);
         coord_sys_src.getCellCenteredMappedCoords(mappedBndryCoord);
         
         RealVect dx = coord_sys_src.dx();

         BoxIterator bit(bndryBox);
         for (bit.begin();bit.ok();++bit) {
            const IntVect& iv = bit();

            //Get mapped cell-centered coord in src block
            RealVect xiSrc;
            for (int dir=0; dir<SpaceDim; ++dir) {
               xiSrc[dir] = mappedBndryCoord(iv,dir);
            }
            
            // Get the global index of the src cell center
            IntVect iv_src;
            for (int dir=0; dir<SpaceDim; ++dir) {
               iv_src[dir] = floor(xiSrc[dir]/dx[dir]);
            }

            // dst block number
            int nDst;

            //Index of the xiSrc in dst block
            IntVect ivDst;
         
            //Second-order interpolation stencil
            Vector<Real> interpStencil(m_sheared_interp_order + 2, -1.0);
            Vector<int> interpStencilOffsets(m_sheared_interp_order + 1, 0);
         
            //Compute ivDst and interpStencil
            //NB: the result iv and interp_stencil should be the same for all ghost cells
            //so, in principle we do not need to recompute it for all toroidal layers
            //however, this is cheap calculation, happens just once, so keep it anyways
            m_coord_sys->toroidalBlockRemapping(ivDst, nDst, interpStencil, interpStencilOffsets, xiSrc, nSrc, a_side);
         
            for (int dir=0; dir<SpaceDim; ++dir) {
               a_remapped_iv[dit](iv,dir) = ivDst[dir];
            }
         
            for (int n=0; n<interpStencil.size(); ++n) {
               a_interp_stencil[dit](iv,n) = interpStencil[n];
            }
            
            for (int n=0; n<interpStencilOffsets.size(); ++n) {
               a_interp_stencil_offsets[dit](iv,n) = interpStencilOffsets[n];
            }

         }
      }
   }
}

Vector<Box>
MagGeom::getToroidalBoundaryBoxes( const int              a_block_number,
                                   const Side::LoHiSide&  a_side ) const
{
   const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(a_block_number);
   const ProblemDomain& domain = block_coord_sys.domain();
   const Box& domain_box = domain.domainBox();
      
   Box bndryBox = adjCellBox(domain_box, m_mb_dir, a_side, 1);

   FArrayBox mappedBndryCoord(bndryBox,SpaceDim);
   block_coord_sys.getCellCenteredMappedCoords(mappedBndryCoord);
         
   RealVect dx = block_coord_sys.dx();

   IntVectSet bndry_cells;

   for ( BoxIterator bit(bndryBox); bit.ok(); ++bit ) {
      const IntVect& iv = bit();

      RealVect xiSrc;
      for (int dir=0; dir<SpaceDim; ++dir) {
         xiSrc[dir] = mappedBndryCoord(iv,dir);
      }
            
      IntVect iv_src;
      for (int dir=0; dir<SpaceDim; ++dir) {
         iv_src[dir] = floor(xiSrc[dir]/dx[dir]);
      }

      int nDst;
      IntVect ivDst;
      Vector<Real> interpStencil(m_sheared_interp_order + 2, 0);
      Vector<int> interpStencilOffsets(m_sheared_interp_order + 1, 0);

      m_coord_sys->toroidalBlockRemapping(ivDst, nDst, interpStencil, interpStencilOffsets, xiSrc, a_block_number, a_side);

      if ( interpStencil[m_sheared_interp_order + 1] >= 0.5 ) {
         bndry_cells |= iv;
      }
   }

   // Convert the IntVectSet into a vector of boxes
   Vector<Box> box_vector = bndry_cells.boxes();

   // Coalesce the IntVectSet into boxes

   Vector<Box> coalesced_box_vector;

   if ( box_vector.size() > 0 ) {
      Box accum_box = box_vector[0];

      for (int n=1; n<box_vector.size(); ++n) {
         Box new_box = accum_box;
         new_box.minBox(box_vector[n]);

         IntVectSet box_ivs(new_box);

         // Check that the new box is full, i.e., bndry_cells
         // contains all of the new box's IntVects
         if ( bndry_cells.contains(box_ivs) ) {
            // Keep accumulating and pushback if no more boxes
            accum_box = new_box;
            if ( n == box_vector.size()-1 ) {
               coalesced_box_vector.push_back(bdryBox(accum_box,m_mb_dir,a_side,1));
            }
         }
         else {
            // Can't accumulate any more.  Push the current accumulation box
            // and start a new one
            coalesced_box_vector.push_back(bdryBox(accum_box,m_mb_dir,a_side,1));
            accum_box = box_vector[n];
         }
      }
   }

   return coalesced_box_vector;
}

#endif

Vector<Box>
MagGeom::getBoundaryBoxes( const int              a_block_number,
                           const int              a_dir,
                           const Side::LoHiSide&  a_side ) const
{
   const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(a_block_number);
   const ProblemDomain& domain = block_coord_sys.domain();
   const Box& domain_box = domain.domainBox();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = m_coord_sys->boundaries();

   Vector<Box> boxes;

#if CFG_DIM==3
   if ( m_sheared_mb_geom && a_dir == m_mb_dir ) {
      Vector<Box> toroidal_boxes = getToroidalBoundaryBoxes(a_block_number, a_side);

      for (int n=0; n<toroidal_boxes.size(); ++n) {
         boxes.push_back(toroidal_boxes[n]);
      }
   }
   else {
      if ( boundaries[a_block_number][a_dir + a_side*SpaceDim].isDomainBoundary() ) {
         boxes.push_back(bdryBox(domain_box, a_dir, a_side, 1));
      }
   }
#else
   if ( boundaries[a_block_number][a_dir + a_side*SpaceDim].isDomainBoundary() ) {
      boxes.push_back(bdryBox(domain_box, a_dir, a_side, 1));
   }
#endif

   return boxes;
}


#include "NamespaceFooter.H"
