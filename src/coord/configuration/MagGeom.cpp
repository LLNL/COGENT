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

#include "FieldSolverF_F.H"
#include "NamespaceHeader.H"

#undef TEST_DIVERGENCE_CLEANING


MagGeom::MagGeom(ParmParse&                 a_pp,
                 const MagCoordSys*         a_coordSysPtr,
                 const DisjointBoxLayout&   a_grids,
                 int                        a_ghosts)
   : MultiBlockLevelGeom(a_coordSysPtr, a_grids, a_ghosts),
     m_coord_sys(a_coordSysPtr),
     m_correct_field(false),
     m_pp_prefix(a_pp.prefix())
{
   if (a_pp.contains("correct_field")) {
     a_pp.get("correct_field", m_correct_field);
   }

   if (a_pp.contains("extrablock_exchange")) {
      a_pp.get("extrablock_exchange", m_extrablock_exchange);
   }
   else {
      m_extrablock_exchange = false;
   }

   if (a_pp.contains("model_geometry")) {
      a_pp.get("model_geometry", m_model_geometry);
   }
   else {
      m_model_geometry = false;
   }

   // Check for full periodicity
   m_fully_periodic = true;
   for (int block=0; block<a_coordSysPtr->numBlocks(); ++block) {
     const ProblemDomain& domain = getBlockCoordSys(block).domain();
     for (int dir=0; dir<SpaceDim; dir++) {
       m_fully_periodic = m_fully_periodic && domain.isPeriodic(dir);
     }
   }

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

   if (a_pp.contains("plot_field_alignment")) {
      bool plot_field_alignment;
      a_pp.get("plot_field_alignment", plot_field_alignment);

      if (plot_field_alignment) {
         // Plot the angle (in degrees) between the poloidal field component and the
         // mapped coordinate system vector in the poloidal direction.
         double time = 0.;
         plotFieldAlignment(time);
      }
   }

   if (a_pp.contains("plot_field_divergence")) {
      bool plot_field_divergence;
      a_pp.get("plot_field_divergence", plot_field_divergence);

      if (plot_field_divergence) {
         double time = 0.;
         plotFieldDivergence(m_BField_fc, time);
      }
   }

   if ( m_coord_sys->numBlocks() > 1 ) {
      m_mblexPtr = new MultiBlockLevelExchangeAverage();

      int spaceOrder(4);
      m_mblexPtr->define( this, m_ghosts, spaceOrder );

      RefCountedPtr<MultiBlockCoordSys> coordSysRCP((MultiBlockCoordSys*)m_coordSysPtr);
      coordSysRCP.neverDelete();  // Workaround for some problem with RefCountedPtr

      m_exchange_transverse_block_register = new BlockRegister(coordSysRCP, a_grids, 1);
   }
   else {
      m_mblexPtr = NULL;
   }
    
    //Initializations for dealignment-correction utilities (for SN geometry only)
    if ( typeid(*m_coord_sys) == typeid(SingleNullCoordSys) ) {
        m_magFS_mapping_cell.define(a_grids, 3, IntVect::Zero);
        m_magFS_mapping_face.define(a_grids, 3, IntVect::Zero);
        computeMagFluxMappingCell( m_magFS_mapping_cell);
        computeMagFluxMappingFace( m_magFS_mapping_face);
    }

}


MagGeom::~MagGeom()
{
   if (m_exchange_transverse_block_register) delete m_exchange_transverse_block_register;
   if (m_mblexPtr) delete m_mblexPtr;
}



void
MagGeom::getMetricTerms( LevelData<FluxBox>& a_N,
                         LevelData<FluxBox>& a_tanGradN ) const
{
   CH_assert(a_N.getBoxes()==a_tanGradN.getBoxes());

   const DisjointBoxLayout& grids = a_N.disjointBoxLayout();

   for (DataIterator dit(a_N.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
      coord_sys.getN(a_N[dit], a_N[dit].box());
   }

   for (DataIterator dit(a_tanGradN.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
      coord_sys.computeTangentialGrad(a_tanGradN[dit], a_N[dit], a_tanGradN[dit].box());
   }

   m_coord_sys->postProcessMetricData(a_N);
   m_coord_sys->postProcessMetricData(a_tanGradN);
}



void
MagGeom::getCellVolumes( LevelData<FArrayBox>& a_volume ) const
{
   DataIterator dit = a_volume.dataIterator();

   if ( !m_cell_volume.isDefined() || !(a_volume.ghostVect() <= m_cell_volume.ghostVect()) ) {

      m_cell_volume.define(a_volume);

      const DisjointBoxLayout& grids = m_cell_volume.disjointBoxLayout();

      for (dit.begin(); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);

         Box grown_box(grow(m_cell_volume[dit].box(),1));
         FluxBox N(grown_box, coord_sys.getNumN());
         coord_sys.getN(N, grown_box);
         coord_sys.cellVol(m_cell_volume[dit], N, m_cell_volume[dit].box());
      }

      m_coord_sys->postProcessMetricData(m_cell_volume);

      m_cell_volume.exchange();
   }


   for (dit.begin(); dit.ok(); ++dit) {
      a_volume[dit].copy(m_cell_volume[dit]);
   }
}



void
MagGeom::getJ( LevelData<FArrayBox>& a_J ) const
{
   if ( !m_J.isDefined() || !(a_J.ghostVect() <= m_J.ghostVect()) ) {

      m_J.define(a_J);

      const DisjointBoxLayout& grids = a_J.disjointBoxLayout();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
         coord_sys.getAvgJ(m_J[dit], m_J[dit].box());
      }

      m_coord_sys->postProcessMetricData(m_J);

      m_J.exchange();
   }

   for (DataIterator dit(a_J.dataIterator()); dit.ok(); ++dit) {
      a_J[dit].copy(m_J[dit]);
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

      m_N_face_centered.exchange();
   }

   for (dit.begin(); dit.ok(); ++dit) {
      a_N[dit].copy(m_N_face_centered[dit]);
   }
}



void
MagGeom::getPointwiseNJinverse( LevelData<FluxBox>& a_NJinverse ) const
{

   DataIterator dit = a_NJinverse.dataIterator();

   if ( !m_NJinverse_face_centered.isDefined() || !(a_NJinverse.ghostVect() <= m_NJinverse_face_centered.ghostVect()) ) {

      m_NJinverse_face_centered.define(a_NJinverse);

      const DisjointBoxLayout& grids = m_NJinverse_face_centered.disjointBoxLayout();

      for (dit.begin(); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
         coord_sys.getPointwiseNJInverse(m_NJinverse_face_centered[dit]);
      }

      m_coord_sys->postProcessMetricData(m_NJinverse_face_centered);

      m_NJinverse_face_centered.exchange();
   }

   for (dit.begin(); dit.ok(); ++dit) {
      a_NJinverse[dit].copy(m_NJinverse_face_centered[dit]);
   }
}


void
MagGeom::divideJonValid( LevelData<FArrayBox>& a_data ) const
{
   // Multiplies the argument by J on valid cells

   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   IntVect grown_vect = IntVect::Unit;
   LevelData<FArrayBox> J(grids, 1, grown_vect);
   getJ(J);

   LevelData<FArrayBox> grown_data(grids, 1, grown_vect);

   for (DataIterator dit(grown_data.dataIterator()); dit.ok(); ++dit) {
      grown_data[dit].copy(a_data[dit]);
   }
   grown_data.exchange();

   // Compute the fourth-order quotient
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
      cellFGToCellF(a_data[dit], grown_data[dit], J[dit], grids[dit], coord_sys.domain(), true);
   }
}


void
MagGeom::multJonValid( LevelData<FArrayBox>& a_data ) const
{
   // Multiplies the argument by J on valid cells

   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   IntVect grown_vect = IntVect::Unit;
   LevelData<FArrayBox> J(grids, 1, grown_vect);
   getJ(J);

   LevelData<FArrayBox> grown_data(grids, 1, grown_vect);

   for (DataIterator dit(grown_data.dataIterator()); dit.ok(); ++dit) {
      grown_data[dit].copy(a_data[dit]);
   }
   grown_data.exchange();

   fourthOrderCellProd(a_data, J, grown_data);

   // Compute the fourth-order product
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
      fourthOrderCellProd(a_data[dit], grown_data[dit], J[dit], grids[dit], coord_sys.domain(), true);
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
   const DisjointBoxLayout& grids = a_product.disjointBoxLayout();

   LevelData<FluxBox> N(grids, SpaceDim*SpaceDim, IntVect::Unit); 
   LevelData<FluxBox> tanGradN(grids, SpaceDim*SpaceDim, IntVect::Unit); 
   getMetricTerms( N, tanGradN );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box& box = grids[dit];
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(box);

      block_coord_sys.computeMetricTermProductAverage(a_product[dit], a_data[dit], N[dit], box, a_fourthOrder);
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
MagGeom::averageAtBlockBoundaries(LevelData<FluxBox>& a_data) const
{
   if ( m_coord_sys->numBlocks() > 1 ) {

      const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

      RefCountedPtr<MultiBlockCoordSys> coordSysRCP((MultiBlockCoordSys*)m_coord_sys);
      coordSysRCP.neverDelete();

      BlockRegister blockRegister(coordSysRCP, grids, 0);

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
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);

     if ( coord_sys.isAxisymmetric() && (!m_model_geometry) ) {

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

   DataIterator dit = a_BField.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
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

         LevelData<FArrayBox> cell_volume(grids, 1, IntVect::Zero);
         getCellVolumes(cell_volume);

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
            uncorrected_divergence[dit] /= cell_volume[dit];
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



void
MagGeom::writeGeometryData( const DisjointBoxLayout& a_grids,
                            const double&            a_time ) const
{
  IntVect geom_data_ghosts = 4*IntVect::Unit;
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


  Box domain_box = a_grids.physDomain().domainBox();
  domain_box.grow(geom_data_ghosts);

  WriteMappedUGHDF5("config_geom_data", a_grids, geom_data, *m_coord_sys, domain_box, a_time);

  LevelData<FArrayBox> volume(a_grids, 1, geom_data_ghosts);
  getCellVolumes(volume);
  WriteMappedUGHDF5("volume", a_grids, volume, *m_coord_sys, domain_box, a_time);


  LevelData<FArrayBox> J(a_grids, 1, geom_data_ghosts);
  getJ(J);
  WriteMappedUGHDF5("J", a_grids, J, *m_coord_sys, domain_box, a_time);
}



void
MagGeom::checkMultiblockMappingConsistency() const
{
   if ( m_coord_sys->numBlocks() > 1 ) {

      const DisjointBoxLayout& grids = gridsFull();

      LevelData<FluxBox> coords(grids, SpaceDim, IntVect::Zero);
      LevelData<FluxBox> metrics(grids, SpaceDim, IntVect::Zero);

      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = m_coord_sys->boundaries();

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

            if ( box.smallEnd(dir) == domain_box.smallEnd(dir) &&
                 boundaries[block_number][dir].isInterface() ) {

               for (BoxIterator bit(this_coord.box());bit.ok();++bit) {
                  IntVect iv = bit();

                  RealVect xi = dx*iv + offset;
                  RealVect X = block_coord_sys.realCoord(xi);
                  this_coord(iv,0) = X[RADIAL_DIR];
                  this_coord(iv,1) = X[POLOIDAL_DIR];
                  // only makes sense to check transverse derivatives
                  this_metric(iv,0) = block_coord_sys.dXdXi(xi,0,1-dir);
                  this_metric(iv,1) = block_coord_sys.dXdXi(xi,1,1-dir);
               }
            }

            if ( box.bigEnd(dir) == domain_box.bigEnd(dir) &&
                 boundaries[block_number][dir+SpaceDim].isInterface() ) {

               for (BoxIterator bit(this_coord.box());bit.ok();++bit) {
                  IntVect iv = bit();

                  RealVect xi = dx*iv + offset;
                  RealVect X = block_coord_sys.realCoord(xi);
                  this_coord(iv,0) = X[RADIAL_DIR];
                  this_coord(iv,1) = X[POLOIDAL_DIR];
                  // only makes sense to check transverse derivatives
                  this_metric(iv,0) = block_coord_sys.dXdXi(xi,0,1-dir);
                  this_metric(iv,1) = block_coord_sys.dXdXi(xi,1,1-dir);
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

      RefCountedPtr<MultiBlockCoordSys> coordSysRCP((MultiBlockCoordSys*)m_coord_sys);
      coordSysRCP.neverDelete();

      BlockRegister blockRegister(coordSysRCP, grids, 0);

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

   LevelData<FArrayBox> cell_volume(grids, 1, IntVect::Zero);
   getCellVolumes(cell_volume);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      divergence[dit] /= cell_volume[dit];
   }
         
   plotCellData( string("field_divergence"), divergence, a_time );
}


void
MagGeom::fillInternalGhosts( LevelData<FArrayBox>& a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   LevelData<FArrayBox> tmp(grids, a_data.nComp(), m_ghosts*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      tmp[dit].copy(a_data[dit]);
   }

   if (m_mblexPtr && !m_extrablock_exchange) {
      m_mblexPtr->interpGhosts(tmp);
   }
   
   else if (m_extrablock_exchange) {
      exchangeExtraBlockGhosts(a_data);
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_data[dit].copy(tmp[dit]);
   }

   a_data.exchange();
}

void
MagGeom::exchangeExtraBlockGhosts( LevelData<FArrayBox>& a_data ) const
{
    /*
     Fills extrablock ghost cells using a pure exchange from valid data
     rather than interpolation (for experimental purposes only).
     */
    
    if ( m_coord_sys->numBlocks() > 1 ) {
        
        RefCountedPtr<MultiBlockCoordSys> coordSysRCP((MultiBlockCoordSys*)m_coord_sys);
        coordSysRCP.neverDelete();
        
        const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
        
        int ncomp = a_data.nComp();
        
        IntVect ghost_vect = a_data.ghostVect();
        int max_ghost = -1;
        for (int dir=0; dir<SpaceDim; ++dir) {
            if (ghost_vect[dir] > max_ghost) max_ghost = ghost_vect[dir];
        }
        
        BlockRegister blockRegister(coordSysRCP, grids, max_ghost);
        
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            for (int dir = 0; dir < SpaceDim; dir++) {
                IntVect grow_vect = ghost_vect;
                grow_vect[dir] = 0;
                int nghost = ghost_vect[dir];
                for (SideIterator sit; sit.ok(); ++sit) {
                    Side::LoHiSide side = sit();
                    if (blockRegister.hasInterface(dit(), dir, side)) {
                        Box fill_box = adjCellBox(grids[dit], dir, side, 1);
                        fill_box.grow(grow_vect);
                        FArrayBox temp(fill_box, ncomp*nghost);
                        for (int n=0; n<nghost; ++n) {
                            temp.shift(dir, -sign(side));
                            temp.copy(a_data[dit],0,n*ncomp,ncomp);
                        }
                        temp.shiftHalf(dir, (2*nghost-1)*sign(side));
                        blockRegister.storeAux(temp, dit(), dir, side);
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
                    }
                } // iterate over dimensions
            } // iterate over sides
        }
    }
}



void
MagGeom::fillTransversePhysicalGhosts( LevelData<FArrayBox>& a_data ) const
{
   CH_assert(a_data.ghostVect() >= IntVect::Unit);

   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys->boundaries();

   a_data.exchange();

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      const ProblemDomain& block_domain = block_coord_sys.domain();
      const int block_number = m_coord_sys->whichBlock(grids[dit]);
      const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];

      for (int dir=0; dir<SpaceDim; ++dir) {
         IntVect grow_vect = a_data.ghostVect();
         grow_vect[dir] = 0;
         Box interior = grow(grids[dit],grow_vect);

         // If the low or high boundaries in this transverse direction is
         // a block interface, then grow the interior box to trick the
         // extrapolation utility into ignoring those sides
         if ( this_block_boundaries[dir].isInterface() ) {
            interior.growLo(dir,1);
         }
         if ( this_block_boundaries[dir + SpaceDim].isInterface() ) {
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
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys->boundaries();

   a_data.exchange();

   IntVect ghostVect = a_data.ghostVect();

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      const ProblemDomain& block_domain = block_coord_sys.domain();
      const int block_number = m_coord_sys->whichBlock(grids[dit]);
      const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];

      for (int dir=0; dir<SpaceDim; ++dir) {
         Box interior = surroundingNodes(grids[dit],dir);
         interior.grow(dir,ghostVect[dir]);

         for (int tdir=0; tdir<SpaceDim; ++tdir) {
            if (tdir != dir) {  // Transverse directions only

               // If the low or high boundaries in this transverse direction is
               // a block interface, then grow the interior box to trick the
               // extrapolation utility into ignoring those sides
               if ( this_block_boundaries[tdir].isInterface() ) {
                  interior.growLo(tdir,1);
               }
               if ( this_block_boundaries[tdir + SpaceDim].isInterface() ) {
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
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys->boundaries();

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      const ProblemDomain& block_domain = block_coord_sys.domain();
      const int block_number = m_coord_sys->whichBlock(grids[dit]);
      const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];

      for (int dir=0; dir<SpaceDim; ++dir) {
         Box interior = surroundingNodes(grids[dit],dir);

         if ( a_do_block_exchange ) {
            for (int tdir=0; tdir<SpaceDim; ++tdir) {
               if (tdir != dir) {  // Transverse directions only

                  // If the low or high boundaries in this transverse direction is
                  // a block interface, then grow the interior box to trick the
                  // extrapolation utility into ignoring those sides
                  if ( this_block_boundaries[tdir].isInterface() ) {
                     interior.growLo(tdir,1);
                  }
                  if ( this_block_boundaries[tdir + SpaceDim].isInterface() ) {
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

void MagGeom::interpolateErFromMagFS(LevelData<FluxBox>&   a_Er_face,
                                     LevelData<FArrayBox>& a_Er_cell) const

{
    //Here, we consider the model geometry, in which Ez along the top cut is the radial component
    //Needs to be generalized for the case of a DIII-D geometry (by extracting the actual radial component)
    //Also needs to be generalized for the 10 block structure
    
    //X-point vicinity, in which we address the dealignemnt
    double Xpt_vicinity = DBL_MAX;
    double Zaxis = 1.75935; //Hardwired value for the model geometry
    double R0 = 1.6; //Hardwired value for the model geometry
    
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
                double delta_R = face_center_coords_dir(iv,0) - ((SingleNullCoordSys*)m_coord_sys)->m_Xpoint[0];
                double delta_Z = face_center_coords_dir(iv,1) - ((SingleNullCoordSys*)m_coord_sys)->m_Xpoint[1];
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
                    
                    if (block_number>RCSOL) this_E_face(iv,2) *= -1.0;
                    
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
            double delta_R = cell_center_coords(iv,0) - ((SingleNullCoordSys*)m_coord_sys)->m_Xpoint[0];
            double delta_Z = cell_center_coords(iv,1) - ((SingleNullCoordSys*)m_coord_sys)->m_Xpoint[1];
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
    const MagBlockCoordSys& lcore_coord_sys = getBlockCoordSys(LCORE);
    RealVect dx_core = lcore_coord_sys.dx();
    int nrad_core = lcore_coord_sys.domain().domainBox().size(RADIAL_DIR);
    
    const MagBlockCoordSys& lcsol_coord_sys = getBlockCoordSys(LCSOL);
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
    const MagBlockCoordSys& lcore_coord_sys = getBlockCoordSys(LCORE);
    RealVect dx_core = lcore_coord_sys.dx();
    int nrad_core = lcore_coord_sys.domain().domainBox().size(RADIAL_DIR);
    
    const MagBlockCoordSys& lcsol_coord_sys = getBlockCoordSys(LCSOL);
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
    
    const MagBlockCoordSys& lcore_coord_sys = getBlockCoordSys(LCORE);
    int nrad_core  = lcore_coord_sys.domain().domainBox().size(RADIAL_DIR);
    int lo_pol_LCORE = lcore_coord_sys.domain().domainBox().smallEnd(POLOIDAL_DIR);
    
    const MagBlockCoordSys& lcsol_coord_sys = getBlockCoordSys(LCSOL);
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
        if (block_number == LCORE) {
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
        
        if (block_number == LCSOL) {
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
    
    const MagBlockCoordSys& lcore_coord_sys = getBlockCoordSys(LCORE);
    int nrad_core  = lcore_coord_sys.domain().domainBox().size(RADIAL_DIR);
    int lo_pol_LCORE = lcore_coord_sys.domain().domainBox().smallEnd(POLOIDAL_DIR);
    
    const MagBlockCoordSys& lcsol_coord_sys = getBlockCoordSys(LCSOL);
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
        if (block_number == LCORE) {
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
        
        if (block_number == LCSOL) {
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
MagGeom::getMagFS( const RealVect& a_xi ) const
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
    const MagBlockCoordSys& lcore_coord_sys =  getBlockCoordSys(LCORE);
    RealVect dx_core = lcore_coord_sys.dx();
    double offset_core = 3.0*dx_core[0];
    RealVect mapped_coord_loZ;
    mapped_coord_loZ[1] = lcore_coord_sys.lowerMappedCoordinate(1);
    mapped_coord_loZ[0] = lcore_coord_sys.lowerMappedCoordinate(0) - offset_core;
    Real lo_Z = lcore_coord_sys.realCoord(mapped_coord_loZ)[1];
    
    const MagBlockCoordSys& lcsol_coord_sys =  getBlockCoordSys(LCSOL);
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
    double magFluxValue = magFlux(a_xi);
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

double MagGeom::magFlux(const RealVect& a_xi) const
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
    
    if (a_xi[1]>Zxpt) {
        val = cos(a*(a_xi[0]-R0)) + b*sin((a_xi[1]-Z0)) - c*(a_xi[1]-Z0);
    }
    else {
        val = cos(a*(a_xi[0]-R0)) + b*sin((2.0*Zxpt-a_xi[1]-Z0)) - c*(2.0*Zxpt-a_xi[1]-Z0);
    }
    
    return val;
}

double MagGeom::magFluxZ(const double a_Z) const
{
    //Presently returns the magnetic flux function for the model geometry case
    //LATER: generalize by making use of the DCT representation of psi
   //    double a=1.2;
    double b=0.9;
    double c=0.7;
    //    double R0=1.6;
    double Z0=1.07967;
    double Zxpt= 0.4;
    double val;
    
    if (a_Z>Zxpt) {
        val = 1.0 + b*sin((a_Z-Z0)) - c*(a_Z-Z0);
    }
    else {
        val = 1.0 + b*sin((2.0*Zxpt-a_Z-Z0)) - c*(2.0*Zxpt-a_Z-Z0);
    }
    
    return val;
}



void
MagGeom::unmapGradient( const LevelData<FArrayBox>& a_mapped_gradient,
                        LevelData<FArrayBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 3);
   CH_assert(a_gradient.nComp() == 3);

   LevelData<FArrayBox> mapped_poloidal_gradient(m_gridsFull, 2, a_mapped_gradient.ghostVect());
   projectPoloidalVector(a_mapped_gradient, mapped_poloidal_gradient);

   LevelData<FArrayBox> poloidal_gradient(m_gridsFull, 2, a_gradient.ghostVect());
   unmapPoloidalGradient(mapped_poloidal_gradient, poloidal_gradient);

   injectPoloidalVector(poloidal_gradient, a_gradient);
}



void
MagGeom::unmapGradient( const LevelData<FluxBox>& a_mapped_gradient,
                        LevelData<FluxBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 3);
   CH_assert(a_gradient.nComp() == 3);

   LevelData<FluxBox> mapped_poloidal_gradient(m_gridsFull, 2, a_mapped_gradient.ghostVect());
   projectPoloidalVector(a_mapped_gradient, mapped_poloidal_gradient);

   LevelData<FluxBox> poloidal_gradient(m_gridsFull, 2, a_gradient.ghostVect());
   unmapPoloidalGradient(mapped_poloidal_gradient, poloidal_gradient);

   injectPoloidalVector(poloidal_gradient, a_gradient);
}



void
MagGeom::unmapPoloidalGradient( const LevelData<FArrayBox>& a_mapped_gradient,
                                LevelData<FArrayBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 2);
   CH_assert(a_gradient.nComp() == 2);

   // Multiply by NJInverse
   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(m_gridsFull[dit]);

      Box box = a_mapped_gradient[dit].box();
      FArrayBox NJInverse(box, 4);
      block_coord_sys.getPointwiseNJInverse(NJInverse);

      FORT_MULT_NJINVERSE(CHF_BOX(box),
                          CHF_CONST_FRA(a_mapped_gradient[dit]),
                          CHF_CONST_FRA(NJInverse),
                          CHF_FRA(a_gradient[dit]));
   }
}



void
MagGeom::unmapPoloidalGradient( const LevelData<FluxBox>& a_mapped_gradient,
                                LevelData<FluxBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 2);
   CH_assert(a_gradient.nComp() == 2);

   // Multiply by NJInverse
   for (DataIterator dit(m_gridsFull); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(m_gridsFull[dit]);

      FluxBox NJInverse(a_mapped_gradient[dit].box(), 4);
      block_coord_sys.getPointwiseNJInverse(NJInverse);

      for (int dir=0; dir<2; dir++) {
         FORT_MULT_NJINVERSE(CHF_BOX(a_mapped_gradient[dit][dir].box()),
                             CHF_CONST_FRA(a_mapped_gradient[dit][dir]),
                             CHF_CONST_FRA(NJInverse[dir]),
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
      for (int dir=0; dir<2; ++dir) {
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
   CH_assert(a_phi.ghostVect() >= 2*IntVect::Unit);
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);
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
         FORT_CELL_CENTERED_GRAD_COMPONENT(CHF_BOX(box),
                                           CHF_CONST_INT(dir),
                                           CHF_CONST_FRA1(a_phi[dit],0),
                                           CHF_CONST_REALVECT(dx),
                                           CHF_CONST_INT(tmp_order),
                                           CHF_FRA1(a_field[dit],dir));
      }
   }
      
   a_field.exchange();
      
   // If fourth-order, recompute the field at valid cell centers to fourth-order
      
   if ( a_order == 4 ) {
         
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
         RealVect dx = block_coord_sys.dx();
            
         for (int dir=0; dir<2; ++dir) {
            FORT_CELL_CENTERED_GRAD_COMPONENT(CHF_BOX(grids[dit]),
                                              CHF_CONST_INT(dir),
                                              CHF_CONST_FRA1(a_phi[dit],0),
                                              CHF_CONST_REALVECT(dx),
                                              CHF_CONST_INT(a_order),
                                              CHF_FRA1(a_field[dit],dir));
         }
      }
   }
}

void
MagGeom::computeMappedPoloidalGradientWithGhosts( const LevelData<FArrayBox>& a_phi,
                                                  LevelData<FluxBox>&         a_field,
                                                  const int                   a_order ) const
{
   CH_assert(a_phi.ghostVect() >= 2*IntVect::Unit);
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);
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
         
         FORT_FACE_INTERPOLATE(CHF_CONST_INT(dir),
                               CHF_BOX(surroundingNodes(box,dir)),
                               CHF_CONST_INT(a_order),
                               CHF_CONST_FRA1(a_phi[dit],0),
                               CHF_FRA1(phi_face[dit][dir],0));
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
         FORT_FACE_CENTERED_GRAD_COMPONENT(CHF_BOX(box_dir),
                                            CHF_CONST_INT(dir),
                                            CHF_CONST_FRA1(a_phi[dit],0),
                                            CHF_CONST_REALVECT(dx),
                                            CHF_CONST_INT(tmp_order),
                                            CHF_FRA1(a_field[dit][dir],dir));
         
         for (int tdir=0; tdir<2; ++tdir) {
            if (tdir != dir) {
               FORT_CELL_CENTERED_GRAD_COMPONENT(CHF_BOX(box_dir),
                                                  CHF_CONST_INT(tdir),
                                                  CHF_CONST_FRA1(phi_face[dit][dir],0),
                                                  CHF_CONST_REALVECT(dx),
                                                  CHF_CONST_INT(tmp_order),
                                                  CHF_FRA1(a_field[dit][dir],tdir));
            }
         }
      }
   }
   a_field.exchange();
   
   // Second-order extrapolate the field to the transverse physical boundary ghosts
   fillTransversePhysicalGhosts(a_field);
   
   // If fourth-order, recompute the field on valid cell faces
   
   if ( a_order == 4 ) {
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
         RealVect dx = block_coord_sys.dx();
         
         for (int dir=0; dir<2; dir++) {
            Box box_dir = surroundingNodes(grids[dit],dir);
            FORT_FACE_CENTERED_GRAD_COMPONENT(CHF_BOX(box_dir),
                                              CHF_CONST_INT(dir),
                                              CHF_CONST_FRA1(a_phi[dit],0),
                                              CHF_CONST_REALVECT(dx),
                                              CHF_CONST_INT(a_order),
                                              CHF_FRA1(a_field[dit][dir],dir));
            
            for (int tdir=0; tdir<2; ++tdir) {
               if (tdir != dir) {
                  FORT_CELL_CENTERED_GRAD_COMPONENT(CHF_BOX(box_dir),
                                                    CHF_CONST_INT(tdir),
                                                    CHF_CONST_FRA1(phi_face[dit][dir],0),
                                                    CHF_CONST_REALVECT(dx),
                                                    CHF_CONST_INT(a_order),
                                                    CHF_FRA1(a_field[dit][dir],tdir));
               }
            }
         }
      }
   }
   
   a_field.exchange();
}
   
#include "NamespaceFooter.H"
