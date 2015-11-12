#include "MagGeom.H"
#include "Poisson.H"
#include "Directions.H"
#include "CONSTANTS.H"
#include "FourthOrderUtil.H"
#include "newMappedGridIO.H"
#include "SimpleDivergence.H"
#include "EdgeToCell.H"
#include "inspect.H"

#include "NamespaceHeader.H"

#undef TEST_DIVERGENCE_CLEANING


MagGeom::MagGeom(ParmParse&                 a_pp,
                 const MagCoordSys*         a_coordSysPtr,
                 const DisjointBoxLayout&   a_grids,
                 int                        a_ghosts)
   : MultiBlockLevelGeom(a_coordSysPtr, a_grids, a_ghosts),
     m_pp_prefix(a_pp.prefix()),
     m_coord_sys(a_coordSysPtr),
     m_correct_field(false)
{
   if (a_pp.contains("correct_field")) {
     a_pp.get("correct_field", m_correct_field);
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
         plotFieldAlignment();
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

      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
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
MagGeom::getNJInverse( LevelData<FluxBox>& a_NJInverse ) const
{
   const DisjointBoxLayout& grids = a_NJInverse.disjointBoxLayout();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);

      Box box(a_NJInverse[dit].box());

      Box grown_box(box);
      grown_box.grow(1);  // Need and extra ghost cell for the transverse gradients in the product formula

      FluxBox N(grown_box, SpaceDim*SpaceDim);
      coord_sys.getN(N, grown_box);

      FluxBox JInverse(grown_box,1);
      coord_sys.getAvgJinverse(JInverse, grown_box);
      coord_sys.computeNJinverse(a_NJInverse[dit], JInverse, N, box);
   }

   m_coord_sys->postProcessMetricData(a_NJInverse);

   a_NJInverse.exchange();
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


void MagGeom::multPointwiseJ( LevelData<FArrayBox>& a_u,
                              const BoundaryBoxLayout& a_bdry_layout ) const
{
   const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
   for (DataIterator dit( bdry_grids.dataIterator() ); dit.ok(); ++dit) {
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
   CH_assert(a_data.ghostVect() >= IntVect::Unit);
   const DisjointBoxLayout& grids = a_product.disjointBoxLayout();

   LevelData<FluxBox> N(grids, SpaceDim*SpaceDim, IntVect::Unit); 
   LevelData<FluxBox> tanGradN(grids, SpaceDim*SpaceDim, IntVect::Unit); 
   getMetricTerms( N, tanGradN );

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
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
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
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

      DataIterator dit = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
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

      for (dit.begin(); dit.ok(); ++dit) {
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
      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys->boundaries();

      for ( int tdir=0; tdir<SpaceDim; ++tdir ) {   // Transverse directions

         m_exchange_transverse_block_register->setToZero( a_data.nComp() );

         for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
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

         for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
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
   DataIterator dit = grids.dataIterator();

   for (dit.begin(); dit.ok(); ++dit) {
     const MagBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);

     if ( coord_sys.isAxisymmetric() ) {

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
   DataIterator dit = a_BField.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
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
      plotFaceData( string("before_correction_BField"), a_BField );
      plotFaceData( string("before_correction_BFieldMag"), a_BFieldMag );
      plotFaceData( string("before_correction_BFieldDir"), a_BFieldDir );
      plotFaceData( string("before_correction_gradB"), a_gradBFieldMag );
      plotFaceData( string("before_correction_curlBFieldDir"), a_curlBFieldDir );
      plotFaceData( string("before_correction_BFieldDirdotcurlBFieldDir"), a_BFieldDirdotcurlBFieldDir );
#endif

      if ( !m_BFieldCorrection.isDefined() ) {
         m_BFieldCorrection.define( grids, 3, IntVect::Zero);

         CH_assert(a_BField.ghostVect() >= IntVect::Unit);
         LevelData<FluxBox> flux( grids, SpaceDim, IntVect::Unit );
         LevelData<FArrayBox> uncorrected_divergence( grids, 1, IntVect::Zero );
#ifdef TEST_DIVERGENCE_CLEANING
         LevelData<FArrayBox> corrected_divergence( grids, 1, IntVect::Zero );
         Box domain_box = grids.physDomain().domainBox();
#endif

         LevelData<FArrayBox> cell_volume(grids, 1, IntVect::Zero);
         getCellVolumes(cell_volume);

         for (dit.begin(); dit.ok(); ++dit) {
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

         fillTransverseGhosts(flux, false);

         // Convert to face-averaged
         fourthOrderAverage(flux);

         fillTransverseGhosts(flux, false);

         computeMappedGridDivergence( flux, uncorrected_divergence, true );

         for (dit.begin(); dit.ok(); ++dit) {
            uncorrected_divergence[dit] /= cell_volume[dit];
         }
         
#ifdef TEST_DIVERGENCE_CLEANING
         WriteMappedUGHDF5("uncorrected_divergence", grids, uncorrected_divergence, *m_coord_sys, domain_box);
#endif

         ParmParse pp( Poisson::pp_name );
         Poisson * poisson = new Poisson( pp, *this );
         poisson->setOperatorCoefficients( m_coord_sys->getDivergenceCleaningBC() );

         LevelData<FArrayBox> phi( grids, 1, 4*IntVect::Unit );
         for (dit.begin(); dit.ok(); ++dit) {
            phi[dit].setVal(0.);
         }

         const PotentialBC& bc = m_coord_sys->getDivergenceCleaningBC();

         poisson->solve( uncorrected_divergence, bc, phi );

#ifdef TEST_DIVERGENCE_CLEANING
         WriteMappedUGHDF5("field_correction_potential", grids, phi, *m_coord_sys, domain_box);
#endif

         poisson->fillGhosts(phi);

         // Compute cell-centered -grad phi
         poisson->computeFaceCenteredField( phi, bc, flux );

         fillTransverseGhosts(flux, false);

         for (dit.begin(); dit.ok(); ++dit) {
            FluxBox& this_BFieldCorrection = m_BFieldCorrection[dit];
            FluxBox& this_flux = flux[dit];

            for (int dir=0; dir<SpaceDim; ++dir) {
               FArrayBox& this_BFieldCorrection_dir = this_BFieldCorrection[dir];
               FArrayBox& this_flux_dir = this_flux[dir];

               this_BFieldCorrection_dir.copy( this_flux_dir, 0, 0, 1 );
               this_BFieldCorrection_dir.setVal(0., 1);
               this_BFieldCorrection_dir.copy( this_flux_dir, 1, 2, 1 );
            }
         }

         for (dit.begin(); dit.ok(); ++dit) {
            a_BField[dit] -= m_BFieldCorrection[dit];
         }

#ifdef TEST_DIVERGENCE_CLEANING

         {
            for (dit.begin(); dit.ok(); ++dit) {
               FluxBox& this_BField = a_BField[dit];
               FluxBox& this_flux = flux[dit];
               for (int dir=0; dir<SpaceDim; ++dir) {
                  FArrayBox& this_BField_dir = this_BField[dir];
                  FArrayBox& this_flux_dir = this_flux[dir];
                  
                  this_flux_dir.copy( this_BField_dir, 0, 0, 1 );
                  this_flux_dir.copy( this_BField_dir, 2, 1, 1 );
               }
            }

            applyAxisymmetricCorrection( flux );

            fillTransverseGhosts(flux, false);

            fourthOrderAverage(flux);

            fillTransverseGhosts(flux, false);

            computeMappedGridDivergence( flux, corrected_divergence, true );

            for (dit.begin(); dit.ok(); ++dit) {
               corrected_divergence[dit] /= cell_volume[dit];
            }

            WriteMappedUGHDF5("corrected_divergence", grids, corrected_divergence, *m_coord_sys, domain_box);
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
         for (dit.begin(); dit.ok(); ++dit) {
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
         for (dit.begin(); dit.ok(); ++dit) {
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

         for (dit.begin(); dit.ok(); ++dit) {
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
         for (dit.begin(); dit.ok(); ++dit) {
            FluxBox& this_change = BFieldMagChange[dit];

            FluxBox denom(this_change.box(), this_change.nComp());
            denom.copy(this_change);
            this_change -= a_BFieldMag[dit];
            for (int dir=0; dir<SpaceDim; ++dir) {
               this_change[dir] /= denom[dir];
            }
         }

         for (dit.begin(); dit.ok(); ++dit) {
            BFieldMagChange[dit] -= a_BFieldMag[dit];
            BFieldDirChange[dit] -= a_BFieldDir[dit];
            gradBFieldMagChange[dit] -= a_gradBFieldMag[dit];
            curlBFieldDirChange[dit] -= a_curlBFieldDir[dit];
            BFieldDirdotcurlBFieldDirChange[dit] -= a_BFieldDirdotcurlBFieldDir[dit];
         }

         plotFaceData( string("BFieldChange"), m_BFieldCorrection );
         plotFaceData( string("BFieldMagChange"), BFieldMagChange );
         plotFaceData( string("BFieldDirChange"), BFieldDirChange );
         plotFaceData( string("gradBFieldMagChange"), gradBFieldMagChange );
         plotFaceData( string("curlBFieldDirChange"), curlBFieldDirChange );
         plotFaceData( string("BFieldDirdotcurlBFieldDirChange"), BFieldDirdotcurlBFieldDirChange );

         plotFaceData( string("after_correction_BField"), a_BField );
         plotFaceData( string("after_correction_BFieldMag"), a_BFieldMag );
         plotFaceData( string("after_correction_BFieldDir"), a_BFieldDir );
         plotFaceData( string("after_correction_gradB"), a_gradBFieldMag );
         plotFaceData( string("after_correction_curlBFieldDir"), a_curlBFieldDir );
         plotFaceData( string("after_correction_BFieldDirdotcurlBFieldDir"), a_BFieldDirdotcurlBFieldDir );
#endif

         if( procID() == 0 ) {
            cout << endl;
         }
      }
   }
}



void MagGeom::plotCellData( const string&               a_file_name,  
                            const LevelData<FArrayBox>& a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   Box domain_box = grids.physDomain().domainBox();
   domain_box.grow(a_data.ghostVect());
   WriteMappedUGHDF5(a_file_name.c_str(), grids, a_data, *m_coord_sys, domain_box);
}



void MagGeom::plotFaceData( const string&             a_file_name,  
                            const LevelData<FluxBox>& a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   // Since the passed data is a LevelData<FluxBox>, we average the face
   // values to cell centers
   LevelData<FArrayBox> tmp(grids, SpaceDim*a_data.nComp(), a_data.ghostVect());
   EdgeToCell(a_data, tmp);

   // Next, average the directions
   LevelData<FArrayBox> data_cell(grids, a_data.nComp(), a_data.ghostVect());

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
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
   WriteMappedUGHDF5(a_file_name.c_str(), grids, data_cell, *m_coord_sys, domain_box);
}



void
MagGeom::writeGeometryData( const DisjointBoxLayout& grids ) const
{
  IntVect geom_data_ghosts = 4*IntVect::Unit;
  LevelData<FArrayBox> geom_data(grids, 6, geom_data_ghosts);

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    const MagBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);

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


  Box domain_box = grids.physDomain().domainBox();
  domain_box.grow(geom_data_ghosts);

  WriteMappedUGHDF5("config_geom_data", grids, geom_data, *m_coord_sys, domain_box);

  LevelData<FArrayBox> volume(grids, 1, geom_data_ghosts);
  getCellVolumes(volume);
  WriteMappedUGHDF5("volume", grids, volume, *m_coord_sys, domain_box);


  LevelData<FArrayBox> J(grids, 1, geom_data_ghosts);
  getJ(J);
  WriteMappedUGHDF5("J", grids, J, *m_coord_sys, domain_box);
}



void
MagGeom::checkMultiblockMappingConsistency() const
{
   if ( m_coord_sys->numBlocks() > 1 ) {

      const DisjointBoxLayout& grids = gridsFull();

      LevelData<FluxBox> coords(grids, SpaceDim, IntVect::Zero);
      LevelData<FluxBox> metrics(grids, SpaceDim, IntVect::Zero);

      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = m_coord_sys->boundaries();

      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
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

      DataIterator dit = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
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

      for (dit.begin(); dit.ok(); ++dit) {
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
MagGeom::plotFieldAlignment() const
{
   const DisjointBoxLayout& grids = gridsFull();
   LevelData<FluxBox> dotprod(grids, 1, IntVect::Zero);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
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

   plotFaceData( string("field_alignment"), dotprod );
}


void
MagGeom::fillInternalGhosts( LevelData<FArrayBox>& a_data ) const
{
   CH_assert(a_data.ghostVect() >= m_ghosts*IntVect::Unit);

   a_data.exchange();

   if (m_mblexPtr) {
     m_mblexPtr->interpGhosts(a_data);
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

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
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



#include "NamespaceFooter.H"
