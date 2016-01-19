#include "MillerCoordSys.H"
#include "MillerBlockCoordSys.H"
#include "FourthOrderUtil.H"
#include "BlockBoundary.H"
#include "SpreadingCopier.H"
#include "ReductionOps.H"
#include "Directions.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

bool MillerCoordSys::isPeriodic( const int a_dir ) const
{
   return m_block_zero_domain.isPeriodic( a_dir );
}

MillerCoordSys::MillerCoordSys( ParmParse&               a_pp,
                                const std::vector<int>&  a_numcells,
                                const std::vector<bool>& a_is_periodic,
                                const std::vector<int>&  a_decomposition )
{
  IntVect lo = IntVect::Zero;
  IntVect hi;
  bool is_periodic[SpaceDim];
  for (int dir=0; dir<SpaceDim; ++dir) {
    hi[dir] = a_numcells[dir] - 1;
    is_periodic[dir] = a_is_periodic[dir];
    m_decomposition[dir] = a_decomposition[dir];
  }

  m_block_zero_domain.define(Box(lo,hi), is_periodic);

  MillerBlockCoordSys* block_coords = new MillerBlockCoordSys( a_pp, m_block_zero_domain );

  Vector<MagBlockCoordSys *> coord_vec;
  coord_vec.push_back(block_coords);

  defineCoordSystemsAndBoundaries(coord_vec);

  m_spread_radially = block_coords->getConstMinorrad() != 0;

  int bc_type = 0;     // Homogeneous Dirichlet
  int bc_value = 0.;
  for (int dir=0; dir<SpaceDim; ++dir) {
     if ( !m_block_zero_domain.isPeriodic(dir) ) {
        for (int side=0; side<2; ++side) {
           m_divergence_cleaning_bc.setBCType(0, dir, side, bc_type);
           m_divergence_cleaning_bc.setBCValue(0, dir, side, bc_value);
        }
     }
  }
}



MillerCoordSys::~MillerCoordSys()
{
}



void
MillerCoordSys::defineBoundaries()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   for (int iblock(0); iblock<numBlocks(); iblock++) {

      const MagBlockCoordSys& coord_sys( *(MagBlockCoordSys*)getCoordSys(iblock) );
      const ProblemDomain& domain( coord_sys.domain() );
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries( m_boundaries[iblock] );
      for (int idir(RADIAL_DIR); idir<SpaceDim; idir++) {
         if ( domain.isPeriodic(idir) ) {
//       if (m_block_zero_domain.isPeriodic(idir)) {
            IntVect shift( m_mappingBlocks[iblock].size(idir) * BASISV(idir) );
            IndicesTransformation it;
            it.defineFromTranslation( shift );
            blockBoundaries[idir].define( it, iblock );
            it.defineFromTranslation( -shift );
            blockBoundaries[idir+SpaceDim].define( it, iblock );
         }
         else {
//          blockBoundaries[idir] = BlockBoundary::BOUNDARY;
//          blockBoundaries[idir+SpaceDim] = BlockBoundary::BOUNDARY;
            blockBoundaries[idir].define(0);
            blockBoundaries[idir+SpaceDim].define(0);
         }
      }
   }
}



void
MillerCoordSys::blockRemapping(RealVect&       a_xi_valid,
                               int&            a_n_valid,
                               const RealVect& a_xiSrc,
                               int             a_nSrc) const
{
   const MillerBlockCoordSys* mag_coord_sys = (MillerBlockCoordSys*)getCoordSys(0);

   double radial_width = mag_coord_sys->getOuterFluxLabel() - mag_coord_sys->getInnerFluxLabel();
   double poloidal_width = 2. * Pi;
#if CFG_DIM==3
   // NEW FOR 5D:
   double toroidal_width = mag_coord_sys->getMaxTorAngle() - mag_coord_sys->getMinTorAngle();
#endif

   a_n_valid = a_nSrc;
   a_xi_valid = a_xiSrc;

   if (a_xiSrc[RADIAL_DIR] < 0.)
      a_xi_valid[RADIAL_DIR] += radial_width;
   else if (a_xiSrc[RADIAL_DIR] > radial_width)
      a_xi_valid[RADIAL_DIR] -= radial_width;

   if (a_xiSrc[POLOIDAL_DIR] < 0.)
      a_xi_valid[POLOIDAL_DIR] += poloidal_width;
   else if (a_xiSrc[POLOIDAL_DIR] > poloidal_width)
      a_xi_valid[POLOIDAL_DIR] -= poloidal_width;

#if CFG_DIM==3
   // NEW FOR 5D:
   if (a_xiSrc[TOROIDAL_DIR] < 0.)
      a_xi_valid[TOROIDAL_DIR] += toroidal_width;
   else if (a_xiSrc[TOROIDAL_DIR] > toroidal_width)
      a_xi_valid[TOROIDAL_DIR] -= toroidal_width;
#endif

}



void
MillerCoordSys::postProcessMetricData( LevelData<FArrayBox>& a_data ) const
{
  if (m_spread_radially) {
    spreadRadially(a_data);
  }
}



void
MillerCoordSys::postProcessMetricData( LevelData<FluxBox>& a_data ) const
{
  if (m_spread_radially) {
    spreadRadially(a_data);
  }
}



void
MillerCoordSys::spreadRadially( LevelData<FArrayBox>& a_data ) const
{
   int radial_dir = 0;
   Box domainBox = m_block_zero_domain.domainBox();
   int rmid = (domainBox.smallEnd(radial_dir) + domainBox.bigEnd(radial_dir) + 1)/2;

   LevelData<FArrayBox> grown_data(a_data.disjointBoxLayout(), a_data.nComp(), a_data.ghostVect()+2*IntVect::Unit);

   DataIterator dit = a_data.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     grown_data[dit].copy(a_data[dit]);
   }
   grown_data.exchange();

   LevelData<FluxBox> face_data(grown_data.disjointBoxLayout(), grown_data.nComp(), a_data.ghostVect());
   fourthOrderCellToFace(face_data, grown_data);

   LevelData<FArrayBox> data2(a_data.disjointBoxLayout(), a_data.nComp(), a_data.ghostVect());
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& radial_face = face_data[dit][radial_dir];
      radial_face.shiftHalf(radial_dir,1);
      data2[dit].copy(radial_face);
   }

   // Make a poloidal slice at the radial midpoint
   //LevelData<FArrayBox> temp;
   const SliceSpec radial_slice(radial_dir,rmid);
   //a_data.degenerate(temp, radial_slice);
   // do this explicitly here
   //   DisjointBoxLayout sliceGrids;
   Box sliceBox(domainBox);
   sliceBox.setSmall(radial_dir,rmid);
   sliceBox.setBig(radial_dir, rmid);
   Vector<Box> sliceBoxes(1,sliceBox);
   Vector<int> bogusProcAssign(1,0);
   DisjointBoxLayout sliceGrids(sliceBoxes, bogusProcAssign, m_block_zero_domain);
#if 0
   const DisjointBoxLayout& grids = a_data.getBoxes();
   grids.degenerate(sliceGrids, radial_slice);
#endif

   LevelData<FArrayBox> temp(sliceGrids, a_data.nComp(), a_data.ghostVect());
   data2.copyTo(temp);

   // Spread the slice radially
   Vector<int> spreadingDirs;
   spreadingDirs.push_back(radial_dir);

   SpreadingCopier spreadCopier(temp.getBoxes(),
                                a_data.getBoxes(),
                                m_block_zero_domain,
                                a_data.ghostVect(),
                                spreadingDirs);
   const SpreadingOp spreadOp(spreadingDirs);

   temp.copyTo(temp.interval(),
               a_data, a_data.interval(),
               spreadCopier,
               spreadOp);
}



void
MillerCoordSys::spreadRadially( LevelData<FluxBox>& a_data ) const
{
   int radial_dir = 0;
   Box domainBox = m_block_zero_domain.domainBox();
   int rmid = (domainBox.smallEnd(radial_dir) + domainBox.bigEnd(radial_dir) + 1)/2;

   LevelData<FluxBox> grown_data(a_data.disjointBoxLayout(), a_data.nComp(), a_data.ghostVect()+2*IntVect::Unit);

   DataIterator dit = a_data.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     grown_data[dit].copy(a_data[dit]);
   }
   grown_data.exchange();

   // Temporaries
   LevelData<FArrayBox> tmp_cell(grown_data.disjointBoxLayout(), grown_data.nComp(), grown_data.ghostVect());
   LevelData<FluxBox> tmp_face(a_data.disjointBoxLayout(), a_data.nComp(), a_data.ghostVect());

   // Average data on non-radial faces.  We do this by shifting to cell
   // centers, apply the 4th-order averaging, then copy the resulting
   // average on the middle flux surface
   for (int dir=0; dir<SpaceDim; ++dir) {
      if (dir != radial_dir) {
         for (int side=-1; side<2; side+=2) {
            for (dit.begin(); dit.ok(); ++dit) {
               FArrayBox& face(grown_data[dit][dir]);
               face.shiftHalf(dir,-side);
               tmp_cell[dit].copy(face);
               face.shiftHalf(dir,side);  // put it back
            }
            tmp_cell.exchange();

            fourthOrderCellToFace(tmp_face, tmp_cell);

            for (dit.begin(); dit.ok(); ++dit) {
               FArrayBox& face(tmp_face[dit][radial_dir]);
               face.shiftHalf(radial_dir,1);
               face.shiftHalf(dir,side);
               a_data[dit][dir].copy(face);
               face.shiftHalf(dir,-side); // put it back
               face.shiftHalf(radial_dir,-1);  // put it back
            }
         }
      }
   }
   a_data.exchange();

   // Make a poloidal slice at the radial midpoint
   //LevelData<FluxBox> temp;
   const SliceSpec radial_slice(radial_dir,rmid);
   //a_data.degenerae(temp, radial_slice);

   // do this explicitly here
#if 0
   DisjointBoxLayout sliceGrids;
   const DisjointBoxLayout& grids = a_data.getBoxes();
   grids.degenerate(sliceGrids, radial_slice);
#endif

   Box sliceBox(domainBox);
   sliceBox.setSmall(radial_dir,rmid);
   sliceBox.setBig(radial_dir, rmid);
   Vector<Box> sliceBoxes(1,sliceBox);
   Vector<int> bogusProcAssign(1,0);
   DisjointBoxLayout sliceGrids(sliceBoxes, bogusProcAssign, m_block_zero_domain);

   LevelData<FluxBox> temp(sliceGrids, a_data.nComp(), a_data.ghostVect());
   a_data.copyTo(temp);

   // Spread the slice radially
   Vector<int> spreadingDirs;
   spreadingDirs.push_back(radial_dir);

   SpreadingCopier spreadCopier(temp.getBoxes(),
                                a_data.getBoxes(),
                                m_block_zero_domain,
                                a_data.ghostVect(),
                                spreadingDirs);
   const FaceSpreadingOp spreadOp(spreadingDirs);

   temp.copyTo(temp.interval(),
               a_data, a_data.interval(),
               spreadCopier,
               spreadOp);
}



#include "NamespaceFooter.H"
