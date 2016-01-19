#include "RectangularTorusCoordSys.H"
#include "RectangularTorusBlockCoordSys.H"
#include "FourthOrderUtil.H"
#include "BlockBoundary.H"
#include "SpreadingCopier.H"
#include "ReductionOps.H"
#include "RectangularTorusDirections.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

bool RectangularTorusCoordSys::isPeriodic( const int a_dir ) const
{
   return m_block_zero_domain.isPeriodic( a_dir );
}

RectangularTorusCoordSys::RectangularTorusCoordSys( ParmParse&               a_pp,
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

  RectangularTorusBlockCoordSys* block_coords = new RectangularTorusBlockCoordSys( a_pp, m_block_zero_domain );

  Vector<MagBlockCoordSys *> coord_vec;
  coord_vec.push_back(block_coords);

  defineCoordSystemsAndBoundaries(coord_vec);

  Tuple<BlockBoundary, 2*SpaceDim>& block_boundaries( m_boundaries[0] );

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



RectangularTorusCoordSys::~RectangularTorusCoordSys()
{
    //Shouldn't I delete block_coords?
}



void
RectangularTorusCoordSys::defineBoundaries()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   for (int iblock(0); iblock<numBlocks(); iblock++) {

      const MagBlockCoordSys& coord_sys( *(MagBlockCoordSys*)getCoordSys(iblock) );
      const ProblemDomain& domain( coord_sys.domain() );
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries( m_boundaries[iblock] );
      for (int idir(R_DIR); idir<SpaceDim; idir++) {
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
RectangularTorusCoordSys::blockRemapping(RealVect&       a_xi_valid,
                               int&            a_n_valid,
                               const RealVect& a_xiSrc,
                               int             a_nSrc) const
{
   const RectangularTorusBlockCoordSys* mag_coord_sys = (RectangularTorusBlockCoordSys*)getCoordSys(0);

   double r_width = mag_coord_sys->getOuterFluxLabel() - mag_coord_sys->getInnerFluxLabel();
   double z_width = mag_coord_sys->getUpperZLabel() - mag_coord_sys->getLowerZLabel();
#if CFG_DIM==3
   // NEW FOR 5D:
   double phi_width = mag_coord_sys->getMaxTorAngle() - mag_coord_sys->getMinTorAngle();
#endif

   a_n_valid = a_nSrc;
   a_xi_valid = a_xiSrc;

   if (a_xiSrc[R_DIR] < 0.)
      a_xi_valid[R_DIR] += r_width;
   else if (a_xiSrc[R_DIR] > r_width)
      a_xi_valid[R_DIR] -= r_width;


#if CFG_DIM==3
   // NEW FOR 5D:
   if (a_xiSrc[PHI_DIR] < 0.)
      a_xi_valid[PHI_DIR] += phi_width;
   else if (a_xiSrc[PHI_DIR] > phi_width)
      a_xi_valid[PHI_DIR] -= phi_width;
   if (a_xiSrc[Z_DIR] < 0.)
      a_xi_valid[Z_DIR] += z_width;
   else if (a_xiSrc[Z_DIR] > z_width)
      a_xi_valid[Z_DIR] -= z_width;
#else //CFG_DIM==2
   if (a_xiSrc[Z_DIR] < 0.)
      a_xi_valid[Z_DIR] += z_width;
   else if (a_xiSrc[Z_DIR] > z_width)
      a_xi_valid[Z_DIR] -= z_width;
#endif

}








#include "NamespaceFooter.H"
