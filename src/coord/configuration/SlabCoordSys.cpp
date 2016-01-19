#include "SlabCoordSys.H"
#include "SlabBlockCoordSys.H"
#include "FourthOrderUtil.H"
#include "BlockBoundary.H"
#include "SpreadingCopier.H"
#include "ReductionOps.H"
#include "Directions.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

bool SlabCoordSys::isPeriodic( const int a_dir ) const
{
   return m_block_zero_domain.isPeriodic( a_dir );
}

SlabCoordSys::SlabCoordSys( ParmParse&               a_pp,
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

  SlabBlockCoordSys* block_coords = new SlabBlockCoordSys( a_pp, m_block_zero_domain );

  Vector<MagBlockCoordSys *> coord_vec;
  coord_vec.push_back(block_coords);

  defineCoordSystemsAndBoundaries(coord_vec);

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



SlabCoordSys::~SlabCoordSys()
{
}



void
SlabCoordSys::defineBoundaries()
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
SlabCoordSys::blockRemapping(RealVect&       a_xi_valid,
                               int&            a_n_valid,
                               const RealVect& a_xiSrc,
                               int             a_nSrc) const
// Don't need to do anything in Slab
// { }
{
   //   double radial_width = m_xmax;
   //   double poloidal_width = m_ymax;
   double radial_width = 1.0;
   double poloidal_width = 1.0;
#if CFG_DIM == 3
   // NEW FOR 5D:
   //double toroidal_width = m_zmax;
   double toroidal_width = 1.0;
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

#include "NamespaceFooter.H"
