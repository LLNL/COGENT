#include "Directions.H"
#include "SlabPhaseCoordSys.H"
#include "BlockBoundary.H"

#include "NamespaceHeader.H"



SlabPhaseCoordSys::SlabPhaseCoordSys( ParmParse&                    a_pp,
                                          const CFG::SlabCoordSys&    a_mag_coords,
                                          const VEL::VelCoordSys&       a_vel_coords,
                                          const Vector<ProblemDomain>&  a_domains )
  : PhaseCoordSys(a_mag_coords, a_vel_coords, a_domains)
{
  for (int dir(RADIAL_DIR); dir<SpaceDim; dir++) {
     if (dir<CFG_DIM) {
        m_is_periodic[dir] = a_mag_coords.isPeriodic( dir );
     }
     else {
        m_is_periodic[dir] = false;
     }
  }

  defineBoundaries();

  initializeBlockTransformations();

  getDecompositionParams( a_pp );
}



void
SlabPhaseCoordSys::defineBoundaries()
{
  CH_assert(gotMappingBlocks());
  m_boundaries.resize(numBlocks());

  for (int block(0); block<numBlocks(); block++) {

//     const PhaseBlockCoordSys& coord_sys( *(PhaseBlockCoordSys*)getCoordSys(block) );
//     const ProblemDomain& domain( coord_sys.domain() );

     Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries( m_boundaries[block] );
     for (int idir(RADIAL_DIR); idir<CFG_DIM; idir++) {
        if (m_is_periodic[idir]) {
           IntVect shift( m_mappingBlocks[block].size(idir) * BASISV(idir) );
           IndicesTransformation it;
           it.defineFromTranslation( shift );
           blockBoundaries[idir].define( it, block );
           it.defineFromTranslation( -shift );
           blockBoundaries[idir+SpaceDim].define( it, block );
        }
        else {
           blockBoundaries[idir].define( 0 );
           blockBoundaries[idir+SpaceDim].define( 0 );
        }
     }
     for (int idir(VPARALLEL_DIR); idir<SpaceDim; idir++) {
        blockBoundaries[idir].define( 0 );
        blockBoundaries[idir + SpaceDim].define( 0 );
     }
  }

  m_gotBoundaries = true;
}



void
SlabPhaseCoordSys::blockRemapping(RealVect&       a_xi_valid,
                                    int&            a_n_valid,
                                    const RealVect& a_xiSrc,
                                    int             a_nSrc) const
{
   const PhaseBlockCoordSys& block_coord_sys = dynamic_cast<const PhaseBlockCoordSys&>(*m_coordSysVect[a_nSrc]);

   // Check to see if the velocity components are within the block
   if ( a_xiSrc[VPARALLEL_DIR] >= block_coord_sys.lowerMappedCoordinate(VPARALLEL_DIR) &&
        a_xiSrc[VPARALLEL_DIR] <= block_coord_sys.upperMappedCoordinate(VPARALLEL_DIR) &&
        a_xiSrc[MU_DIR]        >= block_coord_sys.lowerMappedCoordinate(MU_DIR)        &&
        a_xiSrc[MU_DIR]        <= block_coord_sys.upperMappedCoordinate(MU_DIR)           ) {
         a_n_valid = a_nSrc;
         a_xi_valid = a_xiSrc;
   }
   else {
      a_n_valid = -1;
   }

   if (a_n_valid >= 0) {

      // All blocks are fully periodic in the configuration coordinates
      for (int dir=0; dir<CFG_DIM; ++dir) {
         double width = block_coord_sys.upperMappedCoordinate(dir)
            - block_coord_sys.lowerMappedCoordinate(dir);

         if (a_xiSrc[dir] < 0.) a_xi_valid[dir] += width;
         else if (a_xiSrc[dir] > width) a_xi_valid[dir] -= width;
      }
   }
}



void
SlabPhaseCoordSys::getDecompositionParams( ParmParse& a_pp )
{
  if (a_pp.contains("phase_decomp")) {
    std::vector<int> decomp( SpaceDim );
    a_pp.queryarr( "phase_decomp", decomp, 0, SpaceDim );

    for (int dir=0; dir<SpaceDim; ++dir) {
      m_decomposition[dir] = decomp[dir];
    }
  }
  else {
    m_decomposition = IntVect::Unit;
  }
}



#include "NamespaceFooter.H"
