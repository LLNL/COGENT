#include "Directions.H"
#include "LogRectPhaseCoordSys.H"
#include "BlockBoundary.H"

#include "NamespaceHeader.H"



LogRectPhaseCoordSys::LogRectPhaseCoordSys(ParmParse&                                  a_pp,
                                           const RefCountedPtr<CFG::LogRectCoordSys>&  a_mag_coords,
                                           const RefCountedPtr<VEL::VelCoordSys>&      a_vel_coords,
                                           const Vector<ProblemDomain>&                a_domains )
  : PhaseCoordSys(a_mag_coords, a_vel_coords, a_domains)
{

#if CFG_DIM == 3
  m_mb_dir = TOROIDAL_DIR;
#else
  //this options has not been fully developed yet
  m_mb_dir = POLOIDAL_DIR;
#endif

  m_mag_geom_type = a_mag_coords->getMagBlockGeomType();
  
  defineBoundaries();

  initializeBlockTransformations();

  getDecompositionParams( a_pp );
}



void
LogRectPhaseCoordSys::defineBoundaries()
{
  CH_assert(gotMappingBlocks());
  m_boundaries.resize(numBlocks());

  int num_blocks = numBlocks();
  
  for (int block(0); block<num_blocks; block++) {

    const CFG::MagBlockCoordSys* mag_block_coord_sys
    = (const CFG::MagBlockCoordSys*)m_mag_coords->getCoordSys(block);
    
    const CFG::ProblemDomain& domain( mag_block_coord_sys->domain() );
    
    Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries( m_boundaries[block] );
    for (int idir(RADIAL_DIR); idir<CFG_DIM; idir++) {
      if ( domain.isPeriodic(idir) ) {
        IntVect shift( m_mappingBlocks[block].size(idir) * BASISV(idir) );
        IndicesTransformation it;
        it.defineFromTranslation( shift );
        blockBoundaries[idir].define( it, block );
        it.defineFromTranslation( -shift );
        blockBoundaries[idir+SpaceDim].define( it, block );
      }
     
      else if ((m_mag_geom_type == "toroidal" || num_blocks > 1) && idir == m_mb_dir) {
        IntVect shift;
        IndicesTransformation it;
        
        if ( num_blocks == 1 ) {
          shift = m_mappingBlocks[block].size(m_mb_dir) * BASISV(m_mb_dir);
          it.defineFromTranslation( shift );
          blockBoundaries[m_mb_dir].define( it, block );
          it.defineFromTranslation( -shift );
          blockBoundaries[m_mb_dir+SpaceDim].define( it, block );
        }
        
        else {
          // Lower face coupling
          if ( block > 0 ) {  // Couple to the block on the lower boundary
            shift = -BLOCK_SEPARATION * BASISV(m_mb_dir);
            it.defineFromTranslation( shift );
            blockBoundaries[m_mb_dir].define( it, block - 1);
          }
          
          else {  // First toroidal block: couple it to the last one
            shift = (num_blocks * (BLOCK_SEPARATION + m_mappingBlocks[block].size(m_mb_dir))
                     - BLOCK_SEPARATION) * BASISV(m_mb_dir);
            it.defineFromTranslation( shift );
            blockBoundaries[m_mb_dir].define( it, num_blocks - 1);
          }
          
          // Upper face coupling
          if ( block < num_blocks-1 ) {  // Couple to the block on the upper toroidal boundary
            shift = BLOCK_SEPARATION * BASISV(m_mb_dir);
            it.defineFromTranslation( shift );
            blockBoundaries[m_mb_dir+SpaceDim].define( it, block + 1);
          }
          
          else {  // Last toroidal block: couple it to the first one
            shift = -(num_blocks * (BLOCK_SEPARATION + m_mappingBlocks[block].size(m_mb_dir))
                      - BLOCK_SEPARATION) * BASISV(m_mb_dir);
            it.defineFromTranslation( shift );
            blockBoundaries[m_mb_dir+SpaceDim].define( it, 0);
          }
        }
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
LogRectPhaseCoordSys::blockRemapping(RealVect&       a_xi_valid,
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
LogRectPhaseCoordSys::getDecompositionParams( ParmParse& a_pp )
{
  if (a_pp.contains("phase_decomp")) {
    std::vector<int> decomp( SpaceDim );
    a_pp.queryarr( "phase_decomp", decomp, 0, SpaceDim );

    for (int dir=0; dir<SpaceDim; ++dir) {
      m_decomposition[dir] = decomp[dir];
    }
     
     //Assume that the velocity space is a single block
     //that is numBlocks() = configuration num blocks

     m_decomposition[m_mb_dir] /= numBlocks();
  }
  else {
    m_decomposition = IntVect::Unit;
  }
}



#include "NamespaceFooter.H"
