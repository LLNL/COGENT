#include "Directions.H"
#include "LogRectPhaseCoordSys.H"
#include "BlockBoundary.H"

#include "NamespaceHeader.H"



LogRectPhaseCoordSys::LogRectPhaseCoordSys(ParmParse&                    a_pp,
                                           const CFG::LogRectCoordSys&   a_mag_coords,
                                           const VEL::VelCoordSys&       a_vel_coords,
                                           const Vector<ProblemDomain>&  a_domains )
  : PhaseCoordSys(a_mag_coords, a_vel_coords, a_domains),
    m_config_coords( a_mag_coords )

{

#if CFG_DIM == 3
  m_mb_dir = TOROIDAL_DIR;
#else
  //this options has not been fully developed yet
  m_mb_dir = POLOIDAL_DIR;
#endif

  defineBoundaries();

  initializeBlockTransformations();

  getDecompositionParams( a_pp );
}



void
LogRectPhaseCoordSys::defineBoundaries()
{
  CH_assert(gotMappingBlocks());
  m_boundaries.resize(numBlocks());

  for (int block(0); block<numBlocks(); block++) {

    const CFG::MagBlockCoordSys* mag_block_coord_sys
    = (const CFG::MagBlockCoordSys*)m_config_coords.getCoordSys(block);
    
    const CFG::ProblemDomain& domain( mag_block_coord_sys->domain() );
    
    Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries( m_boundaries[block] );
    for (int idir(RADIAL_DIR); idir<CFG_DIM; idir++) {
      if ( domain.isPeriodic(idir) ) {
	//We have to put something at the periodic domain boundary.
	//The periodicity is used by exchange(), which is informed by the
	//isPeriodic member of the ProblemDomain object. Here, we basically
	//indicate that the boundary is not the one where physical BCs are
	//applied. The value of the "shift" variable in the code below
	//does not matter. It is however would be  important for the
	//multiblock exchange operations requiring the info about the connectivities.
	//NB: blockBoundaries[idir] = BlockBoundary::CONFORMAL did not work as
	// a workaround. M.A.D. 03/20/18
	IntVect shift( m_mappingBlocks[block].size(idir) * BASISV(idir) );
	IndicesTransformation it;
	it.defineFromTranslation( shift );
	blockBoundaries[idir].define( it, block );
	it.defineFromTranslation( -shift );
	blockBoundaries[idir+SpaceDim].define( it, block );
      }
      else if (idir == m_mb_dir && numBlocks()>1) {
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
