#include "SNCorePhaseCoordSys.H"
#include "Directions.H"

#include "NamespaceHeader.H"



SNCorePhaseCoordSys::SNCorePhaseCoordSys( ParmParse&                   a_pp,
                                          const CFG::SNCoreCoordSys&   a_mag_coords,
                                          const VEL::VelCoordSys&      a_vel_coords,
                                          const Vector<ProblemDomain>& a_domains)
   : PhaseCoordSys(a_mag_coords, a_vel_coords, a_domains),
     m_config_coords( a_mag_coords )
{
   defineBoundaries();

   initializeBlockTransformations();

   getDecompositionParams( a_pp );
}



void
SNCorePhaseCoordSys::defineBoundaries()
{
   if ( numBlocks() == 2 ) {
      defineBoundaries2();
   }
   else {
      defineBoundaries3();
   }
}


void
SNCorePhaseCoordSys::defineBoundaries2()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   int bc_tag = 0;  // define this later

   for (int block = 0; block < numBlocks(); block++) {

      IndicesTransformation it;
      IntVect shift;
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block];
      const CFG::SNCoreBlockCoordSys* mag_block_coord_sys
         = (const CFG::SNCoreBlockCoordSys*)m_config_coords.getCoordSys(block);
      int block_type = mag_block_coord_sys->blockType();

      if( block_type == L_CORE ) {

        blockBoundaries[RADIAL_DIR].define(bc_tag);
        blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

        shift = -BLOCK_SEPARATION*(m_mappingBlocks[L_CORE].size(POLOIDAL_DIR) +
                                   m_mappingBlocks[R_CORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
        it.defineFromTranslation(shift);
        blockBoundaries[POLOIDAL_DIR].define(it,R_CORE);

        shift = -(BLOCK_SEPARATION+1)*(m_mappingBlocks[L_CORE].size(POLOIDAL_DIR)
                    + m_mappingBlocks[R_CORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
        it.defineFromTranslation(shift);
        blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,R_CORE);
      }
      else if( block_type == R_CORE ) {

        blockBoundaries[RADIAL_DIR].define(bc_tag);
        blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

        shift = (BLOCK_SEPARATION+1)*(m_mappingBlocks[L_CORE].size(POLOIDAL_DIR)
                   + m_mappingBlocks[R_CORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
        it.defineFromTranslation(shift);
        blockBoundaries[POLOIDAL_DIR].define(it,L_CORE);

        shift = BLOCK_SEPARATION*(m_mappingBlocks[L_CORE].size(POLOIDAL_DIR) +
                                  m_mappingBlocks[R_CORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
        it.defineFromTranslation(shift);
        blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,L_CORE);
      }
      else {
            MayDay::Error("SNCorePhaseCoordSys::defineBoundaries(): case not implemented");
      }

      // Block boundaries in the velocity directions are always physical
      blockBoundaries[VPARALLEL_DIR].define(bc_tag);
      blockBoundaries[VPARALLEL_DIR+SpaceDim].define(bc_tag);
      blockBoundaries[MU_DIR].define(bc_tag);
      blockBoundaries[MU_DIR+SpaceDim].define(bc_tag);
   }

   m_gotBoundaries = true;
}


void
SNCorePhaseCoordSys::defineBoundaries3()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   int numcells_mcore_poloidal = m_mappingBlocks[M_CORE].size(POLOIDAL_DIR);
   int numcells_lcore_poloidal = m_mappingBlocks[L_CORE].size(POLOIDAL_DIR);
   int numcells_rcore_poloidal = m_mappingBlocks[R_CORE].size(POLOIDAL_DIR);

   int lL1 = numcells_mcore_poloidal/2 + numcells_lcore_poloidal;
   int rL1 = numcells_mcore_poloidal/2 + numcells_rcore_poloidal;

   int bc_tag = 0;  // define this later

   for (int block = 0; block < numBlocks(); block++) {

      IndicesTransformation it;
      IntVect shift;
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block];
      const CFG::SingleNullBlockCoordSys* mag_block_coord_sys
         = (const CFG::SingleNullBlockCoordSys*)m_config_coords.getCoordSys(block);
      int block_type = mag_block_coord_sys->blockType();

      if( block_type == M_CORE ) {

         blockBoundaries[RADIAL_DIR].define(bc_tag);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         shift = -(BLOCK_SEPARATION * rL1 + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,R_CORE);

         shift =  (BLOCK_SEPARATION * lL1 + numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,L_CORE);
      }
      else if( block_type == L_CORE ) {

         blockBoundaries[RADIAL_DIR].define(bc_tag);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         shift = -(BLOCK_SEPARATION * lL1 + numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,M_CORE);

         shift = -((BLOCK_SEPARATION+1) * (lL1 + rL1)
                   + numcells_lcore_poloidal + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,R_CORE);
      }
      else if( block_type == R_CORE ) {

         blockBoundaries[RADIAL_DIR].define(bc_tag);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         shift = ((BLOCK_SEPARATION+1) * (lL1 + rL1)
                  + numcells_lcore_poloidal + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,L_CORE);

         shift = (BLOCK_SEPARATION * rL1 + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,M_CORE);
      }
      else {
            MayDay::Error("SNCorePhaseCoordSys::defineBoundaries(): case not implemented");
      }

      // Block boundaries in the velocity directions are always physical
      blockBoundaries[VPARALLEL_DIR].define(bc_tag);
      blockBoundaries[VPARALLEL_DIR+SpaceDim].define(bc_tag);
      blockBoundaries[MU_DIR].define(bc_tag);
      blockBoundaries[MU_DIR+SpaceDim].define(bc_tag);
   }

   m_gotBoundaries = true;
}



void
SNCorePhaseCoordSys::blockRemapping(RealVect&       a_xi_valid,
                                    int&            a_n_valid,
                                    const RealVect& a_xiSrc,
                                    int             a_nSrc) const
{
   // Get the valid coordinate and block number from the configuration coordinate system

   CFG::RealVect xiSrc_config;
   for (int dir=0; dir<CFG_DIM; ++dir) {
      xiSrc_config[dir] = a_xiSrc[dir];
   }

   CFG::RealVect xi_valid_config;
   int n_valid_config;
   m_config_coords.blockRemapping(xi_valid_config, n_valid_config, xiSrc_config, a_nSrc);

   a_n_valid = -1;

   if (n_valid_config >= 0) {

      // If the point is valid with respect to the configuration coordinates, make sure its
      // velocity components are also valid
      const PhaseBlockCoordSys& valid_cs = dynamic_cast<const PhaseBlockCoordSys&>(*m_coordSysVect[n_valid_config]);

      if ( a_xiSrc[VPARALLEL_DIR] >= valid_cs.lowerMappedCoordinate(VPARALLEL_DIR) &&
           a_xiSrc[VPARALLEL_DIR] <= valid_cs.upperMappedCoordinate(VPARALLEL_DIR) &&
           a_xiSrc[MU_DIR]        >= valid_cs.lowerMappedCoordinate(MU_DIR)        &&
           a_xiSrc[MU_DIR]        <= valid_cs.upperMappedCoordinate(MU_DIR)           ) {
         a_n_valid = n_valid_config;
         for (int dir=0; dir<CFG_DIM; ++dir) {
            a_xi_valid[dir] = xi_valid_config[dir];
         }
         for (int dir=CFG_DIM; dir<SpaceDim; ++dir) {
            a_xi_valid[dir] = a_xiSrc[dir-CFG_DIM];
         }
      }
   }

   // If the input point is invalid, just copy it to the outpout point
   if (a_n_valid < 0) {
      a_xi_valid = a_xiSrc;
   }
}



RealVect
SNCorePhaseCoordSys::blockRemappingGeneral(int             a_nDst,
                                           const RealVect& a_xiSrc,
                                           int             a_nSrc) const
{
   CFG::RealVect xiSrc_config;
   for (int dir=0; dir<CFG_DIM; ++dir) {
      xiSrc_config[dir] = a_xiSrc[dir];
   }

   CFG::RealVect xiDst_config = m_config_coords.blockRemappingGeneral(a_nDst, xiSrc_config, a_nSrc);

   RealVect xiDst;
   for (int dir=0; dir<CFG_DIM; ++dir) {
      xiDst[dir] = xiDst_config[dir];
   }
   xiDst[VPARALLEL_DIR] = a_xiSrc[VPARALLEL_DIR];
   xiDst[MU_DIR] = a_xiSrc[MU_DIR];

   return xiDst;
}



IntVect
SNCorePhaseCoordSys::getDecomposition( int a_block ) const
{
   IntVect decomp;

   switch (a_block)
      {
      case M_CORE:
         decomp = m_decomp_mcore;
         break;
      case L_CORE:
         decomp = m_decomp_lcore;
         break;
      case R_CORE:
         decomp = m_decomp_rcore;
         break;
      default:
         MayDay::Error("SingleNullCoordSys::getDecomposition(): Invalid block number");
      }

   return decomp;
}



void
SNCorePhaseCoordSys::getDecompositionParams( ParmParse& a_pp )
{
   std::vector<int> decomp_mcore( SpaceDim );
   std::vector<int> decomp_lcore( SpaceDim );
   std::vector<int> decomp_rcore( SpaceDim );

   if ( a_pp.contains("core.phase") ) {
      // For backwards compatibility with the old way of specifying
      // an 8-block single null geometry

      std::vector<int> decomp_core( SpaceDim );
      a_pp.queryarr( "core.phase", decomp_core, 0, SpaceDim );

      for (int dir=0; dir<SpaceDim; ++dir) {
         decomp_mcore[dir] = 0;
         decomp_lcore[dir] = decomp_core[dir];
         decomp_rcore[dir] = decomp_core[dir];
      }
      if ( decomp_core[POLOIDAL_DIR]%2 != 0) {
         MayDay::Error("SingleNullPhaseCoordSys::getDecompositionParams(): Poloidal core decomposition must be even");
      }
      else {
         decomp_lcore[POLOIDAL_DIR] /= 2;
         decomp_rcore[POLOIDAL_DIR] /= 2;
      }
   }
   else {
      a_pp.queryarr( "mcore.phase", decomp_mcore, 0, SpaceDim );
      a_pp.queryarr( "lcore.phase", decomp_lcore, 0, SpaceDim );
      a_pp.queryarr( "rcore.phase", decomp_rcore, 0, SpaceDim );
   }

   if (procID()==0) {
      cout << endl << "SNCore phase decomposition parameters" << endl << endl;

      cout << "decomp.mcore.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_mcore[i] << " ";
      cout << endl;

      cout << "decomp.lcore.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_lcore[i] << " ";
      cout << endl;

      cout << "decomp.rcore.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_rcore[i] << " ";
      cout << endl;
   }

   for (int dir=0; dir<SpaceDim; ++dir) {
      m_decomp_mcore[dir] = decomp_mcore[dir];
      m_decomp_lcore[dir] = decomp_lcore[dir];
      m_decomp_rcore[dir] = decomp_rcore[dir];
   }
}



#include "NamespaceFooter.H"
