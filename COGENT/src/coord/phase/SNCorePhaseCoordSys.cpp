#include "SNCorePhaseCoordSys.H"
#include "Directions.H"

#include "NamespaceHeader.H"



SNCorePhaseCoordSys::SNCorePhaseCoordSys( ParmParse&                                 a_pp,
                                          const RefCountedPtr<CFG::SNCoreCoordSys>&  a_mag_coords,
                                          const RefCountedPtr<VEL::VelCoordSys>&     a_vel_coords,
                                          const Vector<ProblemDomain>&               a_domains)
   : PhaseCoordSys(a_mag_coords, a_vel_coords, a_domains)
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

   int num_poloidal_blocks = ((RefCountedPtr<CFG::SNCoreCoordSys>)m_mag_coords)->numPoloidalBlocks();
   int num_toroidal_sectors = ((RefCountedPtr<CFG::SNCoreCoordSys>)m_mag_coords)->numToroidalSectors();

   for ( int block_number=0, toroidal_sector=0; toroidal_sector<num_toroidal_sectors; ++toroidal_sector ) {
      for ( int poloidal_block=0; poloidal_block<num_poloidal_blocks; ++poloidal_block, ++block_number ) {

         IndicesTransformation it;
         IntVect shift;
         Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block_number];
         const CFG::SNCoreBlockCoordSys* mag_block_coord_sys
            = (const CFG::SNCoreBlockCoordSys*)m_mag_coords->getCoordSys(block_number);
         int block_type = mag_block_coord_sys->poloidalBlock();

         if( block_type == CFG::SNCoreBlockCoordSys::LCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -POLOIDAL_BLOCK_SEP*(m_mappingBlocks[CFG::SNCoreBlockCoordSys::LCORE].size(POLOIDAL_DIR) +
                                         m_mappingBlocks[CFG::SNCoreBlockCoordSys::RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,CFG::SNCoreBlockCoordSys::RCORE);

            shift = -(POLOIDAL_BLOCK_SEP+1)*(m_mappingBlocks[CFG::SNCoreBlockCoordSys::LCORE].size(POLOIDAL_DIR)
                                             + m_mappingBlocks[CFG::SNCoreBlockCoordSys::RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,CFG::SNCoreBlockCoordSys::RCORE);
         }
         else if( block_type == CFG::SNCoreBlockCoordSys::RCORE ) {
            
            blockBoundaries[RADIAL_DIR].define(bc_tag);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = (POLOIDAL_BLOCK_SEP+1)*(m_mappingBlocks[CFG::SNCoreBlockCoordSys::LCORE].size(POLOIDAL_DIR)
                                            + m_mappingBlocks[CFG::SNCoreBlockCoordSys::RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,CFG::SNCoreBlockCoordSys::LCORE);

            shift = POLOIDAL_BLOCK_SEP*(m_mappingBlocks[CFG::SNCoreBlockCoordSys::LCORE].size(POLOIDAL_DIR) +
                                        m_mappingBlocks[CFG::SNCoreBlockCoordSys::RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,CFG::SNCoreBlockCoordSys::LCORE);
         }
         else {
            MayDay::Error("SNCorePhaseCoordSys::defineBoundaries(): case not implemented");
         }

#if CFG_DIM==3
         if ( num_toroidal_sectors == 1 ) {  // Assuming periodic coupling
            shift = m_mappingBlocks[block_number].size(TOROIDAL_DIR) * BASISV(TOROIDAL_DIR);
            it.defineFromTranslation( shift );
            blockBoundaries[TOROIDAL_DIR].define( it, block_number );
            it.defineFromTranslation( -shift );
            blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, block_number );
         }
         else {
            MayDay::Error("SNCorePhaseCoordSys::defineBoundaries2(): Toroidal block coupling is not defined");
         }
#endif

         // Block boundaries in the velocity directions are always physical
         blockBoundaries[VPARALLEL_DIR].define(bc_tag);
         blockBoundaries[VPARALLEL_DIR+SpaceDim].define(bc_tag);
         blockBoundaries[MU_DIR].define(bc_tag);
         blockBoundaries[MU_DIR+SpaceDim].define(bc_tag);
      }
   }

   m_gotBoundaries = true;
}


void
SNCorePhaseCoordSys::defineBoundaries3()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   int numcells_mcore_poloidal = m_mappingBlocks[CFG::SNCoreBlockCoordSys::MCORE].size(POLOIDAL_DIR);
   int numcells_lcore_poloidal = m_mappingBlocks[CFG::SNCoreBlockCoordSys::LCORE].size(POLOIDAL_DIR);
   int numcells_rcore_poloidal = m_mappingBlocks[CFG::SNCoreBlockCoordSys::RCORE].size(POLOIDAL_DIR);

   int lL1 = numcells_mcore_poloidal/2 + numcells_lcore_poloidal;
   int rL1 = numcells_mcore_poloidal/2 + numcells_rcore_poloidal;

   int bc_tag = 0;  // define this later

   int num_poloidal_blocks = ((RefCountedPtr<CFG::SNCoreCoordSys>)m_mag_coords)->numPoloidalBlocks();
   int num_toroidal_sectors = ((RefCountedPtr<CFG::SNCoreCoordSys>)m_mag_coords)->numToroidalSectors();

   for ( int block_number=0, toroidal_sector=0; toroidal_sector<num_toroidal_sectors; ++toroidal_sector ) {
      for ( int poloidal_block=0; poloidal_block<num_poloidal_blocks; ++poloidal_block, ++block_number ) {

         IndicesTransformation it;
         IntVect shift;
         Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block_number];
         const CFG::SingleNullBlockCoordSys* mag_block_coord_sys
            = (const CFG::SingleNullBlockCoordSys*)m_mag_coords->getCoordSys(block_number);
         int block_type = mag_block_coord_sys->poloidalBlock();

         if( block_type == CFG::SNCoreBlockCoordSys::MCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -(POLOIDAL_BLOCK_SEP * rL1 + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,CFG::SNCoreBlockCoordSys::RCORE);

            shift =  (POLOIDAL_BLOCK_SEP * lL1 + numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,CFG::SNCoreBlockCoordSys::LCORE);
         }
         else if( block_type == CFG::SNCoreBlockCoordSys::LCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -(POLOIDAL_BLOCK_SEP * lL1 + numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,CFG::SNCoreBlockCoordSys::MCORE);

            shift = -((POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1)
                      + numcells_lcore_poloidal + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,CFG::SNCoreBlockCoordSys::RCORE);
         }
         else if( block_type == CFG::SNCoreBlockCoordSys::RCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = ((POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1)
                     + numcells_lcore_poloidal + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,CFG::SNCoreBlockCoordSys::LCORE);

            shift = (POLOIDAL_BLOCK_SEP * rL1 + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,CFG::SNCoreBlockCoordSys::MCORE);
         }
         else {
            MayDay::Error("SNCorePhaseCoordSys::defineBoundaries(): case not implemented");
         }

#if CFG_DIM==3
         if ( num_toroidal_sectors == 1 ) {  // Assuming periodic coupling
            shift = m_mappingBlocks[block_number].size(TOROIDAL_DIR) * BASISV(TOROIDAL_DIR);
            it.defineFromTranslation( shift );
            blockBoundaries[TOROIDAL_DIR].define( it, block_number );
            it.defineFromTranslation( -shift );
            blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, block_number );
         }
         else {
            MayDay::Error("SNCorePhaseCoordSys::defineBoundaries3(): Toroidal block coupling is not defined");
         }
#endif

         // Block boundaries in the velocity directions are always physical
         blockBoundaries[VPARALLEL_DIR].define(bc_tag);
         blockBoundaries[VPARALLEL_DIR+SpaceDim].define(bc_tag);
         blockBoundaries[MU_DIR].define(bc_tag);
         blockBoundaries[MU_DIR+SpaceDim].define(bc_tag);
      }
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
   m_mag_coords->blockRemapping(xi_valid_config, n_valid_config, xiSrc_config, a_nSrc);

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

   CFG::RealVect xiDst_config = m_mag_coords->blockRemappingGeneral(a_nDst, xiSrc_config, a_nSrc);

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

   switch ( poloidalBlock(a_block) )
      {
      case CFG::SNCoreBlockCoordSys::MCORE:
         decomp = m_decomp_mcore;
         break;
      case CFG::SNCoreBlockCoordSys::LCORE:
         decomp = m_decomp_lcore;
         break;
      case CFG::SNCoreBlockCoordSys::RCORE:
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
         MayDay::Error("SNCorePhaseCoordSys::getDecompositionParams(): Poloidal core decomposition must be even");
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

#if CFG_DIM==3
   int num_toroidal_sectors = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->numToroidalSectors();

   if ( decomp_mcore[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("SNCorePhaseCoordSys::getDecompositionParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in MCORE");
   }
   else {
      decomp_mcore[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
   if ( decomp_lcore[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("SNCorePhaseCoordSys::getDecompositionParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in LCORE");
   }
   else {
      decomp_lcore[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
   if ( decomp_rcore[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("SNCorePhaseCoordSys::getDecompositionParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in RCORE");
   }
   else {
      decomp_rcore[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
#endif

   for (int dir=0; dir<SpaceDim; ++dir) {
      m_decomp_mcore[dir] = decomp_mcore[dir];
      m_decomp_lcore[dir] = decomp_lcore[dir];
      m_decomp_rcore[dir] = decomp_rcore[dir];
   }
}


Vector<RealVect>
SNCorePhaseCoordSys::displacements(const Vector<RealVect>&   a_dstCoords,
                                   const Vector<int>&        a_dstBlocks,
                                   const RealVect&           a_srcCoords,
                                   int                       a_srcBlock) const
{
   Vector<RealVect> disps =
      MultiBlockCoordSys::displacements(a_dstCoords, a_dstBlocks,
                                        a_srcCoords, a_srcBlock);

#if CFG_DIM==3

   Vector<CFG::RealVect> dstCoords_config(a_dstCoords.size());
   for (int i=0; i<a_dstCoords.size(); ++i) {
      for (int n=0; n<CFG_DIM; ++n) {
         dstCoords_config[i][n] = a_dstCoords[i][n];
      }
   }

   CFG::RealVect srcCoords_config;
   for (int n=0; n<CFG_DIM; ++n) {
      srcCoords_config[n] = a_srcCoords[n];
   }

   Vector<CFG::RealVect> dips_config = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->displacements(dstCoords_config,
                                                                                                             a_dstBlocks,
                                                                                                             srcCoords_config,
                                                                                                             a_srcBlock);
   for (int i=0; i<disps.size(); ++i) {
      for (int n=0; n<CFG_DIM; ++n) {
         disps[i][n] = dips_config[i][n];
      }
    }

#endif

   return disps;
}


bool
SNCorePhaseCoordSys::containsPhysicalBoundary( int                    a_block_number,
                                               int                    a_dir,
                                               const Side::LoHiSide&  a_side ) const
{
   bool contains_boundary = false;

#if CFG_DIM==2
   contains_boundary = PhaseCoordSys::containsPhysicalBoundary(a_block_number, a_dir, a_side);
#endif

#if CFG_DIM==3
   if ( a_dir == TOROIDAL_DIR ) {

      bool sheared_geom;
      if ( ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->isModelGeom() ) {
         const CFG::SNCoreBlockCoordSysModel* mag_coord_sys = (const CFG::SNCoreBlockCoordSysModel*)m_coordSysVect[a_block_number];
         sheared_geom = mag_coord_sys->isFieldAlignedMapping();
      }
      else {
         const CFG::SNCoreBlockCoordSys* mag_coord_sys = (const CFG::SNCoreBlockCoordSys*)m_coordSysVect[a_block_number];
         sheared_geom = mag_coord_sys->isFieldAlignedMapping();
      }

      if ( sheared_geom ) {
         // Set all toroidal block interfaces to physical boundaries;
         // This is done to deal with saw-tooth divertor BCs
         // Internal part of the block interface is handled by fillInternalGhosts
         contains_boundary = true;
      }
      else {
         contains_boundary = PhaseCoordSys::containsPhysicalBoundary(a_block_number, a_dir, a_side);
      }
   }
   else {
      contains_boundary = PhaseCoordSys::containsPhysicalBoundary(a_block_number, a_dir, a_side);
   }
#endif
   
   return contains_boundary;
}


#include "NamespaceFooter.H"
