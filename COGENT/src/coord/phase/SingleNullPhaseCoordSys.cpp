#include "Directions.H"
#include "SingleNullPhaseCoordSys.H"
#include "BlockBoundary.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "SingleNullBlockCoordSys.H"
#include "SingleNullBlockCoordSysModel.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM


#include "NamespaceHeader.H"



SingleNullPhaseCoordSys::SingleNullPhaseCoordSys( ParmParse&                                     a_pp,
                                                  const RefCountedPtr<CFG::SingleNullCoordSys>&  a_mag_coords,
                                                  const RefCountedPtr<VEL::VelCoordSys>&         a_vel_coords,
                                                  const Vector<ProblemDomain>&                   a_domains)
   : PhaseCoordSys(a_mag_coords, a_vel_coords, a_domains)
{
  defineBoundaries();

  initializeBlockTransformations();

  getDecompositionParams( a_pp );
}



void
SingleNullPhaseCoordSys::defineBoundaries()
{
   if ( ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->numPoloidalBlocks() == 8 ) {
      defineBoundaries8();
   }
   else {
      defineBoundaries10();
   }
}



void
SingleNullPhaseCoordSys::defineBoundaries8()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   int bc_tag = 0;  // define this later

   int num_poloidal_blocks = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->numPoloidalBlocks();
   int num_toroidal_sectors = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->numToroidalSectors();

   for ( int block_number=0, toroidal_sector=0; toroidal_sector<num_toroidal_sectors; ++toroidal_sector ) {
      for ( int poloidal_block=0; poloidal_block<num_poloidal_blocks; ++poloidal_block, ++block_number ) {

         IndicesTransformation it;
         IntVect shift;
         Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block_number];
         const CFG::SingleNullBlockCoordSys* mag_block_coord_sys
            = (const CFG::SingleNullBlockCoordSys*)m_mag_coords->getCoordSys(block_number);
         int block_type = mag_block_coord_sys->poloidalBlock();

         int toroidal_sector = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->toroidalBlockNumber(block_number);
         int LCORE = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::LCORE, toroidal_sector);
         int RCORE = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::RCORE, toroidal_sector);
         int LCSOL = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::LCSOL, toroidal_sector);
         int RCSOL = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::RCSOL, toroidal_sector);
         int LSOL  = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::LSOL,  toroidal_sector);
         int RSOL  = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::RSOL,  toroidal_sector);
         int LPF   = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::LPF,   toroidal_sector);
         int RPF   = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::RPF,   toroidal_sector);
         
         if( block_type == CFG::SingleNullBlockCoordSys::LCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            shift = POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LCSOL);
            
            if ( mag_block_coord_sys->truncated() ) {
               blockBoundaries[POLOIDAL_DIR].define(bc_tag);
            }
            else {
               shift = -POLOIDAL_BLOCK_SEP *
                  (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                   m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               it.defineFromTranslation(shift);
               blockBoundaries[POLOIDAL_DIR].define(it,RCORE);
            }
            shift = -(POLOIDAL_BLOCK_SEP+1) *
                 (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCORE);
         }
         else if( block_type == CFG::SingleNullBlockCoordSys::RCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            shift = POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RCSOL);

            shift = (POLOIDAL_BLOCK_SEP+1) *
                 (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,LCORE);

            if ( mag_block_coord_sys->truncated() ) {
               blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
            }
            else {
               shift = POLOIDAL_BLOCK_SEP *
                  (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                   m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               it.defineFromTranslation(shift);
               blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCORE);
            }
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::LCSOL ) {

            shift = - POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,LCORE);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            if ( mag_block_coord_sys->truncated() ) {
               blockBoundaries[POLOIDAL_DIR].define(bc_tag);
            }
            else {
               shift = -POLOIDAL_BLOCK_SEP *
                  (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                   m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               it.defineFromTranslation(shift);
               blockBoundaries[POLOIDAL_DIR].define(it,RCSOL);
            }
            shift = POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LSOL);
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::RCSOL ) {

            shift = - POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,RCORE);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RSOL);

            if ( mag_block_coord_sys->truncated() ) {
               blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
            }
            else {
               shift = POLOIDAL_BLOCK_SEP *
                  (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                   m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               it.defineFromTranslation(shift);
               blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCSOL);
            }
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::LSOL ) {

            shift = - POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,LPF);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,LCSOL);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::RSOL ) {

            shift = - POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,RPF);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            blockBoundaries[POLOIDAL_DIR].define(bc_tag);
            shift = POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCSOL);
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::LPF ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            shift = POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LSOL);

            shift = -(2*POLOIDAL_BLOCK_SEP+1) *
                 (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RPF);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::RPF ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            shift = POLOIDAL_BLOCK_SEP *
               m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RSOL);

            blockBoundaries[POLOIDAL_DIR].define(bc_tag);
            shift = (2*POLOIDAL_BLOCK_SEP+1) *
                 (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LPF);
         }
         else {
            MayDay::Error("SingleNullPhaseCoordSys::defineBoundaries(): case not implemented");
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
            // Lower face coupling
            if ( toroidal_sector > 0 ) {  // Couple to the block on the lower boundary
               shift = -TOROIDAL_BLOCK_SEP * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR].define( it, block_number - num_poloidal_blocks);
            }
            else {  // First toroidal block: couple it to the last one
               shift = (num_toroidal_sectors * (TOROIDAL_BLOCK_SEP + m_mappingBlocks[block_number].size(TOROIDAL_DIR)) - TOROIDAL_BLOCK_SEP) * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR].define( it, block_number + num_poloidal_blocks * (num_toroidal_sectors-1));
            }

            // Upper face coupling
            if ( toroidal_sector < num_toroidal_sectors-1 ) {  // Couple to the block on the upper toroidal boundary
               shift = TOROIDAL_BLOCK_SEP * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, block_number + num_poloidal_blocks);
            }
            else {  // Last toroidal block: couple it to the first one
               shift = -(num_toroidal_sectors * (TOROIDAL_BLOCK_SEP + m_mappingBlocks[block_number].size(TOROIDAL_DIR)) - TOROIDAL_BLOCK_SEP) * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, block_number - num_poloidal_blocks * (num_toroidal_sectors-1));
            }
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
SingleNullPhaseCoordSys::defineBoundaries10()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   int numcells_mcore_poloidal = m_mappingBlocks[CFG::SingleNullBlockCoordSys::MCORE].size(POLOIDAL_DIR);
   int numcells_lcore_poloidal = m_mappingBlocks[CFG::SingleNullBlockCoordSys::LCORE].size(POLOIDAL_DIR);
   int numcells_rcore_poloidal = m_mappingBlocks[CFG::SingleNullBlockCoordSys::RCORE].size(POLOIDAL_DIR);

   int L0 = m_mappingBlocks[CFG::SingleNullBlockCoordSys::MCORE].size(RADIAL_DIR);
   int lL1 = numcells_mcore_poloidal/2 + numcells_lcore_poloidal;
   int rL1 = numcells_mcore_poloidal/2 + numcells_rcore_poloidal;

   int bc_tag = 0;  // define this later

   int num_poloidal_blocks = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->numPoloidalBlocks();
   int num_toroidal_sectors = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->numToroidalSectors();

   for ( int block_number=0, toroidal_sector=0; toroidal_sector<num_toroidal_sectors; ++toroidal_sector ) {
      for ( int poloidal_block=0; poloidal_block<num_poloidal_blocks; ++poloidal_block, ++block_number ) {

         IndicesTransformation it;
         IntVect shift;
         Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block_number];
         const CFG::SingleNullBlockCoordSys* mag_block_coord_sys
            = (const CFG::SingleNullBlockCoordSys*)m_mag_coords->getCoordSys(block_number);
         int block_type = mag_block_coord_sys->poloidalBlock();

         int toroidal_sector = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->toroidalBlockNumber(block_number);
         int LCORE = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::LCORE, toroidal_sector);
         int RCORE = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::RCORE, toroidal_sector);
         int LCSOL = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::LCSOL, toroidal_sector);
         int RCSOL = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::RCSOL, toroidal_sector);
         int LSOL  = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::LSOL,  toroidal_sector);
         int RSOL  = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::RSOL,  toroidal_sector);
         int LPF   = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::LPF,   toroidal_sector);
         int RPF   = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::RPF,   toroidal_sector);
         int MCORE = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::MCORE, toroidal_sector);
         int MCSOL = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->blockNumber(CFG::SingleNullBlockCoordSys::MCSOL, toroidal_sector);
         
         if( block_type == CFG::SingleNullBlockCoordSys::MCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);

            shift =  (POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,MCSOL);

            shift = -(POLOIDAL_BLOCK_SEP * rL1 + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RCORE);

            shift =  (POLOIDAL_BLOCK_SEP * lL1 + numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCORE);
         }
         else if( block_type == CFG::SingleNullBlockCoordSys::LCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);

            shift = (POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LCSOL);

            shift = -(POLOIDAL_BLOCK_SEP * lL1 + numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,MCORE);

            shift = -((POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1)
                      + numcells_lcore_poloidal + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCORE);
         }
         else if( block_type == CFG::SingleNullBlockCoordSys::RCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);

            shift = (POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RCSOL);

            shift = ((POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1)
                     + numcells_lcore_poloidal + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,LCORE);

            shift = (POLOIDAL_BLOCK_SEP * rL1 + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,MCORE);
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::MCSOL ) {

            shift = -(POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,MCORE);

            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -(POLOIDAL_BLOCK_SEP * rL1 + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RCSOL);

            shift =  (POLOIDAL_BLOCK_SEP * lL1 + numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCSOL);
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::LCSOL ) {

            shift = -(POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,LCORE);

            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -(POLOIDAL_BLOCK_SEP * lL1 + numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,MCSOL);

            shift =  ((POLOIDAL_BLOCK_SEP+1) * lL1 - numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LSOL);
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::RCSOL ) {

            shift = -(POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,RCORE);

            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -((POLOIDAL_BLOCK_SEP+1) * rL1 - numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RSOL);

            shift = (POLOIDAL_BLOCK_SEP * rL1 + numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,MCSOL);
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::LPF ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);

            shift = (POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LSOL);

            shift = -(2*(POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1))*BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RPF);

            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::RPF ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);

            shift = (POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RSOL);

            blockBoundaries[POLOIDAL_DIR].define(bc_tag);

            shift =  (2*(POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1))*BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LPF);
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::LSOL ) {

            shift = -(POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,LPF);

            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -((POLOIDAL_BLOCK_SEP+1) * lL1 - numcells_lcore_poloidal ) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,LCSOL);

            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
         }
         else if ( block_type == CFG::SingleNullBlockCoordSys::RSOL ) {

            shift = -(POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,RPF);

            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            blockBoundaries[POLOIDAL_DIR].define(bc_tag);

            shift = ((POLOIDAL_BLOCK_SEP+1) * rL1 - numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCSOL);
         }
         else {
            MayDay::Error("SingleNullPhaseCoordSys::defineBoundaries(): case not implemented");
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
            // Lower face coupling
            if ( toroidal_sector > 0 ) {  // Couple to the block on the lower boundary
               shift = -TOROIDAL_BLOCK_SEP * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR].define( it, block_number - num_poloidal_blocks);
            }
            else {  // First toroidal block: couple it to the last one
               shift = (num_toroidal_sectors * (TOROIDAL_BLOCK_SEP + m_mappingBlocks[block_number].size(TOROIDAL_DIR)) - TOROIDAL_BLOCK_SEP) * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR].define( it, block_number + num_poloidal_blocks * (num_toroidal_sectors-1));
            }

            // Upper face coupling
            if ( toroidal_sector < num_toroidal_sectors-1 ) {  // Couple to the block on the upper toroidal boundary
               shift = TOROIDAL_BLOCK_SEP * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, block_number + num_poloidal_blocks);
            }
            else {  // Last toroidal block: couple it to the first one
               shift = -(num_toroidal_sectors * (TOROIDAL_BLOCK_SEP + m_mappingBlocks[block_number].size(TOROIDAL_DIR)) - TOROIDAL_BLOCK_SEP) * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, block_number - num_poloidal_blocks * (num_toroidal_sectors-1));
            }
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
SingleNullPhaseCoordSys::blockRemapping(RealVect&       a_xi_valid,
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
SingleNullPhaseCoordSys::blockRemappingGeneral( int              a_nDst,
                                                const RealVect&  a_xiSrc,
                                                int              a_nSrc ) const
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
SingleNullPhaseCoordSys::getDecomposition(int block) const
{
  IntVect decomp;

  switch ( poloidalBlock(block) )
    {
    case CFG::SingleNullBlockCoordSys::MCORE:
      decomp = m_decomp_mcore;
      break;
    case CFG::SingleNullBlockCoordSys::LCORE:
      decomp = m_decomp_lcore;
      break;
    case CFG::SingleNullBlockCoordSys::RCORE:
      decomp = m_decomp_rcore;
      break;
    case CFG::SingleNullBlockCoordSys::LPF:
      decomp = m_decomp_lpf;
      break;
    case CFG::SingleNullBlockCoordSys::RPF:
      decomp = m_decomp_rpf;
      break;
    case CFG::SingleNullBlockCoordSys::MCSOL:
      decomp = m_decomp_mcsol;
      break;
    case CFG::SingleNullBlockCoordSys::LCSOL:
      decomp = m_decomp_lcsol;
      break;
    case CFG::SingleNullBlockCoordSys::RCSOL:
      decomp = m_decomp_rcsol;
      break;
    case CFG::SingleNullBlockCoordSys::LSOL:
      decomp = m_decomp_lsol;
      break;
    case CFG::SingleNullBlockCoordSys::RSOL:
      decomp = m_decomp_rsol;
      break;
    default:
      MayDay::Error("SingleNullCoordSys::getDecomposition(): Invalid block number");
    }

  return decomp;
}



void
SingleNullPhaseCoordSys::getDecompositionParams( ParmParse& a_pp )
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

   std::vector<int> decomp_mcsol( SpaceDim );
   std::vector<int> decomp_lcsol( SpaceDim );
   std::vector<int> decomp_rcsol( SpaceDim );

   if ( a_pp.contains("csol.phase") ) {
      // For backwards compatibility with the old way of specifying
      // an 8-block single null geometry

      std::vector<int> decomp_csol( SpaceDim );
      a_pp.queryarr( "csol.phase", decomp_csol, 0, SpaceDim );

      for (int dir=0; dir<SpaceDim; ++dir) {
         decomp_mcsol[dir] = 0;
         decomp_lcsol[dir] = decomp_csol[dir];
         decomp_rcsol[dir] = decomp_csol[dir];
      }
      if ( decomp_csol[POLOIDAL_DIR]%2 != 0) {
         MayDay::Error("SingleNullPhaseCoordSys::getDecompositionParams(): Poloidal csol decomposition must be even");
      }
      else {
         decomp_lcsol[POLOIDAL_DIR] /= 2;
         decomp_rcsol[POLOIDAL_DIR] /= 2;
      }
   }
   else {
      a_pp.queryarr( "mcsol.phase", decomp_mcsol, 0, SpaceDim );
      a_pp.queryarr( "lcsol.phase", decomp_lcsol, 0, SpaceDim );
      a_pp.queryarr( "rcsol.phase", decomp_rcsol, 0, SpaceDim );
   }

   std::vector<int> decomp_lpf( SpaceDim );
   a_pp.queryarr( "lpf.phase", decomp_lpf, 0, SpaceDim );

   std::vector<int> decomp_rpf( SpaceDim );
   a_pp.queryarr( "rpf.phase", decomp_rpf, 0, SpaceDim );

   std::vector<int> decomp_lsol( SpaceDim );
   a_pp.queryarr( "lsol.phase", decomp_lsol, 0, SpaceDim );

   std::vector<int> decomp_rsol( SpaceDim );
   a_pp.queryarr( "rsol.phase", decomp_rsol, 0, SpaceDim );

   if (procID()==0) {
      cout << endl << "Single Null phase decomposition parameters" << endl << endl;

      cout << "decomp.mcore.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_mcore[i] << " ";
      cout << endl;

      cout << "decomp.lcore.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_lcore[i] << " ";
      cout << endl;

      cout << "decomp.rcore.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_rcore[i] << " ";
      cout << endl;

      cout << "decomp.lpf.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_lpf[i] << " ";
      cout << endl;

      cout << "decomp.rpf.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_rpf[i] << " ";
      cout << endl;

      cout << "decomp.mcsol.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_mcsol[i] << " ";
      cout << endl;

      cout << "decomp.lcsol.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_lcsol[i] << " ";
      cout << endl;

      cout << "decomp.rcsol.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_rcsol[i] << " ";
      cout << endl;

      cout << "decomp.lsol.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_lsol[i] << " ";
      cout << endl;

      cout << "decomp.rsol.phase = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_rsol[i] << " ";
      cout << endl;

      cout << endl;
   }

#if CFG_DIM==3
   int num_toroidal_sectors = ((RefCountedPtr<CFG::SingleNullCoordSys>)m_mag_coords)->numToroidalSectors();

   if ( decomp_mcore[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("Toroidal decomposition is not a multiple of the number of toroidal sectors in MCORE");
   }
   else {
      decomp_mcore[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
   if ( decomp_lcore[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("Toroidal decomposition is not a multiple of the number of toroidal sectors in LCORE");
   }
   else {
      decomp_lcore[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
   if ( decomp_rcore[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("Toroidal decomposition is not a multiple of the number of toroidal sectors in RCORE");
   }
   else {
      decomp_rcore[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
   if ( decomp_lpf[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("Toroidal decomposition is not a multiple of the number of toroidal sectors in LPF");
   }
   else {
      decomp_lpf[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
   if ( decomp_rpf[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("Toroidal decomposition is not a multiple of the number of toroidal sectors in RPF");
   }
   else {
      decomp_rpf[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
   if ( decomp_mcsol[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("Toroidal decomposition is not a multiple of the number of toroidal sectors in MCSOL");
   }
   else {
      decomp_mcsol[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
   if ( decomp_lcsol[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("Toroidal decomposition is not a multiple of the number of toroidal sectors in LCSOL");
   }
   else {
      decomp_lcsol[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
   if ( decomp_rcsol[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("Toroidal decomposition is not a multiple of the number of toroidal sectors in RCSOL");
   }
   else {
      decomp_rcsol[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
   if ( decomp_lsol[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("Toroidal decomposition is not a multiple of the number of toroidal sectors in LSOL");
   }
   else {
      decomp_lsol[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
   if ( decomp_rsol[TOROIDAL_DIR] % num_toroidal_sectors != 0 ) {
      MayDay::Error("Toroidal decomposition is not a multiple of the number of toroidal sectors in RSOL");
   }
   else {
      decomp_rsol[TOROIDAL_DIR] /= num_toroidal_sectors;
   }
#endif
   
   for (int dir=0; dir<SpaceDim; ++dir) {
      m_decomp_mcore[dir] = decomp_mcore[dir];
      m_decomp_lcore[dir] = decomp_lcore[dir];
      m_decomp_rcore[dir] = decomp_rcore[dir];
      m_decomp_lpf[dir] = decomp_lpf[dir];
      m_decomp_rpf[dir] = decomp_rpf[dir];
      m_decomp_mcsol[dir] = decomp_mcsol[dir];
      m_decomp_lcsol[dir] = decomp_lcsol[dir];
      m_decomp_rcsol[dir] = decomp_rcsol[dir];
      m_decomp_lsol[dir] = decomp_lsol[dir];
      m_decomp_rsol[dir] = decomp_rsol[dir];
   }
}


Vector<RealVect>
SingleNullPhaseCoordSys::displacements(const Vector<RealVect>&   a_dstCoords,
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
SingleNullPhaseCoordSys::containsPhysicalBoundary( int                    a_block_number,
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
         const CFG::SingleNullBlockCoordSysModel* mag_block_coord_sys
            = (const CFG::SingleNullBlockCoordSysModel*)m_mag_coords->getCoordSys(a_block_number);
         sheared_geom = mag_block_coord_sys->isFieldAlignedMapping();
      }
      else {
         const CFG::SingleNullBlockCoordSys* mag_block_coord_sys
            = (const CFG::SingleNullBlockCoordSys*)m_mag_coords->getCoordSys(a_block_number);
         sheared_geom = mag_block_coord_sys->isFieldAlignedMapping();
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
