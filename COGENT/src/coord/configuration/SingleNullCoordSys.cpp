#include "SingleNullCoordSys.H"
#include "SingleNullBlockCoordSysModel.H"
#include "BlockBoundary.H"
#include "DataArray.H"
#include "Directions.H"
#include "CONSTANTS.H"
#include <float.h>

#include "NamespaceHeader.H"

const std::string SingleNullCoordSys::pp_name = "singlenull";

SingleNullCoordSys::SingleNullCoordSys( ParmParse& a_pp_grid,
                                        ParmParse& a_pp_geom )
   : m_model_geometry(false),
     m_num_toroidal_sectors(1)
{
   readParams(a_pp_grid, a_pp_geom);

   m_num_poloidal_blocks = m_original_eight_blocks? 8: SingleNullBlockCoordSys::NUM_POLOIDAL_BLOCKS;

   int num_blocks = m_num_poloidal_blocks * m_num_toroidal_sectors;
   
   Vector<Box> domain_boxes(num_blocks);
   for ( int block_number=0; block_number<num_blocks; ++block_number ) {
      domain_boxes[block_number] = Box(lo_mapped_index(block_number), hi_mapped_index(block_number));
   }

   /*
     Compute the computational mesh size.  We arbitrarily set the radial and poloidal
     width of the core computational domain to 1.  The toroidal width is specified in the 
     input as a fraction of two pi.  This, together with the input number of cells across the
     core radial, poloidal and toroidal widths sets the uniform cell size for all of the blocks.
     The computational domain widths of the blocks other than the core blocks are therefore
     determined by the specifed number of cells in those blocks.
   */
   int radial_width;
   int core_poloidal_width;
   if ( m_num_poloidal_blocks == 8 ) {
      radial_width = domain_boxes[SingleNullBlockCoordSys::LCORE].size(RADIAL_DIR);
      core_poloidal_width =  domain_boxes[SingleNullBlockCoordSys::LCORE].size(POLOIDAL_DIR)
         + domain_boxes[SingleNullBlockCoordSys::RCORE].size(POLOIDAL_DIR);
   }
   else {
      radial_width = domain_boxes[SingleNullBlockCoordSys::MCORE].size(RADIAL_DIR);
      core_poloidal_width = domain_boxes[SingleNullBlockCoordSys::MCORE].size(POLOIDAL_DIR)
         + domain_boxes[SingleNullBlockCoordSys::LCORE].size(POLOIDAL_DIR)
         + domain_boxes[SingleNullBlockCoordSys::RCORE].size(POLOIDAL_DIR);
   }

   RealVect dx;
   dx[RADIAL_DIR]   = 1./(double)radial_width;
   dx[POLOIDAL_DIR] = 1./(double)core_poloidal_width;

   bool is_periodic[SpaceDim];
   for (int dir=0; dir<SpaceDim; ++dir) {
      is_periodic[dir] = false;
   }

#if CFG_DIM==3
   dx[TOROIDAL_DIR] = m_toroidal_width / (double)m_numcells_toroidal;
#endif

   for ( int block_number=0, toroidal_sector=0; toroidal_sector<m_num_toroidal_sectors; ++toroidal_sector ) {
      for ( int poloidal_block=0; poloidal_block<m_num_poloidal_blocks; ++poloidal_block, ++block_number ) {
         if ( !m_model_geometry ) {
            SingleNullBlockCoordSys* geom
               = new SingleNullBlockCoordSys(a_pp_geom,
                                             ProblemDomain(domain_boxes[block_number], is_periodic),
                                             dx,
                                             poloidal_block,
                                             toroidal_sector,
                                             TOROIDAL_BLOCK_SEP);
            geom->readFiles(a_pp_geom);
            m_coord_vec.push_back(geom);
         }
         else {
            SingleNullBlockCoordSysModel* geom 
               = new SingleNullBlockCoordSysModel(a_pp_geom,
                                                  ProblemDomain(domain_boxes[block_number], is_periodic),
                                                  dx,
                                                  poloidal_block,
                                                  toroidal_sector,
                                                  TOROIDAL_BLOCK_SEP);
            m_coord_vec.push_back(geom);
         }
      }
   }

   defineCoordSystemsAndBoundaries(m_coord_vec);

   setXPointNeighborhood();

   // Define the boundary conditions for divergence cleaning (whether or not they're used)
   m_divergence_cleaning_bc.setPoloidalBlocksPerSector(numPoloidalBlocks());
   for (int block = 0; block < num_blocks; block++) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         for (int side=0; side<2; ++side) {
            if ( containsPhysicalBoundary(block, dir, (side==0? Side::LoHiSide::Lo: Side::LoHiSide::Hi)) ) {
               double bc_value = 0.;
               int bc_type = EllipticOpBC::DIRICHLET;     // Homogeneous Dirichlet
               int poloidal_block = poloidalBlockNumber(block);
               m_divergence_cleaning_bc.setBCType(poloidal_block, dir, side, bc_type);
               m_divergence_cleaning_bc.setBCValue(poloidal_block, dir, side, bc_value);
            }
         }
      }
   }

   //Checks if the grid cell sizes are continious across block boundaries (assumed for the model geoemtry)
   if ( m_model_geometry ) checkGridConsistency();

   // Check for sheared geometry (by checking block 0)
   if ( m_model_geometry ) {
      const SingleNullBlockCoordSysModel* mag_coord_sys = (const SingleNullBlockCoordSysModel*)m_coordSysVect[0];
      m_sheared_geometry = mag_coord_sys->isFieldAlignedMapping();
   }
   else {
      const SingleNullBlockCoordSys* mag_coord_sys = (const SingleNullBlockCoordSys*)m_coordSysVect[0];
      m_sheared_geometry = mag_coord_sys->isFieldAlignedMapping();
   }
}



SingleNullCoordSys::~SingleNullCoordSys()
{
   for (int i=0; i<m_coord_vec.size(); ++i) {
      delete m_coord_vec[i];
   }
}



void
SingleNullCoordSys::defineBoundaries()
{
   if ( m_original_eight_blocks ) {
      defineBoundaries8();
   }
   else {
      defineBoundaries10();
   }
}



void
SingleNullCoordSys::defineBoundaries8()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   int bc_tag = 0;

   for ( int block_number=0, toroidal_sector=0; toroidal_sector<m_num_toroidal_sectors; ++toroidal_sector ) {
      for ( int poloidal_block=0; poloidal_block<m_num_poloidal_blocks; ++poloidal_block, ++block_number ) {

         IndicesTransformation it;
         IntVect shift;
         Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block_number];
      
         bool poloidally_truncated;
      
         if ( m_model_geometry ) {
            const SingleNullBlockCoordSysModel* coord_sys = (const SingleNullBlockCoordSysModel*)m_coordSysVect[block_number];
            poloidally_truncated = coord_sys->truncated();
         }
         else {
            const SingleNullBlockCoordSys* coord_sys = (const SingleNullBlockCoordSys*)m_coordSysVect[block_number];
            poloidally_truncated = coord_sys->truncated();
         }
     
         int toroidal_block_number = toroidalBlockNumber(block_number);
         int LCORE = blockNumber(SingleNullBlockCoordSys::LCORE, toroidal_block_number);
         int RCORE = blockNumber(SingleNullBlockCoordSys::RCORE, toroidal_block_number);
         int LCSOL = blockNumber(SingleNullBlockCoordSys::LCSOL, toroidal_block_number);
         int RCSOL = blockNumber(SingleNullBlockCoordSys::RCSOL, toroidal_block_number);
         int LSOL  = blockNumber(SingleNullBlockCoordSys::LSOL,  toroidal_block_number);
         int RSOL  = blockNumber(SingleNullBlockCoordSys::RSOL,  toroidal_block_number);
         int LPF   = blockNumber(SingleNullBlockCoordSys::LPF,   toroidal_block_number);
         int RPF   = blockNumber(SingleNullBlockCoordSys::RPF,   toroidal_block_number);
         
         if( poloidal_block == SingleNullBlockCoordSys::LCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            shift = (POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LCSOL);

            if ( poloidally_truncated ) {
               blockBoundaries[POLOIDAL_DIR].define(bc_tag);
            }
            else {
               shift = -POLOIDAL_BLOCK_SEP*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                            m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               it.defineFromTranslation(shift);
               blockBoundaries[POLOIDAL_DIR].define(it,RCORE);
            }
            shift = -(POLOIDAL_BLOCK_SEP+1)*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                           + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCORE);
         }
         else if( poloidal_block == SingleNullBlockCoordSys::RCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            shift = (POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RCSOL);

            shift = (POLOIDAL_BLOCK_SEP+1)*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                          + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,LCORE);

            if ( poloidally_truncated ) {
               blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
            }
            else {
               shift = POLOIDAL_BLOCK_SEP*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                           m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               it.defineFromTranslation(shift);
               blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCORE);
            }
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::LCSOL ) {

            shift = - (POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,LCORE);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            if ( poloidally_truncated ) {
               blockBoundaries[POLOIDAL_DIR].define(bc_tag);
            }
            else {
               shift = -POLOIDAL_BLOCK_SEP*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                            m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               it.defineFromTranslation(shift);
               blockBoundaries[POLOIDAL_DIR].define(it,RCSOL);
            }
            shift = POLOIDAL_BLOCK_SEP*m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LSOL);
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::RCSOL ) {

            shift = - (POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,RCORE);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -POLOIDAL_BLOCK_SEP*m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RSOL);

            if ( poloidally_truncated ) {
               blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
            }
            else {
               shift = POLOIDAL_BLOCK_SEP*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                           m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               it.defineFromTranslation(shift);
               blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCSOL);
            }
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::LSOL ) {

            shift = - (POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,LPF);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -POLOIDAL_BLOCK_SEP*m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,LCSOL);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::RSOL ) {
            
            shift = - (POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,RPF);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            blockBoundaries[POLOIDAL_DIR].define(bc_tag);
            shift = POLOIDAL_BLOCK_SEP*m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCSOL);
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::LPF ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            shift = (POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LSOL);

            shift = -(2*POLOIDAL_BLOCK_SEP+1)*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                             + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RPF);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::RPF ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            shift = (POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RSOL);

            blockBoundaries[POLOIDAL_DIR].define(bc_tag);
            shift = (2*POLOIDAL_BLOCK_SEP+1)*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                            + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LPF);
         }
         else {
            MayDay::Error("SingleNullCoordSys::defineBoundaries(): case not implemented");
         }

#if CFG_DIM==3
         
         // Set toroidal block boundaries

         if ( m_num_toroidal_sectors == 1 ) {  // Assuming periodic coupling
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
               blockBoundaries[TOROIDAL_DIR].define( it, block_number - m_num_poloidal_blocks);
            }
            else {  // First toroidal block: couple it to the last one
               shift = (m_num_toroidal_sectors * (TOROIDAL_BLOCK_SEP + m_mappingBlocks[block_number].size(TOROIDAL_DIR))
                        - TOROIDAL_BLOCK_SEP) * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR].define( it, block_number + m_num_poloidal_blocks * (m_num_toroidal_sectors-1));
            }

            // Upper face coupling
            if ( toroidal_sector < m_num_toroidal_sectors-1 ) {  // Couple to the block on the upper toroidal boundary
               shift = TOROIDAL_BLOCK_SEP * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, block_number + m_num_poloidal_blocks);
            }
            else {  // Last toroidal block: couple it to the first one
               shift = -(m_num_toroidal_sectors * (TOROIDAL_BLOCK_SEP + m_mappingBlocks[block_number].size(TOROIDAL_DIR))
                         - TOROIDAL_BLOCK_SEP) * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, block_number - m_num_poloidal_blocks * (m_num_toroidal_sectors-1));
            }
         }
#endif
      }
   }
}



void
SingleNullCoordSys::defineBoundaries10()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   int L0 = m_numcells_core_radial;
   int lL1 = m_numcells_mcore_poloidal/2 + m_numcells_lcore_poloidal;
   int rL1 = m_numcells_mcore_poloidal/2 + m_numcells_rcore_poloidal;

   int bc_tag = 0;

   for ( int block_number=0, toroidal_sector=0; toroidal_sector<m_num_toroidal_sectors; ++toroidal_sector ) {
      for ( int poloidal_block=0; poloidal_block<m_num_poloidal_blocks; ++poloidal_block, ++block_number ) {

         IndicesTransformation it;
         IntVect shift;
         Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block_number];

         int toroidal_block_number = toroidalBlockNumber(block_number);
         int LCORE = blockNumber(SingleNullBlockCoordSys::LCORE, toroidal_block_number);
         int RCORE = blockNumber(SingleNullBlockCoordSys::RCORE, toroidal_block_number);
         int LCSOL = blockNumber(SingleNullBlockCoordSys::LCSOL, toroidal_block_number);
         int RCSOL = blockNumber(SingleNullBlockCoordSys::RCSOL, toroidal_block_number);
         int LSOL  = blockNumber(SingleNullBlockCoordSys::LSOL,  toroidal_block_number);
         int RSOL  = blockNumber(SingleNullBlockCoordSys::RSOL,  toroidal_block_number);
         int LPF   = blockNumber(SingleNullBlockCoordSys::LPF,   toroidal_block_number);
         int RPF   = blockNumber(SingleNullBlockCoordSys::RPF,   toroidal_block_number);
         int MCORE = blockNumber(SingleNullBlockCoordSys::MCORE, toroidal_block_number);
         int MCSOL = blockNumber(SingleNullBlockCoordSys::MCSOL, toroidal_block_number);
         
         if( poloidal_block == SingleNullBlockCoordSys::MCORE ) {
         
            blockBoundaries[RADIAL_DIR].define(bc_tag);

            shift =  (POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,MCSOL);

            shift = -(POLOIDAL_BLOCK_SEP * rL1 + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RCORE);

            shift =  (POLOIDAL_BLOCK_SEP * lL1 + m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCORE);
         }
         else if( poloidal_block == SingleNullBlockCoordSys::LCORE ) {
            
            blockBoundaries[RADIAL_DIR].define(bc_tag);

            shift = (POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LCSOL);

            shift = -(POLOIDAL_BLOCK_SEP * lL1 + m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,MCORE);

            shift = -((POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1)
                      + m_numcells_lcore_poloidal + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCORE);
         }
         else if( poloidal_block == SingleNullBlockCoordSys::RCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);

            shift = (POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RCSOL);

            shift = ((POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1)
                     + m_numcells_lcore_poloidal + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,LCORE);

            shift = (POLOIDAL_BLOCK_SEP * rL1 + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,MCORE);
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::MCSOL ) {

            shift = -(POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,MCORE);

            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -(POLOIDAL_BLOCK_SEP * rL1 + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RCSOL);

            shift =  (POLOIDAL_BLOCK_SEP * lL1 + m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCSOL);
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::LCSOL ) {

            shift = -(POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,LCORE);

            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -(POLOIDAL_BLOCK_SEP * lL1 + m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,MCSOL);

            shift =  ((POLOIDAL_BLOCK_SEP+1) * lL1 - m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LSOL);
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::RCSOL ) {

            shift = -(POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,RCORE);

            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -((POLOIDAL_BLOCK_SEP+1) * rL1 - m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RSOL);

            shift = (POLOIDAL_BLOCK_SEP * rL1 + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,MCSOL);
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::LPF ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);

            shift = (POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LSOL);

            shift = -(2*(POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1))*BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,RPF);

            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::RPF ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);

            shift = (POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RSOL);

            blockBoundaries[POLOIDAL_DIR].define(bc_tag);

            shift =  (2*(POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1))*BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LPF);
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::LSOL ) {
            
            shift = -(POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,LPF);

            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -((POLOIDAL_BLOCK_SEP+1) * lL1 - m_numcells_lcore_poloidal ) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,LCSOL);

            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
         }
         else if ( poloidal_block == SingleNullBlockCoordSys::RSOL ) {

            shift = -(POLOIDAL_BLOCK_SEP * L0) * BASISV(RADIAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[RADIAL_DIR].define(it,RPF);
            
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            blockBoundaries[POLOIDAL_DIR].define(bc_tag);

            shift = ((POLOIDAL_BLOCK_SEP+1) * rL1 - m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCSOL);
         }
         else {
            MayDay::Error("SingleNullCoordSys::defineBoundaries(): case not implemented");
         }

#if CFG_DIM==3

         // Set toroidal block boundaries

         if ( m_num_toroidal_sectors == 1 ) {  // Assuming periodic coupling
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
               blockBoundaries[TOROIDAL_DIR].define( it, block_number - m_num_poloidal_blocks);
            }
            else {  // First toroidal block: couple it to the last one
               shift = (m_num_toroidal_sectors * (TOROIDAL_BLOCK_SEP + m_mappingBlocks[block_number].size(TOROIDAL_DIR)) 
                        - TOROIDAL_BLOCK_SEP) * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR].define( it, block_number + m_num_poloidal_blocks * (m_num_toroidal_sectors-1));
            }

            // Upper face coupling
            if ( toroidal_sector < m_num_toroidal_sectors-1 ) {  // Couple to the block on the upper toroidal boundary
               shift = TOROIDAL_BLOCK_SEP * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, block_number + m_num_poloidal_blocks);
            }
            else {  // Last toroidal block: couple it to the first one
               shift = -(m_num_toroidal_sectors * (TOROIDAL_BLOCK_SEP + m_mappingBlocks[block_number].size(TOROIDAL_DIR)) 
                         - TOROIDAL_BLOCK_SEP) * BASISV(TOROIDAL_DIR);
               it.defineFromTranslation( shift );
               blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, block_number - m_num_poloidal_blocks * (m_num_toroidal_sectors-1));
            }
         }
#endif
      }
   }
}



void
SingleNullCoordSys::blockRemapping(RealVect&       a_xi_valid,
                                   int&            a_n_valid,
                                   const RealVect& a_xiSrc,
                                   int             a_nSrc) const
{
   if ( m_original_eight_blocks ) {
      blockRemapping8(a_xi_valid, a_n_valid, a_xiSrc, a_nSrc);
   }
   else {
      blockRemapping10(a_xi_valid, a_n_valid, a_xiSrc, a_nSrc);
   }
}


void
SingleNullCoordSys::blockRemapping8(RealVect&       a_xi_valid,
                                    int&            a_n_valid,
                                    const RealVect& a_xiSrc,
                                    int             a_nSrc) const
{
   const SingleNullBlockCoordSys* src_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(a_nSrc);

   if ( src_coord_sys->isValid(a_xiSrc, true) ) {
      a_n_valid = a_nSrc;
      a_xi_valid = a_xiSrc;
   }
   else {

      RealVect X = src_coord_sys->realCoord(a_xiSrc);
      POL::RealVect X_pol = src_coord_sys->restrictPhysCoordToPoloidal(X);

      POL::RealVect XminusXpt = X_pol - m_Xpoint;
      double distance_to_Xpt = XminusXpt.vectorLength();

      int n_valid_poloidal;

      if ( distance_to_Xpt < m_xpoint_radius ) {
         n_valid_poloidal = findBlockNearXpt(X_pol);
      }
      else {

         double lo_mapped_radial   = src_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
         double hi_mapped_radial   = src_coord_sys->upperMappedCoordinate(RADIAL_DIR);
         double lo_mapped_poloidal = src_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);
         double hi_mapped_poloidal = src_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);
         bool poloidally_truncated = src_coord_sys->truncated();

         switch ( poloidalBlockNumber(a_nSrc) )
            {
            case SingleNullBlockCoordSys::LCORE:
               if ( a_xiSrc[RADIAL_DIR] < lo_mapped_radial ||
                    (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal && poloidally_truncated) ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] > hi_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LSOL;
                  }
                  else if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCSOL;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCSOL;
                  }
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::RCORE;
               }
               break;
            case SingleNullBlockCoordSys::RCORE:
               if ( a_xiSrc[RADIAL_DIR] < lo_mapped_radial ||
                    (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal && poloidally_truncated) ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] > hi_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCSOL;
                  }
                  else if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RSOL;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCSOL;
                  }
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::LCORE;
               }
               break;
            case SingleNullBlockCoordSys::LCSOL:
               if ( a_xiSrc[RADIAL_DIR] > hi_mapped_radial ||
                    (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal && poloidally_truncated) ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] < lo_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LPF;
                  }
                  else if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCORE;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCORE;
                  }
               }
               else if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                  n_valid_poloidal = SingleNullBlockCoordSys::LSOL;
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::RCSOL;
               }
               break;
            case SingleNullBlockCoordSys::RCSOL:
               if ( a_xiSrc[RADIAL_DIR] > hi_mapped_radial ||
                    (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal && poloidally_truncated) ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] < lo_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCORE;
                  }
                  else if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RPF;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCORE;
                  }
               }
               else if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                  n_valid_poloidal = SingleNullBlockCoordSys::LCSOL;
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::RSOL;
               }
               break;
            case SingleNullBlockCoordSys::LSOL:
               if ( a_xiSrc[RADIAL_DIR] > hi_mapped_radial ||
                    a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] < lo_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCORE;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::LPF;
                  }
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::LCSOL;
               }
               break;
            case SingleNullBlockCoordSys::RSOL:
               if ( a_xiSrc[RADIAL_DIR] > hi_mapped_radial ||
                    a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] < lo_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCORE;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::RPF;
                  }
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::RCSOL;
               }
               break;
            case SingleNullBlockCoordSys::LPF:
               if ( a_xiSrc[RADIAL_DIR] < lo_mapped_radial ||
                    a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] > hi_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCSOL;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::LSOL;
                  }
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::RPF;
               }
               break;
            case SingleNullBlockCoordSys::RPF:
               if ( a_xiSrc[RADIAL_DIR] < lo_mapped_radial ||
                    a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] > hi_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCSOL;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::RSOL;
                  }
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::LPF;
               }
               break;
            default:
               MayDay::Error("SingleNullCoordSys::blockRemapping(): Unknown block number");
            }
      }

      if (n_valid_poloidal != -1) {

         a_n_valid = toroidalBlockNumber(a_nSrc) * numPoloidalBlocks() + n_valid_poloidal;

         const SingleNullBlockCoordSys* valid_cs = (SingleNullBlockCoordSys*)getCoordSys(a_n_valid);

         a_xi_valid = valid_cs->mappedCoord(X);

         // Check that the point is actually valid
         if ( !valid_cs->isValid(a_xi_valid, true) ) {

            if (distance_to_Xpt < m_xpoint_radius) {

               // If we've arrived here, it's because the separatrix isn't perfectly straight
               // in the X point neighborhood, as was assumed in findBlock().  We therefore
               // have to deal with the individual cases.

               int n_valid_poloidal_new = n_valid_poloidal;

               for (int poloidal_block=0; poloidal_block<numPoloidalBlocks(); ++poloidal_block) {
                  const SingleNullBlockCoordSys* cs = (SingleNullBlockCoordSys*)getCoordSys(poloidal_block);

                  if ( cs->containsPoloidalPoint(X_pol, false) ) {
                     n_valid_poloidal_new = poloidal_block;
                     break;
                  }
               }

#if 1
               if (n_valid_poloidal_new != n_valid_poloidal) {
                  
                  a_n_valid = toroidalBlockNumber(a_nSrc) * numPoloidalBlocks() + n_valid_poloidal_new;

                  const SingleNullBlockCoordSys* new_valid_cs = (SingleNullBlockCoordSys*)getCoordSys(a_n_valid);

                  a_xi_valid = new_valid_cs->mappedCoord(X);

                  if ( !new_valid_cs->isValid(a_xi_valid, true) ) {
#ifndef MODEL_GEOM
		    printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, new_valid_cs);
#endif
                  }
               }
#endif
            }
            else {
	            printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, valid_cs);
            }
         }
      }

#if CFG_DIM==3
#if 0
      // Not needed unless using MultiBlockLevelExchanges to fill toroidal extrablock ghosts
      while ( a_xi_valid[TOROIDAL_DIR] < 0. ) {
         a_xi_valid[TOROIDAL_DIR] += m_toroidal_width;
      }
      while ( a_xi_valid[TOROIDAL_DIR] > m_toroidal_width ) {
         a_xi_valid[TOROIDAL_DIR] -= m_toroidal_width;
      }
#endif
#endif

   }  // end of invalid point
}



void
SingleNullCoordSys::blockRemapping10( RealVect&       a_xi_valid,
                                      int&            a_n_valid,
                                      const RealVect& a_xiSrc,
                                      int             a_nSrc ) const
{
   const SingleNullBlockCoordSys* src_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(a_nSrc);

   if ( src_coord_sys->isValid(a_xiSrc, true) ) {
      a_n_valid = a_nSrc;
      a_xi_valid = a_xiSrc;
   }
   else {

      RealVect X = src_coord_sys->realCoord(a_xiSrc);
      POL::RealVect X_pol = src_coord_sys->restrictPhysCoordToPoloidal(X);

      POL::RealVect XminusXpt = X_pol - m_Xpoint;
      double distance_to_Xpt = XminusXpt.vectorLength();

      int n_valid_poloidal;

      if ( distance_to_Xpt < m_xpoint_radius ) {
         n_valid_poloidal = findBlockNearXpt(X_pol);
      }
      else {

         double lo_mapped_radial   = src_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
         double hi_mapped_radial   = src_coord_sys->upperMappedCoordinate(RADIAL_DIR);
         double lo_mapped_poloidal = src_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);
         double hi_mapped_poloidal = src_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);
         bool poloidally_truncated = src_coord_sys->truncated();

         switch ( poloidalBlockNumber(a_nSrc) )
            {
            case SingleNullBlockCoordSys::MCORE:
               if ( a_xiSrc[RADIAL_DIR] < lo_mapped_radial ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] > hi_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCSOL;
                  }
                  else if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCSOL;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::MCSOL;
                  }
               }
               else {
                  if ( a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal ) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCORE;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCORE;
                  }
               }
               break;
            case SingleNullBlockCoordSys::LCORE:
               if ( a_xiSrc[RADIAL_DIR] < lo_mapped_radial ||
                    (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal && poloidally_truncated) ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] > hi_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LSOL;
                  }
                  else if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::MCSOL;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCSOL;
                  }
               }
               else {
                  if ( a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal ) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCORE;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::MCORE;
                  }
               }
               break;
            case SingleNullBlockCoordSys::RCORE:
               if ( a_xiSrc[RADIAL_DIR] < lo_mapped_radial ||
                    (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal && poloidally_truncated) ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] > hi_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::MCSOL;
                  }
                  else if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RSOL;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCSOL;
                  }
               }
               else {
                  if ( a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal ) {
                     n_valid_poloidal = SingleNullBlockCoordSys::MCORE;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCORE;
                  }
               }
               break;
            case SingleNullBlockCoordSys::MCSOL:
               if ( a_xiSrc[RADIAL_DIR] > hi_mapped_radial ||
                    (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal && poloidally_truncated) ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] < lo_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCORE;
                  }
                  else if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCORE;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::MCORE;
                  }
               }
               else if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                  n_valid_poloidal = SingleNullBlockCoordSys::LCSOL;
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::RCSOL;
               }
               break;
            case SingleNullBlockCoordSys::LCSOL:
               if ( a_xiSrc[RADIAL_DIR] > hi_mapped_radial ||
                    (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal && poloidally_truncated) ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] < lo_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LPF;
                  }
                  else if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::MCORE;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCORE;
                  }
               }
               else if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                  n_valid_poloidal = SingleNullBlockCoordSys::LSOL;
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::MCSOL;
               }
               break;
            case SingleNullBlockCoordSys::RCSOL:
               if ( a_xiSrc[RADIAL_DIR] > hi_mapped_radial ||
                    (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal && poloidally_truncated) ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] < lo_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::MCORE;
                  }
                  else if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RPF;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCORE;
                  }
               }
               else if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                  n_valid_poloidal = SingleNullBlockCoordSys::MCSOL;
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::RSOL;
               }
               break;
            case SingleNullBlockCoordSys::LSOL:
               if ( a_xiSrc[RADIAL_DIR] > hi_mapped_radial ||
                    a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] < lo_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCORE;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::LPF;
                  }
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::LCSOL;
               }
               break;
            case SingleNullBlockCoordSys::RSOL:
               if ( a_xiSrc[RADIAL_DIR] > hi_mapped_radial ||
                    a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] < lo_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCORE;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::RPF;
                  }
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::RCSOL;
               }
               break;
            case SingleNullBlockCoordSys::LPF:
               if ( a_xiSrc[RADIAL_DIR] < lo_mapped_radial ||
                    a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] > hi_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::LCSOL;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::LSOL;
                  }
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::RPF;
               }
               break;
            case SingleNullBlockCoordSys::RPF:
               if ( a_xiSrc[RADIAL_DIR] < lo_mapped_radial ||
                    a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal ) {
                  n_valid_poloidal = -1;
               }
               else if (a_xiSrc[RADIAL_DIR] > hi_mapped_radial) {
                  if (a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal) {
                     n_valid_poloidal = SingleNullBlockCoordSys::RCSOL;
                  }
                  else {
                     n_valid_poloidal = SingleNullBlockCoordSys::RSOL;
                  }
               }
               else {
                  n_valid_poloidal = SingleNullBlockCoordSys::LPF;
               }
               break;
            default:
               MayDay::Error("SingleNullCoordSys::blockRemapping(): Unknown block number");
            }
      }

      if (n_valid_poloidal != -1) {

         a_n_valid = toroidalBlockNumber(a_nSrc) * numPoloidalBlocks() + n_valid_poloidal;

         const SingleNullBlockCoordSys* valid_cs = (SingleNullBlockCoordSys*)getCoordSys(a_n_valid);

         a_xi_valid = valid_cs->mappedCoord(X);

         // Check that the point is actually valid
         if ( !valid_cs->isValid(a_xi_valid, true) ) {

            if (distance_to_Xpt < m_xpoint_radius) {

               // If we've arrived here, it's because the separatrix isn't perfectly straight
               // in the X point neighborhood, as was assumed in findBlockNearXpt().  We therefore
               // have to deal with the special cases.

               int n_valid_poloidal_new = n_valid_poloidal;

               for (int poloidal_block=0; poloidal_block<numPoloidalBlocks(); ++poloidal_block) {
                  const SingleNullBlockCoordSys* cs = (SingleNullBlockCoordSys*)getCoordSys(poloidal_block);

                  if ( cs->containsPoloidalPoint(X_pol, false) ) {
                     n_valid_poloidal_new = poloidal_block;
                     break;
                  }
               }

#if 1
               if (n_valid_poloidal_new != n_valid_poloidal) {

                  a_n_valid = toroidalBlockNumber(a_nSrc) * numPoloidalBlocks() + n_valid_poloidal_new;

                  const SingleNullBlockCoordSys* new_valid_cs = (SingleNullBlockCoordSys*)getCoordSys(a_n_valid);

                  a_xi_valid = new_valid_cs->mappedCoord(X);

                  if ( !new_valid_cs->isValid(a_xi_valid, true) ) {
#ifndef MODEL_GEOM
		      printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, new_valid_cs);
#endif
                  }
               }
#endif
            }
            else {
#if 1
#ifndef MODEL_GEOM
                     printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, valid_cs);
#endif
#endif
            }
         }
      }

#if CFG_DIM==3
#if 0
      // Not needed unless using MultiBlockLevelExchanges to fill toroidal extrablock ghosts
      while ( a_xi_valid[TOROIDAL_DIR] < 0. ) {
         a_xi_valid[TOROIDAL_DIR] += m_toroidal_width;
      }
      while ( a_xi_valid[TOROIDAL_DIR] > m_toroidal_width ) {
         a_xi_valid[TOROIDAL_DIR] -= m_toroidal_width;
      }
#endif
#endif

   }
}



void
SingleNullCoordSys::printInvalidPointDiagnostics( const int                      a_nSrc,
                                                  const RealVect&                a_xiSrc,
                                                  const RealVect&                a_X,
                                                  const RealVect&                a_xi_valid,
                                                  const int                      a_n_valid,
                                                  const SingleNullBlockCoordSys* a_src_cs,
                                                  const SingleNullBlockCoordSys* a_valid_cs ) const
{
#if 0
   cout << "xi not valid, nSrc = " << a_nSrc << " a_xiSrc = " << a_xiSrc << " X = " << a_X
        << " xi_valid = " << a_xi_valid << " nvalid = " << a_n_valid << endl;
   cout << "src block: "
        << a_src_cs->lowerMappedCoordinate(RADIAL_DIR) << " "
        << a_src_cs->upperMappedCoordinate(RADIAL_DIR) << " "
#if CFG_DIM==3
        << a_src_cs->lowerMappedCoordinate(TOROIDAL_DIR) << " "
        << a_src_cs->upperMappedCoordinate(TOROIDAL_DIR) << " "
#endif
        << a_src_cs->lowerMappedCoordinate(POLOIDAL_DIR) << " "
        << a_src_cs->upperMappedCoordinate(POLOIDAL_DIR) << endl;
   cout << "valid block: "
        << a_valid_cs->lowerMappedCoordinate(RADIAL_DIR) << " "
        << a_valid_cs->upperMappedCoordinate(RADIAL_DIR) << " "
#if CFG_DIM==3
        << a_valid_cs->lowerMappedCoordinate(TOROIDAL_DIR) << " "
        << a_valid_cs->upperMappedCoordinate(TOROIDAL_DIR) << " "
#endif
        << a_valid_cs->lowerMappedCoordinate(POLOIDAL_DIR) << " "
        << a_valid_cs->upperMappedCoordinate(POLOIDAL_DIR) << endl;
#endif
}



void
SingleNullCoordSys::setXPointNeighborhood()
{
#if CFG_DIM==2

   {
      // Block boundary between LCORE and RCORE

      const SingleNullBlockCoordSys* lcore_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::LCORE);

      Box boundary_box = surroundingNodes(bdryHi(lcore_coord_sys->domain().domainBox(), POLOIDAL_DIR, 1), RADIAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      lcore_coord_sys->getNodeCenteredRealCoords(X);
      
      m_Xpoint = RealVect(X(hi,0), X(hi,1));

      // Initialize the ray to lie on the block boundary, pointing away from X point

      m_ray[0] = RealVect(X(lo,0), X(lo,1)) - m_Xpoint;

      // Set the radius of the X point neighborhood, which we currently take to as the
      // length of the block interface between the LCORE and RCORE.  If there is any
      // need to reduce the size of the neighborhood, here is where that should be done.  

      m_xpoint_radius = m_ray[0].vectorLength();

      // This loop isn't actually necessary unless the preceding m_xpoint_radius
      // definition has somehow been manually overridden to a smaller value
      for( IntVect iv(hi); iv[RADIAL_DIR] >= lo[RADIAL_DIR]; iv[RADIAL_DIR]--) {

         RealVect diff = RealVect(X(iv,0), X(iv,1)) - m_Xpoint;


         if (diff.vectorLength() > m_xpoint_radius) {
            break;
         }
         else {
            m_ray[0] = diff;
         }
      }
   }

   {
      // Block boundary between RCORE and RCSOL

      const SingleNullBlockCoordSys* rcore_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::RCORE);

      Box boundary_box = surroundingNodes(bdryHi(rcore_coord_sys->domain().domainBox(), RADIAL_DIR, 1), POLOIDAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      rcore_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point

      m_ray[1] = RealVect(X(hi,0), X(hi,1)) - m_Xpoint;

      for( IntVect iv(lo); iv[POLOIDAL_DIR] <= hi[POLOIDAL_DIR]; iv[POLOIDAL_DIR]++) {

         RealVect diff = RealVect(X(iv,0), X(iv,1)) - m_Xpoint;

         if (diff.vectorLength() > m_xpoint_radius) {
            break;
         }
         else {
            m_ray[1] = diff;
         }
      }
   }

   {
      // Block boundary between RCSOL and RSOL

      const SingleNullBlockCoordSys* rcsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::RCSOL);

      Box boundary_box = surroundingNodes(bdryLo(rcsol_coord_sys->domain().domainBox(), POLOIDAL_DIR, 1), RADIAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      rcsol_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point

      m_ray[2] = RealVect(X(hi,0), X(hi,1)) - m_Xpoint;

      for( IntVect iv(lo); iv[RADIAL_DIR] <= hi[RADIAL_DIR]; iv[RADIAL_DIR]++) {

         RealVect diff = RealVect(X(iv,0), X(iv,1)) - m_Xpoint;

         if (diff.vectorLength() > m_xpoint_radius) {
            break;
         }
         else {
            m_ray[2] = diff;
         }
      }
   }

   {
      // Block boundary between RSOL and RPF

      const SingleNullBlockCoordSys* rsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::RSOL);

      Box boundary_box = surroundingNodes(bdryLo(rsol_coord_sys->domain().domainBox(), RADIAL_DIR, 1), POLOIDAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      rsol_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point

      m_ray[3] = RealVect(X(lo,0), X(lo,1)) - m_Xpoint;

      for( IntVect iv(hi); iv[POLOIDAL_DIR] >= lo[POLOIDAL_DIR]; iv[POLOIDAL_DIR]--) {

         RealVect diff = RealVect(X(iv,0), X(iv,1)) - m_Xpoint;

         if (diff.vectorLength() > m_xpoint_radius) {
            break;
         }
         else {
            m_ray[3] = diff;
         }
      }
   }

   {
      // Block boundary between RPF and LPF

      const SingleNullBlockCoordSys* rpf_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::RPF);

      Box boundary_box = surroundingNodes(bdryHi(rpf_coord_sys->domain().domainBox(), POLOIDAL_DIR, 1), RADIAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      rpf_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point

      m_ray[4] = RealVect(X(lo,0), X(lo,1)) - m_Xpoint;

      for( IntVect iv(hi); iv[RADIAL_DIR] >= lo[RADIAL_DIR]; iv[RADIAL_DIR]--) {

         RealVect diff = RealVect(X(iv,0), X(iv,1)) - m_Xpoint;

         if (diff.vectorLength() > m_xpoint_radius) {
            break;
         }
         else {
            m_ray[4] = diff;
         }
      }
   }

   {
      // Block boundary between LPF and LSOL

      const SingleNullBlockCoordSys* lpf_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::LPF);

      Box boundary_box = surroundingNodes(bdryHi(lpf_coord_sys->domain().domainBox(), RADIAL_DIR, 1), POLOIDAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      lpf_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point

      m_ray[5] = RealVect(X(hi,0), X(hi,1)) - m_Xpoint;

      for( IntVect iv(lo); iv[POLOIDAL_DIR] <= hi[POLOIDAL_DIR]; iv[POLOIDAL_DIR]++) {

         RealVect diff = RealVect(X(iv,0), X(iv,1)) - m_Xpoint;

         if (diff.vectorLength() > m_xpoint_radius) {
            break;
         }
         else {
            m_ray[5] = diff;
         }
      }
   }

   {
      // Block boundary between LSOL and LCSOL

      const SingleNullBlockCoordSys* lsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::LSOL);

      Box boundary_box = surroundingNodes(bdryLo(lsol_coord_sys->domain().domainBox(), POLOIDAL_DIR, 1), RADIAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      lsol_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point

      m_ray[6] = RealVect(X(hi,0), X(hi,1)) - m_Xpoint;

      for( IntVect iv(lo); iv[RADIAL_DIR] <= hi[RADIAL_DIR]; iv[RADIAL_DIR]++) {

         RealVect diff = RealVect(X(iv,0), X(iv,1)) - m_Xpoint;

         if (diff.vectorLength() > m_xpoint_radius) {
            break;
         }
         else {
            m_ray[6] = diff;
         }
      }
   }

   {
      // Block boundary between LCSOL and LCORE

      const SingleNullBlockCoordSys* lcsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::LCSOL);

      Box boundary_box = surroundingNodes(bdryLo(lcsol_coord_sys->domain().domainBox(), RADIAL_DIR, 1), POLOIDAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      lcsol_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point

      m_ray[7] = RealVect(X(lo,0), X(lo,1)) - m_Xpoint;

      for( IntVect iv(hi); iv[POLOIDAL_DIR] >= lo[POLOIDAL_DIR]; iv[POLOIDAL_DIR]--) {

         RealVect diff = RealVect(X(iv,0), X(iv,1)) - m_Xpoint;

         if (diff.vectorLength() > m_xpoint_radius) {
            break;
         }
         else {
            m_ray[7] = diff;
         }
      }
   }

   for (int ray=0; ray<8; ++ray) {
      m_ray[ray] /= m_ray[ray].vectorLength();
   }

   if (procID()==0) {
      cout << "X point located at " << m_Xpoint << ", neighborhood radius = " << m_xpoint_radius << endl;
   }

   m_ray_angle[0] = 0.;
   for (int block=1; block<8; ++block) {
      m_ray_angle[block] = m_ray_angle[block-1] + acos(m_ray[block].dotProduct(m_ray[block-1]));
   }
   m_ray_angle[8] = 2.*Pi;

#endif

#if CFG_DIM==3

   double cut_frac = 1.;
   double cut_frac2 = 0.15;

   const SingleNullBlockCoordSys* lcore_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::LCORE);

   RealVect Xpoint_mapped;
   for (int n=0; n<SpaceDim; ++n) {
      Xpoint_mapped[n] = lcore_coord_sys->upperMappedCoordinate(n);
   }

   m_Xpoint = lcore_coord_sys->restrictPhysCoordToPoloidal(lcore_coord_sys->realCoord(Xpoint_mapped));

   double rad_mapped;
   double tor_mapped;
   double pol_mapped;
   double cut_frac_tmp;

   cut_frac_tmp = ( lcore_coord_sys->truncated() || numPoloidalBlocks() == 10 )? cut_frac: cut_frac2;

   rad_mapped = (1.-cut_frac_tmp)*lcore_coord_sys->upperMappedCoordinate(RADIAL_DIR)
                   + cut_frac_tmp*lcore_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
   tor_mapped = lcore_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR);
   pol_mapped = lcore_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);

   RealVect mapped_coord = RealVect(rad_mapped, tor_mapped, pol_mapped);

   m_ray[0] = lcore_coord_sys->restrictPhysCoordToPoloidal(lcore_coord_sys->realCoord(mapped_coord)) - m_Xpoint;

   const SingleNullBlockCoordSys* rcore_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::RCORE);

   cut_frac_tmp = ( rcore_coord_sys->truncated() || numPoloidalBlocks() == 10 )? cut_frac: cut_frac2;

   rad_mapped = rcore_coord_sys->upperMappedCoordinate(RADIAL_DIR);
   tor_mapped = rcore_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR);
   pol_mapped = (1.-cut_frac_tmp)*rcore_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR)
                   + cut_frac_tmp*rcore_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);

   mapped_coord = RealVect(rad_mapped, tor_mapped, pol_mapped);

   m_ray[1] = rcore_coord_sys->restrictPhysCoordToPoloidal(rcore_coord_sys->realCoord(mapped_coord)) - m_Xpoint;

   const SingleNullBlockCoordSys* rcsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::RCSOL);

   cut_frac_tmp = ( rcsol_coord_sys->truncated() || numPoloidalBlocks() == 10 )? cut_frac: cut_frac2;

   rad_mapped = (1.-cut_frac_tmp)*rcsol_coord_sys->lowerMappedCoordinate(RADIAL_DIR)
                   + cut_frac_tmp*rcsol_coord_sys->upperMappedCoordinate(RADIAL_DIR);
   tor_mapped = rcsol_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR);
   pol_mapped = rcsol_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);

   mapped_coord = RealVect(rad_mapped, tor_mapped, pol_mapped);

   m_ray[2] = rcsol_coord_sys->restrictPhysCoordToPoloidal(rcsol_coord_sys->realCoord(mapped_coord)) - m_Xpoint;

   const SingleNullBlockCoordSys* rsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::RSOL);

   rad_mapped = rsol_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
   tor_mapped = rsol_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR);
   pol_mapped = (1.-cut_frac)*rsol_coord_sys->upperMappedCoordinate(POLOIDAL_DIR)
                   + cut_frac*rsol_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);

   mapped_coord = RealVect(rad_mapped, tor_mapped, pol_mapped);

   m_ray[3] = rsol_coord_sys->restrictPhysCoordToPoloidal(rsol_coord_sys->realCoord(mapped_coord)) - m_Xpoint;

   const SingleNullBlockCoordSys* rpf_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::RPF);

   rad_mapped = (1.-cut_frac)*rpf_coord_sys->upperMappedCoordinate(RADIAL_DIR)
                   + cut_frac*rpf_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
   tor_mapped = rpf_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR);
   pol_mapped = rpf_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);

   mapped_coord = RealVect(rad_mapped, tor_mapped, pol_mapped);

   m_ray[4] = rpf_coord_sys->restrictPhysCoordToPoloidal(rpf_coord_sys->realCoord(mapped_coord)) - m_Xpoint;

   const SingleNullBlockCoordSys* lpf_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::LPF);

   rad_mapped = lpf_coord_sys->upperMappedCoordinate(RADIAL_DIR);
   tor_mapped = lpf_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR);
   pol_mapped = (1.-cut_frac)*lpf_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR)
                   + cut_frac*lpf_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);

   mapped_coord = RealVect(rad_mapped, tor_mapped, pol_mapped);

   m_ray[5] = lpf_coord_sys->restrictPhysCoordToPoloidal(lpf_coord_sys->realCoord(mapped_coord)) - m_Xpoint;

   const SingleNullBlockCoordSys* lsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::LSOL);

   rad_mapped = (1.-cut_frac)*lsol_coord_sys->lowerMappedCoordinate(RADIAL_DIR)
                   + cut_frac*lsol_coord_sys->upperMappedCoordinate(RADIAL_DIR);
   tor_mapped = lsol_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR);
   pol_mapped = lsol_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);

   mapped_coord = RealVect(rad_mapped, tor_mapped, pol_mapped);

   m_ray[6] = lsol_coord_sys->restrictPhysCoordToPoloidal(lsol_coord_sys->realCoord(mapped_coord)) - m_Xpoint;

   const SingleNullBlockCoordSys* lcsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(SingleNullBlockCoordSys::LCSOL);

   cut_frac_tmp = ( lcsol_coord_sys->truncated() || numPoloidalBlocks() == 10 )? cut_frac: cut_frac2;

   rad_mapped = lcsol_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
   tor_mapped = lcsol_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR);
   pol_mapped = (1.-cut_frac_tmp)*lcsol_coord_sys->upperMappedCoordinate(POLOIDAL_DIR)
                   + cut_frac_tmp*lcsol_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);

   mapped_coord = RealVect(rad_mapped, tor_mapped, pol_mapped);

   m_ray[7] = lcsol_coord_sys->restrictPhysCoordToPoloidal(lcsol_coord_sys->realCoord(mapped_coord)) - m_Xpoint;

   m_xpoint_radius = m_ray[0].vectorLength();
   m_ray[0] /= m_xpoint_radius;
   for (int block=1; block<8; ++block) {
      double ray_length = m_ray[block].vectorLength();
      if (ray_length < m_xpoint_radius) m_xpoint_radius = ray_length;
      m_ray[block] /= ray_length;
   }

   if (procID()==0) {
      cout << "X point neighborhood radius = " << m_xpoint_radius << endl;
   }

   m_ray_angle[0] = 0.;
   for (int block=1; block<8; ++block) {
      m_ray_angle[block] = m_ray_angle[block-1] + acos(m_ray[block].dotProduct(m_ray[block-1]));
   }
   m_ray_angle[8] = 2.*Pi;

#endif

   m_sector_block[0] = SingleNullBlockCoordSys::RCORE;
   m_sector_block[1] = SingleNullBlockCoordSys::RCSOL;
   m_sector_block[2] = SingleNullBlockCoordSys::RSOL;
   m_sector_block[3] = SingleNullBlockCoordSys::RPF;
   m_sector_block[4] = SingleNullBlockCoordSys::LPF;
   m_sector_block[5] = SingleNullBlockCoordSys::LSOL;
   m_sector_block[6] = SingleNullBlockCoordSys::LCSOL;
   m_sector_block[7] = SingleNullBlockCoordSys::LCORE;
}


int
SingleNullCoordSys::findBlockNearXpt( const POL::RealVect& a_X_pol ) const
{
   POL::RealVect local_X_pol = a_X_pol - m_Xpoint;

   POL::RealVect ray0_normal(m_ray[0][1],-m_ray[0][0]);  // m_ray[0] rotated 90 degrees clockwise

   double angle = acos( local_X_pol.dotProduct(m_ray[0]) / local_X_pol.vectorLength() );

   if (local_X_pol.dotProduct(ray0_normal) < 0.) {
      angle = 2.*Pi - angle;
   }

   int block = -1;

   if ( angle >= m_ray_angle[0] && angle <= m_ray_angle[1] ) {
      block = m_sector_block[0];
   }
   else {
      for (int i=1; i<8; ++i) {
         if ( angle > m_ray_angle[i] && angle <= m_ray_angle[i+1] ) {
            block = m_sector_block[i];
            break;
         }
      }
   }

   if (block == -1) {
      cout << "Looking for block of " << a_X_pol << endl;
      MayDay::Error("SingleNullCoordSys::findBlockNearXpt(): Couldn't find block");
   }

   return block;
}



void
SingleNullCoordSys::defineEllipticOpBC( EllipticOpBC& a_bc ) const
{
   SingleNullEllipticOpBC& bc = (SingleNullEllipticOpBC&)a_bc;

   // Radial core boundary
   RefCountedPtr<GridFunction> radial_core_function = bc.getCoreFunction();
   if ( !radial_core_function.isNull() ) {
      if ( typeid(*radial_core_function) == typeid(DataArray) ) {
         RefCountedPtr<DataArray> data_array = (RefCountedPtr<DataArray>)radial_core_function;
         CH_assert( !data_array.isNull() );

         const Box& lcore_domain_box = getCoordSys(SingleNullBlockCoordSys::LCORE)->domain().domainBox();
         const Box& rcore_domain_box = getCoordSys(SingleNullBlockCoordSys::RCORE)->domain().domainBox();

         Box fab_box = minBox( bdryLo(lcore_domain_box, RADIAL_DIR, 1),
                               bdryLo(rcore_domain_box, RADIAL_DIR, 1) );

         FArrayBox fab(fab_box, 1);
         data_array->setData(fab, false);
      }
   }

   // Radial scrape-off boundary
   RefCountedPtr<GridFunction> radial_sol_function = bc.getSOLFunction();
   if ( !radial_sol_function.isNull() ) {
      if ( typeid(*radial_sol_function) == typeid(DataArray) ) {
         RefCountedPtr<DataArray> data_array = (RefCountedPtr<DataArray>)radial_sol_function;
         CH_assert( !data_array.isNull() );

         const Box& lcsol_domain_box = getCoordSys(SingleNullBlockCoordSys::LCSOL)->domain().domainBox();
         const Box& rcsol_domain_box = getCoordSys(SingleNullBlockCoordSys::RCSOL)->domain().domainBox();
         const Box& lsol_domain_box  = getCoordSys(SingleNullBlockCoordSys::LSOL)->domain().domainBox();
         const Box& rsol_domain_box  = getCoordSys(SingleNullBlockCoordSys::RSOL)->domain().domainBox();

         Box fab_box = minBox( bdryHi(lcsol_domain_box, RADIAL_DIR, 1),
                               bdryHi(rcsol_domain_box, RADIAL_DIR, 1) );
         fab_box = minBox( bdryHi(lsol_domain_box, RADIAL_DIR, 1), fab_box );
         fab_box = minBox( bdryHi(rsol_domain_box, RADIAL_DIR, 1), fab_box );

         FArrayBox fab(fab_box, 1);
         data_array->setData(fab, false);
      }
   }

   // Radial private flux boundary
   RefCountedPtr<GridFunction> radial_pf_function = bc.getPFFunction();
   if ( !radial_pf_function.isNull() ) {
      if ( typeid(*radial_pf_function) == typeid(DataArray) ) {
         RefCountedPtr<DataArray> data_array = (RefCountedPtr<DataArray>)radial_pf_function;
         CH_assert( !data_array.isNull() );

         const Box& lpf_domain_box = getCoordSys(SingleNullBlockCoordSys::LPF)->domain().domainBox();
         const Box& rpf_domain_box = getCoordSys(SingleNullBlockCoordSys::RPF)->domain().domainBox();

         Box fab_box = minBox( bdryLo(lpf_domain_box, RADIAL_DIR, 1),
                               bdryLo(rpf_domain_box, RADIAL_DIR, 1) );

         FArrayBox fab(fab_box, 1);
         data_array->setData(fab, false);
      }
   }

   // Poloidal left divertor boundary
   RefCountedPtr<GridFunction> radial_left_function = bc.getInnerDivertorFunction();
   if ( !radial_left_function.isNull() ) {
      if ( typeid(*radial_left_function) == typeid(DataArray) ) {
         RefCountedPtr<DataArray> data_array = (RefCountedPtr<DataArray>)radial_left_function;
         CH_assert( !data_array.isNull() );

         const Box& lsol_domain_box = getCoordSys(SingleNullBlockCoordSys::LSOL)->domain().domainBox();
         const Box& lpf_domain_box = getCoordSys(SingleNullBlockCoordSys::LPF)->domain().domainBox();

         Box fab_box = minBox( bdryHi(lsol_domain_box, POLOIDAL_DIR, 1),
                               bdryHi(lpf_domain_box, POLOIDAL_DIR, 1) );

         FArrayBox fab(fab_box, 1);
         data_array->setData(fab, false);
      }
   }

   // Poloidal right divertor boundary
   RefCountedPtr<GridFunction> radial_right_function = bc.getOuterDivertorFunction();
   if ( !radial_right_function.isNull() ) {
      if ( typeid(*radial_right_function) == typeid(DataArray) ) {
         RefCountedPtr<DataArray> data_array = (RefCountedPtr<DataArray>)radial_right_function;
         CH_assert( !data_array.isNull() );

         const Box& rsol_domain_box = getCoordSys(SingleNullBlockCoordSys::RSOL)->domain().domainBox();
         const Box& rpf_domain_box = getCoordSys(SingleNullBlockCoordSys::RPF)->domain().domainBox();

         Box fab_box = minBox( bdryLo(rsol_domain_box, POLOIDAL_DIR, 1),
                               bdryLo(rpf_domain_box, POLOIDAL_DIR, 1) );

         FArrayBox fab(fab_box, 1);
         data_array->setData(fab, false);
      }
   }
}


IntVect
SingleNullCoordSys::lo_mapped_index(int a_block_number) const
{
  IntVect index;

  if ( m_original_eight_blocks ) {

     int numcells_core_poloidal = m_numcells_lcore_poloidal + m_numcells_rcore_poloidal;

     switch ( poloidalBlockNumber(a_block_number) )
        {
        case SingleNullBlockCoordSys::LCORE:
           index[RADIAL_DIR] = 0;
           index[POLOIDAL_DIR] = POLOIDAL_BLOCK_SEP*numcells_core_poloidal/2;
           break;
        case SingleNullBlockCoordSys::RCORE:
           index[RADIAL_DIR] = 0;
           index[POLOIDAL_DIR] = -(POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2;
           break;
        case SingleNullBlockCoordSys::LCSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*m_numcells_core_radial;
           index[POLOIDAL_DIR] = POLOIDAL_BLOCK_SEP*numcells_core_poloidal/2;
           break;
        case SingleNullBlockCoordSys::RCSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*m_numcells_core_radial;
           index[POLOIDAL_DIR] = -(POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2;
           break;
        case SingleNullBlockCoordSys::LSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*m_numcells_core_radial;      // same radial extent as CSOL
           index[POLOIDAL_DIR] = (2*POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2;
           break;
        case SingleNullBlockCoordSys::RSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*m_numcells_core_radial;      // same radial extent as CSOL
           index[POLOIDAL_DIR] = -(2*POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2 - m_numcells_rpf_poloidal;
           break;
        case SingleNullBlockCoordSys::LPF:
           index[RADIAL_DIR] = m_numcells_core_radial - m_numcells_pf_radial;
           index[POLOIDAL_DIR] = (2*POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2;
           break;
        case SingleNullBlockCoordSys::RPF:
           index[RADIAL_DIR] = m_numcells_core_radial - m_numcells_pf_radial;
           index[POLOIDAL_DIR] = -(2*POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2 - m_numcells_rpf_poloidal;
           break;
        default:
           MayDay::Error("SingleNullCoordSys::lo_mapped_index(): Invalid poloidal block number");
        }
  }
  else {

     int L0 = m_numcells_core_radial;
     int lL1 = m_numcells_mcore_poloidal/2 + m_numcells_lcore_poloidal;
     int rL1 = m_numcells_mcore_poloidal/2 + m_numcells_rcore_poloidal;

     switch ( poloidalBlockNumber(a_block_number) )
        {
        case SingleNullBlockCoordSys::MCORE:
           index[RADIAL_DIR] = 0;
           index[POLOIDAL_DIR] = -m_numcells_mcore_poloidal/2;
           break;
        case SingleNullBlockCoordSys::LCORE:
           index[RADIAL_DIR] = 0;
           index[POLOIDAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*lL1;
           break;
        case SingleNullBlockCoordSys::RCORE:
           index[RADIAL_DIR] = 0;
           index[POLOIDAL_DIR] = -(POLOIDAL_BLOCK_SEP+1)*rL1 - m_numcells_rcore_poloidal;
           break;
        case SingleNullBlockCoordSys::LPF:
           index[RADIAL_DIR] = L0 - m_numcells_pf_radial;
           index[POLOIDAL_DIR] = 2*(POLOIDAL_BLOCK_SEP+1)*lL1;
           break;
        case SingleNullBlockCoordSys::RPF:
           index[RADIAL_DIR] = L0 - m_numcells_pf_radial;
           index[POLOIDAL_DIR] = -2*(POLOIDAL_BLOCK_SEP+1)*rL1 - m_numcells_rpf_poloidal;
           break;
        case SingleNullBlockCoordSys::MCSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*L0;
           index[POLOIDAL_DIR] = -m_numcells_mcore_poloidal/2;
           break;
        case SingleNullBlockCoordSys::LCSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*L0;
           index[POLOIDAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*lL1;
           break;
        case SingleNullBlockCoordSys::RCSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*L0;
           index[POLOIDAL_DIR] = -(POLOIDAL_BLOCK_SEP+1)*rL1 - m_numcells_rcore_poloidal;
           break;
        case SingleNullBlockCoordSys::LSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*L0;
           index[POLOIDAL_DIR] = 2*(POLOIDAL_BLOCK_SEP+1)*lL1;
           break;
        case SingleNullBlockCoordSys::RSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*L0;
           index[POLOIDAL_DIR] = -2*(POLOIDAL_BLOCK_SEP+1)*rL1 - m_numcells_rpf_poloidal;
           break;
        default:
           MayDay::Error("SingleNullCoordSys::lo_mapped_index(): Invalid poloidal block number");
        }
  }

#if CFG_DIM==3
  int num_block_toroidal_cells = m_numcells_toroidal / m_num_toroidal_sectors;
  
  index[TOROIDAL_DIR] = toroidalBlockNumber(a_block_number) * (num_block_toroidal_cells + TOROIDAL_BLOCK_SEP);
#endif

  return index;
}


IntVect
SingleNullCoordSys::hi_mapped_index(int a_block_number) const
{
  IntVect index;

  if ( m_original_eight_blocks ) {

     int numcells_core_poloidal = m_numcells_lcore_poloidal + m_numcells_rcore_poloidal;

     switch ( poloidalBlockNumber(a_block_number) )
        {
        case SingleNullBlockCoordSys::LCORE:
           index[RADIAL_DIR] = m_numcells_core_radial - 1;
           index[POLOIDAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2 - 1;
           break;
        case SingleNullBlockCoordSys::RCORE:
           index[RADIAL_DIR] = m_numcells_core_radial - 1;
           index[POLOIDAL_DIR] = -POLOIDAL_BLOCK_SEP*numcells_core_poloidal/2 - 1;
           break;
        case SingleNullBlockCoordSys::LCSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*m_numcells_core_radial + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2 - 1;
           break;
        case SingleNullBlockCoordSys::RCSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*m_numcells_core_radial + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = -POLOIDAL_BLOCK_SEP*numcells_core_poloidal/2 - 1;
           break;
        case SingleNullBlockCoordSys::LSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*m_numcells_core_radial + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = (2*POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2 + m_numcells_lpf_poloidal - 1;
           break;
        case SingleNullBlockCoordSys::RSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*m_numcells_core_radial + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = -(2*POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2 - 1;  // poloidal domain shifted downward
           break;
        case SingleNullBlockCoordSys::LPF:
           index[RADIAL_DIR] = m_numcells_core_radial - 1;                                // same hi radial boundary as CORE
           index[POLOIDAL_DIR] = (2*POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2 + m_numcells_lpf_poloidal - 1;
           break;
        case SingleNullBlockCoordSys::RPF:
           index[RADIAL_DIR] = m_numcells_core_radial - 1;                           // same hi radial boundary as CORE
           index[POLOIDAL_DIR] = -(2*POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2 - 1;  // poloidal domain shifted downward
           break;
        default:
           MayDay::Error("SingleNullCoordSys::lo_mapped_index(): Invalid block number");
        }
  }
  else {

     int L0 = m_numcells_core_radial;
     int lL1 = m_numcells_mcore_poloidal/2 + m_numcells_lcore_poloidal;
     int rL1 = m_numcells_mcore_poloidal/2 + m_numcells_rcore_poloidal;

     switch ( poloidalBlockNumber(a_block_number) )
        {
        case SingleNullBlockCoordSys::MCORE:
           index[RADIAL_DIR] = L0 - 1;
           index[POLOIDAL_DIR] = m_numcells_mcore_poloidal/2 - 1;
           break;
        case SingleNullBlockCoordSys::LCORE:
           index[RADIAL_DIR] = L0 - 1;
           index[POLOIDAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*lL1 + m_numcells_lcore_poloidal - 1;
           break;
        case SingleNullBlockCoordSys::RCORE:
           index[RADIAL_DIR] = L0 - 1;
           index[POLOIDAL_DIR] = -(POLOIDAL_BLOCK_SEP+1)*rL1 - 1;
           break;
        case SingleNullBlockCoordSys::LPF:
           index[RADIAL_DIR] = L0 - 1;
           index[POLOIDAL_DIR] = 2*(POLOIDAL_BLOCK_SEP+1)*lL1 + m_numcells_lpf_poloidal - 1; 
           break;
        case SingleNullBlockCoordSys::RPF:
           index[RADIAL_DIR] = L0 - 1;
           index[POLOIDAL_DIR] = -2*(POLOIDAL_BLOCK_SEP+1)*rL1 - 1;
           break;
        case SingleNullBlockCoordSys::MCSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*L0 + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = m_numcells_mcore_poloidal/2 - 1;
           break;
        case SingleNullBlockCoordSys::LCSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*L0 + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*lL1 + m_numcells_lcore_poloidal - 1;
           break;
        case SingleNullBlockCoordSys::RCSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*L0 + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = -(POLOIDAL_BLOCK_SEP+1)*rL1 - 1;
           break;
        case SingleNullBlockCoordSys::LSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*L0 + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = 2*(POLOIDAL_BLOCK_SEP+1)*lL1 + m_numcells_lpf_poloidal - 1;
           break;
        case SingleNullBlockCoordSys::RSOL:
           index[RADIAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*L0 + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = -2*(POLOIDAL_BLOCK_SEP+1)*rL1 - 1;
           break;
        default:
           MayDay::Error("SingleNullCoordSys::lo_mapped_index(): Invalid block number");
        }
  }

#if CFG_DIM==3
  int num_block_toroidal_cells = m_numcells_toroidal / m_num_toroidal_sectors;
  
  index[TOROIDAL_DIR] = toroidalBlockNumber(a_block_number) * (num_block_toroidal_cells + TOROIDAL_BLOCK_SEP)
                        + num_block_toroidal_cells - 1;
#endif

  return index;
}



IntVect
SingleNullCoordSys::getDecomposition( int a_block_number ) const
{
  IntVect decomp;

  switch ( poloidalBlockNumber(a_block_number) )
    {
    case SingleNullBlockCoordSys::MCORE:
      decomp = m_decomp_mcore;
      break;
    case SingleNullBlockCoordSys::LCORE:
      decomp = m_decomp_lcore;
      break;
    case SingleNullBlockCoordSys::RCORE:
      decomp = m_decomp_rcore;
      break;
    case SingleNullBlockCoordSys::MCSOL:
      decomp = m_decomp_mcsol;
      break;
    case SingleNullBlockCoordSys::LCSOL:
      decomp = m_decomp_lcsol;
      break;
    case SingleNullBlockCoordSys::RCSOL:
      decomp = m_decomp_rcsol;
      break;
    case SingleNullBlockCoordSys::LSOL:
      decomp = m_decomp_lsol;
      break;
    case SingleNullBlockCoordSys::RSOL:
      decomp = m_decomp_rsol;
      break;
    case SingleNullBlockCoordSys::LPF:
      decomp = m_decomp_lpf;
      break;
    case SingleNullBlockCoordSys::RPF:
      decomp = m_decomp_rpf;
      break;
    default:
      MayDay::Error("SingleNullCoordSys::getDecomposition(): Invalid block number");
    }

  return decomp;
}



void
SingleNullCoordSys::readParams( ParmParse& a_pp_grid,
                                ParmParse& a_pp_geom )
{
#if CFG_DIM==3
   if ( a_pp_geom.contains("num_toroidal_sectors") ) {
      a_pp_geom.query( "num_toroidal_sectors", m_num_toroidal_sectors);
   }

   if ( a_pp_geom.contains("interpolation_scheme") ) {
      a_pp_geom.query( "interpolation_scheme", m_interpolation_scheme);
   }
   else {
     m_interpolation_scheme = "mapped";
   }
   
   a_pp_geom.get( "toroidal_width_over_2pi", m_toroidal_width);
   CH_assert(m_toroidal_width > 0. && m_toroidal_width <= 1.);
   m_toroidal_width *= 2.*Pi;
#endif

   std::vector<int> decomp_mcore( SpaceDim );
   std::vector<int> decomp_lcore( SpaceDim );
   std::vector<int> decomp_rcore( SpaceDim );
   std::vector<int> decomp_lpf( SpaceDim );
   std::vector<int> decomp_rpf( SpaceDim );
   std::vector<int> decomp_mcsol( SpaceDim );
   std::vector<int> decomp_lcsol( SpaceDim );
   std::vector<int> decomp_rcsol( SpaceDim );
   std::vector<int> decomp_lsol( SpaceDim );
   std::vector<int> decomp_rsol( SpaceDim );

   // If the input specifies the mesh and decomposition for the
   // MCORE block, we assume a 10 poloidal block geometry.  Otherwise,
   // we assume the original 8 poloidal block geometry
   if ( a_pp_grid.contains("numcells.mcore_poloidal") &&
        a_pp_grid.contains("decomp.mcore.configuration") &&
        a_pp_grid.contains("decomp.mcsol.configuration") ) {

      m_original_eight_blocks = false;

      a_pp_grid.query( "numcells.mcore_poloidal", m_numcells_mcore_poloidal );

      // For now, we require that the MCORE has an even number of cells in the 
      // poloidal direction so that it can be located symmetrically about
      // zero in index space.  This assumption can be removed later if needed.
      if ( m_numcells_mcore_poloidal%2 != 0 ) {
         MayDay::Error("SingleNullCoordSys::readGridParams(): Number of mcore cells in poloidal direction must be even");
      }

      a_pp_grid.queryarr( "decomp.mcore.configuration", decomp_mcore, 0, SpaceDim );
      a_pp_grid.queryarr( "decomp.mcsol.configuration", decomp_mcsol, 0, SpaceDim );
   }
   else {
      m_original_eight_blocks = true;
   }

   if ( m_original_eight_blocks && a_pp_grid.contains("numcells.core" ) ) {
      // For backward compatibility with the old way of specifying an eight
      // block single null geometry:
      std::vector<int> numcells_core(SpaceDim);
      a_pp_grid.queryarr( "numcells.core", numcells_core, 0, SpaceDim );

      m_numcells_core_radial = numcells_core[RADIAL_DIR];
      if ( numcells_core[POLOIDAL_DIR]%2 != 0 ) {
         MayDay::Error("SingleNullCoordSys::readGridParams(): Number of core cells in poloidal direction must be even");
      }
      else {
         m_numcells_lcore_poloidal = numcells_core[POLOIDAL_DIR]/2;
         m_numcells_rcore_poloidal = numcells_core[POLOIDAL_DIR]/2;
      }

#if CFG_DIM==3
      m_numcells_toroidal = numcells_core[TOROIDAL_DIR];
#endif
   }
   else {
      a_pp_grid.query( "numcells.core_radial", m_numcells_core_radial );
      a_pp_grid.query( "numcells.lcore_poloidal", m_numcells_lcore_poloidal );
      a_pp_grid.query( "numcells.rcore_poloidal", m_numcells_rcore_poloidal );

#if CFG_DIM==3
      a_pp_grid.query( "numcells.toroidal", m_numcells_toroidal );
#endif
   }

   a_pp_grid.query( "numcells.pf_radial", m_numcells_pf_radial );
   a_pp_grid.query( "numcells.lpf_poloidal", m_numcells_lpf_poloidal );
   a_pp_grid.query( "numcells.rpf_poloidal", m_numcells_rpf_poloidal );
   a_pp_grid.query( "numcells.sol_radial", m_numcells_sol_radial );

   if ( m_original_eight_blocks && a_pp_grid.contains("decomp.core.configuration" ) ) {
      // For backward compatibility with the old way of specifying an eight
      // block single null geometry:
      std::vector<int> decomp_core( SpaceDim );
      a_pp_grid.queryarr( "decomp.core.configuration", decomp_core, 0, SpaceDim );

      for (int dir=0; dir<SpaceDim; ++dir) {
         decomp_lcore[dir] = decomp_core[dir];
         decomp_rcore[dir] = decomp_core[dir];
      }
      if ( decomp_core[POLOIDAL_DIR]%2 != 0 ) {
         MayDay::Error("SingleNullCoordSys::readGridParams(): Core decomposition in poloidal direction must be even");
      }
      else {
         decomp_lcore[POLOIDAL_DIR] /= 2;
         decomp_rcore[POLOIDAL_DIR] /= 2;
      }
   }
   else {
      a_pp_grid.queryarr( "decomp.lcore.configuration", decomp_lcore, 0, SpaceDim );
      a_pp_grid.queryarr( "decomp.rcore.configuration", decomp_rcore, 0, SpaceDim );
   }

   if ( m_original_eight_blocks && a_pp_grid.contains("decomp.csol.configuration" ) ) {
      // For backward compatibility with the old way of specifying an eight
      // block single null geometry:
      std::vector<int> decomp_csol( SpaceDim );
      a_pp_grid.queryarr( "decomp.csol.configuration", decomp_csol, 0, SpaceDim );

      for (int dir=0; dir<SpaceDim; ++dir) {
         decomp_lcsol[dir] = decomp_csol[dir];
         decomp_rcsol[dir] = decomp_csol[dir];
      }
      if ( decomp_csol[POLOIDAL_DIR]%2 != 0 ) {
         MayDay::Error("SingleNullCoordSys::readGridParams(): CSOL decomposition in poloidal direction must be even");
      }
      else {
         decomp_lcsol[POLOIDAL_DIR] /= 2;
         decomp_rcsol[POLOIDAL_DIR] /= 2;
      }
   }
   else {
      a_pp_grid.queryarr( "decomp.lcsol.configuration", decomp_lcsol, 0, SpaceDim );
      a_pp_grid.queryarr( "decomp.rcsol.configuration", decomp_rcsol, 0, SpaceDim );
   }

   a_pp_grid.queryarr( "decomp.lpf.configuration", decomp_lpf, 0, SpaceDim );
   a_pp_grid.queryarr( "decomp.rpf.configuration", decomp_rpf, 0, SpaceDim );
   a_pp_grid.queryarr( "decomp.lsol.configuration", decomp_lsol, 0, SpaceDim );
   a_pp_grid.queryarr( "decomp.rsol.configuration", decomp_rsol, 0, SpaceDim );

   if ( procID() == 0 ) {
      cout << endl << "SingleNull grid and domain decomposition parameters:" << endl << endl;

      cout << "numcells.core_radial = " << m_numcells_core_radial << endl;
      if ( !m_original_eight_blocks ) cout << "numcells.mcore_poloidal = " << m_numcells_mcore_poloidal << endl;
      cout << "numcells.lcore_poloidal = " << m_numcells_lcore_poloidal << endl;
      cout << "numcells.rcore_poloidal = " << m_numcells_rcore_poloidal << endl;
      cout << "numcells.pf_radial = " << m_numcells_pf_radial << endl;
      cout << "numcells.lpf_poloidal = " << m_numcells_lpf_poloidal <<endl;
      cout << "numcells.rpf_poloidal = " << m_numcells_rpf_poloidal << endl;
      cout << "numcells.sol_radial = " << m_numcells_sol_radial << endl;

      if ( !m_original_eight_blocks ) {
         cout << "decomp.mcore.configuration = ";
         for (int i=0; i<SpaceDim; i++) cout << decomp_mcore[i] << " ";
         cout << endl;
      }

      cout << "decomp.lcore.configuration = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_lcore[i] << " ";
      cout << endl;

      cout << "decomp.rcore.configuration = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_rcore[i] << " ";
      cout << endl;

      cout << "decomp.lpf.configuration = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_lpf[i] << " ";
      cout << endl;

      cout << "decomp.rpf.configuration = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_rpf[i] << " ";
      cout << endl;

      if ( !m_original_eight_blocks ) {
         cout << "decomp.mcsol.configuration = ";
         for (int i=0; i<SpaceDim; i++) cout << decomp_mcsol[i] << " ";
         cout << endl;
      }

      cout << "decomp.lcsol.configuration = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_lcsol[i] << " ";
      cout << endl;

      cout << "decomp.rcsol.configuration = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_rcsol[i] << " ";
      cout << endl;

      cout << "decomp.lsol.configuration = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_lsol[i] << " ";
      cout << endl;

      cout << "decomp.rsol.configuration = ";
      for (int i=0; i<SpaceDim; i++) cout << decomp_rsol[i] << " ";
      cout << endl;

      cout << endl;
   }

   // Check decomposition divisibilities

   if ( !m_original_eight_blocks && m_numcells_core_radial % decomp_mcore[RADIAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Radial decomposition does not divide number of cells in MCORE");
   }
   if ( !m_original_eight_blocks && m_numcells_mcore_poloidal % decomp_mcore[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Poloidal decomposition does not divide number of cells in MCORE");
   }

   if ( m_numcells_core_radial % decomp_lcore[RADIAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Radial decomposition does not divide number of cells in LCORE");
   }
   if ( m_numcells_lcore_poloidal % decomp_lcore[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Poloidal decomposition does not divide number of cells in LCORE");
   }

   if ( m_numcells_core_radial % decomp_rcore[RADIAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Radial decomposition does not divide number of cells in RCORE");
   }
   if ( m_numcells_rcore_poloidal % decomp_rcore[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Poloidal decomposition does not divide number of cells in RCORE");
   }

   if ( m_numcells_pf_radial % decomp_lpf[RADIAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Radial decomposition does not divide number of cells in LPF");
   }
   if ( m_numcells_lpf_poloidal % decomp_lpf[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Poloidal decomposition does not divide number of cells in LPF");
   }

   if ( m_numcells_pf_radial % decomp_rpf[RADIAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Radial decomposition does not divide number of cells in RPF");
   }
   if ( m_numcells_rpf_poloidal % decomp_rpf[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Poloidal decomposition does not divide number of cells in RPF");
   }

   if ( !m_original_eight_blocks && m_numcells_sol_radial % decomp_mcsol[RADIAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Radial decomposition does not divide number of cells in MCSOL");
   }
   if ( !m_original_eight_blocks && m_numcells_mcore_poloidal % decomp_mcsol[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Poloidal decomposition does not divide number of cells in MCSOL");
   }

   if ( m_numcells_sol_radial % decomp_lcsol[RADIAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Radial decomposition does not divide number of cells in LCSOL");
   }
   if ( m_numcells_lcore_poloidal % decomp_lcsol[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Poloidal decomposition does not divide number of cells in LCSOL");
   }

   if ( m_numcells_sol_radial % decomp_rcsol[RADIAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Radial decomposition does not divide number of cells in RCSOL");
   }
   if ( m_numcells_rcore_poloidal % decomp_rcsol[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Poloidal decomposition does not divide number of cells in RCSOL");
   }

   if ( m_numcells_sol_radial % decomp_lsol[RADIAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Radial decomposition does not divide number of cells in LSOL");
   }
   if ( m_numcells_lpf_poloidal % decomp_lsol[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Poloidal decomposition does not divide number of cells in LSOL");
   }

   if ( m_numcells_sol_radial % decomp_rsol[RADIAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Radial decomposition does not divide number of cells in RSOL");
   }
   if ( m_numcells_rpf_poloidal % decomp_rsol[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Poloidal decomposition does not divide number of cells in RSOL");
   }

#if CFG_DIM==3
   if ( decomp_mcore[TOROIDAL_DIR] % m_num_toroidal_sectors != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in MCORE");
   }
   else {
      decomp_mcore[TOROIDAL_DIR] /= m_num_toroidal_sectors;
   }
   if ( decomp_lcore[TOROIDAL_DIR] % m_num_toroidal_sectors != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in LCORE");
   }
   else {
      decomp_lcore[TOROIDAL_DIR] /= m_num_toroidal_sectors;
   }
   if ( decomp_rcore[TOROIDAL_DIR] % m_num_toroidal_sectors != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in RCORE");
   }
   else {
      decomp_rcore[TOROIDAL_DIR] /= m_num_toroidal_sectors;
   }
   if ( decomp_lpf[TOROIDAL_DIR] % m_num_toroidal_sectors != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in LPF");
   }
   else {
      decomp_lpf[TOROIDAL_DIR] /= m_num_toroidal_sectors;
   }
   if ( decomp_rpf[TOROIDAL_DIR] % m_num_toroidal_sectors != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in RPF");
   }
   else {
      decomp_rpf[TOROIDAL_DIR] /= m_num_toroidal_sectors;
   }
   if ( decomp_mcsol[TOROIDAL_DIR] % m_num_toroidal_sectors != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in MCSOL");
   }
   else {
      decomp_mcsol[TOROIDAL_DIR] /= m_num_toroidal_sectors;
   }
   if ( decomp_lcsol[TOROIDAL_DIR] % m_num_toroidal_sectors != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in LCSOL");
   }
   else {
      decomp_lcsol[TOROIDAL_DIR] /= m_num_toroidal_sectors;
   }
   if ( decomp_rcsol[TOROIDAL_DIR] % m_num_toroidal_sectors != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in RCSOL");
   }
   else {
      decomp_rcsol[TOROIDAL_DIR] /= m_num_toroidal_sectors;
   }
   if ( decomp_lsol[TOROIDAL_DIR] % m_num_toroidal_sectors != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in LSOL");
   }
   else {
      decomp_lsol[TOROIDAL_DIR] /= m_num_toroidal_sectors;
   }
   if ( decomp_rsol[TOROIDAL_DIR] % m_num_toroidal_sectors != 0 ) {
      MayDay::Error("SingleNullCoordSys::readGridParams(): Toroidal decomposition is not a multiple of the number of toroidal sectors in RSOL");
   }
   else {
      decomp_rsol[TOROIDAL_DIR] /= m_num_toroidal_sectors;
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

   a_pp_geom.query( "model_geometry", m_model_geometry);

   m_provides_flux = a_pp_geom.contains("field_coefficients_file");
}


void
accumulateTuples( const Box&                  a_box,
                  const IntVect&              a_shift_vec,
                  const int                   a_block_number,
                  Vector< Tuple<IntVect,2> >& a_iv_pairs,
                  Vector<int>&                a_block_numbers,
                  IntVectSet&                 a_ghost_ivs )
{
   for (BoxIterator bit(a_box); bit.ok(); ++bit) {
      Tuple<IntVect,2> t;
      t[0] = bit();
      t[1] = bit() - a_shift_vec;
      a_iv_pairs.push_back(t);
      a_block_numbers.push_back(a_block_number);
      a_ghost_ivs |= t[1];
   }
}


void
SingleNullCoordSys::defineStencilsUe( const DisjointBoxLayout&  a_grids,
                                     const int                 a_nghost )
{
   int num_blocks = m_original_eight_blocks? 8: SingleNullBlockCoordSys::NUM_POLOIDAL_BLOCKS;
   
   if ( num_blocks == 8 ) {
      defineStencilsUe8(a_grids, a_nghost);
   }
   else {
      MayDay::Error("SingleNullCoordSys::defineStencilsUe() is only implemented for 8 blocks");
   }
}


void
SingleNullCoordSys::defineStencilsUe8( const DisjointBoxLayout&  a_grids,
                                      const int                 a_nghost )
{
   m_stencil_ue.define(a_grids);
   m_ghostCells_ue.define(a_grids);
   
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      int block_number = whichBlock(a_grids[dit]);
      const MagBlockCoordSys* coord_sys = getCoordSys(block_number);
      const Box& domain_box = coord_sys->domain().domainBox();
      
      Vector< Tuple<IntVect,2> > iv_pairs;
      Vector<int> block_numbers;
      IntVectSet& ghost_ivs = m_ghostCells_ue[dit];
      
      Box grown_box = grow(a_grids[dit], a_nghost);
      
      if ( !domain_box.contains(grown_box) ) {
         
         Box rad_shifted_box;
         Box pol_shifted_box;
         IntVect shift_r, shift_p;
#if CFG_DIM==3
         Box tor_shifted_box;
         IntVect shift_t;
#endif

         int toroidal_block_number = toroidalBlockNumber(block_number);
         int LCORE = blockNumber(SingleNullBlockCoordSys::LCORE, toroidal_block_number);
         int RCORE = blockNumber(SingleNullBlockCoordSys::RCORE, toroidal_block_number);
         int LCSOL = blockNumber(SingleNullBlockCoordSys::LCSOL, toroidal_block_number);
         int RCSOL = blockNumber(SingleNullBlockCoordSys::RCSOL, toroidal_block_number);
         int LSOL  = blockNumber(SingleNullBlockCoordSys::LSOL,  toroidal_block_number);
         int RSOL  = blockNumber(SingleNullBlockCoordSys::RSOL,  toroidal_block_number);
         int LPF   = blockNumber(SingleNullBlockCoordSys::LPF,   toroidal_block_number);
         int RPF   = blockNumber(SingleNullBlockCoordSys::RPF,   toroidal_block_number);
         
         switch ( poloidalBlockNumber(block_number) )
         {
            case SingleNullBlockCoordSys::LCORE:
               
               // RCORE
               shift_p = -POLOIDAL_BLOCK_SEP * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                                m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCORE)->domain().domainBox(), shift_p, RCORE,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               shift_p = -(POLOIDAL_BLOCK_SEP+1) * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                                    m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCORE)->domain().domainBox(), shift_p, RCORE,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // LCSOL
               shift_r = POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(LCSOL)->domain().domainBox(), shift_r, LCSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // RCSOL
               shift_p = -POLOIDAL_BLOCK_SEP * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                                m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCSOL)->domain().domainBox(), shift_r + shift_p, RCSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // LSOL
               shift_p = POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LSOL)->domain().domainBox(), shift_r + shift_p, LSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case SingleNullBlockCoordSys::RCORE:
               
               // LCORE
               shift_p = (POLOIDAL_BLOCK_SEP+1) * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                                 + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCORE)->domain().domainBox(), shift_p, LCORE,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               shift_p = POLOIDAL_BLOCK_SEP * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                             + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCORE)->domain().domainBox(), shift_p, LCORE,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // RCSOL
               shift_r = POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(RCSOL)->domain().domainBox(), shift_r, RCSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // LCSOL
               shift_p = POLOIDAL_BLOCK_SEP * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                               m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCSOL)->domain().domainBox(), shift_r + shift_p, LCSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // RSOL
               shift_p = -POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RSOL)->domain().domainBox(), shift_r + shift_p, RSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               break;
         case SingleNullBlockCoordSys::LCSOL:
               
               // LSOL
               shift_p = POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LSOL)->domain().domainBox(), shift_p, LSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // RCSOL
               shift_p = -POLOIDAL_BLOCK_SEP * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                                m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCSOL)->domain().domainBox(), shift_p, RCSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // LCORE
               shift_r = -POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(LCORE)->domain().domainBox(), shift_r, LCORE,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // RCORE
               shift_p = -POLOIDAL_BLOCK_SEP*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                              m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCORE)->domain().domainBox(), shift_r + shift_p, RCORE,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // LPF
               shift_p = POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LPF)->domain().domainBox(), shift_r + shift_p, LPF,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case SingleNullBlockCoordSys::RCSOL:
               
               // RSOL
               shift_p = -POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RSOL)->domain().domainBox(), shift_p, RSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // LCSOL
               shift_p = POLOIDAL_BLOCK_SEP * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                               m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCSOL)->domain().domainBox(), shift_p, LCSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // RCORE
               shift_r = -POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(RCORE)->domain().domainBox(), shift_r, RCORE,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // LCORE
               shift_p = POLOIDAL_BLOCK_SEP * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                             + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCORE)->domain().domainBox(), shift_r + shift_p, LCORE,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // RPF
               shift_p = -POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RPF)->domain().domainBox(), shift_r + shift_p, RPF,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case SingleNullBlockCoordSys::LSOL:
               
               // LCSOL
               shift_p = -POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCSOL)->domain().domainBox(), shift_p, LCSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // LPF
               shift_r = -POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(LPF)->domain().domainBox(), shift_r, LPF,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // LCORE
               shift_p = -POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCORE)->domain().domainBox(), shift_r + shift_p, LCORE,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case SingleNullBlockCoordSys::RSOL:
               
               // RCSOL
               shift_p = POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCSOL)->domain().domainBox(), shift_p, RCSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // RPF
               shift_r = -POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(RPF)->domain().domainBox(), shift_r, RPF,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // RCORE
               shift_p = POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCORE)->domain().domainBox(), shift_r + shift_p, RCORE,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case SingleNullBlockCoordSys::LPF:
               
               // RPF
               shift_p = -(2*POLOIDAL_BLOCK_SEP+1)*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                                  + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RPF)->domain().domainBox(), shift_p, RPF,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // LSOL
               shift_r = POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(LSOL)->domain().domainBox(), shift_r,  LSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // LCSOL
               shift_p = -POLOIDAL_BLOCK_SEP * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCSOL)->domain().domainBox(), shift_r + shift_p, LCSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case SingleNullBlockCoordSys::RPF:
               
               // LPF
               shift_p = (2*POLOIDAL_BLOCK_SEP+1) * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                                   + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LPF)->domain().domainBox(), shift_p, LPF,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // RSOL
               shift_r = POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(RSOL)->domain().domainBox(), shift_r, RSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               // RCSOL
               shift_p = POLOIDAL_BLOCK_SEP * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCSOL)->domain().domainBox(), shift_r + shift_p, RCSOL,
                                 iv_pairs, block_numbers, ghost_ivs );
               
               break;
            default:
               MayDay::Error("SingleNullCoordSys::blockRemapping(): Unknown poloidal block number");
         }

#if CFG_DIM==3
         if (m_num_toroidal_sectors == 1) {
            shift_t = m_mappingBlocks[block_number].size(TOROIDAL_DIR) * BASISV(TOROIDAL_DIR);
         }
         else {
            shift_t = TOROIDAL_BLOCK_SEP * m_mappingBlocks[block_number].size(TOROIDAL_DIR) * BASISV(TOROIDAL_DIR);
         }

         int other_tor_block_number = toroidal_block_number;

         if ( m_num_toroidal_sectors > 1 && toroidal_block_number < m_num_toroidal_sectors-1 ) {
            other_tor_block_number++;
         }
         
         tor_shifted_box = grown_box + shift_t;
         for (int other_pol_block_number=0; other_pol_block_number<m_num_poloidal_blocks; ++other_pol_block_number) {
            int other_pol_block = blockNumber(other_pol_block_number, other_tor_block_number);
            accumulateTuples( tor_shifted_box & getCoordSys(other_pol_block)->domain().domainBox(), shift_t,
                              other_pol_block, iv_pairs, block_numbers, ghost_ivs );
         }

         if (m_num_toroidal_sectors == 1) {
            shift_t = - m_mappingBlocks[block_number].size(TOROIDAL_DIR) * BASISV(TOROIDAL_DIR);
         }
         else {
            shift_t = -TOROIDAL_BLOCK_SEP * m_mappingBlocks[block_number].size(TOROIDAL_DIR) * BASISV(TOROIDAL_DIR);
         }

         if ( m_num_toroidal_sectors > 1 && toroidal_block_number > 0 ) {
            other_tor_block_number = toroidal_block_number-1;
         }

         tor_shifted_box = grown_box + shift_t;
         for (int other_pol_block_number=0; other_pol_block_number<m_num_poloidal_blocks; ++other_pol_block_number) {
            int other_pol_block = blockNumber(other_pol_block_number, other_tor_block_number);
            accumulateTuples( tor_shifted_box & getCoordSys(other_pol_block)->domain().domainBox(), shift_t,
                              other_pol_block, iv_pairs, block_numbers, ghost_ivs );
         }
#endif
         
         // Now that we have all of the extrablock ghost cells, we can allocate the IVSFAB<MBStencil>
         m_stencil_ue[dit] = RefCountedPtr< IVSFAB<MBStencil> >(new IVSFAB<MBStencil>(ghost_ivs,1));
         IVSFAB<MBStencil>& this_ivs_fab = *m_stencil_ue[dit];
         
         // Make a single element stencil on each extrablock ghost cell
         for (IVSIterator it(ghost_ivs); it.ok(); ++it) {
            IntVect extra_block_iv = it();
            
            Vector<MBStencilElement> stencilVec(1);
            
            bool found_iv = false;
            for (int n=0; n<iv_pairs.size(); ++n) {
               if ( iv_pairs[n][1] == extra_block_iv ) {
                  stencilVec[0].define(iv_pairs[n][0], block_numbers[n], 1.);
                  found_iv = true;
                  break;
               }
            }
            CH_assert(found_iv);
            
            this_ivs_fab(extra_block_iv,0) = RefCountedPtr< Vector<MBStencilElement> >
            (new Vector<MBStencilElement>(stencilVec));
         }
      }
   }
}


bool
SingleNullCoordSys::isCORE(const int a_block_number) const
{
   bool result(false);

   int poloidal_block = poloidalBlockNumber(a_block_number);
   
   if (poloidal_block == SingleNullBlockCoordSys::LCORE ||
       poloidal_block == SingleNullBlockCoordSys::RCORE ||
       poloidal_block == SingleNullBlockCoordSys::MCORE) {

      result = true;
   }
   
   return result;
}

bool
SingleNullCoordSys::isSOL(const int a_block_number) const
{
   bool result(false);

   int poloidal_block = poloidalBlockNumber(a_block_number);
   
   if (poloidal_block == SingleNullBlockCoordSys::LCSOL ||
       poloidal_block == SingleNullBlockCoordSys::RCSOL ||
       poloidal_block == SingleNullBlockCoordSys::MCSOL ||
       poloidal_block == SingleNullBlockCoordSys::LSOL ||
       poloidal_block == SingleNullBlockCoordSys::RSOL ) {

      result = true;
   }
   
   return result;
}

bool
SingleNullCoordSys::isPF(const int a_block_number) const
{
   bool result(false);

   int poloidal_block = poloidalBlockNumber(a_block_number);
   
   if (poloidal_block == SingleNullBlockCoordSys::LPF ||
       poloidal_block == SingleNullBlockCoordSys::RPF ) {

      result = true;
   }
   
   return result;
}

void SingleNullCoordSys::extrapolateEfield( LevelData<FluxBox>&   a_Er_average_face,
                                            LevelData<FArrayBox>& a_Er_average_cell ) const
{
    //Old extrapolation function. Now use something else, but this is left for the record.
    
    const DisjointBoxLayout& grids = a_Er_average_cell.disjointBoxLayout();

    const MagBlockCoordSys& block0_coord_sys = (const MagBlockCoordSys&)(*(getCoordSys(0)));
    int npol_LCORE = block0_coord_sys.domain().domainBox().size(POLOIDAL_DIR);
    int lo_pol_LCORE = block0_coord_sys.domain().domainBox().smallEnd(POLOIDAL_DIR);
    int hi_radial_index = block0_coord_sys.domain().domainBox().bigEnd(RADIAL_DIR);
    
    const MagBlockCoordSys& block1_coord_sys = (const MagBlockCoordSys&)(*(getCoordSys(1)));
    int lo_pol_RCORE = block1_coord_sys.domain().domainBox().smallEnd(POLOIDAL_DIR);
    
    //Initialize 1D arrays to contain boundary data (i.e., physical E-field on the outer core boundary)
    int nghost = 0;
    int size = npol_LCORE + nghost;
    
    double *ER_LCORE_bnd     = new double[size];
    double *ER_LCORE_bnd_loc = new double[size];
    
    double *ER_RCORE_bnd     = new double[size];
    double *ER_RCORE_bnd_loc = new double[size];
    
    double *EZ_LCORE_bnd     = new double[size];
    double *EZ_LCORE_bnd_loc = new double[size];
    
    double *EZ_RCORE_bnd     = new double[size];
    double *EZ_RCORE_bnd_loc = new double[size];
    
    //Initilizing everything with zeros
    for (int i=0; i<size; ++i) {
        ER_LCORE_bnd[i] = 0.0;
        EZ_LCORE_bnd[i] = 0.0;
        ER_RCORE_bnd[i] = 0.0;
        EZ_RCORE_bnd[i] = 0.0;
        
        ER_LCORE_bnd_loc[i] = 0.0;
        EZ_LCORE_bnd_loc[i] = 0.0;
        ER_RCORE_bnd_loc[i] = 0.0;
        EZ_RCORE_bnd_loc[i] = 0.0;
    }
    
    // Fill boundary data
    DataIterator dit( a_Er_average_face.dataIterator() );
    for (dit.begin(); dit.ok(); ++dit) {
        int block_number( whichBlock( grids[dit] ) );
        if (block_number<2) {
            
            FArrayBox& this_E_dir0 = a_Er_average_face[dit][0];
            
            Box box( grids[dit] );
            BoxIterator bit(box);
            for (bit.begin(); bit.ok(); ++bit) {
                
                IntVect iv = bit();
                IntVect iv_face(iv);
                if (iv[0]==hi_radial_index) {
                    iv_face[0]=iv[0]+1;
                    if (block_number==0) {
                        int npol = iv[1]-lo_pol_LCORE;
                        ER_LCORE_bnd_loc[npol] = this_E_dir0(iv_face,0);
                        EZ_LCORE_bnd_loc[npol] = this_E_dir0(iv_face,2);
                    }
                    else {
                        int npol = iv[1]-lo_pol_RCORE;
                        ER_RCORE_bnd_loc[npol] = this_E_dir0(iv_face,0);
                        EZ_RCORE_bnd_loc[npol] = this_E_dir0(iv_face,2);
                    }
                }
            }
        }
    }
    
    MPI_Allreduce(ER_RCORE_bnd_loc, ER_RCORE_bnd, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(EZ_RCORE_bnd_loc, EZ_RCORE_bnd, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    MPI_Allreduce(ER_LCORE_bnd_loc, ER_LCORE_bnd, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(EZ_LCORE_bnd_loc, EZ_LCORE_bnd, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    
    // Copy boundary data to SOL E-field
    for (dit.begin(); dit.ok(); ++dit) {
        int block_number( whichBlock( grids[dit] ) );
        if (block_number==2 || block_number==3) {
            for (int dir = 0; dir < 2; dir++) {
                
                FArrayBox& this_Er_face_dir = a_Er_average_face[dit][dir];
                
                Box box( this_Er_face_dir.box() );
                BoxIterator bit(box);
                for (bit.begin(); bit.ok(); ++bit) {
                    
                    IntVect iv = bit();
                    int npol;
                    if (block_number==2) {
                        if (iv[1]<lo_pol_LCORE) { npol = 0; }
                        else if (iv[1]>lo_pol_LCORE+size-1) { npol = size-1; }
                        else { npol = iv[1]-lo_pol_LCORE; }
                        
                        this_Er_face_dir(iv,0) = ER_LCORE_bnd[npol];
                        this_Er_face_dir(iv,2) = EZ_LCORE_bnd[npol];
                        
                    }
                    
                    if (block_number==3) {
                        if (iv[1]<lo_pol_RCORE) { npol = 0; }
                        else if (iv[1]>lo_pol_RCORE+size-1) { npol = size-1; }
                        else { npol = iv[1]-lo_pol_RCORE; }
                        
                        this_Er_face_dir(iv,0) = ER_RCORE_bnd[npol];
                        this_Er_face_dir(iv,2) = EZ_RCORE_bnd[npol];
                        
                    }
                    
                }
            }
            
            FArrayBox& this_Er_cell = a_Er_average_cell[dit];
            
            Box cellbox( this_Er_cell.box() );
            BoxIterator bit(cellbox);
            for (bit.begin(); bit.ok(); ++bit) {
                
                IntVect iv = bit();
                int npol;
                if (block_number==2) {
                    if (iv[1]<lo_pol_LCORE) { npol = 0; }
                    else if (iv[1]>lo_pol_LCORE+size-1) { npol = size-1; }
                    else { npol = iv[1]-lo_pol_LCORE; }
                    
                    
                    this_Er_cell(iv,0) = ER_LCORE_bnd[npol];
                    this_Er_cell(iv,2) = EZ_LCORE_bnd[npol];
                }
                
                if (block_number==3) {
                    if (iv[1]<lo_pol_RCORE) { npol = 0; }
                    else if (iv[1]>lo_pol_RCORE+size-1) { npol = size-1; }
                    else { npol = iv[1]-lo_pol_RCORE; }
                    
                    this_Er_cell(iv,0) = ER_RCORE_bnd[npol];
                    this_Er_cell(iv,2) = EZ_RCORE_bnd[npol];
                }
                
            }
            
        }
    }
    
    
    delete [] ER_LCORE_bnd_loc;
    delete [] ER_RCORE_bnd_loc;
    delete [] EZ_LCORE_bnd_loc;
    delete [] EZ_RCORE_bnd_loc;
    
    delete [] ER_LCORE_bnd;
    delete [] ER_RCORE_bnd;
    delete [] EZ_LCORE_bnd;
    delete [] EZ_RCORE_bnd;
    
    
}

void SingleNullCoordSys::checkGridConsistency() const
{
   double epsilon = 1.0e-9; //Here, just an arbitrary very small number
   
   const SingleNullBlockCoordSysModel& lcore_coord_sys = (const SingleNullBlockCoordSysModel&)(*(getCoordSys(SingleNullBlockCoordSys::LCORE)));
   const SingleNullBlockCoordSysModel& lsol_coord_sys = (const SingleNullBlockCoordSysModel&)(*(getCoordSys(SingleNullBlockCoordSys::LSOL)));
   const SingleNullBlockCoordSysModel& lpf_coord_sys = (const SingleNullBlockCoordSysModel&)(*(getCoordSys(SingleNullBlockCoordSys::LPF)));
   const SingleNullBlockCoordSysModel& rpf_coord_sys = (const SingleNullBlockCoordSysModel&)(*(getCoordSys(SingleNullBlockCoordSys::RPF)));
   
   int nrad_CORE = lcore_coord_sys.domain().domainBox().size(RADIAL_DIR);
   int npol_CORE = lcore_coord_sys.domain().domainBox().size(POLOIDAL_DIR);
   int nrad_SOL = lsol_coord_sys.domain().domainBox().size(RADIAL_DIR);
   int nrad_PF = lpf_coord_sys.domain().domainBox().size(RADIAL_DIR);
   int npol_LPF = lpf_coord_sys.domain().domainBox().size(POLOIDAL_DIR);
   int npol_RPF = rpf_coord_sys.domain().domainBox().size(POLOIDAL_DIR);
   
   if (fabs(lcore_coord_sys.getCoreWidth()/(double)nrad_CORE - lcore_coord_sys.getSolWidth()/(double)nrad_SOL)>epsilon ) {
      MayDay::Error("SingleNullCoordSysModel::checkGridConsistency(): model geometry requires the same radial cell size in the SOL and CORE blocks");
   }

   if (fabs(lcore_coord_sys.getCoreWidth()/(double)nrad_CORE - lcore_coord_sys.getPFWidth()/(double)nrad_PF)>epsilon ) {
      MayDay::Error("SingleNullCoordSysModel::checkGridConsistency(): model geometry requires the same radial cell size in the PF and CORE blocks");
   }

   if (fabs(lcore_coord_sys.getCoreArcLength()/(double)npol_CORE - lcore_coord_sys.getPFArcLength()/(double)npol_LPF)>epsilon ) {
      MayDay::Error("SingleNullCoordSysModel::checkGridConsistency(): model geometry requires the same poloidal cell size in the LPF and CORE blocks");
   }

   if (fabs(lcore_coord_sys.getCoreArcLength()/(double)npol_CORE - lcore_coord_sys.getPFArcLength()/(double)npol_RPF)>epsilon ) {
      MayDay::Error("SingleNullCoordSysModel::checkGridConsistency(): model geometry requires the same poloidal cell size in the RPF and CORE blocks");
   }

}


#if CFG_DIM == 3
void
SingleNullCoordSys::toroidalBlockRemapping( IntVect&               a_ivDst,
                                            int&                   a_nDst,
                                            Vector<Real>&          a_interpStecil,
                                            Vector<int>&           a_interpStencilOffsets,
                                            const RealVect&        a_xiSrc,
                                            const int              a_nSrc,
                                            const Side::LoHiSide&  a_side ) const
{
   CH_TIME("SingleNullCoordSys::toroidalBlockRemapping");
   
   int src_toroidal_sector = toroidalBlockNumber(a_nSrc);
   
   // Determine the toroidal sector containing the input point, which is assumed to
   // lie in a ghost cell on the toroidal a_side of src_toroidal_sector.  Periodicity
   // is assumed in the toroidal direction.
   int dst_toroidal_sector;
   if (a_side == Side::LoHiSide::Lo) dst_toroidal_sector = src_toroidal_sector - 1;
   if (a_side == Side::LoHiSide::Hi) dst_toroidal_sector = src_toroidal_sector + 1;
   if ( dst_toroidal_sector == m_num_toroidal_sectors ) {
     dst_toroidal_sector = 0;
   }
   if ( dst_toroidal_sector == -1 ) {
     dst_toroidal_sector = m_num_toroidal_sectors - 1;
   }
   
   const SingleNullBlockCoordSys* src_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(a_nSrc);
   double delta_Phi_lo = src_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR) - a_xiSrc[TOROIDAL_DIR];
   double delta_Phi_hi = a_xiSrc[TOROIDAL_DIR] - src_coord_sys->upperMappedCoordinate(TOROIDAL_DIR);

   RealVect X = src_coord_sys->realCoord(a_xiSrc);

   // Get the poloidal (R,Z) coordinates of the input point
   POL::RealVect X_src_pol = src_coord_sys->restrictPhysCoordToPoloidal(X);

   // Determine the poloidal block of dst_toroidal_sector containing the point.  Determine the
   // corresponding mapped coordinate with respect to that block (which only differs from the
   // src mapped coordinate in the toroidal component).
   
   RealVect X_dst_cent;    //Physical coord of the dst cell center
   RealVect X_dst_pol_hi;  //Physical coord of the dst cell hi-poloidal face center
   RealVect X_dst_pol_lo;  //Physical coord of the dst cell lo-poloidal face center
   RealVect xi_dst;        //Mapped coord of the dst point
   RealVect xi_dst_cent;   //Mapped coord of the dst cell center
   
   // Assume flux-aligned mapping
   xi_dst[RADIAL_DIR] = a_xiSrc[RADIAL_DIR];
   
   double R_X = sqrt(X[0]*X[0] + X[1]*X[1]);
   double Z_X = X[2];

   // Figure out in which poloidal block a_xiSrc is valid by determining which poloidal plane it belongs
   // to and then doing a brute force search of that plane for the closest cell center

   double d_min = DBL_MAX;
   int dst_poloidal_block;
   int i_tor_saved;
     
   for (int poloidal_block=0; poloidal_block<m_num_poloidal_blocks; ++poloidal_block) {

      int src_poloidal_block = poloidalBlockNumber(a_nSrc);
	
      if (sameMagFluxSurf(src_poloidal_block, poloidal_block)) {

         int block_number = blockNumber(poloidal_block, dst_toroidal_sector);
         const SingleNullBlockCoordSys* coord_sys = (SingleNullBlockCoordSys*)getCoordSys(block_number);
         const RealVect& dx = coord_sys->getMappedCellSize();

         int i_tor;
         if (a_side == Side::LoHiSide::Lo) {
            xi_dst[TOROIDAL_DIR] = coord_sys->upperMappedCoordinate(TOROIDAL_DIR) - delta_Phi_lo;
            i_tor = floor( (coord_sys->upperMappedCoordinate(TOROIDAL_DIR) - delta_Phi_lo) / dx[TOROIDAL_DIR] );
         }
         else if (a_side == Side::LoHiSide::Hi) {
            xi_dst[TOROIDAL_DIR] = coord_sys->lowerMappedCoordinate(TOROIDAL_DIR) + delta_Phi_hi;
            i_tor = floor( (coord_sys->lowerMappedCoordinate(TOROIDAL_DIR) + delta_Phi_hi) / dx[TOROIDAL_DIR] );
         }

         const Box& domain_box = coord_sys->domain().domainBox();
         IntVect lo = domain_box.smallEnd();
         IntVect hi = domain_box.bigEnd();
         lo[TOROIDAL_DIR] = hi[TOROIDAL_DIR] = i_tor;
         lo[RADIAL_DIR] = hi[RADIAL_DIR] = floor(a_xiSrc[RADIAL_DIR] / dx[RADIAL_DIR]);
	 
         Box search_box(lo,hi);
         CH_assert(search_box.ok());
         FArrayBox phys_coords(search_box,SpaceDim);
         coord_sys->getCellCenteredRealCoords(phys_coords);

         for (BoxIterator bit(search_box); bit.ok(); ++bit) {
            IntVect iv = bit();
	   
            RealVect X_phys;
            for (int n=0; n<SpaceDim; ++n) {
               X_phys[n] = phys_coords(iv,n);
            }
            double this_R_diff = sqrt(X_phys[0]*X_phys[0] + X_phys[1]*X_phys[1]) - R_X;
            double this_Z_diff = X_phys[2] - Z_X;
            double d = sqrt(this_R_diff*this_R_diff + this_Z_diff*this_Z_diff);
            if ( d < d_min ) {
               d_min = d;
               i_tor_saved = i_tor;
               dst_poloidal_block = poloidal_block;
               X_dst_cent = X_phys;

               xi_dst_cent = xi_dst;
               xi_dst_cent[POLOIDAL_DIR] = (iv[POLOIDAL_DIR] + 0.5)*dx[POLOIDAL_DIR];
	     
               RealVect xi_dst_pol_hi(xi_dst);
               xi_dst_pol_hi[POLOIDAL_DIR] = (iv[POLOIDAL_DIR] + 1.0)*dx[POLOIDAL_DIR];
               X_dst_pol_hi = coord_sys -> realCoord(xi_dst_pol_hi);

               RealVect xi_dst_pol_lo(xi_dst);
               xi_dst_pol_lo[POLOIDAL_DIR] = (iv[POLOIDAL_DIR] )*dx[POLOIDAL_DIR];
               X_dst_pol_lo = coord_sys -> realCoord(xi_dst_pol_lo);
            }
         }
      }
   }
 
   a_nDst = blockNumber(dst_poloidal_block, dst_toroidal_sector);
   const SingleNullBlockCoordSys* trial_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(a_nDst);
   const RealVect& trial_dx = trial_coord_sys->getMappedCellSize();

   // Get the poloidal mapped coordinate
   applyToroidalPeriodicity(X, a_side);
   xi_dst[POLOIDAL_DIR] = getInterpMappedPoloidalCoord(X, X_dst_cent, X_dst_pol_hi, X_dst_pol_lo, xi_dst_cent, trial_dx);
   
   // Get the global index of the dst cell center
   for (int dir=0; dir<SpaceDim; ++dir) {
     a_ivDst[dir] = floor(xi_dst[dir]/trial_dx[dir]);
   }

   IntVect iv_saved = a_ivDst;
   RealVect saved_X_dst_cent = X_dst_cent;
   
   if ( fixPoloidalBlock( a_nDst, a_ivDst ) ) {

      const SingleNullBlockCoordSys* coord_sys = (SingleNullBlockCoordSys*)getCoordSys(a_nDst);
      const RealVect& dx = coord_sys->getMappedCellSize();
      const Box& domain_box = coord_sys->domain().domainBox();
      IntVect lo = domain_box.smallEnd();
      IntVect hi = domain_box.bigEnd();
      lo[TOROIDAL_DIR] = hi[TOROIDAL_DIR] = i_tor_saved;
      lo[RADIAL_DIR] = hi[RADIAL_DIR] = floor(a_xiSrc[RADIAL_DIR] / dx[RADIAL_DIR]);
	 
      Box search_box(lo,hi);
      CH_assert(search_box.ok());
      FArrayBox phys_coords(search_box,SpaceDim);
      coord_sys->getCellCenteredRealCoords(phys_coords);

      RealVect X_phys;
      for (int n=0; n<SpaceDim; ++n) {
         X_phys[n] = phys_coords(a_ivDst,n);
      }

      X_dst_cent = X_phys;

      xi_dst_cent = xi_dst;
      xi_dst_cent[POLOIDAL_DIR] = (a_ivDst[POLOIDAL_DIR] + 0.5)*dx[POLOIDAL_DIR];
	     
      RealVect xi_dst_pol_hi(xi_dst);
      xi_dst_pol_hi[POLOIDAL_DIR] = (a_ivDst[POLOIDAL_DIR] + 1.0)*dx[POLOIDAL_DIR];
      X_dst_pol_hi = coord_sys -> realCoord(xi_dst_pol_hi);

      RealVect xi_dst_pol_lo(xi_dst);
      xi_dst_pol_lo[POLOIDAL_DIR] = (a_ivDst[POLOIDAL_DIR] )*dx[POLOIDAL_DIR];
      X_dst_pol_lo = coord_sys -> realCoord(xi_dst_pol_lo);

      // Get the poloidal mapped coordinate
      //      applyToroidalPeriodicity(X);
      xi_dst[POLOIDAL_DIR] = getInterpMappedPoloidalCoord(X, X_dst_cent, X_dst_pol_hi, X_dst_pol_lo, xi_dst_cent, trial_dx);
   
#if 0
      // Get the global index of the dst cell center
      IntVect iv_test;
      for (int dir=0; dir<SpaceDim; ++dir) {
         iv_test[dir] = floor(xi_dst[dir]/trial_dx[dir]);
      }
      
      if ( iv_test != a_ivDst ) {
         MayDay::Error("SingleNullCoordSys::toroidalBlockRemapping(): fixPoloidalBlock didn't work, need to investigate");
      }

#endif      
   }
   
   const SingleNullBlockCoordSys* dst_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(a_nDst);
   const RealVect& dx_dst = dst_coord_sys->getMappedCellSize();
   
   // If this warning comes out it means that the dst point is not really
   // within the cell we though it should be (from the closest center search)
   // This can potentially happen if the dst point is very near poloidal cell
   // boundary (and can potentially be in a different block, if near a block bndry).
   // All this is perhaps OK, since we should have sufficient number of ghost cell
   // in mapping and data available.
#if 0
   IntVect iv_test;
   for (int dir=0; dir<SpaceDim; ++dir) {
      iv_test[dir] = floor(xi_dst_cent[dir]/dx_dst[dir]);
   }
   if ( procID()==0 && a_ivDst != iv_test ) {
      MayDay::Warning("SingleNullCoordSys::toroidalBlockRemapping: closest center search failed");
   }
#endif

   /*
    Get interpolation coefficients (interpolation stencil) and offsets
    The last argument of interpolation_stencil controls saw tooth BC case
    (a) value is less than zero -- treat that cell as an interior cell;
    however still shift an interpolation stencil to avoid the info from poloildal BCs
    (b) Value is between one and zero -- transition cell that is not entirely in the tooth
    (c) Value is larger than one -- cell is entirely inside the tooth
   
    For now we assume that the poloidal size of the cells at the block interafce is not
    more than a factor of 2 different. Otherwise, there may be a situation, where
    a_ivDst is far from the domain boundary, but the lower/upper corner of that cell
    belongs to saw tooth.
   */

   a_interpStencilOffsets[0] = -1.;
   a_interpStencilOffsets[1] =  0.;
   a_interpStencilOffsets[2] =  1.;

   // Set default value to a negative (designates inner cells)
   int stenc_size = a_interpStecil.size();
   a_interpStecil[stenc_size-1] = -1.0;
   
   //Handle the saw-tooth boundary stencil coefficeints
   Box dst_domain_box = dst_coord_sys->domain().domainBox();
   
   if ( a_ivDst[POLOIDAL_DIR] <= dst_domain_box.smallEnd(POLOIDAL_DIR) &&
        containsPhysicalBoundary(a_nDst, POLOIDAL_DIR, Side::LoHiSide::Lo) ) {
      int shift = dst_domain_box.smallEnd(POLOIDAL_DIR) - a_ivDst[POLOIDAL_DIR] + 1;
      for (int n=0; n<3; ++n) {
         a_interpStencilOffsets[n] += shift;
      }

      const RealVect& dx_src = src_coord_sys->getMappedCellSize();

      RealVect xiSrc_lo_pol(a_xiSrc);
      xiSrc_lo_pol[POLOIDAL_DIR] -= 0.5 * dx_src[POLOIDAL_DIR];
      RealVect X_lo_pol = src_coord_sys->realCoord(xiSrc_lo_pol);
      applyToroidalPeriodicity(X_lo_pol, a_side);

      RealVect xiSrc_hi_pol(a_xiSrc);
      xiSrc_hi_pol[POLOIDAL_DIR] += 0.5 * dx_src[POLOIDAL_DIR];
      RealVect X_hi_pol = src_coord_sys->realCoord(xiSrc_hi_pol);
      applyToroidalPeriodicity(X_hi_pol, a_side);

      RealVect xi_dst_bnd(xi_dst);
      xi_dst_bnd[POLOIDAL_DIR] = dst_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);
      RealVect X_dst_bnd = dst_coord_sys->realCoord(xi_dst_bnd);

      RealVect X_hi_to_lo(X_hi_pol);
      X_hi_to_lo -= X_lo_pol;
      RealVect X_bnd_to_low(X_dst_bnd);
      X_bnd_to_low -= X_lo_pol;

      Real tooth_height_fac = X_hi_to_lo.dotProduct(X_bnd_to_low);
      int stenc_size = a_interpStecil.size();
      a_interpStecil[stenc_size-1] = tooth_height_fac / X_hi_to_lo.radSquared();
   }

   if ( a_ivDst[POLOIDAL_DIR] >= dst_domain_box.bigEnd(POLOIDAL_DIR) &&
        containsPhysicalBoundary(a_nDst, POLOIDAL_DIR, Side::LoHiSide::Hi) ) {
      int shift = a_ivDst[POLOIDAL_DIR] - dst_domain_box.bigEnd(POLOIDAL_DIR) + 1;
      for (int n=0; n<3; ++n) {
         a_interpStencilOffsets[n] -= shift;
      }
      
      const RealVect& dx_src = src_coord_sys->getMappedCellSize();

      RealVect xiSrc_lo_pol(a_xiSrc);
      xiSrc_lo_pol[POLOIDAL_DIR] -= 0.5 * dx_src[POLOIDAL_DIR];
      RealVect X_lo_pol = src_coord_sys->realCoord(xiSrc_lo_pol);
      applyToroidalPeriodicity(X_lo_pol, a_side);

      RealVect xiSrc_hi_pol(a_xiSrc);
      xiSrc_hi_pol[POLOIDAL_DIR] += 0.5 * dx_src[POLOIDAL_DIR];
      RealVect X_hi_pol = src_coord_sys->realCoord(xiSrc_hi_pol);
      applyToroidalPeriodicity(X_hi_pol, a_side);

      RealVect xi_dst_bnd(xi_dst);
      xi_dst_bnd[POLOIDAL_DIR] = dst_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);
      RealVect X_dst_bnd = dst_coord_sys->realCoord(xi_dst_bnd);

      RealVect X_hi_to_lo(X_hi_pol);
      X_hi_to_lo -= X_lo_pol;
      RealVect X_hi_to_bnd(X_hi_pol);
      X_hi_to_bnd -= X_dst_bnd;
     
      Real tooth_height_fac = X_hi_to_lo.dotProduct(X_hi_to_bnd);
      int stenc_size = a_interpStecil.size();
      a_interpStecil[stenc_size-1] = tooth_height_fac / X_hi_to_lo.radSquared();
   }

   // We use currently three cells to interpolate;
   // the interpolation error is h^3.
   int order = 3;

   if (m_interpolation_scheme == "mapped") {
     getInterpolationCoefficients(a_interpStecil, a_interpStencilOffsets, xi_dst, a_ivDst, dx_dst, order);
   }
   else if (m_interpolation_scheme == "physical") {
     getInterpolationCoefficients(a_interpStecil, a_interpStencilOffsets, X, a_ivDst, *dst_coord_sys);
   }
   else {
     MayDay::Error("SingleNullCoordSys::toroidalBlockRemapping(): Invalid interpolation scheme type");
   }

#if 0
   if ( a_ivDst[POLOIDAL_DIR] <= dst_domain_box.smallEnd(POLOIDAL_DIR) &&
        containsPhysicalBoundary(a_nDst, POLOIDAL_DIR, Side::LoHiSide::Lo) ||
        a_ivDst[POLOIDAL_DIR] >= dst_domain_box.bigEnd(POLOIDAL_DIR) &&
        containsPhysicalBoundary(a_nDst, POLOIDAL_DIR, Side::LoHiSide::Hi) ) {
      for (int n=0; n<a_interpStecil.size(); ++n) a_interpStecil[n] = 0.;
   }
#endif

#if 0
   if ( a_ivDst[POLOIDAL_DIR] <= dst_domain_box.smallEnd(POLOIDAL_DIR) &&
        containsPhysicalBoundary(a_nDst, POLOIDAL_DIR, Side::LoHiSide::Lo) ) {
      cout << a_ivDst[POLOIDAL_DIR] << " " << dst_domain_box.smallEnd(POLOIDAL_DIR) << " " << a_interpStencilOffsets[0] << " " << a_interpStencilOffsets[1] << " " << a_interpStencilOffsets[2] << endl;
   }
#endif

#if 0
   if ( a_ivDst[POLOIDAL_DIR] >= dst_domain_box.bigEnd(POLOIDAL_DIR) &&
        containsPhysicalBoundary(a_nDst, POLOIDAL_DIR, Side::LoHiSide::Hi) ) {
      cout << a_ivDst[POLOIDAL_DIR] << " " << dst_domain_box.bigEnd(POLOIDAL_DIR) << " " << a_interpStencilOffsets[0] << " " << a_interpStencilOffsets[1] << " " << a_interpStencilOffsets[2] << endl;
   }
#endif

}


void
SingleNullCoordSys::getInterpolationCoefficients( Vector<Real>&    a_coeff,
                                                  Vector<int>&     a_offsets,
                                                  const RealVect&  a_xi0_dst,
                                                  const IntVect&   a_iv0_dst,
                                                  const RealVect&  a_dx_dst,
                                                  const int        a_order ) const

{
   Real h = a_dx_dst[POLOIDAL_DIR];
   
   Real cent = (a_iv0_dst[POLOIDAL_DIR] + a_offsets[1] + 0.5) * h;
   Real lo = cent - h;
   Real hi = cent + h;

   Real coeffLo(0.0);
   Real coeffCent(1.0);
   Real coeffHi(0.0);
   
   Real point = a_xi0_dst[POLOIDAL_DIR];

   CH_assert(a_order == 3);

   if (a_order == 3) {

      Real d0 = lo - point;
      Real d1 = cent - point;
      Real d2 = hi - point;
   
      coeffLo = d1*d2/(d0-d1)/(d0-d2);
      coeffCent = d2*d0/(d1-d2)/(d1-d0);
      coeffHi = d0*d1/(d2-d0)/(d2-d1);
   }
   
   else if (a_order == 2) {
      if (point > cent) {
         coeffLo = 0.0;
         coeffCent = 1.0 - (point - cent)/h;
         coeffHi = (point - cent)/h;
      }
      else {
         coeffLo = (cent - point)/h;
         coeffCent =  1.0 - (cent - point)/h;
         coeffHi = 0.0;
      }
   }
   
   else {
      MayDay::Error("getInterpolationCoefficients:: only order = 2 or 3 is currently supported");
   }
   
   a_coeff[0] = coeffLo;
   a_coeff[1] = coeffCent;
   a_coeff[2] = coeffHi;
}


void
SingleNullCoordSys::getInterpolationCoefficients( Vector<Real>&                  a_coeff,
                                                  Vector<int>&                   a_offsets,
                                                  const RealVect&                a_X,
                                                  const IntVect&                 a_iv0_dst,
                                                  const SingleNullBlockCoordSys& a_coords_dst) const

{
   //Get physical coordinate of the dst cell-center
   RealVect dx = a_coords_dst.dx();
   RealVect xi_dest_cent;
   for (int dir=0; dir<SpaceDim; ++dir) {
      if (dir == POLOIDAL_DIR) {
         xi_dest_cent[dir] = (a_iv0_dst[dir] + a_offsets[1] + 0.5) * dx[dir];
      }
      else {
         xi_dest_cent[dir] = (a_iv0_dst[dir] + 0.5) * dx[dir];
      }
   }
   RealVect X_cent = a_coords_dst.realCoord(xi_dest_cent);
   Real R_cent = sqrt(X_cent[0]*X_cent[0] + X_cent[1]*X_cent[1]);
   Real Z_cent = X_cent[2];

   //Get physical coordinate of the hi cell-center
   RealVect xi_dest_hi(xi_dest_cent);
   xi_dest_hi[POLOIDAL_DIR] += dx[POLOIDAL_DIR];
   RealVect X_hi = a_coords_dst.realCoord(xi_dest_hi);
   Real R_hi = sqrt(X_hi[0]*X_hi[0] + X_hi[1]*X_hi[1]);
   Real Z_hi = X_hi[2];

   //Get physical coordinate of the lo cell-center
   RealVect xi_dest_lo(xi_dest_cent);
   xi_dest_lo[POLOIDAL_DIR] -= dx[POLOIDAL_DIR];
   RealVect X_lo = a_coords_dst.realCoord(xi_dest_lo);
   Real R_lo = sqrt(X_lo[0]*X_lo[0] + X_lo[1]*X_lo[1]);
   Real Z_lo = X_lo[2];
   
   //Get (R,Z) coordinate of the destination point
   Real R_point = sqrt(a_X[0]*a_X[0] + a_X[1]*a_X[1]);
   Real Z_point = a_X[2];

   Real d0;
   Real d1;
   Real d2;
   
   // Check if radial points are monotonically increasing/decreasing
   if ((R_hi - R_cent)*(R_cent - R_lo) > 0 && abs(R_hi-R_lo)>abs(Z_hi-Z_lo)) {
      d0 = R_lo - R_point;
      d1 = R_cent - R_point;
      d2 = R_hi - R_point;
   }
   // else use Z points
   else if ( (Z_hi - Z_cent)*(Z_cent - Z_lo) > 0 && abs(Z_hi-Z_lo)>abs(R_hi-R_lo)) {
      d0 = Z_lo - Z_point;
      d1 = Z_cent - Z_point;
      d2 = Z_hi - Z_point;
   }
   
   else {
      MayDay::Error("getInterpolationCoefficients:: flux-surface appears to be too curved");
   }
   
   //Get inerpolation coefficients
   Real coeffLo = d1*d2/(d0-d1)/(d0-d2);
   Real coeffCent = d2*d0/(d1-d2)/(d1-d0);
   Real coeffHi = d0*d1/(d2-d0)/(d2-d1);

   //Assign interpolation coefficients
   a_coeff[0] = coeffLo;
   a_coeff[1] = coeffCent;
   a_coeff[2] = coeffHi;
}

void
SingleNullCoordSys::applyToroidalPeriodicity(RealVect& a_x,
					     const Side::LoHiSide&  a_side) const
{

  double phi_max = m_toroidal_width;
  double phi = atan2(a_x[1], a_x[0]);
  double r = sqrt(pow(a_x[0],2) + pow(a_x[1],2));

  // Correct for atan2 limits, since valid
  // toroidal angle runs from 0 to phi_max
  if (a_side == Side::LoHiSide::Hi && phi < 0.) {
    phi += 2. * Pi;
  }

  if (phi < 0.0) {
    double phi_reflected = phi + phi_max;
    a_x[0] = r * cos(phi_reflected);
    a_x[1] = r * sin(phi_reflected);
  }

  if (phi > phi_max) {
    double phi_reflected = phi - phi_max;
    a_x[0] = r * cos(phi_reflected);
    a_x[1] = r * sin(phi_reflected);
  }
}


Real
SingleNullCoordSys::getInterpMappedPoloidalCoord(const RealVect& a_X,
						 const RealVect& a_X_dst_pol_cent,
						 const RealVect& a_X_dst_pol_hi,
						 const RealVect& a_X_dst_pol_lo,
						 const RealVect& a_xi_dst_cent,
						 const RealVect& a_dx) const
{
  /*
    This function assumes that all vectors lie in the same toroidal plane
    Further we assume that mapping is accurate enough to assure that all vectors
    lie on the same flux surface (i.e., radial coordinate surface)
   */

  Real result;
  
  RealVect tmp(a_X);
  tmp -= a_X_dst_pol_hi;
  Real dst_to_hi_end = tmp.vectorLength();
  tmp = a_X;
  tmp -= a_X_dst_pol_lo;
  Real dst_to_lo_end = tmp.vectorLength();

  RealVect X_to_Xcent(a_X);
  X_to_Xcent -= a_X_dst_pol_cent;
  
  if (dst_to_hi_end < dst_to_lo_end) {

    RealVect Xhi_to_Xcent(a_X_dst_pol_hi);
    Xhi_to_Xcent -= a_X_dst_pol_cent;

    result = Xhi_to_Xcent.dotProduct(X_to_Xcent);
    result /= Xhi_to_Xcent.radSquared();
    result *= 0.5 * a_dx[POLOIDAL_DIR];
    result += a_xi_dst_cent[POLOIDAL_DIR];
    
  }
  else {

    RealVect Xlo_to_Xcent(a_X_dst_pol_lo);
    Xlo_to_Xcent -= a_X_dst_pol_cent;

    result = Xlo_to_Xcent.dotProduct(X_to_Xcent);
    result /= Xlo_to_Xcent.radSquared();
    result *= -0.5 * a_dx[POLOIDAL_DIR];
    result += a_xi_dst_cent[POLOIDAL_DIR];

  }

  return result;
}

#endif


Vector<RealVect>
SingleNullCoordSys::displacements(const Vector<RealVect>&   a_dstCoords,
                                  const Vector<int>&        a_dstBlocks,
                                  const RealVect&           a_srcCoords,
                                  int                       a_srcBlock) const
{
#if CFG_DIM==2
   return MultiBlockCoordSys::displacements(a_dstCoords, a_dstBlocks, a_srcCoords, a_srcBlock);
#endif
#if CFG_DIM==3
   Vector<RealVect> dstCoords(a_dstCoords.size());

   for (int i=0; i<a_dstCoords.size(); ++i) {
      dstCoords[i] = a_dstCoords[i];

      double toroidal_mapped_coord = a_dstCoords[i][TOROIDAL_DIR];

      if ( toroidal_mapped_coord < a_srcCoords[TOROIDAL_DIR] ) {
         double shifted_coord = toroidal_mapped_coord + m_toroidal_width;
         if ( fabs(shifted_coord - a_srcCoords[TOROIDAL_DIR]) < fabs(toroidal_mapped_coord - a_srcCoords[TOROIDAL_DIR]) ) {
            dstCoords[i][TOROIDAL_DIR] = shifted_coord;
         }
      }
      else if ( toroidal_mapped_coord > a_srcCoords[TOROIDAL_DIR] ) {
         double shifted_coord = toroidal_mapped_coord - m_toroidal_width;
         if ( fabs(shifted_coord - a_srcCoords[TOROIDAL_DIR]) < fabs(toroidal_mapped_coord - a_srcCoords[TOROIDAL_DIR]) ) {
            dstCoords[i][TOROIDAL_DIR] = shifted_coord;
         }
      }
   }

   return MultiBlockCoordSys::displacements(dstCoords, a_dstBlocks, a_srcCoords, a_srcBlock);
#endif
}

bool
SingleNullCoordSys::sameMagFluxSurf(const int& a_pol_block_src,
                                    int&       a_pol_block) const
{
  bool result(false);
  switch ( a_pol_block_src )
    {
    case SingleNullBlockCoordSys::LCORE:
      if (a_pol_block == SingleNullBlockCoordSys::LCORE ||
          a_pol_block == SingleNullBlockCoordSys::RCORE) {
         result = true;
      }
      break;

    case SingleNullBlockCoordSys::RCORE:
      if (a_pol_block == SingleNullBlockCoordSys::RCORE ||
          a_pol_block == SingleNullBlockCoordSys::LCORE) {
         result = true;
      }
      break;

    case SingleNullBlockCoordSys::LCSOL:
      if (a_pol_block == SingleNullBlockCoordSys::LCSOL ||
          a_pol_block == SingleNullBlockCoordSys::RCSOL ||
          a_pol_block == SingleNullBlockCoordSys::LSOL  ||
          a_pol_block == SingleNullBlockCoordSys::RSOL ) {
         result = true;
      }
      break;

    case SingleNullBlockCoordSys::RCSOL:
      if (a_pol_block == SingleNullBlockCoordSys::RCSOL ||
          a_pol_block == SingleNullBlockCoordSys::LCSOL ||
          a_pol_block == SingleNullBlockCoordSys::RSOL  ||
          a_pol_block == SingleNullBlockCoordSys::LSOL ) {
         result = true;
      }
      break;

    case SingleNullBlockCoordSys::LSOL:
      if (a_pol_block == SingleNullBlockCoordSys::LSOL ||
          a_pol_block == SingleNullBlockCoordSys::RCSOL ||
          a_pol_block == SingleNullBlockCoordSys::LCSOL  ||
          a_pol_block == SingleNullBlockCoordSys::RSOL ) {
         result = true;
      }
      break;
         
    case SingleNullBlockCoordSys::RSOL:
      if (a_pol_block == SingleNullBlockCoordSys::RSOL ||
          a_pol_block == SingleNullBlockCoordSys::LCSOL ||
          a_pol_block == SingleNullBlockCoordSys::RCSOL  ||
          a_pol_block == SingleNullBlockCoordSys::LSOL ) {
         result = true;
      }
      break;

    case SingleNullBlockCoordSys::LPF:
      if (a_pol_block == SingleNullBlockCoordSys::LPF ||
          a_pol_block == SingleNullBlockCoordSys::RPF ) {
         result = true;
      }
      break;

    case SingleNullBlockCoordSys::RPF:
      if (a_pol_block == SingleNullBlockCoordSys::RPF ||
          a_pol_block == SingleNullBlockCoordSys::LPF ) {
         result = true;
      }
      break;
    }
   
  return result;
}


bool
SingleNullCoordSys::containsPhysicalBoundary( int                    a_block_number,
                                              int                    a_dir,
                                              const Side::LoHiSide&  a_side ) const
{
   bool contains_boundary = false;

#if CFG_DIM==2
   contains_boundary = MagCoordSys::containsPhysicalBoundary(a_block_number, a_dir, a_side);
#endif

#if CFG_DIM==3
   if ( a_dir == TOROIDAL_DIR ) {

      if ( m_sheared_geometry ) {
         // Set all toroidal block interfaces to physical boundaries;
         // This is done to deal with saw-tooth divertor BCs
         // Internal part of the block interface is handled by fillInternalGhosts
         contains_boundary = true;
      }
      else {
         contains_boundary = MagCoordSys::containsPhysicalBoundary(a_block_number, a_dir, a_side);
      }
   }
   else {
      contains_boundary = MagCoordSys::containsPhysicalBoundary(a_block_number, a_dir, a_side);
   }
#endif
   
   return contains_boundary;
}


#if CFG_DIM==3
Vector<Box>
SingleNullCoordSys::getToroidalBoundaryBoxes( const int              a_block_number,
                                              const Side::LoHiSide&  a_side,
                                              const int              a_shear_interp_order ) const
{
   CH_TIME("SingleNullCoordSys::getToroidalBoundaryBoxes");

   const MagBlockCoordSys& block_coord_sys = *(getCoordSys(a_block_number));
   const ProblemDomain& domain = block_coord_sys.domain();
   const Box& domain_box = domain.domainBox();

   int shear_dir = TOROIDAL_DIR;

   Box bndryBox = adjCellBox(domain_box, shear_dir, a_side, 1);

   FArrayBox mappedBndryCoord(bndryBox,SpaceDim);
   block_coord_sys.getCellCenteredMappedCoords(mappedBndryCoord);
         
   RealVect dx = block_coord_sys.dx();

   IntVectSet bndry_cells;

   for ( BoxIterator bit(bndryBox); bit.ok(); ++bit ) {
      const IntVect& iv = bit();

      RealVect xiSrc;
      for (int dir=0; dir<SpaceDim; ++dir) {
         xiSrc[dir] = mappedBndryCoord(iv,dir);
      }
            
      IntVect iv_src;
      for (int dir=0; dir<SpaceDim; ++dir) {
         iv_src[dir] = floor(xiSrc[dir]/dx[dir]);
      }

      int nDst;
      IntVect ivDst;
      Vector<Real> interpStencil(a_shear_interp_order + 2, 0);
      Vector<int> interpStencilOffsets(a_shear_interp_order + 1, 0);

      toroidalBlockRemapping(ivDst, nDst, interpStencil, interpStencilOffsets, xiSrc, a_block_number, a_side);

      if ( interpStencil[a_shear_interp_order + 1] >= 0.5 ) {
         bndry_cells |= iv;
      }
   }

   // Convert the IntVectSet into a vector of boxes
   Vector<Box> box_vector = bndry_cells.boxes();

   // Coalesce the IntVectSet into boxes

   Vector<Box> coalesced_box_vector;

   if ( box_vector.size() > 0 ) {
      Box accum_box = box_vector[0];

      for (int n=1; n<box_vector.size(); ++n) {
         Box new_box = accum_box;
         new_box.minBox(box_vector[n]);

         IntVectSet box_ivs(new_box);

         // Check that the new box is full, i.e., bndry_cells
         // contains all of the new box's IntVects
         if ( bndry_cells.contains(box_ivs) ) {
            // Keep accumulating and pushback if no more boxes
            accum_box = new_box;
            if ( n == box_vector.size()-1 ) {
               coalesced_box_vector.push_back(bdryBox(accum_box,shear_dir,a_side,1));
            }
         }
         else {
            // Can't accumulate any more.  Push the current accumulation box
            // and start a new one
            coalesced_box_vector.push_back(bdryBox(accum_box,shear_dir,a_side,1));
            accum_box = box_vector[n];
         }
      }
   }

   return coalesced_box_vector;
}
#endif


void
SingleNullCoordSys::setBoundaryBoxes()
{
   m_boundary_boxes.resize(numBlocks() * 2 * SpaceDim);

#if CFG_DIM==3
   int shear_dir = TOROIDAL_DIR;
   int shear_interp_order = 2;
#endif

   for (int block_number=0; block_number<numBlocks(); ++block_number) {
   
      const MagBlockCoordSys& block_coord_sys = *(getCoordSys(block_number));
      const ProblemDomain& domain = block_coord_sys.domain();
      const Box& domain_box = domain.domainBox();

      for (int dir=0; dir<SpaceDim; ++dir) {
         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
         
            Vector<Box> boxes;

#if CFG_DIM==3
            if ( m_sheared_geometry && dir == shear_dir ) {
               Vector<Box> toroidal_boxes = getToroidalBoundaryBoxes(block_number, side, shear_interp_order);

               for (int n=0; n<toroidal_boxes.size(); ++n) {
                  boxes.push_back(toroidal_boxes[n]);
               }
            }
            else {
               if ( m_boundaries[block_number][dir + side*SpaceDim].isDomainBoundary() ) {
                  boxes.push_back(bdryBox(domain_box, dir, side, 1));
               }
            }
#else
            if ( m_boundaries[block_number][dir + side*SpaceDim].isDomainBoundary() ) {
               boxes.push_back(bdryBox(domain_box, dir, side, 1));
            }
#endif

            m_boundary_boxes[block_number*2*SpaceDim + dir + side*SpaceDim] = boxes;
         }
      }
   }
}

#if CFG_DIM==3
bool
SingleNullCoordSys::fixPoloidalBlock( int&      a_num_block,
                                      IntVect&  a_iv ) const
{
   const SingleNullBlockCoordSys* coord_sys = (SingleNullBlockCoordSys*)getCoordSys(a_num_block);
   const ProblemDomain& domain = coord_sys->domain();
   const Box& domain_box = domain.domainBox();

   // Make sure a_iv at least lies within the radial and toroidal block boundaries
   for (int dir=0; dir<SpaceDim; ++dir) {
      if ( dir != POLOIDAL_DIR ) {
         if ( a_iv[dir] < domain_box.smallEnd(dir) || a_iv[dir] > domain_box.bigEnd(dir) ) {
            MayDay::Error("SingleNullCoordSys::fixPoloidalBlock(): a_iv is outside the radial and/or toroidal boundary");

         }
      }
   }

   int poloidal_block = poloidalBlockNumber(a_num_block);
   int toroidal_sector = toroidalBlockNumber(a_num_block);

   int valid_poloidal_block = -1;
   bool block_fixed = false;

   if ( a_iv[POLOIDAL_DIR] < domain_box.smallEnd(POLOIDAL_DIR) ) {

      switch( poloidal_block )
         {
         case SingleNullBlockCoordSys::LCORE:
            valid_poloidal_block = SingleNullBlockCoordSys::RCORE;
            break;
         case SingleNullBlockCoordSys::RCORE:
            valid_poloidal_block = SingleNullBlockCoordSys::LCORE;
            break;
         case SingleNullBlockCoordSys::LPF:
            valid_poloidal_block = SingleNullBlockCoordSys::RPF;
            break;
         case SingleNullBlockCoordSys::RPF:
            valid_poloidal_block = -1;
            break;
         case SingleNullBlockCoordSys::LSOL:
            valid_poloidal_block = SingleNullBlockCoordSys::LCSOL;
            break;
         case SingleNullBlockCoordSys::LCSOL:
            valid_poloidal_block = SingleNullBlockCoordSys::RCSOL;
            break;
         case SingleNullBlockCoordSys::RCSOL:
            valid_poloidal_block = SingleNullBlockCoordSys::RSOL;
            break;
         case SingleNullBlockCoordSys::RSOL:
            valid_poloidal_block = -1;
            break;
         default:
            MayDay::Error("SingleNullCoordSys::checkNonPoloidalIndices(): Unrecogninzed block");
         }

      if ( valid_poloidal_block != -1 ) { 
         a_num_block = blockNumber( valid_poloidal_block, toroidal_sector );

         const SingleNullBlockCoordSys* valid_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(a_num_block);
         const ProblemDomain& valid_domain = valid_coord_sys->domain();
         const Box& valid_domain_box = valid_domain.domainBox();

         a_iv[POLOIDAL_DIR] = valid_domain_box.bigEnd(POLOIDAL_DIR) - (domain_box.smallEnd(POLOIDAL_DIR) - a_iv[POLOIDAL_DIR]) + 1;
         block_fixed = true;
      }
   }
   else if ( a_iv[POLOIDAL_DIR] > domain_box.bigEnd(POLOIDAL_DIR) ) {

      switch( poloidal_block )
         {
         case SingleNullBlockCoordSys::LCORE:
            valid_poloidal_block = SingleNullBlockCoordSys::RCORE;
            break;
         case SingleNullBlockCoordSys::RCORE:
            valid_poloidal_block = SingleNullBlockCoordSys::LCORE;
            break;
         case SingleNullBlockCoordSys::LPF:
            valid_poloidal_block = -1;
            break;
         case SingleNullBlockCoordSys::RPF:
            valid_poloidal_block = SingleNullBlockCoordSys::LPF;
            break;
         case SingleNullBlockCoordSys::LSOL:
            valid_poloidal_block = -1;
            break;
         case SingleNullBlockCoordSys::LCSOL:
            valid_poloidal_block = SingleNullBlockCoordSys::LSOL;
            break;
         case SingleNullBlockCoordSys::RCSOL:
            valid_poloidal_block = SingleNullBlockCoordSys::LCSOL;
            break;
         case SingleNullBlockCoordSys::RSOL:
            valid_poloidal_block = SingleNullBlockCoordSys::RCSOL;
            break;
         default:
            MayDay::Error("SingleNullCoordSys::checkNonPoloidalIndices(): Unrecogninzed block");
         }

      if ( valid_poloidal_block != -1 ) { 
         a_num_block = blockNumber( valid_poloidal_block, toroidal_sector );

         const SingleNullBlockCoordSys* valid_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(a_num_block);
         const ProblemDomain& valid_domain = valid_coord_sys->domain();
         const Box& valid_domain_box = valid_domain.domainBox();

         a_iv[POLOIDAL_DIR] = valid_domain_box.smallEnd(POLOIDAL_DIR) + (a_iv[POLOIDAL_DIR] - domain_box.bigEnd(POLOIDAL_DIR)) - 1;
         block_fixed = true;
      }
   }

   return block_fixed;
}

#endif


#include "NamespaceFooter.H"
