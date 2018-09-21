#include "SingleNullCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "SingleNullBlockCoordSysModel.H"
#include "BlockBoundary.H"
#include "DataArray.H"
#include "Directions.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

const std::string SingleNullCoordSys::pp_name = "singlenull";

SingleNullCoordSys::SingleNullCoordSys( ParmParse& a_pp_grid,
                                        ParmParse& a_pp_geom )
{
   readGridParams( a_pp_grid );
   
   m_model_geometry = false;
   a_pp_geom.query( "model_geometry", m_model_geometry);

   int num_blocks = m_original_eight_blocks? 8: NUM_SINGLE_NULL_BLOCKS;

   Vector<Box> domain_boxes(num_blocks);
   for ( int block_number=0; block_number<num_blocks; ++block_number ) {
      domain_boxes[block_number] = Box(lo_mapped_index(block_number), hi_mapped_index(block_number));
   }

   /*
     Compute the computational mesh size.  We arbitrarily set the radial and poloidal width
     of the core computational domain to 1.  This, together with the input number of
     cells across the core radial and poloidal widths sets the uniform cell size for
     all of the blocks.  The dimension of the blocks other than the core blocks are thereby
     determined by the specifed number of cells in those blocks.
   */
   int radial_width;
   int core_poloidal_width;
   if ( num_blocks == 8 ) {
      radial_width = domain_boxes[LCORE].size(RADIAL_DIR);
      core_poloidal_width =  domain_boxes[LCORE].size(POLOIDAL_DIR)
         + domain_boxes[RCORE].size(POLOIDAL_DIR);
   }
   else {
      radial_width = domain_boxes[MCORE].size(RADIAL_DIR);
      core_poloidal_width = domain_boxes[MCORE].size(POLOIDAL_DIR)
         + domain_boxes[LCORE].size(POLOIDAL_DIR)
         + domain_boxes[RCORE].size(POLOIDAL_DIR);
   }

   RealVect dx;
   dx[RADIAL_DIR]   = 1./(double)radial_width;
   dx[POLOIDAL_DIR] = 1./(double)core_poloidal_width;

   Vector<MagBlockCoordSys *> coord_vec;
   for ( int block_number = 0; block_number < num_blocks; ++block_number ) {
      if (!m_model_geometry) {
         MagBlockCoordSys* geom = new SingleNullBlockCoordSys( a_pp_geom, ProblemDomain(domain_boxes[block_number]), dx, block_number );
         coord_vec.push_back(geom);
      }
      else {
         MagBlockCoordSys* geom = new SingleNullBlockCoordSysModel( a_pp_geom, ProblemDomain(domain_boxes[block_number]), dx, block_number );
         coord_vec.push_back(geom);
      }

   }

   defineCoordSystemsAndBoundaries(coord_vec);

   //if (!m_model_geometry) setXPointNeighborhood();
   setXPointNeighborhood();

   // Define the boundary conditions for divergence cleaning (whether or not they're used)
   for (int block = 0; block < num_blocks; block++) {
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block];

      for (int dir=0; dir<SpaceDim; ++dir) {
         for (int side=0; side<2; ++side) {
            if ( blockBoundaries[dir + side*SpaceDim].isDomainBoundary() ) {
               double bc_value = 0.;
               int bc_type = EllipticOpBC::DIRICHLET;     // Homogeneous Dirichlet
               m_divergence_cleaning_bc.setBCType(block, dir, side, bc_type);
               m_divergence_cleaning_bc.setBCValue(block, dir, side, bc_value);
            }
         }
      }
   }

   //Checks if the grid cell sizes are continious across block boundaries (assumed for the model geoemtry)
   if (m_model_geometry) checkGridConsistency();
   
   
}



SingleNullCoordSys::~SingleNullCoordSys()
{
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

   for (int iblock = 0; iblock < numBlocks(); iblock++) {

      IndicesTransformation it;
      IntVect shift;
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[iblock];
      
      
      int block_type;
      bool poloidally_truncated;
      
      if (m_model_geometry) {
         const SingleNullBlockCoordSysModel* coord_sys = (const SingleNullBlockCoordSysModel*)m_coordSysVect[iblock];
         block_type = coord_sys->blockType();
         poloidally_truncated = coord_sys->truncated();
      }

      else {
         const SingleNullBlockCoordSys* coord_sys = (const SingleNullBlockCoordSys*)m_coordSysVect[iblock];
         block_type = coord_sys->blockType();
         poloidally_truncated = coord_sys->truncated();
      }
     
      
      if( block_type == LCORE ) {

         blockBoundaries[RADIAL_DIR].define(bc_tag);
         shift = (BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LCSOL);

         if ( poloidally_truncated ) {
           blockBoundaries[POLOIDAL_DIR].define(bc_tag);
         }
         else {
           shift = -BLOCK_SEPARATION*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                     m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
           it.defineFromTranslation(shift);
           blockBoundaries[POLOIDAL_DIR].define(it,RCORE);
         }
         shift = -(BLOCK_SEPARATION+1)*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                  + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCORE);
      }
      else if( block_type == RCORE ) {

         blockBoundaries[RADIAL_DIR].define(bc_tag);
         shift = (BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RCSOL);

         shift = (BLOCK_SEPARATION+1)*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,LCORE);

         if ( poloidally_truncated ) {
           blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
         }
         else {
           shift = BLOCK_SEPARATION*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                     m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
           it.defineFromTranslation(shift);
           blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCORE);
         }
      }
      else if ( block_type == LCSOL ) {

         shift = - (BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR].define(it,LCORE);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         if ( poloidally_truncated ) {
           blockBoundaries[POLOIDAL_DIR].define(bc_tag);
         }
         else {
           shift = -BLOCK_SEPARATION*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                     m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
           it.defineFromTranslation(shift);
           blockBoundaries[POLOIDAL_DIR].define(it,RCSOL);
         }
         shift = BLOCK_SEPARATION*m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LSOL);
      }
      else if ( block_type == RCSOL ) {

         shift = - (BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR].define(it,RCORE);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         shift = -BLOCK_SEPARATION*m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,RSOL);

         if ( poloidally_truncated ) {
           blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
         }
         else {
           shift = BLOCK_SEPARATION*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                     m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
           it.defineFromTranslation(shift);
           blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCSOL);
         }
      }
      else if ( block_type == LSOL ) {

         shift = - (BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR].define(it,LPF);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         shift = -BLOCK_SEPARATION*m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,LCSOL);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
      }
      else if ( block_type == RSOL ) {

         shift = - (BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR].define(it,RPF);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         blockBoundaries[POLOIDAL_DIR].define(bc_tag);
         shift = BLOCK_SEPARATION*m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCSOL);
      }
      else if ( block_type == LPF ) {

         blockBoundaries[RADIAL_DIR].define(bc_tag);
         shift = (BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LSOL);

         shift = -(2*BLOCK_SEPARATION+1)*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,RPF);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
      }
      else if ( block_type == RPF ) {

         blockBoundaries[RADIAL_DIR].define(bc_tag);
         shift = (BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(RADIAL_DIR)) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RSOL);

         blockBoundaries[POLOIDAL_DIR].define(bc_tag);
         shift = (2*BLOCK_SEPARATION+1)*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LPF);
      }
      else {
            MayDay::Error("SingleNullCoordSys::defineBoundaries(): case not implemented");
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

   for (int iblock = 0; iblock < numBlocks(); iblock++) {

      IndicesTransformation it;
      IntVect shift;
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[iblock];

      int block_type;
      
      if (m_model_geometry) {
         const SingleNullBlockCoordSysModel* coord_sys = (const SingleNullBlockCoordSysModel*)m_coordSysVect[iblock];
         block_type = coord_sys->blockType();
      }
      
      else {
         const SingleNullBlockCoordSys* coord_sys = (const SingleNullBlockCoordSys*)m_coordSysVect[iblock];
         block_type = coord_sys->blockType();
       }

      if( block_type == MCORE ) {
         
         blockBoundaries[RADIAL_DIR].define(bc_tag);

         shift =  (BLOCK_SEPARATION * L0) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(it,MCSOL);

         shift = -(BLOCK_SEPARATION * rL1 + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,RCORE);

         shift =  (BLOCK_SEPARATION * lL1 + m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCORE);
      }
      else if( block_type == LCORE ) {
            
         blockBoundaries[RADIAL_DIR].define(bc_tag);

         shift = (BLOCK_SEPARATION * L0) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LCSOL);

         shift = -(BLOCK_SEPARATION * lL1 + m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,MCORE);

         shift = -((BLOCK_SEPARATION+1) * (lL1 + rL1)
                   + m_numcells_lcore_poloidal + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCORE);
      }
      else if( block_type == RCORE ) {

         blockBoundaries[RADIAL_DIR].define(bc_tag);

         shift = (BLOCK_SEPARATION * L0) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RCSOL);

         shift = ((BLOCK_SEPARATION+1) * (lL1 + rL1)
                  + m_numcells_lcore_poloidal + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,LCORE);

         shift = (BLOCK_SEPARATION * rL1 + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,MCORE);
      }
      else if ( block_type == MCSOL ) {

         shift = -(BLOCK_SEPARATION * L0) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR].define(it,MCORE);

         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         shift = -(BLOCK_SEPARATION * rL1 + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,RCSOL);

         shift =  (BLOCK_SEPARATION * lL1 + m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LCSOL);
      }
      else if ( block_type == LCSOL ) {

         shift = -(BLOCK_SEPARATION * L0) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR].define(it,LCORE);

         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         shift = -(BLOCK_SEPARATION * lL1 + m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,MCSOL);

         shift =  ((BLOCK_SEPARATION+1) * lL1 - m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LSOL);
      }
      else if ( block_type == RCSOL ) {

         shift = -(BLOCK_SEPARATION * L0) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR].define(it,RCORE);

         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         shift = -((BLOCK_SEPARATION+1) * rL1 - m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,RSOL);

         shift = (BLOCK_SEPARATION * rL1 + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,MCSOL);
      }
      else if ( block_type == LPF ) {

         blockBoundaries[RADIAL_DIR].define(bc_tag);

         shift = (BLOCK_SEPARATION * L0) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(it,LSOL);

         shift = -(2*(BLOCK_SEPARATION+1) * (lL1 + rL1))*BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,RPF);

         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
      }
      else if ( block_type == RPF ) {

         blockBoundaries[RADIAL_DIR].define(bc_tag);

         shift = (BLOCK_SEPARATION * L0) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR + SpaceDim].define(it,RSOL);

         blockBoundaries[POLOIDAL_DIR].define(bc_tag);

         shift =  (2*(BLOCK_SEPARATION+1) * (lL1 + rL1))*BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,LPF);
      }
      else if ( block_type == LSOL ) {
            
         shift = -(BLOCK_SEPARATION * L0) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR].define(it,LPF);

         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         shift = -((BLOCK_SEPARATION+1) * lL1 - m_numcells_lcore_poloidal ) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR].define(it,LCSOL);

         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(bc_tag);
      }
      else if ( block_type == RSOL ) {

         shift = -(BLOCK_SEPARATION * L0) * BASISV(RADIAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[RADIAL_DIR].define(it,RPF);
         
         blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

         blockBoundaries[POLOIDAL_DIR].define(bc_tag);

         shift = ((BLOCK_SEPARATION+1) * rL1 - m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
         it.defineFromTranslation(shift);
         blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,RCSOL);
      }
      else {
         MayDay::Error("SingleNullCoordSys::defineBoundaries(): case not implemented");
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

   if ( src_coord_sys->isValid(a_xiSrc) ) {
      a_n_valid = a_nSrc;
      a_xi_valid = a_xiSrc;
   }
   else {

      RealVect X = src_coord_sys->realCoord(a_xiSrc);

      RealVect XminusXpt = X - m_Xpoint;
      double distance_to_Xpt = XminusXpt.vectorLength();

      if ( distance_to_Xpt < m_xpoint_radius ) {
         a_n_valid = findBlockNearXpt(X);
      }
      else {

         double lo_mapped_radial   = src_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
         double hi_mapped_radial   = src_coord_sys->upperMappedCoordinate(RADIAL_DIR);
         double lo_mapped_poloidal = src_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);
         double hi_mapped_poloidal = src_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);
         bool poloidally_truncated = src_coord_sys->truncated();

         switch (a_nSrc)
            {
            case LCORE:
               if ( a_xiSrc[0] < lo_mapped_radial ||
                    (a_xiSrc[1] < lo_mapped_poloidal && poloidally_truncated) ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] > hi_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = LSOL;
                  }
                  else if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = RCSOL;
                  }
                  else {
                     a_n_valid = LCSOL;
                  }
               }
               else {
                  a_n_valid = RCORE;
               }
               break;
            case RCORE:
               if ( a_xiSrc[0] < lo_mapped_radial ||
                    (a_xiSrc[1] > hi_mapped_poloidal && poloidally_truncated) ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] > hi_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = LCSOL;
                  }
                  else if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = RSOL;
                  }
                  else {
                     a_n_valid = RCSOL;
                  }
               }
               else {
                  a_n_valid = LCORE;
               }
               break;
            case LCSOL:
               if ( a_xiSrc[0] > hi_mapped_radial ||
                    (a_xiSrc[1] < lo_mapped_poloidal && poloidally_truncated) ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] < lo_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = LPF;
                  }
                  else if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = RCORE;
                  }
                  else {
                     a_n_valid = LCORE;
                  }
               }
               else if (a_xiSrc[1] > hi_mapped_poloidal) {
                  a_n_valid = LSOL;
               }
               else {
                  a_n_valid = RCSOL;
               }
               break;
            case RCSOL:
               if ( a_xiSrc[0] > hi_mapped_radial ||
                    (a_xiSrc[1] > hi_mapped_poloidal && poloidally_truncated) ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] < lo_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = LCORE;
                  }
                  else if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = RPF;
                  }
                  else {
                     a_n_valid = RCORE;
                  }
               }
               else if (a_xiSrc[1] > hi_mapped_poloidal) {
                  a_n_valid = LCSOL;
               }
               else if (a_xiSrc[1] < lo_mapped_poloidal) {
                  a_n_valid = RSOL;
               }
               break;
            case LSOL:
               if ( a_xiSrc[0] > hi_mapped_radial ||
                    a_xiSrc[1] > hi_mapped_poloidal ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] < lo_mapped_radial) {
                  if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = LCORE;
                  }
                  else {
                     a_n_valid = LPF;
                  }
               }
               else {
                  a_n_valid = LCSOL;
               }
               break;
            case RSOL:
               if ( a_xiSrc[0] > hi_mapped_radial ||
                    a_xiSrc[1] < lo_mapped_poloidal ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] < lo_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = RCORE;
                  }
                  else {
                     a_n_valid = RPF;
                  }
               }
               else {
                  a_n_valid = RCSOL;
               }
               break;
            case LPF:
               if ( a_xiSrc[0] < lo_mapped_radial ||
                    a_xiSrc[1] > hi_mapped_poloidal ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] > hi_mapped_radial) {
                  if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = LCSOL;
                  }
                  else {
                     a_n_valid = LSOL;
                  }
               }
               else {
                  a_n_valid = RPF;
               }
               break;
            case RPF:
               if ( a_xiSrc[0] < lo_mapped_radial ||
                    a_xiSrc[1] < lo_mapped_poloidal ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] > hi_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = RCSOL;
                  }
                  else {
                     a_n_valid = RSOL;
                  }
               }
               else {
                  a_n_valid = LPF;
               }
               break;
            default:
               MayDay::Error("SingleNullCoordSys::blockRemapping(): Unknown block number");
            }
      }

      if (a_n_valid != -1) {

         const SingleNullBlockCoordSys* valid_cs = (SingleNullBlockCoordSys*)getCoordSys(a_n_valid);

         a_xi_valid = valid_cs->mappedCoord(X);

         // Check that the point is actually valid
         if ( !valid_cs->isValid(a_xi_valid) ) {

            if (distance_to_Xpt < m_xpoint_radius) {

               // If we've arrived here, it's because the separatrix isn't perfectly straight
               // in the X point neighborhood, as was assumed in findBlock().  We therefore
               // have to deal with the special cases.

               int n_valid_new = a_n_valid;

               switch ( a_n_valid )
                  {
                  case LCORE:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = LCSOL;
                     }
                     break;
                  case LCSOL:
                     if (a_xi_valid[RADIAL_DIR] < valid_cs->lowerMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = LCORE;
                     }
                     break;
                  case RCORE:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = RCSOL;
                     }
                     break;
                  case RCSOL:
                     if (a_xi_valid[RADIAL_DIR] < valid_cs->lowerMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = RCORE;
                     }
                     break;
                  case LPF:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = LSOL;
                     }
                     break;
                  case LSOL:
                     if (a_xi_valid[RADIAL_DIR] < valid_cs->lowerMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = LPF;
                     }
                     break;
                  case RPF:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = RSOL;
                     }
                     break;
                  case RSOL:
                     if (a_xi_valid[RADIAL_DIR] < valid_cs->lowerMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = RPF;
                     }
                     break;
                  default:
                     MayDay::Error("SingleNullCoordSys::blockRemapping(): bad block number");
                  }

#if 1
               if (n_valid_new != a_n_valid) {
                  
                  const SingleNullBlockCoordSys* new_valid_cs = (SingleNullBlockCoordSys*)getCoordSys(n_valid_new);

                  a_xi_valid = new_valid_cs->mappedCoord(X);
                  a_n_valid = n_valid_new;

                  if ( !new_valid_cs->isValid(a_xi_valid) ) {
#ifndef MODEL_GEOM
		    printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, new_valid_cs);
                    //                    exit(1);
#endif
                  }
               }
#endif
            }
            else {
#if 0
	            printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, valid_cs);
                    //                    exit(1);
#endif
        }
      }
    }
  }
}



void
SingleNullCoordSys::blockRemapping10( RealVect&       a_xi_valid,
                                      int&            a_n_valid,
                                      const RealVect& a_xiSrc,
                                      int             a_nSrc ) const
{
   const SingleNullBlockCoordSys* src_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(a_nSrc);

   if ( src_coord_sys->isValid(a_xiSrc) ) {
      a_n_valid = a_nSrc;
      a_xi_valid = a_xiSrc;
   }
   else {

      RealVect X = src_coord_sys->realCoord(a_xiSrc);

      RealVect XminusXpt = X - m_Xpoint;
      double distance_to_Xpt = XminusXpt.vectorLength();

      if ( distance_to_Xpt < m_xpoint_radius ) {
         a_n_valid = findBlockNearXpt(X);
      }
      else {

         double lo_mapped_radial   = src_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
         double hi_mapped_radial   = src_coord_sys->upperMappedCoordinate(RADIAL_DIR);
         double lo_mapped_poloidal = src_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);
         double hi_mapped_poloidal = src_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);
         bool poloidally_truncated = src_coord_sys->truncated();

         switch (a_nSrc)
            {
            case MCORE:
               if ( a_xiSrc[0] < lo_mapped_radial ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] > hi_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = LCSOL;
                  }
                  else if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = RCSOL;
                  }
                  else {
                     a_n_valid = MCSOL;
                  }
               }
               else {
                  if ( a_xiSrc[1] > hi_mapped_poloidal ) {
                     a_n_valid = LCORE;
                  }
                  else {
                     a_n_valid = RCORE;
                  }
               }
               break;
            case LCORE:
               if ( a_xiSrc[0] < lo_mapped_radial ||
                    (a_xiSrc[1] < lo_mapped_poloidal && poloidally_truncated) ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] > hi_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = LSOL;
                  }
                  else if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = MCSOL;
                  }
                  else {
                     a_n_valid = LCSOL;
                  }
               }
               else {
                  if ( a_xiSrc[1] > hi_mapped_poloidal ) {
                     a_n_valid = RCORE;
                  }
                  else {
                     a_n_valid = MCORE;
                  }
               }
               break;
            case RCORE:
               if ( a_xiSrc[0] < lo_mapped_radial ||
                    (a_xiSrc[1] > hi_mapped_poloidal && poloidally_truncated) ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] > hi_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = MCSOL;
                  }
                  else if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = RSOL;
                  }
                  else {
                     a_n_valid = RCSOL;
                  }
               }
               else {
                  if ( a_xiSrc[1] > hi_mapped_poloidal ) {
                     a_n_valid = MCORE;
                  }
                  else {
                     a_n_valid = LCORE;
                  }
               }
               break;
            case MCSOL:
               if ( a_xiSrc[0] > hi_mapped_radial ||
                    (a_xiSrc[1] < lo_mapped_poloidal && poloidally_truncated) ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] < lo_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = LCORE;
                  }
                  else if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = RCORE;
                  }
                  else {
                     a_n_valid = MCORE;
                  }
               }
               else if (a_xiSrc[1] > hi_mapped_poloidal) {
                  a_n_valid = LCSOL;
               }
               else {
                  a_n_valid = RCSOL;
               }
               break;
            case LCSOL:
               if ( a_xiSrc[0] > hi_mapped_radial ||
                    (a_xiSrc[1] < lo_mapped_poloidal && poloidally_truncated) ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] < lo_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = LPF;
                  }
                  else if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = MCORE;
                  }
                  else {
                     a_n_valid = LCORE;
                  }
               }
               else if (a_xiSrc[1] > hi_mapped_poloidal) {
                  a_n_valid = LSOL;
               }
               else {
                  a_n_valid = MCSOL;
               }
               break;
            case RCSOL:
               if ( a_xiSrc[0] > hi_mapped_radial ||
                    (a_xiSrc[1] > hi_mapped_poloidal && poloidally_truncated) ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] < lo_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = MCORE;
                  }
                  else if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = RPF;
                  }
                  else {
                     a_n_valid = RCORE;
                  }
               }
               else if (a_xiSrc[1] > hi_mapped_poloidal) {
                  a_n_valid = MCSOL;
               }
               else if (a_xiSrc[1] < lo_mapped_poloidal) {
                  a_n_valid = RSOL;
               }
               break;
            case LSOL:
               if ( a_xiSrc[0] > hi_mapped_radial ||
                    a_xiSrc[1] > hi_mapped_poloidal ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] < lo_mapped_radial) {
                  if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = LCORE;
                  }
                  else {
                     a_n_valid = LPF;
                  }
               }
               else {
                  a_n_valid = LCSOL;
               }
               break;
            case RSOL:
               if ( a_xiSrc[0] > hi_mapped_radial ||
                    a_xiSrc[1] < lo_mapped_poloidal ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] < lo_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = RCORE;
                  }
                  else {
                     a_n_valid = RPF;
                  }
               }
               else {
                  a_n_valid = RCSOL;
               }
               break;
            case LPF:
               if ( a_xiSrc[0] < lo_mapped_radial ||
                    a_xiSrc[1] > hi_mapped_poloidal ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] > hi_mapped_radial) {
                  if (a_xiSrc[1] < lo_mapped_poloidal) {
                     a_n_valid = LCSOL;
                  }
                  else {
                     a_n_valid = LSOL;
                  }
               }
               else {
                  a_n_valid = RPF;
               }
               break;
            case RPF:
               if ( a_xiSrc[0] < lo_mapped_radial ||
                    a_xiSrc[1] < lo_mapped_poloidal ) {
                  a_n_valid = -1;
               }
               else if (a_xiSrc[0] > hi_mapped_radial) {
                  if (a_xiSrc[1] > hi_mapped_poloidal) {
                     a_n_valid = RCSOL;
                  }
                  else {
                     a_n_valid = RSOL;
                  }
               }
               else {
                  a_n_valid = LPF;
               }
               break;
            default:
               MayDay::Error("SingleNullCoordSys::blockRemapping(): Unknown block number");
            }
      }

      if (a_n_valid != -1) {

         const SingleNullBlockCoordSys* valid_cs = (SingleNullBlockCoordSys*)getCoordSys(a_n_valid);

         a_xi_valid = valid_cs->mappedCoord(X);

         // Check that the point is actually valid
         if ( !valid_cs->isValid(a_xi_valid) ) {

            if (a_n_valid == MCORE) {
               cout << "Invalid MCORE point = " << a_xi_valid << endl;
               exit(1);
            }

            if (a_n_valid == MCSOL) {
               cout << "Invalid MCSOL point = " << a_xi_valid << endl;
               exit(1);
            }

            if (distance_to_Xpt < m_xpoint_radius) {

               // If we've arrived here, it's because the separatrix isn't perfectly straight
               // in the X point neighborhood, as was assumed in findBlockNearXpt().  We therefore
               // have to deal with the special cases.

               int n_valid_new = a_n_valid;

               switch ( a_n_valid )
                  {
#if 0
                  case MCORE:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = MCSOL;
                     }
                     break;
                  case MCSOL:
                     if (a_xi_valid[RADIAL_DIR] < valid_cs->lowerMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = MCORE;
                     }
                     break;
#endif
                  case LCORE:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = LCSOL;
                     }
                     break;
                  case LCSOL:
                     if (a_xi_valid[RADIAL_DIR] < valid_cs->lowerMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = LCORE;
                     }
                     break;
                  case RCORE:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = RCSOL;
                     }
                     break;
                  case RCSOL:
                     if (a_xi_valid[RADIAL_DIR] < valid_cs->lowerMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = RCORE;
                     }
                     break;
                  case LPF:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = LSOL;
                     }
                     break;
                  case LSOL:
                     if (a_xi_valid[RADIAL_DIR] < valid_cs->lowerMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = LPF;
                     }
                     break;
                  case RPF:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = RSOL;
                     }
                     break;
                  case RSOL:
                     if (a_xi_valid[RADIAL_DIR] < valid_cs->lowerMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = RPF;
                     }
                     break;
                  default:
                     MayDay::Error("SingleNullCoordSys::blockRemapping(): bad block number");
                  }

#if 1
               if (n_valid_new != a_n_valid) {

                  const SingleNullBlockCoordSys* new_valid_cs = (SingleNullBlockCoordSys*)getCoordSys(n_valid_new);

                  a_xi_valid = new_valid_cs->mappedCoord(X);
                  a_n_valid = n_valid_new;

                  if ( !new_valid_cs->isValid(a_xi_valid) ) {
#ifndef MODEL_GEOM
		      printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, new_valid_cs);
                     //              exit(1);
#endif
                  }
               }
#endif
            }
            else {
#if 1
#ifndef MODEL_GEOM
                     printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, valid_cs);
               //      exit(1);
#endif
#endif
            }
         }
      }
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
  cout << "xi not valid, nSrc = " << a_nSrc << " a_xiSrc = " << a_xiSrc << " X = " << a_X
       << " xi_valid = " << a_xi_valid << " nvalid = " << a_n_valid << endl;
  cout << "src block: "
       << a_src_cs->lowerMappedCoordinate(RADIAL_DIR) << " "
       << a_src_cs->upperMappedCoordinate(RADIAL_DIR) << " "
       << a_src_cs->lowerMappedCoordinate(POLOIDAL_DIR) << " "
       << a_src_cs->upperMappedCoordinate(POLOIDAL_DIR) << endl;
  cout << "valid block: "
       << a_valid_cs->lowerMappedCoordinate(RADIAL_DIR) << " "
       << a_valid_cs->upperMappedCoordinate(RADIAL_DIR) << " "
       << a_valid_cs->lowerMappedCoordinate(POLOIDAL_DIR) << " "
       << a_valid_cs->upperMappedCoordinate(POLOIDAL_DIR) << endl;
}



void
SingleNullCoordSys::setXPointNeighborhood()
#if CFG_DIM==3
{
  double cut_frac = 1.;
  double cut_frac2 = 0.15;

  const SingleNullBlockCoordSys* lcore_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(LCORE);

#if CFG_DIM==3
  RealVect Xpoint_mapped(lcore_coord_sys->upperMappedCoordinate(RADIAL_DIR),lcore_coord_sys->upperMappedCoordinate(POLOIDAL_DIR),0.);
#else
  RealVect Xpoint_mapped(lcore_coord_sys->upperMappedCoordinate(RADIAL_DIR),lcore_coord_sys->upperMappedCoordinate(POLOIDAL_DIR));
#endif
  m_Xpoint = lcore_coord_sys->realCoord(Xpoint_mapped);

  double rad_mapped;
  double pol_mapped;
  double cut_frac_tmp;

  cut_frac_tmp = ( lcore_coord_sys->truncated() || numBlocks() == 10 )? cut_frac: cut_frac2;

  rad_mapped = (1.-cut_frac_tmp)*lcore_coord_sys->upperMappedCoordinate(RADIAL_DIR)
                  + cut_frac_tmp*lcore_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
  pol_mapped = lcore_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);

#if CFG_DIM==3
  m_ray[0] = lcore_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped,0.)) - m_Xpoint;
#else
  m_ray[0] = lcore_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped)) - m_Xpoint;
#endif

  const SingleNullBlockCoordSys* rcore_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(RCORE);

  cut_frac_tmp = ( rcore_coord_sys->truncated() || numBlocks() == 10 )? cut_frac: cut_frac2;

  rad_mapped = rcore_coord_sys->upperMappedCoordinate(RADIAL_DIR);
  pol_mapped = (1.-cut_frac_tmp)*rcore_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR)
    + cut_frac_tmp*rcore_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);

#if CFG_DIM==3
  m_ray[1] = rcore_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped,0.)) - m_Xpoint;
#else
  m_ray[1] = rcore_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped)) - m_Xpoint;
#endif

  const SingleNullBlockCoordSys* rcsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(RCSOL);

  cut_frac_tmp = ( rcsol_coord_sys->truncated() || numBlocks() == 10 )? cut_frac: cut_frac2;

  rad_mapped = (1.-cut_frac_tmp)*rcsol_coord_sys->lowerMappedCoordinate(RADIAL_DIR)
    + cut_frac_tmp*rcsol_coord_sys->upperMappedCoordinate(RADIAL_DIR);
  pol_mapped = rcsol_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);

#if CFG_DIM==3
  m_ray[2] = rcsol_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped,0.)) - m_Xpoint;
#else
  m_ray[2] = rcsol_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped)) - m_Xpoint;
#endif

  const SingleNullBlockCoordSys* rsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(RSOL);

  rad_mapped = rsol_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
  pol_mapped = (1.-cut_frac)*rsol_coord_sys->upperMappedCoordinate(POLOIDAL_DIR)
    + cut_frac*rsol_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);

#if CFG_DIM==3
  m_ray[3] = rsol_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped,0.)) - m_Xpoint;
#else
  m_ray[3] = rsol_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped)) - m_Xpoint;
#endif

  const SingleNullBlockCoordSys* rpf_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(RPF);

  rad_mapped = (1.-cut_frac)*rpf_coord_sys->upperMappedCoordinate(RADIAL_DIR)
    + cut_frac*rpf_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
  pol_mapped = rpf_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);

#if CFG_DIM==3
  m_ray[4] = rpf_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped,0.)) - m_Xpoint;
#else
  m_ray[4] = rpf_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped)) - m_Xpoint;
#endif

  const SingleNullBlockCoordSys* lpf_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(LPF);

  rad_mapped = lpf_coord_sys->upperMappedCoordinate(RADIAL_DIR);
  pol_mapped = (1.-cut_frac)*lpf_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR)
    + cut_frac*lpf_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);

#if CFG_DIM==3
  m_ray[5] = lpf_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped,0.)) - m_Xpoint;
#else
  m_ray[5] = lpf_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped)) - m_Xpoint;
#endif

  const SingleNullBlockCoordSys* lsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(LSOL);

  rad_mapped = (1.-cut_frac)*lsol_coord_sys->lowerMappedCoordinate(RADIAL_DIR)
    + cut_frac*lsol_coord_sys->upperMappedCoordinate(RADIAL_DIR);
  pol_mapped = lsol_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);

#if CFG_DIM==3
  m_ray[6] = lsol_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped,0.)) - m_Xpoint;
#else
  m_ray[6] = lsol_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped)) - m_Xpoint;
#endif

  const SingleNullBlockCoordSys* lcsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(LCSOL);

  cut_frac_tmp = ( lcsol_coord_sys->truncated() || numBlocks() == 10 )? cut_frac: cut_frac2;

  rad_mapped = lcsol_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
  pol_mapped = (1.-cut_frac_tmp)*lcsol_coord_sys->upperMappedCoordinate(POLOIDAL_DIR)
    + cut_frac_tmp*lcsol_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);

#if CFG_DIM==3
  m_ray[7] = lcsol_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped,0.)) - m_Xpoint;
#else
  m_ray[7] = lcsol_coord_sys->realCoord(RealVect(rad_mapped,pol_mapped)) - m_Xpoint;
#endif

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

  m_sector_block[0] = RCORE;
  m_sector_block[1] = RCSOL;
  m_sector_block[2] = RSOL;
  m_sector_block[3] = RPF;
  m_sector_block[4] = LPF;
  m_sector_block[5] = LSOL;
  m_sector_block[6] = LCSOL;
  m_sector_block[7] = LCORE;
}
#else
{
   {
      // Block boundary between LCORE and RCORE

      const SingleNullBlockCoordSys* lcore_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(LCORE);

      Box boundary_box = surroundingNodes(bdryHi(lcore_coord_sys->domain().domainBox(), POLOIDAL_DIR, 1), RADIAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      lcore_coord_sys->getNodeCenteredRealCoords(X);
     
      m_Xpoint = RealVect(X(hi,0),X(hi,1));

      // Initialize the ray to lie on the block boundary, pointing away from X point
      m_ray[0] = RealVect(X(lo,0),X(lo,1)) - m_Xpoint;

      // Set the radius of the X point neighborhood, which we currently take to as the
      // length of the block interface between the LCORE and RCORE.  If there is any
      // need to reduce the size of the neighborhood, here is where that should be done.  

      m_xpoint_radius = m_ray[0].vectorLength();

      // This loop isn't actually necessary unless the preceding m_xpoint_radius
      // definition has somehow been manually overridden to a smaller value
      for( IntVect iv(hi); iv[RADIAL_DIR] >= lo[RADIAL_DIR]; iv[RADIAL_DIR]--) {
         RealVect diff = RealVect(X(iv,0),X(iv,1)) - m_Xpoint;

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

      const SingleNullBlockCoordSys* rcore_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(RCORE);

      Box boundary_box = surroundingNodes(bdryHi(rcore_coord_sys->domain().domainBox(), RADIAL_DIR, 1), POLOIDAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      rcore_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point
      m_ray[1] = RealVect(X(hi,0),X(hi,1)) - m_Xpoint;

      for( IntVect iv(lo); iv[POLOIDAL_DIR] <= hi[POLOIDAL_DIR]; iv[POLOIDAL_DIR]++) {
         RealVect diff = RealVect(X(iv,0),X(iv,1)) - m_Xpoint;

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

      const SingleNullBlockCoordSys* rcsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(RCSOL);

      Box boundary_box = surroundingNodes(bdryLo(rcsol_coord_sys->domain().domainBox(), POLOIDAL_DIR, 1), RADIAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      rcsol_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point
      m_ray[2] = RealVect(X(hi,0),X(hi,1)) - m_Xpoint;

      for( IntVect iv(lo); iv[RADIAL_DIR] <= hi[RADIAL_DIR]; iv[RADIAL_DIR]++) {
         RealVect diff = RealVect(X(iv,0),X(iv,1)) - m_Xpoint;

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

      const SingleNullBlockCoordSys* rsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(RSOL);

      Box boundary_box = surroundingNodes(bdryLo(rsol_coord_sys->domain().domainBox(), RADIAL_DIR, 1), POLOIDAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      rsol_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point
      m_ray[3] = RealVect(X(lo,0),X(lo,1)) - m_Xpoint;

      for( IntVect iv(hi); iv[POLOIDAL_DIR] >= lo[POLOIDAL_DIR]; iv[POLOIDAL_DIR]--) {
         RealVect diff = RealVect(X(iv,0),X(iv,1)) - m_Xpoint;

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

      const SingleNullBlockCoordSys* rpf_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(RPF);

      Box boundary_box = surroundingNodes(bdryHi(rpf_coord_sys->domain().domainBox(), POLOIDAL_DIR, 1), RADIAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      rpf_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point
      m_ray[4] = RealVect(X(lo,0),X(lo,1)) - m_Xpoint;

      for( IntVect iv(hi); iv[RADIAL_DIR] >= lo[RADIAL_DIR]; iv[RADIAL_DIR]--) {
         RealVect diff = RealVect(X(iv,0),X(iv,1)) - m_Xpoint;

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

      const SingleNullBlockCoordSys* lpf_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(LPF);

      Box boundary_box = surroundingNodes(bdryHi(lpf_coord_sys->domain().domainBox(), RADIAL_DIR, 1), POLOIDAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      lpf_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point
      m_ray[5] = RealVect(X(hi,0),X(hi,1)) - m_Xpoint;

      for( IntVect iv(lo); iv[POLOIDAL_DIR] <= hi[POLOIDAL_DIR]; iv[POLOIDAL_DIR]++) {
         RealVect diff = RealVect(X(iv,0),X(iv,1)) - m_Xpoint;

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

      const SingleNullBlockCoordSys* lsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(LSOL);

      Box boundary_box = surroundingNodes(bdryLo(lsol_coord_sys->domain().domainBox(), POLOIDAL_DIR, 1), RADIAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      lsol_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point
      m_ray[6] = RealVect(X(hi,0),X(hi,1)) - m_Xpoint;

      for( IntVect iv(lo); iv[RADIAL_DIR] <= hi[RADIAL_DIR]; iv[RADIAL_DIR]++) {
         RealVect diff = RealVect(X(iv,0),X(iv,1)) - m_Xpoint;

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

      const SingleNullBlockCoordSys* lcsol_coord_sys = (SingleNullBlockCoordSys*)getCoordSys(LCSOL);

      Box boundary_box = surroundingNodes(bdryLo(lcsol_coord_sys->domain().domainBox(), RADIAL_DIR, 1), POLOIDAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      lcsol_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point
      m_ray[7] = RealVect(X(lo,0),X(lo,1)) - m_Xpoint;

      for( IntVect iv(hi); iv[POLOIDAL_DIR] >= lo[POLOIDAL_DIR]; iv[POLOIDAL_DIR]--) {
         RealVect diff = RealVect(X(iv,0),X(iv,1)) - m_Xpoint;

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

   m_sector_block[0] = RCORE;
   m_sector_block[1] = RCSOL;
   m_sector_block[2] = RSOL;
   m_sector_block[3] = RPF;
   m_sector_block[4] = LPF;
   m_sector_block[5] = LSOL;
   m_sector_block[6] = LCSOL;
   m_sector_block[7] = LCORE;
}
#endif



int
SingleNullCoordSys::findBlockNearXpt(const RealVect& X) const
{
   RealVect local_X = X - m_Xpoint;

#if CFG_DIM==3
   RealVect ray0_normal(m_ray[0][1],-m_ray[0][0],0.);  // m_ray[0] rotated 90 degrees clockwise
#else
   RealVect ray0_normal(m_ray[0][1],-m_ray[0][0]);  // m_ray[0] rotated 90 degrees clockwise
#endif

   double angle = acos( local_X.dotProduct(m_ray[0]) / local_X.vectorLength() );

   if (local_X.dotProduct(ray0_normal) < 0.) {
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
      cout << "Looking for block of " << X << endl;
      MayDay::Error("SingleNullCoordSys::findBlock(): Couldn't find block");
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

         const Box& lcore_domain_box = getCoordSys(LCORE)->domain().domainBox();
         const Box& rcore_domain_box = getCoordSys(RCORE)->domain().domainBox();

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

         const Box& lcsol_domain_box = getCoordSys(LCSOL)->domain().domainBox();
         const Box& rcsol_domain_box = getCoordSys(RCSOL)->domain().domainBox();
         const Box& lsol_domain_box  = getCoordSys(LSOL)->domain().domainBox();
         const Box& rsol_domain_box  = getCoordSys(RSOL)->domain().domainBox();

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

         const Box& lpf_domain_box = getCoordSys(LPF)->domain().domainBox();
         const Box& rpf_domain_box = getCoordSys(RPF)->domain().domainBox();

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

         const Box& lsol_domain_box = getCoordSys(LSOL)->domain().domainBox();
         const Box& lpf_domain_box = getCoordSys(LPF)->domain().domainBox();

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

         const Box& rsol_domain_box = getCoordSys(RSOL)->domain().domainBox();
         const Box& rpf_domain_box = getCoordSys(RPF)->domain().domainBox();

         Box fab_box = minBox( bdryLo(rsol_domain_box, POLOIDAL_DIR, 1),
                               bdryLo(rpf_domain_box, POLOIDAL_DIR, 1) );

         FArrayBox fab(fab_box, 1);
         data_array->setData(fab, false);
      }
   }
}


IntVect
SingleNullCoordSys::lo_mapped_index(int block) const
{
  IntVect index;

  if ( m_original_eight_blocks ) {

     int numcells_core_poloidal = m_numcells_lcore_poloidal + m_numcells_rcore_poloidal;

     switch ( block )
        {
        case LCORE:
           index[RADIAL_DIR] = 0;
           index[POLOIDAL_DIR] = BLOCK_SEPARATION*numcells_core_poloidal/2;
           break;
        case RCORE:
           index[RADIAL_DIR] = 0;
           index[POLOIDAL_DIR] = -(BLOCK_SEPARATION+1)*numcells_core_poloidal/2;
           break;
        case LCSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*m_numcells_core_radial;
           index[POLOIDAL_DIR] = BLOCK_SEPARATION*numcells_core_poloidal/2;
           break;
        case RCSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*m_numcells_core_radial;
           index[POLOIDAL_DIR] = -(BLOCK_SEPARATION+1)*numcells_core_poloidal/2;
           break;
        case LSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*m_numcells_core_radial;      // same radial extent as CSOL
           index[POLOIDAL_DIR] = (2*BLOCK_SEPARATION+1)*numcells_core_poloidal/2;
           break;
        case RSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*m_numcells_core_radial;      // same radial extent as CSOL
           index[POLOIDAL_DIR] = -(2*BLOCK_SEPARATION+1)*numcells_core_poloidal/2 - m_numcells_rpf_poloidal;
           break;
        case LPF:
           index[RADIAL_DIR] = m_numcells_core_radial - m_numcells_pf_radial;
           index[POLOIDAL_DIR] = (2*BLOCK_SEPARATION+1)*numcells_core_poloidal/2;
           break;
        case RPF:
           index[RADIAL_DIR] = m_numcells_core_radial - m_numcells_pf_radial;
           index[POLOIDAL_DIR] = -(2*BLOCK_SEPARATION+1)*numcells_core_poloidal/2 - m_numcells_rpf_poloidal;
           break;
        default:
           MayDay::Error("SingleNullCoordSys::lo_mapped_index(): Invalid block number");
        }
  }
  else {

     int L0 = m_numcells_core_radial;
     int lL1 = m_numcells_mcore_poloidal/2 + m_numcells_lcore_poloidal;
     int rL1 = m_numcells_mcore_poloidal/2 + m_numcells_rcore_poloidal;

     switch ( block )
        {
        case MCORE:
           index[RADIAL_DIR] = 0;
           index[POLOIDAL_DIR] = -m_numcells_mcore_poloidal/2;
           break;
        case LCORE:
           index[RADIAL_DIR] = 0;
           index[POLOIDAL_DIR] = (BLOCK_SEPARATION+1)*lL1;
           break;
        case RCORE:
           index[RADIAL_DIR] = 0;
           index[POLOIDAL_DIR] = -(BLOCK_SEPARATION+1)*rL1 - m_numcells_rcore_poloidal;
           break;
        case LPF:
           index[RADIAL_DIR] = L0 - m_numcells_pf_radial;
           index[POLOIDAL_DIR] = 2*(BLOCK_SEPARATION+1)*lL1;
           break;
        case RPF:
           index[RADIAL_DIR] = L0 - m_numcells_pf_radial;
           index[POLOIDAL_DIR] = -2*(BLOCK_SEPARATION+1)*rL1 - m_numcells_rpf_poloidal;
           break;
        case MCSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*L0;
           index[POLOIDAL_DIR] = -m_numcells_mcore_poloidal/2;
           break;
        case LCSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*L0;
           index[POLOIDAL_DIR] = (BLOCK_SEPARATION+1)*lL1;
           break;
        case RCSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*L0;
           index[POLOIDAL_DIR] = -(BLOCK_SEPARATION+1)*rL1 - m_numcells_rcore_poloidal;
           break;
        case LSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*L0;
           index[POLOIDAL_DIR] = 2*(BLOCK_SEPARATION+1)*lL1;
           break;
        case RSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*L0;
           index[POLOIDAL_DIR] = -2*(BLOCK_SEPARATION+1)*rL1 - m_numcells_rpf_poloidal;
           break;
        default:
           MayDay::Error("SingleNullCoordSys::lo_mapped_index(): Invalid block number");
        }
  }

  return index;
}


IntVect
SingleNullCoordSys::hi_mapped_index(int block) const
{
  IntVect index;

  if ( m_original_eight_blocks ) {

     int numcells_core_poloidal = m_numcells_lcore_poloidal + m_numcells_rcore_poloidal;

     switch ( block )
        {
        case LCORE:
           index[RADIAL_DIR] = m_numcells_core_radial - 1;
           index[POLOIDAL_DIR] = (BLOCK_SEPARATION+1)*numcells_core_poloidal/2 - 1;
           break;
        case RCORE:
           index[RADIAL_DIR] = m_numcells_core_radial - 1;
           index[POLOIDAL_DIR] = -BLOCK_SEPARATION*numcells_core_poloidal/2 - 1;
           break;
        case LCSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*m_numcells_core_radial + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = (BLOCK_SEPARATION+1)*numcells_core_poloidal/2 - 1;
           break;
        case RCSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*m_numcells_core_radial + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = -BLOCK_SEPARATION*numcells_core_poloidal/2 - 1;
           break;
        case LSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*m_numcells_core_radial + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = (2*BLOCK_SEPARATION+1)*numcells_core_poloidal/2 + m_numcells_lpf_poloidal - 1;
           break;
        case RSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*m_numcells_core_radial + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = -(2*BLOCK_SEPARATION+1)*numcells_core_poloidal/2 - 1;  // poloidal domain shifted downward
           break;
        case LPF:
           index[RADIAL_DIR] = m_numcells_core_radial - 1;                                // same hi radial boundary as CORE
           index[POLOIDAL_DIR] = (2*BLOCK_SEPARATION+1)*numcells_core_poloidal/2 + m_numcells_lpf_poloidal - 1;
           break;
        case RPF:
           index[RADIAL_DIR] = m_numcells_core_radial - 1;                           // same hi radial boundary as CORE
           index[POLOIDAL_DIR] = -(2*BLOCK_SEPARATION+1)*numcells_core_poloidal/2 - 1;  // poloidal domain shifted downward
           break;
        default:
           MayDay::Error("SingleNullCoordSys::lo_mapped_index(): Invalid block number");
        }
  }
  else {

     int L0 = m_numcells_core_radial;
     int lL1 = m_numcells_mcore_poloidal/2 + m_numcells_lcore_poloidal;
     int rL1 = m_numcells_mcore_poloidal/2 + m_numcells_rcore_poloidal;

     switch ( block )
        {
        case MCORE:
           index[RADIAL_DIR] = L0 - 1;
           index[POLOIDAL_DIR] = m_numcells_mcore_poloidal/2 - 1;
           break;
        case LCORE:
           index[RADIAL_DIR] = L0 - 1;
           index[POLOIDAL_DIR] = (BLOCK_SEPARATION+1)*lL1 + m_numcells_lcore_poloidal - 1;
           break;
        case RCORE:
           index[RADIAL_DIR] = L0 - 1;
           index[POLOIDAL_DIR] = -(BLOCK_SEPARATION+1)*rL1 - 1;
           break;
        case LPF:
           index[RADIAL_DIR] = L0 - 1;
           index[POLOIDAL_DIR] = 2*(BLOCK_SEPARATION+1)*lL1 + m_numcells_lpf_poloidal - 1; 
           break;
        case RPF:
           index[RADIAL_DIR] = L0 - 1;
           index[POLOIDAL_DIR] = -2*(BLOCK_SEPARATION+1)*rL1 - 1;
           break;
        case MCSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*L0 + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = m_numcells_mcore_poloidal/2 - 1;
           break;
        case LCSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*L0 + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = (BLOCK_SEPARATION+1)*lL1 + m_numcells_lcore_poloidal - 1;
           break;
        case RCSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*L0 + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = -(BLOCK_SEPARATION+1)*rL1 - 1;
           break;
        case LSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*L0 + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = 2*(BLOCK_SEPARATION+1)*lL1 + m_numcells_lpf_poloidal - 1;
           break;
        case RSOL:
           index[RADIAL_DIR] = (BLOCK_SEPARATION+1)*L0 + m_numcells_sol_radial - 1;
           index[POLOIDAL_DIR] = -2*(BLOCK_SEPARATION+1)*rL1 - 1;
           break;
        default:
           MayDay::Error("SingleNullCoordSys::lo_mapped_index(): Invalid block number");
        }
  }

  return index;
}



IntVect
SingleNullCoordSys::getDecomposition(int block) const
{
  IntVect decomp;

  switch (block)
    {
    case MCORE:
      decomp = m_decomp_mcore;
      break;
    case LCORE:
      decomp = m_decomp_lcore;
      break;
    case RCORE:
      decomp = m_decomp_rcore;
      break;
    case MCSOL:
      decomp = m_decomp_mcsol;
      break;
    case LCSOL:
      decomp = m_decomp_lcsol;
      break;
    case RCSOL:
      decomp = m_decomp_rcsol;
      break;
    case LSOL:
      decomp = m_decomp_lsol;
      break;
    case RSOL:
      decomp = m_decomp_rsol;
      break;
    case LPF:
      decomp = m_decomp_lpf;
      break;
    case RPF:
      decomp = m_decomp_rpf;
      break;
    default:
      MayDay::Error("SingleNullCoordSys::getDecomposition(): Invalid block number");
    }

  return decomp;
}



void
SingleNullCoordSys::readGridParams( ParmParse& a_pp )
{
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
   // MCORE block, we assume a 10 block geometry.  Otherwise, we
   // assume the original 8 block geometry
   if ( a_pp.contains("numcells.mcore_poloidal") &&
        a_pp.contains("decomp.mcore.configuration") &&
        a_pp.contains("decomp.mcsol.configuration") ) {

      m_original_eight_blocks = false;

      a_pp.query( "numcells.mcore_poloidal", m_numcells_mcore_poloidal );

      // For now, we require that the MCORE has an even number of cells in the 
      // poloidal direction so that it can be located symmetrically about
      // zero in index space.  This assumption can be removed if needed.
      if ( m_numcells_mcore_poloidal%2 != 0 ) {
         MayDay::Error("Number of mcore cells in poloidal direction must be even");
      }

      a_pp.queryarr( "decomp.mcore.configuration", decomp_mcore, 0, SpaceDim );
      a_pp.queryarr( "decomp.mcsol.configuration", decomp_mcsol, 0, SpaceDim );
   }
   else {
      m_original_eight_blocks = true;
   }

   if ( m_original_eight_blocks && a_pp.contains("numcells.core" ) ) {
      // For backward compatibility with the old way of specifying an eight
      // block single null geometry:
      std::vector<int> numcells_core(SpaceDim);
      a_pp.queryarr( "numcells.core", numcells_core, 0, SpaceDim );

      m_numcells_core_radial = numcells_core[RADIAL_DIR];
      if ( numcells_core[POLOIDAL_DIR]%2 != 0 ) {
         MayDay::Error("Number of core cells in poloidal direction must be even");
      }
      else {
         m_numcells_lcore_poloidal = numcells_core[POLOIDAL_DIR]/2;
         m_numcells_rcore_poloidal = numcells_core[POLOIDAL_DIR]/2;
      }
   }
   else {
      a_pp.query( "numcells.core_radial", m_numcells_core_radial );
      a_pp.query( "numcells.lcore_poloidal", m_numcells_lcore_poloidal );
      a_pp.query( "numcells.rcore_poloidal", m_numcells_rcore_poloidal );
   }

   a_pp.query( "numcells.pf_radial", m_numcells_pf_radial );
   a_pp.query( "numcells.lpf_poloidal", m_numcells_lpf_poloidal );
   a_pp.query( "numcells.rpf_poloidal", m_numcells_rpf_poloidal );
   a_pp.query( "numcells.sol_radial", m_numcells_sol_radial );

   if ( m_original_eight_blocks && a_pp.contains("decomp.core.configuration" ) ) {
      // For backward compatibility with the old way of specifying an eight
      // block single null geometry:
      std::vector<int> decomp_core( SpaceDim );
      a_pp.queryarr( "decomp.core.configuration", decomp_core, 0, SpaceDim );

      for (int dir=0; dir<SpaceDim; ++dir) {
         decomp_lcore[dir] = decomp_core[dir];
         decomp_rcore[dir] = decomp_core[dir];
      }
      if ( decomp_core[POLOIDAL_DIR]%2 != 0 ) {
         MayDay::Error("Core decomposition in poloidal direction must be even");
      }
      else {
         decomp_lcore[POLOIDAL_DIR] /= 2;
         decomp_rcore[POLOIDAL_DIR] /= 2;
      }
   }
   else {
      a_pp.queryarr( "decomp.lcore.configuration", decomp_lcore, 0, SpaceDim );
      a_pp.queryarr( "decomp.rcore.configuration", decomp_rcore, 0, SpaceDim );
   }

   if ( m_original_eight_blocks && a_pp.contains("decomp.csol.configuration" ) ) {
      // For backward compatibility with the old way of specifying an eight
      // block single null geometry:
      std::vector<int> decomp_csol( SpaceDim );
      a_pp.queryarr( "decomp.csol.configuration", decomp_csol, 0, SpaceDim );

      for (int dir=0; dir<SpaceDim; ++dir) {
         decomp_lcsol[dir] = decomp_csol[dir];
         decomp_rcsol[dir] = decomp_csol[dir];
      }
      if ( decomp_csol[POLOIDAL_DIR]%2 != 0 ) {
         MayDay::Error("CSOL decomposition in poloidal direction must be even");
      }
      else {
         decomp_lcsol[POLOIDAL_DIR] /= 2;
         decomp_rcsol[POLOIDAL_DIR] /= 2;
      }
   }
   else {
      a_pp.queryarr( "decomp.lcsol.configuration", decomp_lcsol, 0, SpaceDim );
      a_pp.queryarr( "decomp.rcsol.configuration", decomp_rcsol, 0, SpaceDim );
   }

   a_pp.queryarr( "decomp.lpf.configuration", decomp_lpf, 0, SpaceDim );
   a_pp.queryarr( "decomp.rpf.configuration", decomp_rpf, 0, SpaceDim );
   a_pp.queryarr( "decomp.lsol.configuration", decomp_lsol, 0, SpaceDim );
   a_pp.queryarr( "decomp.rsol.configuration", decomp_rsol, 0, SpaceDim );

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
      MayDay::Error("Radial decomposition does not divide number of cells in MCORE");
   }
   if ( !m_original_eight_blocks && m_numcells_mcore_poloidal % decomp_mcore[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("Poloidal decomposition does not divide number of cells in MCORE");
   }

   if ( m_numcells_core_radial % decomp_lcore[RADIAL_DIR] != 0 ) {
      MayDay::Error("Radial decomposition does not divide number of cells in LCORE");
   }
   if ( m_numcells_lcore_poloidal % decomp_lcore[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("Poloidal decomposition does not divide number of cells in LCORE");
   }

   if ( m_numcells_core_radial % decomp_rcore[RADIAL_DIR] != 0 ) {
      MayDay::Error("Radial decomposition does not divide number of cells in RCORE");
   }
   if ( m_numcells_rcore_poloidal % decomp_rcore[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("Poloidal decomposition does not divide number of cells in RCORE");
   }

   if ( m_numcells_pf_radial % decomp_lpf[RADIAL_DIR] != 0 ) {
      MayDay::Error("Radial decomposition does not divide number of cells in LPF");
   }
   if ( m_numcells_lpf_poloidal % decomp_lpf[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("Poloidal decomposition does not divide number of cells in LPF");
   }

   if ( m_numcells_pf_radial % decomp_rpf[RADIAL_DIR] != 0 ) {
      MayDay::Error("Radial decomposition does not divide number of cells in RPF");
   }
   if ( m_numcells_rpf_poloidal % decomp_rpf[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("Poloidal decomposition does not divide number of cells in RPF");
   }

   if ( !m_original_eight_blocks && m_numcells_sol_radial % decomp_mcsol[RADIAL_DIR] != 0 ) {
      MayDay::Error("Radial decomposition does not divide number of cells in MCSOL");
   }
   if ( !m_original_eight_blocks && m_numcells_mcore_poloidal % decomp_mcsol[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("Poloidal decomposition does not divide number of cells in MCSOL");
   }

   if ( m_numcells_sol_radial % decomp_lcsol[RADIAL_DIR] != 0 ) {
      MayDay::Error("Radial decomposition does not divide number of cells in LCSOL");
   }
   if ( m_numcells_lcore_poloidal % decomp_lcsol[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("Poloidal decomposition does not divide number of cells in LCSOL");
   }

   if ( m_numcells_sol_radial % decomp_rcsol[RADIAL_DIR] != 0 ) {
      MayDay::Error("Radial decomposition does not divide number of cells in RCSOL");
   }
   if ( m_numcells_rcore_poloidal % decomp_rcsol[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("Poloidal decomposition does not divide number of cells in RCSOL");
   }

   if ( m_numcells_sol_radial % decomp_lsol[RADIAL_DIR] != 0 ) {
      MayDay::Error("Radial decomposition does not divide number of cells in LSOL");
   }
   if ( m_numcells_lpf_poloidal % decomp_lsol[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("Poloidal decomposition does not divide number of cells in LSOL");
   }

   if ( m_numcells_sol_radial % decomp_rsol[RADIAL_DIR] != 0 ) {
      MayDay::Error("Radial decomposition does not divide number of cells in RSOL");
   }
   if ( m_numcells_rpf_poloidal % decomp_rsol[POLOIDAL_DIR] != 0 ) {
      MayDay::Error("Poloidal decomposition does not divide number of cells in RSOL");
   }

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
   int num_blocks = m_original_eight_blocks? 8: NUM_SINGLE_NULL_BLOCKS;
   
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
         
         switch ( block_number )
         {
            case LCORE:
               
               // RCORE
               shift_p = -BLOCK_SEPARATION * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                              m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCORE)->domain().domainBox(), shift_p, RCORE, iv_pairs, block_numbers, ghost_ivs );
               
               shift_p = -(BLOCK_SEPARATION+1) * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                                  m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCORE)->domain().domainBox(), shift_p, RCORE, iv_pairs, block_numbers, ghost_ivs );
               
               // LCSOL
               shift_r = BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(LCSOL)->domain().domainBox(), shift_r, LCSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // RCSOL
               shift_p = -BLOCK_SEPARATION * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                              m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCSOL)->domain().domainBox(), shift_r + shift_p, RCSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // LSOL
               shift_p = BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LSOL)->domain().domainBox(), shift_r + shift_p, LSOL, iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case RCORE:
               
               // LCORE
               shift_p = (BLOCK_SEPARATION+1) * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                                 + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCORE)->domain().domainBox(), shift_p, LCORE, iv_pairs, block_numbers, ghost_ivs );
               
               shift_p = BLOCK_SEPARATION * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                             + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCORE)->domain().domainBox(), shift_p, LCORE, iv_pairs, block_numbers, ghost_ivs );
               
               // RCSOL
               shift_r = BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(RCSOL)->domain().domainBox(), shift_r, RCSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // LCSOL
               shift_p = BLOCK_SEPARATION * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                             m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCSOL)->domain().domainBox(), shift_r + shift_p, LCSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // RSOL
               shift_p = -BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RSOL)->domain().domainBox(), shift_r + shift_p, RSOL, iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case LCSOL:
               
               // LSOL
               shift_p = BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LSOL)->domain().domainBox(), shift_p, LSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // RCSOL
               shift_p = -BLOCK_SEPARATION * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                              m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCSOL)->domain().domainBox(), shift_p, RCSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // LCORE
               shift_r = -BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(LCORE)->domain().domainBox(), shift_r, LCORE, iv_pairs, block_numbers, ghost_ivs );
               
               // RCORE
               shift_p = -BLOCK_SEPARATION*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                            m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCORE)->domain().domainBox(), shift_r + shift_p, RCORE, iv_pairs, block_numbers, ghost_ivs );
               
               // LPF
               shift_p = BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LPF)->domain().domainBox(), shift_r + shift_p, LPF, iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case RCSOL:
               
               // RSOL
               shift_p = -BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RSOL)->domain().domainBox(), shift_p, RSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // LCSOL
               shift_p = BLOCK_SEPARATION * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR) +
                                             m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCSOL)->domain().domainBox(), shift_p, LCSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // RCORE
               shift_r = -BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(RCORE)->domain().domainBox(), shift_r, RCORE, iv_pairs, block_numbers, ghost_ivs );
               
               // LCORE
               shift_p = BLOCK_SEPARATION * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                             + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCORE)->domain().domainBox(), shift_r + shift_p, LCORE, iv_pairs, block_numbers, ghost_ivs );
               
               // RPF
               shift_p = -BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RPF)->domain().domainBox(), shift_r + shift_p, RPF, iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case LSOL:
               
               // LCSOL
               shift_p = -BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCSOL)->domain().domainBox(), shift_p, LCSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // LPF
               shift_r = -BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(LPF)->domain().domainBox(), shift_r, LPF, iv_pairs, block_numbers, ghost_ivs );
               
               // LCORE
               shift_p = -BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCORE)->domain().domainBox(), shift_r + shift_p, LCORE, iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case RSOL:
               
               // RCSOL
               shift_p = BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCSOL)->domain().domainBox(), shift_p, RCSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // RPF
               shift_r = -BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(RPF)->domain().domainBox(), shift_r, RPF, iv_pairs, block_numbers, ghost_ivs );
               
               // RCORE
               shift_p = BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCORE)->domain().domainBox(), shift_r + shift_p, RCORE, iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case LPF:
               
               // RPF
               shift_p = -(2*BLOCK_SEPARATION+1)*(m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                                  + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RPF)->domain().domainBox(), shift_p, RPF, iv_pairs, block_numbers, ghost_ivs );
               
               // LSOL
               shift_r = BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(LSOL)->domain().domainBox(), shift_r, LSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // LCSOL
               shift_p = -BLOCK_SEPARATION * m_mappingBlocks[LCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LCSOL)->domain().domainBox(), shift_r + shift_p, LCSOL, iv_pairs, block_numbers, ghost_ivs );
               
               break;
            case RPF:
               
               // LPF
               shift_p = (2*BLOCK_SEPARATION+1) * (m_mappingBlocks[LCORE].size(POLOIDAL_DIR)
                                                   + m_mappingBlocks[RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = grown_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(LPF)->domain().domainBox(), shift_p, LPF, iv_pairs, block_numbers, ghost_ivs );
               
               // RSOL
               shift_r = BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(RADIAL_DIR) * BASISV(RADIAL_DIR);
               rad_shifted_box = grown_box + shift_r;
               accumulateTuples( rad_shifted_box & getCoordSys(RSOL)->domain().domainBox(), shift_r, RSOL, iv_pairs, block_numbers, ghost_ivs );
               
               // RCSOL
               shift_p = BLOCK_SEPARATION * m_mappingBlocks[RCORE].size(POLOIDAL_DIR) * BASISV(POLOIDAL_DIR);
               pol_shifted_box = rad_shifted_box + shift_p;
               accumulateTuples( pol_shifted_box & getCoordSys(RCSOL)->domain().domainBox(), shift_r + shift_p, RCSOL, iv_pairs, block_numbers, ghost_ivs );
               
               break;
            default:
               MayDay::Error("SingleNullCoordSys::blockRemapping(): Unknown block number");
         }
         
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
   
   const SingleNullBlockCoordSysModel& lcore_coord_sys = (const SingleNullBlockCoordSysModel&)(*(getCoordSys(LCORE)));
   const SingleNullBlockCoordSysModel& lsol_coord_sys = (const SingleNullBlockCoordSysModel&)(*(getCoordSys(LSOL)));
   const SingleNullBlockCoordSysModel& lpf_coord_sys = (const SingleNullBlockCoordSysModel&)(*(getCoordSys(LPF)));
   const SingleNullBlockCoordSysModel& rpf_coord_sys = (const SingleNullBlockCoordSysModel&)(*(getCoordSys(RPF)));
   
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
SingleNullCoordSys::toroidalBlockRemapping(IntVect& a_ivDst,
                                           Vector<Real>& a_interpStecil,
                                           const RealVect& a_xiSrc,
                                           const int a_nSrc,
                                           const Side::LoHiSide& a_side) const
{
   MayDay::Error("SingleNullCoordSys::toroidalBlockRemapping is not implemented!!!");
}
#endif
      
#include "NamespaceFooter.H"
