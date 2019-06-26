#include "SNCoreCoordSys.H"
#include "Directions.H"
#include "CONSTANTS.H"
#include "DataArray.H"

#include "NamespaceHeader.H"

const std::string SNCoreCoordSys::pp_name = "sncore";

SNCoreCoordSys::SNCoreCoordSys( ParmParse& a_pp_grid,
                                ParmParse& a_pp_geom )
   : m_num_toroidal_blocks(1)
{
   readGridParams( a_pp_grid );

   m_model_geometry = false;
   a_pp_geom.query( "model_geometry", m_model_geometry);

#if CFG_DIM==3
   a_pp_geom.get( "toroidal_width_over_2pi", m_toroidal_width);
   CH_assert(m_toroidal_width > 0. && m_toroidal_width <= 1.);
   m_toroidal_width *= 2.*Pi;
#endif

   define( a_pp_geom );
}


SNCoreCoordSys::SNCoreCoordSys( ParmParse&                a_pp_geom,
                                const SingleNullCoordSys& a_single_null_coord_sys,
                                const DisjointBoxLayout&  a_single_null_dbl,
                                DisjointBoxLayout&        a_dbl )
   : m_num_toroidal_blocks(1)
{
   m_numcells_core_radial    = a_single_null_coord_sys.numCellsCoreRadial();
   m_numcells_mcore_poloidal = a_single_null_coord_sys.numCellsMcorePoloidal();
   m_numcells_lcore_poloidal = a_single_null_coord_sys.numCellsLcorePoloidal();
   m_numcells_rcore_poloidal = a_single_null_coord_sys.numCellsRcorePoloidal();
   m_numcells_toroidal       = a_single_null_coord_sys.numCellsToroidal();
   
#if CFG_DIM==3
   m_toroidal_width = a_single_null_coord_sys.toroidalWidth();
#endif
   
   m_decomp_mcore = a_single_null_coord_sys.getDecomposition(SNCoreBlockCoordSys::MCORE);
   m_decomp_lcore = a_single_null_coord_sys.getDecomposition(SNCoreBlockCoordSys::LCORE);
   m_decomp_rcore = a_single_null_coord_sys.getDecomposition(SNCoreBlockCoordSys::RCORE);

   if ( a_single_null_coord_sys.numBlocks() == 8 ) {
      m_original_two_blocks = true;
   }
   else {
      m_original_two_blocks = false;
   }
   
   m_model_geometry = a_single_null_coord_sys.isModelGeom();

   define( a_pp_geom );

   // Define a new disjointBoxLayout from the passed in full single null
   // disjointBoxLayout by including only those boxes contained in the
   // core blocks and preserving the processor assignments.

   CH_assert( !a_dbl.isClosed() );

   Vector<Box> boxes;
   Vector<int> procMap;

   const Vector<Box>& single_null_box_array = a_single_null_dbl.boxArray();
   const Vector<int>& single_null_procIDs = a_single_null_dbl.procIDs();

   for (int n=0; n<single_null_box_array.size(); ++n) {
      const Box& box = single_null_box_array[n];
      const int proc = single_null_procIDs[n];

      const int block_number = a_single_null_coord_sys.whichBlock(box);

      if (block_number == SNCoreBlockCoordSys::LCORE || block_number == SNCoreBlockCoordSys::RCORE) {
         boxes.push_back(box);
         procMap.push_back(proc);
      }

      if ( !m_original_two_blocks && block_number == SNCoreBlockCoordSys::MCORE ) {
         boxes.push_back(box);
         procMap.push_back(proc);
      }
   }

   bool is_periodic[SpaceDim];
   for (int dir=0; dir<SpaceDim; ++dir) {
      is_periodic[dir] = a_single_null_dbl.physDomain().isPeriodic(dir);
   }

   Box bounding_box;
   for (int n=0; n<boxes.size(); n++) {
      bounding_box = minBox(bounding_box, boxes[n]);
   }
   ProblemDomain prob_domain = ProblemDomain(bounding_box, is_periodic);

   a_dbl.define( boxes, procMap, prob_domain );
   a_dbl.close();
}



SNCoreCoordSys::~SNCoreCoordSys()
{
   for (int i=0; i<m_coord_vec.size(); ++i) {
      delete m_coord_vec[i];
   }
}


void
SNCoreCoordSys::define( ParmParse& a_pp_geom )
{
   m_num_poloidal_blocks = m_original_two_blocks? 2: SNCoreBlockCoordSys::NUM_SNCORE_BLOCKS;

   int num_blocks = m_num_poloidal_blocks * m_num_toroidal_blocks;

   Vector<Box> domain_boxes(num_blocks);
   for ( int block_number=0; block_number<num_blocks; ++block_number ) {
      domain_boxes[block_number] = Box(lo_mapped_index(block_number), hi_mapped_index(block_number));
   }

   /*
     Compute the computational mesh size.  We arbitrarily set the radial and poloidal width
     of the left core computational domain to 1.  This, together with the input number of
     cells across the left core radial and poloidal widths sets the uniform cell size for
     all of the blocks.  The dimension of the blocks other than the left core are thereby
     determined by the specifed number of cells in that block.
   */
   int radial_width;
   int core_poloidal_width;

   if ( m_original_two_blocks ) {
      radial_width = domain_boxes[SNCoreBlockCoordSys::LCORE].size(RADIAL_DIR);
      core_poloidal_width =  domain_boxes[SNCoreBlockCoordSys::LCORE].size(POLOIDAL_DIR)
                           + domain_boxes[SNCoreBlockCoordSys::RCORE].size(POLOIDAL_DIR);
   }
   else {
      radial_width = domain_boxes[SNCoreBlockCoordSys::MCORE].size(RADIAL_DIR);
      core_poloidal_width = domain_boxes[SNCoreBlockCoordSys::MCORE].size(POLOIDAL_DIR)
                          + domain_boxes[SNCoreBlockCoordSys::LCORE].size(POLOIDAL_DIR)
                          + domain_boxes[SNCoreBlockCoordSys::RCORE].size(POLOIDAL_DIR);
   }

   RealVect dx;
   dx[RADIAL_DIR]   = 1./(double)radial_width;
   dx[POLOIDAL_DIR] = 1./(double)core_poloidal_width;

   bool is_periodic[SpaceDim];
   for (int dir=0; dir<SpaceDim; ++dir) {
      is_periodic[dir] = false;
   }

#if CFG_DIM==3
   is_periodic[TOROIDAL_DIR] = true;
   dx[TOROIDAL_DIR] = m_toroidal_width / (double)m_numcells_toroidal;
#endif

   for ( int block_number=0, toroidal_index=0; toroidal_index<m_num_toroidal_blocks; ++toroidal_index ) {
      for ( int poloidal_index=0; poloidal_index<m_num_poloidal_blocks; ++poloidal_index, ++block_number ) {
         if (!m_model_geometry) {
            SNCoreBlockCoordSys* geom
               = new SNCoreBlockCoordSys( a_pp_geom, ProblemDomain(domain_boxes[block_number], is_periodic), dx, block_number, 0 );
            geom->readFiles(a_pp_geom);
            m_coord_vec.push_back(geom);
         }
         else {
            MagBlockCoordSys* geom = new SNCoreBlockCoordSysModel( a_pp_geom, ProblemDomain(domain_boxes[block_number], is_periodic), dx, block_number, 0 );
            m_coord_vec.push_back(geom);
         }
      }
   }

   defineCoordSystemsAndBoundaries(m_coord_vec);

   setXPointNeighborhood();

   bool correct_field;
   if (a_pp_geom.contains("correct_field")) {
     a_pp_geom.get("correct_field", correct_field);
   }
   else {
      correct_field = false;
   }

   if ( correct_field ) {

      int bc_type = 0;     // Homogeneous Dirichlet
      int bc_value = 0.;
      for (int side=0; side<2; ++side) {
         m_divergence_cleaning_bc.setBCType(0, RADIAL_DIR, side, bc_type);
         m_divergence_cleaning_bc.setBCValue(0, RADIAL_DIR, side, bc_value);
      }
   }

   m_provides_flux = a_pp_geom.contains("field_coefficients_file");
}



void
SNCoreCoordSys::defineBoundaries()
{
   if ( m_original_two_blocks ) {
      defineBoundaries2();
   }
   else {
      defineBoundaries3();
   }
}



void
SNCoreCoordSys::defineBoundaries2()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   int bc_tag = 0;

   for ( int block_number=0, toroidal_index=0; toroidal_index<m_num_toroidal_blocks; ++toroidal_index ) {
      for ( int poloidal_index=0; poloidal_index<m_num_poloidal_blocks; ++poloidal_index, ++block_number ) {

         IndicesTransformation it;
         IntVect shift;
         Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block_number];

         int block_type;
         if (m_model_geometry) {
            const SNCoreBlockCoordSysModel* coord_sys = (const SNCoreBlockCoordSysModel*)m_coordSysVect[block_number];
            block_type = coord_sys->poloidalIndex();
         }
         else {
            const SNCoreBlockCoordSysModel* coord_sys = (const SNCoreBlockCoordSysModel*)m_coordSysVect[block_number];
            block_type = coord_sys->poloidalIndex();
         }

         if( block_type == SNCoreBlockCoordSys::LCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);
            
            shift = -POLOIDAL_BLOCK_SEP*(m_mappingBlocks[SNCoreBlockCoordSys::LCORE].size(POLOIDAL_DIR) +
                                         m_mappingBlocks[SNCoreBlockCoordSys::RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,SNCoreBlockCoordSys::RCORE);
            
            shift = -(POLOIDAL_BLOCK_SEP+1)*(m_mappingBlocks[SNCoreBlockCoordSys::LCORE].size(POLOIDAL_DIR)
                                           + m_mappingBlocks[SNCoreBlockCoordSys::RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,SNCoreBlockCoordSys::RCORE);
         }
         else if( block_type == SNCoreBlockCoordSys::RCORE ) {
            
            blockBoundaries[RADIAL_DIR].define(bc_tag);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);
            
            shift = (POLOIDAL_BLOCK_SEP+1)*(m_mappingBlocks[SNCoreBlockCoordSys::LCORE].size(POLOIDAL_DIR)
                                           + m_mappingBlocks[SNCoreBlockCoordSys::RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,SNCoreBlockCoordSys::LCORE);
            
            shift = POLOIDAL_BLOCK_SEP*(m_mappingBlocks[SNCoreBlockCoordSys::LCORE].size(POLOIDAL_DIR) +
                                        m_mappingBlocks[SNCoreBlockCoordSys::RCORE].size(POLOIDAL_DIR)) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,SNCoreBlockCoordSys::LCORE);
         }
         else {
            MayDay::Error("SNCoreCoordSys::defineBoundaries(): case not implemented");
         }
      }
   }
}



void
SNCoreCoordSys::defineBoundaries3()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   int lL1 = m_numcells_mcore_poloidal/2 + m_numcells_lcore_poloidal;
   int rL1 = m_numcells_mcore_poloidal/2 + m_numcells_rcore_poloidal;

   int bc_tag = 0;

   for ( int block_number=0, toroidal_index=0; toroidal_index<m_num_toroidal_blocks; ++toroidal_index ) {
      for ( int poloidal_index=0; poloidal_index<m_num_poloidal_blocks; ++poloidal_index, ++block_number ) {

         IndicesTransformation it;
         IntVect shift;
         Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block_number];

         int block_type;
         if (m_model_geometry) {
            const SNCoreBlockCoordSysModel* coord_sys = (const SNCoreBlockCoordSysModel*)m_coordSysVect[block_number];
            block_type = coord_sys->poloidalIndex();
         }
         else {
            const SNCoreBlockCoordSysModel* coord_sys = (const SNCoreBlockCoordSysModel*)m_coordSysVect[block_number];
            block_type = coord_sys->poloidalIndex();
         }

         if( block_type == SNCoreBlockCoordSys::MCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -(POLOIDAL_BLOCK_SEP * rL1 + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,SNCoreBlockCoordSys::RCORE);

            shift =  (POLOIDAL_BLOCK_SEP * lL1 + m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,SNCoreBlockCoordSys::LCORE);
         }
         else if( block_type == SNCoreBlockCoordSys::LCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = -(POLOIDAL_BLOCK_SEP * lL1 + m_numcells_lcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,SNCoreBlockCoordSys::MCORE);
            
            shift = -((POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1)
                      + m_numcells_lcore_poloidal + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,SNCoreBlockCoordSys::RCORE);
         }
         else if( block_type == SNCoreBlockCoordSys::RCORE ) {

            blockBoundaries[RADIAL_DIR].define(bc_tag);
            blockBoundaries[RADIAL_DIR + SpaceDim].define(bc_tag);

            shift = ((POLOIDAL_BLOCK_SEP+1) * (lL1 + rL1)
                     + m_numcells_lcore_poloidal + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR].define(it,SNCoreBlockCoordSys::LCORE);

            shift = (POLOIDAL_BLOCK_SEP * rL1 + m_numcells_rcore_poloidal) * BASISV(POLOIDAL_DIR);
            it.defineFromTranslation(shift);
            blockBoundaries[POLOIDAL_DIR + SpaceDim].define(it,SNCoreBlockCoordSys::MCORE);
         }
         else {
            MayDay::Error("SNCoreCoordSys::defineBoundaries(): case not implemented");
         }
      }
   }
}



void
SNCoreCoordSys::blockRemapping(RealVect&       a_xi_valid,
                               int&            a_n_valid,
                               const RealVect& a_xiSrc,
                               int             a_nSrc) const
{
   if ( m_original_two_blocks ) {
      blockRemapping2(a_xi_valid, a_n_valid, a_xiSrc, a_nSrc);
   }
   else {
      blockRemapping3(a_xi_valid, a_n_valid, a_xiSrc, a_nSrc);
   }
}


void
SNCoreCoordSys::blockRemapping2(RealVect&       a_xi_valid,
                                int&            a_n_valid,
                                const RealVect& a_xiSrc,
                                int             a_nSrc) const
{
   const SNCoreBlockCoordSys* src_coord_sys = (SNCoreBlockCoordSys*)getCoordSys(a_nSrc);

   if ( src_coord_sys->isValid(a_xiSrc, true) ) {
      a_n_valid = a_nSrc;
      a_xi_valid = a_xiSrc;
   }
   else {

      RealVect X = src_coord_sys->realCoord(a_xiSrc);
      POL::RealVect X_pol = src_coord_sys->restrictPhysCoordToPoloidal(X);
    
      POL::RealVect XminusXpt = X_pol - m_Xpoint;
      double distance_to_Xpt = XminusXpt.vectorLength();

      if ( distance_to_Xpt < m_xpoint_radius ) {
         a_n_valid = findBlockNearXpt(X_pol);
      }
      else {

         double lo_mapped_radial   = src_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
         double hi_mapped_radial   = src_coord_sys->upperMappedCoordinate(RADIAL_DIR);

         if ( a_xiSrc[RADIAL_DIR] < lo_mapped_radial ||
              a_xiSrc[RADIAL_DIR] > hi_mapped_radial ) {
            a_n_valid = -1;
         }
         else {

            switch (a_nSrc)
               {
               case SNCoreBlockCoordSys::LCORE:
                  a_n_valid = SNCoreBlockCoordSys::RCORE;
                  break;
               case SNCoreBlockCoordSys::RCORE:
                  a_n_valid = SNCoreBlockCoordSys::LCORE;
                  break;
               default:
                  MayDay::Error("SNCoreCoordSys::blockRemapping2(): block number does not correspond to LCORE or RCORE");
               }
         }
      }

      if (a_n_valid != -1) {

         const SNCoreBlockCoordSys* valid_cs = (SNCoreBlockCoordSys*)getCoordSys(a_n_valid);

         a_xi_valid = valid_cs->mappedCoord(X);

         // Check that the point is actually valid
         if ( !valid_cs->isValid(a_xi_valid, true) ) {

            if (distance_to_Xpt < m_xpoint_radius) {

               // If we've arrived here, it's because the separatrix isn't perfectly straight
               // in the X point neighborhood, as was assumed in findBlockNearXpt().  We therefore
               // have to deal with the special cases.

               int n_valid_new = a_n_valid;

               switch ( a_n_valid )
                  {
                  case SNCoreBlockCoordSys::LCORE:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = -1;
                     }
                     break;
                  case SNCoreBlockCoordSys::RCORE:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = -1;
                     }
                     break;
                  default:
                     MayDay::Error("SNCoreCoordSys::blockRemapping(): bad block number");
                  }

               if (n_valid_new != -1 && n_valid_new != a_n_valid) {

                  const SNCoreBlockCoordSys* new_valid_cs = (SNCoreBlockCoordSys*)getCoordSys(n_valid_new);

                  a_xi_valid = new_valid_cs->mappedCoord(X);
                  a_n_valid = n_valid_new;

                  if ( !new_valid_cs->isValid(a_xi_valid, true) ) {
#ifndef MODEL_GEOM
                     printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, new_valid_cs);
#endif
                  }
               }
            }
            else {
               printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, valid_cs);
            }
         }
      }

#if CFG_DIM==3
      while ( a_xi_valid[TOROIDAL_DIR] < 0. ) {
         a_xi_valid[TOROIDAL_DIR] += m_toroidal_width;
      }
      while ( a_xi_valid[TOROIDAL_DIR] > m_toroidal_width ) {
         a_xi_valid[TOROIDAL_DIR] -= m_toroidal_width;
      }
#endif

   }  // end of invalid point
}



void
SNCoreCoordSys::blockRemapping3(RealVect&       a_xi_valid,
                                int&            a_n_valid,
                                const RealVect& a_xiSrc,
                                int             a_nSrc) const
{
   const SNCoreBlockCoordSys* src_coord_sys = (SNCoreBlockCoordSys*)getCoordSys(a_nSrc);

   if ( src_coord_sys->isValid(a_xiSrc, true) ) {
      a_n_valid = a_nSrc;
      a_xi_valid = a_xiSrc;
   }
   else {

      RealVect X = src_coord_sys->realCoord(a_xiSrc);
      POL::RealVect X_pol = src_coord_sys->restrictPhysCoordToPoloidal(X);
      
      POL::RealVect XminusXpt = X_pol - m_Xpoint;
      double distance_to_Xpt = XminusXpt.vectorLength();

      if ( distance_to_Xpt < m_xpoint_radius ) {
         a_n_valid = findBlockNearXpt(X_pol);
      }
      else {

         double lo_mapped_radial   = src_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
         double hi_mapped_radial   = src_coord_sys->upperMappedCoordinate(RADIAL_DIR);

         double lo_mapped_poloidal = src_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);
         double hi_mapped_poloidal = src_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);

         if ( a_xiSrc[RADIAL_DIR] < lo_mapped_radial ||
              a_xiSrc[RADIAL_DIR] > hi_mapped_radial ) {
            a_n_valid = -1;
         }
         else {

            switch (a_nSrc)
               {
               case SNCoreBlockCoordSys::MCORE:
                  if ( a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal ) {
                     a_n_valid = SNCoreBlockCoordSys::LCORE;
                  }
                  else if ( a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal ) {
                     a_n_valid = SNCoreBlockCoordSys::RCORE;
                  }
                  else {
                     MayDay::Error("SNCoreCoordSys::blockRemapping3(): Shouldn't be here");
                  }
                  break;
               case SNCoreBlockCoordSys::LCORE:
                  if ( a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal ) {
                     a_n_valid = SNCoreBlockCoordSys::RCORE;
                  }
                  else if ( a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal ) {
                     a_n_valid = SNCoreBlockCoordSys::MCORE;
                  }
                  else {
                     MayDay::Error("SNCoreCoordSys::blockRemapping3(): Shouldn't be here");
                  }
                  break;
               case SNCoreBlockCoordSys::RCORE:
                  if ( a_xiSrc[POLOIDAL_DIR] > hi_mapped_poloidal ) {
                     a_n_valid = SNCoreBlockCoordSys::MCORE;
                  }
                  else if ( a_xiSrc[POLOIDAL_DIR] < lo_mapped_poloidal ) {
                     a_n_valid = SNCoreBlockCoordSys::LCORE;
                  }
                  else {
                     MayDay::Error("SNCoreCoordSys::blockRemapping3(): Shouldn't be here");
                  }
                  break;
               default:
                  MayDay::Error("SNCoreCoordSys::blockRemapping(): block number does not correspond to LCORE or RCORE");
               }
         }
      }

      if (a_n_valid != -1) {

         const SNCoreBlockCoordSys* valid_cs = (SNCoreBlockCoordSys*)getCoordSys(a_n_valid);

         a_xi_valid = valid_cs->mappedCoord(X);

         // Check that the point is actually valid
         if ( !valid_cs->isValid(a_xi_valid, true) ) {

            if (distance_to_Xpt < m_xpoint_radius) {

               // If we've arrived here, it's because the separatrix isn't perfectly straight
               // in the X point neighborhood, as was assumed in findBlockNearXpt().  We therefore
               // have to deal with the special cases.

               int n_valid_new = a_n_valid;

               switch ( a_n_valid )
                  {
                  case SNCoreBlockCoordSys::MCORE:
                  case SNCoreBlockCoordSys::LCORE:
                  case SNCoreBlockCoordSys::RCORE:
                     if (a_xi_valid[RADIAL_DIR] > valid_cs->upperMappedCoordinate(RADIAL_DIR)) {
                        n_valid_new = -1;
                     }
                     break;
                  default:
                     MayDay::Error("SNCoreCoordSys::blockRemapping(): bad block number");
                  }

               if (n_valid_new != -1 && n_valid_new != a_n_valid) {
                  
                  const SNCoreBlockCoordSys* new_valid_cs = (SNCoreBlockCoordSys*)getCoordSys(n_valid_new);

                  a_xi_valid = new_valid_cs->mappedCoord(X);
                  a_n_valid = n_valid_new;

                  if ( !new_valid_cs->isValid(a_xi_valid, true) ) {
#ifndef MODEL_GEOM
                     printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, new_valid_cs);
#endif
                  }
               }
            }
            else {
               printInvalidPointDiagnostics(a_nSrc, a_xiSrc, X, a_xi_valid, a_n_valid, src_coord_sys, valid_cs);
            }
         }
      }

#if CFG_DIM==3
      while ( a_xi_valid[TOROIDAL_DIR] < 0. ) {
         a_xi_valid[TOROIDAL_DIR] += m_toroidal_width;
      }
      while ( a_xi_valid[TOROIDAL_DIR] > m_toroidal_width ) {
         a_xi_valid[TOROIDAL_DIR] -= m_toroidal_width;
      }
#endif

   }
}



void
SNCoreCoordSys::printInvalidPointDiagnostics( const int                  a_nSrc,
                                              const RealVect&            a_xiSrc,
                                              const RealVect&            a_X,
                                              const RealVect&            a_xi_valid,
                                              const int                  a_n_valid,
                                              const SNCoreBlockCoordSys* a_src_cs,
                                              const SNCoreBlockCoordSys* a_valid_cs ) const
{
#if 0
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
#endif
}



void
SNCoreCoordSys::setXPointNeighborhood()
{
#if CFG_DIM==2
   {
      // Block boundary between LCORE and RCORE

      const SNCoreBlockCoordSys* lcore_coord_sys = (SNCoreBlockCoordSys*)getCoordSys(SNCoreBlockCoordSys::LCORE);

      Box boundary_box = surroundingNodes(bdryHi(lcore_coord_sys->domain().domainBox(), POLOIDAL_DIR, 1), RADIAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      lcore_coord_sys->getNodeCenteredRealCoords(X);
     
      m_Xpoint = RealVect(X(hi,0),X(hi,1));

      // Initialize the ray to lie on the block boundary, pointing away from X point
      m_ray[1] = RealVect(X(lo,0),X(lo,1)) - m_Xpoint;

      // Set the radius of the X point neighborhood, which we currently take to as the
      // length of the block interface between the LCORE and RCORE.  If there is any
      // need to reduce the size of the neighborhood, here is where that should be done.  
      
      m_xpoint_radius = m_ray[1].vectorLength();
      
      // This loop isn't actually necessary unless the preceding m_xpoint_radius
      // definition has somehow been manually overridden to a smaller value
      for( IntVect iv(hi); iv[RADIAL_DIR] >= lo[RADIAL_DIR]; iv[RADIAL_DIR]--) {
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
      // LCORE block boundary on separatrix

      const SNCoreBlockCoordSys* lcore_coord_sys = (SNCoreBlockCoordSys*)getCoordSys(SNCoreBlockCoordSys::LCORE);

      Box boundary_box = surroundingNodes(bdryHi(lcore_coord_sys->domain().domainBox(), RADIAL_DIR, 1), POLOIDAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      lcore_coord_sys->getNodeCenteredRealCoords(X);
     
      // Initialize the ray to lie on the block boundary, pointing away from X point
      m_ray[0] = RealVect(X(lo,0),X(lo,1)) - m_Xpoint;

      for( IntVect iv(hi); iv[POLOIDAL_DIR] >= lo[POLOIDAL_DIR]; iv[POLOIDAL_DIR]--) {
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
      // RCORE block boundary on separtrix

      const SNCoreBlockCoordSys* rcore_coord_sys = (SNCoreBlockCoordSys*)getCoordSys(SNCoreBlockCoordSys::RCORE);

      Box boundary_box = surroundingNodes(bdryHi(rcore_coord_sys->domain().domainBox(), RADIAL_DIR, 1), POLOIDAL_DIR);
      IntVect lo = boundary_box.smallEnd();
      IntVect hi = boundary_box.bigEnd();

      FArrayBox X(boundary_box,SpaceDim);
      rcore_coord_sys->getNodeCenteredRealCoords(X);

      // Initialize the ray to lie on the block boundary, pointing away from X point
      m_ray[2] = RealVect(X(hi,0),X(hi,1)) - m_Xpoint;

      for( IntVect iv(lo); iv[POLOIDAL_DIR] <= hi[POLOIDAL_DIR]; iv[POLOIDAL_DIR]++) {
         RealVect diff = RealVect(X(iv,0),X(iv,1)) - m_Xpoint;

         if (diff.vectorLength() > m_xpoint_radius) {
            break;
         }
         else {
            m_ray[2] = diff;
         }
      }
   }

   for (int ray=0; ray<3; ++ray) {
      m_ray[ray] /= m_ray[ray].vectorLength();
   }

   if (procID()==0) {
      cout << "X point located at " << m_Xpoint << ", neighborhood radius = " << m_xpoint_radius << endl;
   }

   m_ray_angle[0] = - acos(m_ray[0].dotProduct(m_ray[1]));
   m_ray_angle[1] = 0.;
   m_ray_angle[2] =   acos(m_ray[2].dotProduct(m_ray[1]));

#endif

#if CFG_DIM==3

   double cut_frac = 1.;
   double cut_frac2 = 0.15;

   const SNCoreBlockCoordSys* lcore_coord_sys = (SNCoreBlockCoordSys*)getCoordSys(SNCoreBlockCoordSys::LCORE);

   RealVect Xpoint_mapped;
   for (int n=0; n<SpaceDim; ++n) {
      Xpoint_mapped[n] = lcore_coord_sys->upperMappedCoordinate(n);
   }

   m_Xpoint = lcore_coord_sys->restrictPhysCoordToPoloidal(lcore_coord_sys->realCoord(Xpoint_mapped));

   double rad_mapped;
   double tor_mapped;
   double pol_mapped;
   double cut_frac_tmp;

   cut_frac_tmp = ( numBlocks() == 3 )? cut_frac: cut_frac2;

   rad_mapped = lcore_coord_sys->upperMappedCoordinate(RADIAL_DIR);
   tor_mapped = lcore_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR);
   pol_mapped = (1.-cut_frac_tmp)*lcore_coord_sys->upperMappedCoordinate(POLOIDAL_DIR)
      + cut_frac_tmp*lcore_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR);

   RealVect mapped_coord = RealVect(rad_mapped, tor_mapped, pol_mapped);

   m_ray[0] = lcore_coord_sys->restrictPhysCoordToPoloidal(lcore_coord_sys->realCoord(mapped_coord)) - m_Xpoint;
   
   rad_mapped = (1.-cut_frac_tmp)*lcore_coord_sys->upperMappedCoordinate(RADIAL_DIR)
      + cut_frac_tmp*lcore_coord_sys->lowerMappedCoordinate(RADIAL_DIR);
   tor_mapped = lcore_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR);
   pol_mapped = lcore_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);

   mapped_coord = RealVect(rad_mapped, tor_mapped, pol_mapped);

   m_ray[1] = lcore_coord_sys->restrictPhysCoordToPoloidal(lcore_coord_sys->realCoord(mapped_coord)) - m_Xpoint;
   
   const SNCoreBlockCoordSys* rcore_coord_sys = (SNCoreBlockCoordSys*)getCoordSys(SNCoreBlockCoordSys::RCORE);

   rad_mapped = rcore_coord_sys->upperMappedCoordinate(RADIAL_DIR);
   tor_mapped = rcore_coord_sys->lowerMappedCoordinate(TOROIDAL_DIR);
   pol_mapped = (1.-cut_frac_tmp)*rcore_coord_sys->lowerMappedCoordinate(POLOIDAL_DIR)
      + cut_frac_tmp*rcore_coord_sys->upperMappedCoordinate(POLOIDAL_DIR);

   mapped_coord = RealVect(rad_mapped, tor_mapped, pol_mapped);

   m_ray[2] = rcore_coord_sys->restrictPhysCoordToPoloidal(rcore_coord_sys->realCoord(mapped_coord)) - m_Xpoint;

   m_xpoint_radius = m_ray[0].vectorLength();
   m_ray[0] /= m_xpoint_radius;
   for (int block=1; block<=2; ++block) {
      double ray_length = m_ray[block].vectorLength();
      if (ray_length < m_xpoint_radius) m_xpoint_radius = ray_length;
      m_ray[block] /= ray_length;
   }

   if (procID()==0) {
      cout << "X point neighborhood radius = " << m_xpoint_radius << endl;
   }

   m_ray_angle[0] = - acos(m_ray[0].dotProduct(m_ray[1]));
   m_ray_angle[1] = 0.;
   m_ray_angle[2] =   acos(m_ray[2].dotProduct(m_ray[1]));

#endif

   m_sector_block[0] = SNCoreBlockCoordSys::LCORE;
   m_sector_block[1] = SNCoreBlockCoordSys::RCORE;

}


int
SNCoreCoordSys::findBlockNearXpt( const POL::RealVect& a_X_pol ) const
{
   POL::RealVect local_X_pol = a_X_pol - m_Xpoint;

   POL::RealVect ray1_normal(m_ray[1][1],-m_ray[1][0]);  // m_ray[1] rotated 90 degrees clockwise

   double angle = acos( local_X_pol.dotProduct(m_ray[1]) / local_X_pol.vectorLength() );

   if (local_X_pol.dotProduct(ray1_normal) < 0.) {
      angle *= -1.;
   }

   return (angle < m_ray_angle[1])? m_sector_block[0]: m_sector_block[1];
}


void
SNCoreCoordSys::defineEllipticOpBC( EllipticOpBC& a_bc ) const
{
   SNCoreEllipticOpBC& bc = (SNCoreEllipticOpBC&)a_bc;

   // Inner radial core boundary
   RefCountedPtr<GridFunction> radial_inner_core_function = bc.getRadialInnerFunction();
   if ( !radial_inner_core_function.isNull() ) {
      if ( typeid(*radial_inner_core_function) == typeid(DataArray) ) {
         RefCountedPtr<DataArray> data_array = (RefCountedPtr<DataArray>)radial_inner_core_function;
         CH_assert( !data_array.isNull() );

         const Box& lcore_domain_box = getCoordSys(SNCoreBlockCoordSys::LCORE)->domain().domainBox();
         const Box& rcore_domain_box = getCoordSys(SNCoreBlockCoordSys::RCORE)->domain().domainBox();

         Box fab_box = minBox( bdryLo(lcore_domain_box, RADIAL_DIR, 1),
                               bdryLo(rcore_domain_box, RADIAL_DIR, 1) );

         FArrayBox fab(fab_box, 1);
         data_array->setData(fab, false);
      }
   }

   // Outer radial core boundary
   RefCountedPtr<GridFunction> radial_outer_core_function = bc.getRadialOuterFunction();
   if ( !radial_outer_core_function.isNull() ) {
      if ( typeid(*radial_outer_core_function) == typeid(DataArray) ) {
         RefCountedPtr<DataArray> data_array = (RefCountedPtr<DataArray>)radial_outer_core_function;
         CH_assert( !data_array.isNull() );

         const Box& lcore_domain_box = getCoordSys(SNCoreBlockCoordSys::LCORE)->domain().domainBox();
         const Box& rcore_domain_box = getCoordSys(SNCoreBlockCoordSys::RCORE)->domain().domainBox();

         Box fab_box = minBox( bdryHi(lcore_domain_box, RADIAL_DIR, 1),
                               bdryHi(rcore_domain_box, RADIAL_DIR, 1) );

         FArrayBox fab(fab_box, 1);
         data_array->setData(fab, false);
      }
   }
}



IntVect
SNCoreCoordSys::lo_mapped_index( int a_block_number ) const
{
   IntVect index;

   if ( m_original_two_blocks ) {

      int numcells_core_poloidal = m_numcells_lcore_poloidal + m_numcells_rcore_poloidal;

      switch ( a_block_number )
         {
         case SNCoreBlockCoordSys::LCORE:
            index[RADIAL_DIR] = 0;
            index[POLOIDAL_DIR] = POLOIDAL_BLOCK_SEP*numcells_core_poloidal/2;
            break;
         case SNCoreBlockCoordSys::RCORE:
            index[RADIAL_DIR] = 0;
            index[POLOIDAL_DIR] = -(POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2;
            break;
         default:
            MayDay::Error("SNCoreCoordSys::lo_mapped_index(): Invalid block number");
         }
   }
   else {

      int lL1 = m_numcells_mcore_poloidal/2 + m_numcells_lcore_poloidal;
      int rL1 = m_numcells_mcore_poloidal/2 + m_numcells_rcore_poloidal;

      switch ( a_block_number )
         {
         case SNCoreBlockCoordSys::MCORE:
            index[RADIAL_DIR] = 0;
            index[POLOIDAL_DIR] = -m_numcells_mcore_poloidal/2;
            break;
         case SNCoreBlockCoordSys::LCORE:
            index[RADIAL_DIR] = 0;
            index[POLOIDAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*lL1;
            break;
         case SNCoreBlockCoordSys::RCORE:
            index[RADIAL_DIR] = 0;
            index[POLOIDAL_DIR] = -(POLOIDAL_BLOCK_SEP+1)*rL1 - m_numcells_rcore_poloidal;
            break;
         default:
            MayDay::Error("SNCoreCoordSys::lo_mapped_index(): Invalid block number");
         }
   }

#if CFG_DIM==3
   int num_block_toroidal_cells = m_numcells_toroidal / m_num_toroidal_blocks;
  
   index[TOROIDAL_DIR] = toroidalBlockNumber(a_block_number) * (num_block_toroidal_cells + TOROIDAL_BLOCK_SEP);
#endif

  return index;
}



IntVect
SNCoreCoordSys::hi_mapped_index( int a_block_number ) const
{
   IntVect index;

   if ( m_original_two_blocks ) {

      int numcells_core_poloidal = m_numcells_lcore_poloidal + m_numcells_rcore_poloidal;

      switch ( a_block_number )
         {
         case SNCoreBlockCoordSys::LCORE:
            index[RADIAL_DIR] = m_numcells_core_radial - 1;
            index[POLOIDAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*numcells_core_poloidal/2 - 1;
            break;
         case SNCoreBlockCoordSys::RCORE:
            index[RADIAL_DIR] = m_numcells_core_radial - 1;
            index[POLOIDAL_DIR] = -POLOIDAL_BLOCK_SEP*numcells_core_poloidal/2 - 1;
            break;
         default:
            MayDay::Error("SNCoreCoordSys::lo_mapped_index(): Invalid block number");
         }
   }
   else {

      int L0 = m_numcells_core_radial;
      int lL1 = m_numcells_mcore_poloidal/2 + m_numcells_lcore_poloidal;
      int rL1 = m_numcells_mcore_poloidal/2 + m_numcells_rcore_poloidal;

      switch ( a_block_number )
         {
         case SNCoreBlockCoordSys::MCORE:
            index[RADIAL_DIR] = L0 - 1;
            index[POLOIDAL_DIR] = m_numcells_mcore_poloidal/2 - 1;
            break;
         case SNCoreBlockCoordSys::LCORE:
            index[RADIAL_DIR] = L0 - 1;
            index[POLOIDAL_DIR] = (POLOIDAL_BLOCK_SEP+1)*lL1 + m_numcells_lcore_poloidal - 1;
            break;
         case SNCoreBlockCoordSys::RCORE:
            index[RADIAL_DIR] = L0 - 1;
            index[POLOIDAL_DIR] = -(POLOIDAL_BLOCK_SEP+1)*rL1 - 1;
            break;
         default:
            MayDay::Error("SNCoreCoordSys::lo_mapped_index(): Invalid block number");
         }
   }

#if CFG_DIM==3
   int num_block_toroidal_cells = m_numcells_toroidal / m_num_toroidal_blocks;
  
   index[TOROIDAL_DIR] = toroidalBlockNumber(a_block_number) * (num_block_toroidal_cells + TOROIDAL_BLOCK_SEP)
                         + num_block_toroidal_cells - 1;
#endif

  return index;
}



IntVect
SNCoreCoordSys::getDecomposition( int a_block ) const
{
  IntVect decomp;

  switch (a_block)
    {
    case SNCoreBlockCoordSys::MCORE:
      decomp = m_decomp_mcore;
      break;
    case SNCoreBlockCoordSys::LCORE:
      decomp = m_decomp_lcore;
      break;
    case SNCoreBlockCoordSys::RCORE:
      decomp = m_decomp_rcore;
      break;
    default:
      MayDay::Error("SNCoreCoordSys::getDecomposition(): Invalid block number");
    }

  return decomp;
}



void
SNCoreCoordSys::readGridParams( ParmParse& a_pp )
{
   std::vector<int> decomp_mcore( SpaceDim );
   std::vector<int> decomp_lcore( SpaceDim );
   std::vector<int> decomp_rcore( SpaceDim );

   // If the input specifies the mesh and decomposition for the
   // MCORE block, we assume a 10 block geometry.  Otherwise, we
   // assume the original 8 block geometry
   if ( a_pp.contains("numcells.mcore_poloidal") &&
        a_pp.contains("decomp.mcore.configuration") ) {

      m_original_two_blocks = false;

      a_pp.query( "numcells.mcore_poloidal", m_numcells_mcore_poloidal );

      // For now, we require that the MCORE has an even number of cells in the 
      // poloidal direction so that it can be located symmetrically about
      // zero in index space.  This assumption can be removed if needed.
      if ( m_numcells_mcore_poloidal%2 != 0 ) {
         MayDay::Error("Number of mcore cells in poloidal direction must be even");
      }

      a_pp.queryarr( "decomp.mcore.configuration", decomp_mcore, 0, SpaceDim );
   }
   else {
      m_original_two_blocks = true;
   }

   if ( m_original_two_blocks && a_pp.contains("numcells.core" ) ) {
      // For backward compatibility with the old way of specifying a two
      // block single null core geometry:
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

#if CFG_DIM==3
      m_numcells_toroidal = numcells_core[TOROIDAL_DIR];
#endif
   }
   else {
      a_pp.query( "numcells.core_radial", m_numcells_core_radial );
      a_pp.query( "numcells.lcore_poloidal", m_numcells_lcore_poloidal );
      a_pp.query( "numcells.rcore_poloidal", m_numcells_rcore_poloidal );

#if CFG_DIM==3
      a_pp.query( "numcells.toroidal", m_numcells_toroidal );
#endif
   }

   if ( m_original_two_blocks && a_pp.contains("decomp.core.configuration" ) ) {
      // For backward compatibility with the old way of specifying a two
      // block single null core geometry:
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

   if ( procID() == 0 ) {
      cout << endl << "SNCore grid and domain decomposition parameters:" << endl << endl;

      cout << "numcells.core_radial = " << m_numcells_core_radial << endl;
      if ( !m_original_two_blocks ) cout << "numcells.mcore_poloidal = " << m_numcells_mcore_poloidal << endl;
      cout << "numcells.lcore_poloidal = " << m_numcells_lcore_poloidal << endl;
      cout << "numcells.rcore_poloidal = " << m_numcells_rcore_poloidal << endl;

      if ( !m_original_two_blocks ) {
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
   }

   // Check decomposition divisibilities

   if ( !m_original_two_blocks ) {
      if ( m_numcells_core_radial % decomp_mcore[RADIAL_DIR] != 0 ) {
         MayDay::Error("Radial decomposition does not divide number of cells in MCORE");
      }
      if ( m_numcells_mcore_poloidal % decomp_mcore[POLOIDAL_DIR] != 0 ) {
         MayDay::Error("Poloidal decomposition does not divide number of cells in MCORE");
      }
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

   for (int dir=0; dir<SpaceDim; ++dir) {
      m_decomp_mcore[dir] = decomp_mcore[dir];
      m_decomp_lcore[dir] = decomp_lcore[dir];
      m_decomp_rcore[dir] = decomp_rcore[dir];
   }
}

#if CFG_DIM == 3
void
SNCoreCoordSys::toroidalBlockRemapping(IntVect& a_ivDst,
                                       Vector<Real>& a_interpStecil,
                                       const RealVect& a_xiSrc,
                                       const int a_nSrc,
                                       const Side::LoHiSide& a_side) const
{
   MayDay::Error("SNCoreCoordSys::toroidalBlockRemapping is not implemented!!!");
}
#endif


Vector<RealVect>
SNCoreCoordSys::displacements(const Vector<RealVect>&   a_dstCoords,
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


#include "NamespaceFooter.H"
