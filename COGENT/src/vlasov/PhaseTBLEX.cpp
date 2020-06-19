#include "PhaseTBLEX.H"
#include "Directions.H"
#include "MBStencilIterator.H"

#define NO_SINGLE_NULL

#ifndef NO_SINGLE_NULL
#include "MagFluxAlignedMBLEXCenter.H"
#endif

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "MagCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM


#include "NamespaceHeader.H"


PhaseTBLEX::PhaseTBLEX( const PhaseGeom&  a_geometry,
                        int               a_ghosts,
                        int               a_order )
   : m_geometry(a_geometry)
{
#if CFG_DIM==3
   const CFG::MagGeom& mag_geom = a_geometry.magGeom();
   
   const RefCountedPtr<CFG::MagCoordSys> mag_coord_sys_ptr = mag_geom.getCoordSys();
#ifdef NO_SINGLE_NULL
   if ( mag_coord_sys_ptr->type() != "LogicallyRectangular" ) {
      MayDay::Error("PhaseTBLEX is currently only implemented for logically rectangular");
   }
#else
   if ( mag_coord_sys_ptr->type() != "LogicallyRectangular" && mag_coord_sys_ptr->type() != "SingleNull" ) {
      MayDay::Error("PhaseTBLEX is currently only implemented for logically rectangular or single null geometry");
   }
#endif

   const DisjointBoxLayout& grids = a_geometry.gridsFull();
   m_ghostCells.define(grids);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_ghostCells[dit] = extraBlockGhosts(grids[dit], a_ghosts);
   }

   m_stencils.define(grids);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const IntVectSet& ghostCellsIVS = m_ghostCells[dit];
      m_stencils[dit] = RefCountedPtr< IVSFAB<MBStencil> >(new IVSFAB<MBStencil>(ghostCellsIVS, 1));
   }

#ifndef NO_SINGLE_NULL
   MagFluxAlignedMBLEXCenter* mfa_mblx_ptr = NULL;
   if ( mag_coord_sys_ptr->type() == "SingleNull" ) {
      mfa_mblx_ptr = new MagFluxAlignedMBLEXCenter;
      mfa_mblx_ptr->define(&a_geometry, a_ghosts, a_order);
   }
#endif

   const LevelData<FArrayBox>& remapped_index = a_geometry.getShearedRemappedIndex();
   const LevelData<FArrayBox>& interp_stencil = a_geometry.getShearedInterpStencil();
   const LevelData<FArrayBox>& sheared_interp_stencil_offsets = a_geometry.getShearedInterpStencilOffsets();

   int order = mag_geom.shearedInterpOrder();
   int num_blocks = mag_coord_sys_ptr->numBlocks();

   const RefCountedPtr<PhaseCoordSys>& coord_sys = a_geometry.phaseCoordSysPtr();
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box box = grids[dit];
      const int block_num = coord_sys->whichBlock(box);
      const ProblemDomain& domain = a_geometry.getBlockCoordSys(box).domain();
      const Box& domain_box = domain.domainBox();
      const FArrayBox& this_remapped_index = remapped_index[dit];
      const FArrayBox& this_interp_stencil = interp_stencil[dit];
      const IntVectSet& ghost_cells = m_ghostCells[dit];
      IVSFAB<MBStencil>& this_stencil_fab = *(m_stencils[dit]);
      
      for (IVSIterator iv_it(ghost_cells); iv_it.ok(); ++iv_it) {
         IntVect iv = iv_it();

         IntVect iv_injected = iv;
         for (int dim=CFG_DIM; dim<SpaceDim; ++dim) {
            iv_injected[dim] = box.smallEnd(dim);
         }

         IntVect iv_remapped_center;
         RefCountedPtr<Vector<MBStencilElement> > stencil_vect_ptr 
            = RefCountedPtr<Vector<MBStencilElement> >(new Vector<MBStencilElement>);

         bool iv_on_lo_toroidal_bndry = iv[TOROIDAL_DIR] < domain_box.smallEnd(TOROIDAL_DIR);
         bool iv_on_hi_toroidal_bndry = iv[TOROIDAL_DIR] > domain_box.bigEnd(TOROIDAL_DIR);

         if ( iv_on_lo_toroidal_bndry || iv_on_hi_toroidal_bndry ) {

            for (int dim=0; dim<CFG_DIM; ++dim) {
               iv_remapped_center[dim] = this_remapped_index(iv_injected,dim);
            }
            for (int dim=CFG_DIM; dim<SpaceDim; ++dim) {
               iv_remapped_center[dim] = iv[dim];
            }
            
            if ( mag_coord_sys_ptr->type() == "LogicallyRectangular" ) {
               CH_assert(domain.isPeriodic(POLOIDAL_DIR));

               int other_block_num;
               if (iv[TOROIDAL_DIR] > domain_box.bigEnd(TOROIDAL_DIR)) {
                  other_block_num = (block_num==num_blocks-1)? 0: block_num+1;
               }
               else {
                  other_block_num = (block_num==0)? num_blocks-1: block_num-1;
               }

               for (int n=0; n<order+1; ++n) {
                  IntVect iv_tmp(iv_remapped_center);
                  iv_tmp[POLOIDAL_DIR] += n - order/2;
                  double weight = this_interp_stencil(iv_injected,n);

                  // Use the assumed periodicity of the poloidal direction to ensure that
                  // iv_tmp is valid.
                  if ( iv_tmp[POLOIDAL_DIR] < domain_box.smallEnd(POLOIDAL_DIR)) {
                     iv_tmp[POLOIDAL_DIR] += domain_box.size(POLOIDAL_DIR);
                  }
                  else if ( iv_tmp[POLOIDAL_DIR] > domain_box.bigEnd(POLOIDAL_DIR)) {
                     iv_tmp[POLOIDAL_DIR] -= domain_box.size(POLOIDAL_DIR);
                  }

                  stencil_vect_ptr->push_back(MBStencilElement(iv_tmp, other_block_num, weight));
               }
            }
            else if ( mag_coord_sys_ptr->type() == "SingleNull" ) {

#ifndef NO_SINGLE_NULL
               int iv_remapped_valid_block = this_remapped_index(iv,SpaceDim);
               CH_assert(iv_remapped_valid_block >= 0);

               double sawtooth_interp_fac = this_interp_stencil(iv,order+1);

               if ( sawtooth_interp_fac < 1. ) {

                  for (int n=0; n<order+1; ++n) {

                     int offset = (int)(sheared_interp_stencil_offsets)[dit](iv_injected,n);

                     // Find the valid index and block corresponding to the offset iv_remapped_index

                     IntVect iv_offset = iv_remapped_center;
                     iv_offset[POLOIDAL_DIR] += offset;

                     Vector<int> valid_blocks;
                     Vector<IntVect> valid_cells;
                     coord_sys->validBlocksCells(valid_blocks, valid_cells, iv_remapped_valid_block, iv_offset);

                     CH_assert (valid_blocks.size() == 1 && valid_cells.size() == 1);
                     int valid_block = valid_blocks[0];
                     IntVect iv_valid = valid_cells[0];
                     CH_assert(valid_block >= 0);

                     double weight = this_interp_stencil(iv,n);
                     if ( sawtooth_interp_fac >= 0. ) weight *= (1. - sawtooth_interp_fac);
                     
                     // Construct the stencil entry and append it
                     stencil_vect_ptr->push_back(MBStencilElement(iv_valid, valid_block, weight));
                  }
               }
#endif

            }
            else {
               MayDay::Error("PhaseTBLEX::PhaseTBLEX(): unknown geometry type");
            }
         }
         else if ( mag_coord_sys_ptr->type() == "LogicallyRectangular" 
                   && (iv[POLOIDAL_DIR] < domain_box.smallEnd(POLOIDAL_DIR) ||
                       iv[POLOIDAL_DIR] > domain_box.bigEnd(POLOIDAL_DIR)) ) {

            iv_remapped_center = iv;
            if (iv[POLOIDAL_DIR] < domain_box.smallEnd(POLOIDAL_DIR)) {
               iv_remapped_center[POLOIDAL_DIR] += domain_box.size(POLOIDAL_DIR);
            }
            else {
               iv_remapped_center[POLOIDAL_DIR] -= domain_box.size(POLOIDAL_DIR);
            }

            stencil_vect_ptr->push_back(MBStencilElement(iv_remapped_center, block_num, 1.));
         }
         else if ( mag_coord_sys_ptr->type() == "SingleNull" ) {

#ifndef NO_SINGLE_NULL
            const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencil = mfa_mblx_ptr->stencils();
            const IVSFAB<MBStencil>& stencil_fab = *stencil[dit];
            
            MBStencilIterator stit(stencil_fab(iv, 0));
            for (stit.begin(); stit.ok(); ++stit) {
               stencil_vect_ptr->push_back(stit());
            }
#endif
         }
         else {
            MayDay::Error("PhaseTBLEX::PhaseTBLEX(): unknown geometry type");
         }

         this_stencil_fab(iv,0) = MBStencil(stencil_vect_ptr);
      }
   }

#ifndef NO_SINGLE_NULL
   if ( mfa_mblx_ptr ) delete mfa_mblx_ptr;
#endif

#else
   MayDay::Error("Attempted to construct a PhaseTBLEX object when not 3D");
#endif
}


void PhaseTBLEX::interpGhosts( LevelData<FArrayBox>& a_data ) const
{
#if CFG_DIM==3
   m_geometry.interpolateFromShearedGhosts(a_data);
#endif
}


IntVectSet PhaseTBLEX::extraBlockGhosts( const Box&  a_box,
                                         int         a_ghosts ) const
{
   // This is a slight modification of MultiBlockUtil::extraBlockGhosts()
   IntVect ghostLayerVect = a_ghosts*IntVect::Unit;
   ghostLayerVect[VPARALLEL_DIR] = 0;  // No extrablock ghosts in VPARALLEL_DIR
   ghostLayerVect[MU_DIR] = 0;         // No extrablock ghosts in MU_DIR

   Box grownBox = grow(a_box, ghostLayerVect);
  
   const RefCountedPtr<PhaseCoordSys> coord_sys_ptr = m_geometry.phaseCoordSysPtr();
   Vector< Tuple<BlockBoundary, 2*SpaceDim> > boundaries = coord_sys_ptr->boundaries();
   int block_number = coord_sys_ptr->whichBlock(a_box);
   const PhaseBlockCoordSys& block_coord_sys 
      = (const PhaseBlockCoordSys&)(*coord_sys_ptr->getCoordSys(a_box));
   const Box& block_box = block_coord_sys.domain().domainBox();

   IntVectSet ivsReturn; // empty

   if ( !block_box.contains(grownBox) ) // otherwise, empty
      {
         // We'll set ivs to extra-block ghost cells of a_box.
         // Start with ivs being the full grownBox.
         DenseIntVectSet ivs = DenseIntVectSet(grownBox, true);
         // Note block_box includes all of a_box
         // so after this, ivs contains ghost cells only.
         ivs -= block_box;
         if ( !ivs.isEmpty() ) {
            // Now remove ghost cells that are outside domain.
            ivs.recalcMinBox();
            for (int idir = 0; idir < CFG_DIM; idir++) {
               // Check idir Lo face.
               int endBlockLo = block_box.smallEnd(idir);
               if (grownBox.smallEnd(idir) < endBlockLo)
                  { // Remove ghost cells from idir Lo face.
                     if (boundaries[block_number][idir].isDomainBoundary())
                        { // Remove ghost cells beyond this face.
                           Box grownBoxFace(grownBox);
                           grownBoxFace.setBig(idir, endBlockLo-1);
                           ivs -= grownBoxFace;
                        }
                  }
               // Check idir Hi face.
               int endBlockHi = block_box.bigEnd(idir);
               if (grownBox.bigEnd(idir) > endBlockHi)
                  { // Remove ghost cells from idir Hi face.
                     if (boundaries[block_number][idir + SpaceDim].isDomainBoundary())
                        { // Remove ghost cells beyond this face.
                           Box grownBoxFace(grownBox);
                           grownBoxFace.setSmall(idir, endBlockHi+1);
                           ivs -= grownBoxFace;
                        }
                  }
            }
            ivs.recalcMinBox();
            // convert from DenseIntVectSet.
            ivsReturn = IntVectSet(ivs);
         }
      }

   return ivsReturn;
}

#include "NamespaceFooter.H"


