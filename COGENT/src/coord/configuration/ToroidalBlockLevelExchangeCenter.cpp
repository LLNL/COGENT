#include "ToroidalBlockLevelExchangeCenter.H"
#include "Directions.H"

#include "NamespaceHeader.H"


ToroidalBlockLevelExchangeCenter::ToroidalBlockLevelExchangeCenter( const MagGeom&  a_geometry,
                                                                    int             a_ghosts,
                                                                    int             a_order )
   : m_geometry(a_geometry)
{
#if CFG_DIM==3
   
   const RefCountedPtr<MagCoordSys> coord_sys_ptr = a_geometry.getCoordSys();
   if ( coord_sys_ptr->type() != "LogicallyRectangular" ) {
      MayDay::Error("ToroidalBlockLevelExchangeCenter is currently only implemented for logically rectangular geometry");
   }

   const DisjointBoxLayout& grids = a_geometry.gridsFull();
   m_ghostCells.define(grids);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_ghostCells[dit] = extraBlockGhosts(grids[dit], a_order/2);
   }

   m_stencils.define(grids);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const IntVectSet& ghostCellsIVS = m_ghostCells[dit];
      m_stencils[dit] = RefCountedPtr< IVSFAB<MBStencil> >
         (new IVSFAB<MBStencil>(ghostCellsIVS, 1));
   }

   const LevelData<FArrayBox>& remapped_index = a_geometry.getShearedRemappedIndex();
   const LevelData<FArrayBox>& interp_stencil = a_geometry.getShearedInterpStencil();
   int order = a_geometry.shearedInterpOrder();
   int num_blocks = coord_sys_ptr->numBlocks();
   
   const RefCountedPtr<MagCoordSys>& coord_sys = a_geometry.getCoordSys();
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box box = grids[dit];
      const int block_num = coord_sys->whichBlock(box);
      const Box& domain_box = a_geometry.getBlockCoordSys(block_num).domain().domainBox();
      const FArrayBox& this_remapped_index = remapped_index[dit];
      const FArrayBox& this_interp_stencil = interp_stencil[dit];
      const IntVectSet& ghost_cells = m_ghostCells[dit];
      IVSFAB<MBStencil>& this_stencil_fab = *(m_stencils[dit]);
      
      for (IVSIterator iv_it(ghost_cells); iv_it.ok(); ++iv_it) {
         IntVect iv = iv_it();
         
         IntVect iv_other_block_center;
         int other_block_num;
         RefCountedPtr<Vector<MBStencilElement> > stencil_vect_ptr 
            = RefCountedPtr<Vector<MBStencilElement> >(new Vector<MBStencilElement>);

         if (iv[TOROIDAL_DIR] < domain_box.smallEnd(TOROIDAL_DIR) ||
             iv[TOROIDAL_DIR] > domain_box.bigEnd(TOROIDAL_DIR) ) {

            if (iv[TOROIDAL_DIR] > domain_box.bigEnd(TOROIDAL_DIR)) {
               other_block_num = (block_num==num_blocks-1)? 0: block_num+1;
            }
            else {
               other_block_num = (block_num==0)? num_blocks-1: block_num-1;
            }

            for (int dim=0; dim<SpaceDim; ++dim) {
               iv_other_block_center[dim] = this_remapped_index(iv,dim);
            }
            
            for (int n=0; n<order+1; ++n) {
               IntVect iv_tmp(iv_other_block_center);
               iv_tmp[POLOIDAL_DIR] += n - order/2;
               double weight = this_interp_stencil(iv,n);

               if ( iv_tmp[POLOIDAL_DIR] < domain_box.smallEnd(POLOIDAL_DIR)) {
                  iv_tmp[POLOIDAL_DIR] += domain_box.size(POLOIDAL_DIR);
               }
               else if ( iv_tmp[POLOIDAL_DIR] > domain_box.bigEnd(POLOIDAL_DIR)) {
                  iv_tmp[POLOIDAL_DIR] -= domain_box.size(POLOIDAL_DIR);
               }

               stencil_vect_ptr->push_back(MBStencilElement(iv_tmp, other_block_num, weight));
            }
         }
         else if (iv[POLOIDAL_DIR] < domain_box.smallEnd(POLOIDAL_DIR) ||
                  iv[POLOIDAL_DIR] > domain_box.bigEnd(POLOIDAL_DIR) ) {

            other_block_num = block_num;

            iv_other_block_center = iv;
            if (iv[POLOIDAL_DIR] < domain_box.smallEnd(POLOIDAL_DIR)) {
               iv_other_block_center[POLOIDAL_DIR] += domain_box.size(POLOIDAL_DIR);
            }
            else {
               iv_other_block_center[POLOIDAL_DIR] -= domain_box.size(POLOIDAL_DIR);
            }

            stencil_vect_ptr->push_back(MBStencilElement(iv_other_block_center, other_block_num, 1.));
         }

         this_stencil_fab(iv,0) = MBStencil(stencil_vect_ptr);
      }
   }

   m_isDefined = true;
#else
   MayDay::Error("Attempted to construct a ToroidalBlockLevelExchangeCenter object when not 3D");
#endif
}


void ToroidalBlockLevelExchangeCenter::interpGhosts( LevelData<FArrayBox>& a_data ) const
{
#if CFG_DIM==3
   m_geometry.interpolateFromShearedGhosts(a_data);
#endif
}


IntVectSet ToroidalBlockLevelExchangeCenter::extraBlockGhosts( const Box&  a_box,
                                                               int         a_ghosts ) const
{
   // This is a slight modification of MultiBlockUtil::extraBlockGhosts()
   IntVect ghostLayerVect = a_ghosts*IntVect::Unit;
   Box grownBox = grow(a_box, ghostLayerVect);
  
   const RefCountedPtr<MagCoordSys> coord_sys_ptr = m_geometry.getCoordSys();
   Vector< Tuple<BlockBoundary, 2*SpaceDim> > boundaries = coord_sys_ptr->boundaries();
   int block_number = coord_sys_ptr->whichBlock(a_box);
   const MagBlockCoordSys* block_coord_sys = coord_sys_ptr->getCoordSys(a_box);
   const Box& block_box = block_coord_sys->domain().domainBox();

   IntVectSet ivsReturn = IntVectSet(); // empty
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
            for (int idir = 0; idir < SpaceDim; idir++) {
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


