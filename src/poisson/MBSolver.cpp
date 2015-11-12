#include "MBSolver.H"
#include "MBSolverF_F.H"
#include "MBStencilIterator.H"
#include "BlockBaseRegister.H"


#include "NamespaceHeader.H"



MBSolver::MBSolver( const MagGeom& a_geom,
                    const int      a_discretization_order )
   : m_geometry(a_geom),
     m_num_potential_ghosts_filled(2),
     m_discretization_order(a_discretization_order),
     m_dropOrder(true)
{
   const DisjointBoxLayout& grids = m_geometry.grids();

   m_volume.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume);

   m_A_diagonal.define(grids, 1, IntVect::Zero);
   m_A_radial.define(grids, 3, IntVect::Zero);
   m_diagonal_increment.define(grids, 1, IntVect::Zero);

   // If there is more than one block, construct the multiblock exchange object
   if ( m_geometry.coordSysPtr()->numBlocks() > 1 ) {
     m_mblex_potential_Ptr = new MultiBlockLevelExchangeCenter();
     int spaceOrder = 4;
     m_mblex_potential_Ptr->define(&a_geom, m_num_potential_ghosts_filled, spaceOrder);
   }
   else {
     m_mblex_potential_Ptr = NULL;
   }
}
      


IntVectSet
MBSolver::getInterBlockCoupledCells( const int  a_block_number,
                                     const int  a_radius,
                                     const Box& a_box ) const
{
  MultiBlockCoordSys* coord_sys_ptr = m_geometry.coordSysPtr();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = coord_sys_ptr->boundaries();

  const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[a_block_number];
  const ProblemDomain& domain = ((const MagBlockCoordSys*)coord_sys_ptr->getCoordSys(a_block_number))->domain();
  const Box& domainBox = domain.domainBox();
  CH_assert(domainBox.contains(a_box));

  IntVectSet ivs;
  for (int dir=0; dir<SpaceDim; ++dir) {
    if ( this_block_boundaries[dir].isInterface() ) {
      int d = a_box.smallEnd(dir) - domainBox.smallEnd(dir);
      if ( d < a_radius ) ivs |= adjCellLo(a_box, dir, -(a_radius - d));
    }
    if ( this_block_boundaries[dir + SpaceDim].isInterface() ) {
      int d = domainBox.bigEnd(dir) - a_box.bigEnd(dir);
      if ( d < a_radius ) ivs |= adjCellHi(a_box, dir, -(a_radius - d));
    }
  }

  return ivs;
}



void
MBSolver::getUnstructuredCouplings(int                                a_radius,
                                   LayoutData< BaseFab<IntVectSet> >& a_unstructured_couplings ) const
{
   if (m_mblex_potential_Ptr) {
  
      LayoutData< Vector< BaseFab<IntVectSet>* > > all_couplings;

      const LayoutData< IntVectSet >& ghostCells = m_mblex_potential_Ptr->ghostCells();
      const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencil = m_mblex_potential_Ptr->stencils();

      MultiBlockCoordSys* coord_sys_ptr = m_geometry.coordSysPtr();
      const DisjointBoxLayout & grids = m_geometry.grids();
    
      all_couplings.define(grids);
      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         const Box& box = grids[dit];
         for (int dir=0; dir<SpaceDim; ++dir) {
            Box dir_box = box;
            dir_box.surroundingNodes(dir);
            BaseFab<IntVectSet>* bi = new BaseFab<IntVectSet>(dir_box,1);
            all_couplings[dit].push_back(bi);
         }
      }
      
      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         const IntVectSet& this_ghosts_fab = ghostCells[dit];
         const IVSFAB<MBStencil>& this_stencil_fab = *stencil[dit];
         int block_number = coord_sys_ptr->whichBlock(grids[dit]);
         const Box& domainBox = (coord_sys_ptr->blockDomainOfBox(grids[dit])).domainBox();
      
         Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
         for (int dir=0; dir<SpaceDim; ++dir) {
            BaseFab<IntVectSet>& this_coupling_dir = *this_coupling[dir];
            const Box& dir_box = this_coupling_dir.box();
            Box valid_domain = domainBox;
            valid_domain.enclosedCells(dir);
        
            Box template_box(IntVect(D_DECL6(-a_radius,-a_radius,-a_radius,-a_radius,-a_radius,-a_radius)),
                             IntVect(D_DECL6(a_radius,a_radius,a_radius,a_radius,a_radius,a_radius)));
            template_box.growHi(dir, -1);
        
            BoxIterator bit(dir_box);
            for (bit.begin();bit.ok();++bit) {
               IntVect iv_face = bit();
               Box stencil_box = template_box + iv_face;

               // Find the cells in the current stencil that are valid in the current block
               //               Box valid_box = dir_box;
               //               valid_box.enclosedCells(dir);
               //               IntVectSet valid_cells(valid_box & stencil_box);
               IntVectSet valid_cells(valid_domain & stencil_box);

               // Find the cells in the current stencil corresponding to extra-block
               // ghost cells filled by interpolation from other blocks
               IntVectSet extra_block_ghosts(this_ghosts_fab & stencil_box);

               // Find the union of the valid cells in the current block and the valid
               // cells used for interpolation into the extra block ghost cells
               IntVectSet union_ivs(valid_cells);
               for (IVSIterator it(extra_block_ghosts); it.ok(); ++it) {
                  MBStencilIterator stit(this_stencil_fab(it(), 0));
                  for (stit.begin(); stit.ok(); ++stit) {
                     const MBStencilElement& stencilElement = stit();
                     const IntVect& stencilCell = stencilElement.cell();

                     union_ivs |= stencilCell;
                  }
               }
               
               this_coupling_dir(iv_face) = union_ivs;
            }
         }
      }

      // For faces on interfaces between blocks, we union the list of unstructured couplings with
      // the list on the corresponding face in the adjacent block.  This is needed to create
      // a system matrix corresponding to a conservative discretization in which fluxes are
      // averaged at multiblock interfaces.

      BlockBaseRegister<BaseFab<SparseCoupling> > blockRegister(coord_sys_ptr, grids, 0);

      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         const Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
         int dst_block_number = coord_sys_ptr->whichBlock(grids[dit]);

         // Find all cells that may be coupled through a block interface
         IntVectSet dst_ivs = getInterBlockCoupledCells(dst_block_number, a_radius, grids[dit]);

         for (int dir = 0; dir < SpaceDim; dir++) {
            BaseFab<IntVectSet>& this_coupling_dir = *this_coupling[dir];
            BaseFab<SparseCoupling> data(surroundingNodes(grids[dit], dir), 1);

            for (SideIterator sit; sit.ok(); ++sit) {
               Side::LoHiSide side = sit();
               if (blockRegister.hasInterface(dit(), dir, side)) {

                  for (IVSIterator it(dst_ivs); it.ok(); ++it) {
                     IntVect iv = it();
                     iv.shift(dir,side);

                     IntVectSet& unstructured_ivs = this_coupling_dir(iv);
                     for (IVSIterator sivsit(unstructured_ivs); sivsit.ok(); ++sivsit) {
                        data(iv).add( sivsit(), 0. );
                     }
                  }

                  blockRegister.store(data, dit(), dir, side);
               }
            }
         }
      }

      blockRegister.exchange();

      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         const Box& baseBox = grids[dit];
         Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
            for (int idir = 0; idir < SpaceDim; idir++) {
               BaseFab<IntVectSet>& this_coupling_dir = *this_coupling[idir];
               if (blockRegister.hasInterface(dit(), idir, side)) {
                  Box faceBox = adjCellBox(baseBox, idir, side, 1);
                  // if Lo, then shift +1; if Hi, then shift -1
                  faceBox.shiftHalf(idir, -sign(side));
                  BaseFab<SparseCoupling> otherFab(faceBox, 1);
                  blockRegister.fill(otherFab, dit(), idir, side, side);
                  BoxIterator bit(faceBox);
                  for (bit.begin(); bit.ok(); ++bit) {
                     IntVect iv = bit();
                     SparseCoupling& coupling = otherFab(iv,0);
                     for (SparseCouplingIterator it(coupling); it.ok(); ++it) {
                        this_coupling_dir(iv) |= coupling[it()];
                     }
                  }
               }
            }
         }
      }

      const Box stencil_offsets(IntVect(D_DECL6(-a_radius,-a_radius,-a_radius,-a_radius,-a_radius,-a_radius)),
                            IntVect(D_DECL6(a_radius,a_radius,a_radius,a_radius,a_radius,a_radius)));

      a_unstructured_couplings.define(grids);
      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         const Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
         int dst_block_number = coord_sys_ptr->whichBlock(grids[dit]);
         Box domain_box = m_geometry.getBlockCoordSys(dst_block_number).domain().domainBox();
         BaseFab<IntVectSet>& this_unstructured_coupling = a_unstructured_couplings[dit];
         
         this_unstructured_coupling.resize(grids[dit], 1);

         // Find all cells that may be coupled through a block interface
         IntVectSet dst_ivs = getInterBlockCoupledCells(dst_block_number, a_radius, grids[dit]);

         for (IVSIterator it(dst_ivs); it.ok(); ++it) {
            IntVect iv = it();
            Box stencil_box = stencil_offsets + iv;
            Box stencil_box_valid = stencil_box & domain_box;

            this_unstructured_coupling(iv) = getUnstructuredCouplingsForCell( this_coupling, stencil_box_valid, iv );
         }
      }
   }
}



IntVectSet
MBSolver::getUnstructuredCouplingsForCell( const Vector< BaseFab<IntVectSet>* >& a_coupling,
                                           const Box&                            a_stencil_box,
                                           const IntVect&                        a_iv ) const
{
   IntVectSet unstructured_ivs;

   for (int dir=0; dir<SpaceDim; ++dir) {
      const BaseFab<IntVectSet>& this_coupling = *a_coupling[dir];
      for (int side=0; side<2; ++side) {
         IntVect index(a_iv);
         index.shift(dir,side);

         const IntVectSet& ivs = this_coupling(index);
                       
         for (IVSIterator ivit(ivs); ivit.ok(); ++ivit) {
            IntVect iv2 = ivit();
            if ( !a_stencil_box.contains(iv2) ) {
               unstructured_ivs |= iv2;
            }
         }
      }
   }

   return unstructured_ivs;
}



void
MBSolver::accumStencilMatrixEntries(const IntVect                           a_index,
                                    const int                               a_dir,
                                    const int                               a_side,
                                    const int                               a_s,
                                    const Tuple<BlockBoundary, 2*SpaceDim>& a_block_boundaries,
                                    const int*                              a_block_bc_types,
                                    const FluxBox&                          a_mapped_coefs,
                                    const RealVect&                         a_dx,
                                    const Box&                              a_domainBox,
                                    const bool                              a_fourthOrder,
                                    const bool                              a_dropOrder,
                                    FArrayBox&                              a_stencil_values) const
{
   int dist2physbdry[] = {2, 2, 2, 2};

   for (int dir=0; dir<SpaceDim; dir++) {
      if ( a_block_boundaries[dir].isDomainBoundary() ) {
         if (a_index[dir]==a_domainBox.smallEnd(dir)) {
            dist2physbdry[dir] = 0;
            dist2physbdry[dir+2] = -1;
         }
         if (a_index[dir]==a_domainBox.smallEnd(dir)+1) dist2physbdry[dir] = -1;
      }
      if ( a_block_boundaries[dir + SpaceDim].isDomainBoundary() ) {
         if (a_index[dir]==a_domainBox.bigEnd(dir)) {
            dist2physbdry[dir+2] = 0;
            dist2physbdry[dir] = 1;
         }
         if (a_index[dir]==a_domainBox.bigEnd(dir)-1) dist2physbdry[dir+2] = 1;
      }
   }

   for (int sprime=0; sprime<SpaceDim; sprime++) {

      // Add stencil contribution from flux on side in this direction 

      if (a_fourthOrder) {
        
         int dropOrder = a_dropOrder? 1: 0;

         FORT_ACCUM_FLUX_STENCIL4(CHF_CONST_INT(a_dir),
                                  CHF_CONST_INT(sprime),
                                  CHF_CONST_INT(a_side),
                                  CHF_CONST_I1D(dist2physbdry,4),
                                  CHF_CONST_REALVECT(a_dx),
                                  CHF_CONST_FRA1(a_mapped_coefs[a_dir],2*a_s+sprime),
                                  CHF_CONST_INTVECT(a_index),
                                  CHF_CONST_INT(dropOrder),
                                  CHF_FRA1(a_stencil_values,0));
      }
      else {

         FORT_ACCUM_FLUX_STENCIL2(CHF_CONST_INT(a_dir),
                                  CHF_CONST_INT(sprime),
                                  CHF_CONST_INT(a_side),
                                  CHF_CONST_I1D(dist2physbdry,4),
                                  CHF_CONST_I1D(a_block_bc_types,4),
                                  CHF_CONST_REALVECT(a_dx),
                                  CHF_CONST_FRA1(a_mapped_coefs[a_dir],2*a_s+sprime),
                                  CHF_CONST_INTVECT(a_index),
                                  CHF_FRA1(a_stencil_values,0));
      }
   }
}



void
MBSolver::setZero( LevelData<FArrayBox>& a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit].setVal(0.);
   }
}



#include "NamespaceFooter.H"
