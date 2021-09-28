#include <math.h>
#include "MBSolver.H"
#include "MBSolverF_F.H"
#include "MBStencilIterator.H"
#include "BlockBaseRegister.H"
#include "BlockRegister.H"
#include "SparseCoupling.H"
#include "newMappedGridIO.H"
#include "MagGeom.H"
#include "DataArray.H"

#include "NamespaceHeader.H"



double
L2Norm( const LevelData<FArrayBox>&  a_data )
{
   const DisjointBoxLayout& grids = a_data.getBoxes();

   double local_sum = 0.;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      local_sum += a_data[dit].dotProduct(a_data[dit],grids[dit]);
   }

   double global_sum;
#ifdef CH_MPI
   MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
   global_sum = local_sum;
#endif

   return sqrt(global_sum);
}


MBSolver::MBSolver( const MultiBlockLevelGeom&      a_geom,
                    const int                       a_discretization_order,
                    MultiBlockLevelExchangeCenter*  a_mblex_ptr )

   : m_geometry(a_geom),
     m_coord_sys_ptr(m_geometry.coordSysPtr()),
     m_mblex_potential_Ptr(a_mblex_ptr),
     m_discretization_order(a_discretization_order),
     m_method_params_set(false),
     m_convergence_params_set(false)
{
   CH_assert(m_discretization_order == 2 || m_discretization_order == 4);

   m_rhs_from_bc.define(a_geom.grids(), 1, IntVect::Zero);
   
   // If there is more than one block and no MultiBlockLevelExchangeCenter object has
   // been provided, construct one
   if ( m_coord_sys_ptr->numBlocks() > 1 && m_mblex_potential_Ptr == NULL) {
      m_mblex_potential_Ptr = new MultiBlockLevelExchangeCenter();
      int num_ghosts_filled = m_discretization_order / 2;
      int spaceOrder = m_discretization_order;

      m_mblex_potential_Ptr->define(&a_geom, num_ghosts_filled, spaceOrder);
      m_mblex_defined_locally = true;
   }
   else {
      m_mblex_defined_locally = false;
   }

   for (int dir=0; dir<SpaceDim; ++dir) {
      m_flux_average[dir] = true;
   }

#if CFG_DIM==3
   if ( ((MagGeom&)m_geometry).shearedMBGeom() ) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         m_flux_average[dir] = false;
      }
   }
#endif
}
      


MBSolver::~MBSolver()
{
   if (m_mblex_defined_locally) delete m_mblex_potential_Ptr;
}



void
MBSolver::setConvergenceParams( const ParmParse&  a_pp )
{
   if ( a_pp.query("tol", m_method_tol) == 0 ) m_method_tol = 0.;
   if ( a_pp.query("max_iter", m_method_max_iter) == 0 ) m_method_max_iter = 1;
   if ( a_pp.query("verbose", m_method_verbose) == 0 ) m_method_verbose = false;

   ParmParse pp_precond( ((string)a_pp.prefix() + ".precond").c_str());
   if ( pp_precond.query("tol", m_precond_tol) == 0 ) m_precond_tol = 0.;
   if ( pp_precond.query("max_iter", m_precond_max_iter) == 0 ) m_precond_max_iter = 0;
   if ( pp_precond.query("verbose", m_precond_verbose) == 0 ) m_precond_verbose = false;

   m_convergence_params_set = true;
}


void
MBSolver::setConvergenceParams( const double  a_method_tol,
                                const int     a_method_max_iter,
                                const bool    a_method_verbose,
                                const double  a_precond_tol,
                                const int     a_precond_max_iter,
                                const bool    a_precond_verbose )
{
   m_method_tol       = a_method_tol;
   m_method_max_iter  = a_method_max_iter;
   m_method_verbose   = a_method_verbose;
   m_precond_tol      = a_precond_tol;
   m_precond_max_iter = a_precond_max_iter;
   m_precond_verbose  = a_precond_verbose;

   m_convergence_params_set = true;
}


void
MBSolver::constructMatrix( LevelData<FluxBox>&  a_tensor_coefficient,
                           const EllipticOpBC&  a_bc )
{
   LevelData<FArrayBox> dummy;
   constructMatrixGeneral(dummy, a_tensor_coefficient, dummy, a_bc );
}



void
MBSolver::constructMatrixBlock( const int            a_block_row,
                                const int            a_block_column,
                                LevelData<FluxBox>&  a_tensor_coefficient,
                                const EllipticOpBC&  a_bc )
{
   LevelData<FArrayBox> dummy;
   constructMatrixBlockGeneral(a_block_row, a_block_column, dummy, a_tensor_coefficient, dummy, a_bc);
}



void
MBSolver::constructMatrix( LevelData<FArrayBox>&  a_alpha_coefficient,
                           LevelData<FluxBox>&    a_tensor_coefficient,
                           const EllipticOpBC&    a_bc )
{
   LevelData<FArrayBox> dummy;
   constructMatrixGeneral(a_alpha_coefficient, a_tensor_coefficient, dummy, a_bc);
}



void
MBSolver::constructMatrixBlock( const int              a_block_row,
                                const int              a_block_column,
                                LevelData<FArrayBox>&  a_alpha_coefficient,
                                LevelData<FluxBox>&    a_tensor_coefficient,
                                const EllipticOpBC&    a_bc )
{
   LevelData<FArrayBox> dummy;
   constructMatrixBlockGeneral(a_block_row, a_block_column, a_alpha_coefficient, a_tensor_coefficient,
                               dummy, a_bc);
}



void
MBSolver::constructMatrix( LevelData<FluxBox>&    a_tensor_coefficient,
                           LevelData<FArrayBox>&  a_beta_coefficient,
                           const EllipticOpBC&    a_bc )
{
   LevelData<FArrayBox> dummy;
   constructMatrixGeneral(dummy, a_tensor_coefficient, a_beta_coefficient, a_bc );
}



void
MBSolver::constructMatrixBlock( const int              a_block_row,
                                const int              a_block_column,
                                LevelData<FluxBox>&    a_tensor_coefficient,
                                LevelData<FArrayBox>&  a_beta_coefficient,
                                const EllipticOpBC&    a_bc )
{
   LevelData<FArrayBox> dummy;
   constructMatrixBlockGeneral(a_block_row, a_block_column, dummy, a_tensor_coefficient, a_beta_coefficient, a_bc );
}



void
MBSolver::constructMatrix( LevelData<FArrayBox>& a_alpha_coefficient,
                           LevelData<FluxBox>&   a_tensor_coefficient,
                           LevelData<FArrayBox>& a_beta_coefficient,
                           const EllipticOpBC&   a_bc )
{
   constructMatrixGeneral( a_alpha_coefficient, a_tensor_coefficient, a_beta_coefficient, a_bc);
}



void
MBSolver::constructMatrixBlock( const int             a_block_row,
                                const int             a_block_column,
                                LevelData<FArrayBox>& a_alpha_coefficient,
                                LevelData<FluxBox>&   a_tensor_coefficient,
                                LevelData<FArrayBox>& a_beta_coefficient,
                                const EllipticOpBC&   a_bc )
{
   constructMatrixBlockGeneral(a_block_row, a_block_column, a_alpha_coefficient, a_tensor_coefficient,
                               a_beta_coefficient, a_bc);
}



void
MBSolver::averageAtBlockBoundaries(LevelData<FluxBox>& a_data) const
{
   if ( m_coord_sys_ptr->numBlocks() > 1 ) {

      const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

      RefCountedPtr<MultiBlockCoordSys> coordSysRCP(m_coord_sys_ptr);
      coordSysRCP.neverDelete();

      BlockRegister blockRegister(coordSysRCP, grids, 0);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int idir = 0; idir < SpaceDim; idir++) {
            if ( m_flux_average[idir] ) {
               for (SideIterator sit; sit.ok(); ++sit) {
                  Side::LoHiSide side = sit();
                  if (blockRegister.hasInterface(dit(), idir, side)) {
                     FArrayBox flux_comp(a_data[dit][idir].box(), a_data.nComp());
                     flux_comp.copy(a_data[dit][idir]);
                     blockRegister.storeFlux(flux_comp, dit(), idir, side);
                  }
               }
            }
         }
      }
      blockRegister.close();

      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = m_coord_sys_ptr->boundaries();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& baseBox = grids[dit];
         int block_number = m_coord_sys_ptr->whichBlock(grids[dit]);
         int faceID = 0;
         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
            for (int idir = 0; idir < SpaceDim; idir++) {
               if ( m_flux_average[idir] && blockRegister.hasInterface(dit(), idir, side) ) {
                  // maybe better if this is done inside BlockRegister
                  const BlockBoundary& bb = boundaries[block_number][faceID];
                  int reorientFace = bb.reorientFace(idir);
                  Box faceBox = adjCellBox(baseBox, idir, side, 1);
                  // if Lo, then shift +1; if Hi, then shift -1
                  faceBox.shiftHalf(idir, -sign(side));
                  Side::LoHiSide sideOther = flip(side);
                  // Need to define these FABs.
                  FArrayBox fluxThisFab(faceBox, a_data.nComp());
                  FArrayBox fluxOtherFab(faceBox, a_data.nComp());
                  blockRegister.getFlux(fluxThisFab, dit(),
                                        idir, side, side);
                  fluxThisFab.mult(reorientFace * 0.5);
                  blockRegister.getFlux(fluxOtherFab, dit(),
                                        idir, side, sideOther);
                  fluxOtherFab.mult(0.5);
                  fluxThisFab += fluxOtherFab;
                  a_data[dit][idir].copy(fluxThisFab);
               }
               faceID++;
            } // iterate over dimensions
         } // iterate over sides
      }
   }
}



IntVectSet
MBSolver::getInterBlockCoupledCells( const int   a_block_number,
                                     const int   a_radius,
                                     const Box&  a_box ) const
{
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys_ptr->boundaries();

  const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[a_block_number];
  const Box& domain_box = (m_coord_sys_ptr->mappingBlocks())[a_block_number];
  CH_assert(domain_box.contains(a_box));

  IntVectSet ivs;
  for (int dir=0; dir<SpaceDim; ++dir) {
    if ( this_block_boundaries[dir].isInterface() ) {
      int d = a_box.smallEnd(dir) - domain_box.smallEnd(dir);
      if ( d < a_radius ) ivs |= adjCellLo(a_box, dir, -(a_radius - d));
    }
    if ( this_block_boundaries[dir + SpaceDim].isInterface() ) {
      int d = domain_box.bigEnd(dir) - a_box.bigEnd(dir);
      if ( d < a_radius ) ivs |= adjCellHi(a_box, dir, -(a_radius - d));
    }
  }

  return ivs;
}



IntVectSet
MBSolver::getBoundaryCoupledCells( const int   a_block_number,
                                   const int   a_radius,
                                   const Box&  a_box ) const
{
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys_ptr->boundaries();

  const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[a_block_number];
  const Box& domain_box = (m_coord_sys_ptr->mappingBlocks())[a_block_number];
  CH_assert(domain_box.contains(a_box));

  IntVectSet ivs;
  for (int dir=0; dir<SpaceDim; ++dir) {
    if ( this_block_boundaries[dir].isDomainBoundary() ) {
      int d = a_box.smallEnd(dir) - domain_box.smallEnd(dir);
      if ( d < a_radius ) ivs |= adjCellLo(a_box, dir, -(a_radius - d));
    }
    if ( this_block_boundaries[dir + SpaceDim].isDomainBoundary() ) {
      int d = domain_box.bigEnd(dir) - a_box.bigEnd(dir);
      if ( d < a_radius ) ivs |= adjCellHi(a_box, dir, -(a_radius - d));
    }
  }

  return ivs;
}



void
MBSolver::getUnstructuredCouplings( int                                 a_radius,
                                    LayoutData< BaseFab<IntVectSet> >&  a_unstructured_couplings ) const
{

  if ( m_coord_sys_ptr->numBlocks() > 1 || ((MagGeom&)m_geometry).shearedMBGeom()) {

      LayoutData< Vector< BaseFab<IntVectSet>* > > all_couplings;

      const LayoutData< IntVectSet >& ghostCells = m_mblex_potential_Ptr->ghostCells();
      const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencil = m_mblex_potential_Ptr->stencils();

      const DisjointBoxLayout & grids = m_geometry.gridsFull();
    
      all_couplings.define(grids);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& box = grids[dit];
         for (int dir=0; dir<SpaceDim; ++dir) {
            Box dir_box = box;
            dir_box.surroundingNodes(dir);
            BaseFab<IntVectSet>* bi = new BaseFab<IntVectSet>(dir_box,1);
            all_couplings[dit].push_back(bi);
         }
      }
      
      const Vector<Box>& mapping_blocks = m_coord_sys_ptr->mappingBlocks();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const IntVectSet& this_ghosts_fab = ghostCells[dit];
         if ( stencil[dit] != NULL  ) {
            const IVSFAB<MBStencil>& this_stencil_fab = *stencil[dit];
            int block_number = m_coord_sys_ptr->whichBlock(grids[dit]);
            const Box& domain_box = mapping_blocks[block_number];
      
            Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
            for (int dir=0; dir<SpaceDim; ++dir) {
               BaseFab<IntVectSet>& this_coupling_dir = *this_coupling[dir];
               const Box& dir_box = this_coupling_dir.box();
        
               Box template_box(IntVect(D_DECL6(-a_radius,-a_radius,-a_radius,-a_radius,-a_radius,-a_radius)),
                                IntVect(D_DECL6(a_radius,a_radius,a_radius,a_radius,a_radius,a_radius)));
               template_box.growLo(dir, 1);
        
               BoxIterator bit(dir_box);
               for (bit.begin();bit.ok();++bit) {
                  IntVect iv_face = bit();
                  Box stencil_box = template_box + iv_face;

                  // Find the cells in the current stencil that are valid in the current block
                  IntVectSet valid_cells(domain_box & stencil_box);

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
      }

      // For faces on interfaces between blocks, we union the list of unstructured couplings with
      // the list on the corresponding face in the adjacent block.  This is needed to create
      // a system matrix corresponding to a conservative discretization in which fluxes are
      // averaged at multiblock interfaces.

      BlockBaseRegister<BaseFab<SparseCoupling> > blockRegister(m_coord_sys_ptr, grids, 0);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
         int dst_block_number = m_coord_sys_ptr->whichBlock(grids[dit]);

         // Find all cells that may be coupled through a block interface
         IntVectSet dst_ivs = getInterBlockCoupledCells(dst_block_number, a_radius, grids[dit]);

         for (int dir = 0; dir < SpaceDim; dir++) {
            if ( m_flux_average[dir] ) {
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
      }

      blockRegister.exchange();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& baseBox = grids[dit];
         Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
            for (int idir = 0; idir < SpaceDim; idir++) {
               if ( m_flux_average[idir] ) {
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
      }

      const Box stencil_offsets(IntVect(D_DECL6(-a_radius,-a_radius,-a_radius,-a_radius,-a_radius,-a_radius)),
                            IntVect(D_DECL6(a_radius,a_radius,a_radius,a_radius,a_radius,a_radius)));

      a_unstructured_couplings.define(grids);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
         int dst_block_number = m_coord_sys_ptr->whichBlock(grids[dit]);
         Box domain_box = mapping_blocks[dst_block_number];
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

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
         for (int i=0; i<this_coupling.size(); ++i) {
            delete this_coupling[i];
         }
      }
   }
}



IntVectSet
MBSolver::getUnstructuredCouplingsForCell( const Vector< BaseFab<IntVectSet>* >&  a_coupling,
                                           const Box&                             a_stencil_box,
                                           const IntVect&                         a_iv ) const
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
MBSolver::accumStencilMatrixEntries(const IntVect    a_index,
                                    const int        a_dir,
                                    const int        a_side,
                                    const int        a_dir2,
                                    const FluxBox&   a_mapped_coefs,
                                    const RealVect&  a_dx,
                                    const bool       a_fourthOrder,
                                    FArrayBox&       a_stencil_values) const
{
   if (a_fourthOrder) {
        
      FORT_ACCUM_FLUX_STENCIL4( CHF_CONST_INT(a_dir),
                                CHF_CONST_INT(a_dir2),
                                CHF_CONST_INT(a_side),
                                CHF_CONST_REALVECT(a_dx),
                                CHF_CONST_FRA1(a_mapped_coefs[a_dir],SpaceDim*a_dir+a_dir2),
                                CHF_CONST_INTVECT(a_index),
                                CHF_FRA1(a_stencil_values,0) );
   }
   else {

      FORT_ACCUM_FLUX_STENCIL2( CHF_CONST_INT(a_dir),
                                CHF_CONST_INT(a_dir2),
                                CHF_CONST_INT(a_side),
                                CHF_CONST_REALVECT(a_dx),
                                CHF_CONST_FRA1(a_mapped_coefs[a_dir],SpaceDim*a_dir+a_dir2),
                                CHF_CONST_INTVECT(a_index),
                                CHF_FRA1(a_stencil_values,0) );
   }
}



void
MBSolver::modifyStencilForBCs( const Vector<CoDim1Stencil>&  a_codim1_stencils,
                               const Vector<CoDim2Stencil>&  a_codim2_stencils,
                               const IntVect&                a_iv,
                               FArrayBox&                    a_stencil_values,
                               FArrayBox&                    a_rhs_from_bc,
                               const bool                    a_rhs_from_bc_only,
                               const bool                    a_force_codim2_condense ) const
{
   bool updating_rhs = a_rhs_from_bc.box().ok();  // FAB defined?

   CH_assert((a_stencil_values.box().size(0)-1)%2 == 0);
   int radius = (a_stencil_values.box().size(0) - 1) / 2;

   Box stencil_box(-radius*IntVect::Unit, radius*IntVect::Unit);
   stencil_box += a_iv;

   IntVect stencil_shift = a_iv - radius*IntVect::Unit;
   IntVectSet condensed_ivs;

   // Condense the codim2 stencils into the codim1 stencils

   for (int i=0; i<a_codim2_stencils.size(); ++i) {
      const CoDim2Stencil& codim2_stencil = a_codim2_stencils[i];

      Box overlap = stencil_box & codim2_stencil.box();

      if (overlap.ok()) {

         for (BoxIterator bit(overlap); bit.ok(); ++bit) {
            IntVect iv_ghost = bit();
            IntVect stencil_iv = iv_ghost - stencil_shift;

            if ( !a_rhs_from_bc_only || a_force_codim2_condense ) {

               double stencil_value = a_stencil_values(stencil_iv,0);

               vector<IntVect> points;
               vector<double> weights;

               codim2_stencil.getStencil(iv_ghost, points, weights);

               for (int j=0; j<points.size(); ++j) {
                  CH_assert(stencil_box.contains(points[j]));
                  a_stencil_values(points[j] - stencil_shift,0) += stencil_value * weights[j];
               }
            }

            // Remember this stencil_iv so that the corresponding stencil value can be zeroed
            // later.  We can't zero it yet, because stencil_iv might be contained in multiple
            // codim2_stencils.

            condensed_ivs |= stencil_iv;
         }
      }
   }

   // Now that all of the codim2 stencils have been condensed out, zero them

   for (IVSIterator it(condensed_ivs); it.ok(); ++it) {
      a_stencil_values(it(),0) = 0.;
   }

   // Condense the codim1 stencils into the valid stencils and accumulate
   // the contribution to the right hand-side due to inhomogeneous boundary
   // values

   for (int i=0; i<a_codim1_stencils.size(); ++i) {
      const CoDim1Stencil& codim1_stencil = a_codim1_stencils[i];

      Box overlap = stencil_box & codim1_stencil.box();

      if (overlap.ok()) {

         for (BoxIterator bit(overlap); bit.ok(); ++bit) {
            IntVect iv_ghost = bit();
            IntVect stencil_iv = iv_ghost - stencil_shift;

            double& stencil_value = a_stencil_values(stencil_iv,0);

            vector<IntVect> points;
            vector<double> weights;
            double bv_contrib;

            codim1_stencil.getStencil(iv_ghost, points, weights, bv_contrib);
                  
            if ( !a_rhs_from_bc_only ) {
               for (int j=0; j<points.size(); ++j) {
                  CH_assert(stencil_box.contains(points[j]));
                  a_stencil_values(points[j] - stencil_shift,0) += stencil_value * weights[j];
               }
            }

            if ( updating_rhs ) {
               a_rhs_from_bc(a_iv,0) -= stencil_value * bv_contrib;
            }

            // We can zero the stencil value here since no other codim1 bc_stencil object
            // will contain the current stencil_iv

            stencil_value = 0.;
         }
      }
   }

   if ( a_rhs_from_bc_only ) {
      a_stencil_values.setVal(0.);
   }
}



void
MBSolver::constructBoundaryStencils( const bool                        a_fourth_order,
                                     const EllipticOpBC&               a_bc,
                                     Vector< Vector<CoDim1Stencil> >&  a_codim1_stencils,
                                     Vector< Vector<CoDim2Stencil> >&  a_codim2_stencils ) const
{
   CH_assert(a_codim1_stencils.size() == 0);
   CH_assert(a_codim2_stencils.size() == 0);

   int order = a_fourth_order? 4: 2;

   int num_blocks = m_coord_sys_ptr->numBlocks();
   a_codim1_stencils.resize(num_blocks);
   a_codim2_stencils.resize(num_blocks);

   const Vector<Box>& mapping_blocks = m_coord_sys_ptr->mappingBlocks();

   for (int block_number=0; block_number<num_blocks; ++block_number) {
      const NewCoordSys* block_coord_sys = m_coord_sys_ptr->getCoordSys(block_number);
      Box domain_box = mapping_blocks[block_number];
      RealVect dx = block_coord_sys->dx();

      Vector<CoDim1Stencil>& codim1_stencil = a_codim1_stencils[block_number];      
      Vector<CoDim2Stencil>& codim2_stencil = a_codim2_stencils[block_number];      

      int num_codim1_neighbors = 0;
      int num_codim2_neighbors = 0;

      for (int dir=0; dir<SpaceDim; ++dir) {
         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();

            Vector<Box> boundary_boxes = ((MagGeom&)m_geometry).getBoundaryBoxes(block_number, dir, side);

            for ( int n=0; n<boundary_boxes.size(); ++n ) {

               codim1_stencil.resize(++num_codim1_neighbors);

               Box box = boundary_boxes[n];
               FArrayBox bv(box,1);

	       Box box_tmp = box;
               box_tmp.shiftHalf(dir,-sign(side));
               FluxBox bv_tmp(box_tmp,1);

               RefCountedPtr<GridFunction> bc_func = a_bc.getBCFunction(block_number, dir, side );
	       if (bc_func && !(typeid(*bc_func) == typeid(DataArray))) {
		  FluxBox real_coords(box_tmp,CFG_DIM);
		  FluxBox norm_flux(box_tmp,1);

                  const MagBlockCoordSys& coord_sys( ((MagGeom&)m_geometry).getBlockCoordSys( block_number ) );

                  coord_sys.getFaceCenteredRealCoords(dir, real_coords[dir]);
		  coord_sys.getNormMagneticFlux(real_coords[dir], norm_flux[dir]);

                  bc_func->assign(bv_tmp, m_geometry, real_coords, norm_flux, block_number, 0., false);
		  
                  for (BoxIterator bit(box); bit.ok(); ++bit) {
                     IntVect iv = bit();
                     iv[dir] = (bv_tmp[dir].box()).sideEnd(side)[dir];
                     bv(bit(),0) = bv_tmp[dir](iv,0);
                  }
               }

	       else if (bc_func && typeid(*bc_func) == typeid(DataArray)) {
		  FArrayBox dummy;
		  bc_func->assign(bv, m_geometry, dummy, dummy, block_number, 0., false);
	       }

	       else {
                  for (BoxIterator bit(box); bit.ok(); ++bit) {
                     bv(bit(),0) = a_bc.getBCValue(block_number, dir, side);
                  }
               }
               
               int bc_type = a_bc.getBCType(block_number, dir, side);

               codim1_stencil[num_codim1_neighbors-1].define(bv, bc_type, dx, dir, side, order);

               Box codim1_box = codim1_stencil[num_codim1_neighbors-1].box();

               CH_assert(order%2 == 0);
               int width = order/2;

               for (int tdir=0; tdir < SpaceDim; ++tdir) {
                  if ( tdir != dir ) {
                     for (SideIterator sit2; sit2.ok(); ++sit2) {
                        Side::LoHiSide side2 = sit2();

                        // Determine if the codim2 stencil we are about to create overlaps a
                        // physical boundary in the tdir direction on side2

                        Vector<Box> transverse_boundary_boxes 
                           = ((MagGeom&)m_geometry).getBoundaryBoxes(block_number, tdir, side2);

                        bool transverse_boundary = transverse_boundary_boxes.size() > 0;

                        if ( transverse_boundary ) {

                           Box grown_codim1_box = codim1_box;
                           grown_codim1_box.grow(tdir,width);

                           for (int n=0; n<transverse_boundary_boxes.size(); ++n) {

                              Box restricted_codim1_box = codim1_box;
                              bool this_transverse_boundary = false;

                              Box tbox = transverse_boundary_boxes[n];
                              tbox.shiftHalf(tdir,sign(side2));
                              tbox.grow(dir,width);
                              Box overlap = tbox & grown_codim1_box;
                              if ( overlap.ok() ) {
                                 for (int tdir2=0; tdir2<SpaceDim; ++tdir2) {
                                    if ( tdir2 != tdir ) {
                                       restricted_codim1_box.setSmall(tdir2, overlap.smallEnd(tdir2));
                                       restricted_codim1_box.setBig(tdir2, overlap.bigEnd(tdir2));
                                    }
                                 }
                                 this_transverse_boundary = true;
                              }

                              codim2_stencil.resize(++num_codim2_neighbors);

                              codim2_stencil[num_codim2_neighbors-1].
                                 define(restricted_codim1_box, bc_type, dx, dir, side, tdir, side2, order, this_transverse_boundary);

                              // Check to see if the new codim2 box overlaps another block

                              const Box& new_codim2_box = codim2_stencil[num_codim2_neighbors-1].box();

                              for (int block_number2=0; block_number2<num_blocks; ++block_number2) {
                                 if (block_number2 != block_number) {
                                    const Box& domain_box2 = mapping_blocks[block_number2];

                                    Box overlap = new_codim2_box & domain_box2;
                                    if ( overlap == new_codim2_box ) {
                                       // Codim2 box is contained in a valid block, so delete it

                                       codim2_stencil.resize(--num_codim2_neighbors);
                                    }
                                    else if ( overlap.ok() ) {
                                       // Codim2 box partially overlaps a valid block.  Need to exit and figure out
                                       // what to do about it.
                                       MayDay::Error("MBSolver::constructBoundaryStencils(): codim2 box partially overlaps a valid block");
                                    }
                                 }
                              }
                           }
                        }
                        else {

                           codim2_stencil.resize(++num_codim2_neighbors);

                           codim2_stencil[num_codim2_neighbors-1].
                              define(codim1_box, bc_type, dx, dir, side, tdir, side2, order, transverse_boundary);

                           // Check to see if the new codim2 box overlaps another block

                           const Box& new_codim2_box = codim2_stencil[num_codim2_neighbors-1].box();

                           for (int block_number2=0; block_number2<num_blocks; ++block_number2) {
                              if (block_number2 != block_number) {
                                 const Box& domain_box2 = mapping_blocks[block_number2];

                                 Box overlap = new_codim2_box & domain_box2;
                                 if ( overlap == new_codim2_box ) {
                                    // Codim2 box is contained in a valid block, so delete it

                                    codim2_stencil.resize(--num_codim2_neighbors);
                                 }
                                 else if ( overlap.ok() ) {
                                    // Codim2 box partially overlaps a valid block.  Need to exit
                                    // and figure out what to do about it.
                                    MayDay::Error("MBSolver::constructBoundaryStencils(): codim2 box partially overlaps a valid block");
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}



void
MBSolver::accumPhysicalGhosts( const Vector< Vector<CoDim1Stencil> >&  a_codim1_stencils,
                               const Vector< Vector<CoDim2Stencil> >&  a_codim2_stencils,
                               const bool                              a_extrapolate_from_interior,
                               const bool                              a_include_bvs,
                               LevelData<FArrayBox>&                   a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = m_coord_sys_ptr->whichBlock(grids[dit]);
      FArrayBox& this_data = a_data[dit];
      const Box& this_data_box = this_data.box();

      // Fill the codim1 ghosts first, which are assumed to depend only upon valid data
      // and boundary values

      const Vector<CoDim1Stencil>& codim1_stencil = a_codim1_stencils[block_number];      

      for (int i=0; i<codim1_stencil.size(); ++i) {
         const CoDim1Stencil& bndry_stencil = codim1_stencil[i];
         const Box overlap = bndry_stencil.box() & this_data_box;

         for (BoxIterator bit(overlap); bit.ok(); ++bit) {
            IntVect iv = bit();

            vector<IntVect> points;
            vector<double> weights;
            double bv_contrib;

            bndry_stencil.getStencil(iv, points, weights, bv_contrib);
                  
            if (a_include_bvs) {
               this_data(iv,0) += bv_contrib;
            }
            if (a_extrapolate_from_interior) {
               for (int j=0; j<points.size(); ++j) {
                  this_data(iv,0) += weights[j] * this_data(points[j],0);
               }
            }
         }
      }

      // Fill the codim2 ghosts next, which are assumed to depend upon valid data
      // and/or codim1 data

      const Vector<CoDim2Stencil>& codim2_stencil = a_codim2_stencils[block_number];      

      for (int i=0; i<codim2_stencil.size(); ++i) {
         const CoDim2Stencil& bndry_stencil = codim2_stencil[i];
         const Box overlap = bndry_stencil.box() & this_data_box;

         for (BoxIterator bit(overlap); bit.ok(); ++bit) {
            IntVect iv = bit();

            vector<IntVect> points;
            vector<double> weights;

            bndry_stencil.getStencil(iv, points, weights);

            for (int j=0; j<points.size(); ++j) {
               this_data(iv,0) += weights[j] * this_data(points[j],0);
            }
         }
      }
   }
}



void
MBSolver::computeFluxNormalFromStencil( const LevelData<FArrayBox>&                        a_in,
                                        LevelData<FluxBox>&                                a_coefs,
                                        const Vector< Tuple<BlockBoundary, 2*SpaceDim> >&  a_block_boundaries,
                                        const Vector< Vector<CoDim1Stencil> >&             a_codim1_stencils,
                                        const Vector< Vector<CoDim2Stencil> >&             a_codim2_stencils,
                                        FArrayBox&                                         a_stencil_values,
                                        const bool                                         a_fourthOrder,
                                        const bool                                         a_extrapolate_from_interior,
                                        const bool                                         a_include_bvs,
                                        LevelData<FluxBox>&                                a_flux_normal ) const
{
   CH_assert(a_flux_normal.ghostVect() == IntVect::Zero);

   CH_assert((a_stencil_values.box().size(0)-1)%2 == 0);
   int radius = (a_stencil_values.box().size(0)-1)/2;

   const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

   LevelData<FArrayBox> tmp(grids, 1, (m_discretization_order/2 + 1)*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_flux_normal[dit].setVal(0.);
   }

   const Vector<Box>& mapping_blocks = m_coord_sys_ptr->mappingBlocks();

   for (int dir_outer=0; dir_outer<SpaceDim; ++dir_outer) {
      for (int dir_inner=0; dir_inner<SpaceDim; ++dir_inner) {

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            tmp[dit].setVal(0.);
            tmp[dit].copy(a_in[dit],grids[dit]);
         }

         // Fill internal ghosts
         if ( m_coord_sys_ptr->numBlocks() > 1 ) {
            ((MagGeom&)m_geometry).fillInternalGhosts(tmp);
         }
         tmp.exchange();

         accumPhysicalGhosts( a_codim1_stencils, a_codim2_stencils, a_extrapolate_from_interior, a_include_bvs, tmp );

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            int block_number = m_coord_sys_ptr->whichBlock(grids[dit]);
            const NewCoordSys* block_coord_sys = m_coord_sys_ptr->getCoordSys(block_number);
            const Box& domain_box = mapping_blocks[block_number];
            RealVect dx = block_coord_sys->dx();
            const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = a_block_boundaries[block_number];

            FArrayBox& this_flux = a_flux_normal[dit][dir_outer];
            const FluxBox& this_coef = a_coefs[dit];

            bool lo_block_interface = this_block_boundaries[dir_outer].isInterface();
            bool hi_block_interface = this_block_boundaries[dir_outer + SpaceDim].isInterface();

            for (BoxIterator bit_face(this_flux.box()); bit_face.ok(); ++bit_face) {
               IntVect iv_face = bit_face();

               bool on_lo_block_interface = (lo_block_interface && iv_face[dir_outer] == domain_box.smallEnd(dir_outer));
               bool on_hi_block_interface = (hi_block_interface && iv_face[dir_outer] == domain_box.bigEnd(dir_outer)+1);

               for (SideIterator sit; sit.ok(); ++sit) {
                  Side::LoHiSide side = sit();

                  IntVect iv = iv_face;
                  iv.shift(dir_outer, -side);

                  a_stencil_values.setVal(0.);

                  accumStencilMatrixEntries(iv, dir_outer, side, dir_inner, this_coef,
                                            dx, a_fourthOrder, a_stencil_values);

                  double flux_contrib = 0.;
                  for (BoxIterator bit_cell(a_stencil_values.box()); bit_cell.ok(); ++bit_cell) {
                     IntVect iv2 = iv + bit_cell() - radius * IntVect::Unit;
                     flux_contrib += a_stencil_values(bit_cell(),0) * tmp[dit](iv2,0);
                  }

                  if ( on_lo_block_interface ) {
                     if ( side == Side::LoHiSide::Lo ) {
                        this_flux(iv_face,0) += (1 - 2*side) * flux_contrib;
                     }
                     else {
                        continue;
                     }
                  }
                  else if ( on_hi_block_interface ) {
                     if ( side == Side::LoHiSide::Hi ) {
                        this_flux(iv_face,0) += (1 - 2*side) * flux_contrib;
                     }
                     else {
                        continue;
                     }
                  }
                  else {
                     this_flux(iv_face,0) += 0.5 * (1 - 2*side) * flux_contrib;
                  }
               }
            }
         }
      }
   }
}



void
MBSolver::testMatrixConstruct( LevelData<FArrayBox>&                              a_alpha_coefficient, 
                               LevelData<FluxBox>&                                a_tensor_coefficient,
                               LevelData<FArrayBox>&                              a_beta_coefficient,
                               const Vector< Tuple<BlockBoundary, 2*SpaceDim> >&  a_block_boundaries,
                               const Vector< Vector<CoDim1Stencil> >&             a_codim1_stencils,
                               const Vector< Vector<CoDim2Stencil> >&             a_codim2_stencils,
                               FArrayBox&                                         a_stencil_values,
                               const bool                                         a_fourthOrder,
                               const LevelData<FArrayBox>&                        a_rhs_from_bc ) const
{
   const DisjointBoxLayout& grids = a_rhs_from_bc.disjointBoxLayout();

   LevelData<FArrayBox> solution(grids, 1, IntVect::Zero);

   srand(time(NULL));
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (BoxIterator bit(grids[dit]); bit.ok(); ++bit) {
         IntVect iv = bit();
         solution[dit](iv,0) = (double)rand() / (double)RAND_MAX;
      }
   }

   bool extrapolate_from_interior, include_bvs;

   // We compute the normal flux in two parts, the first omitting boundary value contributions and the
   // second containing only boundary contributions.  This allows the flux with boundary value contributions
   // omitted to be averaged at block interfaces.

   LevelData<FluxBox> flux_normal_without_bvs(grids, 1, IntVect::Zero);
   extrapolate_from_interior = true;
   include_bvs = false;
   computeFluxNormalFromStencil( solution, a_tensor_coefficient, a_block_boundaries, a_codim1_stencils,
                                 a_codim2_stencils, a_stencil_values, a_fourthOrder,
                                 extrapolate_from_interior, include_bvs, flux_normal_without_bvs );
   averageAtBlockBoundaries(flux_normal_without_bvs);

   LevelData<FArrayBox> zero(grids, 1, IntVect::Zero);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      zero[dit].setVal(0.);
   }

   LevelData<FluxBox> flux_normal_only_bvs(grids, 1, IntVect::Zero);
   extrapolate_from_interior = false;
   include_bvs = true;
   computeFluxNormalFromStencil( zero, a_tensor_coefficient, a_block_boundaries, a_codim1_stencils,
                                 a_codim2_stencils, a_stencil_values, a_fourthOrder,
                                 extrapolate_from_interior, include_bvs, flux_normal_only_bvs );

   LevelData<FArrayBox> operator_eval(grids, 1, IntVect::Zero);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      FArrayBox & this_operator_eval = operator_eval[dit];
      FluxBox & this_flux_normal = flux_normal_without_bvs[dit];

      this_flux_normal += flux_normal_only_bvs[dit];

      const FArrayBox* alpha = NULL;
      if ( a_alpha_coefficient.isDefined() ) {
         alpha = &a_alpha_coefficient[dit];
      }

      const FArrayBox* beta = NULL;
      if ( a_beta_coefficient.isDefined() ) {
         beta = &a_beta_coefficient[dit];
      }

      for (BoxIterator bit(grids[dit]); bit.ok(); ++bit) {
         IntVect iv = bit();

         this_operator_eval(iv,0) = 0.;
         for (int dir=0; dir<SpaceDim; ++dir) {
            for (SideIterator sit; sit.ok(); ++sit) {
               Side::LoHiSide side = sit();

               IntVect iv_face = iv;
               iv_face.shift(dir,side);

               this_operator_eval(iv,0) += (1 - 2*side)*this_flux_normal[dir](iv_face,0);
            }
         }

         if ( alpha ) {
            this_operator_eval(iv,0) *= alpha->operator()(iv);
         }

         if ( beta ) {
            this_operator_eval(iv,0) += beta->operator()(iv) * solution[dit](iv,0);
         }
      }
   }

   plot("operator_eval", operator_eval, 0.);

   LevelData<FArrayBox> matvec(grids, 1, IntVect::Zero);
   multiplyMatrix(solution, matvec);
   plot("matvec", matvec, 0.);

   LevelData<FArrayBox> diff(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      diff[dit].copy(matvec[dit]);
      diff[dit] -= a_rhs_from_bc[dit];
      diff[dit] -= operator_eval[dit];
   }
   plot("diff", diff, 0.);

   double diff_norm = L2Norm(diff);

   if (procID()==0) {
      cout << "Matrix diff = " << diff_norm << " (using random input)" << endl;
   }

   IntVect max_iv;
   double local_max_diff = -DBL_MAX;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (BoxIterator bit(grids[dit]); bit.ok(); ++bit) {
         IntVect iv = bit();

         if ( fabs(diff[dit](iv,0)) > local_max_diff ) {
            local_max_diff = fabs(diff[dit](iv,0));
            max_iv = iv;
         }
      }
   }

   double global_max_diff;
#ifdef CH_MPI
   MPI_Allreduce(&local_max_diff, &global_max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
   global_max_diff = local_max_diff;
#endif

#if 1
   if ( global_max_diff == local_max_diff ) {
      cout << "Max diff at " << max_iv << " = " << global_max_diff << endl;
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (BoxIterator bit(grids[dit]); bit.ok(); ++bit) {
         IntVect iv = bit();
         if ( fabs(diff[dit](iv,0)) > 1.e-14)  cout << iv << " " << fabs(diff[dit](iv,0)) << endl;
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);
   exit(1);
#endif
}



void MBSolver::plot( const string&               a_file_name,  
                     const LevelData<FArrayBox>& a_data,
                     const double&               a_time ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   Box domain_box = grids.physDomain().domainBox();
   domain_box.grow(a_data.ghostVect());
   WriteMappedUGHDF5(a_file_name.c_str(), grids, a_data, *m_coord_sys_ptr, domain_box, a_time);
}



void
MBSolver::setZero( LevelData<FArrayBox>&  a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit].setVal(0.);
   }
}



#include "NamespaceFooter.H"
