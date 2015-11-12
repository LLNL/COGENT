#include "MBHypreSolver.H"
#include "MBStencilIterator.H"
#include "BlockRegister.H"
#include "SparseCoupling.H"
#include "FluxSurface.H"
#include "inspect.H"

#include "NamespaceHeader.H"



MBHypreSolver::MBHypreSolver( const MagGeom& a_geom,
                              const int      a_discretization_order )
   : MBSolver(a_geom, a_discretization_order),
     m_hypre_allocated(false),
     m_A(NULL)
{
   createHypreData();
}
      


MBHypreSolver::~MBHypreSolver()
{
   destroyHypreData();
}



void
MBHypreSolver::multiplyMatrix( const LevelData<FArrayBox>& a_in,
                               LevelData<FArrayBox>&       a_out )
{
   if (!m_hypre_allocated) {
      MayDay::Error( "GKPoisson::applyOperator(): Operator has not yet been initialized!" );
   }

   /* Create temporary Hypre vector objects */

   HYPRE_SStructVector in_vector;
   HYPRE_SStructVector out_vector;

   HYPRE_SStructVectorCreate(MPI_COMM_WORLD, m_grid, &in_vector);
   HYPRE_SStructVectorCreate(MPI_COMM_WORLD, m_grid, &out_vector);

   HYPRE_SStructVectorSetObjectType(in_vector, m_hypre_object_type);
   HYPRE_SStructVectorSetObjectType(out_vector, m_hypre_object_type);

   HYPRE_SStructVectorInitialize(in_vector);
   HYPRE_SStructVectorInitialize(out_vector);

   // Copy the input LevelData to its Hypre vector

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   int var = 0;

   const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = mag_coord_sys->whichBlock(grids[dit]);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      FArrayBox tmp(box,1);
      tmp.copy(a_in[dit]);

      HYPRE_SStructVectorSetBoxValues(in_vector, block_number, lower.dataPtr(), upper.dataPtr(), var, tmp.dataPtr());
   }

   HYPRE_SStructVectorAssemble(in_vector);
   HYPRE_SStructVectorAssemble(out_vector);

   // Do the matrix-vector multiply

   void *matvec_vdata;
   hypre_SStructMatvecCreate( &matvec_vdata );

   hypre_SStructMatvecSetup(matvec_vdata, m_A, in_vector);

   hypre_SStructMatvecCompute(matvec_vdata, 1., m_A, in_vector, 1., out_vector);

   hypre_SStructMatvecDestroy (matvec_vdata);

   if (m_hypre_object_type == HYPRE_PARCSR) {
     HYPRE_SStructVectorGather(out_vector);
   }

   // Copy the output Hypre vector to its LevelData

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = mag_coord_sys->whichBlock(grids[dit]);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      FArrayBox tmp(box,1);

      HYPRE_SStructVectorGetBoxValues(out_vector, block_number, lower.dataPtr(), upper.dataPtr(), var, tmp.dataPtr());

      a_out[dit].copy(tmp);
   }

   HYPRE_SStructVectorDestroy(out_vector);
   HYPRE_SStructVectorDestroy(in_vector);
}



void
MBHypreSolver::updateMatrixDiagonals( const LevelData<FArrayBox>& a_diagonal_increment )
{
   /*
     Add the passed-in data to the diagonals of the system matrix A
   */

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   int var = 0;

   const DisjointBoxLayout & grids = a_diagonal_increment.disjointBoxLayout();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box& box = grids[dit];
      int block_number = mag_coord_sys->whichBlock(grids[dit]);
      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      FArrayBox new_diagonal(box,1);

      new_diagonal.copy(m_A_diagonal[dit]);
      new_diagonal.plus(a_diagonal_increment[dit]);

      HYPRE_SStructMatrixSetBoxValues(m_A, block_number, lower.dataPtr(), upper.dataPtr(), var, 1,
                                      &m_A_diagonal_offset, new_diagonal.dataPtr());

      m_diagonal_increment[dit].copy(a_diagonal_increment[dit]);
   }

   /* This is a collective call finalizing the matrix assembly.
      The matrices are now ``ready to be used'' */
   HYPRE_SStructMatrixAssemble(m_A);
}



void
MBHypreSolver::constructMatrix( LevelData<FluxBox>&   a_diffusion_coefficient,
                                const PotentialBC&    a_bc,
                                const bool            a_poloidally_average )
{
   LevelData<FArrayBox> dummy;
   constructMatrix(dummy, a_diffusion_coefficient, dummy, a_bc, a_poloidally_average);
}



void
MBHypreSolver::constructMatrix( LevelData<FArrayBox>& a_alpha_coefficient,
                                LevelData<FluxBox>&   a_diffusion_coefficient,
                                const PotentialBC&    a_bc,
                                const bool            a_poloidally_average )
{
   LevelData<FArrayBox> dummy;
   constructMatrix(a_alpha_coefficient, a_diffusion_coefficient, dummy, a_bc, a_poloidally_average);
}



void
MBHypreSolver::constructMatrix( LevelData<FluxBox>&   a_diffusion_coefficient,
                                LevelData<FArrayBox>& a_beta_coefficient,
                                const PotentialBC&    a_bc,
                                const bool            a_poloidally_average )
{
   LevelData<FArrayBox> dummy;
   constructMatrix(dummy, a_diffusion_coefficient, a_beta_coefficient, a_bc, a_poloidally_average);
}



void
MBHypreSolver::constructMatrix( LevelData<FArrayBox>& a_alpha_coefficient,
                                LevelData<FluxBox>&   a_diffusion_coefficient,
                                LevelData<FArrayBox>& a_beta_coefficient,
                                const PotentialBC&    a_bc,
                                const bool            a_poloidally_average )
{
   bool fourthOrder = (m_discretization_order == 4);

   constructMatrix(a_alpha_coefficient, a_diffusion_coefficient, a_beta_coefficient, a_bc, m_A_graph,
                   m_A_stencil_size, m_A_stencil_values, m_A_radial, m_A_diagonal, m_A_diagonal_offset,
                   m_A_unstructured_coupling, a_poloidally_average, fourthOrder, m_dropOrder, m_A);

   setZero(m_diagonal_increment);
}



void
MBHypreSolver::solveWithMultigrid( const LevelData<FArrayBox>& a_rhs,
                                   const double                a_tol,
                                   const int                   a_max_iter,
                                   const bool                  a_verbose,
                                   LevelData<FArrayBox>&       a_solution )
{
   copyToHypreVector(a_rhs, m_b);
   copyToHypreVector(a_solution, m_x);

   /* This is a collective call finalizing the vector assembly.
      The vectors are now ``ready to be used'' */
   HYPRE_SStructVectorAssemble(m_b);
   HYPRE_SStructVectorAssemble(m_x);

   MultigridSolve_Hypre( m_A, m_b, a_tol, a_max_iter, a_verbose, m_x );

   copyFromHypreVector(m_x, a_solution);

   a_solution.exchange();
}



int
MBHypreSolver::solveWithGMRES( const LevelData<FArrayBox>& a_rhs,
                               const double                a_tol,
                               const int                   a_max_iter,
                               const double                a_amg_tol,
                               const int                   a_amg_max_iter,
                               const bool                  a_verbose, 
                               LevelData<FArrayBox>&       a_solution )
{
   copyToHypreVector(a_rhs, m_b);
   copyToHypreVector(a_solution, m_x);

   /* This is a collective call finalizing the vector assembly.
      The vectors are now ``ready to be used'' */
   HYPRE_SStructVectorAssemble(m_b);
   HYPRE_SStructVectorAssemble(m_x);

   int num_iters = GMRESSolve_Hypre( m_A, m_A, m_b, a_tol, a_max_iter, a_amg_tol, a_amg_max_iter, a_verbose, m_x );

   copyFromHypreVector(m_x, a_solution);

   a_solution.exchange();

   return num_iters;
}


void
MBHypreSolver::dumpMatrix(const string& a_file_name) const
{
   HYPRE_SStructMatrixPrint(a_file_name.c_str(), m_A, 0);
}


void
MBHypreSolver::createHypreData()
{
   if (m_hypre_allocated) {
      MayDay::Error( "GKPoisson:createHypreData:: data already allocated; deallocate first!" );
   }

   m_hypre_object_type = HYPRE_PARCSR;

   IntVect stencil_box_lo(IntVect::Zero);
   IntVect stencil_box_hi;
   if (m_discretization_order == 4) {
      if (m_dropOrder) {
         stencil_box_hi = 4*IntVect::Unit;
      }
      else {
         stencil_box_hi = 6*IntVect::Unit;
      }
   }
   else {
      stencil_box_hi = 2*IntVect::Unit;
   }
   Box stencil_box(stencil_box_lo, stencil_box_hi);
          
   m_A_stencil_values.define(stencil_box, 1);

   const DisjointBoxLayout & grids = m_geometry.grids();

   MultiBlockCoordSys* coord_sys_ptr = m_geometry.coordSysPtr();
   int num_blocks = coord_sys_ptr->numBlocks();

   int nvar = 1;
   int var = 0;

   // Create the grid

   {
      HYPRE_SStructGridCreate(MPI_COMM_WORLD, SpaceDim, num_blocks, &m_grid);

      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         const Box & box = grids[dit];
         int block_number = coord_sys_ptr->whichBlock(box);
         IntVect lower(box.loVect());
         IntVect upper(box.hiVect());

         HYPRE_SStructGridSetExtents(m_grid, block_number, lower.dataPtr(), upper.dataPtr());
      }

      HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};

      for (int block=0; block<num_blocks; ++block) {
         HYPRE_SStructGridSetVariables(m_grid, block, nvar, vartypes);

         // Set periodic shift vector and check whether the domain is fully periodic
         const ProblemDomain& domain = ((const MagBlockCoordSys*)coord_sys_ptr->getCoordSys(block))->domain();
         IntVect periodic = domain.size();

         for (int dir=0; dir<SpaceDim; dir++) {
            bool this_dir_periodic = domain.isPeriodic(dir);
            if ( !this_dir_periodic ) {
               periodic[dir] = 0;
            }
            else {
               // This direction is periodic.  Make sure that the number of
               // cells in this direction is a power of 2, which the Hypre
               // multigrid solver requires.
               int ncells = periodic[dir];
               while (ncells%2 == 0) {
                  ncells /= 2;
               }
               if (ncells != 1) {
                  if (procID()==0) {
                     MayDay::Error( "GKPoisson: Periodic direction must be a factor of 2!" );
                  }
               }
            }

            HYPRE_SStructGridSetPeriodic(m_grid, block, periodic.dataPtr());
         }
      }

      if (m_discretization_order == 4) {
         if (m_dropOrder) {
            m_A_radius = 2;
         }
         else {
            m_A_radius = 3;
         }
      }
      else {
         m_A_radius = 1;
      }

      int num_ghost[] = {m_A_radius, m_A_radius, m_A_radius, m_A_radius, m_A_radius, m_A_radius};

      HYPRE_SStructGridSetNumGhost(m_grid, num_ghost);

      HYPRE_SStructGridAssemble(m_grid);
   }

   // Define the discretization stencil for the system matrix A

   {

      /* Define the geometry of the stencil. Each represents a
         relative offset (in the index space). */

      if (m_discretization_order == 4) {
         if (m_dropOrder) {
            m_A_stencil_size = 25;
         }
         else {
            m_A_stencil_size = 49;
         }
      }
      else {
         m_A_stencil_size = 9;
      }

      int offsets2[][2] = {{-1,-1}, {0,-1}, {1,-1},
                           {-1,0},  {0,0},  {1,0},
                           {-1,1},  {0,1},  {1,1}};

      int offsets3[][2] = {{-2,-2}, {-1,-2}, {0,-2}, {1,-2}, {2,-2},
                           {-2,-1}, {-1,-1}, {0,-1}, {1,-1}, {2,-1},
                           {-2, 0}, {-1, 0}, {0, 0}, {1, 0}, {2, 0},
                           {-2, 1}, {-1, 1}, {0, 1}, {1, 1}, {2, 1},
                           {-2, 2}, {-1, 2}, {0, 2}, {1, 2}, {2, 2}};

      int offsets4[][2] = {{-3,-3}, {-2,-3}, {-1,-3}, {0,-3}, {1,-3}, {2,-3}, {3,-3},
                           {-3,-2}, {-2,-2}, {-1,-2}, {0,-2}, {1,-2}, {2,-2}, {3,-2},
                           {-3,-1}, {-2,-1}, {-1,-1}, {0,-1}, {1,-1}, {2,-1}, {3,-1},
                           {-3, 0}, {-2, 0}, {-1, 0}, {0, 0}, {1, 0}, {2, 0}, {3, 0},
                           {-3, 1}, {-2, 1}, {-1, 1}, {0, 1}, {1, 1}, {2, 1}, {3, 1},
                           {-3, 2}, {-2, 2}, {-1, 2}, {0, 2}, {1, 2}, {2, 2}, {3, 2},
                           {-3, 3}, {-2, 3}, {-1, 3}, {0, 3}, {1, 3}, {2, 3}, {3, 3}};

      int offsets[m_A_stencil_size][2];
      for (int n=0; n<m_A_stencil_size; ++n) {
         for (int m=0; m<2; ++m) {
            if (m_discretization_order == 4) {
               if (m_dropOrder) {
                  offsets[n][m] = offsets3[n][m];
               }
               else {
                  offsets[n][m] = offsets4[n][m];
               }
            }
            else {
               offsets[n][m] = offsets2[n][m];
            }
         }
      }

      m_A_diagonal_offset = (m_A_stencil_size - 1) / 2;

      /* Create an empty 2D stencil object */
      HYPRE_SStructStencilCreate(SpaceDim, m_A_stencil_size, &m_A_stencil);

      /* Assign each of the stencil entries */
      int entry;
      for (entry = 0; entry < m_A_stencil_size; entry++)
        HYPRE_SStructStencilSetEntry(m_A_stencil, entry, offsets[entry], var);
   }

   // Set up the graph

   {
      HYPRE_SStructGraphCreate(MPI_COMM_WORLD, m_grid, &m_A_graph);

      HYPRE_SStructGraphSetObjectType(m_A_graph, m_hypre_object_type);

      for (int block=0; block<num_blocks; ++block) {
         HYPRE_SStructGraphSetStencil(m_A_graph, block, var, m_A_stencil);
      }

      getUnstructuredCouplings(m_A_radius, m_A_unstructured_coupling);
      addUnstructuredGraphEntries(m_A_radius, m_A_unstructured_coupling, m_A_graph);

      HYPRE_SStructGraphAssemble(m_A_graph);
   }

   // Set up vectors for b and x

   {
      /* Create an empty vector object */
      HYPRE_SStructVectorCreate(MPI_COMM_WORLD, m_grid, &m_b);
      HYPRE_SStructVectorCreate(MPI_COMM_WORLD, m_grid, &m_x);

      HYPRE_SStructVectorSetObjectType(m_b, m_hypre_object_type);
      HYPRE_SStructVectorSetObjectType(m_x, m_hypre_object_type);

      HYPRE_SStructVectorInitialize(m_b);
      HYPRE_SStructVectorInitialize(m_x);
   }

   m_hypre_allocated = true;
}



void
MBHypreSolver::destroyHypreData()
{
   if (m_hypre_allocated) {
      HYPRE_SStructGridDestroy(m_grid);
      HYPRE_SStructStencilDestroy(m_A_stencil);
      if (m_A) {
         HYPRE_SStructMatrixDestroy(m_A);
         m_A = NULL;
      }
      HYPRE_SStructVectorDestroy(m_b);
      HYPRE_SStructVectorDestroy(m_x);

      m_hypre_allocated = false;
   }
}



void
MBHypreSolver::addUnstructuredGraphEntries( const int                                a_radius,
                                            const LayoutData< BaseFab<IntVectSet> >& a_unstructured_coupling,
                                            HYPRE_SStructGraph&                      a_graph ) const
{
   if (m_mblex_potential_Ptr) {

      const Box stencil_offsets(IntVect(D_DECL6(-a_radius,-a_radius,-a_radius,-a_radius,-a_radius,-a_radius)),
                                IntVect(D_DECL6(a_radius,a_radius,a_radius,a_radius,a_radius,a_radius)));

      const MagCoordSys* coord_sys_ptr = m_geometry.getCoordSys();
      const DisjointBoxLayout & grids = m_geometry.grids();

      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         int dst_block_number = coord_sys_ptr->whichBlock(grids[dit]);

         // Find all cells that may be coupled through a block interface
         IntVectSet dst_ivs = getInterBlockCoupledCells(dst_block_number, a_radius, grids[dit]);

         for (IVSIterator it(dst_ivs); it.ok(); ++it) {
            IntVect iv = it();
            Box stencil_box = stencil_offsets + iv;

            const IntVectSet& unstructured_couplings = a_unstructured_coupling[dit](iv);
            for (IVSIterator sivsit(unstructured_couplings); sivsit.ok(); ++sivsit) {
               IntVect siv = sivsit();
               int src_block_number = coord_sys_ptr->whichBlock(siv);
               CH_assert(src_block_number >= 0 && src_block_number < m_geometry.coordSysPtr()->numBlocks());
               HYPRE_SStructGraphAddEntries(a_graph, dst_block_number, iv.dataPtr(), 0,
                                            src_block_number, siv.dataPtr(), 0);
            }
         }
      }
   }
}



int
MBHypreSolver::findHypreEntry(const Box&        a_stencil_box,
                              const IntVectSet& a_unstructured_ivs,
                              const IntVect&    a_iv) const
{
   /*
     The Hypre matrix coefficients corresponding to the equation centered at cell a_iv are
     identified by an entry number, which is returned by this function.  Regular "stencil"
     couplings to cells in a_stencil_box are enumerated first (in the order with which a 
     BoxIterator built on a_stencil_box traverses cells), followed by the unstructured couplings
     to cells in other blocks contained in a_unstructured_ivs.
   */

   bool found_entry = false;

   int entry = 0;
   BoxIterator bit(a_stencil_box);
   for (bit.begin(); bit.ok(); ++bit) {
      if (bit() == a_iv) {
         found_entry = true;
         break;
      }
      entry++;
   }

   if (!found_entry) {
      for (IVSIterator ivit(a_unstructured_ivs); ivit.ok(); ++ivit) {
         if (a_iv == ivit()) {
            found_entry = true;
            break;
         }
         entry++;
      }
      CH_assert(found_entry);
   }

   return entry;
}



void
MBHypreSolver::addUnstructuredMatrixEntries(const LevelData<FArrayBox>&              a_alpha_coefficient,
                                            const LevelData<FluxBox>&                a_diffusion_coefficient,
                                            const PotentialBC&                       a_bc,
                                            FArrayBox&                               a_stencil_values,
                                            const bool                               a_fourthOrder,
                                            const bool                               a_dropOrder,
                                            const LayoutData< BaseFab<IntVectSet> >& a_unstructured_coupling,
                                            HYPRE_SStructMatrix&                     a_matrix) const
{
   if (m_mblex_potential_Ptr) {

      Box stencil_offsets = a_stencil_values.box();
      int stencil_size = stencil_offsets.numPts();
      int radius = (stencil_offsets.size(0)-1)/2;
      stencil_offsets -= radius*IntVect::Unit;

      const DisjointBoxLayout & grids = m_geometry.grids();
      const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = mag_coord_sys->boundaries();

      const LayoutData< IntVectSet >& ghosts = m_mblex_potential_Ptr->ghostCells();
      const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencil = m_mblex_potential_Ptr->stencils();

      BlockBaseRegister<BaseFab<SparseCoupling> > blockRegister(mag_coord_sys, grids, 0);

      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         const Box & box = grids[dit];

         int block_number = mag_coord_sys->whichBlock(box);

         const MagBlockCoordSys* block_coord_sys = mag_coord_sys->getCoordSys(block_number);
         RealVect dx = block_coord_sys->dx();

         const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
         bool lo_block_interface[SpaceDim];
         bool hi_block_interface[SpaceDim];
         for (int dir=0; dir<SpaceDim; ++dir) {
            lo_block_interface[dir] = this_block_boundaries[dir].isInterface();
            hi_block_interface[dir] = this_block_boundaries[dir + SpaceDim].isInterface();
         }

         int block_bc_types[2*SpaceDim];
         for (int dir=0; dir<SpaceDim; ++dir) {
            for (int side=0; side<2; ++side) {
               if ( this_block_boundaries[dir+SpaceDim*side].isDomainBoundary() ) {
                  block_bc_types[dir+2*side] = a_bc.getBCType(block_number, dir, side);
               }
               else {
                  block_bc_types[dir+SpaceDim*side] = -1;
               }
            }
         }

         const FArrayBox* alphaPtr = NULL;
         if ( a_alpha_coefficient.isDefined() ) {
            alphaPtr = &a_alpha_coefficient[dit];
         }

         const Box& domainBox = block_coord_sys->domain().domainBox();
         const IntVectSet& this_ghosts_fab = ghosts[dit];
         const IVSFAB<MBStencil>& this_stencil = *stencil[dit];
         const FArrayBox& this_volume = m_volume[dit];

         BaseFab<SparseCoupling> data[SpaceDim];
         for( int dir = 0; dir < SpaceDim; ++dir ) {
            data[dir].define(surroundingNodes(box, dir), 1);
         }

         IntVectSet interblock_ivs = getInterBlockCoupledCells(block_number, radius, box);
         for (IVSIterator ivsit(interblock_ivs); ivsit.ok(); ++ivsit) {
            IntVect iv = ivsit();
            Box stencil_box = stencil_offsets + iv;
            Box stencil_box_valid = stencil_box & domainBox;
            double this_cell_volume = this_volume(iv);
            IntVectSet extra_block_ghosts = stencil_box & this_ghosts_fab;

            double alpha = alphaPtr? alphaPtr->operator()(iv): 1.;

            // Get all of the indices that are sparsely coupled to this cell
            const IntVectSet& unstructured_couplings = a_unstructured_coupling[dit](iv);

            for (int dir=0; dir<SpaceDim; ++dir) {
               bool at_lo_block_interface = lo_block_interface[dir] && (iv[dir] == domainBox.smallEnd(dir));
               bool at_hi_block_interface = hi_block_interface[dir] && (iv[dir] == domainBox.bigEnd(dir));
               
               for (int side=0; side<2; ++side) {

                  IntVect iv_face;

                  // Halve the contributions on face contained in block interfaces to effect averaging
                  double fac;
                  if ( (at_lo_block_interface && side == 0) || (at_hi_block_interface && side == 1) ) {
                     fac = 0.5;
                     iv_face = iv;
                     iv_face.shift(dir,side);
                  }
                  else {
                     fac = 1.;
                  }
               
                  // Get the matrix couplings to the cell iv associated with the current side
                  a_stencil_values.setVal(0.);
                  accumStencilMatrixEntries(iv, dir, side, dir, this_block_boundaries, block_bc_types, a_diffusion_coefficient[dit],
                                            dx, domainBox, a_fourthOrder, a_dropOrder, a_stencil_values);

                  BoxIterator bit(stencil_box);
                  for (bit.begin(); bit.ok(); ++bit) {
                     IntVect stencil_box_iv = bit();

                     double s = fac * a_stencil_values( stencil_box_iv - iv + radius*IntVect::Unit );

                     if ( s != 0. ) {

                        if ( stencil_box_valid.contains(stencil_box_iv) ) {
                           // If cell iv is on a block interface, add the stencil contribution to the BlockRegister data.
                           // The minus sign accounts for the normal component, which is negated on the adacent cell.
                           if ( (at_lo_block_interface && side == 0) || (at_hi_block_interface && side == 1) ) {
                              data[dir](iv_face).add(stencil_box_iv, -s);

                              int entry = findHypreEntry(stencil_box, unstructured_couplings, stencil_box_iv);
                              double value = s * alpha / this_cell_volume;

                              int num_entries = 1;
                              HYPRE_SStructMatrixAddToValues(a_matrix, block_number, iv.dataPtr(),
                                                             0, num_entries, &entry, &value);

                           }
                        }
                        else if ( extra_block_ghosts.contains(stencil_box_iv) ) {

                           const MBStencil& elements = this_stencil(stencil_box_iv, 0);
                           int num_unstructured_entries = elements.size();
                           double * values = new double[num_unstructured_entries];
                           int * entries = new int[num_unstructured_entries];

                           MBStencilIterator stit(elements);
                           int k = 0;
                           for (stit.begin(); stit.ok(); ++stit) {
                              const MBStencilElement& stencilElement = stit();
                              IntVect interp_cell = stencilElement.cell();

                              entries[k] = findHypreEntry(stencil_box, unstructured_couplings, interp_cell);
                              values[k] = s * stencilElement.weight();

                              // If we are on a block interface, add the stencil contribution to the BlockRegister data.
                              // The minus sign accounts for the normal component, which is negated on the adacent cell.
                              if ( (at_lo_block_interface && side == 0) || (at_hi_block_interface && side == 1) ) {
                                 data[dir](iv_face).add(interp_cell, -values[k]);
                              }

                              values[k] *= (alpha / this_cell_volume);

                              k++;
                           }

                           HYPRE_SStructMatrixAddToValues(a_matrix, block_number, iv.dataPtr(),
                                                          0, num_unstructured_entries, entries, values);

                           delete [] entries;
                           delete [] values;

                        }
                        else {
                           MayDay::Error("Problem adding extra block matrix entries");
                        }
                     }
                  }
               }
            }
         }

         // Store the contributions at block interfaces
         for (int dir = 0; dir < SpaceDim; dir++) {
            for (SideIterator sit; sit.ok(); ++sit) {
               Side::LoHiSide side = sit();
               if (blockRegister.hasInterface(dit(), dir, side)) {
                  blockRegister.store(data[dir], dit(), dir, side);
               }
            }
         }
      }

      // Exchange the block interface contributions and accumulate them to the matrix

      blockRegister.exchange();

      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         int block_number = mag_coord_sys->whichBlock(grids[dit]);
         const FArrayBox& this_volume = m_volume[dit];

         const FArrayBox* alphaPtr = NULL;
         if ( a_alpha_coefficient.isDefined() ) {
            alphaPtr = &a_alpha_coefficient[dit];
         }

         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
            for (int dir = 0; dir < SpaceDim; dir++) {
               if (blockRegister.hasInterface(dit(), dir, side)) {
                  Box faceBox = adjCellBox(grids[dit], dir, side, 1);
                  // if Lo, then shift +1; if Hi, then shift -1
                  faceBox.shiftHalf(dir, -sign(side));
                  BaseFab<SparseCoupling> otherFab(faceBox, 1);
                  blockRegister.fill(otherFab, dit(), dir, side, side);
                  BoxIterator bit(faceBox);
                  for (bit.begin(); bit.ok(); ++bit) {
                     IntVect iv_face = bit();
                     IntVect iv(iv_face);
                     iv.shift(dir,-side);
                     Box stencil_box = stencil_offsets + iv;
                     const IntVectSet& unstructured_couplings = a_unstructured_coupling[dit](iv);

                     double alpha = alphaPtr? alphaPtr->operator()(iv): 1.;

                     SparseCoupling& coupling = otherFab(iv_face,0);
                     int num_entries = coupling.size();
                     double * values = new double[num_entries];
                     int * entries = new int[num_entries];

                     int k = 0;
                     for (SparseCouplingIterator it(coupling); it.ok(); ++it) {
                        IntVect index = coupling[it()];
                        entries[k] = findHypreEntry(stencil_box, unstructured_couplings, index);
                        values[k] = coupling.weight(index) * alpha / this_volume(iv);
                        k++;
                     }  // it

                     HYPRE_SStructMatrixAddToValues(a_matrix, block_number, iv.dataPtr(),
                                                    0, num_entries, entries, values);

                     delete [] entries;
                     delete [] values;

                  }  // bit
               }  // if interface
            }  // dir
         }  // sit
      }  // dit
   }
}



void
MBHypreSolver::constructMatrix( LevelData<FArrayBox>&              a_alpha_coefficient, 
                                LevelData<FluxBox>&                a_diffusion_coefficient,
                                LevelData<FArrayBox>&              a_beta_coefficient,
                                const PotentialBC&                 a_bc,
                                HYPRE_SStructGraph&                a_graph,
                                const int                          a_stencil_size,
                                FArrayBox&                         a_stencil_values,
                                LevelData<FArrayBox>&              a_radial,
                                LevelData<FArrayBox>&              a_diagonal,
                                const int                          a_diagonal_offset,
                                LayoutData< BaseFab<IntVectSet> >& a_unstructured_coupling,
                                const bool                         a_average_poloidally,
                                const bool                         a_fourthOrder,
                                const bool                         a_dropOrder,
                                HYPRE_SStructMatrix&               a_matrix ) const
{
   if (a_matrix == NULL) {
      HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, a_graph, &a_matrix);
      HYPRE_SStructMatrixSetObjectType(a_matrix, m_hypre_object_type);
   }

   HYPRE_SStructMatrixInitialize(a_matrix);

   int var = 0;

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = mag_coord_sys->boundaries();

   // Make sure the coefficient ghost cell values are set.
   // The ghost cell values at physical boundaries aren't used.
   a_diffusion_coefficient.exchange();

   /* Set the matrix coefficients.  Each processor assigns coefficients
      for the boxes in the grid that it owns. Note that the coefficients
      associated with each stencil entry may vary from grid point to grid
      point if desired.  Here, we first set the same stencil entries for
      each grid point.  Then we make modifications to grid points near
      the boundary. */

   int * entries = new int[a_stencil_size];
   for (int n=0; n<a_stencil_size; n++) entries[n] = n;
    
   double * stencil_values = a_stencil_values.dataPtr(0);

   const DisjointBoxLayout & grids = a_diffusion_coefficient.disjointBoxLayout();

   LevelData<FArrayBox> structured_values(grids, a_stencil_size, IntVect::Zero);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     const Box & box = grids[dit];

     int block_number = mag_coord_sys->whichBlock(box);
     const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
     RealVect dx = block_coord_sys.dx();

     const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
     bool lo_block_interface[SpaceDim];
     bool hi_block_interface[SpaceDim];
     if (m_mblex_potential_Ptr) {
        for (int dir=0; dir<SpaceDim; ++dir) {
           lo_block_interface[dir] = this_block_boundaries[dir].isInterface();
           hi_block_interface[dir] = this_block_boundaries[dir + SpaceDim].isInterface();
        }
     }
     else {
        for (int dir=0; dir<SpaceDim; ++dir) {
           lo_block_interface[dir] = false;
           hi_block_interface[dir] = false;
        }
     }

     int block_bc_types[2*SpaceDim];
     for (int dir=0; dir<SpaceDim; ++dir) {
        for (int side=0; side<2; ++side) {
           if ( this_block_boundaries[dir+side*SpaceDim].isDomainBoundary() ) {
              block_bc_types[dir+2*side] = a_bc.getBCType(block_number, dir, side);
           }
        }
     }

     const Box& domainBox = block_coord_sys.domain().domainBox();

     FluxBox& this_mapped_coefs = a_diffusion_coefficient[dit];
     const FArrayBox & this_volume = m_volume[dit];

     const FArrayBox* alpha = NULL;
     if ( a_alpha_coefficient.isDefined() ) {
        alpha = &a_alpha_coefficient[dit];
     }

     const FArrayBox* beta = NULL;
     if ( a_beta_coefficient.isDefined() ) {
        beta = &a_beta_coefficient[dit];
     }

     // Set the intra-block matrix entries

     BoxIterator bit(box);
     for (bit.begin();bit.ok();++bit) {
       IntVect iv = bit();

       a_stencil_values.setVal(0.);

       for (int dir=0; dir<SpaceDim; ++dir) {
          bool at_lo_interface = (lo_block_interface[dir] && iv[dir] == domainBox.smallEnd(dir));
          bool at_hi_interface = (hi_block_interface[dir] && iv[dir] == domainBox.bigEnd(dir));
          for (int side=0; side<2; ++side) {
             if( !(at_lo_interface && side == 0) && !(at_hi_interface && side == 1) ) {
                accumStencilMatrixEntries(iv, dir, side, dir, this_block_boundaries, block_bc_types, this_mapped_coefs,
                                          dx, domainBox, a_fourthOrder, a_dropOrder, a_stencil_values);
             }
          }
       }

       a_stencil_values /= this_volume(iv);

       if ( alpha ) {
          a_stencil_values *= alpha->operator()(iv);
       }

       if ( beta ) {
          stencil_values[a_diagonal_offset] += beta->operator()(iv);
       }

       for (int jj=0; jj<a_stencil_size; jj++) {
          structured_values[dit](iv,jj) = stencil_values[jj];
       }
     }
   }

   if ( a_average_poloidally ) {
      FluxSurface flux_surface_op(m_geometry);

      LevelData<FArrayBox> fs_average(flux_surface_op.grids(), a_stencil_size, IntVect::Zero);
      flux_surface_op.average(structured_values, fs_average);
      flux_surface_op.spread(fs_average, structured_values);
   }

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     const Box & box = grids[dit];
     IntVect lower(box.loVect());
     IntVect upper(box.hiVect());

     int block_number = mag_coord_sys->whichBlock(box);

     // Set the intra-block matrix entries

     int nvalues  = a_stencil_size * box.numPts();
     double * values = new double[nvalues];

     double * radial_ldiag = a_radial[dit].dataPtr(0);
     double * radial_diag  = a_radial[dit].dataPtr(1);
     double * radial_udiag = a_radial[dit].dataPtr(2);

     double * diagonal_values = a_diagonal[dit].dataPtr(0);

     int ii = 0;
     int diagonal_index = 0;

     BoxIterator bit(box);
     for (bit.begin();bit.ok();++bit) {
       IntVect iv = bit();

       for (int jj=0; jj<a_stencil_size; jj++) {
          stencil_values[jj] = structured_values[dit](iv,jj);
       }

       if ( a_average_poloidally ) {
         for (int jj=0; jj<a_stencil_size; jj++) {
           values[ii+jj] = 0.;
         }

#if 0
         values[ii+3] = radial_ldiag[diagonal_index] = stencil_values[3];
         values[ii+4] = radial_diag[diagonal_index]  = stencil_values[4];
         values[ii+5] = radial_udiag[diagonal_index] = stencil_values[5];
#else

         values[ii+3] = radial_ldiag[diagonal_index]
            = stencil_values[0] + stencil_values[3] + stencil_values[6];
         values[ii+4] = radial_diag[diagonal_index]
            = stencil_values[1] + stencil_values[4] + stencil_values[7];
         values[ii+5] = radial_udiag[diagonal_index]
            = stencil_values[2] + stencil_values[5] + stencil_values[8];
#endif
       }
       else {
         for (int jj=0; jj<a_stencil_size; jj++) {
            values[ii+jj] = stencil_values[jj];
         }
       }

       diagonal_values[diagonal_index] = stencil_values[a_diagonal_offset];

       //       cout << iv << " " << radial_ldiag[diagonal_index] << " " << radial_diag[diagonal_index] << " " << radial_udiag[diagonal_index] << endl;

       ii += a_stencil_size;
       diagonal_index++;

     }

     HYPRE_SStructMatrixSetBoxValues(a_matrix, block_number, lower.dataPtr(), upper.dataPtr(),
                                     var, a_stencil_size, entries, values);

     delete [] values;
   }

   delete [] entries;

   addUnstructuredMatrixEntries(a_alpha_coefficient, a_diffusion_coefficient, a_bc, a_stencil_values,
                                a_fourthOrder, a_dropOrder, a_unstructured_coupling, a_matrix);

   /* This is a collective call finalizing the matrix assembly.
      The matrix is now ``ready to be used'' */
   HYPRE_SStructMatrixAssemble(a_matrix);

   //   HYPRE_SStructMatrixPrint("struct.out.A", a_matrix, 0);
}



void
MBHypreSolver::MultigridSolve_Hypre( const HYPRE_SStructMatrix& a_matrix,
                                     const HYPRE_SStructVector& a_b,
                                     const double               a_tol,
                                     const int                  a_max_iter,
                                     const bool                 a_verbose,
                                     const HYPRE_SStructVector& a_x ) const
{
   HYPRE_Solver          par_solver;
   HYPRE_ParCSRMatrix    par_A;
   HYPRE_ParVector       par_b;
   HYPRE_ParVector       par_x;

   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);
   HYPRE_SStructVectorGetObject(a_b, (void **) &par_b);
   HYPRE_SStructVectorGetObject(a_x, (void **) &par_x);

   HYPRE_BoomerAMGCreate(&par_solver);
   HYPRE_BoomerAMGSetTol(par_solver, a_tol);
   HYPRE_BoomerAMGSetMaxIter(par_solver, a_max_iter);
   HYPRE_BoomerAMGSetPrintLevel(par_solver, 0);

   // Algorithm options
   HYPRE_BoomerAMGSetStrongThreshold(par_solver, 0.25);
   HYPRE_BoomerAMGSetCoarsenType(par_solver, 6);  // Falgout coarsening
   //   HYPRE_BoomerAMGSetRelaxType(par_solver, 3);  // hybrid Gauss-Seidel
   //   HYPRE_BoomerAMGSetCycleRelaxType(par_solver, 3, 3);  // hybrid Gauss-Seidel on coarsest level
   //   HYPRE_BoomerAMGSetCycleRelaxType(par_solver, 9, 3);  // Gaussian elimination on coarsest level
   //   HYPRE_BoomerAMGSetCycleType(par_solver, 1);  // V = 1, W = 2
   //   HYPRE_BoomerAMGSetCycleNumSweeps(par_solver, 1, 1);  // 1 sweep on down cycle
   //   HYPRE_BoomerAMGSetCycleNumSweeps(par_solver, 1, 2);  // 1 sweep on up cycle
   //   HYPRE_BoomerAMGSetCycleNumSweeps(par_solver, 1, 3);  // 1 sweeps on coarsest level

   HYPRE_BoomerAMGSetup(par_solver, par_A, par_b, par_x);
   HYPRE_BoomerAMGSolve(par_solver, par_A, par_b, par_x);

   int num_iterations;
   HYPRE_BoomerAMGGetNumIterations(par_solver, &num_iterations);
   double final_res_norm;
   HYPRE_BoomerAMGGetFinalRelativeResidualNorm(par_solver, &final_res_norm);
   if (a_verbose && procID()==0) {
      cout << "      AMG solver residual = " << final_res_norm << " after " << num_iterations << " iterations" << endl;
   }

   HYPRE_BoomerAMGDestroy(par_solver);

   if (m_hypre_object_type == HYPRE_PARCSR) {
      HYPRE_SStructVectorGather(a_x);
   }
}



int
MBHypreSolver::GMRESSolve_Hypre( const HYPRE_SStructMatrix&  a_matrix,
                                 const HYPRE_SStructMatrix&  a_precond,
                                 const HYPRE_SStructVector&  a_b,
                                 const double                a_tol,
                                 const int                   a_max_iter,
                                 const double                a_amg_tol,
                                 const int                   a_amg_max_iter,
                                 const bool                  a_verbose,
                                 const HYPRE_SStructVector&  a_x ) const
{
   HYPRE_Solver          par_solver;
   HYPRE_Solver          par_precond;
   HYPRE_ParCSRMatrix    par_A;
   HYPRE_ParCSRMatrix    par_P;
   HYPRE_ParVector       par_b;
   HYPRE_ParVector       par_x;

   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);
   HYPRE_SStructMatrixGetObject(a_precond, (void **) &par_P);
   HYPRE_SStructVectorGetObject(a_b, (void **) &par_b);
   HYPRE_SStructVectorGetObject(a_x, (void **) &par_x);

   HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &par_solver);

   /* set the GMRES paramaters */
   HYPRE_GMRESSetMaxIter(par_solver, a_max_iter);
   HYPRE_GMRESSetTol(par_solver, a_tol);
   HYPRE_GMRESSetPrintLevel(par_solver, 2);
   HYPRE_GMRESSetLogging(par_solver, 1);

   /* use BoomerAMG as preconditioner */
   HYPRE_BoomerAMGCreate(&par_precond);
   HYPRE_BoomerAMGSetCoarsenType(par_precond, 6);
   HYPRE_BoomerAMGSetStrongThreshold(par_precond, 0.25);
   HYPRE_BoomerAMGSetTol(par_precond, a_amg_tol);
   HYPRE_BoomerAMGSetPrintLevel(par_precond, 0);
   HYPRE_BoomerAMGSetPrintFileName(par_precond, "ex9.out.log");
   HYPRE_BoomerAMGSetMaxIter(par_precond, a_amg_max_iter);

   /* set the preconditioner */
   HYPRE_ParCSRGMRESSetPrecond(par_solver,
                               HYPRE_BoomerAMGSolve,
                               HYPRE_BoomerAMGSetup,
                               par_precond);

   HYPRE_GMRESSetup(par_solver, (HYPRE_Matrix)par_P, (HYPRE_Vector)par_b, (HYPRE_Vector)par_x);
   HYPRE_GMRESSolve(par_solver, (HYPRE_Matrix)par_A, (HYPRE_Vector)par_b, (HYPRE_Vector)par_x);

   // Get number of preconditioner iterations used
   int num_final_iterations;
   HYPRE_GMRESGetNumIterations(par_solver, &num_final_iterations);

   // Get final preconditioner residual norm
   double final_norm2;
   HYPRE_GMRESGetFinalRelativeResidualNorm(par_solver, &final_norm2);

   if (a_verbose && procID()==0) {
      cout << "      GMRES residual = " << final_norm2 << " after " << num_final_iterations << " iterations" << endl;
   }

   if (m_hypre_object_type == HYPRE_PARCSR) {
      HYPRE_SStructVectorGather(a_x);
   }

   return 0;
}



void
MBHypreSolver::copyToHypreVector( const LevelData<FArrayBox>& a_in,
                                  HYPRE_SStructVector&        a_out )
{
   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();

   int var = 0;

   const DisjointBoxLayout & grids = a_in.disjointBoxLayout();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = mag_coord_sys->whichBlock(grids[dit]);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      FArrayBox tmp(box,1);
      tmp.copy(a_in[dit]);

      HYPRE_SStructVectorSetBoxValues(a_out, block_number, lower.dataPtr(), upper.dataPtr(), var, tmp.dataPtr());
   }
} 



void
MBHypreSolver::copyFromHypreVector( const HYPRE_SStructVector& a_in,
                                    LevelData<FArrayBox>&      a_out ) const
{
   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();

   int var = 0;

   const DisjointBoxLayout & grids = a_out.disjointBoxLayout();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = mag_coord_sys->whichBlock(grids[dit]);

      FArrayBox tmp(box,1);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      HYPRE_SStructVectorGetBoxValues(a_in, block_number, lower.dataPtr(), upper.dataPtr(), var, tmp.dataPtr());

      a_out[dit].copy(tmp);
   }
} 



#include "NamespaceFooter.H"
