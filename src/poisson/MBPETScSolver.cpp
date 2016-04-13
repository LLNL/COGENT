#ifdef with_petsc

#include "MBPETScSolver.H"
#include "MBStencilIterator.H"
#include "BlockRegister.H"
#include "SparseCoupling.H"
#include "FluxSurface.H"

#include "NamespaceHeader.H"



MBPETScSolver::MBPETScSolver( const MultiBlockLevelGeom&  a_geom,
                              const int                   a_discretization_order )
   : MBSolver(a_geom, a_discretization_order),
     m_petsc_allocated(false),
     m_amg_solver_allocated(false),
     m_gmres_solver_allocated(false)
{
   createPETScData();
}
      


MBPETScSolver::~MBPETScSolver()
{
   destroyPETScData();
   destroySolver();
}



void
MBPETScSolver::constructMatrixGeneral( LevelData<FArrayBox>& a_alpha_coefficient,
                                       LevelData<FluxBox>&   a_tensor_coefficient,
                                       LevelData<FArrayBox>& a_beta_coefficient,
                                       const PotentialBC&    a_bc )
{
   bool fourthOrder = (m_discretization_order == 4);

   constructPETScMatrix(a_alpha_coefficient, a_tensor_coefficient, a_beta_coefficient, a_bc,
                        m_A_stencil_values, m_A_diagonal_offset,
                        m_A_unstructured_coupling, fourthOrder, m_dropOrder, m_A, m_rhs_from_bc);
}



void
MBPETScSolver::multiplyMatrix( const LevelData<FArrayBox>& a_in,
                               LevelData<FArrayBox>&       a_out ) const
{
   if (!m_petsc_allocated) {
      MayDay::Error( "MBPETScSolver::multiplyMatrix(): PETSc objects have not yet been initialized!" );
   }

   copyToPETScVector(a_in, m_x);

   MatMult(m_A, m_x, m_b);

   copyFromPETScVector(m_b, a_out);
}



void
MBPETScSolver::setParams( const string& a_method,
                          const double  a_method_tol,
                          const int     a_method_max_iter,
                          const bool    a_method_verbose,
                          const string& a_precond_method,
                          const double  a_precond_tol,
                          const int     a_precond_max_iter,
                          const bool    a_precond_verbose )
{
   m_method           = a_method;
   m_method_tol       = a_method_tol;
   m_method_max_iter  = a_method_max_iter;
   m_method_verbose   = a_method_verbose;
   m_precond_method   = a_precond_method;
   m_precond_tol      = a_precond_tol;
   m_precond_max_iter = a_precond_max_iter;
   m_precond_verbose  = a_precond_verbose;

   m_params_set = true;
}



void
MBPETScSolver::solve( const LevelData<FArrayBox>&  a_rhs,
                      LevelData<FArrayBox>&        a_solution,
                      bool                         a_homogeneous_bcs )
{
   if ( !m_params_set ) {
      MayDay::Error("MBPETScSolver::solve(): solver parameters have not been set");
   }

   if ( a_homogeneous_bcs ) {
      copyToPETScVector(a_rhs, m_b);
   }
   else {
      LevelData<FArrayBox> rhs;
      rhs.define(a_rhs);

      for (DataIterator dit(rhs.dataIterator()); dit.ok(); ++dit) {
         rhs[dit] += m_rhs_from_bc[dit];
      }

      copyToPETScVector(rhs, m_b);
   }

   copyToPETScVector(a_solution, m_x);

   if ( m_method == "GMRES" ) {

      if ( !m_gmres_solver_allocated ) {
         createGMRESSolver();
      }

      AMG_Preconditioned_GMRES( m_A, m_b, m_method_tol, m_method_max_iter, m_precond_tol, m_precond_max_iter, m_method_verbose, m_x );
   }
   else if ( m_method == "AMG" ) {

      if ( !m_amg_solver_allocated ) {
         createAMGSolver();
      }

      AMG( m_A, m_b, m_method_tol, m_method_max_iter, m_method_verbose, m_x );
   }
   else {
      MayDay::Error("MBPETScSolver::solve(): Unknown method, only GMRES or AMG recognized");
   }

   copyFromPETScVector(m_x, a_solution);
}



void
MBPETScSolver::dumpMatrix(const string& a_file_name) const
{
   PetscViewer viewer;
   PetscViewerCreate(MPI_COMM_WORLD, &viewer);
   
   PetscViewerASCIIOpen(MPI_COMM_WORLD, a_file_name.c_str(), &viewer);

   MatView(m_A, viewer);

   PetscViewerDestroy(&viewer);
}



void
MBPETScSolver::createPETScData()
{
   if (m_petsc_allocated) {
      MayDay::Error( "GKPoisson:createPETScData:: data already allocated; deallocate first!" );
   }

   // Set the stencil size data

   IntVect stencil_box_lo(IntVect::Zero);
   IntVect stencil_box_hi;
   int radius;
   if (m_discretization_order == 4) {
      if (m_dropOrder) {
         stencil_box_hi = 4*IntVect::Unit;
         radius = 2;
      }
      else {
         stencil_box_hi = 6*IntVect::Unit;
         radius = 3;
      }
   }
   else {
      stencil_box_hi = 2*IntVect::Unit;
      radius = 1;
   }
   Box stencil_box(stencil_box_lo, stencil_box_hi);
   m_A_stencil_values.define(stencil_box, 1);
          
   /* Define the geometry of the stencil. Each represents a
      relative offset (in the index space). */

   int stencil_size = stencil_box.numPts();

   m_A_diagonal_offset = (stencil_size - 1) / 2;

   /*
     PETSc uses global ids (i.e., equation/unknown numbers) to identify matrix and
     vector entries.  Moreover, it is required that global ids are assigned
     contiguously to processors.  In addition to the relatively easy task of
     assigning cell IntVects to global ids consistent with that convention, we
     require a mechanism for determining the inverse of that assignment, namely,
     finding the global id of an arbitrary IntVect.  For this purpose, it is not
     sufficient to simply have each processor store the global ids of cells
     contained in only the boxes assigned to it, or even the global ids of the
     ghost cells of those local boxes filled from other processors using an
     exchange operation.  This is due to the fact that the multiblock interpolation
     used to fill extrablock ghost cells needs to know the global ids of cells
     comprising the interpolation stencil, which in general can be owned by any
     number of neighbor blocks and processors.

     For scalability and efficiency, we obviously don't want every processor
     to have to store the global id of every cell in the entire problem.  We
     therefore adopt the strategy of having every processor maintain a list of
     the boxes contained in each block, together with a corresponding list of
     the starting global id assigned to the box.  Given an arbitray IntVect, the
     block to which it belongs is easily determined.  A search of the corresponding
     box list obtains the box to which the IntVect belongs, and its global id is
     then obtained as an offset to the box's starting global id.

   */

   const DisjointBoxLayout & grids = m_geometry.gridsFull();
   const int num_blocks = m_coord_sys_ptr->numBlocks();

   // For each block, make a list of the boxes it contains, the
   // corresponding processor numbers and the starting global ids

   m_block_boxes.resize(num_blocks);
   m_block_boxes_starting_gid.resize(num_blocks);

   const Vector<Box>& global_boxes = grids.boxArray();
   const Vector<int>& pids = grids.procIDs();

   Vector<int> num_gids(numProc(),0);
   Vector< Vector<int> > block_boxes_pid(num_blocks);

   for (int block=0; block<num_blocks; ++block) {
      const Box& domain_box = (m_coord_sys_ptr->mappingBlocks())[block];

      for (int i=0; i<global_boxes.size(); ++i) {

         Box overlap = global_boxes[i] & domain_box;
         if ( overlap.ok() ) {
            m_block_boxes[block].push_back(overlap);
            block_boxes_pid[block].push_back(pids[i]);
            num_gids[pids[i]] += overlap.numPts();
         }
      }
   }

   // At this point, num_gids[i] has the number of gids allocated
   // to processor i.

   Vector<int> proc_starting_gid(numProc());
   proc_starting_gid[0] = 0;
   for (int i=1; i<num_gids.size(); ++i) {
      proc_starting_gid[i] = proc_starting_gid[i-1] + num_gids[i-1];
   }

   m_starting_gid = proc_starting_gid[procID()];
   m_ending_gid = m_starting_gid + num_gids[procID()] - 1;

   // proc_starting_gid[i] is now the starting gid for the gids to
   // be assigned to processor i.  Next, use it to successively store the
   // starting gid of each box assigned to processor i.

   for (int block=0; block<num_blocks; ++block) {
      const Vector<Box>& block_boxes = m_block_boxes[block];
      const Vector<int>& box_pid = block_boxes_pid[block];
      Vector<int>& starting_gid = m_block_boxes_starting_gid[block];
      for (int i=0; i<block_boxes.size(); ++i) {
         int pid = box_pid[i];
         starting_gid.push_back(proc_starting_gid[pid]);
         proc_starting_gid[pid] += block_boxes[i].numPts();
      }
   }

   // Now that the starting global ids of all boxes have been determined,
   // use them to set the global ids of all of the cells in the boxes
   // owned by this processor.

   m_gids.define(grids, 1, radius*IntVect::Unit);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box& this_box = grids[dit];
      BaseFab<int>& this_gids = m_gids[dit];

      int block_number = m_coord_sys_ptr->whichBlock(this_box);
      const Vector<Box>& block_boxes = m_block_boxes[block_number];
      const Vector<int>& starting_gid = m_block_boxes_starting_gid[block_number];

      for (int i=0; i<block_boxes.size(); ++i) {
         const Box& block_box = block_boxes[i];
         if (block_box == this_box) {
            int gid = starting_gid[i];
            for (BoxIterator bit(this_box); bit.ok(); ++bit) {
               IntVect iv = bit();
               this_gids(iv) = gid++;
            }
         }
      }
   }
   m_gids.exchange();

   // Set up the matrix

   getUnstructuredCouplings(radius, m_A_unstructured_coupling);

   int* d_nnz;   // number of intraprocessor nonzeros
   int* o_nnz;   // number of interprocessor nonzeros

   int num_local_unknowns = numLocalGIDs();

   getNumNonzeros(radius, m_A_unstructured_coupling, &d_nnz, &o_nnz );

   MatCreateAIJ(MPI_COMM_WORLD, num_local_unknowns, num_local_unknowns,
                PETSC_DECIDE, PETSC_DECIDE, 0, d_nnz, 0, o_nnz, &m_A);

   delete [] o_nnz;
   delete [] d_nnz;

   // Set up the right-hand side b, the solution vector x and two
   // temporary vectors for use in multiplyMatrix() (to avoid
   // creating and destroying them each time).

   VecCreate( MPI_COMM_WORLD, &m_b );
   VecSetFromOptions( m_b );
   VecSetSizes( m_b, num_local_unknowns, PETSC_DECIDE );
   VecDuplicate( m_b, &m_x );
   VecDuplicate( m_b, &m_in_vector );
   VecDuplicate( m_b, &m_out_vector );

   m_petsc_allocated = true;
}



void
MBPETScSolver::destroyPETScData()
{
   if (m_petsc_allocated) {
      MatDestroy(&m_A);
      VecDestroy(&m_b);
      VecDestroy(&m_x);

      m_petsc_allocated = false;
   }
}



void
MBPETScSolver::createAMGSolver()
{
   if ( m_amg_solver_allocated ) {
      MayDay::Error("MBPETScSolver::createAMGSolver(): solver is already allocated.");
   }

   KSPCreate(MPI_COMM_WORLD, &m_amg_ksp);
   KSPSetType(m_amg_ksp, KSPRICHARDSON);

   //   PetscOptionsSetValue(NULL, "-ksp_monitor", NULL);
   KSPSetFromOptions(m_amg_ksp);

   KSPGetPC(m_amg_ksp, &m_pc);

   PCSetType(m_pc, PCHYPRE);
   PCHYPRESetType(m_pc, "boomeramg");

   PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_coarsen_type", "Falgout");
   PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_relax_type_all", "Jacobi");
   PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_max_iter", "5");

   PCSetFromOptions(m_pc);

   KSPSetOperators(m_amg_ksp, m_A, m_A);

   m_amg_solver_allocated = true;
}



void
MBPETScSolver::createGMRESSolver()
{
   if ( m_gmres_solver_allocated ) {
      MayDay::Error("MBPETScSolver::createAMGSolver(): solver is already allocated.");
   }

   KSPCreate(MPI_COMM_WORLD, &m_gmres_ksp);
   KSPSetType(m_gmres_ksp, KSPGMRES);

   KSPGetPC(m_gmres_ksp, &m_pc);

   PCSetType(m_pc, PCHYPRE);
   PCHYPRESetType(m_pc, "boomeramg");

   PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_coarsen_type", "Falgout");
   PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_relax_type_all", "Jacobi");
   PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_max_iter", "5");

   PCSetFromOptions(m_pc);

   KSPSetOperators(m_gmres_ksp, m_A, m_A);

   m_gmres_solver_allocated = true;
}



void
MBPETScSolver::destroySolver()
{
   if (m_amg_solver_allocated) {
      KSPDestroy(&m_amg_ksp);
      m_amg_solver_allocated = false;
   }
   if (m_gmres_solver_allocated) {
      KSPDestroy(&m_gmres_ksp);
      m_gmres_solver_allocated = false;
   }
}



int
MBPETScSolver::getGID( const IntVect& a_iv ) const
{
   return getGID( m_coord_sys_ptr->whichBlock(a_iv), a_iv);
}



int
MBPETScSolver::getGID( const int      a_block_number,
                       const IntVect& a_iv ) const
{
   int gid = -1;

   const Vector<Box>& block_boxes = m_block_boxes[a_block_number];
   const Vector<int>& starting_gid = m_block_boxes_starting_gid[a_block_number];

   bool found_gid = false;
   for (int i=0; i<block_boxes.size(); ++i) {
      const Box& block_box = block_boxes[i];
      if ( block_box.contains(a_iv) ) {
         gid = starting_gid[i];
         for (BoxIterator bit(block_box); bit.ok(); ++bit) {
            if ( a_iv == bit() ) {
               found_gid = true;
               break;
            }
            else {
               gid++;
            }
         }
         break;
      }
   }
   
   if ( !found_gid ) {
      MayDay::Error("MBPETScSolver::getGID(): Could not find gid");
   }

   return gid;
}



void
MBPETScSolver::getNumNonzeros( const int                                a_radius,
                               const LayoutData< BaseFab<IntVectSet> >& a_unstructured_couplings,
                               int**                                    a_d_nnz,
                               int**                                    a_o_nnz ) const
{
   const DisjointBoxLayout & grids = m_geometry.gridsFull();

   int num_local_unknowns = numLocalGIDs();

   *a_d_nnz = new int[num_local_unknowns];
   *a_o_nnz = new int[num_local_unknowns];

   for (int i=0; i<num_local_unknowns; ++i) {
      (*a_d_nnz)[i] = (*a_o_nnz)[i] = 0;
   }

   Box template_box(IntVect(D_DECL6(-a_radius,-a_radius,-a_radius,-a_radius,-a_radius,-a_radius)),
                   IntVect(D_DECL6(a_radius,a_radius,a_radius,a_radius,a_radius,a_radius)));

   IntVect periodic = checkPeriodic();

   const Vector<Box>& mapping_blocks = m_coord_sys_ptr->mappingBlocks();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box& box = grids[dit];
      int block_number = m_coord_sys_ptr->whichBlock(box);
      const BaseFab<int>& this_gids = m_gids[dit];
         
      Box domain_box = mapping_blocks[block_number];

      for (int dir=0; dir<SpaceDim; ++dir) {
         if ( periodic[dir] != 0 ) {
            domain_box.grow(dir,a_radius);
         }
      }

      // Count the structured stencil entries
      for (BoxIterator bit(box); bit.ok(); ++bit) {
         IntVect iv = bit();
         Box valid_box = (template_box + iv) & domain_box;
         int index = this_gids(iv) - m_starting_gid;

         for (BoxIterator vit(valid_box); vit.ok(); ++vit) {
            IntVect iv2 = vit();
            int gid = this_gids(iv2);

            if ( gid >= m_starting_gid && gid <= m_ending_gid ) {
               (*a_d_nnz)[index]++;
            }
            else {
               (*a_o_nnz)[index]++;
            }
         }
      }

      if (m_mblex_potential_Ptr) {

         // Find all cells that may be coupled through a block interface
         IntVectSet ivs = getInterBlockCoupledCells(block_number, a_radius, grids[dit]);

         const BaseFab<IntVectSet>& this_unstructured_coupling = a_unstructured_couplings[dit];

         // Count the unstructured stencil entries
         for (IVSIterator it(ivs); it.ok(); ++it) {
            IntVect iv = it();
            int index = this_gids(iv) - m_starting_gid;

            for (IVSIterator ivit(this_unstructured_coupling(iv)); ivit.ok(); ++ivit) {
               int gid = getGID(ivit());

               if ( gid >= m_starting_gid && gid <= m_ending_gid ) {
                  (*a_d_nnz)[index]++;
               }
               else {
                  (*a_o_nnz)[index]++;
               }
            }
         }
      }
   }
}



void
MBPETScSolver::addUnstructuredMatrixEntries(const LevelData<FArrayBox>&              a_alpha_coefficient,
                                            const LevelData<FluxBox>&                a_tensor_coefficient,
                                            const PotentialBC&                       a_bc,
                                            FArrayBox&                               a_stencil_values,
                                            const bool                               a_fourthOrder,
                                            const bool                               a_dropOrder,
                                            const LayoutData< BaseFab<IntVectSet> >& a_unstructured_coupling,
                                            Vector< Vector<CoDim1Stencil> >&         a_codim1_stencils,
                                            Vector< Vector<CoDim2Stencil> >&         a_codim2_stencils,
                                            Mat&                                     a_matrix) const
{
   if (m_mblex_potential_Ptr) {

      Box stencil_offsets = a_stencil_values.box();
      int radius = (stencil_offsets.size(0)-1)/2;
      stencil_offsets -= radius*IntVect::Unit;

      const DisjointBoxLayout & grids = m_geometry.gridsFull();
      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys_ptr->boundaries();

      const LayoutData< IntVectSet >& ghosts = m_mblex_potential_Ptr->ghostCells();
      const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencil = m_mblex_potential_Ptr->stencils();

      BlockBaseRegister<BaseFab<SparseCoupling> > blockRegister(m_coord_sys_ptr, grids, 0);

      const Vector<Box>& mapping_blocks = m_coord_sys_ptr->mappingBlocks();

      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         const Box & box = grids[dit];
         const BaseFab<int>& this_gids = m_gids[dit];

         int block_number = m_coord_sys_ptr->whichBlock(box);

         const NewCoordSys* block_coord_sys = m_coord_sys_ptr->getCoordSys(block_number);
         RealVect dx = block_coord_sys->dx();

         const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
         bool lo_block_interface[SpaceDim];
         bool hi_block_interface[SpaceDim];
         for (int dir=0; dir<SpaceDim; ++dir) {
            lo_block_interface[dir] = this_block_boundaries[dir].isInterface();
            hi_block_interface[dir] = this_block_boundaries[dir + SpaceDim].isInterface();
         }

         const FArrayBox* alphaPtr = NULL;
         if ( a_alpha_coefficient.isDefined() ) {
            alphaPtr = &a_alpha_coefficient[dit];
         }

         const Box& domain_box = mapping_blocks[block_number];

         const IntVectSet& this_ghosts_fab = ghosts[dit];
         const IVSFAB<MBStencil>& this_stencil = *stencil[dit];

         BaseFab<SparseCoupling> data[SpaceDim];
         for( int dir = 0; dir < SpaceDim; ++dir ) {
            data[dir].define(surroundingNodes(box, dir), 1);
         }

         IntVectSet interblock_ivs = getInterBlockCoupledCells(block_number, radius, box);
         for (IVSIterator ivsit(interblock_ivs); ivsit.ok(); ++ivsit) {
            IntVect iv = ivsit();
            Box stencil_box = stencil_offsets + iv;
            Box stencil_box_valid = stencil_box & domain_box;
            IntVectSet extra_block_ghosts = stencil_box & this_ghosts_fab;

            double alpha = alphaPtr? alphaPtr->operator()(iv): 1.;

            for (int dir=0; dir<SpaceDim; ++dir) {
               bool at_lo_block_interface = lo_block_interface[dir] && (iv[dir] == domain_box.smallEnd(dir));
               bool at_hi_block_interface = hi_block_interface[dir] && (iv[dir] == domain_box.bigEnd(dir));
               bool force_codim2_condense = lo_block_interface[dir] || hi_block_interface[dir];
               
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
               
                  for (int dir2=0; dir2<SpaceDim; ++dir2) {

                     a_stencil_values.setVal(0.);

                     accumStencilMatrixEntries(iv, dir, side, dir2, a_tensor_coefficient[dit],
                                               dx, a_fourthOrder, a_dropOrder, a_stencil_values);

                     FArrayBox dummy;
                     modifyStencilForBCs( a_codim1_stencils[block_number], a_codim2_stencils[block_number],
                                          iv, a_stencil_values, dummy, false, force_codim2_condense );
                     
                     BoxIterator bit(stencil_box);
                     for (bit.begin(); bit.ok(); ++bit) {
                        IntVect stencil_box_iv = bit();
                        
                        double s = fac * a_stencil_values( stencil_box_iv - iv + radius*IntVect::Unit );

                        if ( s != 0. ) {
                           
                           if ( stencil_box_valid.contains(stencil_box_iv) ) {

                              // If the current side of cell iv is on a block interface, add the stencil contribution to the BlockRegister data.
                              // The minus sign accounts for the normal component, which is negated on the adacent cell.  Also, since we have not
                              // yet accounted for the regular stencil contributions, do that now.
                              if ( (at_lo_block_interface && side == 0) || (at_hi_block_interface && side == 1) ) {
                                 data[dir](iv_face).add(stencil_box_iv, -s);
                                 
                                 double value = s * alpha;
                                 
                                 int row = this_gids(iv);
                                 int col = getGID(block_number, stencil_box_iv);
                                 MatSetValue(a_matrix, row, col, value, ADD_VALUES);
                              }
                           }
                           else if ( extra_block_ghosts.contains(stencil_box_iv) ) {

                              const MBStencil& elements = this_stencil(stencil_box_iv, 0);
                              MBStencilIterator stit(elements);
                              for (stit.begin(); stit.ok(); ++stit) {
                                 const MBStencilElement& stencilElement = stit();
                                 IntVect interp_cell = stencilElement.cell();

                                 double value = s * stencilElement.weight();

                                 // If we are on a block interface, add the stencil contribution to the BlockRegister data.
                                 // The minus sign accounts for the normal component, which is negated on the adacent cell.
                                 if ( (at_lo_block_interface && side == 0) || (at_hi_block_interface && side == 1) ) {
                                    data[dir](iv_face).add(interp_cell, -value);
                                 }

                                 value *= alpha;
                                 
                                 int row = this_gids(iv);
                                 int col = getGID(interp_cell);
                                 
                                 MatSetValue(a_matrix, row, col, value, ADD_VALUES);
                              }
                           }
                           else {

                              // If we've landed here, then stencil_box_iv is neither a valid cell nor an extrablock
                              // ghost cell (relative to the center cell iv).  The only remaining possibility is that
                              // it's a ghost cell at a codim = 1 or codim = 2 physical boundary, in which case 
                              // the stencil modification to account for boundary conditions should have zeroed
                              // out the corresponding stencil entry (i.e., s = 0).  Something has therefore
                              // gone wrong and needs to be investigated.

                              MayDay::Error("MBPETScSolver::addUnstructuredMatrixEntries(): Invalid stencil_box_iv");
                           }
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
         const BaseFab<int>& this_gids = m_gids[dit];

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

                     double alpha = alphaPtr? alphaPtr->operator()(iv): 1.;

                     SparseCoupling& coupling = otherFab(iv_face,0);
                     for (SparseCouplingIterator it(coupling); it.ok(); ++it) {
                        IntVect index = coupling[it()];
                        double value = coupling.weight(index) * alpha;

                        int row = this_gids(iv);
                        int col = getGID(index);

                        MatSetValue(a_matrix, row, col, value, ADD_VALUES);
                     }  // it

                  }  // bit
               }  // if interface
            }  // dir
         }  // sit
      }  // dit
   }
}



void
MBPETScSolver::constructPETScMatrix( LevelData<FArrayBox>&              a_alpha_coefficient, 
                                     LevelData<FluxBox>&                a_tensor_coefficient,
                                     LevelData<FArrayBox>&              a_beta_coefficient,
                                     const PotentialBC&                 a_bc,
                                     FArrayBox&                         a_stencil_values,
                                     const int                          a_diagonal_offset,
                                     LayoutData< BaseFab<IntVectSet> >& a_unstructured_coupling,
                                     const bool                         a_fourth_order,
                                     const bool                         a_dropOrder,
                                     Mat&                               a_matrix,
                                     LevelData<FArrayBox>&              a_rhs_from_bc ) const
{
   const DisjointBoxLayout & grids = a_tensor_coefficient.disjointBoxLayout();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys_ptr->boundaries();

   MatZeroEntries(a_matrix);

   // Make sure the coefficient ghost cell values are set.
   // The ghost cell values at physical boundaries aren't used.
   a_tensor_coefficient.exchange();

   /* Set the matrix coefficients.  Each processor assigns coefficients
      for the boxes in the grid that it owns. Note that the coefficients
      associated with each stencil entry may vary from grid point to grid
      point if desired.  Here, we first set the same stencil entries for
      each grid point.  Then we make modifications to grid points near
      the boundary. */

   int stencil_size = a_stencil_values.box().numPts();

   double * stencil_values = a_stencil_values.dataPtr(0);

   CH_assert((a_stencil_values.box().size(0)-1)%2 == 0);
   int radius = (a_stencil_values.box().size(0)-1)/2;

   LevelData<FArrayBox> structured_values(grids, stencil_size, IntVect::Unit);
   FArrayBox tmp_stencil_values(a_stencil_values.box(), a_stencil_values.nComp());

   Vector< Vector<CoDim1Stencil> > codim1_stencils;
   Vector< Vector<CoDim2Stencil> > codim2_stencils;

   constructBoundaryStencils(a_fourth_order, a_bc, codim1_stencils, codim2_stencils );

   const Vector<Box>& mapping_blocks = m_coord_sys_ptr->mappingBlocks();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     const Box & box = grids[dit];

     int block_number = m_coord_sys_ptr->whichBlock(box);
     const NewCoordSys& block_coord_sys = *(m_coord_sys_ptr->getCoordSys(block_number));
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

     const Box& domain_box = mapping_blocks[block_number];

     FluxBox& this_coefs = a_tensor_coefficient[dit];

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

       // Compute the stencil couplings to cell iv without regard to boundary conditions, skipping
       // the contributions due to cell faces on block interfaces, which are handled later by
       // addUnstructuredMatrixEntries().

       for (int dir=0; dir<SpaceDim; ++dir) {
          bool at_lo_interface = (lo_block_interface[dir] && iv[dir] == domain_box.smallEnd(dir));
          bool at_hi_interface = (hi_block_interface[dir] && iv[dir] == domain_box.bigEnd(dir));
          bool force_codim2_condense = lo_block_interface[dir] || hi_block_interface[dir];

          for (SideIterator sit; sit.ok(); ++sit) {
             Side::LoHiSide side = sit();

             bool update_rhs_from_bc_only = (at_lo_interface && side == Side::LoHiSide::Lo)
                || (at_hi_interface && side == Side::LoHiSide::Hi);

             for (int dir2=0; dir2<SpaceDim; ++dir2) {

                tmp_stencil_values.setVal(0.);

                accumStencilMatrixEntries(iv, dir, side, dir2, this_coefs, dx,
                                          a_fourth_order, a_dropOrder, tmp_stencil_values);

                modifyStencilForBCs( codim1_stencils[block_number], codim2_stencils[block_number],
                                     iv, tmp_stencil_values, a_rhs_from_bc[dit],
                                     update_rhs_from_bc_only, force_codim2_condense );

                a_stencil_values += tmp_stencil_values;
             }
          }
       }

       if ( alpha ) {
          a_stencil_values *= alpha->operator()(iv);
       }

       if ( beta ) {
          stencil_values[a_diagonal_offset] += beta->operator()(iv);
       }

       for (int jj=0; jj<stencil_size; jj++) {
          structured_values[dit](iv,jj) = stencil_values[jj];
       }
     }
   }

   structured_values.exchange();

   Box tbox = a_stencil_values.box();
   IntVect bshift = radius*IntVect::Unit;
   tbox.shift(-bshift);

   int* indxn = new int[stencil_size];
   double* values = new double[stencil_size];

   IntVect periodic = checkPeriodic();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     const Box & box = grids[dit];
     const BaseFab<int>& this_gids = m_gids[dit];
     const FArrayBox& this_structured_values = structured_values[dit];

     int block_number = m_coord_sys_ptr->whichBlock(box);

     Box domain_box = mapping_blocks[block_number];

     for (int dir=0; dir<SpaceDim; ++dir) {
        if ( periodic[dir] != 0 ) {
           domain_box.grow(dir,radius);
        }
     }

     // Set the intra-block matrix entries

     BoxIterator bit(box);
     for (bit.begin();bit.ok();++bit) {
       IntVect iv = bit();

       for (int jj=0; jj<stencil_size; jj++) {
          stencil_values[jj] = this_structured_values(iv,jj);
       }

       Box loop_box = tbox + iv;
       Box valid_box = loop_box & domain_box;

       // Setting the matrix entries in batches is more efficient than
       // individually, so we load them into temporaries and make a
       // single call.
       int n = valid_box.numPts();
       int indxm = this_gids(iv);
       int k=0;
       int m=0;
       for (BoxIterator lit(loop_box); lit.ok(); ++lit,++m) {
          if ( valid_box.contains(lit()) ) {
             indxn[k] = this_gids(lit());
             values[k] = stencil_values[m];
             k++;
          }
       }

       MatSetValues(m_A, 1, &indxm, n, indxn, values, ADD_VALUES);
     }
   }

   delete [] values;
   delete [] indxn;

   addUnstructuredMatrixEntries(a_alpha_coefficient, a_tensor_coefficient, a_bc, a_stencil_values,
                                a_fourth_order, a_dropOrder, a_unstructured_coupling,
                                codim1_stencils, codim2_stencils, a_matrix);

   MatAssemblyBegin(a_matrix, MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(a_matrix, MAT_FINAL_ASSEMBLY);

   if (m_amg_solver_allocated) {
      KSPSetUp(m_amg_ksp);
   }
   if (m_gmres_solver_allocated) {
      KSPSetUp(m_gmres_ksp);
   }

#if 0

   testMatrixConstruct( a_alpha_coefficient, 
                        a_tensor_coefficient,
                        a_beta_coefficient,
                        block_boundaries,
                        codim1_stencils,
                        codim2_stencils,
                        a_stencil_values,
                        a_fourth_order,
                        a_dropOrder,
                        a_rhs_from_bc );
#endif
}



void
MBPETScSolver::AMG( const Mat&   a_matrix,
                    const Vec&   a_b,
                    const double a_tol,
                    const int    a_max_iter,
                    const bool   a_verbose,
                    const Vec&   a_x ) const
{
   KSPSetTolerances(m_amg_ksp, a_tol, 0., PETSC_DEFAULT, a_max_iter);

   KSPSolve(m_amg_ksp, a_b, a_x);

   double final_res_norm;
   KSPGetResidualNorm(m_amg_ksp, &final_res_norm);

   int num_iterations;
   KSPGetIterationNumber(m_amg_ksp, &num_iterations);

   if (a_verbose && procID()==0) {
      cout << "      AMG solver residual = " << final_res_norm << " after " << num_iterations << " iterations" << endl;
   }
}



void
MBPETScSolver::AMG_Preconditioned_GMRES( const Mat&    a_matrix,
                                         const Vec&    a_b,
                                         const double  a_tol,
                                         const int     a_max_iter,
                                         const double  a_amg_tol,
                                         const int     a_amg_max_iter,
                                         const bool    a_verbose,
                                         const Vec&    a_x ) const
{
   KSPSetTolerances(m_gmres_ksp, a_tol, 0., PETSC_DEFAULT, a_max_iter);

   KSPSolve(m_gmres_ksp, a_b, a_x);

   double final_res_norm;
   KSPGetResidualNorm(m_gmres_ksp, &final_res_norm);

   int num_iterations;
   KSPGetIterationNumber(m_gmres_ksp, &num_iterations);

   if (a_verbose && procID()==0) {
      cout << "      GMRES solver residual = " << final_res_norm << " after " << num_iterations << " iterations" << endl;
   }
}



IntVect
MBPETScSolver::checkPeriodic() const
{
   // Check for periodicity, which is only applicable for single block geometries, since there
   // doesn't seem to be any support for individually periodic blocks in MultiBlockLevelGeom.

   IntVect periodic;

   if (m_coord_sys_ptr->numBlocks() == 1) {
      const ProblemDomain& domain = m_geometry.gridsFull().physDomain();

      for (int dir=0; dir<SpaceDim; ++dir) {
         periodic[dir] = domain.isPeriodic(dir)? 1: 0;
      }
   }
   else {
      for (int dir=0; dir<SpaceDim; ++dir) {
         periodic[dir] = 0;
      }
   }

   return periodic;
}



void
MBPETScSolver::copyToPETScVector( const LevelData<FArrayBox>& a_in,
                                  Vec&                        a_out ) const
{
   const DisjointBoxLayout & grids = a_in.disjointBoxLayout();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      const BaseFab<int>& this_gids = m_gids[dit];
      const FArrayBox& values = a_in[dit];

      for (BoxIterator bit(box); bit.ok(); ++bit) {
         IntVect iv = bit();
         VecSetValues(a_out, 1, &this_gids(iv), &values(iv), INSERT_VALUES);
      }
   }

   VecAssemblyBegin(a_out);
   VecAssemblyEnd(a_out);
} 



void
MBPETScSolver::copyFromPETScVector( const Vec&            a_in,
                                    LevelData<FArrayBox>& a_out ) const
{
   const DisjointBoxLayout & grids = a_out.disjointBoxLayout();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      const BaseFab<int>& this_gids = m_gids[dit];
      FArrayBox& values = a_out[dit];

      for (BoxIterator bit(box); bit.ok(); ++bit) {
         IntVect iv = bit();
         VecGetValues(a_in, 1, &this_gids(iv), &values(iv));
      }
   }
} 



#include "NamespaceFooter.H"
#endif
