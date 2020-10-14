#include "GKVlasovAMG.H"

#include "Directions.H"
#include "MBStencilIterator.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "MagGeom.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"


GKVlasovAMG::GKVlasovAMG( const ParmParse&  a_pp,
                          const PhaseGeom&  a_geometry )
   : m_geometry(a_geometry),
     m_block_register(NULL),
     m_mb_coupling(NULL),
     m_hypre_allocated(false),
     m_A(NULL),
     m_AMG_solver_allocated(false)
{
   if ( a_pp.contains("restrict_type") ) {
      a_pp.get("restrict_type", m_restrict_type);
   }
   else {
      m_restrict_type = 2;
   }

   double tol;
   if ( a_pp.contains("tol") ) {
      a_pp.get("tol", tol);
   }
   else {
      tol = 0.;
   }

   int max_iter;
   if ( a_pp.contains("max_iter") ) {
      a_pp.get("max_iter", max_iter);
   }
   else {
      max_iter = 2;
   }

   bool verbose;
   if ( a_pp.contains("verbose") ) {
      a_pp.get("verbose", verbose);
   }
   else {
      verbose = false;
   }

   setConvergenceParams(tol, max_iter, verbose);

   for (int dir=0; dir<CFG_DIM; ++dir) {
      m_flux_average[dir] = true;
   }

#if CFG_DIM==3
   if ( ((CFG::MagGeom&)m_geometry.magGeom()).shearedMBGeom()) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         m_flux_average[dir] = false;
      }
   }
#endif

   // If we are doing flux averaging in any configuration space direction, then create
   // a BlockRegister.
   for (int dir=0; dir<CFG_DIM; ++dir) {
      if (m_flux_average[dir]) {
         m_block_register = new BlockBaseRegister<BaseFab<SparseCoupling> >(m_geometry.coordSysPtr(),
                                                                            m_geometry.gridsFull(), 0);
         break;
      }
   }

   // Save cell volumes
   m_volume.define(m_geometry.gridsFull(), 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume);
}

GKVlasovAMG::~GKVlasovAMG()
{
   if (m_mb_coupling) delete m_mb_coupling;

   for (int n=0; n<m_stencil_size; ++n) delete [] m_offsets[n];
   delete m_offsets;

   if (m_block_register) delete m_block_register;
}


void GKVlasovAMG::createHypreData()
{
   CH_TIMERS("GKVlasovAMG::createHypreData()");
   CH_TIMER("get_unstructured_graph",t_get_unstructured_graph);
   CH_TIMER("add_unstructured_graph",t_add_unstructured_graph);
   if (m_hypre_allocated) {
      MayDay::Error( "GKVlasovAMG:createHypreData:: data already allocated; deallocate first!" );
   }

   m_hypre_object_type = HYPRE_PARCSR;

   const RefCountedPtr<PhaseCoordSys>& coord_sys_ptr = m_geometry.phaseCoordSysPtr();

   int num_blocks = coord_sys_ptr->numBlocks();

   // Create the grid

   {
      HYPRE_SStructGridCreate(MPI_COMM_WORLD, SpaceDim, num_blocks, &m_grid);

      const DisjointBoxLayout& grids = m_geometry.gridsFull();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& box = grids[dit];
         int block_number = coord_sys_ptr->whichBlock(box);
         IntVect lower(box.loVect());
         IntVect upper(box.hiVect());

         HYPRE_SStructGridSetExtents(m_grid, block_number, lower.dataPtr(), upper.dataPtr());
      }

      HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};

      for (int block=0; block<num_blocks; ++block) {
         HYPRE_SStructGridSetVariables(m_grid, block, 1, vartypes);
      }

      // Specify periodicity.  This is only applicable for single block geometry,
      // since there doesn't seem to be any support for individually periodic blocks
      // in MultiBlockLevelGeom.

      if (num_blocks == 1) {

         const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = coord_sys_ptr->boundaries();

         const ProblemDomain& domain = grids.physDomain();
         IntVect periodic = domain.size();
         int block_number = 0;

         for (int dir=0; dir<SpaceDim; dir++) {
            bool this_dir_periodic = domain.isPeriodic(dir);

            if ( !this_dir_periodic ||
                 (block_boundaries[block_number][dir].isInterface() && m_mb_coupling != NULL) ||
                 (block_boundaries[block_number][dir + SpaceDim].isInterface() && m_mb_coupling != NULL) ) {
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
                     MayDay::Error( "GKVlasovAMG::createHypreData(): Periodic direction must be a power of 2!" );
                  }
               }
            }
         }

         HYPRE_SStructGridSetPeriodic(m_grid, block_number, periodic.dataPtr());
      }

      int num_ghost[CH_SPACEDIM*2];
      for (int n=0; n<2*(CFG_DIM+1); ++n) {
         num_ghost[n] = m_stencil_radius;
      }
      num_ghost[2*(CFG_DIM+1)] = 0;
      num_ghost[2*(CFG_DIM+1)+1] = 0;

      HYPRE_SStructGridSetNumGhost(m_grid, num_ghost);

      HYPRE_SStructGridAssemble(m_grid);
   }

   // Define the discretization stencil for the system matrix A

   {
      /* Define the geometry of the stencil. Each represents a
         relative offset (in the index space). */

      m_stencil_size = 2*(CFG_DIM+1)*m_stencil_radius + 1;
      m_offsets = new int*[m_stencil_size];

      // Stencil center offset
      int n = 0;
      m_offsets[n] = new int[CH_SPACEDIM];
      for (int dir=0; dir<CH_SPACEDIM; ++dir) {
         m_offsets[n][dir] = 0;
      }

      for (int dir=0; dir<CFG_DIM+1; ++dir) {

         for (int m=1; m<=m_stencil_radius; ++m) {
            // Stencil offsets along postive leg
            n++;
            m_offsets[n] = new int[CH_SPACEDIM];
            for (int tdir=0; tdir<CH_SPACEDIM; ++tdir) {
               if (tdir == dir) {
                  m_offsets[n][tdir] = m;
               }
               else {
                  m_offsets[n][tdir] = 0;
               }
            }
            // Stencil offsets along negative leg
            n++;
            m_offsets[n] = new int[CH_SPACEDIM];
            for (int tdir=0; tdir<CH_SPACEDIM; ++tdir) {
               if (tdir == dir) {
                  m_offsets[n][tdir] = -m;
               }
               else {
                  m_offsets[n][tdir] = 0;
               }
            }
         }
                                             
      }

      /* Create an empty stencil object */
      HYPRE_SStructStencilCreate(SpaceDim, m_stencil_size, &m_A_stencil);

      /* Assign each of the stencil entries */
      for (int entry = 0; entry < m_stencil_size; entry++)
         HYPRE_SStructStencilSetEntry(m_A_stencil, entry, m_offsets[entry], 0);
   }

   // Set up the graph

   {
      HYPRE_SStructGraphCreate(MPI_COMM_WORLD, m_grid, &m_A_graph);

      HYPRE_SStructGraphSetObjectType(m_A_graph, m_hypre_object_type);

      for (int block=0; block<num_blocks; ++block) {
         HYPRE_SStructGraphSetStencil(m_A_graph, block, 0, m_A_stencil);
      }

      CH_START(t_get_unstructured_graph);
      getUnstructuredCouplings(m_stencil_radius, m_A_unstructured_coupling);
      CH_STOP(t_get_unstructured_graph);
      CH_START(t_add_unstructured_graph);
      addUnstructuredGraphEntries(m_stencil_radius, m_A_unstructured_coupling, m_A_graph);
      CH_STOP(t_add_unstructured_graph);

      HYPRE_SStructGraphAssemble(m_A_graph);
   }

   // Set up b and x vectors

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


void GKVlasovAMG::AMGSetup( const HYPRE_SStructMatrix& a_matrix )
{
   CH_TIME("GKVlasovAMG::AMGSetup");

   if ( m_AMG_solver_allocated ) {
      HYPRE_BoomerAMGDestroy(m_par_AMG_solver);
   }
   
   HYPRE_BoomerAMGCreate(&m_par_AMG_solver);

   /*
     Instructions from Ruipeng Li (11/6/2019):

     To use the AIR (approximate ideal restriction) algorithm in hypre, you need to set
     HYPRE_BoomerAMGSetRestriction(amg_solver, restri_type) with restri_type = 1, or 2 for distance-1
     and -2 AIR respectively. Distance-2 is more expensive to build and leads to a denser restriction
     operator than the distance-1 but should approximate the ideal restriction better. Do not use aggressive
     coarsening with AIR ( i.e., HYPRE_ParCSRHybridSetAggNumLevels(amg, 0) ). 

     For hyperbolic problems like advection dominant, highly nonsymmetric, using very simple interpolation
     like 1-pt interpolation was often found enough, set HYPRE_BoomerAMGSetInterpType(amg, 100) for that.
     Also,  C/F relaxations can be helpful, in the paper, we used FFC, which can be set by
     HYPRE_BoomerAMGSetGridRelaxPoints(amg_solver, grid_relax_points); with          
         grid_relax_points[2][0] = -1; // F
         grid_relax_points[2][1] = -1; // F
         grid_relax_points[2][2] =  1; // C

     These options can be found in hypre, src/test/ij.c, search for air or AIR to see the default settings.
   */

   HYPRE_BoomerAMGSetRestriction(m_par_AMG_solver, m_restrict_type);

   HYPRE_Int **grid_relax_points = NULL;
   HYPRE_Int ns_down = 1;
   HYPRE_Int ns_up = 3;
   HYPRE_Int ns_coarse = 1;
   /* this is a 2-D 4-by-k array using Double pointers */
   grid_relax_points = hypre_CTAlloc(HYPRE_Int*, 4, HYPRE_MEMORY_HOST);
   grid_relax_points[0] = NULL;
   grid_relax_points[1] = hypre_CTAlloc(HYPRE_Int, ns_down, HYPRE_MEMORY_HOST);
   grid_relax_points[2] = hypre_CTAlloc(HYPRE_Int, ns_up, HYPRE_MEMORY_HOST);
   grid_relax_points[3] = hypre_CTAlloc(HYPRE_Int, ns_coarse, HYPRE_MEMORY_HOST);
   /* down cycle: C */
   for (int i=0; i<ns_down; i++) {
      grid_relax_points[1][i] = 0;//1;
   }
   /* up cycle: F */
   //for (int i=0; i<ns_up; i++)
   //{
   if (ns_up == 3) {
      grid_relax_points[2][0] = -1; // F
      grid_relax_points[2][1] = -1; // F
      grid_relax_points[2][2] =  1; // C
   }
   else if (ns_up == 2) {
      grid_relax_points[2][0] = -1; // F
      grid_relax_points[2][1] = -1; // F
   }
   //}
   /* coarse: all */
   for (int i=0; i<ns_coarse; i++) {
      grid_relax_points[3][i] = 0;
   }
   HYPRE_BoomerAMGSetGridRelaxPoints(m_par_AMG_solver, grid_relax_points);

   HYPRE_BoomerAMGSetRelaxType(m_par_AMG_solver, 0);  // Jacobi

   HYPRE_Int interp_type = 100;  // 1-pt
   HYPRE_BoomerAMGSetInterpType(m_par_AMG_solver, interp_type);
   
   HYPRE_BoomerAMGSetAggNumLevels(m_par_AMG_solver, 0);

   HYPRE_BoomerAMGSetMaxCoarseSize(m_par_AMG_solver, 0);

   //   HYPRE_BoomerAMGSetStrongThreshold(m_par_AMG_solver, 0.25);
   //   HYPRE_BoomerAMGSetStrongThresholdR(m_par_AMG_solver, 0.25);
   //   HYPRE_BoomerAMGSetCoarsenType(m_par_AMG_solver, 6);  // Falgout coarsening
   //   HYPRE_BoomerAMGSetCycleRelaxType(m_par_AMG_solver, 3, 3); // hybrid Gauss-Seidel on coarsest level
   //   HYPRE_BoomerAMGSetCycleRelaxType(m_par_AMG_solver, 9, 3); // Gaussian elimination on coarsest level
   //   HYPRE_BoomerAMGSetCycleType(m_par_AMG_solver, 1);  // V = 1, W = 2
   //   HYPRE_BoomerAMGSetCycleNumSweeps(m_par_AMG_solver, 1, 1);  // 1 sweep on down cycle
   //   HYPRE_BoomerAMGSetCycleNumSweeps(m_par_AMG_solver, 1, 2);  // 1 sweep on up cycle
   //   HYPRE_BoomerAMGSetCycleNumSweeps(m_par_AMG_solver, 1, 3);  // 1 sweeps on coarsest level

   m_AMG_solver_allocated = true;

   HYPRE_ParCSRMatrix par_A;
   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);

   HYPRE_ParVector par_b;  // not used, but needs to be passed to HYPRE_BoomerAMGSetup
   HYPRE_ParVector par_x;  // not used, but needs to be passed to HYPRE_BoomerAMGSetup

   HYPRE_BoomerAMGSetup(m_par_AMG_solver, par_A, par_b, par_x);
}


void GKVlasovAMG::constructMatrix( const string&                                 a_face_avg_type,
                                   const double                                  a_shift,
                                   const LevelData<BaseFab<Vector<IntVect> > >&  a_structured_couplings,
                                   const LevelData<BaseFab<Vector<Real> > >&     a_structured_weights )
{
   CH_TIMERS("GKVlasovAMG::constructMatrix()");
   CH_TIMER("create_hypre_data",t_create_hypre_data);
   CH_TIMER("add_struct_matrix",t_add_struct_matrix);
   CH_TIMER("add_unstruct_matrix",t_add_unstruct_matrix);
   CH_TIMER("amg_setup",t_amg_setup);
   CH_TIMER("hypre_initialize",t_hypre_initialize);
   CH_TIMER("hypre_assemble",t_hypre_assemble);
   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   const RefCountedPtr<PhaseCoordSys>& coord_sys_ptr = m_geometry.phaseCoordSysPtr();

   if ( !m_hypre_allocated ) {

      if ( a_face_avg_type == "uw1" ) {
         m_stencil_radius = 1;
      }
      else if ( a_face_avg_type == "uw3" ) {
         m_stencil_radius = 2;
      }
      else {
         MayDay::Error("GKVlasovAMG::constructMatrix(): Unrecognized face_avg_type");
      }

      if ( coord_sys_ptr->numBlocks() > 1 || ((CFG::MagGeom&)m_geometry.magGeom()).shearedMBGeom()) {
         int ghosts = m_stencil_radius;
         int order = 2;
         m_mb_coupling = new GKVlasovMBCoupling(m_geometry, ghosts, order);
      }
      
      CH_START(t_create_hypre_data);
      createHypreData();
      CH_STOP(t_create_hypre_data);
   }

   CH_START(t_hypre_initialize);
   HYPRE_SStructMatrixDestroy(m_A);
   HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, m_A_graph, &m_A);
   HYPRE_SStructMatrixSetObjectType(m_A, m_hypre_object_type);

   HYPRE_SStructMatrixInitialize(m_A);
   CH_STOP(t_hypre_initialize);

   bool lo_block_interface[CFG_DIM+1];
   bool hi_block_interface[CFG_DIM+1];
   bool lo_block_boundary[CFG_DIM+1];
   bool hi_block_boundary[CFG_DIM+1];

   CH_START(t_add_struct_matrix);
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = coord_sys_ptr->boundaries();
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box& box = grids[dit];
      int block_number = coord_sys_ptr->whichBlock(box);
      const Box& block_box = (coord_sys_ptr->mappingBlocks())[block_number];

      const FArrayBox& this_volume = m_volume[dit];
      const PhaseBlockCoordSys* block_coord_sys = coord_sys_ptr->getCoordSys(block_number);
      double mapped_volume = block_coord_sys->getMappedCellVolume();

      const BaseFab<Vector<IntVect> >& this_couplings = a_structured_couplings[dit];
      const BaseFab<Vector<Real> >& this_weights = a_structured_weights[dit];

      const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
      for (int dir=0; dir<CFG_DIM+1; ++dir) {
         lo_block_interface[dir] = this_block_boundaries[dir].isInterface();
         hi_block_interface[dir] = this_block_boundaries[dir + SpaceDim].isInterface();
         lo_block_boundary[dir] = this_block_boundaries[dir].isDomainBoundary();
         hi_block_boundary[dir] = this_block_boundaries[dir + SpaceDim].isDomainBoundary();
      }

      // Set the intra-block matrix entries

      for (BoxIterator bit(box); bit.ok(); ++bit) {
         IntVect iv = bit();

         // Initialize the matrix with a_shift times the identity.  The J (= physical cell volume
         // divided by mapped cell volume) factor is needed because the matrix being constructed
         // here is assumed to act on vectors in the physical frame, but the time-integration
         // shift is done in the computational frame
         int center_entry = findStructuredEntry(iv, iv);
         double center_value = a_shift * this_volume(iv,0) / mapped_volume;
         HYPRE_SStructMatrixSetValues(m_A, block_number, iv.dataPtr(), 0, 1, &center_entry, 
                                      &center_value);

         for (int dir=0; dir<CFG_DIM+1; ++dir) {
            bool at_lo_interface = lo_block_interface[dir] && (iv[dir] == block_box.smallEnd(dir));
            bool at_hi_interface = hi_block_interface[dir] && (iv[dir] == block_box.bigEnd(dir));

            for (SideIterator sit; sit.ok(); ++sit) {
               Side::LoHiSide side = sit();

               // Couplings corresponding to faces at block interfaces are handled later
               if ( m_mb_coupling == NULL || 
                    (!(side == Side::LoHiSide::Lo && at_lo_interface) &&
                     !(side == Side::LoHiSide::Hi && at_hi_interface) ) ) {

                  int coupling_comp = dir + ((sign(side)+1)/2)*(CFG_DIM+1);
                  for (int n=0; n<this_couplings(iv,coupling_comp).size(); ++n) {
                     IntVect coupling = this_couplings(iv,coupling_comp)[n];
                     double value = this_weights(iv,coupling_comp)[n];

                     // Zero the couplings across physical boundaries
                     if ( (coupling[dir] < block_box.smallEnd(dir) && lo_block_boundary[dir]) ||
                          (coupling[dir] > block_box.bigEnd(dir)   && hi_block_boundary[dir]) ) {
                        value = 0.;
                     }

                     int entry = findStructuredEntry(iv, coupling);
                     HYPRE_SStructMatrixAddToValues(m_A, block_number, iv.dataPtr(), 0, 
                                                    1, &entry, &value);
                  }
               }
            }
         }
      }
   }
   CH_STOP(t_add_struct_matrix);

   CH_START(t_add_unstruct_matrix);
   addUnstructuredMatrixEntries(a_structured_couplings, a_structured_weights,
                                m_A_unstructured_coupling, m_A);
   CH_STOP(t_add_unstruct_matrix);
   
   // Finalize the matrix assembly
   CH_START(t_hypre_assemble);
   HYPRE_SStructMatrixAssemble(m_A);
   CH_STOP(t_hypre_assemble);

   //   HYPRE_SStructMatrixPrint("HYPRE_MATRIX.A", m_A, 0);

   CH_START(t_amg_setup);
   AMGSetup(m_A);
   CH_STOP(t_amg_setup);
}


void
GKVlasovAMG::addUnstructuredMatrixEntries(
               const LevelData<BaseFab<Vector<IntVect> > >&  a_structured_couplings,
               const LevelData<BaseFab<Vector<Real> > >&     a_structured_weights,
               const LayoutData< BaseFab<IntVectSet> >&      a_unstructured_couplings,
               HYPRE_SStructMatrix&                          a_matrix ) const
{
   CH_TIMERS("GKVlasovAMG::addUnstructuredMatrixEntries()");
   CH_TIMER("add_stencil_contrib",t_add_stencil_contrib);
   CH_TIMER("add_extrablock_contrib",t_add_extrablock_contrib);

   if ( m_mb_coupling ) {

      const RefCountedPtr<PhaseCoordSys>& coord_sys_ptr = m_geometry.phaseCoordSysPtr();
      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = coord_sys_ptr->boundaries();
      const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencil = m_mb_coupling->stencils();

      const DisjointBoxLayout & grids = m_geometry.gridsFull();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box & box = grids[dit];

         int block_number = coord_sys_ptr->whichBlock(box);

         const NewCoordSys* block_coord_sys = coord_sys_ptr->getCoordSys(block_number);
         RealVect dx = block_coord_sys->dx();

         const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
         bool lo_block_interface[CFG_DIM+1];
         bool hi_block_interface[CFG_DIM+1];
         for (int dir=0; dir<CFG_DIM; ++dir) {
            lo_block_interface[dir] = this_block_boundaries[dir].isInterface();
            hi_block_interface[dir] = this_block_boundaries[dir + SpaceDim].isInterface();
         }

         const BaseFab<Vector<IntVect> >& this_structured_couplings = a_structured_couplings[dit];
         const BaseFab<Vector<Real> >& this_structured_weights = a_structured_weights[dit];

         const Box& domainBox = (coord_sys_ptr->mappingBlocks())[block_number];
         if ( stencil[dit] != NULL ) {
            const IVSFAB<MBStencil>& this_stencil = *stencil[dit];

            BaseFab<SparseCoupling> data[CFG_DIM+1];
            for(int dir=0; dir<CFG_DIM; ++dir ) {
               if (m_flux_average[dir]) {
                  data[dir].define(surroundingNodes(box, dir), 1);
               }
            }

            // Here, we loop over all cells whose stencil is affected by a block boundary
            IntVectSet interblock_ivs = getInterBlockCoupledCells(block_number, m_stencil_radius, box);
            for (IVSIterator ivsit(interblock_ivs); ivsit.ok(); ++ivsit) {
               IntVect iv = ivsit();

               // Get all of the indices that are sparsely coupled to this cell
               const IntVectSet& unstructured_couplings = a_unstructured_couplings[dit](iv);

               for (int dir=0; dir<CFG_DIM; ++dir) {
                  bool at_lo_interface = lo_block_interface[dir] && (iv[dir] == domainBox.smallEnd(dir));
                  bool at_hi_interface = hi_block_interface[dir] && (iv[dir] == domainBox.bigEnd(dir));
               
                  Box stencil_box(iv,iv);
                  stencil_box.grow(dir, m_stencil_radius);
                  Box stencil_box_valid = stencil_box & domainBox;

                  for (SideIterator sit; sit.ok(); ++sit) {
                     Side::LoHiSide side = sit();

                     IntVect iv_face;

                     // Halve the contributions on face contained in block interfaces to effect averaging
                     double fac;
                     if ( (at_lo_interface && side == Side::LoHiSide::Lo) ||
                          (at_hi_interface && side == Side::LoHiSide::Hi) ) {
                        fac = m_flux_average[dir]? 0.5: 1.0;
                        iv_face = iv;
                        iv_face.shift(dir,side);
                     }
                     else {
                        fac = 1.;
                     }

                     int coupling_comp = dir + ((sign(side)+1)/2)*(CFG_DIM+1);
                     for (int m=0; m<this_structured_couplings(iv,coupling_comp).size(); ++m) {
                        IntVect stencil_box_iv = this_structured_couplings(iv,coupling_comp)[m];

                        double s = fac * this_structured_weights(iv,coupling_comp)[m];

                        if ( s != 0. ) {

                           if ( stencil_box_valid.contains(stencil_box_iv) ) {

                              CH_START(t_add_stencil_contrib);
                              // If the current side of cell iv is on a block interface, add the stencil
                              // contribution to the BlockRegister data.
                              // The minus sign accounts for the normal component, which is negated on
                              // the adacent cell.  Also, since we have not
                              // yet accounted for the regular stencil contributions, do that now.
                              if ( (at_lo_interface && side == Side::LoHiSide::Lo) || 
                                   (at_hi_interface && side == Side::LoHiSide::Hi) ) {
                                 if (m_flux_average[dir]) {
                                    data[dir](iv_face).add(stencil_box_iv, -s);
                                 }

                                 int entry = findEntry(iv, unstructured_couplings, stencil_box_iv);
                                 double value = s;

                                 int num_entries = 1;
                                 HYPRE_SStructMatrixAddToValues(a_matrix, block_number, iv.dataPtr(),
                                                                0, num_entries, &entry, &value);
                              }
                              CH_STOP(t_add_stencil_contrib);
                           }
                           else if ( (stencil_box_iv[dir] < domainBox.smallEnd(dir) 
                                      && lo_block_interface[dir]) ||
                                     (stencil_box_iv[dir] > domainBox.bigEnd(dir)
                                      && hi_block_interface[dir]) ) {

                              CH_START(t_add_extrablock_contrib);

                              const MBStencil& elements = this_stencil(stencil_box_iv, 0);
                              int num_unstructured_entries = elements.size();
                              double * values = new double[num_unstructured_entries];
                              int * entries = new int[num_unstructured_entries];

                              MBStencilIterator stit(elements);
                              int k = 0;
                              for (stit.begin(); stit.ok(); ++stit) {
                                 const MBStencilElement& stencilElement = stit();
                                 IntVect interp_cell = stencilElement.cell();

                                 entries[k] = findEntry(iv, unstructured_couplings, interp_cell);
                                 values[k] = s * stencilElement.weight();

                                 // If we are on a block interface, add the stencil contribution to the
                                 // BlockRegister data.  The minus sign accounts for the normal component,
                                 // which is negated on the adacent cell.
                                 if ( m_flux_average[dir] &&
                                      ((at_lo_interface && side == Side::LoHiSide::Lo) ||
                                       (at_hi_interface && side == Side::LoHiSide::Hi)) ) {
                                    data[dir](iv_face).add(interp_cell, -values[k]);
                                 }

                                 k++;
                              }

                              HYPRE_SStructMatrixAddToValues(a_matrix, block_number, iv.dataPtr(),
                                                             0, num_unstructured_entries, entries, values);

                              delete [] entries;
                              delete [] values;

                              CH_STOP(t_add_extrablock_contrib);
                           }
                        }
                     }
                  }
               }
            }

            // Store the contributions at block interfaces
            for (int dir=0; dir<CFG_DIM; dir++) {
               if ( m_flux_average[dir] ) {
                  for (SideIterator sit; sit.ok(); ++sit) {
                     Side::LoHiSide side = sit();
                     if (m_block_register->hasInterface(dit(), dir, side)) {
                         m_block_register->store(data[dir], dit(), dir, side);
                     }
                  }
               }
            }
         }
      }

      // Exchange the block interface contributions and accumulate them to the matrix

      if ( m_block_register ) {
         m_block_register->exchange();
      }

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         int block_number = coord_sys_ptr->whichBlock(grids[dit]);

         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
            for (int dir=0; dir<CFG_DIM; dir++) {
               if ( m_flux_average[dir] && m_block_register->hasInterface(dit(), dir, side) ) {
                  Box faceBox = adjCellBox(grids[dit], dir, side, 1);
                  // if Lo, then shift +1; if Hi, then shift -1
                  faceBox.shiftHalf(dir, -sign(side));
                  BaseFab<SparseCoupling> otherFab(faceBox, 1);
                  m_block_register->fill(otherFab, dit(), dir, side, side);
                  BoxIterator bit(faceBox);
                  for (bit.begin(); bit.ok(); ++bit) {
                     IntVect iv_face = bit();
                     IntVect iv(iv_face);
                     iv.shift(dir,-side);

                     const IntVectSet& unstructured_couplings = a_unstructured_couplings[dit](iv);

                     SparseCoupling& coupling = otherFab(iv_face,0);
                     int num_entries = coupling.size();
                     double * values = new double[num_entries];
                     int * entries = new int[num_entries];

                     int k = 0;
                     for (SparseCouplingIterator it(coupling); it.ok(); ++it) {
                        IntVect index = coupling[it()];
                        entries[k] = findEntry(iv, unstructured_couplings, index);
                        values[k] = coupling.weight(index);
                        k++;
                     }

                     HYPRE_SStructMatrixAddToValues(a_matrix, block_number, iv.dataPtr(),
                                                    0, num_entries, entries, values);

                     delete [] entries;
                     delete [] values;

                  }
               }
            }
         }
      }
   }
}


int GKVlasovAMG::findStructuredEntry( const IntVect&  a_iv_center,
                                      const IntVect&  a_iv ) const
{
   /*
     The Hypre matrix coefficients corresponding to the equation centered at cell a_iv are
     identified by an entry number, which is returned by this function.  Regular "stencil"
     couplings to cells are enumerated first, followed by the unstructured couplings to
     cells in other blocks contained in a_unstructured_ivs.
   */

   bool found_entry = false;

   int entry;
   for (entry=0; entry<m_stencil_size; ++entry) {
      IntVect offset(m_offsets[entry]);
      offset += a_iv_center;
      if (offset == a_iv) {
         found_entry = true;
         break;
      }
   }

   CH_assert(found_entry);

   return entry;
}


int GKVlasovAMG::findEntry( const IntVect&     a_iv_center,
                            const IntVectSet&  a_unstructured_ivs,
                            const IntVect&     a_iv ) const
{
   CH_TIME("GKVlasovAMG::findEntry()");
   /*
     The Hypre matrix coefficients corresponding to the equation centered at cell a_iv are
     identified by an entry number, which is returned by this function.  Regular stencil
     couplings are enumerated first, followed by the unstructured couplings to cells in
     other blocks contained in a_unstructured_ivs.
   */

   bool found_entry = false;

   int entry;
   for (entry=0; entry<m_stencil_size; ++entry) {
      IntVect offset(m_offsets[entry]);
      offset += a_iv_center;
      if (offset == a_iv) {
         found_entry = true;
         break;
      }
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


IntVectSet
GKVlasovAMG::getInterBlockCoupledCells( const int   a_block_number,
                                        const int   a_radius,
                                        const Box&  a_box ) const
{
   const RefCountedPtr<PhaseCoordSys>& coord_sys_ptr = m_geometry.phaseCoordSysPtr();

   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = coord_sys_ptr->boundaries();

   const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[a_block_number];
   const Box& domain_box = (coord_sys_ptr->mappingBlocks())[a_block_number];
   CH_assert(domain_box.contains(a_box));

   IntVectSet ivs;
   for (int dir=0; dir<CFG_DIM+1; ++dir) {
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


void
GKVlasovAMG::getUnstructuredCouplings( int                                 a_radius,
                                       LayoutData< BaseFab<IntVectSet> >&  a_unstructured_couplings ) const
{
   if ( m_mb_coupling ) {

      const RefCountedPtr<PhaseCoordSys>& coord_sys_ptr = m_geometry.phaseCoordSysPtr();

      LayoutData< Vector< BaseFab<IntVectSet>* > > all_couplings;

      //const LayoutData< IntVectSet >& ghostCells = m_mb_coupling->ghostCells();
      const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencil = m_mb_coupling->stencils();

      const DisjointBoxLayout & grids = m_geometry.gridsFull();
    
      all_couplings.define(grids);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& box = grids[dit];
         for (int dir=0; dir<CFG_DIM+1; ++dir) {
            Box dir_box = box;
            dir_box.surroundingNodes(dir);
            BaseFab<IntVectSet>* bi = new BaseFab<IntVectSet>(dir_box,1);
            all_couplings[dit].push_back(bi);
         }
      }
      
      const Vector<Box>& mapping_blocks = coord_sys_ptr->mappingBlocks();
      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = coord_sys_ptr->boundaries();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         //const IntVectSet& this_ghosts_fab = ghostCells[dit];
         if ( stencil[dit] != NULL  ) {
            const IVSFAB<MBStencil>& this_stencil_fab = *stencil[dit];
            int block_number = coord_sys_ptr->whichBlock(grids[dit]);
            const Box& domain_box = mapping_blocks[block_number];
            const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
      
            Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
            for (int dir=0; dir<CFG_DIM+1; ++dir) {
               bool lo_block_interface = false;
               bool hi_block_interface = false;
               if ( dir < CFG_DIM ) {
                  lo_block_interface = this_block_boundaries[dir].isInterface();
                  hi_block_interface = this_block_boundaries[dir + SpaceDim].isInterface();
               }

               BaseFab<IntVectSet>& this_coupling_dir = *this_coupling[dir];
               const Box& dir_box = this_coupling_dir.box();
        
               Box template_box(IntVect::Zero, IntVect::Zero);
               template_box.grow(dir,a_radius-1);
               template_box.growLo(dir, 1);
        
               BoxIterator bit(dir_box);
               for (bit.begin();bit.ok();++bit) {
                  IntVect iv_face = bit();
                  Box stencil_box = template_box + iv_face;

                  // Find the cells in the current stencil that are valid in the current block
                  IntVectSet valid_cells(domain_box & stencil_box);

                  // Find the cells in the current stencil corresponding to extra-block
                  // ghost cells filled by interpolation from other blocks
                  IntVectSet extra_block_ghosts;
                  for (BoxIterator bit2(stencil_box); bit2.ok(); ++bit2) {
                     IntVect iv2 = bit2();
                     if ( (iv2[dir] < domain_box.smallEnd(dir) && lo_block_interface) ||
                          (iv2[dir] > domain_box.bigEnd(dir) && hi_block_interface) ) {
                        extra_block_ghosts |= iv2;
                     }
                  }

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

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
         int dst_block_number = coord_sys_ptr->whichBlock(grids[dit]);

         // Find all cells that may be coupled through a block interface
         IntVectSet dst_ivs = getInterBlockCoupledCells(dst_block_number, a_radius, grids[dit]);

         for (int dir=0; dir<CFG_DIM; dir++) {
            if ( m_flux_average[dir] ) {
               BaseFab<IntVectSet>& this_coupling_dir = *this_coupling[dir];
               BaseFab<SparseCoupling> data(surroundingNodes(grids[dit], dir), 1);

               for (SideIterator sit; sit.ok(); ++sit) {
                  Side::LoHiSide side = sit();
                  if (m_block_register->hasInterface(dit(), dir, side)) {

                     for (IVSIterator it(dst_ivs); it.ok(); ++it) {
                        IntVect iv = it();
                        iv.shift(dir,side);

                        IntVectSet& unstructured_ivs = this_coupling_dir(iv);
                        for (IVSIterator sivsit(unstructured_ivs); sivsit.ok(); ++sivsit) {
                           data(iv).add( sivsit(), 0. );
                        }
                     }

                     m_block_register->store(data, dit(), dir, side);
                  }
               }
            }
         }
      }

      if ( m_block_register ) {
         m_block_register->exchange();
      }

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& baseBox = grids[dit];
         Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
            for (int dir=0; dir<CFG_DIM; dir++) {
               if ( m_flux_average[dir] ) {
                  BaseFab<IntVectSet>& this_coupling_dir = *this_coupling[dir];
                  if (m_block_register->hasInterface(dit(), dir, side)) {
                     Box faceBox = adjCellBox(baseBox, dir, side, 1);
                     // if Lo, then shift +1; if Hi, then shift -1
                     faceBox.shiftHalf(dir, -sign(side));
                     BaseFab<SparseCoupling> otherFab(faceBox, 1);
                     m_block_register->fill(otherFab, dit(), dir, side, side);
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

      const Box stencil_offsets(-a_radius*IntVect::Unit, a_radius*IntVect::Unit);

      a_unstructured_couplings.define(grids);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Vector< BaseFab<IntVectSet>* >& this_coupling = all_couplings[dit];
         int dst_block_number = coord_sys_ptr->whichBlock(grids[dit]);
         Box domain_box = mapping_blocks[dst_block_number];
         BaseFab<IntVectSet>& this_unstructured_coupling = a_unstructured_couplings[dit];
         
         this_unstructured_coupling.resize(grids[dit], 1);

         // Find all cells that may be coupled through a block interface
         IntVectSet dst_ivs = getInterBlockCoupledCells(dst_block_number, a_radius, grids[dit]);

         for (IVSIterator it(dst_ivs); it.ok(); ++it) {
            IntVect iv = it();
            
            IntVectSet unstructured_ivs;
            for (int dir=0; dir<CFG_DIM+1; ++dir) {

               Box stencil_offsets_dir(IntVect::Zero,IntVect::Zero);
               stencil_offsets_dir.grow(dir, a_radius);
               Box stencil_box = stencil_offsets_dir + iv;
               Box stencil_box_valid = stencil_box & domain_box;

               unstructured_ivs |= 
                  getUnstructuredCouplingsForCell( this_coupling, stencil_box_valid, dir, iv );

            }
            this_unstructured_coupling(iv) = unstructured_ivs;
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
GKVlasovAMG::getUnstructuredCouplingsForCell( const Vector< BaseFab<IntVectSet>* >&  a_coupling,
                                              const Box&                             a_stencil_box,
                                              const int                              a_dir,
                                              const IntVect&                         a_iv ) const
{
   IntVectSet unstructured_ivs;

   const BaseFab<IntVectSet>& this_coupling = *a_coupling[a_dir];
   for (int side=0; side<2; ++side) {
      IntVect index(a_iv);
      index.shift(a_dir,side);

      const IntVectSet& ivs = this_coupling(index);
                       
      for (IVSIterator ivit(ivs); ivit.ok(); ++ivit) {
         IntVect iv2 = ivit();
         if ( !a_stencil_box.contains(iv2) ) {
            unstructured_ivs |= iv2;
         }
      }
   }

   return unstructured_ivs;
}


void
GKVlasovAMG::addUnstructuredGraphEntries( const int                                 a_radius,
                                          const LayoutData< BaseFab<IntVectSet> >&  a_unstructured_coupling,
                                          HYPRE_SStructGraph&                       a_graph ) const
{
   if ( m_mb_coupling ) {

      const DisjointBoxLayout& grids = m_geometry.gridsFull();
      const RefCountedPtr<PhaseCoordSys>& coord_sys_ptr = m_geometry.phaseCoordSysPtr();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         int dst_block_number = coord_sys_ptr->whichBlock(grids[dit]);

         // Find all cells that may be coupled through a block interface
         IntVectSet dst_ivs = getInterBlockCoupledCells(dst_block_number, a_radius, grids[dit]);

         for (IVSIterator it(dst_ivs); it.ok(); ++it) {
            IntVect iv = it();

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


void GKVlasovAMG::multiplyMatrix( const LevelData<FArrayBox>&  a_in,
                                  LevelData<FArrayBox>&        a_out ) const
{
   if (!m_hypre_allocated) {
      MayDay::Error( "GKVlasovAMG::multiplyMatrix(): Operator has not yet been initialized!" );
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

   const DisjointBoxLayout& grids = a_in.disjointBoxLayout();
   const RefCountedPtr<PhaseCoordSys>& coord_sys_ptr = m_geometry.phaseCoordSysPtr();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = coord_sys_ptr->whichBlock(grids[dit]);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      FArrayBox tmp(box,1);
      tmp.copy(a_in[dit]);

      HYPRE_SStructVectorSetBoxValues(in_vector, block_number, lower.dataPtr(), upper.dataPtr(),
                                      0, tmp.dataPtr());
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

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = coord_sys_ptr->whichBlock(grids[dit]);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      FArrayBox tmp(box,1);

      HYPRE_SStructVectorGetBoxValues(out_vector, block_number,
                                      lower.dataPtr(), upper.dataPtr(), 0, tmp.dataPtr());

      a_out[dit].copy(tmp);
   }

   HYPRE_SStructVectorDestroy(out_vector);
   HYPRE_SStructVectorDestroy(in_vector);
}


void GKVlasovAMG::setConvergenceParams( const double  a_tol,
                                        const int     a_max_iter,
                                        const bool    a_verbose )
{
   m_tol       = a_tol;
   m_max_iter  = a_max_iter;
   m_verbose   = a_verbose;

   m_convergence_params_set = true;
}


void GKVlasovAMG::solve( const LevelData<FArrayBox>&  a_rhs,
                         LevelData<FArrayBox>&        a_solution )
{
   CH_TIME("GKVlasovAMG::solve");
   if ( !m_convergence_params_set ) {
      MayDay::Error("GKVlasovAMG::solve(): solver convergence parameters have not been set");
   }

   copyToHypreVector(a_rhs, m_b);
   HYPRE_SStructVectorAssemble(m_b);

   copyToHypreVector(a_solution, m_x);
   HYPRE_SStructVectorAssemble(m_x);

   AMG( m_A, m_b, m_tol, m_max_iter, m_verbose, m_x );

   copyFromHypreVector(m_x, a_solution);

   // The matrix of the linear system just solved acts on vectors in the physical
   // frame, so convert the solution back to the computational frame
   m_geometry.multJonValid(a_solution);
}


void GKVlasovAMG::AMG( const HYPRE_SStructMatrix&  a_matrix,
                       const HYPRE_SStructVector&  a_b,
                       const double                a_tol,
                       const int                   a_max_iter,
                       const bool                  a_verbose,
                       const HYPRE_SStructVector&  a_x )
{
   CH_TIMERS("GKVlasovAMG::AMG");
   CH_TIMER("AMG_solve",t_AMG_solve);

   HYPRE_ParCSRMatrix par_A;
   HYPRE_ParVector    par_b;
   HYPRE_ParVector    par_x;

   HYPRE_SStructMatrixGetObject(a_matrix, (void **) &par_A);
   HYPRE_SStructVectorGetObject(a_b, (void **) &par_b);
   HYPRE_SStructVectorGetObject(a_x, (void **) &par_x);

   HYPRE_BoomerAMGSetTol(m_par_AMG_solver, a_tol);
   HYPRE_BoomerAMGSetMaxIter(m_par_AMG_solver, a_max_iter);
   HYPRE_BoomerAMGSetPrintLevel(m_par_AMG_solver, 1);

   CH_START(t_AMG_solve);
   HYPRE_BoomerAMGSolve(m_par_AMG_solver, par_A, par_b, par_x);
   CH_STOP(t_AMG_solve);

   int num_iterations;
   HYPRE_BoomerAMGGetNumIterations(m_par_AMG_solver, &num_iterations);
   double final_res_norm;
   HYPRE_BoomerAMGGetFinalRelativeResidualNorm(m_par_AMG_solver, &final_res_norm);
   if (a_verbose && procID()==0) {
      cout << "        --> GKVlasovAMG solver residual = " << final_res_norm << " after " 
           << num_iterations << " iterations" << endl;
   }

   if (m_hypre_object_type == HYPRE_PARCSR) {
      HYPRE_SStructVectorGather(a_x);
   }
}

   
void GKVlasovAMG::copyToHypreVector( const LevelData<FArrayBox>&  a_in,
                                     HYPRE_SStructVector&         a_out ) const
{
   const DisjointBoxLayout & grids = a_in.disjointBoxLayout();
   const RefCountedPtr<PhaseCoordSys>& coord_sys_ptr = m_geometry.phaseCoordSysPtr();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = coord_sys_ptr->whichBlock(grids[dit]);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      FArrayBox tmp(box,1);
      tmp.copy(a_in[dit]);

      HYPRE_SStructVectorSetBoxValues(a_out, block_number, lower.dataPtr(), upper.dataPtr(),
                                      0, tmp.dataPtr());
   }
} 



void GKVlasovAMG::copyFromHypreVector( const HYPRE_SStructVector&  a_in,
                                       LevelData<FArrayBox>&       a_out ) const
{
   const DisjointBoxLayout & grids = a_out.disjointBoxLayout();
   const RefCountedPtr<PhaseCoordSys>& coord_sys_ptr = m_geometry.phaseCoordSysPtr();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box & box = grids[dit];
      int block_number = coord_sys_ptr->whichBlock(grids[dit]);

      FArrayBox tmp(box,1);

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      HYPRE_SStructVectorGetBoxValues(a_in, block_number, lower.dataPtr(), upper.dataPtr(),
                                      0, tmp.dataPtr());

      a_out[dit].copy(tmp);
   }
} 


#include "NamespaceFooter.H"

