#include <sstream>
#include "Directions.H"
#include "PhaseGrid.H"
#include "BoxIterator.H"
#include "DataIterator.H"
#include "SPACE.H"

#include "NamespaceHeader.H"


PhaseGrid::PhaseGrid(const Vector<ProblemDomain>&      a_domains,
                     const Vector<IntVect>&            a_decomps,
                     const vector<int>&                a_legacy_decomp,
                     const string&                     a_mag_geom_type)
{
   Vector<Box> phase_boxes;

   for (int block=0; block<a_domains.size(); ++block) {

      const Box& domain_box = a_domains[block].domainBox();

      IntVect decomp = a_decomps[block];

      //        Chop up the domain box in the configuration space directions over the current
      //        number of processors.

      int nproc = 1;
      for (int dir=0; dir<PDIM; ++dir) {
         nproc *= decomp[dir];
      }

      int n_loc[PDIM];
      for (int dir=0; dir<PDIM; ++dir) {
         n_loc[dir] = 0;
      }

      if ( a_mag_geom_type != "SingleNull" && a_mag_geom_type != "SNCore" ) {
         if (a_legacy_decomp.size() > 0) {
            int nproc_phase = 1;
            for (int dir=0; dir<PDIM; ++dir) {
               nproc_phase *= a_legacy_decomp[dir];
            }
            for (int dir=0; dir<PDIM; ++dir) {
               if (domain_box.size(dir)%a_legacy_decomp[dir] != 0) {
                  stringstream msg("Decomposition in phase direction ", ios_base::out|ios_base::ate);
                  msg << dir << " does not evenly divide domain dimension";
                  MayDay::Error( msg.str().c_str() );
               }
               else {
                  n_loc[dir] = domain_box.size(dir) / a_legacy_decomp[dir];
               }
            }
         }
         else {

            /*
              Look for factors p and q of nproc that evenly divide
              the box configuration dimensions.  If there are
              multiple (p,q) pairs that qualify, use the one for which
              p and q are the closest.
            */
            int q;

            for (q=1; q<=(int)sqrt(1.0*nproc); q++) {
               if (nproc%q == 0) {    // will always true at least for q = 1
                  int p = nproc / q;
                  if (domain_box.size(0)%q == 0 && domain_box.size(1)%p == 0) {
                     n_loc[POLOIDAL_DIR] = domain_box.size(1) / p;
                     n_loc[RADIAL_DIR] = domain_box.size(0) / q;
                  }
               }
            }

            n_loc[VPARALLEL_DIR] = domain_box.size(2);
            n_loc[MU_DIR] = domain_box.size(3);
         }
      }
      else {
         int nproc_phase = 1;
         for (int dir=0; dir<PDIM; ++dir) {
            nproc_phase *= decomp[dir];
         }
         for (int dir=0; dir<PDIM; ++dir) {
            int decomp_dir = decomp[dir];
            if (domain_box.size(dir)%decomp_dir != 0) {
               stringstream msg("Decomposition in phase direction ", ios_base::out|ios_base::ate);
               msg << dir << " does not evenly divide domain dimension";
               MayDay::Error( msg.str().c_str() );
            }
            else {
               n_loc[dir] = domain_box.size(dir) / decomp_dir;
            }
         }
      }

      if (n_loc[RADIAL_DIR] != 0 && n_loc[POLOIDAL_DIR] != 0 && n_loc[VPARALLEL_DIR] != 0 && n_loc[MU_DIR] != 0) {

         IntVect box_size(D_DECL6(n_loc[RADIAL_DIR],n_loc[POLOIDAL_DIR],n_loc[VPARALLEL_DIR],n_loc[MU_DIR],1,1));

         for (int dir=0; dir<PDIM; ++dir) {
            if (box_size[dir] < 4) {
               MayDay::Error( "Phase space box is less than 4 cells wide" );
            }
         }

         Box patch(domain_box.smallEnd(), domain_box.smallEnd() + box_size - IntVect::Unit);
         Box skeleton(IntVect::Zero, IntVect(D_DECL6(domain_box.size(0)/n_loc[RADIAL_DIR]-1,
                                                     domain_box.size(1)/n_loc[POLOIDAL_DIR]-1,
                                                     domain_box.size(2)/n_loc[VPARALLEL_DIR]-1,
                                                     domain_box.size(3)/n_loc[MU_DIR]-1,
                                                     0,0)));
         BoxIterator bit(skeleton);
         for (bit.begin();bit.ok();++bit) {
            Box thisBox = patch + bit()*box_size;
            phase_boxes.push_back(thisBox);
         }
      }
      else {
         MayDay::Error( "Phase space domain box cannot be load balanced" );
      }
   }

   List<CFG::Box> config_boxes;
   getConfigBoxes(phase_boxes, config_boxes);

   // Make the layout.  This is where boxes are assigned to processes.
   Vector<int> procMap;
#if CH_MPI
   assignPhaseDecompositionToProcs(phase_boxes, config_boxes, procMap, m_local_config_boxes, m_local_comm);
#else
   assignPhaseDecompositionToProcs(phase_boxes, config_boxes, procMap, m_local_config_boxes);
#endif

   ProblemDomain prob_domain;

   if ( a_mag_geom_type == "Miller" || a_mag_geom_type == "Slab" ) {
      prob_domain = ProblemDomain(a_domains[0]);
   }
   else if ( a_mag_geom_type == "SingleNull" || a_mag_geom_type == "SNCore" ) {
      Box bounding_box;
      for (int n=0; n<phase_boxes.size(); n++) {
         bounding_box = minBox(bounding_box, phase_boxes[n]);
      }
      prob_domain = ProblemDomain(bounding_box);
   }
   else {
      MayDay::Error("Invalid magnetic geometry type");
   }

   m_dbl.define(phase_boxes, procMap, prob_domain);
   m_dbl.close();

   getVelocitySlices(m_dbl, m_local_config_boxes, m_local_velocity_slices);

#if 0
   for (int i=0; i<m_local_velocity_slices.size(); ++i) {
      cout << procID() << " " << m_local_config_boxes[i] << ": ";
      List<VEL::Box>& this_list = m_local_velocity_slices[i];
      ListIterator<VEL::Box> it(this_list);
      for (it.begin(); it.ok(); ++it) {
         cout << " " << it();
      }
      cout << endl;
   }
#endif
}



PhaseGrid::~PhaseGrid()
{
}



void
PhaseGrid::getConfigBoxes(const Vector<Box>& a_phase_boxes,
                          List<CFG::Box>&    a_config_boxes) const
{
   // Make a list of the configuration space boxes corresponding
   // to the input phase space boxes.

   for (int i=0; i<a_phase_boxes.size(); ++i) {
      CFG::Box cbox;
      projectPhaseToConfiguration(a_phase_boxes[i], cbox);
      if ( !a_config_boxes.includes(cbox) ) {
         a_config_boxes.add(cbox);
      }
   }
}


#if CH_MPI
void
PhaseGrid::assignPhaseDecompositionToProcs(const Vector<Box>& a_phase_boxes,
                                           List<CFG::Box>&    a_config_boxes,
                                           Vector<int>&       a_procMap,
                                           Vector<CFG::Box>&  a_local_config_boxes,
                                           Vector<MPI_Comm>&  a_local_comm) const
{
   int nproc = numProc();

   MPI_Group world_group;
   MPI_Comm_group(MPI_COMM_WORLD, &world_group);

   int num_groups = a_config_boxes.length();
   MPI_Group* groups = new MPI_Group[num_groups];

   a_procMap.resize(a_phase_boxes.size());

   if ( nproc < num_groups ) {

      // First, assign processors to configuration boxes
      int* config_procs = new int[num_groups];

      for (int group_number=0; group_number<num_groups; ++group_number) {
         config_procs[group_number] = group_number % nproc;
      }

      // Then assign each phase space box to the same processor as its
      // corresponding configuration box and construct the MPI groups
      int group_number = 0;
      ListIterator<CFG::Box> it(a_config_boxes);
      for (it.begin(); it.ok(); ++it) {
         CFG::Box& config_box = it();

         for (int k=0; k<a_phase_boxes.size(); ++k) {
            CFG::Box cbox;
            projectPhaseToConfiguration(a_phase_boxes[k], cbox);

            if (cbox == config_box) {
               a_procMap[k] = config_procs[group_number];
            }
         }

         MPI_Group_incl(world_group, 1, &config_procs[group_number], &groups[group_number]);
         
         group_number++;
      }

      delete [] config_procs;
   }
   else {

      // Figure out how many processors will be allocated to each group
      int* num_group_processors = new int[num_groups];

      int min_procs_per_group = floor(nproc / num_groups);
      for (int group_number=0; group_number<num_groups; ++group_number) {
         num_group_processors[group_number] = min_procs_per_group;
      }
      int remaining_procs = nproc - min_procs_per_group * num_groups;
      int group_number = 0;
      while ( remaining_procs-- > 0 ) {
         num_group_processors[group_number++]++;
      }

      int max_procs_per_group = ceil(nproc / num_groups);
      int* ranks = new int[max_procs_per_group];
      int proc = 0;
      group_number = 0;
      ListIterator<CFG::Box> it(a_config_boxes);
      for (it.begin(); it.ok(); ++it) {
         CFG::Box& config_box = it();

         int group_size = 0;
         int num_boxes_found = 0;
         for (int k=0; k<a_phase_boxes.size(); ++k) {
            CFG::Box cbox;
            projectPhaseToConfiguration(a_phase_boxes[k], cbox);

            if (cbox == config_box) {
               int local_proc = num_boxes_found % num_group_processors[group_number];
               a_procMap[k] = proc + local_proc;
               num_boxes_found++;
               if (num_boxes_found <= num_group_processors[group_number]) {
                  ranks[group_size] = a_procMap[k];
                  group_size++;
               }
            }
         }

         MPI_Group_incl(world_group, group_size, ranks, &groups[group_number]);

         proc += num_group_processors[group_number];
         group_number++;
      }

      delete [] ranks;
      delete [] num_group_processors;
   }

   List<CFG::Box> box_list;
   List<MPI_Comm> comm_list;

   int group_number = 0;
   ListIterator<CFG::Box> it(a_config_boxes);
   for (it.begin(); it.ok(); ++it) {
      CFG::Box& config_box = it();

      MPI_Comm newcomm;
      MPI_Comm_create(MPI_COMM_WORLD, groups[group_number], &newcomm);

      // A valid communicator is only returned by the preceding call on processors
      // contained in the current group.  We exploit this to construct the list
      // of groups this processor belongs to and the corresponding communicators.
      if (newcomm != MPI_COMM_NULL) {
         box_list.add(config_box);
         comm_list.add(newcomm);
      }

      group_number++;
   }

   // Now that we know how big the lists are, copy them to vectors for easier indexing

   a_local_config_boxes.resize(box_list.length());
   a_local_comm.resize(comm_list.length());

   int index = 0;
   for (ListIterator<CFG::Box> it(box_list); it.ok(); ++it) {
      a_local_config_boxes[index++] = it();
   }

   index = 0;
   for (ListIterator<MPI_Comm> it(comm_list); it.ok(); ++it) {
      a_local_comm[index++] = it();
   }

   delete [] groups;
}
#else

void
PhaseGrid::assignPhaseDecompositionToProcs(const Vector<Box>& a_phase_boxes,
                                           List<CFG::Box>&    a_config_boxes,
                                           Vector<int>&       a_procMap,
                                           Vector<CFG::Box>&  a_local_config_boxes) const
{
   int num_phase_boxes = a_phase_boxes.size();

   a_procMap.resize(num_phase_boxes);
   for (int i=0; i<num_phase_boxes; ++i) {
      a_procMap[i] = 0;
   }

   a_local_config_boxes.resize(a_config_boxes.length());

   int index = 0;
   ListIterator<CFG::Box> it(a_config_boxes);
   for (it.begin(); it.ok(); ++it) {
      a_local_config_boxes[index++] = it();
   }
}

#endif



void
PhaseGrid::getVelocitySlices(const DisjointBoxLayout&  a_dbl,
                             const Vector<CFG::Box>&   a_config_boxes,
                             Vector< List<VEL::Box> >& a_velocity_slices) const
{
   int num_config_boxes = a_config_boxes.size();

   a_velocity_slices.resize(num_config_boxes);

   for (int i=0; i<num_config_boxes; ++i) {
      const CFG::Box& config_box = a_config_boxes[i];

      List<VEL::Box> this_list;

      for (DataIterator dit(a_dbl); dit.ok(); ++dit) {
         const Box& box = a_dbl[dit];

         CFG::Box cbox;
         projectPhaseToConfiguration(box, cbox);

         if (cbox == config_box) {
            VEL::Box vbox;
            projectPhaseToVelocity(box, vbox);
            this_list.add(vbox);
         }
      }

      a_velocity_slices[i] = this_list;
   }
}



void
PhaseGrid::projectPhaseToConfiguration( const Box& a_src,
                                        CFG::Box&  a_dst ) const
{
  a_dst.define(config_restrict(a_src.smallEnd()),
               config_restrict(a_src.bigEnd()),
               config_restrict(a_src.type()));
}



void
PhaseGrid::projectPhaseToVelocity( const Box& a_src,
                                   VEL::Box&  a_dst ) const
{
  a_dst.define(vel_restrict(a_src.smallEnd()),
               vel_restrict(a_src.bigEnd()),
               vel_restrict(a_src.type()));
}



CFG::IntVect
PhaseGrid::config_restrict( const IntVect& a_x ) const
{
   CFG::IntVect restrict_x;

   for (int i=0; i<CFG_DIM; i++) {
      restrict_x[i] = a_x[i];
   }

   return restrict_x;
}



VEL::IntVect
PhaseGrid::vel_restrict( const IntVect& a_x ) const
{
   VEL::IntVect restrict_x;

   for (int i=VPARALLEL_DIR; i<=MU_DIR; i++) {
      restrict_x[i-VPARALLEL_DIR] = a_x[i];
   }

   return restrict_x;
}



void
PhaseGrid::print(const IntVect& a_ghostVect) const
{
   for (DataIterator dit(m_dbl); dit.ok(); ++dit) {
      const Box& box = m_dbl[dit];
      Box grown_box = box;
      grown_box.grow(a_ghostVect);
      int num_valid_cells = box.numPts();
      int num_ghost_cells = grown_box.numPts() - num_valid_cells;
      cout << "   Phase space box " << box << " is assigned to process " << m_dbl.procID(dit())
           << ", has " << num_valid_cells << " valid cells and " << num_ghost_cells
           << " ghost cells" << endl;
   }
}

#include "NamespaceFooter.H"
