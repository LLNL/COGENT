#include "SNCoreSystem.H"
#include "SNCoreSystemF_F.H"
#include "CONSTANTS.H"

#define TIME_MULTIBLOCK

#define MULTI_BLOCK

#include <fstream>
#include <sstream>
#if 1  // warning, OS dependencies, will not work on all platforms
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "CH_HDF5.H"
#include "HDF5Portable.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "MagBlockCoordSys.H"
#include "SNCoreBlockCoordSys.H"
#include "SNCoreCoordSys.H"
#include "newMappedGridIO.H"
#include "FourthOrderUtil.H"
#include "Directions.H"
#include "inspect.H"

#include "MultiBlockLevelExchangeCenter.H"

#include "LocSpatialIBCF_F.H"

#define NORADDECOMP
#undef NORADDECOMP


#include "NamespaceHeader.H"

#ifdef MULTI_BLOCK
MultiBlockLevelExchangeAverage * mblexPtr = NULL;
#endif


SNCoreSystem::SNCoreSystem( ParmParse& a_pp )
   : m_verbosity(0),
     m_hdf_density(true),
     m_hdf_pressure(false),
     m_hdf_vpartheta(false),
     m_hdf_frtheta(false),
     m_hdf_vparmu(false),
     m_ghostVect(4*IntVect::Unit),
     m_ppgksys("gksystem"),
     m_enforce_stage_positivity(false),
     m_enforce_step_positivity(false)
{
   // Read the input
   parseParameters( m_ppgksys );

   /*
     Create the configuration space (magnetic geometry) domains,
     data layout and coordinate system.
   */

   createConfigurationSpace();

   m_operator = new Advect( ParmParse(Advect::pp_name) );

   m_solution.define( m_mag_geom->grids(), 1, IntVect::Zero);

   ParmParse pp_ibc("IBC");

   m_ibc = new LocSpatialIBC( pp_ibc, *m_mag_geom, 0);
}



void
SNCoreSystem::initialize(const int cur_step)
{
   if (procID()==0) cout << "Entering initialize" << endl;

   // Set the initial conditions if this is the first step
   if (cur_step == 0) {
      setInitialConditions( m_solution );
   }

   // Set the boundary data.  m_solution is not modified by this call (since
   // it now contains the initial condition or the solution obtain from a restart),
   // but simply provides a template for the boundary condition object.
   setBoundaryData( m_solution );

   // Make the solution prototype to hold old mapped and real data.
   m_solution_mapped.define(m_solution.disjointBoxLayout(), m_solution.nComp(), m_ghostVect);
   m_solution_mapped_old.define(m_solution.disjointBoxLayout(), m_solution.nComp(), m_ghostVect);

   DataIterator dit = m_solution.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      m_solution_mapped[dit].copy(m_solution[dit]);
      m_solution_mapped_old[dit].copy(m_solution[dit]);
   }

   divideJ(m_solution, m_solution);

#if 0
   double l1_error, l2_error, max_error;
   computeError(0., cur_step, m_solution, l1_error, l2_error, max_error);
   if (procID()==0) cout << "l1, l2, max error at t = 0. is " << l1_error << " " << l2_error << " " << max_error << endl;
#endif

   // Create the integrator
   m_integrator = new RK4Integrator<LevelData<FArrayBox> >(*this, m_solution_mapped);

   MPI_Barrier(MPI_COMM_WORLD);
   if (procID()==0) cout << "Exiting initialize" << endl;
}



SNCoreSystem::~SNCoreSystem()
{
   delete m_mag_geom;
   delete m_mag_geom_coords;
   delete m_operator;
   delete m_integrator;
}



void SNCoreSystem::createConfigurationSpace()
{
   ParmParse pp_grid("sncore");

   string prefix = string("gksystem.magnetic_geometry_mapping.")
     + string(SNCoreBlockCoordSys::pp_name);
   ParmParse pp_geom( prefix.c_str() );

   m_mag_geom_coords = new SNCoreCoordSys(pp_grid, pp_geom);

   DisjointBoxLayout grids;
   createDisjointBoxLayout( grids );

#if 0
   // For testing of single null mapping inverse

   MPI_Barrier(MPI_COMM_WORLD);
   if (procID()==0) cout << "Testing single null mapping inverse" << endl;

   double start_time = MPI_Wtime();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = m_mag_geom_coords->whichBlock(grids[dit]);
      const SNCoreBlockCoordSys* block_coord_sys
         = dynamic_cast<const SNCoreBlockCoordSys*>(m_mag_geom_coords->getCoordSys(block_number));

      RealVect dx = block_coord_sys->dx();
      RealVect offset = 0.5*RealVect::Unit;
      offset *= dx;
      RealVect offset2 = 0.3*RealVect::Unit;
      offset2 *= dx;

      BoxIterator bit(grids[dit]);
      for (bit.begin();bit.ok();++bit) {
         IntVect iv = bit();
         RealVect xi = dx*iv + offset2;
         RealVect X = block_coord_sys->realCoord(xi);
         RealVect xi2 = block_coord_sys->mappedCoord(X);
         RealVect diff = xi2 - xi;
         if (diff.vectorLength() > 1.e-8) cout << block_number << ", iv = " << iv << ", X = " << X << ", xi = " << xi << ", xi2 = " << xi2 << " " << diff.vectorLength() << endl;
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);

   if (procID()==0) cout << "Time = " << MPI_Wtime() - start_time << endl;

   exit(1);
#endif

#if 0
   // For testing of single null block remapping member

   MPI_Barrier(MPI_COMM_WORLD);
   if (procID()==0) cout << "Testing single null block remapping" << endl;

   DataIterator dit2 = grids.dataIterator();
   for (dit2.begin(); dit2.ok(); ++dit2) {
      int block_number = m_mag_geom_coords->whichBlock(grids[dit2]);
      const SNCoreBlockCoordSys* block_coord_sys
         = dynamic_cast<const SNCoreBlockCoordSys*>(m_mag_geom_coords->getCoordSys(block_number));
      RealVect dx = block_coord_sys->dx();
      RealVect offset = 0.3*RealVect::Unit;
      offset *= dx;

      BoxIterator bit(grids[dit2]);
      for (bit.begin();bit.ok();++bit) {
         IntVect iv = bit();
         RealVect xi = dx*iv + offset;
         RealVect xi_valid;
         int n_valid;
         m_mag_geom_coords->blockRemapping(xi_valid, n_valid, xi, block_number);
         if (n_valid != block_number) {
            cout << "block = " << block_number << ", n_valid = " << n_valid << endl;
         }
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);

   exit(1);
#endif

#if 0

   MPI_Barrier(MPI_COMM_WORLD);
   if (procID()==0) cout << "Writing geom_data plot files" << endl;

   // For testing of metrics
   IntVect geom_data_ghosts = 48*IntVect::Unit;
   LevelData<FArrayBox> geom_data(grids, 6, geom_data_ghosts);

   MPI_Barrier(MPI_COMM_WORLD);
   if (procID()==0) cout << "Constructing metrics" << endl;
   double start_time = MPI_Wtime();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = m_mag_geom_coords->whichBlock(grids[dit]);
      const MagBlockCoordSys* block_coord_sys
         = dynamic_cast<const MagBlockCoordSys*>(m_mag_geom_coords->getCoordSys(block_number));

      RealVect dx = block_coord_sys->dx();
      RealVect offset = 0.5*RealVect::Unit;
      offset *= dx;

      FArrayBox& this_geom_data = geom_data[dit];
      BoxIterator bit(geom_data[dit].box());
      for (bit.begin();bit.ok();++bit) {
         IntVect iv = bit();
         RealVect xi = dx*iv + offset;
         RealVect X = block_coord_sys->realCoord(xi);
         this_geom_data(iv,0) = X[RADIAL_DIR];
         this_geom_data(iv,1) = X[POLOIDAL_DIR];

         this_geom_data(iv,2) = block_coord_sys->dXdXi(xi,0,0);
         this_geom_data(iv,3) = block_coord_sys->dXdXi(xi,0,1);
         this_geom_data(iv,4) = block_coord_sys->dXdXi(xi,1,0);
         this_geom_data(iv,5) = block_coord_sys->dXdXi(xi,1,1);
      }
      inspect(this_geom_data);
   }

   MPI_Barrier(MPI_COMM_WORLD);
   if (procID()==0) cout << "Done constructing metrics" << endl;
   if (procID()==0) cout << "Time = " << MPI_Wtime() - start_time << endl;

   Box plot_box = m_domain.domainBox();
   plot_box.grow(geom_data_ghosts);

   WriteMappedUGHDF5("geom_data", grids, geom_data, *m_mag_geom_coords, plot_box);
   //   exit(1);
#endif

#if 0

   MPI_Barrier(MPI_COMM_WORLD);
   if (procID()==0) cout << "Writing Jdata plot files" << endl;

   IntVect Jghosts = 4*IntVect::Unit;

   LevelData<FArrayBox> Jdata(grids, 1, Jghosts);

   m_mag_geom_coords->getJ(Jdata);

   Box plot_box = m_domain.domainBox();
   plot_box.grow(Jghosts);

   WriteMappedUGHDF5("Jdata", grids, Jdata, *m_mag_geom_coords, plot_box);
   exit(1);
#endif

   // Construct the magnetic geometry
   if (procID()==0) cout << "Constructing MultiBlockLevelGeom for magnetic geometry" << endl;
   int ghosts = 0;
   for (int n=0; n<CFG_DIM; ++n) {
      if (m_ghostVect[n] > ghosts) ghosts = m_ghostVect[n];
   }
   int spaceOrder = 4;   // = 2*L with L as in Peter's notes
   m_mag_geom = new MultiBlockLevelGeom();
   m_mag_geom->define(m_mag_geom_coords, grids, ghosts, spaceOrder);
   if (procID()==0) cout << "Done constructing MultiBlockLevelGeom for magnetic geometry" << endl;

   MPI_Barrier(MPI_COMM_WORLD);

#if 0
   // For testing of multiblock exchange
   //   MultiBlockLevelExchangeAverage* mblexPtr = new MultiBlockLevelExchangeAverage();
   MultiBlockLevelExchangeCenter* mblexPtr = new MultiBlockLevelExchangeCenter();
   mblexPtr->define(m_mag_geom, ghosts, spaceOrder);

   LevelData<FArrayBox> tmp_data(grids, 1, 4*IntVect::Unit);

   DataIterator dit2 = grids.dataIterator();
   for (dit2.begin(); dit2.ok(); ++dit2) {
      int block_number = m_mag_geom_coords->whichBlock(grids[dit2]);
      FArrayBox& this_tmp_data = tmp_data[dit2];

      const MagBlockCoordSys* block_coord_sys
         = dynamic_cast<const MagBlockCoordSys*>(m_mag_geom_coords->getCoordSys(block_number));

      this_tmp_data.setVal(0.);

      RealVect dx = block_coord_sys->dx();
      RealVect offset = 0.5*RealVect::Unit;
      offset *= dx;

      BoxIterator bit(grids[dit2]);
      for (bit.begin();bit.ok();++bit) {
         IntVect iv = bit();
         RealVect xi = dx*iv + offset;
         RealVect X = block_coord_sys->realCoord(xi);
         this_tmp_data(iv,0) = X[POLOIDAL_DIR]*X[POLOIDAL_DIR]*X[POLOIDAL_DIR]*X[POLOIDAL_DIR];
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);

   tmp_data.exchange();
   mblexPtr->interpGhosts(tmp_data);

   MPI_Barrier(MPI_COMM_WORLD);

   Box plot_box2(m_domain.domainBox());
   //   plot_box2.grow(tmp_data.ghostVect());

#if 0
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = m_mag_geom_coords->whichBlock(grids[dit]);
      FArrayBox& this_tmp_data = tmp_data[dit];
      const MagBlockCoordSys* block_coord_sys
         = dynamic_cast<const MagBlockCoordSys*>(m_mag_geom_coords->getCoordSys(block_number));

      RealVect dx = block_coord_sys->dx();
      RealVect offset = 0.5*RealVect::Unit;
      offset *= dx;

      Box valid_box = grids[dit];
      BoxIterator bit(this_tmp_data.box());
      for (bit.begin();bit.ok();++bit) {
         IntVect iv = bit();
         if (grids[dit].contains(iv)) {
           this_tmp_data(iv,0) = 0.;
         }
         else {
           RealVect xi = dx*iv + offset;
           RealVect X = block_coord_sys->realCoord(xi);
           double val = fabs( this_tmp_data(iv,0) - X[POLOIDAL_DIR]*X[POLOIDAL_DIR]*X[POLOIDAL_DIR]*X[POLOIDAL_DIR] );
#if 0
           if ( !plot_box2.contains(iv) ) val = 0.;
           if ( (block_number == 0 || block_number == 2) && iv[1] < valid_box.smallEnd(1) ) {
             val = 0.;
           }
           if ( (block_number == 1 || block_number == 3) && iv[1] > valid_box.bigEnd(1) ) {
             val = 0.;
           }
#endif
           this_tmp_data(iv,0) = val;
         }
      }
   }
#endif

   MPI_Barrier(MPI_COMM_WORLD);

   WriteMappedUGHDF5("tmp_data", grids, tmp_data, *m_mag_geom_coords, plot_box2);

   exit(1);
#endif
}



void SNCoreSystem::createDisjointBoxLayout( DisjointBoxLayout& grids )
{
   Vector<Box> boxes;
   for (int block=0; block<m_mag_geom_coords->numBlocks(); ++block) {

      const MagBlockCoordSys* mag_block_coords
       = (MagBlockCoordSys *)m_mag_geom_coords->getCoordSys(block);

      const ProblemDomain& domain = mag_block_coords->domain();

      const IntVect decomp = m_mag_geom_coords->getDecomposition(block);

      const Box& domain_box = domain.domainBox();

      /*
        Chop up the configuration space domain box over the number of processors specified
        for this block.  At this point, we insist that the box decomposes uniformly, or an
        error is thrown.
      */
      int nproc = 1;
      for (int dir=0; dir<CFG_DIM; ++dir) {
         nproc *= decomp[dir];
      }

      int n_loc[CFG_DIM];
      for (int dir=0; dir<CFG_DIM; ++dir) {
         n_loc[dir] = 0;
      }

#if CFG_DIM==2

      for (int dir=0; dir<CFG_DIM; ++dir) {
         int decomp_dir = decomp[dir];
         if (domain_box.size(dir)%decomp_dir != 0) {
            stringstream msg("Decomposition in configuration direction ", ios_base::out|ios_base::ate);
            msg << dir << " does not evenly divide domain dimension";
            MayDay::Error( msg.str().c_str() );
         }
         else {
            n_loc[dir] = domain_box.size(dir) / decomp_dir;
         }
      }

      if (n_loc[RADIAL_DIR] != 0 && n_loc[POLOIDAL_DIR] != 0) {

         IntVect box_size(n_loc[RADIAL_DIR],n_loc[POLOIDAL_DIR]);
         Box patch(domain_box.smallEnd(), domain_box.smallEnd() + box_size-1);
         Box skeleton(IntVect::Zero, IntVect(domain_box.size(0)/n_loc[RADIAL_DIR]-1,domain_box.size(1)/n_loc[POLOIDAL_DIR]-1));
         BoxIterator bit(skeleton);
         for (bit.begin();bit.ok();++bit) {
            Box thisBox = patch + bit()*box_size;
            boxes.push_back(thisBox);
         }
      }
      else {
         MayDay::Error( "Configuration domain box cannot be load balanced" );
      }
#else
      MayDay::Error( "SNCoreSystem::createConfigurationSpace has only been implemented for 2D" );
#endif
   }

   // Make the layout.  This is where boxes are assigned to processes.
   Vector<int> procMap;
   LoadBalance( procMap, boxes );

   Box bounding_box;
   for (int n=0; n<boxes.size(); n++) {
      bounding_box = minBox(bounding_box, boxes[n]);
   }
   m_domain = ProblemDomain(bounding_box);

   grids.define( boxes, procMap, m_domain );
   grids.close();

#ifdef CH_MPI
   if (procID()==0) {
#endif
     if (m_verbosity>0) {
       for (int n=0; n<boxes.size(); n++) {
         Box grown_box = boxes[n];
         grown_box.grow(m_ghostVect);
         int num_valid_cells = boxes[n].numPts();
         int num_ghost_cells = grown_box.numPts() - num_valid_cells;
         cout << "   Box " << boxes[n] << " is assigned to process " << procMap[n]
              << ", has " << num_valid_cells << " valid cells and " << num_ghost_cells
              << " ghost cells" << endl;
       }
     }
#ifdef CH_MPI
   }
#endif
}



void SNCoreSystem::defineRHSData( LevelData<FArrayBox>&       a_rhs,
                              const LevelData<FArrayBox>& a_prototype )
{
   if ( a_rhs.isDefined() ) {
      MayDay::Error("SNCoreSystem::defineRHSData(): a_rhs is already defined");
   }

   a_rhs.define( a_prototype.disjointBoxLayout(),
                 a_prototype.nComp(),
                 IntVect::Zero );

   DataIterator dit = a_prototype.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_rhs[dit].copy(a_prototype[dit]);
   }
}



void SNCoreSystem::defineSolnData( LevelData<FArrayBox>&       a_soln,
                               const LevelData<FArrayBox>& a_prototype )
{
   if ( a_soln.isDefined() ) {
      MayDay::Error("SNCoreSystem::defineSolnData(): a_soln is already defined");
   }

   a_soln.define( a_prototype.disjointBoxLayout(),
                  a_prototype.nComp(),
                  IntVect::Zero );

   DataIterator dit = a_prototype.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_soln[dit].copy(a_prototype[dit]);
   }
}



bool SNCoreSystem::validSolnData( const LevelData<FArrayBox>& a_soln,
                              const LevelData<FArrayBox>& a_protoSoln )
{
   return (a_protoSoln.disjointBoxLayout() == a_soln.disjointBoxLayout()) &&
          (a_protoSoln.nComp()             == a_soln.nComp())             &&
          (a_protoSoln.ghostVect()         == a_soln.ghostVect());
}



bool SNCoreSystem::validRHSData( const LevelData<FArrayBox>& a_rhs,
                             const LevelData<FArrayBox>& a_protoRHS )
{
   return (a_protoRHS.disjointBoxLayout() == a_rhs.disjointBoxLayout()) &&
          (a_protoRHS.nComp()             == a_rhs.nComp());

}



void SNCoreSystem::zeroSolnData( LevelData<FArrayBox>& a_soln )
{
   DataIterator dit = a_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_soln[dit].setVal(0.);
   }
}



void SNCoreSystem::addSolnData( LevelData<FArrayBox>&       a_soln,
                            const LevelData<FArrayBox>& a_increment,
                            const Real                  a_scale )
{
   DataIterator dit = a_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_soln[dit].plus(a_increment[dit], a_scale);
   }
}



void SNCoreSystem::copySolnData( LevelData<FArrayBox>&       a_dstSoln,
                             const LevelData<FArrayBox>& a_srcSoln )
{
  //   CH_assert( validSolnData(a_dstSoln, a_srcSoln) );

   DataIterator dit = a_srcSoln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_dstSoln[dit].copy(a_srcSoln[dit]);
   }
}



void SNCoreSystem::evalRHS( LevelData<FArrayBox>&       a_rhs,
                        const LevelData<FArrayBox>& a_soln,
                        const int                   a_step_number,
                        const Real                  a_time,
                        const int                   a_stage)
{
#ifdef MULTI_BLOCK
   if ( mblexPtr == NULL ) {
      mblexPtr = new MultiBlockLevelExchangeAverage();
      int numGhost = 4;
      int spaceOrder = 4;

      mblexPtr->setGetConditionNumber(true);

#ifdef TIME_MULTIBLOCK
      MPI_Barrier(MPI_COMM_WORLD);
      double define_time = MPI_Wtime()/60.;
#endif

      if (procID()==0) cout << "Beginning multiblock exchange define in evalRHS" << endl;

      mblexPtr->define(m_mag_geom, numGhost, spaceOrder);

#ifdef TIME_MULTIBLOCK
      define_time = MPI_Wtime()/60. - define_time;

      cout << "Exchange define time on proc " << procID() << " = " << define_time << endl;
#endif

      getConditionNumbers(mblexPtr, *m_mag_geom);
   }
#endif

   // Make a non-const copy of the solution so that we can
   // set boundary values and ghost values (via exchange and
   // multiblock interpolation).
   CH_assert( a_soln.ghostVect() == m_ghostVect );
   const DisjointBoxLayout& grids = a_soln.disjointBoxLayout();

   LevelData<FArrayBox> tmp_soln(grids, a_soln.nComp(), m_ghostVect);
   DataIterator dit = tmp_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     tmp_soln[dit].copy(a_soln[dit], grids[dit]);
   }

   divideJ(tmp_soln, tmp_soln, true);

   tmp_soln.exchange();

#ifdef MULTI_BLOCK
   mblexPtr->interpGhosts(tmp_soln);
#endif

   m_ibc->ghostCellBC(tmp_soln, *m_mag_geom, 0., a_time);

   m_operator->evalRHS( *m_mag_geom, a_rhs, tmp_soln, a_time );
}



void
SNCoreSystem::getConditionNumbers(const MultiBlockLevelExchangeAverage* a_mblexPtr,
                              const MultiBlockLevelGeom&            a_geom ) const
{
   const LayoutData< IntVectSet >& ghostCells = a_mblexPtr->ghostCells();
   const LayoutData< IVSFAB<Real>* >& conditionNumber = a_mblexPtr->conditionNumber();

   // Make a LevelData to store the condition number for plotting.  Condition numbers
   // are only in ghost cells, but we can only plot LevelDatas, so we zero out
   // the valid regions.

   const DisjointBoxLayout& grids = a_geom.grids();

   LevelData<FArrayBox> log_plot(grids, 1, m_ghostVect);

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     log_plot[dit].setVal(0.);
   }

   // Find the largest condition number and its location

   const MultiBlockCoordSys* coords = a_geom.coordSysPtr();
   int num_blocks = coords->numBlocks();
   double max_condition_number[num_blocks];
   IntVect max_condition_number_iv[num_blocks];
   for (int block_number=0; block_number<num_blocks; ++block_number) {
     max_condition_number[block_number] = 0.;
     max_condition_number_iv[block_number] = IntVect::Zero;
   }

   for (dit.begin(); dit.ok(); ++dit) {
     int block_number = coords->whichBlock(grids[dit]);
     const IntVectSet& ghostCellsIVS = ghostCells[dit];

     const IVSFAB<Real>& ivsfab = *conditionNumber[dit];
     FArrayBox& this_plot = log_plot[dit];

     for (IVSIterator ivsit(ghostCellsIVS); ivsit.ok(); ++ivsit) {
       IntVect thisGhostCell = ivsit();
       double condition_number = ivsfab(thisGhostCell,0);

       this_plot(thisGhostCell) = log10(condition_number);

       if (condition_number > max_condition_number[block_number]) {
         max_condition_number[block_number] = condition_number;
         max_condition_number_iv[block_number] = thisGhostCell;
       }
     }
   }

#if CH_MPI
   struct {
     double val;
     int proc;
   } local[num_blocks], global[num_blocks];

   for (int block_number=0; block_number<num_blocks; ++block_number) {
     local[block_number].val = max_condition_number[block_number];
     local[block_number].proc = procID();
   }

   MPI_Allreduce( local, global, num_blocks, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
#endif

   for (int block_number=0; block_number<num_blocks; ++block_number) {
     double max_block_condition_number;
#if CH_MPI
     if (procID() == global[block_number].proc) {
       max_block_condition_number = global[block_number].val;
#else
       max_block_condition_number = max_condition_number[block_number];
#endif
       cout << "Max condition number on block " << block_number << " = " << max_block_condition_number
            << ", which occurred at " << max_condition_number_iv[block_number] << endl;
#if CH_MPI
     }
#endif
   }

   Box plot_box = grids.physDomain().domainBox();
   plot_box.grow(log_plot.ghostVect());

   WriteMappedUGHDF5("condition_number_log", grids, log_plot, *coords, plot_box);
}



Real SNCoreSystem::stableDt(const int a_step_number)
{
   return m_operator->computeDt(*m_mag_geom, m_solution_mapped);
}



void SNCoreSystem::enforcePositivity( LevelData<FArrayBox>& a_soln )
{
   m_positivity_post_processor.enforce( a_soln, 1.0 );
}



Real SNCoreSystem::advance( const Real a_cur_time,
                        const Real a_dt,
                        const int  a_step_number)
{
   double new_time = m_integrator->advance( m_solution_mapped,
                                            m_solution_mapped_old,
                                            a_step_number,
                                            a_cur_time,
                                            a_dt );

   if (m_enforce_step_positivity) {
     enforcePositivity( m_solution_mapped );
   }

   copySolnData( m_solution_mapped_old, m_solution_mapped );

   //   divideJ(m_solution_mapped, m_solution);
   divideJ(m_solution_mapped, m_solution, true);
   //   copySolnData(m_solution, m_solution_mapped);

#if 0
   double l1_error, l2_error, max_error;
   computeError(new_time, a_step_number, m_solution, l1_error, l2_error, max_error);
   if (procID()==0) cout << "l1, l2, max error at t = " << new_time << " is " << l1_error << " " << l2_error << " " << max_error << endl;
#endif

   return new_time;
}



void SNCoreSystem::setInitialConditions( LevelData<FArrayBox>& a_soln )
{
   m_ibc->initialize(a_soln, *m_mag_geom, 0., 0.);
}



void SNCoreSystem::setBoundaryData( LevelData<FArrayBox>& a_soln )
{
   m_ibc->setBoundaryData(a_soln, *m_mag_geom, 0., 0.);
}



void SNCoreSystem::multJ( const LevelData<FArrayBox>& a_soln_physical,
                      LevelData<FArrayBox>&       a_soln_mapped)
{
   // Need to grow the soln factor since we need transverse gradients for the fourth-order product
   const DisjointBoxLayout& grids = a_soln_physical.disjointBoxLayout();
   IntVect grown_vect = a_soln_mapped.ghostVect() + IntVect::Unit;
   LevelData<FArrayBox> grown_soln(grids, 1, grown_vect);

   LevelData<FArrayBox> J(grids, 1, grown_vect);
   m_mag_geom_coords->getJ(J);

   DataIterator dit = grown_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_grown_soln = grown_soln[dit];
      Box this_box = a_soln_mapped[dit].box();

      this_grown_soln.copy(a_soln_physical[dit]);

      for (int dir=0; dir<SpaceDim; ++dir) {
         secondOrderTransExtrap(this_grown_soln, dir, this_box);
      }
   }
   grown_soln.exchange();

   // Compute fourth-order product
   fourthOrderCellProd(a_soln_mapped, grown_soln, J);
}



void SNCoreSystem::divideJ( const LevelData<FArrayBox>& a_soln_mapped,
                        LevelData<FArrayBox>&       a_soln_physical,
                        bool                        a_restrict_to_valid)
{
   CH_assert( a_soln_mapped.ghostVect() >= a_soln_physical.ghostVect() );

   const DisjointBoxLayout& grids = a_soln_mapped.disjointBoxLayout();

   IntVect grown_vect = IntVect::Unit;

   LevelData<FArrayBox> J(grids, 1, grown_vect);
   m_mag_geom_coords->getJ(J);

   LevelData<FArrayBox> grown_soln(grids, 1, grown_vect);

   DataIterator dit = a_soln_physical.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     grown_soln[dit].copy(a_soln_mapped[dit]);
   }
   grown_soln.exchange();

   // We now have soln and J with one layer of ghost cells that have been filled at least
   // on interior box sides.  Compute fourth-order quotient using one-sided differences at
   // physical boundaries.

   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = m_mag_geom_coords->whichBlock(grids[dit]);
      const SNCoreBlockCoordSys* coord_sys
         = dynamic_cast<const SNCoreBlockCoordSys*>(m_mag_geom_coords->getCoordSys(block_number));

      cellFGToCellF(a_soln_physical[dit], grown_soln[dit], J[dit], grids[dit], coord_sys->domain(), true);
   }

   a_soln_physical.exchange();
}



void SNCoreSystem::writePlotFile(const char *prefix, const int cur_step )
{
   if (procID()==0) cout << "Writing plotfile" << endl;
   char iter_str[100], dir_prefix[100];

   stringstream filename;

   sprintf(dir_prefix, "%s", prefix);

#if 1  // warning, OS dependencies, will not work on all platforms
   // place density plots in their own directory;
   sprintf(iter_str, "%s_density_plots", prefix);
#ifdef CH_MPI
   if (procID() == 0) {
#endif
      mkdir(iter_str, 0777); // only works the first time, subsequent failure is normal and expected
#ifdef CH_MPI
      }
#endif
   sprintf(dir_prefix, "%s//%s", iter_str, prefix);
#endif

   stringstream filename1;
   sprintf(iter_str, "%s.density%04d.%d.%s.", dir_prefix, cur_step, 0, "hydrogen");
   filename1 << iter_str;

   const DisjointBoxLayout& grids = m_solution.disjointBoxLayout();

   WriteMappedUGHDF5(filename1.str().c_str(), grids, m_solution, *m_mag_geom_coords, m_domain.domainBox());
   if (procID()==0) cout << "Done writing plotfile" << endl;
}



void SNCoreSystem::writeCheckpointFile( HDF5Handle&  handle,
                                    const int    cur_step,
                                    const double cur_time,
                                    const double cur_dt )
{
}



void SNCoreSystem::readCheckpointFile( HDF5Handle& handle, int* cur_step, double* cur_time, double* cur_dt )
{
}


inline void SNCoreSystem::printParameters()
{
   if (procID() == 0 && m_verbosity) {
      cout << "enforce_positivity = " << (m_enforce_step_positivity||m_enforce_stage_positivity) << endl;
      std::string ptype("stage");
      if (m_enforce_step_positivity)
        ptype = "step";
      cout << "enforce_positivity_type = " << ptype << endl;
   }
}



void SNCoreSystem::writeFieldHistory(int cur_step, double cur_time, bool startup_flag)
{
}




void SNCoreSystem::parseParameters( ParmParse&         a_ppgksys )
{
   // This determines the amount of diagnositic output generated
   a_ppgksys.query( "verbosity", m_verbosity );
   CH_assert( m_verbosity >= 0 );

   // Should we make hdf files for density?
   a_ppgksys.query("hdf_density",m_hdf_density);
   // Should we make hdf files for f versus vparallel and poloidal angle?
   a_ppgksys.query("hdf_vpartheta",m_hdf_vpartheta);
   // Should we make hdf files for f versus radius and poloidal angle at a specified vpar, mu?
   a_ppgksys.query("hdf_frtheta",m_hdf_frtheta);
   // Should we make hdf files for pressure?
   a_ppgksys.query("hdf_pressure",m_hdf_pressure);
   // Should we make hdf files for vparallel-mu?
   a_ppgksys.query("hdf_vparmu",m_hdf_vparmu);
   // At what fixed phase space indices should I plot?  (Indices plotted against in a given plot
   //   are ignored.  Specify in 5D; toroidal index ignored in 4D and set to zero in arguments
   //   of hdf write methods.
#if 0
   m_fixed_plotindices.resize( 5 );
   for (int i=0; i<5; ++i) m_fixed_plotindices[i]=0;
   a_ppgksys.queryarr("fixed_plot_indices",m_fixed_plotindices,0,5);
   for (int block=0; block<a_config_blocks.size(); ++block) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         CH_assert( a_config_blocks[block].loMappedIndex()[dir] <= m_fixed_plotindices[dir] );
         CH_assert( a_config_blocks[block].hiMappedIndex()[dir] >  m_fixed_plotindices[dir] );
      }
   }
#endif

   bool enforce_positivity(false);
   if (a_ppgksys.contains("enforce_positivity")) {
     a_ppgksys.get("enforce_positivity", enforce_positivity);
   }

   if (enforce_positivity) {

     std::string ptype("stage");
     m_enforce_stage_positivity = true;
     if (a_ppgksys.contains("enforce_positivity_type")) {
       a_ppgksys.get("enforce_positivity_type", ptype);
       if (ptype=="step") {
         m_enforce_step_positivity = true;
         m_enforce_stage_positivity = false;
       }
       else if (ptype!="stage") {
         MayDay::Error("Invalid positivity enforcement type");
       }
     }

     int n_iter(5);
     if (a_ppgksys.contains("max_positivity_iter")) {
       a_ppgksys.get("max_positivity_iter", n_iter);
     }

     bool verbose(false);
     if (a_ppgksys.contains("positivity_verbose_output")) {
       a_ppgksys.get("positivity_verbose_output", verbose);
     }

     int width(2);
     if (m_enforce_step_positivity) width++;
     IntVect halo( width*IntVect::Unit );
     m_positivity_post_processor.define( halo, n_iter, verbose );
   }

   if (m_verbosity) {
      printParameters();
   }
}

void SNCoreSystem::getSNCoreParams( SNCoreParamsStruct& a_params ) const
{
   ParmParse a_singlenull( "sncore" );

   a_params.numcells_core.resize( 2 );
   a_singlenull.queryarr( "numcells.core",a_params.numcells_core,0,2 );

   a_params.decomp_core_configuration.resize( 2 );
   a_params.decomp_core_velocity.resize( 2 );
   a_params.decomp_core_phase.resize( 4 );
   a_singlenull.queryarr( "decomp.core.configuration",a_params.decomp_core_configuration,0,2 );
   a_singlenull.queryarr( "decomp.core.velocity",a_params.decomp_core_velocity,0,2 );
   a_singlenull.queryarr( "decomp.core.phase",a_params.decomp_core_phase,0,4 );

   if (procID()==0) {
      cout << endl << "SNCore Parameters" << endl << endl;
      int i;
      cout << "numcells.core = ";
      for (int i=0; i<2; i++) cout << a_params.numcells_core[i] << " ";
      cout << endl;

      cout << "decomp.core.configuration = ";
      for (int i=0; i<2; i++) cout << a_params.decomp_core_configuration[i] << " ";
      cout << endl;
      cout << "decomp.core.velocity = ";
      for (int i=0; i<2; i++) cout << a_params.decomp_core_velocity[i] << " ";
      cout << endl;
      cout << "decomp.core.phase = ";
      for (int i=0; i<4; i++) cout << a_params.decomp_core_phase[i] << " ";
      cout << endl;
   }
}



void SNCoreSystem::printDiagnostics()
{
}



void SNCoreSystem::postStageAdvance( LevelData<FArrayBox>& a_soln )
{
   if (m_enforce_stage_positivity) {
      enforcePositivity( a_soln );
   }
}


static bool dir_created = false;

void SNCoreSystem::computeError(const double                a_time,
                            const int                   a_step,
                            const LevelData<FArrayBox>& a_solution,
                            double&                     a_l1_error,
                            double&                     a_l2_error,
                            double&                     a_max_error)
{
  const DisjointBoxLayout& grids = m_mag_geom->grids();

  LevelData<FArrayBox> exact(m_mag_geom->grids(), 1, IntVect::Unit);

  double amp, R0, delR, Z0, delZ;

  ParmParse ppic("IBC");

  if (ppic.contains("amp") &&
      ppic.contains("R_0") && ppic.contains("delR") &&
      ppic.contains("Z_0") && ppic.contains("delZ") ) {

    ppic.query("amp", amp);

    ppic.query("R_0",R0);
    ppic.query("delR",delR);

    ppic.query("Z_0",Z0);
    ppic.query("delZ",delZ);
  }
  else {
    MayDay::Error("GKSystem::computeError(): Invalid input");
  }

  RealVect center(R0, Z0);
  RealVect width(delR, delZ);

  ParmParse ppadv("Advect");
  Vector<Real> velocity_array;
  ppadv.queryarr("velocity", velocity_array, 0, 2);

  RealVect velocity(velocity_array);

  center += a_time * velocity;

  MagCoordSys* mag_coord_sys = static_cast<MagCoordSys*>(m_mag_geom->coordSysPtr());

  LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
  mag_coord_sys->getCellVolumes(volume);

  // Iterate over patches
  DataIterator dit = a_solution.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    int block_number = mag_coord_sys->whichBlock(grids[dit]);
    const MagBlockCoordSys* block_coord_sys = static_cast<const MagBlockCoordSys*>(mag_coord_sys->getCoordSys(block_number));

    Box box(grids[dit]);
    box.grow(1);
    FArrayBox cell_center_coords( box, SpaceDim );
    block_coord_sys->getCellCenteredRealCoords( cell_center_coords );


    FORT_LOC_SPATIAL_IBC(CHF_BOX(box),
                         CHF_CONST_REAL(amp),
                         CHF_CONST_REALVECT(center),
                         CHF_CONST_REALVECT(width),
                         CHF_CONST_FRA(cell_center_coords),
                         CHF_FRA1(exact[dit],0));

    fourthOrderAverageCell(exact[dit]);
  }

  for (dit.begin(); dit.ok(); ++dit) {
    exact[dit] -= a_solution[dit];
    exact[dit].abs();
  }

#if 1
  if (a_step%20 == 0) {

    //  Plot max norm error

    if(procID()==0 && !dir_created) {
      mkdir("error_plots", 0777);
      dir_created = true;
    }

    char iter_str[100];
    sprintf(iter_str, "error_plots/error.%04d.", a_step);
    stringstream filename;
    filename << iter_str;

    Box plot_box = m_domain.domainBox();
    WriteMappedUGHDF5(filename.str().c_str(), grids, exact, *mag_coord_sys, plot_box);
  }
#endif

  // Compute max norm
  double local_max = 0.;
  for (dit.begin(); dit.ok(); ++dit) {
    double box_max = exact[dit].max(grids[dit]);
    if (box_max > local_max) local_max = box_max;
  }

  MPI_Allreduce(&local_max, &a_max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  // Compute L1 norm
  double local_sum1 = 0.;
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox tmp(grids[dit],1);
    tmp.copy(exact[dit]);
    tmp *= volume[dit];
    local_sum1 += tmp.sum(grids[dit],0,1);
  }

  MPI_Allreduce(&local_sum1, &a_l1_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // Compute L2 norm
  double local_sum2 = 0.;
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox tmp(grids[dit],1);
    tmp.copy(exact[dit]);
    tmp *= exact[dit];
    tmp *= volume[dit];
    local_sum2 += tmp.sum(grids[dit],0,1);
  }

  MPI_Allreduce(&local_sum2, &a_l2_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  a_l2_error = sqrt(a_l2_error);

}


#include "NamespaceFooter.H"
