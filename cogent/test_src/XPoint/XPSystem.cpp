#include "XPSystem.H"
#include "XPSystemF_F.H"


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
#include "XPointBlockCoordSys.H"
#include "XPointCoordSys.H"
#include "newMappedGridIO.H"
#include "FourthOrderUtil.H"
#include "Directions.H"
#include "inspect.H"
#include "MBStencilIterator.H"

#define NORADDECOMP
#undef NORADDECOMP


#include "NamespaceHeader.H"

#ifdef MULTI_BLOCK
MultiBlockLevelExchangeAverage * mblexPtr = NULL;
#endif


XPSystem::XPSystem( ParmParse& a_pp )
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
   Vector<CFGBlock> blocks;

   // Read the input
   parseParameters( m_ppgksys, blocks );

   /*
     Create the configuration space (magnetic geometry) domains,
     data layout and coordinate system.
   */

   Vector<ProblemDomain> cs_domains;
   createConfigurationSpace( blocks, cs_domains );

   m_operator = new Advect( ParmParse(Advect::pp_name) );

   m_solution.define( m_mag_geom->grids(), 1, IntVect::Zero);
   m_initial_condition.define( m_mag_geom->grids(), 1, 4*IntVect::Unit);

   ParmParse pp_ibc("IBC");
   m_ibc = new LocSpatialIBC( pp_ibc, *m_mag_geom, 0);
}



void
XPSystem::initialize(const int cur_step)
{
   // Set the initial conditions if this is the first step
   if (cur_step == 0) {
      setInitialConditions( m_solution );
   }

   // Set the boundary data.  m_solution is not modified by this call (since
   // it now contains the initial condition or the solution obtain from a restart),
   // but simply provides a template for the boundary condition object.
   setBoundaryData( m_solution );

   // Clone the solution prototype to hold old mapped and real data.  The cloning
   // operation includes a copy of the current solution to the newly created
   // vectors.
   defineSolnData(m_solution_mapped, m_solution);
   defineSolnData(m_solution_mapped_old, m_solution);

   divideJ(m_solution, m_solution);

   // Create the integrator
   m_integrator = new RK4Integrator<LevelData<FArrayBox> >(*this, m_solution_mapped);
}



XPSystem::~XPSystem()
{
   delete m_mag_geom;
   delete m_mag_geom_coords;
   delete m_operator;
   delete m_integrator;
}



void XPSystem::createBlocks( Vector<CFGBlock>& a_blocks ) const
{
   ParmParse xppp( "xpsystem" );

   Vector<int> block_size(CFG_DIM);
   if (xppp.contains("block_size")) {
      xppp.getarr( "block_size", block_size, 0, CFG_DIM );
   }
   else {
      MayDay::Error("XPSystem::createBlocks(): xpsystem.block_size not found in input file");
   }

   Vector<int> block_decomp(CFG_DIM);
   if (xppp.contains("block_decomp")) {
      xppp.getarr( "block_decomp", block_decomp, 0, CFG_DIM );
   }
   else {
      MayDay::Error("XPSystem::createBlocks(): xpsystem.block_decomp not found in input file");
   }

   for (int block=0; block<8; ++block) {

      IntVect lo_mapped_index, hi_mapped_index, is_periodic, cfg_decomp;
      lo_mapped_index[0] = 2 * block * block_size[0];
      lo_mapped_index[1] = 0;
      hi_mapped_index[0] = lo_mapped_index[0] + block_size[0] - 1;
      hi_mapped_index[1] = block_size[1] - 1;

      for (int dir=0; dir<CFG_DIM; ++dir) {
         is_periodic[dir] = 0;
         cfg_decomp[dir] = block_decomp[dir];
      }

      CFGBlock block(lo_mapped_index, hi_mapped_index, is_periodic, cfg_decomp);

      a_blocks.push_back(block);

   }
}



void XPSystem::createConfigurationSpace(const Vector<CFGBlock>&        a_blocks,
                                        Vector<ProblemDomain>& a_domains)
{
   int num_blocks = a_blocks.size();

   Vector<Box> boxes;
   for (int block=0; block<num_blocks; ++block) {

      const CFGBlock& cfg_block = a_blocks[block];

      const ProblemDomain domain = cfg_block.getConfigurationDomain();
      a_domains.push_back( domain );

      const IntVect decomp = cfg_block.cfgDecomp();

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
      MayDay::Error( "XPSystem::createConfigurationSpace has only been implemented for 2D" );
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

   DisjointBoxLayout grids( boxes, procMap, m_domain );
   grids.close();

#ifdef CH_MPI
   if (procID()==0) {
#endif
     if (m_verbosity>0) {
       for (int n=0; n<boxes.size(); n++) {
         const Box& local_box = boxes[n];
           cout << "   Configuration space box " << local_box << " is assigned to process " << n << endl;
       }
     }
#ifdef CH_MPI
   }
#endif

   stringstream prefix("gksystem.magnetic_geometry_mapping.", ios_base::out|ios_base::ate);

   Vector<MagBlockCoordSys *> coord_vec;

   prefix << string(XPointBlockCoordSys::pp_name);

   ParmParse XPointPP(prefix.str().c_str());

   for ( int block_number = 0; block_number < a_domains.size(); ++block_number ) {
      MagBlockCoordSys* geom
         = new XPointBlockCoordSys( XPointPP, a_domains[block_number], block_number );
      coord_vec.push_back(geom);
   }

   // Construct the multiblock magnetic coordinate system
   XPointCoordSys* coords = new XPointCoordSys();
   coords->define(coord_vec);

   m_mag_geom_coords = coords;

   // Construct the magnetic geometry
   int ghosts = 0;
   for (int n=0; n<CFG_DIM; ++n) {
      if (m_ghostVect[n] > ghosts) ghosts = m_ghostVect[n];
   }

   if (procID()==0) cout << "Constructing magnetic geometry" << endl;
   m_mag_geom = new MultiBlockLevelGeom(m_mag_geom_coords, grids, ghosts);
   if (procID()==0) cout << "Done constructing magnetic geometry" << endl;

#if 0
   // For testing of metrics
   IntVect geom_data_ghosts = 4*IntVect::Unit;
   LevelData<FArrayBox> geom_data(grids, 6, geom_data_ghosts);

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
      //      inspect(this_geom_data);
   }

   Box plot_box = m_domain.domainBox();
   plot_box.grow(geom_data_ghosts);

   WriteMappedUGHDF5("geom_data", grids, geom_data, *m_mag_geom_coords, plot_box);
   exit(1);
#endif

#if 0
   // For testing of multiblock exchange
   MultiBlockLevelExchangeAverage* mblexPtr = new MultiBlockLevelExchangeAverage();
   int numGhost = 4;
   int spaceOrder = 4;   // = 2*L with L as in Peter's notes
   mblexPtr->define(m_mag_geom, numGhost, spaceOrder);

   LevelData<FArrayBox> tmp_data(grids, 1, numGhost*IntVect::Unit);

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = m_mag_geom_coords->whichBlock(grids[dit]);
      FArrayBox& this_tmp_data = tmp_data[dit];
#if 0
      if (block_number == 0) {
         this_tmp_data.setVal(1.);
      }
      else {
         this_tmp_data.setVal(0.);
      }
#else
      const MagBlockCoordSys* block_coord_sys
         = dynamic_cast<const MagBlockCoordSys*>(m_mag_geom_coords->getCoordSys(block_number));

      this_tmp_data.setVal(0.);

      RealVect dx = block_coord_sys->dx();
      RealVect offset = 0.5*RealVect::Unit;
      offset *= dx;

      BoxIterator bit(grids[dit]);
      for (bit.begin();bit.ok();++bit) {
         IntVect iv = bit();
         RealVect xi = dx*iv + offset;
         RealVect X = block_coord_sys->realCoord(xi);
         this_tmp_data(iv,0) = X[RADIAL_DIR];
      }
#endif
   }

   //   tmp_data.exchange();
   mblexPtr->interpGhosts(tmp_data);

   Box plot_box(m_domain.domainBox());
   plot_box.grow(tmp_data.ghostVect());

#if 0
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
         RealVect xi = dx*iv + offset;
         RealVect X = block_coord_sys->realCoord(xi);
         double val = this_tmp_data(iv,0) - X[RADIAL_DIR];
         this_tmp_data(iv,0) = val;
      }
   }
#endif


   WriteMappedUGHDF5("tmp_data", grids, tmp_data, *m_mag_geom_coords, plot_box);

   exit(1);
#endif
}


void XPSystem::defineRHSData( LevelData<FArrayBox>&       a_rhs,
                              const LevelData<FArrayBox>& a_prototype )
{
   if ( a_rhs.isDefined() ) {
      MayDay::Error("XPSystem::defineRHSData(): a_rhs is already defined");
   }

   a_rhs.define( a_prototype.disjointBoxLayout(),
                 a_prototype.nComp(),
                 a_prototype.ghostVect() );

   DataIterator dit = a_prototype.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_rhs[dit].copy(a_prototype[dit]);
   }
}



void XPSystem::defineSolnData( LevelData<FArrayBox>&       a_soln,
                               const LevelData<FArrayBox>& a_prototype )
{
   if ( a_soln.isDefined() ) {
      MayDay::Error("XPSystem::defineSolnData(): a_soln is already defined");
   }

   a_soln.define( a_prototype.disjointBoxLayout(),
                  a_prototype.nComp(),
                  a_prototype.ghostVect() );

   DataIterator dit = a_prototype.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_soln[dit].copy(a_prototype[dit]);
   }
}



bool XPSystem::validSolnData( const LevelData<FArrayBox>& a_soln,
                              const LevelData<FArrayBox>& a_protoSoln )
{
   return (a_protoSoln.disjointBoxLayout() == a_soln.disjointBoxLayout()) &&
          (a_protoSoln.nComp()             == a_soln.nComp())             &&
          (a_protoSoln.ghostVect()         == a_soln.ghostVect());
}



bool XPSystem::validRHSData( const LevelData<FArrayBox>& a_rhs,
                             const LevelData<FArrayBox>& a_protoRHS )
{
   return (a_protoRHS.disjointBoxLayout() == a_rhs.disjointBoxLayout()) &&
          (a_protoRHS.nComp()             == a_rhs.nComp())             &&
          (a_protoRHS.ghostVect()         == a_rhs.ghostVect());
}



void XPSystem::zeroSolnData( LevelData<FArrayBox>& a_soln )
{
   DataIterator dit = a_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_soln[dit].setVal(0.);
   }
}



void XPSystem::addSolnData( LevelData<FArrayBox>&       a_soln,
                            const LevelData<FArrayBox>& a_increment,
                            const Real                  a_scale )
{
   CH_assert( validSolnData(a_soln, a_increment) );

   DataIterator dit = a_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_soln[dit].plus(a_increment[dit], a_scale);
   }
}



void XPSystem::copySolnData( LevelData<FArrayBox>&       a_dstSoln,
                             const LevelData<FArrayBox>& a_srcSoln )
{
   CH_assert( validSolnData(a_dstSoln, a_srcSoln) );

   DataIterator dit = a_srcSoln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_dstSoln[dit].copy(a_srcSoln[dit]);
   }
}



void XPSystem::evalRHS( LevelData<FArrayBox>&       a_rhs,
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
   //   CH_assert( a_soln.ghostVect() == m_ghostVect );
   const DisjointBoxLayout& grids = a_soln.disjointBoxLayout();

   LevelData<FArrayBox> tmp_soln(grids, a_soln.nComp(), m_ghostVect);
   DataIterator dit = tmp_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     tmp_soln[dit].copy(a_soln[dit], grids[dit]);
   }

   divideJ(tmp_soln, tmp_soln);

#ifdef MULTI_BLOCK
   mblexPtr->interpGhosts(tmp_soln);
#endif

   m_ibc->ghostCellBC(tmp_soln, *m_mag_geom, 0., a_time);

   m_operator->evalRHS( *m_mag_geom, a_rhs, tmp_soln, a_time );
}



void
XPSystem::getConditionNumbers(const MultiBlockLevelExchangeAverage* a_mblexPtr,
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

   // const LayoutData< IVSFAB< Vector<IntVect>* >* >& stencils = a_mblexPtr->allStencils();
   // const LayoutData< IVSFAB< Vector<Real>* >* >& weights = a_mblexPtr->allWeights();
  const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& stencils =
    a_mblexPtr->stencils();

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

#if 1
       if (block_number == 0) {

         for (dit.begin(); dit.ok(); ++dit) {
           int block_number2 = coords->whichBlock(grids[dit]);
           if (block_number == block_number2) {

             const IVSFAB<MBStencil>& this_stencil_fab = *stencils[dit];
             // const IVSFAB< Vector<IntVect>* >& this_ivsfab = *(stencils[dit]);
             // const IVSFAB< Vector<Real>* >& this_fab = *(weights[dit]);

             //      IntVect iv = max_condition_number_iv[block_number];
             IntVect iv = IntVect(0,-2);

             // const Vector<IntVect>& ivect = *(this_ivsfab(iv,0));
             // const Vector<Real>& rvect = *(this_fab(iv,0));
             MBStencilIterator sit(this_stencil_fab(iv, 0));
             for (sit.begin(); sit.ok(); ++sit)
               {
                 const MBStencilElement& stencilElement = sit();
                 cout << stencilElement.cell()
                      << " " << stencilElement.weight() << endl;
               }
           }
         }
       }
#endif

#if CH_MPI
     }
#endif
   }

   Box plot_box = grids.physDomain().domainBox();
   plot_box.grow(log_plot.ghostVect());

   WriteMappedUGHDF5("condition_number_log", grids, log_plot, *coords, plot_box);
}




Real XPSystem::stableDt(const int a_step_number)
{
   return m_operator->computeDt(*m_mag_geom, m_solution_mapped);
}



void XPSystem::enforcePositivity( LevelData<FArrayBox>& a_soln )
{
   m_positivity_post_processor.enforce( a_soln, 1.0 );
}



Real XPSystem::advance( const Real a_cur_time,
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
   divideJ(m_solution_mapped, m_solution);
   return new_time;
}



void XPSystem::setInitialConditions( LevelData<FArrayBox>& a_soln )
{
#if 0
   const DisjointBoxLayout& grids = a_soln.disjointBoxLayout();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = m_mag_geom_coords->whichBlock(grids[dit]);
      const MagBlockCoordSys* block_coord_sys
         = dynamic_cast<const MagBlockCoordSys*>(m_mag_geom_coords->getCoordSys(block_number));

      RealVect dx = block_coord_sys->dx();

      FORT_BLOB_INITIAL_CONDITION(CHF_BOX(m_initial_condition[dit].box()),
                                  CHF_CONST_REALVECT(dx),
                                  CHF_FRA1(m_initial_condition[dit],0));
      a_soln[dit].copy(m_initial_condition[dit]);
   }
#else
   m_ibc->initialize(a_soln, *m_mag_geom, 0., 0.);
#endif
}



void XPSystem::setBoundaryData( LevelData<FArrayBox>& a_soln )
{
   m_ibc->setBoundaryData(a_soln, *m_mag_geom, 0., 0.);
}



void XPSystem::multJ( const LevelData<FArrayBox>& a_soln_physical,
                        LevelData<FArrayBox>&     a_soln_mapped)
{
   // Need to grow the soln factor since we need transverse gradients for the fourth-order product
   const DisjointBoxLayout& grids = a_soln_physical.disjointBoxLayout();
   IntVect grown_vect = a_soln_physical.ghostVect() + IntVect::Unit;
   LevelData<FArrayBox> grown_soln(grids, 1, grown_vect);

   LevelData<FArrayBox> J(grids, 1, grown_vect);
   m_mag_geom_coords->getJ(J);

   DataIterator dit = grown_soln.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const FArrayBox& this_soln = a_soln_physical[dit];
      FArrayBox& this_grown_soln = grown_soln[dit];
      Box src_box = this_soln.box();

      this_grown_soln.copy(this_soln);
      for (int dir=0; dir<SpaceDim; ++dir) {
         secondOrderTransExtrap(this_grown_soln, dir, src_box);
      }
   }
   grown_soln.exchange();

   // Compute fourth-order product
   fourthOrderCellProd(a_soln_mapped, grown_soln, J);
}



void XPSystem::divideJ( const LevelData<FArrayBox>& a_soln_mapped,
                        LevelData<FArrayBox>&       a_soln_physical)
{
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

   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = m_mag_geom_coords->whichBlock(grids[dit]);
      const XPointBlockCoordSys* coord_sys
         = dynamic_cast<const XPointBlockCoordSys*>(m_mag_geom_coords->getCoordSys(block_number));

      cellFGToCellF(a_soln_physical[dit], grown_soln[dit], J[dit], grids[dit], coord_sys->domain(), true);
   }

   a_soln_physical.exchange();
}



void XPSystem::writePlotFile(const char *prefix, const int cur_step )
{
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

   WriteMappedUGHDF5(filename1.str().c_str(), m_solution.disjointBoxLayout(), m_solution, *m_mag_geom_coords, m_domain.domainBox());
}



void XPSystem::writeCheckpointFile( HDF5Handle&  handle,
                                    const int    cur_step,
                                    const double cur_time,
                                    const double cur_dt )
{
}



void XPSystem::readCheckpointFile( HDF5Handle& handle, int* cur_step, double* cur_time, double* cur_dt )
{
}


inline void XPSystem::printParameters(const Vector<CFGBlock>& a_blocks)
{
   if (procID() == 0 && m_verbosity) {

      for (int n=0; n<a_blocks.size(); ++n) {
         cout << endl << "Block " << n << " data:" << endl;
         cout << "   lo_mapped_index      = " << a_blocks[n].loMappedIndex() << endl;
         cout << "   hi_mapped_index      = " << a_blocks[n].hiMappedIndex() << endl;
         cout << "   is_periodic          = " << a_blocks[n].isPeriodic() << endl;
         cout << "   configuration_decomp = " << a_blocks[n].cfgDecomp() << endl;
      }
      cout << endl;
   }
}



void XPSystem::writeFieldHistory(int cur_step, double cur_time, bool startup_flag)
{
}




void XPSystem::parseParameters( ParmParse&         a_ppgksys,
                                Vector<CFGBlock>&  a_config_blocks)
{
   // This determines the amount of diagnositic output generated
   a_ppgksys.query( "verbosity", m_verbosity );
   CH_assert( m_verbosity >= 0 );

   createBlocks(a_config_blocks);

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

   if (m_verbosity) {
      printParameters(a_config_blocks);
   }
}


void XPSystem::printDiagnostics()
{
}


void XPSystem::postStageAdvance( LevelData<FArrayBox>& a_soln )
{
   if (m_enforce_stage_positivity) {
      enforcePositivity( a_soln );
   }
}


#include "NamespaceFooter.H"
