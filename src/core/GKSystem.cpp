#include "GKSystem.H"
#include "CONSTANTS.H"
#include "FORT_PROTO.H"
#include "inspect.H"

#include <fstream>
#include <sstream>

#include "MillerPhaseCoordSys.H"
#include "SlabPhaseCoordSys.H"
#include "SingleNullPhaseCoordSys.H"
#include "SNCorePhaseCoordSys.H"

#include "LocalizedF_F.H"

#if 1  // warning, OS dependencies, will not work on all platforms
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "CH_HDF5.H"
#include "HDF5Portable.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "MagBlockCoordSys.H"
#include "MillerBlockCoordSys.H"
#include "MillerCoordSys.H"
#include "SlabBlockCoordSys.H"
#include "SlabCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "SNCoreBlockCoordSys.H"
#include "newMappedGridIO.H"
#include "inspect.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM 3
#include "AMRIO.H"
#include "Slicing.H.transdim"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "EdgeToCell.H"
#include "FourthOrderUtil.H"
#include "Directions.H"

#include "inspect.H"

#include "NamespaceHeader.H"

using namespace CH_MultiDim;

static const bool SKIP_SOLN_DATA = false;
static const IntVect NO_GHOST_CELLS = IntVect::Zero;

GKSystem::GKSystem( ParmParse& a_pp )
   :
     m_using_electrons(false),
     m_enforce_stage_positivity(false),
     m_enforce_step_positivity(false),
     m_max_grid_size(0),
     m_initial_conditions(NULL),
     m_boundary_conditions(NULL),
     m_gk_ops(NULL),
     m_hdf_potential(false),
     m_hdf_efield(false),
     m_hdf_dfn(false),
     m_hdf_deltaF(false),
     m_hdf_dfn_at_mu(false),
     m_hdf_vpartheta(false),
     m_hdf_rtheta(false),
     m_hdf_density(false),
     m_hdf_total_density(false),
     m_hdf_pressure(false),
     m_hdf_temperature(false),
     m_hdf_fourthMoment(false),
     m_hdf_ParallelMomentum(false),
     m_hdf_PoloidalMomentum(false),
     m_hdf_ParticleFlux(false),
     m_hdf_HeatFlux(false),
     m_hdf_vparmu(false),
     m_hdf_vlasov_divergence(false),
     m_ghostVect(4*IntVect::Unit),
     m_state( GKState(m_ghostVect) ),
     m_verbosity(0),
     ti_class("rk"),
     ti_method("4")
{
   ParmParse ppgksys("gksystem");

   m_units = new GKUnits( ppgksys );
   if (procID()==0) {
      m_units->print(cout);
   }
   
   parseParameters( ppgksys );

   createConfigurationSpace();

   createVelocitySpace();

   createPhaseSpace( ppgksys );

   createSpecies();

   m_state.define( m_kinetic_species, m_phase_geom );

   const Real BASE_DT( 1.0 );
   if (ti_class == "rk") m_integrator = new TiRK<GKState, GKRHSData, GKOps>;
   else MayDay::Error("Unrecognized input for ti_class.");
   m_integrator->define(ti_method, m_state, BASE_DT);
   m_gk_ops = &(m_integrator->getOperators());
   m_integrator->setTimeStep(BASE_DT);
   
   if ( m_using_electrons && m_gk_ops->usingBoltzmannElectrons() ) {
      MayDay::Error( "GKSystem::createSpecies():  Electrons input as both kinetic and Boltzmann" );
   }
   else {
      m_using_electrons |= m_gk_ops->usingBoltzmannElectrons();
   }

   bool fixed_efield( m_gk_ops->fixedEField() );
   if (fixed_efield && m_using_electrons)
      MayDay::Warning( "Using electrons with fixed E field" );
      
   if (!fixed_efield && !m_using_electrons)
      MayDay::Warning( "Not using electrons with dynamic E field" );
   
   setupFieldHistories();
}



void
GKSystem::initialize(const int a_cur_step)
{
   if ( a_cur_step == 0 ) {
      // If this is the first step, then set the initial conditions for the
      // full GK system (i.e., distribution functions and potential)
      m_gk_ops->applyInitialConditions( m_kinetic_species, 0.0 );
   }
   else if ( m_gk_ops->fixedEField() ) {
      // If this is a restart and the fixed_efield option is true, just
      // set the potential to its initial condition.  If the fixed_efield
      // option is false, the potential will be set in the below call
      // to m_gk_ops->initialize(), which will compute both the potential
      // and associated field.
      m_gk_ops->initializePotential( 0.0 );
   }

   // Clone the species prototype to hold old mapped and real data.  The cloning
   // operation includes a copy of the current solution to the newly created
   // vectors.

   // Create the integrator
   KineticSpeciesPtrVect& state_data( m_state.data() );
   const int num_kinetic_species( m_kinetic_species.size() );
   for (int species(0); species<num_kinetic_species; species++) {
      state_data[species]->copy( *(m_kinetic_species[species]) );
   }
   m_gk_ops->divideJ( state_data, m_kinetic_species );

   // Initialize the electric field:
   // a.  If the fixed_efield option is true, then the field is calculated
   //     from the potential set in the prior call to m_gk_ops->applyInitialConditions()
   //     or (if restarting) m_gk_ops->initializePotential().
   // b.  If the fixed_efield option is false, then both the potential and
   //     associated field are computed.
   m_gk_ops->initialize( m_kinetic_species, a_cur_step );
   
#ifdef BLOB_ERROR_TEST
   double l1_error, l2_error, max_error;
   computeError(0., a_cur_step, m_kinetic_species, l1_error, l2_error, max_error);
   if (procID()==0) cout << "l1, l2, max error at t = 0. is " << l1_error << " " << l2_error << " " << max_error << endl;
#endif
}



GKSystem::~GKSystem()
{
   delete m_phase_geom;
   delete m_phase_grid;
   delete m_phase_coords;
   delete m_velocity_coords;
   delete m_mag_geom;
   delete m_mag_geom_coords;
   delete m_units;
   delete m_integrator;
}



VEL::ProblemDomain GKSystem::getVelocityDomain() const
{
   VEL::IntVect lo, hi;
   bool isPeriodic[VEL_DIM];

   /* For now, the number of cells in the vparallel direction must be even.
      This is due to the assumed convention used to compute the physical velocity
      directly from the global index and mesh size.  We intend to remove this assumption,
      and the associated limitation tested for here, in the future.
   */
   if ( m_mag_geom_type != "SingleNull" && m_mag_geom_type != "SNCore" ) {
      if ( m_num_cells[VPARALLEL_DIR]%2 != 0 ) {
         MayDay::Warning( "vparallel dimension must be even" );
      }

      int half_vp_dim = m_num_cells[VPARALLEL_DIR] / 2;
      lo = VEL::IntVect(-half_vp_dim,0);
      hi = VEL::IntVect(half_vp_dim-1, m_num_cells[MU_DIR]-1);

      for (int n=0; n<VEL_DIM; n++) {
         isPeriodic[n] = m_is_periodic[n+CFG_DIM];
      }
   }
   else {

      if ( m_num_velocity_cells[0]%2 != 0 ) {
         MayDay::Warning( "vparallel dimension must be even" );
      }

      int half_vp_dim = m_num_velocity_cells[0] / 2;
      lo = VEL::IntVect(-half_vp_dim,0);
      hi = VEL::IntVect(half_vp_dim-1, m_num_velocity_cells[1]-1);

      for (int n=0; n<VEL_DIM; n++) {
         isPeriodic[n] = 0;
      }
   }

   return VEL::ProblemDomain(VEL::Box(lo,hi), isPeriodic);
}



void GKSystem::createConfigurationSpace()
{
  string mag_geom_prefix = string("gksystem.magnetic_geometry_mapping");

  if ( m_mag_geom_type == "SingleNull" ) {

    ParmParse pp_grid("singlenull");

    string prefix = mag_geom_prefix + string(".")
       + string(CFG::SingleNullBlockCoordSys::pp_name);
    ParmParse pp_geom( prefix.c_str() );

    m_mag_geom_coords = new CFG::SingleNullCoordSys(pp_grid, pp_geom);
  }
  else if ( m_mag_geom_type == "SNCore"  ) {

    ParmParse pp_grid("sncore");

    string prefix = mag_geom_prefix + string(".")
       + string(CFG::SNCoreBlockCoordSys::pp_name);
    ParmParse pp_geom( prefix.c_str() );

    m_mag_geom_coords = new CFG::SNCoreCoordSys(pp_grid, pp_geom);
  }
  else if ( m_mag_geom_type == "Miller" ) {

    string prefix = mag_geom_prefix + string(".")
       + string(CFG::MillerBlockCoordSys::pp_name);
    ParmParse pp( prefix.c_str() );

    m_mag_geom_coords
      = new CFG::MillerCoordSys(pp,
                                m_num_cells,
                                m_is_periodic,
                                m_configuration_decomposition);
  }
  else if ( m_mag_geom_type == "Slab" ) {

    string prefix = mag_geom_prefix + string(".")
       + string(CFG::SlabBlockCoordSys::pp_name);
    ParmParse pp( prefix.c_str() );

    m_mag_geom_coords
      = new CFG::SlabCoordSys(pp,
                                m_num_cells,
                                m_is_periodic,
                                m_configuration_decomposition);
  }

  // Construct the phase space DisjointBoxLayout

  CFG::DisjointBoxLayout grids;

  getConfigurationSpaceDisjointBoxLayout( grids );

  // Construct the magnetic geometry

  if (procID()==0) cout << "Constructing magnetic geometry" << endl;

  ParmParse pp_mag_geom( mag_geom_prefix.c_str() );

  int ghosts = 0;
  for (int n=0; n<CFG_DIM; ++n) {
    if (m_ghostVect[n] > ghosts) ghosts = m_ghostVect[n];
  }

  m_mag_geom = new CFG::MagGeom(pp_mag_geom, m_mag_geom_coords, grids, ghosts);

  if (procID()==0) cout << "Done constructing magnetic geometry" << endl;

#if 0
  // For testing of metric data
  m_mag_geom->writeGeometryData(grids);
  exit(1);
#endif
}



void
GKSystem::getConfigurationSpaceDisjointBoxLayout( CFG::DisjointBoxLayout& grids )
{
  Vector<CFG::Box> boxes;
  for (int block=0; block<m_mag_geom_coords->numBlocks(); ++block) {

    const CFG::MagBlockCoordSys* mag_block_coords
      = (CFG::MagBlockCoordSys *)m_mag_geom_coords->getCoordSys(block);

    const CFG::ProblemDomain& domain = mag_block_coords->domain();
    const CFG::Box& domain_box = domain.domainBox();

    const CFG::IntVect decomp = m_mag_geom_coords->getDecomposition(block);

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

    if ( m_mag_geom_type != "SingleNull" && m_mag_geom_type != "SNCore" ) {
      if (m_configuration_decomposition.size() > 0) {
        int nproc_cfg = 1;
        for (int dir=0; dir<CFG_DIM; ++dir) {
          nproc_cfg *= m_configuration_decomposition[dir];
        }
        //      if (nproc_cfg > nproc) {
        //        MayDay::Error( "Requested configuration space decomposition exceeds available processors" );
        //      }
        for (int dir=0; dir<CFG_DIM; ++dir) {
          if (domain_box.size(dir)%m_configuration_decomposition[dir] != 0) {
            stringstream msg("Decomposition in configuration direction ", ios_base::out|ios_base::ate);
            msg << dir << " does not evenly divide domain dimension";
            MayDay::Error( msg.str().c_str() );
          }
          else {
            n_loc[dir] = domain_box.size(dir) / m_configuration_decomposition[dir];
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
            if (domain_box.size(0)%p == 0 && domain_box.size(1)%q == 0) {
              n_loc[RADIAL_DIR] = domain_box.size(0) / p;
              n_loc[POLOIDAL_DIR] = domain_box.size(1) / q;
            }
          }
        }
      }
    }
    else {

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
    }

    if (n_loc[RADIAL_DIR] != 0 && n_loc[POLOIDAL_DIR] != 0) {

      CFG::IntVect box_size(n_loc[RADIAL_DIR],n_loc[POLOIDAL_DIR]);
      CFG::Box patch(domain_box.smallEnd(), domain_box.smallEnd() + box_size-1);
      CFG::Box skeleton(CFG::IntVect::Zero, CFG::IntVect(domain_box.size(0)/n_loc[RADIAL_DIR]-1,domain_box.size(1)/n_loc[POLOIDAL_DIR]-1));
      CFG::BoxIterator bit(skeleton);
      for (bit.begin();bit.ok();++bit) {
        CFG::Box thisBox = patch + bit()*box_size;
        boxes.push_back(thisBox);
      }
    }
    else {
      MayDay::Error( "Configuration domain box cannot be load balanced" );
    }
#else
    MayDay::Error( "GKSystem::createConfigurationSpace has only been implemented for 2D" );
#endif
  }

  // Make the layout.  This is where boxes are assigned to processes.
  Vector<int> procMap;
  CFG::LoadBalance( procMap, boxes );

  CFG::ProblemDomain prob_domain;

  if ( m_mag_geom_type == "Miller" || m_mag_geom_type == "Slab" ) {
    const CFG::MagBlockCoordSys* mag_block_coords
      = (CFG::MagBlockCoordSys *)m_mag_geom_coords->getCoordSys(0);
    prob_domain = mag_block_coords->domain();
  }
  else if ( m_mag_geom_type == "SingleNull" || m_mag_geom_type == "SNCore" ) {
    CFG::Box bounding_box;
    for (int n=0; n<boxes.size(); n++) {
      bounding_box = minBox(bounding_box, boxes[n]);
    }
    prob_domain = CFG::ProblemDomain(bounding_box);
  }
  else {
    MayDay::Error("Invalid magnetic geometry type");
  }

  grids.define( boxes, procMap, prob_domain );
  grids.close();

#ifdef CH_MPI
  if (procID()==0) {
#endif
    if (m_verbosity>0) {
      for (int n=0; n<boxes.size(); n++) {
        const CFG::Box& local_box = boxes[n];
        cout << "   Configuration space box " << local_box << " is assigned to process " << procMap[n] << endl;
      }
    }
#ifdef CH_MPI
  }
#endif
}



void
GKSystem::createVelocitySpace()
{
   const VEL::ProblemDomain& domain = getVelocityDomain();

   VEL::DisjointBoxLayout grids;
   getVelocitySpaceDisjointBoxLayout( domain, grids );

   ParmParse pppsm( "phase_space_mapping" );

   Real v_parallel_max;
   pppsm.get("v_parallel_max", v_parallel_max);
   // v_parallel runs from -v_parallel_max to + v_parallel_max
   Real dv_parallel = 2*v_parallel_max / domain.size(0);

   Real mu_max;
   pppsm.get("mu_max", mu_max);
   Real dmu = mu_max / domain.size(1);

   const VEL::RealVect dv(dv_parallel,dmu);

   string prefix = string("gksystem.") + string(VEL::VelCoordSys::pp_name);
   ParmParse pp_vel(prefix.c_str());

   m_velocity_coords = new VEL::VelCoordSys(pp_vel, grids, domain, dv);
}



void
GKSystem::getVelocitySpaceDisjointBoxLayout( const VEL::ProblemDomain& domain,
                                             VEL::DisjointBoxLayout& grids )
{
   VEL::Box domain_box = domain.domainBox();

   /*
     Chop up the velocity space domain box over the current number of processors.
     At this point, we insist that the box decomposes uniformly, or an error is thrown.
   */
   int nproc = numProc();
   int n_loc[VEL_DIM];

   for (int dir=0; dir<VEL_DIM; ++dir) {
      n_loc[dir] = 0;
   }

   if ( m_mag_geom_type != "SingleNull" && m_mag_geom_type != "SNCore" ) {

      if (m_velocity_decomposition.size() > 0) {
         int nproc_vel = 1;
         for (int dir=0; dir<VEL_DIM; ++dir) {
            nproc_vel *= m_velocity_decomposition[dir];
         }
         //      if (nproc_vel > nproc) {
         //         MayDay::Error( "Requested velocity space decomposition exceeds available processors" );
         //      }
         for (int dir=0; dir<VEL_DIM; ++dir) {
            if (domain_box.size(dir)%m_velocity_decomposition[dir] != 0) {
               stringstream msg("Decomposition in velocity direction ", ios_base::out|ios_base::ate);
               msg << dir << " does not evenly divide domain dimension";
               MayDay::Error( msg.str().c_str() );
            }
            else {
               n_loc[dir] = domain_box.size(dir) / m_velocity_decomposition[dir];
            }
         }
      }
      else {

         /*
           Look for factors p and q of nproc that evenly divide
           the box velocity dimensions.  If there are
           multiple (p,q) pairs that qualify, use the one for which
           p and q are the closest.
         */
         int q;
         for (q=1; q<=(int)sqrt(1.0*nproc); q++) {
            if (nproc%q == 0) {    // will always true at least for q = 1
               int p = nproc / q;
               if (domain_box.size(0)%p == 0 && domain_box.size(1)%q == 0) {
                  n_loc[0] = domain_box.size(0) / p;
                  n_loc[1] = domain_box.size(1) / q;
               }
            }
         }
      }
   }
   else {

      int nproc_vel = 1;
      for (int dir=0; dir<VEL_DIM; ++dir) {
         nproc_vel *= m_velocity_decomposition[dir];
      }
      //      if (nproc_vel > nproc) {
      //         MayDay::Error( "Requested velocity space decomposition exceeds available processors" );
      //      }
      for (int dir=0; dir<VEL_DIM; ++dir) {
         if (domain_box.size(dir)%m_velocity_decomposition[dir] != 0) {
            stringstream msg("Decomposition in velocity direction ", ios_base::out|ios_base::ate);
            msg << dir << " does not evenly divide domain dimension";
            MayDay::Error( msg.str().c_str() );
         }
         else {
            n_loc[dir] = domain_box.size(dir) / m_velocity_decomposition[dir];
         }
      }
   }

   Vector<VEL::Box> boxes;

   if (n_loc[0] != 0 && n_loc[1] != 0) {

      VEL::IntVect box_size(n_loc[0],n_loc[1]);
      VEL::Box patch(domain_box.smallEnd(), domain_box.smallEnd() + box_size-1);
      VEL::Box skeleton(VEL::IntVect::Zero, VEL::IntVect(domain_box.size(0)/n_loc[0]-1,domain_box.size(1)/n_loc[1]-1));
      VEL::BoxIterator bit(skeleton);
      for (bit.begin();bit.ok();++bit) {
         VEL::Box thisBox = patch + bit()*box_size;
         boxes.push_back(thisBox);
      }
   }
   else {
      MayDay::Error( "Velocity domain box cannot be load balanced" );
   }

   // Make the layout.  This is where boxes are assigned to processes.
   Vector<int> procMap;
   VEL::LoadBalance( procMap, boxes );
   grids.define( boxes, procMap, domain );
   grids.close();

#ifdef CH_MPI
   if (procID()==0) {
#endif
     if (m_verbosity>0) {
       for (int n=0; n<boxes.size(); n++) {
         const CFG::Box& local_box = boxes[n];
         cout << "   Velocity space box " << local_box << " is assigned to process " << procMap[n] << endl;
       }
     }
#ifdef CH_MPI
   }
#endif
}



void
GKSystem::createPhaseSpace( ParmParse& a_ppgksys )
{
  // Build a vector of ProblemDomains, one per configuration space block

  const VEL::ProblemDomain& vel_domain = getVelocityDomain();

  Vector<ProblemDomain> domains;
  for (int block=0; block<m_mag_geom_coords->numBlocks(); ++block) {

    // Construct the phase space block domain from the configuration and velocity domains

    const CFG::MagBlockCoordSys* mag_block_coords
      = (CFG::MagBlockCoordSys *)m_mag_geom_coords->getCoordSys(block);

    const CFG::ProblemDomain& cfg_domain = mag_block_coords->domain();
    const CFG::Box& cfg_box = cfg_domain.domainBox();

    IntVect loVect, hiVect;
    bool isPeriodic[PDIM];

    for (int n=0; n<CFG_DIM; n++) {
      loVect[n] = cfg_box.smallEnd(n);
      hiVect[n] = cfg_box.bigEnd(n);
      isPeriodic[n] = cfg_domain.isPeriodic(n);
    }
    const VEL::Box & vel_box = vel_domain.domainBox();
    for (int n=CFG_DIM; n<PDIM; n++) {
      loVect[n] = vel_box.smallEnd(n-CFG_DIM);
      hiVect[n] = vel_box.bigEnd(n-CFG_DIM);
      isPeriodic[n] = vel_domain.isPeriodic(n-CFG_DIM);
    }

    Box domain_box(loVect,hiVect);
    ProblemDomain domain( domain_box, isPeriodic );
    domains.push_back( domain );
  }

  // Construct the multiblock phase space coordinate system

  if ( m_mag_geom_type == "Miller" ) {
    m_phase_coords = new MillerPhaseCoordSys( a_ppgksys,
                                              *(CFG::MillerCoordSys*)m_mag_geom_coords,
                                              *m_velocity_coords,
                                              domains );
  }
  else if ( m_mag_geom_type == "Slab" ) {
    m_phase_coords = new SlabPhaseCoordSys( a_ppgksys,
                                              *(CFG::SlabCoordSys*)m_mag_geom_coords,
                                              *m_velocity_coords,
                                              domains );
  }
  else if ( m_mag_geom_type == "SingleNull" ) {
    ParmParse pp("singlenull.decomp");
    m_phase_coords = new SingleNullPhaseCoordSys( pp,
                                                  *(CFG::SingleNullCoordSys*)m_mag_geom_coords,
                                                  *m_velocity_coords,
                                                  domains );
  }
  else if ( m_mag_geom_type == "SNCore" ) {
    ParmParse pp("sncore.decomp");
    m_phase_coords = new SNCorePhaseCoordSys( pp,
                                              *(CFG::SNCoreCoordSys*)m_mag_geom_coords,
                                              *m_velocity_coords,
                                              domains );
  }
  else {
    MayDay::Error("Invalid magnetic geometry type");
  }

  Vector<IntVect> decomps;
  for (int block=0; block<m_mag_geom_coords->numBlocks(); ++block) {
    decomps.push_back( m_phase_coords->getDecomposition(block) );
  }

  // Construct the phase space grid

  m_phase_grid = new PhaseGrid(domains, decomps, m_phase_decomposition, m_mag_geom_type);

  if (m_verbosity>0) {
     m_phase_grid->print(m_ghostVect);
  }

  // Construct the phase space geometry
  int ghosts = 0;
  for (int n=0; n<SpaceDim; ++n) {
    if (m_ghostVect[n] > ghosts) ghosts = m_ghostVect[n];
  }

  m_phase_geom =
     RefCountedPtr<PhaseGeom>( new PhaseGeom( a_ppgksys,
                                              m_phase_coords,
                                              *m_phase_grid,
                                              *m_mag_geom,
                                              *m_velocity_coords,
                                              ghosts,
                                              m_units->larmorNumber() ) );

  m_phase_geom.neverDelete();  // workaround for some problem with RefCountedPtr
}


void
GKSystem::createSpecies()
{
   /*
     Create the vector of species model (pointers), and when a kinetic
     species is created, add it to the kinetic species vector
     number species.
   */

   if ( m_verbosity && procID() == 0 ) {
      cout << "Adding species and constructing gyrocenter coordinates..." << endl;
   }

   int ghosts = 0;
   for (int n=0; n<SpaceDim; ++n) {
      if (m_ghostVect[n] > ghosts) ghosts = m_ghostVect[n];
   }

   bool more_kinetic_species = true;

   // run through the loop below for all ion species (but skipping electrons); then go through
   // loop one more time for electrons.  Coding here assumes there is a single electron species.
   int i_off(1);
   bool reading_kinetic_electrons = false;
   int i_electron(0);
   string electronname;

   while ( more_kinetic_species || reading_kinetic_electrons ) {

      bool species_is_complete = true;

      // look for data specifying another kinetic species
      stringstream s;
      if ( !reading_kinetic_electrons ){
        s << "kinetic_species." << m_kinetic_species.size()+i_off;
      }
      else{
        // Re-create string for parsing electron data on the (second) electron pass through this loop
        s << "kinetic_species." << i_electron;
      }
      ParmParse ppspecies( s.str().c_str() );

      string name;
      if (ppspecies.contains("name")) {
         ppspecies.get("name", name);
      }
      else species_is_complete = false;

      // check if we are doing electrons;
      if (!reading_kinetic_electrons){
        if ( (name == "electron") || (name == "electrons") ){
           // recognize either "electron" or "electrons", but must be consistent
           i_electron = m_kinetic_species.size() + 1;
           // increment i_off so remaining ion species are gotten starting on next pass
           i_off = 2;
           m_using_electrons = true;
           electronname = name;
           // end this iteration of the "while (more kinetic_species)" loop
           continue;
        }
      }
      Real mass;
      if ( ppspecies.contains("mass") ) {
         ppspecies.get("mass", mass);
      }
      else species_is_complete = false;

      Real charge;
      if ( (name != "electron") && (name != "electrons") ){
        if (ppspecies.contains("charge")) {
          ppspecies.get("charge", charge);
        }
        else species_is_complete = false;
      }
      else {
        // electrons must have charge -1.0
        charge = -1.0;
      }

      if ( species_is_complete ) {
         if ( procID() == 0 && m_verbosity ) {
            cout << "   " << name << ": kinetic species with mass = " << mass << ", charge = "
                 << charge << endl;
         }

         // Get the species geometry
         PhaseGeom* species_geom = new PhaseGeom(*m_phase_geom, mass, charge);

         // Create the species object
         KineticSpecies* kin_spec = new KineticSpecies( name, mass, charge, *species_geom );
         kin_spec->distributionFunction().define(m_phase_grid->disjointBoxLayout(), 1, IntVect::Zero);

         // Add the new species to the solution vector
         KineticSpeciesPtr species(kin_spec);
         m_species.push_back(species);
         m_kinetic_species.push_back(species);

      }
      else {
         more_kinetic_species = false;
      }

      // If all species have been input and kinetic electrons have been included,
      // then we need to go back and read them.
      if ( !more_kinetic_species && m_using_electrons ) {
         reading_kinetic_electrons = !reading_kinetic_electrons;
      }
   }

   if ( procID() == 0 && m_verbosity ) {
      cout << "Done adding species and constructing coordinates" << endl;
   }
}



void GKSystem::defineRHSData( KineticSpeciesPtrVect&       a_rhs,
                              const KineticSpeciesPtrVect& a_prototype )
{
   a_rhs.resize( a_prototype.size() );
   for (int i=0; i<a_prototype.size(); i++) {
      a_rhs[i] = a_prototype[i]->clone( NO_GHOST_CELLS, SKIP_SOLN_DATA );
   }
}



void GKSystem::defineSolnData( KineticSpeciesPtrVect&       a_soln,
                               const KineticSpeciesPtrVect& a_prototype )
{
   a_soln.resize( a_prototype.size() );
   for (int i=0; i<a_prototype.size(); i++) {
      a_soln[i] = a_prototype[i]->clone( NO_GHOST_CELLS, SKIP_SOLN_DATA );
   }
}



bool GKSystem::validSolnData( const KineticSpeciesPtrVect& a_soln,
                              const KineticSpeciesPtrVect& a_protoSoln )
{
   bool valid = ( a_soln.size() == a_protoSoln.size() );
   if (valid) {
      for (int i=0; i<a_protoSoln.size(); i++) {
         valid &= a_soln[i]->conformsTo( *(a_protoSoln[i]) );
      }
   }
   return valid;
}



bool GKSystem::validRHSData( const KineticSpeciesPtrVect& a_rhs,
                             const KineticSpeciesPtrVect& a_protoSoln )
{
   bool valid = ( a_rhs.size() == a_protoSoln.size() );
   if (valid) {
      for (int i=0; i<a_protoSoln.size(); i++) {
         valid &= a_rhs[i]->conformsTo( *(a_protoSoln[i]), false );
      }
   }
   return valid;
}



void GKSystem::zeroSolnData( KineticSpeciesPtrVect& a_soln )
{
   for (int i=0; i<a_soln.size(); i++) {
      a_soln[i]->zeroData();
   }
}



void GKSystem::addSolnData( KineticSpeciesPtrVect&       a_soln,
                            const KineticSpeciesPtrVect& a_increment,
                            const Real            a_scale )
{
   CH_assert( a_soln.size() == a_increment.size() );
   for (int i=0; i<a_soln.size(); i++) {
      CH_assert( a_increment[i]->conformsTo( *(a_soln[i]), false ) );
      a_soln[i]->addData( *(a_increment[i]), a_scale );
   }
}



void GKSystem::copySolnData( KineticSpeciesPtrVect&       a_dstSoln,
                             const KineticSpeciesPtrVect& a_srcSoln )
{
   CH_assert( a_dstSoln.size() == a_srcSoln.size() );
   for (int i=0; i<a_srcSoln.size(); i++) {
      (a_dstSoln[i])->copy(*(a_srcSoln[i]));
   }
}



Real GKSystem::stableDt_Vlasov(const int a_step_number)
{
   return m_gk_ops->stableDt_Vlasov( m_state.data(), a_step_number );
}



Real GKSystem::stableDt_Collisions(const int a_step_number)
{
   return m_gk_ops->stableDt_Collisions( m_state.data(), a_step_number );
}



Real GKSystem::stableDt(const int a_step_number)
{
   return m_gk_ops->stableDt( m_state.data(), a_step_number );
}



void GKSystem::enforcePositivity( KineticSpeciesPtrVect& a_soln )
{
   for (int n(0); n<a_soln.size(); n++) {
      KineticSpecies* kinetic_species( a_soln[n].operator->() );
      m_positivity_post_processor.enforce( kinetic_species->distributionFunction(),
                                           kinetic_species->maxValue() );
   }
}



Real GKSystem::advance( const Real a_cur_time,
                        const Real a_dt,
                        const int  a_step_number)
{
   m_integrator->setCurrentTime( a_cur_time );

   m_integrator->setTimeStep( a_dt );
   m_integrator->advance( a_cur_time, m_state );

   Real new_time;
   m_integrator->getCurrentTime(&new_time);
   
   if (m_enforce_step_positivity) {
      enforcePositivity( m_state.data() );
   }
   m_gk_ops->divideJ( m_state.data(), m_kinetic_species );


#ifdef BLOB_ERROR_TEST
   double l1_error, l2_error, max_error;
   computeError(new_time, a_step_number, m_kinetic_species, l1_error, l2_error, max_error);
   if (procID()==0) cout << "l1, l2, max error at t = " << new_time << " is " << l1_error << " " << l2_error << " " << max_error << endl;
#endif

   return new_time;
}


double GKSystem::sumDfn( const LevelData<FArrayBox>& dfn )
{
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();

   double local_sum = 0.;
   for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      local_sum += dfn[dit].sum(grids[dit],0,1);
   }

   double sum;
#ifdef CH_MPI
   MPI_Allreduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
   sum = local_sum;
#endif

   return sum;
}


inline
std::string dirPrefix( const std::string& a_prefix,
                       const std::string& a_diag_name )
{
   std::string dir_prefix( a_prefix );
#if 1  // warning, OS dependencies, will not work on all platforms
   std::string iter_str( a_prefix + "_" + a_diag_name + "_plots" );
#ifdef CH_MPI
   if (procID() == 0) {
#endif
      // only works the first time, subsequent failure is normal and expected
      mkdir( iter_str.c_str(), 0777 ); 
#ifdef CH_MPI
   }
#endif
   dir_prefix = std::string( iter_str + "/" );
#endif
   return dir_prefix;
}


inline
std::string plotFileName( const std::string& a_prefix,
                          const std::string& a_diag_name,
                          const int a_cur_step )
{
   std::string dir_prefix( dirPrefix( a_prefix, a_diag_name ) );
   char buffer[100];
   sprintf( buffer, "%04d.", a_cur_step );
   std::string filename( dir_prefix + a_prefix + "." + a_diag_name + buffer );
   return filename;
}


inline
std::string plotFileName( const std::string& a_prefix,
                          const std::string& a_diag_name,
                          const std::string& a_species_name,
                          const int a_cur_step,
                          const int a_species_index )
{
   std::string dir_prefix( dirPrefix( a_prefix, a_diag_name ) );
   char buffer0[10];
   char buffer[20];
   sprintf( buffer0, ".%d.", a_species_index);
   sprintf( buffer, "%04d", a_cur_step );
   std::string filename( dir_prefix + a_prefix + buffer0 + a_species_name + "." + a_diag_name + buffer + "." );
   
   return filename;
}


void GKSystem::writePlotFile(const char *prefix, const int cur_step )
{
   // If the efield and potential are fixed, only consider plotting them at step 0
   if ( !m_gk_ops->fixedEField() || cur_step == 0 ) {
   
      if (m_hdf_potential) {
         std::string filename( plotFileName( prefix,
                                             "potential",
                                             cur_step ) );
         m_gk_ops->plotPotential( filename );
      }

      if (m_hdf_efield) {
         std::string filename = plotFileName( prefix,
                                              "efield",
                                              cur_step );
         m_gk_ops->plotEField( filename );
      }
   }

   // Get indices for plots at various specified points
   int radial_index (m_fixed_plotindices[0]);
   int poloidal_index (m_fixed_plotindices[1]);
#if CFG_DIM==2
   int toroidal_index (0);
#endif
#if CFG_DIM==3
   int toroidal_index (m_fixed_plotindices[2]);
#endif
   int vpar_index (m_fixed_plotindices[3]);
   int mu_index (m_fixed_plotindices[4]);

   for (int species(0); species<m_kinetic_species.size(); species++) {
      const KineticSpecies& soln_species( *(m_kinetic_species[species]) );

      // Distribution function

      if (m_hdf_dfn) {
         std::string filename( plotFileName( prefix,
                                             "dfn",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotDistributionFunction( filename, soln_species );

         // Write out BStarParallel at first step
         if (cur_step==0) {
            std::string filename( plotFileName( prefix,
                                                "bstar_par",
                                                soln_species.name(),
                                                cur_step,
                                                species + 1) );

            m_gk_ops->plotBStarParallel( filename, soln_species );
         }
      }

      // Distribution function with subtracted Maxwellian 

      if (m_hdf_deltaF) {
         std::string filename( plotFileName( prefix,
                                             "deltaF",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotDeltaF( filename, soln_species );
      }

      // Distribution function at a specified mu index

      if (m_hdf_dfn_at_mu) {
         std::string filename( plotFileName( prefix,
                                             "dfn_at_mu",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotDistributionFunctionAtMu( filename, soln_species, mu_index );
      }

      // Distribution function vparallel versus theta at specified configuration space point

      if (m_hdf_vpartheta) {
         std::string filename( plotFileName( prefix,
                                             "vpar_poloidal",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotVParallelTheta( filename, soln_species, radial_index, toroidal_index, mu_index );
      }

      // Distribution function r versus theta at specified velocity space point

      if (m_hdf_rtheta) {
         std::string filename( plotFileName( prefix,
                                             "f_rtheta",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotRTheta( filename, soln_species, vpar_index, mu_index );
      }


      // Distribution function vparallel versus mu at fixed configuration coordinate

      if (m_hdf_vparmu) {
         std::string filename( plotFileName( prefix,
                                             "vpar_mu",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotVParallelMu( filename, soln_species, radial_index, poloidal_index );
      }

      // Charge density

      if (m_hdf_density) {
         std::string filename( plotFileName( prefix,
                                             "density",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotChargeDensity( filename, soln_species );
      }

      // Parallel momentum moment

      if (m_hdf_ParallelMomentum) {
         std::string filename( plotFileName( prefix,
                                             "ParallelMomentum",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotParallelMomentum( filename, soln_species );
      }

      // Poloidal momentum moment

      if (m_hdf_PoloidalMomentum) {
         std::string filename( plotFileName( prefix,
                                             "PoloidalMomentum",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotPoloidalMomentum( filename, soln_species );
      }

      // Pressure

      if (m_hdf_pressure) {
         std::string filename( plotFileName( prefix,
                                             "pressure",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotPressure( filename, soln_species );
      }

      // Temperature

      if (m_hdf_temperature) {
         std::string filename( plotFileName( prefix,
                                             "temperature",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotTemperature( filename, soln_species );
      }

      // Fourth momemnt

      if (m_hdf_fourthMoment) {
         std::string filename( plotFileName( prefix,
                                             "fourthMoment",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotFourthMoment( filename, soln_species );
      }

      // Particle flux moment

      if (m_hdf_ParticleFlux) {
         std::string filename( plotFileName( prefix,
                                             "ParticleFlux",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotParticleFlux( filename, soln_species );
      }

      // Heat flux moment

      if (m_hdf_HeatFlux) {
         std::string filename( plotFileName( prefix,
                                             "HeatFlux",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotHeatFlux( filename, soln_species );
      }

      // Vlasov divergence
      if (m_hdf_vlasov_divergence) {
         std::string filename( plotFileName( prefix,
                                             "vlasov_divergence",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotVlasovDivergence( filename, soln_species );
      }
   }

   // Total charge density

   if (m_hdf_total_density) {
      std::string filename( plotFileName( prefix,
                                          "charge_density",
                                          cur_step ) );

      m_gk_ops->plotChargeDensity( filename, m_kinetic_species );
   }
}

      

void GKSystem::writeCheckpointFile( HDF5Handle&  a_handle,
                                    const int    a_cur_step,
                                    const double a_cur_time,
                                    const double a_cur_dt )
{
   pout() << "writing checkpoint file" << endl;

   MPI_Barrier(MPI_COMM_WORLD);
   if (procID()==0) {
      cout << "Writing checkpoint file" << endl;
   }

   HDF5HeaderData header;
   const int RESTART_VERSION(2); // to distinguish from future versions
   header.m_int ["cur_step"]        = a_cur_step;
   header.m_real["cur_time"]        = a_cur_time;
   header.m_real["cur_dt"]          = a_cur_dt;
   header.m_real["Er_lo"]           = m_gk_ops->getLoRadialField();
   header.m_real["Er_hi"]           = m_gk_ops->getHiRadialField();
   header.m_int ["restart_version"] = RESTART_VERSION;
   header.writeToFile( a_handle );

   if ( m_gk_ops->usingAmpereLaw() ) {
      // save the averaged radial E-field in 4D injected form
      LevelData<FArrayBox> Er_cell_injected;
      m_phase_geom->injectConfigurationToPhase( m_gk_ops->getErAverageCell(), Er_cell_injected);
      a_handle.setGroup( "Er_cell" );
      write(a_handle,Er_cell_injected.boxLayout());
      write(a_handle,Er_cell_injected,"data",Er_cell_injected.ghostVect());
 
      LevelData<FluxBox> Er_face_injected;
      m_phase_geom->injectConfigurationToPhase( m_gk_ops->getErAverageFace(), Er_face_injected);
      a_handle.setGroup( "Er_face" );
      write(a_handle,Er_face_injected.boxLayout());
      write(a_handle,Er_face_injected,"data",Er_face_injected.ghostVect());

      LevelData<FArrayBox> E_tilde_cell_injected;
      m_phase_geom->injectConfigurationToPhase( m_gk_ops->getETildeCell(), E_tilde_cell_injected);
      a_handle.setGroup( "E_tilde_cell" );
      write(a_handle,E_tilde_cell_injected.boxLayout());
      write(a_handle,E_tilde_cell_injected,"data",E_tilde_cell_injected.ghostVect());
 
      LevelData<FluxBox> E_tilde_face_injected;
      m_phase_geom->injectConfigurationToPhase( m_gk_ops->getETildeFace(), E_tilde_face_injected);
      a_handle.setGroup( "E_tilde_face" );
      write(a_handle,E_tilde_face_injected.boxLayout());
      write(a_handle,E_tilde_face_injected,"data",E_tilde_face_injected.ghostVect());
   }
   
   for (int species(0); species<m_kinetic_species.size(); species++) {

      // Get solution distribution function for the current species
      KineticSpecies& soln_species( *((m_state.data())[species]) );
      LevelData<FArrayBox> & soln_dfn = soln_species.distributionFunction();
      char buff[100];
      sprintf( buff, "dfn_%d", species + 1 );
      a_handle.setGroup( buff );
      write( a_handle, soln_dfn.boxLayout() );
      write( a_handle, soln_dfn, "data" );
   }

   MPI_Barrier(MPI_COMM_WORLD);
   if (procID()==0) {
      cout << "Writing history file" << endl;
   }

   m_gk_ops->writeCheckpointFile( a_handle );
}


void GKSystem::readCheckpointFile( HDF5Handle& a_handle,
                                   int&        a_cur_step,
                                   double&     a_cur_time,
                                   double&     a_cur_dt )
{
   //pout() << "reading checkpoint file" << endl;
   HDF5HeaderData header;
   header.readFromFile( a_handle );

   a_cur_step = header.m_int ["cur_step"];
   a_cur_time = header.m_real["cur_time"];
   a_cur_dt   = header.m_real["cur_dt"];
   m_gk_ops->setLoRadialField(header.m_real["Er_lo"]);
   m_gk_ops->setHiRadialField(header.m_real["Er_hi"]);
   const int restart_version = header.m_int ["restart_version"];

   if ( m_gk_ops->usingAmpereLaw() ) {
      a_handle.setGroup("Er_cell");
      CFG::LevelData<CFG::FArrayBox> Er_cell_config(m_mag_geom->grids(), 3, 2*CFG::IntVect::Unit);
      LevelData<FArrayBox> Er_cell_injected;
      m_phase_geom->injectConfigurationToPhase(Er_cell_config, Er_cell_injected);
      read( a_handle, Er_cell_injected, "data", Er_cell_injected.disjointBoxLayout() );
      m_gk_ops->setErAverage( Er_cell_injected );
    
      a_handle.setGroup("Er_face");
      CFG::LevelData<CFG::FluxBox> Er_face_config(m_mag_geom->grids(), 3, 2*CFG::IntVect::Unit);
      LevelData<FluxBox> Er_face_injected;
      m_phase_geom->injectConfigurationToPhase(Er_face_config, Er_face_injected);
      read( a_handle, Er_face_injected, "data", Er_face_injected.disjointBoxLayout() );
      m_gk_ops->setErAverage( Er_face_injected );

      a_handle.setGroup("E_tilde_cell");
      CFG::LevelData<CFG::FArrayBox> E_tilde_cell_config(m_mag_geom->grids(), 3, 2*CFG::IntVect::Unit);
      LevelData<FArrayBox> E_tilde_cell_injected;
      m_phase_geom->injectConfigurationToPhase(E_tilde_cell_config, E_tilde_cell_injected);
      read( a_handle, E_tilde_cell_injected, "data", E_tilde_cell_injected.disjointBoxLayout() );
      m_gk_ops->setETilde( E_tilde_cell_injected );
    
      a_handle.setGroup("E_tilde_face");
      CFG::LevelData<CFG::FluxBox> E_tilde_face_config(m_mag_geom->grids(), 3, 2*CFG::IntVect::Unit);
      LevelData<FluxBox> E_tilde_face_injected;
      m_phase_geom->injectConfigurationToPhase(E_tilde_face_config, E_tilde_face_injected);
      read( a_handle, E_tilde_face_injected, "data", E_tilde_face_injected.disjointBoxLayout() );
      m_gk_ops->setETilde( E_tilde_face_injected );
   }

   for (int species(0); species<m_kinetic_species.size(); species++) {

      // Get solution distribution function for the current species
      KineticSpecies& soln_species( *(m_kinetic_species[species]) );
      LevelData<FArrayBox>& soln_dfn = soln_species.distributionFunction();
      char buff[100];
      sprintf( buff, "dfn_%d", species + 1 );
      a_handle.setGroup( buff );
      read( a_handle, soln_dfn, "data", soln_dfn.disjointBoxLayout(), Interval(0,soln_dfn.nComp()-1), false );
   }

   m_gk_ops->readCheckpointFile( a_handle, a_cur_step );
}


void GKSystem::setupFieldHistories()

{
   ParmParse ppsim( "simulation" );
   m_gk_ops->setupFieldHistories( ppsim );
}


void GKSystem::writeFieldHistory(int cur_step, double cur_time, bool startup_flag)
{
  m_gk_ops->writeFieldHistory( cur_step, cur_time, startup_flag );
}


void GKSystem::printDiagnostics()
{
   pout() << "  Distribution Function Extrema:" << std::endl;
   if (procID()==0)
      cout << "  Distribution Function Extrema:" << std::endl;
   for (int n(0); n<m_kinetic_species.size(); n++) {
      const KineticSpecies& kinetic_species( *((m_state.data())[n]) );
      Real maximum( kinetic_species.maxValue() );
      Real minimum( kinetic_species.minValue() );
      pout() << "    Species " << n << ":\t" << maximum << " [max]\t" << minimum << " [min]" << std::endl;
      if (procID()==0)
         cout << "    Species " << n << ":\t" << maximum << " [max]\t" << minimum << " [min]" << std::endl;
   }
   pout() << std::endl;
   if (procID()==0)
      cout << std::endl;
}

void GKSystem::printTimeStep(int cur_step, Real dt)
{
  if (procID() == 0) {
    Real dt_Vlasov     = stableDt_Vlasov(cur_step);
    Real dt_Collisions = stableDt_Collisions(cur_step);
    cout << "  ----\n";
    cout << "  dt: " << dt << std::endl;
    cout << "    dt scales:  ";
    cout << dt_Vlasov << " (Vlasov)  ";
    cout << dt_Collisions << " (Collisions)  ";
    cout << std::endl;
    cout << "    CFL:  ";
    cout << dt/dt_Vlasov << " (Vlasov)  ";
    cout << dt/dt_Collisions << " (Collisions)";
    cout << std::endl;
    cout << "  ----\n";
  }
}


inline void GKSystem::printParameters() const
{
   if (procID() == 0 && m_verbosity) {

      if ( m_mag_geom_type != "SingleNull" && m_mag_geom_type != "SNCore" ) {
         cout << "num_cells = ";
         for (int i=0; i<SpaceDim; i++) cout << m_num_cells[i] << " ";
         cout << endl;
         cout << "is_periodic = ";
         for (int i=0; i<SpaceDim; i++) cout << m_is_periodic[i] << " ";
         cout << endl;
         if (m_configuration_decomposition.size() > 0) {
            cout << "configuration_decomposition = ";
            for (int i=0; i<m_configuration_decomposition.size(); i++)
               cout << m_configuration_decomposition[i] << " ";
            cout << endl;
         }
         if (m_velocity_decomposition.size() > 0) {
            cout << "velocity_decomposition = ";
            for (int i=0; i<m_velocity_decomposition.size(); i++)
               cout << m_velocity_decomposition[i] << " ";
            cout << endl;
         }
         if (m_phase_decomposition.size() > 0) {
            cout << "phase_decomposition = ";
            for (int i=0; i<m_phase_decomposition.size(); i++)
               cout << m_phase_decomposition[i] << " ";
            cout << endl;
         }
      }
      else {
#if 0
         for (int n=0; n<a_blocks.size(); ++n) {
            cout << endl << "Block " << n << " data:" << endl;
            cout << "   configuration_decomp = " << a_blocks[n].cfgDecomp() << endl;
            cout << "   velocity_decomp      = " << a_blocks[n].velDecomp() << endl;
            cout << "   phase_decomp         = " << a_blocks[n].phaseDecomp() << endl;
         }
         cout << endl;
#endif
      }

      cout << "enforce_positivity = " << (m_enforce_step_positivity||m_enforce_stage_positivity) << endl;
      std::string ptype("stage");
      if (m_enforce_step_positivity)
         ptype = "step";
      cout << "enforce_positivity_type = " << ptype << endl;

   }
}



void GKSystem::parseParameters( ParmParse&         a_ppgksys )
{
   // Get magnetic geometry type
   a_ppgksys.get( "magnetic_geometry_mapping", m_mag_geom_type );

   // This determines the amount of diagnositic output generated
   a_ppgksys.query( "verbosity", m_verbosity );
   CH_assert( m_verbosity >= 0 );

   if ( m_mag_geom_type != "SingleNull" && m_mag_geom_type != "SNCore") {
      // Set the grid size
      m_num_cells.resize( PDIM );
      for (int i=0; i<PDIM; ++i) m_num_cells[i] = 0;
      a_ppgksys.getarr( "num_cells", m_num_cells, 0, PDIM );
      for (int i=0; i<PDIM; ++i) CH_assert( m_num_cells[i]>0 );

      // Determine which spatial directions are periodic
      m_is_periodic.resize(PDIM);
      vector<int> isPeriodic( PDIM ); // why should I have to do this?
      a_ppgksys.getarr( "is_periodic", isPeriodic, 0, PDIM );
      for (int dim=0; dim<SpaceDim; dim++)  {
         m_is_periodic[dim] = (isPeriodic[dim] == 1);
      }

      // Get the domain decomposition parameters
      if (a_ppgksys.contains("configuration_decomp")) {
         m_configuration_decomposition.resize( CFG_DIM );
         for (int i=0; i<CFG_DIM; ++i) m_configuration_decomposition[i] = 0;
         a_ppgksys.getarr( "configuration_decomp", m_configuration_decomposition, 0, CFG_DIM );
         for (int i=0; i<CFG_DIM; ++i) CH_assert( m_configuration_decomposition[i]>0 );
      }
      if (a_ppgksys.contains("phase_decomp")) {
         m_phase_decomposition.resize( PDIM );
         for (int i=0; i<PDIM; ++i) m_phase_decomposition[i] = 0;
         a_ppgksys.getarr( "phase_decomp", m_phase_decomposition, 0, PDIM );
         for (int i=0; i<PDIM; ++i) CH_assert( m_phase_decomposition[i]>0 );
      }
   }
   else {
      // Set the velocity space grid size
      m_num_velocity_cells.resize( VEL_DIM );
      for (int i=0; i<VEL_DIM; ++i) m_num_velocity_cells[i] = 0;
      a_ppgksys.getarr( "num_velocity_cells", m_num_velocity_cells, 0, VEL_DIM );
      for (int i=0; i<VEL_DIM; ++i) CH_assert( m_num_velocity_cells[i]>0 );
   }

   if (a_ppgksys.contains("velocity_decomp")) {
      m_velocity_decomposition.resize( VEL_DIM );
      for (int i=0; i<VEL_DIM; ++i) m_velocity_decomposition[i] = 0;
      a_ppgksys.getarr( "velocity_decomp", m_velocity_decomposition, 0, VEL_DIM );
      for (int i=0; i<VEL_DIM; ++i) CH_assert( m_velocity_decomposition[i]>0 );
   }

   // Should we make an hdf file for the potential?
   a_ppgksys.query("hdf_potential",m_hdf_potential);

   // Should we make an hdf file for the electric field?
   a_ppgksys.query("hdf_efield",m_hdf_efield);

   // Should we make an hdf file for the distribution function?
   a_ppgksys.query("hdf_dfn",m_hdf_dfn);

   // Should we make an hdf file for the function minus Maxwellian?
   a_ppgksys.query("hdf_deltaF",m_hdf_deltaF);

   // Should we make an hdf file for the distribution function at a specified mu?
   a_ppgksys.query("hdf_dfn_at_mu",m_hdf_dfn_at_mu);

   // Should we make hdf files for f versus vparallel and poloidal angle?
   a_ppgksys.query("hdf_vpartheta",m_hdf_vpartheta);

   // Should we make hdf files for f versus radius and poloidal angle at a specified vpar, mu?
   a_ppgksys.query("hdf_rtheta",m_hdf_rtheta);

   // Should we make hdf files for vparallel-mu at a specified configuration space point?
   a_ppgksys.query("hdf_vparmu",m_hdf_vparmu);

   // Should we make hdf files for charge density?
   a_ppgksys.query("hdf_density",m_hdf_density);

   // Should we make an hdf file for the total charge density?
   a_ppgksys.query("hdf_total_density",m_hdf_total_density);

   // Should we make hdf files for ParallelMomentum?
   a_ppgksys.query("hdf_ParallelMomentum",m_hdf_ParallelMomentum);

   // Should we make hdf files for PoloidalMomentum?
   a_ppgksys.query("hdf_PoloidalMomentum",m_hdf_PoloidalMomentum);

   // Should we make hdf files for pressure?
   a_ppgksys.query("hdf_pressure",m_hdf_pressure);

   // Should we make hdf files for temperature?
   a_ppgksys.query("hdf_temperature",m_hdf_temperature);

   // Should we make hdf files for fourthMoment?
   a_ppgksys.query("hdf_fourthMoment",m_hdf_fourthMoment);

   // Should we make hdf files for particle flux?
   a_ppgksys.query("hdf_ParticleFlux",m_hdf_ParticleFlux);

   // Should we make hdf files for heat flux?
   a_ppgksys.query("hdf_HeatFlux",m_hdf_HeatFlux);

   // Should we make an hdf file for the Vlasov divergence?
   a_ppgksys.query("hdf_vlasov_divergence",m_hdf_vlasov_divergence);

   // At what fixed phase space indices should I plot?  (Indices plotted against in a given plot
   //   are ignored.  Specify in 5D; toroidal index ignored in 4D and set to zero in arguments
   //   of hdf write methods.
   m_fixed_plotindices.resize( 5 );
   for (int i=0; i<5; ++i) m_fixed_plotindices[i]=0;
   a_ppgksys.queryarr("fixed_plot_indices",m_fixed_plotindices,0,5);
   if ( m_mag_geom_type != "SingleNull" && m_mag_geom_type != "SNCore" ) {       // FIX THIS TO WORK WITH MULTIBLOCK
      // check to make sure these are not out of bounds.
      for (int i=0;i<2; ++i)
         {
            CH_assert( m_fixed_plotindices[i] >= 0 );
            CH_assert( m_fixed_plotindices[i] < m_num_cells[i] );
         }
#if CFG_DIM == 3
      CH_assert( m_fixed_plotindices[2] >= 0 );
      CH_assert( m_fixed_plotindices[2] < m_num_cells[3] );
#endif
      CH_assert( m_fixed_plotindices[3] >= -m_num_cells[VPARALLEL_DIR]/2 );
      CH_assert( m_fixed_plotindices[3] < m_num_cells[VPARALLEL_DIR]/2 );
      CH_assert( m_fixed_plotindices[4] >= 0 );
      CH_assert( m_fixed_plotindices[4] < m_num_cells[MU_DIR] );
   }
   else {
#if 0
      for (int block=0; block<a_config_blocks.size(); ++block) {
         for (int dir=0; dir<CFG_DIM; ++dir) {
            CH_assert( a_config_blocks[block].loMappedIndex()[dir] <= m_fixed_plotindices[dir] );
            CH_assert( a_config_blocks[block].hiMappedIndex()[dir] >  m_fixed_plotindices[dir] );
         }
      }
      CH_assert( m_fixed_plotindices[3] >= -m_num_velocity_cells[0]/2 );
      CH_assert( m_fixed_plotindices[3] < m_num_velocity_cells[0]/2 );
      CH_assert( m_fixed_plotindices[4] >= 0 );
      CH_assert( m_fixed_plotindices[4] < m_num_velocity_cells[1] );
#endif
   }

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

}



void GKSystem::postStageAdvance( KineticSpeciesPtrVect& a_soln )
{
   if (m_enforce_stage_positivity) {
      enforcePositivity( a_soln );
   }
}



#ifdef BLOB_ERROR_TEST

static bool dir_created = false;

void GKSystem::computeError(const double                 a_time,
                            const int                    a_step,
                            const KineticSpeciesPtrVect& a_species,
                            double&                      a_l1_error,
                            double&                      a_l2_error,
                            double&                      a_max_error)
{

  const LevelData<FArrayBox>& solution = a_species[0]->distributionFunction();
  const PhaseGeom& geom = a_species[0]->phaseSpaceGeometry();
  const PhaseCoordSys& coords = geom.phaseCoordSys();

  const DisjointBoxLayout& grids = solution.disjointBoxLayout();

  LevelData<FArrayBox> exact(grids, 1, IntVect::Unit);

  ParmParse pp("kinetic_function_library.blob");

  double amplitude;
  pp.get( "amplitude", amplitude );

  Vector<Real> temp( PDIM );
  temp.assign( 0.0 );
  pp.queryarr( "location", temp, 0, PDIM );
  RealVect location = RealVect( temp );

  temp.assign( 1.0 );
  pp.queryarr( "width", temp, 0, PDIM );
  RealVect width = RealVect( temp );

  temp.assign( 0.0 );
  ParmParse pp_miller("gksystem.magnetic_geometry_mapping.miller");
  pp_miller.queryarr( "origin", temp, 0, CFG_DIM );
  RealVect origin = RealVect( temp );

  double d[CFG_DIM];
  for (int dir=0; dir<CFG_DIM; ++dir) {
     d[dir] = location[dir] - origin[dir];
  }
  double r = sqrt( d[0]*d[0] + d[1]*d[1] );
  double theta_0 = atan2(d[1],d[0]);

  location[0] = r*cos(a_time + theta_0) + origin[0];
  location[1] = r*sin(a_time + theta_0) + origin[1];

  LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
  geom.getCellVolumes(volume);

  // Iterate over patches
  DataIterator dit = solution.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    const PhaseBlockCoordSys& block_coord_sys = geom.getBlockCoordSys(grids[dit]);

    Box box(grids[dit]);
    box.grow(1);
    FArrayBox cell_center_coords( box, SpaceDim );
    block_coord_sys.getCellCenteredRealCoords( cell_center_coords );

    FORT_SET_LOCALIZED(
                       CHF_FRA(exact[dit]),
                       CHF_BOX(box),
                       CHF_CONST_FRA(cell_center_coords),
                       CHF_CONST_REAL(amplitude),
                       CHF_CONST_REALVECT(location),
                       CHF_CONST_REALVECT(width));

    fourthOrderAverageCell(exact[dit]);

    exact[dit] -= solution[dit];
    exact[dit].abs();
  }

  if (a_step%1 == 0) {

    //  Plot max norm error

    if(procID()==0 && !dir_created) {
      mkdir("error_plots", 0777);
      dir_created = true;
    }

    char iter_str[100];
    sprintf(iter_str, "error_plots/error.%04d.", a_step);
    stringstream filename;
    filename << iter_str;

    geom.plotAtVelocityIndex(filename.str(), VEL::IntVect::Zero, exact);
  }

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
#endif


#include "NamespaceFooter.H"
