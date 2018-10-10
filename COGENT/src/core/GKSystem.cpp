#include "GKSystem.H"

#if 1  // warning, OS dependencies, will not work on all platforms
#include <sys/stat.h>
#endif

#include "LogRectPhaseCoordSys.H"
#include "SingleNullPhaseCoordSys.H"
#include "SNCorePhaseCoordSys.H"
#include "Directions.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "LoadBalance.H"
#include "FluidOpVarFactory.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM VEL_DIM
#include "BoxIterator.H"
#include "LoadBalance.H"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM-1
#include "Slicing.H.transdim"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM


#include "NamespaceHeader.H"

using namespace CH_MultiDim;


GKSystem::GKSystem( ParmParse& a_pp, bool a_use_external_TI )
   :
     m_using_electrons(false),
     m_enforce_stage_positivity(false),
     m_enforce_step_positivity(false),
     m_max_grid_size(0),
     m_ghostVect(4*IntVect::Unit),
     m_ti_class("rk"),
     m_ti_method("4"),
     m_gk_ops(NULL),
     m_state_comp( GKState(m_ghostVect) ),
     m_hdf_potential(false),
     m_hdf_efield(false),
     m_hdf_density(false),
     m_hdf_momentum(false),
     m_hdf_vpartheta(false),
     m_hdf_rtheta(false),
     m_hdf_pressure(false),
     m_hdf_parallelHeatFlux(false),
     m_hdf_temperature(false),
     m_hdf_fourthMoment(false),
     m_hdf_ParallelMomentum(false),
     m_hdf_PoloidalMomentum(false),
     m_hdf_ParticleFlux(false),
     m_hdf_HeatFlux(false),
     m_hdf_vparmu(false),
     m_hdf_AmpereErIncrement(false),
     m_hdf_total_density(false),
     m_hdf_dfn(false),
     m_hdf_deltaF(false),
     m_hdf_dfn_at_mu(false),
     m_hdf_fluids(false),
     m_verbosity(0),
     m_use_native_time_integrator( !a_use_external_TI ),
     m_compute_op_matrices(false),
     m_op_matrices_tolerance(1e-6)
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

   createState();

   createGlobalDOF();
   m_state_comp.setGlobalDOF(&m_global_dof);

   const Real BASE_DT( 1.0 );
   if (!m_use_native_time_integrator) {
      m_gk_ops = new GKOps;
      m_gk_ops->define( m_state_comp, BASE_DT );
      m_rhs.define( m_state_comp );
   }
   else {
      m_serialized_vector.define(m_state_comp);
      if (m_ti_class == _TI_RK_) {
         m_integrator = new TiRK<GKVector, GKOps>;
      }
      else if (m_ti_class == _TI_ARK_) {
         m_integrator = new TiARK<GKVector, GKOps>;
      } 
      else {
         MayDay::Error("Unrecognized input for m_ti_class.");
      }
      m_integrator->define( a_pp, m_ti_method, m_serialized_vector, BASE_DT );
      m_gk_ops = &( m_integrator->getOperators() );
   }
   
   if ( m_using_electrons && m_gk_ops->usingBoltzmannElectrons() ) {
      MayDay::Error( "GKSystem::createSpecies():  Electrons input as both kinetic and Boltzmann" );
   }
   else {
      m_using_electrons |= m_gk_ops->usingBoltzmannElectrons();
   }

   bool fixed_efield( m_gk_ops->fixedEField() );
   if (fixed_efield && m_using_electrons)
      MayDay::Warning( "Using electrons with fixed E field" );
      
   //   if (!fixed_efield && !m_using_electrons)
   //      MayDay::Warning( "Not using electrons with dynamic E field" );
   
   setupFieldHistories();

   if (m_compute_op_matrices && m_use_native_time_integrator) {
     m_op_matrices.define(m_gk_ops, m_serialized_vector, m_op_matrices_tolerance);
   }
}



void
GKSystem::initialize( const int a_cur_step )
{
   if ( a_cur_step == 0 ) {
      // If this is the first step, then set the initial conditions for the
      // full system
      m_gk_ops->initializeState( m_state_comp, 0.0 );
   }
   else if ( m_gk_ops->fixedEField() ) {
      // If this is a restart and the fixed_efield option is true, just
      // set the potential to its initial condition.  If the fixed_efield
      // option is false, the potential will be set in the below call
      // to m_gk_ops->initialize(), which will compute both the potential
      // and associated field.
      m_gk_ops->initializePotential( 0.0 );
   }

   // Initialize the physical state variables
   m_gk_ops->convertToPhysical( m_state_comp, m_state_phys );
   
   // Initialize the electric field:
   // a.  If the fixed_efield option is true, then the field is calculated
   //     from the potential set in prior call to m_gk_ops->initializeState()
   //     or (if restarting) m_gk_ops->initializePotential().
   // b.  If the fixed_efield option is false, then both the potential and
   //     associated field are computed.
   m_gk_ops->initializeElectricField( m_state_phys, a_cur_step );
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
#ifdef with_petsc
   if (!m_use_native_time_integrator) {
      delete m_gk_ops;
   }
#endif
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
       + string(CFG::SingleNullCoordSys::pp_name);
    ParmParse pp_geom( prefix.c_str() );

    m_mag_geom_coords = new CFG::SingleNullCoordSys(pp_grid, pp_geom);
  }
  else if ( m_mag_geom_type == "SNCore"  ) {

    ParmParse pp_grid("sncore");

    string prefix = mag_geom_prefix + string(".")
       + string(CFG::SNCoreCoordSys::pp_name);
    ParmParse pp_geom( prefix.c_str() );

    m_mag_geom_coords = new CFG::SNCoreCoordSys(pp_grid, pp_geom);
  }
  else if (   m_mag_geom_type == "miller"
           || m_mag_geom_type == "slab"
           || m_mag_geom_type == "cylindrical" ) {
     
    string prefix = mag_geom_prefix + string(".") + m_mag_geom_type;
     
    ParmParse pp( prefix.c_str() );

    m_mag_geom_coords
      = new CFG::LogRectCoordSys(pp,
                                 m_mag_geom_type,
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
  m_mag_geom->writeGeometryData(grids,0);
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

    CFG::IntVect n_loc = CFG::IntVect::Zero;

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

    int box_cell_num(1);
    for (int dir=0; dir<CFG_DIM; ++dir) {
      box_cell_num *= n_loc[dir];
    }
    
    if (box_cell_num > 0) {
      CFG::IntVect box_size(n_loc);
      CFG::Box patch(domain_box.smallEnd(), domain_box.smallEnd() + box_size-1);
      CFG::IntVect lo = CFG::IntVect::Zero;
      CFG::IntVect hi;
      for (int dir=0; dir<CFG_DIM; ++dir) {
        hi[dir] = domain_box.size(dir)/n_loc[dir]-1;
      }
      CFG::Box skeleton(lo, hi);
      CFG::BoxIterator bit(skeleton);
      for (bit.begin();bit.ok();++bit) {
        CFG::Box thisBox = patch + bit()*box_size;
        boxes.push_back(thisBox);
      }
    }
    else {
      MayDay::Error( "Configuration domain box cannot be load balanced" );
    }
  }

  // Make the layout.  This is where boxes are assigned to processes.
  Vector<int> procMap;
  CFG::LoadBalance( procMap, boxes );

  CFG::ProblemDomain prob_domain;

  if (   m_mag_geom_type == "miller" || m_mag_geom_type == "slab"
      || m_mag_geom_type == "cylindrical" ) {
 
     CFG::Box bounding_box;
     for (int n=0; n<boxes.size(); n++) {
        bounding_box = minBox(bounding_box, boxes[n]);
     }
     bool is_periodic[SpaceDim];
     for (int dir=0; dir<SpaceDim; ++dir) {
        is_periodic[dir] = m_is_periodic[dir];
     }
     prob_domain = CFG::ProblemDomain(bounding_box, is_periodic);
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

   m_dv[0] = dv_parallel;
   m_dv[1] = dmu;
   
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
   
   int n_loc[VEL_DIM];

   for (int dir=0; dir<VEL_DIM; ++dir) {
      n_loc[dir] = 0;
   }

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
         const VEL::Box& local_box = boxes[n];
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
    m_domains.push_back( domain );
  }

  // Construct the multiblock phase space coordinate system

  if (   m_mag_geom_type == "miller" || m_mag_geom_type == "slab"
      || m_mag_geom_type == "cylindrical" ) {
     
    m_phase_coords = new LogRectPhaseCoordSys(a_ppgksys,
                                              *(CFG::LogRectCoordSys*)m_mag_geom_coords,
                                              *m_velocity_coords,
                                              m_domains );
  }
  else if ( m_mag_geom_type == "SingleNull" ) {
    ParmParse pp("singlenull.decomp");
    m_phase_coords = new SingleNullPhaseCoordSys( pp,
                                                  *(CFG::SingleNullCoordSys*)m_mag_geom_coords,
                                                  *m_velocity_coords,
                                                  m_domains );
  }
  else if ( m_mag_geom_type == "SNCore" ) {
    ParmParse pp("sncore.decomp");
    m_phase_coords = new SNCorePhaseCoordSys( pp,
                                              *(CFG::SNCoreCoordSys*)m_mag_geom_coords,
                                              *m_velocity_coords,
                                              m_domains );
  }
  else {
    MayDay::Error("Invalid magnetic geometry type");
  }

  Vector<IntVect> decomps;
  for (int block=0; block<m_mag_geom_coords->numBlocks(); ++block) {
    decomps.push_back( m_phase_coords->getDecomposition(block) );
  }

  // Construct the phase space grid

  m_phase_grid = new PhaseGrid(m_domains, decomps, m_mag_geom_type);

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

inline
void GKSystem::createState()
{
   KineticSpeciesPtrVect kinetic_species;
   createKineticSpecies( kinetic_species );

   CFG::FluidSpeciesPtrVect fluid_species;
   createFluidSpecies( fluid_species );

   ScalarPtrVect scalars;
   createScalars( scalars );

   // Define the computational state object by cloning from the argument vectors.  Ghost cells are added
   // for the kinetic and fluid species.
   m_state_comp.define( kinetic_species, fluid_species, scalars, m_phase_geom );

   // Define the physical state object by cloning from the argument vectors.  No ghost cells are included, since
   // this state object is generally used for output.
   m_state_phys.define( m_state_comp, IntVect::Zero );
}

void
GKSystem::createScalars( ScalarPtrVect& a_scalars )
{
   /*
     Create the vector of time-evolving Real scalars.
   */

   if ( m_verbosity && procID() == 0 ) {
      cout << "Adding scalar variables..." << endl;
   }

   bool more_scalars(true);

   int scalar_number = 1;

   while ( more_scalars ) {

      bool scalars_is_complete(true);

      // look for data specifying another field
      stringstream s;
      s << "scalar." << scalar_number;
      ParmParse ppscalar( s.str().c_str() );
      
      string name;
      if (ppscalar.contains("name")) {
         ppscalar.get("name", name);
      }
      else {
         scalars_is_complete = false;
      }
      
      string op_type;
      if (ppscalar.contains("operator_type")) {
         ppscalar.get("operator_type", op_type);
      }
      else {
         scalars_is_complete = false;
      }

      if ( scalars_is_complete ) {
         
         int depth;
         if ( op_type == "SelfConsistentBCOp" ) {
            depth = 2;
         }
         else {
            MayDay::Error("GKSystem::createScalars(): Unrecognized scalar operator type");
         }

         // Create the scalar object
         Scalar* scalar = new Scalar( name, depth );

         // Add the new scalar to the solution vector
         a_scalars.push_back(ScalarPtr(scalar));

         scalar_number++;
      }
      else {
         more_scalars = false;
      }
   }

   // Check that the required scalar state variables have been defined to support the self-consistent
   // potential boundary condition or Ampere options.

   ParmParse ppgksys("gksystem");

   bool consistent_potential_bcs = false;
   if (ppgksys.contains("consistent_potential_bcs")) {
      ppgksys.query( "consistent_potential_bcs", consistent_potential_bcs );
   }

   bool ampere_law = false;
   if (ppgksys.contains("ampere_law")) {
      ppgksys.query( "ampere_law", ampere_law );
   }

   if ( consistent_potential_bcs || ampere_law ) {

      bool self_consist_defined = false;
      for (int s=0; s<a_scalars.size(); ++s) {
         ScalarPtr scalar_ptr(a_scalars[s]);
         self_consist_defined = ( scalar_ptr->name() == "Er_boundary" );
         if (self_consist_defined ) break;
      }
      if ( !self_consist_defined ) {
         MayDay::Error("GKSystem::createScalars(): Er_boundary must be defined as a scalar state variable when using the consistent_potential_bcs or Ampere option");
      }
   }

   if ( procID() == 0 && m_verbosity ) {
      cout << "Done adding scalar variables" << endl;
   }
}



void
GKSystem::createFluidSpecies( CFG::FluidSpeciesPtrVect& a_fluid_species )
{
   /*
     Create the vector of configuration space variable (pointers)
   */

   if ( m_verbosity && procID() == 0 ) {
      cout << "Adding fluid species..." << endl;
   }

   bool more_vars(true);
   CFG::FluidOpVarFactory var_factory;

   while ( more_vars ) {

      // Look for fluid species.  At a minimum, a name and an operator type
      // must be provided, otherwise the search ends.  Any other information
      // is extracted from the remaining ParmParsed input.

      stringstream s;
      s << "fluid_species." << a_fluid_species.size() + 1;

      ParmParse ppspecies( s.str().c_str() );
      
      string name;
      if (ppspecies.contains("name")) {
         ppspecies.get("name", name);
      }
      else {
         more_vars = false;
      }
      
      string op_type;
      if (ppspecies.contains("operator_type")) {
         ppspecies.get("operator_type", op_type);
      }
      else {
         more_vars = false;
      }
      
      if ( more_vars ) {
         CFG::CFGVars* cfg_var = var_factory.create( s.str(), name, op_type, *m_mag_geom, CFG::IntVect::Zero );
         a_fluid_species.push_back( CFG::FluidSpeciesPtr(cfg_var) );
      }
   }

   ParmParse ppgksys("gksystem");
   bool ampere_law = false;
   if (ppgksys.contains("ampere_law")) {
      ppgksys.query( "ampere_law", ampere_law );
   }

   if ( ampere_law ) {

      bool er_flux_surfaces_defined = false;
      for (int s=0; s<a_fluid_species.size(); ++s) {
         CFG::FluidSpeciesPtr fluid_ptr(a_fluid_species[s]);
         er_flux_surfaces_defined = ( fluid_ptr->name() == "Er_flux_surfaces" );
         if (er_flux_surfaces_defined ) break;
      }
      if ( !er_flux_surfaces_defined ) {
         MayDay::Error("GKSystem::createFluidSpecies(): Er_flux_surfaces must be defined as a fluid state variable when using the Ampere option");
      }
   }

   if ( procID() == 0 && m_verbosity ) {
      cout << "Done adding fluid species" << endl;
   }
}



void
GKSystem::createKineticSpecies( KineticSpeciesPtrVect& a_kinetic_species )
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
        s << "kinetic_species." << a_kinetic_species.size()+i_off;
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
           i_electron = a_kinetic_species.size() + 1;
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

      bool include_velocity_renormalization(false);
      if ( ppspecies.contains("velocity_renormalization") ) {
	ppspecies.get("velocity_renormalization", include_velocity_renormalization);
      }

      if ( species_is_complete ) {
         if ( procID() == 0 && m_verbosity ) {
            cout << "   " << name << ": kinetic species with mass = " << mass << ", charge = "
                 << charge << ", velocity renormalization = " << include_velocity_renormalization <<endl;
         }

         //Construct species-depended velocity space
         string prefix = string("gksystem.") + string(VEL::VelCoordSys::pp_name);
         ParmParse pp_vel(prefix.c_str());

         const VEL::ProblemDomain& domain = getVelocityDomain();
        
         VEL::DisjointBoxLayout grids;
         getVelocitySpaceDisjointBoxLayout( domain, grids );

         VEL::RealVect dv(m_dv);
         if (include_velocity_renormalization) dv[0] = m_dv[0] / sqrt(mass);

         VEL::VelCoordSys* velocity_coords = new VEL::VelCoordSys(pp_vel, grids, domain, dv);

        
         
         // Construct the multiblock species-dependent phase space coordinate system
         
         PhaseCoordSys* phase_coords;
         if (   m_mag_geom_type == "miller" || m_mag_geom_type == "slab"
             || m_mag_geom_type == "cylindrical" ) {
            ParmParse pp("gksystem");
            phase_coords = new LogRectPhaseCoordSys(pp,
                                                    *(CFG::LogRectCoordSys*)m_mag_geom_coords,
                                                    *velocity_coords,
                                                    m_domains );
         }
         else if ( m_mag_geom_type == "SingleNull" ) {
            ParmParse pp("singlenull.decomp");
            phase_coords = new SingleNullPhaseCoordSys(pp,
                                                       *(CFG::SingleNullCoordSys*)m_mag_geom_coords,
                                                       *velocity_coords,
                                                       m_domains );
         }
         else if ( m_mag_geom_type == "SNCore" ) {
            ParmParse pp("sncore.decomp");
            phase_coords = new SNCorePhaseCoordSys(pp,
                                                   *(CFG::SNCoreCoordSys*)m_mag_geom_coords,
                                                   *velocity_coords,
                                                   m_domains );
      
         }
         
         // Get the species geometry
         PhaseGeom* species_geom = new PhaseGeom(*m_phase_geom, phase_coords, *velocity_coords, mass, charge);

         // Create the species object
         KineticSpecies* kin_spec = new KineticSpecies( name, mass, charge, *species_geom );
         kin_spec->distributionFunction().define(m_phase_grid->disjointBoxLayout(), 1, IntVect::Zero);

         // Add the new species to the solution vector
         KineticSpeciesPtr species(kin_spec);
         a_kinetic_species.push_back(species);

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



Real GKSystem::stableDt( const int a_step_number )
{
  if (m_use_native_time_integrator) {
    if ( m_integrator->isExplicit() ) {
      return m_gk_ops->stableDtExpl( m_state_comp, a_step_number );
    } else if ( m_integrator->isImEx() ) {
      return m_gk_ops->stableDtImEx( m_state_comp, a_step_number );
    } else {
      return (1000.0*m_gk_ops->stableDtExpl( m_state_comp, a_step_number ));
    }
  } else {
    return m_gk_ops->stableDtExpl( m_state_comp, a_step_number );
  }
  return DBL_MAX;
}



void GKSystem::enforcePositivity( KineticSpeciesPtrVect& a_soln )
{
   for (int n(0); n<a_soln.size(); n++) {
      KineticSpecies* kinetic_species( a_soln[n].operator->() );
      m_positivity_post_processor.enforce( kinetic_species->distributionFunction(),
                                           kinetic_species->maxValue() );
   }
}



void GKSystem::advance( Real& a_cur_time,
                        Real& a_dt,
                        int&  a_step_number)
{
   CH_assert(m_use_native_time_integrator);
   m_integrator->setTimeStepSize( a_dt );

   m_state_comp.copyTo  ( m_serialized_vector.data() );
   if (m_compute_op_matrices) {
     m_op_matrices.computeMatrices(m_serialized_vector,a_cur_time,0);
   }
   m_integrator->advance( a_cur_time, m_serialized_vector );
   m_state_comp.copyFrom( m_serialized_vector.data() );

   m_integrator->getCurrentTime( a_cur_time );
   m_integrator->getTimeStep( a_step_number );
   m_integrator->getTimeStepSize( a_dt );

   if (m_enforce_step_positivity) {
      enforcePositivity( m_state_comp.dataKinetic() );
   }
   m_gk_ops->convertToPhysical( m_state_comp, m_state_phys );
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


void GKSystem::writePlotFile(const char    *prefix,
                             const int     cur_step,
                             const double& cur_time )
{
   // If the efield and potential are fixed, only consider plotting them at step 0
   if ( !m_gk_ops->fixedEField() || cur_step == 0 ) {
   
      if (m_hdf_potential) {
         std::string filename( plotFileName( prefix,
                                             "potential",
                                             cur_step ) );
         m_gk_ops->plotPotential( filename, cur_time );
      }

      if (m_hdf_efield) {
         std::string filename = plotFileName( prefix,
                                              "efield",
                                              cur_step );
         m_gk_ops->plotEField( filename, cur_time );
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

   const KineticSpeciesPtrVect& kinetic_species( m_state_phys.dataKinetic() );
   for (int species(0); species<kinetic_species.size(); species++) {
      const KineticSpecies& soln_species( *(kinetic_species[species]) );

      // Distribution function

      if (m_hdf_dfn) {
         std::string filename( plotFileName( prefix,
                                             "dfn",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotDistributionFunction( filename, soln_species, cur_time );

         // Write out BStarParallel at first step
         if (cur_step==0) {
            std::string filename( plotFileName( prefix,
                                                "bstar_par",
                                                soln_species.name(),
                                                cur_step,
                                                species + 1) );

            m_gk_ops->plotBStarParallel( filename, soln_species, cur_time );
         }
      }

      // Distribution function with subtracted Maxwellian 

      if (m_hdf_deltaF) {
         std::string filename( plotFileName( prefix,
                                             "deltaF",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotDeltaF( filename, soln_species, cur_time );
      }

      // Distribution function at a specified mu index

      if (m_hdf_dfn_at_mu) {
         std::string filename( plotFileName( prefix,
                                             "dfn_at_mu",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotDistributionFunctionAtMu( filename, soln_species, mu_index, cur_time );
      }

      // Distribution function vparallel versus theta at specified configuration space point

      if (m_hdf_vpartheta) {
         std::string filename( plotFileName( prefix,
                                             "vpar_poloidal",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotVParallelTheta( filename, soln_species, radial_index, toroidal_index, mu_index, cur_time );
      }

      // Distribution function r versus theta at specified velocity space point

      if (m_hdf_rtheta) {
         std::string filename( plotFileName( prefix,
                                             "f_rtheta",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotRTheta( filename, soln_species, vpar_index, mu_index, cur_time );
      }


      // Distribution function vparallel versus mu at fixed configuration coordinate

      if (m_hdf_vparmu) {
         std::string filename( plotFileName( prefix,
                                             "vpar_mu",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotVParallelMu( filename, soln_species, radial_index, poloidal_index, toroidal_index, cur_time );
      }

      // Charge density

      if (m_hdf_density) {
         std::string filename( plotFileName( prefix,
                                             "density",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotChargeDensity( filename, soln_species, cur_time );
      }

      // Momentum density                                                                                                                                                                                        

      if (m_hdf_momentum) {
	std::string filename( plotFileName( prefix,
					    "momentum",
					    soln_species.name(),
					    cur_step,
					    species + 1) );

	m_gk_ops->plotMomentum( filename, soln_species, cur_time );
      }

      // Parallel momentum moment

      if (m_hdf_ParallelMomentum) {
         std::string filename( plotFileName( prefix,
                                             "ParallelMomentum",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotParallelMomentum( filename, soln_species, cur_time );
      }

      // Poloidal momentum moment

      if (m_hdf_PoloidalMomentum) {
         std::string filename( plotFileName( prefix,
                                             "PoloidalMomentum",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotPoloidalMomentum( filename, soln_species, cur_time );
      }

      // Pressure

      if (m_hdf_pressure) {
         std::string filename( plotFileName( prefix,
                                             "pressure",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotPressure( filename, soln_species, cur_time );
      }

      // Parallel heat flux
      
      if (m_hdf_parallelHeatFlux) {
         std::string filename( plotFileName( prefix,
                                            "ParallelHeatFlux",
                                            soln_species.name(),
                                            cur_step,
                                            species + 1) );
         
         m_gk_ops->plotParallelHeatFlux( filename, soln_species, cur_time );
      }

      
      // Temperature

      if (m_hdf_temperature) {
         std::string filename( plotFileName( prefix,
                                             "temperature",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotTemperature( filename, soln_species, cur_time );
      }

      // Fourth momemnt

      if (m_hdf_fourthMoment) {
         std::string filename( plotFileName( prefix,
                                             "fourthMoment",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );

         m_gk_ops->plotFourthMoment( filename, soln_species, cur_time );
      }

      // Increment of Er in the Ampere model                                                                                                                           
      if (m_hdf_AmpereErIncrement) {
         std::string filename( plotFileName( prefix,
					    "AmpereErIncrement",
					    soln_species.name(),
					    cur_step,
					    species + 1) );
         m_gk_ops->plotAmpereErIncrement( filename, m_state_phys.dataFluid(), m_state_phys.dataScalar(), cur_time );
      }


      // Particle flux moment

      if (m_hdf_ParticleFlux) {
         std::string filename( plotFileName( prefix,
                                             "ParticleFlux",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotParticleFlux( filename, soln_species, cur_time );
      }

      // Heat flux moment

      if (m_hdf_HeatFlux) {
         std::string filename( plotFileName( prefix,
                                             "HeatFlux",
                                             soln_species.name(),
                                             cur_step,
                                             species + 1) );
         m_gk_ops->plotHeatFlux( filename, soln_species, cur_time );
      }
   }

   // Total charge density

   if (m_hdf_total_density) {
      std::string filename( plotFileName( prefix,
                                          "charge_density",
                                          cur_step ) );

      m_gk_ops->plotChargeDensity( filename, kinetic_species, cur_time );
   }

   // Fluid species variables

   if (m_hdf_fluids) {

      const CFG::FluidSpeciesPtrVect& fluids( m_state_phys.dataFluid() );
      for (int species(0); species<fluids.size(); species++) {
         const CFG::FluidSpecies& fluid_species( static_cast<const CFG::FluidSpecies&>(*(fluids[species])) );
      
         if ( fluid_species.num_cell_vars() == 1 &&
              fluid_species.name() == fluid_species.cell_var_name(0) ) {
            // Don't repeat the variable name if there's only one and it's the same as the species name
            std::string filename (plotFileName( prefix,
                                                fluid_species.name(),
                                                cur_step ));
            m_gk_ops->plotFluid( filename, fluid_species, fluid_species.cell_var_name(0), cur_time );
         }
         else {
            for (int n=0; n<fluid_species.num_cell_vars(); ++n) {
               std::string filename (plotFileName( prefix,
                                                   fluid_species.cell_var_name(n),
                                                   fluid_species.name(),
                                                   cur_step,
                                                   species + 1));
         
               m_gk_ops->plotFluid( filename, fluid_species, fluid_species.cell_var_name(n), cur_time );
            }
         }
      }
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
   header.m_real["Er_lo"]           = m_gk_ops->getLoRadialField(m_state_comp);
   header.m_real["Er_hi"]           = m_gk_ops->getHiRadialField(m_state_comp);
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
   
   if ( m_gk_ops->stateContainsPotential() ) {
      LevelData<FArrayBox> phi_injected;
      m_phase_geom->injectConfigurationToPhase( m_gk_ops->getPhi(), phi_injected);
      a_handle.setGroup( "potential" );
      write(a_handle,phi_injected.boxLayout());
      write(a_handle,phi_injected,"data",phi_injected.ghostVect());
   }
   
   const KineticSpeciesPtrVect& kinetic_species( m_state_comp.dataKinetic() );
   for (int species(0); species<kinetic_species.size(); species++) {

      // Get solution distribution function for the current species
      KineticSpecies& soln_species( *(kinetic_species[species]) );
      LevelData<FArrayBox> & soln_dfn = soln_species.distributionFunction();
      char buff[100];
      sprintf( buff, "dfn_%d", species + 1 );
      a_handle.setGroup( buff );
      write( a_handle, soln_dfn.boxLayout() );
      write( a_handle, soln_dfn, "data" );
   }

   const CFG::FluidSpeciesPtrVect& fluid_species( m_state_comp.dataFluid() );
   for (int species(0); species<fluid_species.size(); species++) {

      CFG::CFGVars& this_fluid_species( *(fluid_species[species]) );

      char buff[100];
      for (int n=0; n<this_fluid_species.num_cell_vars(); ++n) {
      
         CFG::LevelData<CFG::FArrayBox> & fluid_data = this_fluid_species.cell_var(n);
         const string& var_name = this_fluid_species.cell_var_name(n);
         sprintf( buff, "fluid_%d_%s", species + 1, var_name.c_str() );
         a_handle.setGroup( buff );

         LevelData<FArrayBox> injected_fluid_data;
         m_phase_geom->injectConfigurationToPhase( fluid_data, injected_fluid_data);

         write( a_handle, injected_fluid_data.boxLayout() );
         write( a_handle, injected_fluid_data, "data", injected_fluid_data.ghostVect() );
      }

      // NEED TO ADD FACE DATA
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
   m_gk_ops->setLoRadialField(m_state_comp, header.m_real["Er_lo"]);
   m_gk_ops->setHiRadialField(m_state_comp, header.m_real["Er_hi"]);
   //   const int restart_version = header.m_int ["restart_version"];

   if ( m_gk_ops->usingAmpereLaw() ) {
      a_handle.setGroup("Er_cell");
      CFG::LevelData<CFG::FArrayBox> Er_cell_config(m_mag_geom->grids(), 3, CFG::IntVect::Unit);
      LevelData<FArrayBox> Er_cell_injected;
      m_phase_geom->injectConfigurationToPhase(Er_cell_config, Er_cell_injected);
      read( a_handle, Er_cell_injected, "data", Er_cell_injected.disjointBoxLayout() );
      m_gk_ops->setErAverage( Er_cell_injected );
    
      a_handle.setGroup("Er_face");
      CFG::LevelData<CFG::FluxBox> Er_face_config(m_mag_geom->grids(), 3, CFG::IntVect::Unit);
      LevelData<FluxBox> Er_face_injected;
      m_phase_geom->injectConfigurationToPhase(Er_face_config, Er_face_injected);
      read( a_handle, Er_face_injected, "data", Er_face_injected.disjointBoxLayout() );
      m_gk_ops->setErAverage( Er_face_injected );

      a_handle.setGroup("E_tilde_cell");
      CFG::LevelData<CFG::FArrayBox> E_tilde_cell_config(m_mag_geom->grids(), 3, CFG::IntVect::Unit);
      LevelData<FArrayBox> E_tilde_cell_injected;
      m_phase_geom->injectConfigurationToPhase(E_tilde_cell_config, E_tilde_cell_injected);
      read( a_handle, E_tilde_cell_injected, "data", E_tilde_cell_injected.disjointBoxLayout() );
      m_gk_ops->setETilde( E_tilde_cell_injected );
    
      a_handle.setGroup("E_tilde_face");
      CFG::LevelData<CFG::FluxBox> E_tilde_face_config(m_mag_geom->grids(), 3, CFG::IntVect::Unit);
      LevelData<FluxBox> E_tilde_face_injected;
      m_phase_geom->injectConfigurationToPhase(E_tilde_face_config, E_tilde_face_injected);
      read( a_handle, E_tilde_face_injected, "data", E_tilde_face_injected.disjointBoxLayout() );
      m_gk_ops->setETilde( E_tilde_face_injected );
   }
   
   if ( m_gk_ops->stateContainsPotential() ) {
      a_handle.setGroup("potential");
      CFG::LevelData<CFG::FArrayBox> phi(m_mag_geom->grids(), 1, CFG::IntVect::Zero);
      LevelData<FArrayBox> phi_injected;
      m_phase_geom->injectConfigurationToPhase(phi, phi_injected);
      read( a_handle, phi_injected, "data", phi_injected.disjointBoxLayout() );
      m_gk_ops->setPhi( phi_injected );
   }

   KineticSpeciesPtrVect& kinetic_species( m_state_comp.dataKinetic() );
   for (int species(0); species<kinetic_species.size(); species++) {

      // Get solution distribution function for the current species
      KineticSpecies& soln_species( *(kinetic_species[species]) );
      LevelData<FArrayBox>& soln_dfn = soln_species.distributionFunction();
      char buff[100];
      sprintf( buff, "dfn_%d", species + 1 );
      a_handle.setGroup( buff );
      read( a_handle, soln_dfn, "data", soln_dfn.disjointBoxLayout(), Interval(0,soln_dfn.nComp()-1), false );
   }

   CFG::FluidSpeciesPtrVect& fluid_species( m_state_comp.dataFluid() );
   for (int species(0); species<fluid_species.size(); species++) {

      CFG::CFGVars& this_fluid_species( *(fluid_species[species]) );

      char buff[100];
      for (int n=0; n<this_fluid_species.num_cell_vars(); ++n) {

         CFG::LevelData<CFG::FArrayBox>& fluid_data = this_fluid_species.cell_var(n);

         LevelData<FArrayBox> injected_fluid_data;
         m_phase_geom->injectConfigurationToPhase(fluid_data, injected_fluid_data);

         const string& var_name = this_fluid_species.cell_var_name(n);
         sprintf( buff, "fluid_%d_%s", species + 1, var_name.c_str() );
         a_handle.setGroup( buff );
         read( a_handle, injected_fluid_data, "data", injected_fluid_data.disjointBoxLayout() );
         m_phase_geom->projectPhaseToConfiguration(injected_fluid_data, fluid_data);
      }

      // NEED TO ADD FACE DATA
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
   if (procID()==0) {
      cout << "  Distribution Function Extrema:" << std::endl;
   }
   const KineticSpeciesPtrVect& kinetic_species( m_state_comp.dataKinetic() );
   for (int n(0); n<kinetic_species.size(); n++) {
      const KineticSpecies& species( *(kinetic_species[n]) );
      Real maximum( species.maxValue() );
      Real minimum( species.minValue() );
      pout() << "    Species " << n << ":\t"
             << maximum << " [max]\t" << minimum << " [min]" << std::endl;
      if (procID()==0) {
         cout << "    Species " << n << ":\t" << maximum << " [max]\t" << minimum << " [min]" << std::endl;
      }
   }
   pout() << std::endl;
   if (procID()==0) {
      cout << std::endl;
   }
}

void GKSystem::preTimeStep(int a_cur_step, Real a_cur_time)
{
   if (m_use_native_time_integrator) {
      m_integrator->setCurrentTime( a_cur_time );
      m_integrator->setTimeStep( a_cur_step );
   }
   m_gk_ops->convertToPhysical( m_state_comp, m_state_phys );
   m_gk_ops->preTimeStep( a_cur_step, a_cur_time, m_state_comp, m_state_phys );
}


inline
void printTimeStep( const Real& a_dt,
                    const std::string& a_name,
                    const Real& a_dt_current  )
{
   if (a_dt>= 0) {
      cout << a_name << a_dt
           << " (time scale), " << a_dt_current / a_dt << " (CFL)\n";
   }
}

void GKSystem::postTimeStep(int a_cur_step, Real a_dt, Real a_cur_time)
{
  m_gk_ops->postTimeStep( a_cur_step, a_dt, a_cur_time, m_state_comp );
  m_gk_ops->convertToPhysical( m_state_comp, m_state_phys );
  if (procID() == 0) {
    cout << "  ----\n";
    cout << "  dt: " << a_dt << std::endl;
    Real dt_vlasov = m_gk_ops->dtScaleVlasov( m_state_comp, a_cur_step );
    printTimeStep( dt_vlasov, "    Vlasov    : ", a_dt ); 
    Real dt_collisions = m_gk_ops->dtScaleCollisions( m_state_comp, a_cur_step);
    printTimeStep( dt_collisions, "    Collisions: ", a_dt ); 
    Real dt_transport = m_gk_ops->dtScaleTransport( m_state_comp, a_cur_step );
    printTimeStep( dt_transport, "    Transport : ", a_dt ); 
    Real dt_neutrals = m_gk_ops->dtScaleNeutrals( m_state_comp, a_cur_step );
    printTimeStep( dt_neutrals, "    Neutrals  : ", a_dt ); 
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

   // time integration method to use 
   a_ppgksys.query("ti_class",m_ti_class);
   a_ppgksys.query("ti_method",m_ti_method);

   // compute linearized operator matrices?
   // Caution: a *very* slow process
   a_ppgksys.query("compute_op_matrices",m_compute_op_matrices);
   a_ppgksys.query("op_matrices_tolerance",m_op_matrices_tolerance);

   // Should we make an hdf file for the potential?
   a_ppgksys.query("hdf_potential",m_hdf_potential);

   // Should we make an hdf file for the electric field?
   a_ppgksys.query("hdf_efield",m_hdf_efield);

   // Should we make an hdf file for the distribution function?
   a_ppgksys.query("hdf_dfn",m_hdf_dfn);

   // Should we make an hdf file for fluids?
   a_ppgksys.query("hdf_fluids",m_hdf_fluids);

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

   // Should we make hdf files for momentum?                                                                                                                                                               
   a_ppgksys.query("hdf_momentum",m_hdf_momentum);

   // Should we make an hdf file for the total charge density?
   a_ppgksys.query("hdf_total_density",m_hdf_total_density);

   // Should we make hdf files for ParallelMomentum?
   a_ppgksys.query("hdf_ParallelMomentum",m_hdf_ParallelMomentum);

   // Should we make hdf files for PoloidalMomentum?
   a_ppgksys.query("hdf_PoloidalMomentum",m_hdf_PoloidalMomentum);

   // Should we make hdf files for pressure?
   a_ppgksys.query("hdf_pressure",m_hdf_pressure);

   // Should we make hdf files for parallel heat flux?
   a_ppgksys.query("hdf_parallelHeatFlux",m_hdf_parallelHeatFlux);

   // Should we make hdf files for temperature?
   a_ppgksys.query("hdf_temperature",m_hdf_temperature);

   // Should we make hdf files for fourthMoment?
   a_ppgksys.query("hdf_fourthMoment",m_hdf_fourthMoment);

   // Should we make hdf files for particle flux?
   a_ppgksys.query("hdf_ParticleFlux",m_hdf_ParticleFlux);

   // Should we make hdf files for heat flux?
   a_ppgksys.query("hdf_HeatFlux",m_hdf_HeatFlux);

   // Should we make hdf files for Ampere Er increment?
   a_ppgksys.query("hdf_AmpereErIncrement",m_hdf_AmpereErIncrement);

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
      CH_assert( m_fixed_plotindices[2] < m_num_cells[2] );
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

#include "NamespaceFooter.H"
