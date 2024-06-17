#include "GKSystem.H"

#if 1  // warning, OS dependencies, will not work on all platforms
#include <sys/stat.h>
#endif

#include "LogRectPhaseCoordSys.H"
#include "SingleNullPhaseCoordSys.H"
#include "SNCorePhaseCoordSys.H"
#include "VelocityNormalization.H"
#include "SpaceUtils.H.multidim"
#include "Directions.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "LoadBalance.H"
#include "FluidOpVarFactory.H"
#include "CFGVars.H"
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


GKSystem::GKSystem( int a_sys_id )
   : m_sys_id(a_sys_id),
     m_using_electrons(false),
     m_enforce_stage_positivity(false),
     m_enforce_step_positivity(false),
     m_enforce_positivity_pointwise(false),
     m_enforce_step_floor(false),
     m_max_grid_size(0),
     m_kinetic_ghosts(4),
     m_fluid_ghosts(4),
     m_gk_ops(NULL),
     m_diagnostics(NULL),
     m_verbosity(0),
     m_ti_class("rk"),
     m_ti_method("4"),
     m_use_scales(false),
     m_scale_tol(1e-12),
     m_dt_collisions(-1),
     m_dt_transport(-1),
     m_dt_neutrals(-1)
{ 
   m_diagnostics_dfn_plots.clear();
   m_diagnostics_cfg_field_vars.clear();
   m_diagnostics_cfg_vars_kin_spec.clear();
   m_diagnostics_cfg_vars_total_kin_spec.clear();
   m_diagnostics_fluids = false;

   char sys_id_str[10];
   sprintf(sys_id_str, "%d", m_sys_id);
   m_sys_id_str = std::string(sys_id_str);
}

GKSystem::GKSystem( ParmParse&, int a_sys_id )
   : GKSystem(a_sys_id)
{
   std::string parse_key = "gksystem";
   GKSystemParameters gkparams;
   /* first read all inputs with "gksystem" */
   gkparams.readParams( parse_key );
   appendSysID(parse_key);
   /* now reread any available input with "gksystem.n" */
   gkparams.readParams( parse_key );

   define( gkparams );
   return;
}

GKSystem::GKSystem( const GKSystemParameters& a_gkparams,
                    int                       a_sys_id )
   : GKSystem(a_sys_id)
{
   define( a_gkparams );
}

void GKSystem::define(const GKSystemParameters& a_params )
{
   setParameters( a_params );
   
   m_units = new GKUnits();
   if (procID()==0) {
      m_units->print(cout);
   }

   createConfigurationSpace( a_params );
   createVelocitySpace();
   createPhaseSpace( a_params );

   createState( m_state_comp,
                a_params,
                m_kinetic_ghosts,
                m_fluid_ghosts );
   m_state_phys.define( m_state_comp, IntVect::Zero );

   createGlobalDOF();
   m_state_comp.setGlobalDOF(&m_global_dof);

   m_gk_ops = new GKOps;
   m_gk_ops->define( m_state_comp, a_params );

   m_diagnostics = new GKDiagnostics(m_gk_ops);

   m_rhs.define( m_state_comp );

   giveSpeciesTheirGyroaverageOps();
   
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
}

void GKSystem::defineShell(const GKSystemParameters& a_params )
{
   setParameters( a_params );
   
   m_units = new GKUnits();
   if (procID()==0) {
      m_units->print(cout);
   }

   createConfigurationSpace( a_params );
   createVelocitySpace();
   createPhaseSpace( a_params );

   createState( m_state_phys, a_params, 0, 0 );

   m_gk_ops = new GKOps;
   m_gk_ops->defineShell( m_state_phys, a_params );

   m_diagnostics = new GKDiagnostics(m_gk_ops);

   setupFieldHistories();
}

void
GKSystem::initialize( const int     a_cur_step,
                      const double  a_cur_time )
{
   CH_TIME("GKSystem::initialize");
   if ( a_cur_step == 0 ) {
      // If this is the first step, then set the initial conditions for the
      // full system
      m_gk_ops->initializeState( m_state_comp, a_cur_step, a_cur_time );
   }
   else if ( m_gk_ops->fixedEField() ) {
      // If this is a restart and the fixed_efield option is true, just
      // set the potential to its initial condition.  If the fixed_efield
      // option is false, the potential will be set in the below call
      // to m_gk_ops->initialize(), which will compute both the potential
      // and associated field.
      m_gk_ops->initializePotential( a_cur_time );
   }

   // Initialize physical kinetic & fluid species used in fluidOp (has ghost cells)
   m_gk_ops->initializeKineticSpeciesPhysical(  m_state_comp.dataKinetic(), 
                                                a_cur_time );
   m_gk_ops->initializeFluidSpeciesPhysical(  m_state_comp.dataFluid(), 
                                              a_cur_time );
   
   // Initialize the physical state variables
   m_gk_ops->convertToPhysical( m_state_comp, m_state_phys, a_cur_time );
   
   // Initialize the electric field:
   // a.  If the fixed_efield option is true, then the field is calculated
   //     from the potential set in prior call to m_gk_ops->initializeState()
   //     or (if restarting) m_gk_ops->initializePotential().
   // b.  If the fixed_efield option is false, then both the potential and
   //     associated field are computed.
   m_gk_ops->initializeElectricField( m_state_phys, a_cur_step, a_cur_time );

   // Apply initialization constraints to fluid models (e.g., MHD equilibirum)
   m_gk_ops->applyFluidInitializationConstraints(m_state_comp,
                                                 m_state_phys,
                                                 a_cur_step,
                                                 a_cur_time );
   
   m_gk_ops->initializeTI(  a_cur_step, 
                            a_cur_time, 
                            m_rhs,
                            m_state_comp, 
                            m_state_phys  );

   m_scale_u.define(m_state_comp); m_scale_u.ones();
   m_scale_rhsop.define(m_state_comp); m_scale_rhsop.ones();
   if (!trivialSolutionOp()) {
     m_scale_lhsop.define(m_state_comp); m_scale_lhsop.ones();
   }

   return;
}

void
GKSystem::initializeShell(  const int     a_cur_step,
                            const double  a_cur_time )
{
   CH_TIME("GKSystem::initialize");
   m_gk_ops->initializeElectricFieldShell(  m_state_phys, 
                                            a_cur_step, 
                                            a_cur_time );

   return;
}



GKSystem::~GKSystem()
{
   delete m_units;
   delete m_gk_ops;
   delete m_diagnostics;
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



void GKSystem::createConfigurationSpace( const GKSystemParameters& a_params )
{
  CH_TIME("GKSystem::createConfigurationSpace");

  string mag_geom_prefix = a_params.parseKey() + ".magnetic_geometry_mapping";

  if ( m_mag_geom_type == "SingleNull" ) {

    string parsekey = "singlenull";
    appendSysID(parsekey);
    ParmParse pp_grid(parsekey.c_str());

    string prefix = mag_geom_prefix + string(".")
       + string(CFG::SingleNullCoordSys::pp_name);
    ParmParse pp_geom( prefix.c_str() );

    m_mag_geom_coords = RefCountedPtr<CFG::MagCoordSys>(new CFG::SingleNullCoordSys(pp_grid, pp_geom));
  }
  else if ( m_mag_geom_type == "SNCore"  ) {

    string parsekey = "sncore";
    appendSysID(parsekey);
    ParmParse pp_grid(parsekey.c_str());

    string prefix = mag_geom_prefix + string(".")
       + string(CFG::SNCoreCoordSys::pp_name);
    ParmParse pp_geom( prefix.c_str() );

    m_mag_geom_coords = RefCountedPtr<CFG::MagCoordSys>(new CFG::SNCoreCoordSys(pp_grid, pp_geom));
  }
  else if (   m_mag_geom_type == "miller"
           || m_mag_geom_type == "slab"
           || m_mag_geom_type == "cylindrical"
           || m_mag_geom_type == "toroidal"
           || m_mag_geom_type == "oneblock" )   {
     
    string prefix = mag_geom_prefix + string(".") + m_mag_geom_type;
    ParmParse pp( prefix.c_str() );

    m_mag_geom_coords = RefCountedPtr<CFG::MagCoordSys>(
                          new CFG::LogRectCoordSys( pp,
                                                    m_mag_geom_type,
                                                    m_num_cells,
                                                    m_is_periodic,
                                                    m_configuration_decomposition));
  }
 
  // Construct the phase space DisjointBoxLayout

  CFG::DisjointBoxLayout grids;

  getConfigurationSpaceDisjointBoxLayout( grids );

  // Construct the magnetic geometry

  if (procID()==0) {
    cout << "Constructing magnetic geometry" << endl;
  }

  ParmParse pp_mag_geom( mag_geom_prefix.c_str() );

  int ghosts = (m_kinetic_ghosts > m_fluid_ghosts)  ? m_kinetic_ghosts 
                                                    : m_fluid_ghosts;
  m_mag_geom = RefCountedPtr<CFG::MagGeom>(new CFG::MagGeom(pp_mag_geom, 
                                                            m_mag_geom_coords,
                                                            grids, 
                                                            ghosts));

  if (procID()==0) {
    cout << "Done constructing magnetic geometry" << endl;
  }

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
       = (CFG::MagBlockCoordSys *)(m_mag_geom_coords->getCoordSys(block));

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
      || m_mag_geom_type == "cylindrical" || m_mag_geom_type == "toroidal"
      || m_mag_geom_type == "oneblock") {
 
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
        cout << "   ";
        cout  << "Configuration space box " 
              << local_box << " is assigned to process " 
              << procMap[n] << endl;
      }
    }
#ifdef CH_MPI
  }
#endif
}



void
GKSystem::createVelocitySpace()
{
   CH_TIME("GKSystem::createVelocitySpace");
   const VEL::ProblemDomain& domain = getVelocityDomain();

   VEL::DisjointBoxLayout grids;
   getVelocitySpaceDisjointBoxLayout( domain, grids );

   string parsekey( "phase_space_mapping" );
   ParmParse pppsm( parsekey.c_str() );

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
   ParmParse pp_vel(VEL::VelCoordSys::pp_name.c_str());
   m_velocity_coords = RefCountedPtr<VEL::VelCoordSys>(
                          new VEL::VelCoordSys( pp_vel, 
                                                grids, 
                                                domain, 
                                                dv  ));
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
         cout << "   ";
         cout << "Velocity space box " << local_box 
              << " is assigned to process " << procMap[n] << endl;
       }
     }
#ifdef CH_MPI
   }
#endif
}



void
GKSystem::createPhaseSpace( const GKSystemParameters& a_gkparams )
{
  // Build a vector of ProblemDomains, one per configuration space block

  CH_TIME("GKSystem::createPhaseSpace");
  const VEL::ProblemDomain& vel_domain = getVelocityDomain();

  for (int block=0; block<m_mag_geom_coords->numBlocks(); ++block) {

    // Construct the phase space block domain from the configuration and velocity domains

    const CFG::MagBlockCoordSys* mag_block_coords
       = (CFG::MagBlockCoordSys *)(m_mag_geom_coords->getCoordSys(block));

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
      || m_mag_geom_type == "cylindrical" || m_mag_geom_type == "toroidal"
      || m_mag_geom_type == "oneblock" ) {
     
     m_phase_coords = RefCountedPtr<PhaseCoordSys>(
                        new LogRectPhaseCoordSys( a_gkparams,
                                                  m_mag_geom_coords,
                                                  m_velocity_coords,
                                                  m_domains ));
  }
  else if ( m_mag_geom_type == "SingleNull" ) {

    string parsekey = "singlenull";
    appendSysID(parsekey);
    parsekey += ".decomp";
    ParmParse pp(parsekey.c_str());

    m_phase_coords = RefCountedPtr<PhaseCoordSys>(
                        new SingleNullPhaseCoordSys(  pp,
                                                      m_mag_geom_coords,
                                                      m_velocity_coords,
                                                      m_domains ));
  }
  else if ( m_mag_geom_type == "SNCore" ) {

    string parsekey = "sncore";
    appendSysID(parsekey);
    parsekey += ".decomp";
    ParmParse pp(parsekey.c_str());

    m_phase_coords = RefCountedPtr<PhaseCoordSys>(
                        new SNCorePhaseCoordSys(  pp,
                                                  m_mag_geom_coords,
                                                  m_velocity_coords,
                                                  m_domains ));
  }
  else {
    MayDay::Error("Invalid magnetic geometry type");
  }

  Vector<IntVect> decomps;
  for (int block=0; block<m_mag_geom_coords->numBlocks(); ++block) {
    decomps.push_back( m_phase_coords->getDecomposition(block) );
  }

  // Construct the phase space grid

  m_phase_grid = RefCountedPtr<PhaseGrid>(
                    new PhaseGrid(m_domains, decomps, m_mag_geom_type ));

  // Construct the phase space geometry

  if (m_verbosity>0) {
     m_phase_grid->print(m_kinetic_ghosts*IntVect::Unit);
  }
   
  m_phase_geom =
     RefCountedPtr<PhaseGeom>( new PhaseGeom( m_phase_coords,
                                              m_phase_grid,
                                              m_mag_geom,
                                              m_velocity_coords,
                                              m_kinetic_ghosts,
                                              m_units->larmorNumber() ) );
}

inline
void GKSystem::createState( GKState&                  a_state,
                            const GKSystemParameters& a_params,
                            const int                 a_gpt_kin,
                            const int                 a_gpt_fluid )
{
   CH_TIME("GKSystem::createState");
   KineticSpeciesPtrVect kinetic_species;
   createKineticSpecies( kinetic_species, a_params );

   CFG::FluidSpeciesPtrVect fluid_species;
   createFluidSpecies( fluid_species, a_params );

   ScalarPtrVect scalars;
   createScalars( scalars, a_params );

   a_state.define( kinetic_species,
                   fluid_species,
                   scalars,
                   m_phase_geom,
                   a_gpt_kin * IntVect::Unit,
                   a_gpt_fluid * CFG::IntVect::Unit );
  }

void
GKSystem::createScalars( ScalarPtrVect& a_scalars,
                         const GKSystemParameters& a_gkparams )
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

   bool consistent_potential_bcs = a_gkparams.consistentPotentialBCs();
   bool ampere_law = a_gkparams.ampereLaw();

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
GKSystem::createFluidSpecies( CFG::FluidSpeciesPtrVect& a_fluid_species,
                              const GKSystemParameters& a_gkparams )
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

   bool ampere_law = a_gkparams.ampereLaw();
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
GKSystem::createKineticSpecies( KineticSpeciesPtrVect& a_kinetic_species,
                                const GKSystemParameters& a_params )
{
   /*
     Create the vector of species model (pointers), and when a kinetic
     species is created, add it to the kinetic species vector
     number species.
   */

   if ( m_verbosity && procID() == 0 ) {
      cout << "Adding species and constructing gyrocenter coordinates..." 
           << endl;
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
      string velocity_normalization_type("None");
      if ( ppspecies.contains("velocity_renormalization") ) {
         ppspecies.get( "velocity_renormalization", 
                        include_velocity_renormalization  );
      }
      // NOTE Retaining boolian include_velocity_renormalization 
      // for backward compatibility of input files.
      if ( ppspecies.contains("velocity_normalization") ) {
         ppspecies.get("velocity_normalization", velocity_normalization_type);
      }

      bool is_gyrokinetic(false);
      ppspecies.query("is_gyrokinetic", is_gyrokinetic);

      if ( species_is_complete ) {
         if ( procID() == 0 && m_verbosity ) {
            cout  << "   " << (is_gyrokinetic?"Gyrokinetic":"Drift-kinetic") 
                  << " species "
                  << name << ": mass = " << mass << ", charge = "
                  << charge << ", velocity renormalization = " 
                  << include_velocity_renormalization
                  << ", velocity normalization = " 
                  << velocity_normalization_type <<endl;
         }

         //Construct species-depended velocity space
         ParmParse pp_vel(VEL::VelCoordSys::pp_name.c_str());

         const VEL::ProblemDomain& domain = getVelocityDomain();
        
         VEL::DisjointBoxLayout grids;
         getVelocitySpaceDisjointBoxLayout( domain, grids );

         VEL::RealVect dv(m_dv);
         if (    include_velocity_renormalization 
              || (velocity_normalization_type == "global")) {
            dv[0] = m_dv[0] / sqrt(mass);
         }
         RefCountedPtr<VEL::VelCoordSys> velocity_coords 
            = RefCountedPtr<VEL::VelCoordSys>(
                new VEL::VelCoordSys( pp_vel,
                                      grids,
                                      domain,
                                      dv  ));
         

         // Construct the species-dependent velocity normalization 
         RefCountedPtr<VelocityNormalization> velocity_normalization;
         velocity_normalization = RefCountedPtr<VelocityNormalization>(
                                    new VelocityNormalization(m_mag_geom,
                                                              velocity_coords,
                                                              ppspecies));

         // Construct the multiblock species-dependent 
         // phase space coordinate system
         RefCountedPtr<PhaseCoordSys> phase_coords;
         if (     m_mag_geom_type == "miller" 
              ||  m_mag_geom_type == "slab"
              ||  m_mag_geom_type == "cylindrical" 
              ||  m_mag_geom_type == "toroidal"
              ||  m_mag_geom_type == "oneblock" ) {

            phase_coords = RefCountedPtr<PhaseCoordSys>(
                              new LogRectPhaseCoordSys( a_params,
                                                        m_mag_geom_coords,
                                                        velocity_coords,
                                                        m_domains ));
         } else if ( m_mag_geom_type == "SingleNull" ) {

            string parsekey("singlenull");
            appendSysID(parsekey);
            parsekey += (".decomp");
            ParmParse pp(parsekey.c_str());

            phase_coords = RefCountedPtr<PhaseCoordSys>(
                              new SingleNullPhaseCoordSys(  pp,
                                                            m_mag_geom_coords,
                                                            velocity_coords,
                                                            m_domains ));
         }
         else if ( m_mag_geom_type == "SNCore" ) {

            string parsekey("sncore");
            appendSysID(parsekey);
            parsekey += (".decomp");
            ParmParse pp(parsekey.c_str());

            phase_coords = RefCountedPtr<PhaseCoordSys>(
                              new SNCorePhaseCoordSys(pp,
                                                      m_mag_geom_coords,
                                                      velocity_coords,
                                                      m_domains ));
         }
         
         // Get the species geometry
         RefCountedPtr<PhaseGeom> species_geom 
            = RefCountedPtr<PhaseGeom>(
                new PhaseGeom(  *m_phase_geom,
                                phase_coords,
                                velocity_coords,
                                velocity_normalization,
                                name,
                                mass,
                                charge,
                                is_gyrokinetic  ));

         // Create the species object
         KineticSpecies* kin_spec = new KineticSpecies( name, 
                                                        mass, 
                                                        charge, 
                                                        species_geom, 
                                                        is_gyrokinetic );

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



Real GKSystem::stableDt(  const int                 a_step_number,
                          const TimeIntegratorType  a_ti_type )
{
  if ( a_ti_type == ti_explicit ) {
    return m_gk_ops->stableDtExpl( m_state_comp, a_step_number );
  } else if ( a_ti_type == ti_imex ) {
    return m_gk_ops->stableDtImEx( m_state_comp, a_step_number );
  } else {
    /* if fully implicit time integration, then there is not really a maximum
     * stable dt. So using 1000 times the explicit stable dt as a limit, 
     * because systems get really stiff to solve for larger time steps.
     * January 2019: this is not really important right now since we are not
     * considering fully implicit time integration in the near future. */
    return (1000.0*m_gk_ops->stableDtExpl( m_state_comp, a_step_number ));
  }
}


void GKSystem::enforcePositivity( KineticSpeciesPtrVect& a_soln )
{
   for (int n(0); n<a_soln.size(); n++) {
      KineticSpecies* kinetic_species( a_soln[n].operator->() );

      if (m_enforce_positivity_pointwise) {
	// Replaces all negative values with zeroes (does not conserve particles)
	LevelData<FArrayBox>& dfn = kinetic_species->distributionFunction();
	SpaceUtils::enforcePositivity(dfn);
      }
      else {
	// Enforces positivity with particle conservation. However,
	// this requires ghost cells (including corners) to be properly filled
	// Presently, we have NaNs in all ghosts after restart.
	// Thus, this capability does not preserve the restart functionality.
	// Also, due to NaNs in ghosts, it is not clear whether it even works properly
	m_positivity_post_processor.enforce( 
					    kinetic_species->distributionFunction(),
					    kinetic_species->maxValue() );
      }
   }
}


static inline
std::string dirPrefix( const std::string& a_prefix,
                       const std::string& a_diag_name,
                       const int a_sys_id,
                       const std::string& a_sys_id_str )
{
   std::string dir_prefix( a_prefix );
#if 1  // warning, OS dependencies, will not work on all platforms
   std::string iter_str;
   if (a_sys_id < 0) {
     iter_str = ( a_prefix 
                  + "_" 
                  + a_diag_name 
                  + "_plots" );
   } else {
     iter_str = ( a_prefix 
                  + "_" 
                  + a_sys_id_str
                  + "_" 
                  + a_diag_name 
                  + "_plots" );
   }
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



static inline
std::string plotFileName( const std::string& a_prefix,
                          const std::string& a_diag_name,
                          const int a_cur_step,
                          const int a_sys_id,
                          const std::string& a_sys_id_str )
{
   std::string dir_prefix( dirPrefix( a_prefix, 
                                      a_diag_name, 
                                      a_sys_id, 
                                      a_sys_id_str ) );
   char buffer[100];
   sprintf( buffer, "%04d.", a_cur_step );
   std::string filename( dir_prefix + a_prefix + "." + a_diag_name + buffer );
   return filename;
}


static inline
std::string plotFileName( const std::string& a_prefix,
                          const std::string& a_diag_name,
                          const std::string& a_species_name,
                          const int a_cur_step,
                          const int a_species_index,
                          const int a_sys_id,
                          const std::string& a_sys_id_str )
{
   std::string dir_prefix( dirPrefix( a_prefix, 
                                      a_diag_name, 
                                      a_sys_id, 
                                      a_sys_id_str ) );
   char buffer0[10];
   char buffer[20];
   sprintf( buffer0, ".%d.", a_species_index);
   sprintf( buffer, "%04d", a_cur_step );
   std::string filename(    dir_prefix 
                          + a_prefix 
                          + buffer0 
                          + a_species_name 
                          + "." 
                          + a_diag_name 
                          + buffer 
                          + "." );

   return filename;
}


void GKSystem::writePlotFile(const char    *prefix,
                             const int     cur_step,
                             const double& cur_time )
{
   CH_TIME("GKSystem::writePlotFile()");
   // If the efield and potential are fixed, 
   // only consider plotting them at step 0
   if ( !m_gk_ops->fixedEField() || cur_step == 0 ) {
     for (int i = 0; i < m_diagnostics_cfg_field_vars.size(); i++) {

       std::string varname = m_diagnostics_cfg_field_vars[i];
       if (varname == "potential_non_zonal") continue;

       bool non_zonal = false;
       if (varname == "potential") {
        for (int j = 0; j < m_diagnostics_cfg_field_vars.size(); j++) {
          if (m_diagnostics_cfg_field_vars[j] == "potential_non_zonal") {
            non_zonal = true;
          }
        }
       }

       std::string filename( plotFileName( prefix,
                                           varname,
                                           cur_step,
                                           m_sys_id,
                                           m_sys_id_str ) );

        CFG::LevelData<CFG::FArrayBox> var;
        m_diagnostics->getCfgVar(var, varname, non_zonal);
        m_diagnostics->plotCfgVar(var, filename, cur_time);
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

      if (m_diagnostics_dfn_plots.size() > 0) {
        if (cur_step == 0) {
          for (int i = 0; i<m_diagnostics_dfn_plots.size(); i++) {
            std::string plotname = m_diagnostics_dfn_plots[i];
            if (plotname == "dfn") {
              m_diagnostics_dfn_plots.push_back("bstar_par");
              break;
            }
          }
        }
        std::vector<std::string> filenames(0);
        for (int i = 0; i<m_diagnostics_dfn_plots.size(); i++) {
          std::string plotname = m_diagnostics_dfn_plots[i];
          filenames.push_back( plotFileName( prefix,
                                             plotname,
                                             soln_species.name(),
                                             cur_step,
                                             species + 1,
                                             m_sys_id,
                                             m_sys_id_str ) ) ;
        }
        m_diagnostics->writeDfnPlots( m_diagnostics_dfn_plots,
                                      filenames,
                                      soln_species,
                                      cur_time,
                                      radial_index,
                                      poloidal_index,
                                      toroidal_index,
                                      vpar_index,
                                      mu_index );
      }

      for (int i(0); i < m_diagnostics_cfg_vars_kin_spec.size(); i++) {
        std::string varname = m_diagnostics_cfg_vars_kin_spec[i];

        std::string filename( plotFileName( prefix,
                                            varname,
                                            soln_species.name(),
                                            cur_step,
                                            species + 1,
                                            m_sys_id,
                                            m_sys_id_str ) );

        CFG::LevelData<CFG::FArrayBox> var;
        m_diagnostics->getCfgVar(var, varname, soln_species);
        m_diagnostics->plotCfgVar( var, filename, cur_time );

      }

   }

   for (int i(0); i < m_diagnostics_cfg_vars_total_kin_spec.size(); i++) {
     std::string varname = m_diagnostics_cfg_vars_total_kin_spec[i];

     std::string filename( plotFileName( prefix,
                                         varname,
                                         cur_step,
                                         m_sys_id,
                                         m_sys_id_str ) );

     CFG::LevelData<CFG::FArrayBox> var;
     m_diagnostics->getCfgVar(var, varname, kinetic_species);
     m_diagnostics->plotCfgVar( var, filename, cur_time );
   }

   // Fluid species variables
   if (m_diagnostics_fluids) {

      const CFG::FluidSpeciesPtrVect& fluids( m_state_phys.dataFluid() );

      for (int species(0); species<fluids.size(); species++) {

         const CFG::CFGVars& fluid_vars = *(fluids[species]);

         for (int n=0; n<fluid_vars.num_cell_vars(); ++n) {

            CFG::LevelData<CFG::FArrayBox> var;
            m_diagnostics->getFluidCellVar( var, 
                                            fluid_vars, 
                                            fluid_vars.cell_var_name(n) );

            std::string filename;
            if (  fluid_vars.num_cell_vars() == 1 &&
                  fluid_vars.name() == fluid_vars.cell_var_name(0) ) {
               filename = plotFileName( prefix,
                                        fluid_vars.name(),
                                        cur_step,
                                        m_sys_id,
                                        m_sys_id_str );
            } else {
               filename = plotFileName( prefix,
                                        fluid_vars.cell_var_name(n),
                                        fluid_vars.name(),
                                        cur_step,
                                        species + 1,
                                        m_sys_id,
                                        m_sys_id_str );
            }

            m_diagnostics->plotCfgVar( var, filename, cur_time );
         }
         
         for (int n=0; n<fluid_vars.num_face_vars(); ++n) {

            CFG::LevelData<CFG::FArrayBox> var;
            m_diagnostics->getFluidFaceVarAtCell( var, 
                                                  fluid_vars, 
                                                  fluid_vars.face_var_name(n) );

            std::string filename (plotFileName( prefix,
                                                fluid_vars.face_var_name(n),
                                                fluid_vars.name(),
                                                cur_step,
                                                species + 1,
                                                m_sys_id,
                                                m_sys_id_str ));

            m_diagnostics->plotCfgVar( var, filename, cur_time );
         }
         
         for (int n=0; n<fluid_vars.num_edge_vars(); ++n) {

            CFG::LevelData<CFG::FArrayBox> var;
            m_diagnostics->getFluidEdgeVarAtCell( var, 
                                                  fluid_vars, 
                                                  fluid_vars.edge_var_name(n) );

            std::string filename (plotFileName( prefix,
                                                fluid_vars.edge_var_name(n),
                                                fluid_vars.name(),
                                                cur_step,
                                                species + 1,
                                                m_sys_id,
                                                m_sys_id_str ));

            m_diagnostics->plotCfgVar( var, filename, cur_time );
         }

         for (int n=0; n<fluid_vars.num_node_vars(); ++n) {

           CFG::LevelData<CFG::FArrayBox> var;
           m_diagnostics->getFluidNodeVarAtCell( var,
                                                 fluid_vars,
                                                 fluid_vars.node_var_name(n) );

           std::string filename (plotFileName( prefix,
                                               fluid_vars.node_var_name(n),
                                               fluid_vars.name(),
                                               cur_step,
                                               species + 1,
                                               m_sys_id,
                                               m_sys_id_str ));

           m_diagnostics->plotCfgVar( var, filename, cur_time );

         }
         
         if ( fluid_vars.m_plotMemberVars == 1 ) {
            std::vector<string> varname = fluid_vars.m_plotMemberVarNames;
            for (int i(0); i<varname.size(); i++) {
               std::string filename (plotFileName( prefix,
                                                   varname[i],
                                                   fluid_vars.name(),
                                                   cur_step,
                                                   species + 1,
                                                   m_sys_id,
                                                   m_sys_id_str ));
               CFG::LevelData<CFG::FArrayBox> var;
               m_diagnostics->getFluidOpMember(  var, fluid_vars, varname[i] );
               m_diagnostics->plotCfgVar( var, filename, cur_time );
            }
         }
      }
   }
}

      

void GKSystem::writeCheckpointFile( const std::string& a_chkpt_prefix,
                                    const int          a_cur_step,
                                    const double       a_cur_time,
                                    const double       a_cur_dt )
{
   CH_TIME("GKSystem::writeCheckpointFile");

#ifdef CH_USE_HDF5
   MPI_Barrier(MPI_COMM_WORLD);

   char iter_str[100];
   if (m_sys_id < 0) {
     sprintf( iter_str, "%s%04d.%dd.hdf5",
              a_chkpt_prefix.c_str(), 
              a_cur_step, 
              SpaceDim );
   } else {
     sprintf( iter_str, "%d_%s%04d.%dd.hdf5",
              m_sys_id,
              a_chkpt_prefix.c_str(), 
              a_cur_step, 
              SpaceDim );
   }

   if (procID()==0) {
      if (m_sys_id >= 0) cout << "System " << m_sys_id << ": ";
      cout  << "Writing checkpoint file " 
            << iter_str
            << endl;
   }


   HDF5Handle handle( iter_str, HDF5Handle::CREATE );

   HDF5HeaderData header;
   const int RESTART_VERSION(2); // to distinguish from future versions
   header.m_int ["cur_step"]        = a_cur_step;
   header.m_real["cur_time"]        = a_cur_time;
   header.m_real["cur_dt"]          = a_cur_dt;
   header.m_real["Er_lo"]           = m_gk_ops->getLoRadialField(m_state_comp);
   header.m_real["Er_hi"]           = m_gk_ops->getHiRadialField(m_state_comp);
   header.m_int ["restart_version"] = RESTART_VERSION;
   header.writeToFile( handle );

   if ( m_gk_ops->usingAmpereLaw() ) {
      // save the averaged radial E-field in 4D injected form
      LevelData<FArrayBox> Er_cell_injected;
      m_phase_geom->injectConfigurationToPhase( m_gk_ops->getErAverageCell(), Er_cell_injected);
      handle.setGroup( "Er_cell" );
      write(handle,Er_cell_injected.boxLayout());
      write(handle,Er_cell_injected,"data",Er_cell_injected.ghostVect());
 
      LevelData<FluxBox> Er_face_injected;
      m_phase_geom->injectConfigurationToPhase( m_gk_ops->getErAverageFace(), Er_face_injected);
      handle.setGroup( "Er_face" );
      write(handle,Er_face_injected.boxLayout());
      write(handle,Er_face_injected,"data",Er_face_injected.ghostVect());

      LevelData<FArrayBox> E_tilde_cell_injected;
      m_phase_geom->injectConfigurationToPhase( m_gk_ops->getETildeCell(), E_tilde_cell_injected);
      handle.setGroup( "E_tilde_cell" );
      write(handle,E_tilde_cell_injected.boxLayout());
      write(handle,E_tilde_cell_injected,"data",E_tilde_cell_injected.ghostVect());
 
      LevelData<FluxBox> E_tilde_face_injected;
      m_phase_geom->injectConfigurationToPhase( m_gk_ops->getETildeFace(), E_tilde_face_injected);
      handle.setGroup( "E_tilde_face" );
      write(handle,E_tilde_face_injected.boxLayout());
      write(handle,E_tilde_face_injected,"data",E_tilde_face_injected.ghostVect());
   }
   
   if ( m_old_vorticity_model ) {
      LevelData<FArrayBox> phi_injected;
      m_phase_geom->injectConfigurationToPhase( m_gk_ops->getPhi(), phi_injected);
      handle.setGroup( "potential" );
      write(handle,phi_injected.boxLayout());
      write(handle,phi_injected,"data",phi_injected.ghostVect());
   }
   
   const KineticSpeciesPtrVect& kinetic_species( m_state_comp.dataKinetic() );
   for (int species(0); species<kinetic_species.size(); species++) {

      // Get solution distribution function for the current species
      KineticSpecies& soln_species( *(kinetic_species[species]) );
      LevelData<FArrayBox> & soln_dfn = soln_species.distributionFunction();
      char buff[100];
      sprintf( buff, "dfn_%d", species + 1 );
      handle.setGroup( buff );
      write( handle, soln_dfn.boxLayout() );
      write( handle, soln_dfn, "data" );
   }

   const CFG::FluidSpeciesPtrVect& fluid_species( m_state_comp.dataFluid() );
   for (int species(0); species<fluid_species.size(); species++) {

      CFG::CFGVars& this_fluid_species( *(fluid_species[species]) );

      char buff[100];
      for (int n=0; n<this_fluid_species.num_cell_vars(); ++n) {
      
         CFG::LevelData<CFG::FArrayBox> & fluid_data = this_fluid_species.cell_var(n);
         const string& var_name = this_fluid_species.cell_var_name(n);
         sprintf( buff, "fluid_%d_%s", species + 1, var_name.c_str() );
         handle.setGroup( buff );

         LevelData<FArrayBox> injected_fluid_data;
         m_phase_geom->injectConfigurationToPhase( fluid_data, injected_fluid_data);

         write( handle, injected_fluid_data.boxLayout() );
         write( handle, injected_fluid_data, "data", injected_fluid_data.ghostVect() );
      }

      // NEED TO ADD FACE DATA
   }

   MPI_Barrier(MPI_COMM_WORLD);
   if (procID()==0) {
      if (m_sys_id >= 0) cout << "System " << m_sys_id << ": ";
      cout << "Writing history file" << endl;
   }

   m_gk_ops->writeCheckpointFile( handle );

   handle.close();
   
   //Save geometry data initialization
   if (m_mag_geom->isWritingInitializaionData()) {
      m_mag_geom->writeInitializationData();
   }
   
#else
   MayDay::Error( "restart only defined with hdf5" );
#endif
}


void GKSystem::readCheckpointFile( const std::string& a_chkpt_fname,
                                   int&               a_cur_step,
                                   double&            a_cur_time,
                                   double&            a_cur_dt )
{
   string fname;
   if (m_sys_id >= 0) {
     fname = m_sys_id_str + "_"  + a_chkpt_fname;
   } else {
     fname = a_chkpt_fname;
   }

#ifdef CH_USE_HDF5
   HDF5Handle handle( fname, HDF5Handle::OPEN_RDONLY );

   HDF5HeaderData header;
   header.readFromFile( handle );

   a_cur_step = header.m_int ["cur_step"];
   a_cur_time = header.m_real["cur_time"];
   a_cur_dt   = header.m_real["cur_dt"];
   m_gk_ops->setLoRadialField(m_state_comp, header.m_real["Er_lo"]);
   m_gk_ops->setHiRadialField(m_state_comp, header.m_real["Er_hi"]);
   //   const int restart_version = header.m_int ["restart_version"];

   if ( m_gk_ops->usingAmpereLaw() ) {
      handle.setGroup("Er_cell");
      CFG::LevelData<CFG::FArrayBox> Er_cell_config(m_mag_geom->grids(), 3, CFG::IntVect::Unit);
      LevelData<FArrayBox> Er_cell_injected;
      m_phase_geom->injectConfigurationToPhase(Er_cell_config, Er_cell_injected);
      read( handle, Er_cell_injected, "data", Er_cell_injected.disjointBoxLayout() );
      m_gk_ops->setErAverage( Er_cell_injected );
    
      handle.setGroup("Er_face");
      CFG::LevelData<CFG::FluxBox> Er_face_config(m_mag_geom->grids(), 3, CFG::IntVect::Unit);
      LevelData<FluxBox> Er_face_injected;
      m_phase_geom->injectConfigurationToPhase(Er_face_config, Er_face_injected);
      read( handle, Er_face_injected, "data", Er_face_injected.disjointBoxLayout() );
      m_gk_ops->setErAverage( Er_face_injected );

      handle.setGroup("E_tilde_cell");
      CFG::LevelData<CFG::FArrayBox> E_tilde_cell_config(m_mag_geom->grids(), 3, CFG::IntVect::Unit);
      LevelData<FArrayBox> E_tilde_cell_injected;
      m_phase_geom->injectConfigurationToPhase(E_tilde_cell_config, E_tilde_cell_injected);
      read( handle, E_tilde_cell_injected, "data", E_tilde_cell_injected.disjointBoxLayout() );
      m_gk_ops->setETilde( E_tilde_cell_injected );
    
      handle.setGroup("E_tilde_face");
      CFG::LevelData<CFG::FluxBox> E_tilde_face_config(m_mag_geom->grids(), 3, CFG::IntVect::Unit);
      LevelData<FluxBox> E_tilde_face_injected;
      m_phase_geom->injectConfigurationToPhase(E_tilde_face_config, E_tilde_face_injected);
      read( handle, E_tilde_face_injected, "data", E_tilde_face_injected.disjointBoxLayout() );
      m_gk_ops->setETilde( E_tilde_face_injected );
   }
   
   if ( m_old_vorticity_model ) {
      handle.setGroup("potential");
      CFG::LevelData<CFG::FArrayBox> phi(m_mag_geom->grids(), 1, CFG::IntVect::Zero);
      LevelData<FArrayBox> phi_injected;
      m_phase_geom->injectConfigurationToPhase(phi, phi_injected);
      read( handle, phi_injected, "data", phi_injected.disjointBoxLayout() );
      m_gk_ops->setInitialPhi( phi_injected );
   }

   KineticSpeciesPtrVect& kinetic_species( m_state_comp.dataKinetic() );
   for (int species(0); species<kinetic_species.size(); species++) {

      // Get solution distribution function for the current species
      KineticSpecies& soln_species( *(kinetic_species[species]) );
      LevelData<FArrayBox>& soln_dfn = soln_species.distributionFunction();
      char buff[100];
      sprintf( buff, "dfn_%d", species + 1 );
      handle.setGroup( buff );
      read( handle, soln_dfn, "data", soln_dfn.disjointBoxLayout(), Interval(0,soln_dfn.nComp()-1), false );
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
         handle.setGroup( buff );
         read( handle, injected_fluid_data, "data", injected_fluid_data.disjointBoxLayout() );
         m_phase_geom->projectPhaseToConfiguration( injected_fluid_data, 
                                                    fluid_data);
      }

      // NEED TO ADD FACE DATA
   }

   m_gk_ops->readCheckpointFile( handle, a_cur_step );
   handle.close();
   
#else
   MayDay::Error("restart only defined with hdf5");
#endif
}


void GKSystem::setupFieldHistories()

{
   ParmParse ppsim( "simulation" );
   m_gk_ops->setupFieldHistories( ppsim );
}


void GKSystem::writeFieldHistory( int     cur_step, 
                                  double  cur_time, 
                                  bool    startup_flag)
{
  if (m_sys_id < 0) {
    m_gk_ops->writeFieldHistory( cur_step, cur_time, startup_flag );
  } else {
    string prefix = m_sys_id_str + "_";
    m_gk_ops->writeFieldHistory( cur_step, cur_time, startup_flag, prefix );
  }
}


void GKSystem::getDiagnostics(  std::vector<Real>& a_minvals,
                                std::vector<Real>& a_maxvals )
{
   const KineticSpeciesPtrVect& kinetic_species( m_state_comp.dataKinetic() );
   for (int n(0); n<kinetic_species.size(); n++) {
      const KineticSpecies& species( *(kinetic_species[n]) );
      a_maxvals.push_back( species.maxValue() );
      a_minvals.push_back( species.minValue() );
   }
   return;
}

void GKSystem::printDiagnostics()
{
   std::vector<Real> a_minvals(0), a_maxvals(0);
   getDiagnostics(a_minvals, a_maxvals);
   CH_TIME("GKSystem::printDiagnostics()");
   if (procID()==0) {
      cout << "  ";
      if (m_sys_id >= 0) cout << "System " << m_sys_id << ": ";
      cout << "Distribution Function Extrema:" << std::endl;
   }
   const KineticSpeciesPtrVect& kinetic_species( m_state_comp.dataKinetic() );
   for (int n(0); n<kinetic_species.size(); n++) {
      if (procID()==0) {
         printf("    Species %d:  %5.3e [max], %5.3e [min]\n",
                n, a_maxvals[n], a_minvals[n] );
      }
   }
   if (procID()==0) printf("\n");
}

void printTimeStep( const Real& a_dt,
                    const std::string& a_name,
                    const Real& a_dt_current  )
{
   if (a_dt>= 0) {
     printf("    %-20s: %1.3e (time scale), %1.4f (CFL)\n", 
            a_name.c_str(), a_dt, a_dt_current/a_dt);
   }
}

void GKSystem::preTimeStep( const ODEVector&  a_vec, 
                            int               a_cur_step, 
                            Real              a_cur_time  )
{
   CH_TIME("GKSystem::preTimeStep()");
   copyStateFromArray( a_vec );
   m_gk_ops->convertToPhysical( m_state_comp, m_state_phys, a_cur_time ); 
   m_gk_ops->preTimeStep( a_cur_step, 
                          a_cur_time, 
                          m_state_comp, 
                          m_state_phys );
}

void GKSystem::postTimeStep(ODEVector&  a_vec,
                            int         a_cur_step, 
                            Real        a_dt, 
                            Real        a_cur_time )
{
  CH_TIME("GKSystem::postTimeStep");

  copyStateFromArray( a_vec );

  if (m_enforce_step_positivity) {
     enforcePositivity( m_state_comp.dataKinetic() );
  }
  
  if (m_enforce_step_floor) {
     m_floor_post_processor.enforce( m_state_comp.dataFluid() );
  }

  m_gk_ops->postTimeStep( a_cur_step, a_dt, a_cur_time, m_state_comp );
  m_gk_ops->convertToPhysical( m_state_comp, m_state_phys, a_cur_time );

  m_dt_vlasov = m_gk_ops->dtScaleVlasov( m_state_comp, a_cur_step );
  m_dt_collisions = m_gk_ops->dtScaleCollisions( m_state_comp, a_cur_step);
  m_dt_transport = m_gk_ops->dtScaleTransport( m_state_comp, a_cur_step );
  m_dt_neutrals = m_gk_ops->dtScaleNeutrals( m_state_comp, a_cur_step );

  if ((procID() == 0) && (m_sys_id < 0)) {
    printf("  dt = %f\n", a_dt);
    for (int i = 0; i < m_dt_vlasov.size(); i++) {
      std::string name_str = "Vlasov (" + m_dt_vlasov[i].first + ")"; 
      printTimeStep( m_dt_vlasov[i].second, name_str, a_dt ); 
    }
    printTimeStep( m_dt_collisions, "Collisions", a_dt ); 
    printTimeStep( m_dt_transport, "Transport", a_dt ); 
    printTimeStep( m_dt_neutrals, "Neutrals", a_dt ); 
  }

  copyStateToArray( a_vec );
}

void GKSystem::postTimeStage( const int   a_cur_step, 
                              const Real  a_dt, 
                              const Real  a_cur_time,
                              const int   a_stage )
{
   m_gk_ops->postTimeStage( a_cur_step, 
                            a_cur_time, 
                            a_dt, 
                            m_state_comp, 
                            a_stage );
}

inline void GKSystem::printParameters() const
{
   if (procID() == 0 && m_verbosity) {

      if (m_sys_id >= 0) {
        printf("System %d:", m_sys_id);
      }
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

      cout  << "enforce_positivity = " 
            << (m_enforce_step_positivity||m_enforce_stage_positivity) 
            << endl;
      std::string ptype("stage");
      if (m_enforce_step_positivity)
         ptype = "step";
      cout << "enforce_positivity_type = " << ptype << endl;
      if (m_sys_id >= 0) {
        printf("--");
      }

   }
}


void GKSystem::setParameters( const GKSystemParameters& a_params )
{   
   m_kinetic_ghosts = a_params.kineticGhosts();
   m_fluid_ghosts = a_params.fluidGhosts();

   m_verbosity = a_params.verbosity();
   m_mag_geom_type = a_params.magGeomType();

   if ( m_mag_geom_type != "SingleNull" && m_mag_geom_type != "SNCore") {
      m_num_cells = a_params.numCells();
      m_is_periodic = a_params.isPeriodic();

      m_configuration_decomposition = a_params.cfgDecomp();
      m_phase_decomposition = a_params.phaseDecomp();
   } else {
      m_num_velocity_cells = a_params.numVelocityCells();
   }

   m_velocity_decomposition = a_params.velDecomp();

   m_ti_class = a_params.tiClass();
   m_ti_method = a_params.tiMethod();

   m_diagnostics_dfn_plots = a_params.diagnosticsDfnPlots();
   m_diagnostics_cfg_field_vars = a_params.diagnosticsCfgFieldVars();
   m_diagnostics_cfg_vars_kin_spec = a_params.diagnosticsCfgVarsKinSpec();
   m_diagnostics_cfg_vars_total_kin_spec = a_params.diagnosticsCfgVarsTotalKinSpec();
   m_diagnostics_fluids = a_params.diagnosticsFluids();

   m_fixed_plotindices = a_params.fixedPlotIndices();

   if ( m_mag_geom_type != "SingleNull" && m_mag_geom_type != "SNCore" ) {
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
   } else {
   }

   m_enforce_stage_positivity = a_params.enforceStagePositivity();
   m_enforce_step_positivity = a_params.enforceStepPositivity();
   m_enforce_positivity_pointwise = a_params.enforcePositivityPointwise();
   if (m_enforce_stage_positivity || m_enforce_step_positivity) {

      int n_iter = a_params.positivityNIter();
      bool verbose = a_params.positivityVerbose();

      int width(2);
      if (m_enforce_step_positivity) width++;
      IntVect halo( width*IntVect::Unit );
      m_positivity_post_processor.define( halo, n_iter, verbose );
   }


   m_enforce_step_floor = a_params.enforceStepFloor();
   if (m_enforce_step_floor) {
      Real floor_value = a_params.floorValue();
      bool absolute_floor = a_params.absoluteFloor();
      m_floor_post_processor.define( floor_value, absolute_floor );
   }
   
   m_old_vorticity_model = a_params.oldVorticityModel();

   return;
}


void GKSystem::giveSpeciesTheirGyroaverageOps()
{
  std::map<std::string,GyroaverageOperator*>& 
                gyroavg_ops = m_gk_ops->getGyroaverageOps();

  KineticSpeciesPtrVect& species_comp(m_state_comp.dataKinetic());
  KineticSpeciesPtrVect& species_phys(m_state_phys.dataKinetic());

  for (int s(0); s < species_comp.size(); s++) {
    KineticSpecies& sp_comp = *(species_comp[s]);
    KineticSpecies& sp_phys = *(species_phys[s]);

    if (sp_comp.isGyrokinetic()) sp_comp.gyroaverageOp(gyroavg_ops[sp_comp.name()]);
    if (sp_phys.isGyrokinetic()) sp_phys.gyroaverageOp(gyroavg_ops[sp_phys.name()]);
  }

  return;
}

#include "NamespaceFooter.H"
