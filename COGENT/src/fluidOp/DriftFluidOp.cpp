#include "DriftFluidOp.H"
#include "EllipticOpBCFactory.H"
#include "FluidVarBCFactory.H"
#include "FourthOrderUtil.H"
#include "MagBlockCoordSys.H"
#include "SpaceUtils.H.multidim"
#include "CONSTANTS.H"

#include "NodeFArrayBox.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "FluidOpF_F.H"

#include "NamespaceHeader.H" 

const char* DriftFluidOp::pp_name = {"DriftFluid_op"};


DriftFluidOp::DriftFluidOp(const string&   a_pp_str,
                           const string&   a_species_name,
                           const MagGeom&  a_geometry,
                           const Real      a_larmor,
                           const int       a_verbosity )
   : m_verbosity(a_verbosity),
     m_larmor(a_larmor),
     m_geometry(a_geometry),
     m_is_time_implicit(false),
     m_first_call(true),
     m_is_linearized(true),
     m_include_hyperviscosity(false),
     m_use_hyperviscosity_bcs(false),
     m_include_boundary_relaxation(false),
     m_advection_scheme("uw3"),
     m_harmonic_filtering(SpaceDim, 0),
     m_species_name(a_species_name),
     m_opt_string(a_pp_str),
     m_my_pc_idx(-1),
     m_max_dt(DBL_MAX),
     m_time_step_diagnostics(false),
     m_n0_func(NULL)

{
   ParmParse pp(a_pp_str.c_str());
   parseParameters( pp );
   if (m_verbosity>0) {
      printParameters();
   }

   // Create FluidOpUtils class that handles hyperviscosity and boundary buffers
   m_fluid_op_utils = new FluidOpUtils(pp, m_geometry, m_larmor, m_verbosity);
   
   string variable_name = "density";

   // Input the initial conditions
   parseInitialConditions(a_species_name, variable_name);

   // Input the boundary conditions
   FluidVarBCFactory fluid_var_bc_factory;
   m_fluid_variable_bc = fluid_var_bc_factory.create(a_species_name,
                                                     variable_name,
                                                     a_geometry.getCoordSys()->type(),
                                                     false);

      
   // Input the hyperviscosity parameters and bcs
   if (m_include_hyperviscosity) {

      const std::string prefixOp( a_pp_str + ".hyperviscosity");
      ParmParse ppOp( prefixOp.c_str() );

      const std::string name( m_species_name + ".hyperviscosity");
      const std::string prefixBC( "BC." + name);
      ParmParse ppOpBC( prefixBC.c_str() );

      EllipticOpBCFactory elliptic_op_bc_factory;
      m_hyperviscosity_op_bcs = RefCountedPtr<EllipticOpBC>(elliptic_op_bc_factory.create(name,
                                                                                          ppOpBC,
                                                                                          *(a_geometry.getCoordSys()),
                                                                                          false ));
      
      m_fluid_op_utils->initializeHyperViscosityOp(ppOp,
                                                   *m_hyperviscosity_op_bcs,
                                                   m_use_hyperviscosity_bcs);
      
   }
   
   // Initialize boundary buffers
   if (m_include_boundary_relaxation) {
      m_fluid_op_utils->initializeBoundaryBuffers(pp, m_include_boundary_relaxation);
   }
}


DriftFluidOp::~DriftFluidOp()
{
   if (m_fluid_variable_bc) delete m_fluid_variable_bc;
   if (m_fluid_op_utils) delete m_fluid_op_utils;
}


void DriftFluidOp::accumulateRHS(FluidSpeciesPtrVect&              a_rhs,
                                const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                const FluidSpeciesPtrVect&         a_fluid_species,
                                const PS::ScalarPtrVect&           a_scalars,
                                const EMFields&                    a_EM_fields,
                                const int                          a_fluidVecComp,
                                const Real                         a_time)
{
   /*
    This class computes the RHS for :
      Linearized model
         Physics model: d(n)/dt = -1/c * ExB/B^2 * nabla(n0)
         COGENT equation: d(J*n)/dt = -J * Larmor * ExB/B^2 * nabla(n0)
      Nonlinear model
         Physics model: d(n)/dt = -1/c div(n * ExB/B^2)
         COGENT equation: d(J*n)/dt = -J * Larmor div(n * ExB/B^2)
   */
   
   // Get fluid rhs
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   LevelData<FArrayBox>& rhs_data( rhs_fluid.cell_var(0) );

   // Get fluid soln
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
   const LevelData<FArrayBox>& soln_data( soln_fluid.cell_var(0) );


   // Get geometry parameters
   const DisjointBoxLayout& grids( soln_data.getBoxes() );

   // Add advection term
   addAdvection(rhs_data, soln_data, a_EM_fields, a_time);

   // Add hyperviscosity
   if (m_include_hyperviscosity) {
      m_fluid_op_utils->addHyperViscosity(rhs_data,
                                          soln_data,
                                          true);
   }
   
   // Add relaxation to a reference solution in boundary buffers
   if (m_include_boundary_relaxation) {
       m_fluid_op_utils->addRelaxationInBoundaryBuffers(rhs_data,
                                                        soln_data,
                                                        true);
    }
   
   // We keep a local copy of Efield to be used in CFL calculations
   // This is a temporary approach;
   m_E_field.define(a_EM_fields.getEFieldFace());

   m_first_call = false;
}


void DriftFluidOp::accumulateExplicitRHS(FluidSpeciesPtrVect&               a_rhs,
                                        const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                        const FluidSpeciesPtrVect&         a_fluid_species,
                                        const PS::ScalarPtrVect&           a_scalars,
                                        const EMFields&                    a_EM_fields,
                                        const int                          a_fluidVecComp,
                                        const Real                         a_time)
{
   if (!m_is_time_implicit) {
      accumulateRHS(a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_EM_fields, a_fluidVecComp, a_time);
   }
}


void DriftFluidOp::accumulateImplicitRHS(FluidSpeciesPtrVect&              a_rhs,
                                        const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                        const FluidSpeciesPtrVect&        a_fluid_species,
                                        const PS::ScalarPtrVect&          a_scalars,
                                        const EMFields&                   a_EM_fields,
                                        const int                         a_fluidVecComp,
                                        const Real                        a_time )
{
   if (m_is_time_implicit) {
      accumulateRHS(a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_EM_fields, a_fluidVecComp, a_time);
   }
}

void DriftFluidOp::addAdvection(LevelData<FArrayBox>&       a_rhs,
                                const LevelData<FArrayBox>& a_soln,
                                const EMFields&             a_EM_fields,
                                const Real                  a_time)
{
   // Get geometry parameters
   const DisjointBoxLayout& grids = m_geometry.grids();

   // Initialize background density
   if (m_is_linearized && m_first_call) {
     m_n0.define(grids, 1, 2*IntVect::Unit);
     m_n0_func->assign( m_n0, m_geometry, a_time);
   }

   // Compute physical ExB advection term = div[V_ExB * n]
   LevelData<FArrayBox> ExB_advection(grids, 1, IntVect::Zero);

   if (m_is_linearized) {
     m_fluid_op_utils->computeExBAdvection(ExB_advection,
                                           m_n0,
                                           a_EM_fields.getEFieldFace(),
                                           m_advection_scheme,
                                           false);
   }
   else {
     m_fluid_op_utils->computeExBAdvection(ExB_advection,
                                           a_soln,
                                           a_EM_fields.getEFieldFace(),
                                           m_advection_scheme,
                                           false);
   }

   // Get mapped advection
   m_geometry.multJonValid(ExB_advection);

   // Update RHS
   for (DataIterator dit(grids); dit.ok(); ++dit) {
     a_rhs[dit] -= ExB_advection[dit];
   }
}


void DriftFluidOp::subtractExBGradN(LevelData<FArrayBox>&  a_rhs,
                                    const LevelData<FArrayBox>& a_soln,
                                    const EMFields&        a_EM_fields,
                                    const Real             a_time)
{
      
   // Get geometry parameters
   const DisjointBoxLayout& grids = m_geometry.grids();
   int order = (m_geometry.secondOrder()) ? 2 : 4;
   
   LevelData<FArrayBox> ExB_grad_n(grids, 1, IntVect::Zero);
   m_fluid_op_utils->computeExBGradData(ExB_grad_n,
                                        a_soln,
                                        a_EM_fields.getEFieldCell(),
                                        "uw2",
                                        false);
   
   // Get mapped ExBgradN
   m_geometry.multJonValid(ExB_grad_n);

   // Update RHS
   for (DataIterator dit(grids); dit.ok(); ++dit) {
     a_rhs[dit] -= ExB_grad_n[dit];
   }
}

void DriftFluidOp::subtractExBGradN0(LevelData<FArrayBox>&  a_rhs,
                                     const EMFields&        a_EM_fields,
                                     const Real             a_time)
{
   /*
    Because this term only represets a "source term" for the density equation
    we compute gradN0 using centered scheme without worrying about stability
    */
   
   // Get geometry parameters
   const DisjointBoxLayout& grids = m_geometry.grids();
   int order = (m_geometry.secondOrder()) ? 2 : 4;
   
   // Get background density gradient
   if (m_first_call) {
      LevelData<FArrayBox> n0(grids, 1, 2*IntVect::Unit);
      m_n0_func->assign( n0, m_geometry, a_time);
      
      LevelData<FArrayBox> n0_gradient_mapped(grids, CFG_DIM, IntVect::Unit);
      if (SpaceDim == 3) {
         m_geometry.computeMappedGradient(n0, n0_gradient_mapped, order);
      }
      if (SpaceDim == 2) {
         m_geometry.computeMappedPoloidalGradientWithGhosts(n0, n0_gradient_mapped, order);
      }

      m_n0_gradient.define(grids, CFG_DIM, IntVect::Unit);
      if (SpaceDim == 3) {
         m_geometry.unmapGradient(n0_gradient_mapped, m_n0_gradient );
      }
      if (SpaceDim == 2) {
         m_geometry.unmapPoloidalGradient(n0_gradient_mapped, m_n0_gradient );
      }
   }
   
   // Compute the ExB drift on cell centers
   LevelData<FArrayBox> ExB_drift(grids, 3, IntVect::Zero);
   m_geometry.computeEXBDrift(a_EM_fields.getEFieldCell(), ExB_drift);
   
   LevelData<FArrayBox> ExB_grad_n0(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      ExB_grad_n0[dit].setVal(0.);
      FArrayBox tmp(grids[dit],1);
      
      if (SpaceDim ==3) {
         for (int n=0; n<SpaceDim; n++) {
            tmp.copy(ExB_drift[dit],n,0,1);
            tmp.mult(m_n0_gradient[dit],n,0,1);
            ExB_grad_n0[dit].plus(tmp);
         }
      }
      
      if (SpaceDim == 2) {
         tmp.copy(ExB_drift[dit],0,0,1);
         tmp.mult(m_n0_gradient[dit],0,0,1);
         ExB_grad_n0[dit].plus(tmp);
         
         tmp.copy(ExB_drift[dit],2,0,1);
         tmp.mult(m_n0_gradient[dit],1,0,1);
         ExB_grad_n0[dit].plus(tmp);
      }
      
      ExB_grad_n0[dit].mult(m_larmor);
   }
   
   m_geometry.multJonValid(ExB_grad_n0);

   // Update RHS
   for (DataIterator dit(grids); dit.ok(); ++dit) {
     a_rhs[dit] -= ExB_grad_n0[dit];
   }
}

void DriftFluidOp::postTimeEval(FluidSpecies&  a_species_comp,
                                FluidSpecies&  a_species_phys,
                                const Real     a_dt,
                                const Real     a_time,
                                const int      a_stage )
{
   LevelData<FArrayBox>& soln_comp = a_species_comp.cell_var(0);
   LevelData<FArrayBox>& soln_phys = a_species_phys.cell_var(0);
   
   for (int dir=0; dir<SpaceDim; dir++) {
      if (m_harmonic_filtering[dir] == 1) {
         SpaceUtils::applyHarmonicFiltering(soln_comp, dir);
         SpaceUtils::applyHarmonicFiltering(soln_phys, dir);
      }
   }
}

void DriftFluidOp::defineBlockPC( std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>& a_pc,
                                 std::vector<PS::DOFList>&                                    a_dof_list,
                                 const PS::ODEVector&                                         a_soln_vec,
                                 void*                                                        a_gkops,
                                 const std::string&                                           a_out_string,
                                 const std::string&                                           a_opt_string,
                                 bool                                                         a_im,
                                 const FluidSpecies&                                          a_fluid_species,
                                 const PS::GlobalDOFFluidSpecies&                             a_global_dofs,
                                 const int                                                    a_species_idx,
                                 const int                                                    a_id )
{
  if (a_im && m_is_time_implicit) {
    CH_assert(a_pc.size() == a_dof_list.size());
  
    if (!procID()) {
      std::cout << "  Fluid Species " << a_species_idx
                << " : "
                << " creating " << _FLUID_OP_PC_ << " preconditioner"
                << " (index = " << a_pc.size() << ").\n";
    }
  
    PS::Preconditioner<PS::ODEVector,PS::AppCtxt> *pc;
    pc = new PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt>;
    pc->setSysID(a_id);
    dynamic_cast<PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt>*>
      (pc)->define(a_soln_vec, a_gkops, *this, m_opt_string, m_opt_string, a_im);
    dynamic_cast<PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt>*>
      (pc)->speciesIndex(a_species_idx);
  
    PS::DOFList dof_list(0);
  
    for (int var=0; var<a_fluid_species.num_cell_vars(); ++var) {

       const LevelData<FArrayBox>& soln_data   (a_fluid_species.cell_var(var));
       const DisjointBoxLayout&    grids       (soln_data.disjointBoxLayout());
       const int                   n_comp      (soln_data.nComp());
       const LevelData<FArrayBox>& pMapping    (a_global_dofs.data(var)); 
  
       for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
          const Box& grid = grids[dit];
          const FArrayBox& pMap = pMapping[dit];
          for (BoxIterator bit(grid); bit.ok(); ++bit) {
             IntVect ic = bit();
             for (int n(0); n < n_comp; n++) {
                dof_list.push_back((int) pMap.get(ic ,n) - a_global_dofs.mpiOffset());
             }
          }
       }
    }
  
    m_my_pc_idx = a_pc.size();

    a_pc.push_back(pc);
    a_dof_list.push_back(dof_list);
  }

  return;
}


void DriftFluidOp::updateBlockPC(std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>&  a_pc,
                                const PS::KineticSpeciesPtrVect&                              a_kin_species_phys,
                                const FluidSpeciesPtrVect&                                    a_fluid_species,
                                const Real                                                    a_time,
                                const int                                                     a_step,
                                const int                                                     a_stage,
                                const Real                                                    a_shift,
                                const bool                                                    a_im,
                                const int                                                     a_species_idx )
{
  if (a_im && m_is_time_implicit) {
    CH_assert(m_my_pc_idx >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx);
  
    if (!procID()) {
      std::cout << "    ==> Updating " << _FLUID_OP_PC_ << " preconditioner " 
                << " for IdealMhDriftFluidf fluid species " << a_species_idx << ".\n";
    }
  
    PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt> *pc 
      = dynamic_cast<PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt>*>(a_pc[m_my_pc_idx]);
    CH_assert(pc != NULL);
    pc->update( a_kin_species_phys, 
                a_fluid_species, 
                a_time, 
                a_step,
                a_stage,
                a_shift, 
                a_im, 
                a_species_idx);
  }
  return;
}


void DriftFluidOp::updatePCImEx(const FluidSpeciesPtrVect&       a_fluid_species,
                               const PS::KineticSpeciesPtrVect& a_kinetic_species,
                               const double                     a_time,
                               const int                        a_step,
                               const int                        a_stage,
                               const double                     a_shift,
                               const int                        a_component,
                               const std::string& )
{
   CH_TIME("DriftFluidOp::updatePCImEx()");
   
   // Update preconditioner for implicit terms here
   
}


void DriftFluidOp::solvePCImEx(FluidSpeciesPtrVect&              a_fluid_species_solution,
                              const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                              const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                              const std::string&,
                              const int                         a_component )
{
   CH_TIME("DriftFluidOp::solvePCImEx()");
   if (!m_is_time_implicit) {
     MayDay::Error("Error in DriftFluidOp::solvePCImEx(): m_is_time_implicit is false. This function shouldn't have been called at all!");
   }
   
   // Solve preconditioner operator for implicit terms here

}

void DriftFluidOp::fillGhostCells(FluidSpecies&  a_species_phys,
                                  const double   a_time )
{
   CH_TIME("DriftFluidOp::fillGhostCells()");
   
   // Do we need it here, or we alredy pre-filled the ghosts
   // Definitely don't need for the ExB*gradn0 case
   
   // Fill ghost cells except for those on physical boundaries
   LevelData<FArrayBox>& fld( a_species_phys.cell_var(0) );
   m_geometry.fillInternalGhosts( fld );

   // Fill ghost cells on physical boundaries
   m_fluid_variable_bc->apply( a_species_phys, a_time );
}


Real DriftFluidOp::computeDtImExTI( const FluidSpeciesPtrVect&  a_fluid_comp )
{
   CH_TIME("DriftFluidOp::computeDtImExTI()");
   
   Real dt = computeCourantTimeStep( a_fluid_comp );
   
   if (dt>m_max_dt) dt=m_max_dt;
   
   return dt;
}

Real DriftFluidOp::computeDtExplicitTI( const FluidSpeciesPtrVect&  a_fluid_comp )
{
   CH_TIME("DriftFluidOp::computeDtExplicit()");
   
   Real dt = computeCourantTimeStep( a_fluid_comp );
   
   if (dt>m_max_dt) dt=m_max_dt;
   
   return dt;
}

Real DriftFluidOp::computeCourantTimeStep( const FluidSpeciesPtrVect&  a_fluid_comp )
{
   
   const DisjointBoxLayout& grids = m_geometry.grids();
   
   // Compute physical ExB velocity drift on face centers
   LevelData<FluxBox> ExB_drift(grids, 3, IntVect::Unit);
   m_geometry.computeEXBDrift(m_E_field, ExB_drift);
      
   LevelData<FluxBox> vel_ExB(grids, SpaceDim, IntVect::Unit);
   for (DataIterator dit(vel_ExB.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         vel_ExB[dit][dir].copy(ExB_drift[dit][dir],0,0);
         if (SpaceDim == 2) {
            vel_ExB[dit][dir].copy(ExB_drift[dit][dir],2,1);
         }
         if (SpaceDim == 3) {
            vel_ExB[dit][dir].copy(ExB_drift[dit][dir],1,1);
            vel_ExB[dit][dir].copy(ExB_drift[dit][dir],2,2);
         }
         vel_ExB[dit][dir] *= m_larmor;
      }
   }
   m_geometry.applyAxisymmetricCorrection(vel_ExB);
   
   // Transform to mapped ExB velocity
   m_geometry.multNTransposePointwise( vel_ExB );
   
   Real dt = m_fluid_op_utils->computeAdvectionDt(vel_ExB,
                                                  m_advection_scheme,
                                                  m_time_step_diagnostics);
   
   return dt;
}


void DriftFluidOp::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "time_implicit", m_is_time_implicit);
   
   a_pp.query( "linearized_model", m_is_linearized);
      
   GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
   std::string grid_function_name;
   
   if (m_is_linearized) {
      if (a_pp.contains("background_density")) {
         a_pp.get("background_density", grid_function_name );
         m_n0_func = grid_library->find( grid_function_name );
      }
      else {
         MayDay::Error("DriftFluidOp::background density must be specified for the linearized model");
      }
   }
   
   a_pp.query("advection_scheme", m_advection_scheme );
   a_pp.query("max_dt", m_max_dt);
   a_pp.query("time_step_diagnostics", m_time_step_diagnostics);
   
   a_pp.query( "include_hyperviscosity", m_include_hyperviscosity);
   a_pp.query( "use_hyperviscosity_bcs", m_use_hyperviscosity_bcs );
      
   a_pp.query( "include_boundary_relaxation", m_include_boundary_relaxation);
   
   a_pp.queryarr( "harmonic_filtering", m_harmonic_filtering, 0, SpaceDim );
}


void DriftFluidOp::printParameters()
{
   if (procID()==0) {
      std::cout << "DriftFluidOp parameters:" << std::endl;
      std::cout << " time_implicit  =  " << m_is_time_implicit << std::endl;
   }
}


#include "NamespaceFooter.H"
