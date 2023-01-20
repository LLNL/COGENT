#include "DriftFluidOp.H"
#include "EllipticOpBCFactory.H"
#include "FluidVarBCFactory.H"
#include "FourthOrderUtil.H"
#include "MagBlockCoordSys.H"
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
     m_include_hyperviscosity(false),
     m_use_hyperviscosity_bcs(false),
     m_species_name(a_species_name),
     m_opt_string(a_pp_str),
     m_my_pc_idx(-1),
     m_n0_func(NULL),
     m_diffusion_op(NULL),
     m_diffusion_op_bcs(NULL)

{
   ParmParse pp(a_pp_str.c_str());
   parseParameters( pp );
   if (m_verbosity>0) {
      printParameters();
   }

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
      m_diffusion_op = new Diffusion(ppOp, a_geometry);

      const std::string name( m_species_name + ".hyperviscosity");
      const std::string prefixBC( "BC." + name);
      ParmParse ppOpBC( prefixBC.c_str() );

      EllipticOpBCFactory elliptic_op_bc_factory;
      m_diffusion_op_bcs = elliptic_op_bc_factory.create(name,
                                                         ppOpBC,
                                                         *(a_geometry.getCoordSys()),
                                                         false );
   }
}


DriftFluidOp::~DriftFluidOp()
{
   if (m_fluid_variable_bc) delete m_fluid_variable_bc;
   if (m_diffusion_op_bcs) delete m_diffusion_op_bcs;
   if (m_diffusion_op) delete m_diffusion_op;
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
      Physics model: d(n)/dt = -1/c * ExB/B^2 * nabla(n0)
      COGENT equation: d(J*n)/dt = -J * Larmor * ExB/B^2 * nabla(n0)
   */
   
   // Get fluid rhs
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   LevelData<FArrayBox>& rhs_data( rhs_fluid.cell_var(0) );

   // Get fluid soln
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
   const LevelData<FArrayBox>& soln_data( soln_fluid.cell_var(0) );


   // Get geoemtry parameters
   int order = (m_geometry.secondOrder()) ? 2 : 4;
   const DisjointBoxLayout& grids( soln_data.getBoxes() );


   // Get backgruond density gradient
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
      m_geometry.unmapGradient(n0_gradient_mapped, m_n0_gradient );
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

   // update RHS
   for (DataIterator dit(grids); dit.ok(); ++dit) {
     rhs_data[dit] -= ExB_grad_n0[dit];
   }

   // Add parallel hyperviscosity
   if (m_include_hyperviscosity) {

      // Get Dshape function
      if (m_first_call) {
         LevelData<FArrayBox> D_par_cc(grids, 1, 2*IntVect::Unit);
         m_D_par_func->assign( D_par_cc, m_geometry, a_time);
         
         m_geometry.fillInternalGhosts( D_par_cc );
         
         m_D_par_fc.define(grids, 1, IntVect::Zero);
         fourthOrderCellToFaceCenters(m_D_par_fc, D_par_cc);
      }
   
      // Compute diffusion coefficient
      LevelData<FluxBox> D_tensor(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
      LevelData<FluxBox> D_tensor_mapped(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);

      // Get parallel elliptic coefficients
      const LevelData<FluxBox>& par_coeff = m_geometry.getEllipticOpParCoeff();
      const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_tensor[dit].copy(par_coeff[dit]);
         D_tensor_mapped[dit].copy(par_coeff_mapped[dit]);
         for (int dir = 0; dir < SpaceDim; dir++) {
            for (int n = 0; n < SpaceDim*SpaceDim; n++) {
               D_tensor[dit][dir].mult(m_D_par_fc[dit][dir],0,n,1);
               D_tensor_mapped[dit][dir].mult(m_D_par_fc[dit][dir],0,n,1);
            }
         }
      }
   
      // Set diffisionOf coefficients
      m_diffusion_op->setOperatorCoefficients(D_tensor, D_tensor_mapped, *m_diffusion_op_bcs);

      // Compute diffusion flux divergence
      LevelData<FArrayBox> flux_div(grids, 1, IntVect::Zero);
      m_diffusion_op->computeFluxDivergence( soln_data, flux_div, false, !m_use_hyperviscosity_bcs);

      LevelData<FArrayBox> hypervisc(grids, 1, IntVect::Zero);
      m_diffusion_op->computeFluxDivergence( flux_div, hypervisc, false, !m_use_hyperviscosity_bcs);
      m_geometry.multJonValid(hypervisc);
   
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs_data[dit] -= flux_div[dit];
      }
   }
   
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

void DriftFluidOp::fillGhostCells( FluidSpecies&  a_species_phys,
                              const double   a_time )
{
#if 0
   CH_TIME("DriftFluidOp::fillGhostCells()");
   
   // Do we need it here, or we alredy pre-filled the ghosts
   // Definitely don't need for the ExB*gradn0 case
   
   // Fill ghost cells except for those on physical boundaries
   LevelData<FArrayBox>& fld( a_species_phys.cell_var(0) );
   m_geometry.fillInternalGhosts( fld );
   // Fill ghost cells on physical boundaries
   m_fluid_variable_bc->apply( a_species_phys, a_time );
#endif
}


void DriftFluidOp::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "time_implicit", m_is_time_implicit);
      
   GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
   std::string grid_function_name;
   a_pp.get("background_density", grid_function_name );
   m_n0_func = grid_library->find( grid_function_name );
   
   a_pp.query( "include_parallel_hyperviscosity", m_include_hyperviscosity);
   a_pp.query( "use_hyperviscosity_bcs", m_use_hyperviscosity_bcs );
   
   if (m_include_hyperviscosity) {
      if (a_pp.contains("D_par_hyperviscosity")) {
         a_pp.get("D_par_hyperviscosity", grid_function_name );
         m_D_par_func = grid_library->find( grid_function_name );
      }
      else {
         MayDay::Error("parallel hyperviscosity coefficient must be specified");
      }
   }
   
}


void DriftFluidOp::printParameters()
{
   if (procID()==0) {
      std::cout << "DriftFluidOp parameters:" << std::endl;
      std::cout << " time_implicit  =  " << m_is_time_implicit << std::endl;
   }
}


#include "NamespaceFooter.H"
