#include "OneFieldOp.H"
#include "EllipticOpBCFactory.H"
#include "FluidVarBCFactory.H"
#include "FourthOrderUtil.H"


#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "GKOps.H"
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H" 

const char* OneFieldOp::pp_name = {"one_field_op"};


OneFieldOp::OneFieldOp( const string&   a_pp_str,
                        const string&   a_species_name,
                        const MagGeom&  a_geometry,
                        const int       a_verbosity )
   : m_verbosity(a_verbosity),
     m_geometry(a_geometry),
     m_is_time_implicit(true),
     m_opt_string(a_pp_str),
     m_my_pc_idx(-1),
     m_D_rad(NULL),
     m_D_par(NULL),
     m_D_perp(NULL)

{
   ParmParse pp(a_pp_str.c_str());
   parseParameters( pp );
   if (m_verbosity>0) {
      printParameters();
   }

   const std::string name("neutrals_diffusion");
   const std::string prefix( "BC." + name );
   ParmParse ppsp( prefix.c_str() );

   m_diffusion_op = new Diffusion(ppsp, a_geometry);

   EllipticOpBCFactory elliptic_op_bc_factory;
   m_bc = elliptic_op_bc_factory.create( name,
                                         ppsp,
                                         a_geometry.getCoordSys()->type(),
                                         false );

   string variable_name = "density";

   // Input the initial conditions

   parseInitialConditions(a_species_name, variable_name);

   // Input the boundary conditions

   FluidVarBCFactory fluid_var_bc_factory;
   m_fluid_variable_bc = fluid_var_bc_factory.create(a_species_name, variable_name, a_geometry.getCoordSys()->type(), false);
}


OneFieldOp::~OneFieldOp()
{
   delete m_fluid_variable_bc;
   delete m_bc;
   delete m_diffusion_op;
}


void OneFieldOp::accumulateRHS(  FluidSpeciesPtrVect&               a_rhs,
                                 const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                 const FluidSpeciesPtrVect&         a_fluid_species,
                                 const PS::ScalarPtrVect&           a_scalars,
                                 const EField&                      a_E_field,
                                 const int                          a_fluidVecComp,
                                 const Real                         a_time)
{
   // Get fluid rhs
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   LevelData<FArrayBox>& rhs_data( rhs_fluid.cell_var(0) );

   // Get fluid soln
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
   const LevelData<FArrayBox>& soln_data( soln_fluid.cell_var(0) );

   const DisjointBoxLayout& grids( rhs_data.getBoxes() );
   
   LevelData<FluxBox> D_tensor(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   LevelData<FluxBox> D_tensor_mapped(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   computeDiffusionCoefficients(D_tensor, D_tensor_mapped, a_time);
   
   LevelData<FArrayBox> flux_div;
   flux_div.define(rhs_data);

   m_diffusion_op->setOperatorCoefficients(D_tensor, D_tensor_mapped, *m_bc);
   m_diffusion_op->computeFluxDivergence( soln_data, flux_div, false, false);

   // The preceding call computes a negated physical flux divergence cell average.
   // Since this function is being used to define the right-hand side of an evolution
   // equation for a mapped unknown (i.e., containing a J factor), we therefore need
   // to multiply by J.  We also flip the sign.
   m_geometry.multJonValid(flux_div);
   
   for (DataIterator dit( flux_div.dataIterator() ); dit.ok(); ++dit) {
      rhs_data[dit] -= flux_div[dit];
   }

   //add source term
   if (m_source != NULL) {
     LevelData<FArrayBox> source_dst(grids, 1, IntVect::Zero);
     m_source->assign( source_dst, m_geometry, a_time);

     m_geometry.multJonValid(source_dst);

     for (DataIterator dit(grids); dit.ok(); ++dit) {
       rhs_data[dit].plus(source_dst[dit]);
     }
   }
}


void OneFieldOp::accumulateExplicitRHS(  FluidSpeciesPtrVect&               a_rhs,
                                         const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                         const FluidSpeciesPtrVect&         a_fluid_species,
                                         const PS::ScalarPtrVect&           a_scalars,
                                         const EField&                      a_E_field,
                                         const int                          a_fluidVecComp,
                                         const Real                         a_time)
{
  if (!m_is_time_implicit) {
    accumulateRHS(a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, a_fluidVecComp, a_time);
  }
}


void OneFieldOp::accumulateImplicitRHS(  FluidSpeciesPtrVect&               a_rhs,
                                         const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                         const FluidSpeciesPtrVect&         a_fluid_species,
                                         const PS::ScalarPtrVect&           a_scalars,
                                         const EField&                      a_E_field,
                                         const int                          a_fluidVecComp,
                                         const Real                         a_time)
{
  if (m_is_time_implicit) {
    accumulateRHS(a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, a_fluidVecComp, a_time);
  }
}


void OneFieldOp::defineBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
                                 std::vector<PS::DOFList>&                                 a_dof_list,
                                 const PS::GKVector&                                       a_soln_vec,
                                 PS::GKOps&                                                a_gkops,
                                 const std::string&                                        a_out_string,
                                 const std::string&                                        a_opt_string,
                                 bool                                                      a_im,
                                 const FluidSpecies&                                       a_fluid_species,
                                 const PS::GlobalDOFFluidSpecies&                          a_global_dofs,
                                 int                                                       a_species_idx )
{
  if (a_im && m_is_time_implicit) {
    CH_assert(a_pc.size() == a_dof_list.size());
  
    if (!procID()) {
      std::cout << "  Fluid Species " << a_species_idx
                << " : "
                << " creating " << _FLUID_OP_PC_ << " preconditioner"
                << " (index = " << a_pc.size() << ").\n";
    }
  
    PS::Preconditioner<PS::GKVector, PS::GKOps> *pc;
    pc = new PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>;
    dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>
      (pc)->define(a_soln_vec, a_gkops, *this, m_opt_string, m_opt_string, a_im);
    dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>
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


void OneFieldOp::updateBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
                                 const PS::KineticSpeciesPtrVect&                          a_kin_species_phys,
                                 const FluidSpeciesPtrVect&                                a_fluid_species,
                                 const Real                                                a_time,
                                 const Real                                                a_shift,
                                 const bool                                                a_im,
                                 const int                                                 a_species_idx )
{
  if (a_im && m_is_time_implicit) {
    CH_assert(m_my_pc_idx >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx);
  
    if (!procID()) {
      std::cout << "    ==> Updating " << _FLUID_OP_PC_ << " preconditioner " 
                << " for OneFieldOp of fluid species " << a_species_idx << ".\n";
    }
  
    PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps> *pc 
      = dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>(a_pc[m_my_pc_idx]);
    CH_assert(pc != NULL);
    pc->update(a_kin_species_phys, a_fluid_species, a_time, a_shift, a_im, a_species_idx);
  }
  return;
}


void OneFieldOp::updatePCImEx( const PS::KineticSpeciesPtrVect& a_kinetic_species,
                               const double                     a_time,
                               const double                     a_shift )
{
   CH_TIME("OneFieldOp::updatePCImEx");

   m_diffusion_op->updateImExPreconditioner( a_shift, *m_bc );
}


void OneFieldOp::solvePCImEx( FluidSpeciesPtrVect&              a_fluid_species_solution,
                              const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                              const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                              const int                         a_component )
{
   CH_TIME("OneFieldOp::solvePCImEx");
   if (!m_is_time_implicit) {
     MayDay::Error("Error in OneFieldOp::solvePCImEx(): m_is_time_implicit is false. This function shouldn't have been called at all!");
   }
   const FluidSpecies& rhs_species = static_cast<const FluidSpecies&>(*a_fluid_species_rhs[a_component]);
   LevelData<FArrayBox> r_phys;
   r_phys.define(rhs_species.cell_var(0));

   // Convert the right-hand side to physical space
   m_geometry.divideJonValid(r_phys);
   
   FluidSpecies& sol_species = static_cast<FluidSpecies&>(*a_fluid_species_solution[a_component]);
   LevelData<FArrayBox>& z = sol_species.cell_var(0);

   for (DataIterator dit(z.dataIterator()); dit.ok(); ++dit) {
      z[dit].setVal(0.);
   }

   double tol = 0.;
   double max_iter = 1;
   double precond_tol = 0.;
   double precond_max_iter = 1;
   m_diffusion_op->setImExPreconditionerConvergenceParams(tol, max_iter, precond_tol, precond_max_iter);

   m_diffusion_op->solveImExPreconditioner(r_phys, z);

   // Convert the solution to mapped space
   m_geometry.multJonValid(z);
}


void OneFieldOp::computeDiffusionCoefficients(LevelData<FluxBox>& a_D_tensor,
                                              LevelData<FluxBox>& a_D_tensor_mapped,
                                              const Real          a_time) const
{
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   LevelData<FluxBox> D_tmp(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   LevelData<FluxBox> D_tmp_mapped(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   
   //Parallel and perpendicular tensor coefficeints
   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& par_coeff = m_geometry.getEllipticOpParCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();
   
   //Temporary to hold scalar diffusion coefficients
   LevelData<FArrayBox> D_cell(grids, 1, 4*IntVect::Unit);
   LevelData<FluxBox> D_face(grids, 1, 2*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     a_D_tensor[dit].setVal(0.0);
     a_D_tensor_mapped[dit].setVal(0.0);
   }   

   //Add perpendicular diffusion term, if any
   if (m_D_perp != NULL) {
      m_D_perp->assign( D_cell, m_geometry, a_time);
      m_geometry.fillInternalGhosts( D_cell );
      fourthOrderCellToFaceCenters(D_face, D_cell);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_tmp[dit].copy(perp_coeff[dit]);
         D_tmp_mapped[dit].copy(perp_coeff_mapped[dit]);
         for (int dir = 0; dir < SpaceDim; dir++) {
            for (int n = 0; n < SpaceDim*SpaceDim; n++) {
               D_tmp[dit][dir].mult(D_face[dit][dir],0,n,1);
               D_tmp_mapped[dit][dir].mult(D_face[dit][dir],0,n,1);
            }
         }
         a_D_tensor[dit]+=D_tmp[dit];
         a_D_tensor_mapped[dit]+=D_tmp_mapped[dit];
      }
   }
   
   //Add parallel diffusion term, if any
   if (m_D_par != NULL) {
      m_D_par->assign( D_cell, m_geometry, a_time);
      m_geometry.fillInternalGhosts( D_cell );
      fourthOrderCellToFaceCenters(D_face, D_cell);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_tmp[dit].copy(par_coeff[dit]);
         D_tmp_mapped[dit].copy(par_coeff_mapped[dit]);
         for (int dir = 0; dir < SpaceDim; dir++) {
            for (int n = 0; n < SpaceDim*SpaceDim; n++) {
               D_tmp[dit][dir].mult(D_face[dit][dir],0,n,1);
               D_tmp_mapped[dit][dir].mult(D_face[dit][dir],0,n,1);
            }
         }
         a_D_tensor[dit]+=D_tmp[dit];
         a_D_tensor_mapped[dit]+=D_tmp_mapped[dit];
      }
   }
   
   //Add radial diffusion term, if any
   if (m_D_rad != NULL) {
      m_D_rad->assign( D_cell, m_geometry, a_time);
      m_geometry.fillInternalGhosts( D_cell );
      fourthOrderCellToFaceCenters(D_face, D_cell);
      
      m_geometry.getEllipticOpRadCoeff(D_tmp);
      m_geometry.getEllipticOpRadCoeffMapped(D_tmp_mapped);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir = 0; dir < SpaceDim; dir++) {
            for (int n = 0; n < SpaceDim*SpaceDim; n++) {
               D_tmp[dit][dir].mult(D_face[dit][dir],0,n,1);
               D_tmp_mapped[dit][dir].mult(D_face[dit][dir],0,n,1);
            }
         }
         a_D_tensor[dit]+=D_tmp[dit];
         a_D_tensor_mapped[dit]+=D_tmp_mapped[dit];
      }
   }
}


void OneFieldOp::fillGhostCells( FluidSpecies&  a_species_phys,
                                 const double   a_time )
{
}


void OneFieldOp::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "time_implicit", m_is_time_implicit);

   if (a_pp.contains("D_rad")) {
     GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
     std::string grid_function_name;
     a_pp.get("D_rad", grid_function_name );
     m_D_rad = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("D_perp")) {
      GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
      std::string grid_function_name;
      a_pp.get("D_perp", grid_function_name );
      m_D_perp = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("D_par")) {
      GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
      std::string grid_function_name;
      a_pp.get("D_par", grid_function_name );
      m_D_par = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("source")) {
     GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
     std::string grid_function_name;
     a_pp.get("source", grid_function_name );
     m_source = grid_library->find( grid_function_name );
   }

}


void OneFieldOp::printParameters()
{
   if (procID()==0) {
      std::cout << "OneFieldOp parameters:" << std::endl;
      std::cout << "  time_implicit  =  " << m_is_time_implicit << std::endl;
   }
}

Real OneFieldOp::computeDt( const FluidSpeciesPtrVect&  a_fluid_species )
{
   return DBL_MAX;
}


#include "NamespaceFooter.H"
