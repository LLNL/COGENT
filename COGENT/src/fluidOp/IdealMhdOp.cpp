#include "IdealMhdOp.H"
#include "EllipticOpBCFactory.H"
#include "FluidVarBCFactory.H"
#include "FourthOrderUtil.H"
#include "MagBlockCoordSys.H"
#include "inspect.H"
#include "CONSTANTS.H"
//#include "FluidOpF_F.H"

#include "NodeFArrayBox.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
//#include "FluidOpF_F.H"
#include "GKOps.H"
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

//#include "SPACE.H"
#include "FluidOpF_F.H"

#include "NamespaceHeader.H" 

const char* IdealMhdOp::pp_name = {"idealmhd_op"};


IdealMhdOp::IdealMhdOp( const string&   a_pp_str,
                        const string&   a_species_name,
                        const MagGeom&  a_geometry,
                        const int       a_verbosity )
   : m_verbosity(a_verbosity),
     m_geometry(a_geometry),
     m_is_time_implicit(true),
     m_advScheme("c2"),
     m_species_name(a_species_name),
     m_courant_time_step(DBL_MAX),
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
   
   //  get initial conditions for cell variables and define boundary conditions instance
   //
   const std::string icsn( "IC." + a_species_name);
   std:: string variable_name;
   ParmParse ppic( icsn.c_str() );
   FluidVarBCFactory fluid_var_bc_factory;
  
   // set IC and create BC for mass density
   //
   variable_name = "density";
   CH_assert( ppic.contains("density.function") );
   parseInitialConditions(a_species_name, variable_name);
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );
 
   // set IC and create BC for momentum density
   //
   variable_name = "momentumDensity";
   CH_assert( ppic.contains("momentumDensity_x.function") );
   CH_assert( ppic.contains("momentumDensity_y.function") );
   parseInitialConditions(a_species_name, "momentumDensity_x");
   parseInitialConditions(a_species_name, "momentumDensity_y");
   if (ppic.contains("momentumDensity_z.function")) {
      parseInitialConditions(a_species_name, "momentumDensity_z");
   }
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );
 
   // set IC and create BC for energy density
   //
   variable_name = "energyDensity";
   CH_assert( ppic.contains("energyDensity.function") );
   parseInitialConditions(a_species_name, variable_name);
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );
   
   
   if (ppic.contains("momentumDensity_virtual.function")) {
      variable_name = "momentumDensity_virtual";
      parseInitialConditions(a_species_name, variable_name);
      m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                            a_geometry.getCoordSys()->type(), false) );
   }
   
   if (ppic.contains("magneticField_virtual.function")) {
      variable_name = "magneticField_virtual";
      parseInitialConditions(a_species_name, variable_name);
      m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                            a_geometry.getCoordSys()->type(), false) );
   }
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   IntVect ghostVect = 2*IntVect::Unit;
   //const IntVect& ghostVect = m_geometry.m_num_ghosts*IntVect::Unit;
   if(m_advScheme == "weno5" || m_advScheme == "bweno" || m_advScheme == "uw5") {
      //CH_assert(m_geometry.m_num_ghosts>=3);
      ghostVect = 4*IntVect::Unit;
   }
   defineLevelDatas( grids, ghostVect );
   //cout << "JRA: IdealMHDOp: ghostVect = " << ghostVect << endl;

}


IdealMhdOp::~IdealMhdOp()
{
   //delete m_fluid_variable_bc;
   for (std::vector<FluidVarBC*>::iterator it = m_fluid_bc.begin(); it != m_fluid_bc.end(); ++it)
   {
      delete (*it);
   }
   m_fluid_bc.clear();
   delete m_bc;
   delete m_diffusion_op;
}


void IdealMhdOp::accumulateRHS(  FluidSpeciesPtrVect&               a_rhs,
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
   LevelData<FArrayBox>& rhs_data1( rhs_fluid.cell_var(1) );

   // Get fluid soln
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
   const LevelData<FArrayBox>& soln_data( soln_fluid.cell_var(0) );
   const LevelData<FArrayBox>& soln_data1( soln_fluid.cell_var(1) );

   const DisjointBoxLayout& grids( rhs_data.getBoxes() );
   
   LevelData<FluxBox> D_tensor(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   LevelData<FluxBox> D_tensor_mapped(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   computeDiffusionCoefficients(D_tensor, D_tensor_mapped, a_time);
   
   LevelData<FArrayBox> flux_div, flux_div1;
   flux_div.define(rhs_data);
   flux_div1.define(rhs_data1);

   m_diffusion_op->setOperatorCoefficients(D_tensor, D_tensor_mapped, *m_bc);
   m_diffusion_op->computeFluxDivergence( soln_data, flux_div, false, false);
   m_diffusion_op->computeFluxDivergence( soln_data1, flux_div1, false, false);

   // The preceding call computes a negated physical flux divergence cell average.
   // Since this function is being used to define the right-hand side of an evolution
   // equation for a mapped unknown (i.e., containing a J factor), we therefore need
   // to multiply by J.  We also flip the sign.
   m_geometry.multJonValid(flux_div);
   m_geometry.multJonValid(flux_div1);
   
   for (DataIterator dit( flux_div.dataIterator() ); dit.ok(); ++dit) {
      rhs_data[dit]  -= flux_div[dit];
      rhs_data1[dit] -= flux_div1[dit];
   }

}


void IdealMhdOp::accumulateExplicitRHS( FluidSpeciesPtrVect&         a_rhs,
                                  const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                  const FluidSpeciesPtrVect&         a_fluid_species,
                                  const PS::ScalarPtrVect&           a_scalars,
                                  const EField&                      a_E_field,
                                  const int                          a_fluidVecComp,
                                  const Real                         a_time)
{
   CH_TIME("IdealMhdOp::accumulateExplicitRHS()");

   // Get RHS fluid species
   //
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   LevelData<FArrayBox>& rhs_data( rhs_fluid.cell_var(0) );
   const DisjointBoxLayout& grids( rhs_data.getBoxes() );
   //const IntVect& ghostVect( rhs_data.ghostVect() ); // caution, this is zero
   //cout << "JRA: rhs.ghostVect" << ghostVect << endl;   
   //cout << grids << endl;   

   // Get solution fluid species
   //
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
   const LevelData<FArrayBox>& soln_data( soln_fluid.cell_var(0) );
   //const IntVect& ghostVect( soln_data.ghostVect() );
   //cout << "JRA: soln_data.ghostVect() = " << ghostVect << endl;   
   
   //CH_assert(soln_fluid.m_evolve_massDensity==1);
   CH_assert(soln_fluid.m_evolve_momentumDensity==1);
   CH_assert(soln_fluid.m_evolve_energyDensity==1);
   
   setCellCenterValues( soln_fluid );

   setCellCenterFluxes( soln_fluid );
   
   // compute flux-freezing speed at cell-center for each direction
   // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/N); i = q0, q1, q2
   //
   const double gamma = soln_fluid.m_gamma;
   setCspeed( m_rhoDen_cc, gamma );
   SpaceUtils::upWindToFaces(m_CspeedR_norm, m_Cspeed_cc, m_CspeedR_norm, "c2"); // 3rd arg not used for "c2"
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_CspeedL_norm[dit].copy(m_CspeedR_norm[dit],m_CspeedL_norm[dit].box());
      m_CspeedL_norm[dit].negate();
   }
   
   setMagneticFieldTerms( soln_fluid );

   setFaceCenteredFluxes( soln_fluid );

   enforceFluxBCs( rhs_fluid, soln_fluid, a_time );

   updateRHSs( rhs_fluid, soln_fluid );

}


void IdealMhdOp::accumulateImplicitRHS( FluidSpeciesPtrVect&        a_rhs,
                                  const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                  const FluidSpeciesPtrVect&        a_fluid_species,
                                  const PS::ScalarPtrVect&          a_scalars,
                                  const EField&                     a_E_field,
                                  const int                         a_fluidVecComp,
                                  const Real                        a_time )
{
  // JRA, turning diffusion off for now
  //
  //if (m_is_time_implicit) {
  //  accumulateRHS(a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, a_fluidVecComp, a_time);
  //}
}


void IdealMhdOp::defineBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
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


void IdealMhdOp::updateBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
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
                << " for IdealMhIdealMhdf fluid species " << a_species_idx << ".\n";
    }
  
    PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps> *pc 
      = dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>(a_pc[m_my_pc_idx]);
    CH_assert(pc != NULL);
    pc->update(a_kin_species_phys, a_fluid_species, a_time, a_shift, a_im, a_species_idx);
  }
  return;
}


void IdealMhdOp::updatePCImEx(const FluidSpeciesPtrVect&       a_fluid_species, 
                              const PS::KineticSpeciesPtrVect& a_kinetic_species,
                              const double                     a_time,
                              const double                     a_shift,
                              const int                        a_component)
{
   CH_TIME("IdealMhdOp::updatePCImEx()");
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   LevelData<FArrayBox> shift(grids, 1, IntVect::Zero);
   for (DataIterator dit(shift.dataIterator()); dit.ok(); ++dit) {
      shift[dit].setVal(a_shift);
   }

   m_diffusion_op->updateImExPreconditioner( shift, *m_bc );
}


void IdealMhdOp::solvePCImEx( FluidSpeciesPtrVect&              a_fluid_species_solution,
                           const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                           const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                           const int                         a_component )
{
   CH_TIME("IdealMhdOp::solvePCImEx()");
   if (!m_is_time_implicit) {
     MayDay::Error("Error in IdealMhdOp::solvePCImEx(): m_is_time_implicit is false. This function shouldn't have been called at all!");
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


void IdealMhdOp::computeDiffusionCoefficients(LevelData<FluxBox>& a_D_tensor,
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

void IdealMhdOp::inspectFArrayBox(const LevelData<FArrayBox>& a_F0,
                               const int                   a_comp) const
{
   CH_assert(a_comp<a_F0.nComp());
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& F0_on_patch = a_F0[dit];   
      const Box& thisbox = F0_on_patch.box();

      FORT_INSPECT_FARRAYBOX( CHF_BOX(thisbox), 
                              CHF_CONST_FRA1(F0_on_patch,a_comp) );
   }
}

void IdealMhdOp::inspectFluxBox(const LevelData<FluxBox>& a_Flux,
                             const int                 a_dir) const
{
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FluxBox& Flux_on_patch = a_Flux[dit]; 
      const Box& thisbox = Flux_on_patch.box();

      //Box face_box( Flux_on_patch.box());
      //face_box.surroundingNodes( a_dir );
      //face_box.grow( a_dir, 1 );

      const FArrayBox& Flux_on_dir = Flux_on_patch[a_dir]; 

      FORT_INSPECT_FLUXBOX( CHF_BOX(thisbox), 
                            CHF_CONST_FRA(Flux_on_dir),
                            CHF_CONST_INT(a_dir) );
   }
}

void IdealMhdOp::inspectEdgeDataBox(const LevelData<EdgeDataBox>& a_Edge,
                                    const int                     a_dir) const
{
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const EdgeDataBox& Edge_on_patch = a_Edge[dit]; 
      const Box& thisbox = Edge_on_patch.box();

      const FArrayBox& Edge_on_dir = Edge_on_patch[a_dir]; 

      FORT_INSPECT_FLUXBOX( CHF_BOX(thisbox), 
                            CHF_CONST_FRA(Edge_on_dir),
                            CHF_CONST_INT(a_dir) );
   }
}

void IdealMhdOp::convertPhysFluxToMappedFlux(LevelData<FArrayBox>& a_mapped_flux,
                                       const LevelData<FArrayBox>& a_physical_flux,
                                       const LevelData<FArrayBox>& a_N) const
{
   const DisjointBoxLayout& grids( m_geometry.grids() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {

      //const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);

      const FArrayBox& phys_on_patch = a_physical_flux[dit];   
      const FArrayBox& N_on_patch = a_N[dit]; 
      FArrayBox& mapped_on_patch = a_mapped_flux[dit];
      
      convertPhysFluxToMappedFlux(mapped_on_patch,phys_on_patch,N_on_patch);
      /*
      const Box& thisbox = mapped_on_patch.box();
      FORT_CONVERT_TO_MAPPED( CHF_BOX(thisbox), 
                              CHF_CONST_FRA(phys_on_patch),
                              CHF_CONST_FRA(N_on_patch),
                              CHF_CONST_INT(SpaceDim),
                              CHF_FRA(mapped_on_patch) );
      */
   }
}

void IdealMhdOp::convertPhysFluxToMappedFlux( FArrayBox&  a_mapped_on_patch,
                                        const FArrayBox&  a_phys_on_patch,
                                        const FArrayBox&  a_N_on_patch ) const
{
   CH_TIME("IdealMhdOp::convertPhysFluxToMappedFlux()");
   const Box& thisbox = a_mapped_on_patch.box();

   FORT_CONVERT_TO_MAPPED( CHF_BOX(thisbox), 
                           CHF_CONST_FRA(a_phys_on_patch),
                           CHF_CONST_FRA(a_N_on_patch),
                           CHF_CONST_INT(SpaceDim),
                           CHF_FRA(a_mapped_on_patch) );

}

void IdealMhdOp::computeNTFfaceArea( LevelData<FluxBox>&    a_Flux_norm,
                               const LevelData<FArrayBox>&  a_Flux_cc,
                               const LevelData<FArrayBox>&  a_Jf_cc,
                               const LevelData<FArrayBox>&  a_Cspeed_cc,
                               const LevelData<FluxBox>&    a_CspeedR_norm,
                               const LevelData<FluxBox>&    a_CspeedL_norm ) const
{
   CH_TIME("IdealMhdOp::computeNTFfaceArea()");
   CH_assert(a_Flux_norm.nComp() == 1);
   CH_assert(a_Flux_cc.nComp() == SpaceDim);
   
   const DisjointBoxLayout& grids( a_Jf_cc.getBoxes() );
   //const DisjointBoxLayout& grids2( a_Flux_norm.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
   
      // convert cell-center flux from physical to mapped 
      //
      const Box& cell_box( a_Jf_cc[dit].box() );
      FArrayBox Flux_cc_mapped(cell_box, a_Flux_cc.nComp());
      convertPhysFluxToMappedFlux(Flux_cc_mapped,a_Flux_cc[dit],m_Nmatrix[dit]);

      // compute Lax splitting on mapped flux and interpolate FluxL and FluxR to face norms
      //
      //   get left and right going flux at cell-center
      //   fluxR = 0.5*(flux + Cspeed*fun),
      //   fluxL = 0.5*(flux - Cspeed*fun),
      //   Cspeed = abs(max(eigenValue of Flux Jacobian))
      //
      /*
      FArrayBox FluxL_cc;
      FArrayBox FluxR_cc;
      SpaceUtils::computeLaxSplitting(FluxR_cc,FluxL_cc,Flux_cc_mapped,a_Cspeed_cc[dit],a_Jf_cc[dit]);
      */

      FArrayBox FluxL_cc(cell_box, SpaceDim);
      FArrayBox FluxR_cc(cell_box, SpaceDim);
      const FArrayBox& Cspeed_on_patch = a_Cspeed_cc[dit];
      const FArrayBox& Jf_on_patch = a_Jf_cc[dit];

      FORT_LAX_FLUX_SPLITTING( CHF_BOX(cell_box), 
                               CHF_CONST_FRA(FluxL_cc),
                               CHF_CONST_FRA(FluxR_cc),
                               CHF_CONST_FRA(Flux_cc_mapped),
                               CHF_CONST_FRA(Cspeed_on_patch),
                               CHF_CONST_FRA1(Jf_on_patch,0) );

      const Box& grid_box( grids[dit] );             // this box has no ghost
      //const Box& grid_box2( grids2[dit] );           // this box has no ghost
      //cout << "grid_box()  = " << grid_box << endl;
      //cout << "grid_box2() = " << grid_box2 << endl;
      const Box& flux_box( a_Flux_norm[dit].box() ); // this box has ghost
      FluxBox thisFluxL_norm( flux_box, 1 );  // need box with ghosts here
      FluxBox thisFluxR_norm( flux_box, 1 );
      const FluxBox& thisCspeedL_norm = a_CspeedL_norm[dit];
      const FluxBox& thisCspeedR_norm = a_CspeedR_norm[dit];
      if(m_advScheme=="weno5") {
         const FArrayBox& thisN = m_Nmatrix[dit]; 
         FArrayBox  thisSmooth(cell_box, SpaceDim);
         
         FORT_MAG_NT_ROW_SQ( CHF_BOX(cell_box), 
                             CHF_CONST_FRA(thisN),
                             CHF_FRA(thisSmooth) );
         SpaceUtils::interpToFacesWENO(thisFluxL_norm, FluxL_cc, thisCspeedL_norm, thisSmooth, grid_box, m_advScheme);
         SpaceUtils::interpToFacesWENO(thisFluxR_norm, FluxR_cc, thisCspeedR_norm, thisSmooth, grid_box, m_advScheme);
      }
      else {
         SpaceUtils::upWindToFaces(thisFluxL_norm, FluxL_cc, thisCspeedL_norm, grid_box, m_advScheme);
         SpaceUtils::upWindToFaces(thisFluxR_norm, FluxR_cc, thisCspeedR_norm, grid_box, m_advScheme);
      }

      // compute total norm flux by adding left and right together 
      // and then multiply by face area
      //
      a_Flux_norm[dit].copy(thisFluxR_norm, flux_box );
      a_Flux_norm[dit] += thisFluxL_norm;
      const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const RealVect& faceArea = coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& this_face = a_Flux_norm[dit][dir];
         this_face.mult(faceArea[dir]);
         //cout << "faceArea[0] " << faceArea[0] << endl;
         //cout << "faceArea[1] " << faceArea[1] << endl;
      }
   }
   a_Flux_norm.exchange();
 
}

void IdealMhdOp::computeIdealEatEdges( LevelData<EdgeDataBox>&  a_Edge_covar,
                                 const LevelData<FArrayBox>&    a_V_phys_cc,
                                 const LevelData<FArrayBox>&    a_B_contr_cc,
                                 const LevelData<FArrayBox>&    a_Cspeed_cc ) const
{
   CH_TIME("IdealMhdOp::computeIdealEatEdges()");
   CH_assert(a_Edge_covar.nComp() == 1);
   CH_assert(a_V_phys_cc.nComp() == SpaceDim);
   CH_assert(a_Cspeed_cc.nComp() == SpaceDim);
  
  
    // interpolate components of Cspeed from cell-center to +/- face norms
   //
   const DisjointBoxLayout& grids( a_B_contr_cc.getBoxes() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
    
      const Box& grid_box( grids[dit] ); // this box has no ghosts
      const Box& cell_box(a_B_contr_cc[dit].box());
      FluxBox CspeedL_norm(cell_box, 1);
      FluxBox CspeedR_norm(cell_box, 1);
      SpaceUtils::upWindToFaces(CspeedR_norm, a_Cspeed_cc[dit], CspeedR_norm, grid_box, "c2"); // 3rd arg not used for "c2"
      CspeedL_norm.copy(CspeedR_norm, cell_box );
      CspeedL_norm.negate();
      
      // convert cell-center Velocity from physical to contravariant 
      // and multiply by contravariant B
      //
      FArrayBox VB_contr_cc(cell_box, a_V_phys_cc.nComp());
      FArrayBox JB_contr_cc(cell_box, a_B_contr_cc.nComp());
      convertPhysFluxToMappedFlux(VB_contr_cc,a_V_phys_cc[dit],m_Nmatrix[dit]);
      for (int n=0; n<SpaceDim; n++) {
         VB_contr_cc.mult(a_B_contr_cc[dit],cell_box,cell_box,0,n,1);
      }
      JB_contr_cc.copy(a_B_contr_cc[dit], cell_box );
      JB_contr_cc.mult(m_Jacobian[dit]);
  
      // compute Lax splitting and interpolate FluxL and FluxR to face norms
      //
      /*
      FArrayBox FluxL_cc;
      FArrayBox FluxR_cc;
      SpaceUtils::computeLaxSplitting(FluxR_cc,FluxL_cc,VB_contr_cc,a_Cspeed_cc[dit],JB_contr_cc);
      */
      
      FArrayBox FluxL_cc(cell_box, SpaceDim);
      FArrayBox FluxR_cc(cell_box, SpaceDim);
      const FArrayBox& Cspeed_on_patch = a_Cspeed_cc[dit];

      FORT_LAX_FLUX_SPLITTING( CHF_BOX(cell_box), 
                               CHF_CONST_FRA(FluxL_cc),
                               CHF_CONST_FRA(FluxR_cc),
                               CHF_CONST_FRA(VB_contr_cc),
                               CHF_CONST_FRA(Cspeed_on_patch),
                               CHF_CONST_FRA1(JB_contr_cc,0) );
     
      const Box& flux_box( a_Edge_covar[dit].box() ); // this box has ghosts
      FluxBox FluxL_norm( flux_box, 1 );
      FluxBox FluxR_norm( flux_box, 1 );
      FluxBox VB_norm( flux_box, 1 );
      if(m_advScheme=="weno5") {
         const FArrayBox& thisN = m_Nmatrix[dit]; 
         FArrayBox thisSmooth( cell_box, SpaceDim );

         FORT_MAG_NT_ROW_SQ( CHF_BOX(cell_box), 
                             CHF_CONST_FRA(thisN),
                             CHF_FRA(thisSmooth) );
         SpaceUtils::interpToFacesWENO(FluxL_norm, FluxL_cc, CspeedL_norm, thisSmooth, grid_box, m_advScheme);
         SpaceUtils::interpToFacesWENO(FluxR_norm, FluxR_cc, CspeedR_norm, thisSmooth, grid_box, m_advScheme);
      } 
      else {
         SpaceUtils::upWindToFaces(FluxL_norm, FluxL_cc, CspeedL_norm, grid_box, m_advScheme);
         SpaceUtils::upWindToFaces(FluxR_norm, FluxR_cc, CspeedR_norm, grid_box, m_advScheme);
      }

      // compute total flux on faces by adding left and right together
      //
      VB_norm.copy(FluxR_norm,VB_norm.box());
      VB_norm += FluxL_norm;
      a_Edge_covar[dit][0].copy(VB_norm[1]);
      a_Edge_covar[dit][1].copy(VB_norm[0]);
      a_Edge_covar[dit][1].negate();
   } 

}

void IdealMhdOp::fillGhostCells( FluidSpecies&  a_species_phys,
                              const double   a_time )
{
   CH_TIME("IdealMhdOp::fillGhostCells()");
   for (int n=0; n<a_species_phys.num_cell_vars(); ++n) {
      LevelData<FArrayBox>& fld( a_species_phys.cell_var(n) );
      // Fill ghost cells except for those on physical boundaries
      m_geometry.fillInternalGhosts( fld );
      m_fluid_bc.at(n)->apply( a_species_phys, a_time );
   }
}


void IdealMhdOp::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "time_implicit", m_is_time_implicit);

   GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
   std::string grid_function_name;
   
   if (a_pp.contains("D_rad")) {
     a_pp.get("D_rad", grid_function_name );
     m_D_rad = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("D_perp")) {
      a_pp.get("D_perp", grid_function_name );
      m_D_perp = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("D_par")) {
      a_pp.get("D_par", grid_function_name );
      m_D_par = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("advScheme")) {
      a_pp.get("advScheme", m_advScheme );
   }
   
}


void IdealMhdOp::printParameters()
{
   if (procID()==0) {
      std::cout << "IdealMhdOp parameters:" << std::endl;
      std::cout << " time_implicit  =  " << m_is_time_implicit << std::endl;
      std::cout << " advScheme  =  " << m_advScheme << std::endl;
   }
}

void IdealMhdOp::setCourantTimeStep( const LevelData<FArrayBox>&  a_Cspeed )
{
   CH_assert(a_Cspeed.nComp() == SpaceDim);
    
   Vector<double> local_Cmax(SpaceDim,-CH_BADVAL);
   RealVect dX_mapped;
   const DisjointBoxLayout& grids( a_Cspeed.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      dX_mapped = coord_sys.getMappedCellSize();
      
      Box box( grids[dit] );
      for (int dir=0; dir<SpaceDim; ++dir) {
         double box_max( a_Cspeed[dit].max( box,dir ) );
         local_Cmax[dir] = Max( local_Cmax[dir], box_max );
      }
   }
   //cout << "local_Cmax[0] = " << local_Cmax[0] << endl;
   //cout << "local_Cmax[1] = " << local_Cmax[1] << endl;
   //cout << "dX_mapped[0] = " << dX_mapped[0] << endl;
   //cout << "dX_mapped[1] = " << dX_mapped[1] << endl;
    
   
   Vector<double> Cmax( local_Cmax );
   double this_dt_max=DBL_MAX;
   for (int dir=0; dir<SpaceDim; ++dir) {
#ifdef CH_MPI
      MPI_Allreduce( &local_Cmax[dir], &Cmax[dir], 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif
      this_dt_max = Min( dX_mapped[dir]/Cmax[dir], this_dt_max);
   }
   //cout << "Cmax[0] = " << Cmax[0] << endl;
   //cout << "Cmax[1] = " << Cmax[1] << endl;
   //cout << "this_dt_max = " << this_dt_max << endl;

   m_courant_time_step = this_dt_max;
   //cout << "m_courant_time_step = " << m_courant_time_step << endl;

}

void IdealMhdOp::setCspeed( const LevelData<FArrayBox>&  a_density,
                            const double                 a_gamma )
{
   CH_TIME("IdealMhdOp::setCspeed()");

   CH_assert(m_Cspeed_cc.nComp() == SpaceDim);
   
   const DisjointBoxLayout& grids( m_Cspeed_cc.getBoxes() ); 
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& rho_on_patch  = a_density[dit];   
      const FArrayBox& V_on_patch    = m_velocity[dit]; 
      const FArrayBox& P_on_patch    = m_pressure[dit]; 
      const FArrayBox& magP_on_patch = m_Bpressure[dit]; 
      const FArrayBox& N_on_patch    = m_Nmatrix[dit]; 
      const FArrayBox& J_on_patch    = m_Jacobian[dit]; 
      FArrayBox& C_on_patch = m_Cspeed_cc[dit];
      const Box& thisbox = C_on_patch.box();
      
      FORT_EVAL_FLUX_FREEZING_SPEED( CHF_BOX(thisbox), 
                                     CHF_CONST_FRA1(rho_on_patch,0),
                                     CHF_CONST_FRA1(P_on_patch,0),
                                     CHF_CONST_FRA1(magP_on_patch,0),
                                     CHF_CONST_FRA(V_on_patch),
                                     CHF_CONST_FRA(N_on_patch),
                                     CHF_CONST_FRA1(J_on_patch,0),
                                     CHF_CONST_REAL(a_gamma),
                                     CHF_CONST_INT(SpaceDim),
                                     CHF_FRA(C_on_patch) );
      
   }

}

Real IdealMhdOp::computeDtExplicitTI( const FluidSpeciesPtrVect&  a_fluid_species )
{
   CH_TIME("Real IdealMhdOp::computeDtExplicitTI()");   
   // get physical values for m_species_name and calculate freezing speed
   // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/N); i = q0, q1, q2
   //
   double gamma;
   //LevelData<FArrayBox> density, Cspeed;
   for (int species(0); species<a_fluid_species.size(); species++) {
      const FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
      const std::string species_name( fluid_species.name() );
      //species_name = fluid_species.name();
      if(species_name==m_species_name) {
         //const LevelData<FArrayBox>& soln_data( fluid_species.cell_var(0) );
         //const DisjointBoxLayout& grids( soln_data.getBoxes() );
         gamma = fluid_species.m_gamma;
         //density.define(grids, 1, IntVect::Zero);
         fluid_species.massDensity(m_rhoDen_cc);
         fluid_species.velocity(m_velocity);  // in-plane velocity vector
         fluid_species.pressure(m_pressure);
         fluid_species.Bpressure(m_Bpressure);
         //Cspeed.define(grids, SpaceDim, IntVect::Zero);
   
         // Calls for physical variables here are from mapped state
         // vector (has not been convertToPhysical )
         // and so need to divide by J
         //
         m_geometry.divideJonValid(m_rhoDen_cc);
         //m_geometry.divideJonValid(m_velocity);
         m_geometry.divideJonValid(m_pressure);
         m_geometry.divideJonValid(m_Bpressure);
         m_geometry.divideJonValid(m_Bpressure);
         
         break;
      }
   }
   setCspeed( m_rhoDen_cc, gamma );
   setCourantTimeStep(m_Cspeed_cc);   
   

   //cout << "m_courant_time_step = " << m_courant_time_step << endl;
   return m_courant_time_step;
   //return DBL_MAX;

}

void IdealMhdOp::defineLevelDatas( const DisjointBoxLayout&  a_grids,
                                   const IntVect&            a_ghostVect )
{
   // get the transformation matrix N and Jacobian at cell-center needed to transform 
   // from physical space to mapped space
   //
   m_Jacobian.define(a_grids, 1, a_ghostVect);
   m_Nmatrix.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);
   m_X_phys.define(a_grids, SpaceDim, a_ghostVect);
   m_g_y.define(a_grids, 1, a_ghostVect); // used only when solving virtual B
   m_cellVol.define(a_grids, 1, IntVect::Zero);
   double twoPi = 2.0*Pi;
   //double mapVol;
   m_geometry.getXphys(m_X_phys);
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(a_grids[dit]);
      m_twoDaxisymm = coord_sys.isAxisymmetric();
      if(m_twoDaxisymm) { // g_y is covariant basis vector in virtual direction
         m_g_y[dit].copy(m_X_phys[dit],0,0,1);
         m_g_y[dit] *= twoPi;
      } 
      else { 
         m_g_y[dit].setVal(1.0);
      }
   }
   m_geometry.getJ(m_Jacobian);
   m_geometry.getN(m_Nmatrix);
   m_geometry.getCellVolumes(m_cellVol);
   //inspectFArrayBox(Nmatrix,0);
   //inspectFArrayBox(Jacobian,0);
   //inspectFArrayBox(m_X_phys,0);
      
   m_rhoDen_cc.define(a_grids, 1, a_ghostVect);
   m_momDen_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_eneDen_cc.define(a_grids, 1, a_ghostVect);
   m_momDen_virtual_cc.define(a_grids, 1, a_ghostVect);
   m_magField_virtual_cc.define(a_grids, 1, a_ghostVect);
   //
   m_pressure.define(a_grids, 1, a_ghostVect);
   m_Bpressure.define(a_grids, 1, a_ghostVect);
   m_velocity.define(a_grids, SpaceDim, a_ghostVect);
   m_velocity_virtual.define(a_grids, 1, a_ghostVect);
   //
   m_rhoFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mxFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_myFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mzFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_enFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mvFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   //
   m_Cspeed_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_CspeedL_norm.define(a_grids, 1, a_ghostVect);
   m_CspeedR_norm.define(a_grids, 1, a_ghostVect);
   //
   m_rhoFlux_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_mxFlux_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_myFlux_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_mzFlux_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_enFlux_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_momFlux_norm.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_mvFlux_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_div_IdentFlux_R.define(a_grids, 1, IntVect::Zero);
   m_div_IdentFlux_Z.define(a_grids, 1, IntVect::Zero);
   LevelData<FluxBox> IdentFlux_phys_R( a_grids, SpaceDim, 1*IntVect::Unit );
   LevelData<FluxBox> IdentFlux_phys_Z( a_grids, SpaceDim, 1*IntVect::Unit );
   if(SpaceDim==2 && m_twoDaxisymm) {
      for (DataIterator dit(a_grids); dit.ok(); ++dit) {
         IdentFlux_phys_R[dit].setVal( 1.0,0,0,1); // val=1, dir=R, comp=R
         IdentFlux_phys_R[dit].setVal( 0.0,0,1,1); // val=0, dir=R, comp=Z
         IdentFlux_phys_R[dit].setVal( 1.0,1,0,1); // val=1, dir=Z, comp=R
         IdentFlux_phys_R[dit].setVal( 0.0,1,1,1); // val=0, dir=Z, comp=Z
         //
         IdentFlux_phys_Z[dit].setVal( 0.0,0,0,1); // val=0, dir=R, comp=R
         IdentFlux_phys_Z[dit].setVal( 1.0,0,1,1); // val=1, dir=R, comp=Z
         IdentFlux_phys_Z[dit].setVal( 0.0,1,0,1); // val=0, dir=Z, comp=R
         IdentFlux_phys_Z[dit].setVal( 1.0,1,1,1); // val=1, dir=Z, comp=Z
      }
      m_geometry.applyAxisymmetricCorrection( IdentFlux_phys_R );
      m_geometry.applyAxisymmetricCorrection( IdentFlux_phys_Z );
      m_geometry.computeMappedGridDivergence( IdentFlux_phys_R, m_div_IdentFlux_R, 0.0);
      m_geometry.computeMappedGridDivergence( IdentFlux_phys_Z, m_div_IdentFlux_Z, 0.0);
      // computeMappedGridDivergence transforms phys flux to mapped (NT*Flux) 
      // then subtracts cell fluxes in each direction and multiplies
      // by face area
   } 
      
   //virtualEleField.define(a_grids, 1, a_ghostVect);
   m_E0_ce.define(a_grids, 1, IntVect::Unit);
   m_J0_ce.define(a_grids, 1, IntVect::Unit);
   m_J0_cc.define(a_grids, SpaceDim, IntVect::Zero);
   m_E0_cc.define(a_grids, SpaceDim, IntVect::Zero);
   m_curlE0.define(a_grids, 1, IntVect::Zero);
   m_By_contr_cc.define(a_grids, 1, a_ghostVect);
   m_By_covar_cc.define(a_grids, 1, a_ghostVect);
   m_JJdotE.define(a_grids, 1, IntVect::Zero); // Jacobian*J\cdotE
   m_JJcrossB.define(a_grids, SpaceDim, IntVect::Zero); // Jacobian*J\timesB
   
   //
   m_rhoFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_mxFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_myFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_mzFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_enFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_mvFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);
   //
   m_EdgeBC_zeros.define(a_grids, 1, 1*IntVect::Unit);
   m_FluxBC_zeros.define(a_grids, 1, 1*IntVect::Unit);
   m_rhoFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_mxFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_myFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_mzFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_enFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_mvFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_momFluxBC_norm.define(a_grids, SpaceDim, 1*IntVect::Unit);
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_EdgeBC_zeros[dit].setVal(0.0,1);
      m_FluxBC_zeros[dit].setVal(0.0);
   }

}

void IdealMhdOp::setCellCenterValues( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setCellCenterValues()");
   
   const LevelData<FArrayBox>& soln_rhoDen( a_soln_fluid.cell_var(0) );
   const LevelData<FArrayBox>& soln_momDen( a_soln_fluid.cell_var("momentumDensity") );
   const LevelData<FArrayBox>& soln_eneDen( a_soln_fluid.cell_var("energyDensity") );
   
   // set values for MHD variables
   //
   const DisjointBoxLayout& grids( soln_rhoDen.getBoxes() ); 
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_rhoDen_cc[dit].copy( soln_rhoDen[dit], m_rhoDen_cc[dit].box() );
      m_momDen_cc[dit].copy( soln_momDen[dit], m_momDen_cc[dit].box() );
      m_eneDen_cc[dit].copy( soln_eneDen[dit], m_eneDen_cc[dit].box() );
   }
   if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
      const LevelData<FArrayBox>& soln_momDen_virtual( a_soln_fluid.cell_var("momentumDensity_virtual") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_momDen_virtual_cc[dit].copy( soln_momDen_virtual[dit] );
      }
   }
   if(a_soln_fluid.m_evolve_magneticField_virtual) {
      const LevelData<FArrayBox>& soln_magField_virtual( a_soln_fluid.cell_var("magneticField_virtual") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_magField_virtual_cc[dit].copy( soln_magField_virtual[dit] );
      }
   }  
   //inspectFArrayBox(m_momDen_cc,0);
   //inspectFArrayBox(m_momDen_cc,1);
   //inspectFArrayBox(m_velocity,0);
   //inspectFArrayBox(m_velocity,1);
   
   // set physical derived variables for fluid species
   //
   a_soln_fluid.pressure(m_pressure);
   a_soln_fluid.Bpressure(m_Bpressure);
   a_soln_fluid.velocity(m_velocity);                 // in-plane velocity vector
   a_soln_fluid.velocity_virtual(m_velocity_virtual); // out-of-plane velocity vector
   //inspectFArrayBox(m_Bpressure,0);
}   

void IdealMhdOp::setCellCenterFluxes( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setCellCenterFluxes()");
   
   // Compute cell-center fluxes for each variable
   //
   // rhoFlux_x = mx,         rhoFlux_y = my,         rhoFlux_z = mz;
   // mxFlux_x  = mx*ux + P,  mxFlux_y  = mx*uy,      mxFlux_z  = mx*uz
   // myFlux_x  = my*ux,      myFlux_y  = my*uy + P,  myFlux_z  = my*uz
   // mzFlux_x  = mz*ux,      mzFlux_y  = mz*uy,      mzFlux_z  = mz*uz + P
   // enFlux_x  = (en+P)*ux,  enFlux_y  = (en+P)*uy,  enFlux_z  = (en+P)*uz
   // 
   
   const DisjointBoxLayout& grids( m_rhoDen_cc.getBoxes() ); 
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      /*
      for (int n=0; n<SpaceDim; ++n) {
         m_rhoFlux_cc[dit].copy( m_momDen_cc[dit],n,n,1 );
         m_mxFlux_cc[dit].copy( m_momDen_cc[dit],0,n,1 );
         m_mxFlux_cc[dit].mult( m_velocity[dit],n,n,1 );
         m_myFlux_cc[dit].copy( m_momDen_cc[dit],1,n,1 );
         m_myFlux_cc[dit].mult( m_velocity[dit],n,n,1 );
         if (SpaceDim==3) {
            m_mzFlux_cc[dit].copy( m_momDen_cc[dit],2,n,1 );
            m_mzFlux_cc[dit].mult( m_velocity[dit],n,n,1 );
         }
         m_enFlux_cc[dit].copy( m_eneDen_cc[dit],0,n,1 );
         m_enFlux_cc[dit].plus( m_pressure[dit],0,n,1 );
         m_enFlux_cc[dit].mult( m_velocity[dit],n,n,1);
      }
      m_mxFlux_cc[dit].plus( m_pressure[dit],0,0,1 );
      m_myFlux_cc[dit].plus( m_pressure[dit],0,1,1 );
      if (SpaceDim==3) {
         m_mzFlux_cc[dit].plus( m_pressure[dit],0,2,1 );
      }
      //if(a_soln_fluid.m_evolve_magneticField_virtual) {
      //   m_mxFlux_cc[dit].plus( m_Bpressure[dit],0,0,1 );
      //   m_myFlux_cc[dit].plus( m_Bpressure[dit],0,1,1 );
      //   if (SpaceDim==3) {
      //      m_mzFlux_cc[dit].plus( m_Bpressure[dit],0,2,1 );
      //   }
      //}
      */
      if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
         for (int n=0; n<SpaceDim; ++n) {
            m_mvFlux_cc[dit].copy( m_momDen_virtual_cc[dit],0,n,1 );
            m_mvFlux_cc[dit].mult( m_velocity[dit],n,n,1);
         }
      }
      
      
      const FArrayBox& eneDen_on_patch = m_eneDen_cc[dit];   
      const FArrayBox& momDen_on_patch = m_momDen_cc[dit];   
      const FArrayBox& V_on_patch = m_velocity[dit];   
      const FArrayBox& P_on_patch = m_pressure[dit];   
      const Box& cellbox = eneDen_on_patch.box();

      FArrayBox& rhoFlux_on_patch = m_rhoFlux_cc[dit];
      FArrayBox& mxFlux_on_patch  = m_mxFlux_cc[dit];
      FArrayBox& myFlux_on_patch  = m_myFlux_cc[dit];
      FArrayBox& mzFlux_on_patch  = m_mzFlux_cc[dit];
      FArrayBox& enFlux_on_patch  = m_enFlux_cc[dit];

      FORT_CELL_CENTER_FLUXES( CHF_BOX(cellbox), 
                               CHF_CONST_FRA1(eneDen_on_patch,0),
                               CHF_CONST_FRA(momDen_on_patch),
                               CHF_CONST_FRA(V_on_patch),
                               CHF_CONST_FRA1(P_on_patch,0),
                               CHF_FRA(rhoFlux_on_patch),
                               CHF_FRA(mxFlux_on_patch),
                               CHF_FRA(myFlux_on_patch),
                               CHF_FRA(mzFlux_on_patch),
                               CHF_FRA(enFlux_on_patch) );
      
   }
}

void IdealMhdOp::setMagneticFieldTerms( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setMagneticFieldTerms()");
   
   //  EXPERIMENTAL STUFF FOR EdgeDataBox AND NodeFArrayBox STUFF
   //  NEEDED FOR MAGNETIC FIELD ON CELL FACES AND ELECTRIC FIELD
   //  AT CELL EDGES
   //
   //
   const DisjointBoxLayout& grids( m_rhoDen_cc.getBoxes() ); 
   if(a_soln_fluid.m_evolve_magneticField_virtual) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& thisbox = m_E0_cc[dit].box();
         for (int n=0; n<SpaceDim; ++n) {
            m_E0_cc[dit].copy( m_magField_virtual_cc[dit],0,n,1 ); 
         }
         m_E0_cc[dit].mult( m_velocity[dit],thisbox,1,0,1 ); // Ex =  Vz*By
         m_E0_cc[dit].mult( m_velocity[dit],thisbox,0,1,1 ); // Ez = -Vx*By
         m_E0_cc[dit].mult( -1.0,thisbox,1,1 );
      }
   
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& thisbox = m_By_covar_cc[dit].box();
         m_By_covar_cc[dit].copy( m_magField_virtual_cc[dit], thisbox); // phys
         m_By_contr_cc[dit].copy( m_magField_virtual_cc[dit], thisbox); // phys
         if(m_twoDaxisymm) {
            m_By_covar_cc[dit].mult(m_g_y[dit]);
            m_By_contr_cc[dit].divide(m_g_y[dit], thisbox, 0,0,1);
         }
      }
 
      // The following function takes covar By and returns 
      // J*contravar(curl(B)) at cell edges
      m_geometry.mappedGridCurlofVirtComp(m_By_covar_cc,m_J0_ce); // Jacobian * curl(By)\cdot g^perp
      SpaceUtils::interpEdgesToCell(m_J0_cc,m_J0_ce,"c2");
      m_geometry.convertPhysToCovar(m_J0_cc,1); // actually Jacobian times physical curl(B) at cell center
      //m_geometry.plotCellData( "m_By_covar_cc", m_By_covar_cc, a_time);
      //m_geometry.plotCellData( "m_J0_cc", m_J0_cc, a_time);
      //inspectFArrayBox(m_magField_virtual_cc,0);

      // calculate source terms for energy density and momentum density equations
      //
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& thisbox = m_J0_cc[dit].box();
         FArrayBox m_JJdotE_comp(grids[dit], SpaceDim);
         m_JJdotE_comp.copy(m_J0_cc[dit], thisbox);
         m_JJdotE_comp.mult(m_E0_cc[dit]);
         m_JJdotE[dit].copy(m_JJdotE_comp,0,0,1);
         m_JJdotE[dit].plus(m_JJdotE_comp,thisbox,1,0,1);
         //
         m_JJcrossB[dit].copy(m_J0_cc[dit],1,0,1);
         m_JJcrossB[dit].mult(m_magField_virtual_cc[dit],thisbox,0,0,1);
         m_JJcrossB[dit].negate(0,1);
         m_JJcrossB[dit].copy(m_J0_cc[dit],0,1,1);
         m_JJcrossB[dit].mult(m_magField_virtual_cc[dit],thisbox,0,1,1);
      }
      //m_geometry.plotCellData( "m_JJcrossB", m_JJcrossB, a_time);
   }
   //
   //
   //
}   

void IdealMhdOp::setFaceCenteredFluxes( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setFaceCenteredFluxes()");
   
   // Get mapped flux normals on each face multiplied by face area (NTFlux*faceArea)
   //
   //
   const DisjointBoxLayout& grids( m_rhoDen_cc.getBoxes() ); 
   LevelData<FArrayBox> thisJf(grids, 1, m_rhoDen_cc.ghostVect());
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      thisJf[dit].copy(m_rhoDen_cc[dit]);
      thisJf[dit].mult(m_Jacobian[dit],0,0,1);
   }
   computeNTFfaceArea(m_rhoFlux_norm,m_rhoFlux_cc,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   
   for (int n=0; n<SpaceDim; ++n) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         thisJf[dit].copy(m_momDen_cc[dit],n,0,1);
         thisJf[dit].mult(m_Jacobian[dit],0,0,1);
      }
      if(n==0) computeNTFfaceArea(m_mxFlux_norm,m_mxFlux_cc,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
      if(n==1) computeNTFfaceArea(m_myFlux_norm,m_myFlux_cc,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
      if(n==2) computeNTFfaceArea(m_mzFlux_norm,m_mzFlux_cc,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   }
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_momFlux_norm[dit].copy( m_mxFlux_norm[dit],0,0,1 );
      m_momFlux_norm[dit].copy( m_myFlux_norm[dit],0,1,1 );
      if(SpaceDim==3) m_momFlux_norm[dit].copy( m_mzFlux_norm[dit],0,2,1 );
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      thisJf[dit].copy(m_eneDen_cc[dit]);
      thisJf[dit].mult(m_Jacobian[dit],0,0,1);
   }
   computeNTFfaceArea(m_enFlux_norm,m_enFlux_cc,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   
   if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         thisJf[dit].copy(m_momDen_virtual_cc[dit]);
         thisJf[dit].mult(m_Jacobian[dit],0,0,1);
      }
      computeNTFfaceArea(m_mvFlux_norm,m_mvFlux_cc,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   }
   
}

void IdealMhdOp::enforceFluxBCs( FluidSpecies&  a_rhs_fluid, 
                           const FluidSpecies&  a_soln_fluid, 
                           const Real           a_time )
{
   CH_TIME("IdealMhdOp::enforceFluxBCs()");
   
   // enforce flux BC for mass density flux
   //
   SpaceUtils::upWindToFaces( m_rhoFlux_ce, m_rhoFlux_cc, m_rhoFlux_ce, "c2" );
   m_geometry.applyAxisymmetricCorrection( m_rhoFlux_ce );
   m_geometry.computeMetricTermProductAverage( m_rhoFluxBC_norm, m_rhoFlux_ce, 0 );
   m_fluid_bc.at(0)->applyFluxBC( a_rhs_fluid, m_rhoFlux_norm, m_rhoFluxBC_norm, a_time );
   
   // enforce flux BC for momentum density flux
   //
   int this_cvc;
   SpaceUtils::upWindToFaces( m_mxFlux_ce, m_mxFlux_cc, m_mxFlux_ce, "c2" );
   SpaceUtils::upWindToFaces( m_myFlux_ce, m_myFlux_cc, m_myFlux_ce, "c2" );
   m_geometry.applyAxisymmetricCorrection( m_mxFlux_ce );
   m_geometry.applyAxisymmetricCorrection( m_myFlux_ce );
   m_geometry.computeMetricTermProductAverage( m_mxFluxBC_norm, m_mxFlux_ce, 0 );
   m_geometry.computeMetricTermProductAverage( m_myFluxBC_norm, m_myFlux_ce, 0 );
   const DisjointBoxLayout& grids( m_rhoDen_cc.getBoxes() ); 
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_momFluxBC_norm[dit].copy( m_mxFluxBC_norm[dit],0,0,1 );
      m_momFluxBC_norm[dit].copy( m_myFluxBC_norm[dit],0,1,1 );
   }
   if(SpaceDim==3) {
      SpaceUtils::upWindToFaces( m_mzFlux_ce, m_mzFlux_cc, m_mzFlux_ce, "c2");
      m_geometry.applyAxisymmetricCorrection( m_mzFlux_ce );
      m_geometry.computeMetricTermProductAverage( m_mzFluxBC_norm, m_mzFlux_ce, 0 );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_momFluxBC_norm[dit].copy( m_mzFluxBC_norm[dit],0,2,1 );
      }
   }
   this_cvc = a_soln_fluid.cell_var_component("momentumDensity");
   m_fluid_bc.at(this_cvc)->applyFluxBC( a_rhs_fluid, m_momFlux_norm, m_momFluxBC_norm, a_time );
  
   // enforce flux BC for energy density flux
   //
   SpaceUtils::upWindToFaces( m_enFlux_ce, m_enFlux_cc, m_enFlux_ce, "c2" );
   m_geometry.applyAxisymmetricCorrection( m_enFlux_ce );
   m_geometry.computeMetricTermProductAverage( m_enFluxBC_norm, m_enFlux_ce, 0 );
   this_cvc = a_soln_fluid.cell_var_component("energyDensity");
   m_fluid_bc.at(this_cvc)->applyFluxBC( a_rhs_fluid, m_enFlux_norm, m_enFluxBC_norm, a_time );

   // enforce flux BC for virtual momentum density flux
   //
   if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
      SpaceUtils::upWindToFaces( m_mvFlux_ce, m_mvFlux_cc, m_mvFlux_ce, "c2" );
      m_geometry.applyAxisymmetricCorrection( m_mvFlux_ce );
      m_geometry.computeMetricTermProductAverage( m_mvFluxBC_norm, m_mvFlux_ce, 0 );
      this_cvc = a_soln_fluid.cell_var_component("momentumDensity_virtual");
      m_fluid_bc.at(this_cvc)->applyFluxBC( a_rhs_fluid, m_mvFlux_norm, m_mvFluxBC_norm, a_time );
   }
   
   // enforce edge BC for virtual magnetic field
   //
   if(a_soln_fluid.m_evolve_magneticField_virtual) {
      computeIdealEatEdges(m_E0_ce,m_velocity,m_By_contr_cc,m_Cspeed_cc);
      //m_geometry.convertPhysToCovar(m_E0_cc,0);
      //SpaceUtils::interpCellToEdges(m_E0_ce, m_E0_cc, m_E0_ce, "c2"); // Need to upwind using Cspeed
      this_cvc = a_soln_fluid.cell_var_component("magneticField_virtual");
      m_fluid_bc.at(this_cvc)->applyEdgeBC( a_rhs_fluid, m_E0_ce, m_EdgeBC_zeros, a_time );
      m_geometry.mappedGridCurlFromEdgeTans(m_E0_ce,2,m_curlE0); // takes covar E and returns J*curl(E)dot(g^y)
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         if(m_twoDaxisymm) {
            m_curlE0[dit].mult(m_g_y[dit]); // g^y*g_y = 1
         }
      }
   }
   //inspectEdgeDataBox(m_E0_ce,1);
   //inspectFluxBox(m_rhoFlux_norm,0);
   //inspectFluxBox(m_momFlux_norm,0);
   //inspectFluxBox(m_enFlux_norm,0);
}   

void IdealMhdOp::updateRHSs( FluidSpecies&  a_rhs_fluid,
                       const FluidSpecies&  a_soln_fluid ) 
{
   CH_TIME("IdealMhdOp::updateRHSs()");
   
   // compute divergence of Fluxes and add (subtract) to RHS 
   // Note that J/cellVol = 1/mapVol
   //
   const DisjointBoxLayout& grids( m_rhoDen_cc.getBoxes() ); 

   // update RHS for mass density
   //
   LevelData<FArrayBox>& rhs_data( a_rhs_fluid.cell_var(0) );
   for (DataIterator dit(rhs_data.dataIterator()); dit.ok(); ++dit) {
      FArrayBox div_rhoFlux(grids[dit], 1);
      m_geometry.mappedGridDivergenceFromFluxNorms(m_rhoFlux_norm[dit], div_rhoFlux);
      div_rhoFlux /= m_cellVol[dit];
      div_rhoFlux.mult(m_Jacobian[dit]);
      rhs_data[dit].minus(div_rhoFlux);
   }
  
   // update RHS for momentum density
   // 
   LevelData<FArrayBox>& rhs_data2( a_rhs_fluid.cell_var("momentumDensity") );
   for (DataIterator dit(rhs_data2.dataIterator()); dit.ok(); ++dit) {
      FArrayBox div_momFlux(grids[dit], SpaceDim);
      m_geometry.mappedGridDivergenceFromFluxNorms(m_momFlux_norm[dit], div_momFlux);
      for (int n=0; n<div_momFlux.nComp(); ++n) {
         if(n==0 && m_twoDaxisymm) {
            //if(a_soln_fluid.m_evolve_magneticField_virtual) {
            //   m_pressure[dit].minus(m_Bpressure[dit],0,0,1);
            //}
            FArrayBox pdivIdent_R(grids[dit],1);
            FArrayBox pdivIdent_Z(grids[dit],1);
            pdivIdent_R.copy( m_div_IdentFlux_R[dit] );
            pdivIdent_Z.copy( m_div_IdentFlux_Z[dit] );
            pdivIdent_R.mult(m_pressure[dit],0,0,1);
            pdivIdent_Z.mult(m_pressure[dit],0,0,1);
            div_momFlux.minus(pdivIdent_R,0,0,1);
            div_momFlux.minus(pdivIdent_Z,0,1,1);
         }
         div_momFlux.divide(m_cellVol[dit],0,n,1);
         div_momFlux.mult(m_Jacobian[dit],0,n,1);
      }
      //div_momFlux.setVal(0.0,0); // JRA, no mr update
      rhs_data2[dit].minus(div_momFlux);
      if(m_twoDaxisymm && a_soln_fluid.m_evolve_momentumDensity_virtual) {
         FArrayBox JaRhoVthsqOverR(grids[dit], 1);
         JaRhoVthsqOverR.copy(m_momDen_virtual_cc[dit]);
         JaRhoVthsqOverR.mult(m_velocity_virtual[dit]);
         JaRhoVthsqOverR.mult(m_Jacobian[dit]); // Jacobian = Ja = 2*pi*R_phys * J_RZ
         JaRhoVthsqOverR.divide(m_X_phys[dit],0,0,1);
         rhs_data2[dit].plus(JaRhoVthsqOverR,0,0,1);            
      }
      if(a_soln_fluid.m_evolve_magneticField_virtual) {
         rhs_data2[dit].plus(m_JJcrossB[dit]);
      }
   }
   //inspectFluxBox(m_momFlux_norm,0);
   
   // update RHS for energy density
   //
   LevelData<FArrayBox>& rhs_data3( a_rhs_fluid.cell_var("energyDensity") );
   for (DataIterator dit(rhs_data3.dataIterator()); dit.ok(); ++dit) {
      FArrayBox div_enFlux(grids[dit], 1);
      m_geometry.mappedGridDivergenceFromFluxNorms(m_enFlux_norm[dit], div_enFlux);
      div_enFlux /= m_cellVol[dit];
      div_enFlux.mult(m_Jacobian[dit]);
      rhs_data3[dit].minus( div_enFlux );
      if(a_soln_fluid.m_evolve_magneticField_virtual) {
         rhs_data3[dit].plus(m_JJdotE[dit]);
      }
   }
   
   // update RHS for momentum density in virtual direction
   //
   if(a_soln_fluid.m_evolve_momentumDensity_virtual){
      LevelData<FArrayBox>& rhs_data4( a_rhs_fluid.cell_var("momentumDensity_virtual") );
      for (DataIterator dit(rhs_data4.dataIterator()); dit.ok(); ++dit) {
         FArrayBox div_mvFlux(grids[dit], 1);
         m_geometry.mappedGridDivergenceFromFluxNorms(m_mvFlux_norm[dit], div_mvFlux);
         div_mvFlux /= m_cellVol[dit];
         div_mvFlux.mult(m_Jacobian[dit]);
         rhs_data4[dit].minus( div_mvFlux );
         if(m_twoDaxisymm) {
            FArrayBox JaRhoVthVrOverR(grids[dit], 1);
            JaRhoVthVrOverR.copy(m_momDen_virtual_cc[dit]);
            JaRhoVthVrOverR.mult(m_velocity[dit],0,0,1);
            JaRhoVthVrOverR.mult(m_Jacobian[dit]); // Jacobian = Ja = 2*pi*R_phys * J_RZ
            JaRhoVthVrOverR.divide(m_X_phys[dit],0,0,1);
            rhs_data4[dit].minus(JaRhoVthVrOverR,0,0,1);            
         }
      }
   }

   // update RHS for virtual magnetic field
   //
   if(a_soln_fluid.m_evolve_magneticField_virtual){
      LevelData<FArrayBox>& rhs_data5( a_rhs_fluid.cell_var("magneticField_virtual") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs_data5[dit].minus( m_curlE0[dit] );
      }
   }

}   

#include "NamespaceFooter.H"
