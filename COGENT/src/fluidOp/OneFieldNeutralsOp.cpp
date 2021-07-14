#include "OneFieldNeutralsOp.H"
#include "EllipticOpBCFactory.H"
#include "FluidVarBCFactory.H"
#include "FourthOrderUtil.H"
#include "SimpleDivergence.H"
#include "inspect.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "PhaseGeom.H"
#include "MomentOp.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NeutralsOpF_F.H"

#include "NamespaceHeader.H" 

const char* OneFieldNeutralsOp::pp_name = {"one_field_neutrals_op"};


OneFieldNeutralsOp::OneFieldNeutralsOp(const string&   a_pp_str,
                                       const string&   a_species_name,
                                       const MagGeom&  a_geometry,
                                       const Real      a_larmor,
                                       const int       a_verbosity )
   : m_Tg_func(NULL),
     m_larmor(a_larmor),
     m_consistent_diffusion(false),
     m_D_rad_func(NULL),
     m_D_par_func(NULL),
     m_D_perp_func(NULL),
     m_ne_func(NULL),
     m_Te_func(NULL),
     m_geometry(a_geometry),
     m_verbosity(a_verbosity),
     m_is_time_implicit(true),
     m_first_call(true),
     m_include_ionization(false),
     m_opt_string(a_pp_str),
     m_my_pc_idx(-1)
{
   ParmParse pp(a_pp_str.c_str());
   parseParameters( pp );
   if (m_verbosity>0) {
      printParameters();
   }

   ParmParse pp_vlasov("gkvlasov");
   m_vlasov = new PS::GKVlasov(pp_vlasov, a_larmor);
   
   const std::string name("neutrals_diffusion");
   ParmParse ppOp( name.c_str() );
   
   const std::string prefix( "BC." + name );
   ParmParse ppOpBC( prefix.c_str() );

   m_diffusion_op = new Diffusion(ppOp, a_geometry);

   EllipticOpBCFactory elliptic_op_bc_factory;
   m_bc = elliptic_op_bc_factory.create( name,
                                         ppOpBC,
                                         *(a_geometry.getCoordSys()),
                                         false );

   string variable_name = "density";

   // Input the initial conditions
   parseInitialConditions(a_species_name, variable_name);

   // Input the boundary conditions
   FluidVarBCFactory fluid_var_bc_factory;
   m_fluid_variable_bc = fluid_var_bc_factory.create(a_species_name, variable_name, a_geometry.getCoordSys()->type(), false);
   m_dens_recycling_bc = m_fluid_variable_bc->isRecyclingBC();
   
   // Initialize objects
   const DisjointBoxLayout& grids = m_geometry.grids();

   m_volume.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume);
   
   if (m_consistent_density || m_consistent_temperature) {
      // Try to optimize expensive moment calculation
      // by computing it in valid cells only and then
      // by using extrapolation routins at phys. boundaries
      m_ghostVect = IntVect::Zero;
   }
   else {
      m_ghostVect = 2*IntVect::Unit;
   }
   
}


OneFieldNeutralsOp::~OneFieldNeutralsOp()
{
   delete m_fluid_variable_bc;
   delete m_bc;
   delete m_diffusion_op;
   delete m_vlasov;
}


void OneFieldNeutralsOp::accumulateRHS(FluidSpeciesPtrVect&               a_rhs,
                                       const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                       const FluidSpeciesPtrVect&         a_fluid_species,
                                       const PS::ScalarPtrVect&           a_scalars,
                                       const EField&                      a_E_field,
                                       const int                          a_fluidVecComp,
                                       const Real                         a_time)
{
   /*
    This class computes the RHS for (NB: rho_g is the mass density):
      drho_g/dt = \nabla( D \nabla(rho_g T_g)) - nu_iz * rho_g + Source
   */
   
   // Get fluid rhs
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   LevelData<FArrayBox>& rhs_data( rhs_fluid.cell_var(0) );

   // Get fluid soln
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
   const LevelData<FArrayBox>& soln_data( soln_fluid.cell_var(0) );

   //Get species mass
   const Real species_mass = soln_fluid.mass();

   // Get geometry
   const DisjointBoxLayout& grids( rhs_data.getBoxes() );

   // Get electron density (n_e), temperature (T_e)
   // and neutral temperature (T_g)
   getPlasmaParameters(a_time);
   
   // Compute diffusion coefficient
   LevelData<FluxBox> D_tensor(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   LevelData<FluxBox> D_tensor_mapped(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   computeDiffusionCoefficients(D_tensor, D_tensor_mapped, species_mass, a_time);

   // Set diffisionOf coefficients
   m_diffusion_op->setOperatorCoefficients(D_tensor, D_tensor_mapped, *m_bc);

   // Get neutral pressure (times mass) P_g = rho_g * T_g
   LevelData<FArrayBox> neutr_pressure(grids, 1, IntVect::Zero);
   for (DataIterator dit( neutr_pressure.dataIterator() ); dit.ok(); ++dit) {
      neutr_pressure[dit].copy(soln_data[dit]);
      neutr_pressure[dit].mult(m_Tg[dit]);
   }
   
   // Compute diffusion flux divergence
   LevelData<FArrayBox> flux_div(grids, 1, IntVect::Zero);
   
   if (!m_dens_recycling_bc) {
      m_diffusion_op->computeFluxDivergence( neutr_pressure, flux_div, false, false);
   }
   else {
      computeDiffusionWithRecyclingBC(flux_div, *m_fluid_variable_bc, neutr_pressure,
                                      *m_diffusion_op, m_ion_normal_flux, a_time);
   }
   
   // The preceding call computes a negated physical flux divergence cell average.
   // Since this function is being used to define the right-hand side of an evolution
   // equation for a mapped unknown (i.e., containing a J factor), we therefore need
   // to multiply by J.  We also flip the sign.
   m_geometry.multJonValid(flux_div);
   
   for (DataIterator dit( flux_div.dataIterator() ); dit.ok(); ++dit) {
      rhs_data[dit] -= flux_div[dit];
   }

   // Add ionization
   if (m_include_ionization) {
      
      LevelData<FArrayBox> ionization_rate(grids, 1, IntVect::Zero);
      computeIonizationRate(ionization_rate, a_time);
      
      LevelData<FArrayBox> ionization_term(grids, 1, IntVect::Zero);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         ionization_term[dit].copy(ionization_rate[dit]);
         ionization_term[dit].mult(soln_data[dit]);
      }

      m_geometry.multJonValid(ionization_term);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs_data[dit].minus(ionization_term[dit]);
      }
   }
   
   //Add fixed source
   if (m_fixed_source_func != NULL) {
     LevelData<FArrayBox> source_dst(grids, 1, IntVect::Zero);
     m_fixed_source_func->assign( source_dst, m_geometry, a_time);

     m_geometry.multJonValid(source_dst);

     for (DataIterator dit(grids); dit.ok(); ++dit) {
       rhs_data[dit].plus(source_dst[dit]);
     }
   }

   m_first_call = false;
}

void OneFieldNeutralsOp::accumulateExplicitRHS(FluidSpeciesPtrVect&               a_rhs,
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


void OneFieldNeutralsOp::accumulateImplicitRHS(FluidSpeciesPtrVect&               a_rhs,
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

void OneFieldNeutralsOp::preOpEval(const PS::KineticSpeciesPtrVect&   a_kinetic_species,
                                   const FluidSpeciesPtrVect&         a_fluid_species,
                                   const PS::ScalarPtrVect&           a_scalars,
                                   const EField&                      a_E_field,
                                   const double                       a_time )
{
   CH_TIME("OneFieldNeutralsOp::preOpEval");
   
   const DisjointBoxLayout& grids = m_geometry.grids();
   
   // The moment calculations used below compute moment including the ghost-cell regions
   // of the kinent distritbuion function, and some of them (e.g., temperature) expect
   // (through the use of intermidiate variables) that cfg objects have the same number of ghosts.
   // This can be fixed by using direct functions from momentOp class if needed for performance.
   
   IntVect ghostVectMom(IntVect::Zero);
   if (m_consistent_density || m_consistent_temperature) {
     const PS::KineticSpecies& master_species( *(a_kinetic_species[0]) );
     const PS::PhaseGeom& phase_geom = master_species.phaseSpaceGeometry();
     const PS::LevelData<PS::FArrayBox>& master_dfn( master_species.distributionFunction() );
     PS::IntVect phase_ghosts = master_dfn.ghostVect();
     ghostVectMom = phase_geom.config_restrict(phase_ghosts);
   }
   
   // Compute electron density (by summing charge density from all speceis)
   if (m_consistent_density) {
      if (!m_ne.isDefined()) m_ne.define(grids, 1, ghostVectMom);
      computeIonChargeDensity( m_ne, a_kinetic_species );
   }
  
   //Compute temperature of the zero species (assume it to be the main species)
   if (m_consistent_temperature) {
      if (!m_Tg.isDefined()) m_Tg.define(grids, 1, ghostVectMom);
      const PS::KineticSpecies& this_species( *(a_kinetic_species[0]) );
      this_species.temperature( m_Tg );
   }
   
   // Compute normal particle flux thourgh a cell face (NT * Gamma_phys * dA_mapped )
   if (m_dens_recycling_bc) {
      if (!m_ion_normal_flux.isDefined()) m_ion_normal_flux.define(grids, 1, IntVect::Zero);
      computeIntegratedNormalIonParticleFlux( m_ion_normal_flux, a_kinetic_species, a_E_field, a_time );
   }
}

void OneFieldNeutralsOp::computeIonChargeDensity(LevelData<FArrayBox>&               a_ion_charge_density,
                                                 const PS::KineticSpeciesPtrVect&    a_species ) const
{
   CH_TIME("OneFieldNeutralsOp::computeIonChargeDensity");
   
   // Container for individual species charge density
   LevelData<FArrayBox> ion_species_charge;
   ion_species_charge.define(a_ion_charge_density);
   
   setZero( a_ion_charge_density );
   
   for (int species(0); species<a_species.size(); species++) {
      
      const PS::KineticSpecies& this_species( *(a_species[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      // Compute the charge density for this species including ghost cells
      this_species.chargeDensity( ion_species_charge );
      
      for (DataIterator dit(a_ion_charge_density.dataIterator()); dit.ok(); ++dit) {
        a_ion_charge_density[dit] += ion_species_charge[dit];
      }
   }
}

void OneFieldNeutralsOp::computeIntegratedNormalIonParticleFlux(LevelData<FluxBox>&                a_ion_particle_flux,
                                                                const PS::KineticSpeciesPtrVect&   a_species,
                                                                const EField&                      a_E_field,
                                                                const Real&                        a_time) const
{

   /*
    Compute normal particle flux thourgh a cell face (NT * Gamma_phys * dA_mapped )
    */
   
   CH_TIME("OneFieldNeutralsOp::computeNormalIonParticleFlux");
   
   const DisjointBoxLayout& grids = m_geometry.grids();

   IntVect ghosts = (m_geometry.secondOrder()) ? IntVect::Zero : IntVect::Unit;
   
   LevelData<FluxBox> ion_species_mapped_flux(grids, 1, ghosts);
   
   setZero( a_ion_particle_flux );
   
   for (int species(0); species<a_species.size(); species++) {
      
      const PS::KineticSpecies& this_species( *(a_species[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      //Gyroaverage calculation requires phi. Pass zero phi for now. Fix later
      LevelData<FArrayBox> phi_tmp(grids, 1, IntVect::Zero);
      setZero(phi_tmp);
      m_vlasov->computeIntegratedMomentFluxNormals(ion_species_mapped_flux, this_species, phi_tmp, a_E_field,
                                                   PS::PhaseGeom::FULL_VELOCITY, "particle", a_time);

      for (DataIterator dit(a_ion_particle_flux.dataIterator()); dit.ok(); ++dit) {
         a_ion_particle_flux[dit] += ion_species_mapped_flux[dit];
      }
   }
}

void OneFieldNeutralsOp::defineBlockPC(std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>& a_pc,
                                       std::vector<PS::DOFList>&                                    a_dof_list,
                                       const PS::ODEVector&                                         a_soln_vec,
                                       void*                                                        a_gkops,
                                       const std::string&                                           a_out_string,
                                       const std::string&                                           a_opt_string,
                                       bool                                                         a_im,
                                       const FluidSpecies&                                          a_fluid_species,
                                       const PS::GlobalDOFFluidSpecies&                             a_global_dofs,
                                       int                                                          a_species_idx )
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


void OneFieldNeutralsOp::updateBlockPC(std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>& a_pc,
                                       const PS::KineticSpeciesPtrVect&                             a_kin_species_phys,
                                       const FluidSpeciesPtrVect&                                   a_fluid_species,
                                       const Real                                                   a_time,
                                       const int                                                    a_step,
                                       const int                                                    a_stage,
                                       const Real                                                   a_shift,
                                       const bool                                                   a_im,
                                       const int                                                    a_species_idx )
{
  if (a_im && m_is_time_implicit) {
    CH_assert(m_my_pc_idx >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx);
  
    if (!procID()) {
      std::cout << "    ==> Updating " << _FLUID_OP_PC_ << " preconditioner " 
                << " for OneFieldNeutralsOp of fluid species " << a_species_idx << ".\n";
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
                a_species_idx );
  }
  return;
}


void OneFieldNeutralsOp::updatePCImEx(const FluidSpeciesPtrVect&       a_fluid_species,
                                      const PS::KineticSpeciesPtrVect& a_kinetic_species,
                                      const double                     a_time,
                                      const int                        a_step,
                                      const int                        a_stage,
                                      const double                     a_shift,
                                      const int                        a_component)
{
   CH_TIME("OneFieldNeutralsOp::updatePCImEx");
   
   //Update linear coefficient beta = nu_iz + 1/(dt * a_ii)
   //in the precond diffusion eqution solved by hypre
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   LevelData<FArrayBox> beta(grids, 1, IntVect::Zero);
   for (DataIterator dit(beta.dataIterator()); dit.ok(); ++dit) {
      beta[dit].setVal(a_shift);
   }

   if (m_include_ionization) {
      LevelData<FArrayBox> ionization_rate(grids, 1, IntVect::Zero);
      computeIonizationRate(ionization_rate, a_time);

      for (DataIterator dit(beta.dataIterator()); dit.ok(); ++dit) {
         beta[dit].plus(ionization_rate[dit]);
      }
   }
   
   for (DataIterator dit(beta.dataIterator()); dit.ok(); ++dit) {
      beta[dit].divide(m_Tg[dit]);
   }
   
   m_diffusion_op->updateImExPreconditioner( beta, *m_bc );
}


void OneFieldNeutralsOp::solvePCImEx(FluidSpeciesPtrVect&              a_fluid_species_solution,
                                     const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                     const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                     const int                         a_component )
{
   CH_TIME("OneFieldNeutralsOp::solvePCImEx");
   if (!m_is_time_implicit) {
     MayDay::Error("Error in OneFieldNeutralsOp::solvePCImEx(): m_is_time_implicit is false. This function shouldn't have been called at all!");
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
   
   m_diffusion_op->solveImExPreconditioner(r_phys, z);
   
   // Convert the elliptic solutions to the state vector components (mapped space)
   for (DataIterator dit(z.dataIterator()); dit.ok(); ++dit) {
      z[dit].divide(m_Tg[dit]);
   }
   m_geometry.multJonValid(z);

}


void OneFieldNeutralsOp::computeDiffusionCoefficients(LevelData<FluxBox>& a_D_tensor,
                                                      LevelData<FluxBox>& a_D_tensor_mapped,
                                                      const Real          a_mass,
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
   LevelData<FArrayBox> D_cell(grids, 1, 2*IntVect::Unit);
   LevelData<FluxBox> D_face(grids, 1, 2*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     a_D_tensor[dit].setVal(0.0);
     a_D_tensor_mapped[dit].setVal(0.0);
   }   

   //Add consistent isotropic diffusion
   //D = 1/(mg * ne * <sigma_cx*Vti>)
   
   if (m_consistent_diffusion) {

      LevelData<FArrayBox> chx_rate(grids, 1, m_ghostVect);
      computeChargeExchangeRate(chx_rate, a_time);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_cell[dit].copy(chx_rate[dit]);
         D_cell[dit].mult(a_mass);
         D_cell[dit].invert(1.0);
      }
      
      if (m_ghostVect == IntVect::Zero) {
         // Chx-exchange rate is only defined in valid region, so we
         // extrapolate to physical ghosts to get D_coeff on boundary faces
         bool fourth_order = (m_geometry.secondOrder()) ? false : true;
         m_geometry.extrapolateToPhysicalGhosts(D_cell, fourth_order);
      }

      // Get diffusion coefficient on faces
      convertCellToFace(D_face, D_cell);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_tmp[dit].copy(perp_coeff[dit]);
         D_tmp[dit]+= par_coeff[dit];
         D_tmp_mapped[dit].copy(perp_coeff_mapped[dit]);
         D_tmp_mapped[dit] += par_coeff_mapped[dit];
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
   
   //Add prescribed perpendicular diffusion term, if any
   
   if (m_D_perp_func != NULL) {
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_cell[dit].copy(m_D_perp[dit]);
      }
      
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
   
   //Add prescribed parallel diffusion term, if any
   
   if (m_D_par_func != NULL) {
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_cell[dit].copy(m_D_par[dit]);
      }
      
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
   
   //Add prescribed radial diffusion term, if any
   
   if (m_D_rad_func != NULL) {
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_cell[dit].copy(m_D_perp[dit]);
      }
      
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

void OneFieldNeutralsOp::getPlasmaParameters(const Real a_time)
{
   const DisjointBoxLayout& grids = m_geometry.grids();
   
   if (!m_ne.isDefined()) m_ne.define(grids, 1, m_ghostVect);
   if (!m_Te.isDefined()) m_Te.define(grids, 1, m_ghostVect);
   if (!m_Tg.isDefined()) m_Tg.define(grids, 1, m_ghostVect);
   
   //Initialize plasma background (if prescribed)
   if (m_ne_func != NULL && m_first_call) {
      m_ne_func->assign( m_ne, m_geometry, a_time);
   }
   
   if (m_Te_func != NULL && m_first_call) {
      m_Te_func->assign( m_Te, m_geometry, a_time);
   }

   if (m_Tg_func != NULL && m_first_call) {
      m_Tg_func->assign( m_Tg, m_geometry, a_time);
   }
   
   //Set prescribed diffusion coefficients
   // 4ghost is likely too many (fix later)
   if (m_D_par_func != NULL && m_first_call) {
      m_D_par.define(grids, 1, 2*IntVect::Unit);
      m_D_par_func->assign( m_D_par, m_geometry, a_time);
   }
   if (m_D_perp_func != NULL && m_first_call) {
      m_D_perp.define(grids, 1, 2*IntVect::Unit);
      m_D_perp_func->assign( m_D_perp, m_geometry, a_time);
   }
   if (m_D_rad_func != NULL && m_first_call) {
      m_D_rad.define(grids, 1, 2*IntVect::Unit);
      m_D_rad_func->assign( m_D_rad, m_geometry, a_time);
   }
}

void OneFieldNeutralsOp::computeIonizationRate(LevelData<FArrayBox>&   a_ionization_rate,
                                               const Real              a_time) const
{
   
   /*
      This function computes ionization
      frequency (in COGENT units)
    */
   
   //Universal constants (in CGS)
   double mp = 1.6726e-24;
   
   //Get normalization parameters (units)
   double N, T, L;
   ParmParse ppunits( "units" );
   ppunits.get("number_density",N);  //[m^{-3}]
   ppunits.get("temperature",T);     //[eV]
   ppunits.get("length",L);          //[m]
   
   double Tcgs = 1.602e-12 * T; //[erg]
   double Lcgs  = 1.0e2 * L;   //[cm]
   
   double time_norm = Lcgs / sqrt(Tcgs/mp); //[s]

   const DisjointBoxLayout& grids( m_geometry.grids() );

   LevelData<FArrayBox> a2(grids, 1, IntVect::Zero);
   double fac = pow(T/10.0,2);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a2[dit].copy(m_Te[dit]);
      a2[dit].mult(m_Te[dit]);
      a2[dit].mult(fac);
   }
   
   //Get <sigmaV> in m^3/s
   LevelData<FArrayBox> sigmaV(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      sigmaV[dit].copy(a2[dit]);
      sigmaV[dit].plus(3.0);
      sigmaV[dit].divide(a2[dit]);
      sigmaV[dit].invert(1.0);
      sigmaV[dit].mult(3.0e-14);
   }
   
   //Compute ionization frequency in COGENT time units
   double norm = time_norm * N;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_ionization_rate[dit].copy(sigmaV[dit]);
      a_ionization_rate[dit].mult(m_ne[dit]);
      a_ionization_rate[dit].mult(norm);
   }
}

void OneFieldNeutralsOp::computeChargeExchangeRate(LevelData<FArrayBox>&   a_chx_rate,
                                                   const Real              a_time) const
{
   
   /*
      This function computes charge-exchange
      frequency (in COGENT units)
   */
   
   //Universal constants (in CGS)
   double mp = 1.6726e-24;
   
   //Get normalization parameters (units)
   double N, T, L;
   ParmParse ppunits( "units" );
   ppunits.get("number_density",N);  //[m^{-3}]
   ppunits.get("temperature",T);     //[eV]
   ppunits.get("length",L);          //[m]
   
   double Tcgs = 1.602e-12 * T; //[erg]
   double Lcgs  = 1.0e2 * L;   //[cm]
   
   double time_norm = Lcgs / sqrt(Tcgs/mp); //[s]

   const DisjointBoxLayout& grids( m_geometry.grids() );

   //Compute charge-exchange normalization factor

   double norm_chx;
   if (m_uedge_chx) {
     //UEDGE expression 
     norm_chx = time_norm * 5.0e-15 * sqrt(T) * N;
   }
   else {
     //Expression (27) from Rozhansky et al NF, 49, 025007 (2009)                                                             
     norm_chx = time_norm * 1.98e-14 * sqrt(T) * N;   
   }
      
   //Get sqrt(Ti); assume T_i=T_g
   LevelData<FArrayBox> sqrtTi(grids, 1, m_ghostVect);
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      FORT_COMPUTE_SQRT_FUNC(CHF_BOX(sqrtTi[dit].box()),
                             CHF_CONST_FRA1(m_Tg[dit],0),
                             CHF_FRA1(sqrtTi[dit],0));
   }
   
   //Compute charge-exchange frequency in COGENT time units
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_chx_rate[dit].copy(m_ne[dit]);
      a_chx_rate[dit].mult(sqrtTi[dit]);
      a_chx_rate[dit].mult(norm_chx);
   }
}


void OneFieldNeutralsOp::computeDiffusionWithRecyclingBC(LevelData<FArrayBox>&        a_out,
                                                         FluidVarBC&                  a_recycl_var_bc,
                                                         const LevelData<FArrayBox>&  a_in,
                                                         const Diffusion&             a_diffusion_op,
                                                         const LevelData<FluxBox>&    a_ion_normal_flux,
                                                         const Real                   a_time) const
{
   // This is a (second-order) low-polution option of EllipticOp::computeFLuxDivergence,
   // but with the option to assign flux BC (from recycling)
   
   CH_TIME("EllipticOp::computeFluxDivergence");
   const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

   LevelData<FArrayBox> phi(grids, 1, 2*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi[dit].copy(a_in[dit]);
   }

   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);

   if (SpaceDim == 3) {
      a_diffusion_op.computeMapped3DFieldWithBCs(phi, flux, false);
   }

   else  {
      a_diffusion_op.computeMappedPoloidalFieldWithBCs(phi, flux, false);
   }

   // Multiply the field by the mapped, face-centered GKP coefficients
   a_diffusion_op.multiplyCoefficients(flux, true);

   m_geometry.fillTransversePhysicalGhosts(flux);

   // Convert to face-averaged
   if ( !(a_diffusion_op.secondOrder()) ) fourthOrderAverage(flux);

   m_geometry.averageAtBlockBoundaries(flux);

   LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
   for (DataIterator dit(flux.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect faceArea = block_coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; ++dir) {
         NTF_normal[dit][dir].copy(flux[dit][dir],dir,0,1);
         NTF_normal[dit][dir].mult(faceArea[dir]);
      }
   }
   
   // Apply recycling BC
   a_recycl_var_bc.applyRecyclingBC(NTF_normal, a_ion_normal_flux, a_time);

   RealVect fakeDx = RealVect::Unit;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      simpleDivergence(a_out[dit], NTF_normal[dit], grids[dit], fakeDx);
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_out[dit] /= m_volume[dit];
   }
}

void OneFieldNeutralsOp::fillGhostCells(FluidSpecies&  a_species_phys,
                                        const double   a_time )
{
   //Alghough diffusion operator has its own methods for BC specification
   // we need to make a call to apply() to create boundary structures
   // (perhaps a single call would be sufficient)
   if (m_dens_recycling_bc) {
      // Fill ghost cells except for those on physical boundaries
      LevelData<FArrayBox>& fld( a_species_phys.cell_var(0) );
      m_geometry.fillInternalGhosts( fld );
      // Fill ghost cells on physical boundaries
      m_fluid_variable_bc->apply( a_species_phys, a_time );
   }
}

void OneFieldNeutralsOp::convertCellToFace(LevelData<FluxBox>&         a_face_data,
                                           const LevelData<FArrayBox>& a_cell_data) const
{
   const DisjointBoxLayout& grids = m_geometry.grids();

   int order = (m_geometry.secondOrder()) ? 2 : 4;

   if (order == 4) {
      CH_assert(a_cell_data.ghostVect() >= 2*IntVect::Unit );
   }
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {

       for (int dir=0; dir<SpaceDim; dir++) {
       
       Box box(grids[dit]);
          
       if (order == 4 && (a_face_data.ghostVect() >= IntVect::Unit)) {
          for (int tdir(0); tdir<SpaceDim; tdir++) {
             if (tdir!=dir) {
                const int TRANSVERSE_GROW(1);
                box.grow( tdir, TRANSVERSE_GROW );
             }
          }
       }
       
       SpaceUtils::faceInterpolate(dir,
                                   surroundingNodes(box,dir),
                                   order,
                                   a_cell_data[dit],
                                   a_face_data[dit][dir] );
       }
    }
}

void OneFieldNeutralsOp::setZero( LevelData<FArrayBox>& a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit].setVal(0.);
   }
}

void OneFieldNeutralsOp::setZero( LevelData<FluxBox>& a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit].setVal(0.);
   }
}

void OneFieldNeutralsOp::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "time_implicit", m_is_time_implicit);

   GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
   std::string grid_function_name;
   
   if (a_pp.contains("D_rad")) {
     a_pp.get("D_rad", grid_function_name );
     m_D_rad_func = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("D_perp")) {
      a_pp.get("D_perp", grid_function_name );
      m_D_perp_func = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("D_par")) {
      a_pp.get("D_par", grid_function_name );
      m_D_par_func = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("fixed_source")) {
     a_pp.get("fixed_source", grid_function_name );
     m_fixed_source_func = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("electron_temperature")) {
      a_pp.get( "electron_temperature", grid_function_name );
      m_Te_func = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("electron_density")) {
      a_pp.get( "electron_density", grid_function_name );
      m_ne_func = grid_library->find( grid_function_name );
      m_consistent_density = false;
   }
   else {
      m_consistent_density = true;
   }

   if (a_pp.contains("neutral_temperature")) {
      a_pp.get( "neutral_temperature", grid_function_name );
      m_Tg_func = grid_library->find( grid_function_name );
      m_consistent_temperature = false;
   }
   else {
      m_consistent_temperature = true;
   }  
   
   a_pp.query("consistent_diffusion", m_consistent_diffusion);
   a_pp.query("include_ionization", m_include_ionization);
   a_pp.query("consistent_density", m_consistent_density);
   
   if (a_pp.contains("use_uedge_chx_coefficient")) {
     a_pp.get( "use_uedge_chx_coefficient", m_uedge_chx );
   }
   else {
     m_uedge_chx = false;
   }
   
   //Check input for consistency
   if (m_include_ionization) {
      if (m_Te_func == NULL) {
         MayDay::Error("Electron temperature must be specified to compute ionization");
      }
      
      if (!m_consistent_density && m_ne_func == NULL) {
         MayDay::Error("Electron density must be specified or set consistent to compute ionization");
      }
   }

   if (m_consistent_diffusion && !m_consistent_density && m_ne_func == NULL) {
     MayDay::Error("Electron density must be specified or set consistent to compute diffusion coefficient");
   }

   if (m_consistent_density && m_ne_func != NULL) {
     MayDay::Error("Please use either consistent or prescribed plasmas density");
   }
}


void OneFieldNeutralsOp::printParameters()
{
   if (procID()==0) {
      std::cout << "OneFieldNeutralsOp parameters:" << std::endl;
      std::cout << "  time_implicit  =  " << m_is_time_implicit << std::endl;
      std::cout << "  consistent_diffusion  =  " << m_consistent_diffusion << std::endl;
      std::cout << "  consistent_temperature  =  " << m_consistent_temperature << std::endl;
      std::cout << "  consistent_density  =  " << m_consistent_density << std::endl;
      std::cout << "  include_ionization  =  " << m_include_ionization << std::endl;
   }
}

Real OneFieldNeutralsOp::computeDtExplicitTI( const FluidSpeciesPtrVect&  a_fluid_species )
{
   // Replace with an estimate of the dt
   return DBL_MAX;
}

Real OneFieldNeutralsOp::computeDtImExTI( const FluidSpeciesPtrVect&  a_fluid_species )
{
   if (m_is_time_implicit) {
     return DBL_MAX;
   } else {
     return computeDtExplicitTI(a_fluid_species);
   }
}


#include "NamespaceFooter.H"
