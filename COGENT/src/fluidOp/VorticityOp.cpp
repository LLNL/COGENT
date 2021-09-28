#include "VorticityOp.H"
#include "Directions.H"

#include "EllipticOpBCFactory.H"
#include "LogRectCoordSys.H"
#include "SNCoreCoordSys.H"
#include "FourthOrderUtil.H"
#include "SimpleDivergence.H"
#include "SpaceUtils.H.multidim"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "PhaseGeom.H"
#include "MomentOp.H"
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#undef OLD_DIVJPERP
#undef FULL_DIAMAG_CURRENT

#include "NamespaceHeader.H" 

const char* VorticityOp::pp_name = {"vorticity_op"};

VorticityOp::VorticityOp( const ParmParse&    a_pp,
                          const MagGeom&      a_geometry,
                          const double        a_larmor,
                          const int           a_verbosity )
   : m_verbosity(a_verbosity),
     m_geometry(a_geometry),
     m_larmor(a_larmor),
     m_my_pc_idx_e(-1),
     m_my_pc_idx_i(-1),
     m_sigma_div_e_coefs_set(false),
     m_reynolds_stress(false),
     m_ExB_current_model(true),
     m_include_pol_den_correction(false),
     m_include_pol_den_correction_precond(false),
     m_include_pol_den_correction_to_pe(false),
     m_second_order(false),
     m_low_pollution(false),
     m_remove_axisymmetric_phi(false),
     m_update_pc_freq_oldmodel(1),
     m_it_counter(0),
     m_electron_temperature_func(NULL),
     m_update_pc_freq_e(1),
     m_update_pc_freq_i(1),
     m_update_pc_skip_stage_e(false),
     m_update_pc_skip_stage_i(false)
{
   const std::string name("potential");
   const std::string prefix( "BC." + name );
   ParmParse ppsp( prefix.c_str() );
   EllipticOpBCFactory elliptic_op_bc_factory;

   // This might be needed for seetin up fluid preconditioner  
   m_opt_string = "gksystem";

   bool old_model = (string(a_pp.prefix()) == "vorticity_op");
   
   // Input the initial conditions
   parseInitialConditions("", "potential");

   // Parse remaining parameters 
   parseParameters(a_pp);
   if (m_verbosity) {
     printParameters();
   }

   // Set up parallel current divergence operator
   m_parallel_current_divergence_op = new GKVorticity(a_pp, a_pp, a_geometry, a_larmor, 0., m_second_order,
                                                      m_low_pollution, m_include_pol_den_correction_precond);
   m_parallel_current_divergence_op->setModel("ParallelCurrent");
   m_parallel_current_divergence_op_bcs = elliptic_op_bc_factory.create( name,
                                                                         ppsp,
                                                                         *(a_geometry.getCoordSys()),
                                                                         false );

   // Set up parallel conductivity divergence operator
   m_par_cond_op = new GKVorticity(a_pp, a_pp, a_geometry, 0., 0., m_second_order, m_low_pollution,
                                   m_include_pol_den_correction_precond);
   m_par_cond_op_bcs = elliptic_op_bc_factory.create(name, ppsp, *(a_geometry.getCoordSys()), false);

   // Set up ImEx preconditioner operator
   if ( old_model ) {
      m_imex_pc_op = new GKVorticity(a_pp, a_pp, a_geometry, a_larmor, 0., m_second_order, m_low_pollution,
                                     m_include_pol_den_correction_precond);
   }
   else {
      string s_pp_base = string(a_pp.prefix()) + ".pc_imex";
      ParmParse pp_base( s_pp_base.c_str() );

      m_imex_pc_op = new GKVorticity(a_pp, pp_base, a_geometry, a_larmor, 0., m_second_order, m_low_pollution,
                                     m_include_pol_den_correction_precond);
   }
   m_imex_pc_op_bcs = elliptic_op_bc_factory.create(name, ppsp, *(a_geometry.getCoordSys()), false);

   // Set up gyrokinetic Poisson operator
   if ( old_model ) {
      m_gyropoisson_op = new GKPoisson(a_pp, a_geometry, a_larmor, 0., true, m_second_order, m_low_pollution);
   }
   else {
      string s_pp_base = string(a_pp.prefix()) + ".soln_op";
      ParmParse pp_base( s_pp_base.c_str() );

      m_gyropoisson_op = new GKPoisson(pp_base, a_geometry, a_larmor, 0., true, m_second_order, m_low_pollution);
   }
   m_gyropoisson_op_bcs = elliptic_op_bc_factory.create(name, ppsp, *(a_geometry.getCoordSys()), false);

   ParmParse pp_vlasov("vorticity_vlasov");
   m_vlasov = new PS::GKVlasov(pp_vlasov, a_larmor);

   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   m_negativeDivJpar.define(grids, 1, IntVect::Zero);

   m_ion_charge_density.define(grids, 1, 2*IntVect::Unit);
   m_ion_mass_density.define(grids, 1, IntVect::Zero);

   m_volume.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume);

   m_divJperp_mag_i.define(grids, 1, IntVect::Zero);
   m_divJperp_mag_e.define(grids, 1, IntVect::Zero);
   
   // Assign electron temperature
   m_electron_temperature.define(grids, 1, 2*IntVect::Unit);
   if (m_electron_temperature_func == NULL) {
     for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
       m_electron_temperature[dit].setVal(1.0);
     }
   }
   else {
     m_electron_temperature_func->assign( m_electron_temperature, m_geometry, 0.0);     
   }


}


VorticityOp::~VorticityOp()
{
   delete m_vlasov;
   delete m_gyropoisson_op;
   delete m_imex_pc_op;
   delete m_par_cond_op;
   delete m_parallel_current_divergence_op;
   delete m_parallel_current_divergence_op_bcs;
   delete m_par_cond_op_bcs;
   delete m_imex_pc_op_bcs;
   delete m_gyropoisson_op_bcs;
}


void VorticityOp::accumulateExplicitRHS( FluidSpeciesPtrVect&               a_rhs,
                                         const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                         const FluidSpeciesPtrVect&         a_fluid_species,
                                         const PS::ScalarPtrVect&           a_scalars,
                                         const EField&                      a_E_field,
                                         const int                          a_fluid_vec_comp,
                                         const Real                         a_time )
{
}


void VorticityOp::accumulateImplicitRHS( FluidSpeciesPtrVect&               a_rhs,
                                         const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                         const FluidSpeciesPtrVect&         a_fluid_species,
                                         const PS::ScalarPtrVect&           a_scalars,
                                         const EField&                      a_E_field,
                                         const int                          a_fluid_vec_comp,
                                         const Real                         a_time )
{
   CH_TIMERS("VorticityOp::accumulateImplicitRHS");
   CH_TIMER("sigma_div_e",t_sigma_div_e);
   CH_TIMER("div_jpar",t_div_jpar);
   CH_TIMER("div_jperp",t_div_jperp);

   CFGVars& rhs_fluid_species = *(a_rhs[a_fluid_vec_comp]);
   LevelData<FArrayBox>& rhs_data = rhs_fluid_species.cell_var(0);

   const CFGVars& sol_fluid_species = *(a_fluid_species[a_fluid_vec_comp]);
   const LevelData<FArrayBox>& sol_data = sol_fluid_species.cell_var(0);

   const DisjointBoxLayout& grids = m_geometry.gridsFull();

   CH_START(t_sigma_div_e);

   LevelData<FArrayBox> zero(grids, 1, IntVect::Zero);
   setZero(zero);
      
   if (a_scalars.size() > 0) {
      const Vector<Real>& scalar_data = a_scalars[0]->data();
      setCoreBC( scalar_data[0], -scalar_data[1], *m_par_cond_op_bcs );
   }

   if ( m_sigma_div_e_coefs_set ) {
      m_par_cond_op->updateBoundaries(*m_par_cond_op_bcs);
   }
   else {
      m_par_cond_op->setOperatorCoefficients( zero, *m_par_cond_op_bcs, false );
      m_sigma_div_e_coefs_set = true;
   }

   LevelData<FArrayBox> sigma_E_div_par(grids, 1, IntVect::Zero);
   m_par_cond_op->computeFluxDivergence(sol_data, sigma_E_div_par, false, false);

   CH_STOP(t_sigma_div_e);

   for (DataIterator dit(rhs_data.dataIterator()); dit.ok(); ++dit) {
      rhs_data[dit] -= sigma_E_div_par[dit];
   }

   CH_START(t_div_jpar);

   m_parallel_current_divergence_op->computeFluxDivergence(m_ion_charge_density, m_negativeDivJpar, false, true);
      
   CH_STOP(t_div_jpar);

   CH_START(t_div_jperp);

   // Compute the divergence of the perpendicular ExB ion current
   LevelData<FArrayBox> divJperp_ExB(grids, 1, IntVect::Zero);
   computeDivPerpIonExBCurrentDensity(divJperp_ExB, a_E_field, a_kinetic_species_phys, m_ion_charge_density, a_time);

   CH_STOP(t_div_jperp);

   // Compute rhs_data that describe time evolution of the negative vorticity
   
   // Add divergence of total parallel and perpendicular ion currents
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      rhs_data[dit] += m_negativeDivJpar[dit];
      rhs_data[dit] -= m_divJperp_mag_i[dit];
   }

   // Add divergence of perpendicular electron current
   if (m_ExB_current_model) {

     CH_START(t_div_jperp);
     // Compute the divergence of the perpendicular ExB ion current
     LevelData<FArrayBox> divJperp_ExB(grids, 1, IntVect::Zero);
     computeDivPerpIonExBCurrentDensity(divJperp_ExB, a_E_field, a_kinetic_species_phys, m_ion_charge_density, a_time);
     CH_STOP(t_div_jperp);

     for (DataIterator dit(grids); dit.ok(); ++dit) {
       rhs_data[dit] -= divJperp_ExB[dit];
     }
   }

   else {
     for (DataIterator dit(grids); dit.ok(); ++dit) {
       rhs_data[dit] -= m_divJperp_mag_e[dit];
     }
   }

   // Add reynolds_stress term
   if (m_reynolds_stress) {

     LevelData<FArrayBox> negative_vorticity(grids, 1, IntVect::Zero);
     m_gyropoisson_op->computeFluxDivergence(sol_data, negative_vorticity, false);

     LevelData<FArrayBox> negative_div_ExB_vorticity(grids, 1, IntVect::Zero);
     computeReynoldsStressTerm(negative_div_ExB_vorticity, negative_vorticity, a_E_field);
     for (DataIterator dit(grids); dit.ok(); ++dit) {
       rhs_data[dit] -= negative_div_ExB_vorticity[dit];
     }
   }

   // Add polarization density correction into electron pressure
   if (m_include_pol_den_correction) {
      
      // Compute negative vorticity
      LevelData<FArrayBox> negative_vorticity(grids, 1, IntVect::Zero);
      m_gyropoisson_op->computeFluxDivergence(sol_data, negative_vorticity, false);
      
      // Compute par_grad * sigma * par_grad of negative_vorticity
      LevelData<FArrayBox> pol_dens_correction(grids, 1, IntVect::Zero);
      //m_par_cond_op->computeFluxDivergence(negative_vorticity, pol_dens_correction, false, true);
      m_parallel_current_divergence_op->computeFluxDivergence(negative_vorticity, pol_dens_correction, false, true);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs_data[dit] -= pol_dens_correction[dit];
      }
   }
}


void VorticityOp::preSolutionOpEval( const PS::KineticSpeciesPtrVect&   a_kinetic_species,
                                     const FluidSpeciesPtrVect&         a_fluid_species,
                                     const PS::ScalarPtrVect&           a_scalars,
                                     const EField&                      a_E_field,
                                     const double                       a_time )
{
   CH_TIME("VorticityOp::preSolutionOpEval");

   computeIonMassDensity( m_ion_mass_density, a_kinetic_species );

   m_gyropoisson_op->setOperatorCoefficients( m_ion_mass_density, *m_gyropoisson_op_bcs, true );
}


void VorticityOp::preOpEval( const PS::KineticSpeciesPtrVect&   a_kinetic_species,
                             const FluidSpeciesPtrVect&         a_fluid_species,
                             const PS::ScalarPtrVect&           a_scalars,
                             const EField&                      a_E_field,
                             const double                       a_time )
{
   CH_TIME("VorticityOp::preOpEval");

   preSolutionOpEval(a_kinetic_species, a_fluid_species, a_scalars, a_E_field, a_time);

   // Make sure there are enough ghost cells to compute the charge density in the ghost cells
   // of m_ion_charge_density
   for (int species(0); species<a_kinetic_species.size(); species++) {
      PS::LevelData<PS::FArrayBox>& dfn = a_kinetic_species[species]->distributionFunction();
      const PS::PhaseGeom& phase_geom = a_kinetic_species[species]->phaseSpaceGeometry();
      CH_assert(phase_geom.config_restrict(dfn.ghostVect()) >= m_ion_charge_density.ghostVect());
   }

   computeIonChargeDensity( m_ion_charge_density, a_kinetic_species );

   m_parallel_current_divergence_op->setOperatorCoefficients( m_ion_charge_density,
                                                              *m_parallel_current_divergence_op_bcs,
                                                              false );

   computeDivPerpIonMagCurrentDensity(m_divJperp_mag_i, a_E_field, a_kinetic_species, a_time);

   if (!m_include_pol_den_correction_to_pe) {
     const DisjointBoxLayout& grids = m_ion_charge_density.disjointBoxLayout();
     LevelData<FArrayBox> zero(grids, 1, IntVect::Zero);
     setZero(zero);
     computeDivPerpElectronMagCurrentDensity(m_divJperp_mag_e, m_ion_charge_density, zero, a_time);
   }
   else {
     MayDay::Error("VorticityOp::preOpEval: polarization density corrections to electron pressure are not yet implemented for the new vorticity model ");
   }
}


void VorticityOp::evalSolutionOp( FluidSpeciesPtrVect&               a_rhs,
                                  const PS::KineticSpeciesPtrVect&   a_kinetic_species_comp,
                                  const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                  const FluidSpeciesPtrVect&         a_fluid_species_comp,
                                  const FluidSpeciesPtrVect&         a_fluid_species_phys,
                                  const PS::ScalarPtrVect&           a_scalars,
                                  const int                          a_component,
                                  const Real                         a_time )
{
   CH_TIME("VorticityOp::evalSolutionOp");

   if (a_scalars.size() > 0) {
      const Vector<Real>& scalar_data = a_scalars[0]->data();
      setCoreBC( scalar_data[0], -scalar_data[1], *m_gyropoisson_op_bcs );
      //      if (procID()==0) cout << "Er_lo = " << scalar_data[0] << ", Er_hi = " << scalar_data[1] << endl;
   }

   m_gyropoisson_op->updateBoundaries(*m_gyropoisson_op_bcs);

   const FluidSpecies& sol_species = static_cast<const FluidSpecies&>(*a_fluid_species_phys[a_component]);
   FluidSpecies& rhs_species = static_cast<FluidSpecies&>(*a_rhs[a_component]);

   const LevelData<FArrayBox>& sol = sol_species.cell_var(0);
   LevelData<FArrayBox>& rhs = rhs_species.cell_var(0);

   m_gyropoisson_op->computeFluxDivergence(sol, rhs, false, false);
}

void VorticityOp::defineBlockPC(  std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>&  a_pc,
                                  std::vector<PS::DOFList>&                                     a_dof_list,
                                  const PS::ODEVector&                                          a_soln_vec,
                                  void*                                                         a_gkops,
                                  const std::string&                                            a_out_string,
                                  const std::string&                                            a_opt_string,
                                  bool                                                          a_im,
                                  const FluidSpecies&                                           a_fluid_species,
                                  const PS::GlobalDOFFluidSpecies&                              a_global_dofs,
                                  const int                                                     a_species_idx,
                                  const int                                                     a_id )
{
  
  CH_assert(a_pc.size() == a_dof_list.size());

  if (!procID()) {
    if (a_im) {
      std::cout << "  Fluid Species " << a_species_idx
                << " : "
                << " creating " << _FLUID_OP_PC_ << " preconditioner RHSOp"
                << " (index = " << a_pc.size() << ").\n";
    } else {
      std::cout << "  Fluid Species " << a_species_idx
                << " : "
                << " creating " << _FLUID_OP_PC_ << " preconditioner for solutionOp"
                << " (index = " << a_pc.size() << ").\n";
    }
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

  if (a_im)  m_my_pc_idx_i = a_pc.size();
  else       m_my_pc_idx_e = a_pc.size();

  a_pc.push_back(pc);
  a_dof_list.push_back(dof_list);

  return;
}

void VorticityOp::updateBlockPC(  std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>&  a_pc,
                                  const PS::KineticSpeciesPtrVect&                              a_kin_species_phys,
                                  const FluidSpeciesPtrVect&                                    a_fluid_species,
                                  const Real                                                    a_time,
                                  const int                                                     a_step,
                                  const int                                                     a_stage,
                                  const Real                                                    a_shift,
                                  const bool                                                    a_im,
                                  const int                                                     a_species_idx )
{
  if (a_im) {
    CH_assert(m_my_pc_idx_i >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx_i);
    /* see if this update can be skipped */
    if (a_step >= 0) {
      if (a_step%m_update_pc_freq_i != 0) return;
      if ((a_stage >= 0) && m_update_pc_skip_stage_i) return; 
      if ((a_stage < 0) && (!m_update_pc_skip_stage_i)) return;
    }
  } else {
    CH_assert(m_my_pc_idx_e >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx_e);
    /* see if this update can be skipped */
    if (a_step >= 0) {
      if (a_step%m_update_pc_freq_e != 0) return;
      if ((a_stage >= 0) && m_update_pc_skip_stage_e) return; 
      if ((a_stage < 0) && (!m_update_pc_skip_stage_e)) return;
    }
  }

  if (!procID()) {
    std::cout << "    ==> Updating " << _FLUID_OP_PC_ << " preconditioner " 
              << " for VorticityOp " << (a_im ? "RHS " : "LHS ")
              << "of fluid species " << a_species_idx 
              << ".\n";
  }

  PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt> *pc 
    = dynamic_cast<PS::FluidOpPreconditioner<PS::ODEVector,PS::AppCtxt>*>
      (a_pc[(a_im ? m_my_pc_idx_i : m_my_pc_idx_e)]);
  CH_assert(pc != NULL);
  pc->update( a_kin_species_phys, 
              a_fluid_species, 
              a_time, 
              a_step,
              a_stage,
              a_shift, 
              a_im, 
              a_species_idx );

  return;
}

void VorticityOp::solveSolutionPC( FluidSpeciesPtrVect&              a_fluid_species_solution,
                                   const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                   const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                   const int                         a_component )
{
   CH_TIME("VorticityOp::solveSolutionPC");
   FluidSpecies& sol_species = static_cast<FluidSpecies&>(*a_fluid_species_solution[a_component]);
   const FluidSpecies& rhs_species = static_cast<const FluidSpecies&>(*a_fluid_species_rhs[a_component]);

   const LevelData<FArrayBox>& r = rhs_species.cell_var(0);
   LevelData<FArrayBox>& z = sol_species.cell_var(0);

   setZero(z);
   m_gyropoisson_op->solvePreconditioner(r, z);
}

void VorticityOp::updatePCImEx( const FluidSpeciesPtrVect&       a_fluid_species,
                                const PS::KineticSpeciesPtrVect& a_kinetic_species,
                                const double                     a_time,
                                const int                        a_step,
                                const int                        a_stage,
                                const double                     a_shift,
                                const int                        a_component)
{
   CH_TIME("VorticityOp::updatePCImEx");

   m_imex_pc_op->m_alpha = a_shift;

   m_imex_pc_op->setVorticityOperatorCoefficients(m_ion_mass_density, m_ion_charge_density,
                                                  *m_imex_pc_op_bcs, true);
}


void VorticityOp::solvePCImEx( FluidSpeciesPtrVect&              a_fluid_species_solution,
                               const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                               const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                               const int                         a_component )
{
   CH_TIME("VorticityOp::solvePCImEx");

   const FluidSpecies& rhs_species = static_cast<const FluidSpecies&>(*a_fluid_species_rhs[a_component]);
   const LevelData<FArrayBox>& r = rhs_species.cell_var(0);

   FluidSpecies& sol_species = static_cast<FluidSpecies&>(*a_fluid_species_solution[a_component]);
   LevelData<FArrayBox>& z = sol_species.cell_var(0);

   setZero(z);
   m_imex_pc_op->solvePreconditioner(r, z);
}


void VorticityOp::fillGhostCells( FluidSpecies&  a_species_phys,
                                  const double   a_time )
{
}


void VorticityOp::computeDivPerpIonMagCurrentDensity( LevelData<FArrayBox>&             a_divJperp_mag,
                                                      const EField&                     a_E_field,
                                                      const PS::KineticSpeciesPtrVect&  a_species_phys,
                                                      const Real&                       a_time )
{
   CH_TIME("VorticityOp::computeDivPerpIonMagCurrentDensity");

   PS::MomentOp& moment_op = PS::MomentOp::instance();

   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   LevelData<FArrayBox> species_div_Jperp(grids, 1, IntVect::Zero);

   setZero(a_divJperp_mag);
      
   for (int species(0); species<a_species_phys.size(); species++) {
      
      const PS::KineticSpecies& this_species( *(a_species_phys[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      const RefCountedPtr<PS::PhaseGeom>& phase_geometry( this_species.phaseSpaceGeometryPtr() );
      PS::KineticSpecies tmp_rhs_species("tmp_rhs", this_species.mass(), this_species.charge(), phase_geometry);
      PS::LevelData<PS::FArrayBox>& rhs_dfn = tmp_rhs_species.distributionFunction();
      rhs_dfn.define( this_species.distributionFunction().getBoxes(), 1, PS::IntVect::Zero);
      
      // Compute the divergence of -Jperp due to magnetic drifts for this species
      
      //Gyroaverage calculation requires phi. Pass zero phi for now. Fix later
      LevelData<FArrayBox> phi_tmp(grids, 1, IntVect::Zero);
      setZero(phi_tmp);
      m_vlasov->evalRHS(tmp_rhs_species, this_species, phi_tmp, a_E_field,
                        PS::PhaseGeom::MAGNETIC_DRIFT_VELOCITY, a_time);

      // Divide by J to get physical cell averages
      phase_geometry->divideJonValid( rhs_dfn );
      
      moment_op.compute( species_div_Jperp, tmp_rhs_species, PS::ChargeDensityKernel<PS::FArrayBox>() );
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_divJperp_mag[dit] -= species_div_Jperp[dit];
      }
   }
}

void VorticityOp::computeDivPerpIonExBCurrentDensity( LevelData<FArrayBox>&             a_div_Jperp,
                                                      const EField&                     a_E_field,
                                                      const PS::KineticSpeciesPtrVect&  a_species_phys,
                                                      const LevelData<FArrayBox>&       a_ion_charge_density,
                                                      const Real&                       a_time)
                                                  
{
   CH_TIME("VorticityOp::computeDivPerpIonExBCurrentDensity");
   
   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   
   // Compute the J_ExB divergence
   
#ifndef OLD_DIVJPERP
   
   // Make sure there are enough ghost cells to compute the face average density below
   for (int species(0); species<a_species_phys.size(); species++) {
      PS::LevelData<PS::FArrayBox>& dfn = a_species_phys[species]->distributionFunction();
      CH_assert(dfn.ghostVect() >= 2*PS::IntVect::Unit);
   }

   // Get GKVelocity option (presently, use the first species)
   const PS::KineticSpecies& species0( *(a_species_phys[0]) );
   const RefCountedPtr<PS::PhaseGeom>& phase_geometry( species0.phaseSpaceGeometryPtr() );
   bool divFreeVelocity = phase_geometry->divFreeVelocity();
   bool fourthOrder = !m_geometry.secondOrder();

   // Compute the current due to ExB drifts
   if (divFreeVelocity) {

      LevelData<FluxBox> J_ExB(grids, 1, IntVect::Unit);
      m_geometry.computeBxEIntegrals(a_E_field.getPhiNode(), true, J_ExB);
      
      LevelData<FluxBox> charge_density_face(grids, 1, IntVect::Zero);
      fourthOrderCellToFace(charge_density_face, a_ion_charge_density);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0, sign=-1; dir<SpaceDim; ++dir, sign*=-1) {
            J_ExB[dit][dir] *= sign * m_larmor;
            J_ExB[dit][dir] *= charge_density_face[dit][dir];
         }
      }
      
      RealVect fakeDx = RealVect::Unit;
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         simpleDivergence(a_div_Jperp[dit], J_ExB[dit], grids[dit], fakeDx);
      }
   }
   
   else {
      
      // Compute the ExB drift on face centers
      LevelData<FluxBox> ExB_drift(grids, 3, IntVect::Unit);
      m_geometry.computeEXBDrift(a_E_field.getFaceCenteredField(), ExB_drift);
      
      LevelData<FluxBox> charge_density_face(grids, 1, IntVect::Unit);
      fourthOrderCellToFace(charge_density_face, a_ion_charge_density);
      
      LevelData<FluxBox> J_ExB(grids, SpaceDim, IntVect::Unit);
      for (DataIterator dit(J_ExB.dataIterator()); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; dir++) {
            for (int nComp=0; nComp<SpaceDim; nComp++) {
               J_ExB[dit][dir].copy(charge_density_face[dit][dir],0,nComp);
            }
            J_ExB[dit][dir].mult(ExB_drift[dit][dir],0,0);
            if (SpaceDim == 2) {
               J_ExB[dit][dir].mult(ExB_drift[dit][dir],2,1);
            }
            if (SpaceDim == 3) {
               J_ExB[dit][dir].mult(ExB_drift[dit][dir],1,1);
               J_ExB[dit][dir].mult(ExB_drift[dit][dir],2,2);
            }
	    J_ExB[dit][dir] *= m_larmor;
         }
      }

      m_geometry.applyAxisymmetricCorrection(J_ExB);
      m_geometry.computeMappedGridDivergence(J_ExB, a_div_Jperp, fourthOrder);
      
   }
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_div_Jperp[dit] /= m_volume[dit];
   }
   
#else
   
   PS::MomentOp& moment_op = PS::MomentOp::instance();
   
   LevelData<FArrayBox> species_div_Jperp(grids, 1, IntVect::Zero);
   
   setZero(a_div_Jperp);
   
   for (int species(0); species<a_species_phys.size(); species++) {
      
      const PS::KineticSpecies& this_species( *(a_species_phys[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      const RefCountedPtr<PS::PhaseGeom>& phase_geometry( this_species.phaseSpaceGeometry() );
      PS::KineticSpecies tmp_rhs_species("tmp_rhs", this_species.mass(), this_species.charge(), phase_geometry);
      PS::LevelData<PS::FArrayBox>& rhs_dfn = tmp_rhs_species.distributionFunction();
      rhs_dfn.define( this_species.distributionFunction().getBoxes(), 1, PS::IntVect::Zero);
      
      // Compute the divergence of -Jperp due to ExB drifts for this species
      bool fourth_order = !m_parallel_current_divergence_op->secondOrder();
      
      //Gyroaverage calculation requires phi. Pass zero phi for now. Fix later
      LevelData<FArrayBox> phi_tmp(grids, 1, IntVect::Zero);
      setZero(phi_tmp);
      m_vlasov->evalRHS(tmp_rhs_species, this_species, phi_tmp, a_E_field,
                        PS::PhaseGeom::MAGNETIC_DRIFT_VELOCITY, a_time);
     
      // Divide by J to get physical cell averages      
      phase_geometry->divideJonValid( rhs_dfn );
      
      moment_op.compute( species_div_Jperp, tmp_rhs_species, PS::ChargeDensityKernel<PS::FArrayBox>() );
      
      for (DataIterator dit(a_div_Jperp.dataIterator()); dit.ok(); ++dit) {
         a_div_Jperp[dit] -= species_div_Jperp[dit];
      }
   }
   
#endif
}


void VorticityOp::computeDivPerpElectronMagCurrentDensity( LevelData<FArrayBox>&             a_div_Jperp,
							   const LevelData<FArrayBox>&       a_ion_charge_density,
							   const LevelData<FArrayBox>&       a_negative_vorticity,
							   const Real&                       a_time)

{
   CH_TIME("VorticityOp::computeDivPerpElectronMagCurrentDensity");

   const DisjointBoxLayout& grids = m_geometry.gridsFull();

   // Get magnetic field data
   const LevelData<FluxBox>& gradB = m_geometry.getFCGradBFieldMag();
   const LevelData<FluxBox>& curlb = m_geometry.getFCCurlBFieldDir();
   const LevelData<FArrayBox>& Bfield = m_geometry.getCCBFieldMag();

   // Compute advection velocity
   LevelData<FluxBox> adv_vel(grids, 3, IntVect::Unit);
   m_geometry.computeEXBDrift(gradB, adv_vel);
   for (DataIterator dit(adv_vel.dataIterator()); dit.ok(); ++dit) {
     for (int dir=0; dir<SpaceDim; dir++) {
       adv_vel[dit][dir].negate();
       adv_vel[dit][dir].plus(curlb[dit][dir]);
     }
   }

#ifdef FULL_DIAMAG_CURRENT
   // Compute -grad(n) [the calculation below extrapolates n in ghosts]
   LevelData<FArrayBox> ne(grids, 1, IntVect::Zero);
   for (DataIterator dit(ne.dataIterator()); dit.ok(); ++dit) {
     ne[dit].copy(a_ion_charge_density[dit]);
     if (m_include_pol_den_correction_to_pe) {
       ne[dit].minus(a_negative_vorticity[dit]);
     }
   }
   LevelData<FluxBox> neg_grad_n(grids, 3, IntVect::Unit);
   m_gyropoisson_op->computeField(ne, neg_grad_n);

   
   // Compute -grad(n)xB/B^2
   m_geometry.computeEXBDrift(neg_grad_n, adv_vel);
#endif

   // Compute scalar advective quantity (u)
   LevelData<FArrayBox> u(grids, 1, 2*IntVect::Unit);
   for (DataIterator dit(u.dataIterator()); dit.ok(); ++dit) {
     u[dit].copy(a_ion_charge_density[dit]);

     if (m_include_pol_den_correction_to_pe) {
       u[dit].minus(a_negative_vorticity[dit]);
     }

     u[dit].mult(m_electron_temperature[dit]);
     u[dit].divide(Bfield[dit]);
     fourthOrderCellExtrap(u[dit],grids[dit]);
     
#ifdef FULL_DIAMAG_CURRENT
     u[dit].setVal(1.0);
#endif     
   }

   m_geometry.fillInternalGhosts(u);

   // Convert cell-averaged vorticity to face-averaged quantity
   // presently, use centered scheme, update to high-order UW later
   LevelData<FluxBox> u_fc(grids, 1, IntVect::Unit);
   fourthOrderCellToFace(u_fc, u);

   // Compute v times u flux
   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);
   for (DataIterator dit(flux.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         for (int nComp=0; nComp<SpaceDim; nComp++) {
            flux[dit][dir].copy(u_fc[dit][dir],0,nComp);
         }
         flux[dit][dir].mult(adv_vel[dit][dir],0,0);
         if (SpaceDim == 2) {
            flux[dit][dir].mult(adv_vel[dit][dir],2,1);
         }
         if (SpaceDim == 3) {
            flux[dit][dir].mult(adv_vel[dit][dir],1,1);
            flux[dit][dir].mult(adv_vel[dit][dir],2,2);
         }
         flux[dit][dir] *= m_larmor;
      }
   }

   m_geometry.applyAxisymmetricCorrection(flux);

   bool fourthOrder = !m_geometry.secondOrder();

   m_geometry.computeMappedGridDivergence(flux, a_div_Jperp, fourthOrder);
      
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_div_Jperp[dit] /= m_volume[dit];
   }   
}


void VorticityOp::computeReynoldsStressTerm( LevelData<FArrayBox>&             a_div_ExB_times_vorticity,
                                             const LevelData<FArrayBox>&       a_vorticity,
                                             const EField&                     a_E_field) const
{
   CH_TIME("VorticityOp::computeReynoldsStressTerm");
   
   // The present implementation supports only the second-order option
   // since vorticity_grown has only 2 layers of ghosts and
   // we use fourthOrderCellToFace. Extend to 4th order later.
   CH_assert(m_geometry.secondOrder());
   
   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   
   // Compute the ExB drift on face centers
   LevelData<FluxBox> ExB_drift(grids, 3, IntVect::Unit);
   m_geometry.computeEXBDrift(a_E_field.getFaceCenteredField(), ExB_drift);

   // Fill ghost cells (currently use 4th order extrapolation)
   LevelData<FArrayBox> vorticity_grown(grids, 1, 2*IntVect::Unit);
   
   for (DataIterator dit(a_vorticity.dataIterator()); dit.ok(); ++dit) {
      vorticity_grown[dit].copy(a_vorticity[dit]);
      fourthOrderCellExtrap(vorticity_grown[dit],grids[dit]);
   }
  
   m_imex_pc_op->fillInternalGhosts(vorticity_grown);

   // Convert cell-averaged vorticity to face-averaged quantity
   // presently, use centered scheme, update to high-order UW later
   LevelData<FluxBox> faceVorticity(grids, 1, IntVect::Unit);
   fourthOrderCellToFace(faceVorticity, vorticity_grown);

#if 1
   // Experimental option to use uw3 method for ExB vorticity advection
   // Might not be useful, since cfl is likely dominated by stiff k^4 term
   LevelData<FluxBox> ExB_drift_mapped(grids, SpaceDim, IntVect::Unit);
   for (DataIterator dit(ExB_drift_mapped.dataIterator()); dit.ok(); ++dit) {
      if (SpaceDim == 3) {
         ExB_drift_mapped[dit].copy(ExB_drift[dit]);
      }
      if (SpaceDim == 2) {
         ExB_drift_mapped[dit].copy(ExB_drift[dit],0,0,1);
         ExB_drift_mapped[dit].copy(ExB_drift[dit],2,1,1);
      }
   }
   m_geometry.multNTransposePointwise(ExB_drift_mapped);
   
   SpaceUtils::upWindToFaces(faceVorticity,
                             vorticity_grown,
                             ExB_drift_mapped,
                             "uw3");
#endif
   
   // Compute ExB times vorticity flux
   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);
   for (DataIterator dit(a_vorticity.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         for (int nComp=0; nComp<SpaceDim; nComp++) {
            flux[dit][dir].copy(faceVorticity[dit][dir],0,nComp);
         }
         flux[dit][dir].mult(ExB_drift[dit][dir],0,0);
         if (SpaceDim == 2) {
            flux[dit][dir].mult(ExB_drift[dit][dir],2,1);
         }
         if (SpaceDim == 3) {
            flux[dit][dir].mult(ExB_drift[dit][dir],1,1);
            flux[dit][dir].mult(ExB_drift[dit][dir],2,2);
         }
         flux[dit][dir] *= m_larmor;
      }
   }

   m_geometry.applyAxisymmetricCorrection(flux);

   bool fourthOrder = !m_geometry.secondOrder();
   m_geometry.computeMappedGridDivergence(flux, a_div_ExB_times_vorticity, fourthOrder);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_div_ExB_times_vorticity[dit] /= m_volume[dit];
   }
}


void VorticityOp::computeIonMassDensity( LevelData<FArrayBox>&             a_mass_density,
                                         const PS::KineticSpeciesPtrVect&  a_species ) const
{
   CH_TIME("VorticityOp::computeIonMassDensity");
      
   // Container for individual species charge density
   LevelData<FArrayBox> species_mass_density;
   species_mass_density.define(a_mass_density);

   setZero( a_mass_density );

   for (int species(0); species<a_species.size(); species++) {
         
      const PS::KineticSpecies& this_species( *(a_species[species]) );
      if ( this_species.charge() < 0.0 ) continue;

      // Compute the mass density for this species including ghost cells
      this_species.massDensity( species_mass_density );
      
      for (DataIterator dit(a_mass_density.dataIterator()); dit.ok(); ++dit) {
         a_mass_density[dit] += species_mass_density[dit];
      }
   }
}


void VorticityOp::computeIonChargeDensity( LevelData<FArrayBox>&               a_ion_charge_density,
                                           const PS::KineticSpeciesPtrVect&    a_species ) const
{
   CH_TIME("VorticityOp::computeIonChargeDensity");
   
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


void VorticityOp::setCoreBC( const double   a_core_inner_bv,
                             const double   a_core_outer_bv,
                             EllipticOpBC&  a_bc ) const 
{
   if ( typeid(*(m_geometry.getCoordSys())) == typeid(LogRectCoordSys) ) {
     if (m_consistent_lower_bc_only) { 
       a_bc.setBCValue(0,RADIAL_DIR,0,a_core_inner_bv);
     }
     else if (m_consistent_upper_bc_only) {
       a_bc.setBCValue(0,RADIAL_DIR,1,a_core_outer_bv);
     }
     else{
       a_bc.setBCValue(0,RADIAL_DIR,0,a_core_inner_bv);
       a_bc.setBCValue(0,RADIAL_DIR,1,a_core_outer_bv);
     }
   }
   else if ( typeid(*(m_geometry.getCoordSys())) == typeid(SNCoreCoordSys) ) {
      a_bc.setBCValue(SNCoreBlockCoordSys::LCORE,RADIAL_DIR,0,a_core_inner_bv);
      a_bc.setBCValue(SNCoreBlockCoordSys::LCORE,RADIAL_DIR,1,a_core_outer_bv);
   }
   else if ( typeid(*(m_geometry.getCoordSys())) == typeid(SingleNullCoordSys) ) {
      a_bc.setBCValue(SingleNullBlockCoordSys::LCORE,RADIAL_DIR,0,a_core_inner_bv);
   }
   else {
      MayDay::Error("EField::setCoreBC(): unknown geometry with consistent bcs");
   }
}


inline
void VorticityOp::setZero( LevelData<FArrayBox>& a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit].setVal(0.);
   }
}


void VorticityOp::initializeOldModel( const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                      const Vector<Real>&               a_scalar_data )
{
   computeIonMassDensity(m_ion_mass_density, a_kinetic_species);

   if (a_scalar_data.size() > 0) {
      setCoreBC( a_scalar_data[0], -a_scalar_data[1], *m_gyropoisson_op_bcs );
   }
   m_gyropoisson_op->setOperatorCoefficients(m_ion_mass_density, *m_gyropoisson_op_bcs, false);
}


void VorticityOp::updatePotentialOldModel( LevelData<FArrayBox>&             a_phi,
                                           EField&                           a_E_field,
                                           const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                           const Vector<Real>&               a_scalar_data,
                                           const Real                        a_dt,
                                           const Real                        a_time )
{
   CH_TIME("VorticityOp::updatePotentialOldModel");
   const DisjointBoxLayout& grids = m_geometry.gridsFull();

   // Compute the polarization density term (at the old time) using the operator and boundary
   // conditions set in the previous call of this function or initializeOldModel().
   // NB: gkPoissonRHS = -vorticity
   LevelData<FArrayBox> gkPoissonRHS(grids, 1, IntVect::Zero);
   m_gyropoisson_op->computeFluxDivergence(a_phi, gkPoissonRHS, false);
   
   // Compute the Reynolds stress term (at the old time)
   // NB: -div(ExB*vorticity) = div(ExB*gkPoissonRHS)
   LevelData<FArrayBox> negative_div_ExB_vorticity(grids, 1, IntVect::Zero);
   if (m_reynolds_stress) {
     computeReynoldsStressTerm(negative_div_ExB_vorticity, gkPoissonRHS, a_E_field);
   }
   
   // Compute the parallel current divergence
   computeIonChargeDensity(m_ion_charge_density, a_kinetic_species);
   m_parallel_current_divergence_op->setImplicitDt(a_dt);
   m_parallel_current_divergence_op->setOperatorCoefficients(m_ion_charge_density, *m_parallel_current_divergence_op_bcs, false);
   m_parallel_current_divergence_op->computeFluxDivergence(m_ion_charge_density, m_negativeDivJpar, false, true);
    
   // Compute the perpendicular current divergence
   LevelData<FArrayBox> div_J_electron_mag(grids, 1, IntVect::Zero);
   if (m_ExB_current_model) {
     computeDivPerpIonExBCurrentDensity(div_J_electron_mag, a_E_field, a_kinetic_species, m_ion_charge_density, a_time);
   }
   else{
     computeDivPerpElectronMagCurrentDensity(div_J_electron_mag, m_ion_charge_density, gkPoissonRHS, a_time);
   }
   
   LevelData<FArrayBox> div_J_ion_mag(grids, 1, IntVect::Zero);
   computeDivPerpIonMagCurrentDensity(div_J_ion_mag, a_E_field, a_kinetic_species, a_time);

   // Get the RHS for the implicit vorticity solve
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      gkPoissonRHS[dit] += m_negativeDivJpar[dit];
      
      div_J_electron_mag[dit] *= a_dt;
      gkPoissonRHS[dit] -= div_J_electron_mag[dit];

      div_J_ion_mag[dit] *= a_dt;
      gkPoissonRHS[dit] -= div_J_ion_mag[dit];

      if (m_reynolds_stress) {
         negative_div_ExB_vorticity[dit] *= a_dt;
         gkPoissonRHS[dit] -= negative_div_ExB_vorticity[dit];
      }
   }
      
   // Do the implicit vorticity solve
   computeIonMassDensity(m_ion_mass_density, a_kinetic_species);
   m_imex_pc_op->setImplicitDt(a_dt);
      
   if (a_scalar_data.size() > 0) setCoreBC( a_scalar_data[0], -a_scalar_data[1], *m_imex_pc_op_bcs );

   bool update_precond = (m_it_counter % m_update_pc_freq_oldmodel == 0) ? true : false;
   m_imex_pc_op->setVorticityOperatorCoefficients(m_ion_mass_density, m_ion_charge_density, *m_imex_pc_op_bcs,
                                                  update_precond);

   m_imex_pc_op->solveWithBCs( a_phi, gkPoissonRHS );
   m_imex_pc_op->fillInternalGhosts( a_phi );

   // Update nodal phi if supporting the calculation of a divergence-free phase velocity
   a_E_field.interpToNodes(a_phi);

   // Set up the polarization density operator for the next step
   if (a_scalar_data.size() > 0) setCoreBC( a_scalar_data[0], -a_scalar_data[1], *m_gyropoisson_op_bcs );
   m_gyropoisson_op->setOperatorCoefficients(m_ion_mass_density, *m_gyropoisson_op_bcs, false);

   // Optional phi post-processing to remove the axisymmetric component
   if (m_remove_axisymmetric_phi) {
     FluxSurface flux_surface(m_geometry);
     LevelData<FArrayBox> phi_fs(grids, 1, IntVect::Zero);
     flux_surface.averageAndSpread(a_phi, phi_fs);

     for (DataIterator dit(grids); dit.ok(); ++dit) {
       a_phi[dit] -= phi_fs[dit];
     }
   }

   m_it_counter++;
}


void VorticityOp::parseParameters( const ParmParse& a_pp )
{
   a_pp.query( "include_reynolds_stress", m_reynolds_stress);
   a_pp.query( "ExB_current_model", m_ExB_current_model );
   a_pp.query( "include_pol_den_correction", m_include_pol_den_correction );

   m_include_pol_den_correction_precond = m_include_pol_den_correction; // default precond value to operator value
   a_pp.query( "include_pol_den_correction_precond", m_include_pol_den_correction_precond );

   a_pp.query( "include_pol_den_correction_to_pe", m_include_pol_den_correction_to_pe );
   a_pp.query( "second_order", m_second_order );
   a_pp.query( "low_pollution", m_low_pollution );
   a_pp.query( "remove_axisymmetric_phi", m_remove_axisymmetric_phi );

   if (a_pp.contains("update_precond_interval")) {
      a_pp.get( "update_precond_interval", m_update_pc_freq_oldmodel );
      m_update_pc_freq_e = m_update_pc_freq_oldmodel;
      m_update_pc_freq_i = m_update_pc_freq_oldmodel;
   } else {
      a_pp.query( "update_precond_interval_LHS", m_update_pc_freq_e);
      a_pp.query( "update_precond_interval_RHS", m_update_pc_freq_i);
   }
   if (a_pp.contains("update_precond_skip_stage")) {
      a_pp.get( "update_precond_skip_stage", m_update_pc_skip_stage_e );
      m_update_pc_skip_stage_i = m_update_pc_skip_stage_e;
   } else {
      a_pp.query( "update_precond_skip_stage_LHS", m_update_pc_skip_stage_e );
      a_pp.query( "update_precond_skip_stage_RHS", m_update_pc_skip_stage_i );
   }
   
   
   if (a_pp.contains("consistent_upper_bc_only")) {
     a_pp.get("consistent_upper_bc_only", m_consistent_upper_bc_only);
   }
   else {
     m_consistent_upper_bc_only = false;
   }

   if (a_pp.contains("consistent_lower_bc_only")) {
     a_pp.get("consistent_lower_bc_only",m_consistent_lower_bc_only);
   }
   else	{
     m_consistent_lower_bc_only = false;
   }
   
   if (a_pp.contains("electron_temperature")) {
      GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
      std::string grid_function_name;
      a_pp.get("electron_temperature", grid_function_name );
      m_electron_temperature_func = grid_library->find( grid_function_name );
   }
}


void VorticityOp::printParameters()
{
   if (procID()==0) {
      std::cout << "Vorticity model parameters:" << std::endl;
      std::cout << "  include_reynolds_stress  =  " << m_reynolds_stress
                << ", ExB_current_model = " << m_ExB_current_model
                << std::endl;

      if (m_electron_temperature_func != NULL) {
         std::cout << "  Electron Temperature Function:" << std::endl;
         m_electron_temperature_func->printParameters();
      }
   }
}

#include "NamespaceFooter.H"
