#include "VorticityOp.H"
#include "Directions.H"

#include "EllipticOpBCFactory.H"
#include "LogRectCoordSys.H"
#include "SNCoreCoordSys.H"
#include "FourthOrderUtil.H"
#include "SimpleDivergence.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "PhaseGeom.H"
#include "MomentOp.H"
#include "GKOps.H"
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#undef OLD_DIVJPERP

#include "NamespaceHeader.H" 

const char* VorticityOp::pp_name = {"vorticity_op"};

VorticityOp::VorticityOp( const std::string&  a_pp_str,
                          const MagGeom&      a_geometry,
                          const double        a_larmor,
                          const int           a_verbosity )
   : m_geometry(a_geometry),
     m_larmor(a_larmor),
     m_opt_string(a_pp_str),
     m_my_pc_idx_e(-1),
     m_my_pc_idx_i(-1),
     m_sigma_div_e_coefs_set(false)
{
   const std::string name("potential");
   const std::string prefix( "BC." + name );
   ParmParse ppsp( prefix.c_str() );
   EllipticOpBCFactory elliptic_op_bc_factory;
   ParmParse pp_vorticity_op("gkpoisson");

   m_parallel_current_divergence_op = new GKPoisson(pp_vorticity_op, a_geometry, a_larmor, 0.);
   m_parallel_current_divergence_op->m_model = "ParallelCurrent";
   m_parallel_current_divergence_op->m_dt_implicit = 1.;
   m_parallel_current_divergence_op_bcs = elliptic_op_bc_factory.create( name,
                                                                         ppsp,
                                                                         "SingleNull",
                                                                         false );
   setCoreBC( 0., 0., *m_parallel_current_divergence_op_bcs );

   m_par_cond_op = new GKPoisson(pp_vorticity_op, a_geometry, 0., 0.);
   m_par_cond_op->m_model = "Vorticity";
   m_par_cond_op->m_dt_implicit = 1.;
   m_par_cond_op_bcs = elliptic_op_bc_factory.create( name,
                                                      ppsp,
                                                      "SingleNull",
                                                      false );

   m_imex_pc_op = new GKPoisson(pp_vorticity_op, a_geometry, a_larmor, 0.);
   m_imex_pc_op->m_model = "Vorticity";
   m_imex_pc_op->m_dt_implicit = 1.;
   m_imex_pc_op_bcs = elliptic_op_bc_factory.create( name,
                                                     ppsp,
                                                     "SingleNull",
                                                     false );

   m_gyropoisson_op = new GKPoisson(pp_vorticity_op, a_geometry, a_larmor, 0.);
   m_gyropoisson_op->m_model = "GyroPoisson";
   m_gyropoisson_op->m_dt_implicit = 1.;
   m_gyropoisson_op_bcs = elliptic_op_bc_factory.create( name,
                                                         ppsp,
                                                         "SingleNull",
                                                         false );

   ParmParse pp_vlasov("vorticity_vlasov");
   m_vlasov = new PS::GKVlasov(pp_vlasov, a_larmor);

   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   m_negativeDivJpar.define(grids, 1, IntVect::Zero);
   m_divJperp.define(grids, 1, IntVect::Zero);

   m_ion_charge_density.define(grids, 1, 2*IntVect::Unit);
   m_ion_mass_density.define(grids, 1, IntVect::Zero);

   m_volume.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume);

   m_divJperp_mag.define(grids, 1, IntVect::Zero);

   // Input the initial conditions
   parseInitialConditions("", "potential");
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
   else {
      setCoreBC( 0., 0., *m_par_cond_op_bcs );
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

   // Compute the divergence of the perpendicular ion current
   computeDivPerpIonCurrentDensity(m_divJperp, a_E_field, a_kinetic_species_phys, m_ion_charge_density, a_time, false);

   CH_STOP(t_div_jperp);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      rhs_data[dit] += m_negativeDivJpar[dit];
      rhs_data[dit] -= m_divJperp[dit];
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

   computeDivPerpIonMagCurrentDensity(m_divJperp_mag, a_E_field, a_kinetic_species, a_time);
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
   else {
      setCoreBC( 0., 0., *m_gyropoisson_op_bcs );
   }

   m_gyropoisson_op->updateBoundaries(*m_gyropoisson_op_bcs);

   const FluidSpecies& sol_species = static_cast<const FluidSpecies&>(*a_fluid_species_phys[a_component]);
   FluidSpecies& rhs_species = static_cast<FluidSpecies&>(*a_rhs[a_component]);

   const LevelData<FArrayBox>& sol = sol_species.cell_var(0);
   LevelData<FArrayBox>& rhs = rhs_species.cell_var(0);

   m_gyropoisson_op->computeFluxDivergence(sol, rhs, false, false);
}

void VorticityOp::defineBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
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

  if (a_im)  m_my_pc_idx_i = a_pc.size();
  else       m_my_pc_idx_e = a_pc.size();

  a_pc.push_back(pc);
  a_dof_list.push_back(dof_list);

  return;
}

void VorticityOp::updateBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
                                  const PS::KineticSpeciesPtrVect&                          a_kin_species_phys,
                                  const FluidSpeciesPtrVect&                                a_fluid_species,
                                  const Real                                                a_time,
                                  const Real                                                a_shift,
                                  const bool                                                a_im,
                                  const int                                                 a_species_idx )
{
  if (a_im) {
    CH_assert(m_my_pc_idx_i >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx_i);
  } else {
    CH_assert(m_my_pc_idx_e >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx_e);
  }

  if (!procID()) {
    std::cout << "    ==> Updating " << _FLUID_OP_PC_ << " preconditioner " 
              << " for VorticityOp " << (a_im ? "RHS " : "LHS ")
              << "of fluid species " << a_species_idx << ".\n";
  }

  PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps> *pc 
    = dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>
      (a_pc[(a_im ? m_my_pc_idx_i : m_my_pc_idx_e)]);
  CH_assert(pc != NULL);
  pc->update(a_kin_species_phys, a_fluid_species, a_time, a_shift, a_im, a_species_idx);

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
   m_gyropoisson_op->setPreconditionerConvergenceParams(0., 1, 0., 1);
   m_gyropoisson_op->solvePreconditioner(r, z);
}

void VorticityOp::updatePCImEx( const PS::KineticSpeciesPtrVect& a_kinetic_species,
                                const double                     a_time,
                                const double                     a_shift )
{
   CH_TIME("VorticityOp::updatePCImEx");

   m_imex_pc_op->m_alpha = a_shift;
   m_imex_pc_op->setOperatorCoefficients( m_ion_mass_density, *m_imex_pc_op_bcs, true );
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
   m_imex_pc_op->setPreconditionerConvergenceParams(0., 1, 0., 1);
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
      bool fourth_order = !m_parallel_current_divergence_op->secondOrder();
      m_vlasov->evalRHS(tmp_rhs_species, this_species, a_E_field.getCellCenteredField(), a_E_field.getPhiNode(),
                        fourth_order, PS::PhaseGeom::MAGNETIC_DRIFT_VELOCITY, a_time);

      // Divide by J to get physical cell averages
      phase_geometry->divideJonValid( rhs_dfn );
      
      moment_op.compute( species_div_Jperp, tmp_rhs_species, PS::ChargeDensityKernel() );
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_divJperp_mag[dit] -= species_div_Jperp[dit];
      }
   }
}

void VorticityOp::computeDivPerpIonCurrentDensity( LevelData<FArrayBox>&             a_div_Jperp,
                                                   const EField&                     a_E_field,
                                                   const PS::KineticSpeciesPtrVect&  a_species_phys,
                                                   const LevelData<FArrayBox>&       a_ion_charge_density,
                                                   const Real&                       a_time,
                                                   const bool                        a_compute_mag_contrib )
{
   CH_TIME("VorticityOp::computeDivPerpIonCurrentDensity");

   const DisjointBoxLayout& grids = m_geometry.gridsFull();

   if ( a_compute_mag_contrib ) {
      // Compute the ion perpendicular ion current divergence due solely to magnetic drifts
      computeDivPerpIonMagCurrentDensity(m_divJperp_mag, a_E_field, a_species_phys, a_time);
   }

   // Compute the J_ExB divergence

#ifndef OLD_DIVJPERP
   
   // Make sure there are enough ghost cells to compute the face average density below
   for (int species(0); species<a_species_phys.size(); species++) {
      PS::LevelData<PS::FArrayBox>& dfn = a_species_phys[species]->distributionFunction();
      CH_assert(dfn.ghostVect() >= 2*PS::IntVect::Unit);
   }

   // Compute the current due to ExB drifts

   LevelData<FluxBox> J_ExB(grids, 1, IntVect::Unit);
   m_geometry.computeBxEIntegrals(a_E_field.getPhiNode(), true, J_ExB);

   //   LevelData<FArrayBox> charge_density(grids, 1, 2*IntVect::Unit);
   //   computeIonChargeDensity(charge_density, a_species_phys);

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
      m_vlasov->evalRHS(tmp_rhs_species, this_species, a_E_field.getCellCenteredField(), a_E_field.getPhiNode(),
                        fourth_order, PS::PhaseGeom::EXB_DRIFT_VELOCITY, a_time);

      // Divide by J to get physical cell averages

      phase_geometry->divideJonValid( rhs_dfn );
      
      moment_op.compute( species_div_Jperp, tmp_rhs_species, PS::ChargeDensityKernel() );

      for (DataIterator dit(a_div_Jperp.dataIterator()); dit.ok(); ++dit) {
         a_div_Jperp[dit] -= species_div_Jperp[dit];
      }
   }

#endif

   // Add the pre-computed divJperp due solely to magnetic drifts
   for (DataIterator dit(a_div_Jperp.dataIterator()); dit.ok(); ++dit) {
      a_div_Jperp[dit] += m_divJperp_mag[dit];
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

      // Compute the mass density for this species
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
      
      // Compute the charge density for this species
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
      a_bc.setBCValue(0,RADIAL_DIR,0,a_core_inner_bv);
      a_bc.setBCValue(0,RADIAL_DIR,1,a_core_outer_bv);
   }
   else if ( typeid(*(m_geometry.getCoordSys())) == typeid(SNCoreCoordSys) ) {
      a_bc.setBCValue(L_CORE,RADIAL_DIR,0,a_core_inner_bv);
      a_bc.setBCValue(L_CORE,RADIAL_DIR,1,a_core_outer_bv);
   }
   else if ( typeid(*(m_geometry.getCoordSys())) == typeid(SingleNullCoordSys) ) {
      a_bc.setBCValue(LCORE,RADIAL_DIR,0,a_core_inner_bv);
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

   setCoreBC( a_scalar_data[0], -a_scalar_data[1], *m_gyropoisson_op_bcs );

   m_gyropoisson_op->setOperatorCoefficients(m_ion_mass_density, *m_gyropoisson_op_bcs, false);
}


void VorticityOp::updatePotentialOldModel( LevelData<FArrayBox>&             a_phi,
                                           EField&                           a_E_field,
                                           const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                           const Vector<Real>&               a_scalar_data,
                                           const Real                        a_dt,
                                           const Real                        a_time )
{
   const DisjointBoxLayout& grids = m_geometry.gridsFull();
      
   // Compute the polarization density term (at the old time) using the operator and boundary
   // conditions set in the previous call of this function or initializeOldModel().
   LevelData<FArrayBox> gkPoissonRHS(grids, 1, IntVect::Zero);
   m_gyropoisson_op->computeFluxDivergence(a_phi, gkPoissonRHS, false);
      
   // Compute the parallel current divergence
   computeIonChargeDensity(m_ion_charge_density, a_kinetic_species);
   m_parallel_current_divergence_op->m_dt_implicit = a_dt;
   m_parallel_current_divergence_op->setOperatorCoefficients(m_ion_charge_density, *m_parallel_current_divergence_op_bcs, true);
   m_parallel_current_divergence_op->computeFluxDivergence(m_ion_charge_density, m_negativeDivJpar, false, true);
      
   // Compute the perpendicular current divergence
   computeDivPerpIonCurrentDensity(m_divJperp, a_E_field, a_kinetic_species, m_ion_charge_density, a_time);

   // Get the RHS for the implicit vorticity solve
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      gkPoissonRHS[dit] += m_negativeDivJpar[dit];
      m_divJperp[dit] *= a_dt;
      gkPoissonRHS[dit] -= m_divJperp[dit];
   }
      
   double Er_lo = a_scalar_data[0];
   double Er_hi = a_scalar_data[1];

   // Do the implicit vorticity solve
   computeIonMassDensity(m_ion_mass_density, a_kinetic_species);
   m_imex_pc_op->m_dt_implicit = a_dt;
   setCoreBC( Er_lo, -Er_hi, *m_imex_pc_op_bcs );
   m_imex_pc_op->setOperatorCoefficients(m_ion_mass_density, *m_imex_pc_op_bcs, true);

   m_imex_pc_op->computePotential( a_phi, gkPoissonRHS );
   m_imex_pc_op->fillInternalGhosts( a_phi );

   // Update nodal phi if supporting the calculation of a divergence-free phase velocity
   a_E_field.interpToNodes(a_phi);

   // Set up the polarization density operator for the next step
   setCoreBC( Er_lo, -Er_hi, *m_gyropoisson_op_bcs );
   m_gyropoisson_op->setOperatorCoefficients(m_ion_mass_density, *m_gyropoisson_op_bcs, false);
}


#include "NamespaceFooter.H"
