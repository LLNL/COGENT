#include "VorticityOp.H"
#include "Directions.H"

#include "EllipticOpBCFactory.H"
#include "LogRectCoordSys.H"
#include "SNCoreCoordSys.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "PhaseGeom.H"
#include "MomentOp.H"
#include "GKOps.H"
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H" 


VorticityOp::VorticityOp( const std::string&  a_pp_str,
                          const MagGeom&      a_geometry,
                          const double        a_larmor,
                          const int           a_verbosity )
   : m_geometry(a_geometry),
     m_sigma_div_e_coefs_set(false),
     m_opt_string(a_pp_str),
     m_my_pc_idx_e(-1),
     m_my_pc_idx_i(-1)
{
   ParmParse pp1(a_pp_str.c_str());
   parseParameters( pp1 );
   if (m_verbosity>0) {
      printParameters();
   }

   const std::string name("potential");
   const std::string prefix( "BC." + name );
   ParmParse ppsp( prefix.c_str() );
   EllipticOpBCFactory elliptic_op_bc_factory;
   ParmParse pp("gkpoisson");

   m_parallel_current_divergence_op = new GKPoisson(pp, a_geometry, a_larmor, 0.);
   m_parallel_current_divergence_op->m_model = "ParallelCurrent";
   m_parallel_current_divergence_op->m_dt_implicit = 1.;
   m_parallel_current_divergence_op_bcs = elliptic_op_bc_factory.create( name,
                                                                         ppsp,
                                                                         "SingleNull",
                                                                         false );
   setCoreBC( 0., 0., *m_parallel_current_divergence_op_bcs );

   m_par_cond_op = new GKPoisson(pp, a_geometry, 0., 0.);
   m_par_cond_op->m_model = "Vorticity";
   m_par_cond_op->m_dt_implicit = 1.;
   m_par_cond_op_bcs = elliptic_op_bc_factory.create( name,
                                                      ppsp,
                                                      "SingleNull",
                                                      false );

   m_imex_pc_op = new GKPoisson(pp, a_geometry, a_larmor, 0.);
   m_imex_pc_op->m_model = "Vorticity";
   m_imex_pc_op->m_dt_implicit = 1.;
   m_imex_pc_op_bcs = elliptic_op_bc_factory.create( name,
                                                     ppsp,
                                                     "SingleNull",
                                                     false );

   m_gyropoisson_op = new GKPoisson(pp, a_geometry, a_larmor, 0.);
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
   CH_TIMERS("VorticityOp::accumulateExplicitRHS");
}


void VorticityOp::accumulateImplicitRHS( FluidSpeciesPtrVect&               a_rhs,
                                         const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                         const FluidSpeciesPtrVect&         a_fluid_species,
                                         const PS::ScalarPtrVect&           a_scalars,
                                         const EField&                      a_E_field,
                                         const int                          a_fluid_vec_comp,
                                         const bool                         a_recompute_kinetic_terms,
                                         const Real                         a_time )
{
   CH_TIMERS("VorticityOp::accumulateImplicitRHS");
   CH_TIMER("sigma_div_e",t_sigma_div_e);
   CH_TIMER("div_jpar",t_div_jpar);
   CH_TIMER("div_jperp",t_div_jperp);

   CFGVar& rhs_fluid_species = *(a_rhs[a_fluid_vec_comp]);
   LevelData<FArrayBox>& rhs_data = rhs_fluid_species.cell_data();

   const CFGVar& sol_fluid_species = *(a_fluid_species[a_fluid_vec_comp]);
   const LevelData<FArrayBox>& sol_data = sol_fluid_species.cell_data();

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

   // Since the evaluation of terms involving the kinetic species is somewhat expensive,
   // the following conditional provides the opportunity to cache them between calls

   if ( a_recompute_kinetic_terms ) {
   
      CH_START(t_div_jpar);

      // Compute the divergence of the parallel ion current
      LevelData<FArrayBox> ion_charge_density( grids, 1, IntVect::Zero );
      computeIonChargeDensity( ion_charge_density, a_kinetic_species_phys );

      m_parallel_current_divergence_op->setOperatorCoefficients( ion_charge_density, *m_parallel_current_divergence_op_bcs, false );
      m_parallel_current_divergence_op->computeFluxDivergence(ion_charge_density, m_negativeDivJpar, false, true);
      
      CH_STOP(t_div_jpar);

      CH_START(t_div_jperp);

      // Compute the divergence of the perpendicular ion current
      computeDivPerpIonCurrentDensity(m_divJperp, a_E_field, a_kinetic_species_phys, a_time);

      CH_STOP(t_div_jperp);
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      rhs_data[dit] += m_negativeDivJpar[dit];
      rhs_data[dit] -= m_divJperp[dit];
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
#if 1
   LevelData<FArrayBox> ion_mass_density(m_geometry.gridsFull(), 1, IntVect::Zero);
   computeIonMassDensity( ion_mass_density, a_kinetic_species_phys );

   if (a_scalars.size() > 0) {
      const Vector<Real>& scalar_data = a_scalars[0]->data();
      setCoreBC( scalar_data[0], -scalar_data[1], *m_gyropoisson_op_bcs );
      //      if (procID()==0) cout << "Er_lo = " << scalar_data[0] << ", Er_hi = " << scalar_data[1] << endl;
   }
   else {
      setCoreBC( 0., 0., *m_gyropoisson_op_bcs );
   }

   m_gyropoisson_op->setOperatorCoefficients( ion_mass_density, *m_gyropoisson_op_bcs, true );
#endif

   const FluidSpecies& sol_species = static_cast<const FluidSpecies&>(*a_fluid_species_phys[a_component]);
   FluidSpecies& rhs_species = static_cast<FluidSpecies&>(*a_rhs[a_component]);

   const LevelData<FArrayBox>& sol = sol_species.cell_data();
   LevelData<FArrayBox>& rhs = rhs_species.cell_data();

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

  const LevelData<FArrayBox>& soln_data   (a_fluid_species.cell_data());
  const DisjointBoxLayout&    grids       (soln_data.disjointBoxLayout());
  const int                   n_comp      (soln_data.nComp());
  const LevelData<FArrayBox>& pMapping    (a_global_dofs.data()); 

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

  if (a_im)  m_my_pc_idx_i = a_pc.size();
  else       m_my_pc_idx_e = a_pc.size();

  a_pc.push_back(pc);
  a_dof_list.push_back(dof_list);

  return;
}

void VorticityOp::updateBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
                                  const PS::KineticSpeciesPtrVect&                          a_kin_species_phys,
                                  const FluidSpeciesPtrVect&                                a_fluid_species,
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
//  if (pc == NULL) {
//    std::cout << "PC idx: " << m_my_pc_idx << ", a_im: " << a_im 
//              << "\n";
//  }
  CH_assert(pc != NULL);
  pc->update(a_kin_species_phys, a_fluid_species, a_shift, a_im, a_species_idx);

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

   const LevelData<FArrayBox>& r = rhs_species.cell_data();
   LevelData<FArrayBox>& z = sol_species.cell_data();

   setZero(z);
   m_gyropoisson_op->setPreconditionerConvergenceParams(0., 1, 0., 1);
   m_gyropoisson_op->solvePreconditioner(r, z);
}

void VorticityOp::updatePCImEx( const PS::KineticSpeciesPtrVect& a_kinetic_species,
                                const double                     a_mshift )
{
   CH_TIME("VorticityOp::updatePCImEx");
   LevelData<FArrayBox> ion_mass_density(m_geometry.gridsFull(), 1, IntVect::Zero);
   computeIonMassDensity( ion_mass_density, a_kinetic_species );

   m_imex_pc_op->m_alpha = a_mshift;
   m_imex_pc_op->setOperatorCoefficients( ion_mass_density, *m_imex_pc_op_bcs, true );
}


void VorticityOp::solvePCImEx( FluidSpeciesPtrVect&              a_fluid_species_solution,
                               const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                               const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                               const int                         a_component )
{
   CH_TIME("VorticityOp::solvePCImEx");
   const FluidSpecies& rhs_species = static_cast<const FluidSpecies&>(*a_fluid_species_rhs[a_component]);
   const LevelData<FArrayBox>& r = rhs_species.cell_data();

   FluidSpecies& sol_species = static_cast<FluidSpecies&>(*a_fluid_species_solution[a_component]);
   LevelData<FArrayBox>& z = sol_species.cell_data();

   setZero(z);
   m_imex_pc_op->setPreconditionerConvergenceParams(0., 1, 0., 1);
   m_imex_pc_op->solvePreconditioner(r, z);
}


void VorticityOp::computeDivPerpIonCurrentDensity( LevelData<FArrayBox>&             a_div_Jperp,
                                                   const EField&                     a_E_field,
                                                   const PS::KineticSpeciesPtrVect&  a_species_phys,
                                                   const Real&                       a_time ) const
{
   CH_TIME("VorticityOp::computeDivPerpIonCurrentDensity");
   LevelData<FArrayBox> species_div_Jperp(a_div_Jperp.disjointBoxLayout(), 1, IntVect::Zero);
   PS::MomentOp& moment_op = PS::MomentOp::instance();

   setZero(a_div_Jperp);
      
   for (int species(0); species<a_species_phys.size(); species++) {
      
      const PS::KineticSpecies& this_species( *(a_species_phys[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      const PS::PhaseGeom& phase_geometry( this_species.phaseSpaceGeometry() );
      PS::KineticSpecies tmp_rhs_species("tmp_rhs", this_species.mass(), this_species.charge(), phase_geometry);
      PS::LevelData<PS::FArrayBox>& rhs_dfn = tmp_rhs_species.distributionFunction();
      rhs_dfn.define( this_species.distributionFunction().getBoxes(), 1, PS::IntVect::Zero);
      
      // Compute the divergence of -Jperp for this species
      bool fourth_order = !m_parallel_current_divergence_op->secondOrder();
      int velocity_option = 1; // Drift terms only
      m_vlasov->evalRHS(tmp_rhs_species, this_species, a_E_field.cell_data(), a_E_field.getPhiNode(), fourth_order, velocity_option, a_time);

      // Divide by J to get physical cell averages
      phase_geometry.divideJonValid( rhs_dfn );
      
      moment_op.compute( species_div_Jperp, tmp_rhs_species, PS::ChargeDensityKernel() );
      
      DataIterator dit( a_div_Jperp.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         a_div_Jperp[dit].minus( species_div_Jperp[dit] );
      }
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


inline
void VorticityOp::parseParameters( ParmParse& a_pp )
{
}


inline
void VorticityOp::printParameters()
{
   if (procID()==0) {
      std::cout << "VorticityOp parameters:" << std::endl;
   }
}


#include "NamespaceFooter.H"
