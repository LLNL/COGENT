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
     m_courant_time_step(DBL_MAX),
     m_geometry(a_geometry),
     //m_is_time_implicit(true),
     m_implicitMhd(false),
     m_implicitViscosity(false),
     m_advScheme("c2"),
     m_species_name(a_species_name),
     m_CGL(0),
     m_initializeBfromVectorPotential(0),
     m_opt_string(a_pp_str),
     m_my_pc_idx(-1),
     m_tau(0.0),
     m_etaMin(0.0)
{
   ParmParse pp(a_pp_str.c_str());
   parseParameters( pp );
   if (m_verbosity>0) {
      printParameters();
   }

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
   m_density_cbc = m_fluid_bc.size();
   parseInitialConditions(a_species_name, variable_name);
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );
 
   // set IC and create BC for momentum density
   //
   variable_name = "momentumDensity";
   CH_assert( ppic.contains("momentumDensity_0.function") );
   CH_assert( ppic.contains("momentumDensity_1.function") );
   parseInitialConditions(a_species_name, "momentumDensity_0");
   parseInitialConditions(a_species_name, "momentumDensity_1");
   if (ppic.contains("momentumDensity_2.function")) {
      parseInitialConditions(a_species_name, "momentumDensity_2");
   }
   m_momentumDensity_cbc = m_fluid_bc.size();
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );
 

   
   // construct elliptic operators for viscosity
   //
   std::string name, prefix, prefix_BC;
   EllipticOpBCFactory elliptic_op_bc_factory;
   
   name = "mom0_visc";
   prefix = name;
   ParmParse ppVisc_mom0( prefix.c_str() );
   m_diffusionOp_mom0 = new Diffusion(ppVisc_mom0, a_geometry);
   //
   prefix_BC = "BC." + name;
   ParmParse ppVisc_mom0_BC( prefix_BC.c_str() );
   m_bc_mom0 = elliptic_op_bc_factory.create( name, ppVisc_mom0_BC,
                                              *(a_geometry.getCoordSys()), false );
   
   name = "mom1_visc";
   prefix = name;
   ParmParse ppVisc_mom1( prefix.c_str() );
   m_diffusionOp_mom1 = new Diffusion(ppVisc_mom1, a_geometry);
   //
   prefix_BC = "BC." + name;
   ParmParse ppVisc_mom1_BC( prefix_BC.c_str() );
   m_bc_mom1 = elliptic_op_bc_factory.create( name, ppVisc_mom1_BC,
                                              *(a_geometry.getCoordSys()), false );
   
   if(SpaceDim==3) {
      name = "mom2_visc";
      prefix = name;
      ParmParse ppVisc_mom2( prefix.c_str() );
      m_diffusionOp_mom2 = new Diffusion(ppVisc_mom2, a_geometry);
      //
      prefix_BC = "BC." + name;
      ParmParse ppVisc_mom2_BC( prefix_BC.c_str() );
      m_bc_mom2 = elliptic_op_bc_factory.create( name, ppVisc_mom2_BC,
                                                 *(a_geometry.getCoordSys()), false );
      //
      //
      name = "Apar_init";
      prefix = name;
      ParmParse ppVisc_Apar( prefix.c_str() );
      m_diffusionOp_Apar = new Poisson(ppVisc_Apar, a_geometry);
      //
      prefix_BC = "BC." + name;
      ParmParse ppVisc_Apar_BC( prefix_BC.c_str() );
      m_bc_Apar = elliptic_op_bc_factory.create( name, ppVisc_Apar_BC,
                                                 *(a_geometry.getCoordSys()), false );
   } 
   else {
      m_diffusionOp_mom2 = NULL;
      m_bc_mom2 = NULL;
      m_diffusionOp_Apar = NULL;
      m_bc_Apar = NULL;
   }

   //
   ///////////////////////


   // set IC and create BC for energy density
   //
   variable_name = "energyDensity";
   CH_assert( ppic.contains("energyDensity.function") 
           || ppic.contains("energyDensity_0.function") );
   if( ppic.contains("energyDensity.function") ) {
      parseInitialConditions(a_species_name, variable_name);
   } 
   else {
      CH_assert( ppic.contains("energyDensity_1.function") );
      parseInitialConditions(a_species_name, "energyDensity_0");
      parseInitialConditions(a_species_name, "energyDensity_1");
      m_CGL = 1;
   }
   m_energyDensity_cbc = m_fluid_bc.size();
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );
   
   // set IC and create BC for in-plane magnetic field
   //
   if (ppic.contains("magneticField_0.function")) {
      variable_name = "magneticField";
      CH_assert( ppic.contains("magneticField_0.function") );
      CH_assert( ppic.contains("magneticField_1.function") );
      parseInitialConditions(a_species_name, "magneticField_0");
      parseInitialConditions(a_species_name, "magneticField_1");
      if (ppic.contains("magneticField_2.function")) {
         parseInitialConditions(a_species_name, "magneticField_2");
      }
      m_magneticField_fbc = m_fluid_face_var_bc.size();
      m_fluid_face_var_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name,
                                     a_geometry.getCoordSys()->type(), false) );
   }

   // set IC and create BC for virtual momentum density
   //
   if (ppic.contains("momentumDensity_virtual.function")) {
      variable_name = "momentumDensity_virtual";
      parseInitialConditions(a_species_name, variable_name);
      m_momentumDensity_virtual_cbc = m_fluid_bc.size();
      m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                            a_geometry.getCoordSys()->type(), false) );
   }
   
   if (ppic.contains("magneticField_virtual.function")) {
      variable_name = "magneticField_virtual";
      parseInitialConditions(a_species_name, variable_name);
      m_magneticField_virtual_cbc = m_fluid_bc.size();
      m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                            a_geometry.getCoordSys()->type(), false) );
   }
   
   // get fluid ghost cells and define member LevelData's
   Real ghosts = 4;
   ParmParse ppgksys("gksystem");
   ppgksys.query( "fluid_ghost_width", ghosts );
   IntVect ghostVect = ghosts*IntVect::Unit;
   if(m_advScheme == "weno5" || m_advScheme == "bweno" || m_advScheme == "uw5") {
      CH_assert(ghosts >= 3);
   }
   const DisjointBoxLayout& grids( m_geometry.grids() );
   defineLevelDatas( grids, ghostVect );
  
}


IdealMhdOp::~IdealMhdOp()
{
   //delete m_fluid_variable_bc;
   for (std::vector<FluidVarBC*>::iterator it = m_fluid_bc.begin(); it != m_fluid_bc.end(); ++it)
   {
      delete (*it);
   }
   m_fluid_bc.clear();
   delete m_bc_mom0; 
   delete m_bc_mom1; 
   delete m_bc_mom2; 
   delete m_bc_Apar; 
   delete m_diffusionOp_mom0; 
   delete m_diffusionOp_mom1; 
   delete m_diffusionOp_mom2; 
   delete m_diffusionOp_Apar; 
   
   for (std::vector<FluidVarBC*>::iterator it = m_fluid_face_var_bc.begin(); it != m_fluid_face_var_bc.end(); ++it)
   {
      delete (*it);
   }
   m_fluid_face_var_bc.clear();
   
}


void IdealMhdOp::accumulateMhdRHS( FluidSpecies&  a_rhs_fluid,
                             const FluidSpecies&  a_soln_fluid,
                             const Real           a_time )
{
   CH_TIME("IdealMhdOp::accumulateMhdRHS()");
   
   // Get RHS fluid species
   //
   LevelData<FArrayBox>& rhs_data( a_rhs_fluid.cell_var(0) );
   const DisjointBoxLayout& grids( rhs_data.getBoxes() );

   // Get solution fluid species
   //
   //const LevelData<FArrayBox>& soln_data( a_soln_fluid.cell_var(0) );
   //const IntVect& ghostVect( soln_data.ghostVect() );
   //cout << "JRA: soln_data.ghostVect() = " << ghostVect << endl;   
   
   //CH_assert(soln_fluid.m_evolve_massDensity==1);
   CH_assert(a_soln_fluid.m_evolve_momentumDensity==1);
   CH_assert(a_soln_fluid.m_evolve_energyDensity==1);
  
   setCellCenterValues( a_soln_fluid );
   
   if(m_CGL) {
      setCellCenterFluxesCGL( a_soln_fluid );
   } 
   else {
      setCellCenterFluxes( a_soln_fluid );
   }

   // compute flux-freezing speed at cell-center for each direction
   // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/N); i = q0, q1, q2
   //
   const double gamma = a_soln_fluid.m_gamma;
   setMappedCspeed( gamma, 1, 1 );
   
   SpaceUtils::upWindToFaces(m_CspeedR_norm, m_Cspeed_cc, m_CspeedR_norm, "c2");
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_CspeedL_norm[dit].copy(m_CspeedR_norm[dit],m_CspeedL_norm[dit].box());
      m_CspeedL_norm[dit].negate();
   }
   
   setMagneticFieldTerms( a_soln_fluid );
   
   setFaceCenteredFluxes( a_soln_fluid );

   enforceFluxBCs( a_soln_fluid, a_time );

   updateRHSs( a_rhs_fluid, a_soln_fluid );
    
}

void IdealMhdOp::accumulateViscRHS( FluidSpecies&  a_rhs_fluid,
                              const FluidSpecies&  a_soln_fluid,
                              const int            a_accumulateMom,
                              const int            a_accumulateEne )
{
   CH_TIME("IdealMhdOp::accumulateViscRHS()");
   
   setCellCenterValues( a_soln_fluid ); 
 
   setViscosityCoefficient( m_etaVisc, m_etaVisc_cf, a_soln_fluid );
 
   updateRHSs_visc( a_rhs_fluid, a_soln_fluid, a_accumulateMom, a_accumulateEne );

}

void IdealMhdOp::accumulateExplicitRHS( FluidSpeciesPtrVect&        a_rhs,
                                  const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                  const FluidSpeciesPtrVect&        a_fluid_species,
                                  const PS::ScalarPtrVect&          a_scalars,
                                  const EField&                     a_E_field,
                                  const int                         a_fluidVecComp,
                                  const Real                        a_time )
{
   CH_TIME("IdealMhdOp::accumulateExplicitRHS()");
   
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
      
   if(!m_implicitMhd) {
      accumulateMhdRHS( rhs_fluid, soln_fluid, a_time );
      if(m_tau>0.0 || m_etaMin>0.0) {
         int accumulateMom = 1;
         if(m_implicitViscosity) accumulateMom = 0;
         accumulateViscRHS( rhs_fluid, soln_fluid, accumulateMom, 1 );
      }
   }

}


void IdealMhdOp::accumulateImplicitRHS( FluidSpeciesPtrVect&        a_rhs,
                                  const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                  const FluidSpeciesPtrVect&        a_fluid_species,
                                  const PS::ScalarPtrVect&          a_scalars,
                                  const EField&                     a_E_field,
                                  const int                         a_fluidVecComp,
                                  const Real                        a_time )
{
   CH_TIME("IdealMhdOp::accumulateImplicitRHS()");

   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );

   if(m_implicitMhd) {
      accumulateMhdRHS( rhs_fluid, soln_fluid, a_time );
   }
   if((m_tau>0.0 || m_etaMin>0.0) && m_implicitViscosity) {
      int accumulateMom = 1;
      int accumulateEne = 0;
      if(m_implicitMhd) accumulateEne = 1;
      accumulateViscRHS( rhs_fluid, soln_fluid, accumulateMom, accumulateEne );
   }
   
}

void IdealMhdOp::getMemberVar( LevelData<FArrayBox>&  a_Var,
                         const FluidSpecies&          a_fluid_species,
                         const string&                a_name ) const
{
   const DisjointBoxLayout& grids( m_geometry.grids() );
   const IntVect ghostVect = 0*IntVect::Unit;
   //const IntVect& ghostVect( soln_data.ghostVect() );
   
   if(a_name=="electricField") {
      a_Var.define(grids, SpaceDim, ghostVect);
      if(a_fluid_species.m_evolve_magneticField_virtual==1) {
         const LevelData<FArrayBox>& magField_phys_cc( a_fluid_species.cell_var("magneticField_virtual") );
         LevelData<FArrayBox>& vel_phys_cc = m_dummyFArray_spaceDim;
         a_fluid_species.velocity(vel_phys_cc);  // in-plane velocity vector
         computePhysicalElectricField( a_Var, vel_phys_cc, magField_phys_cc, 1 );
      }
      else if(a_fluid_species.m_evolve_magneticField==1) {
         const LevelData<FluxBox>& magField_contra_cf( a_fluid_species.face_var("magneticField") );
         LevelData<FArrayBox>& vel_phys_cc = m_dummyFArray_spaceDim;
         a_fluid_species.velocity(vel_phys_cc);  // in-plane velocity vector
         LevelData<FArrayBox> magField_phys_cc;
         magField_phys_cc.define(grids, SpaceDim, ghostVect);
         SpaceUtils::interpFaceVectorToCell(magField_phys_cc,magField_contra_cf,"c2");
         computePhysicalElectricField( a_Var, vel_phys_cc, magField_phys_cc, 0 );
      }
      //for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
      //   a_Var[dit].copy(m_E0_cc[dit]);
      //}
   }
   if(a_name=="electricField_upwinded") { // need ghost cells to calc here, so one stage out of sync 
      a_Var.define(grids, SpaceDim, ghostVect);
      SpaceUtils::interpEdgesToCell(a_Var,m_E0_ce,"c2"); // covar at cell center
      m_geometry.convertPhysToCovar(a_Var,1);            // phys at cell center
   }
   if(a_name=="currentDensity") {
      a_Var.define(grids, SpaceDim, ghostVect);
      const LevelData<FluxBox>& magField_contra_cf( a_fluid_species.face_var("magneticField") );
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         m_dummyFlux_oneComp[dit].copy(magField_contra_cf[dit]);
      }
      SpaceUtils::exchangeFluxBox(m_dummyFlux_oneComp);
      const double dummy_time = 0.0;
      int this_fbc = m_magneticField_fbc;
      m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( a_fluid_species, m_dummyFlux_oneComp, dummy_time );
      computeMappedCurrentDensity( a_Var, m_dummyFlux_oneComp ); // Ja*curl(B)
      //for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
      //   a_Var[dit].copy(m_JaJ0_cc[dit]); // off by time stage
      //}
      m_geometry.divideJonValid(a_Var);
   }
   if(a_name=="divB") {
      a_Var.define(grids, 1, IntVect::Zero);
      CH_assert(a_fluid_species.m_evolve_magneticField==1);
      const LevelData<FluxBox>& magField_contra_cf( a_fluid_species.face_var("magneticField") );
      computePhysicalDivergence( a_Var, magField_contra_cf, 0);
   }
   if(a_name=="divCurlE") {
      a_Var.define(grids, 1, IntVect::Zero);
      computePhysicalDivergence( a_Var, m_curlE0_cf, 1);
   }
   if(a_name=="magPressure") {
      a_Var.define(grids, 1, IntVect::Zero);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_Bpressure[dit]);
      }
   }
   if(a_name=="Apar") {
      a_Var.define(grids, m_Apar_cc.nComp(), IntVect::Zero);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_Apar_cc[dit]);
      }
   }
   if(a_name=="Cspeed") { // mapped Cspeed
      a_Var.define(grids, m_Cspeed_cc.nComp(), IntVect::Zero);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_Cspeed_cc[dit]);
      }
   }
   if(a_name=="JaJcrossB") {
      a_Var.define(grids, m_JaJcrossB.nComp(), IntVect::Zero);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_JaJcrossB[dit]);
      }
   }
   if(a_name=="JaJdotE") {
      a_Var.define(grids, m_JaJdotE.nComp(), IntVect::Zero);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_JaJdotE[dit]);
      }
   }

}

void IdealMhdOp::defineBlockPC(  std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>& a_pc,
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
  if (a_im && (m_implicitViscosity || m_implicitMhd)) {
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


void IdealMhdOp::updateBlockPC(  std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>& a_pc,
                              const PS::KineticSpeciesPtrVect&                                a_kin_species_phys,
                              const FluidSpeciesPtrVect&                                      a_fluid_species,
                              const Real                                                      a_time,
                              const int                                                       a_step,
                              const int                                                       a_stage,
                              const Real                                                      a_shift,
                              const bool                                                      a_im,
                              const int                                                       a_species_idx )
{
  if (a_im && (m_implicitViscosity || m_implicitMhd)) {
    CH_assert(m_my_pc_idx >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx);
  
    if (!procID()) {
      std::cout << "    ==> Updating " << _FLUID_OP_PC_ << " preconditioner " 
                << " for IdealMhd fluid species " << a_species_idx << ".\n";
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


void IdealMhdOp::updatePCImEx(const FluidSpeciesPtrVect&       a_fluid_species, 
                              const PS::KineticSpeciesPtrVect& a_kinetic_species,
                              const double                     a_time,
                              const int                        a_step,
                              const int                        a_stage,
                              const double                     a_shift,
                              const int                        a_component)
{
   CH_TIME("IdealMhdOp::updatePCImEx()");
   
   // Get computational fluid species data
   //
   const FluidSpecies& fluid_species = static_cast<const FluidSpecies&>(*a_fluid_species[a_component]);
   const LevelData<FArrayBox>& Jrho( fluid_species.cell_var("density") );
   m_alpha_PC = a_shift;   

   for (DataIterator dit(m_shift_PC.dataIterator()); dit.ok(); ++dit) {
      m_Jrho_PC[dit].copy(Jrho[dit]);
      m_shift_PC[dit].setVal(m_alpha_PC);
      m_shift_PC[dit].mult(m_Jrho_PC[dit]);
   }
   m_geometry.divideJonValid(m_shift_PC);

   m_diffusionOp_mom0->updateImExPreconditioner( m_shift_PC, *m_bc_mom0 );
   m_diffusionOp_mom1->updateImExPreconditioner( m_shift_PC, *m_bc_mom1 );
   if(SpaceDim==3) m_diffusionOp_mom2->updateImExPreconditioner( m_shift_PC, *m_bc_mom2 );
   
   // Create elliptic coefficients
   // Need to update m_etaVisc_cf ? careful, m_pressure above is mapped
   //
   const double a_dummy_time=0.0;
   computeDiffusionCoefficients(m_ellip_coeff, m_ellip_coeff_mapped, m_etaVisc_cf, a_dummy_time);
   m_diffusionOp_mom0->setOperatorCoefficients(m_ellip_coeff, m_ellip_coeff_mapped, *m_bc_mom0);
   m_diffusionOp_mom1->setOperatorCoefficients(m_ellip_coeff, m_ellip_coeff_mapped, *m_bc_mom1);
   if(SpaceDim==3) m_diffusionOp_mom2->setOperatorCoefficients(m_ellip_coeff, m_ellip_coeff_mapped, *m_bc_mom2);
   //SpaceUtils::inspectFluxBox(m_ellip_coeff,1);
   //m_geometry.plotFaceData( "m_ellip_coeff_mapped", m_ellip_coeff_mapped, a_dummy_time);
   
}


void IdealMhdOp::solvePCImEx( FluidSpeciesPtrVect&        a_fluid_species_solution,
                        const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                        const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                        const int                         a_component )
{
   CH_TIME("IdealMhdOp::solvePCImEx()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // This function computes solution to the preconditioned problem.
   // treating Fi(y) as linear operator, the equation solved here is
   // alpha*y - Fi(y) = RHS, y is mapped state vector
   //
   // momentum density and energy density have implicit parts
   // 

   if (!m_implicitViscosity) {
     MayDay::Error("Error in IdealMhdOp::solvePCImEx(): m_implicitViscosity is false." 
                   "This function shouldn't have been called at all!");
   }

   // Get RHS data for each implicit variable
   //
   const FluidSpecies& rhs_species = static_cast<const FluidSpecies&>(*a_fluid_species_rhs[a_component]);
   const LevelData<FArrayBox>& rhs_mom( rhs_species.cell_var("momentumDensity") );
 
   // Get initial solution state for each implicit variable (y)
   //
   FluidSpecies& soln_species = static_cast<FluidSpecies&>(*a_fluid_species_solution[a_component]);
   LevelData<FArrayBox>& z_mom = soln_species.cell_var("momentumDensity");

   // solve modified linear precondition problem for momentum:
   // alpha*rho*z_mom - nabla\cdot(D*nabla(z_mom)) = rhs_mom/J
   // with 
   // z_mom = y_mom/J/rho (velocity) 
   //
   LevelData<FArrayBox> this_rhs_mom(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> this_z_mom(grids, 1, IntVect::Zero);

   for (int n=0; n<SpaceDim; n++) {

      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         this_rhs_mom[dit].copy(rhs_mom[dit],n,0,1);
         this_z_mom[dit].copy(z_mom[dit],n,0,1);
         this_z_mom[dit].divide(m_Jrho_PC[dit],0,0,1);
         this_z_mom[dit] /= m_alpha_PC; // initial guess is alph*y = rhs
      }
      m_geometry.divideJonValid(this_rhs_mom);

      // Get convergence parameters for diffusion MBSolver
      //
      double tol = 0.;
      double max_iter = 1;
      double precond_tol = 0.;
      double precond_max_iter = 1;

      if(n==0) { // Solve for z_mom0
         tol = m_diffusionOp_mom0->getPrecondTol();
         max_iter = m_diffusionOp_mom0->getPrecondMaxIter();
         m_diffusionOp_mom0->setImExPreconditionerConvergenceParams(tol, max_iter, precond_tol, precond_max_iter);
         m_diffusionOp_mom0->solveImExPreconditioner(this_rhs_mom, this_z_mom);
      }
      if(n==1) { // Solve for z_mom1
         tol = m_diffusionOp_mom1->getPrecondTol();
         max_iter = m_diffusionOp_mom1->getPrecondMaxIter();
         m_diffusionOp_mom1->setImExPreconditionerConvergenceParams(tol, max_iter, precond_tol, precond_max_iter);
         m_diffusionOp_mom1->solveImExPreconditioner(this_rhs_mom, this_z_mom);
      }
      if(n==2) { // Solve for z_mom2
         tol = m_diffusionOp_mom2->getPrecondTol();
         max_iter = m_diffusionOp_mom2->getPrecondMaxIter();
         m_diffusionOp_mom2->setImExPreconditionerConvergenceParams(tol, max_iter, precond_tol, precond_max_iter);
         m_diffusionOp_mom2->solveImExPreconditioner(this_rhs_mom, this_z_mom);
      }
   
      // Convert the elliptic solution back to the state vector components
      //
      for(DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         z_mom[dit].copy(this_z_mom[dit],0,n,1);
         z_mom[dit].mult(m_Jrho_PC[dit],0,n,1);
      }

   }
 
}


void IdealMhdOp::computeDiffusionCoefficients( LevelData<FluxBox>&  a_ellip_coeff,
                                               LevelData<FluxBox>&  a_ellip_coeff_mapped,
                                         const LevelData<FluxBox>&  a_D_fc,
                                         const Real                 a_time) const
{
   CH_TIME("IdealMhdOp::computeDiffusionCoefficients()");
  
   // compute isotropic diffusion coefficient needed for elliptic solver
   // D = 4/3*etaVisc*I
   //

   // parallel and perpendicular tensor coefficients
   //
   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& par_coeff  = m_geometry.getEllipticOpParCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();
   
   
   // multiply identity tensor by diffusion coefficient
   //
   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_ellip_coeff[dit].copy(perp_coeff[dit]);
      a_ellip_coeff[dit] += par_coeff[dit];
      a_ellip_coeff_mapped[dit].copy(perp_coeff_mapped[dit]);
      a_ellip_coeff_mapped[dit] += par_coeff_mapped[dit];
      for (int dir = 0; dir < SpaceDim; dir++) {
         for (int n = 0; n < SpaceDim*SpaceDim; n++) {
            a_ellip_coeff[dit][dir].mult(a_D_fc[dit][dir],0,n,1);
            a_ellip_coeff_mapped[dit][dir].mult(a_D_fc[dit][dir],0,n,1);
         }
      }
      a_ellip_coeff[dit] *= 4.0/3.0;
      a_ellip_coeff_mapped[dit] *= 4.0/3.0;
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
                               const LevelData<FArrayBox>&  a_Flux_mapped_cc,
                               const LevelData<FArrayBox>&  a_Jf_cc,
                               const LevelData<FArrayBox>&  a_Cspeed_cc,
                               const LevelData<FluxBox>&    a_CspeedR_norm,
                               const LevelData<FluxBox>&    a_CspeedL_norm ) const
{
   CH_TIME("IdealMhdOp::computeNTFfaceArea()");
   CH_assert(a_Flux_norm.nComp() == 1);
   CH_assert(a_Flux_mapped_cc.nComp() == SpaceDim);
   
   const DisjointBoxLayout& grids( a_Jf_cc.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
   
      // compute Lax splitting on mapped flux and interpolate FluxL and FluxR to face norms
      //
      //   get left and right going flux at cell-center
      //   fluxR = 0.5*(flux + Cspeed*fun),
      //   fluxL = 0.5*(flux - Cspeed*fun),
      //   Cspeed = abs(max(eigenValue of Flux Jacobian))
      //

      const Box& cell_box( a_Jf_cc[dit].box() );
      FArrayBox FluxL_cc(cell_box, SpaceDim);
      FArrayBox FluxR_cc(cell_box, SpaceDim);
      const FArrayBox& Flux_on_patch= a_Flux_mapped_cc[dit];
      const FArrayBox& Cspeed_on_patch = a_Cspeed_cc[dit];
      const FArrayBox& Jf_on_patch = a_Jf_cc[dit];

      FORT_LAX_FLUX_SPLITTING( CHF_BOX(cell_box), 
                               CHF_CONST_FRA(FluxL_cc),
                               CHF_CONST_FRA(FluxR_cc),
                               CHF_CONST_FRA(Flux_on_patch),
                               CHF_CONST_FRA(Cspeed_on_patch),
                               CHF_CONST_FRA1(Jf_on_patch,0) );

      const Box& grid_box( grids[dit] );             // this box has no ghost
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
         SpaceUtils::interpToFacesWENO(thisFluxL_norm, FluxL_cc, thisCspeedL_norm, 
                                       thisSmooth, grid_box, m_advScheme);
         SpaceUtils::interpToFacesWENO(thisFluxR_norm, FluxR_cc, thisCspeedR_norm, 
                                       thisSmooth, grid_box, m_advScheme);
      }
      else {
         SpaceUtils::upWindToFaces(thisFluxL_norm, FluxL_cc, thisCspeedL_norm, 
                                   grid_box, m_advScheme);
         SpaceUtils::upWindToFaces(thisFluxR_norm, FluxR_cc, thisCspeedR_norm, 
                                   grid_box, m_advScheme);
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
      }
   }
   //a_Flux_norm.exchange();
 
}

void IdealMhdOp::computeIdealEatEdges( LevelData<EdgeDataBox>&  a_Edge_covar,
                                 const LevelData<FArrayBox>&    a_V_phys_cc,
                                 const LevelData<FluxBox>&      a_B_contra_cf,
                                 const LevelData<FArrayBox>&    a_Cspeed_cc ) const
{
   CH_TIME("IdealMhdOp::computeIdealEatEdges()");
   CH_assert(a_Edge_covar.nComp() == 1);
   CH_assert(a_B_contra_cf.nComp() == 1);
   CH_assert(a_V_phys_cc.nComp() == SpaceDim);
   CH_assert(a_Cspeed_cc.nComp() == SpaceDim);
  
   // multiply contravariant B on faces by Jacobian
   //
   const DisjointBoxLayout& grids( m_geometry.grids() );
   LevelData<FluxBox>& JaB_contra_cf(m_dummyFlux_oneComp);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      JaB_contra_cf[dit].copy(a_B_contra_cf[dit]);
   } 
   m_geometry.multJonFaces(JaB_contra_cf);
 
   // compute constrained covariant E on cell edges
   //
   string this_advScheme = m_advScheme;
   SpaceUtils::constrainedCrossProductOnEdges( a_Edge_covar, JaB_contra_cf,
                                               a_V_phys_cc, a_Cspeed_cc, this_advScheme );
   
   // redefine dummy to infty for safety
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      JaB_contra_cf[dit].setVal(1.0/0.0);
   } 

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
      FArrayBox JaVB_covar_cc(cell_box, a_V_phys_cc.nComp());
      FArrayBox JaB_contr_cc(cell_box, a_B_contr_cc.nComp());
      convertPhysFluxToMappedFlux(JaVB_covar_cc,a_V_phys_cc[dit],m_Nmatrix[dit]);
      for (int n=0; n<SpaceDim; n++) {
         JaVB_covar_cc.mult(a_B_contr_cc[dit],cell_box,cell_box,0,n,1);
      }
      JaB_contr_cc.copy(a_B_contr_cc[dit], cell_box );
      JaB_contr_cc.mult(m_Jacobian[dit]);
  
      // compute Lax splitting and interpolate FluxL and FluxR to face norms
      //
      FArrayBox FluxL_cc(cell_box, SpaceDim);
      FArrayBox FluxR_cc(cell_box, SpaceDim);
      const FArrayBox& Cspeed_on_patch = a_Cspeed_cc[dit];

      FORT_LAX_FLUX_SPLITTING( CHF_BOX(cell_box), 
                               CHF_FRA(FluxL_cc),
                               CHF_FRA(FluxR_cc),
                               CHF_CONST_FRA(JaVB_covar_cc),
                               CHF_CONST_FRA(Cspeed_on_patch),
                               CHF_CONST_FRA1(JaB_contr_cc,0) );
     
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
                           const double         a_time )
{
   CH_TIME("IdealMhdOp::fillGhostCells()");
   for (int n=0; n<a_species_phys.num_cell_vars(); ++n) {
      m_fluid_bc.at(n)->apply( a_species_phys, a_time );
   }
   
   for (int n=0; n<a_species_phys.num_face_vars(); ++n) {
      const std::string& this_var_name = a_species_phys.face_var_name(n);
      CH_assert(m_fluid_face_var_bc.at(n)->isForVariable(this_var_name));
      LevelData<FluxBox>& faceVar( a_species_phys.face_var(n) );
      m_fluid_face_var_bc.at(n)->applyFluxBC( a_species_phys, faceVar, a_time );
   }
}


void IdealMhdOp::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "implicitViscosity", m_implicitViscosity);
   a_pp.query( "implicitMhd", m_implicitMhd);

   //GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
   //std::string grid_function_name;
   
   if (a_pp.contains("advScheme")) {
      a_pp.get("advScheme", m_advScheme );
   }
   
   if (a_pp.contains("tau")) {
      a_pp.get("tau", m_tau );
   }
   
   if (a_pp.contains("etaMin")) {
      a_pp.get("etaMin", m_etaMin );
   }

   a_pp.query("initializeBfromVectorPotential", m_initializeBfromVectorPotential);
   
}


void IdealMhdOp::printParameters()
{
   if (procID()==0) {
      std::cout << "IdealMhdOp parameters:" << std::endl;
      std::cout << " advScheme  =  " << m_advScheme << std::endl;
      std::cout << " tau     =  " << m_tau << std::endl;
      std::cout << " etaMin  =  " << m_etaMin << std::endl;
      std::cout << " implicitMhd  =  " << m_implicitMhd << std::endl;
      std::cout << " implicitVisosity  =  " << m_implicitViscosity << std::endl;
      std::cout << " initializeBfromVectorPotential = " << m_initializeBfromVectorPotential << std::endl;
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
   
   Vector<double> Cmax( local_Cmax );
   double this_dt_max=DBL_MAX;
   for (int dir=0; dir<SpaceDim; ++dir) {
#ifdef CH_MPI
      MPI_Allreduce( &local_Cmax[dir], &Cmax[dir], 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif
      this_dt_max = Min( dX_mapped[dir]/Cmax[dir], this_dt_max);
      //if(!procID()) {
      //   cout << "dir = " << dir << endl;
      //   cout << "dX_mapped[dir]  = " << dX_mapped[dir] << endl;
      //   cout << "Cmax[dir]       = " << Cmax[dir] << endl;
      //   cout << "this_dt_max     = " << this_dt_max << endl;
      //}
   }

   m_courant_time_step = this_dt_max;
   //cout << "m_courant_time_step = " << m_courant_time_step << endl;

}

void IdealMhdOp::setMappedCspeed( const double  a_gamma,
                                  const int     a_Ptherm,
                                  const int     a_Pmag )
{
   CH_TIME("IdealMhdOp::setMappedCspeed()");

   // compute flux-freezing speed at cell-center for each direction
   // Cspeed_i = |NTVdotqihot| + |N_i|sqrt(gamma*P/N + B^2/N); i = q0, q1, q2

   CH_assert(m_Cspeed_cc.nComp() == SpaceDim);
   
   const DisjointBoxLayout& grids( m_Cspeed_cc.getBoxes() ); 
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& rho_on_patch  = m_rhoDen_cc[dit];   
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
                                     CHF_CONST_INT(a_Ptherm),
                                     CHF_CONST_INT(a_Pmag),
                                     CHF_FRA(C_on_patch) );
      
   }

}

Real IdealMhdOp::computeDtImExTI( const FluidSpeciesPtrVect&  a_fluid_comp )
{
   CH_TIME("Real IdealMhdOp::computeDtImExTI()");   
   
   double thisdt=DBL_MAX;
   if(m_courant_time_step==DBL_MAX) {
      thisdt = computeDtMhd( a_fluid_comp );
      return thisdt;
   }
   else {
      return m_courant_time_step;
   }

}

Real IdealMhdOp::computeDtExplicitTI( const FluidSpeciesPtrVect&  a_fluid_comp )
{
   CH_TIME("Real IdealMhdOp::computeDtExplicit()");   
   double thisdt=DBL_MAX;
   thisdt = computeDtMhd( a_fluid_comp );
   return thisdt;
}

Real IdealMhdOp::computeDtMhd( const FluidSpeciesPtrVect&  a_fluid_comp )
{
   CH_TIME("Real IdealMhdOp::computeDtMhd()");   
   // get physical values for m_species_name and calculate freezing speed
   // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/rho); i = q0, q1, q2
   //
   double gamma;
   for (int species(0); species<a_fluid_comp.size(); species++) {
      const FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_comp[species])) );
      const std::string species_name( fluid_species.name() );
      //species_name = fluid_species.name();
      if(species_name==m_species_name) {
         gamma = fluid_species.m_gamma;
         fluid_species.massDensity(m_rhoDen_cc);
         fluid_species.velocity(m_velocity);  // in-plane velocity vector
         fluid_species.pressure(m_pressure);
         fluid_species.Bpressure(m_Bpressure);
         
         // Calls for physical variables abve are from mapped state
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

   setMappedCspeed( gamma, 1, 1 );
   setCourantTimeStep(m_Cspeed_cc);   
   
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
   m_Xphys.define(a_grids, SpaceDim, a_ghostVect);
   m_g_y.define(a_grids, 1, a_ghostVect); // used only when solving virtual B
   m_cellVol.define(a_grids, 1, IntVect::Zero);
   double twoPi = 2.0*Pi;
   //double mapVol;
   m_geometry.getXphys(m_Xphys);
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(a_grids[dit]);
      m_twoDaxisymm = coord_sys.isAxisymmetric();
      if(m_twoDaxisymm) { // g_y is covariant basis vector in virtual direction
         m_g_y[dit].copy(m_Xphys[dit],0,0,1);
         m_g_y[dit] *= twoPi;
      } 
      else { 
         m_g_y[dit].setVal(1.0);
      }
   }
   m_geometry.getJ(m_Jacobian);
   m_geometry.getN(m_Nmatrix);
   m_geometry.getCellVolumes(m_cellVol);
   
   // containers for viscosity tensor
   //
   m_etaVisc.define(a_grids, 1, a_ghostVect);
   m_momVisc_source.define(a_grids, 1, IntVect::Zero);
   m_etaVisc_cf.define(a_grids, 1, 1*IntVect::Unit);
   m_JaW_cf.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);
   m_m0JaFluxVisc_cf.define(a_grids, SpaceDim, a_ghostVect);
   m_m1JaFluxVisc_cf.define(a_grids, SpaceDim, a_ghostVect);
   m_m2JaFluxVisc_cf.define(a_grids, SpaceDim, a_ghostVect);
   m_enJaFluxVisc_cf.define(a_grids, SpaceDim, a_ghostVect);
   m_m0JaFluxVisc_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_m1JaFluxVisc_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_m2JaFluxVisc_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_momJaFluxVisc_norm.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_enJaFluxVisc_norm.define(a_grids, 1, 1*IntVect::Unit);
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_momVisc_source[dit].setVal(0.0);
   }
   
   // preconditioner stuff
   //
   m_Jrho_PC.define(a_grids, 1, IntVect::Zero);
   m_shift_PC.define(a_grids, 1, IntVect::Zero);
   m_ellip_coeff.define(a_grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   m_ellip_coeff_mapped.define(a_grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   
   // state variables
   //
   m_rhoDen_cc.define(a_grids, 1, a_ghostVect);
   m_momDen_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_eneDen_cc.define(a_grids, 1+m_CGL, a_ghostVect);
   m_momDen_virtual_cc.define(a_grids, 1, a_ghostVect);
   m_magField_virtual_cc.define(a_grids, 1, a_ghostVect);
   //
   m_pressure.define(a_grids, 1+m_CGL, a_ghostVect);
   m_Bpressure.define(a_grids, 1, a_ghostVect);
   m_velocity.define(a_grids, SpaceDim, a_ghostVect);
   m_velocity_cf.define(a_grids, SpaceDim, a_ghostVect);
   m_velocity_virtual.define(a_grids, 1, a_ghostVect);
   //
   m_rhoFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mxFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_myFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mzFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_enFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_enParFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mvFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   //
   m_Cspeed_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_CspeedL_norm.define(a_grids, 1, a_ghostVect);
   m_CspeedR_norm.define(a_grids, 1, a_ghostVect);
   //
   m_rhoFlux_norm.define(a_grids, 1, 0*IntVect::Unit);
   m_mxFlux_norm.define(a_grids, 1, 0*a_ghostVect);
   m_myFlux_norm.define(a_grids, 1, 0*a_ghostVect);
   m_mzFlux_norm.define(a_grids, 1, 0*a_ghostVect);
   m_enFlux_norm.define(a_grids, 1, 0*IntVect::Unit);
   m_enParFlux_norm.define(a_grids, 1, 0*IntVect::Unit);
   m_momFlux_norm.define(a_grids, SpaceDim, 0*a_ghostVect);
   m_mvFlux_norm.define(a_grids, 1, 0*IntVect::Unit);
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
      
   m_E0_ce.define(a_grids, 1, a_ghostVect);
   m_JaJ0_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_Apar_cc.define(a_grids, 1, a_ghostVect);
   m_E0_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_curlE0_virtual.define(a_grids, 1, IntVect::Zero);
   m_curlE0_cf.define(a_grids, 1, IntVect::Zero);
   m_By_contr_cc.define(a_grids, 1, a_ghostVect);
   m_B_contr_cf.define(a_grids, 1, a_ghostVect);
   m_JaJdotE.define(a_grids, 1, a_ghostVect); // Jacobian*J\cdotE
   m_JaJcrossB.define(a_grids, SpaceDim, a_ghostVect); // Jacobian*J\timesB
   //
   m_EdgeBC_zeros.define(a_grids, 1, a_ghostVect);
   m_FluxBC_zeros.define(a_grids, 1, 0*a_ghostVect);
   m_rhoFluxBC_norm.define(a_grids, 1, 0*a_ghostVect);
   m_mxFluxBC_norm.define(a_grids, 1, 0*a_ghostVect);
   m_myFluxBC_norm.define(a_grids, 1, 0*a_ghostVect);
   m_mzFluxBC_norm.define(a_grids, 1, 0*a_ghostVect);
   m_enFluxBC_norm.define(a_grids, 1, 0*IntVect::Unit);
   m_enParFluxBC_norm.define(a_grids, 1, 0*IntVect::Unit);
   m_mvFluxBC_norm.define(a_grids, 1, 0*IntVect::Unit);
   m_momFluxBC_norm.define(a_grids, SpaceDim, 0*a_ghostVect);
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_EdgeBC_zeros[dit].setVal(0.0);
      m_FluxBC_zeros[dit].setVal(0.0);
   }
   
   m_dummyDiv.define(a_grids, 1, IntVect::Zero);
   m_dummyFArray_oneComp.define(a_grids, 1, a_ghostVect);
   m_dummyFArray_spaceDim.define(a_grids, SpaceDim, a_ghostVect);
   m_dummyDiv_mom.define(a_grids, SpaceDim, IntVect::Zero);
   m_dummyFlux_oneComp.define(a_grids, 1, a_ghostVect);
   m_dummyFlux_spaceDim.define(a_grids, SpaceDim, a_ghostVect);
   m_dummyEdge_oneComp.define(a_grids, 1, a_ghostVect);

}

void IdealMhdOp::setCellCenterValues( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setCellCenterValues()");
   
   const LevelData<FArrayBox>& soln_rhoDen( a_soln_fluid.cell_var(0) );
   const LevelData<FArrayBox>& soln_momDen( a_soln_fluid.cell_var("momentumDensity") );
   const LevelData<FArrayBox>& soln_eneDen( a_soln_fluid.cell_var("energyDensity") );
   
   // set values for MHD variables
   //
   const DisjointBoxLayout& grids( m_geometry.grids() );
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
   
   // set physical derived variables for fluid species
   //
   a_soln_fluid.pressure(m_pressure);
   a_soln_fluid.Bpressure(m_Bpressure);
   a_soln_fluid.velocity(m_velocity);                 // in-plane velocity vector
   a_soln_fluid.velocity_virtual(m_velocity_virtual); // out-of-plane velocity vector
  
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
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
     
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

void IdealMhdOp::setCellCenterFluxesCGL( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setCellCenterFluxesCGL()");
   
   // Compute cell-center fluxes for each variable when using 2-pressure CGL model
   // (assumes B is in third direction or virtual direction for 2D)
   //
   // rhoFlux_x = mx,             rhoFlux_y = my,             rhoFlux_z = mz
   // mxFlux_x  = mx*ux + Pperp,  mxFlux_y  = mx*uy,          mxFlux_z  = mx*uz
   // myFlux_x  = my*ux,          myFlux_y  = my*uy + Pperp,  myFlux_z  = my*uz
   // mzFlux_x  = mz*ux,          mzFlux_y  = mz*uy,          mzFlux_z  = mz*uz + Ppar
   // enPerpFlux_x = (enPerp+Pperp)*ux,  enPerpFlux_y = (enPerp+Pperp)*uy,  enPerpFlux_z = enPerp*uz
   // enParFlux_x  = enPar*ux,           enParFlux_y  = enPar*uy,           enParFlux_z  = (enPar+Ppar)*uz
   // 
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
   
      if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
         for (int n=0; n<SpaceDim; ++n) {
            m_mvFlux_cc[dit].copy( m_momDen_virtual_cc[dit],0,n,1 );
            m_mvFlux_cc[dit].mult( m_velocity[dit],n,n,1);
         }
      }
      
      const FArrayBox& eneDen_on_patch = m_eneDen_cc[dit];  // has two comps now
      const FArrayBox& momDen_on_patch = m_momDen_cc[dit];   
      const FArrayBox& V_on_patch = m_velocity[dit];   
      const FArrayBox& P_on_patch = m_pressure[dit];  // has two comps now 
      const Box& cellbox = eneDen_on_patch.box();

      FArrayBox& rhoFlux_on_patch = m_rhoFlux_cc[dit];
      FArrayBox& mxFlux_on_patch  = m_mxFlux_cc[dit];
      FArrayBox& myFlux_on_patch  = m_myFlux_cc[dit];
      FArrayBox& mzFlux_on_patch  = m_mzFlux_cc[dit];
      FArrayBox& enPerpFlux_on_patch = m_enFlux_cc[dit];    // perp energy flux
      FArrayBox& enParFlux_on_patch  = m_enParFlux_cc[dit]; // par energy flux

      FORT_CELL_CENTER_FLUXES_CGL( CHF_BOX(cellbox), 
                                   CHF_CONST_FRA1(eneDen_on_patch,0),
                                   CHF_CONST_FRA1(eneDen_on_patch,1),
                                   CHF_CONST_FRA(momDen_on_patch),
                                   CHF_CONST_FRA(V_on_patch),
                                   CHF_CONST_FRA1(P_on_patch,0),
                                   CHF_CONST_FRA1(P_on_patch,1),
                                   CHF_FRA(rhoFlux_on_patch),
                                   CHF_FRA(mxFlux_on_patch),
                                   CHF_FRA(myFlux_on_patch),
                                   CHF_FRA(mzFlux_on_patch),
                                   CHF_FRA(enPerpFlux_on_patch),
                                   CHF_FRA(enParFlux_on_patch) );
      
   }
      
}

void IdealMhdOp::setMagneticFieldTerms( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setMagneticFieldTerms()");
   
   //  NEEDED FOR MAGNETIC FIELD ON CELL FACES AND ELECTRIC FIELD
   //  AT CELL EDGES
   //
   //
   const DisjointBoxLayout& grids( m_geometry.grids() );
   if(a_soln_fluid.m_evolve_magneticField_virtual) {
      
      const LevelData<FArrayBox>& magField_virtual_cc( a_soln_fluid.cell_var("magneticField_virtual") );

      // compute physical E=BxV
      //
      computePhysicalElectricField( m_E0_cc, m_velocity, magField_virtual_cc, 1 );

      // compute Ja*curl(B)_phys (aka mapped)
      // 
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& thisbox = m_By_contr_cc[dit].box();
         m_By_contr_cc[dit].copy( magField_virtual_cc[dit], thisbox); // phys
         if(m_twoDaxisymm) {
            m_By_contr_cc[dit].divide(m_g_y[dit], thisbox, 0,0,1);
         }
      }
      computeMappedCurrentDensity(m_JaJ0_cc,m_By_contr_cc);

      // calculate source terms for energy density and momentum density equations
      //
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& thisbox = m_JaJ0_cc[dit].box();
         FArrayBox m_JaJdotE_comp(thisbox, SpaceDim);
         m_JaJdotE_comp.copy(m_JaJ0_cc[dit], thisbox);
         m_JaJdotE_comp.mult(m_E0_cc[dit]);
         m_JaJdotE[dit].copy(m_JaJdotE_comp,0,0,1);
         m_JaJdotE[dit].plus(m_JaJdotE_comp,thisbox,1,0,1);
         //
         m_JaJcrossB[dit].copy(m_JaJ0_cc[dit],1,0,1);
         m_JaJcrossB[dit].mult(magField_virtual_cc[dit],thisbox,0,0,1);
         m_JaJcrossB[dit].negate(0,1);
         m_JaJcrossB[dit].copy(m_JaJ0_cc[dit],0,1,1);
         m_JaJcrossB[dit].mult(magField_virtual_cc[dit],thisbox,0,1,1);
      }
   }

   // if using in-plane magnetic field
   //
   if(a_soln_fluid.m_evolve_magneticField) {
      CH_assert(!a_soln_fluid.m_evolve_magneticField_virtual); // don't work together yet !!!
      const LevelData<FluxBox>& soln_magField_cf( a_soln_fluid.face_var("magneticField") );
      CH_assert(SpaceDim==3); // only works in 2D for now !!!
 
      // compute Ja*curl(B)_phys (aka mapped)
      //
      computeMappedCurrentDensity(m_JaJ0_cc,soln_magField_cf);
  
      // compute physical B at cell center
      //
      LevelData<FArrayBox>& B_phys_cc = m_dummyFArray_spaceDim;
      SpaceUtils::interpFaceVectorToCell(B_phys_cc,soln_magField_cf,"c2");
      m_geometry.convertPhysToContravar(B_phys_cc,1);
      
      // compute physical E=BxV at cell center
      //
      computePhysicalElectricField( m_E0_cc, m_velocity, B_phys_cc, 0 );

      // compute Ja*JxB and Ja*JdotE at cell center 
      //
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         
         const FArrayBox& JaJ_on_patch     = m_JaJ0_cc[dit];
         const FArrayBox& E_on_patch       = m_E0_cc[dit];
         const FArrayBox& B_on_patch       = B_phys_cc[dit];
               FArrayBox& JaJdotE_on_patch = m_JaJdotE[dit];
               FArrayBox& JaJxB_on_patch = m_JaJcrossB[dit];

         const Box& cellbox = JaJdotE_on_patch.box();
         FORT_EVAL_DOT_PRODUCT( CHF_BOX(cellbox),
                                CHF_CONST_FRA(JaJ_on_patch),
                                CHF_CONST_FRA(E_on_patch),
                                CHF_FRA1(JaJdotE_on_patch,0) );
         
         FORT_EVAL_CROSS_PRODUCT( CHF_BOX(cellbox),
                                  CHF_CONST_FRA(JaJ_on_patch),
                                  CHF_CONST_FRA(B_on_patch),
                                  CHF_FRA(JaJxB_on_patch) );

         // reset dummy containers to infty
         //
         m_dummyFArray_spaceDim[dit].setVal(1./0.);
      }
      
   }

}   

void IdealMhdOp::setFaceCenteredFluxes( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setFaceCenteredFluxes()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // Get mapped flux normals on each face multiplied by face area (NTFlux*faceArea)
   //
   LevelData<FArrayBox>& thisJf = m_dummyFArray_oneComp;
   LevelData<FArrayBox>& thisFlux_mapped = m_dummyFArray_spaceDim;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      thisJf[dit].copy(m_rhoDen_cc[dit]);
      thisJf[dit].mult(m_Jacobian[dit],0,0,1);
   }
   m_geometry.multiplyNTranspose(thisFlux_mapped, m_rhoFlux_cc); // convert phys flux to mapped at cc
   computeNTFfaceArea(m_rhoFlux_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   
   for (int n=0; n<SpaceDim; ++n) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         thisJf[dit].copy(m_momDen_cc[dit],n,0,1);
         thisJf[dit].mult(m_Jacobian[dit],0,0,1);
         if(n==0) m_geometry.multiplyNTranspose(thisFlux_mapped, m_mxFlux_cc);
         if(n==1) m_geometry.multiplyNTranspose(thisFlux_mapped, m_myFlux_cc);
         if(n==2) m_geometry.multiplyNTranspose(thisFlux_mapped, m_mzFlux_cc);
      }
      if(n==0) computeNTFfaceArea(m_mxFlux_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
      if(n==1) computeNTFfaceArea(m_myFlux_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
      if(n==2) computeNTFfaceArea(m_mzFlux_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   }
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_momFlux_norm[dit].copy( m_mxFlux_norm[dit],0,0,1 );
      m_momFlux_norm[dit].copy( m_myFlux_norm[dit],0,1,1 );
      if(SpaceDim==3) m_momFlux_norm[dit].copy( m_mzFlux_norm[dit],0,2,1 );
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      thisJf[dit].copy(m_eneDen_cc[dit],0,0,1);
      thisJf[dit].mult(m_Jacobian[dit],0,0,1);
   }
   m_geometry.multiplyNTranspose(thisFlux_mapped, m_enFlux_cc); // convert phys flux to mapped at cc
   computeNTFfaceArea(m_enFlux_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   
   if(m_CGL) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         thisJf[dit].copy(m_eneDen_cc[dit],1,0,1);
         thisJf[dit].mult(m_Jacobian[dit],0,0,1);
      }
      m_geometry.multiplyNTranspose(thisFlux_mapped, m_enParFlux_cc); // convert phys flux to mapped at cc
      computeNTFfaceArea(m_enParFlux_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   }

   if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         thisJf[dit].copy(m_momDen_virtual_cc[dit]);
         thisJf[dit].mult(m_Jacobian[dit],0,0,1);
      }
      m_geometry.multiplyNTranspose(thisFlux_mapped, m_mvFlux_cc); // convert phys flux to mapped at cc
      computeNTFfaceArea(m_mvFlux_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   }
     
   // set m_dummyFArray containers to infinity at end of function for safety
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_oneComp[dit].setVal(1.0/0.0);
      m_dummyFArray_spaceDim[dit].setVal(1.0/0.0);
   }

}

void IdealMhdOp::setFaceCenteredFlux( LevelData<FluxBox>&    a_Flux_norm,
                                const LevelData<FArrayBox>&  a_Flux_phys_cc )
{
   CH_TIME("IdealMhdOp::setFaceCenteredFlux()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   const int method = 2;

   // original way
   //
   if(method==0) {
      SpaceUtils::upWindToFaces( m_dummyFlux_spaceDim, a_Flux_phys_cc, m_dummyFlux_spaceDim, "c2" );
      m_geometry.applyAxisymmetricCorrection( m_dummyFlux_spaceDim );
      m_geometry.computeMetricTermProductAverage( a_Flux_norm, m_dummyFlux_spaceDim, 0 );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_dummyFlux_spaceDim[dit].setVal(1.0/0.0);
      }
   }
   
   // convert cell-center flux from physical to mapped and interp to faces 
   //
   if(method==1) {
      m_geometry.multiplyNTranspose(m_dummyFArray_spaceDim, a_Flux_phys_cc);
      SpaceUtils::upWindToFaces( a_Flux_norm, m_dummyFArray_spaceDim, a_Flux_norm, "c2" );
   
      for (DataIterator dit(grids); dit.ok(); ++dit) { // mult by face area
         const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         const RealVect& faceArea = coord_sys.getMappedFaceArea();
         for (int dir=0; dir<SpaceDim; dir++) {
            FArrayBox& this_face = a_Flux_norm[dit][dir];
            this_face.mult(faceArea[dir]);
         }
         m_dummyFArray_spaceDim[dit].setVal(1.0/0.0); // reset dummy to infinty for safety
      }
   }
   
   // new method to be consistent with how symmetry BCs are applied
   //
   if(method==2) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_dummyFArray_spaceDim[dit].copy(a_Flux_phys_cc[dit]);
      }
      m_geometry.convertPhysToContravar(m_dummyFArray_spaceDim, 0);
      SpaceUtils::upWindToFaces( a_Flux_norm, m_dummyFArray_spaceDim, a_Flux_norm, "c2" );
      m_geometry.multJonFaces(a_Flux_norm);  
 
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         const RealVect& faceArea = coord_sys.getMappedFaceArea();
         for (int dir=0; dir<SpaceDim; dir++) {
            FArrayBox& this_face = a_Flux_norm[dit][dir];
            this_face.mult(faceArea[dir]);
         }
         m_dummyFArray_spaceDim[dit].setVal(1.0/0.0); // reset dummy for safety
      }
   }

}

void IdealMhdOp::enforceFluxBCs( const FluidSpecies&  a_soln_fluid, 
                                 const Real           a_time )
{
   CH_TIME("IdealMhdOp::enforceFluxBCs()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // enforce flux BC for mass density flux
   //
   setFaceCenteredFlux(m_rhoFluxBC_norm, m_rhoFlux_cc);
   m_fluid_bc.at(0)->setFluxBC( a_soln_fluid, m_rhoFlux_norm, m_rhoFluxBC_norm, a_time );
   //m_fluid_bc.at(0)->setFluxBC( a_soln_fluid, m_rhoFlux_norm, m_FluxBC_zeros, a_time );
   
   // enforce flux BC for momentum density flux
   //
   int this_cbc;
   setFaceCenteredFlux(m_mxFluxBC_norm, m_mxFlux_cc);
   setFaceCenteredFlux(m_myFluxBC_norm, m_myFlux_cc);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_momFluxBC_norm[dit].copy( m_mxFluxBC_norm[dit],0,0,1 );
      m_momFluxBC_norm[dit].copy( m_myFluxBC_norm[dit],0,1,1 );
   }
   if(SpaceDim==3) {
      setFaceCenteredFlux(m_mzFluxBC_norm, m_mzFlux_cc);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_momFluxBC_norm[dit].copy( m_mzFluxBC_norm[dit],0,2,1 );
      }
   }
   this_cbc = m_momentumDensity_cbc;
   m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_momFlux_norm, m_momFluxBC_norm, a_time );
  
   // enforce flux BC for energy density flux
   //
   setFaceCenteredFlux(m_enFluxBC_norm, m_enFlux_cc);
   this_cbc = m_energyDensity_cbc;
   m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_enFlux_norm, m_enFluxBC_norm, a_time );
   //m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_enFlux_norm, m_FluxBC_zeros, a_time );
   
   // enforce flux BC for parallel energy density flux if using CGL
   //
   if(m_CGL) {
      setFaceCenteredFlux(m_enParFluxBC_norm, m_enParFlux_cc);
      this_cbc = m_energyDensity_cbc;
      m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_enParFlux_norm, m_enParFluxBC_norm, a_time );
   }

   // enforce flux BC for virtual momentum density flux
   //
   if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
      setFaceCenteredFlux(m_mvFluxBC_norm, m_mvFlux_cc);
      this_cbc = m_momentumDensity_virtual_cbc;
      m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_mvFlux_norm, m_mvFluxBC_norm, a_time );
   }
   
   // enforce edge BC for virtual magnetic field
   //
   if(a_soln_fluid.m_evolve_magneticField_virtual) {
      computeIdealEatEdges(m_E0_ce,m_velocity,m_By_contr_cc,m_Cspeed_cc);
      this_cbc = m_magneticField_virtual_cbc;
      m_fluid_bc.at(this_cbc)->setEdgeBC( a_soln_fluid, m_E0_ce, m_EdgeBC_zeros, a_time );
      m_geometry.mappedGridCurlFromEdgeTans(m_E0_ce,2,m_curlE0_virtual); // takes covar E and returns Ja*curl(E)dot(g^y)
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         if(m_twoDaxisymm) {
            m_curlE0_virtual[dit].mult(m_g_y[dit]); // g^y*g_y = 1 (Ja*curlE_phys)
         }
      }
   }
   
   // compute E0 (should do this elsewhere) and enforce boundary conditions
   //
   if(a_soln_fluid.m_evolve_magneticField) {
      const LevelData<FluxBox>& soln_magField_cf( a_soln_fluid.face_var("magneticField") );
      computeIdealEatEdges(m_E0_ce,m_velocity,soln_magField_cf,m_Cspeed_cc);
      int this_fbc = m_magneticField_fbc;
      m_fluid_face_var_bc.at(this_fbc)->setEdgeBC( a_soln_fluid, m_E0_ce, m_EdgeBC_zeros, a_time );
      SpaceUtils::exchangeEdgeDataBox(m_E0_ce); // call exchange after setBC !!!!!
      m_geometry.mappedGridCurl3D(m_curlE0_cf,m_E0_ce);
   }

}   

void IdealMhdOp::setViscosityCoefficient( LevelData<FArrayBox>&  a_etaVisc,
                                          LevelData<FluxBox>&    a_etaVisc_cf,
                                    const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setViscosityCoefficient()");

   // set the ion viscosity coefficients
   // eta = 0.96*tau*Pi + etaMin
   //
   const double factor = 0.96*0.5; // 1/2 for only ion pressure
   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_etaVisc[dit].copy(m_pressure[dit],0,0,1);
      a_etaVisc[dit].mult(factor*m_tau);
      a_etaVisc[dit].plus(m_etaMin);
   }
   SpaceUtils::upWindToFaces(a_etaVisc_cf,a_etaVisc,a_etaVisc_cf,"c2"); 
   

   // Create elliptic coefficients
   //
   /*
   LevelData<FluxBox> ellip_coeff(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   LevelData<FluxBox> ellip_coeff_mapped(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   const double a_dummy_time=0.0;
   computeDiffusionCoefficients(ellip_coeff, ellip_coeff_mapped, a_etaVisc_cf, a_dummy_time);
   m_diffusionOp_mom0->setOperatorCoefficients(ellip_coeff, ellip_coeff_mapped, *m_bc_mom0);
   m_diffusionOp_mom1->setOperatorCoefficients(ellip_coeff, ellip_coeff_mapped, *m_bc_mom1);
   //SpaceUtils::inspectFluxBox(ellip_coeff,1);
   //m_geometry.plotFaceData( "ellip_coeff_mapped", ellip_coeff_mapped, a_dummy_time);
   //SpaceUtils::inspectFluxBox(ellip_coeff_mapped,1);
   */

}

void IdealMhdOp::computeViscSourceAxisymm( LevelData<FArrayBox>&  a_Wthth,
                                     const LevelData<FArrayBox>&  a_velocity_phys,
                                     const LevelData<FluxBox>&    a_velocity_phys_cf )
{
   CH_TIME("IdealMhdOp::computeViscSourceAxisymm()");

   // calculate div stress tensor source terms at cell center
   // Wthth = 2*Vr/r - 2/3*divV
   // Wthr  = d(Vth)/dr + d(Vr)/dtheta/r
   //
   
   //   compute divergence of physical velocity field
   //
   const DisjointBoxLayout& grids( m_geometry.grids() );
   LevelData<FluxBox> velocity_norm(grids, 1, IntVect::Unit);
   //m_geometry.applyAxisymmetricCorrection( a_velocity_phys_cf );
   m_geometry.computeMetricTermProductAverage( velocity_norm, a_velocity_phys_cf, 0 );
   m_geometry.mappedGridDivergenceFromFluxNorms( velocity_norm, m_dummyDiv );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyDiv[dit].divide(m_cellVol[dit]);
   }
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_Wthth[dit].copy(a_velocity_phys[dit],0,0,1);
      a_Wthth[dit].divide(m_Xphys[dit],0,0,1);
      a_Wthth[dit].mult(2.0);
      m_dummyDiv[dit].mult(-2.0/3.0);
      a_Wthth[dit].plus(m_dummyDiv[dit],0,0,1);
   }

}

void IdealMhdOp::computeViscosityEnergyFlux( LevelData<FluxBox>&    a_UJaPi_Flux_phys,
                                       const LevelData<FluxBox>&    a_velocity_phys_cf )
{
   CH_TIME("IdealMhdOp::computeViscosityEnergyFlux()");

   // calculate physical viscosity fluxes for total energy equation
   // on each face (multiplied by Jacobian)
   //
   // Ja*Flux_0 = U0*JaPi00 + U1*JaPi10 + U2*JaPi20
   // Ja*Flux_1 = U0*JaPi01 + U1*JaPi11 + U2*JaPi21
   // Ja*Flux_2 = U0*JaPi02 + U1*JaPi12 + U2*JaPi22
   //

   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; ++dir) {

         const FArrayBox& vel_on_dir   = a_velocity_phys_cf[dit][dir];
         const FArrayBox& JaPi0_on_dir = m_m0JaFluxVisc_cf[dit][dir];
               FArrayBox& UJaPi_on_dir = a_UJaPi_Flux_phys[dit][dir];

         const Box& thisbox = UJaPi_on_dir.box();
         FORT_EVAL_DOT_PRODUCT( CHF_BOX(thisbox),
                                CHF_CONST_FRA(vel_on_dir),
                                CHF_CONST_FRA(JaPi0_on_dir),
                                CHF_FRA1(UJaPi_on_dir,0) );

         const FArrayBox& JaPi1_on_dir = m_m1JaFluxVisc_cf[dit][dir];
         FORT_EVAL_DOT_PRODUCT( CHF_BOX(thisbox),
                                CHF_CONST_FRA(vel_on_dir),
                                CHF_CONST_FRA(JaPi1_on_dir),
                                CHF_FRA1(UJaPi_on_dir,1) );
         if(SpaceDim==3) {
            const FArrayBox& JaPi2_on_dir = m_m2JaFluxVisc_cf[dit][dir];
            FORT_EVAL_DOT_PRODUCT( CHF_BOX(thisbox),
                                   CHF_CONST_FRA(vel_on_dir),
                                   CHF_CONST_FRA(JaPi2_on_dir),
                                   CHF_FRA1(UJaPi_on_dir,2) );
         }

      }

   }

}

void IdealMhdOp::updateRHSs_visc( FluidSpecies&  a_rhs_fluid,
                            const FluidSpecies&  a_soln_fluid,
                            const int            a_accumulateMom,
                            const int            a_accumulateEne )
{
   CH_TIME("IdealMhdOp::updateRHSs_visc()");

   // compute physical viscosity fluxes on cell faces
   // for each momentum equation (multiplied by Jacobian)
   // where viscosity tensor is Pi = -etaVis*W with 
   // rate of strain tensor W defined as
   // W = \nabla(V) + (\nabla(V))^T - 2/3*I\nabla\cdot V
   //
   m_geometry.computeJaStrainTensorPhys( m_JaW_cf,
                                         m_velocity );
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int n=0; n<SpaceDim; ++n) {
         m_m0JaFluxVisc_cf[dit].copy( m_JaW_cf[dit],n,n,1 ); 
         m_m1JaFluxVisc_cf[dit].copy( m_JaW_cf[dit],SpaceDim+n,n,1 );
         if(SpaceDim==3) m_m2JaFluxVisc_cf[dit].copy( m_JaW_cf[dit],2*SpaceDim+n,n,1 );
      } 
      const Box& cell_box( m_etaVisc_cf[dit].box() );
      for (int n=0; n<SpaceDim; ++n) {
         m_m0JaFluxVisc_cf[dit].mult( m_etaVisc_cf[dit],cell_box,0,n,1 ); 
         m_m1JaFluxVisc_cf[dit].mult( m_etaVisc_cf[dit],cell_box,0,n,1 ); 
         if(SpaceDim==3) m_m2JaFluxVisc_cf[dit].mult( m_etaVisc_cf[dit],cell_box,0,n,1 ); 
      }
      m_m0JaFluxVisc_cf[dit].negate(); 
      m_m1JaFluxVisc_cf[dit].negate(); 
      if(SpaceDim==3) m_m2JaFluxVisc_cf[dit].negate(); 
   }

   // calculate viscosity flux for energy density equation
   // and source term for momentum flux for axisymmetric
   //
   SpaceUtils::upWindToFaces(m_velocity_cf, m_velocity, m_velocity_cf, "c2");
   if(m_twoDaxisymm) {
      computeViscSourceAxisymm( m_momVisc_source, m_velocity, m_velocity_cf );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_momVisc_source[dit].mult( m_etaVisc[dit] ); 
         m_momVisc_source[dit].negate(); 
      }
   } 
   
   computeViscosityEnergyFlux( m_enJaFluxVisc_cf,
                               m_velocity_cf );

   //   update RHSs for momentum equation
   //
   if(a_accumulateMom) {

      // convert physical Ja*Flux to contravar and mult by faceArea
      // (Result is same as applying N^T*Flux*faceArea as is done by 
      // computeMetricTermProductAverage), which is not done here 
      // since Flux is already multiplied by Jacobian
      //
      m_geometry.applyAxisymmetricCorrection( m_m0JaFluxVisc_cf );
      m_geometry.applyAxisymmetricCorrection( m_m1JaFluxVisc_cf );
      if(SpaceDim==3) m_geometry.applyAxisymmetricCorrection( m_m2JaFluxVisc_cf );
      m_geometry.computedxidXProductAverage(m_m0JaFluxVisc_norm, m_m0JaFluxVisc_cf, 0);
      m_geometry.computedxidXProductAverage(m_m1JaFluxVisc_norm, m_m1JaFluxVisc_cf, 0);
      if(SpaceDim==3) m_geometry.computedxidXProductAverage(m_m2JaFluxVisc_norm, m_m2JaFluxVisc_cf, 0);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_momJaFluxVisc_norm[dit].copy( m_m0JaFluxVisc_norm[dit],0,0,1 ); 
         m_momJaFluxVisc_norm[dit].copy( m_m1JaFluxVisc_norm[dit],0,1,1 );
         if(SpaceDim==3) m_momJaFluxVisc_norm[dit].copy( m_m2JaFluxVisc_norm[dit],0,2,1 );
      }
   
      LevelData<FArrayBox>& rhs_mom( a_rhs_fluid.cell_var("momentumDensity") );
      m_geometry.mappedGridDivergenceFromFluxNorms(m_momJaFluxVisc_norm, m_dummyDiv_mom);
      for (DataIterator dit(rhs_mom.dataIterator()); dit.ok(); ++dit) {
         for (int n=0; n<SpaceDim; ++n) {
            m_dummyDiv_mom[dit].divide(m_cellVol[dit],0,n,1);
            m_dummyDiv_mom[dit].mult(m_Jacobian[dit],0,n,1);
            if(n==0 && m_twoDaxisymm) {
               FArrayBox JaPiyyOverR(grids[dit], 1);
               JaPiyyOverR.copy(m_momVisc_source[dit]);
               JaPiyyOverR.mult(m_Jacobian[dit]); // Jacobian = Ja = 2*pi*R_phys * J_RZ
               JaPiyyOverR.divide(m_Xphys[dit],0,0,1);
               m_dummyDiv_mom[dit].minus(JaPiyyOverR,0,0,1);
            }
         }
         //rhs_mom[dit].minus(m_dummyDiv_mom[dit],1,1,1);
         rhs_mom[dit].minus(m_dummyDiv_mom[dit]);
      }
   
   }   

   //   update RHSs for energy density equation
   //
   if(a_accumulateEne) {

      // convert physical Ja*Flux to contravar and mult by faceArea
      // (Result is same as applying N^T*Flux*faceArea as is done by 
      // computeMetricTermProductAverage), which is not done here 
      // since Flux is already multiplied by Jacobian
      //
      m_geometry.applyAxisymmetricCorrection( m_enJaFluxVisc_cf );
      m_geometry.computedxidXProductAverage(m_enJaFluxVisc_norm, m_enJaFluxVisc_cf, 0);

      LevelData<FArrayBox>& rhs_ene( a_rhs_fluid.cell_var("energyDensity") );
      m_geometry.mappedGridDivergenceFromFluxNorms(m_enJaFluxVisc_norm, m_dummyDiv);
      for (DataIterator dit(rhs_ene.dataIterator()); dit.ok(); ++dit) {
         m_dummyDiv[dit].divide(m_cellVol[dit]);
         m_dummyDiv[dit].mult(m_Jacobian[dit]);
         rhs_ene[dit].minus( m_dummyDiv[dit],0,0,1 ); // specify comp since can have two (CGL?)
      }

   }

}


void IdealMhdOp::updateRHSs( FluidSpecies&  a_rhs_fluid,
                       const FluidSpecies&  a_soln_fluid ) 
{
   CH_TIME("IdealMhdOp::updateRHSs()");
   
   // compute divergence of Fluxes and add (subtract) to RHS 
   // Note that Ja/cellVol = 1/mapVol
   //

   // update RHS for mass density
   //
   LevelData<FArrayBox>& rhs_rho( a_rhs_fluid.cell_var(0) );
   m_geometry.mappedGridDivergenceFromFluxNorms(m_rhoFlux_norm, m_dummyDiv);
   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (DataIterator dit(rhs_rho.dataIterator()); dit.ok(); ++dit) {
      m_dummyDiv[dit].divide(m_cellVol[dit]);
      m_dummyDiv[dit].mult(m_Jacobian[dit]);
      rhs_rho[dit].minus(m_dummyDiv[dit]);
   }
  
   // update RHS for momentum density
   // 
   LevelData<FArrayBox>& rhs_mom( a_rhs_fluid.cell_var("momentumDensity") );
   m_geometry.mappedGridDivergenceFromFluxNorms(m_momFlux_norm, m_dummyDiv_mom);
   for (DataIterator dit(rhs_mom.dataIterator()); dit.ok(); ++dit) {
      for (int n=0; n<SpaceDim; ++n) {
         if(n==0 && m_twoDaxisymm) {
            //if(a_soln_fluid.m_evolve_magneticField_virtual) {
            //   m_pressure[dit].minus(m_Bpressure[dit],0,0,1);
            //}
            FArrayBox pdivIdent_R(grids[dit],1);
            FArrayBox pdivIdent_Z(grids[dit],1);
            pdivIdent_R.copy( m_div_IdentFlux_R[dit] );
            pdivIdent_Z.copy( m_div_IdentFlux_Z[dit] );
            pdivIdent_R.mult(m_pressure[dit],0+m_CGL,0,1);
            pdivIdent_Z.mult(m_pressure[dit],0+m_CGL,0,1);
            m_dummyDiv_mom[dit].minus(pdivIdent_R,0,0,1);
            m_dummyDiv_mom[dit].minus(pdivIdent_Z,0,1,1);
         }
         m_dummyDiv_mom[dit].divide(m_cellVol[dit],0,n,1);
         m_dummyDiv_mom[dit].mult(m_Jacobian[dit],0,n,1);
      }
      //m_dummyDiv_mom[dit].setVal(0.0,0); // JRA, no mr update
      rhs_mom[dit].minus(m_dummyDiv_mom[dit]);
      if(m_twoDaxisymm && a_soln_fluid.m_evolve_momentumDensity_virtual) {
         FArrayBox JaRhoVthsqOverR(grids[dit], 1);
         JaRhoVthsqOverR.copy(m_momDen_virtual_cc[dit]);
         JaRhoVthsqOverR.mult(m_velocity_virtual[dit]);
         JaRhoVthsqOverR.mult(m_Jacobian[dit]); // Jacobian = Ja = 2*pi*R_phys * J_RZ
         JaRhoVthsqOverR.divide(m_Xphys[dit],0,0,1);
         rhs_mom[dit].plus(JaRhoVthsqOverR,0,0,1);            
      }
      if(a_soln_fluid.m_evolve_magneticField_virtual ||
         a_soln_fluid.m_evolve_magneticField) {
         rhs_mom[dit].plus(m_JaJcrossB[dit]); // JRA !!!!!!!!!!!!!!!!
      }
   }
   
   // update RHS for energy density
   //
   LevelData<FArrayBox>& rhs_ene( a_rhs_fluid.cell_var("energyDensity") );
   m_geometry.mappedGridDivergenceFromFluxNorms(m_enFlux_norm, m_dummyDiv);
   for (DataIterator dit(rhs_ene.dataIterator()); dit.ok(); ++dit) {
      m_dummyDiv[dit].divide(m_cellVol[dit]);
      m_dummyDiv[dit].mult(m_Jacobian[dit]);
      rhs_ene[dit].minus( m_dummyDiv[dit],0,0,1 );
      if(a_soln_fluid.m_evolve_magneticField_virtual ||
         a_soln_fluid.m_evolve_magneticField) {
         rhs_ene[dit].plus( m_JaJdotE[dit],0,0,1 ); //JRA !!!!!!!!!!!!!!!!
      }
   }
   
   // update RHS for parallel energy density if using CGL model
   //
   if(m_CGL) {
      m_geometry.mappedGridDivergenceFromFluxNorms(m_enParFlux_norm, m_dummyDiv);
      for (DataIterator dit(rhs_ene.dataIterator()); dit.ok(); ++dit) {
         m_dummyDiv[dit].divide(m_cellVol[dit]);
         m_dummyDiv[dit].mult(m_Jacobian[dit]);
         rhs_ene[dit].minus( m_dummyDiv[dit],0,1,1 );
         if(m_twoDaxisymm) {
            m_dummyFArray_oneComp[dit].copy(m_velocity[dit],0,0,1);
            m_dummyFArray_oneComp[dit].mult(m_eneDen_cc[dit],1,0,1);
            m_dummyFArray_oneComp[dit].mult(m_Jacobian[dit],0,0,1);
            m_dummyFArray_oneComp[dit].divide(m_Xphys[dit],0,0,1);
            m_dummyFArray_oneComp[dit].mult(2.0); // Ja*2*Ur*Epar/r
            rhs_ene[dit].plus(m_dummyFArray_oneComp[dit],0,0,1);            
            rhs_ene[dit].minus(m_dummyFArray_oneComp[dit],0,1,1);
         }            
      }
   }
   
   // set m_dummyFArray containers to infinity after use for safety
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_oneComp[dit].setVal(1.0/0.0);
   }
 
   /*
   LevelData<FluxBox> constFlux(grids, 1, 0*IntVect::Unit); 
   LevelData<FArrayBox> div_constFlux(grids, 1, 0*IntVect::Unit); 
   for (DataIterator dit(rhs_data3.dataIterator()); dit.ok(); ++dit) {
      constFlux[dit].setVal(1.0);
   }
   m_geometry.multJonFaces(constFlux);
   m_geometry.mappedGridDivergenceFromFluxNorms(constFlux, div_constFlux);
   */

   // update RHS for momentum density in virtual direction
   //
   if(a_soln_fluid.m_evolve_momentumDensity_virtual){
      LevelData<FArrayBox>& rhs_momV( a_rhs_fluid.cell_var("momentumDensity_virtual") );
      m_geometry.mappedGridDivergenceFromFluxNorms(m_mvFlux_norm, m_dummyDiv);
      for (DataIterator dit(rhs_momV.dataIterator()); dit.ok(); ++dit) {
         FArrayBox& div_mvFlux(m_dummyDiv[dit]);
         div_mvFlux /= m_cellVol[dit];
         div_mvFlux.mult(m_Jacobian[dit]);
         rhs_momV[dit].minus( div_mvFlux );
         if(m_twoDaxisymm) {
            FArrayBox JaRhoVthVrOverR(grids[dit], 1);
            JaRhoVthVrOverR.copy(m_momDen_virtual_cc[dit]);
            JaRhoVthVrOverR.mult(m_velocity[dit],0,0,1);
            JaRhoVthVrOverR.mult(m_Jacobian[dit]); // Jacobian = Ja = 2*pi*R_phys * J_RZ
            JaRhoVthVrOverR.divide(m_Xphys[dit],0,0,1);
            rhs_momV[dit].minus(JaRhoVthVrOverR,0,0,1);            
         }
      }
   }

   // update RHS for virtual magnetic field
   //
   if(a_soln_fluid.m_evolve_magneticField_virtual){
      LevelData<FArrayBox>& rhs_magV( a_rhs_fluid.cell_var("magneticField_virtual") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs_magV[dit].minus( m_curlE0_virtual[dit] );
      }
   }
   
   // update RHS for magnetic field
   //
   if(a_soln_fluid.m_evolve_magneticField){
      LevelData<FluxBox>& rhs_magField( a_rhs_fluid.face_var("magneticField") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs_magField[dit].minus( m_curlE0_cf[dit],rhs_magField[dit].box(),0,0,1 );
      }
   }

}   

void IdealMhdOp::initialize( CFGVars&    a_species, 
                             const Real  a_time )
{
   CH_TIME("IdealMhdOp::initialize()");
  
   // initilize cell center variables
   //
   //m_CGL = a_species.m_CGL; // CFG var here rather than FluidSpecies...dynamic cast?
   for (int n=0; n<a_species.num_cell_vars(); ++n) {

      const int num_comps = a_species.cell_var(n).nComp();
      if (num_comps>1) {
        string this_cell_var_comp_name; 
        for (int m=0; m<num_comps; m++) {
            if(m==0) this_cell_var_comp_name = a_species.cell_var_name(n)+"_0";
            if(m==1) this_cell_var_comp_name = a_species.cell_var_name(n)+"_1";
            if(m==2) this_cell_var_comp_name = a_species.cell_var_name(n)+"_2";
            const GridFunction& fsic( fluidSpeciesIC( this_cell_var_comp_name ) );
            LevelData<FArrayBox> this_cell_var_comp( a_species.cell_var(n).getBoxes(), 1,
                                                     a_species.cell_var(n).ghostVect() );
            fsic.assign( this_cell_var_comp, a_species.configurationSpaceGeometry(), a_time );
            LevelData<FArrayBox>& cell_var_ref( a_species.cell_var(n) );
            for (DataIterator dit( this_cell_var_comp.dataIterator() ); dit.ok(); ++dit) {
               cell_var_ref[dit].copy(this_cell_var_comp[dit],0,m,1);
            }
         }
      }
      else {
         const GridFunction& fsic( fluidSpeciesIC( a_species.cell_var_name(n) ) );
         fsic.assign( a_species.cell_var(n), a_species.configurationSpaceGeometry(), a_time );
      }

   }
   
   // initilize face center variables
   //
   for (int n=0; n<a_species.num_face_vars(); ++n) { // initialize contravar on faces from phys on cell
      string this_face_var_comp_name;
      LevelData<FArrayBox> this_face_var_on_cell( a_species.face_var(n).getBoxes(), SpaceDim,
                                                  a_species.face_var(n).ghostVect() );
      for(int m=0; m<SpaceDim; ++m) {
         if(m==0) this_face_var_comp_name = a_species.face_var_name(n)+"_0";
         if(m==1) this_face_var_comp_name = a_species.face_var_name(n)+"_1";
         if(m==2) this_face_var_comp_name = a_species.face_var_name(n)+"_2";
         const GridFunction& fsic( fluidSpeciesIC( this_face_var_comp_name ) );
         LevelData<FArrayBox> this_face_var_on_cell_comp( a_species.face_var(n).getBoxes(), 1,
                                                          a_species.face_var(n).ghostVect() );
         fsic.assign( this_face_var_on_cell_comp, a_species.configurationSpaceGeometry(), a_time );
         for (DataIterator dit( this_face_var_on_cell_comp.dataIterator() ); dit.ok(); ++dit) {
            this_face_var_on_cell[dit].copy(this_face_var_on_cell_comp[dit],0,m,1);
         }
      }
      const MagGeom& geometry = a_species.configurationSpaceGeometry();
      geometry.convertPhysToContravar(this_face_var_on_cell,0);
      SpaceUtils::upWindToFaces(a_species.face_var(n),this_face_var_on_cell,a_species.face_var(n),"c2");
   }

   // initilize edge center variables
   // 
   for (int n=0; n<a_species.num_edge_vars(); ++n) { // initialize covar on edges from phys on cell
      string this_edge_var_comp_name;
      LevelData<FArrayBox> this_edge_var_on_cell( a_species.edge_var(n).getBoxes(), SpaceDim,
                                                  a_species.edge_var(n).ghostVect() );
      for(int m=0; m<SpaceDim; ++m) {
         if(m==0) this_edge_var_comp_name = a_species.edge_var_name(n)+"_0";
         if(m==1) this_edge_var_comp_name = a_species.edge_var_name(n)+"_1";
         if(m==2) this_edge_var_comp_name = a_species.edge_var_name(n)+"_2";
         const GridFunction& fsic( fluidSpeciesIC( this_edge_var_comp_name ) );
         LevelData<FArrayBox> this_edge_var_on_cell_comp( a_species.edge_var(n).getBoxes(), 1,
                                                          a_species.edge_var(n).ghostVect() );
         fsic.assign( this_edge_var_on_cell_comp, a_species.configurationSpaceGeometry(), a_time );
         for (DataIterator dit( this_edge_var_on_cell_comp.dataIterator() ); dit.ok(); ++dit) {
            this_edge_var_on_cell[dit].copy(this_edge_var_on_cell_comp[dit],0,m,1);
         }
      }
      const MagGeom& geometry = a_species.configurationSpaceGeometry();
      geometry.convertPhysToCovar(this_edge_var_on_cell,0);
      SpaceUtils::interpCellToEdges(a_species.edge_var(n),this_edge_var_on_cell,a_species.edge_var(n),"c2");
   }

   if ( typeid(a_species) == typeid(FluidSpecies&) ) {
      (static_cast<FluidSpecies&>(a_species)).convertFromPhysical(); // this just multiplies by Jacobian
   }

}

void IdealMhdOp::initializeWithBC( FluidSpecies&  a_species_comp,
                                   FluidSpecies&  a_species_phys,
                             const double         a_time )
{
   CH_TIME("IdealMhdOp::initializeWithBC()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   if(a_species_comp.m_evolve_magneticField && m_initializeBfromVectorPotential) {
         
      LevelData<FluxBox>& contravar_B_withBCs( a_species_phys.face_var("magneticField") );
      int this_fbc = m_magneticField_fbc;
      
      // compute Ja*curlB_phys at cell center using input B field
      computeMappedCurrentDensity( m_JaJ0_cc, contravar_B_withBCs );
      
      //   B = nabla x (Apar z_hat) ==> nabla^2 Apar = -Jz
      if(SpaceDim==3) { // Solve for Apar

         // set mapped coefficients and BCs for Poisson solver
         m_diffusionOp_Apar->setOperatorCoefficients(*m_bc_Apar);
     
         // fill dummy with Jz at cell center
         LevelData<FArrayBox>& Jz_cc = m_dummyFArray_oneComp;
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            Jz_cc[dit].copy(m_JaJ0_cc[dit],2,0,1); 
            m_Apar_cc[dit].setVal(0.0); // saw this in TEST_DIVERGENCE_CLEA MagGeom.cpp
         }
         m_geometry.divideJonValid(Jz_cc); // physical Jz at cell center

         // solve nabla^2 Az = -Jz
         m_diffusionOp_Apar->solvePreconditioner(Jz_cc, m_Apar_cc); // physical Apar ?

         // interpolate Az to edges and set BCs there
         LevelData<FArrayBox> A_cc;
         A_cc.define(grids,SpaceDim,2*IntVect::Unit);
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            A_cc[dit].setVal(0.0);
            A_cc[dit].copy(m_Apar_cc[dit],0,2,1); 
         }
         m_geometry.convertPhysToCovar(A_cc,0); // covariant Apar at cc
         LevelData<EdgeDataBox> Az_ce;
         Az_ce.define(grids,1,m_dummyFArray_oneComp.ghostVect());
         SpaceUtils::interpCellToEdges(Az_ce,A_cc,Az_ce,"c2"); // covariant Apar at ce
         m_fluid_face_var_bc.at(this_fbc)->setEdgeBC( a_species_phys, Az_ce, m_EdgeBC_zeros, a_time );
         SpaceUtils::exchangeEdgeDataBox(Az_ce);
         //
         //

         // compute B = curlA (divergence free and satisfies BCs)
         m_geometry.mappedGridCurl3D(contravar_B_withBCs,Az_ce);
         m_geometry.divideJonFaces(contravar_B_withBCs);
         SpaceUtils::exchangeFluxBox(contravar_B_withBCs); // must exchange before applyBC !!!!
         m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( a_species_phys, contravar_B_withBCs, a_time );
         LevelData<FluxBox>& mapped_contravar_B_noBCs( a_species_comp.face_var("magneticField") );
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            mapped_contravar_B_noBCs[dit].copy(contravar_B_withBCs[dit]);
         }
         m_geometry.multJonFaces(mapped_contravar_B_noBCs);         
   
      }
      //
      //
      ////////////////////////////////////////////////

      // compute Ja time contravariant curlE0 on faces
      setCellCenterValues( a_species_phys );
      const double gamma = a_species_phys.m_gamma;
      setMappedCspeed( gamma, 1, 1 );
      
      computeIdealEatEdges(m_E0_ce,m_velocity,contravar_B_withBCs,m_Cspeed_cc); // Ja*E0_covar
      m_fluid_face_var_bc.at(this_fbc)->setEdgeBC( a_species_phys, m_E0_ce, m_EdgeBC_zeros, a_time );
      SpaceUtils::exchangeEdgeDataBox(m_E0_ce); // call exchange after setBC !!!!!
      m_geometry.mappedGridCurl3D(m_curlE0_cf,m_E0_ce);

      // compute physical E0 at cell center
      SpaceUtils::interpEdgesToCell(m_E0_cc,m_E0_ce,"c2"); // covar at cell center
      m_geometry.convertPhysToCovar(m_E0_cc,1);            // phys at cell center
      
      //
      //
      ////////////////////////////////////////////////
      
   }
   
}     

void IdealMhdOp::computePhysicalDivergence( LevelData<FArrayBox>&  a_divF_phys,
                                      const LevelData<FluxBox>&    a_F_mapped,
                                      const int                    a_mapped ) const
{
   CH_TIME("IdealMhdOp::computePhysicalDivergence()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // compute physical divergence of F at cell center
   // 
   // if a_mapped=1, then a_F_mapped is Ja*F_contravariant on faces
   // else a_F_mapped is F_contravariant on faces
   //
   LevelData<FluxBox>& F_norm = m_dummyFlux_oneComp;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      F_norm[dit].copy(a_F_mapped[dit]);
   }
   if(a_mapped==0) m_geometry.multJonFaces(F_norm);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const RealVect& faceArea = coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& this_face = F_norm[dit][dir];
         this_face.mult(faceArea[dir]);
      }
   }
   m_geometry.mappedGridDivergenceFromFluxNorms( F_norm, a_divF_phys );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_divF_phys[dit].divide(m_cellVol[dit]);
         F_norm[dit].setVal(1./0.);
   }

}     

void IdealMhdOp::computeMappedCurrentDensity( LevelData<FArrayBox>&  a_JaJ0_phys,
                                        const LevelData<FArrayBox>&  a_By_contra ) const
{
   CH_TIME("IdealMhdOp::computeMappedCurrentDensity() (virtual B)");
   const DisjointBoxLayout& grids( m_geometry.grids() );
 
   // This function takes contravar By at cell center and returns 
   // Ja*curl(B)_phys at cell center (aka mapped)
   //
   LevelData<FArrayBox>& By_covar = m_dummyFArray_oneComp;
   LevelData<EdgeDataBox>& JaJ0_ce = m_dummyEdge_oneComp;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      By_covar[dit].copy(a_By_contra[dit]);
      if(m_twoDaxisymm) {
         By_covar[dit].mult(m_g_y[dit]);
         By_covar[dit].mult(m_g_y[dit]);
      }
   }
 
   m_geometry.mappedGridCurlofVirtComp(JaJ0_ce,By_covar); // Jacobian * curl(By)\cdot g^perp
   SpaceUtils::interpEdgesToCell(a_JaJ0_phys,JaJ0_ce,"c2");
   m_geometry.convertPhysToContravar(a_JaJ0_phys,1); // Ja*curlB_phys at cell center
   
   // reset dummys to infinity
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      By_covar[dit].setVal(1./0.);
      JaJ0_ce[dit].setVal(1./0.);
   }
   
}     

void IdealMhdOp::computeMappedCurrentDensity( LevelData<FArrayBox>&  a_JaJ0_phys,
                                        const LevelData<FluxBox>&    a_B_contra ) const
{
   CH_TIME("IdealMhdOp::computeMappedCurrentDensity()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
 
   // This function takes contravar B on faces and returns 
   // Ja*curl(B)_phys at cell center (aka mapped)
   //
   LevelData<FluxBox>& B_covar = m_dummyFlux_oneComp;
   LevelData<EdgeDataBox>& JaJ0_ce = m_dummyEdge_oneComp;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      B_covar[dit].copy(a_B_contra[dit]);
   }
   m_geometry.convertContravarToCovar(B_covar,0);
   m_geometry.mappedGridCurl3D(JaJ0_ce,B_covar); // Ja*curl(B)\cdot g^perp (contravariant)
   SpaceUtils::exchangeEdgeDataBox(JaJ0_ce);
      
   // compute physical Ja*curlB at cell center (aka mapped)
   //
   SpaceUtils::interpEdgesToCell(a_JaJ0_phys,JaJ0_ce,"c2");
   m_geometry.convertPhysToContravar(a_JaJ0_phys,1);
   
   // reset dummys to infinity
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      B_covar[dit].setVal(1./0.);
      JaJ0_ce[dit].setVal(1./0.);
   }
   
}     

void IdealMhdOp::computePhysicalElectricField( LevelData<FArrayBox>&  a_E_phys,
                                         const LevelData<FArrayBox>&  a_V_phys, 
                                         const LevelData<FArrayBox>&  a_B_phys,
                                         const int                    a_virtual ) const
{
   CH_TIME("IdealMhdOp::computePhysicalElectricField()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
 
   // compute physical E=BxV at cell center
   //

   if(a_virtual==1) {   
      CH_assert(a_B_phys.nComp()==1);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& thisbox = a_E_phys[dit].box();
         for (int n=0; n<SpaceDim; ++n) {
            a_E_phys[dit].copy( a_B_phys[dit],0,n,1 ); 
         }
         a_E_phys[dit].mult( a_V_phys[dit],thisbox,1,0,1 ); // Ex =  Vz*By
         a_E_phys[dit].mult( a_V_phys[dit],thisbox,0,1,1 ); // Ez = -Vx*By
         a_E_phys[dit].mult( -1.0,thisbox,1,1 );
      }
   }
   else {
      CH_assert(a_B_phys.nComp()==SpaceDim);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         
         const FArrayBox& B_on_patch = a_B_phys[dit];
         const FArrayBox& V_on_patch = a_V_phys[dit];
               FArrayBox& E_on_patch = a_E_phys[dit];

         const Box& cellbox = E_on_patch.box();
         FORT_EVAL_CROSS_PRODUCT( CHF_BOX(cellbox),
                                  CHF_CONST_FRA(B_on_patch),
                                  CHF_CONST_FRA(V_on_patch),
                                  CHF_FRA(E_on_patch) );
      }
   }
   
}     

#include "NamespaceFooter.H"
