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
     m_energy_is_thermal(false),
     m_energy_is_thermal_2(false),
     m_isothermal(false),
     m_enforceFluxBCs(false),
     m_use_perturbative_form(false),
     m_free_stream_preserving(false),
     m_advScheme("c2"),
     m_species_name(a_species_name),
     m_CGL(0),
     m_magP(1),
     m_blankN(false),
     m_blankB(false),
     m_magPforE0(1),
     m_thermPforE0(1),
     m_initializeBfromVectorPotential(0),
     m_initializeProfilesFromPSI(0),
     m_initializeBackgroundPlusPert(0), 
     m_addEnergyPerturbation(0),
     m_perturbDensity(0),
     m_perturbBtor(0),
     m_rbtor_fact(1.0),
     m_bpol_fact(1.0),
     m_opt_string(a_pp_str),
     m_my_pc_idx(-1),
     m_relative_density_floor(0.0),
     m_tau(0.0),
     m_etaMin(0.0)
{
   ParmParse pp(a_pp_str.c_str());
   parseParameters( pp );
   if (m_verbosity>0) {
      printParameters();
   }

   //  get initial conditions for cell variables and define boundary conditions instance
   const std::string icsn( "IC." + a_species_name);
   std:: string variable_name;
   ParmParse ppic( icsn.c_str() );
   FluidVarBCFactory fluid_var_bc_factory;
  
   // set IC and create BC for mass density
   variable_name = "density";
   CH_assert( ppic.contains("density.function") );
   m_density_cbc = m_fluid_bc.size();
   parseInitialConditions(a_species_name, variable_name);
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );
 
   // set IC and create BC for momentum density
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
   std::string name, prefix, prefix_BC;
   EllipticOpBCFactory elliptic_op_bc_factory;
  
   m_diffusionOp_mom0 = NULL;
   m_diffusionOp_mom1 = NULL;
   m_diffusionOp_mom2 = NULL;
   m_bc_mom0 = NULL;
   m_bc_mom1 = NULL;
   m_bc_mom2 = NULL;
   if(m_implicitViscosity) {
 
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
      }

   }
   
   m_diffusionOp_Apar = NULL;
   m_bc_Apar = NULL;
   if(SpaceDim==3 && m_initializeBfromVectorPotential) {
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

   // set IC and create BC for energy density
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
   
   if (m_initializeProfilesFromPSI) {

      variable_name = "RBTOR";
      parseInitialConditions(a_species_name, variable_name);

      variable_name = "PSI";
      parseInitialConditions(a_species_name, variable_name);
     
      if(m_addEnergyPerturbation) { 
         variable_name = "deltaE";
         parseInitialConditions(a_species_name, variable_name);
      }
      
      if(m_perturbBtor) { 
         variable_name = "deltaB";
         parseInitialConditions(a_species_name, variable_name);
      }

   }
   
   if (m_initializeBackgroundPlusPert) {

      variable_name = "RHO0";
      parseInitialConditions(a_species_name, variable_name);

      variable_name = "T0";
      parseInitialConditions(a_species_name, variable_name);
       
      variable_name = "BY0";
      parseInitialConditions(a_species_name, variable_name);
      
      variable_name = "deltaRHO";
      parseInitialConditions(a_species_name, variable_name);
     
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
   const DisjointBoxLayout& grids( m_geometry.grids() );
      
   setCellCenterValues( a_soln_fluid );
   
   if(m_advScheme=="ZIP") {
 
      setFaceCenteredFluxes_ZIP( a_soln_fluid );
      setMagneticFieldTerms( a_soln_fluid );
      setCurlE( a_soln_fluid, a_time );
  
   }
   else {
   
      if(m_CGL) setCellCenterFluxesCGL( a_soln_fluid );
      else setCellCenterFluxes( a_soln_fluid );

      // compute flux-freezing speed at cell-center for each direction
      // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/N); i = q0, q1, q2
      const double gamma = a_soln_fluid.m_gamma;
      setMappedCspeed( gamma, 1, m_magP );
   
      SpaceUtils::upWindToFaces(m_CspeedR_norm, m_Cspeed_cc, m_CspeedR_norm, "c2");
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_CspeedL_norm[dit].copy(m_CspeedR_norm[dit],m_CspeedL_norm[dit].box());
         m_CspeedL_norm[dit].negate();
      }

      setMagneticFieldTerms( a_soln_fluid );
   
      if(m_free_stream_preserving) setFaceCenteredFluxes_freestreamPreserving( a_soln_fluid );
      else setFaceCenteredFluxes( a_soln_fluid );

      enforceFluxBCs( a_soln_fluid, a_time );
   
      setMappedCspeed( gamma, m_thermPforE0, m_magPforE0 );
      setCurlE( a_soln_fluid, a_time );

   }

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
                                  const EMFields&                   a_EM_fields,
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
                                  const EMFields&                   a_EM_field,
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
                               const CFGVars&         a_fluid_vars,
                               const string&          a_name ) const
{
   const FluidSpecies& fluid_species = static_cast<const FluidSpecies&>(a_fluid_vars);

   const DisjointBoxLayout& grids( m_geometry.grids() );
   const IntVect ghostVect = 2*IntVect::Unit;
   
   if(a_name=="electricField") {
      a_Var.define(grids, SpaceDim, 0*IntVect::Unit);
   
      switch(fluid_species.m_evolveB_type) {
         case FluidSpecies::EVOLVE_B_TYPE::NONE : 
            break;
         case FluidSpecies::EVOLVE_B_TYPE::TWOD_VIRTUAL_B : 
            getE0_2Dvirtual( a_Var, fluid_species );
            break;
         case FluidSpecies::EVOLVE_B_TYPE::TWOD_INPLANE_B : 
            getE0_2DinPlane( a_Var, fluid_species );
            break;
         case FluidSpecies::EVOLVE_B_TYPE::TWOD_FULL_B : 
            getE0_2Dfull( a_Var, fluid_species );
            break;
         case FluidSpecies::EVOLVE_B_TYPE::THREED : 
            getE0_3D( a_Var, fluid_species );
            break;
         default :
            MayDay::Error("Invalid evolve B type in getMemberVar");
      }
   }
   if(a_name=="electricField0") { // off by time stage from "electricField" above
      a_Var.define(grids, SpaceDim, 2*IntVect::Unit);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_E0_cc[dit]);
      }
   }
   if(a_name=="electricField_upwinded") { // need ghost cells to calc here, so one stage out of sync 
      a_Var.define(grids, SpaceDim, 2*IntVect::Unit);
      SpaceUtils::interpEdgesToCell(a_Var,m_E0_ce,"c2"); // covar at cell center
      m_geometry.convertPhysToCovar(a_Var,1);            // phys at cell center
   }
   if(a_name=="electricField0_virtual") { // off by time stage from "electricField" above
      a_Var.define(grids, 1, ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_E0_virtual_cc[dit]);
      }
   }
   if(a_name=="electricField0_virtual_nc") { // off by time stage from "electricField" above
      a_Var.define(grids, 1, ghostVect);
      SpaceUtils::interpNodesToCells( a_Var, m_E0_virtual_nc, "c2" );
      for (auto dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         if(m_twoDaxisymm) {
            a_Var[dit].divide(m_g_y[dit],0,0,1); // convert to physical
         }
      }
      a_Var.exchange();
   }
   if(a_name=="curlE_virtual") { // off by time stage from "electricField" above
      a_Var.define(grids, 1, ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_curlE0_virtual[dit]);
      }
   }
   if(a_name=="curlE") { // off by time stage
      a_Var.define(grids, SpaceDim, IntVect::Zero);
      SpaceUtils::interpFaceVectorToCell( a_Var, m_curlE0_cf, "c2" );
      m_geometry.convertPhysToContravar(a_Var,1);
      m_geometry.divideJonValid(a_Var); // physical curlE0
   }
   if(a_name=="currentDensity") {
      //a_Var.define(grids, SpaceDim, ghostVect);
      a_Var.define(grids, SpaceDim, 0*IntVect::Unit);
      
      switch(fluid_species.m_evolveB_type) {
         case FluidSpecies::EVOLVE_B_TYPE::NONE : 
            break;
         case FluidSpecies::EVOLVE_B_TYPE::TWOD_VIRTUAL_B : 
            getJ0_2Dvirtual( a_Var, fluid_species );
            break;
         case FluidSpecies::EVOLVE_B_TYPE::TWOD_INPLANE_B : 
            getJ0_2DinPlane( a_Var, fluid_species );
            break;
         case FluidSpecies::EVOLVE_B_TYPE::TWOD_FULL_B : 
            getJ0_2Dfull( a_Var, fluid_species );
            break;
         case FluidSpecies::EVOLVE_B_TYPE::THREED : 
            getJ0_3D( a_Var, fluid_species );
            break;
         default :
            MayDay::Error("Invalid evolve B type in getMemberVar");
      }
      //for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
      //   a_Var[dit].copy(m_JaJ0_cc[dit]); // off by time stage
      //}
      //m_geometry.divideJonValid(a_Var);
   }
   if(a_name=="currentDensity2") { // off by time stage from "electricField" above
      a_Var.define(grids, SpaceDim, IntVect::Zero);
      getJ0_3D2( a_Var, fluid_species );
   }
   if(a_name=="currentDensity_virtual") {
      a_Var.define(grids, 1, ghostVect);
   
      const LevelData<FluxBox>& B_contra_cf( fluid_species.face_var("magneticField") );
      LevelData<FluxBox>& B_contra_withBCs = m_dummyFlux_oneComp;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         B_contra_withBCs[dit].copy(B_contra_cf[dit]);
      }
      SpaceUtils::exchangeFluxBox(B_contra_withBCs);
   
      const double dummy_time = 0.0;
      int this_fbc = m_magneticField_fbc;
      m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( fluid_species, B_contra_withBCs, dummy_time );
   
      LevelData<FArrayBox>& B_covar_cc = m_dummyFArray_spaceDim;
      SpaceUtils::interpFaceVectorToCell(B_covar_cc, B_contra_withBCs,"c2");
      m_geometry.convertPhysToContravar(B_covar_cc,1);
      m_geometry.convertPhysToCovar(B_covar_cc,0);
      
      // JRA, fix this. Need correct BCs for B_covar_cc. What are they?
      //int this_cbc = m_energyDensity_cbc;
      //m_fluid_bc.at(this_cbc)->applyCellBC( fluid_species, B_covar_cc, dummy_time );

      m_geometry.mappedGridCurl( a_Var, B_covar_cc ); // Ja*curl(B)_phys
      m_geometry.divideJonValid(a_Var);
            
      // rebase dummy container
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         B_contra_withBCs[dit].setVal(1./0.);
         B_covar_cc[dit].setVal(1./0.);
      }
  
   }
   if(a_name=="density_withBCs") {
      a_Var.define(grids, 1, ghostVect);
      const LevelData<FArrayBox>& N( fluid_species.cell_var(0) );
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(N[dit]);
      }
      a_Var.exchange();
      int this_cbc = m_density_cbc;
      const double dummy_time = 0.0;
      m_fluid_bc.at(this_cbc)->applyCellBC( fluid_species, a_Var, dummy_time );
   }
   if(a_name=="momentumDensity_withBCs") {
      a_Var.define(grids, SpaceDim, ghostVect);
      const LevelData<FArrayBox>& rhoV( fluid_species.cell_var("momentumDensity") );
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(rhoV[dit]);
      }
      m_geometry.convertPhysToContravar(a_Var,0);
      a_Var.exchange();
      int this_cbc = m_momentumDensity_cbc;
      const double dummy_time = 0.0;
      m_fluid_bc.at(this_cbc)->applyCellBC( fluid_species, a_Var, dummy_time );
      //m_geometry.convertContravarToCovar(a_Var,0); // output covar
      m_geometry.convertPhysToContravar(a_Var,1); // output phys
   }
   if(a_name=="divB") {
      a_Var.define(grids, 1, IntVect::Zero);
      CH_assert(fluid_species.m_evolve_magneticField==1);
      const LevelData<FluxBox>& magField_contra_cf( fluid_species.face_var("magneticField") );
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
   if(a_name=="pressure") {
      a_Var.define(grids, 1+m_CGL, IntVect::Zero);
      const double gamma = fluid_species.m_gamma;
      if(m_energy_is_thermal || m_energy_is_thermal_2) {
         const LevelData<FArrayBox>& eneDen( fluid_species.cell_var("energyDensity") );
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            a_Var[dit].copy(eneDen[dit],0,0,1+m_CGL);
            a_Var[dit].mult(gamma-1.0);
         }         
      }
      else if(m_isothermal) {
         const LevelData<FArrayBox>& rhoDen( fluid_species.cell_var(0) );
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            a_Var[dit].copy(rhoDen[dit],0,0,1);
            a_Var[dit].mult(m_T0[dit],0,0,1);
            a_Var[dit].mult(2.0);
         }
      }
      else {
         fluid_species.pressure(a_Var);
      }
   }
   if(a_name=="P0") {
      a_Var.define(grids, 1, 2*IntVect::Unit);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_P0[dit]);
      }
   }
   if(a_name=="N0") {
      a_Var.define(grids, 1, 2*IntVect::Unit);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_N0[dit]);
      }
   }
   if(a_name=="B0") {
      a_Var.define(grids, SpaceDim, IntVect::Zero);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_B0_phys[dit]);
      }
   }
   if(a_name=="J00") {
      a_Var.define(grids, SpaceDim, IntVect::Zero);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_JaJ00[dit]);
      }
      m_geometry.divideJonValid(a_Var);
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
   if(a_name=="PSI") {
      a_Var.define(grids, m_PSI.nComp(), IntVect::Zero);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_PSI[dit]);
      }
   }
   if(a_name=="gradPressure") {

      a_Var.define(grids, SpaceDim, IntVect::Unit);

      LevelData<FArrayBox>& pressure = m_dummyFArray_oneComp;
      fluid_species.pressure(pressure);
      pressure.exchange();

      const double dummy_time = 0.0;
      int this_cbc = m_energyDensity_cbc;
      m_fluid_bc.at(this_cbc)->applyCellBC( fluid_species, pressure, dummy_time );

      //m_geometry.computeMappedPoloidalGradientWithGhosts( pressure, a_Var, 2);
      m_geometry.computeMappedGradient( pressure, a_Var, 2);
      m_geometry.convertPhysToCovar(a_Var,1);

   }
   if(a_name=="forceDensity") {
      a_Var.define(grids, SpaceDim, IntVect::Zero);
      for (auto dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_forceDensity[dit]);
      }
      m_geometry.divideJonValid(a_Var);
   }
   if(a_name=="rhs_rho") {
      a_Var.define(grids, 1, IntVect::Zero);
      for (auto dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_rhs_rho[dit]);
      }
      m_geometry.divideJonValid(a_Var);
   }
   if(a_name=="mxFlux") {
      a_Var.define(grids, SpaceDim, 2*IntVect::Unit);
      m_geometry.multiplyNTranspose(a_Var, m_mxFlux_cc);
   }
   if(a_name=="myFlux") {
      a_Var.define(grids, SpaceDim, 2*IntVect::Unit);
      m_geometry.multiplyNTranspose(a_Var, m_myFlux_cc);
   }
   if(a_name=="mzFlux") {
      a_Var.define(grids, SpaceDim, 2*IntVect::Unit);
      m_geometry.multiplyNTranspose(a_Var, m_mzFlux_cc);
   }
   if(a_name=="contravarB") {
      a_Var.define(grids, SpaceDim, 2*IntVect::Unit);
      const LevelData<FluxBox>& B_contra_cf( fluid_species.face_var("magneticField") );
      LevelData<FluxBox>& B_contra_withBCs = m_dummyFlux_oneComp;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         B_contra_withBCs[dit].copy(B_contra_cf[dit]);
      }
      SpaceUtils::exchangeFluxBox(B_contra_withBCs);
      const double dummy_time = 0.0;
      int this_fbc = m_magneticField_fbc;
      m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( fluid_species, B_contra_withBCs, dummy_time );
      SpaceUtils::interpFaceVectorToCell(a_Var, B_contra_withBCs,"c2");
   }
   if(a_name=="covarB") {
      a_Var.define(grids, SpaceDim, 0*IntVect::Unit);
      const LevelData<FluxBox>& B_contra_cf( fluid_species.face_var("magneticField") );
      LevelData<FluxBox>& B_contra_withBCs = m_dummyFlux_oneComp;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         B_contra_withBCs[dit].copy(B_contra_cf[dit]);
      }
      SpaceUtils::exchangeFluxBox(B_contra_withBCs);
      const double dummy_time = 0.0;
      int this_fbc = m_magneticField_fbc;
      m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( fluid_species, B_contra_withBCs, dummy_time );
      
      SpaceUtils::interpFaceVectorToCell(a_Var, B_contra_withBCs, "c2");
      m_geometry.convertContravarToCovar(a_Var,0);
   }
   if(a_name=="covarB2_0") {
      a_Var.define(grids, SpaceDim, 0*IntVect::Unit);
      const LevelData<FluxBox>& B_contra_cf( fluid_species.face_var("magneticField") );
      LevelData<FluxBox>& B_contra_withBCs = m_dummyFlux_oneComp;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         B_contra_withBCs[dit].copy(B_contra_cf[dit]);
      }
      SpaceUtils::exchangeFluxBox(B_contra_withBCs);
      const double dummy_time = 0.0;
      int this_fbc = m_magneticField_fbc;
      m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( fluid_species, B_contra_withBCs, dummy_time );
      
      LevelData<FluxBox>& Bvec_covar = m_dummyFlux_spaceDim;
      SpaceUtils::interpFacesToFaces(Bvec_covar, B_contra_withBCs);
      m_geometry.convertContravarToCovar(Bvec_covar,0);
      SpaceUtils::interpFaceScalarToCell(a_Var, Bvec_covar, 0);
   }
   if(a_name=="covarB2_1") {
      a_Var.define(grids, SpaceDim, 0*IntVect::Unit);
      const LevelData<FluxBox>& B_contra_cf( fluid_species.face_var("magneticField") );
      LevelData<FluxBox>& B_contra_withBCs = m_dummyFlux_oneComp;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         B_contra_withBCs[dit].copy(B_contra_cf[dit]);
      }
      SpaceUtils::exchangeFluxBox(B_contra_withBCs);
      const double dummy_time = 0.0;
      int this_fbc = m_magneticField_fbc;
      m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( fluid_species, B_contra_withBCs, dummy_time );
      
      LevelData<FluxBox>& Bvec_covar = m_dummyFlux_spaceDim;
      SpaceUtils::interpFacesToFaces(Bvec_covar, B_contra_withBCs);
      m_geometry.convertContravarToCovar(Bvec_covar,0);
      SpaceUtils::interpFaceScalarToCell(a_Var, Bvec_covar, 1);
   }
   if(a_name=="covarB2_2") {
      a_Var.define(grids, SpaceDim, 0*IntVect::Unit);
      const LevelData<FluxBox>& B_contra_cf( fluid_species.face_var("magneticField") );
      LevelData<FluxBox>& B_contra_withBCs = m_dummyFlux_oneComp;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         B_contra_withBCs[dit].copy(B_contra_cf[dit]);
      }
      SpaceUtils::exchangeFluxBox(B_contra_withBCs);
      const double dummy_time = 0.0;
      int this_fbc = m_magneticField_fbc;
      m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( fluid_species, B_contra_withBCs, dummy_time );
      
      LevelData<FluxBox>& Bvec_covar = m_dummyFlux_spaceDim;
      SpaceUtils::interpFacesToFaces(Bvec_covar, B_contra_withBCs);
      m_geometry.convertContravarToCovar(Bvec_covar,0);
      SpaceUtils::interpFaceScalarToCell(a_Var, Bvec_covar, 2);
   }
   if(a_name=="curlConst") {
      a_Var.define(grids, SpaceDim, 0*IntVect::Unit);

      // create constant vector on cells
      LevelData<FArrayBox>& const_vect = m_dummyFArray_spaceDim;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         const_vect[dit].setVal(1.0);
      }
      
      // convert to covariant
      m_geometry.convertPhysToCovar(const_vect,0);

      // compute physical curl
      m_geometry.mappedGridCurl(a_Var,const_vect);
      m_geometry.convertPhysToContravar(a_Var,1);
      m_geometry.divideJonValid(a_Var); // physical curl
   }
   if(a_name=="curlConst2") {
      a_Var.define(grids, SpaceDim, 0*IntVect::Unit);

      // create constant vector on faces
      LevelData<FluxBox>& const_vect = m_dummyFlux_spaceDim;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         const_vect[dit].setVal(1.0);
      }
      
      // convert to covariant
      m_geometry.convertPhysToCovar(const_vect,0);

      // compute physical curl
      m_geometry.mappedGridCurl(a_Var,const_vect);
      m_geometry.convertPhysToContravar(a_Var,1);
      m_geometry.divideJonValid(a_Var); // physical curl
   }
   if(a_name=="contravarV") {
      a_Var.define(grids, SpaceDim, 2*IntVect::Unit);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
        a_Var[dit].copy(m_velocity[dit]);
      }
      m_geometry.convertPhysToContravar(a_Var,0);
   } 


}

void IdealMhdOp::getE0_2Dvirtual( LevelData<FArrayBox>&  a_E0,
                            const FluidSpecies&          a_fluid_species ) const
{

   const LevelData<FArrayBox>& B_virtual( a_fluid_species.cell_var("magneticField_virtual") );
   LevelData<FArrayBox>& V_inPlane = m_dummyFArray_spaceDim;
   a_fluid_species.velocity(V_inPlane);
   computeCrossProduct( a_E0, B_virtual, V_inPlane, 2 );
      
   // rebase dummy container
   for (DataIterator dit(a_E0.dataIterator()); dit.ok(); ++dit) {
      V_inPlane[dit].setVal(1./0.);
   }

}

void IdealMhdOp::getE0_2DinPlane( LevelData<FArrayBox>&  a_E0,
                            const FluidSpecies&          a_fluid_species ) const
{
            
   const LevelData<FluxBox>& B_contra_cf( a_fluid_species.face_var("magneticField") );
   LevelData<FArrayBox>& V_virtual = m_dummyFArray_oneComp;
   a_fluid_species.velocity_virtual(V_virtual);
   LevelData<FArrayBox>& B_inPlane = m_dummyFArray_spaceDim;
   SpaceUtils::interpFaceVectorToCell(B_inPlane,B_contra_cf,"c2");
   m_geometry.convertPhysToContravar(B_inPlane,1);
   computeCrossProduct( a_E0, B_inPlane, V_virtual, 3 );

   // rebase dummy container
   for (DataIterator dit(a_E0.dataIterator()); dit.ok(); ++dit) {
      B_inPlane[dit].setVal(1./0.);
      V_virtual[dit].setVal(1./0.);
   }

}

void IdealMhdOp::getE0_2Dfull( LevelData<FArrayBox>&  a_E0,
                         const FluidSpecies&          a_fluid_species ) const
{
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   const IntVect ghostVect = a_E0.ghostVect();
            
   const LevelData<FArrayBox>& B_virtual( a_fluid_species.cell_var("magneticField_virtual") );
   const LevelData<FluxBox>& B_contra_cf( a_fluid_species.face_var("magneticField") );
   LevelData<FArrayBox>& V_inPlane = m_dummyFArray_spaceDim;
   LevelData<FArrayBox>& V_virtual = m_dummyFArray_oneComp;
   a_fluid_species.velocity(V_inPlane);
   a_fluid_species.velocity_virtual(V_virtual);

   LevelData<FArrayBox> B_inPlane;
   B_inPlane.define(grids, SpaceDim, ghostVect);
   SpaceUtils::interpFaceVectorToCell(B_inPlane,B_contra_cf,"c2");
   m_geometry.convertPhysToContravar(B_inPlane,1);

   computeCrossProduct( a_E0, B_inPlane, V_inPlane,  
                              B_virtual, V_virtual );

   // rebase dummy container
   for (DataIterator dit(a_E0.dataIterator()); dit.ok(); ++dit) {
      B_inPlane[dit].setVal(1./0.);
      V_virtual[dit].setVal(1./0.);
   }

}

void IdealMhdOp::getE0_3D( LevelData<FArrayBox>&  a_E0,
                     const FluidSpecies&          a_fluid_species ) const
{
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   const IntVect ghostVect = a_E0.ghostVect();
   
   const LevelData<FluxBox>& B_contra_cf( a_fluid_species.face_var("magneticField") );

   LevelData<FArrayBox>& V_inPlane = m_dummyFArray_spaceDim;
   a_fluid_species.velocity(V_inPlane);

   LevelData<FArrayBox> B_inPlane;
   B_inPlane.define(grids, SpaceDim, ghostVect);
   SpaceUtils::interpFaceVectorToCell(B_inPlane,B_contra_cf,"c2");
   m_geometry.convertPhysToContravar(B_inPlane,1);

   computeCrossProduct( a_E0, B_inPlane, V_inPlane, 0 ); 

   // rebase dummy container
   for (DataIterator dit(a_E0.dataIterator()); dit.ok(); ++dit) {
      //a_E0[dit].copy(B_inPlane[dit]); // testing....
      //a_E0[dit].copy(V_inPlane[dit]); // testing....
      V_inPlane[dit].setVal(1./0.);
   }

}

void IdealMhdOp::getJ0_2Dvirtual( LevelData<FArrayBox>&  a_J0,
                            const FluidSpecies&          a_fluid_species ) const
{

   const LevelData<FArrayBox>& B_virtual( a_fluid_species.cell_var("magneticField_virtual") );
   LevelData<FArrayBox>& By_covar = m_dummyFArray_oneComp;
   for (DataIterator dit(a_J0.dataIterator()); dit.ok(); ++dit) {
      By_covar[dit].copy(B_virtual[dit]);
   }
   By_covar.exchange();
   
   const double dummy_time = 0.0;
   int this_cbc = m_magneticField_virtual_cbc;
   m_fluid_bc.at(this_cbc)->applyCellBC( a_fluid_species, By_covar, dummy_time );
   for (DataIterator dit(a_J0.dataIterator()); dit.ok(); ++dit) {
      if(m_twoDaxisymm) {
         By_covar[dit].mult(m_g_y[dit],0,0,1);
      }
   }
   m_geometry.mappedGridCurlofVirtComp(a_J0,By_covar); // Ja*curl(By)\cdot g^perp
   m_geometry.convertPhysToContravar(a_J0,1);          // Ja*curlB_phys at cell center
   m_geometry.divideJonValid(a_J0);
   
   // rebase dummy container
   for (DataIterator dit(a_J0.dataIterator()); dit.ok(); ++dit) {
      By_covar[dit].setVal(1./0.);
   }

}

void IdealMhdOp::getJ0_2DinPlane( LevelData<FArrayBox>&  a_J0,
                            const FluidSpecies&          a_fluid_species ) const
{

   for (DataIterator dit(a_J0.dataIterator()); dit.ok(); ++dit) {
      a_J0[dit].setVal(0.0);
   }
           
}

void IdealMhdOp::getJ0_2Dfull( LevelData<FArrayBox>&  a_J0,
                         const FluidSpecies&          a_fluid_species ) const
{
   getJ0_2Dvirtual (a_J0, a_fluid_species );   
}

void IdealMhdOp::getJ0_3D( LevelData<FArrayBox>&  a_J0,
                     const FluidSpecies&          a_fluid_species ) const
{
   
   const LevelData<FluxBox>& B_contra_cf( a_fluid_species.face_var("magneticField") );
   LevelData<FluxBox>& B_contra_withBCs = m_dummyFlux_oneComp;
   for (DataIterator dit(a_J0.dataIterator()); dit.ok(); ++dit) {
      B_contra_withBCs[dit].copy(B_contra_cf[dit]);
   }
   SpaceUtils::exchangeFluxBox(B_contra_withBCs);
   const double dummy_time = 0.0;
   int this_fbc = m_magneticField_fbc;
   m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( a_fluid_species, B_contra_withBCs, dummy_time );
  
   LevelData<FArrayBox>& B_covar_cc = m_dummyFArray_spaceDim;
   SpaceUtils::interpFaceVectorToCell(B_covar_cc, B_contra_withBCs,"c2");
   m_geometry.convertPhysToContravar(B_covar_cc,1);
   m_geometry.convertPhysToCovar(B_covar_cc,0);
   m_geometry.mappedGridCurl( a_J0, B_covar_cc); // Ja*curl(B)_contra
   
   m_geometry.convertPhysToContravar(a_J0,1);
   m_geometry.divideJonValid(a_J0);

}  

void IdealMhdOp::getJ0_3D2( LevelData<FArrayBox>&  a_J0,
                     const FluidSpecies&          a_fluid_species ) const
{
   
   const LevelData<FluxBox>& B_contra_cf( a_fluid_species.face_var("magneticField") );
   LevelData<FluxBox>& B_contra_withBCs = m_dummyFlux_oneComp;
   for (DataIterator dit(a_J0.dataIterator()); dit.ok(); ++dit) {
      B_contra_withBCs[dit].copy(B_contra_cf[dit]);
   }
   SpaceUtils::exchangeFluxBox(B_contra_withBCs);
   const double dummy_time = 0.0;
   int this_fbc = m_magneticField_fbc;
   m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( a_fluid_species, B_contra_withBCs, dummy_time );
  
   LevelData<FluxBox>& Bvec_covar = m_dummyFlux_spaceDim;
   SpaceUtils::interpFacesToFaces(Bvec_covar, B_contra_withBCs);
   m_geometry.convertContravarToCovar(Bvec_covar,0);
   m_geometry.mappedGridCurl( a_J0, Bvec_covar ); // Ja*curl(B)_contra

   m_geometry.convertPhysToContravar(a_J0,1);
   m_geometry.divideJonValid(a_J0);

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
                                 const int                                                    a_species_idx,
                                 const int                                                    a_id )
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
                              const int                        a_component,
                              const std::string& )
{
   CH_TIME("IdealMhdOp::updatePCImEx()");
   
   if (!m_implicitViscosity) {
     MayDay::Error("Error in IdealMhdOp::updatePCImEx(): m_implicitViscosity is false." 
                   "This function shouldn't have been called at all!");
   }
   
   // Get computational fluid species data
   const FluidSpecies& fluid_species = static_cast<const FluidSpecies&>(*a_fluid_species[a_component]);
   const LevelData<FArrayBox>& Jrho( fluid_species.cell_var("density") );
   m_alpha_PC = a_shift;   

   for (DataIterator dit(m_shift_PC.dataIterator()); dit.ok(); ++dit) {
      m_Jrho_PC[dit].copy(Jrho[dit]);
      m_shift_PC[dit].setVal(m_alpha_PC);
      m_shift_PC[dit].mult(m_Jrho_PC[dit]);
   }
   m_geometry.divideJonValid(m_shift_PC);

   // Create elliptic coefficients
   // Need to update m_etaVisc_cf ? careful, m_pressure above is mapped
   const double a_dummy_time=0.0;
   computeDiffusionCoefficients(m_ellip_coeff, m_ellip_coeff_mapped, m_etaVisc_cf, a_dummy_time);
   m_diffusionOp_mom0->setOperatorCoefficients(m_ellip_coeff, m_ellip_coeff_mapped, *m_bc_mom0);
   m_diffusionOp_mom1->setOperatorCoefficients(m_ellip_coeff, m_ellip_coeff_mapped, *m_bc_mom1);
   if(SpaceDim==3) m_diffusionOp_mom2->setOperatorCoefficients(m_ellip_coeff, m_ellip_coeff_mapped, *m_bc_mom2);
   //SpaceUtils::inspectFluxBox(m_ellip_coeff,1);
   //m_geometry.plotFaceData( "m_ellip_coeff_mapped", m_ellip_coeff_mapped, a_dummy_time);
   
   m_diffusionOp_mom0->updateImExPreconditioner( m_shift_PC, *m_bc_mom0 );
   m_diffusionOp_mom1->updateImExPreconditioner( m_shift_PC, *m_bc_mom1 );
   if(SpaceDim==3) m_diffusionOp_mom2->updateImExPreconditioner( m_shift_PC, *m_bc_mom2 );
}


void IdealMhdOp::solvePCImEx( FluidSpeciesPtrVect&        a_fluid_species_solution,
                        const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                        const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                        const std::string&,
                        const int                         a_component )
{
   CH_TIME("IdealMhdOp::solvePCImEx()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // This function computes solution to the preconditioned problem.
   // treating Fi(y) as linear operator, the equation solved here is
   // alpha*y - Fi(y) = RHS, y is mapped state vector
   //
   // momentum density and energy density have implicit parts

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
 
}

void IdealMhdOp::upwindFluxVector( LevelData<FluxBox>&    a_Flux_phys_fc,
                             const LevelData<FArrayBox>&  a_Flux_phys_cc,
                             const LevelData<FArrayBox>&  a_f_cc,
                             const LevelData<FArrayBox>&  a_Cspeed_cc,
                             const LevelData<FluxBox>&    a_CspeedR_norm,
                             const LevelData<FluxBox>&    a_CspeedL_norm ) const
{
   CH_TIME("IdealMhdOp::upwindFluxVector()");
   CH_assert(a_Flux_phys_fc.nComp() == SpaceDim);
   CH_assert(a_Flux_phys_cc.nComp() == SpaceDim);
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      // compute Lax splitting on mapped flux and interpolate FluxL and FluxR to face norms
      //
      //   get left and right going flux at cell-center
      //   fluxR = 0.5*(flux + Cspeed*fun),
      //   fluxL = 0.5*(flux - Cspeed*fun),
      //   Cspeed = abs(max(eigenValue of Flux Jacobian))

      const Box& cell_box( a_f_cc[dit].box() );
      FArrayBox FluxL_cc(cell_box, SpaceDim);
      FArrayBox FluxR_cc(cell_box, SpaceDim);
      const FArrayBox& Flux_on_patch = a_Flux_phys_cc[dit];
      const FArrayBox& Cspeed_on_patch = a_Cspeed_cc[dit];
      const FArrayBox& f_on_patch = a_f_cc[dit];
      for (int n=0; n<SpaceDim; n++) {
         FORT_LAX_FLUX_SPLITTING_NEW( CHF_BOX(cell_box), 
                                  CHF_FRA(FluxL_cc),
                                  CHF_FRA(FluxR_cc),
                                  CHF_CONST_FRA1(Flux_on_patch,n),
                                  CHF_CONST_FRA(Cspeed_on_patch),
                                  CHF_CONST_FRA1(f_on_patch,0) );

         const Box& grid_box( grids[dit] );             // this box has no ghost
         const Box& flux_box( a_Flux_phys_fc[dit].box() ); // this box has ghost
         FluxBox thisFluxL( flux_box, 1 );  // need box with ghosts here
         FluxBox thisFluxR( flux_box, 1 );
         const FluxBox& thisCspeedL_norm = a_CspeedL_norm[dit];
         const FluxBox& thisCspeedR_norm = a_CspeedR_norm[dit];
         SpaceUtils::upWindToFaces(thisFluxL, FluxL_cc, thisCspeedL_norm, grid_box, m_advScheme);
         SpaceUtils::upWindToFaces(thisFluxR, FluxR_cc, thisCspeedR_norm, grid_box, m_advScheme);

         // compute total norm flux by adding left and right together 
         thisFluxL += thisFluxR;
         a_Flux_phys_fc[dit].copy(thisFluxL, 0, n, 1 );
      }
      
   }

}

void IdealMhdOp::computeIdealEatEdges( LevelData<EdgeDataBox>&  a_E_covar,
                                 const LevelData<FArrayBox>&    a_V_phys_cc,
                                 const LevelData<FluxBox>&      a_B_contra,
                                 const LevelData<FArrayBox>&    a_Cspeed_cc ) const
{
   CH_TIME("IdealMhdOp::computeIdealEatEdges()");
   CH_assert(a_E_covar.nComp() == 1);
   CH_assert(a_B_contra.nComp() == 1);
   CH_assert(a_V_phys_cc.nComp() == SpaceDim);
   CH_assert(a_Cspeed_cc.nComp() == SpaceDim);
   
   //
   // compute covar E0 = -VxB on edges:
   // E_i = Ja*(B^jV^k - B^kV^j)
   //
  
   const DisjointBoxLayout& grids( m_geometry.grids() );
  
   // multiply contravariant B on faces by Jacobian and convert
   // physical velocity to contravariant
   LevelData<FArrayBox>& V_contra_cc = m_dummyFArray_spaceDim;
   LevelData<FluxBox>& JaB_contra_cf(m_dummyFlux_oneComp);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      JaB_contra_cf[dit].copy(a_B_contra[dit]);
      const Box& cell_box(a_V_phys_cc[dit].box());
      SpaceUtils::copy(V_contra_cc[dit],a_V_phys_cc[dit],cell_box);
   } 
   m_geometry.multJonFaces(JaB_contra_cf);
   m_geometry.convertPhysToContravar(V_contra_cc,0);

   // compute constrained covariant E on cell edges
   SpaceUtils::constrainedCrossProductOnEdges( a_E_covar, JaB_contra_cf,
                                               V_contra_cc, a_Cspeed_cc, m_advScheme_B );
   
   // rebase dummies
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      V_contra_cc[dit].setVal(1.0/0.0);
      JaB_contra_cf[dit].setVal(1.0/0.0);
   } 

}

void IdealMhdOp::computeIdealEatEdges( LevelData<EdgeDataBox>&  a_E_covar,
                                 const LevelData<FArrayBox>&    a_V_phys_cc,
                                 const LevelData<FArrayBox>&    a_B_contr_cc,
                                 const LevelData<FArrayBox>&    a_Cspeed_cc ) const
{
   CH_TIME("IdealMhdOp::computeIdealEatEdges()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   CH_assert(a_E_covar.nComp() == 1);
   CH_assert(a_V_phys_cc.nComp() == SpaceDim);
   CH_assert(a_Cspeed_cc.nComp() == SpaceDim);
  
   // compute covar E0=-VxB at cell center
   // E_x =  Ja*V^z*B^y
   // E_z = -Ja*V^x*B^y
   LevelData<FArrayBox>& E0_covar_cc = m_dummyFArray_spaceDim;
   LevelData<FArrayBox>& JaB_contr_cc = m_dummyFArray_oneComp;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box& cell_box(a_V_phys_cc[dit].box());
      SpaceUtils::copy(E0_covar_cc[dit],a_V_phys_cc[dit],cell_box);
   }
   m_geometry.convertPhysToContravar(E0_covar_cc,0); // contains contravar V
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box& cell_box(a_V_phys_cc[dit].box());
      SpaceUtils::copy(JaB_contr_cc[dit],a_B_contr_cc[dit],cell_box);
      JaB_contr_cc[dit].mult(m_Jacobian[dit]);
      for (int n=0; n<SpaceDim; n++) {
         E0_covar_cc[dit].mult(JaB_contr_cc[dit],cell_box,cell_box,0,n,1); // [Ja*V^x*B^y, Ja*V^z*B^y]
      }
   }

   // upwind E0_covar_cc appropriately to cell edges
   for (DataIterator dit(grids); dit.ok(); ++dit) {
    
      // Compute Left/Right Cspeeds
      const Box& grid_box( grids[dit] ); // this box has no ghosts
      const Box& cell_box(a_B_contr_cc[dit].box());
      FluxBox CspeedL_norm(cell_box, 1);
      FluxBox CspeedR_norm(cell_box, 1);
      SpaceUtils::upWindToFaces(CspeedR_norm, a_Cspeed_cc[dit], CspeedR_norm, grid_box, "c2"); // 3rd arg not used for "c2"
      CspeedL_norm.copy(CspeedR_norm, cell_box );
      CspeedL_norm.negate();
      
      // Do Lax flux splitting
      FArrayBox FluxL_cc(cell_box, SpaceDim);
      FArrayBox FluxR_cc(cell_box, SpaceDim);
      const FArrayBox& Cspeed_on_patch = a_Cspeed_cc[dit];

      FORT_LAX_FLUX_SPLITTING( CHF_BOX(cell_box),
                               CHF_FRA(FluxL_cc),
                               CHF_FRA(FluxR_cc),
                               CHF_CONST_FRA(E0_covar_cc[dit]),
                               CHF_CONST_FRA(Cspeed_on_patch),
                               CHF_CONST_FRA1(JaB_contr_cc[dit],0) );
     
      // Upwind
      const Box& flux_box( a_E_covar[dit].box() ); // this box has ghosts
      FluxBox FluxL_norm( flux_box, 1 );
      FluxBox FluxR_norm( flux_box, 1 );
      FluxBox VB_norm( flux_box, 1 );
      if(m_advScheme_B=="weno5") {
         cout << "JRA: weno5 not working here !!!! " << endl;
         //const FArrayBox& thisSmooth = m_NTrowSq[dit];
         //SpaceUtils::interpToFacesWENO(FluxL_norm, FluxL_cc, CspeedL_norm, thisSmooth, grid_box, m_advScheme_B);
         //SpaceUtils::interpToFacesWENO(FluxR_norm, FluxR_cc, CspeedR_norm, thisSmooth, grid_box, m_advScheme_B);
      } 
      else {
         SpaceUtils::upWindToFaces(FluxL_norm, FluxL_cc, CspeedL_norm, grid_box, m_advScheme_B);
         SpaceUtils::upWindToFaces(FluxR_norm, FluxR_cc, CspeedR_norm, grid_box, m_advScheme_B);
      }

      // compute total flux on faces by adding left and right together
      VB_norm.copy(FluxR_norm,VB_norm.box());
      VB_norm += FluxL_norm;
      a_E_covar[dit][0].copy(VB_norm[1]);
      a_E_covar[dit][1].copy(VB_norm[0]);
      a_E_covar[dit][1].negate();
   } 

   // rebase all dummys
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_oneComp[dit].setVal(1./0.);
      m_dummyFArray_spaceDim[dit].setVal(1./0.);
   }


}

void IdealMhdOp::computeIdealEatNodes( LevelData<NodeFArrayBox>&  a_E_covar,
                                 const LevelData<FArrayBox>&      a_V_phys_cc,
                                 const LevelData<FluxBox>&        a_B_contra,
                                 const LevelData<FArrayBox>&      a_Cspeed_cc ) const
{
   CH_TIME("IdealMhdOp::computeIdealEatNodes() ");
   CH_assert(a_E_covar.nComp() == 1);
   CH_assert(a_B_contra.nComp() == 1);
   CH_assert(a_V_phys_cc.nComp() == SpaceDim);
   CH_assert(a_Cspeed_cc.nComp() == SpaceDim);
  
   //
   // compute covar E0 = -VxB at nodes:
   // E_i = Ja*(B^jV^k - B^kV^j)
   //
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
  
   // multiply contravariant B on faces by Jacobian and convert
   // physical velocity to contravariatn
   LevelData<FArrayBox>& V_contra_cc = m_dummyFArray_spaceDim;
   LevelData<FluxBox>& JaB_contra_cf(m_dummyFlux_oneComp);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      JaB_contra_cf[dit].copy(a_B_contra[dit]);
      const Box& cell_box(a_V_phys_cc[dit].box());
      SpaceUtils::copy(V_contra_cc[dit],a_V_phys_cc[dit],cell_box);
   } 
   m_geometry.multJonFaces(JaB_contra_cf);
   m_geometry.convertPhysToContravar(V_contra_cc,0);
 
   // compute constrained covariant E on cell edges
   SpaceUtils::constrainedCrossProductOnNodes( a_E_covar, JaB_contra_cf,
                                               V_contra_cc, a_Cspeed_cc, m_advScheme_B );
   
   // rebase dummies
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      V_contra_cc[dit].setVal(1.0/0.0);
      JaB_contra_cf[dit].setVal(1.0/0.0);
   } 

}

void IdealMhdOp::fillGhostCells( FluidSpecies&  a_species_phys,
                           const double         a_time )
{
   CH_TIME("IdealMhdOp::fillGhostCells()");
   
   for (int n=0; n<a_species_phys.num_cell_vars(); ++n) {
      LevelData<FArrayBox>& cellVar( a_species_phys.cell_var(n) );
      const std::string& this_var_name = a_species_phys.cell_var_name(n);
      if(this_var_name=="momentumDensity") m_geometry.convertPhysToContravar(cellVar,0);
      m_geometry.fillInternalGhosts( cellVar );
      m_fluid_bc.at(n)->apply( a_species_phys, a_time );
      if(this_var_name=="momentumDensity") m_geometry.convertPhysToContravar(cellVar,1);
   }
   
   for (int n=0; n<a_species_phys.num_face_vars(); ++n) {
      const std::string& this_var_name = a_species_phys.face_var_name(n);
      CH_assert(m_fluid_face_var_bc.at(n)->isForVariable(this_var_name));
      LevelData<FluxBox>& faceVar( a_species_phys.face_var(n) );
      SpaceUtils::exchangeFluxBox( faceVar ); 
      m_fluid_face_var_bc.at(n)->applyFluxBC( a_species_phys, faceVar, a_time );
   }

}

void IdealMhdOp::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "verbosity", m_verbosity);

   a_pp.query( "implicitViscosity", m_implicitViscosity);
   a_pp.query( "implicitMhd", m_implicitMhd);
   
   a_pp.query( "energy_is_thermal", m_energy_is_thermal);
   if(!m_energy_is_thermal) a_pp.query( "energy_is_thermal_2", m_energy_is_thermal_2);
   a_pp.query("isothermal", m_isothermal);
   if(m_isothermal) {
      m_energy_is_thermal = false;
      m_energy_is_thermal_2 = false;
   }
   a_pp.query("enforceFluxBCs", m_enforceFluxBCs);
   a_pp.query( "use_perturbative_form", m_use_perturbative_form);
   a_pp.query( "free_stream_preserving", m_free_stream_preserving);

   if (a_pp.contains("advScheme"))  a_pp.get("advScheme", m_advScheme );
   m_advScheme_B = m_advScheme;
   a_pp.query("advScheme_B", m_advScheme_B );
   
   if(m_advScheme=="ZIP") m_advScheme_B = "c2";
   
   a_pp.query( "magP_inCspeed", m_magP);
   a_pp.query( "magPforE0_inCspeed", m_magPforE0);
   a_pp.query( "thermPforE0_inCspeed", m_thermPforE0);
   m_magPforE0 *= m_magP;
   a_pp.query( "blankN", m_blankN);
   a_pp.query( "blankB", m_blankB);
   
   if (a_pp.contains("relative_density_floor")) {
      a_pp.get("relative_density_floor", m_relative_density_floor );
   }
   
   if (a_pp.contains("tau")) {
      a_pp.get("tau", m_tau );
   }
   
   if (a_pp.contains("etaMin")) {
      a_pp.get("etaMin", m_etaMin );
   }

   a_pp.query("initializeBfromVectorPotential", m_initializeBfromVectorPotential);
   a_pp.query("initializeProfilesFromPSI", m_initializeProfilesFromPSI);
   a_pp.query("initializeBackgroundPlusPert", m_initializeBackgroundPlusPert);
   a_pp.query("addEnergyPerturbation", m_addEnergyPerturbation);
   a_pp.query("perturbDensity", m_perturbDensity);
   a_pp.query("perturbBtor", m_perturbBtor);
   if(m_addEnergyPerturbation && m_isothermal) CH_assert(m_perturbDensity);
   a_pp.query("RBTOR_FACT", m_rbtor_fact);
   a_pp.query("BPOL_FACT", m_bpol_fact);

   int init_profiles = m_initializeBfromVectorPotential*m_initializeProfilesFromPSI;
   init_profiles *= m_initializeBackgroundPlusPert;
   CH_assert(init_profiles==0);
  
}


void IdealMhdOp::printParameters()
{
   if (procID()==0) {
      std::cout << "IdealMhdOp parameters:" << std::endl;
      std::cout << " advScheme    =  " << m_advScheme << std::endl;
      std::cout << " advScheme_B  =  " << m_advScheme_B << std::endl;
      std::cout << " blankN  =  " << m_blankN << std::endl;
      std::cout << " blankB  =  " << m_blankB << std::endl;
      std::cout << " magP_inCspeed =  " << m_magP << std::endl;
      std::cout << " magPforE0_inCspeed =  " << m_magPforE0 << std::endl;
      std::cout << " thermPforE0_inCspeed =  " << m_thermPforE0 << std::endl;
      std::cout << " rel density floor  =  " << m_relative_density_floor << std::endl;
      std::cout << " tau     =  " << m_tau << std::endl;
      std::cout << " etaMin  =  " << m_etaMin << std::endl;
      std::cout << " energy_is_thermal  =  " << m_energy_is_thermal << std::endl;
      std::cout << " energy_is_thermal_2  =  " << m_energy_is_thermal_2 << std::endl;
      std::cout << " isothermal = " << m_isothermal << std::endl;
      std::cout << " enforceFluxBCs = " << m_enforceFluxBCs << std::endl;
      std::cout << " use_perturbative_form  =  " << m_use_perturbative_form << std::endl;
      std::cout << " free_stream_preserving  =  " << m_free_stream_preserving << std::endl;
      std::cout << " implicitMhd  =  " << m_implicitMhd << std::endl;
      std::cout << " implicitVisosity  =  " << m_implicitViscosity << std::endl;
      std::cout << " initializeBfromVectorPotential = " << m_initializeBfromVectorPotential << std::endl;
      std::cout << " initializeProfilesFromPSI = " << m_initializeProfilesFromPSI << std::endl;
      std::cout << " initializeBackgroundPlusPert = " << m_initializeBackgroundPlusPert << std::endl;
      std::cout << " addEnergyPerturbation = " << m_addEnergyPerturbation << std::endl;
      std::cout << " rbtor_fact = " << m_rbtor_fact << std::endl;
      std::cout << " bpol_fact = " << m_bpol_fact << std::endl;
      std::cout << " perturbBtor = " << m_perturbBtor << std::endl;
      std::cout << " perturbDensity = " << m_perturbDensity << std::endl << std::endl;
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
   // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/N + B^2/N); i = q0, q1, q2

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
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // get physical values for m_species_name and calculate freezing speed
   // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/rho); i = q0, q1, q2
   
   double gamma;
   LevelData<FArrayBox>& velocity = m_dummyFArray_spaceDim;
   for (int species(0); species<a_fluid_comp.size(); species++) {
      const FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_comp[species])) );
      const std::string species_name( fluid_species.name() );
      if(species_name==m_species_name) {
         gamma = fluid_species.m_gamma;
         fluid_species.massDensity(m_rhoDen_cc);
         fluid_species.velocity(velocity);  // in-plane velocity vector
         if(m_energy_is_thermal || m_energy_is_thermal_2) {
            for (DataIterator dit(grids); dit.ok(); ++dit) {
               m_pressure[dit].copy(m_eneDen_cc[dit],0,0,1);
               m_pressure[dit].mult(gamma-1.0);
            }         
         }
         else if(m_isothermal) {
            for (DataIterator dit(grids); dit.ok(); ++dit) {
               m_pressure[dit].copy(m_rhoDen_cc[dit],0,0,1);
               m_pressure[dit].mult(m_T0[dit]);
               m_pressure[dit].mult(2.0);
            }
         }
         else {
            fluid_species.pressure(m_pressure);
         }

         for (DataIterator dit(grids); dit.ok(); ++dit) {
            m_temperature[dit].copy(m_pressure[dit]);
            m_temperature[dit].divide(m_rhoDen_cc[dit],0,0,1);
            m_temperature[dit].divide(2.0);
         }

         fluid_species.Bpressure(m_Bpressure,true);
         
         // Calls for physical variables abve are from mapped state
         // vector (has not been convertToPhysical )
         // and so need to divide by J
         m_geometry.divideJonValid(m_rhoDen_cc);
         m_geometry.divideJonValid(m_pressure);
         
         break;
      }
   }

   //setMappedCspeed( gamma, 1, m_magP );
   setMappedCspeed( gamma, 1, 1 );
   setCourantTimeStep(m_Cspeed_cc);   
   
   return m_courant_time_step;

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
   double twoPi = 2.0*Pi;
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
   m_By_phys_cc.define(a_grids, 1, a_ghostVect);
   //
   m_By_contr_cc.define(a_grids, 1, a_ghostVect);
   m_By_covar_cc.define(a_grids, 1, a_ghostVect);
   m_B_covar_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_B0_phys.define(a_grids, SpaceDim, a_ghostVect);
   m_B1_phys.define(a_grids, SpaceDim, a_ghostVect);
   //
   m_pressure.define(a_grids, 1+m_CGL, a_ghostVect);
   m_temperature.define(a_grids, 1+m_CGL, a_ghostVect);
   m_P0.define(a_grids, 1+m_CGL, a_ghostVect);
   m_N0.define(a_grids, 1, a_ghostVect);
   m_T0.define(a_grids, 1+m_CGL, a_ghostVect);
   m_deltaN.define(a_grids, 1, a_ghostVect);
   m_deltaT.define(a_grids, 1+m_CGL, a_ghostVect);
   m_Bpressure.define(a_grids, 1, a_ghostVect);
   m_velocity.define(a_grids, SpaceDim, a_ghostVect);
   m_velocity_cf.define(a_grids, SpaceDim, a_ghostVect);
   m_velocity_virtual.define(a_grids, 1, a_ghostVect);
   //
   m_velocity_contra.define(a_grids, SpaceDim, a_ghostVect);
   m_momDen_contra.define(a_grids, SpaceDim, a_ghostVect);
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
   m_JaJ00.define(a_grids, SpaceDim, a_ghostVect);
   m_JaJ01.define(a_grids, SpaceDim, a_ghostVect);
   m_Apar_cc.define(a_grids, 1, a_ghostVect);
   m_E0_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_curlE0_virtual.define(a_grids, 1, IntVect::Zero);
   m_curlE0_cf.define(a_grids, 1, IntVect::Zero);
   m_B_contr_cf.define(a_grids, 1, a_ghostVect);
   m_JaUdotGradP.define(a_grids, 1, IntVect::Zero); // Jacobian*U\cdot\grad(P)
   m_JaPdivU.define(a_grids, 1, IntVect::Zero); // Jacobian*U\cdot\grad(P)
   m_JaJdotE.define(a_grids, 1, IntVect::Zero); // Jacobian*J\cdotE
   m_JaJcrossB.define(a_grids, SpaceDim, a_ghostVect); // Jacobian*J\timesB
   m_JaJcrossB_10.define(a_grids, SpaceDim, a_ghostVect); // Jacobian*J\timesB
   //
   m_JaJ0_virtual_cc.define(a_grids,1,a_ghostVect);
   m_JaJcrossB_virtual.define(a_grids,1,a_ghostVect);
   m_E0_virtual_cc.define(a_grids,1,a_ghostVect);
   m_E0_virtual_nc.define(a_grids,1,a_ghostVect);
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_JaJ0_virtual_cc[dit].setVal(0.0);
      m_JaJcrossB_virtual[dit].setVal(0.0);
      m_velocity_virtual[dit].setVal(0.0);
      m_E0_virtual_cc[dit].setVal(0.0);
   }
   //
   m_NodeBC_zeros.define(a_grids, 1, a_ghostVect);
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
      m_NodeBC_zeros[dit].setVal(0.0);
      m_EdgeBC_zeros[dit].setVal(0.0);
      m_FluxBC_zeros[dit].setVal(0.0);
   }
   
   m_dummyDiv.define(a_grids, 1, IntVect::Zero);
   m_dummyFArray_oneComp.define(a_grids, 1, a_ghostVect);
   m_dummyFArray_spaceDim.define(a_grids, SpaceDim, a_ghostVect);
   m_dummyDiv_mom.define(a_grids, SpaceDim, IntVect::Zero);
   m_forceDensity.define(a_grids, SpaceDim, IntVect::Zero);
   m_rhs_rho.define(a_grids, 1, IntVect::Zero);
   m_dummyFlux_oneComp.define(a_grids, 1, a_ghostVect);
   m_dummyFlux_spaceDim.define(a_grids, SpaceDim, a_ghostVect);
   m_dummyEdge_oneComp.define(a_grids, 1, a_ghostVect);
   m_dummyNode_oneComp.define(a_grids, 1, a_ghostVect);
   
   m_PSI.define(a_grids, 1, a_ghostVect);

}

void IdealMhdOp::setCellCenterValues( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setCellCenterValues()");
   
   const LevelData<FArrayBox>& soln_rhoDen( a_soln_fluid.cell_var(0) );
   const LevelData<FArrayBox>& soln_momDen( a_soln_fluid.cell_var("momentumDensity") );
   const LevelData<FArrayBox>& soln_eneDen( a_soln_fluid.cell_var("energyDensity") );
   
   // set values for MHD variables
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
         m_By_phys_cc[dit].copy( soln_magField_virtual[dit] );
         m_By_contr_cc[dit].copy( soln_magField_virtual[dit] );
         m_By_covar_cc[dit].copy( soln_magField_virtual[dit] );
         if(m_twoDaxisymm) m_By_contr_cc[dit].divide(m_g_y[dit], 0,0,1);
         if(m_twoDaxisymm) m_By_covar_cc[dit].mult(m_g_y[dit], 0,0,1);
      }
   }  
   if(a_soln_fluid.m_evolve_magneticField) {
      const LevelData<FluxBox>& soln_magField( a_soln_fluid.face_var("magneticField") );
      SpaceUtils::interpFaceVectorToCell(m_B_covar_cc,soln_magField,"c2");
      m_geometry.convertPhysToContravar(m_B_covar_cc,1);
      m_geometry.convertPhysToCovar(m_B_covar_cc,0);
   }
   
   a_soln_fluid.Bpressure(m_Bpressure);
   a_soln_fluid.velocity(m_velocity);
   if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
      a_soln_fluid.velocity_virtual(m_velocity_virtual);
   }
 
   // set physical derived variables for fluid species
   const double gamma = a_soln_fluid.m_gamma;
   if(m_energy_is_thermal || m_energy_is_thermal_2) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_pressure[dit].copy(m_eneDen_cc[dit],0,0,1);
         m_pressure[dit].mult(gamma-1.0);
      }         
   }
   else if(m_isothermal) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_pressure[dit].copy(m_rhoDen_cc[dit],0,0,1);
         m_pressure[dit].mult(m_T0[dit],0,0,1);
         m_pressure[dit].mult(2.0);
      }
   }
   else {
      a_soln_fluid.pressure(m_pressure);
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_temperature[dit].copy(m_pressure[dit]);
      m_temperature[dit].divide(m_rhoDen_cc[dit],0,0,1);
      m_temperature[dit].divide(2.0);
   }
 
   if(m_energy_is_thermal) { // compute Ja*U\cdotGradP_phys
      LevelData<FArrayBox>& gradP_phys = m_dummyFArray_spaceDim;
      m_geometry.computeMappedGradient( m_pressure, gradP_phys, 2);
      m_geometry.convertPhysToCovar(gradP_phys,1);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
               FArrayBox& JaUdotGradP_on_patch = m_JaUdotGradP[dit];
         const FArrayBox& gradP_on_patch       = gradP_phys[dit];
         const FArrayBox& U_on_patch           = m_velocity[dit];

         const Box& cellbox = JaUdotGradP_on_patch.box();
         FORT_EVAL_DOT_PRODUCT( CHF_BOX(cellbox),
                                CHF_CONST_FRA(gradP_on_patch),
                                CHF_CONST_FRA(U_on_patch),
                                CHF_FRA1(JaUdotGradP_on_patch,0) );
      }
      m_geometry.multJonValid(m_JaUdotGradP);
   }
   else if(m_energy_is_thermal_2) { // compute Ja*P*Div(U)
      LevelData<FluxBox> velocity_norm(grids, 1, IntVect::Unit);
      SpaceUtils::upWindToFaces(m_velocity_cf, m_velocity, m_velocity_cf, "c2");
      m_geometry.computeMetricTermProductAverage( velocity_norm, m_velocity_cf, 0 );
      m_geometry.mappedGridDivergenceFromFluxNorms( velocity_norm, m_JaPdivU );
      for (DataIterator dit(m_JaPdivU.dataIterator()); dit.ok(); ++dit) {
         const MagBlockCoordSys& coords( m_geometry.getBlockCoordSys( grids[dit] ) );
         m_JaPdivU[dit].mult( 1.0/coords.getMappedCellVolume() );
         m_JaPdivU[dit].mult(m_pressure[dit]);
      }

   }

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

      if(m_use_perturbative_form) {
         const FArrayBox& P0 = m_P0[dit];
         mxFlux_on_patch.minus(P0,0,0,1);
         myFlux_on_patch.minus(P0,0,1,1);
         mzFlux_on_patch.minus(P0,0,2,1);
      }
 
      if(m_energy_is_thermal_2) {
         for (int n=0; n<SpaceDim; ++n) {
            m_enFlux_cc[dit].copy( eneDen_on_patch,0,n,1 );
            m_enFlux_cc[dit].mult( m_velocity[dit],n,n,1);
         }
      }
      
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
   
   //  set the magnetic field source terms for momentum density and 
   //  energy density equations:
   //  momentum source: Ja * J_phys x B_phys
   //  energy source:   Ja * J_phys dot E_phys

   switch(a_soln_fluid.m_evolveB_type) {
      case FluidSpecies::EVOLVE_B_TYPE::NONE :
          break;
      case FluidSpecies::EVOLVE_B_TYPE::TWOD_VIRTUAL_B :
          setMagneticFieldTerms_2Dvirtual( a_soln_fluid );
          break;
      case FluidSpecies::EVOLVE_B_TYPE::TWOD_INPLANE_B :
          setMagneticFieldTerms_2DinPlane( a_soln_fluid );
          break;
      case FluidSpecies::EVOLVE_B_TYPE::TWOD_FULL_B :
          setMagneticFieldTerms_2Dfull( a_soln_fluid );
          break;
      case FluidSpecies::EVOLVE_B_TYPE::THREED :
          setMagneticFieldTerms_3D( a_soln_fluid );
          break;
      default :
          cout << "m_evolveB_type = " << a_soln_fluid.m_evolveB_type << endl;
          MayDay::Error("Invalid evolve B type in setMagneticFieldTerms() ");
   }

}   

void IdealMhdOp::setMagneticFieldTerms_2Dvirtual( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setMagneticFieldTerms_2Dvirtual()");
   
   //  E=BxV:   x: (By*Vz        )
   //           y: (             ) (virtual)
   //           z: (      - By*Vx)
   //  J=curlB: x: (       - dBy/dz)
   //           y: (               ) (virtual)
   //           z: (dBy/dx         )
   //  JxB:     x: (      - Jz*By)
   //           y: (             ) (virtual)
   //           z: (Jx*By        )
   //  JdotE:  (Jx*Ex + Jz*Ez)
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   const LevelData<FArrayBox>& B_virtual( a_soln_fluid.cell_var("magneticField_virtual") );

   // compute E = BxV
   computeCrossProduct( m_E0_cc, B_virtual, m_velocity, 2 );

   // compute Ja*JxB      
   m_geometry.mappedGridCurlofVirtComp(m_JaJ0_cc,m_By_covar_cc); // Ja*curl(By)\cdot g^perp
   m_geometry.convertPhysToContravar(m_JaJ0_cc,1);          // Ja*curlB_phys at cell center
   computeCrossProduct( m_JaJcrossB, m_JaJ0_cc, B_virtual, 3 );

   // compute Ja*JdotE
   for (DataIterator dit(grids); dit.ok(); ++dit) {
         
      const FArrayBox& JaJ_on_patch     = m_JaJ0_cc[dit];
      const FArrayBox& E_on_patch       = m_E0_cc[dit];
            FArrayBox& JaJdotE_on_patch = m_JaJdotE[dit];

      const Box& cellbox = JaJdotE_on_patch.box();
      FORT_EVAL_DOT_PRODUCT( CHF_BOX(cellbox),
                             CHF_CONST_FRA(JaJ_on_patch),
                             CHF_CONST_FRA(E_on_patch),
                             CHF_FRA1(JaJdotE_on_patch,0) );
         
      // reset dummy containers to infty
      m_dummyFArray_spaceDim[dit].setVal(1./0.);
      
   }

}

void IdealMhdOp::setMagneticFieldTerms_2DinPlane( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setMagneticFieldTerms_2DinPlane()");
   
   //  E=BxV:   x: (      - Bz*Vy)
   //           y: (Bz*Vx - Bx*Vz) (virtual)
   //           z: (Bx*Vy -      )
   //  J=curlB: x: (               )
   //           y: (dBx/dz - dBz/dx) (virtual)
   //           z: (               )
   //  JxB:     x: (Jy*Bz        )
   //           y: (             ) (virtual)
   //           z: (      - Jy*Bx)
   //  JdotE:  (Jy*Ey)
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   const LevelData<FluxBox>& soln_magField_cf( a_soln_fluid.face_var("magneticField") );
 
   // compute physical B at cell center
   LevelData<FArrayBox>& B_inPlane = m_dummyFArray_spaceDim;
   SpaceUtils::interpFaceVectorToCell(B_inPlane,soln_magField_cf,"c2");
   m_geometry.convertPhysToContravar(B_inPlane,1);
      
   // compute E = BxV
   computeCrossProduct( m_E0_cc, B_inPlane, m_velocity_virtual, 3 );
   computeCrossProduct( m_E0_virtual_cc, B_inPlane, m_velocity, 1 );
   
   // compute Ja*JxB
   m_geometry.mappedGridCurl( m_JaJ0_virtual_cc, m_B_covar_cc ); // Ja*curl(B)_phys
   computeCrossProduct( m_JaJcrossB, m_JaJ0_virtual_cc, B_inPlane, 2 );
   
   // compute Ja*JdotE
   for (DataIterator dit(grids); dit.ok(); ++dit) {
         
            FArrayBox& JaJdotE_on_patch = m_JaJdotE[dit];
      const FArrayBox& JaJv_on_patch = m_JaJ0_virtual_cc[dit];
      const FArrayBox& Ev_on_patch   = m_E0_virtual_cc[dit];
         
      const Box& cellbox = JaJdotE_on_patch.box();
      Real JaJ0v, E0v, JaJdotE;
      for (BoxIterator bit(cellbox); bit.ok(); ++bit) {
          IntVect ic = bit();
          JaJ0v = JaJv_on_patch.get(ic,0);
          E0v   = Ev_on_patch.get(ic,0);
          JaJdotE = JaJ0v*E0v;
          JaJdotE_on_patch.set(ic,0,JaJdotE);
      }

      // reset dummy containers to infty
      m_dummyFArray_spaceDim[dit].setVal(1./0.);

   }

}

void IdealMhdOp::setMagneticFieldTerms_2Dfull( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setMagneticFieldTerms_2Dfull()");
   
   //  E=BxV:   x: (By*Vz - Bz*Vy)
   //           y: (Bz*Vx - Bx*Vz) (virtual)
   //           z: (Bx*Vy - By*Vx)
   //  J=curlB: x: (       - dBy/dz)
   //           y: (dBx/dz - dBz/dx) (virtual)
   //           z: (dBy/dx         )
   //  JxB:     x: (Jy*Bz - Jz*By)
   //           y: (Jz*Bx - Jx*Bz) (virtual)
   //           z: (Jx*By - Jy*Bx)
   //  JdotE:  (Jx*Ex + Jz*Ez + Jy*Ey)
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   const LevelData<FluxBox>& soln_magField_cf( a_soln_fluid.face_var("magneticField") );
   const LevelData<FArrayBox>& B_virtual( a_soln_fluid.cell_var("magneticField_virtual") );
      
   // compute B at cell center
   LevelData<FArrayBox>& B_inPlane = m_dummyFArray_spaceDim;
   SpaceUtils::interpFaceVectorToCell(B_inPlane,soln_magField_cf,"c2");
   m_geometry.convertPhysToContravar(B_inPlane,1);
      
   // compute E = BxV
   computeCrossProduct( m_E0_cc, B_inPlane, m_velocity,
                        B_virtual, m_velocity_virtual );
   computeCrossProduct( m_E0_virtual_cc, B_inPlane, m_velocity, 1);

   // compute Ja*JxB
   m_geometry.mappedGridCurlofVirtComp(m_JaJ0_cc,m_By_covar_cc); // Ja*curl(By)\cdot g^perp
   m_geometry.convertPhysToContravar(m_JaJ0_cc,1);          // Ja*curlB_phys at cell center
   m_geometry.mappedGridCurl( m_JaJ0_virtual_cc, m_B_covar_cc ); // Ja*curl(B)_phys
   computeCrossProduct( m_JaJcrossB, m_JaJ0_cc, B_inPlane, 
                        m_JaJ0_virtual_cc, B_virtual );
   computeCrossProduct( m_JaJcrossB_virtual, m_JaJ0_cc, B_inPlane, 1 );
   
   // compute Ja*JdotE
   for (DataIterator dit(grids); dit.ok(); ++dit) {
         
      const FArrayBox& JaJ_on_patch     = m_JaJ0_cc[dit];
      const FArrayBox& E_on_patch       = m_E0_cc[dit];
            FArrayBox& JaJdotE_on_patch = m_JaJdotE[dit];

      const Box& cellbox = JaJdotE_on_patch.box();
      FORT_EVAL_DOT_PRODUCT( CHF_BOX(cellbox),
                             CHF_CONST_FRA(JaJ_on_patch),
                             CHF_CONST_FRA(E_on_patch),
                             CHF_FRA1(JaJdotE_on_patch,0) );
     
      // add virtual contribution to Ja*JdotE 
      const FArrayBox& JaJv_on_patch = m_JaJ0_virtual_cc[dit];
      const FArrayBox& Ev_on_patch   = m_E0_virtual_cc[dit];
      
      Real JaJ0v, E0v, JaJdotE;
      for (BoxIterator bit(cellbox); bit.ok(); ++bit) {
          IntVect ic = bit();
          JaJ0v = JaJv_on_patch.get(ic,0);
          E0v   = Ev_on_patch.get(ic,0);
          JaJdotE  = JaJdotE_on_patch.get(ic,0);
          JaJdotE += JaJ0v*E0v;
          JaJdotE_on_patch.set(ic,0,JaJdotE);
      }   
         
      // reset dummy containers to infty
      m_dummyFArray_spaceDim[dit].setVal(1./0.);
      
   }
   
}

void IdealMhdOp::setMagneticFieldTerms_3D( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setMagneticFieldTerms_3D()");
   
   //  E=BxV:   x: (By*Vz - Bz*Vy)
   //           y: (Bz*Vx - Bx*Vz)
   //           z: (Bx*Vy - By*Vx)
   //  J=curlB: x: (dBz/dy - dBy/dz)
   //           y: (dBx/dz - dBz/dx)
   //           z: (dBy/dx - dBx/dy)
   //  JxB:     x: (Jy*Bz - Jz*By)
   //           y: (Jz*Bx - Jx*Bz)
   //           z: (Jx*By - Jy*Bx)
   //  JdotE:  (Jx*Ex + Jy*Ey + Jz*Ez)
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   const LevelData<FluxBox>& soln_magField_cf( a_soln_fluid.face_var("magneticField") );
 
   // compute physical B at cell center
   LevelData<FArrayBox>& B_inPlane = m_dummyFArray_spaceDim;
   SpaceUtils::interpFaceVectorToCell(B_inPlane,soln_magField_cf,"c2");
   m_geometry.convertPhysToContravar(B_inPlane,1);
      
   // compute E = BxV
   computeCrossProduct( m_E0_cc, B_inPlane, m_velocity, 0 );
   
   // compute Ja*JxB
   //LevelData<FluxBox>& Bvec_covar = m_dummyFlux_spaceDim;
   //SpaceUtils::interpFacesToFaces(Bvec_covar, soln_magField_cf);
   //m_geometry.convertContravarToCovar(Bvec_covar,0);
   //m_geometry.mappedGridCurl( m_JaJ0_cc, Bvec_covar ); // Ja*curl(B)_contra
   m_geometry.mappedGridCurl( m_JaJ0_cc, m_B_covar_cc ); // Ja*curl(B)_contra
   m_geometry.convertPhysToContravar(m_JaJ0_cc,1);

   if(m_use_perturbative_form) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_JaJ01[dit].copy(m_JaJ0_cc[dit]);
         m_JaJ01[dit].minus(m_JaJ00[dit]);
         m_B1_phys[dit].copy(B_inPlane[dit]);
         m_B1_phys[dit].minus(m_B0_phys[dit]);
      }
      computeCrossProduct( m_JaJcrossB, m_JaJ00, m_B1_phys, 0 );
      computeCrossProduct( m_JaJcrossB_10, m_JaJ01, m_B0_phys, 0 );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_JaJcrossB[dit].plus(m_JaJcrossB_10[dit]);
      }
   }  
   else {
      computeCrossProduct( m_JaJcrossB, m_JaJ0_cc, B_inPlane, 0 );
   }

   // compute Ja*JdotE
   for (DataIterator dit(grids); dit.ok(); ++dit) {
         
      const FArrayBox& JaJ_on_patch     = m_JaJ0_cc[dit];
      const FArrayBox& E_on_patch       = m_E0_cc[dit];
            FArrayBox& JaJdotE_on_patch = m_JaJdotE[dit];

      const Box& cellbox = JaJdotE_on_patch.box();
      FORT_EVAL_DOT_PRODUCT( CHF_BOX(cellbox),
                             CHF_CONST_FRA(JaJ_on_patch),
                             CHF_CONST_FRA(E_on_patch),
                             CHF_FRA1(JaJdotE_on_patch,0) );
         
      // reset dummy containers to infty
      m_dummyFArray_spaceDim[dit].setVal(1./0.);
      m_dummyFlux_spaceDim[dit].setVal(1./0.);
      
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
#if CH_SPACEDIM==3
      m_momFlux_norm[dit].copy( m_mzFlux_norm[dit],0,2,1 );
#endif
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

void IdealMhdOp::setFaceCenteredFluxes_freestreamPreserving( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setFaceCenteredFluxes_freestreamPreserving()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   LevelData<FluxBox>& thisFlux_phys_fc = m_dummyFlux_spaceDim;
   LevelData<FArrayBox>& thisf = m_dummyFArray_oneComp;
  
   // compute normal flux for mass density
   upwindFluxVector(thisFlux_phys_fc,m_rhoFlux_cc,m_rhoDen_cc,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   m_geometry.computeMetricTermProductAverage(m_rhoFlux_norm, thisFlux_phys_fc, false);

   // compute normal flux for momentum density vector 
   for (int n=0; n<SpaceDim; ++n) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         thisf[dit].copy(m_momDen_cc[dit],n,0,1);
      }
      if(n==0) {
         upwindFluxVector(thisFlux_phys_fc,m_mxFlux_cc,thisf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
         m_geometry.computeMetricTermProductAverage(m_mxFlux_norm, thisFlux_phys_fc, false);
      }
      if(n==1) {
         upwindFluxVector(thisFlux_phys_fc,m_myFlux_cc,thisf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
         m_geometry.computeMetricTermProductAverage(m_myFlux_norm, thisFlux_phys_fc, false);
      }
      if(n==2) {
         upwindFluxVector(thisFlux_phys_fc,m_mzFlux_cc,thisf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
         m_geometry.computeMetricTermProductAverage(m_mzFlux_norm, thisFlux_phys_fc, false);
      }
   }
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_momFlux_norm[dit].copy( m_mxFlux_norm[dit],0,0,1 );
      m_momFlux_norm[dit].copy( m_myFlux_norm[dit],0,1,1 );
#if CH_SPACEDIM==3
      m_momFlux_norm[dit].copy( m_mzFlux_norm[dit],0,2,1 );
#endif
   }
   
   // compute normal flux for ion energy density
   upwindFluxVector(thisFlux_phys_fc,m_enFlux_cc,m_eneDen_cc,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   m_geometry.computeMetricTermProductAverage(m_enFlux_norm, thisFlux_phys_fc, false);

}

void IdealMhdOp::setFaceCenteredFluxes_ZIP( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("IdealMhdOp::setFaceCenteredFluxes_ZIP()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
 
   const double dummy_time = 0.0;

   // compute contravariant velocity and momentum density
   for (DataIterator dit(grids); dit.ok(); ++dit) {
       m_velocity_contra[dit].copy(m_velocity[dit]);
       m_momDen_contra[dit].copy(m_momDen_cc[dit]);
   }
   m_geometry.convertPhysToContravar(m_velocity_contra, 0);
   m_geometry.convertPhysToContravar(m_momDen_contra, 0);
 
   // compute contraviant fluxes on cell faces for density and energy density
   SpaceUtils::interpToFacesZIP( m_rhoFlux_norm, m_rhoDen_cc, m_velocity_contra );
   
   if(m_energy_is_thermal_2) {
      SpaceUtils::interpToFacesZIP( m_enFlux_norm, m_eneDen_cc, m_velocity_contra );
   }
   else if(m_isothermal) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_enFlux_norm[dit].setVal(0.0);
      }
   }
   else {
      LevelData<FArrayBox>& eneP  = m_dummyFArray_oneComp;
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         eneP[dit].copy(m_eneDen_cc[dit]);
         eneP[dit].plus(m_pressure[dit]);
      }
      SpaceUtils::interpToFacesZIP( m_enFlux_norm, eneP, m_velocity_contra );
   }
 
   // compute flux for momentum density vector 
   //if(m_use_perturbative_form) {
   //   for (DataIterator dit(grids); dit.ok(); ++dit) {
   //      m_mxFlux_norm[dit].setVal(0.0);
   //      m_myFlux_norm[dit].setVal(0.0);
   //      m_mzFlux_norm[dit].setVal(0.0);
   //   }
   //} 
   //else {
      for (int dir=0; dir<SpaceDim; ++dir) {
         if(dir==0) SpaceUtils::interpToFacesZIP( m_mxFlux_norm, m_velocity, m_momDen_cc, 
                                                  m_velocity_contra, m_momDen_contra, dir );
         if(dir==1) SpaceUtils::interpToFacesZIP( m_myFlux_norm, m_velocity, m_momDen_cc, 
                                                  m_velocity_contra, m_momDen_contra, dir );
         if(dir==2) SpaceUtils::interpToFacesZIP( m_mzFlux_norm, m_velocity, m_momDen_cc, 
                                                  m_velocity_contra, m_momDen_contra, dir );
      }

      if(m_enforceFluxBCs) {
         m_fluid_bc.at(0)->setFluxBC( a_soln_fluid, m_rhoFlux_norm, m_FluxBC_zeros, dummy_time );
         m_fluid_bc.at(0)->setFluxBC( a_soln_fluid, m_enFlux_norm, m_FluxBC_zeros, dummy_time );
         m_fluid_bc.at(0)->setFluxBC( a_soln_fluid, m_mxFlux_norm, m_FluxBC_zeros, dummy_time );
         m_fluid_bc.at(0)->setFluxBC( a_soln_fluid, m_myFlux_norm, m_FluxBC_zeros, dummy_time );
#if CH_SPACEDIM==3
         m_fluid_bc.at(0)->setFluxBC( a_soln_fluid, m_mzFlux_norm, m_FluxBC_zeros, dummy_time );
#endif
      }
   //}

   // interpolate pressure to faces
   LevelData<FluxBox>& pressure_fc = m_dummyFlux_oneComp;
   if(m_use_perturbative_form) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_deltaN[dit].copy(m_rhoDen_cc[dit]);
         m_deltaN[dit].minus(m_N0[dit]);
         m_deltaT[dit].copy(m_temperature[dit]);
         m_deltaT[dit].minus(m_T0[dit]);
      }
      LevelData<FluxBox> pressure_fc2(grids,1,m_rhoDen_cc.ghostVect());
      SpaceUtils::interpToFacesZIP( pressure_fc, m_N0, m_deltaT );
      SpaceUtils::interpToFacesZIP( pressure_fc2, m_deltaN, m_T0 );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for(int dir=0; dir<SpaceDim; dir++) {
            pressure_fc[dit][dir].plus(pressure_fc2[dit][dir]);
         }
      }
   }
   else {
      SpaceUtils::interpToFacesZIP( pressure_fc, m_rhoDen_cc, m_temperature );
   }
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for(int dir=0; dir<SpaceDim; dir++) pressure_fc[dit][dir].mult(2.0);
   }
   LevelData<FluxBox>& Pcontra = m_dummyFlux_spaceDim;

   // convert pressure flux to contravariant for x-momentum
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      Pcontra[dit].setVal(0.0);
      for (int dir=0; dir<SpaceDim; dir++) {
         Pcontra[dit][dir].copy(pressure_fc[dit][dir],0,0,1 );
      }
   }
   m_geometry.convertPhysToContravar(Pcontra, 0);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         m_mxFlux_norm[dit][dir].plus( Pcontra[dit][dir],dir,0,1 );
      }
   }
   
   // convert pressure flux to contravariant for y-momentum
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      Pcontra[dit].setVal(0.0);
      for (int dir=0; dir<SpaceDim; dir++) {
         Pcontra[dit][dir].copy(pressure_fc[dit][dir],0,1,1 );
      }
   }
   m_geometry.convertPhysToContravar(Pcontra, 0);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         m_myFlux_norm[dit][dir].plus( Pcontra[dit][dir],dir,0,1 );
      }
   }
   
#if CH_SPACEDIM==3
   // convert pressure flux to contravariant for z-momentum
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      Pcontra[dit].setVal(0.0);
      for (int dir=0; dir<SpaceDim; dir++) {
         Pcontra[dit][dir].copy(pressure_fc[dit][dir],0,2,1 );
      }
   }
   m_geometry.convertPhysToContravar(Pcontra, 0);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         m_mzFlux_norm[dit][dir].plus( Pcontra[dit][dir],dir,0,1 );
      }
   }
#endif

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         m_momFlux_norm[dit][dir].copy( m_mxFlux_norm[dit][dir],0,0,1 );
         m_momFlux_norm[dit][dir].copy( m_myFlux_norm[dit][dir],0,1,1 );
#if CH_SPACEDIM==3
         m_momFlux_norm[dit][dir].copy( m_mzFlux_norm[dit][dir],0,2,1 );
#endif
      }
   }
   
   // multiply fluxes by face area and Jacobian
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const RealVect& faceArea = coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; dir++) {
         m_rhoFlux_norm[dit][dir].mult(faceArea[dir]);
         m_momFlux_norm[dit][dir].mult(faceArea[dir]);
         m_enFlux_norm[dit][dir].mult(faceArea[dir]);
      }
   }
   m_geometry.multJonFaces( m_rhoFlux_norm );
   m_geometry.multJonFaces( m_momFlux_norm );
   m_geometry.multJonFaces( m_enFlux_norm );
   
}

void IdealMhdOp::setFaceCenteredFlux( LevelData<FluxBox>&    a_Flux_norm,
                                const LevelData<FArrayBox>&  a_Flux_phys_cc )
{
   CH_TIME("IdealMhdOp::setFaceCenteredFlux()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   int method = 2;
   if(m_free_stream_preserving) method = 0;

   if(method==0) {
      SpaceUtils::upWindToFaces( m_dummyFlux_spaceDim, a_Flux_phys_cc, m_dummyFlux_spaceDim, "c2" );
      m_geometry.applyAxisymmetricCorrection( m_dummyFlux_spaceDim );
      m_geometry.computeMetricTermProductAverage( a_Flux_norm, m_dummyFlux_spaceDim, 0 );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_dummyFlux_spaceDim[dit].setVal(1.0/0.0);
      }
   }
   
   // convert cell-center flux from physical to mapped and interp to faces 
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
   setFaceCenteredFlux(m_rhoFluxBC_norm, m_rhoFlux_cc);
   m_fluid_bc.at(0)->setFluxBC( a_soln_fluid, m_rhoFlux_norm, m_rhoFluxBC_norm, a_time );
   //m_fluid_bc.at(0)->setFluxBC( a_soln_fluid, m_rhoFlux_norm, m_FluxBC_zeros, a_time );
   
   // enforce flux BC for momentum density flux
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
   setFaceCenteredFlux(m_enFluxBC_norm, m_enFlux_cc);
   this_cbc = m_energyDensity_cbc;
   m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_enFlux_norm, m_enFluxBC_norm, a_time );
   //m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_enFlux_norm, m_FluxBC_zeros, a_time );
   
   // enforce flux BC for parallel energy density flux if using CGL
   if(m_CGL) {
      setFaceCenteredFlux(m_enParFluxBC_norm, m_enParFlux_cc);
      this_cbc = m_energyDensity_cbc;
      m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_enParFlux_norm, m_enParFluxBC_norm, a_time );
   }

   // enforce flux BC for virtual momentum density flux
   if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
      setFaceCenteredFlux(m_mvFluxBC_norm, m_mvFlux_cc);
      this_cbc = m_momentumDensity_virtual_cbc;
      m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_mvFlux_norm, m_mvFluxBC_norm, a_time );
   }
}

void IdealMhdOp::setCurlE( const FluidSpecies&  a_soln_fluid, 
                           const Real           a_time )
{
   CH_TIME("IdealMhdOp::computeCurlE()");
   const DisjointBoxLayout& grids( m_geometry.grids() );  
   
   if(a_soln_fluid.m_evolve_magneticField_virtual) {
      computeIdealEatEdges(m_E0_ce,m_velocity,m_By_contr_cc,m_Cspeed_cc);
      int this_cbc = m_magneticField_virtual_cbc;
      m_fluid_bc.at(this_cbc)->setEdgeBC( a_soln_fluid, m_E0_ce, m_EdgeBC_zeros, a_time );
      m_geometry.mappedGridCurl2D(m_curlE0_virtual,m_E0_ce); // takes covar E and returns Ja*curl(E)dot(g^y)
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         if(m_twoDaxisymm) m_curlE0_virtual[dit].mult(m_g_y[dit]); // g^y*g_y = 1 (Ja*curlE_phys)
      }
   }
   
   if(a_soln_fluid.m_evolve_magneticField) {
      const LevelData<FluxBox>& soln_magField_cf( a_soln_fluid.face_var("magneticField") );
#if CH_SPACEDIM==3
      computeIdealEatEdges(m_E0_ce,m_velocity,soln_magField_cf,m_Cspeed_cc);
      int this_fbc = m_magneticField_fbc;
      m_fluid_face_var_bc.at(this_fbc)->setEdgeBC( a_soln_fluid, m_E0_ce, m_EdgeBC_zeros, a_time );
      SpaceUtils::exchangeEdgeDataBox(m_E0_ce); // call exchange after setBC !!!!!
      m_geometry.mappedGridCurl3D(m_curlE0_cf,m_E0_ce);
#else
      // for 2D, only virtual E0 enters into curl for in-plane B
      computeIdealEatNodes(m_E0_virtual_nc,m_velocity,soln_magField_cf,m_Cspeed_cc);
      int this_fbc = m_magneticField_fbc;
      m_fluid_face_var_bc.at(this_fbc)->setNodeBC( a_soln_fluid, m_E0_virtual_nc, m_NodeBC_zeros, a_time );
      m_E0_virtual_nc.exchange();
      m_geometry.mappedGridCurl2D(m_curlE0_cf,m_E0_virtual_nc);
#endif
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
   m_geometry.divideCellVolume( m_dummyDiv );
   
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

   SpaceUtils::upWindToFaces(m_velocity_cf, m_velocity, m_velocity_cf, "c2");
   if(m_twoDaxisymm) {
      computeViscSourceAxisymm( m_momVisc_source, m_velocity, m_velocity_cf );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_momVisc_source[dit].mult( m_etaVisc[dit] ); 
         m_momVisc_source[dit].negate(); 
      }
   } 
   
   //   update RHSs for momentum equation
   if(a_accumulateMom || m_energy_is_thermal || m_energy_is_thermal_2) {

      // convert physical Ja*Flux to contravar and mult by faceArea
      // (Result is same as applying N^T*Flux*faceArea as is done by 
      // computeMetricTermProductAverage), which is not done here 
      // since Flux is already multiplied by Jacobian

      if(SpaceDim<3) {
         m_geometry.applyAxisymmetricCorrection( m_m0JaFluxVisc_cf );
         m_geometry.applyAxisymmetricCorrection( m_m1JaFluxVisc_cf );
      }
      m_geometry.computedxidXProductNorm(m_m0JaFluxVisc_norm, m_m0JaFluxVisc_cf);
      m_geometry.computedxidXProductNorm(m_m1JaFluxVisc_norm, m_m1JaFluxVisc_cf);
      if(SpaceDim==3) m_geometry.computedxidXProductNorm(m_m2JaFluxVisc_norm, m_m2JaFluxVisc_cf);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_momJaFluxVisc_norm[dit].copy( m_m0JaFluxVisc_norm[dit],0,0,1 ); 
         m_momJaFluxVisc_norm[dit].copy( m_m1JaFluxVisc_norm[dit],0,1,1 );
         if(SpaceDim==3) m_momJaFluxVisc_norm[dit].copy( m_m2JaFluxVisc_norm[dit],0,2,1 );
      }
   
      LevelData<FArrayBox>& rhs_mom( a_rhs_fluid.cell_var("momentumDensity") );
      m_geometry.mappedGridDivergenceFromFluxNorms(m_momJaFluxVisc_norm, m_dummyDiv_mom);
      for (DataIterator dit(rhs_mom.dataIterator()); dit.ok(); ++dit) {
         for (int n=0; n<SpaceDim; ++n) {
            const MagBlockCoordSys& coords( m_geometry.getBlockCoordSys( grids[dit] ) );
            m_dummyDiv_mom[dit].mult( 1.0/coords.getMappedCellVolume(),n,1 );
            if(n==0 && m_twoDaxisymm) {
               FArrayBox JaPiyyOverR(grids[dit], 1);
               JaPiyyOverR.copy(m_momVisc_source[dit]);
               JaPiyyOverR.mult(m_Jacobian[dit]); // Jacobian = Ja = 2*pi*R_phys * J_RZ
               JaPiyyOverR.divide(m_Xphys[dit],0,0,1);
               m_dummyDiv_mom[dit].minus(JaPiyyOverR,0,0,1);
            }
         }
         if(a_accumulateMom) rhs_mom[dit].minus(m_dummyDiv_mom[dit]);
      }
   
   }   

   //   update RHSs for energy density equation
   if(a_accumulateEne) {
 
      computeViscosityEnergyFlux( m_enJaFluxVisc_cf,
                                  m_velocity_cf );

      // convert physical Ja*Flux to contravar and mult by faceArea
      // (Result is same as applying N^T*Flux*faceArea as is done by 
      // computeMetricTermProductAverage), which is not done here 
      // since Flux is already multiplied by Jacobian
      
      m_geometry.applyAxisymmetricCorrection( m_enJaFluxVisc_cf );
      m_geometry.computedxidXProductNorm(m_enJaFluxVisc_norm, m_enJaFluxVisc_cf);

      LevelData<FArrayBox>& rhs_ene( a_rhs_fluid.cell_var("energyDensity") );
      m_geometry.mappedGridDivergenceFromFluxNorms(m_enJaFluxVisc_norm, m_dummyDiv);
      for (DataIterator dit(rhs_ene.dataIterator()); dit.ok(); ++dit) {
         const MagBlockCoordSys& coords( m_geometry.getBlockCoordSys( grids[dit] ) );
         m_dummyDiv[dit].mult( 1.0/coords.getMappedCellVolume() );
         rhs_ene[dit].minus( m_dummyDiv[dit],0,0,1 ); // specify comp since can have two (CGL?)
      }
      
      if(m_energy_is_thermal || m_energy_is_thermal_2) {

         LevelData<FArrayBox>& UdotDivMomFlux = m_dummyFArray_oneComp;
         for (DataIterator dit(grids); dit.ok(); ++dit) {

            const FArrayBox& this_velocity  = m_velocity[dit];
            const FArrayBox& this_momDiv = m_dummyDiv_mom[dit];
            FArrayBox& this_UdotMomDiv = UdotDivMomFlux[dit];

            const Box& thisbox = this_UdotMomDiv.box();
            FORT_EVAL_DOT_PRODUCT( CHF_BOX(thisbox),
                                   CHF_CONST_FRA(this_velocity),
                                   CHF_CONST_FRA(this_momDiv),
                                   CHF_FRA1(this_UdotMomDiv,0) );

            rhs_ene[dit].plus(this_UdotMomDiv,0,0,1); 

         }

      }

   }

}


void IdealMhdOp::updateRHSs( FluidSpecies&  a_rhs_fluid,
                       const FluidSpecies&  a_soln_fluid ) 
{
   CH_TIME("IdealMhdOp::updateRHSs()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // compute divergence of Fluxes and add (subtract) to RHS 
   // Note that Ja/cellVol = 1/mapVol

   // update RHS for mass density
   if(!m_blankN){
   LevelData<FArrayBox>& rhs_rho( a_rhs_fluid.cell_var(0) );
   m_geometry.mappedGridDivergenceFromFluxNorms(m_rhoFlux_norm, m_dummyDiv);
   for (DataIterator dit(rhs_rho.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& coords( m_geometry.getBlockCoordSys( grids[dit] ) );
      m_dummyDiv[dit].mult( 1.0/coords.getMappedCellVolume() );
      rhs_rho[dit].minus(m_dummyDiv[dit]);
      m_rhs_rho[dit].copy(rhs_rho[dit]);
   }
   }
  
   // update RHS for momentum density
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
         const MagBlockCoordSys& coords( m_geometry.getBlockCoordSys( grids[dit] ) );
         m_dummyDiv_mom[dit].mult( 1.0/coords.getMappedCellVolume(),n,1 );
      }
      //m_dummyDiv_mom[dit].setVal(0.0,0); // JRA, no mr update
      rhs_mom[dit].minus(m_dummyDiv_mom[dit]);
      if(m_twoDaxisymm && a_soln_fluid.m_evolve_momentumDensity_virtual) {
         FArrayBox JaRhoVthsqOverR(grids[dit], 1);
         JaRhoVthsqOverR.copy(m_momDen_virtual_cc[dit]);
         JaRhoVthsqOverR.mult(m_velocity_virtual[dit]);
         JaRhoVthsqOverR.mult(m_Jacobian[dit]); // Jacobian = 2*pi*R_phys * J_RZ
         JaRhoVthsqOverR.divide(m_Xphys[dit],0,0,1);
         rhs_mom[dit].plus(JaRhoVthsqOverR,0,0,1);            
      }
      if(a_soln_fluid.m_evolve_magneticField_virtual ||
         a_soln_fluid.m_evolve_magneticField) {
         rhs_mom[dit].plus(m_JaJcrossB[dit]);
      }
      m_forceDensity[dit].copy(rhs_mom[dit]);
   }
   
   // update RHS for energy density
   LevelData<FArrayBox>& rhs_ene( a_rhs_fluid.cell_var("energyDensity") );
   m_geometry.mappedGridDivergenceFromFluxNorms(m_enFlux_norm, m_dummyDiv);
   for (DataIterator dit(rhs_ene.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& coords( m_geometry.getBlockCoordSys( grids[dit] ) );
      m_dummyDiv[dit].mult( 1.0/coords.getMappedCellVolume() );
      rhs_ene[dit].minus( m_dummyDiv[dit],0,0,1 );
      if(m_energy_is_thermal) {
         rhs_ene[dit].plus( m_JaUdotGradP[dit],0,0,1 );
      }
      else if(m_energy_is_thermal_2) {
         rhs_ene[dit].minus( m_JaPdivU[dit],0,0,1 );
      }
      else if(a_soln_fluid.m_evolve_magneticField_virtual ||
              a_soln_fluid.m_evolve_magneticField) {
         rhs_ene[dit].plus( m_JaJdotE[dit],0,0,1 );
      }
   }
   
   // update RHS for parallel energy density if using CGL model
   if(m_CGL) {
      m_geometry.mappedGridDivergenceFromFluxNorms(m_enParFlux_norm, m_dummyDiv);
      for (DataIterator dit(rhs_ene.dataIterator()); dit.ok(); ++dit) {
         const MagBlockCoordSys& coords( m_geometry.getBlockCoordSys( grids[dit] ) );
         m_dummyDiv[dit].mult( 1.0/coords.getMappedCellVolume() );
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
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_oneComp[dit].setVal(1.0/0.0);
   }
 
   // update RHS for momentum density in virtual direction
   if(a_soln_fluid.m_evolve_momentumDensity_virtual){
      LevelData<FArrayBox>& rhs_momV( a_rhs_fluid.cell_var("momentumDensity_virtual") );
      m_geometry.mappedGridDivergenceFromFluxNorms(m_mvFlux_norm, m_dummyDiv);
      for (DataIterator dit(rhs_momV.dataIterator()); dit.ok(); ++dit) {
         FArrayBox& div_mvFlux(m_dummyDiv[dit]);
         const MagBlockCoordSys& coords( m_geometry.getBlockCoordSys( grids[dit] ) );
         div_mvFlux.mult( 1.0/coords.getMappedCellVolume() );
         rhs_momV[dit].minus( div_mvFlux );
         if(m_twoDaxisymm) {
            FArrayBox JaRhoVthVrOverR(grids[dit], 1);
            JaRhoVthVrOverR.copy(m_momDen_virtual_cc[dit]);
            JaRhoVthVrOverR.mult(m_velocity[dit],0,0,1);
            JaRhoVthVrOverR.mult(m_Jacobian[dit]); // Jacobian = 2*pi*R_phys * J_RZ
            JaRhoVthVrOverR.divide(m_Xphys[dit],0,0,1);
            rhs_momV[dit].minus(JaRhoVthVrOverR,0,0,1);            
         }
         rhs_momV[dit].plus(m_JaJcrossB_virtual[dit]);
      }
   }

   // update RHS for virtual magnetic field
   if(a_soln_fluid.m_evolve_magneticField_virtual){
      LevelData<FArrayBox>& rhs_magV( a_rhs_fluid.cell_var("magneticField_virtual") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs_magV[dit].minus( m_curlE0_virtual[dit] );
      }
   }
   
   // update RHS for magnetic field
   if(a_soln_fluid.m_evolve_magneticField && !m_blankB){
      LevelData<FluxBox>& rhs_magField( a_rhs_fluid.face_var("magneticField") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; dir++) {
            rhs_magField[dit][dir].minus( m_curlE0_cf[dit][dir] );
         }
      }
   }

}   

void IdealMhdOp::initialize( CFGVars&    a_species, 
                             const Real  a_time )
{
   CH_TIME("IdealMhdOp::initialize()");
  
   // initilize cell center variables
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

bool IdealMhdOp::isInitializationConstrained(const FluidSpecies& a_fluid_phys)
{
   return true;
}

void IdealMhdOp::applyInitializationConstraints(FluidSpeciesPtrVect&               a_fluid_comp,
                                                FluidSpeciesPtrVect&               a_fluid_phys,
                                                const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                                const EMFields&                    a_EM_fields,
                                                const int                          a_component,
                                                const double                       a_time )
{
   CH_TIME("IdealMhdOp::applyInitializationConstraints()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
     
   FluidSpecies& species_phys( static_cast<FluidSpecies&>(*(a_fluid_phys[a_component])) );
   FluidSpecies& species_comp( static_cast<FluidSpecies&>(*(a_fluid_comp[a_component])) );
 
   if(m_initializeBackgroundPlusPert) { 

      // set background profiles for fluid
      const GridFunction& fsic_rho0( fluidSpeciesIC( "RHO0" ) );
      fsic_rho0.assign( m_N0, species_phys.configurationSpaceGeometry(), a_time );
      m_N0.exchange();
     
      const GridFunction& fsic_T0( fluidSpeciesIC( "T0" ) );
      fsic_T0.assign( m_T0, species_phys.configurationSpaceGeometry(), a_time );
      m_T0.exchange();
      
      for (auto dit(m_P0.dataIterator()); dit.ok(); ++dit) {
         m_P0[dit].copy(m_N0[dit]);
         m_P0[dit].mult(m_T0[dit]);
         m_P0[dit].mult(2.0);
      }
      
      // get the perturbation for density
      const GridFunction& fsic_deltaN( fluidSpeciesIC( "deltaRHO" ) );
      LevelData<FArrayBox>& deltaN = m_dummyFArray_oneComp;
      fsic_deltaN.assign( deltaN, species_phys.configurationSpaceGeometry(), a_time );
      deltaN.exchange();

      // now set the phys and comp vars
      LevelData<FArrayBox>& rhoDen( species_phys.cell_var(0) );
      LevelData<FArrayBox>& eneDen( species_phys.cell_var("energyDensity") );
      LevelData<FArrayBox>& rhoDen_comp( species_comp.cell_var(0) );
      LevelData<FArrayBox>& eneDen_comp( species_comp.cell_var("energyDensity") );
      const double gamma = species_phys.m_gamma;
      for (auto dit(eneDen.dataIterator()); dit.ok(); ++dit) {
         rhoDen[dit].copy(m_N0[dit]);
         rhoDen[dit].mult(deltaN[dit]);
         eneDen[dit].copy(m_P0[dit]);
         eneDen[dit].mult(deltaN[dit]);
         eneDen[dit].divide(gamma-1.0);
         rhoDen_comp[dit].copy(rhoDen[dit]);
         eneDen_comp[dit].copy(eneDen[dit]);
      }
      m_geometry.multJonValid(rhoDen_comp);
      m_geometry.multJonValid(eneDen_comp);

      //
      //   initialize magnetic field
      //

      // get By0 from input file
      const GridFunction& fsic_By0( fluidSpeciesIC( "BY0" ) );
      LevelData<FArrayBox>& By0 = m_dummyFArray_oneComp;
      fsic_By0.assign( By0, species_phys.configurationSpaceGeometry(), a_time );
      By0.exchange();
      
      // copy By0 to m_B0_phys tor_dir component
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const int tor_dir = 1;
         m_B0_phys[dit].setVal(0.0);
         m_B0_phys[dit].copy(By0[dit],0,tor_dir,1);
      }
   
      // convert to contravar and put it on cell faces
      LevelData<FArrayBox>& B0_contra = m_dummyFArray_spaceDim;
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         B0_contra[dit].copy(m_B0_phys[dit]);
      }
      m_geometry.convertPhysToContravar(B0_contra,0);
      LevelData<FluxBox>& B_withBCs(species_phys.face_var("magneticField"));
      SpaceUtils::upWindToFaces(B_withBCs, B0_contra, B_withBCs, "c2");
      
      // apply BCs   
      int this_fbc = m_magneticField_fbc;
      SpaceUtils::exchangeFluxBox(B_withBCs);
      m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( species_phys, B_withBCs, a_time );
      
      // compute Ja*J00_phys
      LevelData<FArrayBox>& B0_covar = m_dummyFArray_spaceDim;
      SpaceUtils::interpFaceVectorToCell(B0_covar,B_withBCs,"c2");
      m_geometry.convertContravarToCovar(B0_covar,0);
      //for (DataIterator dit(grids); dit.ok(); ++dit) {
      //   B0_covar[dit].copy(m_B0_phys[dit]);
      //}
      //m_geometry.convertPhysToCovar(B0_covar,0);
      m_geometry.mappedGridCurl( m_JaJ00, B0_covar ); // Ja*curl(B0)_contra
      m_geometry.convertPhysToContravar(m_JaJ00,1);
      
      // set the comp B = Jacobian * B^q
      LevelData<FluxBox>& B_noBCs( species_comp.face_var("magneticField") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            B_noBCs[dit][dir].copy(B_withBCs[dit][dir]);
         }       
      }
      m_geometry.multJonFaces(B_noBCs);         
   
    
   }

#if CH_SPACEDIM==3 
   if(m_initializeProfilesFromPSI && species_phys.m_evolve_magneticField) {
 
      double rhoDen_floor = 0.0;   
      double eneDen_floor = 0.0;   
   
      // do this before adding perturbation 
      species_phys.pressure(m_P0);
      m_P0.exchange();
      
      LevelData<FArrayBox>& rhoDen( species_phys.cell_var(0) );
      LevelData<FArrayBox>& eneDen( species_phys.cell_var("energyDensity") );
      for (auto dit(rhoDen.dataIterator()); dit.ok(); ++dit) {
         m_N0[dit].copy(rhoDen[dit]);
         m_T0[dit].copy(m_P0[dit]);
         m_T0[dit].divide(m_N0[dit]);
         m_T0[dit].divide(2.0);
      }
   
      if(m_relative_density_floor > 0.0) {
         double local_Nmax, local_Emax;
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            Box box( grids[dit] );
            double box_Nmax( rhoDen[dit].max( box ) );
            double box_Emax( eneDen[dit].max( box ) );
            local_Nmax = Max( local_Nmax, box_Nmax );
            local_Emax = Max( local_Emax, box_Emax );
         }
   
         double Nmax = local_Nmax;
         double Emax = local_Emax;
#ifdef CH_MPI
         MPI_Allreduce( &local_Nmax, &Nmax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
         MPI_Allreduce( &local_Emax, &Emax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif
         rhoDen_floor = m_relative_density_floor*Nmax;
         eneDen_floor = m_relative_density_floor*Emax;
      
         // adjust background density and pressure to include density floor
         for (auto dit(rhoDen.dataIterator()); dit.ok(); ++dit) {
            m_N0[dit].plus(rhoDen_floor);
            m_P0[dit].copy(m_N0[dit]);
            m_P0[dit].mult(m_T0[dit]);
            m_P0[dit].mult(2.0);
         }
      }
 
      LevelData<FArrayBox>& deltaE = m_dummyFArray_oneComp;
      if(m_addEnergyPerturbation) { 
         const GridFunction& fsic_deltaE( fluidSpeciesIC( "deltaE" ) );
         fsic_deltaE.assign( deltaE, species_phys.configurationSpaceGeometry(), a_time );
         deltaE.exchange();
      }
      else {
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            deltaE[dit].setVal(1.0);
         }
      }

      LevelData<FArrayBox>& eneDen_comp( species_comp.cell_var("energyDensity") );
      for (auto dit(eneDen.dataIterator()); dit.ok(); ++dit) {
         eneDen[dit].mult(deltaE[dit]);
         eneDen[dit].plus(eneDen_floor); // do this after mult by perturbation
         eneDen_comp[dit].copy(eneDen[dit]);
      }
      m_geometry.multJonValid(eneDen_comp);

      LevelData<FArrayBox>& rhoDen_comp( species_comp.cell_var(0) );
      for (auto dit(eneDen.dataIterator()); dit.ok(); ++dit) {
         if(m_perturbDensity) rhoDen[dit].mult(deltaE[dit]);
         rhoDen[dit].plus(rhoDen_floor); // do this after mult by perturbation
         rhoDen_comp[dit].copy(rhoDen[dit]);
      }
      m_geometry.multJonValid(rhoDen_comp);

      //
      //
      //

      const GridFunction& fsic_RBTOR( fluidSpeciesIC( "RBTOR" ) );
      LevelData<FArrayBox>& RBTOR = m_dummyFArray_oneComp;
      fsic_RBTOR.assign( RBTOR, species_phys.configurationSpaceGeometry(), a_time );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         RBTOR[dit].mult(m_rbtor_fact);
      }
      LevelData<FluxBox>& RBTOR_fc = m_dummyFlux_oneComp;
      SpaceUtils::upWindToFaces(RBTOR_fc, RBTOR, RBTOR_fc, "c2");

      const GridFunction& fsic_PSI( fluidSpeciesIC( "PSI" ) );
      LevelData<FArrayBox>& PSI = m_PSI;
      fsic_PSI.assign( PSI, species_phys.configurationSpaceGeometry(), a_time );
      PSI.exchange();

      // Need dxi0/dZ and dxi0/dR as in axisymmetric to compute Bcyl
      m_geometry.setPointwiseNJInverseOnFaces(grids, m_Xphys.ghostVect());
      const LevelData<FluxBox>& Xphys = m_geometry.getFaceCenteredRealCoords();
      const LevelData<FluxBox>& dxidX = m_geometry.getFaceCentereddxidX();
      LevelData<FluxBox> dxi0dR(grids,1,dxidX.ghostVect());
      LevelData<FluxBox> dxi0dZ(grids,1,dxidX.ghostVect());
      LevelData<FluxBox> Rcyl(grids,1,dxidX.ghostVect());
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         Real x, y, R, costh, sinth, dxi0dX, dxi0dY;
         for (int dir = 0; dir<SpaceDim; dir++) {
            const FArrayBox& Xcar = Xphys[dit][dir];
            const Box& box = dxi0dR[dit][dir].box();
            for (BoxIterator bit(box); bit.ok(); ++bit) {
               IntVect iv = bit();

               x = Xcar(iv,0);
               y = Xcar(iv,1);
               R = sqrt(x*x + y*y);
               costh = x/R;
               sinth = y/R;
 
               dxi0dX = dxidX[dit][dir](iv,0);
               dxi0dY = dxidX[dit][dir](iv,3);
               dxi0dR[dit][dir](iv,0) = dxi0dX*costh + dxi0dY*sinth;
               dxi0dZ[dit][dir](iv,0) = dxidX[dit][dir](iv,6);
               Rcyl[dit][dir](iv,0) = R;

            }
         }
      }
       
      // compute Bfield in cyl coords: 
      // B_p = dPSi/dxi0/R
      // B_R = -B_p*dxi0/dZ, B_Z = Bp*dxi0/dR
      LevelData<FluxBox>& Bfield = m_dummyFlux_spaceDim;
      m_geometry.computeMappedGradient( PSI, Bfield, 2);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir = 0; dir<SpaceDim; dir++) {
            Bfield[dit][dir].mult(m_bpol_fact);
         }
      }
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir = 0; dir<SpaceDim; dir++) {
            Bfield[dit][dir].copy(Bfield[dit][dir],0,2,1);
            Bfield[dit][dir].mult(dxi0dR[dit][dir],0,2,1);
            Bfield[dit][dir].divide(Rcyl[dit][dir],0,2,1);
            //
            Bfield[dit][dir].copy(RBTOR_fc[dit][dir],0,1,1);
            Bfield[dit][dir].divide(Rcyl[dit][dir],0,1,1);
            //
            Bfield[dit][dir].mult(dxi0dZ[dit][dir],0,0,1);
            Bfield[dit][dir].divide(Rcyl[dit][dir],0,0,1);
            Bfield[dit][dir].negate(0);
         }
      }
       
      // convert Bfield from cyl to car
      for (DataIterator dit(grids); dit.ok(); ++dit) {

         Real x, y, R, costh, sinth, B_R, B_T;
         for (int dir = 0; dir<SpaceDim; dir++) {
            const FArrayBox& Xcar = Xphys[dit][dir];
            const Box& box = Bfield[dit][dir].box();
            for (BoxIterator bit(box); bit.ok(); ++bit) {
               IntVect iv = bit();
               x = Xcar(iv,0);
               y = Xcar(iv,1);
               R = sqrt(x*x + y*y); 
               costh = x/R;
               sinth = y/R;
               B_R = Bfield[dit][dir](iv,0);
               B_T = Bfield[dit][dir](iv,1);
               Bfield[dit][dir](iv,0) = costh * B_R - sinth * B_T;
               Bfield[dit][dir](iv,1) = sinth * B_R + costh * B_T;
            }
         }

      }
       
      // convert Bield from car to contravariant
      m_geometry.convertPhysToContravar(Bfield,0);
       
      // copy normal component on each dir to phys variables
      LevelData<FluxBox>& B_withBCs(species_phys.face_var("magneticField"));
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         B_withBCs[dit][0].copy(Bfield[dit][0],0,0,1);
         B_withBCs[dit][1].copy(Bfield[dit][1],1,0,1);
         B_withBCs[dit][2].copy(Bfield[dit][2],2,0,1);
      }
      
      // apply BCs   
      int this_fbc = m_magneticField_fbc;
      SpaceUtils::exchangeFluxBox(B_withBCs);
      m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( species_phys, B_withBCs, a_time );
   

      // compute Ja*J00_phys and B0_phys
      LevelData<FArrayBox>& B0_covar = m_dummyFArray_spaceDim;
      SpaceUtils::interpFaceVectorToCell(B0_covar,B_withBCs,"c2");
      m_geometry.convertPhysToContravar(B0_covar,1);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_B0_phys[dit].copy(B0_covar[dit]);
      }
      m_geometry.convertPhysToCovar(B0_covar,0);          
      m_geometry.mappedGridCurl( m_JaJ00, B0_covar ); // Ja*curl(B0)_contra
    
      //LevelData<FluxBox>& Bvec_covar = m_dummyFlux_spaceDim;
      //SpaceUtils::interpFacesToFaces(Bvec_covar, B_withBCs);
      //m_geometry.convertContravarToCovar(Bvec_covar,0);
      //m_geometry.mappedGridCurl( m_JaJ00, Bvec_covar ); // Ja*curl(B)_contra
   
      m_geometry.convertPhysToContravar(m_JaJ00,1);
      
      if(m_perturbBtor) {

         const GridFunction& fsic_deltaB( fluidSpeciesIC( "deltaB" ) );
         LevelData<FArrayBox>& deltaB = m_dummyFArray_oneComp;
         fsic_deltaB.assign( deltaB, species_phys.configurationSpaceGeometry(), a_time );
         deltaB.exchange();
         
         LevelData<FluxBox>& deltaB_fc = m_dummyFlux_oneComp;
         SpaceUtils::upWindToFaces(deltaB_fc, deltaB, deltaB_fc, "c2");
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            const int tor_dir = 1;
            B_withBCs[dit][tor_dir].mult(deltaB_fc[dit][tor_dir]);
         }

      }
 
      // set the comp B = Jacobian * B^q
      LevelData<FluxBox>& B_noBCs( species_comp.face_var("magneticField") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            B_noBCs[dit][dir].copy(B_withBCs[dit][dir]);
         }       
      }
      m_geometry.multJonFaces(B_noBCs);         
   

   }
#else 
   if(m_initializeProfilesFromPSI && species_phys.m_evolve_magneticField_virtual) {

      // do this before adding perturbation 
      species_phys.pressure(m_P0);
      m_P0.exchange();
      
      if(m_addEnergyPerturbation) { 

         const GridFunction& fsic_deltaE( fluidSpeciesIC( "deltaE" ) );
         LevelData<FArrayBox>& deltaE = m_dummyFArray_oneComp;
         fsic_deltaE.assign( deltaE, species_phys.configurationSpaceGeometry(), a_time );
         deltaE.exchange();

         LevelData<FArrayBox>& eneDen( species_phys.cell_var("energyDensity") );
         LevelData<FArrayBox>& eneDen_comp( species_comp.cell_var("energyDensity") );
         for (auto dit(eneDen.dataIterator()); dit.ok(); ++dit) {
            eneDen[dit].mult(deltaE[dit]);
            eneDen_comp[dit].copy(eneDen[dit]);
         }
         m_geometry.multJonValid(eneDen_comp);
         
         if(m_perturbDensity) {
            LevelData<FArrayBox>& den( species_phys.cell_var(0) );
            LevelData<FArrayBox>& den_comp( species_comp.cell_var(0) );
            for (auto dit(eneDen.dataIterator()); dit.ok(); ++dit) {
               den[dit].mult(deltaE[dit]);
               den_comp[dit].copy(den[dit]);
            }
            m_geometry.multJonValid(den_comp);
         }

      }

      // get the toroidal magnetic field profile from R*BT
      const GridFunction& fsic_RBTOR( fluidSpeciesIC( "RBTOR" ) );
      LevelData<FArrayBox>& BT_withBCs(species_phys.cell_var("magneticField_virtual"));
      fsic_RBTOR.assign( BT_withBCs, species_phys.configurationSpaceGeometry(), a_time );
      LevelData<FArrayBox>& BT_noBCs( species_comp.cell_var("magneticField_virtual") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         BT_withBCs[dit].divide(m_Xphys[dit],0,0,1);
         BT_noBCs[dit].copy(BT_withBCs[dit]);
      }
      m_geometry.multJonValid(BT_noBCs);         
       
      // compute the contravariant poloidal magnetic field profile from PSI: 
      // B^1 = d(PSI)/dxi0 / J_2D*R = 2*pi*d(PSI)/dxi0/Ja

      const GridFunction& fsic_PSI( fluidSpeciesIC( "PSI" ) );
      LevelData<FArrayBox>& PSI = m_PSI;
      fsic_PSI.assign( PSI, species_phys.configurationSpaceGeometry(), a_time );
      PSI.exchange();

      LevelData<FluxBox>& Bp_tmp = m_dummyFlux_spaceDim;
      m_geometry.computeMappedGradient( PSI, Bp_tmp, 2);

      LevelData<FluxBox>& B_withBCs(species_phys.face_var("magneticField"));
      LevelData<FluxBox>& B_noBCs( species_comp.face_var("magneticField") );
      double twoPi = 2.0*Pi;
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            if(dir==1) {
               B_withBCs[dit][dir].copy(Bp_tmp[dit][dir],0,0,1);
               B_withBCs[dit][dir].mult(twoPi);
            } 
            else {
               B_withBCs[dit][dir].setVal(0.0);
            }
            B_noBCs[dit][dir].copy(B_withBCs[dit][dir]);
         }
      }
      m_geometry.divideJonFaces(B_withBCs);         
       
   }
#endif 
 
   if(species_comp.m_evolve_magneticField && m_initializeBfromVectorPotential) {
         
      LevelData<FluxBox>& contravar_B_withBCs( species_phys.face_var("magneticField") );
      int this_fbc = m_magneticField_fbc;
     
      // compute Ja*curlB_phys at cell center using input B field
      SpaceUtils::interpFaceVectorToCell(m_B_covar_cc,contravar_B_withBCs,"c2");
      m_geometry.convertPhysToContravar(m_B_covar_cc,1);
      m_geometry.convertPhysToCovar(m_B_covar_cc,0);
      m_geometry.mappedGridCurl( m_JaJ0_cc, m_B_covar_cc ); // Ja*curl(B)_contra
      m_geometry.convertPhysToContravar(m_JaJ0_cc,1);

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
         m_fluid_face_var_bc.at(this_fbc)->setEdgeBC( species_phys, Az_ce, m_EdgeBC_zeros, a_time );
         SpaceUtils::exchangeEdgeDataBox(Az_ce);

         // compute B = curlA (divergence free and satisfies BCs)
         m_geometry.mappedGridCurl3D(contravar_B_withBCs,Az_ce);
         m_geometry.divideJonFaces(contravar_B_withBCs);
         SpaceUtils::exchangeFluxBox(contravar_B_withBCs); // must exchange before applyBC !!!!
         m_fluid_face_var_bc.at(this_fbc)->applyFluxBC( species_phys, contravar_B_withBCs, a_time );
         LevelData<FluxBox>& mapped_contravar_B_noBCs( species_comp.face_var("magneticField") );
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            mapped_contravar_B_noBCs[dit].copy(contravar_B_withBCs[dit]);
         }
         m_geometry.multJonFaces(mapped_contravar_B_noBCs);         
   
      }

      // compute Ja time contravariant curlE0 on faces
      setCellCenterValues( species_phys );
      const double gamma = species_phys.m_gamma;
      setMappedCspeed( gamma, m_thermPforE0, m_magPforE0 ); // JRA
      
      computeIdealEatEdges(m_E0_ce,m_velocity,contravar_B_withBCs,m_Cspeed_cc); // Ja*E0_covar
      m_fluid_face_var_bc.at(this_fbc)->setEdgeBC( species_phys, m_E0_ce, m_EdgeBC_zeros, a_time );
      SpaceUtils::exchangeEdgeDataBox(m_E0_ce); // call exchange after setBC !!!!! 
      m_geometry.mappedGridCurl3D(m_curlE0_cf,m_E0_ce);

      // compute physical E0 at cell center
      SpaceUtils::interpEdgesToCell(m_E0_cc,m_E0_ce,"c2"); // covar at cell center
      m_geometry.convertPhysToCovar(m_E0_cc,1);            // phys at cell center
      
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
   m_geometry.divideCellVolume( a_divF_phys );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      F_norm[dit].setVal(1./0.);
   }

}     

void IdealMhdOp::computeCrossProduct( LevelData<FArrayBox>&  a_A,
                                const LevelData<FArrayBox>&  a_B,
                                const LevelData<FArrayBox>&  a_C, 
                                const int                    a_virtual_case ) const
{
   CH_TIME("IdealMhdOp::computeCrossProduct()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
 
   // compute A = BxC
   
   if(a_virtual_case==1) { // A is virtual: Ay = Bz*Cx - Bx*Cz
      CH_assert(a_A.nComp()==1);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
               FArrayBox& A_on_patch = a_A[dit];
         const FArrayBox& B_on_patch = a_B[dit];
         const FArrayBox& C_on_patch = a_C[dit];

         const Box& cellbox = A_on_patch.box();
         Real Ay, Bx, Bz, Cx, Cz;
         for (BoxIterator bit(cellbox); bit.ok(); ++bit) {
             IntVect ic = bit();
             Bx = B_on_patch.get(ic,0);
             Bz = B_on_patch.get(ic,1);
             Cx = C_on_patch.get(ic,0);
             Cz = C_on_patch.get(ic,1);
             Ay = Bz*Cx - Bx*Cz;
             A_on_patch.set(ic,0,Ay);
         }
      }
   }
   else if(a_virtual_case==2) { // B is virtual
      CH_assert(a_B.nComp()==1);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& thisbox = a_A[dit].box();
         for (int n=0; n<SpaceDim; ++n) a_A[dit].copy( a_B[dit],0,n,1 ); 
         a_A[dit].mult( a_C[dit],thisbox,1,0,1 ); // Ax =  By*Cz
         a_A[dit].mult( a_C[dit],thisbox,0,1,1 ); // Az = -By*Cx
         a_A[dit].mult( -1.0,thisbox,1,1 );
      }
   }
   else if(a_virtual_case==3) { // C is virtual   
      CH_assert(a_C.nComp()==1);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& thisbox = a_A[dit].box();
         for (int n=0; n<SpaceDim; ++n) a_A[dit].copy( a_C[dit],0,n,1 ); 
         a_A[dit].mult( a_B[dit],thisbox,1,0,1 ); // Ax = -Bz*Cy
         a_A[dit].mult( a_B[dit],thisbox,0,1,1 ); // Az =  Bx*Cy
         a_A[dit].mult( -1.0,thisbox,0,1 );
      }
   }
   else { // nothing is virtual
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         
               FArrayBox& A_on_patch = a_A[dit];
         const FArrayBox& B_on_patch = a_B[dit];
         const FArrayBox& C_on_patch = a_C[dit];

         const Box& cellbox = A_on_patch.box();
         FORT_EVAL_CROSS_PRODUCT( CHF_BOX(cellbox),
                                  CHF_CONST_FRA(B_on_patch),
                                  CHF_CONST_FRA(C_on_patch),
                                  CHF_FRA(A_on_patch) );
      }
   }
   
}     

void IdealMhdOp::computeCrossProduct( LevelData<FArrayBox>&  a_A_inPlane,
                                const LevelData<FArrayBox>&  a_B_inPlane, 
                                const LevelData<FArrayBox>&  a_C_inPlane,
                                const LevelData<FArrayBox>&  a_B_virtual, 
                                const LevelData<FArrayBox>&  a_C_virtual ) const
{
   CH_TIME("IdealMhdOp::computeCrossProduct()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   CH_assert(SpaceDim==2);
 
   // general 2D cross product: A = BxC
   // Ax = By*Cz - Bz*Cy
   // Az = Bx*Cy - By*Cx

   for (DataIterator dit(grids); dit.ok(); ++dit) {

            FArrayBox& A_on_patch = a_A_inPlane[dit];
      const FArrayBox& B_on_patch = a_B_inPlane[dit];
      const FArrayBox& C_on_patch = a_C_inPlane[dit];
      const FArrayBox& Bv_on_patch = a_B_virtual[dit];
      const FArrayBox& Cv_on_patch = a_C_virtual[dit];

      const Box& cellbox = A_on_patch.box();
      FORT_EVAL_CROSS_PRODUCT_2D( CHF_BOX(cellbox),
                                  CHF_CONST_FRA(B_on_patch),
                                  CHF_CONST_FRA(C_on_patch),
                                  CHF_CONST_FRA1(Bv_on_patch,0),
                                  CHF_CONST_FRA1(Cv_on_patch,0),
                                  CHF_FRA(A_on_patch) );
   }
   
}     

#include "NamespaceFooter.H"
