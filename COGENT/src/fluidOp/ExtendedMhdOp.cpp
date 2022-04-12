#include "ExtendedMhdOp.H"
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
#include "ExtendedMhdOpF_F.H"

#include "NamespaceHeader.H" 

const char* ExtendedMhdOp::pp_name = {"extendedmhd_op"};


ExtendedMhdOp::ExtendedMhdOp( const string&   a_pp_str,
                              const string&   a_species_name,
                              const MagGeom&  a_geometry,
                              const double    a_Nscale,
                              const double    a_Tscale,
                              const double    a_Xscale,
                              const double    a_Mscale,
                              const int       a_verbosity )
   : m_verbosity(a_verbosity),
     m_geometry(a_geometry),
     m_updateHydroVars(true),
     m_useRelaxedHeatFlux(false), 
     m_useRelaxedViscosity(false), 
     m_useRelaxedOhmsLaw(true), 
     m_include_divJStress(false),
     m_By_insulator_bc(false),
     m_tauiMin(1.0e-4),
     m_taueMin(1.0e-4),
     m_etaResMin(1.0e-4),
     m_etaAnomFactor(0.0),
     m_NvacC(0.0e-1),
     m_NvacP(4.0),
     m_Nmin(1.0e-4),
     m_Tmin(1.0e-1),
     m_Tmax(1.0e5),
     m_etaMin(1.0e-3),
     m_etaVis_ele(0.0),
     m_viscLimiter(1.0e20),
     m_gyroViscLimiter(1.0e20),
     m_gyroHeatViscCoeff(0.0),
     m_kappaMin(1.0e-3),
     m_heatLimiter(1.0e20),
     m_heatLimiter_ele(1.0e20),
     m_gyroHeatLimiter(1.0e20),
     m_xionFactor(1.0),
     m_xeleFactor(1.0),
     m_advScheme("c2"),
     m_meRel(1.0),
     m_epRel(1.0),
     m_useJ0forJcc(false),
     m_fullyImplicitJ(false),
     m_iterHallMax(1),
     m_HallResTol(1.0e-6),
     m_advSchemeHall("c2"),
     m_advSchemeJStress("c2"),
     m_modelDriftTerms(false),
     m_initializeJfromAmperesLaw(false),
     m_initializeEfromOhmsLaw(false),
     m_species_name(a_species_name),
     m_courant_time_step(DBL_MAX),
     m_opt_string(a_pp_str),
     m_my_pc_idx(-1)

{

   ParmParse pp(a_pp_str.c_str());
   parseParameters( pp, 
                    a_Nscale, a_Tscale,
                    a_Xscale, a_Mscale );
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
   parseInitialConditions(a_species_name, variable_name);
   m_density_cbc = m_fluid_bc.size();
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );

   // set IC and create BC for electron density
   //
   if (ppic.contains("electronDensity.function")) { // JRA, need way to not go here if not evolving
      variable_name = "electronDensity";
      parseInitialConditions(a_species_name, variable_name);
      m_electronDensity_cbc = m_fluid_bc.size();
      m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                            a_geometry.getCoordSys()->type(), true) );
   }

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
 
   // set IC and create BC for ion energy density
   //
   variable_name = "energyDensity";
   CH_assert( ppic.contains("energyDensity.function") );
   parseInitialConditions(a_species_name, variable_name);
   m_energyDensity_cbc = m_fluid_bc.size();
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );
   
   // set IC and create BC for electron energy density
   //
   variable_name = "energyDensity_ele";
   CH_assert( ppic.contains("energyDensity_ele.function") );
   parseInitialConditions(a_species_name, variable_name);
   m_energyDensity_ele_cbc = m_fluid_bc.size();
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );
   
   // set IC and create BC for out-of-plane momentum density
   //
   if (ppic.contains("momentumDensity_virtual.function")) {
      variable_name = "momentumDensity_virtual";
      parseInitialConditions(a_species_name, variable_name);
      m_momentumDensity_virtual_cbc = m_fluid_bc.size();
      m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                            a_geometry.getCoordSys()->type(), false) );
   }
   
   // set IC and create BC for out-of-plane magnetic field
   //
   if (ppic.contains("magneticField_virtual.function")) {
      variable_name = "magneticField_virtual";
      parseInitialConditions(a_species_name, variable_name);
      m_magneticField_virtual_cbc = m_fluid_bc.size();
      m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                            a_geometry.getCoordSys()->type(), false) );
   
      //Check for insulator BCs for virtual magnetic field
      m_By_insulator_bc = m_fluid_bc.back()->isInsulatorConductorBC();
      //m_By_insulator_bc = m_fluid_bc.at(m_magneticField_virtual_cbc)->isInsulatorConductorBC();
   }
 
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

   // set IC and create BC for electric field
   //
   if (ppic.contains("electricField_0.function")) {
      variable_name = "electricField";
      CH_assert( ppic.contains("electricField_0.function") );
      CH_assert( ppic.contains("electricField_1.function") );
      parseInitialConditions(a_species_name, "electricField_0");
      parseInitialConditions(a_species_name, "electricField_1");
      if (ppic.contains("electricField_2.function")) {
         parseInitialConditions(a_species_name, "electricField_2");
      }
      m_electricField_ebc = m_fluid_edge_var_bc.size();
      m_fluid_edge_var_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name,
                                     a_geometry.getCoordSys()->type(), false) );
   }

   // set IC and create BC for current density
   //
   if (ppic.contains("currentDensity_0.function")) {
      variable_name = "currentDensity";
      CH_assert( ppic.contains("currentDensity_0.function") );
      CH_assert( ppic.contains("currentDensity_1.function") );
      parseInitialConditions(a_species_name, "currentDensity_0");
      parseInitialConditions(a_species_name, "currentDensity_1");
      if (ppic.contains("currentDensity_2.function")) {
         parseInitialConditions(a_species_name, "currentDensity_2");
      }
      m_currentDensity_ebc = m_fluid_edge_var_bc.size();
      m_fluid_edge_var_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name,
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
   setMappedCLight(); // only need to do this once for static grids

}


ExtendedMhdOp::~ExtendedMhdOp()
{
   for (std::vector<FluidVarBC*>::iterator it = m_fluid_bc.begin(); it != m_fluid_bc.end(); ++it)
   {
      delete (*it);
   }
   m_fluid_bc.clear();
   
   for (std::vector<FluidVarBC*>::iterator it = m_fluid_face_var_bc.begin(); it != m_fluid_face_var_bc.end(); ++it)
   {
      delete (*it);
   }
   m_fluid_face_var_bc.clear();
   
   for (std::vector<FluidVarBC*>::iterator it = m_fluid_edge_var_bc.begin(); it != m_fluid_edge_var_bc.end(); ++it)
   {
      delete (*it);
   }
   m_fluid_edge_var_bc.clear();
   
}


void ExtendedMhdOp::accumulateMhdRHS( FluidSpecies&  a_rhs_fluid,
                                const FluidSpecies&  a_soln_fluid,
                                const Real           a_time )
{
   CH_TIME("ExtendedMhdOp::accumulateMhdRHS()");

   //CH_assert(soln_fluid.m_evolve_massDensity==1);
   CH_assert(a_soln_fluid.m_evolve_momentumDensity==1);
   CH_assert(a_soln_fluid.m_evolve_energyDensity==1);
   //CH_assert(a_soln_fluid.m_evolve_energyDensity_ele==1);
   if(m_modelDriftTerms) CH_assert(a_soln_fluid.m_evolve_currentDensity);
   
   setCellCenterValues( a_soln_fluid );

   setCellCenterFluxes( a_soln_fluid );
   
   // compute flux-freezing speed at cell-center for each direction
   // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/N); i = q0, q1, q2
   //
   const double gamma = a_soln_fluid.m_gamma;
   setMappedCspeed( gamma, 1, 1 );

   updateCollisionTerms( a_soln_fluid, m_pressure_ele, m_pressure_ion ); 

   setMagneticFieldTerms( a_soln_fluid );

   setFaceCenteredFluxes( a_soln_fluid );
   
   enforceFluxBCs( a_soln_fluid, a_time );

   updateRHSs( a_rhs_fluid, a_soln_fluid );
    
}

void ExtendedMhdOp::accumulateHeatRHS( FluidSpecies&  a_rhs_fluid,
                                 const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("ExtendedMhdOp::accumulateHeatRHS()");

   const DisjointBoxLayout& grids( m_geometry.grids() );
  
   if(!m_useRelaxedHeatFlux) {
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_qe_cf[dit].setVal( 0.0 );
         m_qi_cf[dit].setVal( 0.0 );
      }
 
   }  

   updateRHSs_heat( a_rhs_fluid, a_soln_fluid );

}

void ExtendedMhdOp::accumulateViscRHS( FluidSpecies&  a_rhs_fluid,
                                 const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("ExtendedMhdOp::accumulateViscRHS()");
  
   if(!m_useRelaxedViscosity) {
      
      updateCollisionTerms( a_soln_fluid, m_pressure_ele, m_pressure_ion ); 
      
      setIonViscosityCoefficients( m_taui, m_xion, m_pressure_ion );
 
      SpaceUtils::upWindToFaces(m_velocity_cf,m_velocity,m_velocity_cf,"c2");
      
      computeViscFluxes( m_m0JaFluxVisc_cf, m_m1JaFluxVisc_cf, 
                         m_momVisc_source,
                         m_velocity_cf, m_velocity,
                         m_etaVisC_cf,  m_etaVisC );
      
      computeViscEnergyFlux( m_enJaFluxVisc_cf,
                             m_m0JaFluxVisc_cf,
                             m_m1JaFluxVisc_cf,
                             m_velocity_cf );
 
   }  
   updateRHSs_visc( a_rhs_fluid, a_soln_fluid );

}

void ExtendedMhdOp::accumulateExplicitRHS( FluidSpeciesPtrVect&        a_rhs,
                                  const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                  const FluidSpeciesPtrVect&        a_fluid_species,
                                  const PS::ScalarPtrVect&          a_scalars,
                                  const EField&                     a_E_field,
                                  const int                         a_fluidVecComp,
                                  const Real                        a_time )
{
   CH_TIME("ExtendedMhdOp::accumulateExplicitRHS()");
  
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );

   accumulateMhdRHS( rhs_fluid, soln_fluid, a_time );
   if(m_viscLimiter>0.0) {
      accumulateViscRHS( rhs_fluid, soln_fluid );
   }
   if(m_heatLimiter>0.0) {
      accumulateHeatRHS( rhs_fluid, soln_fluid );
   }

}

void ExtendedMhdOp::accumulateImplicitRHS( FluidSpeciesPtrVect&        a_rhs,
                                  const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                  const FluidSpeciesPtrVect&        a_fluid_species,
                                  const PS::ScalarPtrVect&          a_scalars,
                                  const EField&                     a_E_field,
                                  const int                         a_fluidVecComp,
                                  const Real                        a_time )
{
   CH_TIME("ExtendedMhdOp::accumulateImplicitRHS()");
}

void ExtendedMhdOp::enforcePositivity( FluidSpecies&  a_species_comp )
{
   CH_TIME("ExtendedMhdOp::enforcePositivity()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   const double gamma = a_species_comp.m_gamma;
         LevelData<FArrayBox>& soln_rhoDen( a_species_comp.cell_var(0) );
   const LevelData<FArrayBox>& soln_momDen( a_species_comp.cell_var("momentumDensity") );
         LevelData<FArrayBox>& soln_eneDen( a_species_comp.cell_var("energyDensity") );
         LevelData<FArrayBox>& soln_eneDen_ele( a_species_comp.cell_var("energyDensity_ele") );
         
   // set mass density limits
   for (DataIterator dit(grids); dit.ok(); ++dit) {
            FArrayBox& JNh_on_patch = soln_rhoDen[dit];
      const FArrayBox& Jac_on_patch  = m_Jacobian[dit];
      const Box& grid_box( grids[dit] );
      FORT_SET_DENSITY_FLOOR( CHF_BOX(grid_box),
                              CHF_CONST_REAL(m_Nmin),
                              CHF_CONST_FRA1(Jac_on_patch,0),
                              CHF_FRA1(JNh_on_patch,0) );
   }         

   if(a_species_comp.m_evolve_electronDensity) {
      // set Zbar limits
      LevelData<FArrayBox>& soln_eleDen( a_species_comp.cell_var("electronDensity") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const FArrayBox& JNh_on_patch = soln_rhoDen[dit];
               FArrayBox& JNe_on_patch = soln_eleDen[dit];
         const Box& grid_box( grids[dit] );  // this box has no ghost
         const double Zmin = a_species_comp.m_Zmin;
         FORT_SET_ZBAR_LIMITS( CHF_BOX(grid_box),
                               CHF_CONST_REAL(Zmin),
                               CHF_CONST_FRA1(JNh_on_patch,0),
                               CHF_CONST_FRA1(JNe_on_patch,0) );
      }

      // set energy density limits using Ne=Zbar*Nh
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const FArrayBox& JNh_on_patch = soln_rhoDen[dit];
         const FArrayBox& JNe_on_patch = soln_eleDen[dit];
         const FArrayBox& Jmom_on_patch = soln_momDen[dit];
         const FArrayBox& Jac_on_patch  = m_Jacobian[dit];
               FArrayBox& Jenei_on_patch = soln_eneDen[dit];
               FArrayBox& Jenee_on_patch = soln_eneDen_ele[dit];

         const Box& grid_box( grids[dit] );  // this box has no ghost
         FORT_SET_ENERGY_FLOOR( CHF_BOX(grid_box),
                                CHF_CONST_REAL(m_Tmin),
                                CHF_CONST_REAL(m_Tmax),
                                CHF_CONST_REAL(gamma),
                                CHF_CONST_FRA1(Jac_on_patch,0),
                                CHF_CONST_FRA1(JNh_on_patch,0),
                                CHF_CONST_FRA1(JNe_on_patch,0),
                                CHF_CONST_FRA(Jmom_on_patch),
                                CHF_FRA1(Jenei_on_patch,0),
                                CHF_FRA1(Jenee_on_patch,0) );
      }
   }
   else {
      // set energy density limits using Ne=Nh
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const FArrayBox& JNh_on_patch = soln_rhoDen[dit];
         const FArrayBox& JNe_on_patch = soln_rhoDen[dit];
         const FArrayBox& Jmom_on_patch = soln_momDen[dit];
         const FArrayBox& Jac_on_patch  = m_Jacobian[dit];
               FArrayBox& Jenei_on_patch = soln_eneDen[dit];
               FArrayBox& Jenee_on_patch = soln_eneDen_ele[dit];

         const Box& grid_box( grids[dit] );  // this box has no ghost
         FORT_SET_ENERGY_FLOOR( CHF_BOX(grid_box),
                                CHF_CONST_REAL(m_Tmin),
                                CHF_CONST_REAL(m_Tmax),
                                CHF_CONST_REAL(gamma),
                                CHF_CONST_FRA1(Jac_on_patch,0),
                                CHF_CONST_FRA1(JNh_on_patch,0),
                                CHF_CONST_FRA1(JNe_on_patch,0),
                                CHF_CONST_FRA(Jmom_on_patch),
                                CHF_FRA1(Jenei_on_patch,0),
                                CHF_FRA1(Jenee_on_patch,0) );
      }
   }

}

void ExtendedMhdOp::setDivMassFluxBinary( LevelData<FArrayBox>&  a_divMassFluxBinary,
                                    const LevelData<FArrayBox>&  a_divMassFlux ) const
{
   CH_TIME("ExtendedMhdOp::setDivMassFluxBinary()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // divMassFluxBinary is set to 1 if mass density is above the user specified
   // minimum value. If the mass density is equal to or below the mininim, then
   // divmassFluxBinary is set to 0 if divMassFlux is positive, which would cause
   // the mass density in a give cell to decrease further.

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FArrayBox& rhoDen_on_patch = m_rhoDen_cc[dit];
      const FArrayBox& divMF_on_patch  = a_divMassFlux[dit];
            FArrayBox& divMFB_on_patch = a_divMassFluxBinary[dit];
      const Box& grid_box = divMFB_on_patch.box();
      FORT_SET_DIV_MASS_FLUX_BINARY( CHF_BOX(grid_box),
                                     CHF_CONST_REAL(m_Nmin),
                                     CHF_CONST_FRA1(rhoDen_on_patch,0),
                                     CHF_CONST_FRA1(divMF_on_patch,0),
                                     CHF_FRA1(divMFB_on_patch,0) );
   }         

}

void ExtendedMhdOp::postStageEval( FluidSpecies&  a_species_comp,
                                   FluidSpecies&  a_species_phys,
                                   const Real     a_dt,
                                   const Real     a_time,
                                   const int      a_stage )
{
   CH_TIME("ExtendedMhdOp::postStageEval()");
   //if (!procID()) {
   //   cout << "JRA: ExtendedMhdOp::postStageEval() a_stage = " << a_stage << endl;
   //   cout << "JRA: ExtendedMhdOp::postStageEval() a_time  = " << a_time << endl;
   //   cout << "JRA: ExtendedMhdOp::postStageEval() a_dt    = " << a_dt << endl;
   //}

   // update relaxed variables for stage>0 
   // (advance for stage=0 done from direct call to postTimeEval())
   //
   if(a_stage>0 || a_time==0) { // a_time==0 is to set initial value
      postTimeEval(a_species_comp, a_species_phys, a_dt, a_time, a_stage);
   }

}

void ExtendedMhdOp::postTimeEval( FluidSpecies&  a_species_comp,
                                  FluidSpecies&  a_species_phys,
                                  const Real     a_dt,
                                  const Real     a_time,
                                  const int      a_stage )
{
   CH_TIME("ExtendedMhdOp::postTimeEval()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // update some cell-center variables and update collision terms
   setCellCenterValues( a_species_phys );
   SpaceUtils::upWindToFaces(m_rhoDen_cf,m_rhoDen_cc,m_rhoDen_cf,"c2");
   SpaceUtils::upWindToFaces(m_eleDen_cf,m_eleDen_cc,m_eleDen_cf,"c2");
   updateCollisionTerms( a_species_phys, m_pressure_ele, m_pressure_ion );

   // set insulator BC to resistivity
   if(a_species_phys.m_evolve_magneticField_virtual && m_By_insulator_bc) {
      m_etaRes_BC.exchange();
      int this_cbc = m_magneticField_virtual_cbc;
      m_fluid_bc.at(this_cbc)->setInsulatorBC( a_species_phys, m_etaRes, m_etaRes_BC, a_time );
   }
   SpaceUtils::interpCellToEdges(m_etaRes_ce, m_etaRes, m_etaRes_ce, "c2");  

   /////////////////

   if(a_species_phys.m_evolve_currentDensity==1) {

      // first-order backward euler advance of J and E
      CH_assert(a_species_phys.m_evolve_currentDensity==1);
      CH_assert(a_species_phys.m_evolve_electricField==1);
      int this_ebc;

      LevelData<EdgeDataBox>& soln_Jcovar( a_species_phys.edge_var("currentDensity") );
      LevelData<EdgeDataBox>& soln_Ecovar( a_species_phys.edge_var("electricField") );
      const LevelData<FArrayBox>& soln_By( a_species_phys.cell_var("magneticField_virtual") );

      // compute covariant gradient of Pe on edges
      m_geometry.computeMappedPoloidalGradientWithGhosts( m_pressure_ele, m_gradPe_covar_cc, 2);
      SpaceUtils::interpCellToEdges(m_gradPe_covar_ce, m_gradPe_covar_cc, m_gradPe_covar_ce, "c2");  
      SpaceUtils::interpCellToEdges(m_eleDen_ce, m_eleDen_cc, m_eleDen_ce,"c2");
      //SpaceUtils::exchangeEdgeDataBox(m_eleDen_ce); 
   
      // compute contra and covar components of virtual magnetic field
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& thisbox = m_By_covar_cc[dit].box();
         m_By_covar_cc[dit].copy( soln_By[dit], thisbox); // phys
         m_By_contr_cc[dit].copy( soln_By[dit], thisbox); // phys
         if(m_twoDaxisymm) {
            m_By_covar_cc[dit].mult(m_g_y[dit]);
            m_By_contr_cc[dit].divide(m_g_y[dit], thisbox, 0,0,1);
         }
      }

      // compute covariant ideal electric field, E0=-VxB, at cell edges
      if(m_updateHydroVars) {
         setMappedCspeed( a_species_phys.m_gamma, 1, 1 );
         computeIdealEatEdges(m_E0_ce,m_velocity,m_By_contr_cc,m_Cspeed_cc);
         this_ebc = m_electricField_ebc;
         m_fluid_edge_var_bc.at(this_ebc)->setEdgeBC( a_species_phys, m_E0_ce, m_EdgeBC_zeros, a_time );
         SpaceUtils::exchangeEdgeDataBox(m_E0_ce);
      }
      else {
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            for (int dir=0; dir<SpaceDim; dir++) m_E0_ce[dit][dir].setVal(0.0);
         }
      }

      // compute covariant J0=curl(By) at cell-edges and apply BC
      computeMappedCovariantCurl( m_J0_ce, m_By_covar_cc ); // Jacobian * curl(By)\cdot g_l (Ja*covar curl)
      m_geometry.divideJonEdges(m_J0_ce);                         // curl(By)_covar
      int this_cbc = m_magneticField_virtual_cbc;
      if(m_twoDaxisymm) m_fluid_bc.at(this_cbc)->setOnAxisCurlBC( m_J0_ce, soln_By, a_species_phys );
      SpaceUtils::exchangeEdgeDataBox(m_J0_ce);
      //
      this_ebc = m_currentDensity_ebc; // this needs to be done for Jmethod = 2  ...
      m_fluid_edge_var_bc.at(this_ebc)->applyEdgeBC( a_species_phys, m_J0_ce, a_time );

      if(a_time==0) { // set initial old values of J and E variables
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            m_Jold[dit].copy( soln_Jcovar[dit],0,0,soln_Jcovar.nComp() );
            m_Eold[dit].copy( soln_Ecovar[dit],0,0,soln_Ecovar.nComp() );
            m_gammaHallold_ce[dit].copy( m_gammaHall_ce[dit],0,0,m_gammaHall_ce.nComp() );
         }
         if(m_etaVis_ele>0.0) {
            computeElectronVelForVisc( m_Ve_ce, 0, soln_Jcovar, m_eleDen_ce, m_eleDen_cc );
            computeStrainTensorStag( m_Pie_diag_cn, m_Pie_xz_cc, m_Pie_yy_ce, m_Ve_ce);
         }
      }

      if(m_etaVis_ele>0.0) { // add div if ele visc tensor to pressure gradient
         computeDivOfElectronVisc( m_divPie_ce, m_Pie_diag_cn, m_Pie_xz_cc, m_Pie_yy_ce );
         addDivPieToGradPele( m_gradPe_covar_ce, m_divPie_ce );
      }

      // compute etaHall and Ehall at cell edges to use for implicit Hall solve
      if(m_modelDriftTerms) {
         computeHallEatEdges( m_etaHall_ce, m_etaHall_op, m_Ehall_ce, m_Ehall_op, 
                              m_J0_ce, m_By_contr_cc, 0 );
         //this_ebc = m_electricField_ebc;
         //m_fluid_edge_var_bc.at(this_ebc)->setEdgeBC( a_species_phys, m_Ehall_ce, m_EdgeBC_zeros, a_time );
      }

      // compute covariant div of stress tensor for Ohm's law on cell edges
      if(m_include_divJStress) {
         SpaceUtils::interpEdgesToCell(m_Jcc,m_J0_ce,"c2"); // curlB_covar at cell center
         computeElectronVelOnEdges( m_Ve_ce, soln_Jcovar, m_eleDen_ce, m_velocity );
         computeDivJStress( m_divJStress_covar, soln_Jcovar, m_Ve_ce, m_eleDen_ce,
                            m_velocity, m_Jcc, m_eleDen_cc );
      }

      // compute g-fun for implicit Hall solve on cell-edges 
      // and interpolate to opposite edges
      computeGfunEdges( m_gfun_covar, a_dt, m_E0_ce, m_gradPe_covar_ce, 
                        m_gammaHallold_ce, m_divJStress_covar,
                        m_J0_ce,  m_eleDen_ce, m_Eold, m_Jold );
      SpaceUtils::exchangeEdgeDataBox(m_gfun_covar); 
      this_ebc = m_currentDensity_ebc;
      m_fluid_edge_var_bc.at(this_ebc)->applyEdgeBC( a_species_phys, m_gfun_covar, a_time );
      SpaceUtils::interpEdgesToEdges(m_gfun_covar_op,m_gfun_covar,"c2");
 
      // advance J using fully implicit scheme or J0 for Ehall
      int Jmethod = 1; 
      if(m_useJ0forJcc || !m_fullyImplicitJ) Jmethod = 0;  
      if(!m_useRelaxedOhmsLaw) Jmethod = 2;
      LevelData<EdgeDataBox>& hallDriftScaleFactor_ce = m_dummyEdge_oneComp;
      for (int dir=0; dir<SpaceDim; dir++) {
         advanceEandJ( soln_Ecovar, soln_Jcovar, hallDriftScaleFactor_ce, m_gammaHall_ce,
                       a_dt, m_E0_ce,
                       m_etaRes_ce, m_Ehall_ce, m_Ehall_op,
                       m_etaHall_ce, m_etaHall_op,
                       m_gfun_covar, m_gfun_covar_op,
                       m_gradPe_covar_ce, m_J0_ce, m_eleDen_ce,
                       m_Eold, m_Jold, m_gammaHallold_ce, dir, Jmethod );
      }
      SpaceUtils::interpEdgesToCell(m_hallDriftScaleFactor,hallDriftScaleFactor_ce,"c2");
      SpaceUtils::exchangeEdgeDataBox(soln_Jcovar); 
      this_ebc = m_currentDensity_ebc;
      m_fluid_edge_var_bc.at(this_ebc)->applyEdgeBC( a_species_phys, soln_Jcovar, a_time );

      //////////////////////////////////////////////////////////////////////////////////////////
      //
      //    start iterative process for Hall term here
      //
      int iterHall = 0;
      Vector<double> local_normJres(SpaceDim,CH_BADVAL);
      Vector<double> normJres0(SpaceDim,1.0);
      double normJresMax(1.0e6);
      while ( normJresMax > m_HallResTol ) {
         if (m_useJ0forJcc || !m_modelDriftTerms) break; 
         if (!m_useRelaxedOhmsLaw) break; 
         if (iterHall==m_iterHallMax) break;

         // interp J from edges to cell center
         if(iterHall==0) {
            SpaceUtils::interpEdgesToCell(m_Jcc,soln_Jcovar,"c2"); // covar at cell center
            for (DataIterator dit(grids); dit.ok(); ++dit) {
               m_Jcc_res[dit].copy(m_Jcc[dit]);
            }
         }

         // compute covariant Ehall = Li0/Ne*JxB = J x etaHall
         computeHallEatEdges( m_etaHall_ce, m_etaHall_op, m_Ehall_ce, m_Ehall_op, 
                              soln_Jcovar, m_By_contr_cc, iterHall );
         //this_ebc = m_electricField_ebc;
         //m_fluid_edge_var_bc.at(this_ebc)->setEdgeBC( a_species_phys, m_Ehall_ce, m_EdgeBC_zeros, a_time );

         // advance E and J on edges using relaxation scheme
         Jmethod = 0; // 0 uses Ehall since trying to achieve self-consistent solution with J living on edges
         LevelData<EdgeDataBox>& hallDriftScaleFactor_ce = m_dummyEdge_oneComp;
         for (int dir=0; dir<SpaceDim; dir++) {
            advanceEandJ( soln_Ecovar, soln_Jcovar, hallDriftScaleFactor_ce, m_gammaHall_ce,
                          a_dt, m_E0_ce,
                          m_etaRes_ce, m_Ehall_ce, m_Ehall_op, 
                          m_etaHall_ce, m_etaHall_op,
                          m_gfun_covar, m_gfun_covar_op,
                          m_gradPe_covar_ce, m_J0_ce, m_eleDen_ce,
                          m_Eold, m_Jold, m_gammaHallold_ce, dir, Jmethod );
         }
         SpaceUtils::interpEdgesToCell(m_hallDriftScaleFactor,hallDriftScaleFactor_ce,"c2");

         // fill ghost cells and perform exchange
         SpaceUtils::exchangeEdgeDataBox(soln_Jcovar); 
         this_ebc = m_currentDensity_ebc;
         m_fluid_edge_var_bc.at(this_ebc)->applyEdgeBC( a_species_phys, soln_Jcovar, a_time );

         // update Jcc with updated soln_Jcovar and compute residual matrix
         SpaceUtils::interpEdgesToCell(m_Jcc,soln_Jcovar,"c2"); // covar at cell center
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            m_Jcc_res[dit].minus(m_Jcc[dit]);
         }
  
         // calculate residual
         for (int dir=0; dir<SpaceDim; dir++) {
            local_normJres[dir]=0.0;
            for (DataIterator dit(grids); dit.ok(); ++dit) {
               local_normJres[dir] += m_Jcc_res[dit].norm(2,dir,1);
            }
         }
   
         Vector<double> normJres( local_normJres );
         normJresMax=0.0;
         for (int dir=0; dir<SpaceDim; ++dir) {
#ifdef CH_MPI
            MPI_Allreduce( &local_normJres[dir], &normJres[dir], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#endif
            if(iterHall==0 && normJres[dir]>1.0) normJres0[dir] = normJres[dir];
            normJres[dir] /= normJres0[dir];
            normJresMax = Max( normJres[dir], normJresMax );
         }
         if(!procID() && iterHall>=10) {
            cout << "iterHall = " << iterHall << endl;
            cout << "normJres[0] = " << normJres[0] << endl;
            cout << "normJres[1] = " << normJres[1] << endl;
         }

         // reset residual matrix for next iteration and update iteration value
         iterHall += 1;
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            m_Jcc_res[dit].copy(m_Jcc[dit]);
         }
  
      }
      //
      //    end iterative process for hall term here
      //
      //////////////////////////////////////////////////////////////////////////////////////////
   
      SpaceUtils::exchangeEdgeDataBox(soln_Ecovar);
      SpaceUtils::exchangeEdgeDataBox(m_gammaHall_ce);
      this_ebc = m_electricField_ebc;
      m_fluid_edge_var_bc.at(this_ebc)->applyEdgeBC( a_species_phys, soln_Ecovar, a_time );

      // set E and J comp and phys (same for edge vars, both covar)
      //
      LevelData<EdgeDataBox>& soln_Jcomp( a_species_comp.edge_var("currentDensity") );
      LevelData<EdgeDataBox>& soln_Ecomp( a_species_comp.edge_var("electricField") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         SpaceUtils::copyEdgeDataBox(soln_Jcomp[dit],soln_Jcovar[dit]);
         SpaceUtils::copyEdgeDataBox(soln_Ecomp[dit],soln_Ecovar[dit]);
         if(a_stage==0) { // update old values
            SpaceUtils::copyEdgeDataBox(m_Jold[dit],soln_Jcomp[dit]);
            SpaceUtils::copyEdgeDataBox(m_Eold[dit],soln_Ecomp[dit]);
            SpaceUtils::copyEdgeDataBox(m_gammaHallold_ce[dit],m_gammaHall_ce[dit]);
         }
      }

      // compute charge density at nodes from electric field
      // and interpolate to cell center for diagnostic
      //
      m_geometry.physicalDivergenceFromEdgeTans(m_divE, soln_Ecovar);
      SpaceUtils::interpNodesToEdges(m_deltaN_ce, m_divE, "c2");
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; dir++) {
            FArrayBox& thisDeltaN( m_deltaN_ce[dit][dir] ); 
            thisDeltaN.mult( m_delta0*m_lambdai0/m_epRel );
         }
      }
      SpaceUtils::interpEdgesToCell(m_deltaN_cc, m_deltaN_ce, "c2");
  
   
      // advance electron viscosity flux 
      // (This needs to be just after J update)
      // 
      if(m_etaVis_ele>0.0) {
         computeElectronVelForVisc( m_Ve_ce, 0, soln_Jcovar, m_eleDen_ce, m_eleDen_cc );
         advanceElectronViscFluxes( m_Pie_diag_cn, m_Pie_xz_cc, m_Pie_yy_ce, 
                                    m_Ve_ce, a_stage, a_time, a_dt ); 
      
      }   

   } // end if statement for J and E advance


   //////////////////////////////////////////////////////////////////
   //
   //  advance collisional transport coefficents 
   //  using relaxation scheme
   //
   //
   if(m_viscLimiter>0.0 && m_useRelaxedViscosity) {

      // compute the relaxed viscosity fluxes on cell faces
      setIonViscosityCoefficients( m_taui, m_xion, m_pressure_ion );
      
      SpaceUtils::upWindToFaces(m_velocity_cf,m_velocity,m_velocity_cf,"c2");
      
      computeViscFluxes( m_m0JaFluxVisc_cf0, m_m1JaFluxVisc_cf0, 
                         m_momVisc_source0, 
                         m_velocity_cf, m_velocity,
                         m_etaVisC_cf,  m_etaVisC );
 
      // set initial and old values if t=0
      if(a_time==0) {
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            //m_m0JaFluxVisc_cf0[dit].setVal(0.0); JRA, not for shear flow
            //m_m1JaFluxVisc_cf0[dit].setVal(0.0); JRA, not for shear flow
            m_m0JaFluxVisc_cf[dit].copy( m_m0JaFluxVisc_cf0[dit],0,0,SpaceDim );
            m_m1JaFluxVisc_cf[dit].copy( m_m1JaFluxVisc_cf0[dit],0,0,SpaceDim );
            m_momVisc_source[dit].copy( m_momVisc_source0[dit] );
            m_m0JaFluxVisc_cf_old[dit].copy( m_m0JaFluxVisc_cf[dit],0,0,SpaceDim );
            m_m1JaFluxVisc_cf_old[dit].copy( m_m1JaFluxVisc_cf[dit],0,0,SpaceDim );
            m_momVisc_source_old[dit].copy( m_momVisc_source[dit] );
         }
      }
     
      // update viscosity fluxes using relaxation scheme
      for (DataIterator dit(grids); dit.ok(); ++dit) {

         for (int dir=0; dir<SpaceDim; dir++) {

            const FArrayBox& this_N_dir     = m_rhoDen_cf[dit][dir];
            const FArrayBox& this_Dx_dir    = m_DiffVx_cf[dit][dir];
            const FArrayBox& this_Dz_dir    = m_DiffVz_cf[dit][dir];
            const FArrayBox& this_Pxold_dir = m_m0JaFluxVisc_cf_old[dit][dir];
            const FArrayBox& this_Pzold_dir = m_m1JaFluxVisc_cf_old[dit][dir];
            const FArrayBox& this_Px0_dir   = m_m0JaFluxVisc_cf0[dit][dir];
            const FArrayBox& this_Pz0_dir   = m_m1JaFluxVisc_cf0[dit][dir];
                  FArrayBox& this_Px_dir    = m_m0JaFluxVisc_cf[dit][dir];
                  FArrayBox& this_Pz_dir    = m_m1JaFluxVisc_cf[dit][dir];

            Box facebox( grids[dit] ); // no ghost 
            facebox.surroundingNodes( dir );
            for (int n=0; n<this_Px_dir.nComp(); n++) {
               FORT_ADVANCE_USING_RELAX( CHF_BOX(facebox),
                                         CHF_CONST_REAL( a_dt ),
                                         CHF_CONST_REAL( m_delta0 ),
                                         CHF_CONST_FRA1( this_N_dir,0 ),
                                         CHF_CONST_FRA1( this_Dx_dir,n ),
                                         CHF_CONST_FRA1( this_Pxold_dir,n ),
                                         CHF_CONST_FRA1( this_Px0_dir,n ),
                                         CHF_FRA1( this_Px_dir,n ) );
               // 
               FORT_ADVANCE_USING_RELAX( CHF_BOX(facebox),
                                         CHF_CONST_REAL( a_dt ),
                                         CHF_CONST_REAL( m_delta0 ),
                                         CHF_CONST_FRA1( this_N_dir,0 ),
                                         CHF_CONST_FRA1( this_Dz_dir,n ),
                                         CHF_CONST_FRA1( this_Pzold_dir,n ),
                                         CHF_CONST_FRA1( this_Pz0_dir,n ),
                                         CHF_FRA1( this_Pz_dir,n ) );
            }
            
         }
         SpaceUtils::exchangeFluxBox(m_m0JaFluxVisc_cf);
         SpaceUtils::exchangeFluxBox(m_m1JaFluxVisc_cf);
         
         if(a_stage==0) { // update old values
            m_m0JaFluxVisc_cf_old[dit].copy( m_m0JaFluxVisc_cf[dit],0,0,SpaceDim );
            m_m1JaFluxVisc_cf_old[dit].copy( m_m1JaFluxVisc_cf[dit],0,0,SpaceDim );
         }
         
         // advance source term for axisymmetric
         if(m_twoDaxisymm) {
            
            const FArrayBox& this_N    = m_rhoDen_cc[dit];
            const FArrayBox& this_Dth  = m_DiffVth[dit];
            const FArrayBox& this_Sold = m_momVisc_source_old[dit];
            const FArrayBox& this_S0   = m_momVisc_source0[dit];

            FArrayBox& this_S = m_momVisc_source[dit];

            const Box& cellbox = this_S.box();
            FORT_ADVANCE_USING_RELAX( CHF_BOX(cellbox),
                                      CHF_CONST_REAL( a_dt ),
                                      CHF_CONST_REAL( m_delta0 ),
                                      CHF_CONST_FRA1( this_N,0 ),
                                      CHF_CONST_FRA1( this_Dth,0 ),
                                      CHF_CONST_FRA1( this_Sold,0 ),
                                      CHF_CONST_FRA1( this_S0,0 ),
                                      CHF_FRA1( this_S,0 ) );
              
            if(a_stage==0) { // update old values
               m_momVisc_source_old[dit].copy( m_momVisc_source[dit] );
            }
            
         }
         
      }

      // compute the viscosity flux for energy density 
      computeViscEnergyFlux( m_enJaFluxVisc_cf,
                             m_m0JaFluxVisc_cf,
                             m_m1JaFluxVisc_cf,
                             m_velocity_cf );

   }
   

   //////////////////////////////////////////////////////////////////
   //
   //  advance collisional heat flux 
   //  using relaxation scheme
   //
   //
   if(m_heatLimiter>0.0 && m_useRelaxedHeatFlux) {
   
         // calculate temperature for each species at cell center and 
         // cell faces
         SpaceUtils::upWindToFaces(m_Tele_cf,m_Tele_cc,m_Tele_cf,"c2");
         SpaceUtils::upWindToFaces(m_Tion_cf,m_Tion_cc,m_Tion_cf,"c2");

         // compute the relaxed heat fluxes on cell faces for each species
         const double gamma = a_species_phys.m_gamma;
         setHeatFluxCoefficients( gamma, 
                                  m_taue, m_taui,
                                  m_xele, m_xion,
                                  m_pressure_ele, m_pressure_ion );
         computeCollHeatFlux( m_qe0, m_Dqe, m_Tele_cc, m_Tele_cf, m_eleDen_cf, m_kappae_cf, 
                             -1, m_taueMin, m_heatLimiter_ele, m_Ve0sq, gamma);
         computeCollHeatFlux( m_qi0, m_Dqi, m_Tion_cc, m_Tion_cf, m_rhoDen_cf, m_kappai_cf, 
                              1, m_tauiMin, m_heatLimiter, m_Vi0sq, gamma);
      
      // set initial and old values if t=0
      if(a_time==0) {
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            m_qi_cf[dit].copy( m_qi0[dit],0,0,SpaceDim );
            m_qe_cf[dit].copy( m_qe0[dit],0,0,SpaceDim );
            m_qi_old[dit].copy( m_qi_cf[dit],0,0,SpaceDim );
            m_qe_old[dit].copy( m_qe_cf[dit],0,0,SpaceDim );
         }
      }
     
      // update coll heat fluxes using relaxation scheme
      for (DataIterator dit(grids); dit.ok(); ++dit) {

         for (int dir=0; dir<SpaceDim; dir++) {

            const FArrayBox& this_N_dir     = m_rhoDen_cf[dit][dir];
            const FArrayBox& this_Ne_dir    = m_eleDen_cf[dit][dir];
            const FArrayBox& this_Dqe_dir   = m_Dqe[dit][dir];
            const FArrayBox& this_Dqi_dir   = m_Dqi[dit][dir];
            const FArrayBox& this_qeold_dir = m_qe_old[dit][dir];
            const FArrayBox& this_qiold_dir = m_qi_old[dit][dir];
            const FArrayBox& this_qe0_dir   = m_qe0[dit][dir];
            const FArrayBox& this_qi0_dir   = m_qi0[dit][dir];
                  FArrayBox& this_qe_dir    = m_qe_cf[dit][dir];
                  FArrayBox& this_qi_dir    = m_qi_cf[dit][dir];

            Box facebox( grids[dit] ); // no ghost 
            facebox.surroundingNodes( dir );
            for (int n=0; n<this_qe_dir.nComp(); n++) {
               FORT_ADVANCE_USING_RELAX( CHF_BOX(facebox),
                                         CHF_CONST_REAL( a_dt ),
                                         CHF_CONST_REAL( m_delta0 ),
                                         CHF_CONST_FRA1( this_Ne_dir,0 ),
                                         CHF_CONST_FRA1( this_Dqe_dir,n ),
                                         CHF_CONST_FRA1( this_qeold_dir,n ),
                                         CHF_CONST_FRA1( this_qe0_dir,n ),
                                         CHF_FRA1( this_qe_dir,n ) );
               // 
               FORT_ADVANCE_USING_RELAX( CHF_BOX(facebox),
                                         CHF_CONST_REAL( a_dt ),
                                         CHF_CONST_REAL( m_delta0 ),
                                         CHF_CONST_FRA1( this_N_dir,0 ),
                                         CHF_CONST_FRA1( this_Dqi_dir,n ),
                                         CHF_CONST_FRA1( this_qiold_dir,n ),
                                         CHF_CONST_FRA1( this_qi0_dir,n ),
                                         CHF_FRA1( this_qi_dir,n ) );
            }
            
         }
         SpaceUtils::exchangeFluxBox(m_qe_cf);
         SpaceUtils::exchangeFluxBox(m_qi_cf);
         
         if(a_stage==0) { // update old values
            m_qe_old[dit].copy( m_qe_cf[dit],0,0,SpaceDim );
            m_qi_old[dit].copy( m_qi_cf[dit],0,0,SpaceDim );
         }
      }

   }         

}

void ExtendedMhdOp::getMemberVar( LevelData<FArrayBox>&  a_Var,
                            const FluidSpecies&          a_fluid_species,
                            const string&                a_name ) const
{
   const DisjointBoxLayout& grids( m_geometry.grids() );
   const IntVect ghostVect = 1*IntVect::Unit;
   //const IntVect& ghostVect( soln_data.ghostVect() );
   
   if(a_name=="chargeDensity") {
      a_Var.define(grids, m_deltaN_cc.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_deltaN_cc[dit]);
      }
   }
   else if(a_name=="nue_izn") {
      a_Var.define(grids, m_nue_izn.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_nue_izn[dit]);
      }
   }
   else if(a_name=="gradVeperp") {
      a_Var.define(grids, m_gradVeperp_cc.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_gradVeperp_cc[dit]);
      }
   }
   else if(a_name=="gradVepara") {
      a_Var.define(grids, m_gradVepara_cc.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_gradVepara_cc[dit]);
      }
   }
   else if(a_name=="kappae_wedge") {
      a_Var.define(grids, m_kappae_wedge.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_kappae_wedge[dit]);
      }
   }
   else if(a_name=="kappai_wedge") {
      a_Var.define(grids, m_kappai_wedge.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_kappai_wedge[dit]);
      }
   }
   else if(a_name=="kappae") {
      a_Var.define(grids, m_kappae.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_kappae[dit]);
      }
   }
   else if(a_name=="kappai") {
      a_Var.define(grids, m_kappai.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_kappai[dit]);
      }
   }
   else if(a_name=="etaVisC") {
      a_Var.define(grids, m_etaVisC.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_etaVisC[dit]);
      }
   }
   else if(a_name=="etaVisG") {
      a_Var.define(grids, m_etaVisG.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_etaVisG[dit]);
      }
   }
   else if(a_name=="qi_wedge") {
      a_Var.define(grids, m_qi_wedge_cc.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_qi_wedge_cc[dit]);
      }
   }
   else if(a_name=="gfun") {
      a_Var.define(grids, SpaceDim, ghostVect);
      SpaceUtils::interpEdgesToCell(a_Var,m_gfun_covar,"c2");
      m_geometry.convertPhysToCovar(a_Var,1);
   }
   else if(a_name=="gfun_op") {
      a_Var.define(grids, SpaceDim, ghostVect);
      SpaceUtils::interpEdgesToCell(a_Var,m_gfun_covar_op,"c2");
      m_geometry.convertPhysToCovar(a_Var,1);
   }
   else if(a_name=="VhallDrift") {
      a_Var.define(grids, SpaceDim, IntVect::Zero);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_VhallDrift_cc[dit]);
      }
   }
   else if(a_name=="hallDriftScaleFactor") {
      a_Var.define(grids, 1, IntVect::Zero);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_hallDriftScaleFactor[dit]);
      }
   }
   else if(a_name=="Ehall") {
      a_Var.define(grids, SpaceDim, ghostVect);
      SpaceUtils::interpEdgesToCell(a_Var,m_Ehall_ce,"c2");
      m_geometry.convertPhysToCovar(a_Var,1);
   }
   else if(a_name=="Eideal") {
      a_Var.define(grids, SpaceDim, ghostVect);
      SpaceUtils::interpEdgesToCell(a_Var,m_E0_ce,"c2");
      m_geometry.convertPhysToCovar(a_Var,1);
   }
   else if(a_name=="curlB") {
      a_Var.define(grids, SpaceDim, ghostVect);
      SpaceUtils::interpEdgesToCell(a_Var,m_J0_ce,"c2"); // covar at cell center
      m_geometry.convertPhysToCovar(a_Var,1);
   }
   else if(a_name=="etaHall") {
      a_Var.define(grids, SpaceDim, ghostVect);
      SpaceUtils::interpEdgesToCell(a_Var,m_etaHall_ce,"c2");
   }
   else if(a_name=="divJStress") {
      a_Var.define(grids, SpaceDim, ghostVect);
      SpaceUtils::interpEdgesToCell(a_Var,m_divJStress_covar,"c2");
      m_geometry.convertPhysToCovar(a_Var,1);
   }
   else if(a_name=="qi_perp") {
      a_Var.define(grids, SpaceDim, ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         m_dummyFlux_oneComp[dit].copy(m_qi_norm[dit]);
         const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         const RealVect& faceArea = coord_sys.getMappedFaceArea();
         for (int dir=0; dir<SpaceDim; dir++) {
            FArrayBox& this_face = m_dummyFlux_oneComp[dit][dir];
            this_face.divide(faceArea[dir]);
         }
      }
      SpaceUtils::interpFaceVectorToCell(a_Var,m_dummyFlux_oneComp,"c2");
   }
   else if(a_name=="qe_perp") {
      a_Var.define(grids, SpaceDim, ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         m_dummyFlux_oneComp[dit].copy(m_qe_norm[dit]);
         const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         const RealVect& faceArea = coord_sys.getMappedFaceArea();
         for (int dir=0; dir<SpaceDim; dir++) {
            FArrayBox& this_face = m_dummyFlux_oneComp[dit][dir];
            this_face.divide(faceArea[dir]);
         }
      }
      SpaceUtils::interpFaceVectorToCell(a_Var,m_dummyFlux_oneComp,"c2");
   }
   else if(a_name=="etaRes") {
      a_Var.define(grids, m_etaRes.nComp(), ghostVect);
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         a_Var[dit].copy(m_etaRes[dit]);
      }
   }
   else if(a_name=="JaW_cf_0") {
      a_Var.define(grids, m_JaW_cf.nComp(), IntVect::Unit);
      const int dir=0;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         //Box& cellbox( a_Var[dit].box() );
         //cellbox.surroundingNodes( dir ); // convert to node type in dir direction
         const FArrayBox& thisJaW = m_JaW_cf[dit][dir];
         //a_Var[dit].copy(thisJaW,cellbox);
         SpaceUtils::copyNodeToCell(a_Var[dit],thisJaW,a_Var[dit].box(),dir);
         //cellbox.enclosedCells( dir ); // convert back to cell type in dir direction
      }
      a_Var.exchange();
   }
   else if(a_name=="JaW_cf_1") {
      a_Var.define(grids, m_JaW_cf.nComp(), IntVect::Unit);
      const int dir = 1;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         //Box& cellbox( a_Var[dit].box() );
         //cellbox.surroundingNodes( dir ); // convert to node type in dir direction
         const FArrayBox& thisJaW = m_JaW_cf[dit][dir];
         //a_Var[dit].copy(thisJaW,cellbox);
         SpaceUtils::copyNodeToCell(a_Var[dit],thisJaW,a_Var[dit].box(),dir);
         //cellbox.enclosedCells( dir ); // convert back to cell type in dir direction
      }
      a_Var.exchange();
   }
   else if(a_name=="divPie_0") {
      a_Var.define(grids, m_divPie_ce.nComp(), IntVect::Unit);
      const int dir = 0;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         const FArrayBox& thisdivPie = m_divPie_ce[dit][dir];
         SpaceUtils::copyNodeToCell(a_Var[dit],thisdivPie,a_Var[dit].box(),1);
      }
      a_Var.exchange();
   }
   else if(a_name=="divPie_1") {
      a_Var.define(grids, m_divPie_ce.nComp(), IntVect::Unit);
      const int dir = 1;
      for (DataIterator dit(a_Var.dataIterator()); dit.ok(); ++dit) {
         const FArrayBox& thisdivPie = m_divPie_ce[dit][dir];
         SpaceUtils::copyNodeToCell(a_Var[dit],thisdivPie,a_Var[dit].box(),0);
      }
      a_Var.exchange();
   }
   else {
      cout << "Trying to getMemberVar() from ExtendedMhdOp for a_name = " <<
              a_name << endl;
      MayDay::Error("Trying to getMemberVar() from ExtendedMhdOp!!!  ");
   }

}


void ExtendedMhdOp::defineBlockPC(std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>&  a_pc,
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
  if (a_im) {
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


void ExtendedMhdOp::updateBlockPC(  std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>&  a_pc,
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
    CH_assert(m_my_pc_idx >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx);
  
    if (!procID()) {
      std::cout << "    ==> Updating " << _FLUID_OP_PC_ << " preconditioner " 
                << " for ExtendedMhd fluid species " << a_species_idx << ".\n";
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


void ExtendedMhdOp::updatePCImEx(const FluidSpeciesPtrVect&     a_fluid_species, 
                              const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                              const double                      a_time,
                              const int                         a_step,
                              const int                         a_stage,
                              const double                      a_shift,
                              const int                         a_component,
                              const std::string& )
{
   CH_TIME("ExtendedMhdOp::updatePCImEx()");
     
   MayDay::Error("Error in ExtendedMhdOp::updatePCImEx()" 
                 "This function shouldn't have been called at all!"
                 "No implicit for ExtendedMhdOp()!");
   
}


void ExtendedMhdOp::solvePCImEx( FluidSpeciesPtrVect&     a_fluid_species_solution,
                        const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                        const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                        const std::string&,
                        const int                         a_component )
{
   CH_TIME("ExtendedMhdOp::solvePCImEx()");
   
   MayDay::Error("Error in ExtendedMhdOp::solvePCImEx()" 
                 "This function shouldn't have been called at all!"
                 "No implicit for ExtendedMhdOp()!");
   
}

void ExtendedMhdOp::computeGfunEdges( LevelData<EdgeDataBox>&  a_gfun_covar,
                                const double                   a_dt,
                                const LevelData<EdgeDataBox>&  a_E0_covar,
                                const LevelData<EdgeDataBox>&  a_gradPe_covar,
                                const LevelData<EdgeDataBox>&  a_gammaHallold_ce,
                                const LevelData<EdgeDataBox>&  a_divJStress_covar,
                                const LevelData<EdgeDataBox>&  a_J0_covar,
                                const LevelData<EdgeDataBox>&  a_eleDen_ce,
                                const LevelData<EdgeDataBox>&  a_Eold,
                                const LevelData<EdgeDataBox>&  a_Jold ) const
{
   CH_TIME("ExtendedMhdOp::computeGfunEdges()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // compute covar g on edges needed to solve the following equation
   // alpha*J = g + wlh x J

   int useFiniteIon = 0;
   if(m_modelDriftTerms) useFiniteIon = 1;
 
   for (DataIterator dit(grids); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; dir++) {

         const FArrayBox& this_Eideal_dir   = a_E0_covar[dit][dir];
         const FArrayBox& this_gradPe_dir   = a_gradPe_covar[dit][dir];
         const FArrayBox& this_gammaHallold_dir = a_gammaHallold_ce[dit][dir];
         const FArrayBox& this_divJStress_dir = a_divJStress_covar[dit][dir];
         const FArrayBox& this_J0_dir     = a_J0_covar[dit][dir];
         const FArrayBox& this_Ne_dir     = a_eleDen_ce[dit][dir];
         const FArrayBox& this_Eold_dir   = a_Eold[dit][dir];
         const FArrayBox& this_Jold_dir   = a_Jold[dit][dir];   
 
         FArrayBox& this_gfun_dir = a_gfun_covar[dit][dir];

         Box edgebox(grids[dit]); // no ghosts
         edgebox.surroundingNodes();
         edgebox.enclosedCells(dir);

         FORT_COMPUTE_GFUN_EDGES( CHF_BOX(edgebox),
                                  CHF_CONST_INT(useFiniteIon),
                                  CHF_CONST_REAL( a_dt ),
                                  CHF_CONST_REAL( m_epsilon0 ),
                                  CHF_CONST_REAL( m_delta0 ),
                                  CHF_CONST_REAL( m_lambdai0 ),
                                  CHF_CONST_FRA1( this_Eideal_dir, 0 ),
                                  CHF_CONST_FRA1( this_gradPe_dir, 0 ),
                                  CHF_CONST_FRA1( this_gammaHallold_dir, 0 ),
                                  CHF_CONST_FRA1( this_divJStress_dir, 0 ),
                                  CHF_CONST_FRA1( this_J0_dir, 0 ),
                                  CHF_CONST_FRA1( this_Ne_dir, 0 ),
                                  CHF_CONST_FRA1( this_Eold_dir, 0 ),
                                  CHF_CONST_FRA1( this_Jold_dir, 0 ),
                                  CHF_FRA1( this_gfun_dir, 0 ) );
      }
   }

}

void ExtendedMhdOp::advanceEandJ( LevelData<EdgeDataBox>&  soln_Ecovar,
                                  LevelData<EdgeDataBox>&  soln_Jcovar,
                                  LevelData<EdgeDataBox>&  a_hallDriftScaleFactor_ce,
                                  LevelData<EdgeDataBox>&  a_gammaHall_ce,
                            const double                   a_dt,
                            const LevelData<EdgeDataBox>&  a_E0_covar,
                            const LevelData<EdgeDataBox>&  a_eta_ce,
                            const LevelData<EdgeDataBox>&  a_Ehall_covar,
                            const LevelData<EdgeDataBox>&  a_Ehall_op,
                            const LevelData<EdgeDataBox>&  a_etaHall_ce,
                            const LevelData<EdgeDataBox>&  a_etaHall_op,
                            const LevelData<EdgeDataBox>&  a_gfun_covar,
                            const LevelData<EdgeDataBox>&  a_gfun_covar_op,
                            const LevelData<EdgeDataBox>&  a_gradPe_covar,
                            const LevelData<EdgeDataBox>&  a_J0_covar,
                            const LevelData<EdgeDataBox>&  a_eleDen_ce,
                            const LevelData<EdgeDataBox>&  a_Eold,
                            const LevelData<EdgeDataBox>&  a_Jold,
                            const LevelData<EdgeDataBox>&  a_gammaHallold_ce,
                            const int                      dir,
                            const int                      a_method ) const
{
   CH_TIME("ExtendedMhdOp::advanceEandJ()");
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // advance E and J on edges using relaxation scheme

   int useFiniteIon = 0;
   if(m_modelDriftTerms) useFiniteIon = 1;
 
   for (DataIterator dit(grids); dit.ok(); ++dit) {
         
      Box edgebox( grids[dit] ); // no ghost 
      edgebox.surroundingNodes( );
      edgebox.enclosedCells( dir );
        
      const FArrayBox& this_gammaHallold_dir = a_gammaHallold_ce[dit][dir];
      const FArrayBox& this_dxidX_dir    = m_dxidX_ce[dit][dir];
 
      const FArrayBox& this_Eideal_dir  = a_E0_covar[dit][dir];
      const FArrayBox& this_eta_dir     = a_eta_ce[dit][dir];
      const FArrayBox& this_Ehall_dir   = a_Ehall_covar[dit][dir];
      const FArrayBox& this_EhallOp_dir = a_Ehall_op[dit][dir];
      const FArrayBox& this_etaHall_dir   = a_etaHall_ce[dit][dir];
      const FArrayBox& this_etaHallOp_dir = a_etaHall_op[dit][dir];
      const FArrayBox& this_gfun_dir    = a_gfun_covar[dit][dir];
      const FArrayBox& this_gfun_op_dir = a_gfun_covar_op[dit][dir];
      const FArrayBox& this_gradPe_dir  = a_gradPe_covar[dit][dir];
      const FArrayBox& this_J0_dir     = a_J0_covar[dit][dir];
      const FArrayBox& this_Ne_dir     = a_eleDen_ce[dit][dir];
      const FArrayBox& this_Eold_dir   = a_Eold[dit][dir];
      const FArrayBox& this_Jold_dir   = a_Jold[dit][dir];

      FArrayBox& this_EonEdges_dir = soln_Ecovar[dit][dir];
      FArrayBox& this_JonEdges_dir = soln_Jcovar[dit][dir];
      FArrayBox& this_hallDriftSF_dir = a_hallDriftScaleFactor_ce[dit][dir];
      FArrayBox& this_gammaHall_dir = a_gammaHall_ce[dit][dir];

      FORT_ADVANCE_EANDJ_USING_RELAX( CHF_BOX(edgebox),
                                CHF_CONST_INT( dir ),
                                CHF_CONST_INT( a_method ),
                                CHF_CONST_INT( useFiniteIon ),
                                CHF_CONST_REAL( a_dt ),
                                CHF_CONST_REAL( m_epsilon0 ),
                                CHF_CONST_REAL( m_delta0 ),
                                CHF_CONST_REAL( m_lambdai0 ),
                                CHF_CONST_FRA( this_dxidX_dir ),
                                CHF_CONST_FRA1( this_Eideal_dir, 0 ),
                                CHF_CONST_FRA1( this_eta_dir, 0 ),
                                CHF_CONST_FRA1( this_Ehall_dir, 0 ),
                                CHF_CONST_FRA1( this_EhallOp_dir, 0 ),
                                CHF_CONST_FRA1( this_etaHall_dir, 0 ),
                                CHF_CONST_FRA1( this_etaHallOp_dir, 0 ),
                                CHF_CONST_FRA1( this_gfun_dir, 0 ),
                                CHF_CONST_FRA1( this_gfun_op_dir, 0 ),
                                CHF_CONST_FRA1( this_gradPe_dir, 0 ),
                                CHF_CONST_FRA1( this_gammaHallold_dir, 0 ),
                                CHF_CONST_FRA1( this_J0_dir, 0 ),
                                CHF_CONST_FRA1( this_Ne_dir, 0 ),
                                CHF_CONST_FRA1( this_Eold_dir, 0 ),
                                CHF_CONST_FRA1( this_Jold_dir, 0 ),
                                CHF_FRA1( this_gammaHall_dir, 0 ),
                                CHF_FRA1( this_hallDriftSF_dir, 0 ),
                                CHF_FRA1( this_EonEdges_dir, 0 ),
                                CHF_FRA1( this_JonEdges_dir, 0 ) );

           
   }

}

void ExtendedMhdOp::computeNTFfaceArea( LevelData<FluxBox>&    a_Flux_norm,
                                  const LevelData<FArrayBox>&  a_Flux_mapped_cc,
                                  const LevelData<FArrayBox>&  a_Jf_cc,
                                  const LevelData<FArrayBox>&  a_Cspeed_cc,
                                  const LevelData<FluxBox>&    a_CspeedR_norm,
                                  const LevelData<FluxBox>&    a_CspeedL_norm ) const
{
   CH_TIME("ExtendedMhdOp::computeNTFfaceArea()");
   CH_assert(a_Flux_norm.nComp() == 1);
   CH_assert(a_Flux_mapped_cc.nComp() == SpaceDim);
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      // compute Lax splitting on mapped flux and interpolate FluxL and FluxR to face norms
      //
      //   get left and right going flux at cell-center
      //   fluxR = 0.5*(flux + Cspeed*fun),
      //   fluxL = 0.5*(flux - Cspeed*fun),
      //   Cspeed = abs(max(eigenValue of Flux Jacobian))

      const Box& cell_box( a_Jf_cc[dit].box() );
      FArrayBox FluxL_cc(cell_box, SpaceDim);
      FArrayBox FluxR_cc(cell_box, SpaceDim);
      const FArrayBox& Flux_on_patch = a_Flux_mapped_cc[dit];
      const FArrayBox& Cspeed_on_patch = a_Cspeed_cc[dit];
      const FArrayBox& Jf_on_patch = a_Jf_cc[dit];

      FORT_LAX_FLUX_SPLITTING( CHF_BOX(cell_box), 
                               CHF_FRA(FluxL_cc),
                               CHF_FRA(FluxR_cc),
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
         const FArrayBox& thisSmooth = m_NTrowSq[dit];
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
      }
   }
 
}

void ExtendedMhdOp::computeIdealEatEdges( LevelData<EdgeDataBox>&  a_Edge_covar,
                                    const LevelData<FArrayBox>&    a_V_phys_cc,
                                    const LevelData<FArrayBox>&    a_B_contr_cc,
                                    const LevelData<FArrayBox>&    a_Cspeed_cc ) const
{
   CH_TIME("ExtendedMhdOp::computeIdealEatEdges()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   CH_assert(a_Edge_covar.nComp() == 1);
   CH_assert(a_V_phys_cc.nComp() == SpaceDim);
   CH_assert(a_Cspeed_cc.nComp() == SpaceDim);

   // compute covar E0=-VxB at cell center
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
         E0_covar_cc[dit].mult(JaB_contr_cc[dit],cell_box,cell_box,0,n,1); // now is covar E0
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
      const Box& flux_box( a_Edge_covar[dit].box() ); // this box has ghosts
      FluxBox FluxL_norm( flux_box, 1 );
      FluxBox FluxR_norm( flux_box, 1 );
      FluxBox VB_norm( flux_box, 1 );
      if(m_advScheme=="weno5") {
         const FArrayBox& thisSmooth = m_NTrowSq[dit];
         SpaceUtils::interpToFacesWENO(FluxL_norm, FluxL_cc, CspeedL_norm, thisSmooth, grid_box, m_advScheme);
         SpaceUtils::interpToFacesWENO(FluxR_norm, FluxR_cc, CspeedR_norm, thisSmooth, grid_box, m_advScheme);
      } 
      else {
         SpaceUtils::upWindToFaces(FluxL_norm, FluxL_cc, CspeedL_norm, grid_box, m_advScheme);
         SpaceUtils::upWindToFaces(FluxR_norm, FluxR_cc, CspeedR_norm, grid_box, m_advScheme);
         //SpaceUtils::upWindToFaces(FluxL_norm, FluxL_cc, CspeedL_norm, grid_box, "c2");
         //SpaceUtils::upWindToFaces(FluxR_norm, FluxR_cc, CspeedR_norm, grid_box, "c2");
      }

      // compute total flux on faces by adding left and right together
      VB_norm.copy(FluxR_norm,VB_norm.box());
      VB_norm += FluxL_norm;
      a_Edge_covar[dit][0].copy(VB_norm[1]);
      a_Edge_covar[dit][1].copy(VB_norm[0]);
      a_Edge_covar[dit][1].negate();
   } 

   // rebase all dummys
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_oneComp[dit].setVal(1./0.);
      m_dummyFArray_spaceDim[dit].setVal(1./0.);
   }

}

void ExtendedMhdOp::computeMappedCovariantCurl( LevelData<EdgeDataBox>&  a_JaCurlF_covar,
                                          const LevelData<FArrayBox>&    a_F_covar ) const
{
   CH_TIME("ExtendedMhdOp::computeMappedCovariantCurl()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   //CH_assert(SpaceDim==2);

   // This function takes covariant F on cells and returns the covariant curl at cell edges
   // multiplied by the Jacobian (aka mapped): Ja*curl(F)_covar

   // compute each component of the covariant gradient on each face (dF/dxi) 
   LevelData<FluxBox>& gradF_covar = m_dummyFlux_spaceDim;
   m_geometry.computeMappedGradient(a_F_covar, gradF_covar, 2);   
   //SpaceUtils::exchangeFluxBox(gradF_covar);

   // compute Ja times the covariant curl on cell edges
   m_geometry.mappedCovarCurlofVirtComp(a_JaCurlF_covar, gradF_covar); // Jacobian * curl(F)\cdot g_l
   //SpaceUtils::exchangeEdgeDataBox(a_JaCurlF_covar);
      
   // reset dummy container to infinity
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFlux_spaceDim[dit].setVal(1./0.);
   }
   
}     

void ExtendedMhdOp::computeDivJStress( LevelData<EdgeDataBox>&  a_divJStress_covar,
                                 const LevelData<EdgeDataBox>&  a_J_covar,
                                 const LevelData<EdgeDataBox>&  a_Ue_covar,
                                 const LevelData<EdgeDataBox>&  a_Ne_ce,
                                 const LevelData<FArrayBox>&    a_U_phys_cc,
                                 const LevelData<FArrayBox>&    a_J_covar_cc,
                                 const LevelData<FArrayBox>&    a_Ne_cc ) const
{
   CH_TIME("ExtendedMhdOp::computeDivJStress()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
      
   // compute covariant div of stress tensor for Ohm's law on cell edges
   //
   // divJStress = div(UJ+JUe) = div(UJ+JU-J*lambdai*J/Ne)
   //            = (divU)J + UdotGrad(J) + JdotGrad(U) + (divJ)Ue - JdotGrad(lambdai*J/Ne)
   //            = (divU)J + Grad(UdotJ) - U x (Grad x J)
   //            + (divJ)Ue - J x (Grad x Ue) - lambdai*Ne*Grad(J^2/Ne^2/2) 

   // initialize a_divJStress_covar to zeros
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_divJStress_covar[dit].setVal( 0.0 );
   }

   /////////////////////////////////////
   //
   // add -lambdai*Grad(J^2/2/Ne^2) term 
   // (method here ensures curl of gradient will be zero)
   //
   /////////////////////////////////////

   // convert J_covar_cc to phys
   LevelData<FArrayBox>& J_phys_cc = m_dummyFArray_spaceDim;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      J_phys_cc[dit].copy( a_J_covar_cc[dit] );
   } 
   m_geometry.convertPhysToCovar(J_phys_cc,1);
   
   // calc phi = -lambdai*(Jphys/Ne)^2/2.0 at cell center 
   LevelData<FArrayBox>& phi = m_dummyFArray_oneComp;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      SpaceUtils::localVectorNorm(phi[dit], J_phys_cc[dit]);
      phi[dit].mult( phi[dit] );
      phi[dit].divide( a_Ne_cc[dit] );
      phi[dit].divide( a_Ne_cc[dit] );
      phi[dit].mult( -m_lambdai0/2.0 );
   } 

   // compute covar gradient at cell center and interpolate to edges
   LevelData<FArrayBox>& gradphi_covar_cc = m_dummyFArray_spaceDim;
   LevelData<EdgeDataBox>& gradphi_covar_ce = m_dummyEdge_oneComp;
   m_geometry.computeMappedPoloidalGradientWithGhosts( phi, gradphi_covar_cc, 2);
   SpaceUtils::interpCellToEdges(gradphi_covar_ce, gradphi_covar_cc, gradphi_covar_ce, "c2");  

   // multiply by density and add result to a_divJStress
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         gradphi_covar_ce[dit][dir].mult( a_Ne_ce[dit][dir] );
         a_divJStress_covar[dit][dir] += gradphi_covar_ce[dit][dir];
      }
   }

   // rebase all dummys
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_oneComp[dit].setVal(1./0.);
      m_dummyFArray_spaceDim[dit].setVal(1./0.);
      m_dummyEdge_oneComp[dit].setVal(1./0.);
   }
   
   /////////////////////////////////////
   //
   // add -J x (Grad x Ue) term
   // ( covar: (J x F)*g_i = Ja*(J^j*F^k - J^k*F^j)
   //                      = Ja*(|g^j|^2*J_j*F^k - |g^k|^2*J_k*F^j) )
   //   contra: F^k = (Grad x Ue)*g^k
   //
   /////////////////////////////////////

   // the following takes covar Ue and returns Ja*curl(Ue)_contra
   LevelData<FArrayBox>& JaCurlUe_contra_cc = m_dummyFArray_oneComp;
   m_geometry.mappedGridCurlFromEdgeTans(a_Ue_covar,2,JaCurlUe_contra_cc);
  
   // upwind Ja*curl(Ue)_contra from cells to faces using Vhall
   // then copy to edgebox
   LevelData<FluxBox>& JaCurlUe_contra_cf = m_dummyFlux_oneComp;
   SpaceUtils::interpToFaces( JaCurlUe_contra_cf,
                              JaCurlUe_contra_cc,
                              JaCurlUe_contra_cc,
                              m_Vhall_cf,
                              m_advSchemeJStress );
   //
   LevelData<EdgeDataBox>& JxCurlUe_covar_ce = m_dummyEdge_oneComp;  
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      JxCurlUe_covar_ce[dit][0].copy( JaCurlUe_contra_cf[dit][1] );
      JxCurlUe_covar_ce[dit][1].copy( JaCurlUe_contra_cf[dit][0] );
   }
 
   // convert J_covar_cc to contra
   LevelData<FArrayBox>& J_contra_cc = m_dummyFArray_spaceDim;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      J_contra_cc[dit].copy( a_J_covar_cc[dit] );
   } 
   m_geometry.convertPhysToCovar(J_contra_cc,1);
   m_geometry.convertPhysToContravar(J_contra_cc,0);
 
   // interp to faces, then copy to edges
   LevelData<FluxBox>& J_contra_cf = m_dummyFlux_oneComp;
   SpaceUtils::upWindToFaces(J_contra_cf, J_contra_cc, J_contra_cf, "c2");

   // take cross product on edges
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      JxCurlUe_covar_ce[dit][0].mult( J_contra_cf[dit][1] );
      JxCurlUe_covar_ce[dit][1].mult( J_contra_cf[dit][0] );
      JxCurlUe_covar_ce[dit][0].negate();
   }

   // add -result to divJStress_covar
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         a_divJStress_covar[dit][dir] -= JxCurlUe_covar_ce[dit][dir];
      }
   }
   
   // rebase all dummys
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_oneComp[dit].setVal(1./0.);
      m_dummyEdge_oneComp[dit].setVal(1./0.);
   }
   
   /////////////////////////////////////
   //
   // add -U x (Grad x J) term
   //
   /////////////////////////////////////

}



void ExtendedMhdOp::computeHallVelocity( LevelData<FArrayBox>&    a_Vhall_cc,
                                   const LevelData<EdgeDataBox>&  a_Jcovar_ce,
                                   const LevelData<EdgeDataBox>&  a_Ne_ce ) const
{
   CH_TIME("ExtendedMhdOp::computeHallVelocity()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // compute physical Vhall = -lambdai*Jphys/Ne at cell center
   LevelData<EdgeDataBox>& Vhall_ce = m_dummyEdge_oneComp;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      Vhall_ce[dit].copy( a_Jcovar_ce[dit] );
      for (int dir=0; dir<SpaceDim; ++dir) {
         Vhall_ce[dit][dir].divide( a_Ne_ce[dit][dir],0,0,1 );
         Vhall_ce[dit][dir].mult( -m_lambdai0 );
      }
   }
   SpaceUtils::interpEdgesToCell(a_Vhall_cc,Vhall_ce,"c2");
   m_geometry.convertPhysToCovar(a_Vhall_cc,1);

   // rebase dummy containers used
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyEdge_oneComp[dit].setVal(1./0.);
   }

}

void ExtendedMhdOp::computeHallEatEdges( LevelData<EdgeDataBox>&  a_etaHall_ce,
                                         LevelData<EdgeDataBox>&  a_etaHall_op,
                                         LevelData<EdgeDataBox>&  a_Ehall_covar,
                                         LevelData<EdgeDataBox>&  a_Ehall_op,
                                   const LevelData<EdgeDataBox>&  a_J_covar_ce,
                                   const LevelData<FArrayBox>&    a_B_contr_cc,
                                   const int                      a_NoUpdateBonFaces )
{
   CH_TIME("ExtendedMhdOp::computeHallEatEdges()");
   CH_assert(a_Ehall_covar.nComp() == 1);
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // Hall electric field: E = lambdi0/Ne*J x B = -Vhall x B = J x etaHall
   // covar: E_i = lambdai0/Ne*Ja*(J^j*B^k - J^k*B^j)
   //            = lambdai0/Ne*Ja*(|g^j|^2*J_j*B^k - |g^k|^2*J_k*B^j)

   // compute contravariant Vhall = -lambdai*J/Ne at cell center
   // and interp to cell faces (yes, faces) using c2
   
   LevelData<FArrayBox>& Vhall_cc = m_dummyFArray_spaceDim;
   SpaceUtils::interpEdgesToCell(Vhall_cc,a_J_covar_ce,"c2");
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      Vhall_cc[dit].negate();
   }
   SpaceUtils::upWindToFaces(m_Vhall_cf, Vhall_cc, m_Vhall_cf, "c2");

   if(!a_NoUpdateBonFaces) {

      // compute Ehall_covar on cell faces by upwiding Ja*B_contravar
      LevelData<FArrayBox>& VhallBy_cc = m_dummyFArray_spaceDim;
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int n=0; n<SpaceDim; ++n) {
            VhallBy_cc[dit].copy( a_B_contr_cc[dit],0,n,1 );
            VhallBy_cc[dit].mult( m_Jacobian[dit],0,n,1 );
         }
      }
      
      LevelData<FluxBox>& etaHall_cf = m_dummyFlux_oneComp;
      SpaceUtils::interpToFaces( etaHall_cf, 
                                 VhallBy_cc, 
                                 VhallBy_cc, 
                                 m_Vhall_cf, 
                                 m_advSchemeHall );

      // set etaHall on faces, copy to edgebox for etaHall_op, 
      // and then convert (multiply by appropriate |g|^2)
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            const Box& thisBox = etaHall_cf[dit][dir].box();
            etaHall_cf[dit][dir].divide(m_eleDen_cf[dit][dir],thisBox,0,0,1);
            etaHall_cf[dit][dir].mult( m_lambdai0 );
         }
         a_etaHall_op[dit][0].copy(etaHall_cf[dit][1]);
         a_etaHall_op[dit][1].copy(etaHall_cf[dit][0]);
      }  
      m_geometry.convertContravarToCovar(a_etaHall_op,1);

      // convert etaHall_cf and then copy to edgebox
      m_geometry.convertContravarToCovar(etaHall_cf,1);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_etaHall_ce[dit][0].copy(etaHall_cf[dit][1]);
         a_etaHall_ce[dit][1].copy(etaHall_cf[dit][0]);
      }

   }

   // compute covariant Ehall = -VhallxB at cell edges
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_Ehall_covar[dit][0].copy(a_etaHall_ce[dit][0]);
      a_Ehall_covar[dit][0].mult(m_Vhall_cf[dit][1]);
      a_Ehall_covar[dit][1].copy(a_etaHall_ce[dit][1]);
      a_Ehall_covar[dit][1].mult(m_Vhall_cf[dit][0]);
      a_Ehall_covar[dit][1].negate();
   }

   // compute covariant Ehall = -VhallxB at oposite cell edges
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_Ehall_op[dit][0].copy(a_J_covar_ce[dit][0]);
      a_Ehall_op[dit][0].copy(a_etaHall_op[dit][0]);
      a_Ehall_op[dit][1].mult(a_J_covar_ce[dit][1]);
      a_Ehall_op[dit][1].mult(a_etaHall_op[dit][1]);
      a_Ehall_op[dit][1].negate();
   }

   // rebase dummy containers used
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_spaceDim[dit].setVal(1./0.);
      m_dummyFlux_oneComp[dit].setVal(1./0.);
   }

}


void ExtendedMhdOp::computeEleVelOnCells( LevelData<FArrayBox>&  a_Ue_cc,
                                    const LevelData<FArrayBox>&  a_J_cc,
                                    const LevelData<FArrayBox>&  a_Ne_cc,
                                    const LevelData<FArrayBox>&  a_U_cc )
{
   CH_TIME("ExtendedMhdOp::computeEleVelOnCells()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // compute electron velocity at cell center
   // Ue = U - lambdai*J/Ne
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_Ue_cc[dit].copy( a_J_cc[dit] );
      for (int n=0; n<SpaceDim; ++n) {
         a_Ue_cc[dit].divide( a_Ne_cc[dit],0,n,1 );
      }
      a_Ue_cc[dit].mult( -m_lambdai0 );
      for (int n=0; n<SpaceDim; ++n) {
      //   a_Ue_cc[dit].plus( a_U_cc[dit],n,n,1 );
      }
   }

}


void ExtendedMhdOp::setHallDriftSpeed( const LevelData<FArrayBox>&  a_Ne,
                                       const LevelData<FArrayBox>&  a_By_covar )
{
   CH_TIME("ExtendedMhdOp::setHallDriftSpeed()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // compute contravariant Vhalldrift speed at cell center for
   // time step calculation
   //
   // Vhd = -lambdai/Ne*grad(Ne) x B = A x B
   // (contra) Vhd^i = (A_jB_k - A_kB_j)/Ja
   // (covar) grad(Ne)_j = d(Ne)/dxi_j
   
   // compute mapped grad(Ne) = d(Ne)/dxi at cell center
   m_geometry.computeMappedPoloidalGradientWithGhosts( m_eleDen_cc, m_dummyFArray_spaceDim, 2);
   
   // set VHallDrift = Li0*BcrossGrad(Ne)/Ne^2 contravariant
   for (DataIterator dit(grids); dit.ok(); ++dit) {
     
      FArrayBox& Vhd_on_patch = m_VhallDrift_cc[dit];
      const FArrayBox& gradNe_on_patch = m_dummyFArray_spaceDim[dit];   
      const FArrayBox& Ne_on_patch     = m_eleDen_cc[dit];   
      const FArrayBox& By_on_patch     = m_By_covar_cc[dit]; 
      const FArrayBox& Jac_on_patch    = m_Jacobian[dit]; 
      const FArrayBox& hdSF_on_patch   = m_hallDriftScaleFactor[dit];

      const Box& thisbox( grids[dit] ); // no ghosts
      FORT_EVAL_HALL_DRIFT_SPEED( CHF_BOX(thisbox), 
                                  CHF_CONST_FRA(gradNe_on_patch),
                                  CHF_CONST_FRA1(Ne_on_patch,0),
                                  CHF_CONST_FRA1(By_on_patch,0),
                                  CHF_CONST_FRA1(Jac_on_patch,0),
                                  CHF_CONST_FRA1(hdSF_on_patch,0),
                                  CHF_CONST_REAL(m_lambdai0),
                                  CHF_FRA(Vhd_on_patch) );
      Vhd_on_patch.abs( ); // only used for time-step
   }
   
}

void ExtendedMhdOp::computeElectronVelOnEdges( LevelData<EdgeDataBox>&  a_Ve_covar,
                                         const LevelData<EdgeDataBox>&  a_J_covar,
                                         const LevelData<EdgeDataBox>&  a_Ne_ce,
                                         const LevelData<FArrayBox>&    a_U_phys_cc ) const
{
   CH_TIME("ExtendedMhdOp::computeElectronVelOnEdges()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
  
   // compute contravariant electron velocity at cell edges
   // Ve = U - lambdai*J/Ne
   // 
   LevelData<FArrayBox>& U_covar_cc = m_dummyFArray_spaceDim;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      U_covar_cc[dit].copy(a_U_phys_cc[dit]);
   }
   m_geometry.convertPhysToCovar(U_covar_cc,0);
   
   LevelData<EdgeDataBox>& U_covar_ce = m_dummyEdge_oneComp;
   SpaceUtils::interpCellToEdges(U_covar_ce,U_covar_cc,U_covar_ce,"c2");
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_Ve_covar[dit].copy( a_J_covar[dit] );
      for(int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_Ve(a_Ve_covar[dit][dir]);
         this_Ve.mult(-m_lambdai0);
         this_Ve.divide(a_Ne_ce[dit][dir]);
         this_Ve.plus(U_covar_ce[dit][dir]);
      }
   }
   SpaceUtils::exchangeEdgeDataBox(a_Ve_covar); // always do this... don't ask why...
  
}

void ExtendedMhdOp::computeElectronVelForVisc( LevelData<EdgeDataBox>&  a_Ve_ce,
                                         const int                      a_method,
                                         const LevelData<EdgeDataBox>&  a_J_ce,
                                         const LevelData<EdgeDataBox>&  a_N_ce,
                                         const LevelData<FArrayBox>&    a_N_cc )
{
   CH_TIME("ExtendedMhdOp::computeElectronVelForVisc()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
  
   // This function computes the effective electron velocity 
   // to be used for viscosity in Ohm's law. The effective diffusion
   // coefficients needed for the relaxation scheme depends on the
   // form used for Ve and are thus also set here
   //
 
   // compute contravariant electron velocity at cell edges
   // Ve = U - lambdai*J/Ne (method = 2) 
   // Ve = - lambdai*J/Ne   (method = 1)
   // Ve = - lambdai*J     (method = 0) (should use this for etaVis_ele=const)
   // 
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_Ve_ce[dit].copy( a_J_ce[dit] );
      for(int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_Ve(a_Ve_ce[dit][dir]);
         this_Ve.mult(-m_lambdai0);
         if(a_method>0) this_Ve.divide(a_N_ce[dit][dir]);
         //if(a_method==2) this_Ve.plus(-m_velocity_ce[dit]);
      }
   }
   SpaceUtils::exchangeEdgeDataBox(a_Ve_ce); // always do this... don't ask why...
  
   
   // set the effective diffusion coefficients for relaxation scheme
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_Pie_diagD_cn[dit].setVal( 2.0*m_etaVis_ele/m_mM/m_meRel );
      m_Pie_xzD_cc[dit].setVal( 2.0*m_etaVis_ele/m_mM/m_meRel );
      m_Pie_yyD_ce[dit].setVal( 2.0*m_etaVis_ele/m_mM/m_meRel );
      if(a_method>0) m_Pie_xzD_cc[dit].divide(a_N_cc[dit]);
   }
   if(a_method>0) {
      SpaceUtils::interpCellToEdges(m_Pie_yyD_ce,m_Pie_xzD_cc,m_Pie_yyD_ce,"c2");
      SpaceUtils::interpCellsToNodes(m_Pie_diagD_cn,m_Pie_xzD_cc,"c2");
   }   

}


void ExtendedMhdOp::computeStrainTensorStag( LevelData<NodeFArrayBox>&  a_Pie_diag_cn,
                                             LevelData<FArrayBox>&      a_Pie_xz_cc,
                                             LevelData<EdgeDataBox>&    a_Pie_yy_ce,
                                       const LevelData<EdgeDataBox>&    a_Ve_ce )
{
   CH_TIME("ExtendedMhdOp::computeStrainTensorStag()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // calculate the gradient of Ve at nodes and cells:
   // m_gradVepara_cn.comp0 = dVex/dx 
   // m_gradVepara_cn.comp1 = dVez/dz
   // m_gradVeperp_cc.comp0 = dVex/dz
   // m_gradVeperp_cc.comp1 = dVez/dx
   //
   m_geometry.ParaGradAtNodesFromEdges( m_gradVepara_cn, a_Ve_ce );
   m_geometry.PerpGradAtCellsFromEdges( m_gradVeperp_cc, a_Ve_ce );
   m_gradVepara_cn.exchange(); 
   m_gradVeperp_cc.exchange(); 
   SpaceUtils::interpNodesToCells(m_gradVepara_cc, m_gradVepara_cn, "c2");

   // compute the divergence of Ve at cell nodes and interpolae to cell edges
   //
   m_geometry.physicalDivergenceFromEdgeTans(m_divVe_cn, a_Ve_ce);
   m_divVe_cn.exchange(); 
   //for (DataIterator dit(grids); dit.ok(); ++dit) {
   //   FArrayBox& thisdivVe(m_divVe_cn[dit].getFab());
   //   thisdivVe.setVal(0.0);
   //}
   SpaceUtils::interpNodesToEdges(m_divVe_ce, m_divVe_cn, "c2");

   // set components of stress tensor on cell edges
   //
   // Pie_diag_cn.comp0 = 2*dVex/dx - 2/3*divVe
   // Pie_diag_cn.comp1 = 2*dVez/dz - 2/3*divVe
   // Pie_yy_ce = 2*Vex/x   - 2/3*divVe
   // Pie_xz_cc = dVez/dx + dVex/dz
   //
   const LevelData<EdgeDataBox>& Xphys_ce = m_geometry.getEdgeCenteredRealCoords();
   for (DataIterator dit(grids); dit.ok(); ++dit) {
     
      FArrayBox& Pie_diag_on_patch = a_Pie_diag_cn[dit].getFab();
      const FArrayBox& gradVe_on_patch = m_gradVepara_cn[dit].getFab();   
      const FArrayBox& divVe_on_patch  = m_divVe_cn[dit].getFab();   
      const Box& nodebox = Pie_diag_on_patch.box();
      FORT_EVAL_PI_ELE_DIAG( CHF_BOX(nodebox), 
                             CHF_CONST_FRA1(gradVe_on_patch,0),
                             CHF_CONST_FRA1(gradVe_on_patch,1),
                             CHF_CONST_FRA1(divVe_on_patch,0),
                             CHF_FRA(Pie_diag_on_patch) );
      //
      //
      //
      int dir=0;
      FArrayBox& Pie_yy_on_patch = a_Pie_yy_ce[dit][dir]; // only need it for dir=0
      const FArrayBox& Ve_on_patch = a_Ve_ce[dit][dir];   
      const FArrayBox& divVe_ce_on_patch  = m_divVe_ce[dit][dir];   
      const FArrayBox& X_on_patch = Xphys_ce[dit][dir];
      const Box& edgebox = Pie_yy_on_patch.box();
      FORT_EVAL_PI_ELE_YY( CHF_BOX(edgebox), 
                           CHF_CONST_FRA1(Ve_on_patch,dir),
                           CHF_CONST_FRA1(X_on_patch,0),
                           CHF_CONST_FRA1(divVe_ce_on_patch,dir),
                           CHF_FRA1(Pie_yy_on_patch,dir) );
      //
      //
      //
      a_Pie_xz_cc[dit].copy(m_gradVeperp_cc[dit],0,0,1);
      a_Pie_xz_cc[dit].plus(m_gradVeperp_cc[dit],1,0,1);
   
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_Pie_diag_cn[dit].mult( -m_etaVis_ele );
      for (int dir=0; dir<SpaceDim; ++dir) {
         a_Pie_yy_ce[dit][dir] *=  -m_etaVis_ele ;
      }
      a_Pie_xz_cc[dit].mult( -m_etaVis_ele );
   }

   a_Pie_diag_cn.exchange();
   a_Pie_xz_cc.exchange();
   a_Pie_yy_ce.exchange();

}

void ExtendedMhdOp::computeDivOfElectronVisc( LevelData<EdgeDataBox>&    a_divPie_ce,
                                        const LevelData<NodeFArrayBox>&  a_Pie_para_cn,
                                        const LevelData<FArrayBox>&      a_Pie_xz_cc,
                                        const LevelData<EdgeDataBox>&    a_Pie_yy_ce )
{
   CH_TIME("ExtendedMhdOp::computeDivOfElectronVisc()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   m_geometry.physicalDivergenceAtEdges( a_divPie_ce, a_Pie_para_cn, a_Pie_xz_cc );
   
   // add viscosity flux source term to divegence for axisymmetric
   //
   if(m_twoDaxisymm) {
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
     
         const int thisdir = 0;      

         const FArrayBox& Xphys_onPatch  = m_Xphys[dit];
         const FArrayBox& Pyy_onPatch    = a_Pie_yy_ce[dit][thisdir];
               FArrayBox& divPie_onPatch = a_divPie_ce[dit][thisdir];
      
         Box edgebox( grids[dit] ); // no ghosts
         edgebox.surroundingNodes(  );
         edgebox.enclosedCells( thisdir );

         FORT_ADD_SOURCE_TO_DIVPIE( CHF_BOX(edgebox),
                                    CHF_CONST_FRA1(Xphys_onPatch,0),
                                    CHF_CONST_FRA1(Pyy_onPatch,0),
                                    CHF_FRA1(divPie_onPatch,0) );
      } 
      
   }
   SpaceUtils::exchangeEdgeDataBox(a_divPie_ce);
 
}

void ExtendedMhdOp::addDivPieToGradPele( LevelData<EdgeDataBox>&  a_gradPe_ce,
                                   const LevelData<EdgeDataBox>&  a_divPie_ce )
{
   CH_TIME("ExtendedMhdOp::addDivPieToGradPele()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         const FArrayBox& thisDivPie = a_divPie_ce[dit][dir];
               FArrayBox& thisGradPe = a_gradPe_ce[dit][dir];
         thisGradPe.plus( thisDivPie );
      } 
   }
  
}

void ExtendedMhdOp::advanceElectronViscFluxes( LevelData<NodeFArrayBox>&  a_Pie_diag_cn,
                                               LevelData<FArrayBox>&      a_Pie_xz_cc,
                                               LevelData<EdgeDataBox>&    a_Pie_yy_ce,
                                         const LevelData<EdgeDataBox>&    a_Ve_ce,
                                         const int                        a_stage,
                                         const double                     a_time,
                                         const double                     a_dt )
{
   CH_TIME("ExtendedMhdOp::advanceElectronViscFluxes()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // compute the relaxed viscosity fluxes on cell faces
   //
   computeStrainTensorStag( m_Pie_diag_cn0, m_Pie_xz_cc0, m_Pie_yy_ce0, a_Ve_ce);  
 
   // set initial and old values if t=0
   //
   if(a_time==0) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_Pie_diag_cn[dit].copy( m_Pie_diag_cn0[dit] );
         a_Pie_xz_cc[dit].copy( m_Pie_xz_cc0[dit],0,0,1 );
         a_Pie_yy_ce[dit][0].copy( m_Pie_yy_ce0[dit][0] );
         m_Pie_diag_cn_old[dit].copy( m_Pie_diag_cn0[dit] );
         m_Pie_xz_cc_old[dit].copy( m_Pie_xz_cc0[dit],0,0,1 );
         m_Pie_yy_ce_old[dit][0].copy( m_Pie_yy_ce0[dit][0] );
      }
   }
      
   // update viscosity fluxes using relaxation scheme
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {

      const FArrayBox& this_D_cn      = m_Pie_diagD_cn[dit].getFab();
      const FArrayBox& this_Pieold_cn = m_Pie_diag_cn_old[dit].getFab();
      const FArrayBox& this_Pie0_cn   = m_Pie_diag_cn0[dit].getFab();
            FArrayBox& this_Pie_cn    = m_Pie_diag_cn[dit].getFab();

      Box nodebox( grids[dit] );   // no ghost 
      nodebox.surroundingNodes( ); // grow hi end by 1 in all dirs
      for (int n=0; n<this_Pie_cn.nComp(); n++) {
         FORT_ADVANCE_USING_RELAX_NEW( CHF_BOX(nodebox),
                                   CHF_CONST_REAL( a_dt ),
                                   CHF_CONST_REAL( m_delta0 ),
                                   CHF_CONST_FRA1( this_D_cn,0 ),
                                   CHF_CONST_FRA1( this_Pieold_cn,n ),
                                   CHF_CONST_FRA1( this_Pie0_cn,n ),
                                   CHF_FRA1( this_Pie_cn,n ) );
      }
      //
      //
      //
      const FArrayBox& this_D_cc      = m_Pie_xzD_cc[dit];
      const FArrayBox& this_Pieold_cc = m_Pie_xz_cc_old[dit];
      const FArrayBox& this_Pie0_cc   = m_Pie_xz_cc0[dit];
            FArrayBox& this_Pie_cc    = m_Pie_xz_cc[dit];

      Box cellbox( grids[dit] );   // no ghost 
      cellbox.grow( 1 );           // add 1 layer of ghost cells
      for (int n=0; n<this_Pie_cc.nComp(); n++) {
         FORT_ADVANCE_USING_RELAX_NEW( CHF_BOX(nodebox),
                                   CHF_CONST_REAL( a_dt ),
                                   CHF_CONST_REAL( m_delta0 ),
                                   CHF_CONST_FRA1( this_D_cc,0 ),
                                   CHF_CONST_FRA1( this_Pieold_cc,n ),
                                   CHF_CONST_FRA1( this_Pie0_cc,n ),
                                   CHF_FRA1( this_Pie_cc,n ) );
      }
      m_Pie_xz_cc.exchange();
         
      if(a_stage==0) { // update old values
         m_Pie_diag_cn_old[dit].copy( m_Pie_diag_cn[dit] );
         m_Pie_xz_cc_old[dit].copy( m_Pie_xz_cc[dit] );
      }
      
      // advance source term for axisymmetric
      //
      if(m_twoDaxisymm) {
         const int thisdir = 0;
         const FArrayBox& this_D_ce      = m_Pie_yyD_ce[dit][thisdir];
         const FArrayBox& this_Pieold_ce = m_Pie_yy_ce_old[dit][thisdir];
         const FArrayBox& this_Pie0_ce   = m_Pie_yy_ce0[dit][thisdir];
               FArrayBox& this_Pie_ce    = m_Pie_yy_ce[dit][thisdir];

         const Box& cellbox = this_Pie_ce.box();
         FORT_ADVANCE_USING_RELAX_NEW( CHF_BOX(cellbox),
                                   CHF_CONST_REAL( a_dt ),
                                   CHF_CONST_REAL( m_delta0 ),
                                   CHF_CONST_FRA1( this_D_ce,0 ),
                                   CHF_CONST_FRA1( this_Pieold_ce,0 ),
                                   CHF_CONST_FRA1( this_Pie0_ce,0 ),
                                   CHF_FRA1( this_Pie_ce,0 ) );
              
         if(a_stage==0) { // update old values
            m_Pie_yy_ce_old[dit][thisdir].copy( m_Pie_yy_ce[dit][thisdir] );
         }
      }
   
   }

}

void ExtendedMhdOp::initializeWithBC( FluidSpecies&  a_species_comp,
                                      FluidSpecies&  a_species_phys,
                                const double         a_time )
{
   CH_TIME("ExtendedMhdOp::initializeWithBC()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   int this_ebc;

   if(m_initializeJfromAmperesLaw) {
    
      CH_assert(a_species_phys.m_evolve_currentDensity==1);

      // compute covar virtual magnetic field
      const LevelData<FArrayBox>& phys_By( a_species_phys.cell_var("magneticField_virtual") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& thisbox = m_By_covar_cc[dit].box();
         m_By_covar_cc[dit].copy( phys_By[dit], thisbox); // phys
         if(m_twoDaxisymm) {
            m_By_covar_cc[dit].mult(m_g_y[dit]); // covar
         }
      }
      
      // compute covariant J0=curl(By) at cell-edges and apply BC
      computeMappedCovariantCurl( m_J0_ce, m_By_covar_cc ); // Jacobian * curl(By)\cdot g_l (Ja*covar curl)
      m_geometry.divideJonEdges(m_J0_ce);                   // curl(By)\cdot g_l
      int this_cbc = m_magneticField_virtual_cbc;
      if(m_twoDaxisymm) m_fluid_bc.at(this_cbc)->setOnAxisCurlBC( m_J0_ce, phys_By, a_species_phys );
      SpaceUtils::exchangeEdgeDataBox(m_J0_ce);
      
      this_ebc = m_currentDensity_ebc;
      m_fluid_edge_var_bc.at(this_ebc)->applyEdgeBC( a_species_phys, m_J0_ce, a_time );

      // copy covar J0 to phys_Jcovar (has ghosts) and comp_Jcovar (no ghosts)
      LevelData<EdgeDataBox>& phys_Jcovar( a_species_phys.edge_var("currentDensity") );
      LevelData<EdgeDataBox>& comp_Jcovar( a_species_comp.edge_var("currentDensity") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
          phys_Jcovar[dit].copy( m_J0_ce[dit],0,0,phys_Jcovar.nComp() );
          comp_Jcovar[dit].copy( m_J0_ce[dit],0,0,comp_Jcovar.nComp() );
      }
         
   }
     
   if(m_initializeEfromOhmsLaw) {
    
      CH_assert(a_species_phys.m_evolve_electricField==1);
    
      const LevelData<FArrayBox>& phys_By( a_species_phys.cell_var("magneticField_virtual") );
      const LevelData<EdgeDataBox>& phys_Jcovar( a_species_phys.edge_var("currentDensity") );

      // update cell-center member variables
      setCellCenterValues( a_species_phys );
      SpaceUtils::upWindToFaces(m_eleDen_cf,m_eleDen_cc,m_eleDen_cf,"c2");
      updateCollisionTerms( a_species_phys, m_pressure_ele, m_pressure_ion );
      SpaceUtils::interpCellToEdges(m_etaRes_ce, m_etaRes, m_etaRes_ce, "c2");
      SpaceUtils::exchangeEdgeDataBox(m_etaRes_ce);

      // compute covariant boltzmann electric field, Eboltz=-lambdai0/Ne*grad(Pe), on edges
      m_geometry.computeMappedPoloidalGradientWithGhosts( m_pressure_ele, m_gradPe_covar_cc, 2);
      SpaceUtils::interpCellToEdges(m_gradPe_covar_ce, m_gradPe_covar_cc, m_gradPe_covar_ce, "c2");
      SpaceUtils::interpCellToEdges(m_eleDen_ce, m_eleDen_cc, m_eleDen_ce,"c2");
      LevelData<EdgeDataBox> Eboltz;
      Eboltz.define(grids, 1, phys_Jcovar.ghostVect());
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         Eboltz[dit].copy( m_gradPe_covar_ce[dit],0,0,1 );
         for (int dir=0; dir<SpaceDim; ++dir) {
            Eboltz[dit][dir].divide( m_eleDen_ce[dit][dir],0,0,1 );
            Eboltz[dit][dir].mult( -m_lambdai0 );
         }
      }
   
      // compute contravariant virtual magnetic field
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_By_contr_cc[dit].copy( phys_By[dit] ); // phys
         if(m_twoDaxisymm) {
            m_By_contr_cc[dit].divide(m_g_y[dit]); // contravar
         }
      }
      
      // compute covariant ideal electric field, E0=-VxB, at cell edges
      setMappedCspeed( a_species_phys.m_gamma, 1, 1 );
      computeIdealEatEdges(m_E0_ce,m_velocity,m_By_contr_cc,m_Cspeed_cc);
      this_ebc = m_electricField_ebc;
      m_fluid_edge_var_bc.at(this_ebc)->setEdgeBC( a_species_phys, m_E0_ce, m_EdgeBC_zeros, a_time );
      SpaceUtils::exchangeEdgeDataBox(m_E0_ce);
 
      // compute covariant hall electric field, Ehall=lambdai0/N*JxB, at cell edges
      if(m_modelDriftTerms) {
         computeHallEatEdges( m_etaHall_ce, m_etaHall_op, m_Ehall_ce, m_Ehall_op, 
                              phys_Jcovar, m_By_contr_cc, 0 );
        // this_ebc = m_electricField_ebc;
        // m_fluid_edge_var_bc.at(this_ebc)->setEdgeBC( a_species_phys, m_Ehall_ce, m_EdgeBC_zeros, a_time );
      }

      // compute covariant div of stress tensor for Ohm's law on cell edges
      if(m_include_divJStress) {
         SpaceUtils::interpEdgesToCell(m_Jcc,phys_Jcovar,"c2"); // covar at cell center
         computeElectronVelOnEdges( m_Ve_ce, phys_Jcovar, m_eleDen_ce, m_velocity );
         computeDivJStress( m_divJStress_covar, phys_Jcovar, m_Ve_ce, m_eleDen_ce,
                            m_velocity, m_Jcc, m_eleDen_cc );
      }
      LevelData<EdgeDataBox>& EdivJStress = m_dummyEdge_oneComp;
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         EdivJStress[dit].copy( m_divJStress_covar[dit],0,0,1 );
         for (int dir=0; dir<SpaceDim; ++dir) {
            EdivJStress[dit][dir].divide( m_eleDen_ce[dit][dir],0,0,1 );
            EdivJStress[dit][dir].mult( m_epsilon0 );
         }
      }


      // copy covar E = Eideal + Ehall + Eboltz + Eres to phys_Ecovar
      LevelData<EdgeDataBox>& phys_Ecovar( a_species_phys.edge_var("electricField") );
      LevelData<EdgeDataBox>& comp_Ecovar( a_species_comp.edge_var("electricField") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         phys_Ecovar[dit].copy( m_E0_ce[dit],0,0,1 );
         if(m_modelDriftTerms) {
            for (int dir=0; dir<SpaceDim; ++dir) {
               phys_Ecovar[dit][dir].plus( m_Ehall_ce[dit][dir],0,0,1 );
               phys_Ecovar[dit][dir].plus( Eboltz[dit][dir],0,0,1 );
               phys_Ecovar[dit][dir].plus( EdivJStress[dit][dir],0,0,1 );
            }
         }
         //phys_Ecovar[dit][dir].plus( Eres_ce[dit][dir],0,0,1 );
      }
      
      // exchange, fill ghost cells, and copy to comp state, which is the same as
      // the phys state for edge vars with expection that phys has ghost cells
      SpaceUtils::exchangeEdgeDataBox(phys_Ecovar); 
      this_ebc = m_electricField_ebc;
      m_fluid_edge_var_bc.at(this_ebc)->applyEdgeBC( a_species_phys, phys_Ecovar, a_time );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         SpaceUtils::copyEdgeDataBox(comp_Ecovar[dit],phys_Ecovar[dit]);
      }
         

      // compute charge density at nodes from electric field
      // and interpolate to cell center for diagnostic
      m_geometry.physicalDivergenceFromEdgeTans(m_divE, comp_Ecovar);
      SpaceUtils::interpNodesToEdges(m_deltaN_ce, m_divE, "c2");
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; dir++) {
            FArrayBox& thisDeltaN( m_deltaN_ce[dit][dir] ); 
            thisDeltaN.mult( m_delta0*m_lambdai0/m_epRel );
         }
      }
      SpaceUtils::interpEdgesToCell(m_deltaN_cc, m_deltaN_ce, "c2");

      // rebase dummys
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_dummyEdge_oneComp[dit].setVal(1./0.);
      }
   
   }

}

void ExtendedMhdOp::fillGhostCells( FluidSpecies&  a_species_phys,
                              const double   a_time )
{
   CH_TIME("ExtendedMhdOp::fillGhostCells()");

   computeSelfConsistentBoundaryMomentum( a_species_phys );

   for (int n=0; n<a_species_phys.num_cell_vars(); ++n) {
      //LevelData<FArrayBox>& fld( a_species_phys.cell_var(n) );
      //m_geometry.fillInternalGhosts( fld );
      const std::string& this_var_name = a_species_phys.cell_var_name(n);
      CH_assert(m_fluid_bc.at(n)->isForVariable(this_var_name));
      m_fluid_bc.at(n)->apply( a_species_phys, a_time );
   }

   for (int n=0; n<a_species_phys.num_face_vars(); ++n) {
      const std::string& this_var_name = a_species_phys.face_var_name(n);
      CH_assert(m_fluid_face_var_bc.at(n)->isForVariable(this_var_name));
      //LevelData<FluxBox>& faceVar( a_species_phys.face_var(n) );
      //m_fluid_face_var_bc.at(n)->applyFluxBC( a_species_phys, faceVar, a_time );
   }
   
   for (int n=0; n<a_species_phys.num_edge_vars(); ++n) {
      const std::string& this_var_name = a_species_phys.edge_var_name(n);
      CH_assert(m_fluid_edge_var_bc.at(n)->isForVariable(this_var_name));
      LevelData<EdgeDataBox>& edgeVar( a_species_phys.edge_var(n) );
      m_fluid_edge_var_bc.at(n)->applyEdgeBC( a_species_phys, edgeVar, a_time );
   }

}

void ExtendedMhdOp::computeSelfConsistentBoundaryMomentum( FluidSpecies&  a_species_phys )
{
   CH_TIME("ExtendedMhdOp::computeSelfConsistentBoundaryMomentum()");

   // this function computes the boundary value of the momentum normal to a conducting
   // surface that is self-consistent with the those for the electric field (E) and 
   // current density (J) when using extended Ohms law with hall and pressure terms
   //
   // mom_norm = lambdai0*(J - b_hat x nabla(Piee)/B)\cdot\normdir 
   // where Piee = Pe + pi_e is total pressure including viscosity
   //
   //
   // Note that solution needs to be self-consistent with the semi-implicit time-advance scheme
   // such that the boundary value of E0 is indeed such that the Ohm's law gives zero value
   // for component tangential to surface
   //

}

void ExtendedMhdOp::parseParameters( ParmParse& a_pp,
                               const double     a_Nscale,
                               const double     a_Tscale,
                               const double     a_Xscale,
                               const double     a_Mscale )
{
   a_pp.query( "updateHydroVars", m_updateHydroVars);

   a_pp.query( "useRelaxedViscosity", m_useRelaxedViscosity);
   
   a_pp.query( "useRelaxedHeatFlux", m_useRelaxedHeatFlux);
   
   a_pp.query( "useRelaxedOhmsLaw", m_useRelaxedOhmsLaw);

   if (a_pp.contains("advScheme")) {
      a_pp.get("advScheme", m_advScheme );
   }
   
   if (a_pp.contains("advSchemeHall")) {
      a_pp.get("advSchemeHall", m_advSchemeHall );
   }

   if (a_pp.contains("advSchemeJStress")) {
      a_pp.get("advSchemeJStress", m_advSchemeJStress );
   }   

   if (a_pp.contains("modelDriftTerms")) {
      a_pp.get("modelDriftTerms", m_modelDriftTerms );
   }
   
   if (a_pp.contains("initializeJfromAmperesLaw")) {
      a_pp.get("initializeJfromAmperesLaw", m_initializeJfromAmperesLaw );
   }
   
   if (a_pp.contains("initializeEfromOhmsLaw")) {
      a_pp.get("initializeEfromOhmsLaw", m_initializeEfromOhmsLaw );
   }
   
   if (a_pp.contains("fullyImplicitJ")) {
      a_pp.get("fullyImplicitJ", m_fullyImplicitJ );
   }
   
   if (a_pp.contains("iterHallMax")) {
      a_pp.get("iterHallMax", m_iterHallMax );
   }
   
   if (a_pp.contains("HallResTol")) {
      a_pp.get("HallResTol", m_HallResTol );
   }
   
   if (a_pp.contains("useJ0forJcc")) {
      a_pp.get("useJ0forJcc", m_useJ0forJcc );
   }
   
   if (a_pp.contains("etaResMin")) {
      a_pp.get("etaResMin", m_etaResMin );
   }
   
   if (a_pp.contains("etaAnomFactor")) {
      a_pp.get("etaAnomFactor", m_etaAnomFactor );
   }

   if (a_pp.contains("NvacC")) {
      a_pp.get("NvacC", m_NvacC );
   }

   if (a_pp.contains("NvacP")) {
      a_pp.get("NvacP", m_NvacP );
   }

   if (a_pp.contains("Nmin")) {
      a_pp.get("Nmin", m_Nmin );
   }
   
   if (a_pp.contains("Tmin")) {
      a_pp.get("Tmin", m_Tmin );
   }
   
   if (a_pp.contains("Tmax")) {
      a_pp.get("Tmax", m_Tmax );
   }

   if (a_pp.contains("etaMin")) {
      a_pp.get("etaMin", m_etaMin );
   }
   
   if (a_pp.contains("etaVis_ele")) {
      a_pp.get("etaVis_ele", m_etaVis_ele );
   }

   if (a_pp.contains("include_divJStress")) {
      a_pp.get("include_divJStress", m_include_divJStress );
   }

   if (a_pp.contains("viscLimiter")) {
      a_pp.get("viscLimiter", m_viscLimiter );
   }
   
   if (a_pp.contains("gyroViscLimiter")) {
      a_pp.get("gyroViscLimiter", m_gyroViscLimiter );
   }
   
   if (a_pp.contains("gyroHeatViscCoeff")) {
      a_pp.get("gyroHeatViscCoeff", m_gyroHeatViscCoeff );
   }
   
   if (a_pp.contains("tauiMin")) {
      a_pp.get("tauiMin", m_tauiMin );
   }
   
   if (a_pp.contains("taueMin")) {
      a_pp.get("taueMin", m_taueMin );
   }
   
   if (a_pp.contains("kappaMin")) {
      a_pp.get("kappaMin", m_kappaMin );
   }
   
   if (a_pp.contains("heatLimiter")) {
      a_pp.get("heatLimiter", m_heatLimiter );
      m_heatLimiter_ele = m_heatLimiter;
   }
   
   if (a_pp.contains("heatLimiter_ele")) {
      a_pp.get("heatLimiter_ele", m_heatLimiter_ele );
   }
   
   if (a_pp.contains("gyroHeatLimiter")) {
      a_pp.get("gyroHeatLimiter", m_gyroHeatLimiter );
   }
   
   if (a_pp.contains("xionFactor")) {
      a_pp.get("xionFactor", m_xionFactor );
   }
   
   if (a_pp.contains("xeleFactor")) {
      a_pp.get("xeleFactor", m_xeleFactor );
   }
   
   if (a_pp.contains("meRel")) {
      a_pp.get("meRel", m_meRel );
   }
   
   if (a_pp.contains("epRel")) {
      a_pp.get("epRel", m_epRel );
   }
      
   // get fundamental constants
   const double cvac = Constants::SPEED_OF_LIGHT;
   const double mu0 = 4.0e-7*Constants::PI;
   const double qe  = Constants::ELEMENTARY_CHARGE;
   const double me  = Constants::MASS_OF_ELECTRON;
   const double amu = 1.660539e-27;  

   m_Nscale_SI = a_Nscale; // density scale 1/m^3

   // derived scales
   const double Mn = a_Mscale*amu;               // neutral mass [kg]
   const double Mi = Mn - me;                    // ion mass [kg]
   const double Pscale = qe*a_Tscale*a_Nscale;   // species pressure [J/m^3]
   const double rhoScale = Mn*a_Nscale;          // mass density scale [kg/m^3]
   const double Vscale = sqrt(Pscale/rhoScale);  // velocity scale [m/s]
   const double Bscale = sqrt(mu0*Pscale);       // magnetic field scale [Tesla]
   const double Jscale = Bscale/a_Xscale/mu0;    // current density scale [A/m^2]
   const double Escale = Vscale*Bscale;          // electric field scale [V/m]
   const double tscale = a_Xscale/Vscale;        // time scale [s]

   m_mM = me/Mi;
   m_wpescale_Hz = 5.64e4*pow(a_Nscale/1.0e6,0.5);
   m_wpiscale_Hz = m_wpescale_Hz*pow(me/Mi,0.5);
   m_wcescale_Hz = qe*Bscale/me;
   m_wciscale_Hz = qe*Bscale/Mi;
   
   const double Lescale_m = cvac/m_wpescale_Hz;
   const double Liscale_m = cvac/m_wpiscale_Hz;

   m_delta0   = m_epRel*pow(Vscale/cvac,2);
   m_Le0sq    = pow(Lescale_m/a_Xscale,2);
   m_epsilon0 = m_meRel*m_Le0sq;
   m_lambdai0 = Liscale_m/a_Xscale;
   
   m_Tscale_eV = a_Tscale;
   m_tscale_s  = tscale;
   m_tauescale1eV_s = 3.44e5/10.0*pow(1.0,1.5)/(a_Nscale/1.0e6);
   m_tauiscale1eV_s = sqrt(2.0/m_mM)*m_tauescale1eV_s;

   m_taue0 = 3.44e5/10.0*pow(m_Tscale_eV,1.5)/(a_Nscale/1.0e6)/m_tscale_s;
   m_taui0 = sqrt(2.0/m_mM)*m_taue0;

   //const double Vte_mps = sqrt(qe/me)*sqrt(m_Tscale_eV); // char ele therm speed [m/s]
   m_Ve0sq = Mn/me; // Vte_mps*Vte_mps/Vscale/Vscale;
   m_Vi0sq = Mn/Mi; // m_mM*m_Ve0sq;
   
   // overwrite Tmax so that free-streaming ele heat flux cant be larger than 
   // effective speed of light ?
   //const double Tmax_fsehf = cvac*cvac/m_epRe./(qe/me)/m_heatLim_ele/m_heatLim_ele;

   // calculate normalized resistivity at Te=Tscale and N=Nscale 
   // just for print out purposes
   //
   //const double etascale = pow(a_Xscale,2)*mu0/m_tscale_s;
   //const double etanorm1 = 1.03e-4*10.0/pow(a_Tscale,1.5)/etascale;
   const double etanorm = m_Le0sq*m_tscale_s/m_tauescale1eV_s/pow(a_Tscale,1.5);
   
   if(!procID() && m_verbosity>0) {
      cout << "fundamental scales in ExtendedMhdOp: " << endl;
      cout << "length scale [cm]              = " << a_Xscale*1e2 << endl;
      cout << "density scale [1/cm^3]         = " << a_Nscale/1e6 << endl;
      cout << "species temperature scale [eV] = " << a_Tscale << endl;
      cout << "atomic mass                    = " << a_Mscale << endl;
      cout << endl;
      cout << "derived scales in ExtendedMhdOp: " << endl;
      cout << "species pressure scale [J/m^3] = " << Pscale << endl;
      cout << "velocity scale [m/s]           = " << Vscale << endl;
      cout << "current density scale [A/m^2]  = " << Jscale << endl;
      cout << "magnetic field scale [T]       = " << Bscale << endl;
      cout << "electric field scale [V/m]     = " << Escale << endl;
      cout << "time scale [s]                 = " << tscale << endl;
      cout << "wpe scale [Hz]                 = " << m_wpescale_Hz << endl;
      cout << endl;
      cout << "dimensionless scales in ExtendedMhdOp: " << endl;
      cout << "m_meRel   = " << m_meRel << endl;
      cout << "m_epRel   = " << m_epRel << endl;
      cout << "m_lambdai0 = Li/X        = " << m_lambdai0 << endl; 
      cout << "m_epsilon0 = (Leeff/X)^2 = " << m_epsilon0 << endl;
      cout << "m_delta0   = (V/ceff)^2  = " << m_delta0 << endl;
      cout << "m_taui0    = taui/tscale = " << m_taui0 << endl;
      cout << "m_taue0    = taue/tscale = " << m_taue0 << endl;
      cout << "etanorm(Tscale,Nscale)  = " << etanorm << endl;
      cout << endl;
   }

}

void ExtendedMhdOp::printParameters()
{
   if (procID()==0) {
      cout << "ExtendedMhdOp parameters:" << std::endl;
      cout << "advScheme   =  " << m_advScheme << std::endl;
      cout << "updateHydroVars      =  " << m_updateHydroVars << std::endl;
      cout << "useRelaxedViscosity  =  " << m_useRelaxedViscosity << std::endl;
      cout << "etaVis_ele        =  " << m_etaVis_ele << std::endl;
      cout << "etaResMin         =  " << m_etaResMin << std::endl;
      cout << "etaAnomFactor     =  " << m_etaAnomFactor << std::endl;
      cout << "NvacC             =  " << m_NvacC << std::endl;
      cout << "NvacP             =  " << m_NvacP << std::endl;
      cout << "Nmin              =  " << m_Nmin << std::endl;
      cout << "Tmin              =  " << m_Tmin << std::endl;
      cout << "Tmax              =  " << m_Tmax << std::endl;
      cout << "etaMin            =  " << m_etaMin << std::endl;
      cout << "viscLimiter       =  " << m_viscLimiter << std::endl;
      cout << "gyroViscLimiter   =  " << m_gyroViscLimiter << std::endl;
      cout << "gyroHeatViscCoeff =  " << m_gyroHeatViscCoeff << std::endl;
      cout << "useRelaxedHeatFlux   =  " << m_useRelaxedHeatFlux << std::endl;
      cout << "tauiMin(code units)  =  " << m_tauiMin << std::endl;
      cout << "taueMin(code units)  =  " << m_taueMin << std::endl;
      cout << "kappaMin          =  " << m_kappaMin << std::endl;
      cout << "heatLimiter       =  " << m_heatLimiter << std::endl;
      cout << "heatLimiter_ele   =  " << m_heatLimiter_ele << std::endl;
      cout << "gyroHeatLimiter   =  " << m_gyroHeatLimiter << std::endl;
      cout << "xionFactor        =  " << m_xionFactor << std::endl;
      cout << "xeleFactor        =  " << m_xeleFactor << std::endl;
      cout << "useRelaxedOhmsLaw   =  " << m_useRelaxedOhmsLaw << std::endl;
      cout << "include_divJStress  =  " << m_include_divJStress << std::endl;
      cout << "useJ0forJcc         =  " << m_useJ0forJcc << std::endl;
      cout << "advSchemeHall    =  " << m_advSchemeHall << std::endl;
      cout << "advSchemeJStress =  " << m_advSchemeJStress << std::endl;
      cout << "modelDriftTerms =  " << m_modelDriftTerms << std::endl;
      cout << "initializeJfromAmperesLaw =  " << m_initializeJfromAmperesLaw << std::endl;
      cout << "initializeEfromOhmsLaw =  " << m_initializeEfromOhmsLaw << std::endl;
      cout << "fullyImplicitJ  =  " << m_fullyImplicitJ << std::endl;
      cout << "iterHallMax     =  " << m_iterHallMax << std::endl;
      cout << "HallResTol      =  " << m_HallResTol << std::endl;
      cout << endl;
   }
}

void ExtendedMhdOp::setCourantTimeStep( const LevelData<FArrayBox>&  a_Cspeed )
{
   CH_TIME("ExtendedMhdOp::setCourantTimeStep()");
   CH_assert(a_Cspeed.nComp() == SpaceDim);
   const DisjointBoxLayout& grids( m_geometry.grids() );
    
   Vector<double> local_CSoundmax(SpaceDim,-CH_BADVAL);
   Vector<double> local_CElectmax(SpaceDim,-CH_BADVAL);
   Vector<double> local_CHallDmax(SpaceDim,-CH_BADVAL);
   Vector<double> local_CLightmax(SpaceDim,-CH_BADVAL);
   double local_Nemin(CH_BADVAL);
   double local_Nemax(CH_BADVAL);
   double local_Nmin(CH_BADVAL);
   double local_Nmax(-CH_BADVAL);
   RealVect dX_mapped;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      dX_mapped = coord_sys.getMappedCellSize();
      
      const Box& box( grids[dit] ); // this box has no ghosts
      for (int dir=0; dir<SpaceDim; ++dir) {
         double box_max( a_Cspeed[dit].max( box,dir ) );
         local_CSoundmax[dir] = Max( local_CSoundmax[dir], box_max );
         //
         box_max = m_CLight_cc[dit].max( box,dir );
         local_CLightmax[dir] = Max( local_CLightmax[dir], box_max );
         //
         box_max = m_Ve_cc[dit].max( box,dir );
         local_CElectmax[dir] = Max( local_CElectmax[dir], box_max );
         //
         box_max = m_VhallDrift_cc[dit].max( box,dir );
         local_CHallDmax[dir] = Max( local_CHallDmax[dir], box_max );
      }
      local_Nemin = m_eleDen_cc[dit].min( box, 0 );
      local_Nemax = m_eleDen_cc[dit].max( box, 0 );
      local_Nmin = m_rhoDen_cc[dit].min( box, 0 );
      local_Nmax = m_rhoDen_cc[dit].max( box, 0 );
   }
   
   Vector<double> CSoundmax( local_CSoundmax );
   Vector<double> CElectmax( local_CElectmax );
   Vector<double> CHallDmax( local_CHallDmax );
   Vector<double> CLightmax( local_CLightmax );
   double Nemin( local_Nemin );
   double Nemax( local_Nemax );
   double Nmin( local_Nmin );
   double Nmax( local_Nmax );
   double dtSound_max=DBL_MAX;
   double dtElect_max=DBL_MAX;
   double dtHallD_max=DBL_MAX;
   double dtLight_max=DBL_MAX;
   double dtwpe_max=DBL_MAX;
   for (int dir=0; dir<SpaceDim; ++dir) {
#ifdef CH_MPI
      MPI_Allreduce( &local_CSoundmax[dir], &CSoundmax[dir], 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
      MPI_Allreduce( &local_CElectmax[dir], &CElectmax[dir], 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
      MPI_Allreduce( &local_CHallDmax[dir], &CHallDmax[dir], 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
      MPI_Allreduce( &local_CLightmax[dir], &CLightmax[dir], 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif
      dtSound_max = Min( dX_mapped[dir]/CSoundmax[dir], dtSound_max);
      dtElect_max = Min( abs(dX_mapped[dir]/CElectmax[dir]), dtElect_max);
      dtHallD_max = Min( abs(dX_mapped[dir]/CHallDmax[dir]), dtHallD_max);
      dtLight_max = Min( dX_mapped[dir]/CLightmax[dir], dtLight_max);
   }
   //if(procID()) cout << "dX_mapped[0] = " << dX_mapped[0] << endl; // 2*pi/n0 for slab
   //if(procID()) cout << "dX_mapped[1] = " << dX_mapped[1] << endl; // 2*pi/n1 for slab

#ifdef CH_MPI
      MPI_Allreduce( &local_Nemin, &Nemin, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
      MPI_Allreduce( &local_Nemax, &Nemax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
      MPI_Allreduce( &local_Nmin, &Nmin, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
      MPI_Allreduce( &local_Nmax, &Nmax, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif
   dtwpe_max = 1.0/(m_wpescale_Hz*m_tscale_s*pow(Nmax/m_meRel,0.5));
   if(!procID()) {
      cout << "dtcfl_sound = " << dtSound_max << endl;
      cout << "dtcfl_Elect = " << dtElect_max << endl;
      cout << "dtcfl_HallD = " << dtHallD_max << endl;
      cout << "dtcfl_light = " << dtLight_max << endl;
      cout << "dt_wpe      = " << dtwpe_max << endl;
      cout << "min(Ne)     = " << Nemin << endl;
      cout << "max(Ne)     = " << Nemax << endl;
      cout << "min(rho)    = " << Nmin << endl;
      cout << "max(rho)    = " << Nmax << endl;
   }
   if(std::isnan(Nmin*Nmax)) MayDay::Error("ExendedMhdOp: bad density!!!  ");

   m_courant_time_step = dtLight_max;
   m_courant_time_step = Min( m_courant_time_step, dtSound_max );
   if(m_modelDriftTerms) {
      m_courant_time_step = Min( m_courant_time_step, dtElect_max );
      m_courant_time_step = Min( m_courant_time_step, dtHallD_max );
   }

}

void ExtendedMhdOp::setMappedCspeed( const double  a_gamma,
                                     const int     a_Ptherm,
                                     const int     a_Pmag )
{
   CH_TIME("ExtendedMhdOp::setMappedCspeed()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // compute flux-freezing speed at cell-center for each direction
   // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/N + B^2/N); i = q0, q1, q2

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_Cspeed_cc[dit].copy(m_velocity[dit]);
   }
   m_geometry.convertPhysToContravar(m_Cspeed_cc,0);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& rho_on_patch  = m_rhoDen_cc[dit];   
      const FArrayBox& P_on_patch    = m_pressure[dit]; 
      const FArrayBox& magP_on_patch = m_Bpressure[dit]; 
      const FArrayBox& NTrowJinverse_on_patch = m_NTrowJinverse[dit]; 
      FArrayBox& C_on_patch = m_Cspeed_cc[dit];
      const Box& thisbox = C_on_patch.box();
      
      FORT_EVAL_FLUX_FREEZING_SPEED_EMHD( CHF_BOX(thisbox), 
                                          CHF_CONST_FRA1(rho_on_patch,0),
                                          CHF_CONST_FRA1(P_on_patch,0),
                                          CHF_CONST_FRA1(magP_on_patch,0),
                                          CHF_CONST_FRA(NTrowJinverse_on_patch),
                                          CHF_CONST_REAL(a_gamma),
                                          CHF_CONST_REAL(m_delta0),
                                          CHF_CONST_INT(a_Ptherm),
                                          CHF_CONST_INT(a_Pmag),
                                          CHF_FRA(C_on_patch) );
   }
      
   // interp to cell faces
   SpaceUtils::upWindToFaces(m_CspeedR_norm, m_Cspeed_cc, m_CspeedR_norm, "c2");
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_CspeedL_norm[dit].copy(m_CspeedR_norm[dit],m_CspeedL_norm[dit].box());
      m_CspeedL_norm[dit].negate();
   }

}

void ExtendedMhdOp::setMappedCLight()
{  
   CH_TIME("ExtendedMhdOp::setMappedCLight()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // comute mapped speed of light at cell-center for each direction
   // CLight_i = |N_i|*1/sqrt(m_delta0); i = q0, q1, q2

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_CLight_cc[dit].setVal(1.0/sqrt(m_delta0));
      for (int dir=0; dir<SpaceDim; ++dir) {
         m_CLight_cc[dit].mult(m_NTrowJinverse[dit],dir,dir,1);
      }
   }

}


Real ExtendedMhdOp::computeDtImExTI( const FluidSpeciesPtrVect&  a_fluid_comp )
{
   CH_TIME("Real ExtendedMhdOp::computeDtImExTI()");   
   
   double thisdt=DBL_MAX;
   if(m_courant_time_step==DBL_MAX) {
      thisdt = computeDtMhd( a_fluid_comp );
      return thisdt;
   }
   else {
      return m_courant_time_step;
   }

}

Real ExtendedMhdOp::computeDtExplicitTI( const FluidSpeciesPtrVect&  a_fluid_comp )
{
   CH_TIME("Real ExtendedMhdOp::computeDtExplicit()");   
   double thisdt=DBL_MAX;
   thisdt = computeDtMhd( a_fluid_comp );
   return thisdt;
}

Real ExtendedMhdOp::computeDtMhd( const FluidSpeciesPtrVect&  a_fluid_comp )
{
   CH_TIME("Real ExtendedMhdOp::computeDtMhd()");   
   
   // get physical values for m_species_name and calculate freezing speed
   // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/rho); i = q0, q1, q2
   
   double gamma;
   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (int species(0); species<a_fluid_comp.size(); species++) {
      const FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_comp[species])) );
      const std::string species_name( fluid_species.name() );
      if(species_name==m_species_name) {
         if(fluid_species.m_evolve_magneticField_virtual) {
            const LevelData<FArrayBox>& soln_By( fluid_species.cell_var("magneticField_virtual") );
            for (DataIterator dit(grids); dit.ok(); ++dit) {
               const Box& thisbox = m_By_covar_cc[dit].box();
               m_By_covar_cc[dit].copy( soln_By[dit], thisbox); // comp
               if(m_twoDaxisymm) {
                  m_By_covar_cc[dit].mult(m_g_y[dit]);
               }
            }
         } else {
            for (DataIterator dit(grids); dit.ok(); ++dit) {
               m_By_covar_cc[dit].setVal(0.0);
            }
         }

         gamma = fluid_species.m_gamma;
         fluid_species.massDensity(m_rhoDen_cc);
         fluid_species.electronDensity(m_eleDen_cc);
         fluid_species.velocity(m_velocity);  // in-plane velocity vector
         fluid_species.pressure(m_pressure_ion);
         fluid_species.pressure_ele(m_pressure_ele);
         for (DataIterator dit(m_pressure.dataIterator()); dit.ok(); ++dit) {
            m_pressure[dit].copy( m_pressure_ion[dit] );
            m_pressure[dit].plus( m_pressure_ele[dit] );
         }
         fluid_species.Bpressure(m_Bpressure);
   
         // Calls for physical variables here are from mapped state
         // vector (has not been convertToPhysical )
         // and so need to divide by J
         m_geometry.divideJonValid(m_eleDen_cc);
         m_geometry.divideJonValid(m_rhoDen_cc);
         m_geometry.divideJonValid(m_pressure);
         m_geometry.divideJonValid(m_Bpressure);
         m_geometry.divideJonValid(m_Bpressure);
         m_geometry.divideJonValid(m_By_covar_cc);
  
         // need to fill ghosts for eleDen for grad calc
         // (comp state passed here has no ghosts)
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
            const ProblemDomain& domain = block_coord_sys.domain();
            fourthOrderCellExtrapAtDomainBdry(m_eleDen_cc[dit], domain, grids[dit]);
            fourthOrderCellExtrapAtDomainBdry(m_By_covar_cc[dit], domain, grids[dit]);
         }
         m_eleDen_cc.exchange();

         // compute contravariant electron velocity magnitude
         if(fluid_species.m_evolve_currentDensity) {
            SpaceUtils::interpCellToEdges(m_eleDen_ce, m_eleDen_cc, m_eleDen_ce,"c2");
            const LevelData<EdgeDataBox>& soln_Jcovar( fluid_species.edge_var("currentDensity") );
            computeHallVelocity( m_Vhall_cc, soln_Jcovar, m_eleDen_ce ); // phys Ve-Vi at cell center
            for (DataIterator dit(grids); dit.ok(); ++dit) {
               m_Ve_cc[dit].copy( m_Vhall_cc[dit] );
            }
            m_geometry.convertPhysToContravar(m_Ve_cc,0);
            for (DataIterator dit(grids); dit.ok(); ++dit) {
               m_Ve_cc[dit].abs( ); // only used for time-step
            }
         }
 
         break;
      }
   }
   setMappedCspeed( gamma, 1, 1 );
   setHallDriftSpeed(m_eleDen_cc, m_By_covar_cc);

   setCourantTimeStep(m_Cspeed_cc);   
   return m_courant_time_step;

}

void ExtendedMhdOp::defineLevelDatas( const DisjointBoxLayout&  a_grids,
                                      const IntVect&            a_ghostVect )
{
   // Initialize exchangeCopier with a_ghostVect
   m_exchangeCopier.define(a_grids, a_grids, a_ghostVect, true);

   // get the transformation matrix N and Jacobian at cell-center needed to transform 
   // from physical space to mapped space
   m_Jacobian.define(a_grids, 1, a_ghostVect);
   m_Jacobian_cf.define(a_grids, 1, a_ghostVect);
   m_Nmatrix.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);
   m_dxidX_ce.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);
   m_Xphys.define(a_grids, SpaceDim, a_ghostVect);
   m_g_y.define(a_grids, 1, a_ghostVect); // used only when solving virtual B
   m_cellVol.define(a_grids, 1, IntVect::Zero);
   double twoPi = 2.0*Constants::PI;
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
   m_geometry.getJ(m_Jacobian_cf);
   m_geometry.getN(m_Nmatrix);
   m_geometry.getPointwiseNJInverse(m_dxidX_ce);
   m_geometry.getCellVolumes(m_cellVol);
   
   // define some metric scale factors for time step calc
   m_NTrowSq.define(a_grids, SpaceDim, a_ghostVect);        // square of magnitude of rows of N=J*dxi/dX transpose
   m_NTrowJinverse.define(a_grids, SpaceDim, a_ghostVect);  // magnitude of rows of dxi/dX transpose
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      const FArrayBox& thisN = m_Nmatrix[dit];
            FArrayBox& thisNTrowSq = m_NTrowSq[dit];

      const Box& cell_box( thisNTrowSq.box() );
      FORT_MAG_NT_ROW_SQ( CHF_BOX(cell_box),
                          CHF_CONST_FRA(thisN),
                          CHF_FRA(thisNTrowSq) );

      BoxIterator bit(cell_box);
      IntVect ib;
      Real this_val;
      FArrayBox& thisNTrowJinverse = m_NTrowJinverse[dit];
      for(bit.begin(); bit.ok(); ++bit) {
         ib = bit();
         for (int n=0; n<SpaceDim; ++n) {
            this_val = sqrt(thisNTrowSq.get(ib,n));
            thisNTrowJinverse.set(ib,n,this_val);
         }
      }
   }
   m_geometry.divideJonValid( m_NTrowJinverse );
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_NTrowJinverse[dit].abs( ); // Ja can be negative in ghost cells at axis (2*pi*R)
   }

   // containers for viscosity tensor
   m_DiffVx_cf.define(a_grids, SpaceDim, a_ghostVect);
   m_DiffVz_cf.define(a_grids, SpaceDim, a_ghostVect);
   m_DiffVth.define(a_grids, SpaceDim, a_ghostVect);
   m_etaVisC.define(a_grids, 3, a_ghostVect);
   m_etaVisG.define(a_grids, 2, a_ghostVect);
   m_JaPion.define(a_grids, 1, a_ghostVect);
   m_momVisc_source.define(a_grids, 1, IntVect::Zero);
   m_etaVisC_cf.define(a_grids, 3, a_ghostVect);
   m_etaVisG_cf.define(a_grids, 2, a_ghostVect);
   m_JaPion_cf.define(a_grids, 1, a_ghostVect);
   m_JaW_cf.define(a_grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_m0JaFluxVisc_cf.define(a_grids, SpaceDim, IntVect::Unit); // 1 ghost needed to do exchange
   m_m1JaFluxVisc_cf.define(a_grids, SpaceDim, IntVect::Unit); // 1 ghost needed to do exchange
   m_enJaFluxVisc_cf.define(a_grids, SpaceDim, a_ghostVect);
   m_m0JaFluxVisc_norm.define(a_grids, 1, a_ghostVect);
   m_m1JaFluxVisc_norm.define(a_grids, 1, a_ghostVect);
   m_momJaFluxVisc_norm.define(a_grids, SpaceDim, a_ghostVect);
   m_enJaFluxVisc_norm.define(a_grids, 1, a_ghostVect);
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_momVisc_source[dit].setVal(0.0);
   }
   //
   m_qi_wedge_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_JaWq_cf.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);   
   
   // containers for viscosity in Ohm's law
   m_gradVepara_cn.define(a_grids, SpaceDim, a_ghostVect);   
   m_gradVepara_cc.define(a_grids, SpaceDim, a_ghostVect);   
   m_gradVeperp_cc.define(a_grids, SpaceDim, a_ghostVect);   
   m_divVe_cn.define(a_grids, 1, a_ghostVect);
   m_divVe_ce.define(a_grids, 1, a_ghostVect);
   m_Pie_diag_cn.define(a_grids, SpaceDim, a_ghostVect);   
   m_Pie_diag_cn_old.define(a_grids, SpaceDim, a_ghostVect);   
   m_Pie_diag_cn0.define(a_grids, SpaceDim, a_ghostVect);   
   m_Pie_diagD_cn.define(a_grids, 1, a_ghostVect);   
   m_Pie_xz_cc.define(a_grids, 1, a_ghostVect);   
   m_Pie_xz_cc0.define(a_grids, 1, a_ghostVect);   
   m_Pie_xz_cc_old.define(a_grids, 1, a_ghostVect);   
   m_Pie_xzD_cc.define(a_grids, 1, a_ghostVect);   
   m_Pie_yy_ce.define(a_grids, 1, a_ghostVect);   
   m_Pie_yy_ce0.define(a_grids, 1, a_ghostVect);   
   m_Pie_yy_ce_old.define(a_grids, 1, a_ghostVect);   
   m_Pie_yyD_ce.define(a_grids, 1, a_ghostVect);   
   m_divPie_ce.define(a_grids, 1, a_ghostVect);   
   m_enEleFluxVisc_cf.define(a_grids, SpaceDim, a_ghostVect);
   m_enEleFluxVisc_norm.define(a_grids, 1, a_ghostVect);
   
   // below are for relaxation scheme
   m_rhoDen_cf.define(a_grids, 1, a_ghostVect);
   m_m0JaFluxVisc_cf_old.define(a_grids, SpaceDim, IntVect::Zero);
   m_m1JaFluxVisc_cf_old.define(a_grids, SpaceDim, IntVect::Zero);
   m_momVisc_source_old.define(a_grids, 1, IntVect::Zero);
   m_m0JaFluxVisc_cf0.define(a_grids, SpaceDim, IntVect::Zero);
   m_m1JaFluxVisc_cf0.define(a_grids, SpaceDim, IntVect::Zero);
   m_momVisc_source0.define(a_grids, 1, IntVect::Zero);
   
   // containers for heat flux
   m_Tele_cc.define(a_grids, 1, a_ghostVect);
   m_Tion_cc.define(a_grids, 1, a_ghostVect);
   m_Tele_cf.define(a_grids, 1, a_ghostVect);
   m_Tion_cf.define(a_grids, 1, a_ghostVect);
   m_Dqe.define(a_grids, SpaceDim, a_ghostVect);
   m_Dqi.define(a_grids, SpaceDim, a_ghostVect);
   //
   m_qe_old.define(a_grids, SpaceDim, IntVect::Unit);
   m_qe_cf.define(a_grids, SpaceDim, IntVect::Unit);
   m_qe_norm.define(a_grids, 1, a_ghostVect);
   m_qe0.define(a_grids, SpaceDim, a_ghostVect);
   //
   m_qi_old.define(a_grids, SpaceDim, IntVect::Unit);
   m_qi_cf.define(a_grids, SpaceDim, IntVect::Unit);
   m_qi_norm.define(a_grids, 1, a_ghostVect);
   m_qi0.define(a_grids, SpaceDim, a_ghostVect);
   //
   m_kappae.define(a_grids, 1, a_ghostVect);
   m_kappai.define(a_grids, 1, a_ghostVect);
   m_kappae_wedge.define(a_grids, 1, a_ghostVect);
   m_kappai_wedge.define(a_grids, 1, a_ghostVect);
   //
   m_kappae_cf.define(a_grids, 1, a_ghostVect);
   m_kappai_cf.define(a_grids, 1, a_ghostVect);
   m_kappae_wedge_cf.define(a_grids, 1, a_ghostVect);
   m_kappai_wedge_cf.define(a_grids, 1, a_ghostVect);

   ////////////////////////////////////////////

   // cell center vars
   m_rhoDen_cc.define(a_grids, 1, a_ghostVect);
   m_momDen_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_eneDeni_cc.define(a_grids, 1, a_ghostVect);
   m_eneDene_cc.define(a_grids, 1, a_ghostVect);
   m_momDen_virtual_cc.define(a_grids, 1, a_ghostVect);
   m_magField_virtual_cc.define(a_grids, 1, a_ghostVect);
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_magField_virtual_cc[dit].setVal(0.0);
   }
   //
   m_pressure.define(a_grids, 1, a_ghostVect);
   m_pressure_ele.define(a_grids, 1, a_ghostVect);
   m_pressure_ion.define(a_grids, 1, a_ghostVect);
   m_Bpressure.define(a_grids, 1, a_ghostVect);
   m_velocity.define(a_grids, SpaceDim, a_ghostVect);
   m_velocity_cf.define(a_grids, SpaceDim, a_ghostVect);
   m_velocity_virtual.define(a_grids, 1, a_ghostVect);
   //
   m_eleDen_cc.define(a_grids, 1, a_ghostVect);
   m_eleDen_cf.define(a_grids, 1, a_ghostVect);
   m_Zbar_cc.define(a_grids, 1, a_ghostVect);
   m_eleFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_eleFlux_norm.define(a_grids, 1, a_ghostVect);
   m_JaSeleDen.define(a_grids, 1, a_ghostVect);
   m_nue_izn.define(a_grids, 1, a_ghostVect);
   m_nue_neu.define(a_grids, 1, a_ghostVect);
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_nue_izn[dit].setVal( 0.0 );
      m_nue_neu[dit].setVal( 0.0 );
   }
   //
   m_rhoFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mxFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_myFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mzFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_enFluxi_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_enFluxe_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mvFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   //
   m_Cspeed_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_CspeedL_norm.define(a_grids, 1, a_ghostVect);
   m_CspeedR_norm.define(a_grids, 1, a_ghostVect);
   m_CLight_cc.define(a_grids, SpaceDim, IntVect::Zero);
   //
   m_rhoFlux_norm.define(a_grids, 1, IntVect::Zero);
   m_mxFlux_norm.define(a_grids, 1, IntVect::Zero);
   m_myFlux_norm.define(a_grids, 1, IntVect::Zero);
   m_mzFlux_norm.define(a_grids, 1, IntVect::Zero);
   m_enFluxi_norm.define(a_grids, 1, IntVect::Zero);
   m_enFluxe_norm.define(a_grids, 1, IntVect::Zero);
   m_enFluxeHall_norm.define(a_grids, 1, IntVect::Zero);
   m_momFlux_norm.define(a_grids, SpaceDim, IntVect::Zero);
   m_mvFlux_norm.define(a_grids, 1, IntVect::Zero);
   m_div_IdentFlux_R.define(a_grids, 1, IntVect::Zero);
   m_div_IdentFlux_Z.define(a_grids, 1, IntVect::Zero);
   LevelData<FluxBox> IdentFlux_phys_R( a_grids, SpaceDim, a_ghostVect );
   LevelData<FluxBox> IdentFlux_phys_Z( a_grids, SpaceDim, a_ghostVect );
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
   
   ////////////////////
   
   //virtualEleField.define(a_grids, 1, a_ghostVect);
   m_divE.define(a_grids, 1, IntVect::Zero);
   m_deltaN_ce.define(a_grids, 1, IntVect::Zero);
   m_deltaN_cc.define(a_grids, 1, IntVect::Zero);
   //
   m_E0_ce.define(a_grids, 1, a_ghostVect);
   m_J0_ce.define(a_grids, 1, a_ghostVect);
   m_J0_cc.define(a_grids, SpaceDim, IntVect::Zero);
   m_E0_cc.define(a_grids, SpaceDim, IntVect::Zero);
   m_curlE.define(a_grids, 1, IntVect::Zero);
   m_By_contr_cc.define(a_grids, 1, a_ghostVect);
   m_By_covar_cc.define(a_grids, 1, a_ghostVect);
   m_JaSeneDeni.define(a_grids, 1, a_ghostVect);   // ion energy den source
   m_JaSeneDene.define(a_grids, 1, a_ghostVect);   // ele energy den source
   m_JaNUdotE.define(a_grids, 1, a_ghostVect);     // Jacobian*N*U\cdotE
   m_JaJdotE.define(a_grids, 1, IntVect::Zero);  // Jacobian*J\cdotE
   m_JaJcrossBv_cc.define(a_grids, SpaceDim, IntVect::Zero); // Ja*JxBvirt at cell center
   
   ///////////////////

   m_EdgeBC_zeros.define(a_grids, 1, a_ghostVect);
   m_FluxBC_zeros.define(a_grids, 1, a_ghostVect);
   m_momFluxBC_norm.define(a_grids, SpaceDim, IntVect::Zero);
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_EdgeBC_zeros[dit].setVal(0.0);
      m_FluxBC_zeros[dit].setVal(0.0);
   }
   
   m_divMassFluxBinary.define(a_grids, 1, IntVect::Zero);
   m_dummyDiv.define(a_grids, 1, IntVect::Zero);
   m_dummyDiv_mom.define(a_grids, SpaceDim, IntVect::Zero);
   m_dummyFArray_oneComp.define(a_grids, 1, a_ghostVect);
   m_dummyFArray_spaceDim.define(a_grids, SpaceDim, a_ghostVect);
   m_dummyFluxBC.define(a_grids, 1, IntVect::Zero);
   m_dummyFlux_oneComp.define(a_grids, 1, a_ghostVect);
   m_dummyFlux_spaceDim.define(a_grids, SpaceDim, a_ghostVect);
   m_dummyEdge_oneComp.define(a_grids, 1, a_ghostVect);

   m_eleDen_ce.define(a_grids, 1, a_ghostVect);
   m_Eold.define(a_grids, 1, IntVect::Zero);
   m_Jold.define(a_grids, 1, IntVect::Zero);
   m_Ecc.define(a_grids, SpaceDim, a_ghostVect);
   m_Jcc.define(a_grids, SpaceDim, a_ghostVect);
   m_Jcc_res.define(a_grids, SpaceDim, IntVect::Zero); // for implicit hall solve

   /////////////////

   m_Ehall_ce.define(a_grids, 1, a_ghostVect);          // Ehall = -VhallxB covar at cell edges
   m_Ehall_op.define(a_grids, 1, IntVect::Zero);          // Ehall = -VhallxB covar at cell edges
   m_gfun_covar.define(a_grids, 1, a_ghostVect);        // Ehall = -VhallxB covar at cell edges
   m_gfun_covar_op.define(a_grids, 1, IntVect::Zero);     // Ehall = -VhallxB covar at cell edges
   m_Vhall_cc.define(a_grids, SpaceDim, a_ghostVect);   // contra Vhall = -lambdai*J/Ne at cell center
   m_Vhall_cf.define(a_grids, 1, a_ghostVect);          // contra Vhall = -lambdai*J/Ne at cell faces
   m_VhallDrift_cc.define(a_grids, SpaceDim, IntVect::Zero); // VhallD = lambdai*Bcrossgrad(Ne)/Ne^2 contravar at cell center
   m_hallDriftScaleFactor.define(a_grids, 1, IntVect::Zero);
   m_gammaHall_ce.define(a_grids, 1, IntVect::Unit);
   m_gammaHallold_ce.define(a_grids, 1, IntVect::Zero);
   m_Ve_cc.define(a_grids, SpaceDim, a_ghostVect);      // Ve = U - lambdai*J/Ne at cell center
   m_Ve_ce.define(a_grids, 1, a_ghostVect);             // Ve = U - lambdai*J/Ne contravar at edge center
   m_Ve_cf.define(a_grids, SpaceDim, a_ghostVect);      // Ve = U - lambdai*J/Ne physical at face center
   m_etaHall_ce.define(a_grids, 1, a_ghostVect);        // Ehall = lambdai/Ne*JxB = J x etaHall
   m_etaHall_op.define(a_grids, 1, IntVect::Zero);
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_Ehall_ce[dit].setVal( 0.0 );
      m_Ehall_op[dit].setVal( 0.0 );
      m_etaHall_ce[dit].setVal( 0.0 );
      m_etaHall_op[dit].setVal( 0.0 );
      m_gammaHall_ce[dit].setVal( 1.0 );
   } 

   m_divJStress_covar.define(a_grids, 1, a_ghostVect);  // covariant div(UJ+JUe)
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_divJStress_covar[dit].setVal( 0.0 );
      m_hallDriftScaleFactor[dit].setVal( 0.0 );
   }

   m_gradPe_covar_cc.define(a_grids, SpaceDim, IntVect::Unit);
   m_gradTeeff_covar_cc.define(a_grids, SpaceDim, IntVect::Unit);
   m_gradT_cc.define(a_grids, SpaceDim, 4*IntVect::Unit);
   m_gradPe_covar_ce.define(a_grids, 1, a_ghostVect);
   
   /////////////////

   m_taui.define(a_grids, 1, a_ghostVect);  // ion collision time
   m_taue.define(a_grids, 1, a_ghostVect);  // electron collision time
   m_etaRes.define(a_grids, 1, a_ghostVect);      // collisional resistivity
   m_etaRes_BC.define(a_grids, 1, a_ghostVect);   // insulator BC resistivity
   for (DataIterator dit(a_grids); dit.ok(); ++dit) {
      m_etaRes_BC[dit].setVal(m_etaResMin); // need to init to zero
   }
   m_etaRes_ce.define(a_grids, 1, a_ghostVect);   // collisional resistivity on edges
   m_Qie.define(a_grids, 1, a_ghostVect);   // thermal energy exchange term
   m_xion.define(a_grids, 1, a_ghostVect);    // wci*taui
   m_xele.define(a_grids, 1, a_ghostVect);    // wce*taue

}

void ExtendedMhdOp::setCellCenterValues( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("ExtendedMhdOp::setCellCenterValues()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   const LevelData<FArrayBox>& soln_rhoDen( a_soln_fluid.cell_var(0) );
   const LevelData<FArrayBox>& soln_momDen( a_soln_fluid.cell_var("momentumDensity") );
   const LevelData<FArrayBox>& soln_eneDeni( a_soln_fluid.cell_var("energyDensity") );
   const LevelData<FArrayBox>& soln_eneDene( a_soln_fluid.cell_var("energyDensity_ele") );
   
   // set values for MHD variables
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_rhoDen_cc[dit].copy( soln_rhoDen[dit], m_rhoDen_cc[dit].box() );
      m_momDen_cc[dit].copy( soln_momDen[dit], m_momDen_cc[dit].box() );
      m_eneDeni_cc[dit].copy( soln_eneDeni[dit], m_eneDeni_cc[dit].box() );
      m_eneDene_cc[dit].copy( soln_eneDene[dit], m_eneDene_cc[dit].box() );
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
   if(a_soln_fluid.m_evolve_electronDensity) {
      const LevelData<FArrayBox>& soln_electronDensity( a_soln_fluid.cell_var("electronDensity") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_eleDen_cc[dit].copy( soln_electronDensity[dit] );
         m_Zbar_cc[dit].copy( soln_electronDensity[dit] );
         m_Zbar_cc[dit].divide( soln_rhoDen[dit] );
      }
   } 
   else {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_eleDen_cc[dit].copy(soln_rhoDen[dit]);
         m_Zbar_cc[dit].setVal( 1.0 );
      }
   }
   
   // set physical derived variables for fluid species
   a_soln_fluid.pressure(m_pressure_ion);
   a_soln_fluid.pressure_ele(m_pressure_ele);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_pressure[dit].copy( m_pressure_ion[dit] );
      m_pressure[dit].plus( m_pressure_ele[dit] );
      //
      m_Tele_cc[dit].copy( m_pressure_ele[dit] );
      m_Tele_cc[dit].divide( m_eleDen_cc[dit] );
      //
      m_Tion_cc[dit].copy( m_pressure_ion[dit] );
      m_Tion_cc[dit].divide( m_rhoDen_cc[dit] );
   }
   a_soln_fluid.Bpressure(m_Bpressure);
   a_soln_fluid.velocity(m_velocity);                 // in-plane velocity vector
   a_soln_fluid.velocity_virtual(m_velocity_virtual); // out-of-plane velocity vector
   
}   

void ExtendedMhdOp::setCellCenterFluxes( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("ExtendedMhdOp::setCellCenterFluxes()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // Compute cell-center fluxes for each variable
   //
   // rhoFlux_x = mx           rhoFlux_y = my           rhoFlux_z = mz
   // eleFlux_x = Zbar*mx      eleFlux_y = Zbar*my      eleFlux_z = Zbar*mz
   // mxFlux_x  = mx*ux + P    mxFlux_y  = mx*uy        mxFlux_z  = mx*uz
   // myFlux_x  = my*ux        myFlux_y  = my*uy + P    myFlux_z  = my*uz
   // mzFlux_x  = mz*ux        mzFlux_y  = mz*uy        mzFlux_z  = mz*uz + P
   // enFluxi_x = (eni+Pi)*ux  enFluxi_y = (eni+Pi)*uy  enFluxi_z = (eni+Pi)*uz
   // enFluxe_x = (ene+Pe)*ux  enFluxe_y = (ene+Pe)*uy  enFluxe_z = (ene+Pe)*uz
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
         for (int n=0; n<SpaceDim; ++n) {
            m_mvFlux_cc[dit].copy( m_momDen_virtual_cc[dit],0,n,1 );
            m_mvFlux_cc[dit].mult( m_velocity[dit],n,n,1);
         }
      }
      
      const FArrayBox& eneDeni_on_patch = m_eneDeni_cc[dit];   
      const FArrayBox& eneDene_on_patch = m_eneDene_cc[dit];   
      const FArrayBox& momDen_on_patch = m_momDen_cc[dit];   
      const FArrayBox& V_on_patch = m_velocity[dit];   
      const FArrayBox& Pi_on_patch = m_pressure_ion[dit];   
      const FArrayBox& Pe_on_patch = m_pressure_ele[dit];   

      FArrayBox& rhoFlux_on_patch = m_rhoFlux_cc[dit];
      FArrayBox& mxFlux_on_patch  = m_mxFlux_cc[dit];
      FArrayBox& myFlux_on_patch  = m_myFlux_cc[dit];
      FArrayBox& mzFlux_on_patch  = m_mzFlux_cc[dit];
      FArrayBox& enFluxi_on_patch  = m_enFluxi_cc[dit];
      FArrayBox& enFluxe_on_patch  = m_enFluxe_cc[dit];

      const Box& cellbox = rhoFlux_on_patch.box();
      FORT_CELL_CENTER_FLUXES_TWOTEMP( CHF_BOX(cellbox), 
                                       CHF_CONST_FRA1(eneDeni_on_patch,0),
                                       CHF_CONST_FRA1(eneDene_on_patch,0),
                                       CHF_CONST_FRA(momDen_on_patch),
                                       CHF_CONST_FRA(V_on_patch),
                                       CHF_CONST_FRA1(Pi_on_patch,0),
                                       CHF_CONST_FRA1(Pe_on_patch,0),
                                       CHF_FRA(rhoFlux_on_patch),
                                       CHF_FRA(mxFlux_on_patch),
                                       CHF_FRA(myFlux_on_patch),
                                       CHF_FRA(mzFlux_on_patch),
                                       CHF_FRA(enFluxi_on_patch),
                                       CHF_FRA(enFluxe_on_patch) );
    
      if(a_soln_fluid.m_evolve_electronDensity) {
         for (int n=0; n<SpaceDim; ++n) {
            m_eleFlux_cc[dit].copy( m_rhoFlux_cc[dit],n,n,1 );
            m_eleFlux_cc[dit].mult( m_Zbar_cc[dit],0,n,1);
         }
      }
 
   }

}

void ExtendedMhdOp::addHallEnergyFlux( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("ExtendedMhdOp::addHallEnergyFlux()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // Compute normal Hall energy flux on faces and add 
   // to norm electron energy flux
   //
   // enFluxeHall = Vhall*(ene+Pe)
   // Vhall = -J/Ne*lambdai (Ve = V+Vhall)
   //
   // enFluxeHall_norm = Ja*Vhall_contra*(ene+Pe)*faceArea 
   // 
  
   // set contra -J (Vhall ~ -J) on cell faces
   // needed for electron energy density flux
   LevelData<FArrayBox>& Vhall_cc = m_dummyFArray_spaceDim;
   if(m_useJ0forJcc) {
      m_geometry.mappedGridCurlofVirtComp(m_J0_ce,m_By_covar_cc); // Jacobian * curl(By)\cdot g^l
      SpaceUtils::interpEdgesToCell(Vhall_cc,m_J0_ce,"c2");     // Ja*curl(By)_contra at cell center
      m_geometry.divideJonValid( Vhall_cc );                    // contra curl(By) at cell center
   }
   else {
      const LevelData<EdgeDataBox>& soln_Jcovar( a_soln_fluid.edge_var("currentDensity") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_J0_ce[dit].copy( soln_Jcovar[dit] );
      }
      m_geometry.multJonEdges(m_J0_ce);
      SpaceUtils::interpEdgesToCell(Vhall_cc,m_J0_ce,"c2"); // Ja*J_covar at cell center
      m_geometry.convertPhysToCovar(Vhall_cc,1);            // Ja*J_phys at cell center
      m_geometry.convertPhysToContravar(Vhall_cc,0);        // Ja*J_contr at cell center
   }
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      Vhall_cc[dit].mult( -m_lambdai0 );
   }
   SpaceUtils::upWindToFaces(m_Vhall_cf, Vhall_cc, m_Vhall_cf, "c2");

   // upwind (ene+Pe)/Ne from cell center to cell face
   // and multiply by Ja (mult Ja before upwinding?) 
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_oneComp[dit].copy( m_eneDene_cc[dit] );
      m_dummyFArray_oneComp[dit].plus( m_pressure_ele[dit] );
      m_dummyFArray_oneComp[dit].divide( m_eleDen_cc[dit] );
   }
   SpaceUtils::upWindToFaces(m_enFluxeHall_norm, m_dummyFArray_oneComp, m_enFluxeHall_norm, "c2");
  
   // compute norm flux on faces (multiplied by face area)
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const RealVect& faceArea = coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; dir++) {
         m_enFluxeHall_norm[dit][dir].mult( m_Vhall_cf[dit][dir] );
         FArrayBox& this_face = m_enFluxeHall_norm[dit][dir];
         this_face.mult( m_Vhall_cf[dit][dir] );
         this_face.mult(faceArea[dir]);
      }
   
      // add norm hall energy flux to electron energy density flux
      m_enFluxe_norm[dit] += m_enFluxeHall_norm[dit];
   }
 
   // rebase dummy containers 
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_oneComp[dit].setVal( 1./0. );
      m_dummyFArray_spaceDim[dit].setVal( 1./0. );
   }
 
}

void ExtendedMhdOp::setMagneticFieldTerms( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("ExtendedMhdOp::setMagneticFieldTerms()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // define source terms for density, momentum, and energy density equations
   // 
   // Sele  = Ja*(nu_izn - nu_recomb)*ne
   // Smom  = Ja*(JxB)_phys
   // SEion = Ja*(Zbar*momDen*(E_phys - eta*J_phys) + Qie)
   // SEele = Ja*( JdotE_phys - Zbar*momDen*(E_phys - eta*J_phys) + Qie ) - Uizn*Sele
   
   if(a_soln_fluid.m_evolve_magneticField_virtual) {

      const LevelData<FArrayBox>& soln_momDen( a_soln_fluid.cell_var("momentumDensity") );
      const LevelData<FArrayBox>& soln_By( a_soln_fluid.cell_var("magneticField_virtual") );

      // compute contravar and covar components of virtual magnetic field
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const Box& thisbox = m_By_covar_cc[dit].box();
         m_By_covar_cc[dit].copy( m_magField_virtual_cc[dit], thisbox); // phys
         m_By_contr_cc[dit].copy( m_magField_virtual_cc[dit], thisbox); // phys
         if(m_twoDaxisymm) {
            m_By_covar_cc[dit].mult(m_g_y[dit]);
            m_By_contr_cc[dit].divide(m_g_y[dit], thisbox, 0,0,1);
         }
      }

      const LevelData<EdgeDataBox>& soln_Ecovar( a_soln_fluid.edge_var("electricField") );
      SpaceUtils::interpEdgesToCell(m_Ecc,soln_Ecovar,"c2");
      m_geometry.convertPhysToCovar(m_Ecc,1);

      if(m_useJ0forJcc) {
         m_geometry.mappedGridCurlofVirtComp(m_J0_ce,m_By_covar_cc); // Jacobian * curl(By)\cdot g^l
         SpaceUtils::interpEdgesToCell(m_Jcc,m_J0_ce,"c2"); // Ja*contravariant at cell center
      }
      else {
         const LevelData<EdgeDataBox>& soln_Jcovar( a_soln_fluid.edge_var("currentDensity") );
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            m_J0_ce[dit].copy( soln_Jcovar[dit] );
         }
         m_geometry.multJonEdges(m_J0_ce);
         SpaceUtils::interpEdgesToCell(m_Jcc,m_J0_ce,"c2"); // Ja*Jcovar at cell center
         m_geometry.convertPhysToCovar(m_Jcc,1);                // Ja*Jphys at cell center
      }

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int n=0; n<SpaceDim; ++n) {
            m_JaJcrossBv_cc[dit].copy( soln_By[dit],0,n,1 );
         }
         m_JaJcrossBv_cc[dit].mult( m_Jcc[dit],1,0,1 );
         m_JaJcrossBv_cc[dit].negate(0);
         m_JaJcrossBv_cc[dit].mult( m_Jcc[dit],0,1,1 );
      }

      // compute source terms for energy equations
      //  
      // SEion = Ja*(NeUdotE + Qie)
      // SEele = Ja*JdotE - SEion
      //
      // if(m_useFiniteIonCorrections) {
      //    m_JaNeUdotE = Ja*Ne(Vz*(Ezcc - eta*Jzcc) + Vx*(Excc - eta*Jxcc));
      //    m_JaNeUdotE /= m_lambdai0;
      // }
      // else {
      //    NeUdotE = Vz*(Jxcc*By - dPedz) - Vx*(Jzcc*By + dPedx);
      // }
      // JdotE = Jzcc*Ezcc + Jxcc*Excc;

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_dummyFArray_oneComp[dit].copy( m_eneDene_cc[dit] );
         m_dummyFArray_oneComp[dit].plus( m_pressure_ele[dit] );
         m_dummyFArray_oneComp[dit].divide( m_eleDen_cc[dit] );
      }
      m_geometry.computeMappedPoloidalGradientWithGhosts( m_pressure_ele, m_gradPe_covar_cc, 2);
      m_geometry.computeMappedPoloidalGradientWithGhosts( m_dummyFArray_oneComp, m_gradTeeff_covar_cc, 2);
      m_geometry.convertPhysToCovar(m_gradPe_covar_cc,1);    // grad(Pe) phys at cell center
      m_geometry.convertPhysToCovar(m_gradTeeff_covar_cc,1); // grad((Ee+Pe)/Ne) phys at cell center
      m_geometry.divideJonValid( m_Jcc );                    // J_phys at cell center
      int useFiniteIon = 0;
      if(m_modelDriftTerms) useFiniteIon = 1;
      for (DataIterator dit(grids); dit.ok(); ++dit) {

         const FArrayBox& this_Nh        = m_rhoDen_cc[dit];
         const FArrayBox& this_Ne        = m_eleDen_cc[dit];
         const FArrayBox& this_JaJxB     = m_JaJcrossBv_cc[dit];
         const FArrayBox& this_gradPe    = m_gradPe_covar_cc[dit];
         const FArrayBox& this_gradTeeff = m_gradTeeff_covar_cc[dit];

         const FArrayBox& this_NU   = soln_momDen[dit];
         const FArrayBox& this_Ja   = m_Jacobian[dit];
         const FArrayBox& this_Ecc  = m_Ecc[dit];
         const FArrayBox& this_Jcc  = m_Jcc[dit];
         const FArrayBox& this_eta  = m_etaRes[dit];
         const FArrayBox& this_Qie  = m_Qie[dit];
               FArrayBox& this_JaSenei  = m_JaSeneDeni[dit];
               FArrayBox& this_JaSenee  = m_JaSeneDene[dit];

         const Box& cellbox = grids[dit]; // no ghosts
         FORT_EVAL_ENERGY_SOURCE_TERMS( CHF_BOX(cellbox),
                                        CHF_CONST_INT(useFiniteIon),
                                        CHF_CONST_REAL(m_lambdai0),
                                        CHF_CONST_FRA1(this_Ja,0),
                                        CHF_CONST_FRA1(this_eta,0),
                                        CHF_CONST_FRA1(this_Qie,0),
                                        CHF_CONST_FRA1(this_Nh,0),
                                        CHF_CONST_FRA1(this_Ne,0),
                                        CHF_CONST_FRA(this_JaJxB),
                                        CHF_CONST_FRA(this_gradPe),
                                        CHF_CONST_FRA(this_gradTeeff),
                                        CHF_CONST_FRA(this_NU),
                                        CHF_CONST_FRA(this_Ecc),
                                        CHF_CONST_FRA(this_Jcc),
                                        CHF_FRA1(this_JaSenei,0),
                                        CHF_FRA1(this_JaSenee,0) );

      }

   } else {

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_JaSeneDeni[dit].setVal(0.0);
         m_JaSeneDene[dit].setVal(0.0);
      }

   }
   
   if(a_soln_fluid.m_evolve_electronDensity) {

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const FArrayBox& this_Ja         = m_Jacobian[dit];
         const FArrayBox& this_Ne         = m_eleDen_cc[dit];
         const FArrayBox& this_nue_izn    = m_nue_izn[dit];
               FArrayBox& this_JaSeleDen  = m_JaSeleDen[dit];
               FArrayBox& this_JaSeneDene = m_JaSeneDene[dit];

         const Box& cellbox = grids[dit]; // no ghosts
         FORT_ADD_IONIZATION_SOURCE_TERMS( CHF_BOX(cellbox),
                                           CHF_CONST_REAL(m_Uizn_eV),
                                           CHF_CONST_REAL(m_Tscale_eV),
                                           CHF_CONST_FRA1(this_Ja,0),
                                           CHF_CONST_FRA1(this_Ne,0),
                                           CHF_CONST_FRA1(this_nue_izn,0),
                                           CHF_FRA1(this_JaSeleDen,0),
                                           CHF_FRA1(this_JaSeneDene,0) );
      }
   } else {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_JaSeleDen[dit].setVal(0.0);
      }
   }

}   

void ExtendedMhdOp::setFaceCenteredFluxes( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("ExtendedMhdOp::setFaceCenteredFluxes()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // Get mapped flux normals on each face multiplied by face area (NTFlux*faceArea)
   LevelData<FArrayBox>& thisJf = m_dummyFArray_oneComp;
   LevelData<FArrayBox>& thisFlux_mapped = m_dummyFArray_spaceDim;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      thisJf[dit].copy(m_rhoDen_cc[dit]);
      thisJf[dit].mult(m_Jacobian[dit],0,0,1);
   }
   m_geometry.multiplyNTranspose(thisFlux_mapped, m_rhoFlux_cc); // convert phys flux to mapped at cc
   computeNTFfaceArea(m_rhoFlux_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);

   if(a_soln_fluid.m_evolve_electronDensity) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         thisJf[dit].copy(m_eleDen_cc[dit]);
         thisJf[dit].mult(m_Jacobian[dit],0,0,1);
      }
      m_geometry.multiplyNTranspose(thisFlux_mapped, m_eleFlux_cc);
      computeNTFfaceArea(m_eleFlux_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   }

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
      thisJf[dit].copy(m_eneDeni_cc[dit]);
      thisJf[dit].mult(m_Jacobian[dit],0,0,1);
   }
   m_geometry.multiplyNTranspose(thisFlux_mapped, m_enFluxi_cc);
   computeNTFfaceArea(m_enFluxi_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      thisJf[dit].copy(m_eneDene_cc[dit]);
      thisJf[dit].mult(m_Jacobian[dit],0,0,1);
   }
   m_geometry.multiplyNTranspose(thisFlux_mapped, m_enFluxe_cc);
   computeNTFfaceArea(m_enFluxe_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
 
   if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         thisJf[dit].copy(m_momDen_virtual_cc[dit]);
         thisJf[dit].mult(m_Jacobian[dit],0,0,1);
      }
      m_geometry.multiplyNTranspose(thisFlux_mapped, m_mvFlux_cc);
      computeNTFfaceArea(m_mvFlux_norm,thisFlux_mapped,thisJf,m_Cspeed_cc,m_CspeedR_norm,m_CspeedL_norm);
   }

   // rebase dummy containers
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_oneComp[dit].setVal(1.0/0.0);
      m_dummyFArray_spaceDim[dit].setVal(1.0/0.0);
   }
   
}

void ExtendedMhdOp::setFaceCenteredFlux( LevelData<FluxBox>&    a_Flux_norm,
                                   const LevelData<FArrayBox>&  a_Flux_phys_cc ) const
{
   CH_TIME("IdealMhdOp::setFaceCenteredFlux()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   const int method = 2;

   // original way
   if(method==0) {
      SpaceUtils::upWindToFaces( m_dummyFlux_spaceDim, a_Flux_phys_cc, m_dummyFlux_spaceDim, "c2" );
      m_geometry.applyAxisymmetricCorrection( m_dummyFlux_spaceDim );
      m_geometry.computeMetricTermProductAverage( a_Flux_norm, m_dummyFlux_spaceDim, 0 );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_dummyFlux_spaceDim[dit].setVal(1./0.);
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
         m_dummyFArray_spaceDim[dit].setVal(1./0.); // reset dummy to infinty for safety
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
         m_dummyFArray_spaceDim[dit].setVal(1./0.); // reset dummy for safety
      }
   }

}

void ExtendedMhdOp::enforceFluxBCs( const FluidSpecies&  a_soln_fluid, 
                                    const Real           a_time )
{
   CH_TIME("ExtendedMhdOp::enforceFluxBCs()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   LevelData<FluxBox>& FluxBC_norm = m_dummyFluxBC;
   int this_cbc;

   // enforce flux BC for mass density flux
   setFaceCenteredFlux(FluxBC_norm, m_rhoFlux_cc);
   m_fluid_bc.at(0)->setFluxBC( a_soln_fluid, m_rhoFlux_norm, FluxBC_norm, a_time );

   // enforce flux BC for electron density flux
   if(a_soln_fluid.m_evolve_electronDensity) {
      setFaceCenteredFlux(FluxBC_norm, m_eleFlux_cc);
      this_cbc = m_electronDensity_cbc;
      m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_eleFlux_norm, FluxBC_norm, a_time );
   }

   // enforce flux BC for momentum density flux
   for(int n=0; n<SpaceDim; ++n) {
      if(n==0) setFaceCenteredFlux(FluxBC_norm, m_mxFlux_cc);
      if(n==1) setFaceCenteredFlux(FluxBC_norm, m_myFlux_cc);
      if(n==2) setFaceCenteredFlux(FluxBC_norm, m_mzFlux_cc);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_momFluxBC_norm[dit].copy( FluxBC_norm[dit],0,n,1 );
      }
   }
   this_cbc = m_momentumDensity_cbc;
   m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_momFlux_norm, m_momFluxBC_norm, a_time );

   // enforce flux BC for energy density flux
   setFaceCenteredFlux(FluxBC_norm, m_enFluxi_cc);
   this_cbc = m_energyDensity_cbc;
   m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_enFluxi_norm, FluxBC_norm, a_time );

   // enforce flux BC for electron energy density flux
   setFaceCenteredFlux(FluxBC_norm, m_enFluxe_cc);
   this_cbc = m_energyDensity_ele_cbc;
   m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_enFluxe_norm, FluxBC_norm, a_time );

   // enforce flux BC for virtual momentum density flux
   if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
      setFaceCenteredFlux(FluxBC_norm, m_mvFlux_cc);
      this_cbc = m_momentumDensity_virtual_cbc;
      m_fluid_bc.at(this_cbc)->setFluxBC( a_soln_fluid, m_mvFlux_norm, FluxBC_norm, a_time );
   }

   // rebase dummy container
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      FluxBC_norm[dit].setVal( 1./0. );
   }
 
}   

void ExtendedMhdOp::updateCollisionTerms( const FluidSpecies&          a_soln_fluid,
                                          const LevelData<FArrayBox>&  a_Pele, 
                                          const LevelData<FArrayBox>&  a_Pion )
{
   CH_TIME("ExtendedMhdOp::updateCollisionTerms()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   const LevelData<FArrayBox>& soln_rhoDen( a_soln_fluid.cell_var(0) );
   const double gamma = a_soln_fluid.m_gamma;
  
   if(a_soln_fluid.m_evolve_electronDensity) {

      // set neutral interaction rates
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const LevelData<FArrayBox>& soln_eleDen( a_soln_fluid.cell_var("electronDensity") );
         const FArrayBox& this_Nh      = soln_rhoDen[dit];
         const FArrayBox& this_Ne      = soln_eleDen[dit];
         const FArrayBox& this_Pe      = a_Pele[dit];
               FArrayBox& this_nue_izn = m_nue_izn[dit];
               FArrayBox& this_nue_neu = m_nue_neu[dit];

         const Box& cellbox = this_nue_izn.box(); // has ghosts
         FORT_EVAL_NEUTRAL_RATES( CHF_BOX(cellbox),
                                  CHF_CONST_REAL(m_Nscale_SI),
                                  CHF_CONST_REAL(m_tscale_s),
                                  CHF_CONST_REAL(m_Uizn_eV),
                                  CHF_CONST_REAL(m_Tscale_eV),
                                  CHF_CONST_FRA1(this_Nh,0),
                                  CHF_CONST_FRA1(this_Ne,0),
                                  CHF_CONST_FRA1(this_Pe,0),
                                  CHF_FRA1(this_nue_izn,0),
                                  CHF_FRA1(this_nue_neu,0) );
      }

   }

   // update collision terms
   //
   // m_taui = m_tauiscale1eV_s*10.0/Clog*Ti_eV^1.5/Nh/m_tscale_s (code norm) (Nh?)
   // m_taue = m_tauescale1eV_s*10.0/Clog*Te_eV^1.5/Ne/m_tscale_s (code norm)
   // m_etaRes = m_Le0sq/Ne*(1.0/taue + 1.0/tauVac + 1.0/tauAnom) (code norm)
   // m_Qie = 2.0/(gamma0-1.0)*m_mM*(Ne/taue)*(Te-Ti) (code norm)
   // m_xion = m_wciscale_Hz*a_B*m_tscale_s*m_taui (dimensionless)
   // m_xele = m_wcescale_Hz*a_B*m_tscale_s*m_taue (dimensionless)
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& this_B    = m_magField_virtual_cc[dit];
      const FArrayBox& this_Nh   = soln_rhoDen[dit];
      const FArrayBox& this_Ne   = m_eleDen_cc[dit];
      const FArrayBox& this_nue_neu = m_nue_neu[dit];
      const FArrayBox& this_Pele = a_Pele[dit];
      const FArrayBox& this_Pion = a_Pion[dit];
            FArrayBox& this_taui = m_taui[dit];
            FArrayBox& this_taue = m_taue[dit];
            FArrayBox& this_eta  = m_etaRes[dit];
            FArrayBox& this_Qie  = m_Qie[dit];
            FArrayBox& this_xi   = m_xion[dit];
            FArrayBox& this_xe   = m_xele[dit];
            
      const Box& cellbox = this_eta.box(); // has ghosts
      FORT_EVAL_COLLISION_TERMS( CHF_BOX(cellbox),
                                 CHF_CONST_REAL(m_tscale_s),
                                 CHF_CONST_REAL(m_Tscale_eV),
                                 CHF_CONST_REAL(m_tauiscale1eV_s),
                                 CHF_CONST_REAL(m_tauescale1eV_s),
                                 CHF_CONST_REAL(m_Le0sq),
                                 CHF_CONST_REAL(m_wciscale_Hz),
                                 CHF_CONST_REAL(m_wcescale_Hz),
                                 CHF_CONST_REAL(m_mM),
                                 CHF_CONST_REAL(gamma),
                                 CHF_CONST_REAL(m_etaResMin),
                                 CHF_CONST_REAL(m_NvacC),
                                 CHF_CONST_REAL(m_NvacP),
                                 CHF_CONST_FRA1(this_B,0),
                                 CHF_CONST_FRA1(this_Nh,0),
                                 CHF_CONST_FRA1(this_Ne,0),
                                 CHF_CONST_FRA1(this_nue_neu,0),
                                 CHF_CONST_FRA1(this_Pele,0),
                                 CHF_CONST_FRA1(this_Pion,0),
                                 CHF_FRA1(this_taui,0),
                                 CHF_FRA1(this_taue,0),
                                 CHF_FRA1(this_eta,0),
                                 CHF_FRA1(this_Qie,0),
                                 CHF_FRA1(this_xi,0),
                                 CHF_FRA1(this_xe,0) );

      if(m_By_insulator_bc) {
         FArrayBox& this_etaBC  = m_etaRes_BC[dit];
         FORT_SET_MAX_VAL( CHF_BOX(cellbox),
                           CHF_CONST_FRA1(this_eta,0),
                           CHF_FRA1(this_etaBC,0) );
      } 
 
   }

   // add anomalous term to resistivity
   if(a_soln_fluid.m_evolve_currentDensity && m_etaAnomFactor>0.0) {
      const LevelData<EdgeDataBox>& soln_Jcovar( a_soln_fluid.edge_var("currentDensity") );
      SpaceUtils::interpCellToEdges(m_eleDen_ce, m_eleDen_cc, m_eleDen_ce,"c2");
      computeHallVelocity( m_Vhall_cc, soln_Jcovar, m_eleDen_ce ); // phys Ve-Vi at cell center

      for (DataIterator dit(grids); dit.ok(); ++dit) {

         const FArrayBox& this_Vhall = m_Vhall_cc[dit];
         const FArrayBox& this_B     = m_magField_virtual_cc[dit];
         const FArrayBox& this_Nh    = soln_rhoDen[dit];
         const FArrayBox& this_Ne    = m_eleDen_cc[dit];
         const FArrayBox& this_Pele  = a_Pele[dit];
         const FArrayBox& this_Pion  = a_Pion[dit];
         const FArrayBox& this_xe    = m_xele[dit];
         const FArrayBox& this_taue  = m_taue[dit];
               FArrayBox& this_eta   = m_etaRes[dit];

         const Box& cellbox = this_eta.box(); // has ghosts
         FORT_ADD_ANOMALOUS_RESISTIVITY( CHF_BOX(cellbox),
                                         CHF_CONST_REAL(m_etaAnomFactor),
                                         CHF_CONST_REAL(m_tscale_s),
                                         CHF_CONST_REAL(gamma),
                                         CHF_CONST_REAL(m_wpiscale_Hz),
                                         CHF_CONST_REAL(m_wciscale_Hz),
                                         CHF_CONST_REAL(m_wcescale_Hz),
                                         CHF_CONST_REAL(m_Le0sq),
                                         CHF_CONST_FRA(this_Vhall),
                                         CHF_CONST_FRA1(this_B,0),
                                         CHF_CONST_FRA1(this_Nh,0),
                                         CHF_CONST_FRA1(this_Ne,0),
                                         CHF_CONST_FRA1(this_Pele,0),
                                         CHF_CONST_FRA1(this_Pion,0),
                                         CHF_CONST_FRA1(this_xe,0),
                                         CHF_CONST_FRA1(this_taue,0),
                                         CHF_FRA1(this_eta,0) );

         if(m_By_insulator_bc) {
            FArrayBox& this_etaBC  = m_etaRes_BC[dit];
            FORT_SET_MAX_VAL( CHF_BOX(cellbox),
                              CHF_CONST_FRA1(this_eta,0),
                              CHF_FRA1(this_etaBC,0) );
         }

      }
      m_etaRes.exchange();

   }

}

void ExtendedMhdOp::setHeatFluxCoefficients( const double                 a_gamma,
                                             const LevelData<FArrayBox>&  a_taue, 
                                             const LevelData<FArrayBox>&  a_taui, 
                                             const LevelData<FArrayBox>&  a_xele, 
                                             const LevelData<FArrayBox>&  a_xion, 
                                             const LevelData<FArrayBox>&  a_Pele, 
                                             const LevelData<FArrayBox>&  a_Pion )
{
   CH_TIME("ExtendedMhdOp::setHeatFluxCoefficients()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& taui_onPatch = a_taui[dit];
      const FArrayBox& xion_onPatch = a_xion[dit];
      const FArrayBox& Pion_onPatch = a_Pion[dit];
            FArrayBox& kappai_onPatch = m_kappai[dit];
            FArrayBox& kappai_wedge_onPatch = m_kappai_wedge[dit];
      
      const Box& cellbox = kappai_onPatch.box();
      FORT_EVAL_ION_HEAT_FLUX_COEFFS( CHF_BOX(cellbox),
                                      CHF_CONST_REAL(m_xionFactor),
                                      CHF_CONST_REAL(a_gamma),
                                      CHF_CONST_REAL(m_Vi0sq),
                                      CHF_CONST_FRA1(taui_onPatch,0),
                                      CHF_CONST_FRA1(xion_onPatch,0),
                                      CHF_CONST_FRA1(Pion_onPatch,0),
                                      CHF_FRA1(kappai_onPatch,0),
                                      CHF_FRA1(kappai_wedge_onPatch,0) );
      //
      //
      //
      const FArrayBox& taue_onPatch = a_taue[dit];
      const FArrayBox& xele_onPatch = a_xele[dit];
      const FArrayBox& Pele_onPatch = a_Pele[dit];
            FArrayBox& kappae_onPatch = m_kappae[dit];
            FArrayBox& kappae_wedge_onPatch = m_kappae_wedge[dit];
      
      FORT_EVAL_ELE_HEAT_FLUX_COEFFS( CHF_BOX(cellbox),
                                      CHF_CONST_REAL(m_xeleFactor),
                                      CHF_CONST_REAL(a_gamma),
                                      CHF_CONST_REAL(m_Ve0sq),
                                      CHF_CONST_FRA1(taue_onPatch,0),
                                      CHF_CONST_FRA1(xele_onPatch,0),
                                      CHF_CONST_FRA1(Pele_onPatch,0),
                                      CHF_FRA1(kappae_onPatch,0),
                                      CHF_FRA1(kappae_wedge_onPatch,0) );

   }
   
   SpaceUtils::upWindToFaces(m_kappai_cf,m_kappai,m_kappai_cf,"c2"); 
   SpaceUtils::upWindToFaces(m_kappai_wedge_cf,m_kappai_wedge,m_kappai_wedge_cf,"c2"); 
   SpaceUtils::upWindToFaces(m_kappae_cf,m_kappae,m_kappae_cf,"c2"); 
   SpaceUtils::upWindToFaces(m_kappae_wedge_cf,m_kappae_wedge,m_kappae_wedge_cf,"c2"); 

}

void ExtendedMhdOp::computeCollHeatFlux( LevelData<FluxBox>&    a_q,
                                         LevelData<FluxBox>&    a_Dq,
                                   const LevelData<FArrayBox>&  a_Temp_cc,
                                   const LevelData<FluxBox>&    a_Temp_cf,
                                   const LevelData<FluxBox>&    a_Dens_cf,
                                   const LevelData<FluxBox>&    a_kappa_cf,
                                   const int                    a_charge,
                                   const double                 a_tauMin,
                                   const double                 a_limiter,
                                   const double                 a_V0sq,
                                   const double                 a_gamma )
{
   CH_TIME("ExtendedMhdOp::computeCollHeatFlux()");

   // compute physical collisional heat fluxes on cell faces
   //
   // q = -kappa*grad(T), kappa is flux limited
   //

   double tauMinEff; 
   if(a_charge==-1) tauMinEff = a_tauMin*11.92/3.77;
   if(a_charge==1) tauMinEff = a_tauMin*2.645/0.677;
   
   m_geometry.computeMappedGradient(a_Temp_cc, m_dummyFlux_spaceDim, 2);   
   m_geometry.convertPhysToCovar(m_dummyFlux_spaceDim, 1);

   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
         
      for (int dir=0; dir<SpaceDim; ++dir) {
         
         const FArrayBox& Temp_on_dir  = a_Temp_cf[dit][dir];
         const FArrayBox& Dens_on_dir  = a_Dens_cf[dit][dir];
         const FArrayBox& kappa_on_dir = a_kappa_cf[dit][dir];
         const FArrayBox& gradT_on_dir = m_dummyFlux_spaceDim[dit][dir];
               FArrayBox& Dq_on_dir    = a_Dq[dit][dir];
               FArrayBox& q_on_dir     = a_q[dit][dir];
      
         Box facebox( grids[dit] ); // no ghost 
         facebox.surroundingNodes( dir );
         FORT_EVAL_COLL_HEAT_FLUX( CHF_BOX(facebox),
                                   CHF_CONST_REAL(a_gamma),
                                   CHF_CONST_REAL(tauMinEff),
                                   CHF_CONST_REAL(m_kappaMin),
                                   CHF_CONST_REAL(a_limiter),
                                   CHF_CONST_REAL(a_V0sq),
                                   CHF_CONST_FRA1(Temp_on_dir,0),
                                   CHF_CONST_FRA1(Dens_on_dir,0),
                                   CHF_CONST_FRA1(kappa_on_dir,0),
                                   CHF_CONST_FRA(gradT_on_dir),
                                   CHF_FRA(Dq_on_dir),
                                   CHF_FRA(q_on_dir) );

      }
   
   }
   
}

void ExtendedMhdOp::setIonViscosityCoefficients( const LevelData<FArrayBox>&  a_taui, 
                                                 const LevelData<FArrayBox>&  a_xion, 
                                                 const LevelData<FArrayBox>&  a_Pion )
{
   CH_TIME("ExtendedMhdOp::setIonViscosityCoefficients()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& taui_onPatch = a_taui[dit];
      const FArrayBox& xion_onPatch = a_xion[dit];
      const FArrayBox& Pion_onPatch = a_Pion[dit];
            FArrayBox& etaVisC_onPatch = m_etaVisC[dit];
            FArrayBox& etaVisG_onPatch = m_etaVisG[dit];
      
      const Box& cellbox = etaVisC_onPatch.box();
      FORT_EVAL_ION_VISC_COEFFS( CHF_BOX(cellbox),
                                 CHF_CONST_REAL(m_xionFactor),
                                 CHF_CONST_REAL(m_etaMin),
                                 CHF_CONST_FRA1(taui_onPatch,0),
                                 CHF_CONST_FRA1(xion_onPatch,0),
                                 CHF_CONST_FRA1(Pion_onPatch,0),
                                 CHF_FRA(etaVisC_onPatch),
                                 CHF_FRA(etaVisG_onPatch) );

   }
   
   SpaceUtils::upWindToFaces(m_etaVisC_cf,m_etaVisC,m_etaVisC_cf,"c2"); 
   SpaceUtils::upWindToFaces(m_etaVisG_cf,m_etaVisG,m_etaVisG_cf,"c2"); 

}

void ExtendedMhdOp::computeViscFluxes( LevelData<FluxBox>&    a_m0JaFluxVisc_cf,
                                       LevelData<FluxBox>&    a_m1JaFluxVisc_cf,
                                       LevelData<FArrayBox>&  a_momVisc_source,
                                 const LevelData<FluxBox>&    a_velocity_phys_cf,
                                 const LevelData<FArrayBox>&  a_velocity_phys,
                                 const LevelData<FluxBox>&    a_etaVisC_cf,
                                 const LevelData<FArrayBox>&  a_etaVisC )
{
   CH_TIME("ExtendedMhdOp::computeViscFluxes()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // compute physical viscosity fluxes on cell faces
   // for each momentum equation (multiplied by Jacobian)
   // where viscosity tensor is Pi = -etaVis*W with 
   // rate of strain tensor W defined as
   // W = \nabla(V) + (\nabla(V))^T - 2/3*I\nabla\cdot V
   //
   m_geometry.computeJaStrainTensorPhys( m_JaW_cf,
                                         a_velocity_phys );
 
   
   // calculate Ja*Pion for viscosity flux limiting
   for (DataIterator dit(grids); dit.ok(); ++dit) {
        m_JaPion[dit].copy( m_pressure_ion[dit] ); 
        m_JaPion[dit].mult( m_Jacobian[dit] ); 
   }
   SpaceUtils::upWindToFaces(m_JaPion_cf,m_JaPion,m_JaPion_cf,"c2"); 
  
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      for (int dir=0; dir<SpaceDim; ++dir) {
         
         const FArrayBox& etaVisC_on_dir = a_etaVisC_cf[dit][dir];
         const FArrayBox& JaPion_on_dir  = m_JaPion_cf[dit][dir];
         const FArrayBox& Ja_on_dir    = m_Jacobian_cf[dit][dir];
         const FArrayBox& JaW_on_dir   = m_JaW_cf[dit][dir];
               FArrayBox& Dx_on_dir    = m_DiffVx_cf[dit][dir];
               FArrayBox& Dz_on_dir    = m_DiffVz_cf[dit][dir];
               FArrayBox& JaPi0_on_dir = a_m0JaFluxVisc_cf[dit][dir];
               FArrayBox& JaPi1_on_dir = a_m1JaFluxVisc_cf[dit][dir];
      
         Box facebox( grids[dit] ); // no ghost 
         facebox.surroundingNodes( dir );
         FORT_EVAL_COLL_VISC_FLUX( CHF_BOX(facebox),
                                   CHF_CONST_REAL(m_tauiMin),
                                   CHF_CONST_REAL(m_etaMin),
                                   CHF_CONST_REAL(m_viscLimiter),
                                   CHF_CONST_FRA(etaVisC_on_dir),
                                   CHF_CONST_FRA1(JaPion_on_dir,0),
                                   CHF_CONST_FRA1(Ja_on_dir,0),
                                   CHF_CONST_FRA(JaW_on_dir),
                                   CHF_FRA(Dx_on_dir),
                                   CHF_FRA(Dz_on_dir),
                                   CHF_FRA(JaPi0_on_dir),
                                   CHF_FRA(JaPi1_on_dir) );

      }
   
   }
   
   // calculate viscosity flux source term for axisymmetric
   if(m_twoDaxisymm) {
      computeViscSourceAxisymm( a_momVisc_source, a_velocity_phys, a_velocity_phys_cf );
      //for (DataIterator dit(grids); dit.ok(); ++dit) {
      //   a_momVisc_source[dit].mult( a_etaVisC[dit],0,0,1 ); 
      //   a_momVisc_source[dit].negate(); 
      //}
      for (DataIterator dit(grids); dit.ok(); ++dit) {
      
         const FArrayBox& etaVisC_onPatch = a_etaVisC[dit];
         const FArrayBox& Pion_onPatch  = m_pressure_ion[dit];
               FArrayBox& Dth_onPatch   = m_DiffVth[dit];
               FArrayBox& Piyy_onPatch  = a_momVisc_source[dit];
      
         const Box& cellbox = Piyy_onPatch.box();
         FORT_EVAL_COLL_VISC_SOURCE( CHF_BOX(cellbox),
                                     CHF_CONST_REAL(m_tauiMin),
                                     CHF_CONST_REAL(m_etaMin),
                                     CHF_CONST_REAL(m_viscLimiter),
                                     CHF_CONST_FRA1(Pion_onPatch,0),
                                     CHF_CONST_FRA1(etaVisC_onPatch,0),
                                     CHF_FRA1(Dth_onPatch,0),
                                     CHF_FRA1(Piyy_onPatch,0) );
      } 
   }

}

void ExtendedMhdOp::computeViscSourceAxisymm( LevelData<FArrayBox>&  a_Wthth,
                                        const LevelData<FArrayBox>&  a_velocity_phys,
                                        const LevelData<FluxBox>&    a_velocity_phys_cf )
{
   CH_TIME("ExtendedMhdOp::computeViscSourceAxisymm()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // calculate div stress tensor source terms at cell center
   // Wthth = 2*Vr/r - 2/3*divV
   // Wthr  = d(Vth)/dr + d(Vr)/dtheta/r
   
   //   compute divergence of physical velocity field
   LevelData<FluxBox> velocity_norm(grids, 1, IntVect::Unit);
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

void ExtendedMhdOp::computeViscEnergyFlux( LevelData<FluxBox>&  a_UJaPi_Flux_phys,
                                     const LevelData<FluxBox>&  a_JaPix_Flux_phys,
                                     const LevelData<FluxBox>&  a_JaPiy_Flux_phys,
                                     const LevelData<FluxBox>&  a_velocity_phys_cf )
{
   CH_TIME("ExtendedMhdOp::computeViscEnergyFlux()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // calculate physical viscosity fluxes for total energy equation
   // on each face (multiplied by Jacobian)
   //
   // Ja*Flux_0 = U0*JaPi_00 + U1*JaPi_10
   // Ja*Flux_1 = U0*JaPi_01 + U1*JaPi_11

   for (DataIterator dit(grids); dit.ok(); ++dit) {

      for (int dir=0; dir<SpaceDim; ++dir) {

         const FArrayBox& vel_on_dir   = a_velocity_phys_cf[dit][dir];
         const FArrayBox& JaPix_on_dir = a_JaPix_Flux_phys[dit][dir];
         const FArrayBox& JaPiy_on_dir = a_JaPiy_Flux_phys[dit][dir];
               FArrayBox& UJaPi_on_dir = a_UJaPi_Flux_phys[dit][dir];

         Box facebox( grids[dit] ); // no ghost 
         facebox.surroundingNodes( dir );
         FORT_EVAL_VISC_ENERGY_FLUX( CHF_BOX(facebox),
                                     CHF_CONST_FRA(vel_on_dir),
                                     CHF_CONST_FRA(JaPix_on_dir),
                                     CHF_CONST_FRA(JaPiy_on_dir),
                                     CHF_FRA(UJaPi_on_dir) );

      }

   }

}

void ExtendedMhdOp::computeEleViscEnergyFlux( LevelData<FluxBox>&        a_UePie_Flux_phys,
                                        const LevelData<FArrayBox>&      a_Pie_xz_cc,
                                        const LevelData<NodeFArrayBox>&  a_Pie_diag_cn,
                                        const LevelData<FluxBox>&        a_velocity_ele_cf )
{
   CH_TIME("ExtendedMhdOp::computeEleViscEnergyFlux()");
   CH_assert( a_Pie_xz_cc.ghostVect() >= 1*IntVect::Unit );
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // calculate physical electron viscosity fluxes for electron 
   // energy equation on each face.
   //
   // Note that ele visc fluxes live at nodes and cells, but energy
   // flux lives at faces. Interpolations all done here inside
   // call to FORTRAN to avoid defining extra containers
   //
   // Flux_0 = Ue0*Pie_00 + Ue1*Pie_10
   // Flux_1 = Ue0*Pie_01 + Ue1*Pie_11
   //

   for (DataIterator dit(grids); dit.ok(); ++dit) {
         
      const FArrayBox& Pie_xz_on_dir = a_Pie_xz_cc[dit];

      for (int dir=0; dir<SpaceDim; ++dir) {

         const FArrayBox& Ue_on_dir       = a_velocity_ele_cf[dit][dir];
         const FArrayBox& Pie_diag_on_dir = a_Pie_diag_cn[dit].getFab();
               FArrayBox& UePie_on_dir    = a_UePie_Flux_phys[dit][dir];

         Box facebox( grids[dit] ); // no ghost 
         facebox.surroundingNodes( dir );
         FORT_EVAL_ELE_VISC_ENERGY_FLUX( CHF_BOX(facebox),
                                         CHF_CONST_INT(dir),
                                         CHF_CONST_FRA(Ue_on_dir),
                                         CHF_CONST_FRA1(Pie_xz_on_dir,0),
                                         CHF_CONST_FRA(Pie_diag_on_dir),
                                         CHF_FRA(UePie_on_dir) );

      }

   }

}

void ExtendedMhdOp::addGyroHeatToSource( LevelData<FArrayBox>&  a_RHS,
                                   const LevelData<FArrayBox>&  a_Temp_cc,
                                   const LevelData<FArrayBox>&  a_Dens_cc,
                                   const LevelData<FArrayBox>&  a_kappaG_cc,
                                   const double                 a_charge,
                                   const double                 a_limiter,
                                   const double                 a_V0sq )
{
   CH_TIME("ExtendedMhdOp::addGyroHeatToSource()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // update the RHS += -JaDivqG to include divergence of
   // the diamagnetic heat flux
   // 
   // qG = +kappaG*bxgrad(T) (for ions)
   // qG = -kappaG*bxgrad(T) (for eles)
   //
   // This is non-conservative formulation that directly
   // exploits Div(Curl(T)) = 0
   //

   // compute physical gradT
   //
   m_geometry.computeMappedPoloidalGradientWithGhosts(a_Temp_cc, m_gradT_cc, 2);   
   m_geometry.convertPhysToCovar(m_gradT_cc, 1);

   // multiply kappaG by charge sign and Jacobian, then compute physical Ja*curl at cell center
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_dummyFArray_oneComp[dit].copy(a_kappaG_cc[dit]);   
      m_dummyFArray_oneComp[dit].mult(a_charge);   
      if(m_twoDaxisymm) {
         m_dummyFArray_oneComp[dit].mult(m_g_y[dit]);
      }
   }
   m_geometry.mappedGridCurlofVirtComp(m_dummyFArray_spaceDim,m_dummyFArray_oneComp); // Ja*curl(A)\cdot g^l
   m_geometry.convertPhysToContravar(m_dummyFArray_spaceDim, 1);
   

   // updated RHS to included -JaDivqG
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {

      const FArrayBox& gradA_onPatch = m_dummyFArray_spaceDim[dit];
      const FArrayBox& gradT_onPatch = m_gradT_cc[dit];
            FArrayBox& RHS_onPatch   = a_RHS[dit];
      
      Box cellbox( grids[dit] ); // no ghost 
      FORT_ADD_GYRO_HEAT_SOURCE( CHF_BOX(cellbox),
                                 CHF_CONST_FRA(gradA_onPatch),
                                 CHF_CONST_FRA(gradT_onPatch),
                                 CHF_FRA1(RHS_onPatch,0) );

   }
   //SpaceUtils::inspectFArrayBox(m_gradT_cc,1);
   //SpaceUtils::inspectFArrayBox(a_kappaG_cc,0);
   //SpaceUtils::inspectFArrayBox(m_dummyFArray_spaceDim,0,0);

}

void ExtendedMhdOp::computeGyroHeatFluxPhys( LevelData<FArrayBox>&  a_q_wedge,
                                       const LevelData<FArrayBox>&  a_Temp,
                                       const LevelData<FArrayBox>&  a_kappa_wedge,
                                       const double                 a_charge )
{
   CH_TIME("ExtendedMhdOp::computeGyroHeatFluxPhys()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // compute physical gyro heat flux at cell center needed for
   // generalized gyro viscosity calculation
   // 
   // q_wedge = +kappa_wedge*bxgrad(T) (for ions)
   // q_wedge = -kappa_wedge*bxgrad(T) (for eles)
   //

   m_geometry.computeMappedPoloidalGradientWithGhosts(a_Temp, m_gradT_cc, 2);   
   // This fills two codim 1 ghost cell layers at all block boundaries
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = block_coord_sys.domain();
      fourthOrderCellExtrapAtDomainBdry(m_gradT_cc[dit], domain, grids[dit]);
      //secondOrderCellExtrapAtDomainBdry(m_gradT_cc[dit], domain, grids[dit]);
   }
   m_geometry.convertPhysToCovar(m_gradT_cc, 1);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
         
         const FArrayBox& thiskappa = a_kappa_wedge[dit];
         const FArrayBox& thisgradT = m_gradT_cc[dit];
               FArrayBox& thisq     = a_q_wedge[dit];
      
         const Box& cellbox = thisq.box(); 
         FORT_COMPUTE_GYRO_HEAT( CHF_BOX(cellbox),
                                 CHF_CONST_REAL(a_charge),
                                 CHF_CONST_FRA1(thiskappa,0),
                                 CHF_CONST_FRA(thisgradT),
                                 CHF_FRA(thisq) );

   }

   // extrapolate to one-layer of ghosts and exchange
   //
   a_q_wedge.exchange();
   
}

void ExtendedMhdOp::addGyroHeatToFlux( LevelData<FluxBox>&    a_q,
                                 const LevelData<FArrayBox>&  a_Temp_cc,
                                 const LevelData<FluxBox>&    a_Temp_cf,
                                 const LevelData<FluxBox>&    a_Dens_cf,
                                 const LevelData<FluxBox>&    a_kappa_cf,
                                 const double                 a_charge,
                                 const double                 a_limiter,
                                 const double                 a_V0sq )
{
   CH_TIME("ExtendedMhdOp::addGyroHeatToFlux()");

   // update the collisional heat flux on cell faces
   // to include the diamagnetic heat flux
   // 
   // q += +kappa_wedge*bxgrad(T) (for ions)
   // q += -kappa_wedge*bxgrad(T) (for eles)
   //

   m_geometry.computeMappedGradient(a_Temp_cc, m_dummyFlux_spaceDim, 2);   
   m_geometry.convertPhysToCovar(m_dummyFlux_spaceDim, 1);

   const DisjointBoxLayout& grids( m_geometry.grids() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
         
      for (int dir=0; dir<SpaceDim; ++dir) {
         
         const FArrayBox& Temp_on_dir  = a_Temp_cf[dit][dir];
         const FArrayBox& Dens_on_dir  = a_Dens_cf[dit][dir];
         const FArrayBox& kappa_on_dir = a_kappa_cf[dit][dir];
         const FArrayBox& gradT_on_dir = m_dummyFlux_spaceDim[dit][dir];
               FArrayBox& q_on_dir     = a_q[dit][dir];
      
         Box facebox( grids[dit] ); // no ghost 
         facebox.surroundingNodes( dir );
         FORT_ADD_GYRO_HEAT( CHF_BOX(facebox),
                             CHF_CONST_REAL(a_charge),
                             CHF_CONST_REAL(a_limiter),
                             CHF_CONST_REAL(a_V0sq),
                             CHF_CONST_FRA1(Temp_on_dir,0),
                             CHF_CONST_FRA1(Dens_on_dir,0),
                             CHF_CONST_FRA1(kappa_on_dir,0),
                             CHF_CONST_FRA(gradT_on_dir),
                             CHF_FRA(q_on_dir) );

      }
   
   }
   
}


void ExtendedMhdOp::addGyroViscToFluxes( LevelData<FluxBox>&    a_m0JaFluxVisc_cf,
                                         LevelData<FluxBox>&    a_m1JaFluxVisc_cf,
                                   const LevelData<FArrayBox>&  a_velocity_phys,
                                   const LevelData<FluxBox>&    a_etaVisG_cf )
{
   CH_TIME("ExtendedMhdOp::addGyroViscToFluxes()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // update physical viscosity fluxes on faces to include
   // gyroviscosity (multiplied by Jacobian)
   //
   m_geometry.computeJaStrainTensorPhys( m_JaW_cf,
                                         a_velocity_phys );
   
   // calculate Ja*Pion for viscosity flux limiting
   //
   for (DataIterator dit(grids); dit.ok(); ++dit) {
        m_JaPion[dit].copy( m_pressure_ion[dit] ); 
        m_JaPion[dit].mult( m_Jacobian[dit] ); 
   }
   SpaceUtils::upWindToFaces(m_JaPion_cf,m_JaPion,m_JaPion_cf,"c2"); 
 
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      for (int dir=0; dir<SpaceDim; ++dir) {
         
         const FArrayBox& etaVisG_on_dir = a_etaVisG_cf[dit][dir];
         const FArrayBox& JaPion_on_dir  = m_JaPion_cf[dit][dir];
         const FArrayBox& JaW_on_dir   = m_JaW_cf[dit][dir];
               FArrayBox& JaPi0_on_dir = a_m0JaFluxVisc_cf[dit][dir];
               FArrayBox& JaPi1_on_dir = a_m1JaFluxVisc_cf[dit][dir];
      
         const Box& facebox = JaW_on_dir.box();
         FORT_ADD_GYRO_VISC( CHF_BOX(facebox),
                             CHF_CONST_REAL(m_gyroViscLimiter),
                             CHF_CONST_FRA1(JaPion_on_dir,0),
                             CHF_CONST_FRA(etaVisG_on_dir),
                             CHF_CONST_FRA(JaW_on_dir),
                             CHF_FRA(JaPi0_on_dir),
                             CHF_FRA(JaPi1_on_dir) );

      }
      
   }

}  

void ExtendedMhdOp::addGyroHeatToViscFluxes( LevelData<FluxBox>&    a_m0JaFluxVisc_cf,
                                             LevelData<FluxBox>&    a_m1JaFluxVisc_cf,
                                       const double                 a_gamma,
                                       const LevelData<FArrayBox>&  a_qwedge_phys,
                                       const LevelData<FluxBox>&    a_etaVisG_cf )
{
   CH_TIME("ExtendedMhdOp::addGyroHeatToViscToFluxes()");
   const DisjointBoxLayout& grids( m_geometry.grids() );

   // update physical viscosity fluxes on faces to include
   // gyroHeatviscosity (multiplied by Jacobian)
   //
   // stress tensor here is computed using physical gyro heat flux in replace of
   // velocity. The stress tensor is the scaled by 2/5/Pion
   //
   // Ja*Wq = Ja*(gamma-1/gamma)/P*nabla(qwedge)
   //
   m_geometry.computeJaStrainTensorPhys( m_JaWq_cf,
                                         a_qwedge_phys );
   
   SpaceUtils::upWindToFaces(m_dummyFlux_oneComp,m_pressure_ion,m_dummyFlux_oneComp,"c2"); 
   const double gmoog = (a_gamma-1.0)/a_gamma;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         const Box& thisBox = m_JaWq_cf[dit][dir].box();
         for (int n=0; n<m_JaWq_cf.nComp(); ++n) {
            m_JaWq_cf[dit][dir].divide(m_dummyFlux_oneComp[dit][dir],thisBox,0,n,1);
         }
         m_JaWq_cf[dit][dir].mult( gmoog );
         m_JaWq_cf[dit][dir].mult( m_gyroHeatViscCoeff );
      }
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      for (int dir=0; dir<SpaceDim; ++dir) {
         
         const FArrayBox& etaVisG_on_dir = a_etaVisG_cf[dit][dir];
         const FArrayBox& JaPion_on_dir  = m_JaPion_cf[dit][dir];
         const FArrayBox& JaWq_on_dir    = m_JaWq_cf[dit][dir];
               FArrayBox& JaPi0_on_dir = a_m0JaFluxVisc_cf[dit][dir];
               FArrayBox& JaPi1_on_dir = a_m1JaFluxVisc_cf[dit][dir];
      
         const Box& facebox = JaWq_on_dir.box();
         FORT_ADD_GYRO_VISC( CHF_BOX(facebox),
                             CHF_CONST_REAL(m_gyroViscLimiter),
                             CHF_CONST_FRA1(JaPion_on_dir,0),
                             CHF_CONST_FRA(etaVisG_on_dir),
                             CHF_CONST_FRA(JaWq_on_dir),
                             CHF_FRA(JaPi0_on_dir),
                             CHF_FRA(JaPi1_on_dir) );

      }
      
   }

}  

void ExtendedMhdOp::updateRHSs_heat( FluidSpecies&  a_rhs_fluid,
                               const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("ExtendedMhdOp::updateRHSs_heat()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
         
   // add diamagnetic heat flux
   //
   /*
   if(m_gyroHeatLimiter>0.0) {
      
      // calculate temperature for each species at cell center and 
      // cell faces
      //
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_Tele_cc[dit].copy( m_pressure_ele[dit] );
         m_Tele_cc[dit].divide( m_eleDen_cc[dit] );
         //
         m_Tion_cc[dit].copy( m_pressure_ion[dit] );
         m_Tion_cc[dit].divide( m_rhoDen_cc[dit] );
      }
      SpaceUtils::upWindToFaces(m_Tele_cf,m_Tele_cc,m_Tele_cf,"c2");
      SpaceUtils::upWindToFaces(m_Tion_cf,m_Tion_cc,m_Tion_cf,"c2");

      // update physical heat flux for each species to include
      // the diamagnetic heat flux
      //
      setHeatFluxCoefficients( a_soln_fluid.m_gamma, 
                               m_taue, m_taui,
                               m_xele, m_xion,
                               m_pressure_ele, m_pressure_ion );
      double thisCharge = -1.0;
      addGyroHeatToFlux( m_qe_cf, 
                         m_Tele_cc, m_Tele_cf, m_eleDen_cf, 
                         m_kappae_wedge_cf, 
                         thisCharge, m_gyroHeatLimiter, m_Ve0sq );
      thisCharge = 1.0;
      addGyroHeatToFlux( m_qi_cf, 
                         m_Tion_cc, m_Tion_cf, m_rhoDen_cf, 
                         m_kappai_wedge_cf, 
                         thisCharge, m_gyroHeatLimiter, m_Vi0sq );
   
   }
   */  
 

   //   update RHSs for ion energy density equation
   //
   m_geometry.applyAxisymmetricCorrection( m_qi_cf );
   m_geometry.computeMetricTermProductAverage(m_qi_norm, m_qi_cf, 0);
   
   LevelData<FArrayBox>& rhs_ene_ion( a_rhs_fluid.cell_var("energyDensity") );
   m_geometry.mappedGridDivergenceFromFluxNorms(m_qi_norm, m_dummyDiv);
   for (DataIterator dit(rhs_ene_ion.dataIterator()); dit.ok(); ++dit) {
      m_dummyDiv[dit].divide(m_cellVol[dit]);
      m_dummyDiv[dit].mult(m_Jacobian[dit]);
      rhs_ene_ion[dit].minus( m_dummyDiv[dit] );
   }
  
   //   update RHSs for electron energy density equation
   //
   m_geometry.applyAxisymmetricCorrection( m_qe_cf );
   m_geometry.computeMetricTermProductAverage(m_qe_norm, m_qe_cf, 0);
   
   LevelData<FArrayBox>& rhs_ene_ele( a_rhs_fluid.cell_var("energyDensity_ele") );
   m_geometry.mappedGridDivergenceFromFluxNorms(m_qe_norm, m_dummyDiv);
   for (DataIterator dit(rhs_ene_ele.dataIterator()); dit.ok(); ++dit) {
      m_dummyDiv[dit].divide(m_cellVol[dit]);
      m_dummyDiv[dit].mult(m_Jacobian[dit]);
      rhs_ene_ele[dit].minus( m_dummyDiv[dit] );
   }
   
   // add (subtract) div of diamagnetic heat flux to RHS
   //
   if(m_gyroHeatLimiter>0.0) {
      
      // calculate temperature for each species 
      //
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_Tele_cc[dit].copy( m_pressure_ele[dit] );
         m_Tele_cc[dit].divide( m_eleDen_cc[dit] );
         //
         m_Tion_cc[dit].copy( m_pressure_ion[dit] );
         m_Tion_cc[dit].divide( m_rhoDen_cc[dit] );
      }

      // update physical heat flux for each species to include
      // the diamagnetic heat flux
      //
      setHeatFluxCoefficients( a_soln_fluid.m_gamma, 
                               m_taue, m_taui,
                               m_xele, m_xion,
                               m_pressure_ele, m_pressure_ion );
      double thisCharge = -1.0;
      addGyroHeatToSource( rhs_ene_ele, 
                           m_Tele_cc, m_eleDen_cc, 
                           m_kappae_wedge, 
                           thisCharge, m_gyroHeatLimiter, m_Ve0sq );
      thisCharge = 1.0;
      addGyroHeatToSource( rhs_ene_ion, 
                           m_Tion_cc, m_rhoDen_cc, 
                           m_kappai_wedge, 
                           thisCharge, m_gyroHeatLimiter, m_Vi0sq );
   
   }
   

  
}

void ExtendedMhdOp::updateRHSs_visc( FluidSpecies&  a_rhs_fluid,
                               const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("ExtendedMhdOp::updateRHSs_visc()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // convert physical Ja*Flux to contravar and mult by faceArea
   // (Result is same as applying N^T*Flux*faceArea as is done by 
   // computeMetricTermProductAverage), which is not done here 
   // since Flux is already multiplied by Jacobian
   //

   // add gyroviscosity to vicsosity momentum flux
   //
   if(m_gyroViscLimiter>0.0) {

      updateCollisionTerms( a_soln_fluid, m_pressure_ele, m_pressure_ion ); 
      setIonViscosityCoefficients( m_taui, m_xion, m_pressure_ion );
      addGyroViscToFluxes( m_m0JaFluxVisc_cf, m_m1JaFluxVisc_cf, m_velocity, m_etaVisG_cf);

      ///////////////////////////////////////////////////////////////////////////////////
      //
      // compute gyro heat flux at cell center for generalized gyro visc
      // with heat flux included in rate of strain tensor
      //
      if(m_gyroHeatViscCoeff>0.0) {
         setHeatFluxCoefficients( a_soln_fluid.m_gamma, 
                                  m_taue, m_taui, m_xele, m_xion,
                                  m_pressure_ele, m_pressure_ion );
         double thisCharge = 1.0;
         computeGyroHeatFluxPhys( m_qi_wedge_cc, m_Tion_cc, m_kappai_wedge, thisCharge ); 

         addGyroHeatToViscFluxes( m_m0JaFluxVisc_cf, m_m1JaFluxVisc_cf, a_soln_fluid.m_gamma,
                                  m_qi_wedge_cc, m_etaVisG_cf);
         // note that this contribution is not traceless and should be included in energy
         // flux
         computeViscEnergyFlux( m_enJaFluxVisc_cf,
                                m_m0JaFluxVisc_cf,
                                m_m1JaFluxVisc_cf,
                                m_velocity_cf );
      }
      //
      //
      ////////////////////////////////////////////////////////////////////////////////////

   }

   //
   m_geometry.applyAxisymmetricCorrection( m_m0JaFluxVisc_cf );
   m_geometry.applyAxisymmetricCorrection( m_m1JaFluxVisc_cf );
   m_geometry.computedxidXProductNorm(m_m0JaFluxVisc_norm, m_m0JaFluxVisc_cf);
   m_geometry.computedxidXProductNorm(m_m1JaFluxVisc_norm, m_m1JaFluxVisc_cf);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_momJaFluxVisc_norm[dit].copy( m_m0JaFluxVisc_norm[dit],0,0,1 ); 
      m_momJaFluxVisc_norm[dit].copy( m_m1JaFluxVisc_norm[dit],0,1,1 );
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
  
 
   //   update RHS for ion energy density equation
   //
   // convert physical Ja*Flux to contravar and mult by faceArea
   // (Result is same as applying N^T*Flux*faceArea as is done by 
   // computeMetricTermProductAverage), which is not done here 
   // since Flux is already multiplied by Jacobian
   //
   m_geometry.applyAxisymmetricCorrection( m_enJaFluxVisc_cf );
   m_geometry.computedxidXProductNorm(m_enJaFluxVisc_norm, m_enJaFluxVisc_cf);

   LevelData<FArrayBox>& rhs_ene_ion( a_rhs_fluid.cell_var("energyDensity") );
   m_geometry.mappedGridDivergenceFromFluxNorms(m_enJaFluxVisc_norm, m_dummyDiv);
   for (DataIterator dit(rhs_ene_ion.dataIterator()); dit.ok(); ++dit) {
      m_dummyDiv[dit].divide(m_cellVol[dit]);
      m_dummyDiv[dit].mult(m_Jacobian[dit]);
      rhs_ene_ion[dit].minus( m_dummyDiv[dit] );
   }
   
   
   //   update RHS for electron energy density equation
   //
   /*
   if(m_etaVis_ele>0.0) {
      
      // calculate ele energy density fluxes
      //
      LevelData<EdgeDataBox>& soln_Jcovar( a_rhs_fluid.edge_var("currentDensity") );
      SpaceUtils::interpEdgesToCell(m_Jcc,soln_Jcovar,"c2"); // covar at cell center
      m_geometry.convertPhysToCovar(m_Jcc,1);                // phys at cell center
      computeEleVelOnCells(m_Ve_cc, m_Jcc, m_eleDen_cc, m_velocity);
      SpaceUtils::upWindToFaces(m_Ve_cf,m_Ve_cc,m_Ve_cf,"c2");
      computeEleViscEnergyFlux( m_enEleFluxVisc_cf,
                                m_Pie_xz_cc,
                                m_Pie_diag_cn,
                                m_Ve_cf );
 
      // compute div of flux and add to RHS
      //
      m_geometry.applyAxisymmetricCorrection( m_enEleFluxVisc_cf );
      m_geometry.computeMetricTermProductAverage(m_enEleFluxVisc_norm, m_enEleFluxVisc_cf, 0);

      LevelData<FArrayBox>& rhs_ene_ele( a_rhs_fluid.cell_var("energyDensity_ele") );
      m_geometry.mappedGridDivergenceFromFluxNorms(m_enEleFluxVisc_norm, m_dummyDiv);
      for (DataIterator dit(rhs_ene_ele.dataIterator()); dit.ok(); ++dit) {
         m_dummyDiv[dit].divide(m_cellVol[dit]);
         m_dummyDiv[dit].mult(m_Jacobian[dit]);
         rhs_ene_ele[dit].minus( m_dummyDiv[dit] );
      }
   
   }
   */
}


void ExtendedMhdOp::updateRHSs( FluidSpecies&  a_rhs_fluid,
                          const FluidSpecies&  a_soln_fluid ) 
{
   CH_TIME("ExtendedMhdOp::updateRHSs()");
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // compute divergence of Fluxes and add (subtract) to RHS 
   // Note that Ja/cellVol = 1/mapVol
   
   if(m_updateHydroVars) {  

   // update RHS for mass density
   LevelData<FArrayBox>& rhs_rho( a_rhs_fluid.cell_var(0) );
   m_geometry.mappedGridDivergenceFromFluxNorms(m_rhoFlux_norm, m_dummyDiv);
   setDivMassFluxBinary(m_divMassFluxBinary,m_dummyDiv);
   for (DataIterator dit(rhs_rho.dataIterator()); dit.ok(); ++dit) {
      m_dummyDiv[dit].divide(m_cellVol[dit]);
      m_dummyDiv[dit].mult(m_Jacobian[dit]);
      m_dummyDiv[dit].mult(m_divMassFluxBinary[dit]);
      rhs_rho[dit].minus(m_dummyDiv[dit]);
   }

   // update RHS for electron density
   if(a_soln_fluid.m_evolve_electronDensity) {
      LevelData<FArrayBox>& rhs_ele( a_rhs_fluid.cell_var("electronDensity") );
      m_geometry.mappedGridDivergenceFromFluxNorms(m_eleFlux_norm, m_dummyDiv);
      for (DataIterator dit(rhs_ele.dataIterator()); dit.ok(); ++dit) {
         m_dummyDiv[dit].divide(m_cellVol[dit]);
         m_dummyDiv[dit].mult(m_Jacobian[dit]);
         m_dummyDiv[dit].mult(m_divMassFluxBinary[dit]);
         rhs_ele[dit].minus(m_dummyDiv[dit]);
         rhs_ele[dit].plus(m_JaSeleDen[dit]);
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
            pdivIdent_R.mult(m_pressure[dit],0,0,1);
            pdivIdent_Z.mult(m_pressure[dit],0,0,1);
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
      if(a_soln_fluid.m_evolve_magneticField_virtual) {
         rhs_mom[dit].plus(m_JaJcrossBv_cc[dit]);
      }
   }
   
   // update RHS for ion energy density
   LevelData<FArrayBox>& rhs_ene_ion( a_rhs_fluid.cell_var("energyDensity") );
   m_geometry.mappedGridDivergenceFromFluxNorms(m_enFluxi_norm, m_dummyDiv);
   for (DataIterator dit(rhs_ene_ion.dataIterator()); dit.ok(); ++dit) {
      m_dummyDiv[dit].divide(m_cellVol[dit]);
      m_dummyDiv[dit].mult(m_Jacobian[dit]);
      rhs_ene_ion[dit].minus( m_dummyDiv[dit] );
      rhs_ene_ion[dit].plus( m_JaSeneDeni[dit] );
   }
   
   // update RHS for ele energy density
   LevelData<FArrayBox>& rhs_ene_ele( a_rhs_fluid.cell_var("energyDensity_ele") );
   m_geometry.mappedGridDivergenceFromFluxNorms(m_enFluxe_norm, m_dummyDiv);
   for (DataIterator dit(rhs_ene_ele.dataIterator()); dit.ok(); ++dit) {
      m_dummyDiv[dit].divide(m_cellVol[dit]);
      m_dummyDiv[dit].mult(m_Jacobian[dit]);
      rhs_ene_ele[dit].minus( m_dummyDiv[dit] );
      rhs_ene_ele[dit].plus( m_JaSeneDene[dit] );
   }
   
   // update RHS for momentum density in virtual direction
   if(a_soln_fluid.m_evolve_momentumDensity_virtual){
      LevelData<FArrayBox>& rhs_momV( a_rhs_fluid.cell_var("momentumDensity_virtual") );
      m_geometry.mappedGridDivergenceFromFluxNorms(m_mvFlux_norm, m_dummyDiv);
      for (DataIterator dit(rhs_momV.dataIterator()); dit.ok(); ++dit) {
         FArrayBox& div_mvFlux = m_dummyDiv[dit];
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
   
   } // end if updateHydroVars

   // update RHS for virtual magnetic field
   if(a_soln_fluid.m_evolve_magneticField_virtual){
      
      // the following takes covar E and returns Jacobian*curl(E)cdot(g^y)
      const LevelData<EdgeDataBox>& soln_Ecovar( a_soln_fluid.edge_var("electricField") );
      m_geometry.mappedGridCurlFromEdgeTans(soln_Ecovar,2,m_curlE);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         if(m_twoDaxisymm) {
            m_curlE[dit].mult(m_g_y[dit]); // g^y*g_y = 1
         }
      } 
      
      LevelData<FArrayBox>& rhs_magV( a_rhs_fluid.cell_var("magneticField_virtual") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs_magV[dit].minus( m_curlE[dit] );
      }

   }

}   

void ExtendedMhdOp::initialize( CFGVars&  a_species, 
                          const Real      a_time )
{
   CH_TIME("ExtendedMhdOp::initialize()");
   
   // initilize cell center variables
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
   for (int n=0; n<a_species.num_face_vars(); ++n) {
      string this_face_var_comp_name;
      for(int m=0; m<SpaceDim; ++m) {
         if(m==0) this_face_var_comp_name = a_species.face_var_name(n)+"_0";
         if(m==1) this_face_var_comp_name = a_species.face_var_name(n)+"_1";
         if(m==2) this_face_var_comp_name = a_species.face_var_name(n)+"_2";
         const GridFunction& fsic( fluidSpeciesIC( this_face_var_comp_name ) );
         LevelData<FluxBox> this_face_var_comp( a_species.face_var(n).getBoxes(), 1,
                                                a_species.face_var(n).ghostVect() );
         fsic.assign( this_face_var_comp, a_species.configurationSpaceGeometry(), a_time );
         LevelData<FluxBox>& face_var_ref( a_species.face_var(n) );
         for (DataIterator dit( this_face_var_comp.dataIterator() ); dit.ok(); ++dit) {
            face_var_ref[dit].copy(this_face_var_comp[dit],m,0,0,1);
         }
      }
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
      if(a_species.edge_var_name(n)=="currentDensity") {
         geometry.convertPhysToCovar(this_edge_var_on_cell,0);
         SpaceUtils::interpCellToEdges(a_species.edge_var(n),this_edge_var_on_cell,a_species.edge_var(n),"c2");
      }

      if(a_species.edge_var_name(n)=="electricField") {

         // overwrite initial E with that from Ohm's law for initial value
         // (currently just adding ideal term to non-ideal term, which is in input file.... HACK!)
         geometry.convertPhysToCovar(this_edge_var_on_cell,0);
         SpaceUtils::interpCellToEdges(a_species.edge_var(n),this_edge_var_on_cell,a_species.edge_var(n),"c2");

         // compute charge density at nodes from electric field
         // and interpolate to cell center for diagnostic
         SpaceUtils::exchangeEdgeDataBox(a_species.edge_var(n));
         const int this_ebc = m_electricField_ebc;
         m_fluid_edge_var_bc.at(this_ebc)->applyEdgeBC( static_cast<FluidSpecies&>(a_species), 
                                                        a_species.edge_var(n), a_time );
      
         const DisjointBoxLayout& grids( a_species.edge_var(n).getBoxes() ); 
         m_geometry.physicalDivergenceFromEdgeTans(m_divE, a_species.edge_var(n));
         SpaceUtils::interpNodesToEdges(m_deltaN_ce, m_divE, "c2");
         for (DataIterator dit(grids); dit.ok(); ++dit) {
            for (int dir=0; dir<SpaceDim; dir++) {
               FArrayBox& thisDeltaN( m_deltaN_ce[dit][dir] ); 
               thisDeltaN.mult( m_delta0*m_lambdai0/m_epRel );
            }
         }
         SpaceUtils::interpEdgesToCell(m_deltaN_cc, m_deltaN_ce, "c2");
      
      }
   
   }

   if ( typeid(a_species) == typeid(FluidSpecies&) ) {
      (static_cast<FluidSpecies&>(a_species)).convertFromPhysical(); // this just multiplies by Jacobian
   }

}

#include "NamespaceFooter.H"
