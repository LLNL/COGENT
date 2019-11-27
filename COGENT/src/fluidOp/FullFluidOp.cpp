#include "FullFluidOp.H"
#include "EllipticOpBCFactory.H"
#include "FluidVarBCFactory.H"
#include "FourthOrderUtil.H"
#include "MagBlockCoordSys.H"
#include "inspect.H"
#include "CONSTANTS.H"

#include "NodeFArrayBox.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "GKOps.H"
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "FluidOpF_F.H"

#include "NamespaceHeader.H" 

const char* FullFluidOp::pp_name = {"fullfluid_op"};


FullFluidOp::FullFluidOp( const string&   a_pp_str,
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
     m_my_pc_idx(-1)

{
   ParmParse pp(a_pp_str.c_str());
   parseParameters( pp );
   if (m_verbosity>0) {
      printParameters();
   }

   // Get initial conditions for cell variables and define boundary conditions instance
   const std::string icsn( "IC." + a_species_name);
   std:: string variable_name;
   ParmParse ppic( icsn.c_str() );
   FluidVarBCFactory fluid_var_bc_factory;
  
   // Set IC and create BC for mass density
   variable_name = "density";
   CH_assert( ppic.contains("density.function") );
   parseInitialConditions(a_species_name, variable_name);
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );
 
   // Set IC and create BC for momentum density
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
 
   // Set IC and create BC for energy density
   variable_name = "energyDensity";
   CH_assert( ppic.contains("energyDensity.function") );
   parseInitialConditions(a_species_name, variable_name);
   m_fluid_bc.push_back( fluid_var_bc_factory.create(a_species_name, variable_name, 
                         a_geometry.getCoordSys()->type(), false) );
   
   // Set IC and create BC for the momentum in the "ignorable" direction (used for 2D)
   if (ppic.contains("momentumDensity_virtual.function")) {
      variable_name = "momentumDensity_virtual";
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

}


FullFluidOp::~FullFluidOp()
{
   //delete m_fluid_variable_bc;
   for (std::vector<FluidVarBC*>::iterator it = m_fluid_bc.begin(); it != m_fluid_bc.end(); ++it)
   {
      delete (*it);
   }
   m_fluid_bc.clear();
}


void FullFluidOp::accumulateRHS(FluidSpeciesPtrVect&               a_rhs,
                                const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                const FluidSpeciesPtrVect&         a_fluid_species,
                                const PS::ScalarPtrVect&           a_scalars,
                                const EField&                      a_E_field,
                                const int                          a_fluidVecComp,
                                const Real                         a_time)
{
   // Placeholder for functional calls that could be
   // used in either Explicit or Implicit RHS's
   
}


void FullFluidOp::accumulateExplicitRHS(FluidSpeciesPtrVect&               a_rhs,
                                        const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                        const FluidSpeciesPtrVect&         a_fluid_species,
                                        const PS::ScalarPtrVect&           a_scalars,
                                        const EField&                      a_E_field,
                                        const int                          a_fluidVecComp,
                                        const Real                         a_time)
{
   CH_TIME("FullFluidOp::accumulateExplicitRHS()");

   // Get RHS fluid species
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   LevelData<FArrayBox>& rhs_data( rhs_fluid.cell_var(0) );
   const DisjointBoxLayout& grids( rhs_data.getBoxes() );

   // Get solution fluid species
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
   //const LevelData<FArrayBox>& soln_data( soln_fluid.cell_var(0) );
   
  
   CH_assert(soln_fluid.m_evolve_momentumDensity == 1);
   CH_assert(soln_fluid.m_evolve_energyDensity == 1);
   
   setCellCenterValues( soln_fluid );

   setCellCenterFluxes( soln_fluid );
   
   // Compute flux-freezing speed at cell-center for each direction
   // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/N); i = q0, q1, q2
   const double gamma = soln_fluid.m_gamma;
   setCspeed( m_rhoDen_cc, gamma );
   SpaceUtils::upWindToFaces(m_CspeedR_norm, m_Cspeed_cc, m_CspeedR_norm, "c2"); // 3rd arg not used for "c2"
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_CspeedL_norm[dit].copy(m_CspeedR_norm[dit],m_CspeedL_norm[dit].box());
      m_CspeedL_norm[dit].negate();
   }
   
   setFaceCenteredFluxes( soln_fluid );

   enforceFluxBCs( rhs_fluid, soln_fluid, a_time );

   const Real mass = soln_fluid.mass();
   const Real charge = soln_fluid.charge();
   
   updateRHSs( rhs_fluid, soln_fluid, a_E_field.getCellCenteredField(), mass, charge );

}


void FullFluidOp::accumulateImplicitRHS(FluidSpeciesPtrVect&              a_rhs,
                                        const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                        const FluidSpeciesPtrVect&        a_fluid_species,
                                        const PS::ScalarPtrVect&          a_scalars,
                                        const EField&                     a_E_field,
                                        const int                         a_fluidVecComp,
                                        const Real                        a_time )
{

}


void FullFluidOp::defineBlockPC( std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
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


void FullFluidOp::updateBlockPC(std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
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
                << " for IdealMhFullFluidf fluid species " << a_species_idx << ".\n";
    }
  
    PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps> *pc 
      = dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>(a_pc[m_my_pc_idx]);
    CH_assert(pc != NULL);
    pc->update(a_kin_species_phys, a_fluid_species, a_time, a_shift, a_im, a_species_idx);
  }
  return;
}


void FullFluidOp::updatePCImEx(const FluidSpeciesPtrVect&       a_fluid_species,
                               const PS::KineticSpeciesPtrVect& a_kinetic_species,
                               const double                     a_time,
                               const double                     a_shift,
                               const int                        a_component)
{
   CH_TIME("FullFluidOp::updatePCImEx()");
   
   // Update preconditioner for implicit terms here
   
}


void FullFluidOp::solvePCImEx(FluidSpeciesPtrVect&              a_fluid_species_solution,
                              const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                              const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                              const int                         a_component )
{
   CH_TIME("FullFluidOp::solvePCImEx()");
   if (!m_is_time_implicit) {
     MayDay::Error("Error in FullFluidOp::solvePCImEx(): m_is_time_implicit is false. This function shouldn't have been called at all!");
   }
   
   // Solve preconditioner operator for implicit terms here

}

void FullFluidOp::convertPhysFluxToMappedFlux(LevelData<FArrayBox>&        a_mapped_flux,
                                              const LevelData<FArrayBox>&  a_physical_flux,
                                              const LevelData<FArrayBox>&  a_N) const
{
   const DisjointBoxLayout& grids( m_geometry.grids() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {

      const FArrayBox& phys_on_patch = a_physical_flux[dit];
      const FArrayBox& N_on_patch = a_N[dit]; 
      FArrayBox& mapped_on_patch = a_mapped_flux[dit];
      
      convertPhysFluxToMappedFlux(mapped_on_patch,phys_on_patch,N_on_patch);

   }
}

void FullFluidOp::convertPhysFluxToMappedFlux(FArrayBox&          a_mapped_on_patch,
                                              const FArrayBox&    a_phys_on_patch,
                                              const FArrayBox&    a_N_on_patch ) const
{
   CH_TIME("FullFluidOp::convertPhysFluxToMappedFlux()");
   const Box& thisbox = a_mapped_on_patch.box();

   FORT_CONVERT_TO_MAPPED( CHF_BOX(thisbox), 
                           CHF_CONST_FRA(a_phys_on_patch),
                           CHF_CONST_FRA(a_N_on_patch),
                           CHF_CONST_INT(SpaceDim),
                           CHF_FRA(a_mapped_on_patch) );

}

void FullFluidOp::computeNTFfaceArea(LevelData<FluxBox>&          a_Flux_norm,
                                     const LevelData<FArrayBox>&  a_Flux_cc,
                                     const LevelData<FArrayBox>&  a_Jf_cc,
                                     const LevelData<FArrayBox>&  a_Cspeed_cc,
                                     const LevelData<FluxBox>&    a_CspeedR_norm,
                                     const LevelData<FluxBox>&    a_CspeedL_norm ) const
{
   CH_TIME("FullFluidOp::computeNTFfaceArea()");
   CH_assert(a_Flux_norm.nComp() == 1);
   CH_assert(a_Flux_cc.nComp() == SpaceDim);
   
   const DisjointBoxLayout& grids( a_Jf_cc.getBoxes() );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
   
      // Convert cell-center flux from physical to mapped
      const Box& cell_box( a_Jf_cc[dit].box() );
      FArrayBox Flux_cc_mapped(cell_box, a_Flux_cc.nComp());
      convertPhysFluxToMappedFlux(Flux_cc_mapped,a_Flux_cc[dit],m_Nmatrix[dit]);

      /* Compute Lax splitting on mapped flux and interpolate FluxL and FluxR to face norms
      //
      //   get left and right going flux at cell-center
      //   fluxR = 0.5*(flux + Cspeed*fun),
      //   fluxL = 0.5*(flux - Cspeed*fun),
      //   Cspeed = abs(max(eigenValue of Flux Jacobian))
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

      const Box& grid_box( grids[dit] );              // this box has no ghost
      const Box& flux_box( a_Flux_norm[dit].box() );  // this box has ghost
      FluxBox thisFluxL_norm( flux_box, 1 );          // need box with ghosts here
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

      // Compute total norm flux by adding left and right together
      // and then multiply by face area
      a_Flux_norm[dit].copy(thisFluxR_norm, flux_box );
      a_Flux_norm[dit] += thisFluxL_norm;
      const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const RealVect& faceArea = coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& this_face = a_Flux_norm[dit][dir];
         this_face.mult(faceArea[dir]);
      }
   }
   a_Flux_norm.exchange();
 
}

void FullFluidOp::fillGhostCells( FluidSpecies&  a_species_phys,
                              const double   a_time )
{
   CH_TIME("FullFluidOp::fillGhostCells()");
   for (int n=0; n<a_species_phys.num_cell_vars(); ++n) {
      LevelData<FArrayBox>& fld( a_species_phys.cell_var(n) );
      // Fill ghost cells except for those on physical boundaries
      m_geometry.fillInternalGhosts( fld );
      m_fluid_bc.at(n)->apply( a_species_phys, a_time );
   }
}


void FullFluidOp::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "time_implicit", m_is_time_implicit);

   GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
   std::string grid_function_name;
   
   if (a_pp.contains("advScheme")) {
      a_pp.get("advScheme", m_advScheme );
   }
   
}


void FullFluidOp::printParameters()
{
   if (procID()==0) {
      std::cout << "FullFluidOp parameters:" << std::endl;
      std::cout << " time_implicit  =  " << m_is_time_implicit << std::endl;
      std::cout << " advScheme  =  " << m_advScheme << std::endl;
   }
}

void FullFluidOp::setCourantTimeStep( const LevelData<FArrayBox>&  a_Cspeed )
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
   }
   m_courant_time_step = this_dt_max;
}

void FullFluidOp::setCspeed(const LevelData<FArrayBox>&  a_density,
                            const double                 a_gamma )
{
   CH_TIME("FullFluidOp::setCspeed()");

   CH_assert(m_Cspeed_cc.nComp() == SpaceDim);
   
   const DisjointBoxLayout& grids( m_Cspeed_cc.getBoxes() ); 
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const FArrayBox& rho_on_patch  = a_density[dit];   
      const FArrayBox& V_on_patch    = m_velocity[dit]; 
      const FArrayBox& P_on_patch    = m_pressure[dit]; 
      const FArrayBox& N_on_patch    = m_Nmatrix[dit];
      const FArrayBox& J_on_patch    = m_Jacobian[dit]; 
      FArrayBox& C_on_patch = m_Cspeed_cc[dit];
      const Box& thisbox = C_on_patch.box();
      
      FArrayBox magP_on_patch(thisbox, 1);
      magP_on_patch.setVal(0.);
      
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

#if 0
Real FullFluidOp::computeDtExplicitTI( const FluidSpeciesPtrVect&  a_fluid_species )
{
   CH_TIME("Real FullFluidOp::computeDtExplicitTI()");

   // Get physical values for m_species_name and calculate freezing speed
   // Cspeed_i = |NTVdotqihat| + |N_i|sqrt(gamma*P/N); i = q0, q1, q2

   double gamma;

   for (int species(0); species<a_fluid_species.size(); species++) {
      const FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
      const std::string species_name( fluid_species.name() );
      if(species_name==m_species_name) {
         gamma = fluid_species.m_gamma;
         fluid_species.massDensity(m_rhoDen_cc);
         fluid_species.velocity(m_velocity);  // in-plane velocity vector
         fluid_species.pressure(m_pressure);
   
         // Calls for physical variables here are from mapped state
         // vector (has not been convertToPhysical )
         // and so need to divide by J
         m_geometry.divideJonValid(m_rhoDen_cc);
         m_geometry.divideJonValid(m_pressure);
         
         break;
      }
   }
   setCspeed( m_rhoDen_cc, gamma );
   setCourantTimeStep(m_Cspeed_cc);   
   

   return m_courant_time_step;
}
#endif

void FullFluidOp::defineLevelDatas( const DisjointBoxLayout&  a_grids,
                                   const IntVect&            a_ghostVect )
{
   // Get the transformation matrix N and Jacobian at cell-center needed to transform
   // from physical space to mapped space
   //
   m_Jacobian.define(a_grids, 1, a_ghostVect);
   m_Nmatrix.define(a_grids, SpaceDim*SpaceDim, a_ghostVect);
   m_X_phys.define(a_grids, SpaceDim, a_ghostVect);
   m_g_y.define(a_grids, 1, a_ghostVect); // used only when solving virtual B
   m_cellVol.define(a_grids, 1, IntVect::Zero);
   double twoPi = 2.0*Pi;
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
   
   m_rhoDen_cc.define(a_grids, 1, a_ghostVect);
   m_momDen_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_eneDen_cc.define(a_grids, 1, a_ghostVect);
   m_momDen_virtual_cc.define(a_grids, 1, a_ghostVect);

   m_pressure.define(a_grids, 1, a_ghostVect);
   m_velocity.define(a_grids, SpaceDim, a_ghostVect);
   m_velocity_virtual.define(a_grids, 1, a_ghostVect);

   m_rhoFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mxFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_myFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mzFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_enFlux_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_mvFlux_cc.define(a_grids, SpaceDim, a_ghostVect);

   m_Cspeed_cc.define(a_grids, SpaceDim, a_ghostVect);
   m_CspeedL_norm.define(a_grids, 1, a_ghostVect);
   m_CspeedR_norm.define(a_grids, 1, a_ghostVect);

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
      
   m_rhoFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_mxFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_myFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_mzFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_enFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);
   m_mvFlux_ce.define(a_grids, SpaceDim, 1*IntVect::Unit);

   m_rhoFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_mxFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_myFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_mzFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_enFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_mvFluxBC_norm.define(a_grids, 1, 1*IntVect::Unit);
   m_momFluxBC_norm.define(a_grids, SpaceDim, 1*IntVect::Unit);
}

void FullFluidOp::setCellCenterValues( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("FullFluidOp::setCellCenterValues()");
   
   const LevelData<FArrayBox>& soln_rhoDen( a_soln_fluid.cell_var(0) );
   const LevelData<FArrayBox>& soln_momDen( a_soln_fluid.cell_var("momentumDensity") );
   const LevelData<FArrayBox>& soln_eneDen( a_soln_fluid.cell_var("energyDensity") );
   
   // Set values for fluid variables
   const DisjointBoxLayout& grids( soln_rhoDen.getBoxes() );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_rhoDen_cc[dit].copy( soln_rhoDen[dit], m_rhoDen_cc[dit].box() );
      m_momDen_cc[dit].copy( soln_momDen[dit], m_momDen_cc[dit].box() );
      m_eneDen_cc[dit].copy( soln_eneDen[dit], m_eneDen_cc[dit].box() );
   }

   if (a_soln_fluid.m_evolve_momentumDensity_virtual) {
      const LevelData<FArrayBox>& soln_momDen_virtual( a_soln_fluid.cell_var("momentumDensity_virtual") );
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_momDen_virtual_cc[dit].copy( soln_momDen_virtual[dit] );
      }
   }

   // Set physical derived variables for fluid species
   a_soln_fluid.pressure(m_pressure);
   a_soln_fluid.velocity(m_velocity);                 // in-plane velocity vector
   a_soln_fluid.velocity_virtual(m_velocity_virtual); // out-of-plane velocity vector
}

void FullFluidOp::setCellCenterFluxes( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("FullFluidOp::setCellCenterFluxes()");
   
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
      */
      if (a_soln_fluid.m_evolve_momentumDensity_virtual) {
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


void FullFluidOp::setFaceCenteredFluxes( const FluidSpecies&  a_soln_fluid )
{
   CH_TIME("FullFluidOp::setFaceCenteredFluxes()");
   
   // Get mapped flux normals on each face multiplied by face area (NTFlux*faceArea)
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

void FullFluidOp::enforceFluxBCs(FluidSpecies&  a_rhs_fluid,
                                 const FluidSpecies&  a_soln_fluid,
                                 const Real           a_time )
{
   CH_TIME("FullFluidOp::enforceFluxBCs()");
   
   // Enforce flux BC for mass density flux
   SpaceUtils::upWindToFaces( m_rhoFlux_ce, m_rhoFlux_cc, m_rhoFlux_ce, "c2" );
   m_geometry.applyAxisymmetricCorrection( m_rhoFlux_ce );
   m_geometry.computeMetricTermProductAverage( m_rhoFluxBC_norm, m_rhoFlux_ce, 0 );
   m_fluid_bc.at(0)->applyFluxBC( a_rhs_fluid, m_rhoFlux_norm, m_rhoFluxBC_norm, a_time );
   
   // Enforce flux BC for momentum density flux
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
  
   // Enforce flux BC for energy density flux
   SpaceUtils::upWindToFaces( m_enFlux_ce, m_enFlux_cc, m_enFlux_ce, "c2" );
   m_geometry.applyAxisymmetricCorrection( m_enFlux_ce );
   m_geometry.computeMetricTermProductAverage( m_enFluxBC_norm, m_enFlux_ce, 0 );
   this_cvc = a_soln_fluid.cell_var_component("energyDensity");
   m_fluid_bc.at(this_cvc)->applyFluxBC( a_rhs_fluid, m_enFlux_norm, m_enFluxBC_norm, a_time );

   // Enforce flux BC for virtual momentum density flux
   if(a_soln_fluid.m_evolve_momentumDensity_virtual) {
      SpaceUtils::upWindToFaces( m_mvFlux_ce, m_mvFlux_cc, m_mvFlux_ce, "c2" );
      m_geometry.applyAxisymmetricCorrection( m_mvFlux_ce );
      m_geometry.computeMetricTermProductAverage( m_mvFluxBC_norm, m_mvFlux_ce, 0 );
      this_cvc = a_soln_fluid.cell_var_component("momentumDensity_virtual");
      m_fluid_bc.at(this_cvc)->applyFluxBC( a_rhs_fluid, m_mvFlux_norm, m_mvFluxBC_norm, a_time );
   }
}

void FullFluidOp::updateRHSs(FluidSpecies&                a_rhs_fluid,
                             const FluidSpecies&          a_soln_fluid,
                             const LevelData<FArrayBox>&  a_E_field,
                             const Real&                  a_mass,
                             const Real&                  a_charge)
{
   CH_TIME("FullFluidOp::updateRHSs()");
   
   // compute divergence of Fluxes and add (subtract) to RHS 
   // Note that J/cellVol = 1/mapVol
   //
   const DisjointBoxLayout& grids( m_rhoDen_cc.getBoxes() ); 

   // Update RHS for mass density
   LevelData<FArrayBox>& rhs_data( a_rhs_fluid.cell_var(0) );
   for (DataIterator dit(rhs_data.dataIterator()); dit.ok(); ++dit) {
      FArrayBox div_rhoFlux(grids[dit], 1);
      m_geometry.mappedGridDivergenceFromFluxNorms(m_rhoFlux_norm[dit], div_rhoFlux);
      div_rhoFlux /= m_cellVol[dit];
      div_rhoFlux.mult(m_Jacobian[dit]);
      rhs_data[dit].minus(div_rhoFlux);
   }
  
   // Update RHS for momentum density
   LevelData<FArrayBox>& rhs_data2( a_rhs_fluid.cell_var("momentumDensity") );
   for (DataIterator dit(rhs_data2.dataIterator()); dit.ok(); ++dit) {
      FArrayBox div_momFlux(grids[dit], SpaceDim);
      m_geometry.mappedGridDivergenceFromFluxNorms(m_momFlux_norm[dit], div_momFlux);
      for (int n=0; n<div_momFlux.nComp(); ++n) {
         if(n==0 && m_twoDaxisymm) {
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

      rhs_data2[dit].minus(div_momFlux);
      if(m_twoDaxisymm && a_soln_fluid.m_evolve_momentumDensity_virtual) {
         FArrayBox JaRhoVthsqOverR(grids[dit], 1);
         JaRhoVthsqOverR.copy(m_momDen_virtual_cc[dit]);
         JaRhoVthsqOverR.mult(m_velocity_virtual[dit]);
         JaRhoVthsqOverR.mult(m_Jacobian[dit]); // Jacobian = Ja = 2*pi*R_phys * J_RZ
         JaRhoVthsqOverR.divide(m_X_phys[dit],0,0,1);
         rhs_data2[dit].plus(JaRhoVthsqOverR,0,0,1);            
      }

      // Add source terms (J*n*q*[vxB], J*q*n*E)
      accumulateMomDensitySources(rhs_data2[dit], a_E_field[dit], m_rhoDen_cc[dit],
                                  m_Jacobian[dit], a_mass, a_charge);
     
   }
   
   // Update RHS for energy density
   LevelData<FArrayBox>& rhs_data3( a_rhs_fluid.cell_var("energyDensity") );
   for (DataIterator dit(rhs_data3.dataIterator()); dit.ok(); ++dit) {
      FArrayBox div_enFlux(grids[dit], 1);
      m_geometry.mappedGridDivergenceFromFluxNorms(m_enFlux_norm[dit], div_enFlux);
      div_enFlux /= m_cellVol[dit];
      div_enFlux.mult(m_Jacobian[dit]);
      rhs_data3[dit].minus( div_enFlux );

      // Add source terms here (J*v*q*n*E)
      accumulateEnergyDensitySources(rhs_data2[dit], a_E_field[dit], m_rhoDen_cc[dit],
                                     m_velocity[dit], m_Jacobian[dit], a_mass, a_charge);
   }
   
   // Update RHS for momentum density in virtual direction
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
}

void FullFluidOp::accumulateMomDensitySources(FArrayBox& a_rhs,
                                              const FArrayBox& a_E_field,
                                              const FArrayBox& a_rhoDen,
                                              const FArrayBox& a_J,
                                              const Real& a_mass,
                                              const Real& a_charge) const
{
   // Compute source = J * (rho/m) * q * E
   const Box& box = a_rhs.box();
   FArrayBox source(box, SpaceDim);
   
   if (SpaceDim == 3) {
      source.copy(a_E_field);
   }
   else {
      source.copy(a_E_field, 0, 0, 1);
      source.copy(a_E_field, 2, 1, 1);
   }

   source.divide(a_mass);
   source.mult(a_charge);
   
   for (int n=0; n<SpaceDim; ++n) {
      source.mult(a_rhoDen, 0, n, 1);
      source.mult(a_J, 0, n, 1);
   }
   
   // Add the source term to rhs
   a_rhs.plus(source);
}

void FullFluidOp::accumulateEnergyDensitySources(FArrayBox& a_rhs,
                                                 const FArrayBox& a_E_field,
                                                 const FArrayBox& a_rhoDen,
                                                 const FArrayBox& a_velocity,
                                                 const FArrayBox& a_J,
                                                 const Real& a_mass,
                                                 const Real& a_charge) const
{
   // Compute J * (rho/m) * q * v * E
   const Box& box = a_rhs.box();
   FArrayBox source(box, SpaceDim);
   
   if (SpaceDim == 3) {
      source.copy(a_E_field);
   }
   else {
      source.copy(a_E_field, 0, 0, 1);
      source.copy(a_E_field, 2, 1, 1);
   }
   
   source.divide(a_mass);
   source.mult(a_charge);
   
   for (int n=0; n<SpaceDim; ++n) {
      source.mult(a_rhoDen, 0, n, 1);
      source.mult(a_velocity, 0, n, 1);
      source.mult(a_J, 0, n, 1);
   }
   
   // Add the source term to rhs
   a_rhs.plus(source);
}

#include "NamespaceFooter.H"
