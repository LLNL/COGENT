#include "GKVlasov.H"
#include "CONSTANTS.H"
#include "Directions.H"
#include "MomentOp.H"
#include "Kernels.H"

#include "FourthOrderUtil.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"

#include "mappedLimiterF_F.H"
#include "altFaceAverages.H"
#include "mappedAdvectionFlux.H"
#include "KineticFunctionLibrary.H"
#include "GKVlasovPreconditioner.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "SimpleDivergence.H"
#include "FluxSurface.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM


#undef TEST_ZERO_DIVERGENCE
#undef TEST_ENERGY_CONSERVATION

#include "NamespaceHeader.H"

using std::max;

const char* GKVlasov::pp_name = {"gkvlasov"};

Real GKVlasov::s_stability_bound[NUM_FLUX] = {2.06,2.7852,1.7453,1.7320,1.7320,1.7453};

Real
MaxNorm( const LevelData<FArrayBox>& a )
{
   const DisjointBoxLayout& grids = a.disjointBoxLayout();

   double local_max = -DBL_MAX;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      double this_max = a[dit].max();
      if (this_max > local_max) local_max = this_max;
   }

   double global_max;
#ifdef CH_MPI
   MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
   global_max = local_max;
#endif

   return global_max;
}

GKVlasov::GKVlasov( ParmParse&                      a_pp,
                    const Real                      a_larmor_number )
  : m_larmor_number(a_larmor_number),
    m_face_avg_type(INVALID),
    m_dt_dim_factor(1.0),
    m_first_call(true),
    m_include_FLR_effects(false),
    m_plot_FLR_switch(false),
    m_update_pc_freq(1),
    m_update_pc_skip_stage(false)
{
   if (a_pp.contains("limiter")) {
      if ( procID()==0 ) MayDay::Warning("GKVlasov: Use of input flag 'limiter' deprecated");
   }

   if (a_pp.contains("verbose") && procID()==0) {
      a_pp.get("verbose", m_verbose);
   }
   else {
      m_verbose = false;
   }

   a_pp.query("include_FLR_effects", m_include_FLR_effects);
   if (m_include_FLR_effects) {
     a_pp.query("FLR_plot_switch_func", m_plot_FLR_switch);
     if (a_pp.contains("FLR_switch_func")) {
       CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
       std::string func_name;
       a_pp.get("FLR_switch_func", func_name);
       m_FLR_switch_func = grid_library->find(func_name);
     }
   }
   a_pp.query("update_precond_interval", m_update_pc_freq);
   if (m_update_pc_freq > 1) {
      m_update_pc_skip_stage = true;
      a_pp.query("update_precond_skip_stage", m_update_pc_skip_stage);
   }

   m_ti_op_type.clear();
   
   if (a_pp.contains("subtract_maxwellian_background")) {
      a_pp.get("subtract_maxwellian_background", m_subtract_maxwellian);
      m_update_maxwellian = true;
   }
   else {
      m_subtract_maxwellian = false;
      m_update_maxwellian = false;
   }
   
   KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   std::string function_name;
   if (a_pp.contains( "ref_function")) {
      a_pp.get( "ref_function", function_name );
      m_ref_func = library->find( function_name );
      m_update_maxwellian = false;
   }

   if (a_pp.contains("time_step_diagnostics")) {
      a_pp.get("time_step_diagnostics", m_time_step_diagnostics);
   }
   else {
      m_time_step_diagnostics = false;
   }

   if (a_pp.contains("compute_cfl")) {
      a_pp.get("compute_cfl", m_compute_cfl);
   }
   else {
      m_compute_cfl = true;
   }
   
   if (a_pp.contains("face_avg_type")) {
      std::string dummy;
      a_pp.get("face_avg_type", dummy);
      if (dummy.compare("uw1")==0) {
         m_face_avg_type = UW1;
      }
      else if (dummy.compare("uw3")==0) {
         m_face_avg_type = UW3;
      }
      else if (dummy.compare("uw5")==0) {
         m_face_avg_type = UW5;
      }
      else if (dummy.compare("weno5")==0) {
         m_face_avg_type = WENO5;
      }
      else if (dummy.compare("bweno")==0) {
         m_face_avg_type = BWENO;
      }
      m_method = dummy;
   }
   else {
      if ( procID()==0 ) MayDay::Warning("Flux type unspecified; using default (BWENO)");
      m_face_avg_type = BWENO;
      m_method = "bweno";
   }
   // In theory, advection is in (PDIM-1) dimensions, so we could relax this
   // a little, but for now, let's be conservative.
   m_dt_dim_factor = (m_face_avg_type>PPM) ? sqrt(PDIM) : 1.0;

   string precond_prefix = string(a_pp.prefix()) + ".precond";

   ParmParse pp_precond(precond_prefix.c_str());

   if ( pp_precond.contains("build_test") ) {
      pp_precond.get("build_test", m_precond_build_test);
   }
   else {
      m_precond_build_test = false;
   }

   if ( pp_precond.contains("face_avg_type") ) {
      pp_precond.get("face_avg_type", m_precond_face_avg_type);
   }
   else {
      m_precond_face_avg_type = "uw1";
   }

	// parameter parsing for sparse grids
   ParmParse ppsg("sparse_grids");
   m_useSG = false;  //Don't use sparse grids by default
   ppsg.query( "useSGstencils", m_useSG );

}


GKVlasov::GKVlasov( ParmParse&                      a_pp,
                    const Real                      a_larmor_number,
                    const std::vector<std::string>& a_name_list )
  : GKVlasov(a_pp, a_larmor_number) 
{
   m_ti_op_type.clear();
   for (int n=0; n<a_name_list.size(); n++) {
     std::string species_name = a_name_list[n];
     std::string key = species_name + ".time_implicit";
     bool flag = false;
     a_pp.query(key.c_str(), flag);
     m_ti_op_type[species_name] = (flag ? _implicit_op_ : _explicit_op_);
   }
}


GKVlasov::~GKVlasov()
{
}


void GKVlasov::accumulateRHS( GKRHSData&                            a_rhs,
                              const KineticSpeciesPtrVect&          a_kinetic_phys,
                              const CFG::LevelData<CFG::FArrayBox>& a_phi,
                              const CFG::EField&                    a_E_field,
                              const bool                            a_implicit,
                              const Real&                           a_time )
{
   /*
     Evaluates the (negated) phase space divergence:
        rhs = - divergence_R ( R_dot soln ) - divergence_v ( v_dot soln )
   */
   KineticSpeciesPtrVect& rhs_kinetic = a_rhs.dataKinetic();

   for (int s(0); s<rhs_kinetic.size(); s++) {
      const KineticSpecies& soln_species( *(a_kinetic_phys[s]) );
      KineticSpecies& rhs_species( *(rhs_kinetic[s]) );
      if (a_implicit) {
         evalRHSImplicit(rhs_species,
                         soln_species,
                         a_phi,
                         a_E_field,
                         PhaseGeom::FULL_VELOCITY,
                         a_time );
         
      } else {
         evalRHSExplicit(rhs_species,
                         soln_species,
                         a_phi,
                         a_E_field,
                         PhaseGeom::FULL_VELOCITY,
                         a_time );
      }
   }
}


void
GKVlasov::evalRHSExplicit( KineticSpecies&                        a_rhs_species,
                           const KineticSpecies&                  a_soln_species,
                           const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                           const CFG::EField&                     a_E_field,
                           const int                              a_velocity_option,
                           const Real                             a_time )
{
  const std::string& name = a_soln_species.name();

  std::map<std::string, tiOpType>::const_iterator it;
  it = m_ti_op_type.find(name);
  if (it == m_ti_op_type.end()) {
    MayDay::Error("in GKVlasov::evalRHSExplicit - m_ti_op_type does not contain species info.");
  }

  if (m_ti_op_type[name] == _explicit_op_) {
    evalRHS(a_rhs_species, a_soln_species, a_phi, a_E_field, a_velocity_option, a_time);
  }
}


void
GKVlasov::evalRHSImplicit( KineticSpecies&                        a_rhs_species,
                           const KineticSpecies&                  a_soln_species,
                           const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                           const CFG::EField&                     a_E_field,
                           const int                              a_velocity_option,
                           const Real                             a_time )
{
  const std::string& name = a_soln_species.name();

  std::map<std::string, tiOpType>::const_iterator it;
  it = m_ti_op_type.find(name);
  if (it == m_ti_op_type.end()) {
    MayDay::Error("in GKVlasov::evalRHSImplicit - m_ti_op_type does not contain species info.");
  }

  if (m_ti_op_type[name] == _implicit_op_) {
    evalRHS(a_rhs_species, a_soln_species, a_phi, a_E_field, a_velocity_option, a_time);
  }
}


void
GKVlasov::evalRHS( KineticSpecies&                        a_rhs_species,
                   const KineticSpecies&                  a_soln_species,
                   const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                   const CFG::EField&                     a_E_field,
                   const int                              a_velocity_option,
                   const Real                             a_time )
{
   CH_TIME("GKVlasov::evalRHS");

   /*
     Evaluates the (negated) phase space divergence:
        rhs = - divergence_R ( R_dot soln ) - divergence_v ( v_dot soln )
   */
   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   
   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& dbl( soln_dfn.getBoxes() );

   const PhaseGeom& geometry( a_rhs_species.phaseSpaceGeometry() );
   IntVect ghostVect = (geometry.secondOrder()) ? IntVect::Zero : IntVect::Unit;

   if ( geometry.divFreeVelocity() ) {
      bool fourth_order_Efield = !a_E_field.secondOrder();
      
      if ( !m_flux_normal.isDefined()) {
         m_flux_normal.define( dbl, 1, ghostVect );
      }
   
      computeIntegratedFluxNormals(m_flux_normal,
                                   a_soln_species,
                                   a_phi,
                                   a_E_field.getCellCenteredField(),
                                   a_E_field.getPhiNode(),
                                   fourth_order_Efield,
                                   a_velocity_option,
                                   a_time);
   
      geometry.mappedGridDivergenceFromIntegratedFluxNormals( rhs_dfn, m_flux_normal );
   }
   
   else {
      if ( !m_flux.isDefined()) {
         m_flux.define( dbl, SpaceDim, ghostVect );
      }
   
      computePhysicalFlux(m_flux,
                          a_soln_species,
                          a_phi,
                          a_E_field,
                          a_velocity_option,
                          a_time,
                          true);
   

      const bool OMIT_NT(false);
      geometry.mappedGridDivergence( rhs_dfn, m_flux, OMIT_NT );
   }

   // Divide by cell volume and negate
   for (DataIterator dit( rhs_dfn.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys&
         block_coord_sys( geometry.getBlockCoordSys( dbl[dit] ) );
      double fac( -1.0 / block_coord_sys.getMappedCellVolume() );
      rhs_dfn[dit].mult( fac );
   }

#ifdef TEST_ZERO_DIVERGENCE
   if ( phase_geom.divFreeVelocity() ) {
      testZeroDivergence(m_velocity_normal, geometry);
   }
   else {
      testZeroDivergence(m_velocity, geometry);
   }
#endif

#ifdef TEST_ENERGY_CONSERVATION
   testEnergyConservation(rhs_dfn, a_soln_species, a_phi, a_E_field);
#endif
}

void
GKVlasov::computePhysicalFlux( LevelData<FluxBox>&                    a_flux,
                               const KineticSpecies&                  a_soln_species,
                               const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                               const CFG::EField&                     a_E_field,
                               const int                              a_velocity_option,
                               const Real                             a_time,
                               const bool                             a_apply_axisymmetric_correction)
{
   CH_TIME("GKVlasov::computePhysicalFlux");

   /*
     Evaluates the physical phase space flux: R_dot soln + v_dot soln
     modified by axisymmetric corrections. NB:: the corrections
     are not applied for the case of pointwise metrics
   */
   
   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& dbl( soln_dfn.getBoxes() );
   const PhaseGeom& geometry( a_soln_species.phaseSpaceGeometry() );

   IntVect ghostVect = (geometry.secondOrder()) ? IntVect::Zero : IntVect::Unit;
   if ( !m_velocity.isDefined()) {
     m_velocity.define( dbl, SpaceDim, ghostVect );
   }

   // Compute full gyrokinetic velocity
   if (a_soln_species.isGyrokinetic() && m_include_FLR_effects) {

     LevelData<FluxBox> gyroaveraged_E_field;
     getEFieldwFLR(gyroaveraged_E_field, a_soln_species, a_phi, a_E_field);
     a_soln_species.computeVelocity(  m_velocity,
                                      gyroaveraged_E_field,
                                      true,
                                      a_velocity_option,
                                      a_time,
                                      a_apply_axisymmetric_correction);

   } else {
     a_soln_species.computeVelocity(  m_velocity,
                                      a_E_field.getInjectedField(),
                                      false,
                                      a_velocity_option,
                                      a_time,
                                      a_apply_axisymmetric_correction);
   }

   // Compute Vlasov_flux[dfn]
   if (!m_subtract_maxwellian) {
      computeAdvectionFluxPhys( soln_dfn, m_velocity, a_flux, geometry );
   }
   
   // Compute Vlasov_flux[dfn-F0] + Vlasov_flux_high_order[F0]
   else {
      
      if ( !m_dfn_no_bstar.isDefined()) {
        m_dfn_no_bstar.define( dbl, soln_dfn.nComp(), soln_dfn.ghostVect() );
      }
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
         m_dfn_no_bstar[dit].copy(soln_dfn[dit]);
      }
      geometry.divideBStarParallel(m_dfn_no_bstar);
      
      
      if (!m_maxwellian_dfn.isDefined()) {
        m_maxwellian_dfn.define(dbl, soln_dfn.nComp(), soln_dfn.ghostVect());
      }
      if (!m_delta_dfn.isDefined()) {
         m_delta_dfn.define(dbl, soln_dfn.nComp(), soln_dfn.ghostVect());
      }
            
      computeDeltaF(a_soln_species, m_dfn_no_bstar, m_delta_dfn, m_maxwellian_dfn);
      geometry.multBStarParallel(m_delta_dfn);
      
      if (!m_delta_flux.isDefined()) {
         m_delta_flux.define( dbl, SpaceDim, ghostVect );
      }
      computeAdvectionFluxPhys( m_delta_dfn, m_velocity, m_delta_flux, geometry );
      
      // Update gk-velocity to include no zero-order terms
      if (a_soln_species.isGyrokinetic() && m_include_FLR_effects) {

        LevelData<FluxBox> gyroaveraged_E_field;
        getEFieldwFLR(gyroaveraged_E_field, a_soln_species, a_phi, a_E_field);
        a_soln_species.computeVelocity( m_velocity,
                                        gyroaveraged_E_field,
                                        true,
                                        PhaseGeom::NO_ZERO_ORDER_TERMS,
                                        a_time );
      
      } else {

        a_soln_species.computeVelocity( m_velocity,
                                        a_E_field.getInjectedField(),
                                        false,
                                        PhaseGeom::NO_ZERO_ORDER_TERMS,
                                        a_time );

      }

      geometry.multBStarParallel(m_maxwellian_dfn);
      computeAdvectionFluxPhys( m_maxwellian_dfn, m_velocity, a_flux, geometry);
  
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
         a_flux[dit] += m_delta_flux[dit];
      }
   }
}

void
GKVlasov::computeIntegratedFluxNormals(LevelData<FluxBox>&                    a_flux_normal,
                                       const KineticSpecies&                  a_soln_species,
                                       const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                                       const CFG::LevelData<CFG::FArrayBox>&  a_Efield_cell,
                                       const CFG::LevelData<CFG::FArrayBox>&  a_phi_node,
                                       const bool                             a_fourth_order_Efield,
                                       const int                              a_velocity_option,
                                       const Real                             a_time )
{
   /*
     Evaluates the normal component of the mapped
     phase space flux integrated over a cell face
     NT * (R_dot soln + v_dot soln) * face_area
   */

   CH_TIME("GKVlasov::computeIntegratedFluxNormals");
   
   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& dbl( soln_dfn.getBoxes() );

   const PhaseGeom& geometry( a_soln_species.phaseSpaceGeometry() );
   const CFG::MagGeom& mag_geom = geometry.magGeom();
   
   if ( !m_dfn_no_bstar.isDefined()) {
     m_dfn_no_bstar.define( dbl, soln_dfn.nComp(), soln_dfn.ghostVect() );
   }
   for (DataIterator dit(dbl); dit.ok(); ++dit) {
      m_dfn_no_bstar[dit].copy(soln_dfn[dit]);
   }
   geometry.divideBStarParallel(m_dfn_no_bstar);
      
   IntVect ghostVect = (geometry.secondOrder()) ? IntVect::Zero : IntVect::Unit;

   if ( !m_velocity_normal.isDefined()) {
     m_velocity_normal.define( dbl, 1, ghostVect );
   }

   a_soln_species.computeMappedVelocityNormals( m_velocity_normal,
                                                a_Efield_cell,
                                                a_phi_node,
                                                a_fourth_order_Efield,
                                                a_velocity_option );
   
     
   // Compute normals of Vlasov_flux[dfn]
   if (!m_subtract_maxwellian) {
      computeAdvectionFluxMapped( m_dfn_no_bstar, m_velocity_normal, a_flux_normal, geometry );
   }
   
   // Compute normals of Vlasov_flux[dfn-F0] + Vlasov_flux_high_order[F0]
   else {

      if (!m_maxwellian_dfn.isDefined()) {
         m_maxwellian_dfn.define(dbl, soln_dfn.nComp(), soln_dfn.ghostVect());
      }
      if (!m_delta_dfn.isDefined()) {
         m_delta_dfn.define(dbl, soln_dfn.nComp(), soln_dfn.ghostVect());
      }
           
      computeDeltaF(a_soln_species, m_dfn_no_bstar, m_delta_dfn, m_maxwellian_dfn);
           
      if (!m_delta_flux.isDefined()) {
         m_delta_flux_normal.define( dbl, 1, ghostVect );
      }
      
      computeAdvectionFluxMapped( m_delta_dfn, m_velocity_normal, m_delta_flux_normal, geometry );
       
      a_soln_species.computeMappedVelocityNormals(  m_velocity_normal,
                                                    a_Efield_cell,
                                                    a_phi_node,
                                                    a_fourth_order_Efield,
                                                    PhaseGeom::NO_ZERO_ORDER_TERMS );

      computeAdvectionFluxMapped( m_maxwellian_dfn, m_velocity_normal, a_flux_normal, geometry);

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
         a_flux_normal[dit] += m_delta_flux_normal[dit];
      }
   }

   // Enforce conservation
   if (!mag_geom.extrablockExchange()) {
      geometry.averageAtBlockBoundaries(a_flux_normal);
   }
}

void
GKVlasov::computeAdvectionFluxPhys( const LevelData<FArrayBox>& a_dist_fn,
                                    const LevelData<FluxBox>&   a_velocity,
                                    LevelData<FluxBox>&         a_flux,
                                    const PhaseGeom&            a_phase_geom )
{
   CH_TIMERS("GKVlasov::computeAdvectionFluxPhys");
   CH_TIMER("compFaceAver", t_face_aver);

   /*
      Compute the phase space flux given the input phase space advection
      velocity.   This is where the hyperbolic stuff connects.
   */

   IntVect ghostVect = (a_phase_geom.secondOrder()) ? IntVect::Zero : IntVect::Unit;
   if ( !m_face_dist.isDefined()) {
      m_face_dist.define( a_dist_fn.getBoxes(), a_dist_fn.nComp(), ghostVect );
   }
   
   CH_START(t_face_aver);
   AltFaceAverages alt_face_average(a_phase_geom);
   
   // If we're limiting the face-centered values, do it here
   bool mapped_velocity = false;
   if (m_face_avg_type==PPM) {
      computeFaceAverages( m_face_dist, a_dist_fn, a_phase_geom.secondOrder() );
      m_face_dist.exchange();
      applyMappedLimiter( m_face_dist, a_dist_fn, a_velocity, a_phase_geom, mapped_velocity );
   }
   else  {
      alt_face_average.computeFaceAverages(m_face_dist,
                                           a_dist_fn,
                                           a_velocity,
                                           a_phase_geom,
                                           m_method,
                                           mapped_velocity );
   }

   CH_STOP(t_face_aver);

   if ( a_phase_geom.secondOrder() ) {

      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxes( a_flux, m_face_dist, a_velocity, false ); // if second-order, don't need to worry about SG stencils
   }
   else {

      const DisjointBoxLayout& grids = a_flux.disjointBoxLayout();
		int num_gpts_in_flux=0;
		if (m_useSG) { num_gpts_in_flux=1; } // if you're using sparse grids, your flux needs a ghost cell
		IntVect gvec = num_gpts_in_flux*IntVect::Unit;
		LevelData<FluxBox> fourth_order_flux(grids, SpaceDim, gvec); 
      
      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxesGen( fourth_order_flux, m_face_dist, a_velocity, true, m_useSG ); // if fourth order, you do

      // Compute the second-order flux in valid plus ghost cell faces,
      // then overwrite with the fourth-order flux on the valid faces.
      CH_assert(a_flux.ghostVect() == IntVect::Unit);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_flux[dit].copy(a_velocity[dit]);

         Box box = grow(grids[dit],1);
         for (int dir=0; dir<SpaceDim; ++dir) {
            a_flux[dit].mult(m_face_dist[dit],box,0,dir,1);
         }
         
         a_flux[dit].copy(fourth_order_flux[dit],grids[dit]);
      }
   }

   if ( !a_phase_geom.secondOrder() ) {
     a_flux.exchange();
   }
}


void
GKVlasov::computeAdvectionFluxMapped( const LevelData<FArrayBox>& a_dist_fn,
                                      const LevelData<FluxBox>&   a_velocity_normal,
                                      LevelData<FluxBox>&         a_flux_normal,
                                      const PhaseGeom&            a_phase_geom )
{

   /*
      Compute the phase space flux normal integrated over a cell face.
      This is where the hyperbolic stuff connects.
   */

   CH_TIME("GKVlasov::computeAdvectionFluxMapped");

   // Construct appropriately accurate face-averages of phi and advVel
   IntVect ghostVect = (a_phase_geom.secondOrder()) ? IntVect::Zero : IntVect::Unit;
   if ( !m_face_dist.isDefined()) {
      m_face_dist.define( a_dist_fn.getBoxes(), a_dist_fn.nComp(), a_dist_fn.ghostVect() );
   }

   AltFaceAverages alt_face_average(a_phase_geom);
   
   // If we're limiting the face-centered values, do it here
   bool mapped_velocity = true;
   if (m_face_avg_type==PPM) {
      computeFaceAverages( m_face_dist, a_dist_fn, a_phase_geom.secondOrder() );
      m_face_dist.exchange();
      applyMappedLimiter( m_face_dist, a_dist_fn, a_velocity_normal, a_phase_geom, mapped_velocity );
   }
   else  {
      alt_face_average.computeFaceAverages(m_face_dist,
                                           a_dist_fn,
                                           a_velocity_normal,
                                           a_phase_geom,
                                           m_method,
                                           mapped_velocity );
   }

   if ( a_phase_geom.secondOrder() ) {

      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxesNormal( a_flux_normal, m_face_dist, a_velocity_normal, false );
   }
   else {

      const DisjointBoxLayout& grids = a_flux_normal.disjointBoxLayout();
      LevelData<FluxBox> fourth_order_flux_normal(grids, 1, IntVect::Zero);
      
      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxesNormal( fourth_order_flux_normal, m_face_dist, a_velocity_normal, true );

      // Compute the second-order flux in valid plus ghost cell faces,
      // then overwrite with the fourth-order flux on the valid faces.
      CH_assert(a_flux_normal.ghostVect() == IntVect::Unit);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_flux_normal[dit].copy(a_velocity_normal[dit]);

         Box box = grow(grids[dit],1);
         a_flux_normal[dit].mult(m_face_dist[dit],box,0,0,1);
         
         a_flux_normal[dit].copy(fourth_order_flux_normal[dit],grids[dit]);
      }
   }

   a_flux_normal.exchange();
}

void
GKVlasov::computeDeltaF(const KineticSpecies&         a_soln_species,
                        const LevelData<FArrayBox>&   a_dfn_no_bstar,
                        LevelData<FArrayBox>&         a_delta_dfn,
                        LevelData<FArrayBox>&         a_maxwellian_dfn)

{
   
   const PhaseGeom& geometry = a_soln_species.phaseSpaceGeometry();

   //Presently implemented only to second order
   CH_assert(geometry.secondOrder());

   if (!m_update_maxwellian) {

      if (m_first_call) {
         const DisjointBoxLayout& dbl( a_dfn_no_bstar.getBoxes() );
         m_F0.define(dbl, a_dfn_no_bstar.nComp(), a_dfn_no_bstar.ghostVect());

         KineticSpeciesPtr ref_species( a_soln_species.clone( IntVect::Unit, false ) );
         m_ref_func->assign( *ref_species, 0.0 );
         LevelData<FArrayBox>& ref_dfn( ref_species->distributionFunction() );
         for (DataIterator dit(m_F0.dataIterator()); dit.ok(); ++dit) {
            m_F0[dit].copy(ref_dfn[dit]);
         }
         geometry.divideJonValid(m_F0);
      
         //Filling the ghosts of m_F0
         for (DataIterator dit(m_F0.dataIterator()); dit.ok(); ++dit) {
            const PhaseBlockCoordSys& coord_sys = geometry.getBlockCoordSys(dbl[dit]);
            const ProblemDomain& domain = coord_sys.domain();
            fourthOrderCellExtrapAtDomainBdry(m_F0[dit], domain, dbl[dit]);
         }

         geometry.fillInternalGhosts(m_F0);
         geometry.divideBStarParallel(m_F0);
      }

      for (DataIterator dit(m_F0.dataIterator()); dit.ok(); ++dit) {
         a_maxwellian_dfn[dit].copy(m_F0[dit]);
      }
   }
   
   else {
      a_soln_species.computeFSavgMaxwellian(a_maxwellian_dfn);
   }
   
   for (DataIterator dit(a_delta_dfn.dataIterator()); dit.ok(); ++dit) {
      a_delta_dfn[dit].copy(a_dfn_no_bstar[dit]);
      a_delta_dfn[dit].minus(a_maxwellian_dfn[dit]);
   }
   
   m_first_call = false;
   
}


Real
GKVlasov::computeDtExplicitTI( const CFG::EField&            a_E_field,
                               const KineticSpeciesPtrVect&  a_species_vect,
                               const Real&                   a_time)
{
   CH_TIME("GKVlasov::computeDtExplicitTI");

   Real dt(BASEFAB_REAL_SETVAL);

   if ( a_E_field.supportsDivFreePhaseVel() ) {

      const CFG::LevelData<CFG::FArrayBox>& Efield_cell = a_E_field.getCellCenteredField();
      const CFG::LevelData<CFG::FArrayBox>& phi_node = a_E_field.getPhiNode();

      for (int s(0); s<a_species_vect.size(); s++) {

         const KineticSpecies& species( *(a_species_vect[s]) );
         const LevelData<FArrayBox>& dfn( species.distributionFunction() );

         LevelData<FluxBox> velocity_normals( dfn.getBoxes(), 1, IntVect::Unit );
         species.computeMappedVelocityNormals(  velocity_normals, 
                                                Efield_cell, 
                                                phi_node, 
                                                false, 
                                                PhaseGeom::FULL_VELOCITY );

         const Real UNIT_CFL(1.0);
         const PhaseGeom& geometry( species.phaseSpaceGeometry() );
         Real speciesDt(  computeMappedDtSpeciesFromNormals(  velocity_normals, 
                                                              geometry, 
                                                              UNIT_CFL ) );
         CH_assert(speciesDt >= 0);

         dt = Min( dt, speciesDt );
      }
   }
   else {

      const LevelData<FluxBox>& injected_E_field = a_E_field.getInjectedField();

      for (int s(0); s<a_species_vect.size(); s++) {

         const KineticSpecies& species( *(a_species_vect[s]) );
         const LevelData<FArrayBox>& dfn( species.distributionFunction() );
         const PhaseGeom& geometry( species.phaseSpaceGeometry() );

         IntVect ghostVect = (geometry.secondOrder()) ? IntVect::Zero : IntVect::Unit;
         if ( !m_velocity.isDefined()) {
            m_velocity.define( dfn.getBoxes(), SpaceDim, ghostVect );
         }
         species.computeMappedVelocity( m_velocity, 
                                        injected_E_field, 
                                        false, 
                                        a_time );

         const Real UNIT_CFL(1.0);
         Real speciesDt( computeMappedDtSpecies( m_velocity, geometry, UNIT_CFL ) );
         CH_assert(speciesDt >= 0);

         dt = Min( dt, speciesDt );
      }
   }

   return dt;
}

Real
GKVlasov::computeDtImExTI( const CFG::EField&            a_E_field,
                           const KineticSpeciesPtrVect&  a_species_vect,
                           const Real&                   a_time)
{
   CH_TIME("GKVlasov::computeDtImExTI");
  
   Real dt(BASEFAB_REAL_SETVAL);

   if ( a_E_field.supportsDivFreePhaseVel() ) {

      const CFG::LevelData<CFG::FArrayBox>& Efield_cell = a_E_field.getCellCenteredField();
      const CFG::LevelData<CFG::FArrayBox>& phi_node = a_E_field.getPhiNode();

      for (int s(0); s<a_species_vect.size(); s++) {

         const KineticSpecies& species( *(a_species_vect[s]) );
         const std::string& name = species.name();

         std::map<std::string, tiOpType>::const_iterator it;
         it = m_ti_op_type.find(name);
         if (it == m_ti_op_type.end()) {
           MayDay::Error("in GKVlasov::computeDtImExTI - m_ti_op_type does not contain species info.");
         }

         if (m_ti_op_type[name] == _explicit_op_) {

           const LevelData<FArrayBox>& dfn( species.distributionFunction() );
  
           LevelData<FluxBox> velocity_normals( dfn.getBoxes(), 1, IntVect::Unit );
           species.computeMappedVelocityNormals(  velocity_normals, 
                                                  Efield_cell, 
                                                  phi_node, 
                                                  false, 
                                                  PhaseGeom::FULL_VELOCITY );
  
           const Real UNIT_CFL(1.0);
           const PhaseGeom& geometry( species.phaseSpaceGeometry() );
           Real speciesDt( computeMappedDtSpeciesFromNormals( velocity_normals, geometry, UNIT_CFL ) );
           CH_assert(speciesDt >= 0);
  
           dt = Min( dt, speciesDt );
         }
      }
   }
   else {

      const LevelData<FluxBox>& injected_E_field = a_E_field.getInjectedField();
      
      for (int s(0); s<a_species_vect.size(); s++) {

         const KineticSpecies& species( *(a_species_vect[s]) );
         const std::string& name = species.name();

         std::map<std::string, tiOpType>::const_iterator it;
         it = m_ti_op_type.find(name);
         if (it == m_ti_op_type.end()) {
           MayDay::Error("in GKVlasov::computeDtImExTI - m_ti_op_type does not contain species info.");
         }

         if (m_ti_op_type[name] == _explicit_op_) {

           const LevelData<FArrayBox>& dfn( species.distributionFunction() );
           const PhaseGeom& geometry( species.phaseSpaceGeometry() );
 
           IntVect ghostVect = (geometry.secondOrder()) ? IntVect::Zero : IntVect::Unit;
           if ( !m_velocity.isDefined()) {
             m_velocity.define( dfn.getBoxes(), SpaceDim, ghostVect );
           }
           species.computeMappedVelocity( m_velocity, 
                                          injected_E_field, 
                                          false, 
                                          a_time );
            
           const Real UNIT_CFL(1.0);
           Real speciesDt( computeMappedDtSpecies( m_velocity, geometry, UNIT_CFL ) );
           CH_assert(speciesDt >= 0);
  
           dt = Min( dt, speciesDt );
         }
      }
   }

   return dt;
}


Real
GKVlasov::computeTimeScale( const CFG::EField&           a_E_field,
                            const KineticSpeciesPtrVect& a_species_vect,
                            const Real&                  a_time)
{
   CH_TIME("GKVlasov::computeTimeScale");

   Real dt(BASEFAB_REAL_SETVAL);

   if ( a_E_field.supportsDivFreePhaseVel() ) {

      const CFG::LevelData<CFG::FArrayBox>& Efield_cell = a_E_field.getCellCenteredField();
      const CFG::LevelData<CFG::FArrayBox>& phi_node = a_E_field.getPhiNode();

      for (int s(0); s<a_species_vect.size(); s++) {

         const KineticSpecies& species( *(a_species_vect[s]) );
         const LevelData<FArrayBox>& dfn( species.distributionFunction() );

         LevelData<FluxBox> velocity_normals( dfn.getBoxes(), 1, IntVect::Unit );
         species.computeMappedVelocityNormals(  velocity_normals, 
                                                Efield_cell, 
                                                phi_node, 
                                                false,
                                                PhaseGeom::FULL_VELOCITY);

         //      const Real UNIT_CFL(1.0);
         const PhaseGeom& geometry( species.phaseSpaceGeometry() );
         Real speciesDt( computeMappedTimeScaleSpeciesFromNormals( velocity_normals, geometry) );
         CH_assert(speciesDt >= 0);

         dt = Min( dt, speciesDt );
      }
   }
   else {

      const LevelData<FluxBox>& injected_E_field = a_E_field.getInjectedField();

      for (int s(0); s<a_species_vect.size(); s++) {

         const KineticSpecies& species( *(a_species_vect[s]) );
         const LevelData<FArrayBox>& dfn( species.distributionFunction() );
         const PhaseGeom& geometry( species.phaseSpaceGeometry() );
         
         IntVect ghostVect = (geometry.secondOrder()) ? IntVect::Zero : IntVect::Unit;
         if ( !m_velocity.isDefined()) {
            m_velocity.define( dfn.getBoxes(), SpaceDim, ghostVect );
         }
         species.computeMappedVelocity( m_velocity, 
                                        injected_E_field, 
                                        false, 
                                        a_time );
         
         //      const Real UNIT_CFL(1.0);
         Real speciesDt( computeMappedTimeScaleSpecies( m_velocity, geometry) );
         CH_assert(speciesDt >= 0);
         
         dt = Min( dt, speciesDt );
      }
   }

   return dt;
}



void
GKVlasov::initialize( KineticSpeciesPtrVect& soln,
                      const Real      time )
{
   MayDay::Error( "not implemented!" );
}



inline void
GKVlasov::computeFaceAverages( LevelData<FluxBox>&         a_face_data,
                               const LevelData<FArrayBox>& a_cell_data,
                               const bool                  a_second_order  ) const
{
   if ( a_second_order) {
      CellToEdge( a_cell_data, a_face_data );
   }
   else {
      fourthOrderCellToFace( a_face_data, a_cell_data );
   }
}


#undef PLOT_STABLEDT

Real
GKVlasov::computeMappedDtSpecies( const LevelData<FluxBox>& a_faceVel,
                                  const PhaseGeom&          a_geom,
                                  Real                      a_cfl )
{
   CH_TIMERS("GKVlasov::computeMappedDtSpecies");
   CH_TIMER("create_face_norm_vel",t_create_face_norm_vel); 
   CH_TIMER("fill_face_norm_vel",t_fill_face_norm_vel);
   CH_TIMER("get_cell_vol",t_get_cell_vol);

   const DisjointBoxLayout grids = a_faceVel.getBoxes();
   CH_assert(grids == a_geom.gridsFull());

   MultiBlockCoordSys* coords = a_geom.coordSysPtr();

   // Get velocities normal to cell faces and compute area-weighted normal
   // velocities -- 2nd order should be good enough.
   // Note: For axisymmetric 2-spatial-D configuration geometries, the area-weighted
   // velocities returned by the following function call contain a factor
   // of 2piR times major radius R; a consequence of the fact that this
   // factor is included in all of the metric factors.  However, the
   // physical cell volume being divided out below also contains
   // a 2piR factor (having been derived from the metrics), and therefore
   // this extra factor has no net effect (to second order).
   // (Of course this is irrrelevant in 3D)

   CH_START(t_create_face_norm_vel);
   LevelData<FluxBox> faceNormalVel(grids, 1);
   CH_STOP(t_create_face_norm_vel);

   CH_START(t_fill_face_norm_vel);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coords->whichBlock(grids[dit]);
      const PhaseBlockCoordSys* block_coord_sys
         = static_cast<const PhaseBlockCoordSys*>(coords->getCoordSys(block_number));
      RealVect face_area = block_coord_sys->getMappedFaceArea();

      const FluxBox& thisFaceVel = a_faceVel[dit];
      FluxBox& thisNormalVel = faceNormalVel[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& thisNormalVel_dir = thisNormalVel[dir];
         thisNormalVel_dir.copy(thisFaceVel[dir],dir,0,1);
         thisNormalVel_dir *= face_area[dir];
      }
   }
   CH_STOP(t_fill_face_norm_vel);

#ifdef PLOT_STABLEDT
   LevelData<FArrayBox> stableDt(grids, 1);
#endif

   // now average back to cell-centers and divide by cell volumes.
   // instead of averaging face->cell we pick
   // the max absolute value on the two opposing faces.

   CH_START(t_get_cell_vol);
   LevelData<FArrayBox> cellVolumes(grids, 1, IntVect::Zero);
   a_geom.getCellVolumes(cellVolumes);
   CH_STOP(t_get_cell_vol);

   struct {
      double val;
      int rank;
   } pair_in, pair_out;

   Real maxVelLoc = -1.;
   int maxblockLoc(-1);
   IntVect maxindLoc;
   int maxDirLoc(-1);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coords->whichBlock(grids[dit]);

      FArrayBox cellVel(grids[dit], 1);
      FArrayBox cellVelDir(grids[dit],SpaceDim);
      cellVel.setVal(0.);
      for (int dir=0; dir<SpaceDim; dir++) {
         // average face velocity on the two faces in this direction
         int faceComp = 0;
         int cellComp = dir;
         EdgeToCell(faceNormalVel[dit], faceComp,
                    cellVelDir, cellComp,
                    dir);
         cellVelDir.abs(dir,1);
         cellVel.plus(cellVelDir, dir, 0, 1);
      }

      // Divide by the cell physical volume, which is J times the computational cell volume
      cellVel.divide(cellVolumes[dit], 0, 0, 1);

      // now compute maxVelLoc on this box
      // note that this is essentially an estimate of max(vel/dx)
      Real thisMax = cellVel.norm(0,0,1);
      if (thisMax > maxVelLoc) {
         maxVelLoc = thisMax;
         if (m_time_step_diagnostics) {
            maxblockLoc = block_number;
            maxindLoc = cellVel.maxIndex();

            // Figure out which direction made the biggest contribution
            double max_cv = 0.;
            for (int dir=0; dir<SpaceDim; ++dir) {
               if (cellVelDir(maxindLoc,dir) > max_cv) {
                  max_cv = cellVelDir(maxindLoc,dir);
                  maxDirLoc = dir;
               }
            }
         }
      }

#ifdef PLOT_STABLEDT
      stableDt[dit].copy(cellVel);
#endif
   }

   RealVect X;
   if (m_time_step_diagnostics) {
      
     if (maxblockLoc < 0) MayDay::Error( "GKVlasov::computeMappedDtSpecies: time step diagnostic calculation failed, check if the number of processes used is consistent with the domain decomposition " );
     
      RealVect dx = coords->getCoordSys(maxblockLoc)->dx();
      RealVect xi = dx*maxindLoc;
      xi += 0.5*dx;
      X = coords->getCoordSys(maxblockLoc)->realCoord(xi);
   }

   Real maxVel = maxVelLoc;
   IntVect maxind = maxindLoc;
   int maxDir = maxDirLoc;
#ifdef CH_MPI
   if (m_time_step_diagnostics) {
      pair_in.val = maxVel;
      pair_in.rank = procID();
      MPI_Allreduce(&pair_in, &pair_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      maxVel = pair_out.val;
      MPI_Bcast(X.dataPtr(), SpaceDim, MPI_DOUBLE, pair_out.rank, MPI_COMM_WORLD);
      MPI_Bcast(maxind.dataPtr(), SpaceDim, MPI_INT, pair_out.rank, MPI_COMM_WORLD);
      MPI_Bcast(&maxDir, 1, MPI_INT, pair_out.rank, MPI_COMM_WORLD);
   }
   else {
      MPI_Allreduce(&maxVelLoc, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   }
#endif

   Real smallVal = 1.0e-15;
   Real newDt = 0;

   if (maxVel > smallVal)
      {
         // 2.06 factor per our JCP paper, equation (80), however approximate
         newDt = a_cfl * s_stability_bound[m_face_avg_type] / (maxVel * m_dt_dim_factor);
         //       newDt = a_cfl * 2.06 / maxVel;
      }

   if (m_time_step_diagnostics && procID()==0) {
#if PDIM==4
      cout << "Vlasov operator time step ("<< newDt <<") was limited by the velocity at (R,Z,vparallel,mu) = " << X << " and mapped coordinate = " << maxind << endl;
#endif
#if PDIM==5
      cout << "Vlasov operator time step ("<< newDt <<") was limited by the velocity at (X,Y,Z,vparallel,mu) = " << X << " and mapped coordinate = " << maxind << endl;
#endif
      cout << "The ";
      switch(maxDir)
         {
         case RADIAL_DIR:
            cout << "radial";
            break;
         case POLOIDAL_DIR:
            cout << "poloidal";
            break;
#if PDIM==5
         case TOROIDAL_DIR:
            cout << "toroidal";
            break;
#endif
         case VPARALLEL_DIR:
            cout << "vparallel";
            break;
         case MU_DIR:
            cout << "mu";
            break;
         }
      cout << " direction makes the largest contribution to the stable dt reciprocal at that point" << endl;
   }

#ifdef PLOT_STABLEDT
   a_geom.plotAtVelocityIndex( "stabledt", a_geom.vel_restrict(maxind), stableDt, 0.0);
   exit(1);
#endif

   return newDt;
}


Real
GKVlasov::computeMappedDtSpeciesFromNormals( const LevelData<FluxBox>& a_velocity_normal,
                                             const PhaseGeom&          a_geom,
                                             Real                      a_cfl )
{
   // N.B.: The velocity being passed in here is assumed to be the normal component of
   // the mapped velocity integrated over cell faces.

   const DisjointBoxLayout grids = a_velocity_normal.getBoxes();
   CH_assert(grids == a_geom.gridsFull());

#ifdef PLOT_STABLEDT
   LevelData<FArrayBox> stableDt(grids, 1);
#endif

   // now average back to cell-centers and divide by cell volumes times
   // BStarParallel.  The physical cell volume is J times the mapped
   // grid cell volume, so this division removes the J*BStarParallel
   // factor carried in the computational unknown.  In the non-divergence-free
   // velocity option, the velocities already include a 1/BStarParallel
   // factor

   // instead of averaging face->cell we pick
   // the max absolute value on the two opposing faces.
   LevelData<FArrayBox> cellVolumesBStarPar(grids, 1, IntVect::Zero);
   a_geom.getCellVolumes(cellVolumesBStarPar);
   a_geom.multBStarParallel(cellVolumesBStarPar);

   MultiBlockCoordSys* coords = a_geom.coordSysPtr();

   struct {
      double val;
      int rank;
   } pair_in, pair_out;

   Real maxVelLoc = 0.;
   int maxblockLoc(-1);
   IntVect maxindLoc;
   int maxDirLoc(-1);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coords->whichBlock(grids[dit]);

      FArrayBox cellVel(grids[dit], 1);
      FArrayBox cellVelDir(grids[dit],SpaceDim);
      cellVel.setVal(0.);
      for (int dir=0; dir<SpaceDim; dir++) {
         // average face velocity on the two faces in this direction
         int faceComp = 0;
         int cellComp = dir;
         EdgeToCell(a_velocity_normal[dit], faceComp,
                    cellVelDir, cellComp,
                    dir);
         cellVelDir.abs(dir,1);
         cellVel.plus(cellVelDir, dir, 0, 1);
      }

      cellVel.divide(cellVolumesBStarPar[dit], 0, 0, 1);

      // now compute maxVelLoc on this box
      // note that this is essentially an estimate of max(vel/dx)
      Real thisMax = cellVel.norm(0,0,1);
      if (thisMax > maxVelLoc) {
         maxVelLoc = thisMax;
         if (m_time_step_diagnostics) {
            maxblockLoc = block_number;
            maxindLoc = cellVel.maxIndex();

            // Figure out which direction made the biggest contribution
            double max_cv = 0.;
            for (int dir=0; dir<SpaceDim; ++dir) {
               if (cellVelDir(maxindLoc,dir) > max_cv) {
                  max_cv = cellVelDir(maxindLoc,dir);
                  maxDirLoc = dir;
               }
            }
         }
      }

#ifdef PLOT_STABLEDT
      stableDt[dit].copy(cellVel);
#endif
   }

   RealVect X;
   if (m_time_step_diagnostics) {

     if (maxblockLoc < 0) MayDay::Error( "GKVlasov::computeMappedDtSpeciesFromNormals: time step diagnostic calculation failed, check if the number of processes used is consistent with the domain decomposition");

      RealVect dx = coords->getCoordSys(maxblockLoc)->dx();
      RealVect xi = dx*maxindLoc;
      xi += 0.5*dx;
      X = coords->getCoordSys(maxblockLoc)->realCoord(xi);
   }

   Real maxVel = maxVelLoc;
   IntVect maxind = maxindLoc;
   int maxDir = maxDirLoc;
#ifdef CH_MPI
   if (m_time_step_diagnostics) {
      pair_in.val = maxVel;
      pair_in.rank = procID();
      MPI_Allreduce(&pair_in, &pair_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      maxVel = pair_out.val;
      MPI_Bcast(X.dataPtr(), SpaceDim, MPI_DOUBLE, pair_out.rank, MPI_COMM_WORLD);
      MPI_Bcast(maxind.dataPtr(), SpaceDim, MPI_INT, pair_out.rank, MPI_COMM_WORLD);
      MPI_Bcast(&maxDir, 1, MPI_INT, pair_out.rank, MPI_COMM_WORLD);
   }
   else {
      MPI_Allreduce(&maxVelLoc, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   }
#endif

   Real smallVal = 1.0e-15;
   Real newDt = 0;

   if (maxVel > smallVal)
      {
         // 2.06 factor per our JCP paper, equation (80), however approximate
         newDt = a_cfl * s_stability_bound[m_face_avg_type] / (maxVel * m_dt_dim_factor);
         //       newDt = a_cfl * 2.06 / maxVel;
      }

   if (m_time_step_diagnostics && procID()==0) {
#if PDIM==4
      cout << "Vlasov operator time step was limited by the velocity at (R,Z,vparallel,mu) = " << X << " and mapped coordinate = " << maxind << endl;
#endif
#if PDIM==5
      cout << "Vlasov operator time step was limited by the velocity at (R,phi,Z,vparallel,mu) = " << X << " and mapped coordinate = " << maxind << endl;
#endif
      cout << "The ";
      switch(maxDir)
         {
         case RADIAL_DIR:
            cout << "radial";
            break;
         case POLOIDAL_DIR:
            cout << "poloidal";
            break;
#if PDIM==5
         case TOROIDAL_DIR:
            cout << "toroidal";
            break;
#endif
         case VPARALLEL_DIR:
            cout << "vparallel";
            break;
         case MU_DIR:
            cout << "mu";
            break;
         }
      cout << " direction makes the largest contribution to the stable dt reciprocal at that point" << endl;
   }

#ifdef PLOT_STABLEDT
   a_geom.plotAtVelocityIndex( "stabledt", a_geom.vel_restrict(maxind), stableDt );
   exit(1);
#endif

   return newDt;
}

Real
GKVlasov::computeMappedTimeScaleSpecies( const LevelData<FluxBox>& a_faceVel,
                                         const PhaseGeom&          a_geom )
{
   const DisjointBoxLayout grids = a_faceVel.getBoxes();
   CH_assert(grids == a_geom.gridsFull());

   MultiBlockCoordSys* coords = a_geom.coordSysPtr();

   // Get velocities normal to cell faces and compute area-weighted normal
   // velocities -- 2nd order should be good enough.
   // Note: For axisymmetric 2-spatial-D configuration geometries, the area-weighted
   // velocities returned by the following function call contain a factor
   // of 2piR times major radius R; a consequence of the fact that this
   // factor is included in all of the metric factors.  However, the
   // physical cell volume being divided out below also contains
   // a 2piR factor (having been derived from the metrics), and therefore
   // this extra factor has no net effect (to second order).
   // (Of course this is irrrelevant in 3D)

   LevelData<FluxBox> faceNormalVel(grids, 1);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coords->whichBlock(grids[dit]);
      const PhaseBlockCoordSys* block_coord_sys
         = static_cast<const PhaseBlockCoordSys*>(coords->getCoordSys(block_number));
      RealVect face_area = block_coord_sys->getMappedFaceArea();

      const FluxBox& thisFaceVel = a_faceVel[dit];
      FluxBox& thisNormalVel = faceNormalVel[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& thisNormalVel_dir = thisNormalVel[dir];
         thisNormalVel_dir.copy(thisFaceVel[dir],dir,0,1);
         thisNormalVel_dir *= face_area[dir];
      }
   }

   // now average back to cell-centers and divide by cell volumes
   // instead of averaging face->cell we pick
   // the max absolute value on the two opposing faces.
   LevelData<FArrayBox> cellVolumes(grids, 1, IntVect::Zero);
   a_geom.getCellVolumes(cellVolumes);

   struct {
      double val;
      int rank;
   } pair_in, pair_out;

   Real maxVelLoc = -1.;
   int maxblockLoc(-1);
   IntVect maxindLoc;
   int maxDirLoc(-1);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coords->whichBlock(grids[dit]);

      FArrayBox cellVel(grids[dit], 1);
      FArrayBox cellVelDir(grids[dit],SpaceDim);
      cellVel.setVal(0.);
      for (int dir=0; dir<SpaceDim; dir++) {
         // average face velocity on the two faces in this direction
         int faceComp = 0;
         int cellComp = dir;
         EdgeToCell(faceNormalVel[dit], faceComp,
                    cellVelDir, cellComp,
                    dir);
         cellVelDir.abs(dir,1);
         cellVel.plus(cellVelDir, dir, 0, 1);
      }

      // Divide by the cell physical volume, which is J times the computational cell volume
      cellVel.divide(cellVolumes[dit], 0, 0, 1);

      // now compute maxVelLoc on this box
      // note that this is essentially an estimate of max(vel/dx)
      Real thisMax = cellVel.norm(0,0,1);
      if (thisMax > maxVelLoc) {
         maxVelLoc = thisMax;
         if (m_time_step_diagnostics) {
            maxblockLoc = block_number;
            maxindLoc = cellVel.maxIndex();

            // Figure out which direction made the biggest contribution
            double max_cv = 0.;
            for (int dir=0; dir<SpaceDim; ++dir) {
               if (cellVelDir(maxindLoc,dir) > max_cv) {
                  max_cv = cellVelDir(maxindLoc,dir);
                  maxDirLoc = dir;
               }
            }
         }
      }
   }

   RealVect X;
   if (m_time_step_diagnostics) {
      RealVect dx = coords->getCoordSys(maxblockLoc)->dx();
      RealVect xi = dx*maxindLoc;
      xi += 0.5*dx;
      X = coords->getCoordSys(maxblockLoc)->realCoord(xi);
   }

   Real maxVel = maxVelLoc;
   IntVect maxind = maxindLoc;
   int maxDir = maxDirLoc;
#ifdef CH_MPI
   if (m_time_step_diagnostics) {
      pair_in.val = maxVel;
      pair_in.rank = procID();
      MPI_Allreduce(&pair_in, &pair_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      maxVel = pair_out.val;
      MPI_Bcast(X.dataPtr(), SpaceDim, MPI_DOUBLE, pair_out.rank, MPI_COMM_WORLD);
      MPI_Bcast(maxind.dataPtr(), SpaceDim, MPI_INT, pair_out.rank, MPI_COMM_WORLD);
      MPI_Bcast(&maxDir, 1, MPI_INT, pair_out.rank, MPI_COMM_WORLD);
   }
   else {
      MPI_Allreduce(&maxVelLoc, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   }
#endif

   Real smallVal = 1.0e-15;
   Real newDt = 0;
   if (maxVel > smallVal) newDt = 1.0 / maxVel;

   return newDt;
}



Real
GKVlasov::computeMappedTimeScaleSpeciesFromNormals( const LevelData<FluxBox>& a_velocity_normal,
                                                    const PhaseGeom&          a_geom )
{
   // N.B.: The velocity being passed in here is assumed to be the normal component of
   // the mapped velocity integrated over cell faces.

   const DisjointBoxLayout grids = a_velocity_normal.getBoxes();
   CH_assert(grids == a_geom.gridsFull());

   // now average back to cell-centers and divide by cell volumes times
   // BStarParallel.  The physical cell volume is J times the mapped
   // grid cell volume, so this division removes the J*BStarParallel
   // factor carried in the computational unknown.  In the non-divergence-free
   // velocity option, the velocities already include a 1/BStarParallel
   // factor

   // instead of averaging face->cell we pick
   // the max absolute value on the two opposing faces.
   LevelData<FArrayBox> cellVolumesBStarPar(grids, 1, IntVect::Zero);
   a_geom.getCellVolumes(cellVolumesBStarPar);
   a_geom.multBStarParallel(cellVolumesBStarPar);

   MultiBlockCoordSys* coords = a_geom.coordSysPtr();

   struct {
      double val;
      int rank;
   } pair_in, pair_out;

   Real maxVelLoc = 0.;
   int maxblockLoc(-1);
   IntVect maxindLoc;
   int maxDirLoc(-1);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coords->whichBlock(grids[dit]);

      FArrayBox cellVel(grids[dit], 1);
      FArrayBox cellVelDir(grids[dit],SpaceDim);
      cellVel.setVal(0.);
      for (int dir=0; dir<SpaceDim; dir++) {
         // average face velocity on the two faces in this direction
         int faceComp = 0;
         int cellComp = dir;
         EdgeToCell(a_velocity_normal[dit], faceComp,
                    cellVelDir, cellComp,
                    dir);
         cellVelDir.abs(dir,1);
         cellVel.plus(cellVelDir, dir, 0, 1);
      }

      cellVel.divide(cellVolumesBStarPar[dit], 0, 0, 1);

      // now compute maxVelLoc on this box
      // note that this is essentially an estimate of max(vel/dx)
      Real thisMax = cellVel.norm(0,0,1);
      if (thisMax > maxVelLoc) {
         maxVelLoc = thisMax;
         if (m_time_step_diagnostics) {
            maxblockLoc = block_number;
            maxindLoc = cellVel.maxIndex();

            // Figure out which direction made the biggest contribution
            double max_cv = 0.;
            for (int dir=0; dir<SpaceDim; ++dir) {
               if (cellVelDir(maxindLoc,dir) > max_cv) {
                  max_cv = cellVelDir(maxindLoc,dir);
                  maxDirLoc = dir;
               }
            }
         }
      }
   }

   RealVect X;
   if (m_time_step_diagnostics) {
      RealVect dx = coords->getCoordSys(maxblockLoc)->dx();
      RealVect xi = dx*maxindLoc;
      xi += 0.5*dx;
      X = coords->getCoordSys(maxblockLoc)->realCoord(xi);
   }

   Real maxVel = maxVelLoc;
   IntVect maxind = maxindLoc;
   int maxDir = maxDirLoc;
#ifdef CH_MPI
   if (m_time_step_diagnostics) {
      pair_in.val = maxVel;
      pair_in.rank = procID();
      MPI_Allreduce(&pair_in, &pair_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      maxVel = pair_out.val;
      MPI_Bcast(X.dataPtr(), SpaceDim, MPI_DOUBLE, pair_out.rank, MPI_COMM_WORLD);
      MPI_Bcast(maxind.dataPtr(), SpaceDim, MPI_INT, pair_out.rank, MPI_COMM_WORLD);
      MPI_Bcast(&maxDir, 1, MPI_INT, pair_out.rank, MPI_COMM_WORLD);
   }
   else {
      MPI_Allreduce(&maxVelLoc, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   }
#endif

   Real smallVal = 1.0e-15;
   Real newDt = 0;
   if (maxVel > smallVal) newDt = 1.0 / maxVel;

   return newDt;
}



void
GKVlasov::applyMappedLimiter( LevelData<FluxBox>&         a_facePhi,
                              const LevelData<FArrayBox>& a_cellPhi,
                              const LevelData<FluxBox>&   a_faceVel,
                              const PhaseGeom&            a_geom,
                              const bool                  a_mapped_velocity )
{
   int nComp = a_facePhi.nComp();

   // this specifies the number of ghost faces we will need in
   // each transverse and normal direction
   int transverseGrow = 2;
   int normalGrow = 0;

   // use C value from Colella and Sekora
   Real limiterConst =  1.25;

   // may need to do exchange on facePhi
   a_facePhi.exchange();

   const DisjointBoxLayout& grids = a_facePhi.getBoxes();

   // in order to do upwinding, need normal velocities in
   // computational space. computing Fourth-order face
   // averages just in order to upwinding is over-working the
   // issue, but it's the most convenient thing to do here.
   LevelData<FluxBox> normalVel(grids, 1, a_faceVel.ghostVect());
   if ( a_mapped_velocity ) {
      // The passed in velocity is in the mapped frame
      int nvcomp = a_faceVel.nComp();
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const PhaseBlockCoordSys& block_coord_sys = a_geom.getBlockCoordSys(grids[dit]);
         RealVect face_area = block_coord_sys.getMappedFaceArea();

         const FluxBox& thisFaceVel( a_faceVel[dit] );
         FluxBox& thisNormalVel( normalVel[dit] );

         for (int dir=0; dir<SpaceDim; ++dir) {
            int src_comp = (nvcomp==1)? 0: dir;
            FArrayBox& thisNormalVel_dir( thisNormalVel[dir] );
            thisNormalVel_dir.copy(thisFaceVel[dir],src_comp,0,1);
            thisNormalVel_dir *= face_area[dir];
         }
      }
   }
   else {
      // The passed in velocity is in the physical frame (default)
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         normalVel[dit].setVal( 0.0 );
      }
      a_geom.computeMetricTermProductAverage( normalVel, a_faceVel, false );
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = a_geom.getBlockCoordSys(grids[dit]);
      const RealVect& dx = block_coord_sys.dx();

      const Box& gridBox = grids[dit];
      Box lapBox(gridBox);
      // need two laplacian ghost cells for this limiter
      int lapBoxGrow = max(normalGrow+2, transverseGrow);
      lapBox.grow(lapBoxGrow);
      FArrayBox ccLaplacians(lapBox, nComp);
      FluxBox& thisFacePhi = a_facePhi[dit];
      const FArrayBox& thisCellPhi = a_cellPhi[dit];
      const FluxBox& thisNormalVel = normalVel[dit];
      // check to be sure that we have enough phi to compute
      // all of these Laplacians
      {
         Box LapPhiBox(lapBox);
         LapPhiBox.grow(1);
         CH_assert(thisCellPhi.box().contains(LapPhiBox));
      }
      
      for (int dir=0; dir<SpaceDim; dir++)
         {
            // box of valid edges for this grid
            Box faceBox(gridBox);
            // for 4th order, we need extra faces in the transverse
            // directions. (for mapped-grid case, we actually need
            // _2_ transverse faces

            // need extra faces in the tangential directions in order
            // to handle 4th-order multiplies
            faceBox.grow(transverseGrow);
            // however, need different number of extra faces in normal dir
            faceBox.grow(dir,normalGrow - transverseGrow);

            faceBox.surroundingNodes(dir);
            FArrayBox& thisFacePhiDir = thisFacePhi[dir];
            const FArrayBox& thisNormalVelDir = thisNormalVel[dir];
            {
               // context for computing limited face values
               Box grownFaceBox(faceBox);
               // need an extra face's worth of the FC laplacians
               grownFaceBox.grow(dir,1);

               FArrayBox centeredLaplacian(grownFaceBox, nComp);

               // compute centered Laplacian
               centeredLaplacian.setVal(0.0);
               FORT_CENTEREDLAPLACIAN(CHF_FRA(centeredLaplacian),
                                      CHF_CONST_FRA(thisFacePhiDir),
                                      CHF_CONST_FRA(thisCellPhi),
                                      CHF_BOX(grownFaceBox),
                                      CHF_CONST_INT(dir),
                                      CHF_CONST_REAL(dx[dir]));

               // compute cell-centered Laplacians
               ccLaplacians.setVal(0.0);
               FORT_CCLAPLACIAN(CHF_FRA(ccLaplacians),
                                CHF_CONST_FRA(thisCellPhi),
                                CHF_BOX(lapBox),
                                CHF_CONST_INT(dir),
                                CHF_CONST_REAL(dx[dir]));

               // now compute limited face value
               FORT_LIMITFACEVALUES(CHF_FRA(thisFacePhiDir),
                                    CHF_CONST_FRA(thisCellPhi),
                                    CHF_CONST_FRA(centeredLaplacian),
                                    CHF_CONST_FRA(ccLaplacians),
                                    CHF_BOX(grownFaceBox),
                                    CHF_CONST_INT(dir),
                                    CHF_CONST_REAL(dx[dir]),
                                    CHF_CONST_REAL(limiterConst));
               // end context for computing limited face values
            } // (centeredLaplacian goes out of scope)

            // now compute parabolic interpolants
            // storage for cell-centered D^2a
            Box lapBoxDir(gridBox);
            // need transverse faces for 4th order
            lapBoxDir.grow(transverseGrow);
            // need this to be grown by one in normal dir
            lapBoxDir.grow(dir,normalGrow-transverseGrow+1);

            FArrayBox D2a(lapBoxDir,nComp);
            FArrayBox D2aLim(lapBoxDir,nComp);

            // initialize D2a to be -a_6/(h^2)
            // first compute a_6...
            FORT_COMPUTEA6(CHF_FRA(D2a),
                           CHF_CONST_FRA(thisCellPhi),
                           CHF_CONST_FRA(thisFacePhiDir),
                           CHF_BOX(lapBoxDir),
                           CHF_CONST_INT(dir));

            // then multiply by -2/(h^2)
            // (dfm 5/5/09 -- was missing the factor of 2)
            Real mult = -2.0/(dx[dir]*dx[dir]);
            D2a *= mult;

            // now limit D2a w/r/t other 2nd-derivatives
            FORT_LIMITD2A(CHF_FRA(D2aLim),
                          CHF_CONST_FRA(D2a),
                          CHF_CONST_FRA(ccLaplacians),
                          CHF_BOX(lapBoxDir),
                          CHF_CONST_INT(dir),
                          CHF_CONST_REAL(limiterConst));

            // storage for left and right interpolants
            //   note that while leftPhi and rightPhi are associated
            //   with face indicies, it is more natural to compute the
            //   PPM limited face values on a per-cell basis.  Thus,
            //   each of these FArrayBoxes is one cell longer than it
            //   need to be
            Box growBox( faceBox );
            growBox.grow(dir,1);
            FArrayBox leftPhi(growBox, nComp);
            FArrayBox rightPhi(growBox, nComp);

            // We operate on the cells of the domain plus one ghost on each
            // end in the current direction
            Box cellBox( growBox );
            cellBox.growHi(dir,-1);
            FORT_LEFTRIGHTSTATES(CHF_FRA(leftPhi),
                                 CHF_FRA(rightPhi),
                                 CHF_CONST_FRA(thisFacePhiDir),
                                 CHF_CONST_FRA(thisCellPhi),
                                 CHF_CONST_FRA(D2aLim),
                                 CHF_CONST_FRA(D2a),
                                 CHF_BOX(cellBox),
                                 CHF_CONST_INT(dir));

            // need to do something about boundary values here

            // now pick upwind state
            // note that normalVel only has a single component, so we
            // use 0 as the component argument for the CHF_FRA1 macro
            // as a debugging check, setThisFacePhiDir to a bogus value
            // first
            thisFacePhiDir.setVal(1000000000);
            Box overlap( faceBox & thisNormalVelDir.box() );
            FORT_SELECTUPWIND(CHF_FRA(thisFacePhiDir),
                              CHF_CONST_FRA(leftPhi),
                              CHF_CONST_FRA(rightPhi),
                              CHF_CONST_FRA1(thisNormalVelDir,0),
                              CHF_BOX(overlap));

         } // end loop over directions

   } // end loop over grids
}


double
GKVlasov::globalMax( const double a_data ) const
{
   double global_max;

#ifdef CH_MPI
   double local_data = a_data;
   MPI_Allreduce(&local_data, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
   global_max = a_data;
#endif

   return global_max;
}


void 
GKVlasov::defineMultiPhysicsPC(std::vector<Preconditioner<ODEVector,AppCtxt>*>& a_pc,
                               std::vector<DOFList>&                            a_dof_list,
                               const KineticSpeciesPtrVect&                     a_soln,
                               const GlobalDOFKineticSpeciesPtrVect&            a_gdofs,
                               const ODEVector&                                 a_x,
                               void*                                            a_ops,
                               const std::string&                               a_out_string,
                               const std::string&                               a_opt_string,
                               bool                                             a_im,
                               const int                                        a_id )
{
  for (int species(0); species<a_soln.size(); species++) {
    const KineticSpecies&           soln_species(*(a_soln[species]));
    const GlobalDOFKineticSpecies&  gdofs_species(*(a_gdofs[species]));
    const std::string               species_name(soln_species.name());

    defineBlockPC(  a_pc, 
                    a_dof_list, 
                    a_x, 
                    a_ops, 
                    a_out_string, 
                    a_opt_string, 
                    a_im, 
                    soln_species, 
                    gdofs_species, 
                    species,
                    a_id );
  }
}


void 
GKVlasov::defineBlockPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>& a_pc,
                         std::vector<DOFList>&                            a_dof_list,
                         const ODEVector&                                 a_X,
                         void*                                            a_ops,
                         const std::string&                               a_out_string,
                         const std::string&                               a_opt_string,
                         bool                                             a_im,
                         const KineticSpecies&                            a_species,
                         const GlobalDOFKineticSpecies&                   a_global_dofs,
                         const int                                        a_species_index,
                         const int                                        a_id )
{
  const std::string& name = a_species.name();

  {
    std::map<std::string, tiOpType>::const_iterator it;
    it = m_ti_op_type.find(name);
    if (it == m_ti_op_type.end()) {
      MayDay::Error("in GKVlasov::defineBlockPC - m_ti_op_type does not contain species info.");
    }
  }

  if ( a_im && (m_ti_op_type[name] == _implicit_op_) ) {
  
    CH_assert(a_pc.size() == a_dof_list.size());

    if (!procID()) {
      std::cout << "  Kinetic Species " << a_species_index
                << " - Vlasov term: "
                << " creating " << _GK_VLASOV_PC_ << " preconditioner "
                << " (index = " << a_pc.size() << ").\n";
    }

    Preconditioner<ODEVector,AppCtxt> *pc;
    pc = new GKVlasovPreconditioner<ODEVector,AppCtxt>;
    pc->setSysID(a_id);

    m_my_pc_idx[name] = a_pc.size();

    std::string out_string = pp_name + std::string(".") + name;
    std::string opt_string = pp_name + std::string(".") + name;
    dynamic_cast<GKVlasovPreconditioner<ODEVector,AppCtxt>*>
       (pc)->define(a_X, a_ops, out_string, opt_string, a_im);
    dynamic_cast<GKVlasovPreconditioner<ODEVector,AppCtxt>*>
       (pc)->speciesIndex(a_species_index);

    std::string precond_prefix = string(pp_name) + ".precond";
    ParmParse pp_precond(precond_prefix.c_str());
    dynamic_cast<GKVlasovPreconditioner<ODEVector,AppCtxt>*>
       (pc)->defineAMGPC(pp_precond, a_species.phaseSpaceGeometry());

    a_pc.push_back(pc);

    {
      DOFList dof_list(0);
  
      const LevelData<FArrayBox>& soln_dfn(a_species.distributionFunction());
      const DisjointBoxLayout&    grids   (soln_dfn.disjointBoxLayout());
      const int                   n_comp  (soln_dfn.nComp());
      const LevelData<FArrayBox>& pMapping(a_global_dofs.data()); 
    
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
      a_dof_list.push_back(dof_list);
    }

  }

  return;
}


void 
GKVlasov::updateMultiPhysicsPC(std::vector<Preconditioner<ODEVector,AppCtxt>*>& a_pc,
                               const KineticSpeciesPtrVect&                     a_soln,
                               const CFG::LevelData<CFG::FArrayBox>&            a_phi,
                               const CFG::EField&                               a_E_field,
                               const GlobalDOFKineticSpeciesPtrVect&            a_gdofs,
                               const Real                                       a_time,
                               const int                                        a_step,
                               const int                                        a_stage,
                               const Real                                       a_shift,
                               const bool                                       a_im )
{
   for (int species(0); species<a_soln.size(); species++) {
      const KineticSpecies&           soln_species(*(a_soln[species]));
      const GlobalDOFKineticSpecies&  gdofs_species(*(a_gdofs[species]));
      const std::string               species_name(soln_species.name());

      updateBlockPC( a_pc, 
                     soln_species, 
                     a_phi, 
                     a_E_field, 
                     gdofs_species, 
                     a_time, 
                     a_step,
                     a_stage,
                     a_shift, 
                     a_im, 
                     species);
   }
}


void 
GKVlasov::updateBlockPC( std::vector<Preconditioner<ODEVector,AppCtxt>*>& a_pc,
                         const KineticSpecies&                            a_species,
                         const CFG::LevelData<CFG::FArrayBox>&            a_phi,
                         const CFG::EField&                               a_E_field,
                         const GlobalDOFKineticSpecies&                   a_global_dofs,
                         const Real                                       a_time,
                         const int                                        a_step,
                         const int                                        a_stage,
                         const Real                                       a_shift,
                         const bool                                       a_im,
                         const int                                        a_species_index )
{
   CH_TIMERS("GKVlasov::updateBlockPC");
   const LevelData<FArrayBox>& dfn( a_species.distributionFunction() );

   const std::string& name = a_species.name();

   {
     std::map<std::string, tiOpType>::const_iterator it;
     it = m_ti_op_type.find(name);
     if (it == m_ti_op_type.end()) {
       MayDay::Error("in GKVlasov::updateBlockPC - m_ti_op_type does not contain species info.");
     }
   }

   {
     std::map<std::string, int>::iterator it;
     it = m_my_pc_idx.find(name);
     if ( (it == m_my_pc_idx.end()) && (m_ti_op_type[name] == _implicit_op_) ) {
        MayDay::Error("in GKVlasov::updateBlockPC - m_my_pc_idx does not contain species info.");
     }
   }

   if ( a_im && (m_ti_op_type[name] == _implicit_op_) ) {

     /* see if this update can be skipped */
     if (a_step >= 0) {
      if (a_step%m_update_pc_freq != 0) return;
      if ((a_stage >= 0) && m_update_pc_skip_stage) return; 
      if ((a_stage < 0) && (!m_update_pc_skip_stage)) return;
     }

     CH_assert(m_my_pc_idx[name] >= 0);
     CH_assert(a_pc.size() > m_my_pc_idx[name]);
 
     if (!procID()) {
       std::cout << "    ==> Updating " << _GK_VLASOV_PC_ << " preconditioner "
                 << " for Vlasov term of kinetic species " << a_species_index << ".\n";
     }

     CH_TIMER("compute_velocity", t_compute_velocity);
     CH_TIMER("construct_matrix", t_construct_matrix);
     const DisjointBoxLayout& dbl( dfn.getBoxes() );
     const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
  
     IntVect ghostVect = (geometry.secondOrder()) ? IntVect::Zero : IntVect::Unit;
     if ( !m_velocity.isDefined()) {
        m_velocity.define( dbl, SpaceDim, ghostVect );
     }
     // Compute full gyrokinetic velocity
     CH_START(t_compute_velocity);
     if (a_species.isGyrokinetic() && m_include_FLR_effects) {
  
        LevelData<FluxBox> gyroaveraged_E_field;
        getEFieldwFLR(gyroaveraged_E_field, a_species, a_phi, a_E_field);
        a_species.computeVelocity(  m_velocity, 
                                    gyroaveraged_E_field, 
                                    true,
                                    PhaseGeom::FULL_VELOCITY, 
                                    0.);
     } else {
        a_species.computeVelocity(  m_velocity, 
                                    a_E_field.getInjectedField(), 
                                    false,
                                    PhaseGeom::FULL_VELOCITY, 
                                    0. );
     }
     CH_STOP(t_compute_velocity);

     //    const MultiBlockCoordSys* coord_sys = geometry.coordSysPtr();
     const DisjointBoxLayout& grids = geometry.gridsFull();
 
     int coupling_comps = 2*(CFG_DIM+1);
     IntVect coupling_ghosts = IntVect::Unit - BASISV(MU_DIR);  // No coupling in MU
     LevelData<BaseFab<Vector<IntVect> > > couplings(grids, coupling_comps, coupling_ghosts);
     LevelData<BaseFab<Vector<Real> > > weights(grids, coupling_comps, coupling_ghosts);
 
     constructStencils( m_velocity, geometry, couplings, weights);
 
     GKVlasovAMG* this_pc = dynamic_cast<GKVlasovPreconditioner<ODEVector,AppCtxt>*>
                              (a_pc[m_my_pc_idx[name]])->getAMGPC();
 
     CH_START(t_construct_matrix);
     this_pc->constructMatrix(m_precond_face_avg_type, a_shift, couplings, weights);
     CH_STOP(t_construct_matrix);
 
     if ( m_precond_build_test ) {
        testPC(this_pc, geometry, m_velocity, a_shift);
     }

   }

   return;
}


void GKVlasov::testPC( const GKVlasovAMG*         a_pc,
                       const PhaseGeom&           a_geometry,
                       const LevelData<FluxBox>&  a_velocity,
                       const Real                 a_shift )
{
#if CFG_DIM == 3
   const CFG::MagGeom& mag_geom = a_geometry.magGeom();
   if ( typeid(*mag_geom.getCoordSys()) == typeid(CFG::SingleNullCoordSys) &&
        m_precond_face_avg_type == "uw3" ) {
      // The extrapolation at outflow boundaries does not work at sawteeth.  A getBoundaryBoxes()
      // member needs to be added to PhaseGeom (like MagGeom's) in order to finish this.
      MayDay::Error("GKVlasov::testPC() is not fully implemented for single null with UW3");
   }
#endif

   if ( a_geometry.divFreeVelocity() ) {
      MayDay::Error("GKVlasov::testPC() is not yet implemented for the divergence-free velocity option");
   }

   if ( m_subtract_maxwellian ) {
      MayDay::Error("GKVlasov::testPC() is not yet implemented for the delta-f option");
   }

   const DisjointBoxLayout& grids = a_geometry.gridsFull();

   // Get the normal component of the mapped velocity for use with UW3
   LevelData<FluxBox> normal_vel(grids, 1, IntVect::Zero);
   if (m_precond_face_avg_type == "uw3") {
      a_geometry.computeMetricTermProductAverage( normal_vel, a_velocity, false );
   }

   LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
   a_geometry.getCellVolumes(volume);

   // Construct a random vector between -1 and 1 in the computational frame

   LevelData<FArrayBox> vec_comp(grids, 1, IntVect::Zero);

   srand(time(NULL));
   for (DataIterator dit(grids); dit.ok(); ++dit) {
#if 1
      for (BoxIterator bit(grids[dit]); bit.ok(); ++bit) {
         IntVect iv = bit();
         vec_comp[dit](iv,0) = (double)rand() / (double)RAND_MAX;
      }
#else
      const PhaseBlockCoordSys& block_coord_sys = a_geometry.getBlockCoordSys(grids[dit]);
      double mapped_volume = block_coord_sys.getMappedCellVolume();

      // Set to J
      vec_comp[dit].copy(volume[dit]);
      vec_comp[dit] /= mapped_volume;
#endif
   }
   a_geometry.multBStarParallel(vec_comp);

   MultiBlockCoordSys* coord_sys_ptr = a_geometry.coordSysPtr();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = coord_sys_ptr->boundaries();
      
   // Convert to physical coordinates with ghosts.  Fill the internal ghosts,
   // set the physical ghosts at inflow boundaries to zero and fill the
   // physical ghosts at outflow boundaries using second-order extrapolation.

   LevelData<FArrayBox> vec_phys(grids, 1, 2*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box& box = grids[dit];

      const PhaseBlockCoordSys& block_coord_sys = a_geometry.getBlockCoordSys(box);
      double mapped_volume = block_coord_sys.getMappedCellVolume();

      vec_phys[dit].setVal(0.);
      vec_phys[dit].copy(vec_comp[dit],box);  // no ghosts copied here; remain zero

      // Divide J
      vec_phys[dit] *= mapped_volume;
      vec_phys[dit] /= volume[dit];

      if (m_precond_face_avg_type == "uw3") {  // Second-order extrapolation on outflow physical boundary
         int block_number = coord_sys_ptr->whichBlock(box);
         const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
         const Box& block_box = block_coord_sys.domain().domainBox();

         for (int dir=0; dir<CFG_DIM+1; ++dir) {
            if ( this_block_boundaries[dir].isDomainBoundary() ) {
               Box lo_ghost_box = adjCellLo(block_box,dir,1);
               Box grown_box = box;
               grown_box.grow(dir,1);
               Box overlap = grown_box & lo_ghost_box;
               
               if ( overlap.ok() ) {
                  for (BoxIterator bit(overlap); bit.ok(); ++bit) {
                     IntVect iv  = bit();
                     IntVect iv1 = iv;  iv1[dir]++;
                     IntVect iv2 = iv1; iv2[dir]++;

                     if ( normal_vel[dit][dir](iv1,0) <= 0. ) {
                        vec_phys[dit](iv,0) = 2.*vec_phys[dit](iv1,0) - vec_phys[dit](iv2,0);
                     }
                  }
               }
            }

            if ( this_block_boundaries[dir + SpaceDim].isDomainBoundary() ) {
               Box hi_ghost_box = adjCellHi(block_box,dir,1);
               Box grown_box = box;
               grown_box.grow(dir,1);
               Box overlap = grown_box & hi_ghost_box;

               if ( overlap.ok() ) {
                  for (BoxIterator bit(overlap); bit.ok(); ++bit) {
                     IntVect iv  = bit();
                     IntVect iv1 = iv;  iv1[dir]--;
                     IntVect iv2 = iv1; iv2[dir]--;

                     if ( normal_vel[dit][dir](iv,0) >= 0. ) {
                        vec_phys[dit](iv,0) = 2.*vec_phys[dit](iv1,0) - vec_phys[dit](iv2,0);
                     }
                  }
               }
            }

         }
      }
   }

   a_geometry.fillInternalGhosts(vec_phys);

   // Compute the matrix vector product.  N.B.: the matrix operates on
   // vectors in the physical frame

   LevelData<FArrayBox> matvec(grids, 1, IntVect::Zero);
   a_pc->multiplyMatrix(vec_phys, matvec);

   // Compute the negative flux divergence divided by the mapped cell volume

   LevelData<FArrayBox> pc_operator(grids, 1, IntVect::Zero);
   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Zero);

   // Compute Vlasov_flux[dfn]
   computeAdvectionFluxPhys( vec_phys, a_velocity, flux, a_geometry );
   
   const bool OMIT_NT(false);
   a_geometry.mappedGridDivergence( pc_operator, flux, OMIT_NT );

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = a_geometry.getBlockCoordSys(grids[dit]);
      double fac( -1.0 / block_coord_sys.getMappedCellVolume() );
      pc_operator[dit].mult( fac );
   }

   // Negate the negative flux divergence and add the shifted identity contribution
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      pc_operator[dit].negate();
      pc_operator[dit].plus(vec_comp[dit], a_shift);
   }

   // Compute the difference between the matrix vector product and the evaluated
   // preconditioner operator

   LevelData<FArrayBox> diff(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      diff[dit].copy(matvec[dit]);
      diff[dit] -= pc_operator[dit];
      diff[dit].abs();
   }

#if 0
   a_geometry.plotAtVelocityIndex("pc_operator", VEL::IntVect::Zero, pc_operator, 0.);
   a_geometry.plotAtVelocityIndex("matvec", VEL::IntVect::Zero, matvec, 0.);
   a_geometry.plotAtVelocityIndex("diff", VEL::IntVect::Zero, diff, 0.);
#endif
      
   double relative_difference = MaxNorm(diff) / MaxNorm(pc_operator);

   if (procID()==0) {
      cout << "GKVlasov::testPC(): relative difference = " << relative_difference << " (using random input)" << endl;
   }
}


void GKVlasov::solvePCImEx( Preconditioner<ODEVector,AppCtxt>* const  a_pc,
                            KineticSpeciesPtrVect&                    a_kinetic_species_solution,
                            const KineticSpeciesPtrVect&              a_kinetic_species_rhs,
                            const CFG::FluidSpeciesPtrVect&           a_fluid_species_rhs,
                            int                                       a_idx ) const
{
   CH_TIME("GKVlasov::solvePCImEx");
   CH_assert(a_idx >= 0);

   KineticSpecies& kinetic_species_soln( static_cast<KineticSpecies&>(*(a_kinetic_species_solution[a_idx])) );
   {
     const std::string& species_name( kinetic_species_soln.name() );
     CH_assert(m_ti_op_type.at(species_name) == _implicit_op_);
   }
   GKVlasovAMG* this_pc = dynamic_cast<GKVlasovPreconditioner<ODEVector,AppCtxt>*>
                            (a_pc)->getAMGPC();
   KineticSpecies& kinetic_species_rhs( static_cast<KineticSpecies&>(*(a_kinetic_species_rhs[a_idx])) );
   this_pc->solve(kinetic_species_rhs.distributionFunction(), kinetic_species_soln.distributionFunction());
}


void
GKVlasov::constructStencils( const LevelData<FluxBox>&               a_physical_velocity,
                             const PhaseGeom&                        a_geometry,
                             LevelData<BaseFab<Vector<IntVect> > >&  a_couplings,
                             LevelData<BaseFab<Vector<Real> > >&     a_weights ) const
{
   CH_TIME("GKVlasov::constructStencils");
   
   const DisjointBoxLayout& grids( a_geometry.gridsFull() );

   // The passed in velocity is assumed to be in the physical frame.  Get the
   // normal component of the mapped velocity

   LevelData<FluxBox> normal_vel(grids, 1, IntVect::Zero);
   a_geometry.computeMetricTermProductAverage( normal_vel, a_physical_velocity, false );

   // The dimensioning of the following temporary arrays may need to be increased if other methods are added
   double method_stencil_weights[3];
   double method_stencil_weights_outflow[3];
   IntVect coupling_iv[3];
   double coupling_weight[3];

   int method_stencil_length;
   if ( m_precond_face_avg_type == "uw1" ) {
      method_stencil_length = 1;
      method_stencil_weights[0] = one;
      method_stencil_weights_outflow[0] = one;      // No changes for extrapolation at outflow boundaries
   }
   else if ( m_precond_face_avg_type == "uw3" ) {
      method_stencil_length = 3;

      method_stencil_weights[0] = third;
      method_stencil_weights[1] = 5.*sixth;
      method_stencil_weights[2] = -sixth;

      method_stencil_weights_outflow[0] = 0.;       // Modifications incorporating second-order
      method_stencil_weights_outflow[1] = 3.*half;  // extrapolation at outflow boundaries
      method_stencil_weights_outflow[2] = -half;
   }
   else {
      MayDay::Error("GKVlasov::computeStencils(): Unknown face_avg_type");
   }

   CH_assert(method_stencil_length <= 3);  // See comment above

   MultiBlockCoordSys* coord_sys_ptr = a_geometry.coordSysPtr();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = coord_sys_ptr->boundaries();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FluxBox& this_normal_vel( normal_vel[dit] );
      BaseFab<Vector<IntVect> >& this_couplings = a_couplings[dit];
      BaseFab<Vector<Real> >& this_weights = a_weights[dit];

      const PhaseBlockCoordSys& block_coord_sys = a_geometry.getBlockCoordSys(grids[dit]);
      double mapped_cell_volume = block_coord_sys.getMappedCellVolume();

      const Box& block_box = block_coord_sys.domain().domainBox();

      int block_number = coord_sys_ptr->whichBlock(grids[dit]);
      const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];

      for (int dir=0; dir<CFG_DIM+1; dir++) {
         Box face_box = surroundingNodes(grids[dit], dir);
         Box block_face_box = surroundingNodes(block_box, dir);

         bool lo_block_boundary = this_block_boundaries[dir].isDomainBoundary();
         bool hi_block_boundary = this_block_boundaries[dir + SpaceDim].isDomainBoundary();

         for (BoxIterator bit(face_box); bit.ok(); ++bit) {
            IntVect iv_face = bit();
            IntVect iv_lo_cell = iv_face; iv_lo_cell[dir]--;
            IntVect iv_hi_cell = iv_face;
            
            double v = this_normal_vel[dir](iv_face,0);

            if ( v >= 0. ) {
               coupling_iv[0] = iv_lo_cell;
            }
            else {
               coupling_iv[0] = iv_hi_cell;
            }

            if ( m_precond_face_avg_type == "uw3" ) {
               if ( v >= 0. ) {
                  coupling_iv[1] = coupling_iv[0];
                  coupling_iv[0][dir]++;
                  coupling_iv[2] = coupling_iv[1]; coupling_iv[2][dir]--;
               }
               else {
                  coupling_iv[1] = coupling_iv[0];
                  coupling_iv[0][dir]--;
                  coupling_iv[2] = coupling_iv[1]; coupling_iv[2][dir]++;
               }
            }

            double fac = v / mapped_cell_volume;

            if ( (lo_block_boundary && iv_face[dir] == block_face_box.smallEnd(dir) && v < 0.) ||
                 (hi_block_boundary && iv_face[dir] == block_face_box.bigEnd(dir) && v >= 0.) ) {
               // Using weights incorporating extrapolation at outflow boundaries
               for (int m=0; m<method_stencil_length; ++m) {
                  coupling_weight[m] = method_stencil_weights_outflow[m] * fac;
               }
            }
            else {
               for (int m=0; m<method_stencil_length; ++m) {
                  coupling_weight[m] = method_stencil_weights[m] * fac;
               }
            }
            
            // Update the couplings on the low side of the current face.
            for (int m=0; m<method_stencil_length; ++m) {
               int num_couplings = this_couplings(iv_lo_cell,dir).size();
               bool new_coupling = true;
               for (int n=0; n<num_couplings; ++n) {
                  if ( (this_couplings(iv_lo_cell,dir))[n] == coupling_iv[m] ) {
                     this_weights(iv_lo_cell,dir)[n] += coupling_weight[m];
                     new_coupling = false;
                     break;
                  }
               }
               if ( new_coupling ) {
                  this_couplings(iv_lo_cell,dir).push_back(coupling_iv[m]);
                  this_weights(iv_lo_cell,dir).push_back(coupling_weight[m]);
               }
            }               

            // Update the couplings on the high side of the current face.
            for (int m=0; m<method_stencil_length; ++m) {
               int num_couplings = this_couplings(iv_hi_cell,dir + CFG_DIM+1).size();
               bool new_coupling = true;
               for (int n=0; n<num_couplings; ++n) {
                  if ( (this_couplings(iv_hi_cell,dir + CFG_DIM+1))[n] == coupling_iv[m] ) {
                     this_weights(iv_hi_cell,dir + CFG_DIM+1)[n] -= coupling_weight[m];
                     new_coupling = false;
                     break;
                  }
               }
               if ( new_coupling ) {
                  this_couplings(iv_hi_cell,dir + CFG_DIM+1).push_back(coupling_iv[m]);
                  this_weights(iv_hi_cell,dir + CFG_DIM+1).push_back(-coupling_weight[m]);
               }
            }
         }
      }

   }
}

void
GKVlasov::computeIntegratedMomentFluxNormals(CFG::LevelData<CFG::FluxBox>&          a_flux_normal,
                                             const KineticSpecies&                  a_soln_species,
                                             const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                                             const CFG::EField&                     a_E_field,
                                             const int                              a_velocity_option,
                                             const string&                          a_moment_type,
                                             const Real                             a_time)
{
   /*
      Computes normal component of a moment flux integrated over a cell face
    */
   
   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& dbl( soln_dfn.getBoxes() );

   const PhaseGeom& geometry( a_soln_species.phaseSpaceGeometry() );
   const CFG::MagGeom& mag_geom = geometry.magGeom();
   const CFG::DisjointBoxLayout& config_grids = a_flux_normal.getBoxes();
   
   IntVect ghostVect = (geometry.secondOrder()) ? IntVect::Zero : IntVect::Unit;
   
   MomentOp& moment_op = MomentOp::instance();
   
   if ( geometry.divFreeVelocity() ) {
      bool fourth_order_Efield = !a_E_field.secondOrder();
      
      if ( !m_flux_normal.isDefined()) {
         m_flux_normal.define( dbl, 1, ghostVect );
      }
   
      computeIntegratedFluxNormals(m_flux_normal,
                                   a_soln_species,
                                   a_phi,
                                   a_E_field.getCellCenteredField(),
                                   a_E_field.getPhiNode(),
                                   fourth_order_Efield,
                                   a_velocity_option,
                                   a_time);

      if ( a_moment_type == "particle" ) {
        moment_op.compute( a_flux_normal, a_soln_species, m_flux_normal, DensityKernel<FluxBox>() );
      }
      else if ( a_moment_type == "energy" ) {
       moment_op.compute( a_flux_normal, a_soln_species, m_flux_normal, EnergyKernel<FluxBox>(a_phi) );
      }
      else {
       MayDay::Error("GKVlasov::computeIntegratedMomentFluxNormals:: unknwn moment option");
      }

      const VEL::VelCoordSys& vel_coords = geometry.velSpaceCoordSys();
      double vel_area = vel_coords.pointwiseJ(VEL::RealVect::Zero) * vel_coords.dx().product();

      for (CFG::DataIterator dit(config_grids); dit.ok(); ++dit) {
         a_flux_normal[dit] *= (1./vel_area);
      }
   }
   
   else {
      if ( !m_flux.isDefined()) {
         m_flux.define( dbl, SpaceDim, ghostVect );
      }
   
      computePhysicalFlux(m_flux,
                          a_soln_species,
                          a_phi,
                          a_E_field,
                          a_velocity_option,
                          a_time,
                          true);
      
      const CFG::IntVect& ghosts = a_flux_normal.ghostVect();
      
      CFG::LevelData<CFG::FluxBox> phys_flux(config_grids, SpaceDim, ghosts);
      if ( a_moment_type == "particle" ) {
        moment_op.compute( phys_flux, a_soln_species, m_flux, DensityKernel<FluxBox>() );
      }
      else if ( a_moment_type == "energy" ) {
        moment_op.compute( phys_flux, a_soln_species, m_flux, EnergyKernel<FluxBox>(a_phi) );
      }
      else {
        MayDay::Error("GKVlasov::computeIntegratedMomentFluxNormals:: unknwn moment option");
      }
      
      CFG::LevelData<CFG::FluxBox> phys_flux_cfg(config_grids, CFG_DIM, ghosts);
      for (CFG::DataIterator dit(config_grids); dit.ok(); ++dit) {
         for (int dir=0; dir<CFG_DIM; ++dir) {
            phys_flux_cfg[dit][dir].copy(phys_flux[dit][dir],0,0,CFG_DIM);
         }
      }

      // Compute NT * phys_flux * face_area
      bool fourth_order = (mag_geom.secondOrder()) ? false : true;
      
      // computeMetricTermProductAverage will only work correctly if
      // the assertions below are satisfied
      if (fourth_order) {
        CH_assert(a_flux_normal.ghostVect() == CFG::IntVect::Unit);
        CH_assert(!geometry.secondOrder());
      }
      mag_geom.computeMetricTermProductAverage(a_flux_normal, phys_flux_cfg, fourth_order);
   }
}

void
GKVlasov::testZeroDivergence(const LevelData<FluxBox>&  a_velocity,
                             const PhaseGeom&           a_geometry) const
{
   const DisjointBoxLayout& dbl = a_velocity.disjointBoxLayout();
   LevelData<FArrayBox> velocity_divergence(dbl, 1, IntVect::Zero);
   
   // Test old velocity calculation
   if (a_velocity.nComp() == SpaceDim) {
      a_geometry.mappedGridDivergence( velocity_divergence, m_velocity, false );
   }
   // Test new velocity (a.k.a. divfree) calculation
   else if (a_velocity.nComp() == 1) {
      a_geometry.mappedGridDivergenceFromIntegratedFluxNormals( velocity_divergence, a_velocity );
   }
   else {
      MayDay::Error("GKVlasov::testZeroDivergence: numer of velocity components must be either 1 or SpaceDim");
   }

   // Divide by mapped cell volume and negate to bring on the same footing as
   // the RHS flux divergence
   for (DataIterator dit( a_velocity.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys&
         block_coord_sys( a_geometry.getBlockCoordSys( dbl[dit] ) );
      double fac( -1.0 / block_coord_sys.getMappedCellVolume() );
      velocity_divergence[dit].mult( fac );
   }

   double veldiv_norm = MaxNorm(velocity_divergence);
   if (procID()==0) cout << "velocity divergence norm = " << veldiv_norm << endl;
}

void
GKVlasov::testEnergyConservation(const LevelData<FArrayBox>&            a_rhs,
                                 const KineticSpecies&                  a_soln_species,
                                 const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                                 const CFG::EField&                     a_E_field)
{
  
  /*
   This function compares the energy moment of the Vlasov RHS (\int d^3 {Q * div_phase(Flux_phase)}),
   with the configuration divergence of the transport power ( div_cfg( \int {d^3 Q * Flux_cfg} )
   */
  
  const PhaseGeom& geometry( a_soln_species.phaseSpaceGeometry() );
  const CFG::MagGeom& mag_geom = geometry.magGeom();
    
  CFG::IntVect ghostVect = (geometry.secondOrder()) ? CFG::IntVect::Zero : CFG::IntVect::Unit;
  const CFG::DisjointBoxLayout& grids_cfg = mag_geom.grids();
  
  // Compute the configuration divergence of the transport power
  CFG::LevelData<CFG::FluxBox> flux_normal(grids_cfg, 1, ghostVect);
  
  computeIntegratedMomentFluxNormals(flux_normal,
                                     a_soln_species,
                                     a_phi,
                                     a_E_field,
                                     PhaseGeom::FULL_VELOCITY,
                                     "energy",
                                     0.);
  
  CFG::LevelData<CFG::FArrayBox> div_transport_power(grids_cfg, 1, CFG::IntVect::Zero);
  CFG::LevelData<CFG::FArrayBox> cell_voll_cfg(grids_cfg, 1, CFG::IntVect::Zero);
  mag_geom.getCellVolumes(cell_voll_cfg);
  
  CFG::RealVect fakeDx = CFG::RealVect::Unit;
  for (CFG::DataIterator dit(grids_cfg); dit.ok(); ++dit) {
    simpleDivergence(div_transport_power[dit], flux_normal[dit], grids_cfg[dit], fakeDx);
    div_transport_power[dit].divide(cell_voll_cfg[dit]);
  }
  
  // Compute the energy moment of the Vlasov RHS
  CFG::LevelData<CFG::FArrayBox> rhs_energy_mom(grids_cfg, 1, CFG::IntVect::Zero);

  MomentOp& moment_op = MomentOp::instance();
  moment_op.compute( rhs_energy_mom, a_soln_species, a_rhs, EnergyKernel<FArrayBox>(a_phi) );
  mag_geom.divideJonValid(rhs_energy_mom);
  for (CFG::DataIterator dit(grids_cfg); dit.ok(); ++dit) {
    rhs_energy_mom[dit].mult(-1.0);
  }

  mag_geom.plotCellData("rhs_energy_mom", rhs_energy_mom, 0.);
  mag_geom.plotCellData("div_transport_power", div_transport_power, 0.);

  // Compare fs-integrals
  CFG::FluxSurface fs_shell(mag_geom);
  
  CFG::LevelData<CFG::FArrayBox> rhs_energy_mom_fs(grids_cfg, 1, CFG::IntVect::Zero);
  fs_shell.averageAndSpread(rhs_energy_mom, rhs_energy_mom_fs);

  CFG::LevelData<CFG::FArrayBox> div_transport_power_fs(grids_cfg, 1, CFG::IntVect::Zero);
  fs_shell.averageAndSpread(div_transport_power, div_transport_power_fs);
  
  mag_geom.plotCellData("rhs_energy_mom_fs", rhs_energy_mom_fs, 0.);
  mag_geom.plotCellData("div_transport_power_fs", div_transport_power_fs, 0.);
  
  // Obtain fs-integrated transport power
  CFG::LevelData<CFG::FArrayBox> face_areas(grids_cfg, CFG_DIM, CFG::IntVect::Zero);
  mag_geom.getPointwiseFaceAreas(face_areas);

  CFG::LevelData<CFG::FArrayBox> flux_normal_cc(grids_cfg, 1, CFG::IntVect::Zero);
  for (CFG::DataIterator dit(grids_cfg); dit.ok(); ++dit) {
    CFG::BoxIterator bit(flux_normal_cc[dit].box());
    for (bit.begin(); bit.ok(); ++bit) {
      CFG::IntVect iv = bit();
      CFG::IntVect ivR(iv);
      ivR[RADIAL_DIR] += 1;
      flux_normal_cc[dit](iv,0) = 0.5*(flux_normal[dit][RADIAL_DIR](iv,0) + flux_normal[dit][RADIAL_DIR](ivR,0));
    }
    flux_normal_cc[dit].divide(face_areas[dit],RADIAL_DIR,0,1);
  }
  
  CFG::FluxSurface fs_face(mag_geom, false);
  CFG::LevelData<CFG::FArrayBox> heat_flux_fs(grids_cfg, 1, CFG::IntVect::Zero);
  fs_face.averageAndSpread(flux_normal_cc, heat_flux_fs);
  mag_geom.plotCellData("heat_flux_fs", heat_flux_fs, 0.);
  
  // Obtain difference between FS-integrated radial heat fluxes
  CFG::LevelData<CFG::FluxBox> flux_normal_rad;
  flux_normal_rad.define(flux_normal);
  for (CFG::DataIterator dit(grids_cfg); dit.ok(); ++dit) {
    for (int dir=0; dir<CFG_DIM; ++dir) {
      if (dir != RADIAL_DIR) flux_normal_rad[dit][dir].setVal(0.);
    }
  }
  CFG::LevelData<CFG::FArrayBox> energy_src_loc(grids_cfg, 1, CFG::IntVect::Zero);
  for (CFG::DataIterator dit(grids_cfg); dit.ok(); ++dit) {
     simpleDivergence(energy_src_loc[dit], flux_normal_rad[dit], grids_cfg[dit], fakeDx);
     energy_src_loc[dit].divide(cell_voll_cfg[dit]);
  }
  CFG::LevelData<CFG::FArrayBox> energy_src_fs(grids_cfg, 1, CFG::IntVect::Zero);
  fs_shell.averageAndSpread(energy_src_loc, energy_src_fs);
  mag_geom.plotCellData("energy_src_fs", energy_src_fs, 0.);

  exit(1);
}


void
GKVlasov::getEFieldwFLR(  LevelData<FluxBox>&                   a_E_wFLR,
                          const KineticSpecies&                 a_species,
                          const CFG::LevelData<CFG::FArrayBox>& a_phi,
                          const CFG::EField&                    a_Efield )
{
  CH_assert(  a_species.isGyrokinetic() && m_include_FLR_effects );
  static const int ncomp_Efield = 3;

  int order (a_Efield.secondOrder() ? 2 : 4);
  const PhaseGeom& phase_geom( a_species.phaseSpaceGeometry() );
  a_species.gyroaveragedEField( a_E_wFLR, a_phi, order );

  if (m_FLR_switch_func != NULL) {


    if (!m_FLR_bc_factor.isDefined()) {
      const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
      CFG::LevelData<CFG::FluxBox> FLR_bc_factor_cfg;
      FLR_bc_factor_cfg.define( mag_geom.grids(), 
                                1, 
                                phase_geom.config_restrict(a_E_wFLR.ghostVect()) );
      m_FLR_switch_func->assign( FLR_bc_factor_cfg, mag_geom, 0.0 );

      if (m_plot_FLR_switch) {
        mag_geom.plotFaceData("FLR_switch_var", FLR_bc_factor_cfg, 0.0);
      }


      phase_geom.injectConfigurationToPhase(  FLR_bc_factor_cfg,
                                              m_FLR_bc_factor );
    }

    const LevelData<FluxBox>& E_field_noFLR = a_Efield.getInjectedField();

    for (DataIterator dit(a_E_wFLR.dataIterator()); dit.ok(); ++dit) {
      for (int dir = 0 ; dir < SpaceDim; dir++) {

        FArrayBox& efield_wFLR = a_E_wFLR[dit][dir];
        const FArrayBox& efield_noFLR = E_field_noFLR[dit][dir];
        const FArrayBox& switch_fac = m_FLR_bc_factor[dit][dir];

        const Box& box = efield_wFLR.box();
        const Box& flat_box = efield_noFLR.box();

        int vpar_idx = flat_box.smallEnd(VPARALLEL_DIR);
        int mu_idx = flat_box.smallEnd(MU_DIR);

        VEL::IntVect vel_idx(vpar_idx, mu_idx);

        for (BoxIterator bit(box); bit.ok(); ++bit) {
          IntVect iv(bit());

          CFG::IntVect iv_cfg = phase_geom.config_restrict(iv);
          IntVect iv_flat = phase_geom.tensorProduct(iv_cfg, vel_idx);

          Real alpha = switch_fac(iv_flat, 0);
          for (int n=0; n < ncomp_Efield; n++) {
            Real e_wflr = efield_wFLR(iv,n);
            Real e_noflr = efield_noFLR(iv_flat,n);
            efield_wFLR(iv,n) = alpha*e_wflr + (1.0-alpha)*e_noflr;
          }
        }
      }
    }

  }

  return;
}

#include "NamespaceFooter.H"
