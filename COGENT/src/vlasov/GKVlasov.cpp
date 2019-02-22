#include "GKVlasov.H"
#include "Directions.H"

#include "FourthOrderUtil.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"

#include "mappedLimiterF_F.H"
#include "altFaceAverages.H"
#include "mappedAdvectionFlux.H"
#include "KineticFunctionLibrary.H"

#undef TEST_ZERO_DIVERGENCE

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



GKVlasov::GKVlasov( ParmParse& a_pp,
                    const Real a_larmor_number )
  : m_larmor_number(a_larmor_number),
    m_face_avg_type(INVALID),
    m_dt_dim_factor(1.0)
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
   }
   else {
      if ( procID()==0 ) MayDay::Warning("Flux type unspecified; using default (BWENO)");
      m_face_avg_type = BWENO;
   }
   // In theory, advection is in (PDIM-1) dimensions, so we could relax this
   // a little, but for now, let's be conservative.
   m_dt_dim_factor = (m_face_avg_type>PPM) ? sqrt(PDIM) : 1.0;
}


void GKVlasov::accumulateRHS( GKRHSData&                            a_rhs,
                              const KineticSpeciesPtrVect&          a_kinetic_phys,
                              const CFG::LevelData<CFG::FArrayBox>& a_phi,
                              const CFG::EField&                    a_E_field,
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
      const PhaseGeom& phase_geom = rhs_species.phaseSpaceGeometry();
      if ( phase_geom.divFreeVelocity() ) {
         bool fourth_order_Efield = !a_E_field.secondOrder();
         evalRHS( rhs_species, 
                  soln_species, 
                  a_phi, 
                  a_E_field.getCellCenteredField(), 
                  a_E_field.getPhiNode(),
                  fourth_order_Efield, 
                  PhaseGeom::FULL_VELOCITY, 
                  a_time );
      }
      else {
         evalRHS( rhs_species, soln_species, a_phi, a_E_field, a_time );
      }
   }
}


void
GKVlasov::evalRHS( KineticSpecies&                        a_rhs_species,
                   const KineticSpecies&                  a_soln_species,
                   const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                   const CFG::EField&                     a_E_field,
                   const Real                             a_time )
{
   /*
     Evaluates the (negated) phase space divergence:
        rhs = - divergence_R ( R_dot soln ) - divergence_v ( v_dot soln )
   */
   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& dbl( soln_dfn.getBoxes() );
   const PhaseGeom& geometry( a_rhs_species.phaseSpaceGeometry() );

   LevelData<FluxBox> velocity( dbl, SpaceDim, IntVect::Unit );

   if (a_soln_species.isGyrokinetic()) {

     LevelData<FluxBox> gyroaveraged_E_field;
     int order (a_E_field.secondOrder() ? 2 : 4);
     a_soln_species.gyroaveragedEField( gyroaveraged_E_field, 
                                        a_phi, 
                                        order );

     a_soln_species.computeVelocity( velocity, gyroaveraged_E_field);

   } else {

     LevelData<FluxBox> injected_E_field;
     geometry.injectConfigurationToPhase( a_E_field.getFaceCenteredField(),
                                          a_E_field.getCellCenteredField(),
                                          injected_E_field );
  
     a_soln_species.computeVelocity( velocity, injected_E_field );

   }

   LevelData<FluxBox> flux( dbl, SpaceDim, IntVect::Unit );
   computeFlux( soln_dfn, velocity, flux, geometry );

   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   const bool OMIT_NT(false);
   geometry.mappedGridDivergence( rhs_dfn, flux, OMIT_NT );

#ifdef TEST_ZERO_DIVERGENCE
   LevelData<FArrayBox> velocity_divergence(dbl, 1, IntVect::Zero);
   geometry.mappedGridDivergence( velocity_divergence, velocity, OMIT_NT );
#endif

   // Divide by cell volume and negate
   for (DataIterator dit( rhs_dfn.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys&
         block_coord_sys( geometry.getBlockCoordSys( dbl[dit] ) );
      double fac( -1.0 / block_coord_sys.getMappedCellVolume() );
      rhs_dfn[dit].mult( fac );
#ifdef TEST_ZERO_DIVERGENCE
      velocity_divergence[dit].mult( fac );
#endif
   }

#ifdef TEST_ZERO_DIVERGENCE
   double veldiv_norm = MaxNorm(velocity_divergence);
   double rhs_norm = MaxNorm(rhs_dfn);
   if (procID()==0) cout << "velocity divergence norm = " << veldiv_norm << ", rhs norm = " << rhs_norm << endl;
#endif
}



void
GKVlasov::evalRHS( KineticSpecies&                        a_rhs_species,
                   const KineticSpecies&                  a_soln_species,
                   const CFG::LevelData<CFG::FArrayBox>&  a_phi,
                   const CFG::LevelData<CFG::FArrayBox>&  a_Efield_cell,
                   const CFG::LevelData<CFG::FArrayBox>&  a_phi_node,
                   const bool                             a_fourth_order_Efield,
                   const int                              a_velocity_option,
                   const Real                             a_time )
{
   /*
     Evaluates the (negated) phase space divergence:
        rhs = - divergence_R ( R_dot soln ) - divergence_v ( v_dot soln )
   */
   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& dbl( soln_dfn.getBoxes() );

   const PhaseGeom& geometry( a_rhs_species.phaseSpaceGeometry() );
   const CFG::MagGeom& mag_geom = geometry.magGeom();
   
   LevelData<FArrayBox> dfn_no_bstar(dbl, soln_dfn.nComp(), soln_dfn.ghostVect());
   for (DataIterator dit(dbl); dit.ok(); ++dit) {
      dfn_no_bstar[dit].copy(soln_dfn[dit]);
   }
   geometry.divideBStarParallel(dfn_no_bstar);

   LevelData<FluxBox> flux_normal( dbl, 1, IntVect::Unit );
   LevelData<FluxBox> velocity_normal( dbl, 1, IntVect::Unit );
   
   // Compute normal flux
   if (!m_subtract_maxwellian) {
      a_soln_species.computeMappedVelocityNormals( velocity_normal, a_Efield_cell, a_phi_node, a_fourth_order_Efield, a_velocity_option );
      computeFluxNormal( dfn_no_bstar, velocity_normal, flux_normal, geometry );
   }
   else {

      LevelData<FArrayBox> maxwellian_dfn(dbl, soln_dfn.nComp(), soln_dfn.ghostVect());
      LevelData<FArrayBox> delta_dfn(dbl, soln_dfn.nComp(), soln_dfn.ghostVect());
      computeDeltaF(a_soln_species, dfn_no_bstar, delta_dfn, maxwellian_dfn);
      
      LevelData<FluxBox> flux_tmp( dbl, 1, IntVect::Unit );
      a_soln_species.computeMappedVelocityNormals( velocity_normal, a_Efield_cell, a_phi_node, a_fourth_order_Efield, a_velocity_option );
      computeFluxNormal( delta_dfn, velocity_normal, flux_tmp, geometry );
 
      a_soln_species.computeMappedVelocityNormals( velocity_normal, a_Efield_cell, a_phi_node, a_fourth_order_Efield,
                                                   PhaseGeom::NO_ZERO_ORDER_PARALLEL_VELOCITY );

      if (!m_update_maxwellian) {
	computeFluxNormal( m_F0, velocity_normal, flux_normal, geometry);
      }
      else {
	computeFluxNormal( maxwellian_dfn, velocity_normal, flux_normal, geometry);
      }

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
         flux_normal[dit] += flux_tmp[dit];
      }

   }

   // Enforce conservation
   if (!mag_geom.extrablockExchange()) {
      geometry.averageAtBlockBoundaries(flux_normal);
   }
   
   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   geometry.mappedGridDivergenceFromIntegratedFluxNormals( rhs_dfn, flux_normal );

#ifdef TEST_ZERO_DIVERGENCE
   LevelData<FArrayBox> velocity_divergence(dbl, 1, IntVect::Zero);
   geometry.mappedGridDivergenceFromIntegratedFluxNormals( velocity_divergence, velocity_normal );
#endif

   // Divide by cell volume and negate
   for (DataIterator dit( rhs_dfn.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys&
         block_coord_sys( geometry.getBlockCoordSys( dbl[dit] ) );
      double fac( -1.0 / block_coord_sys.getMappedCellVolume() );
      rhs_dfn[dit].mult( fac );
#ifdef TEST_ZERO_DIVERGENCE
      velocity_divergence[dit].mult( fac );
#endif
   }

#ifdef TEST_ZERO_DIVERGENCE
   double veldiv_norm = MaxNorm(velocity_divergence);
   double rhs_norm = MaxNorm(rhs_dfn);
   if (procID()==0) cout << "velocity divergence norm = " << veldiv_norm << ", rhs norm = " << rhs_norm << endl;
#endif

}

void
GKVlasov::evalRHS( KineticSpecies&                        a_rhs_species,
                   const KineticSpecies&                  a_soln_species,
                   const CFG::LevelData<CFG::FArrayBox>&  a_Efield_cell,
                   const CFG::LevelData<CFG::FArrayBox>&  a_phi_node,
                   const bool                             a_fourth_order_Efield,
                   const int                              a_velocity_option,
                   const Real                             a_time )
{
   /*
     Evaluates the (negated) phase space divergence:
        rhs = - divergence_R ( R_dot soln ) - divergence_v ( v_dot soln )
   */
   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& dbl( soln_dfn.getBoxes() );

   const PhaseGeom& geometry( a_rhs_species.phaseSpaceGeometry() );
   const CFG::MagGeom& mag_geom = geometry.magGeom();
   
   LevelData<FArrayBox> dfn_no_bstar(dbl, soln_dfn.nComp(), soln_dfn.ghostVect());
   for (DataIterator dit(dbl); dit.ok(); ++dit) {
      dfn_no_bstar[dit].copy(soln_dfn[dit]);
   }
   geometry.divideBStarParallel(dfn_no_bstar);

   LevelData<FluxBox> flux_normal( dbl, 1, IntVect::Unit );
   LevelData<FluxBox> velocity_normal( dbl, 1, IntVect::Unit );
   
   // Compute normal flux
   if (!m_subtract_maxwellian) {
      a_soln_species.computeMappedVelocityNormals( velocity_normal, a_Efield_cell, a_phi_node, a_fourth_order_Efield, a_velocity_option );
      computeFluxNormal( dfn_no_bstar, velocity_normal, flux_normal, geometry );
   }
   else {

      LevelData<FArrayBox> maxwellian_dfn(dbl, soln_dfn.nComp(), soln_dfn.ghostVect());
      LevelData<FArrayBox> delta_dfn(dbl, soln_dfn.nComp(), soln_dfn.ghostVect());
      computeDeltaF(a_soln_species, dfn_no_bstar, delta_dfn, maxwellian_dfn);
      
      LevelData<FluxBox> flux_tmp( dbl, 1, IntVect::Unit );
      a_soln_species.computeMappedVelocityNormals( velocity_normal, a_Efield_cell, a_phi_node, a_fourth_order_Efield, a_velocity_option );
      computeFluxNormal( delta_dfn, velocity_normal, flux_tmp, geometry );
 
      a_soln_species.computeMappedVelocityNormals( velocity_normal, a_Efield_cell, a_phi_node, a_fourth_order_Efield,
                                                   PhaseGeom::NO_ZERO_ORDER_PARALLEL_VELOCITY );

      if (!m_update_maxwellian) {
	computeFluxNormal( m_F0, velocity_normal, flux_normal, geometry);
      }
      else {
	computeFluxNormal( maxwellian_dfn, velocity_normal, flux_normal, geometry);
      }

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
         flux_normal[dit] += flux_tmp[dit];
      }

   }

   // Enforce conservation
   if (!mag_geom.extrablockExchange()) {
      geometry.averageAtBlockBoundaries(flux_normal);
   }
   
   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   geometry.mappedGridDivergenceFromIntegratedFluxNormals( rhs_dfn, flux_normal );

#ifdef TEST_ZERO_DIVERGENCE
   LevelData<FArrayBox> velocity_divergence(dbl, 1, IntVect::Zero);
   geometry.mappedGridDivergenceFromIntegratedFluxNormals( velocity_divergence, velocity_normal );
#endif

   // Divide by cell volume and negate
   for (DataIterator dit( rhs_dfn.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys&
         block_coord_sys( geometry.getBlockCoordSys( dbl[dit] ) );
      double fac( -1.0 / block_coord_sys.getMappedCellVolume() );
      rhs_dfn[dit].mult( fac );
#ifdef TEST_ZERO_DIVERGENCE
      velocity_divergence[dit].mult( fac );
#endif
   }

#ifdef TEST_ZERO_DIVERGENCE
   double veldiv_norm = MaxNorm(velocity_divergence);
   double rhs_norm = MaxNorm(rhs_dfn);
   if (procID()==0) cout << "velocity divergence norm = " << veldiv_norm << ", rhs norm = " << rhs_norm << endl;
#endif

}


void
GKVlasov::computeFlux( const LevelData<FArrayBox>& a_dist_fn,
                       const LevelData<FluxBox>&   a_velocity,
                       LevelData<FluxBox>&         a_flux,
                       const PhaseGeom&            a_phase_geom )
{

   /*
      Compute the phase space flux given the input phase space advection
      velocity.   This is where the hyperbolic stuff connects.
   */

   // Construct appropriately accurate face-averages of phi and advVel
   LevelData<FluxBox> faceDist(a_dist_fn.getBoxes(), a_dist_fn.nComp(), a_dist_fn.ghostVect() );

   // If we're limiting the face-centered values, do it here
   bool mapped_velocity = false;
   if (m_face_avg_type==PPM) {
      computeFaceAverages( faceDist, a_dist_fn, a_phase_geom.secondOrder() );
      faceDist.exchange();
      applyMappedLimiter( faceDist, a_dist_fn, a_velocity, a_phase_geom, mapped_velocity );
   }
   else if (m_face_avg_type==UW1) {
      uw1FaceAverages( faceDist, a_dist_fn, a_velocity, a_phase_geom, mapped_velocity );
   }
   else if (m_face_avg_type==UW3) {
      uw3FaceAverages( faceDist, a_dist_fn, a_velocity, a_phase_geom, mapped_velocity );
   }
   else if (m_face_avg_type==UW5) {
      uw5FaceAverages( faceDist, a_dist_fn, a_velocity, a_phase_geom, mapped_velocity );
   }
   else if (m_face_avg_type==WENO5) {
      weno5FaceAverages( faceDist, a_dist_fn, a_velocity, a_phase_geom, mapped_velocity );
   }
   else if (m_face_avg_type==BWENO) {
      bwenoFaceAverages( faceDist, a_dist_fn, a_velocity, a_phase_geom, mapped_velocity );
   }

   if ( a_phase_geom.secondOrder() ) {

      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxes( a_flux, faceDist, a_velocity, false );
   }
   else {

      const DisjointBoxLayout& grids = a_flux.disjointBoxLayout();
      LevelData<FluxBox> fourth_order_flux(grids, SpaceDim, IntVect::Zero);
      
      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxes( fourth_order_flux, faceDist, a_velocity, true );

      // Compute the second-order flux in valid plus ghost cell faces,
      // then overwrite with the fourth-order flux on the valid faces.
      CH_assert(a_flux.ghostVect() == IntVect::Unit);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_flux[dit].copy(a_velocity[dit]);

         Box box = grow(grids[dit],1);
         for (int dir=0; dir<SpaceDim; ++dir) {
            a_flux[dit].mult(faceDist[dit],box,0,dir,1);
         }
         
         a_flux[dit].copy(fourth_order_flux[dit],grids[dit]);
      }
   }

   a_flux.exchange();
}


void
GKVlasov::computeFluxNormal( const LevelData<FArrayBox>& a_dist_fn,
                             const LevelData<FluxBox>&   a_velocity_normal,
                             LevelData<FluxBox>&         a_flux_normal,
                             const PhaseGeom&            a_phase_geom )
{

   /*
      Compute the phase space flux normal given the input phase space advection
      velocity normal.   This is where the hyperbolic stuff connects.
   */

   // Construct appropriately accurate face-averages of phi and advVel
   LevelData<FluxBox> faceDist(a_dist_fn.getBoxes(), a_dist_fn.nComp(), a_dist_fn.ghostVect() );

   // If we're limiting the face-centered values, do it here
   bool mapped_velocity = true;
   if (m_face_avg_type==PPM) {
      computeFaceAverages( faceDist, a_dist_fn, a_phase_geom.secondOrder() );
      faceDist.exchange();
      applyMappedLimiter( faceDist, a_dist_fn, a_velocity_normal, a_phase_geom, mapped_velocity );
   }
   else if (m_face_avg_type==UW1) {
      uw1FaceAverages( faceDist, a_dist_fn, a_velocity_normal, a_phase_geom, mapped_velocity );
   }
   else if (m_face_avg_type==UW3) {
      uw3FaceAverages( faceDist, a_dist_fn, a_velocity_normal, a_phase_geom, mapped_velocity );
   }
   else if (m_face_avg_type==UW5) {
      uw5FaceAverages( faceDist, a_dist_fn, a_velocity_normal, a_phase_geom, mapped_velocity );
   }
   else if (m_face_avg_type==WENO5) {
      weno5FaceAverages( faceDist, a_dist_fn, a_velocity_normal, a_phase_geom, mapped_velocity );
   }
   else if (m_face_avg_type==BWENO) {
      bwenoFaceAverages( faceDist, a_dist_fn, a_velocity_normal, a_phase_geom, mapped_velocity );
   }

   if ( a_phase_geom.secondOrder() ) {

      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxesNormal( a_flux_normal, faceDist, a_velocity_normal, false );
   }
   else {

      const DisjointBoxLayout& grids = a_flux_normal.disjointBoxLayout();
      LevelData<FluxBox> fourth_order_flux_normal(grids, 1, IntVect::Zero);
      
      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxesNormal( fourth_order_flux_normal, faceDist, a_velocity_normal, true );

      // Compute the second-order flux in valid plus ghost cell faces,
      // then overwrite with the fourth-order flux on the valid faces.
      CH_assert(a_flux_normal.ghostVect() == IntVect::Unit);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_flux_normal[dit].copy(a_velocity_normal[dit]);

         Box box = grow(grids[dit],1);
         a_flux_normal[dit].mult(faceDist[dit],box,0,0,1);
         
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
   //Presently works only to second order
   const PhaseGeom& geometry = a_soln_species.phaseSpaceGeometry();

   static bool first_call = true;
   if (!m_update_maxwellian) {

      if (first_call) {
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
         a_delta_dfn[dit].copy(a_dfn_no_bstar[dit]);
         a_delta_dfn[dit].minus(m_F0[dit]);
      }
   }
   
   else {

      a_soln_species.computeFSavgMaxwellian(a_maxwellian_dfn);

      for (DataIterator dit(a_delta_dfn.dataIterator()); dit.ok(); ++dit) {
         a_delta_dfn[dit].copy(a_dfn_no_bstar[dit]);
         a_delta_dfn[dit].minus(a_maxwellian_dfn[dit]);
      }
   }
   
   first_call = false;
   
}


Real
GKVlasov::computeDt( const CFG::EField&            a_E_field,
                     const KineticSpeciesPtrVect&  a_species_vect )
{
   Real dt(BASEFAB_REAL_SETVAL);

   if ( a_E_field.supportsDivFreePhaseVel() ) {

      const CFG::LevelData<CFG::FArrayBox>& Efield_cell = a_E_field.getCellCenteredField();
      const CFG::LevelData<CFG::FArrayBox>& phi_node = a_E_field.getPhiNode();

      for (int s(0); s<a_species_vect.size(); s++) {

         const KineticSpecies& species( *(a_species_vect[s]) );
         const LevelData<FArrayBox>& dfn( species.distributionFunction() );

         LevelData<FluxBox> velocity_normals( dfn.getBoxes(), 1, IntVect::Unit );
         species.computeMappedVelocityNormals( velocity_normals, Efield_cell, phi_node, false, PhaseGeom::FULL_VELOCITY );

         const Real UNIT_CFL(1.0);
         const PhaseGeom& geometry( species.phaseSpaceGeometry() );
         Real speciesDt( computeMappedDtSpeciesFromNormals( velocity_normals, geometry, UNIT_CFL ) );
         CH_assert(speciesDt >= 0);

         dt = Min( dt, speciesDt );
      }
   }
   else {

      LevelData<FluxBox> injected_E_field;

      for (int s(0); s<a_species_vect.size(); s++) {

         const KineticSpecies& species( *(a_species_vect[s]) );
         const LevelData<FArrayBox>& dfn( species.distributionFunction() );
         const PhaseGeom& geometry( species.phaseSpaceGeometry() );

         if ( !injected_E_field.isDefined() ) {
            geometry.injectConfigurationToPhase( a_E_field.getFaceCenteredField(),
                                                 a_E_field.getCellCenteredField(),
                                                 injected_E_field );
         }

         LevelData<FluxBox> velocity( dfn.getBoxes(), SpaceDim, IntVect::Unit );
         species.computeMappedVelocity( velocity, injected_E_field );

         const Real UNIT_CFL(1.0);
         Real speciesDt( computeMappedDtSpecies( velocity, geometry, UNIT_CFL ) );
         CH_assert(speciesDt >= 0);

         dt = Min( dt, speciesDt );
      }
   }

   return dt;
}


Real
GKVlasov::computeTimeScale( const CFG::EField&           a_E_field,
                            const KineticSpeciesPtrVect& a_species_vect )
{
   Real dt(BASEFAB_REAL_SETVAL);

   if ( a_E_field.supportsDivFreePhaseVel() ) {

      const CFG::LevelData<CFG::FArrayBox>& Efield_cell = a_E_field.getCellCenteredField();
      const CFG::LevelData<CFG::FArrayBox>& phi_node = a_E_field.getPhiNode();

      for (int s(0); s<a_species_vect.size(); s++) {

         const KineticSpecies& species( *(a_species_vect[s]) );
         const LevelData<FArrayBox>& dfn( species.distributionFunction() );

         LevelData<FluxBox> velocity_normals( dfn.getBoxes(), 1, IntVect::Unit );
         species.computeMappedVelocityNormals( velocity_normals, Efield_cell, phi_node, false,
                                               PhaseGeom::FULL_VELOCITY);

         //      const Real UNIT_CFL(1.0);
         const PhaseGeom& geometry( species.phaseSpaceGeometry() );
         Real speciesDt( computeMappedTimeScaleSpeciesFromNormals( velocity_normals, geometry) );
         CH_assert(speciesDt >= 0);

         dt = Min( dt, speciesDt );
      }
   }
   else {

      LevelData<FluxBox> injected_E_field;

      for (int s(0); s<a_species_vect.size(); s++) {

         const KineticSpecies& species( *(a_species_vect[s]) );
         const LevelData<FArrayBox>& dfn( species.distributionFunction() );
         const PhaseGeom& geometry( species.phaseSpaceGeometry() );

         if ( !injected_E_field.isDefined() ) {
            geometry.injectConfigurationToPhase( a_E_field.getFaceCenteredField(),
                                                 a_E_field.getCellCenteredField(),
                                                 injected_E_field );
         }

         LevelData<FluxBox> velocity( dfn.getBoxes(), SpaceDim, IntVect::Unit );
         species.computeMappedVelocity( velocity, injected_E_field );
         
         //      const Real UNIT_CFL(1.0);
         Real speciesDt( computeMappedTimeScaleSpecies( velocity, geometry) );
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

#ifdef PLOT_STABLEDT
   LevelData<FArrayBox> stableDt(grids, 1);
#endif

   // now average back to cell-centers and divide by cell volumes.
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
      cout << "The (";
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
      cout << ") direction makes the largest contribution to the stable dt reciprocal at that point" << endl;
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



#include "NamespaceFooter.H"
