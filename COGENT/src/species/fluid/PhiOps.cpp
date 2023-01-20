#include "PhiOps.H"
#include "NewGKPoissonBoltzmann.H"
#include "Directions.H"
#include "LogRectCoordSys.H"
#include "ConstFact.H"
#include "SpaceUtils.H.multidim"

#include "NamespaceHeader.H"


PhiOps::PhiOps( const string&   a_pp_prefix,
                const MagGeom&  a_geometry )
   : m_geometry(a_geometry),
     m_poisson(NULL),
     m_FLR_charge_density(false),
     m_defined(false)
{
}


void PhiOps::define( const double                      a_larmor,
                     const double                      a_debye,
                     const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                     BoltzmannElectron*                a_boltzmann_electron,
                     const bool                        a_fixed_efield,
                     const bool                        a_support_divfree_phase_vel,
                     const int                         a_potential_ghosts,
                     const int                         a_cur_step )
{
   if ( m_defined ) {
      MayDay::Error("PhiOps::define(): Efield has already been created");
   }

   m_boltzmann_electron = a_boltzmann_electron;
   m_fixed_efield       = a_fixed_efield;
   m_include_pol_dens_i = true;
   m_include_pol_dens_e = false;
   
   const DisjointBoxLayout& grids = m_geometry.gridsFull();

   if ( m_boltzmann_electron ) {

      LevelData<FArrayBox> ion_charge_density( grids, 1, IntVect::Zero);
      computeIonChargeDensity( ion_charge_density, a_kinetic_species );

      ParmParse pp( GKPoissonBoltzmann::pp_name );
      parseParameters(pp);
      if ( typeid(*(m_geometry.getCoordSys())) == typeid(SingleNullCoordSys) ) {
         m_poisson = new NewGKPoissonBoltzmann( pp, m_geometry, a_larmor, a_debye, ion_charge_density );
      }
      else {
         m_poisson = new GKPoissonBoltzmann( pp, m_geometry, a_larmor, a_debye, ion_charge_density );
      }
   }
   else {
      ParmParse pp( GKPoisson::pp_name );
      parseParameters(pp);
      m_poisson = new GKPoisson(  pp, 
                                  m_geometry, 
                                  a_kinetic_species, 
                                  a_larmor, 
                                  a_debye, 
                                  true );
   }

   m_defined = true;
}


void PhiOps::defineEField( EMFields& a_EM_fields ) const
{
   LevelData<FluxBox>& E_field_face = a_EM_fields.getEFieldFace();
   E_field_face.define(m_geometry.gridsFull(), 3, IntVect::Unit);

   LevelData<FArrayBox>& E_field_cell = a_EM_fields.getEFieldCell();
   E_field_cell.define(m_geometry.gridsFull(), 3, IntVect::Unit);

   a_EM_fields.setSecondOrder(m_poisson->secondOrder());
}


void PhiOps::definePotential( const LevelData<FArrayBox>&  a_phi,
                              const bool                   a_include_phi_node,
                              EMFields&                    a_EM_fields ) const
{
   LevelData<FArrayBox>& phi = a_EM_fields.getPhi();
   phi.define(m_geometry.gridsFull(), 1, a_phi.ghostVect());

   if ( a_include_phi_node ) {
      LevelData<FArrayBox>& phi_node = a_EM_fields.getPhiNode();
      phi_node.define(m_geometry.gridsFull(), 1, 2*IntVect::Unit);
   }
}


void PhiOps::updateEField( EMFields&                         a_EM_fields,
                           const bool                        a_update_potential,
                           const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                           const FluidSpeciesPtrVect&        a_fluid_species,
                           const PS::ScalarPtrVect&          a_scalars,
                           EllipticOpBC&                     a_bc,
                           const bool                        a_inject_Efield,
                           const bool                        a_initial_time )
{
   LevelData<FArrayBox>& E_field_cell = a_EM_fields.getEFieldCell();
   CH_assert(E_field_cell.ghostVect() == IntVect::Unit);
   LevelData<FluxBox>& E_field_face = a_EM_fields.getEFieldFace();
   CH_assert(E_field_face.ghostVect() == IntVect::Unit);
   LevelData<FArrayBox>& phi = a_EM_fields.getPhi();

   if ( m_defined ) {

      if ( !m_fixed_efield ) {

         const DisjointBoxLayout& grids = m_geometry.gridsFull();

         bool single_null = typeid(*(m_geometry.getCoordSys())) == typeid(SingleNullCoordSys);

         // Update the potential and field, if not fixed_efield
         if ( a_update_potential ) {

            // ion_mass_density has two components: first one is this quantity summed
            // over *all* species; second one is this quantity summed over non-FLR species
            LevelData<FArrayBox> ion_mass_density( grids, 2, IntVect::Zero );
            computePolarizationMassDensity( ion_mass_density, a_kinetic_species, a_fluid_species );
            m_poisson->setOperatorCoefficients( a_kinetic_species, ion_mass_density, a_bc, true );
            
            if (m_boltzmann_electron == NULL) {
            
               LevelData<FArrayBox> total_charge_density( grids, 1, IntVect::Zero );
               computeTotalChargeDensity( total_charge_density, a_kinetic_species, a_fluid_species );
               
               if (!m_fixed_krho2) {
                  if (m_poisson->isBoussinesq()) {
                     applyBoussinesqCorrections(total_charge_density,ion_mass_density);
                  }
                  for (DataIterator dit(grids); dit.ok(); ++dit) {
                     phi[dit].setVal(0.);
                  }
                  m_poisson->solveWithBCs( phi, total_charge_density );
               }
               else {
                  for (DataIterator dit(phi.dataIterator()); dit.ok(); ++dit) {
                     phi[dit].setVal(0.);
                     phi[dit].copy(total_charge_density[dit]);
                     phi[dit].divide(m_krho2_value);
                  }
               }
               fillInternalGhosts(phi);
            }
            
            else {
               
               // Boltzmann electron model
               
               LevelData<FArrayBox> ion_charge_density( grids, 1, IntVect::Zero );
               computeIonChargeDensity( ion_charge_density, a_kinetic_species );
               
               if ( single_null ) {

                  LevelData<FArrayBox> ion_parallel_current_density( grids, 1, IntVect::Zero );
                  computeIonParallelCurrentDensity( ion_parallel_current_density, a_kinetic_species );

                  ((NewGKPoissonBoltzmann*)m_poisson)
                     ->setDivertorBVs( ion_charge_density, ion_parallel_current_density, a_bc );
              
                  ((NewGKPoissonBoltzmann*)m_poisson)
                     ->computePotentialAndElectronDensity( phi,
                                                           *m_boltzmann_electron,
                                                           ion_charge_density,
                                                           a_kinetic_species,
                                                           a_bc,
                                                           a_initial_time );
               }
               else {
                  ((GKPoissonBoltzmann*)m_poisson)
                     ->computePotentialAndElectronDensity( phi,
                                                           *m_boltzmann_electron,
                                                           ion_charge_density,
                                                           a_kinetic_species,
                                                           a_bc,
                                                           a_initial_time );
               }
            }

            // Filter out all harmonics except for the primary mode
            if (m_apply_harm_filtering) {
               SpaceUtils::applyHarmonicFiltering(phi, m_harm_filtering_dir);
            }
         }
      }

      if ( !m_fixed_efield || a_initial_time ) {
         fillInternalGhosts(phi);
         m_poisson->computeField( phi, E_field_cell );
         m_poisson->computeField( phi, E_field_face );
         
         // Update nodal phi if supporting the calculation of a divergence-free phase velocity
         interpToNodes(phi, a_EM_fields.getPhiNode());
      }

   }
   else {
      // Arrival here means that the define() member was never called.  Returning zero.
      for (DataIterator dit( E_field_face.dataIterator() ); dit.ok(); ++dit) {
         E_field_cell[dit].setVal(0.);
         E_field_face[dit].setVal(0.);
      }
   }

   // Inject the E-field
   if ( a_inject_Efield ) {
      inject(a_EM_fields.getEField(), a_kinetic_species, E_field_face, E_field_cell);
   }
}


void PhiOps::inject( PS::LevelData<PS::FluxBox>&       a_injected_data,
                     const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                     const LevelData<FluxBox>&         a_data_face,
                     const LevelData<FArrayBox>&       a_data_cell )
{
  if (a_kinetic_species.size() > 0 ) {
      const PS::PhaseGeom& phaseGeom( a_kinetic_species[0]->phaseSpaceGeometry() );

#if 0
      phaseGeom.injectConfigurationToPhase(a_data_face,
                                              a_data_cell,
                                              a_injected_data);
      
#else
      //This code provides some performance optimization
      if (!phaseGeom.secondOrder()) {
         phaseGeom.injectConfigurationToPhase(a_data_face,
                                              a_data_cell,
                                              a_injected_data);
      }
      else {
         //Strip ghost since they are not needed for 2nd order
	 //Work with assertions in field calculations later 
         const DisjointBoxLayout& grids = m_geometry.gridsFull();
         LevelData<FArrayBox> cell_tmp(grids, 3, IntVect::Zero);
         LevelData<FluxBox> face_tmp(grids, 3, IntVect::Zero);
         for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
            cell_tmp[dit].copy(a_data_cell[dit]);
            face_tmp[dit].copy(a_data_face[dit]);
         }

         phaseGeom.injectConfigurationToPhase(face_tmp,
                                              cell_tmp,
                                              a_injected_data,
                                              false );

      }
#endif
   }
}


void PhiOps::computeIonChargeDensity( LevelData<FArrayBox>&               a_ion_charge_density,
                                      const PS::KineticSpeciesPtrVect&    a_kinetic_species_phys,
                                      const FluidSpeciesPtrVect&          a_fluid_species_comp) const
{
   CH_TIME("VorticityOp::computeIonChargeDensity");
   
   // Container for individual species charge density
   LevelData<FArrayBox> species_charge_density;
   species_charge_density.define(a_ion_charge_density);
   
   for (DataIterator dit(a_ion_charge_density.dataIterator()); dit.ok(); ++dit) {
      a_ion_charge_density[dit].setVal(0.);
   }
   
   // Accumulate contribution from kinetic species
   for (int species(0); species<a_kinetic_species_phys.size(); species++) {
      
      const PS::KineticSpecies& this_species( *(a_kinetic_species_phys[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      // Compute the charge density for this species including ghost cells
      this_species.chargeDensity( species_charge_density );
      
      for (DataIterator dit(a_ion_charge_density.dataIterator()); dit.ok(); ++dit) {
         a_ion_charge_density[dit] += species_charge_density[dit];
      }
   }
   
   // Accumulate constribution from fluid species
   for (int species(0); species<a_fluid_species_comp.size(); species++) {

      if (typeid(*(a_fluid_species_comp[species])) == typeid(FluidSpecies)) {
         const FluidSpecies& this_species( static_cast<FluidSpecies&>(*(a_fluid_species_comp[species])) );

         if ( this_species.charge() < 0.0 ) continue;
         
         this_species.chargeDensity( species_charge_density );
         m_geometry.divideJonValid(species_charge_density);
      
         DataIterator dit( a_ion_charge_density.dataIterator() );
         for (dit.begin(); dit.ok(); ++dit) {
            a_ion_charge_density[dit].plus( species_charge_density[dit] );
         }
      }
   }
}


void PhiOps::setCoreBC( const double   a_core_inner_bv,
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
      MayDay::Error("PhiOps::setCoreBC(): unknown geometry with consistent bcs");
   }
}


void PhiOps::computeQuasiNeutralElectronDensity( LevelData<FArrayBox>&        a_quasineutral_density,
                                                 LevelData<FArrayBox>&        a_potential,
                                                 EllipticOpBC&                a_bc, 
                                                 const LevelData<FArrayBox>&  a_ion_density) const
{
   a_ion_density.copyTo( a_quasineutral_density );

   m_poisson->setOperatorCoefficients( a_ion_density, a_bc, true );
   const DisjointBoxLayout& grids( a_ion_density.disjointBoxLayout() );
   LevelData<FArrayBox> polarization_density( grids, 1, a_quasineutral_density.ghostVect() );

   m_poisson->computeFluxDivergence( a_potential, polarization_density, false );
   for (DataIterator cdit( grids.dataIterator() ); cdit.ok(); ++cdit) {
      a_quasineutral_density[cdit] -= polarization_density[cdit];
   }
}


void PhiOps::computeIonChargeDensity( LevelData<FArrayBox>&             a_ion_charge_density,
                                      const PS::KineticSpeciesPtrVect&  a_species) const
{
   // Container for individual species charge density
   const DisjointBoxLayout& mag_grids = m_geometry.gridsFull();
   LevelData<FArrayBox> species_charge_density( mag_grids, 1, IntVect::Zero );

   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      a_ion_charge_density[dit].setVal(0.);
   }

   // Accumulate contribution from kinetic species; add fluid species contribution later
   for (int species(0); species<a_species.size(); species++) {

      const PS::KineticSpecies& this_species( *(a_species[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      // Compute the charge density for this species
      if (this_species.isGyrokinetic() && m_FLR_charge_density) {
        this_species.gyroaveragedChargeDensity( species_charge_density );
      } else {
        this_species.chargeDensity( species_charge_density );
      }
      
      DataIterator dit( a_ion_charge_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         a_ion_charge_density[dit].plus( species_charge_density[dit] );
      }
   }
}


void PhiOps::computePolarizationMassDensity( LevelData<FArrayBox>&             a_mass_density,
                                             const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                             const FluidSpeciesPtrVect&        a_fluid_species ) const
{
   /* If a_mass_density has 2 components,
    * GKPoisson with FLR effects needs both the total mass density over all species
    * and a total mass density for all non-gyrokinetic species; the 2nd component
    * will contain the latter. */
   int nc = a_mass_density.nComp();

   // Container for individual species charge density
   const DisjointBoxLayout& mag_grids = m_geometry.gridsFull();
   LevelData<FArrayBox> species_mass_density( mag_grids, 1, IntVect::Zero );

   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      a_mass_density[dit].setVal(0.);
   }
   
   // Accumulate contribution from kinetic species
   for (int species(0); species<a_kinetic_species.size(); species++) {
         
      const PS::KineticSpecies& this_species( *(a_kinetic_species[species]) );
      if (( !m_include_pol_dens_e && this_species.charge() < 0.0 ) ||
          ( !m_include_pol_dens_i && this_species.charge() > 0.0 )) continue;

      // Compute the charge density for this species
      this_species.massDensity( species_mass_density );
      
      DataIterator dit( a_mass_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         a_mass_density[dit].plus( species_mass_density[dit], 0, 0 );
         if ((!this_species.isGyrokinetic()) && (nc == 2)) {
          a_mass_density[dit].plus( species_mass_density[dit], 0, 1 );
         }
      }
   }

   // Accumulate contribution from fluid species
   for (int species(0); species<a_fluid_species.size(); species++) {
      if (typeid(*(a_fluid_species[species])) == typeid(FluidSpecies)) {
         
         const FluidSpecies& this_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
         if (( !m_include_pol_dens_e && this_species.charge() < 0.0 ) ||
             ( !m_include_pol_dens_i && this_species.charge() > 0.0 )) continue;
      
         // Compute the charge density for this species
         this_species.massDensity( species_mass_density );
      
         DataIterator dit( a_mass_density.dataIterator() );
         for (dit.begin(); dit.ok(); ++dit) {
            a_mass_density[dit].plus( species_mass_density[dit], 0, 0 );
            if (nc == 2) {
              a_mass_density[dit].plus( species_mass_density[dit], 0, 1 );
            }
         }
      }
   }
}


void PhiOps::computeTotalChargeDensity( LevelData<FArrayBox>&             a_charge_density,
                                        const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                        const FluidSpeciesPtrVect&        a_fluid_species ) const
{
   // Container for individual species charge density
   const DisjointBoxLayout& mag_grids = m_geometry.gridsFull();
   LevelData<FArrayBox> species_charge_density( mag_grids, 1, IntVect::Zero );

   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      a_charge_density[dit].setVal(0.);
   }

   // Accumulate contribution from kinetic species
   for (int species(0); species<a_kinetic_species.size(); species++) {
      const PS::KineticSpecies& this_species( *(a_kinetic_species[species]) );

      // Compute the charge density for this species
      if (this_species.isGyrokinetic() && m_FLR_charge_density) {
        this_species.gyroaveragedChargeDensity( species_charge_density );
      } else {
        this_species.chargeDensity( species_charge_density );
      }

      DataIterator dit( a_charge_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         a_charge_density[dit].plus( species_charge_density[dit] );
      }
   }

   // Accumulate constribution from fluid species
   for (int species(0); species<a_fluid_species.size(); species++) {

      if (typeid(*(a_fluid_species[species])) == typeid(FluidSpecies)) {
         const FluidSpecies& this_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
         this_species.chargeDensity( species_charge_density );
      
         DataIterator dit( a_charge_density.dataIterator() );
         for (dit.begin(); dit.ok(); ++dit) {
            a_charge_density[dit].plus( species_charge_density[dit] );
         }
      }
   }
}


void PhiOps::computeIonParallelCurrentDensity( LevelData<FArrayBox>&             a_ion_current_density,
                                               const PS::KineticSpeciesPtrVect&  a_species ) const
{
   // Container for individual species charge density
   const DisjointBoxLayout& mag_grids = m_geometry.gridsFull();
   LevelData<FArrayBox> species_current_density( mag_grids, 1, IntVect::Zero );

   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      a_ion_current_density[dit].setVal(0.);
   }

   for (int species(0); species<a_species.size(); species++) {

      const PS::KineticSpecies& this_species( *(a_species[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      // Compute the charge density for this species
      this_species.parallelParticleFlux( species_current_density );
      
      DataIterator dit( a_ion_current_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         species_current_density[dit].mult( this_species.charge() );
         a_ion_current_density[dit].plus( species_current_density[dit] );
      }
   }
}

void
PhiOps::applyBoussinesqCorrections( LevelData<FArrayBox>& a_total_charge_density,
                                    const LevelData<FArrayBox>& a_ion_mass_density ) const
{
   
   const LevelData<FArrayBox>& BfieldMag = m_geometry.getCCBFieldMag();
   
   const DisjointBoxLayout& grids = a_ion_mass_density.getBoxes();
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_total_charge_density[dit].divide(a_ion_mass_density[dit]);
      a_total_charge_density[dit].mult(BfieldMag[dit]);
      a_total_charge_density[dit].mult(BfieldMag[dit]);
   }
}


void
PhiOps::parseParameters( ParmParse& a_pp)
{
  a_pp.query( "include_ion_polarization_density", m_include_pol_dens_i);
  a_pp.query( "include_electron_polarization_density", m_include_pol_dens_e);
  a_pp.query( "FLR_charge_density", m_FLR_charge_density);
   
  if (a_pp.contains("harmonic_filtering")) {
     a_pp.get("harmonic_filtering", m_apply_harm_filtering);
  }
  else {
     m_apply_harm_filtering = false;
  }
  
  if (m_apply_harm_filtering) {
     a_pp.get("harmonic_filtering_dir", m_harm_filtering_dir);
  }

  if (a_pp.contains("fixed_krho2")) {
    a_pp.get("fixed_krho2", m_fixed_krho2);
  }
  else {
    m_fixed_krho2 = false;
  }
   
  if (m_fixed_krho2) {
    a_pp.get("krho2_value", m_krho2_value);
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
  else {
    m_consistent_lower_bc_only = false;
  }
}


#include "NamespaceFooter.H"
