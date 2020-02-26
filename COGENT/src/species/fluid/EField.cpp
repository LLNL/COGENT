#include "EField.H"
#include "NewGKPoissonBoltzmann.H"
#include "Directions.H"
#include "LogRectCoordSys.H"
#include "ConstFact.H"

#include "NamespaceHeader.H"


void EField::define( const double                      a_larmor,
                     const double                      a_debye,
                     const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                     BoltzmannElectron*                a_boltzmann_electron,
                     const bool                        a_fixed_efield,
                     const bool                        a_support_divfree_phase_vel,
                     const int                         a_cur_step )
{
   if ( m_defined ) {
      MayDay::Error("EField::define(): GKPoisson has already been created");
   }

   m_boltzmann_electron = a_boltzmann_electron;
   m_fixed_efield       = a_fixed_efield;
   m_include_pol_dens_i = true;
   m_include_pol_dens_e = false;
   
   const MagGeom& mag_geom = configurationSpaceGeometry();
   const DisjointBoxLayout& grids = mag_geom.gridsFull();

   if ( m_boltzmann_electron ) {

      LevelData<FArrayBox> ion_charge_density( grids, 1, IntVect::Zero);
      computeIonChargeDensity( ion_charge_density, a_kinetic_species );

      ParmParse pp( GKPoissonBoltzmann::pp_name );
      parseParameters(pp);
      if ( typeid(*(mag_geom.getCoordSys())) == typeid(SingleNullCoordSys) ) {
         m_poisson = new NewGKPoissonBoltzmann( pp, mag_geom, a_larmor, a_debye, ion_charge_density );
      }
      else {
         m_poisson = new GKPoissonBoltzmann( pp, mag_geom, a_larmor, a_debye, ion_charge_density );
      }
   }
   else {
      ParmParse pp( GKPoisson::pp_name );
      parseParameters(pp);
      m_poisson = new GKPoisson( pp, mag_geom, a_larmor, a_debye );
      dynamic_cast<GKPoisson*>(m_poisson)->setPhaseSpaceObjs(a_kinetic_species);
   }

   if ( a_support_divfree_phase_vel ) {
      m_phi_node.define(grids, 1, 2*IntVect::Unit);
   }

   m_defined = true;
}


void EField::computeEField( const PS::GKState&                a_state,
                            const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                            const FluidSpeciesPtrVect&        a_fluid_species,
                            const PS::ScalarPtrVect&          a_scalars,
                            LevelData<FArrayBox>&             a_phi,
                            EllipticOpBC&                     a_bc,
                            const bool                        a_update_potential,
                            const bool                        a_inject_field,
                            const bool                        a_initial_time )
{
   LevelData<FArrayBox>& E_field_cell = cell_var(0);
   CH_assert(E_field_cell.ghostVect() == IntVect::Unit);
   LevelData<FluxBox>& E_field_face = face_var(0);
   CH_assert(E_field_face.ghostVect() == IntVect::Unit);

   if ( m_defined ) {

      if ( !m_fixed_efield ) {

         const MagGeom& mag_geom = configurationSpaceGeometry();
         const DisjointBoxLayout& grids = mag_geom.gridsFull();

         bool single_null = typeid(*(mag_geom.getCoordSys())) == typeid(SingleNullCoordSys);

         // Update the potential and field, if not fixed_efield
         if ( a_update_potential ) {

            LevelData<FArrayBox> ion_mass_density( grids, 1, IntVect::Zero );
            computePolarizationMassDensity( ion_mass_density, a_kinetic_species, a_fluid_species );
            m_poisson->setOperatorCoefficients( a_kinetic_species, ion_mass_density, a_bc, true );
            
            if (m_boltzmann_electron == NULL) {
            
               LevelData<FArrayBox> total_charge_density( grids, 1, IntVect::Zero );
               computeTotalChargeDensity( total_charge_density, a_kinetic_species, a_fluid_species );
               m_poisson->computePotential( a_phi, total_charge_density );
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
               }
               
               if (single_null) {
                  ((NewGKPoissonBoltzmann*)m_poisson)
                     ->computePotentialAndElectronDensity( a_phi,
                                                           *m_boltzmann_electron,
                                                           ion_charge_density,
                                                           a_bc,
                                                           a_initial_time );
               }
               else {
                  ((GKPoissonBoltzmann*)m_poisson)
                     ->computePotentialAndElectronDensity( a_phi,
                                                           *m_boltzmann_electron,
                                                           ion_charge_density,
                                                           a_bc,
                                                           a_initial_time );
               }
            }
         }

         // Filter out all harmonics except for the primary mode
         if (m_apply_harm_filtering) {
            applyHarmonicFiltering(a_phi, m_harm_filtering_dir);
         }
      }

      if ( !m_fixed_efield || a_initial_time ) {
         fillInternalGhosts(a_phi);
         m_poisson->computeField( a_phi, E_field_cell );
         m_poisson->computeField( a_phi, E_field_face );
         
         // Update nodal phi if supporting the calculation of a divergence-free phase velocity
         interpToNodes(a_phi);
      }
   }
   else {
      // Arrival here means that the define() member was never called.  Returning zero.
      for (DataIterator dit( E_field_face.dataIterator() ); dit.ok(); ++dit) {
         E_field_cell[dit].setVal(0.);
         E_field_face[dit].setVal(0.);
      }
   }
   
   // Compute injected E-field
   if (a_inject_field) {
      computeInjectedField(m_injected_E_field, a_kinetic_species, E_field_face, E_field_cell);
   }
}


void EField::computeInjectedField(PS::LevelData<PS::FluxBox>& a_injected_E_field,
                                  const PS::KineticSpeciesPtrVect& a_kinetic_species,
                                  const LevelData<FluxBox>& a_E_field_face,
                                  const LevelData<FArrayBox>& a_E_field_cell)
{
  if (a_kinetic_species.size() > 0 ) {
      const PS::PhaseGeom& phaseGeom( a_kinetic_species[0]->phaseSpaceGeometry() );

#if 0
      phaseGeom.injectConfigurationToPhase(a_E_field_face,
                                              a_E_field_cell,
                                              a_injected_E_field);
      
#else
      //This code provides some performance optimization
      if (!phaseGeom.secondOrder()) {
         phaseGeom.injectConfigurationToPhase(a_E_field_face,
                                              a_E_field_cell,
                                              a_injected_E_field);
      }
      else {
         //Strip ghost since they are not needed for 2nd order
	 //Work with assertions in field calculations later 
         const MagGeom& mag_geom = configurationSpaceGeometry();
         const DisjointBoxLayout& grids = mag_geom.gridsFull();
         LevelData<FArrayBox> E_cell_tmp(grids, 3, IntVect::Zero);
         LevelData<FluxBox> E_face_tmp(grids, 3, IntVect::Zero);
         for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
            E_cell_tmp[dit].copy(a_E_field_cell[dit]);
            E_face_tmp[dit].copy(a_E_field_face[dit]);
         }

         phaseGeom.injectConfigurationToPhase(E_face_tmp,
                                              E_cell_tmp,
                                              a_injected_E_field,
                                              false );

      }
#endif
   }
}

void EField::setCoreBC( const double   a_core_inner_bv,
                        const double   a_core_outer_bv,
                        EllipticOpBC&  a_bc ) const 
{
   const MagGeom& mag_geom = configurationSpaceGeometry();

   if ( typeid(*(mag_geom.getCoordSys())) == typeid(LogRectCoordSys) ) {
      a_bc.setBCValue(0,RADIAL_DIR,0,a_core_inner_bv);
      a_bc.setBCValue(0,RADIAL_DIR,1,a_core_outer_bv);
   }
   else if ( typeid(*(mag_geom.getCoordSys())) == typeid(SNCoreCoordSys) ) {
      a_bc.setBCValue(SNCoreBlockCoordSys::LCORE,RADIAL_DIR,0,a_core_inner_bv);
      a_bc.setBCValue(SNCoreBlockCoordSys::LCORE,RADIAL_DIR,1,a_core_outer_bv);
   }
   else if ( typeid(*(mag_geom.getCoordSys())) == typeid(SingleNullCoordSys) ) {
      a_bc.setBCValue(SingleNullBlockCoordSys::LCORE,RADIAL_DIR,0,a_core_inner_bv);
   }
   else {
      MayDay::Error("EField::setCoreBC(): unknown geometry with consistent bcs");
   }
}


RefCountedPtr<CFGVars>
EField::clone( const bool a_copy_data ) const
{
   EField* efield = new EField( pp_prefix(), name(), configurationSpaceGeometry(), IntVect::Unit );

   if (a_copy_data) {
      LevelData<FArrayBox>& dst_cell_data = efield->cell_var(0);
      const LevelData<FArrayBox>& src_cell_data = cell_var(0);
      for (DataIterator dit( src_cell_data.dataIterator() ); dit.ok(); ++dit) {
         dst_cell_data[dit].copy( src_cell_data[dit] );
      }

      LevelData<FluxBox>& dst_face_data = efield->face_var(0);
      const LevelData<FluxBox>& src_face_data = face_var(0);

      for (DataIterator dit( src_face_data.dataIterator() ); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            dst_face_data[dit][dir].copy( src_face_data[dit][dir] );
         }
      }
   }

   LevelData<FArrayBox>& dst_phi_node = efield->m_phi_node;
   const LevelData<FArrayBox>& src_phi_node = m_phi_node;
   if ( src_phi_node.isDefined() ) {
      dst_phi_node.define( src_phi_node.disjointBoxLayout(),
                           src_phi_node.nComp(),
                           src_phi_node.ghostVect() );
   }

   efield->m_boltzmann_electron = m_boltzmann_electron;
   efield->m_fixed_efield = m_fixed_efield;
   efield->m_poisson = m_poisson;
   efield->m_defined = true;

   return RefCountedPtr<CFGVars>(efield);
}


void EField::computeQuasiNeutralElectronDensity( LevelData<FArrayBox>&        a_quasineutral_density,
                                                 LevelData<FArrayBox>&        a_potential,
                                                 const EllipticOpBC&          a_bc, 
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


void EField::computeIonChargeDensity( LevelData<FArrayBox>&             a_ion_charge_density,
                                      const PS::KineticSpeciesPtrVect&  a_species) const
{
   // Container for individual species charge density
   const MagGeom& mag_geom = configurationSpaceGeometry();
   const DisjointBoxLayout& mag_grids = mag_geom.gridsFull();
   LevelData<FArrayBox> species_charge_density( mag_grids, 1, IntVect::Zero );

   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      a_ion_charge_density[dit].setVal(0.);
   }

   // Accumulate contribution from kinetic species; add fluid species contribution later
   for (int species(0); species<a_species.size(); species++) {

      const PS::KineticSpecies& this_species( *(a_species[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      // Compute the charge density for this species
      if (this_species.isGyrokinetic()) {
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


void EField::computePolarizationMassDensity(LevelData<FArrayBox>&             a_mass_density,
                                            const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                            const FluidSpeciesPtrVect&        a_fluid_species) const
{
   // Container for individual species charge density
   const MagGeom& mag_geom = configurationSpaceGeometry();
   const DisjointBoxLayout& mag_grids = mag_geom.gridsFull();
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
         a_mass_density[dit].plus( species_mass_density[dit] );
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
            a_mass_density[dit].plus( species_mass_density[dit] );
         }
      }
   }
}


void EField::computeTotalChargeDensity( LevelData<FArrayBox>&             a_charge_density,
                                        const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                        const FluidSpeciesPtrVect&        a_fluid_species) const
{
   // Container for individual species charge density
   const MagGeom& mag_geom = configurationSpaceGeometry();
   const DisjointBoxLayout& mag_grids = mag_geom.gridsFull();
   LevelData<FArrayBox> species_charge_density( mag_grids, 1, IntVect::Zero );

   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      a_charge_density[dit].setVal(0.);
   }

   // Accumulate contribution from kinetic species
   for (int species(0); species<a_kinetic_species.size(); species++) {
      const PS::KineticSpecies& this_species( *(a_kinetic_species[species]) );

      // Compute the charge density for this species
      if (this_species.isGyrokinetic()) {
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


void EField::computeIonParallelCurrentDensity( LevelData<FArrayBox>&             a_ion_current_density,
                                               const PS::KineticSpeciesPtrVect&  a_species ) const
{
   // Container for individual species charge density
   const MagGeom& mag_geom = configurationSpaceGeometry();
   const DisjointBoxLayout& mag_grids = mag_geom.gridsFull();
   LevelData<FArrayBox> species_current_density( mag_grids, 1, IntVect::Zero );

   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      a_ion_current_density[dit].setVal(0.);
   }

   for (int species(0); species<a_species.size(); species++) {

      const PS::KineticSpecies& this_species( *(a_species[species]) );
      if ( this_species.charge() < 0.0 ) continue;
      
      // Compute the charge density for this species
      this_species.ParallelMomentum( species_current_density );
      
      DataIterator dit( a_ion_current_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         species_current_density[dit].mult( this_species.charge() );
         a_ion_current_density[dit].plus( species_current_density[dit] );
      }
   }
}

void
EField::applyHarmonicFiltering(LevelData<FArrayBox>& a_phi,
                               const int& a_dir) const
{
  CH_TIME("EField::applyHarmonicFiltering");
  
   Real pi = Constants::PI;
   const DisjointBoxLayout& grids = a_phi.getBoxes();
   
   LevelData<FArrayBox> tmp;
   tmp.define(a_phi);
   
   Box domain_box = grids.physDomain().domainBox();
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      int Npts = grids[dit].size(a_dir);
      
      if (Npts != domain_box.size(a_dir)) {
         MayDay::Error("EField::applyHarmonicFiltering(): the direction of filtering must be single-block and has no domain decomposition");
      }
      
      for (BoxIterator bit( grids[dit] ); bit.ok(); ++bit) {
         IntVect iv( bit() );
         
         double A = 0.0;
         double B = 0.0;
         for (int n=0; n<Npts; ++n) {
            IntVect ivSum(iv);
            ivSum[a_dir] = n;
            double phase = 2.0 * pi * (n+0.5)/Npts;
            A += (2.0/Npts) * tmp[dit](ivSum,0) * cos(phase);
            B += (2.0/Npts) * tmp[dit](ivSum,0) * sin(phase);
         }
         double phase = 2.0 * pi * (iv[a_dir]+0.5)/Npts;
         a_phi[dit](iv,0) = A*cos(phase)+B*sin(phase);
      }
   }
}


void
EField::parseParameters( ParmParse& a_pp)
{
  a_pp.query( "include_ion_polarization_density", m_include_pol_dens_i);
  a_pp.query( "include_electron_polarization_density", m_include_pol_dens_e);
   
  if (a_pp.contains("harmonic_filtering")) {
    a_pp.get("harmonic_filtering", m_apply_harm_filtering);
  }
  else {
    m_apply_harm_filtering = false;
  }
  
  if (m_apply_harm_filtering) {
    a_pp.get("harmonic_filtering_dir", m_harm_filtering_dir);
  }
}

#include "NamespaceFooter.H"
