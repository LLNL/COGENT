#include "EFieldAmpere.H"
#include "NewGKPoissonBoltzmann.H"
#include "Directions.H"

#include "NamespaceHeader.H"


EFieldAmpere::EFieldAmpere( const string&      a_pp_prefix,
                            const std::string& a_name,
                            const MagGeom&     a_geometry,
                            const IntVect&     a_ghost_vect )
   : EField(a_pp_prefix, a_name, a_geometry, a_ghost_vect)
{
   ParmParse ppgksys("gksystem");

   if (ppgksys.contains("ampere_cold_electrons")) {
      ppgksys.query( "ampere_cold_electrons", m_cold_electrons );
   }
   else {
      m_cold_electrons = false;
   }

   m_radial_gkp_divergence_average.define(a_geometry.gridsFull(), 1, IntVect::Zero);
}


void EFieldAmpere::define( const double                      a_larmor,
                           const double                      a_debye,
                           const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                           BoltzmannElectron*                a_boltzmann_electron,
                           const bool                        a_fixed_efield,
                           const bool                        a_support_divfree_phase_vel,
                           const int                         a_cur_step )
{
   EField::define(a_larmor, a_debye, a_kinetic_species, a_boltzmann_electron, a_fixed_efield,
                  a_support_divfree_phase_vel, a_cur_step);

   const MagGeom& mag_geom = configurationSpaceGeometry();
   const DisjointBoxLayout& grids = mag_geom.gridsFull();

   m_Er_average_face.define( grids, 3, IntVect::Unit );
   m_Er_average_cell.define( grids, 3, IntVect::Unit );
   m_E_tilde_face.define( grids, 3, IntVect::Unit );
   m_E_tilde_cell.define( grids, 3, IntVect::Unit );

   if ( a_cur_step == 0 ) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         m_Er_average_face[dit].setVal(0.);
         m_Er_average_cell[dit].setVal(0.);
         m_E_tilde_face[dit].setVal(0.);
         m_E_tilde_cell[dit].setVal(0.);
      }
   }
}


void EFieldAmpere::computeEField( const PS::GKState&                a_state,
                                  const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                  const FluidSpeciesPtrVect&        a_fluid_species,
                                  const PS::ScalarPtrVect&          a_scalars,
                                  LevelData<FArrayBox>&             a_phi,
                                  PotentialBC&                      a_bc,
                                  const bool                        a_update_potential,
                                  const bool                        a_initial_time )
{
   LevelData<FArrayBox>& E_field_cell = cell_data();
   CH_assert(E_field_cell.ghostVect() == IntVect::Unit);
   LevelData<FluxBox>& E_field_face = face_data();
   CH_assert(E_field_face.ghostVect() == IntVect::Unit);

   if (m_poisson) {

      if ( !m_fixed_efield ) {

         const MagGeom& mag_geom = configurationSpaceGeometry();
         const DisjointBoxLayout& grids = mag_geom.gridsFull();

         bool single_null = typeid(*(mag_geom.getCoordSys())) == typeid(SingleNullCoordSys);

         // Update the potential and field, if not fixed_efield

         LevelData<FArrayBox> ion_mass_density( grids, 1, IntVect::Zero );
         computeIonMassDensity( ion_mass_density, a_kinetic_species );

         const Vector<Real>& scalar_data = a_scalars[a_state.getScalarComponent("Er_boundary")]->data();
         double Er_lo = scalar_data[0];
         double Er_hi = scalar_data[1];

         setCoreBC( Er_lo, -Er_hi, a_bc );

         if ( m_boltzmann_electron ) {
            if ( single_null ) {
               double core_outer_bv = -Er_hi;
               
               ((NewGKPoissonBoltzmann*)m_poisson)
                  ->setOperatorCoefficients( ion_mass_density, a_bc, core_outer_bv, m_radial_gkp_divergence_average_lo,
                                             m_radial_gkp_divergence_average_hi, m_radial_gkp_divergence_average );
            }
            else {
               ((GKPoissonBoltzmann*)m_poisson)
                  ->setOperatorCoefficients( ion_mass_density, a_bc, m_radial_gkp_divergence_average_lo,
                                             m_radial_gkp_divergence_average_hi, m_radial_gkp_divergence_average );
            }
         }
         else {
            m_poisson->setOperatorCoefficients( ion_mass_density, a_bc, true, m_radial_gkp_divergence_average_lo,
                                                m_radial_gkp_divergence_average_hi, m_radial_gkp_divergence_average );
         }

         if ( a_update_potential ) {

            if (m_boltzmann_electron == NULL) {

               LevelData<FArrayBox> gkPoissonRHS( grids, 1, IntVect::Zero );
               computeTotalChargeDensity( gkPoissonRHS, a_kinetic_species );

               m_poisson->computePotential( a_phi, gkPoissonRHS );
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
            }
         }
      }

      if ( !m_fixed_efield || a_initial_time ) {

         fillInternalGhosts(a_phi);

         for (DataIterator dit(E_field_cell.dataIterator()); dit.ok(); ++dit) {
            E_field_face[dit].copy(m_Er_average_face[dit]);
            E_field_face[dit] += m_E_tilde_face[dit];
            E_field_cell[dit].copy(m_Er_average_cell[dit]);
            E_field_cell[dit] += m_E_tilde_cell[dit];
         }

         // Update nodal phi if supporting the calculation of a divergence-free phase velocity
         interpToNodes(a_phi);
      }
   }
   else {
      for (DataIterator dit( E_field_face.dataIterator() ); dit.ok(); ++dit) {
         E_field_cell[dit].setVal(0.);
         E_field_face[dit].setVal(0.);
      }
   }
}


void EFieldAmpere::updateAveragedEfield( LevelData<FluxBox>&          a_Er_average_face,
                                         LevelData<FArrayBox>&        a_Er_average_cell,
                                         const LevelData<FArrayBox>&  a_Er,
                                         const double                 a_Er_lo,
                                         const double                 a_Er_hi ) const
{
   CH_assert(a_Er_average_face.ghostVect() == IntVect::Unit);
   CH_assert(a_Er_average_cell.ghostVect() == IntVect::Unit);

   bool m_Esol_extrapolation = true;
   bool m_dealignment_corrections = false;

   //Geometry parameters
   const MagGeom& mag_geom( configurationSpaceGeometry() );
   const MagCoordSys& coords = *mag_geom.getCoordSys();
   const DisjointBoxLayout& mag_grids = mag_geom.gridsFull();
    
   //Create face-centered increment for the mapped E-field (- grad Phi )
   LevelData<FluxBox> Er_mapped_face(mag_grids, 3, IntVect::Unit );
    
   //Creating tmp arrays with extra layers of ghost cells
   LevelData<FArrayBox> flux_divergence_tmp( mag_grids, 1, 2*IntVect::Unit );
    
   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      flux_divergence_tmp[dit].setVal(0.);
   }
    
   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      flux_divergence_tmp[dit].copy(a_Er[dit]);
   }
    
   flux_divergence_tmp.exchange();
    
   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
        
      Er_mapped_face[dit].setVal(0.0);
        
      FArrayBox& this_flux_divergence = flux_divergence_tmp[dit];
        
      int block_number( coords.whichBlock(mag_grids[dit]) );
      const MagBlockCoordSys& block_coord_sys = (const MagBlockCoordSys&)(*(coords.getCoordSys(block_number)));
      int lo_radial_index = block_coord_sys.domain().domainBox().smallEnd(RADIAL_DIR);
      int hi_radial_index = block_coord_sys.domain().domainBox().bigEnd(RADIAL_DIR);
        
      //This is the main calculation of Er on closed flux surfaces (within the core).
      if ( block_number <= RCORE || block_number == MCORE ) {
            
         for (int dir = 0; dir < 2; dir++) {
                
            FArrayBox& this_Er_mapped_dir = Er_mapped_face[dit][dir];
                
            Box box( this_Er_mapped_dir.box() );
            BoxIterator bit(box);
            for (bit.begin(); bit.ok(); ++bit) {
                    
               IntVect iv = bit();
                    
               if (iv[0]<lo_radial_index) {
                  this_Er_mapped_dir(iv,0) = a_Er_lo;
               }
               else if (iv[0]>hi_radial_index) {
                  this_Er_mapped_dir(iv,0) = a_Er_hi;
               }
               else {
                        
                  IntVect iv_tmp = bit();
                  iv_tmp[1] = mag_grids[dit].smallEnd(POLOIDAL_DIR);
                        
                  if ( dir == RADIAL_DIR ) {
                     this_Er_mapped_dir(iv,0) =   this_flux_divergence(iv_tmp,0);
                  }
                        
                  if ( dir == POLOIDAL_DIR ) {
                     IntVect iv_incr = iv_tmp;
                     iv_incr[0] = iv_tmp[0] + 1;
                     if (iv[0]<hi_radial_index) {
                        this_Er_mapped_dir(iv,0) =   0.5 * ( this_flux_divergence(iv_tmp,0) + this_flux_divergence(iv_incr,0) );
                     }
                     else  {
                        this_Er_mapped_dir(iv,0) =   0.5 * ( this_flux_divergence(iv_tmp,0) + a_Er_hi );
                     }
                  }
               }
            }
         }
      }
        
      //Linear extrapolation of Er into the SOL blocks (here, assumed that the cut is vertical, i.e., dPsi/dz = dPsi/dr)
      //To provide continuity, the mapped field is scaled by |grad_Psi|, i.e., Emapped_sol = Emapped_core * gradPsi_core / gradPsi_sol
      //Fix later for 10 block and DIII-D
      else if (m_Esol_extrapolation && ((block_number < LPF) || (block_number == MCORE) || (block_number == MCSOL))) {
         const MagBlockCoordSys& lcore_coord_sys = (const MagBlockCoordSys&)(*(coords.getCoordSys(LCORE)));
         const MagBlockCoordSys& lsol_coord_sys = (const MagBlockCoordSys&)(*(coords.getCoordSys(LCSOL)));
         double lo_rad_end = block_coord_sys.lowerMappedCoordinate(0);
         double hi_rad_end = block_coord_sys.upperMappedCoordinate(0);
            
         RealVect topSep_core;
         topSep_core[0] = lcore_coord_sys.upperMappedCoordinate(0);
         topSep_core[1] = lcore_coord_sys.lowerMappedCoordinate(1);
            
         RealVect topSep_sol;
         topSep_sol[0] = lsol_coord_sys.lowerMappedCoordinate(0);
         topSep_sol[1] = lsol_coord_sys.lowerMappedCoordinate(1);
            
         Real dZdR_core = lcore_coord_sys.dXdXi(topSep_core,1,1)/lcore_coord_sys.dXdXi(topSep_core,0,1);
         Real dPsidZ_core = 1.0/(lcore_coord_sys.dXdXi(topSep_core,1,0) - lcore_coord_sys.dXdXi(topSep_core,0,0)*dZdR_core);
            
         Real dZdR_sol = lsol_coord_sys.dXdXi(topSep_sol,1,1)/lsol_coord_sys.dXdXi(topSep_sol,0,1);
         Real dPsidZ_sol = 1.0/(lsol_coord_sys.dXdXi(topSep_sol,1,0) - lsol_coord_sys.dXdXi(topSep_sol,0,0)*dZdR_sol);
            
         for (int dir = 0; dir < 2; dir++) {
            FArrayBox& this_Er_mapped_dir = Er_mapped_face[dit][dir];
                
            Box box( this_Er_mapped_dir.box() );
            FArrayBox xi(box,2);
            block_coord_sys.getFaceCenteredMappedCoords(dir, xi);
                
            BoxIterator bit(box);
            for (bit.begin(); bit.ok(); ++bit) {
               IntVect iv = bit();
               double amplitude = a_Er_hi;
               amplitude *= dPsidZ_core / dPsidZ_sol;
               this_Er_mapped_dir(iv,0) =  amplitude  * (1.0 - (xi(iv,0) - lo_rad_end)/(hi_rad_end - lo_rad_end) );
            }
         }
      }
      else {
         double E_open = 0.0;
         Er_mapped_face[dit].setVal(E_open,RADIAL_DIR,0,1);
         Er_mapped_face[dit].setVal(E_open,POLOIDAL_DIR,0,1);
      }
   }

   //Compute cell-centered increment for the mapped E-field
   LevelData<FArrayBox> Er_mapped_cell(mag_grids, 3, IntVect::Unit );
   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      Box box( Er_mapped_cell[dit].box() );
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); ++bit) {
         IntVect iv = bit();
         IntVect iv_r = bit();
         iv_r[0] = iv[0] + 1;
         Er_mapped_cell[dit](iv,0) = 0.5*(Er_mapped_face[dit][0](iv,0)+Er_mapped_face[dit][0](iv_r,0));
         Er_mapped_cell[dit](iv,1) = 0.5*(Er_mapped_face[dit][0](iv,1)+Er_mapped_face[dit][0](iv_r,1));
         Er_mapped_cell[dit](iv,2) = 0.5*(Er_mapped_face[dit][0](iv,2)+Er_mapped_face[dit][0](iv_r,2));
      }
   }
    
   //Multiply by NJtranspose (E_fs_average = - NJInverse * grad Phi)
   LevelData<FluxBox> Er_face_incr( mag_grids, 3, IntVect::Unit );
   mag_geom.unmapGradient(Er_mapped_face, Er_face_incr);
    
   LevelData<FArrayBox> Er_cell_incr( mag_grids, 3, IntVect::Unit );
   mag_geom.unmapGradient(Er_mapped_cell, Er_cell_incr);
    
   //Update E-field
   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      a_Er_average_face[dit].copy(Er_face_incr[dit]);
      a_Er_average_cell[dit].copy(Er_cell_incr[dit]);
   }
    
   //Improve Er field calculation to take into accout the dealigment between the grid and magnetic surfaces
   if ( (typeid(coords) == typeid(SingleNullCoordSys)) && (m_dealignment_corrections)) {
        
      mag_geom.interpolateErFromMagFS(a_Er_average_face, a_Er_average_cell);
        
      //If not doing extrapolation, zero out efield in the PF region,
      //which appeares in the present version of interpolateFromMagFS
      //as symmetric reflection of E in the core region
      if (!m_Esol_extrapolation) {
         for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
            int block_number( coords.whichBlock( mag_grids[dit] ) );
            if ((block_number == RPF) || (block_number == LPF)) {
               a_Er_average_face[dit].setVal(0.0);
               a_Er_average_cell[dit].setVal(0.0);
            }
         }
      }
   }
}


void EFieldAmpere::updateEfieldPoloidalVariation( const LevelData<FluxBox>&   a_E_tilde_mapped_face,
                                                  const LevelData<FArrayBox>& a_E_tilde_mapped_cell,
                                                  LevelData<FluxBox>&         a_E_tilde_phys_face,
                                                  LevelData<FArrayBox>&       a_E_tilde_phys_cell,
                                                  LevelData<FArrayBox>&       a_phi_tilde_fs_average,
                                                  double&                     a_phi_tilde_fs_average_lo,
                                                  double&                     a_phi_tilde_fs_average_hi,
                                                  const LevelData<FArrayBox>& a_radial_gkp_divergence_average,
                                                  const double                a_radial_gkp_divergence_average_lo,
                                                  const double                a_radial_gkp_divergence_average_hi ) const
{
   CH_assert(a_E_tilde_phys_face.ghostVect() == IntVect::Unit);
   CH_assert(a_E_tilde_phys_cell.ghostVect() == IntVect::Unit);

   const MagGeom& mag_geom( configurationSpaceGeometry() );
   const DisjointBoxLayout& mag_grids = mag_geom.gridsFull();

   for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
      a_E_tilde_phys_face[dit].setVal(0.0);
      a_E_tilde_phys_cell[dit].setVal(0.0);
   }

   if ( !m_cold_electrons ) {
      CH_assert(a_E_tilde_mapped_face.ghostVect() == IntVect::Unit);
      CH_assert(a_E_tilde_mapped_cell.ghostVect() == IntVect::Unit);
      CH_assert(a_phi_tilde_fs_average.ghostVect() == IntVect::Zero);

      mag_geom.unmapGradient(a_E_tilde_mapped_face, a_E_tilde_phys_face);
      mag_geom.unmapGradient(a_E_tilde_mapped_cell, a_E_tilde_phys_cell);
      
      LevelData<FluxBox> Er_average_face_pol_contrib(mag_grids, 3, IntVect::Unit);
      LevelData<FArrayBox> Er_average_cell_pol_contrib(mag_grids, 3, IntVect::Unit);
 
      LevelData<FArrayBox> flux_div_ratio(mag_grids, 1, IntVect::Zero);

      for ( DataIterator dit(mag_grids); dit.ok(); ++dit ) {
         flux_div_ratio[dit].copy(a_phi_tilde_fs_average[dit]);
         flux_div_ratio[dit] /= a_radial_gkp_divergence_average[dit];
      }

      double lo_ratio = a_phi_tilde_fs_average_lo / a_radial_gkp_divergence_average_lo;
      double hi_ratio = a_phi_tilde_fs_average_hi / a_radial_gkp_divergence_average_hi;

      updateAveragedEfield( Er_average_face_pol_contrib, Er_average_cell_pol_contrib, flux_div_ratio, lo_ratio, hi_ratio );
       
      for (DataIterator dit(mag_grids); dit.ok(); ++dit) {
         a_E_tilde_phys_face[dit] -= Er_average_face_pol_contrib[dit];
         a_E_tilde_phys_cell[dit] -= Er_average_cell_pol_contrib[dit];
      }
   }
}


void EFieldAmpere::updateAverageAndPerturbation( const Vector<Real>&               a_scalar_data,
                                                 const LevelData<FArrayBox>&       a_Er,
                                                 const PS::KineticSpeciesPtrVect&  a_species_phys,
                                                 const Real                        a_time )
{
   double Er_lo = a_scalar_data[0];
   double Er_hi = a_scalar_data[1];
   if (procID()==0) cout << "  Er_lo = " << Er_lo << ", Er_hi = " << Er_hi << endl;

   LevelData<FluxBox> E_tilde_mapped_face;
   LevelData<FArrayBox> E_tilde_mapped_cell;
   LevelData<FArrayBox> phi_tilde_fs_average;
   double phi_tilde_fs_average_lo;
   double phi_tilde_fs_average_hi;

   if ( !m_cold_electrons ) {

      const MagGeom& mag_geom = configurationSpaceGeometry();
      const DisjointBoxLayout& mag_grids = mag_geom.gridsFull();

      // Compute the poloidal variation phi_tilde.  The poloidal variation of the
      // current potential is used as the initial guess for the interative procedure.
      int num_phi_ghosts_filled = m_poisson->numPotentialGhosts();
      LevelData<FArrayBox> phi_tilde(mag_grids, 1, num_phi_ghosts_filled*IntVect::Unit);

      computePhiTilde( a_species_phys, *m_boltzmann_electron, phi_tilde );

      // Fill the three ghost cell layers at block boundaries as required by the
      // subsequent field calculations
      fillInternalGhosts( phi_tilde );

      E_tilde_mapped_face.define(mag_grids, 3, IntVect::Unit);
      m_poisson->computeMappedField( phi_tilde, E_tilde_mapped_face );

      E_tilde_mapped_cell.define(mag_grids, 3, IntVect::Unit);
      m_poisson->computeMappedField( phi_tilde, E_tilde_mapped_cell );

      LevelData<FluxBox> Er_tilde_mapped(mag_grids, 1, IntVect::Zero);
      m_poisson->extractNormalComponent(E_tilde_mapped_face, RADIAL_DIR, Er_tilde_mapped);

      phi_tilde_fs_average.define(mag_grids, 1, IntVect::Zero);
      m_poisson->computeRadialFSAverage( Er_tilde_mapped, phi_tilde_fs_average_lo, phi_tilde_fs_average_hi, phi_tilde_fs_average );
   }

   updateAveragedEfield( m_Er_average_face, m_Er_average_cell, a_Er, Er_lo, Er_hi );
          
   updateEfieldPoloidalVariation( E_tilde_mapped_face, E_tilde_mapped_cell, m_E_tilde_face, m_E_tilde_cell,
                                  phi_tilde_fs_average, phi_tilde_fs_average_lo, phi_tilde_fs_average_hi,
                                  m_radial_gkp_divergence_average, m_radial_gkp_divergence_average_lo, m_radial_gkp_divergence_average_hi );
}


void EFieldAmpere::computePhiTilde( const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                    const BoltzmannElectron&          a_ne,
                                    LevelData<FArrayBox>&             a_phi_tilde ) const
{
   const MagGeom& mag_geom = configurationSpaceGeometry();
   const DisjointBoxLayout& mag_grids = mag_geom.gridsFull();
   LevelData<FArrayBox> Zni(mag_grids, 1, IntVect::Zero);
   computeIonChargeDensity( Zni, a_kinetic_species );

   if ( GKPoissonBoltzmann* gkpb = dynamic_cast<GKPoissonBoltzmann*>(m_poisson) ) {
      gkpb->getPhiTilde( Zni, a_ne, a_phi_tilde );
   }
   else if ( NewGKPoissonBoltzmann* new_gkpb = dynamic_cast<NewGKPoissonBoltzmann*>(m_poisson) ) {
      new_gkpb->getPhiTilde( Zni, a_ne, a_phi_tilde );
   }
   else {
      MayDay::Error("GKOps::computeGKPPhiTildeDivergence() can only be called with Boltzmann electrons");
   }
}


void EFieldAmpere::updateImplicitPotential(LevelData<FArrayBox>&             a_phi,
                                           const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                           const Vector<Real>&               a_scalar_data,
                                           LevelData<FArrayBox>&             a_divJperp,
                                           PotentialBC&                      a_bc,
                                           const Real                        a_dt )
{
   double Er_lo = a_scalar_data[0];
   double Er_hi = a_scalar_data[1];

   setCoreBC( Er_lo, -Er_hi, a_bc );

   EField::updateImplicitPotential(a_phi, a_kinetic_species, a_scalar_data, a_divJperp, a_bc, a_dt );
}


RefCountedPtr<CFGVar>
EFieldAmpere::clone( const bool a_copy_data ) const
{
   EFieldAmpere* efield = new EFieldAmpere( pp_prefix(), name(), configurationSpaceGeometry(), IntVect::Unit );

   if ( a_copy_data ) {
      LevelData<FArrayBox>& dst_cell_data = efield->cell_data();
      const LevelData<FArrayBox>& src_cell_data = cell_data();
      for (DataIterator dit( src_cell_data.dataIterator() ); dit.ok(); ++dit) {
         dst_cell_data[dit].copy( src_cell_data[dit] );
      }

      LevelData<FluxBox>& dst_face_data = efield->face_data();
      const LevelData<FluxBox>& src_face_data = face_data();
      for (DataIterator dit( src_face_data.dataIterator() ); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            dst_face_data[dit][dir].copy( src_face_data[dit][dir] );
         }
      }

      efield->m_fixed_efield = m_fixed_efield;
      efield->m_cold_electrons = m_cold_electrons;
   }

   LevelData<FArrayBox>& dst_phi_node = efield->m_phi_node;
   const LevelData<FArrayBox>& src_phi_node = m_phi_node;
   if ( src_phi_node.isDefined() ) {
      dst_phi_node.define( src_phi_node.disjointBoxLayout(),
                           src_phi_node.nComp(),
                           src_phi_node.ghostVect() );
   }

   LevelData<FArrayBox>& dst_rad_div = efield->m_radial_gkp_divergence_average;
   const LevelData<FArrayBox>& src_rad_div = m_radial_gkp_divergence_average;
   dst_rad_div.define( src_rad_div.disjointBoxLayout(),
                       src_rad_div.nComp(),
                       src_rad_div.ghostVect() );

   if (m_Er_average_cell.isDefined()) efield->m_Er_average_cell.define(m_Er_average_cell);
   if (m_Er_average_face.isDefined()) efield->m_Er_average_face.define(m_Er_average_face);
   if (m_E_tilde_cell.isDefined()) efield->m_E_tilde_cell.define(m_E_tilde_cell);
   if (m_E_tilde_face.isDefined()) efield->m_E_tilde_face.define(m_E_tilde_face);

   efield->m_boltzmann_electron = m_boltzmann_electron;
   efield->m_poisson = m_poisson;
   efield->m_defined = true;
   
   return RefCountedPtr<CFGVar>(efield);
}


#include "NamespaceFooter.H"
