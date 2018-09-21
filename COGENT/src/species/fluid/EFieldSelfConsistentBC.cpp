#include "EFieldSelfConsistentBC.H"
#include "GKPoissonBoltzmann.H"
#include "NewGKPoissonBoltzmann.H"

#include "NamespaceHeader.H"


EFieldSelfConsistentBC::EFieldSelfConsistentBC( const string&      a_pp_prefix,
                                                const std::string& a_name,
                                                const MagGeom&     a_geometry,
                                                const IntVect&     a_ghost_vect )
   : EField(a_pp_prefix, a_name, a_geometry, a_ghost_vect)
{
   // This is actually just a scratch variable expected to be filled by calls to the
   // solver setOperatorCoefficients() members
   m_radial_gkp_divergence_average.define(a_geometry.gridsFull(), 1, IntVect::Zero);
}


void EFieldSelfConsistentBC::computeEField( const PS::GKState&                a_state,
                                            const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                            const FluidSpeciesPtrVect&        a_fluid_species,
                                            const PS::ScalarPtrVect&          a_scalars,
                                            LevelData<FArrayBox>&             a_phi,
                                            EllipticOpBC&                     a_bc,
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
      for (DataIterator dit( E_field_face.dataIterator() ); dit.ok(); ++dit) {
         E_field_cell[dit].setVal(0.);
         E_field_face[dit].setVal(0.);
      }
   }
}


void EFieldSelfConsistentBC::updateImplicitPotential(LevelData<FArrayBox>&             a_phi,
                                                     const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                                     const Vector<Real>&               a_scalar_data,
                                                     LevelData<FArrayBox>&             a_divJperp,
                                                     EllipticOpBC&                     a_bc,
                                                     const Real                        a_dt )
{
   setCoreBC(a_scalar_data[0], -a_scalar_data[1], a_bc);

   EField::updateImplicitPotential(a_phi, a_kinetic_species, a_scalar_data, a_divJperp, a_bc, a_dt );
}


RefCountedPtr<CFGVar>
EFieldSelfConsistentBC::clone( const bool a_copy_data ) const
{
   EFieldSelfConsistentBC* efield = new EFieldSelfConsistentBC( pp_prefix(), name(), configurationSpaceGeometry(), IntVect::Unit );

   if (a_copy_data) {
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

   efield->m_boltzmann_electron = m_boltzmann_electron;
   efield->m_fixed_efield = m_fixed_efield;
   efield->m_poisson = m_poisson;
   efield->m_radial_gkp_divergence_average_lo = m_radial_gkp_divergence_average_lo;
   efield->m_radial_gkp_divergence_average_hi = m_radial_gkp_divergence_average_hi;
   efield->m_defined = true;

   return RefCountedPtr<CFGVar>(efield);
}


#include "NamespaceFooter.H"
