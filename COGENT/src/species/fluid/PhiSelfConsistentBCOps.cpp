#include "PhiSelfConsistentBCOps.H"
#include "GKPoissonBoltzmann.H"
#include "NewGKPoissonBoltzmann.H"

#include "NamespaceHeader.H"


PhiSelfConsistentBCOps::PhiSelfConsistentBCOps( const string&      a_pp_prefix,
                                                const MagGeom&     a_geometry )
   : PhiOps(a_pp_prefix, a_geometry)
{
   // This is actually just a scratch variable expected to be filled by calls to the
   // solver setOperatorCoefficients() members
   m_radial_gkp_divergence_average.define(a_geometry.gridsFull(), 1, IntVect::Zero);
}


void PhiSelfConsistentBCOps::updateEField( EMFields&                         a_EM_fields,
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

   if (m_poisson) {

      if ( !m_fixed_efield ) {

         const DisjointBoxLayout& grids = m_geometry.gridsFull();

         bool single_null = typeid(*(m_geometry.getCoordSys())) == typeid(SingleNullCoordSys);

         // Update the potential and field, if not fixed_efield

         LevelData<FArrayBox> ion_mass_density( grids, 1, IntVect::Zero );
         computePolarizationMassDensity( ion_mass_density, a_kinetic_species, a_fluid_species );

         int component = -1;
         for (int i=0; i<a_scalars.size(); ++i) {
            if (a_scalars[i]->name() == "Er_boundary") {
               component = i;
               break;
            }
         }
         CH_assert(component >= 0);
         const Vector<Real>& scalar_data = a_scalars[component]->data();

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
            m_poisson->setOperatorCoefficients( ion_mass_density, a_bc, a_update_potential,
                                                m_radial_gkp_divergence_average_lo,
                                                m_radial_gkp_divergence_average_hi,
                                                m_radial_gkp_divergence_average );
         }

         if ( a_update_potential ) {
         
            if (m_boltzmann_electron == NULL) {
            
               LevelData<FArrayBox> gkPoissonRHS( grids, 1, IntVect::Zero );
               computeTotalChargeDensity( gkPoissonRHS, a_kinetic_species, a_fluid_species );

               for (DataIterator dit(grids); dit.ok(); ++dit) {
                  phi[dit].setVal(0.);
               }
               m_poisson->solveWithBCs( phi, gkPoissonRHS );
            }
            else {

               // Boltzmann electron model

               LevelData<FArrayBox> ion_charge_density( grids, 1, IntVect::Zero );
               computeIonChargeDensity( ion_charge_density, a_kinetic_species );

               for (DataIterator dit(grids); dit.ok(); ++dit) {
                  phi[dit].setVal(0.);
               }

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


#include "NamespaceFooter.H"
