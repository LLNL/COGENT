#include "VorticityBEOp.H"
#include "EllipticOpBCFactory.H"

#include "NamespaceHeader.H" 

const char* VorticityBEOp::pp_name = {"vorticityBE_op"};

VorticityBEOp::VorticityBEOp(const ParmParse&    a_pp,
                             const MagGeom&      a_geometry,
                             const double        a_larmor,
                             const int           a_verbosity )
   : VorticityOp(a_pp, a_geometry, a_larmor, a_verbosity),
     m_split_op(false),
     m_update_pc_freq(1)

{
   // Parse parameters
   parseParameters(a_pp);
   
   // Print parameters
   if (m_verbosity) {
     printParameters();
   }
   
   // Set up GKVorticity operator
   m_vorticity_op = new GKVorticityBE(a_pp,
                                      a_pp,
                                      a_geometry,
                                      a_larmor,
                                      0.,
                                      m_second_order,
                                      m_low_pollution,
                                      m_include_pol_den_correction,
                                      m_include_diffusion);
   
   m_vorticity_op_bcs = m_potential_bcs->clone();
   
   if (m_include_diffusion) {
      m_vorticity_op->setDiffusionOperatorCoeff(m_pol_diffusion_face,
                                                m_perp_diffusion_face);
   }
   
   if (m_use_vorticity_bcs) {
      m_vorticity_op->setVorticityBc(*m_vorticity_bcs);
   }
   
   
   // Set up GKPoisson operator; N.B. we cannot use
   // the base class member, because we pass a
   // different parsing object here
   m_gyropoisson_op = new GKPoisson(a_pp,
                                    a_geometry,
                                    a_larmor,
                                    0.,
                                    true,
                                    m_second_order,
                                    m_low_pollution );
   
   m_gyropoisson_op_bcs = m_potential_bcs->clone();
   
}


VorticityBEOp::~VorticityBEOp()
{
   delete m_gyropoisson_op;
   delete m_vorticity_op;
}


void VorticityBEOp::initializeBEModel(const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                      const FluidSpeciesPtrVect&        a_fluid_species_phys,
                                      const Vector<Real>&               a_scalar_data )
{
   
   // Need to convert phys fluid species to comp,
   // becasue the rest of the vorticity model
   // (e.g., computeIonMassDensity()) assumes computational fluid species
   FluidSpeciesPtrVect fluid_species_comp;
   int num_fluid_species = a_fluid_species_phys.size();
   fluid_species_comp.resize(num_fluid_species);
   for (int species(0); species<num_fluid_species; species++) {
      fluid_species_comp[species] = a_fluid_species_phys[species]->clone(IntVect::Zero, true);
      fluid_species_comp[species]->convertFromPhysical();
   }
   computeIonMassDensity(m_ion_mass_density, a_kinetic_species_phys, fluid_species_comp);
   fluid_species_comp.clear();
   
   if (a_scalar_data.size() > 0) {
      setCoreBC( a_scalar_data[0], -a_scalar_data[1], *m_gyropoisson_op_bcs );
   }
   
   m_gyropoisson_op->setOperatorCoefficients(m_ion_mass_density, *m_gyropoisson_op_bcs, false);
}


void VorticityBEOp::updatePotentialBEModel(EMFields&                         a_EM_fields,
                                           PhiOps&                           a_phi_ops,
                                           bool                              a_cbc_or_ampere,
                                           const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                           const FluidSpeciesPtrVect&        a_fluid_species_comp,
                                           const PS::ScalarPtrVect&          a_scalar_data,
                                           EllipticOpBC&                     a_potential_bc,
                                           const Real                        a_dt,
                                           const Real                        a_time )
{
   CH_TIME("VorticityBEOp::updatePotentialOldModel");
   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   
   LevelData<FArrayBox> zero(grids, 1, IntVect::Zero);
   setZero(zero);

   const LevelData<FArrayBox>& current_phi = a_EM_fields.getPhi();
   
   // Compute the polarization density term (at the old time) using the operator and boundary
   // conditions set in the previous call of this function or initializeOldModel().
   // NB: gkPoissonRHS = -vorticity
   LevelData<FArrayBox> gkPoissonRHS(grids, 1, IntVect::Zero);
   m_gyropoisson_op->computeFluxDivergence(current_phi, gkPoissonRHS, false);

   // Compute the Reynolds stress term (at the old time)
   // NB: -div(ExB*vorticity) = div(ExB*gkPoissonRHS)
   LevelData<FArrayBox> negative_div_ExB_vorticity(grids, 1, IntVect::Zero);
   if (m_reynolds_stress) {
     m_fluid_op_utils->computeExBAdvection(negative_div_ExB_vorticity,
                                           gkPoissonRHS,
                                           a_EM_fields.getEFieldFace(),
                                           m_advection_scheme,
                                           true);
   }
   
   // Compute the parallel current divergence
   computeIonChargeDensity(m_ion_charge_density, a_kinetic_species_phys, a_fluid_species_comp);
   m_parallel_current_divergence_op->setImplicitDt(a_dt);
   m_parallel_current_divergence_op->setVorticityOperatorCoefficients(zero,
                                                                      m_ion_charge_density,
                                                                      m_electron_temperature,
                                                                      *m_parallel_current_divergence_op_bcs,
                                                                      false);
   m_parallel_current_divergence_op->computeFluxDivergence(m_ion_charge_density, m_negativeDivJpar, false, true);
    
   // Compute the perpendicular current divergence
   LevelData<FArrayBox> div_J_electron_mag(grids, 1, IntVect::Zero);
   if (m_ExB_current_model) {
      computeDivPerpIonExBCurrentDensity(div_J_electron_mag, a_EM_fields, a_kinetic_species_phys, m_ion_charge_density, a_time);
      
   }
   else{
     computeDivPerpElectronMagCurrentDensity(div_J_electron_mag, m_ion_charge_density, gkPoissonRHS, a_time);
   }
   
   LevelData<FArrayBox> div_J_ion_mag(grids, 1, IntVect::Zero);
   computeDivPerpIonMagCurrentDensity(div_J_ion_mag, a_EM_fields, a_kinetic_species_phys, a_time);

   // Get the RHS for the implicit vorticity solve
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      gkPoissonRHS[dit] += m_negativeDivJpar[dit];
      
      div_J_electron_mag[dit] *= a_dt;
      gkPoissonRHS[dit] -= div_J_electron_mag[dit];

      div_J_ion_mag[dit] *= a_dt;
      gkPoissonRHS[dit] -= div_J_ion_mag[dit];

      if (m_reynolds_stress) {
         negative_div_ExB_vorticity[dit] *= a_dt;
         gkPoissonRHS[dit] -= negative_div_ExB_vorticity[dit];
      }
   }

   // Add stablization terms
   if (m_include_stabilization) {
      LevelData<FArrayBox> negative_vorticity_times_dt(grids, 1, IntVect::Zero);
      m_gyropoisson_op->computeFluxDivergence(current_phi, negative_vorticity_times_dt, false);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         negative_vorticity_times_dt[dit].mult(a_dt);
      }
      addStabilizationTerms(gkPoissonRHS, negative_vorticity_times_dt);
   }
   
   // Do the implicit vorticity solve
   computeIonMassDensity(m_ion_mass_density, a_kinetic_species_phys, a_fluid_species_comp);
   m_vorticity_op->setImplicitDt(a_dt);
      
   Vector<Real> Er_scalar_vect;
   if ( a_cbc_or_ampere ) {
      int component = -1;
      for (int i=0; i<a_scalar_data.size(); ++i) {
         if (a_scalar_data[i]->name() == "Er_boundary") {
            component = i;
            break;
         }
      }
      CH_assert(component >= 0);
      Er_scalar_vect = a_scalar_data[component]->data();
   }

   if (Er_scalar_vect.size() > 0) setCoreBC( Er_scalar_vect[0], -Er_scalar_vect[1], *m_vorticity_op_bcs );

   bool update_precond = (m_it_counter % m_update_pc_freq == 0) ? true : false;
   m_vorticity_op->setVorticityOperatorCoefficients(m_ion_mass_density,
                                                    m_ion_charge_density,
                                                    m_electron_temperature,
                                                    *m_vorticity_op_bcs,
                                                    update_precond);

   LevelData<FArrayBox>& new_phi = a_EM_fields.getPhi();
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      new_phi[dit].copy(current_phi[dit]);
   }

   m_vorticity_op->solveWithBCs( new_phi, gkPoissonRHS );
   
   if (m_split_op) {
      LevelData<FArrayBox> negativeVortOld(grids, 1, IntVect::Zero);
      m_gyropoisson_op->computeFluxDivergence(new_phi, negativeVortOld, false);
   
      //Advance the high-order operator to update negative vorticity
      LevelData<FArrayBox> diffusionRHS(grids, 1, IntVect::Zero);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         diffusionRHS[dit].copy(negativeVortOld[dit]);
         diffusionRHS[dit].divide(a_dt);
      }
   
      m_vorticity_diffusion_op->setImplicitDt(a_dt);
      m_vorticity_diffusion_op->setBeta(1.0/a_dt);
      m_vorticity_diffusion_op->setVorticityOperatorCoefficients(zero,
                                                                 m_ion_charge_density,
                                                                 m_electron_temperature,
                                                                 *m_vorticity_bcs,
                                                                 update_precond);
   
      LevelData<FArrayBox> negativeVortNew(grids, 1, IntVect::Zero);
      m_vorticity_diffusion_op->solveWithBCs( negativeVortNew, diffusionRHS );
   
      //Obtain potential from the vorticity at the end of the step
      m_gyropoisson_op->setOperatorCoefficients(m_ion_mass_density, *m_gyropoisson_op_bcs, update_precond);
      m_gyropoisson_op->solveWithBCs( new_phi, negativeVortNew );
   }
   
   m_vorticity_op->fillInternalGhosts( new_phi );

   // Update nodal phi if supporting the calculation of a divergence-free phase velocity
   a_phi_ops.interpToNodes(new_phi, a_EM_fields.getPhiNode());

   // Set up the polarization density operator for the next step
   if (Er_scalar_vect.size() > 0) setCoreBC( Er_scalar_vect[0], -Er_scalar_vect[1], *m_gyropoisson_op_bcs );
   m_gyropoisson_op->setOperatorCoefficients(m_ion_mass_density, *m_gyropoisson_op_bcs, false);
   
   // Optional phi post-processing to remove the axisymmetric component
   if (m_remove_axisymmetric_phi) {
     FluxSurface flux_surface(m_geometry);
     LevelData<FArrayBox> phi_fs(grids, 1, IntVect::Zero);
     flux_surface.averageAndSpread(new_phi, phi_fs);

     for (DataIterator dit(grids); dit.ok(); ++dit) {
       new_phi[dit] -= phi_fs[dit];
     }
   }

   m_it_counter++;
}


void VorticityBEOp::parseParameters( const ParmParse& a_pp )
{
   a_pp.query( "update_precond_interval", m_update_pc_freq);
   a_pp.query( "use_operator_split", m_split_op);
}


void VorticityBEOp::printParameters()
{

}

#include "NamespaceFooter.H"
