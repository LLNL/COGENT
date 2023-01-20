#include "EMVorticityBEOp.H"
#include "EllipticOpBCFactory.H"

#include "NamespaceHeader.H" 

const char* EMVorticityBEOp::pp_name = {"EMvorticityBE_op"};

EMVorticityBEOp::EMVorticityBEOp(const ParmParse&    a_pp,
                                 const MagGeom&      a_geometry,
                                 const double        a_larmor,
                                 const double        a_ion_skin_depth,
                                 const double        a_electron_skin_depth,
                                 const int           a_verbosity )
   : EMVorticityOp(a_pp,
                   a_geometry,
                   a_larmor,
                   a_ion_skin_depth,
                   a_electron_skin_depth,
                   a_verbosity),
     m_dt_sub_cycle(DBL_MAX),
     m_update_pc_freq(1)

{
   // Parse parameters
   parseParameters(a_pp);
   
   // Print parameters
   if (m_verbosity) {
     printParameters();
   }
   
   
   // Set up GKPoisson and Asoln_op operator;
   // N.B. we cannot use the base class members
   // because we pass a different ParmParse object
 
   // Set up gyrokinetic Poisson operator
   string s_pp_base = string(a_pp.prefix()) + ".phi_soln_op";
   ParmParse pp_base = s_pp_base.c_str();

   m_gyropoisson_op = new GKPoisson(pp_base,
                                    a_geometry,
                                    a_larmor,
                                    0.,
                                    true,
                                    m_second_order,
                                    m_low_pollution );
   
   m_gyropoisson_op_bcs = m_potential_bcs->clone();
   
   // Set up perpendicular Lplacian operator
   s_pp_base = string(a_pp.prefix()) + ".Apar_soln_op";
   pp_base = s_pp_base.c_str();

   m_Apar_soln_op = new GKEMVorticity(a_pp,
                                      pp_base,
                                      a_geometry,
                                      a_larmor,
                                      m_ion_skin_depth,
                                      m_electron_skin_depth,
                                      m_second_order,
                                      m_low_pollution,
                                      m_include_pol_den_correction,
                                      m_include_diffusion,
                                      "AparSolnOp");
   
   m_Apar_soln_op_bcs = m_A_parallel_bcs->clone();
 
   // Create m_Apar; presently it is internal to the model
   // Should be taken outise this class, when we include Apar
   // corrections to the GK velocity
   const DisjointBoxLayout& grids = m_geometry.grids();
   m_Apar.define(grids, 1, IntVect::Zero);
   
   //TEMPORARY
   // Becasue m_Apar is currently not visible from outside this class
   // we have to initialize it here;
   // need to fix it for restart purposes though!!!
   setZero(m_Apar);
}


EMVorticityBEOp::~EMVorticityBEOp()
{
   delete m_gyropoisson_op;
   delete m_Apar_soln_op;
}


void EMVorticityBEOp::initializeBEModel(const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                        const FluidSpeciesPtrVect&        a_fluid_species_phys,
                                        const Vector<Real>&               a_scalar_data )
{
   
}

void EMVorticityBEOp::updateBEModelData(const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                        const FluidSpeciesPtrVect&        a_fluid_species_comp,
                                        const Vector<Real>&               a_scalar_data )
{
   if (a_scalar_data.size() > 0) {
      setCoreBC( a_scalar_data[0], -a_scalar_data[1], *m_gyropoisson_op_bcs );
   }
   
   m_gyropoisson_op->setOperatorCoefficients(m_ion_mass_density,
                                             *m_gyropoisson_op_bcs,
                                             true);
   
   
   m_Apar_soln_op->setVorticityOperatorCoefficients(m_ion_mass_density,
                                                    m_ion_charge_density,
                                                    m_electron_temperature,
                                                    *m_Apar_soln_op_bcs,
                                                    true);
   
}


void EMVorticityBEOp::updatePotentialBEModel(EMFields&                         a_EM_fields,
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
 
   PS::ScalarPtrVect dummy_scalars;
   preOpEval(a_kinetic_species_phys, a_fluid_species_comp, dummy_scalars, a_EM_fields, a_time);
   
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

   updateBEModelData(a_kinetic_species_phys, a_fluid_species_comp, Er_scalar_vect);
   
   const DisjointBoxLayout& grids = m_geometry.grids();
   LevelData<FArrayBox> rhs_phi(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> rhs_Apar(grids, 1, IntVect::Zero);

   LevelData<FArrayBox> phi_half(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> Apar_half(grids, 1, IntVect::Zero);

   const LevelData<FArrayBox>& current_phi = a_EM_fields.getPhi();
   LevelData<FArrayBox> new_phi;
   new_phi.define(current_phi);
   
   // Do sub-cycling time integration with RK2
   Real time = a_time;
   while (time < a_time + a_dt) {

      Real dt;

      if (time + m_dt_sub_cycle < a_time + a_dt) {
         dt = m_dt_sub_cycle;
      }
      else {
         dt = a_time + a_dt - time;
      }
      
      // Get RHS for phi and Apar half step update
      setZero(rhs_phi);
      setZero(rhs_Apar);
      EMVorticityOp::accumulateRHS(rhs_phi,
                                   rhs_Apar,
                                   new_phi,
                                   m_Apar,
                                   a_kinetic_species_phys,
                                   a_EM_fields,
                                   time);
      
      // Get phi and Apar half step values
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         phi_half[dit].copy(new_phi[dit]);
         Apar_half[dit].copy(m_Apar[dit]);
      }
      Real half_dt = 0.5*dt;
      advance(phi_half, Apar_half, rhs_phi, rhs_Apar, half_dt);
      
      // Get RHS for phi and Apar full step update
      setZero(rhs_phi);
      setZero(rhs_Apar);
      Real time_half = time + half_dt;
      EMVorticityOp::accumulateRHS(rhs_phi,
                                   rhs_Apar,
                                   phi_half,
                                   Apar_half,
                                   a_kinetic_species_phys,
                                   a_EM_fields,
                                   time_half);

      // Get phi and Apar end step values
      advance(new_phi, m_Apar, rhs_phi, rhs_Apar, dt);

      time = time + dt;
   }

   // Now that have the new phi, update it in a_EM_fields and compute the corresponding EField

   a_EM_fields.setPhi(new_phi);
   a_phi_ops.updateEField( a_EM_fields,
                           false,
                           a_kinetic_species_phys,
                           a_fluid_species_comp,
                           a_scalar_data,
                           a_potential_bc,
                           true,
                           false );

   m_it_counter++;
}


void EMVorticityBEOp::advance(LevelData<FArrayBox>&   a_phi,
                              LevelData<FArrayBox>&   a_Apar,
                              LevelData<FArrayBox>&   a_rhs_phi,
                              LevelData<FArrayBox>&   a_rhs_Apar,
                              const Real              a_dt)
{
   const DisjointBoxLayout& grids = m_geometry.grids();
   
   // Multiply by dt
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_rhs_phi[dit].mult(a_dt);
      a_rhs_Apar[dit].mult(a_dt);
   }
   
   // Add M_phi[phi_old]=-neg_vort_old to the RHS for phi
   LevelData<FArrayBox> soln_op_old(grids, 1, IntVect::Zero);
   m_gyropoisson_op->computeFluxDivergence(a_phi, soln_op_old, false);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_rhs_phi[dit].plus(soln_op_old[dit]);
    }

   // Add M_Apar[Apar_old] to the RHS for Apar
   m_Apar_soln_op->applyOp(soln_op_old, a_Apar, false);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_rhs_Apar[dit].plus(soln_op_old[dit]);
    }

   // Solve phi and Apar
   m_gyropoisson_op->solveWithBCs( a_phi, a_rhs_phi );
   m_Apar_soln_op->solveWithBCs( a_Apar, a_rhs_Apar );
}

void EMVorticityBEOp::parseParameters( const ParmParse& a_pp )
{
   a_pp.query( "update_precond_interval", m_update_pc_freq);
   a_pp.query( "dt_subcycle", m_dt_sub_cycle);
}


void EMVorticityBEOp::printParameters()
{

}

#include "NamespaceFooter.H"
