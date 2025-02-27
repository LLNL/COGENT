#include "EMVorticityOp.H"
#include "Directions.H"
#include "SimpleDivergence.H"

#include "EllipticOpBCFactory.H"
#include "SpaceUtils.H.multidim"

#undef VERIFY_MATRIX

#undef ALT_PAR_GRAD_PHI
#define ALT_DIV_JPAR

#include "NamespaceHeader.H" 

const char* EMVorticityOp::pp_name = {"EMvorticity_op"};

EMVorticityOp::EMVorticityOp(const ParmParse&    a_pp,
                             const MagGeom&      a_geometry,
                             const double        a_larmor,
                             const double        a_ion_skin_depth,
                             const double        a_electron_skin_depth,
                             const int           a_verbosity )
   : VorticityOp(a_pp, a_geometry, a_larmor, a_verbosity),
     m_geometry(a_geometry),
     m_is_time_implicit(false),
     m_ion_skin_depth(a_ion_skin_depth),
     m_electron_skin_depth(a_electron_skin_depth)
{
      
   // Input the initial conditions
   parseInitialConditions("", "potential");
   parseInitialConditions("", "A_parallel");

   // Parse remaining parameters 
   parseParameters(a_pp);
   if (m_verbosity) {
     printParameters();
   }
       
   // Create Apar BCs
   const std::string name("A_parallel");
   const std::string prefix( "BC." + name );
   ParmParse ppsp( prefix.c_str() );
   EllipticOpBCFactory elliptic_op_bc_factory;

   m_A_parallel_bcs = elliptic_op_bc_factory.create(name,
                                                    ppsp,
                                                    *(a_geometry.getCoordSys()),
                                                    false);
    
   
   // Set up gyrokinetic Poisson operator
   string s_pp_base = string(a_pp.prefix()) + ".em_soln_op";
   ParmParse pp_base = s_pp_base.c_str();

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
   
   // Set up ImEx preconditioner operator
   s_pp_base = string(a_pp.prefix()) + ".pc_imex";
   pp_base =s_pp_base.c_str();
   string pc_model;
   pp_base.query("pc_model", pc_model);
   m_imex_pc_op = new GKEMVorticity(a_pp,
                                    pp_base,
                                    a_geometry,
                                    a_larmor,
                                    m_ion_skin_depth,
                                    m_electron_skin_depth,
                                    m_second_order,
                                    m_low_pollution,
                                    m_include_pol_den_correction,
                                    m_include_diffusion,
                                    pc_model);
   
   // Get jparallel BCs
   const std::string name_high_order("jparallel");
   const std::string prefix_high_order( "BC." + name_high_order );
   ParmParse pp_high_order( prefix_high_order.c_str() );
   m_jpar_bcs = RefCountedPtr<EllipticOpBC>(elliptic_op_bc_factory.create(name_high_order,
                                                                          pp_high_order,
                                                                          *(a_geometry.getCoordSys()),
                                                                          false));
   m_imex_pc_op->setJparallelBc(*m_jpar_bcs);
   
   // Set up negative perpendicular Laplacian operator for initialization and ALT_DIV_J option
   const DisjointBoxLayout& grids = m_geometry.grids();
   
   s_pp_base = string(a_pp.prefix()) + ".nlp_Apar_op";
   pp_base =s_pp_base.c_str();
   m_nlp_Apar_op = new Diffusion(pp_base, m_geometry);
   
   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   LevelData<FluxBox> D_tensor(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   LevelData<FluxBox> D_tensor_mapped(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   double fac = pow(m_ion_skin_depth, 2);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      D_tensor[dit].copy(perp_coeff[dit]);
      D_tensor_mapped[dit].copy(perp_coeff_mapped[dit]);
      D_tensor[dit] *= fac;
      D_tensor_mapped[dit] *= fac;
   }
   m_nlp_Apar_op->setOperatorCoefficients(D_tensor, D_tensor_mapped, *m_A_parallel_bcs);
   m_nlp_Apar_op->updatePreconditioner(*m_A_parallel_bcs);
   
   
   // Set up div(b...) operator for initialization and ALT_DIV_J option
   s_pp_base = string(a_pp.prefix()) + ".div_jpar_op";
   pp_base =s_pp_base.c_str();
   m_div_jpar_op = new Advection(pp_base, m_geometry);

   const LevelData<FluxBox>& bunit = m_geometry.getFCBFieldDir();
   m_div_jpar_op->setOperatorCoefficients(bunit);
   m_div_jpar_op->updatePreconditioner(*m_jpar_bcs);
   
   
   // Set up vorticity diffusion; we cannot inherit this class from VorticityOp
   // because the base class member could have polarization density corrections
   // this class handles polarization density corrections separately from this diffusion op
   if (m_include_diffusion) {
       
      s_pp_base = string(a_pp.prefix()) + ".diffusion_op";
      pp_base = s_pp_base.c_str();
       
      m_vorticity_diffusion_op = new GKVorticity(a_pp,
                                                 pp_base,
                                                 a_geometry,
                                                 a_larmor,
                                                 0.,
                                                 m_second_order,
                                                 m_low_pollution,
                                                 false,
                                                 m_include_diffusion,
                                                 "VorticityDiffusion" );

      m_vorticity_diffusion_op->setDiffusionOperatorCoeff(m_pol_diffusion_face,
                                                          m_perp_diffusion_face);
   }
   else {
      m_vorticity_diffusion_op = NULL;
   }
   
   m_par_cond.define(grids, 1, IntVect::Zero);
   
   m_volume.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume);
}


EMVorticityOp::~EMVorticityOp()
{
   if (m_imex_pc_op) delete m_imex_pc_op;
   if (m_Apar_soln_op) delete m_Apar_soln_op;
   if (m_nlp_Apar_op) delete m_nlp_Apar_op;
   if (m_div_jpar_op) delete m_div_jpar_op;
   if (m_vorticity_diffusion_op) delete m_vorticity_diffusion_op;
   if (m_A_parallel_bcs) delete m_A_parallel_bcs;
 }


void EMVorticityOp::accumulateExplicitRHS(FluidSpeciesPtrVect&               a_rhs,
                                          const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                          const FluidSpeciesPtrVect&         a_fluid_species,
                                          const PS::ScalarPtrVect&           a_scalars,
                                          const EMFields&                    a_EM_fields,
                                          const int                          a_fluid_vec_comp,
                                          const Real                         a_time )
{
  
   // Get soln data
   const CFGVars& sol_fluid_species = *(a_fluid_species[a_fluid_vec_comp]);
   const LevelData<FArrayBox>& soln_phi = sol_fluid_species.cell_var("potential");
   const LevelData<FArrayBox>& soln_Apar = sol_fluid_species.cell_var("A_parallel");

   // Get RHS data
   CFGVars& rhs_fluid_species = *(a_rhs[a_fluid_vec_comp]);
   LevelData<FArrayBox>& rhs_phi = rhs_fluid_species.cell_var("potential");
   LevelData<FArrayBox>& rhs_Apar = rhs_fluid_species.cell_var("A_parallel");
  
   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   
#ifndef ALFVEN_WAVE_TEST
   // Add divergence of the perpendicular ion current
   // N.B. if we want to experiment with placing divJperp,i into implicitOpImEx
   // we need to make sure that kinetic ghosts are filled at chkpt_stage_func_0
   // default optimization implementation does not fill kinetic ghosts at chkpt_stage_func_0
   computeDivPerpIonMagCurrentDensity(m_divJperp_mag_i, a_EM_fields, a_kinetic_species_phys, a_time);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
       rhs_phi[dit] -= m_divJperp_mag_i[dit];
   }
#endif
   
   // Add divergence of the perpendicular electron current
   CH_assert(!m_include_pol_den_correction_to_pe);
   LevelData<FArrayBox> zero(grids, 1, IntVect::Zero);
   setZero(zero);
   // N.B. m_ion_charge used ensity here comes from the previous call to preOpEval.
   // For the optimzied time integration that call corresponds to chkpt_stage_func_0,
   // which is fine for the case where the evolution of dfn is explicit and thus density
   // does not change within the implicit stage solve.
   computeDivPerpElectronMagCurrentDensity(m_divJperp_mag_e, m_ion_charge_density, zero, a_time);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      rhs_phi[dit] -= m_divJperp_mag_e[dit];
   }
   
   //Subtract initial value of the negative divJperp
   if (m_apply_initialization_constraints) {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         rhs_phi[dit].minus(m_negative_divJperp_0[dit]);
      }
   }
   
   if (!m_is_time_implicit) {
      accumulateRHS(rhs_phi, rhs_Apar, soln_phi, soln_Apar, a_kinetic_species_phys, a_EM_fields, a_time);
   }
}


void EMVorticityOp::accumulateImplicitRHS(FluidSpeciesPtrVect&               a_rhs,
                                          const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                          const FluidSpeciesPtrVect&         a_fluid_species,
                                          const PS::ScalarPtrVect&           a_scalars,
                                          const EMFields&                    a_EM_fields,
                                          const int                          a_fluid_vec_comp,
                                          const Real                         a_time )
{
   if (m_is_time_implicit) {
      // Get soln data
      const CFGVars& sol_fluid_species = *(a_fluid_species[a_fluid_vec_comp]);
      const LevelData<FArrayBox>& soln_phi = sol_fluid_species.cell_var("potential");
      const LevelData<FArrayBox>& soln_Apar = sol_fluid_species.cell_var("A_parallel");

      // Get RHS data
      CFGVars& rhs_fluid_species = *(a_rhs[a_fluid_vec_comp]);
      LevelData<FArrayBox>& rhs_phi = rhs_fluid_species.cell_var("potential");
      LevelData<FArrayBox>& rhs_Apar = rhs_fluid_species.cell_var("A_parallel");
     
      accumulateRHS(rhs_phi, rhs_Apar, soln_phi, soln_Apar, a_kinetic_species_phys, a_EM_fields, a_time);
   }
   
}


void EMVorticityOp::accumulateRHS(LevelData<FArrayBox>&              a_rhs_phi,
                                  LevelData<FArrayBox>&              a_rhs_Apar,
                                  const LevelData<FArrayBox>&        a_soln_phi,
                                  const LevelData<FArrayBox>&        a_soln_Apar,
                                  const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                  const EMFields&                    a_EM_fields,
                                  const Real                         a_time )
{
   accumulatePhiRHS(a_rhs_phi, a_soln_phi, a_soln_Apar, a_kinetic_species_phys, a_EM_fields, a_time);
   accumulateAparRHS(a_rhs_Apar, a_soln_phi, a_soln_Apar, a_kinetic_species_phys, a_EM_fields, a_time);
}


void EMVorticityOp::accumulatePhiRHS(LevelData<FArrayBox>&              a_rhs_phi,
                                     const LevelData<FArrayBox>&        a_soln_phi,
                                     const LevelData<FArrayBox>&        a_soln_Apar,
                                     const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                     const EMFields&                    a_EM_fields,
                                     const Real                         a_time )
{
   /*
    Compute phi RHS
   */

   // Get boxes
   const DisjointBoxLayout& grids( a_soln_phi.getBoxes() );

   // Compute jpar
   LevelData<FArrayBox> j_par(grids, 1, IntVect::Zero);

   //((EllipticOp*)m_Apar_soln_op)->computeFluxDivergence(a_soln_Apar, j_par, false);
   m_Apar_soln_op->EllipticOp::computeFluxDivergence(a_soln_Apar, j_par, false);
   
   double fac = pow(m_ion_skin_depth, 2);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      j_par[dit].mult(fac);
   }

   
   // Add divergence of the parallel current
   addDivJpar(a_rhs_phi, a_soln_Apar, j_par);

   // Add Reynolds Stress term
   if (m_reynolds_stress) {
      addReynoldsStressTerm(a_rhs_phi, a_soln_phi, a_EM_fields);
   }
   
   // Add vorticity diffusion
   if (m_include_diffusion) {

      // Compute negative vorticity
      LevelData<FArrayBox> negative_vorticity(grids, 1, IntVect::Zero);
      m_gyropoisson_op->computeFluxDivergence(a_soln_phi, negative_vorticity, false);
       
      // Compute -perp_grad * D * perp_grad of negative_vorticity
      LevelData<FArrayBox> pol_dens_correction(grids, 1, IntVect::Zero);
      m_vorticity_diffusion_op->computeFluxDivergence(negative_vorticity, pol_dens_correction, false, !m_use_vorticity_bcs);
       
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_rhs_phi[dit] -= pol_dens_correction[dit];
      }
   }

   // Add diamagnetic corrections
   if (m_include_diam_correction) {
      addDiamagneticCorrection(a_rhs_phi, a_soln_phi, a_EM_fields);
   }
   
   // Add stabilization terms
   if (m_include_stabilization) {
      
      LevelData<FArrayBox> negative_vorticity(grids, 1, IntVect::Zero);
      m_gyropoisson_op->computeFluxDivergence(a_soln_phi, negative_vorticity, false);

      addStabilizationTerms(a_rhs_phi, negative_vorticity);
   }
   
}


void EMVorticityOp::accumulateAparRHS(LevelData<FArrayBox>&              a_rhs_Apar,
                                      const LevelData<FArrayBox>&        a_soln_phi,
                                      const LevelData<FArrayBox>&        a_soln_Apar,
                                      const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                      const EMFields&                    a_EM_fields,
                                      const Real                         a_time )
{
   
   // Get boxes
   const DisjointBoxLayout& grids( a_soln_phi.getBoxes() );

   /*
    Compute A_par RHS
   */
   
   // Compute -nabla_phi with BCs
   LevelData<FArrayBox> neg_grad_phi(grids, SpaceDim, IntVect::Unit);
   m_gyropoisson_op->compute3DFieldWithBCs(a_soln_phi, neg_grad_phi, false);
   
   // Compute gradPe/ne term; begin with getting electron density
   LevelData<FArrayBox> electron_pressure(grids, 1, IntVect::Zero);
   for (DataIterator dit(electron_pressure.dataIterator()); dit.ok(); ++dit) {
      electron_pressure[dit].copy(m_ion_charge_density[dit]);
   }

   // Add polarization density corrections
   if (m_include_pol_den_correction) {
      
      // Compute negative_vorticity with phi BCs
      LevelData<FArrayBox> neg_vort(grids, 1, IntVect::Zero);
      m_gyropoisson_op->computeFluxDivergence(a_soln_phi, neg_vort, false);
      
      // Add polarization density correction
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         electron_pressure[dit].minus(neg_vort[dit]);
      }
   }

   // Compute corrected Pe
   for (DataIterator dit(electron_pressure.dataIterator()); dit.ok(); ++dit) {
      electron_pressure[dit].mult(m_electron_temperature[dit]);
   }
   
   // Compute -gradPe using extrapolation (since we don't have BCs for
   // Te * vorticity quantity anyways )
   LevelData<FArrayBox> grad_pe_over_ne(grids, SpaceDim, IntVect::Unit);
   m_gyropoisson_op->computeField(electron_pressure, grad_pe_over_ne);
   
   // Compute gradPe/ne
   for (DataIterator dit(electron_pressure.dataIterator()); dit.ok(); ++dit) {
      for (int n = 0; n < SpaceDim; n++) {
	grad_pe_over_ne[dit].divide(m_ion_charge_density[dit],0,n,1);
      }
      grad_pe_over_ne[dit].mult(-1.0);
   }
   
   // Compute 0.71*gradTe
   LevelData<FArrayBox> grad_Te_term(grids, SpaceDim, IntVect::Unit);
   m_gyropoisson_op->computeField(m_electron_temperature, grad_Te_term);
   for (DataIterator dit(electron_pressure.dataIterator()); dit.ok(); ++dit) {
      grad_Te_term[dit].mult(-0.71);
   }
   
   // Compute all parallel gradient terms
   LevelData<FArrayBox> par_grad_terms(grids, 1, IntVect::Zero);
#ifdef ALT_PAR_GRAD_PHI
   m_imex_pc_op->computeAltParGrad(a_soln_phi, par_grad_terms);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      par_grad_terms[dit].negate();
   }

#ifndef ALFVEN_WAVE_TEST
   LevelData<FArrayBox> sum_grad_terms(grids, SpaceDim, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      sum_grad_terms[dit].copy(grad_pe_over_ne[dit]);
      sum_grad_terms[dit].plus(grad_Te_term[dit]);
   }
   LevelData<FArrayBox> par_grad_terms2(grids, 1, IntVect::Zero);
   m_geometry.computeParallelProjection(par_grad_terms2, sum_grad_terms);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      par_grad_terms[dit] += par_grad_terms2;
   }
#endif

#else
   LevelData<FArrayBox> sum_grad_terms(grids, SpaceDim, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      sum_grad_terms[dit].copy(neg_grad_phi[dit]);
#ifndef ALFVEN_WAVE_TEST
      sum_grad_terms[dit].plus(grad_pe_over_ne[dit]);
      sum_grad_terms[dit].plus(grad_Te_term[dit]);
#endif
   }
   m_geometry.computeParallelProjection(par_grad_terms, sum_grad_terms);
#endif
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_rhs_Apar[dit].plus(par_grad_terms[dit]);
   }

   // Add -j_par/sigma, first compute j_par
   LevelData<FArrayBox> j_par(grids, 1, IntVect::Zero);

   //((EllipticOp*)m_Apar_soln_op)->computeFluxDivergence(a_soln_Apar, j_par, false);
   m_Apar_soln_op->EllipticOp::computeFluxDivergence(a_soln_Apar, j_par, false);
   
   double fac = pow(m_ion_skin_depth, 2);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      j_par[dit].mult(fac);
   }

#ifndef ALFVEN_WAVE_TEST
   LevelData<FArrayBox> j_par_over_sigma(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      j_par_over_sigma[dit].copy(j_par[dit]);
      j_par_over_sigma[dit].divide(m_par_cond[dit]);
      a_rhs_Apar[dit].minus(j_par_over_sigma[dit]);
   }
#endif
   
   // Add stabilization terms
   if (m_include_stabilization) {
      addStabilizationTerms(a_rhs_Apar, a_soln_Apar);
   }
}


void EMVorticityOp::addDivJpar(LevelData<FArrayBox>&        a_rhs,
                               const LevelData<FArrayBox>&  a_Apar,
                               const LevelData<FArrayBox>&  a_jpar)
{

   const DisjointBoxLayout& grids = m_geometry.grids();

   bool fourth_order = (m_geometry.secondOrder()) ? false : true;

   int ghosts = (fourth_order) ? 2 : 1;
   IntVect ghostVect = ghosts*IntVect::Unit;

   // Extrapolate parallel current to physical ghosts
   LevelData<FArrayBox> jpar_grown(grids, 1, ghostVect);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      jpar_grown[dit].copy(a_jpar[dit]);
   }
   m_geometry.extrapolateToPhysicalGhosts(jpar_grown, fourth_order);
   
   // Fill internal ghosts
   m_geometry.fillInternalGhosts(jpar_grown);
   
   // Get parallel velocity on faces
   LevelData<FluxBox> jpar_on_faces( grids, 1, (ghosts-1)*IntVect::Unit );
  
   //Put this in a geometry or space utils class
   int order = (fourth_order) ? 4 : 2;
   SpaceUtils::interpToFaces(jpar_on_faces, jpar_grown, order);

   // Force zero flux on radial and poloidal faces for a field-aligned grid
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      if (block_coord_sys.isFieldAligned()) {
         jpar_on_faces[dit][RADIAL_DIR].setVal(0.0);
         jpar_on_faces[dit][POLOIDAL_DIR].setVal(0.0);
      }
   }
   // Get unit b vector
   const LevelData<FluxBox>& b_dir = m_geometry.getFCBFieldDir();
   
   // Compute b*j_par (physical flux)
   LevelData<FluxBox> flux( grids, SpaceDim, (ghosts-1)*IntVect::Unit);
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      flux[dit].copy(b_dir[dit]);
      for (int dir = 0; dir < SpaceDim; dir++) {
         for (int n = 0; n < SpaceDim; n++) {
           flux[dit][dir].mult(jpar_on_faces[dit][dir],0,n,1);
         }
      }
   }
   
   // Compute integrated normal flux (b*jpar)
   LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
   m_geometry.computeMetricTermProductAverage(NTF_normal, flux, fourth_order);
   
   // Compute div(b*jpar)
   LevelData<FArrayBox> div_jpar(grids, 1, IntVect::Zero);

#ifndef ALT_DIV_JPAR
   RealVect fakeDx = RealVect::Unit;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      simpleDivergence(div_jpar[dit], NTF_normal[dit], grids[dit], fakeDx);
      div_jpar[dit] /= m_volume[dit];
   }
#else
   LevelData<FArrayBox> nlpApar(grids, 1, IntVect::Zero);
   m_nlp_Apar_op->applyPCOp(a_Apar, nlpApar);

   m_div_jpar_op->applyPCOp(nlpApar, div_jpar);
#endif

#ifndef ALFVEN_WAVE_TEST
   if (m_geometry.shearedMBGeom()) {

      //Subtract FA averaged parallel flux divergence on closed field-lines;
      //N.B. this function takes mapped flux as input, without mapped_face_area factor
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         RealVect faceArea = block_coord_sys.getMappedFaceArea();
         for (int dir=0; dir<SpaceDim; ++dir) {
            NTF_normal[dit][dir].divide(faceArea[dir]);
         }
      }
      m_parallel_current_divergence_op->subtractFSAverParDiv(div_jpar, NTF_normal);
   }
#endif

   // Add -div(b*jpar) to phi rhs
   for (DataIterator dit( div_jpar.dataIterator() ); dit.ok(); ++dit) {
      a_rhs[dit].minus(div_jpar[dit]);
   }

}

void EMVorticityOp::setAparDerivatives( const PS::KineticSpeciesPtrVect&  a_kinetic_species_phys,
                                        const LevelData<FArrayBox>&       a_Apar,
                                        EMFields&                         a_EM_fields )
{
   const DisjointBoxLayout& grids = m_geometry.grids();
   
   // Compute the gradient of Apar

   LevelData<FArrayBox> grad_Apar(grids, 3, IntVect::Unit);
   m_gyropoisson_op->compute3DFieldWithBCs(a_Apar, grad_Apar, false);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      grad_Apar[dit].negate();  // Get rid of the minus sign included in the field calculation
   }

   // Compute the time derivative of Apar

   LevelData<FArrayBox> rhs(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      rhs[dit].setVal(0.);
   }
   accumulateAparRHS(rhs, a_EM_fields.getPhi(), a_Apar, a_kinetic_species_phys, a_EM_fields, 0.);

   LevelData<FArrayBox> dApardt(grids, 1, IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      dApardt[dit].setVal(0.);  // Initial guess
   }
   m_Apar_soln_op->solveWithBCs(dApardt, rhs);

   // Cell-centered Apar derivatives

   LevelData<FArrayBox>& Apar_derivs_cell = a_EM_fields.getAparDerivsCell();
   if ( !Apar_derivs_cell.isDefined() ) {
      Apar_derivs_cell.define(grids, 4, IntVect::Unit);
   }

   // Load up the derivatives into the components, gradient first then the time derivative
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      Apar_derivs_cell[dit].copy(grad_Apar[dit],0,0,3);
      Apar_derivs_cell[dit].copy(dApardt[dit],0,3,1);
   }

   // Face-centered Apar derivatives

   bool fourth_order = false;
   m_geometry.extrapolateToPhysicalGhosts(Apar_derivs_cell, fourth_order);

   LevelData<FluxBox>& Apar_derivs_face = a_EM_fields.getAparDerivsFace();
   if ( !Apar_derivs_face.isDefined() ) {
      Apar_derivs_face.define(grids, 4, IntVect::Unit);
   }

   int order = 2;
   SpaceUtils::interpToFaces(Apar_derivs_face, Apar_derivs_cell, order);
}


void EMVorticityOp::preSolutionOpEval(const PS::KineticSpeciesPtrVect&   a_kinetic_species,
                                      const FluidSpeciesPtrVect&         a_fluid_species,
                                      const PS::ScalarPtrVect&           a_scalars,
                                      const EMFields&                    a_EM_fields,
                                      const double                       a_time )
{
   CH_TIME("EMVorticityOp::preSolutionOpEval");
   
   VorticityOp::preSolutionOpEval(a_kinetic_species, a_fluid_species, a_scalars, a_EM_fields, a_time);
   
   computeIonChargeDensity( m_ion_charge_density, a_kinetic_species, a_fluid_species );
   
   m_Apar_soln_op->setVorticityOperatorCoefficients(m_ion_mass_density,
                                                    m_ion_charge_density,
                                                    m_electron_temperature,
                                                    *m_A_parallel_bcs,
                                                    true);
}

void EMVorticityOp::preOpEval(const PS::KineticSpeciesPtrVect&   a_kinetic_species,
                              const FluidSpeciesPtrVect&         a_fluid_species,
                              const PS::ScalarPtrVect&           a_scalars,
                              const EMFields&                    a_EM_fields,
                              const double                       a_time )
{

   CH_TIME("EMVorticityOp::preOpEval");

   preSolutionOpEval(a_kinetic_species, a_fluid_species, a_scalars, a_EM_fields, a_time);
      
   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   LevelData<FArrayBox> zero(grids, 1, IntVect::Zero);
   setZero(zero);

   // Update parallel current divergence coefficients,
   // do we actually need this?? -- perhasp not!
   m_parallel_current_divergence_op->setVorticityOperatorCoefficients(zero,
                                                                      m_ion_charge_density,
                                                                      m_electron_temperature,
                                                                      *m_parallel_current_divergence_op_bcs,
                                                                      false );

   // Get parallel conductivity
   m_parallel_current_divergence_op->updateParallelConductivity(m_electron_temperature, m_par_cond);

   // N.B. Diffusion operator only applies perpendicular diffision for
   // the EM case; parallel high-order terms are included elswehere
   if (m_include_diffusion) {
     m_vorticity_diffusion_op->setVorticityOperatorCoefficients(zero,
                                                                zero,
                                                                zero,
                                                                *m_vorticity_bcs,
                                                                false);
   }
   
}


void EMVorticityOp::evalSolutionOp(FluidSpeciesPtrVect&               a_rhs,
                                   const PS::KineticSpeciesPtrVect&   a_kinetic_species_comp,
                                   const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                   const FluidSpeciesPtrVect&         a_fluid_species_comp,
                                   const FluidSpeciesPtrVect&         a_fluid_species_phys,
                                   const PS::ScalarPtrVect&           a_scalars,
                                   const int                          a_component,
                                   const Real                         a_time )
{
   CH_TIME("EMVorticityOp::evalSolutionOp");

   // This will evaluate solution operator for potential
   VorticityOp::evalSolutionOp(a_rhs,
                               a_kinetic_species_comp,
                               a_kinetic_species_phys,
                               a_fluid_species_comp,
                               a_fluid_species_phys,
                               a_scalars,
                               a_component,
                               a_time);

   // Evalute solution operator for Aparallel
   // Don't need to update bndry data, becasue that was done earlier and
   // scalar BCs are not changing for A_parallel (in contrast to phi)
   //m_Apar_soln_op->updateBoundaryDataAddingCoefs(*m_A_parallel_bcs);
   
   const FluidSpecies& sol_species = static_cast<const FluidSpecies&>(*a_fluid_species_phys[a_component]);
   FluidSpecies& rhs_species = static_cast<FluidSpecies&>(*a_rhs[a_component]);

   const LevelData<FArrayBox>& sol = sol_species.cell_var("A_parallel");
   LevelData<FArrayBox>& rhs = rhs_species.cell_var("A_parallel");

   m_Apar_soln_op->applyOp(rhs, sol, false);
}

void EMVorticityOp::solveSolutionPC(FluidSpeciesPtrVect&              a_fluid_species_solution,
                                    const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                    const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                    const int                         a_component )
{

   VorticityOp::solveSolutionPC(a_fluid_species_solution,
                                a_kinetic_species_rhs,
                                a_fluid_species_rhs,
                                a_component);
   
   FluidSpecies& sol_species = static_cast<FluidSpecies&>(*a_fluid_species_solution[a_component]);
   const FluidSpecies& rhs_species = static_cast<const FluidSpecies&>(*a_fluid_species_rhs[a_component]);

   const LevelData<FArrayBox>& r = rhs_species.cell_var("A_parallel");
   LevelData<FArrayBox>& z = sol_species.cell_var("A_parallel");

   setZero(z);
   m_Apar_soln_op->solvePreconditioner(r, z);
   
}

void EMVorticityOp::updatePCImEx(const FluidSpeciesPtrVect&       a_fluid_species,
                                 const PS::KineticSpeciesPtrVect& a_kinetic_species,
                                 const double                     a_time,
                                 const int                        a_step,
                                 const int                        a_stage,
                                 const double                     a_shift,
                                 const int                        a_component,
                                 const std::string&               a_op_name )
{
   CH_TIME("EMVorticityOp::updatePCImEx");
   
   m_imex_pc_op->m_alpha = a_shift;
   m_imex_pc_op->setVorticityOperatorCoefficients(m_ion_mass_density,
                                                  m_ion_charge_density,
                                                  m_electron_temperature,
                                                  *m_potential_bcs,
                                                  *m_A_parallel_bcs,
                                                  true);
   
#ifdef VERIFY_MATRIX
   verifyMatrix(a_shift);
#endif

}


void EMVorticityOp::solvePCImEx(FluidSpeciesPtrVect&              a_fluid_species_solution,
                                const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                const std::string&                a_op_name,
                                const int                         a_component )
{
   CH_TIME("EMVorticityOp::solvePCImEx");
   
   // Get rhs data
   const FluidSpecies& rhs_species
        = static_cast<const FluidSpecies&>(*a_fluid_species_rhs[a_component]);
   const LevelData<FArrayBox>& r_phi = rhs_species.cell_var("potential");
   const LevelData<FArrayBox>& r_Apar = rhs_species.cell_var("A_parallel");

   const DisjointBoxLayout& grids = m_geometry.grids();
   LevelData<FArrayBox> r(grids, 2, IntVect::Zero);
   for (DataIterator dit( r.dataIterator() ); dit.ok(); ++dit) {
      r[dit].copy(r_phi[dit], 0, 0, 1);
      r[dit].copy(r_Apar[dit], 0, 1, 1);
   }
   
   // Solve preconditioner
   LevelData<FArrayBox> z(grids, 2, IntVect::Zero);
   setZero(z);
   m_imex_pc_op->solvePreconditioner(r, z);
   
   // Copy pc vector solution to phi and Apar components
   FluidSpecies& sol_species = static_cast<FluidSpecies&>(*a_fluid_species_solution[a_component]);
   LevelData<FArrayBox>& z_phi = sol_species.cell_var("potential");
   LevelData<FArrayBox>& z_Apar = sol_species.cell_var("A_parallel");

   for (DataIterator dit( r.dataIterator() ); dit.ok(); ++dit) {
      z_phi[dit].copy(z[dit], 0, 0, 1);
      z_Apar[dit].copy(z[dit], 1, 0, 1);
    }
}


void EMVorticityOp::fillGhostCells(FluidSpecies&  a_species_phys,
                                   const double   a_time )
{
}

void EMVorticityOp::postTimeEval(FluidSpecies&  a_species_comp,
                                 FluidSpecies&  a_species_phys,
                                 const Real     a_dt,
                                 const Real     a_time,
                                 const int      a_stage )
{
   VorticityOp::postTimeEval(a_species_comp, a_species_phys, a_dt, a_time, a_stage);
   
   LevelData<FArrayBox>& Apar_comp = a_species_comp.cell_var("A_parallel");
   LevelData<FArrayBox>& Apar_phys = a_species_phys.cell_var("A_parallel");
   
   for (int dir=0; dir<SpaceDim; dir++) {
      if (m_harmonic_filtering[dir] == 1) {
         SpaceUtils::applyHarmonicFiltering(Apar_comp, dir);
         SpaceUtils::applyHarmonicFiltering(Apar_phys, dir);
      }
   }
}


bool EMVorticityOp::isInitializationConstrained(const FluidSpecies& a_fluid_phys,
						const int           a_step)
{
   return m_apply_initialization_constraints;
}

void EMVorticityOp::applyInitializationConstraints(FluidSpeciesPtrVect&               a_fluid_comp,
                                                   FluidSpeciesPtrVect&               a_fluid_phys,
                                                   const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                                   const EMFields&                    a_EM_fields,
                                                   const int                          a_component,
                                                   const double                       a_time )
{
   const DisjointBoxLayout& grids( m_geometry.grids() );

   //Computes the initial value of negative divJperp
   //m_negative_divJperp_0, which is subtracted from phi RHS
   VorticityOp::applyInitializationConstraints(a_fluid_comp,
                                               a_fluid_phys,
                                               a_kinetic_species_phys,
                                               a_EM_fields,
                                               a_component,
                                               a_time);
   
#if 0
   
   // This is a development version to initialize with proper jpar
   
   // Create RHS for j_par equation
   PS::ScalarPtrVect  scalar_species_dummy; //zero lenght dummy
   
   VorticityOp::preOpEval(a_kinetic_species_phys,
                          a_fluid_comp,
                          scalar_species_dummy,
                          a_EM_fields,
                          a_time);

   LevelData<FArrayBox> neg_divJperp_mag(grids, 1, IntVect::Zero);
   for (DataIterator dit( neg_divJperp_mag.dataIterator() ); dit.ok(); ++dit) {
      neg_divJperp_mag[dit].copy(m_divJperp_mag_e[dit]);
      neg_divJperp_mag[dit].plus(m_divJperp_mag_i[dit]);
      neg_divJperp_mag[dit].mult(-1.0);
   }
   
   // Subtract FS average to satisfy solvability constraints
   FluxSurface flux_surface(m_geometry);
   LevelData<FArrayBox> jperp_fs(grids, 1, IntVect::Zero);
   flux_surface.averageAndSpread(neg_divJperp_mag, jperp_fs);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      neg_divJperp_mag[dit] -= jperp_fs[dit];
   }
      
   // Solve for j_par
   LevelData<FArrayBox> jpar(grids, 1, IntVect::Zero);
   m_div_jpar_op->solveWithBCs(jpar, neg_divJperp_mag);
      
   // Solve for A_par
   LevelData<FArrayBox> Apar(grids, 1, IntVect::Zero);
   m_nlp_Apar_op->solveWithBCs(Apar, jpar);

   // Revise the initial condition for Apar
   CFGVars& sol_fluid_comp = *(a_fluid_comp[a_component]);
   CFGVars& sol_fluid_phys = *(a_fluid_phys[a_component]);

   LevelData<FArrayBox>& soln_Apar_comp = sol_fluid_comp.cell_var("A_parallel");
   LevelData<FArrayBox>& soln_Apar_phys = sol_fluid_phys.cell_var("A_parallel");

   setZero(soln_Apar_comp);
   setZero(soln_Apar_phys);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      soln_Apar_comp[dit].copy(Apar[dit]);
      soln_Apar_phys[dit].copy(Apar[dit]);
   }
#endif
}

void EMVorticityOp::parseParameters( const ParmParse& a_pp )
{

   a_pp.query( "implicit", m_is_time_implicit);
   
}


void EMVorticityOp::printParameters()
{
   if (procID()==0) {

   }
}

void EMVorticityOp::verifyMatrix(const double   a_shift )
{
 
   EMFields* EM_fields_dummy = new EMFields;
   PS::KineticSpeciesPtrVect  kinetic_species_dummy;

   const DisjointBoxLayout& grids = m_geometry.grids();
   LevelData<FArrayBox> phi(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> Apar(grids, 1, IntVect::Zero);
   
   //set phi and Apar here
   for (DataIterator dit( phi.dataIterator() ); dit.ok(); ++dit) {
      phi[dit].setVal(-1.0);
      phi[dit].plus(m_ion_mass_density[dit]);
      Apar[dit].setVal(-1.0);
      Apar[dit].plus(m_ion_mass_density[dit]);
   }
   
   LevelData<FArrayBox> vorticity(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> jpar(grids, 1, IntVect::Zero);

   m_gyropoisson_op->computeFluxDivergence(phi, vorticity, true);
   m_Apar_soln_op->computeFluxDivergence(Apar, jpar, true);
   
   for (DataIterator dit( phi.dataIterator() ); dit.ok(); ++dit) {
      vorticity[dit].mult(-1.0);
      jpar[dit].mult(m_ion_skin_depth);
      jpar[dit].mult(m_ion_skin_depth);
   }
   
   
   // Compute alpha * M
   LevelData<FArrayBox> alpha_M_phi(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> alpha_M_Apar(grids, 1, IntVect::Zero);
   m_gyropoisson_op->computeFluxDivergence(phi, alpha_M_phi, true);
   m_Apar_soln_op->applyOp(alpha_M_Apar, Apar, true);
   
   for (DataIterator dit( phi.dataIterator() ); dit.ok(); ++dit) {
      alpha_M_phi[dit].mult(a_shift);
      alpha_M_Apar[dit].mult(a_shift);
   }

   // Compute F_implicit(phi, Apar)
   LevelData<FArrayBox> implicitOp_phi(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> implicitOp_Apar(grids, 1, IntVect::Zero);
   setZero(implicitOp_Apar);
   setZero(implicitOp_phi);
   
   accumulateRHS(implicitOp_phi, implicitOp_Apar, phi, Apar, kinetic_species_dummy, *EM_fields_dummy, 0.);

   // Compute action of the physical op alpha*M - F_i
   LevelData<FArrayBox> physOp_phi(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> physOp_Apar(grids, 1, IntVect::Zero);
   for (DataIterator dit( phi.dataIterator() ); dit.ok(); ++dit) {
      physOp_phi[dit].copy(alpha_M_phi[dit]);
      physOp_phi[dit].minus(implicitOp_phi[dit]);
      physOp_Apar[dit].copy(alpha_M_Apar[dit]);
      physOp_Apar[dit].minus(implicitOp_Apar[dit]);
   }
   
   // Now compute action of the preconditioner op
   LevelData<FArrayBox> z(grids, 4, IntVect::Zero);
   for (DataIterator dit( phi.dataIterator() ); dit.ok(); ++dit) {
      z[dit].copy(phi[dit], 0, 0, 1);
      z[dit].copy(vorticity[dit], 0, 1, 1);
      z[dit].copy(Apar[dit], 0, 2, 1);
      z[dit].copy(jpar[dit], 0, 3, 1);
   }
   
   LevelData<FArrayBox> r(grids, 4, IntVect::Zero);
   m_imex_pc_op->applyPCOp(z, r);
   
   LevelData<FArrayBox> pcOp_phi(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> pcOp_Apar(grids, 1, IntVect::Zero);
   for (DataIterator dit( phi.dataIterator() ); dit.ok(); ++dit) {
      pcOp_phi[dit].copy(r[dit],0,0,1);
      pcOp_Apar[dit].copy(r[dit],2,0,1);
   }
   
   m_geometry.plotCellData("phi", phi, 0.);
   
   m_geometry.plotCellData("pcOp_phi", pcOp_phi, 0.);
   m_geometry.plotCellData("pcOp_Apar", pcOp_Apar, 0.);
   m_geometry.plotCellData("physOp_phi", physOp_phi, 0.);
   m_geometry.plotCellData("physOp_Apar", physOp_Apar, 0.);
   
   delete EM_fields_dummy;
   exit(1);
}

#include "NamespaceFooter.H"
