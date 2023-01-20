#include "GKVorticityBE.H"
#include "GKPoissonF_F.H"
#include "FourthOrderUtil.H"
#include "SimpleDivergence.H"
#include "SpaceUtils.H.multidim"
#include "Directions.H"
#include "ToroidalBlockLevelExchangeCenter.H"

//#define VERIFY_MATRIX

#include "NamespaceHeader.H"
 
const char* GKVorticityBE::pp_name = {"GKVorticityBE"};

GKVorticityBE::GKVorticityBE(const ParmParse&   a_pp,
                             const ParmParse&   a_pp_base,
                             const MagGeom&     a_geom,
                             const Real         a_larmor_number,
                             const Real         a_debye_number,
                             const bool         a_second_order,
                             const bool         a_low_pollution,
                             const bool         a_include_pol_den_correction,
                             const bool         a_include_diffusion,
                             const std::string& a_model)
   : GKVorticity(a_pp,
                 a_pp_base,
                 a_geom,
                 a_larmor_number,
                 a_debye_number,
                 a_second_order,
                 a_low_pollution,
                 a_include_pol_den_correction,
                 a_include_diffusion,
                 a_model)

{
   parseParameters( a_pp );

   if (m_verbosity>0) {
      printParameters();
   }
}
      
GKVorticityBE::~GKVorticityBE()
{
}

void
GKVorticityBE::setVorticityOperatorCoefficients(const LevelData<FArrayBox>&  a_ion_mass_density,
                                                const LevelData<FArrayBox>&  a_ion_charge_density,
                                                const LevelData<FArrayBox>&  a_electron_temperature,
                                                EllipticOpBC&                a_bc,
                                                const bool                   a_update_preconditioner )
{
   
   /*
    GKVorticityBE class solves both the preconditioner and
    the physical operator problems. When the high-order terms are
    present and we use bcs for high-order vorticy variable,
    we need to maintain a local copy of low-order potential variable
    bcs (passed as the a_bc argument). More details are included in
    the comments to the applyHighOrderCorrectionOp() function
    */
   
   // Update tensor coefficeints, low-order phi-bc relevant data and a preconditioner
   GKVorticity::setVorticityOperatorCoefficients(a_ion_mass_density,
                                                 a_ion_charge_density,
                                                 a_electron_temperature,
                                                 a_bc,
                                                 a_update_preconditioner);
      
   // Save a local copy of potential BCs
   m_potential_bc = a_bc.clone();

   // Now compute inhomogeneous contribution
   // needed to solve the physical-operator problem
   computeBcDivergence(m_bc_divergence);
      
#ifdef VERIFY_MATRIX
   if (m_include_high_order_corr) {
      verifyMatrix(a_preconditioner);
   }
#endif
}


void
GKVorticityBE::verifyMatrix( const MBSolver*  a_matrix )
{
   CH_assert(m_include_high_order_corr);

   // Check the difference between the input block matrix and the
   // operator for which it's intended to server as a preconditioner.

   const DisjointBoxLayout& grids = m_geometry.grids();

   LevelData<FArrayBox> u(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> in_vec(grids, 2, IntVect::Zero);
   srand(time(NULL));
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (BoxIterator bit(grids[dit]); bit.ok(); ++bit) {
         IntVect iv = bit();

         u[dit](iv,0) = (double)rand() / (double)RAND_MAX;
         in_vec[dit](iv,0) = u[dit](iv,0);
      }
   }

   // Compute Mu
   LevelData<FArrayBox> Mu(grids, 1, IntVect::Zero);
   if ( !m_low_pollution ) {
      computeFluxDivergenceWithCoeff(Mu, u, m_M_unmapped, true, false);
   }
   else {
      computeFluxDivergenceWithCoeff(Mu, u, m_M_mapped, true, false);
   }
      
   if ( m_symmetrized_preconditioner ) {

      LevelData<FArrayBox> Dv(grids, 1, IntVect::Zero);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         Dv[dit].copy(Mu[dit]);
         Dv[dit] /= m_ne_over_Te[dit];
         Dv[dit].negate();

         // Set the vector (ubar,Dv) = (u - Dv, Dv)
         in_vec[dit].minus(Dv[dit],0,0,1);
         in_vec[dit].copy(Dv[dit],0,1,1);
      }
   }
   else {

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         Mu[dit].negate();
         in_vec[dit].copy(Mu[dit],0,1,1);
      }
   }

   LevelData<FArrayBox> matvec(grids, 2, IntVect::Zero);
   a_matrix->multiplyMatrix(in_vec, matvec);

   LevelData<FArrayBox> block_equation_0(grids, 1, IntVect::Zero);
   applyOp(block_equation_0, u, true);

   LevelData<FArrayBox> block_equation_1(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      block_equation_0[dit].minus(matvec[dit],0,0,1);
      block_equation_0[dit].abs();
      block_equation_1[dit].copy(matvec[dit],1,0,1);
      block_equation_1[dit].abs();
   }

   double norm0 = SpaceUtils::MaxNorm(block_equation_0);
   double norm1 = SpaceUtils::MaxNorm(block_equation_1);
   
   if (procID()==0) {
      cout << "Block equation 0 diff norm = " << norm0 
           << ", Block equation diff norm = " << norm1 << endl;
   }

   m_geometry.plotCellData("block_0_diff", block_equation_0, 0.);
   m_geometry.plotCellData("block_1_diff", block_equation_1, 0.);

   exit(1);
}


void
GKVorticityBE::computeFluxDivergenceWithCoeff(LevelData<FArrayBox>&       a_out,
                                              const LevelData<FArrayBox>& a_in,
                                              const LevelData<FluxBox>&   a_coeff,
                                              const bool                  a_homogeneous_bcs,
                                              const bool                  a_subtract_fs_par_div,
                                              const bool                  a_extrap_to_ghosts)
{
  /*
   Computes -\nabla(a_coeff * \nabla a_in)
   */
  CH_TIME("GKVorticityBE::computeFluxDivergenceWithCoeff");
   
  const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

  LevelData<FArrayBox> phi(grids, 1, 3*IntVect::Unit);

  for (DataIterator dit(grids); dit.ok(); ++dit) {
     phi[dit].copy(a_in[dit]);
  }

  LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);

  if ( !m_low_pollution ) {

     if (SpaceDim == 3) {
        if (a_extrap_to_ghosts) computeField(phi, flux);
        else compute3DFieldWithBCs(phi, flux, a_homogeneous_bcs);

     }
     else {
        if (a_extrap_to_ghosts) computePoloidalField(phi, flux);
        else computePoloidalFieldWithBCs(phi, flux, a_homogeneous_bcs);
     }

     // Multiply the field by the unmapped, face-centered GKP coefficients
     m_geometry.multiplyMatrix(flux, a_coeff);

     m_geometry.fillTransversePhysicalGhosts(flux);

     m_geometry.applyAxisymmetricCorrection(flux);

     // Convert to face-averaged
     if (!m_second_order) fourthOrderAverage(flux);

     m_geometry.computeMappedGridDivergence(flux, a_out, !m_second_order);

     for (DataIterator dit(grids); dit.ok(); ++dit) {
        a_out[dit] /= m_volume[dit];
     }
  
  }
  else {
  
     if (SpaceDim == 3) {
        if (a_extrap_to_ghosts) computeMapped3DField(phi, flux);
        else computeMapped3DFieldWithBCs(phi, flux, a_homogeneous_bcs);
     }

     else  {
        if (a_extrap_to_ghosts) computeMappedPoloidalField(phi, flux);
        else computeMappedPoloidalFieldWithBCs(phi, flux, a_homogeneous_bcs);
     }

     // Multiply the field by the mapped, face-centered GKP coefficients
     m_geometry.multiplyMatrix(flux, a_coeff);

     m_geometry.fillTransversePhysicalGhosts(flux);

     // Convert to face-averaged
     if (!m_second_order) fourthOrderAverage(flux);

     m_geometry.averageAtBlockBoundaries(flux);

     LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
     for (DataIterator dit(flux.dataIterator()); dit.ok(); ++dit) {
        const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
        RealVect faceArea = block_coord_sys.getMappedFaceArea();
        for (int dir=0; dir<SpaceDim; ++dir) {
           NTF_normal[dit][dir].copy(flux[dit][dir],dir,0,1);
           NTF_normal[dit][dir].mult(faceArea[dir]);
        }
     }

     RealVect fakeDx = RealVect::Unit;
     for (DataIterator dit(grids); dit.ok(); ++dit) {
        simpleDivergence(a_out[dit], NTF_normal[dit], grids[dit], fakeDx);
     }

     for (DataIterator dit(grids); dit.ok(); ++dit) {
        a_out[dit] /= m_volume[dit];
     }
  
     if (a_subtract_fs_par_div) {
        subtractFSAverParDiv(a_out, flux);
     }
  }  
}


void
GKVorticityBE::applyOp(LevelData<FArrayBox>&       a_out,
                       const LevelData<FArrayBox>& a_in,
                       bool                        a_homogeneous )
{
   const DisjointBoxLayout& grids = a_out.disjointBoxLayout();

   computeFluxDivergence(a_in, a_out, a_homogeneous);

   // Add high-order corrections (vorticity parallel and perp diffusion)
   if (m_include_high_order_corr) {
      
      LevelData<FArrayBox> pol_den_correction(grids, 1, IntVect::Zero);
      applyHighOrderCorrectionOp(pol_den_correction, a_in, a_homogeneous);

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_out[dit] += pol_den_correction[dit];
      }
   }
   
   // Add linear term from the first-order Backward Euler solve
   // used in oldVorticy model if operatorSplitMethod is used
   if (m_model == "VorticityDiffusion" ) {

     for (DataIterator dit(grids); dit.ok(); ++dit) {
         FArrayBox shift(grids[dit],1);
         shift.copy(a_in[dit]);
         shift.mult(m_beta);
         a_out[dit] += shift;
      }
   }
}

void
GKVorticityBE::applyHighOrderCorrectionOp(LevelData<FArrayBox>&       a_out,
                                          const LevelData<FArrayBox>& a_in,
                                          bool                        a_homogeneous )
{
   /*
    This function computes high-order corrections by first
    computing the vorticity quantity, where potential bcs are used.
    Then elliptic operator is applied to the vorticity quantity where
    either vorticity bcs or extrapolation bcs are used
    */
   
   const DisjointBoxLayout& grids = m_geometry.grids();

   // Compute negative vorticity; boundary data for the elliptic
   // operator was updated with phi BCs as part of the
   // earlier call to setVorticityOperatorCoefficients()
   LevelData<FArrayBox> negative_vorticity(grids, 1, IntVect::Zero);
   if ( !m_low_pollution ) {
      computeFluxDivergenceWithCoeff(negative_vorticity,
                                     a_in,
                                     m_M_unmapped,
                                     a_homogeneous,
                                     false);
   }
   else {
      computeFluxDivergenceWithCoeff(negative_vorticity,
                                     a_in,
                                     m_M_mapped,
                                     a_homogeneous,
                                     false);
   }

  
   
   // Update boundary data using vorticity bcs; presently only homogeneous
   // vorticity BCs are supported; overwise the code will stuck
   // in infinite recursion (updateBoundaryData->computeBcDivergence->
   // applyOp->applyHighOrderCorrectionOp->updateBoundaryData loop)
   if (m_use_vorticity_bcs) {
      bool homogeneous_bcs = true;
      updateBoundaryData(m_N2_unmapped, *m_vorticity_bc, homogeneous_bcs = true);
   }
   
    // Compute -(par_grad * sigma * par_grad + perp_grad D perp_grad) of negative_vorticity
   if ( !m_low_pollution ) {

      computeFluxDivergenceWithCoeff(a_out,
                                     negative_vorticity,
                                     m_N2_unmapped,
                                     a_homogeneous,
                                     m_subtract_fs_par_div,
                                     !m_use_vorticity_bcs);
   }
   else {

      computeFluxDivergenceWithCoeff(a_out,
                                     negative_vorticity,
                                     m_N2_mapped,
                                     a_homogeneous,
                                     m_subtract_fs_par_div,
                                     !m_use_vorticity_bcs);
   }

   // Return boundary data to be consistent with low-order potential bcs.
   // Do not update m_bc_divergence, since we did that earlier
   // (as part of the setVorticityOperatorCoefficients() call)
   if (m_use_vorticity_bcs) {
      bool homogeneous_bcs = true;
      updateBoundaryData(m_unmapped_coefficients, *m_potential_bc, homogeneous_bcs = true);
    }

}

void
GKVorticityBE::computeBcDivergence( LevelData<FArrayBox>& a_out )
{

   // Can't we use base class function for this -- no we can't!
   // First, base class calls computeFluxDivergence, not applyOp
   // Second, the case class computeBcDivergence() function seems
   // to call the base class computeFluxDivergence(), and not
   // its overwritten derived class version
   const DisjointBoxLayout& grids = a_out.disjointBoxLayout();
   LevelData<FArrayBox> phi(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi[dit].setVal(0.);
   }

   applyOp(a_out, phi, false);
}

void
GKVorticityBE::parseParameters( const ParmParse&   a_pp )
{
}


void
GKVorticityBE::printParameters()
{
}

#include "NamespaceFooter.H"
