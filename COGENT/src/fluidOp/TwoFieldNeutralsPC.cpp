#include "TwoFieldNeutralsPC.H"
#include "GKPoissonF_F.H"
#include "FourthOrderUtil.H"
#include "SimpleDivergence.H"
#include "Directions.H"
#include "ToroidalBlockLevelExchangeCenter.H"
#include "AMRIO.H"

//#define VERIFY_MATRIX


#include "NamespaceHeader.H"
 

TwoFieldNeutralsPC::TwoFieldNeutralsPC(const ParmParse&   a_pp,
                                       const MagGeom&     a_geom )
   : EllipticOp(a_pp, a_geom),
     m_diffusion_coeffs_defined(false),
     m_viscosity_coeffs_defined(false),
     m_bottomLeft_coeffs_defined(false),
     m_imex_preconditioner(NULL)
{



   const DisjointBoxLayout& grids = a_geom.grids();

   m_volume_reciprocal.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume_reciprocal);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_volume_reciprocal[dit].invert(1.0);
   }

   m_topLeft_mapped_tensor.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_bottomLeft_mapped_tensor.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_bottomRight_mapped_tensor.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_topLeft_unmapped_tensor.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_bottomLeft_unmapped_tensor.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_bottomRight_unmapped_tensor.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_topLeft_beta.define(grids, 1, IntVect::Zero);
   m_bottomLeft_beta.define(grids, 1, IntVect::Zero);
   m_bottomRight_beta.define(grids, 1, IntVect::Zero);

   m_verbosity = false;

   m_precond_verbose = true;
   m_precond_precond_verbose = true;

   int preconditioner_order = 2;

   if ( a_geom.extrablockExchange() ) {
      m_mblx_ptr = new MagFluxAlignedMBLEXCenter;
      m_mblx_ptr->define(&a_geom, preconditioner_order, preconditioner_order);
   }
   else {
      m_mblx_ptr = NULL;
   }

   if ( m_imex_preconditioner == NULL ) {
      m_imex_preconditioner = new MBHypreSolver(m_geometry, 2, preconditioner_order, m_mblx_ptr);
      ParmParse pp_precond( ((string)a_pp.prefix() + ".linear_solver.precond").c_str());
      m_imex_preconditioner->setMethodParams(pp_precond);
      m_imex_preconditioner->setConvergenceParams(pp_precond);
   }
   else {
      MayDay::Error("TwoFieldNeutralsPC::init(): m_imex_preconditioner is already allocated");
   }

}

TwoFieldNeutralsPC::~TwoFieldNeutralsPC()
{
   if ( m_imex_preconditioner ) delete m_imex_preconditioner;
   if ( m_mblx_ptr ) delete m_mblx_ptr;
}

void
TwoFieldNeutralsPC::updatePreconditioner( const EllipticOpBC& a_diffusion_bc,
                                          const EllipticOpBC& a_viscosity_bc )
{
   CH_TIME("TwoFieldNeutralsPC::updatePreconditioner");

   const DisjointBoxLayout& grids = m_geometry.grids();

   LevelData<FluxBox> zeroTensor(grids, SpaceDim*SpaceDim, IntVect::Unit);
   LevelData<FArrayBox> zeroBeta(grids, 1, IntVect::Zero);
   LevelData<FArrayBox> negVolRec(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      zeroTensor[dit].setVal(0.);
      zeroBeta[dit].setVal(0.);
      negVolRec[dit].copy(m_volume_reciprocal[dit]);
      negVolRec[dit].negate();
   }

   m_imex_preconditioner->initializeMatrix();

   // Upper-left coefficient
   m_imex_preconditioner->constructMatrixBlock(0, 0, m_volume_reciprocal, m_topLeft_mapped_tensor, m_topLeft_beta, a_diffusion_bc);

   // Upper-right coefficient
   m_imex_preconditioner->constructMatrixBlock(0, 1, zeroBeta, zeroTensor, zeroBeta, a_viscosity_bc);

   // Lower-left coefficient
   m_imex_preconditioner->constructMatrixBlock(1, 0, m_volume_reciprocal, m_bottomLeft_mapped_tensor, m_bottomLeft_beta, a_diffusion_bc);
   // m_imex_preconditioner->constructMatrixBlock(1, 0, negVolRec, zeroTensor, zeroBeta, a_diffusion_bc);

   // Lower-right coefficient
   m_imex_preconditioner->constructMatrixBlock(1, 1, m_volume_reciprocal, m_bottomRight_mapped_tensor, m_bottomRight_beta, a_viscosity_bc);

   m_imex_preconditioner->finalizeMatrix();
}

void TwoFieldNeutralsPC::solvePreconditioner( const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out )
{
  // Don't call this function. Call the other solvePreconditioner
}

void
TwoFieldNeutralsPC::solvePreconditioner( const LevelData<FArrayBox>& a_rhs_density,
                                         const LevelData<FArrayBox>& a_rhs_par_mom,
                                         LevelData<FArrayBox>&       a_z_density,
                                         LevelData<FArrayBox>&       a_z_par_mom )
{
  CH_TIME("TwoFieldNeutralsPC::solvePreconditioner");

  const DisjointBoxLayout& grids = m_geometry.grids();

  LevelData<FArrayBox> block_in(grids, 2, IntVect::Zero);
  LevelData<FArrayBox> block_out(grids, 2, IntVect::Zero);

  for (DataIterator dit(grids); dit.ok(); ++dit) {
    block_in[dit].copy(a_rhs_density[dit],0,0,1);
    block_in[dit].copy(a_rhs_par_mom[dit],0,1,1);
    block_out[dit].copy(a_z_density[dit],0,0,1);
    block_out[dit].copy(a_z_par_mom[dit],0,1,1);
  }

  m_imex_preconditioner->solve(block_in, block_out, true);

  for (DataIterator dit(grids); dit.ok(); ++dit) {
    a_z_density[dit].copy(block_out[dit],0,0,1);
    a_z_par_mom[dit].copy(block_out[dit],1,0,1);
  }
}

void
TwoFieldNeutralsPC::updateSourceTermComponents( const LevelData<FArrayBox>& a_topLeft_beta,
                                                const LevelData<FArrayBox>& a_bottomLeft_beta,
                                                const LevelData<FArrayBox>& a_bottomRight_beta )
{
   CH_TIME("TwoFieldNeutralsPC::updateSourceTermComponents");

   a_topLeft_beta.copyTo(m_topLeft_beta);
   a_bottomLeft_beta.copyTo(m_bottomLeft_beta);
   a_bottomRight_beta.copyTo(m_bottomRight_beta);
}

void
TwoFieldNeutralsPC::updateBottomRightCoefficients( const LevelData<FluxBox>& a_unmapped_coeff,
                                                   const LevelData<FluxBox>& a_mapped_coeff)
{
   CH_TIME("TwoFieldNeutralsPC::updateBottomRightCoefficients");

   const DisjointBoxLayout & grids = m_geometry.grids();
   
   const IntVect grown_ghosts = m_bottomRight_mapped_tensor.ghostVect() + IntVect::Unit;
   CH_assert(grown_ghosts == 2*IntVect::Unit);
   CH_assert(a_mapped_coeff.ghostVect() >= 2*IntVect::Unit);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      m_bottomRight_mapped_tensor[dit].copy(a_mapped_coeff[dit]);
      m_bottomRight_unmapped_tensor[dit].copy(a_unmapped_coeff[dit]);
   }

   
   // The mapped coefficients must now be converted to face averages, which
   // requires a layer of transverse ghost faces that we don't have at the
   // physical boundary, so we need to extrapolate them. DO WE NEED THIS?
   m_geometry.fillTransverseGhosts(m_bottomRight_mapped_tensor, false);
   m_geometry.fillTransverseGhosts(m_bottomRight_unmapped_tensor, false);

   m_viscosity_coeffs_defined = true;
}

void
TwoFieldNeutralsPC::updateBottomLeftCoefficients( const LevelData<FluxBox>& a_unmapped_coeff,
                                                  const LevelData<FluxBox>& a_mapped_coeff)
{
   CH_TIME("TwoFieldNeutralsPC::updateBottomLeftCoefficients");

   const DisjointBoxLayout & grids = m_geometry.grids();
   
   const IntVect grown_ghosts = m_bottomLeft_mapped_tensor.ghostVect() + IntVect::Unit;
   CH_assert(grown_ghosts == 2*IntVect::Unit);
   CH_assert(a_mapped_coeff.ghostVect() >= 2*IntVect::Unit);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      m_bottomLeft_mapped_tensor[dit].copy(a_mapped_coeff[dit]);
      m_bottomLeft_unmapped_tensor[dit].copy(a_unmapped_coeff[dit]);
   }

   // The mapped coefficients must now be converted to face averages, which
   // requires a layer of transverse ghost faces that we don't have at the
   // physical boundary, so we need to extrapolate them. DO WE NEED THIS?
   m_geometry.fillTransverseGhosts(m_bottomLeft_mapped_tensor, false);
   m_geometry.fillTransverseGhosts(m_bottomLeft_unmapped_tensor, false);

   m_bottomLeft_coeffs_defined = true;
}

void
TwoFieldNeutralsPC::updateTopLeftCoefficients( const LevelData<FluxBox>& a_unmapped_coeff,
                                               const LevelData<FluxBox>& a_mapped_coeff)
{
   CH_TIME("TwoFieldNeutralsPC::updateDiffusionCoefficients");

   const DisjointBoxLayout & grids = m_geometry.grids();
   
   const IntVect grown_ghosts = m_topLeft_mapped_tensor.ghostVect() + IntVect::Unit;
   CH_assert(grown_ghosts == 2*IntVect::Unit);
   CH_assert(a_mapped_coeff.ghostVect() >= 2*IntVect::Unit);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      m_topLeft_mapped_tensor[dit].copy(a_mapped_coeff[dit]);
      m_topLeft_unmapped_tensor[dit].copy(a_unmapped_coeff[dit]);
   }
   
   // The mapped coefficients must now be converted to face averages, which
   // requires a layer of transverse ghost faces that we don't have at the
   // physical boundary, so we need to extrapolate them. DO WE NEED THIS?
   m_geometry.fillTransverseGhosts(m_topLeft_mapped_tensor, false);
   m_geometry.fillTransverseGhosts(m_topLeft_unmapped_tensor, false);

   m_diffusion_coeffs_defined = true;
}

void
TwoFieldNeutralsPC::setImExPreconditionerConvergenceParams( const double a_tol,
                                                            const int    a_max_iter,
                                                            const double a_precond_tol,
                                                            const int    a_precond_max_iter )
{
  m_imex_preconditioner->setConvergenceParams(a_tol, a_max_iter, m_precond_verbose,
                                              a_precond_tol, a_precond_max_iter,
                                              m_precond_precond_verbose);
}

void
TwoFieldNeutralsPC::multiplyCoefficients( LevelData<FluxBox>& a_data,
                                          const bool          a_mapped_coeff ) const
{

}

#include "NamespaceFooter.H"
