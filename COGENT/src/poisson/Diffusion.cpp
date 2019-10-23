#include "Diffusion.H"
#include "DiffusionF_F.H"
#include "FourthOrderUtil.H"

#include "NamespaceHeader.H"

const char* Diffusion::pp_name = {"diffusion"};


Diffusion::Diffusion( const ParmParse&  a_pp,
                      const MagGeom&    a_geom )
   : EllipticOp(a_pp, a_geom),
     m_coefficients_defined(false)
{
   // We give the coefficients one ghost cell layer so that the
   // second-order centered difference formula can be used to compute
   // the transverse gradients needed for the fourth-order formulas even
   // at box boundaries.
   const DisjointBoxLayout& grids = a_geom.grids();
   m_mapped_coefficients.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_unmapped_coefficients.define(grids, SpaceDim*SpaceDim, IntVect::Unit);

   m_volume_reciprocal.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume_reciprocal);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_volume_reciprocal[dit].invert(1.0);
   }
   
   int preconditioner_order = 2;

   if ( a_geom.extrablockExchange() ) {
      m_mblx_ptr = new MagFluxAlignedMBLEXCenter;
      m_mblx_ptr->define(&a_geom, preconditioner_order, preconditioner_order);
   }
   else {
      m_mblx_ptr = NULL;
   }

#ifdef with_petsc
   m_preconditioner      = new MBPETScSolver(a_geom, preconditioner_order, m_mblx_ptr);
   m_imex_preconditioner = new MBPETScSolver(a_geom, preconditioner_order, m_mblx_ptr);
#else
   m_preconditioner      = new MBHypreSolver(a_geom, preconditioner_order, m_mblx_ptr);
   m_imex_preconditioner = new MBHypreSolver(a_geom, preconditioner_order, m_mblx_ptr);
#endif

   m_preconditioner->setMethodParams(m_precond_method, m_precond_precond_method);
   m_imex_preconditioner->setMethodParams(m_precond_method, m_precond_precond_method);
}
      


Diffusion::~Diffusion()
{
   if (m_mblx_ptr) delete m_mblx_ptr;
   if (m_preconditioner) delete m_preconditioner;
   if (m_imex_preconditioner) delete m_imex_preconditioner;
}


void
Diffusion::updateBoundaries( const EllipticOpBC&  a_bc )
{
   CH_TIME("Diffusion::updateBoundaries");
   CH_assert(m_coefficients_defined);
   
   setBc(a_bc);
   
   computeBcDivergence( m_bc_divergence );
}


void
Diffusion::setOperatorCoefficients( const LevelData<FluxBox>&           a_unmapped_coefficients,
                                    const LevelData<FluxBox>&           a_mapped_coefficients,
                                    const EllipticOpBC&                 a_bc )
{
   CH_TIME("Diffusion::setOperatorCoefficients");
   
   const DisjointBoxLayout & grids = m_geometry.grids();
   
   const IntVect grown_ghosts = m_mapped_coefficients.ghostVect() + IntVect::Unit;
   CH_assert(grown_ghosts == 2*IntVect::Unit);
   CH_assert(a_mapped_coefficients.ghostVect() >= 2*IntVect::Unit);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      m_unmapped_coefficients[dit].copy(a_unmapped_coefficients[dit]);
      m_mapped_coefficients[dit].copy(a_mapped_coefficients[dit]);
   }
   
   // The mapped coefficients must now be converted to face averages, which
   // requires a layer of transverse ghost faces that we don't have at the
   // physical boundary, so we need to extrapolate them. DO WE NEED THIS?
   m_geometry.fillTransverseGhosts(m_mapped_coefficients, false);
   m_geometry.fillTransverseGhosts(m_unmapped_coefficients, false);

   // Convert the mapped coefficients from face-centered to face-averaged
   if (!m_second_order) fourthOrderAverage(m_mapped_coefficients);
   
   m_coefficients_defined = true;

   updateBoundaries(a_bc);
}


void
Diffusion::updateImExPreconditioner( const LevelData<FArrayBox>&  a_mshift,
                                     const EllipticOpBC&          a_bc )
{
   LevelData<FArrayBox> beta(m_geometry.grids(), 1, IntVect::Zero);

   for (DataIterator dit(beta.dataIterator()); dit.ok(); ++dit ) {
      beta[dit].copy(a_mshift[dit]);
   }

   m_imex_preconditioner->constructMatrix(m_volume_reciprocal, m_mapped_coefficients, beta, a_bc);

}

void
Diffusion::setPreconditionerConvergenceParams( const double a_tol,
                                               const int    a_max_iter,
                                               const double a_precond_tol,
                                               const int    a_precond_max_iter )
{
   m_preconditioner->setConvergenceParams(a_tol, a_max_iter, m_precond_verbose,
                                          a_precond_tol, a_precond_max_iter,
                                          m_precond_precond_verbose);
}


void
Diffusion::setImExPreconditionerConvergenceParams( const double a_tol,
                                                   const int    a_max_iter,
                                                   const double a_precond_tol,
                                                   const int    a_precond_max_iter )
{
   m_imex_preconditioner->setConvergenceParams(a_tol, a_max_iter, m_precond_verbose,
                                               a_precond_tol, a_precond_max_iter,
                                               m_precond_precond_verbose);
}



void
Diffusion::solvePreconditioner( const LevelData<FArrayBox>& a_in,
                                LevelData<FArrayBox>&       a_out )
{
   m_preconditioner->solve(a_in, a_out, true);
}



void
Diffusion::solveImExPreconditioner( const LevelData<FArrayBox>& a_in,
                                    LevelData<FArrayBox>&       a_out )
{
   m_imex_preconditioner->solve(a_in, a_out, true);
}



void
Diffusion::multiplyCoefficients( LevelData<FluxBox>& a_data,
                                 const bool          a_mapped_coeff ) const
{
   CH_assert(a_data.ghostVect() <= m_unmapped_coefficients.ghostVect());
   CH_assert(m_coefficients_defined);

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      FluxBox& this_data = a_data[dit];
      FluxBox saved_data(this_data.box(),SpaceDim);
      saved_data.copy(this_data);
      this_data.setVal(0.);
      for (int dir=0; dir<SpaceDim; ++dir) { 
         FArrayBox& this_data_dir = this_data[dir];
         FArrayBox& this_saved_data_dir = saved_data[dir];
         FArrayBox tmp(this_data_dir.box(),1);
         for (int i=0; i<SpaceDim; ++i) {
            for (int j=0; j<SpaceDim; ++j) {
               if (a_mapped_coeff) {
                  tmp.copy(m_mapped_coefficients[dit][dir], SpaceDim*i+j, 0, 1);
               }
               else {
                  tmp.copy(m_unmapped_coefficients[dit][dir], SpaceDim*i+j, 0, 1);
               }
               tmp.mult(this_saved_data_dir, j, 0, 1);
               this_data_dir.plus(tmp, 0, i, 1);
            }
         }
      }
   }

   a_data.exchange();
}


#include "NamespaceFooter.H"
