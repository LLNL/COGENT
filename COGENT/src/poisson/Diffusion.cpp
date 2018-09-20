#include "Diffusion.H"
#include "DiffusionF_F.H"
#include "FourthOrderUtil.H"

#include "NamespaceHeader.H"

const char* Diffusion::pp_name = {"diffusion"};


Diffusion::Diffusion( const ParmParse&  a_pp,
                      const MagGeom&    a_geom )
   : FieldSolver(a_pp, a_geom)
{
   // We give the coefficients one ghost cell layer so that the
   // second-order centered difference formula can be used to compute
   // the transverse gradients needed for the fourth-order formulas even
   // at box boundaries.
   const DisjointBoxLayout& grids = a_geom.grids();
   m_mapped_coefficients.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_unmapped_coefficients.define(grids, SpaceDim*SpaceDim, IntVect::Unit);

   m_mapped_volume_reciprocal.define(grids, SpaceDim*SpaceDim, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = a_geom.getBlockCoordSys(grids[dit]);
      m_mapped_volume_reciprocal[dit].setVal(1.0 / block_coord_sys.getMappedCellVolume());
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
   m_preconditioner = new MBPETScSolver(a_geom, preconditioner_order, m_mblx_ptr);
#else
   m_preconditioner = new MBHypreSolver(a_geom, preconditioner_order, m_mblx_ptr);
#endif
}
      


Diffusion::~Diffusion()
{
   if (m_preconditioner) delete m_preconditioner;
}


void
Diffusion::updateBoundaries( const PotentialBC&  a_bc )
{
   // N.B.: The operator coefficient must computed prior to calling this function
   CH_TIME("Diffusion::updateBoundaries");
   setBc(a_bc);
   
   computeBcDivergence( m_bc_divergence );
}


void
Diffusion::setOperatorCoefficients( const double                       a_scalar_factor,
                                    const RefCountedPtr<GridFunction>  a_coefficient_function,
                                    const PotentialBC&                 a_bc,
                                    const bool                         a_update_preconditioner )
{
   CH_TIME("Diffusion::setOperatorCoefficients");

   computeCoefficients(a_scalar_factor, a_coefficient_function, m_unmapped_coefficients, m_mapped_coefficients);

   updateBoundaries(a_bc);

   if ( a_update_preconditioner ) {
      m_preconditioner->constructMatrix(m_mapped_volume_reciprocal, m_mapped_coefficients, a_bc);
   }
}


void
Diffusion::computeCoefficients( const double                       a_scalar_factor,
                                const RefCountedPtr<GridFunction>  a_coefficient_function,
                                LevelData<FluxBox>&                a_unmapped_coefficients,
                                LevelData<FluxBox>&                a_mapped_coefficients )
{
   const DisjointBoxLayout & grids = m_geometry.grids();

   const IntVect grown_ghosts = a_mapped_coefficients.ghostVect() + IntVect::Unit;
   CH_assert(grown_ghosts == 2*IntVect::Unit);

   // Compute the coefficients at the face centers of a grown (by one ghost cell)
   // LevelData so that they can then be face averaged
   LevelData<FluxBox> grown_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);

   LevelData<FluxBox> N(grids, SpaceDim*SpaceDim, grown_ghosts);
   m_geometry.getPointwiseN(N);

   LevelData<FluxBox> NJinverse(grids, SpaceDim*SpaceDim, grown_ghosts);
   m_geometry.getPointwiseNJinverse(NJinverse);

   LevelData<FArrayBox> D_cell( grids, 1, grown_ghosts + IntVect::Unit );
   if ( !a_coefficient_function.isNull() ) {
      a_coefficient_function->assign( D_cell, m_geometry, 0.0);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_cell[dit] *= a_scalar_factor;
      }
   }
   else {
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_cell[dit].setVal(a_scalar_factor);
      }
   }

   m_geometry.fillInternalGhosts( D_cell );

   LevelData<FluxBox> D_face( grids, 1, grown_ghosts );
   fourthOrderCellToFaceCenters(D_face, D_cell);

   LevelData<FluxBox> D_tensor(grids, SpaceDim*SpaceDim, grown_ghosts);

   // Tensor is diagonal for now
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_D =        D_face[dit][dir];
         FArrayBox& this_D_tensor = D_tensor[dit][dir];
         this_D_tensor.setVal(0.);
         for (int n=0; n<SpaceDim; ++n) {
            this_D_tensor.copy(this_D,0,n*(SpaceDim+1),1);
         }
      }
   }      

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     FluxBox& this_coef = grown_coefficients[dit];

     for (int dir=0; dir<SpaceDim; ++dir) {
       FArrayBox& this_coef_dir = this_coef[dir];
       const Box& box_dir = this_coef_dir.box();

       FORT_COMPUTE_MAPPED_DIFFUSION_COEFFICIENTS(CHF_BOX(box_dir),
                                                  CHF_CONST_FRA(N[dit][dir]),
                                                  CHF_CONST_FRA(D_tensor[dit][dir]),
                                                  CHF_CONST_FRA(NJinverse[dit][dir]),
                                                  CHF_FRA(this_coef_dir));
     }
   }

   // The mapped coefficients must now be converted to face averages, which
   // requires a layer of transverse ghost faces that we don't have at the
   // physical boundary, so we need to extrapolate them.
   m_geometry.fillTransverseGhosts(grown_coefficients, false);

   // Convert the mapped coefficients from face-centered to face-averaged
   fourthOrderAverage(grown_coefficients);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     a_unmapped_coefficients[dit].copy(D_tensor[dit]);
     a_mapped_coefficients[dit].copy(grown_coefficients[dit]);
   }
}


void
Diffusion::setPreconditionerConvergenceParams( const double a_tol,
                                               const int    a_max_iter,
                                               const double a_precond_tol,
                                               const int    a_precond_max_iter )
{
   m_preconditioner->setParams(m_precond_method, a_tol, a_max_iter, m_precond_verbose,
                               m_precond_precond_method, a_precond_tol, a_precond_max_iter,
                               m_precond_precond_verbose);
}



void
Diffusion::solvePreconditioner( const LevelData<FArrayBox>& a_in,
                                LevelData<FArrayBox>&       a_out )
{
   m_preconditioner->solve(a_in, a_out, true);
}



void
Diffusion::multiplyCoefficients( LevelData<FluxBox>& a_data,
                                 const bool          a_mapped_coeff ) const
{
   CH_assert(a_data.ghostVect() <= m_unmapped_coefficients.ghostVect());

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
