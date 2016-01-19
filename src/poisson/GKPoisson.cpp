#include "GKPoisson.H"
#include "GKPoissonF_F.H"
#include "FourthOrderUtil.H"
#include "FluxSurface.H"
#include "Directions.H"

#include "NamespaceHeader.H"

const char* GKPoisson::pp_name = {"gkpoisson"};


GKPoisson::GKPoisson( const ParmParse&   a_pp,
                      const MagGeom&     a_geom,
                      const Real         a_larmor_number,
                      const Real         a_debye_number )
   : FieldSolver(a_pp, a_geom),
     m_solver(a_geom, 2),
     m_larmor_number2(a_larmor_number*a_larmor_number),
     m_debye_number2(a_debye_number*a_debye_number)
{
   // We give the coefficients one ghost cell layer so that the
   // second-order centered difference formula can be used to compute
   // the transverse gradients needed for the fourth-order formulas even
   // at box boundaries.
   const DisjointBoxLayout& grids = a_geom.grids();
   m_mapped_coefficients.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_unmapped_coefficients.define(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);

   // If there is more than one block, construct the multiblock exchange object
   if ( m_geometry.coordSysPtr()->numBlocks() > 1 ) {
      m_mblex_density_Ptr = new MultiBlockLevelExchangeAverage();
      m_density_interp_ghosts = 2;
      int spaceOrder = 4;
      m_mblex_density_Ptr->define(&a_geom, m_density_interp_ghosts, spaceOrder);
   }
   else {
      m_density_interp_ghosts = 0;
      m_mblex_density_Ptr = NULL;
   }
}
      


GKPoisson::~GKPoisson()
{
   if (m_mblex_density_Ptr) delete m_mblex_density_Ptr;
}



void
GKPoisson::setOperatorCoefficients( const LevelData<FArrayBox>& a_ion_mass_density,
                                    const PotentialBC&          a_bc )
{
   computeCoefficients( a_ion_mass_density, m_mapped_coefficients, m_unmapped_coefficients );

   setBcDivergence( a_bc, m_bc_divergence );

   m_solver.constructMatrix(m_mapped_coefficients, a_bc, false, false);
}



void
GKPoisson::computeCoefficients( const LevelData<FArrayBox>& a_ion_mass_density,
                                LevelData<FluxBox>&         a_mapped_coefficients,
                                LevelData<FluxBox>&         a_unmapped_coefficients )
{
   const DisjointBoxLayout& grids( m_geometry.grids() );

   LevelData<FluxBox> density_sum_face( grids, 1, 2*IntVect::Unit );
   LevelData<FArrayBox> density_sum_cell( grids, 1, 4*IntVect::Unit );

   a_ion_mass_density.copyTo( density_sum_cell );

   // This fills two codim 1 ghost cell layers at all block boundaries
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
     const ProblemDomain& domain = block_coord_sys.domain();

     fourthOrderCellExtrapAtDomainBdry(density_sum_cell[dit], domain, grids[dit]);
   }

   // This fills extrablock ghosts and performs an exchange.  We now have 4 layers
   // of ghost cells filled, except at the physical boundary were we only have 2.
   m_geometry.fillInternalGhosts( density_sum_cell );

   // Convert from cell averages to face centers.  We lose 2 ghost cells in this
   // process, so we now only have data on 2 ghost cell layers, except at the
   // physical boundary where we have none.
   fourthOrderCellToFaceCenters(density_sum_face, density_sum_cell);

   const IntVect grown_ghosts( m_mapped_coefficients.ghostVect() + IntVect::Unit );
   CH_assert(grown_ghosts == 2*IntVect::Unit);
   const LevelData<FluxBox>& BFieldMag = m_geometry.getFCBFieldMag();
   CH_assert(BFieldMag.ghostVect()>=grown_ghosts);
   const LevelData<FluxBox>& BFieldDir = m_geometry.getFCBFieldDir();
   CH_assert(BFieldDir.ghostVect()>=grown_ghosts);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     FluxBox& this_density_sum = density_sum_face[dit];
     FluxBox& this_coef = a_unmapped_coefficients[dit];
     const FluxBox& this_BMag = BFieldMag[dit];
     const FluxBox& this_bunit = BFieldDir[dit];

     for (int dir=0; dir<SpaceDim; ++dir) {
       FArrayBox& this_coef_dir = this_coef[dir];
       const Box& box_dir = this_coef_dir.box();

       FORT_COMPUTE_GKP_COEFFICIENTS(CHF_BOX(box_dir),
                                     CHF_CONST_REAL(m_debye_number2),
                                     CHF_CONST_REAL(m_larmor_number2),
                                     CHF_CONST_FRA(this_bunit[dir]),
                                     CHF_CONST_FRA1(this_BMag[dir],0),
                                     CHF_CONST_FRA1(this_density_sum[dir],0),
                                     CHF_FRA(this_coef_dir));
     }
   }

   // Compute the mapped coefficients

   LevelData<FluxBox> grown_mapped_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);

   LevelData<FluxBox> N(grids, SpaceDim*SpaceDim, grown_ghosts);
   m_geometry.getPointwiseN(N);

   LevelData<FluxBox> NJinverse(grids, SpaceDim*SpaceDim, grown_ghosts);
   m_geometry.getPointwiseNJinverse(NJinverse);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     FluxBox& this_unmapped_coef = a_unmapped_coefficients[dit];
     FluxBox& this_mapped_coef = grown_mapped_coefficients[dit];
     FluxBox& this_N = N[dit];
     FluxBox& this_NJinverse = NJinverse[dit];

     for (int dir=0; dir<SpaceDim; ++dir) {
       FArrayBox& this_mapped_coef_dir = this_mapped_coef[dir];
       const Box& box_dir = this_mapped_coef_dir.box();

       FORT_COMPUTE_MAPPED_GKP_COEFFICIENTS(CHF_BOX(box_dir),
                                            CHF_CONST_FRA(this_unmapped_coef[dir]),
                                            CHF_CONST_FRA(this_N[dir]),
                                            CHF_CONST_FRA(this_NJinverse[dir]),
                                            CHF_FRA(this_mapped_coef_dir));
     }
   }

   a_unmapped_coefficients.exchange();
   grown_mapped_coefficients.exchange();

   // The mapped coefficients must now be converted to face averages, which
   // requires a layer of transverse ghost faces that we don't have at the
   // physical boundary, so we need to extrapolate them.
   m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);

   // Convert the mapped coefficients from face-centered to face-averaged
   fourthOrderAverage(grown_mapped_coefficients);

   m_geometry.fillTransverseGhosts(a_unmapped_coefficients, false);
   m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
     a_mapped_coefficients[dit].copy(grown_mapped_coefficients[dit]);
   }
}


void
GKPoisson::fillDensityGhosts(LevelData<FArrayBox>& a_density) const
{
   a_density.exchange();

   if (m_mblex_density_Ptr) {
      CH_assert(a_density.ghostVect() <= IntVect(D_DECL(m_density_interp_ghosts,
                                                        m_density_interp_ghosts,
                                                        m_density_interp_ghosts)));
      m_mblex_density_Ptr->interpGhosts(a_density);
   }
}



void
GKPoisson::getMinMax(LevelData<FArrayBox>& a_data, double& a_min, double& a_max) const
{
   const DisjointBoxLayout & grids = a_data.disjointBoxLayout();

   a_min =  DBL_MAX;
   a_max = -DBL_MAX;
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      FArrayBox& this_data = a_data[dit];
      Box box(grids[dit]);
      if (this_data.min(box) < a_min) a_min = this_data.min(box);
      if (this_data.max(box) > a_max) a_max = this_data.max(box);
   }   

#ifdef CH_MPI
   double local_min = a_min;
   MPI_Allreduce(&local_min, &a_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   double local_max = a_max;
   MPI_Allreduce(&local_max, &a_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}



void
GKPoisson::getMinMax(LevelData<FluxBox>& a_data, double& a_min, double& a_max) const
{
   const DisjointBoxLayout & grids = a_data.disjointBoxLayout();

   a_min =  DBL_MAX;
   a_max = -DBL_MAX;
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         Box box(surroundingNodes(grids[dit],dir));
         FArrayBox& this_data = a_data[dit][dir];
         if (this_data.min(box) < a_min) a_min = this_data.min(box);
         if (this_data.max(box) > a_max) a_max = this_data.max(box);
      }
   }   

#ifdef CH_MPI
   double local_min = a_min;
   MPI_Allreduce(&local_min, &a_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   double local_max = a_max;
   MPI_Allreduce(&local_max, &a_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}



void
GKPoisson::solvePreconditioner( const LevelData<FArrayBox>& a_in,
                                LevelData<FArrayBox>&       a_out )
{
   m_solver.solveWithMultigrid(a_in, m_amg_tol, m_amg_max_iter, m_verbose, a_out);
}



void
GKPoisson::multiplyUnmappedCoefficients( LevelData<FluxBox>& a_data ) const
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
         const FArrayBox& this_coef_dir = m_unmapped_coefficients[dit][dir];

         FArrayBox tmp(this_data_dir.box(),1);
         for (int i=0; i<SpaceDim; ++i) {
            for (int j=0; j<SpaceDim; ++j) {
               tmp.copy(this_coef_dir, SpaceDim*i+j, 0, 1);
               tmp.mult(this_saved_data_dir, j, 0, 1);
               this_data_dir.plus(tmp, 0, i, 1);
            }
         }
      }
   }

   a_data.exchange();
}



#include "NamespaceFooter.H"
