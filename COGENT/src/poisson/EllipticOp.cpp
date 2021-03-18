#include "EllipticOp.H"
#include "BlockRegister.H"
#include "EdgeToCell.H"
#include "FourthOrderUtil.H"
#include "Directions.H"
#include "FluxSurface.H"
#include "BiCGStabSolver.H"
#include "GMRESSolver.H"
#include "SingleNullCoordSys.H"
#include "MBSolverF_F.H"
#include "SimpleDivergence.H"
#include "CoDimCopyManager.H"
#include "DataArray.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM POL_DIM
#include "BilinearInterp.H"
#include "BicubicInterp.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H"

#undef INTERP_TO_NODES_WITH_BCS


EllipticOp::EllipticOp( const ParmParse& a_pp,
                        const MagGeom&   a_geometry )
   : m_geometry(a_geometry)
{
   if ( a_pp.query("second_order", m_second_order) == 0 ) m_second_order = false;
   if ( a_pp.query("low_pollution", m_low_pollution) == 0 ) m_low_pollution = false;

   init(a_pp);
}
      

EllipticOp::EllipticOp( const ParmParse&  a_pp,
                        const MagGeom&    a_geometry,
                        const bool        a_second_order,
                        const bool        a_low_pollution )
   : m_geometry(a_geometry),
     m_second_order(a_second_order),
     m_low_pollution(a_low_pollution)
{
   init(a_pp);
}
      

EllipticOp::~EllipticOp()
{
   if (m_Chombo_solver) delete m_Chombo_solver;
   if (m_mblex_potential_Ptr) delete m_mblex_potential_Ptr;
}


void
EllipticOp::init( const ParmParse& a_pp )
{
   const DisjointBoxLayout& grids = m_geometry.grids();

   // Check minimum box size.  Needed for boundary extrapolation.
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         if(grids[dit].size(dir) < 3) {
            MayDay::Error( "Poisson: Minimum box dimension is 3!" );
         }
      }
   }

   if (m_second_order) {
      m_num_potential_ghosts = 2;
   }
   else {
      m_num_potential_ghosts = 3;
   }

   m_Chombo_solver = NULL;
   ParmParse pp_linear_solver( ((string)a_pp.prefix() + ".linear_solver").c_str());
   defineSolver( pp_linear_solver );

   // If there is more than one block, construct the multiblock exchange object
   if ( (m_geometry.coordSysPtr()->numBlocks() > 1) && !m_geometry.extrablockExchange() && !m_geometry.shearedMBGeom()) {
     m_mblex_potential_Ptr = new MultiBlockLevelExchangeCenter();
     int spaceOrder = 4;
     m_mblex_potential_Ptr->define(&m_geometry, m_num_potential_ghosts, spaceOrder);
   }
   else {
     m_mblex_potential_Ptr = NULL;
   }

   m_volume.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume);

   m_volume_reciprocal.define(m_volume);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_volume_reciprocal[dit].invert(1.);
   }

   m_bc_divergence.define(grids, 1, IntVect::Zero);
}
      

void
EllipticOp::defineSolver( const ParmParse& a_pp )
{
   // Defines the linear solver from ParmParse data or defaults

   string method;
   if ( a_pp.query("method", method) == 0 ) method = "BiCGStab";

   double tol;
   if ( a_pp.query("tol", tol) == 0 ) tol = 1.e-6;

   int max_iter;
   if ( a_pp.query("max_iter", max_iter) == 0 ) max_iter = 20;

   bool verbose;
   if ( a_pp.query("verbose", verbose) == 0 ) verbose = false;

   defineSolver(method, tol, max_iter, verbose);
}


void
EllipticOp::defineSolver( const string&  a_method,
                          const double   a_tol,
                          const int      a_max_iter,
                          const bool     a_verbose )
{
   if ( m_Chombo_solver ) delete m_Chombo_solver;

   m_method = a_method;
   m_tol = a_tol;
   m_max_iter = a_max_iter;
   m_verbose = a_verbose;

   if ( m_method == "BiCGStab" ) {
      m_Chombo_solver = new BiCGStabSolver< LevelData<FArrayBox> >;
   }
   else {
      m_Chombo_solver = new GMRESSolver< LevelData<FArrayBox> >;
   }

   m_Chombo_solver->define(this, true);
}


void
EllipticOp::solveWithBCs( LevelData<FArrayBox>&        a_solution,
                          const LevelData<FArrayBox>&  a_rhs )
{
   LevelData<FArrayBox> rhs_with_bcs;
   rhs_with_bcs.define(a_rhs);

   // Add any inhomogeneous boundary values to the right-hand side
   subtractBcDivergence(rhs_with_bcs);

   solve( rhs_with_bcs, a_solution );
}


void
EllipticOp::solve( const LevelData<FArrayBox>&  a_rhs,
                   LevelData<FArrayBox>&        a_solution )
{
   if ( m_method == "BiCGStab" ) {
      ((BiCGStabSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_eps = m_tol;
      ((BiCGStabSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_imax = m_max_iter;
      ((BiCGStabSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_verbosity = m_verbose? 5: 0;
   }
   else {
      ((GMRESSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_eps = m_tol;
      ((GMRESSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_imax = m_max_iter;
      ((GMRESSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_verbosity = m_verbose? 5: 0;
   }

   setToZero(a_solution);
   m_Chombo_solver->solve(a_solution, a_rhs);

   if (m_verbose && procID() == 0) {
      if ( m_method == "BiCGStab" ) {
         int exit_status = ((BiCGStabSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_exitStatus;
         if ( exit_status == 1 ) {
            cout << "      --> BiCGStab converged to residual norm " 
                 << ((BiCGStabSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_residual << " in " 
                 << ((BiCGStabSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_iter << " iterations" 
                 << endl;
         }
         else {
            cout << "      --> BiCGStab solver did not converge, returning " << exit_status << endl;
         }
      }
      else if ( m_method == "GMRES" ) {
         int exit_status = ((GMRESSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_exitStatus;
         if ( exit_status == 1 ) {
            cout << "      --> GMRES converged to residual norm " 
                 << ((GMRESSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_residual << " in " 
                 << ((GMRESSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_iter << " iterations" 
                 << endl;
         }
         else {
            cout << "      --> GMRES solver did not converge, returning " << exit_status << endl;
         }
      }
   }
}



void
EllipticOp::compute3DFieldWithBCs( const LevelData<FArrayBox>&  a_phi,
                                   LevelData<FluxBox>&          a_field,
                                   const bool                   a_homogeneousBCs ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FluxBox> mapped_field(m_geometry.grids(), 3, a_field.ghostVect());
   computeMapped3DFieldWithBCs(a_phi, mapped_field, a_homogeneousBCs);

   m_geometry.unmap3DGradient(mapped_field, a_field);
}


void
EllipticOp::computePoloidalFieldWithBCs( const LevelData<FArrayBox>&  a_phi,
                                         LevelData<FluxBox>&          a_field,
                                         const bool                   a_homogeneousBCs ) const
{
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FluxBox> mapped_field(m_geometry.grids(), 2, a_field.ghostVect());
   computeMappedPoloidalFieldWithBCs(a_phi, mapped_field, a_homogeneousBCs);

   m_geometry.unmapPoloidalGradient(mapped_field, a_field);
}


void
EllipticOp::computeMapped3DFieldWithBCs( const LevelData<FArrayBox>&  a_phi,
                                         LevelData<FluxBox>&          a_field,
                                         const bool                   a_homogeneousBCs ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   // Make a temporary with ghost cells and copy the potential on valid cells

   const DisjointBoxLayout& grids = m_geometry.grids();
   CH_assert(m_num_potential_ghosts >=2);
   LevelData<FArrayBox> phi_cell(grids, 1, m_num_potential_ghosts*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi_cell[dit].setVal(0.);
      phi_cell[dit].copy(a_phi[dit],grids[dit]);
   }

   fillInternalGhosts(phi_cell);

   // Fill the physical boundary ghost cells
   bool extrapolate_from_interior = true;
   bool include_bvs = !a_homogeneousBCs;
   accumPhysicalGhosts( m_codim1_stencils, m_codim2_stencils, extrapolate_from_interior, include_bvs, phi_cell );

   int order = m_second_order? 2: 4;

   computeMapped3DFieldWithGhosts(phi_cell, a_field, order);

}


void
EllipticOp::computeMappedPoloidalFieldWithBCs( const LevelData<FArrayBox>&  a_phi,
                                               LevelData<FluxBox>&          a_field,
                                               const bool                   a_homogeneousBCs ) const
{
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   // Make a temporary with ghost cells and copy the potential on valid cells

   const DisjointBoxLayout& grids = m_geometry.grids();
   CH_assert(m_num_potential_ghosts >=2);
   LevelData<FArrayBox> phi_cell(grids, 1, m_num_potential_ghosts*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi_cell[dit].setVal(0.);
      phi_cell[dit].copy(a_phi[dit],grids[dit]);
   }

   fillInternalGhosts(phi_cell);

   // Fill the physical boundary ghost cells
   bool extrapolate_from_interior = true;
   bool include_bvs = !a_homogeneousBCs;
   accumPhysicalGhosts( m_codim1_stencils, m_codim2_stencils, extrapolate_from_interior, include_bvs, phi_cell );
         
   int order = m_second_order? 2: 4;

   computeMappedPoloidalFieldWithGhosts(phi_cell, a_field, order);
}


void
EllipticOp::computePoloidalField( const LevelData<FArrayBox>&  a_phi,
                                  LevelData<FArrayBox>&        a_field ) const
{
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FArrayBox> mapped_field(m_geometry.grids(), 2, IntVect::Unit);
   computeMappedPoloidalField(a_phi, mapped_field);

   m_geometry.unmapPoloidalGradient(mapped_field, a_field);
}


void
EllipticOp::computePoloidalField( const LevelData<FArrayBox>&  a_phi,
                                  LevelData<FluxBox>&          a_field ) const
{
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FluxBox> mapped_field(m_geometry.grids(), 2, IntVect::Unit);
   computeMappedPoloidalField(a_phi, mapped_field);

   m_geometry.unmapPoloidalGradient(mapped_field, a_field);
}


void
EllipticOp::computeMapped3DField( const LevelData<FArrayBox>&  a_phi,
                                  LevelData<FArrayBox>&        a_field ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   // Make a temporary with ghost cells and copy the potential on valid cells

   const DisjointBoxLayout & grids = m_geometry.grids();
   CH_assert(m_num_potential_ghosts >= 2);
   LevelData<FArrayBox> phi_cell(grids, 1, m_num_potential_ghosts*IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi_cell[dit].setVal(0.);
      phi_cell[dit].copy(a_phi[dit],grids[dit]);
   }

   fillInternalGhosts(phi_cell);

   // Extrapolate the potential at physical boundaries

   int order = m_second_order? 2: 4;

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   //const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = mag_coord_sys->boundaries();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = mag_coord_sys->whichBlock(grids[dit]);
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
      const ProblemDomain& domain = block_coord_sys.domain();
      
      for (int dir = 0; dir < SpaceDim; dir++) {
	for (SideIterator sit; sit.ok(); ++sit) {
          Side::LoHiSide side = sit();
	  
          if (mag_coord_sys->containsPhysicalBoundary(block_number, dir, side)) {
            IntVect grow_vec = 2*IntVect::Unit;
            grow_vec[dir] = 0;
            Box interior_box = grow(grids[dit],grow_vec);
            Box domain_box = grow(domain.domainBox(),grow_vec);
	    
            SpaceUtils::extrapBoundaryGhostsForCC(phi_cell[dit], interior_box, domain_box, dir, order, sign(side));
          }
        }
      }
   }

   if (m_geometry.mixedBoundaries()) {
     fillInternalGhosts(phi_cell);
   }
   
   computeMapped3DFieldWithGhosts(phi_cell, a_field, order);
}


void
EllipticOp::computeMappedPoloidalField( const LevelData<FArrayBox>&  a_phi,
                                        LevelData<FArrayBox>&        a_field ) const
{
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   // Make a temporary with ghost cells and copy the potential on valid cells

   const DisjointBoxLayout & grids = m_geometry.grids();
   CH_assert(m_num_potential_ghosts >= 2);
   LevelData<FArrayBox> phi_cell(grids, 1, m_num_potential_ghosts*IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi_cell[dit].setVal(0.);
      phi_cell[dit].copy(a_phi[dit],grids[dit]);
   }

   fillInternalGhosts(phi_cell);

   // Extrapolate the potential at physical boundaries

   int order = m_second_order? 2: 4;

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   //const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = mag_coord_sys->boundaries();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = mag_coord_sys->whichBlock(grids[dit]);
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
      const ProblemDomain& domain = block_coord_sys.domain();

      for (int dir = 0; dir < SpaceDim; dir++) {
        for (SideIterator sit; sit.ok(); ++sit) {
          Side::LoHiSide side = sit();

          if (mag_coord_sys->containsPhysicalBoundary(block_number, dir, side)) {
            IntVect grow_vec = 2*IntVect::Unit;
            grow_vec[dir] = 0;
            Box interior_box = grow(grids[dit],grow_vec);
            Box domain_box = grow(domain.domainBox(),grow_vec);

            SpaceUtils::extrapBoundaryGhostsForCC(phi_cell[dit], interior_box, domain_box, dir, order, sign(side));
          }
	}
      }
   }

   if (m_geometry.mixedBoundaries()) {
     fillInternalGhosts(phi_cell);
   }
     
   
   computeMappedPoloidalFieldWithGhosts(phi_cell, a_field, order);
}


void
EllipticOp::computeMappedPoloidalField( const LevelData<FArrayBox>&  a_phi,
                                        LevelData<FluxBox>&          a_field ) const
{
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   // Make a temporary with ghost cells and copy the potential on valid cells

   const DisjointBoxLayout & grids = m_geometry.grids();
   CH_assert(m_num_potential_ghosts >= 2);
   LevelData<FArrayBox> phi_cell(grids, 1, m_num_potential_ghosts*IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi_cell[dit].setVal(0.);
      phi_cell[dit].copy(a_phi[dit],grids[dit]);
   }

   fillInternalGhosts(phi_cell);

   // Extrapolate the potential to two cells at physical boundaries

   int order = m_second_order? 2: 4;

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   //const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = mag_coord_sys->boundaries();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = mag_coord_sys->whichBlock(grids[dit]);
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
      const ProblemDomain& domain = block_coord_sys.domain();

      for (int dir = 0; dir < SpaceDim; dir++) {
        for (SideIterator sit; sit.ok(); ++sit) {
          Side::LoHiSide side = sit();

          if (mag_coord_sys->containsPhysicalBoundary(block_number, dir, side)) {
            IntVect grow_vec = phi_cell.ghostVect();
            grow_vec[dir] = 0;
            Box interior_box = grow(grids[dit],grow_vec);
            Box domain_box = grow(domain.domainBox(),grow_vec);
      
            SpaceUtils::extrapBoundaryGhostsForFC(phi_cell[dit], interior_box, domain_box, dir, order, sign(side));
	  }
	}
      }
   }

   if (m_geometry.mixedBoundaries()) {
     fillInternalGhosts(phi_cell);
   }
      
   computeMappedPoloidalFieldWithGhosts(phi_cell, a_field, order);
}


void
EllipticOp::computeMapped3DField( const LevelData<FArrayBox>&  a_phi,
                                  LevelData<FluxBox>&          a_field ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   // Make a temporary with ghost cells and copy the potential on valid cells

   const DisjointBoxLayout & grids = m_geometry.grids();
   CH_assert(m_num_potential_ghosts >= 2);
   LevelData<FArrayBox> phi_cell(grids, 1, m_num_potential_ghosts*IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi_cell[dit].setVal(0.);
      phi_cell[dit].copy(a_phi[dit],grids[dit]);
   }

   fillInternalGhosts(phi_cell);

   // Extrapolate the potential to two cells at physical boundaries

   int order = m_second_order? 2: 4;

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   //const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = mag_coord_sys->boundaries();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = mag_coord_sys->whichBlock(grids[dit]);
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
      const ProblemDomain& domain = block_coord_sys.domain();

      for (int dir = 0; dir < SpaceDim; dir++) {
        for (SideIterator sit; sit.ok(); ++sit) {
          Side::LoHiSide side = sit();

          if (mag_coord_sys->containsPhysicalBoundary(block_number, dir, side)) {
            IntVect grow_vec = phi_cell.ghostVect();
            grow_vec[dir] = 0;
            Box interior_box = grow(grids[dit],grow_vec);
            Box domain_box = grow(domain.domainBox(),grow_vec);

            SpaceUtils::extrapBoundaryGhostsForFC(phi_cell[dit], interior_box, domain_box, dir, order, sign(side));
          }
	}
      }
   }

   if (m_geometry.mixedBoundaries()) {
     fillInternalGhosts(phi_cell);
   }

   computeMapped3DFieldWithGhosts(phi_cell, a_field, order);
}


void
EllipticOp::computeField( const LevelData<FArrayBox>&  a_phi,
                          LevelData<FArrayBox>&        a_field ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FArrayBox> mapped_field(a_field.disjointBoxLayout(), 3, a_field.ghostVect());
   computeMappedField( a_phi, mapped_field );

   m_geometry.unmapGradient(mapped_field, a_field);
}



void
EllipticOp::computeField( const LevelData<FArrayBox>&  a_phi,
                          LevelData<FluxBox>&          a_field ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FluxBox> mapped_field(a_field.disjointBoxLayout(), 3, a_field.ghostVect());
   computeMappedField( a_phi, mapped_field );

   m_geometry.unmapGradient(mapped_field, a_field);
}


void
EllipticOp::computeMappedField( const LevelData<FArrayBox>&  a_phi,
                                LevelData<FArrayBox>&        a_field ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

#if CFG_DIM==3
   computeMapped3DField(a_phi, a_field);
#else

   LevelData<FArrayBox> poloidal_field(a_field.disjointBoxLayout(), 2, IntVect::Unit);
   computeMappedPoloidalField(a_phi, poloidal_field);

   m_geometry.injectPoloidalVector(poloidal_field, a_field);
#endif

}


void
EllipticOp::computeMappedField( const LevelData<FArrayBox>&  a_phi,
                                LevelData<FluxBox>&          a_field ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

#if CFG_DIM==3
   computeMapped3DField(a_phi, a_field);
#else

   LevelData<FluxBox> poloidal_field(a_field.disjointBoxLayout(), 2, IntVect::Unit);
   computeMappedPoloidalField(a_phi, poloidal_field);

   m_geometry.injectPoloidalVector(poloidal_field, a_field);
#endif
}


void
EllipticOp::computeFluxDivergence( const LevelData<FArrayBox>&  a_in,
                                   LevelData<FArrayBox>&        a_out,
                                   const bool                   a_homogeneous_bcs,
                                   const bool                   a_extrap_to_ghosts )
{
   CH_TIME("EllipticOp::computeFluxDivergence");
   const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

   LevelData<FArrayBox> phi(grids, 1, 3*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi[dit].copy(a_in[dit]);
   }

   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);

   if (!m_low_pollution) {

     if (SpaceDim == 3) {
       if (a_extrap_to_ghosts) computeField(phi, flux);
       else compute3DFieldWithBCs(phi, flux, a_homogeneous_bcs);

     }
     else {
       if (a_extrap_to_ghosts) computePoloidalField(phi, flux);
       else computePoloidalFieldWithBCs(phi, flux, a_homogeneous_bcs);
     }

     // Multiply the field by the unmapped, face-centered GKP coefficients
     multiplyCoefficients(flux, false);

     m_geometry.fillTransversePhysicalGhosts(flux);

     m_geometry.applyAxisymmetricCorrection(flux);

     // Convert to face-averaged
     if (!m_second_order) fourthOrderAverage(flux);

     m_geometry.computeMappedGridDivergence(flux, a_out, !m_second_order);
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
     multiplyCoefficients(flux, true);

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

   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_out[dit] /= m_volume[dit];
   }
}


void
EllipticOp::computeBcDivergence( LevelData<FArrayBox>& a_out ) 
{
   const DisjointBoxLayout& grids = a_out.disjointBoxLayout();
   LevelData<FArrayBox> phi(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi[dit].setVal(0.);
   }

   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);

   if (!m_low_pollution) {

#if CFG_DIM==3
      LevelData<FluxBox> mapped_field(grids, CFG_DIM, IntVect::Unit);
      computeMapped3DFieldWithBCs( phi, mapped_field, false );
      m_geometry.unmap3DGradient(mapped_field, flux);
#else
      computePoloidalFieldWithBCs(phi, flux, false);
#endif
      
      // Multiply the field by the unmapped, face-centered GKP coefficients
      multiplyCoefficients(flux, false);
      m_geometry.applyAxisymmetricCorrection(flux);
   }
   else {
#if CFG_DIM==3
      computeMapped3DFieldWithBCs(phi, flux, false);
#else
      computeMappedPoloidalFieldWithBCs(phi, flux, false);
#endif
      
      // Multiply the field by the mapped, face-centered GKP coefficients
      multiplyCoefficients(flux, true);
   }

   // Convert to face-averaged if we're fourth-order
   if (!m_second_order) fourthOrderAverage(flux);

   if (!m_low_pollution) {
      m_geometry.computeMappedGridDivergence(flux, a_out, !m_second_order);
   }

   else {
      LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
      for (DataIterator dit(flux.dataIterator()); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         RealVect dx = block_coord_sys.dx();
         for (int dir=0; dir<SpaceDim; ++dir) {
            Real transverse_area_mapped(1.0);
            for (int tdir = 0; tdir<SpaceDim; tdir++) {
               if (tdir != dir) {
                  transverse_area_mapped *= dx[tdir];
               }
            }
            NTF_normal[dit][dir].copy(flux[dit][dir],dir,0,1);
            NTF_normal[dit][dir].mult(transverse_area_mapped);
         }
      }
      RealVect fakeDx = RealVect::Unit;
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         simpleDivergence(a_out[dit], NTF_normal[dit], grids[dit], fakeDx);
      }
   }
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_out[dit] /= m_volume[dit];
   }
}


void
EllipticOp::setBc( const EllipticOpBC& a_bc )
{
   m_codim1_stencils.clear();
   m_codim2_stencils.clear();
   constructBoundaryStencils(!m_second_order, a_bc, m_codim1_stencils, m_codim2_stencils );
}


void
EllipticOp::addBcDivergence( LevelData<FArrayBox>& a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit] += m_bc_divergence[dit];
   }
}


void
EllipticOp::subtractBcDivergence( LevelData<FArrayBox>& a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit] -= m_bc_divergence[dit];
   }
}


void
EllipticOp::computeRadialFSAverage( const LevelData<FluxBox>&  a_in,
                                    double&                    a_lo_value,
                                    double&                    a_hi_value,
                                    LevelData<FArrayBox>&      a_out ) const
{
   CH_assert(a_in.nComp() == 1);
   CH_assert(a_in.ghostVect() == IntVect::Zero);
   CH_assert(a_out.ghostVect() == IntVect::Zero);

   const DisjointBoxLayout& grids = m_geometry.grids();
   const MagCoordSys& coords = *m_geometry.getCoordSys();
   
   LevelData<FluxBox> mapped_flux_even(grids, SpaceDim, IntVect::Unit);
   LevelData<FluxBox> mapped_flux_odd(grids, SpaceDim, IntVect::Unit);
    
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
            
         const FArrayBox& this_flux_input_dir = a_in[dit][dir];
         FArrayBox& this_flux_even_dir = mapped_flux_even[dit][dir];
         FArrayBox& this_flux_odd_dir = mapped_flux_odd[dit][dir];
            
         Box box( this_flux_input_dir.box() );
            
         this_flux_even_dir.setVal(0.);
         this_flux_odd_dir.setVal(0.);
            
         int block_number = coords.whichBlock(grids[dit]);
            
         if ( dir == RADIAL_DIR ) {
            
            if ((typeid(coords) != typeid(SingleNullCoordSys)) ||
               ((const SingleNullCoordSys&)coords).isCORE(block_number))  {
               
               BoxIterator bit(box);
               for (bit.begin(); bit.ok(); ++bit) {
                  IntVect iv = bit();
                  if (iv[0]%2 == 0 ) {
                     this_flux_even_dir(iv,0) = this_flux_input_dir(iv,0);
                  }
                  else {
                     this_flux_odd_dir(iv,0) = this_flux_input_dir(iv,0);
                  }
               }
            }
         }
      }
   }
   mapped_flux_even.exchange();
   mapped_flux_odd.exchange();
    
   LevelData<FluxBox> flux_even(grids, SpaceDim, IntVect::Unit);
   LevelData<FluxBox> flux_odd(grids, SpaceDim, IntVect::Unit);

   if (SpaceDim == 2) {
      m_geometry.unmapPoloidalGradient(mapped_flux_even, flux_even);
      m_geometry.unmapPoloidalGradient(mapped_flux_odd, flux_odd);
   }

   else {
      m_geometry.unmap3DGradient(mapped_flux_even, flux_even);
      m_geometry.unmap3DGradient(mapped_flux_odd, flux_odd);
   }

   // Multiply the flux by the unmapped, face-centered coefficients
   multiplyCoefficients(flux_even, false);
   multiplyCoefficients(flux_odd, false);
   
   m_geometry.applyAxisymmetricCorrection(flux_even);
   m_geometry.applyAxisymmetricCorrection(flux_odd);
 
   // At this point, we have face-centered fluxes including one
   // ghost cell layer.

   if (!m_second_order) {
     // Convert to face-averaged
     fourthOrderAverage(flux_even);
     fourthOrderAverage(flux_odd);
   }
    
   // Now we have face-averaged fluxes on valid faces.  In one
   // layer of ghost cells, we still have the face-centered values,
   // which approximate face-averaged values to second-order.

   LevelData<FArrayBox> divergence_even(grids, 1, IntVect::Zero);
   m_geometry.computeMappedGridDivergence(flux_even, divergence_even, !m_second_order);
    
   LevelData<FArrayBox> divergence_odd(grids, 1, IntVect::Zero);
   m_geometry.computeMappedGridDivergence(flux_odd, divergence_odd, !m_second_order);
    
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      divergence_even[dit] /= m_volume[dit];
      divergence_odd[dit] /= m_volume[dit];
   }

   FluxSurface fs(m_geometry);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_out[dit].copy(divergence_even[dit]);  // using a_out as a temp here
   }
   fs.averageAndSpread(a_out, divergence_even);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_out[dit].copy(divergence_odd[dit]);  // using a_out as a temp here
   }
   fs.averageAndSpread(a_out, divergence_odd);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (BoxIterator bit(grids[dit]); bit.ok(); ++bit) {
         IntVect iv = bit();
         if (iv[0]%2 == 0 ) {
            a_out[dit](iv,0) = -divergence_even[dit](iv,0);
         }
         else { 
            a_out[dit](iv,0) = -divergence_odd[dit](iv,0);
         }
      }
   }

   const ProblemDomain& domain0 = m_geometry.getBlockCoordSys(0).domain();
   const Box& domain_box = domain0.domainBox();

   a_lo_value = a_hi_value = -DBL_MAX;
   const MagCoordSys& coord_sys( *(m_geometry.getCoordSys()) );
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {

      const Box& box = grids[dit];
      int block_number = coord_sys.whichBlock(box);

      if ((typeid(coords) != typeid(SingleNullCoordSys)) ||
          ((const SingleNullCoordSys&)coords).isCORE(block_number))  {
         
         IntVect iv = box.smallEnd();
         if ( iv[RADIAL_DIR] == domain_box.smallEnd(RADIAL_DIR) ) {
            if (iv[0]%2 == 0 ) {
               a_lo_value = -divergence_even[dit](iv,0);
            }
            else { 
               a_lo_value = -divergence_odd[dit](iv,0);
            }
         }
         iv = box.bigEnd();
         if ( iv[RADIAL_DIR] == domain_box.bigEnd(RADIAL_DIR) ) {
            if (iv[0]%2 == 0 ) {
               a_hi_value = divergence_odd[dit](iv,0);
            }
            else { 
               a_hi_value = divergence_even[dit](iv,0);
            }
         }
      }
   }
    
   a_lo_value = globalMax(a_lo_value);
   a_hi_value = globalMax(a_hi_value);
}


void
EllipticOp::fillInternalGhosts( LevelData<FArrayBox>& a_phi ) const
{
  // Because mulliBlockExchange objects may have diferent (e.g. fewer) number of ghosts 
  // this branch is called separately from its counterpart in m_geometry.fillInternalGhosts()
  if (m_mblex_potential_Ptr && !m_geometry.extrablockExchange() && !m_geometry.shearedMBGeom()) {

    const IntVect& nghost = a_phi.ghostVect();    

    if (nghost < m_num_potential_ghosts * IntVect::Unit) {
      // interpGhosts() can't seem to handle a smaller number of ghost cells than
      // where specified when creating the MBLevelExchange object.
      const DisjointBoxLayout& grids = m_geometry.grids();
      LevelData<FArrayBox> tmp(grids, a_phi.nComp(), m_num_potential_ghosts * IntVect::Unit);
      for (DataIterator dit(a_phi.dataIterator()); dit.ok(); ++dit) {
	tmp[dit].copy(a_phi[dit]);
      }
      
      m_mblex_potential_Ptr->interpGhosts(tmp);
      for (DataIterator dit(a_phi.dataIterator()); dit.ok(); ++dit) {
	a_phi[dit].copy(tmp[dit]);
      }
    }
    else {
      m_mblex_potential_Ptr->interpGhosts(a_phi);
    }
    a_phi.exchange();
  }

  else {
    m_geometry.fillInternalGhosts(a_phi);
  }
}


void
EllipticOp::extractNormalComponent( const LevelData<FluxBox>&  a_in,
                                    const int                  a_dir,
                                    LevelData<FluxBox>&        a_out ) const
{
   CH_assert(a_out.nComp() == 1);
   CH_assert(a_in.nComp() >= a_dir);

   for (DataIterator dit(a_out.dataIterator()); dit.ok(); ++dit) {
      FluxBox& this_fb = a_out[dit];
      this_fb.setVal(0.);
      this_fb[a_dir].copy(a_in[dit][a_dir],a_dir,0,1);
   }
}


// ----> Begin LinOp virtuals


void
EllipticOp::preCond( LevelData<FArrayBox>&        a_cor,
                     const LevelData<FArrayBox>&  a_residual )
{
   CH_TIME("EllipticOp::preCond");
   setToZero(a_cor);
   solvePreconditioner(a_residual, a_cor);
}


void
EllipticOp::applyOp( LevelData<FArrayBox>&        a_out,
                     const LevelData<FArrayBox>&  a_in,
                     bool                         a_homogeneous )
{
   CH_TIME("EllipticOp::applyOp");
   // We've already accounted for inhomogeneous boundary values, so if
   // a_homogeneous is ever passed in as false we need to know why.
   CH_assert( a_homogeneous );

   computeFluxDivergence(a_in, a_out, a_homogeneous);
}


void
EllipticOp::residual(  LevelData<FArrayBox>&        a_lhs,
                       const LevelData<FArrayBox>&  a_phi,
                       const LevelData<FArrayBox>&  a_rhs,
                       bool                         a_homogeneous)
{
   // We've already accounted for inhomogeneous boundary values, so if
   // a_homogeneous is ever passed in as false we need to know why.
   CH_assert( a_homogeneous );

   applyOp(a_lhs, a_phi, a_homogeneous );

   for (DataIterator dit(a_lhs.dataIterator()); dit.ok(); ++dit) {
      a_lhs[dit].negate();
      a_lhs[dit].plus(a_rhs[dit]);
   }
}


void
EllipticOp::create( LevelData<FArrayBox>&        a_lhs,
                    const LevelData<FArrayBox>&  a_rhs )
{
   a_lhs.define(a_rhs.disjointBoxLayout(), a_rhs.nComp(), a_rhs.ghostVect());
}


void
EllipticOp::assign( LevelData<FArrayBox>&        a_lhs,
                    const LevelData<FArrayBox>&  a_rhs )
{
   for (DataIterator dit(a_lhs.dataIterator()); dit.ok(); ++dit) {
      a_lhs[dit].copy(a_rhs[dit]);
   }
}


Real
EllipticOp::dotProduct( const LevelData<FArrayBox>&  a_1,
                        const LevelData<FArrayBox>&  a_2 )
{
   CH_TIME("EllipticOp::dotProduct");
   const DisjointBoxLayout& grids = a_1.disjointBoxLayout();

   double local_sum = 0.;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      local_sum += a_1[dit].dotProduct(a_2[dit],grids[dit]);
   }

   double global_sum;
#ifdef CH_MPI
   MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
   global_sum = local_sum;
#endif

   return global_sum;
}


void
EllipticOp::incr( LevelData<FArrayBox>&        a_lhs,
                  const LevelData<FArrayBox>&  a_x,
                  Real                         a_scale )
{
   for (DataIterator dit(a_lhs.dataIterator()); dit.ok(); ++dit) {
      a_lhs[dit].plus(a_x[dit], a_scale);
   }
}


void
EllipticOp::axby( LevelData<FArrayBox>&        a_lhs,
                  const LevelData<FArrayBox>&  a_x,
                  const LevelData<FArrayBox>&  a_y,
                  Real                         a_a,
                  Real                         a_b )
{
   for (DataIterator dit(a_lhs.dataIterator()); dit.ok(); ++dit) {
      FArrayBox& this_fab = a_lhs[dit];
      this_fab.copy(a_x[dit]);
      this_fab *= a_a;
      this_fab.plus(a_y[dit], a_b);
   }
}


void
EllipticOp::scale( LevelData<FArrayBox>&  a_lhs,
                   const Real&            a_scale )
{
   for (DataIterator dit(a_lhs.dataIterator()); dit.ok(); ++dit) {
      a_lhs[dit] *= a_scale;
   }
}


Real
EllipticOp::norm( const LevelData<FArrayBox>&  a_rhs,
                  int                          a_ord )
{
   CH_TIME("EllipticOp::norm");
   return CH_XD::norm(a_rhs, a_rhs.interval(), a_ord);
}


void
EllipticOp::setToZero( LevelData<FArrayBox>& a_lhs )
{
   for (DataIterator dit(a_lhs.dataIterator()); dit.ok(); ++dit) {
      a_lhs[dit].setVal(0.);
   }
}


// <---- End LinOp virtuals


void
EllipticOp::computeMapped3DFieldWithGhosts( const LevelData<FArrayBox>&  a_phi,
                                            LevelData<FArrayBox>&        a_field,
                                            const int                    a_order ) const
{
   CH_assert(a_phi.ghostVect() >= 2*IntVect::Unit);
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);
   CH_assert(a_order == 2 || a_order == 4);

   // The following assumes that we have potential values in at least two layers of ghost cells

   // Compute the field to second_order including one layer of ghost cells

   const DisjointBoxLayout& grids = m_geometry.grids();

   int tmp_order = 2;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();
      Box box = grow(grids[dit],1);

      // Initializing to NaN to ensure that the field is actually computed everywhere it's used
      a_field[dit].setVal(1./0.);

      for (int dir=0; dir<SpaceDim; dir++) {
         SpaceUtils::cellCenteredGradientComponent( box,
                                                    dir,
                                                    a_phi[dit],
                                                    dx,
                                                    tmp_order,
                                                    a_field[dit] );
      }
      a_field[dit].negate();

   }
   a_field.exchange();

   // If fourth-order, recompute the field at valid cell centers to fourth-order

   if ( a_order == 4 ) {

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         RealVect dx = block_coord_sys.dx();

         for (int dir=0; dir<SpaceDim; ++dir) {
            SpaceUtils::cellCenteredGradientComponent(  grids[dit],
                                                        dir,
                                                        a_phi[dit],
                                                        dx,
                                                        a_order,
                                                        a_field[dit] );
         }
         a_field[dit].negate(grids[dit], 0, SpaceDim);
      }

      a_field.exchange();
   }
}


void
EllipticOp::computeMappedPoloidalFieldWithGhosts( const LevelData<FArrayBox>&  a_phi,
                                                  LevelData<FArrayBox>&        a_field,
                                                  const int                    a_order ) const
{
   CH_assert(a_phi.ghostVect() >= 2*IntVect::Unit);
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);
   CH_assert(a_order == 2 || a_order == 4);

   // The following assumes that we have potential values in at least two layers of ghost cells

   // Compute the field to second_order including one layer of ghost cells

   const DisjointBoxLayout& grids = m_geometry.grids();

   int tmp_order = 2;

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();
      Box box = grow(grids[dit],1);

      // Initializing to NaN to ensure that the field is actually computed everywhere it's used
      a_field[dit].setVal(1./0.);

      for (int dir=0; dir<SpaceDim; dir++) {
         SpaceUtils::cellCenteredGradientComponent( box,
                                                    dir,
                                                    a_phi[dit],
                                                    dx,
                                                    tmp_order,
                                                    a_field[dit] );
      }
      a_field[dit].negate();
   }

   a_field.exchange();

   // If fourth-order, recompute the field at valid cell centers to fourth-order

   if ( a_order == 4 ) {

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         RealVect dx = block_coord_sys.dx();

         for (int dir=0; dir<SpaceDim; ++dir) {
            SpaceUtils::cellCenteredGradientComponent(  grids[dit],
                                                        dir,
                                                        a_phi[dit],
                                                        dx,
                                                        a_order,
                                                        a_field[dit]);
         }
         a_field[dit].negate( grids[dit], 0, SpaceDim );
      }

      a_field.exchange();
   }
}


void
EllipticOp::computeMapped3DFieldWithGhosts( const LevelData<FArrayBox>&  a_phi,
                                            LevelData<FluxBox>&          a_field,
                                            const int                    a_order ) const
{
   CH_assert(a_phi.ghostVect() >= 2*IntVect::Unit);
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);
   CH_assert(a_order == 2 || a_order == 4);

   // The following assumes that we have potential values in at least two layers of ghost cells at
   // all block boundaries.

   const DisjointBoxLayout& grids = m_geometry.grids();
   LevelData<FluxBox> phi_face(grids, 1, 2*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         IntVect grow_vect = 2*IntVect::Unit;
         grow_vect[dir] = 0;
         Box box = grow(grids[dit],grow_vect);

         SpaceUtils::faceInterpolate(dir,
                                     surroundingNodes(box,dir),
                                     a_order,
                                     a_phi[dit],
                                     phi_face[dit][dir] );
      }
   }
   phi_face.exchange();

   // Compute the field to second-order including one layer of transverse faces

   int tmp_order = 2;

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();
      Box box = grow(grids[dit],1);

      // Initializing to NaN to ensure that the field is actually computed everywhere it's used
      a_field[dit].setVal(1./0.);

      for (int dir=0; dir<SpaceDim; dir++) {
         Box box_dir = surroundingNodes(box,dir);
         SpaceUtils::faceCenteredGradientComponent( box_dir,
                                                    dir,
                                                    a_phi[dit],
                                                    dx,
                                                    tmp_order,
                                                    a_field[dit][dir] );

         for (int tdir=0; tdir<SpaceDim; ++tdir) {
            if (tdir != dir) {
               SpaceUtils::cellCenteredGradientComponent( box_dir,
                                                          tdir,
                                                          phi_face[dit][dir],
                                                          dx,
                                                          tmp_order,
                                                          a_field[dit][dir] );
            }
         }
         a_field[dit][dir].negate();
      }
   }
   a_field.exchange();

   // Second-order extrapolate the field to the transverse physical boundary ghosts
   m_geometry.fillTransversePhysicalGhosts(a_field);

   // If fourth-order, recompute the field on valid cell faces

   if ( a_order == 4 ) {

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         RealVect dx = block_coord_sys.dx();


         for (int dir=0; dir<SpaceDim; dir++) {
            Box box_dir = surroundingNodes(grids[dit],dir);
            SpaceUtils::faceCenteredGradientComponent(  box_dir,
                                                        dir,
                                                        a_phi[dit],
                                                        dx,
                                                        a_order,
                                                        a_field[dit][dir] );

            a_field[dit][dir].negate(box_dir, dir);
            for (int tdir=0; tdir<SpaceDim; ++tdir) {
               if (tdir != dir) {
                  SpaceUtils::cellCenteredGradientComponent(  box_dir,
                                                              tdir,
                                                              phi_face[dit][dir],
                                                              dx,
                                                              a_order,
                                                              a_field[dit][dir] );
                  a_field[dit][dir].negate(box_dir, tdir);
               }
            }
         }
      }

      a_field.exchange();
   }

}


void
EllipticOp::computeMappedPoloidalFieldWithGhosts( const LevelData<FArrayBox>&  a_phi,
                                                  LevelData<FluxBox>&          a_field,
                                                  const int                    a_order ) const
{
   CH_assert(a_phi.ghostVect() >= 2*IntVect::Unit);
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);
   CH_assert(a_order == 2 || a_order == 4);

   // The following assumes that we have potential values in at least two layers of ghost cells at
   // all block boundaries.

   const DisjointBoxLayout& grids = m_geometry.grids();
   LevelData<FluxBox> phi_face(grids, 1, 2*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         IntVect grow_vect = 2*IntVect::Unit;
         grow_vect[dir] = 0;
         Box box = grow(grids[dit],grow_vect);

         SpaceUtils::faceInterpolate( dir,
                                      surroundingNodes(box,dir),
                                      a_order,
                                      a_phi[dit],
                                      phi_face[dit][dir] );
      }
   }
   phi_face.exchange();

   // Compute the field to second-order including one layer of transverse faces

   int tmp_order = 2;

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();
      Box box = grow(grids[dit],1);

      // Initializing to NaN to ensure that the field is actually computed everywhere it's used
      a_field[dit].setVal(1./0.);

      for (int dir=0; dir<SpaceDim; dir++) {
         Box box_dir = surroundingNodes(box,dir);
         SpaceUtils::faceCenteredGradientComponent( box_dir,
                                                    dir,
                                                    a_phi[dit],
                                                    dx,
                                                    tmp_order,
                                                    a_field[dit][dir] );

         for (int tdir=0; tdir<SpaceDim; ++tdir) {
            if (tdir != dir) {
               SpaceUtils::cellCenteredGradientComponent( box_dir,
                                                          tdir,
                                                          phi_face[dit][dir],
                                                          dx,
                                                          tmp_order,
                                                          a_field[dit][dir] );
            }
         }
         a_field[dit][dir].negate();
      }
   }
   a_field.exchange();

   // Second-order extrapolate the field to the transverse physical boundary ghosts
   m_geometry.fillTransversePhysicalGhosts(a_field);

   // If fourth-order, recompute the field on valid cell faces

   if ( a_order == 4 ) {

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         RealVect dx = block_coord_sys.dx();

         for (int dir=0; dir<SpaceDim; dir++) {
            Box box_dir = surroundingNodes(grids[dit],dir);
            SpaceUtils::faceCenteredGradientComponent(  box_dir,
                                                        dir,
                                                        a_phi[dit],
                                                        dx,
                                                        a_order,
                                                        a_field[dit][dir] );

            a_field[dit][dir].negate(box_dir, dir);
            for (int tdir=0; tdir<SpaceDim; ++tdir) {
               if (tdir != dir) {
                  SpaceUtils::cellCenteredGradientComponent(  box_dir,
                                                              tdir,
                                                              phi_face[dit][dir],
                                                              dx,
                                                              a_order,
                                                              a_field[dit][dir] );
                  a_field[dit][dir].negate(box_dir, tdir);
               }
            }
         }
      }

      a_field.exchange();
   }
}


void
EllipticOp::constructBoundaryStencils( const bool                         a_fourth_order,
                                       const EllipticOpBC&                a_bc,
                                       Vector< Vector<CoDim1Stencil> >&   a_codim1_stencils,
                                       Vector< Vector<CoDim2Stencil> >&   a_codim2_stencils ) const
{
   CH_assert(a_codim1_stencils.size() == 0);
   CH_assert(a_codim2_stencils.size() == 0);

   int order = a_fourth_order? 4: 2;

   const RefCountedPtr<MagCoordSys>& coord_sys_ptr = m_geometry.getCoordSys();

   int num_blocks = coord_sys_ptr->numBlocks();
   a_codim1_stencils.resize(num_blocks);
   a_codim2_stencils.resize(num_blocks);

   const Vector<Box>& mapping_blocks = coord_sys_ptr->mappingBlocks();

   for (int block_number=0; block_number<num_blocks; ++block_number) {
      const NewCoordSys* block_coord_sys = coord_sys_ptr->getCoordSys(block_number);
      Box domain_box = mapping_blocks[block_number];
      RealVect dx = block_coord_sys->dx();

      Vector<CoDim1Stencil>& codim1_stencil = a_codim1_stencils[block_number];      
      Vector<CoDim2Stencil>& codim2_stencil = a_codim2_stencils[block_number];      

      int num_codim1_neighbors = 0;
      int num_codim2_neighbors = 0;

      for (int dir=0; dir<SpaceDim; ++dir) {
         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();

            Vector<Box> boundary_boxes = ((MagGeom&)m_geometry).getBoundaryBoxes(block_number, dir, side);

            for ( int n=0; n<boundary_boxes.size(); ++n ) {

               codim1_stencil.resize(++num_codim1_neighbors);

               Box box = boundary_boxes[n];
               FArrayBox bv(box,1);

	       Box box_tmp = box;
               box_tmp.shiftHalf(dir,-sign(side));
               FluxBox bv_tmp(box_tmp,1);
               
               RefCountedPtr<GridFunction> bc_func = a_bc.getBCFunction(block_number, dir, side );
               if (bc_func && !(typeid(*bc_func) == typeid(DataArray))) {
 		  FluxBox real_coords(box_tmp,CFG_DIM);
		  FluxBox norm_flux(box_tmp,1);

		  const MagBlockCoordSys& coord_sys( m_geometry.getBlockCoordSys( block_number ) );
		  
		  coord_sys.getFaceCenteredRealCoords(dir, real_coords[dir]);
		  coord_sys.getNormMagneticFlux(real_coords[dir], norm_flux[dir]);
		  
                  bc_func->assign(bv_tmp, m_geometry, real_coords, norm_flux, block_number, 0., false);

                  for (BoxIterator bit(box); bit.ok(); ++bit) {
                     IntVect iv = bit();
                     iv[dir] = (bv_tmp[dir].box()).sideEnd(side)[dir];
                     bv(bit(),0) = bv_tmp[dir](iv,0);
                  }
               }

	       else if (bc_func && (typeid(*bc_func) == typeid(DataArray))) {
                 FArrayBox dummy;
                 bc_func->assign(bv, m_geometry, dummy, dummy, block_number, 0., false);
	       }
	       
               else {
                  for (BoxIterator bit(box); bit.ok(); ++bit) {
                     bv(bit(),0) = a_bc.getBCValue(block_number, dir, side);
                  }
               }

               int bc_type = a_bc.getBCType(block_number, dir, side);

               codim1_stencil[num_codim1_neighbors-1].define(bv, bc_type, dx, dir, side, order);

               Box codim1_box = codim1_stencil[num_codim1_neighbors-1].box();

               CH_assert(order%2 == 0);
               int width = order/2;

               for (int tdir=0; tdir < SpaceDim; ++tdir) {
                  if ( tdir != dir ) {
                     for (SideIterator sit2; sit2.ok(); ++sit2) {
                        Side::LoHiSide side2 = sit2();

                        // Determine if the codim2 stencil we are about to create overlaps a physical boundary in the tdir direction on side2

                        Box grown_codim1_box = codim1_box;
                        grown_codim1_box.grow(tdir,width);

                        Vector<Box> transverse_boundary_boxes = ((MagGeom&)m_geometry).getBoundaryBoxes(block_number, tdir, side2);

                        bool transverse_boundary = transverse_boundary_boxes.size() > 0;

                        if ( transverse_boundary ) {
                           // The transverse block boundary contains some part of the physical boundary,
                           // indicated by the transverse_boundary_boxes vector.

                           for (int n=0; n<transverse_boundary_boxes.size(); ++n) {

                              Box restricted_codim1_box = codim1_box;
                              bool this_transverse_boundary = false;

                              Box tbox = transverse_boundary_boxes[n];
                              tbox.shiftHalf(tdir,sign(side2));
                              tbox.grow(dir,width);
                              Box overlap = tbox & grown_codim1_box;
                              if ( overlap.ok() ) {
                                 for (int tdir2=0; tdir2<SpaceDim; ++tdir2) {
                                    if ( tdir2 != tdir ) {
                                       restricted_codim1_box.setSmall(tdir2, overlap.smallEnd(tdir2));
                                       restricted_codim1_box.setBig(tdir2, overlap.bigEnd(tdir2));
                                    }
                                 }
                                 this_transverse_boundary = true;
                              }

                              codim2_stencil.resize(++num_codim2_neighbors);

                              codim2_stencil[num_codim2_neighbors-1].
                                 define(restricted_codim1_box, dx, dir, side, tdir, side2, order, this_transverse_boundary);

                              // Check to see if the new codim2 box overlaps another block

                              const Box& new_codim2_box = codim2_stencil[num_codim2_neighbors-1].box();

                              for (int block_number2=0; block_number2<num_blocks; ++block_number2) {
                                 if (block_number2 != block_number) {
                                    const Box& domain_box2 = mapping_blocks[block_number2];

                                    Box overlap = new_codim2_box & domain_box2;
                                    if ( overlap == new_codim2_box ) {
                                       // Codim2 box is contained in a valid block, so delete it

                                       codim2_stencil.resize(--num_codim2_neighbors);
                                    }
                                    else if ( overlap.ok() ) {
                                       // Codim2 box partially overlaps a valid block.  Need to exit and figure out
                                       // what to do about it.
                                       MayDay::Error("MBSolver::constructBoundaryStencils(): codim2 box partially overlaps a valid block");
                                    }
                                 }
                              }
                           }
                        }
                        else {
                           // The transverse block boundary contains no part of the physical boundary,

                           codim2_stencil.resize(++num_codim2_neighbors);

                           codim2_stencil[num_codim2_neighbors-1].
                              define(codim1_box, dx, dir, side, tdir, side2, order, transverse_boundary);

                           // Check to see if the new codim2 box overlaps another block

                           const Box& new_codim2_box = codim2_stencil[num_codim2_neighbors-1].box();

                           for (int block_number2=0; block_number2<num_blocks; ++block_number2) {
                              if (block_number2 != block_number) {
                                 const Box& domain_box2 = mapping_blocks[block_number2];

                                 Box overlap = new_codim2_box & domain_box2;
                                 if ( overlap == new_codim2_box ) {
                                    // Codim2 box is contained in a valid block, so delete it

                                    codim2_stencil.resize(--num_codim2_neighbors);
                                 }
                                 else if ( overlap.ok() ) {
                                    // Codim2 box partially overlaps a valid block.  Need to exit and figure out
                                    // what to do about it.
                                    MayDay::Error("MBSolver::constructBoundaryStencils(): codim2 box partially overlaps a valid block");
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}


void
EllipticOp::accumPhysicalGhosts( const Vector< Vector<CoDim1Stencil> >&   a_codim1_stencils,
                                 const Vector< Vector<CoDim2Stencil> >&   a_codim2_stencils,
                                 const bool                               a_extrapolate_from_interior,
                                 const bool                               a_include_bvs,
                                 LevelData<FArrayBox>&                    a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   const MultiBlockCoordSys* coord_sys_ptr = m_geometry.coordSysPtr();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coord_sys_ptr->whichBlock(grids[dit]);
      FArrayBox& this_data = a_data[dit];
      const Box& this_data_box = this_data.box();

      // Fill the codim1 ghosts first, which are assumed to depend only upon valid data
      // and boundary values

      const Vector<CoDim1Stencil>& codim1_stencil = a_codim1_stencils[block_number];      

      for (int i=0; i<codim1_stencil.size(); ++i) {
         const CoDim1Stencil& bndry_stencil = codim1_stencil[i];
         const Box overlap = bndry_stencil.box() & this_data_box;

         for (BoxIterator bit(overlap); bit.ok(); ++bit) {
            IntVect iv = bit();

            vector<IntVect> points;
            vector<double> weights;
            double bv_contrib;

            bndry_stencil.getStencil(iv, points, weights, bv_contrib);
                  
            if (a_include_bvs) {
               this_data(iv,0) += bv_contrib;
            }
            if (a_extrapolate_from_interior) {
               for (int j=0; j<points.size(); ++j) {
                  this_data(iv,0) += weights[j] * this_data(points[j],0);
               }
            }
         }
      }

      // Fill the codim2 ghosts next, which are assumed to depend upon valid data
      // and/or codim1 data

      const Vector<CoDim2Stencil>& codim2_stencil = a_codim2_stencils[block_number];      

      for (int i=0; i<codim2_stencil.size(); ++i) {
         const CoDim2Stencil& bndry_stencil = codim2_stencil[i];
         const Box overlap = bndry_stencil.box() & this_data_box;

         for (BoxIterator bit(overlap); bit.ok(); ++bit) {
            IntVect iv = bit();

            vector<IntVect> points;
            vector<double> weights;

            bndry_stencil.getStencil(iv, points, weights);

            for (int j=0; j<points.size(); ++j) {
               this_data(iv,0) += weights[j] * this_data(points[j],0);
            }
         }
      }
   }
}


void
EllipticOp::interpToNodes( const LevelData<FArrayBox>&  a_phi,
                           LevelData<FArrayBox>&        a_phi_node ) const
{
#if CFG_DIM==3
   MayDay::Error("EllipticOp::interpToNodes() not yet implemented");
#else
   
   CH_assert(a_phi_node.ghostVect() >= 2*IntVect::Unit);

   // Make a temporary with ghost cells and copy the potential on valid cells

   const DisjointBoxLayout & grids = m_geometry.grids();
   CH_assert(m_num_potential_ghosts >= 2);
   LevelData<FArrayBox> phi_cell(grids, 1, m_num_potential_ghosts*IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi_cell[dit].setVal(0.);
      phi_cell[dit].copy(a_phi[dit],grids[dit]);
   }

   fillInternalGhosts(phi_cell);

#ifdef INTERP_TO_NODES_WITH_BCS
   
   CH_assert(m_second_order);
   
   // Fill the physical boundary ghost cells
   bool extrapolate_from_interior = true;
   bool include_bvs = true;
   accumPhysicalGhosts( m_codim1_stencils, m_codim2_stencils, extrapolate_from_interior, include_bvs, phi_cell );

#else

   // If fourth-order, extrapolate the potential to three cells at physical boundaries
   // If second-order, extrapolate the potential to two cells at physical boundaries

   int order = m_second_order? 2: 4;

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   //const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = mag_coord_sys->boundaries();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = mag_coord_sys->whichBlock(grids[dit]);
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
      const ProblemDomain& domain = block_coord_sys.domain();

      for (int dir = 0; dir < SpaceDim; dir++) {
        for (SideIterator sit; sit.ok(); ++sit) {
          Side::LoHiSide side = sit();

          if (mag_coord_sys->containsPhysicalBoundary(block_number, dir, side)) {
            IntVect grow_vec = phi_cell.ghostVect();
            grow_vec[dir] = 0;
            Box interior_box = grow(grids[dit],grow_vec);
            Box domain_box = grow(domain.domainBox(),grow_vec);

            SpaceUtils::extrapBoundaryGhostsForCC(phi_cell[dit], interior_box, domain_box, dir, order, sign(side));
          }
	}
      }
   }

   if (m_geometry.mixedBoundaries()) {
     fillInternalGhosts(phi_cell);
   }
     
#endif
   
   ParmParse pp;

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      FArrayBox& this_phi_cell = phi_cell[dit];

      const Box& phi_box = this_phi_cell.box();
      FArrayBox interp_cell_coords(phi_box, SpaceDim);
      block_coord_sys.MagBlockCoordSys::getCellCenteredMappedCoords(interp_cell_coords);

      Interp* interp_to_nodes;
      if ( order == 2 ) {
         interp_to_nodes = new BilinearInterp(pp, interp_cell_coords, this_phi_cell);
      }
      else {
         interp_to_nodes = new BicubicInterp(pp, interp_cell_coords, this_phi_cell);
      }

      Box node_box = grow(surroundingNodes(grids[dit]), 1);
      FArrayBox& this_phi_node = a_phi_node[dit];
      FArrayBox xi_node(node_box,SpaceDim);
      block_coord_sys.MagBlockCoordSys::getNodeCenteredMappedCoords(xi_node);

      interp_to_nodes->interpolate(0, 0, 0, xi_node, node_box, this_phi_node);
      delete interp_to_nodes;
   }

   a_phi_node.exchange();
#endif
}


double
EllipticOp::L2Norm( const LevelData<FArrayBox>& a_data )
{
   return sqrt(dotProduct(a_data, a_data));
}


double
EllipticOp::globalMax(const double a_data) const
{
   double global_max;

#ifdef CH_MPI
   double local_data = a_data;
   MPI_Allreduce(&local_data, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
   global_max = a_data;
#endif

   return global_max;
}


#include "NamespaceFooter.H"
