#include "FieldSolver.H"
#include "FieldSolverF_F.H"
#include "BlockRegister.H"
#include "EdgeToCell.H"
#include "FourthOrderUtil.H"
#include "Directions.H"
#include "FluxSurface.H"
#include "BiCGStabSolver.H"
#include "GMRESSolver.H"
#include "MBSolverF_F.H"


#include "NamespaceHeader.H"


FieldSolver::FieldSolver( const ParmParse& a_pp,
                          const MagGeom&   a_geom )
   : m_geometry(a_geom),
     m_num_potential_ghosts(3)
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

   if (a_pp.contains("second_order")) {
     a_pp.get("second_order", m_second_order);
   }
   else {
     m_second_order = false;
   }

   ParmParse pp_linear_solver( ((string)a_pp.prefix() + ".linear_solver").c_str());
   defineLinearSolver( pp_linear_solver );

   // If there is more than one block, construct the multiblock exchange object
   if ( m_geometry.coordSysPtr()->numBlocks() > 1 ) {
     m_mblex_potential_Ptr = new MultiBlockLevelExchangeCenter();
     int spaceOrder = 4;
     m_mblex_potential_Ptr->define(&a_geom, m_num_potential_ghosts, spaceOrder);
   }
   else {
     m_mblex_potential_Ptr = NULL;
   }

   m_volume.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume);
   m_bc_divergence.define(grids, 1, IntVect::Zero);
}
      


FieldSolver::~FieldSolver()
{
   if (m_Chombo_solver) delete m_Chombo_solver;
   if (m_mblex_potential_Ptr) delete m_mblex_potential_Ptr;
}



void
FieldSolver::defineLinearSolver( const ParmParse& a_pp )
{
   // Defaults; over-ridden if present in ParmParse object
   m_method   = "BiCGStab";
   m_tol      = 1.e-6;
   m_max_iter = 20;
   m_verbose  = false;
   parseMethodAndParams(a_pp, m_method, m_tol, m_max_iter, m_verbose);

   ParmParse pp_precond( ((string)a_pp.prefix() + ".precond").c_str());

   // Defaults; over-ridden if present in ParmParse object
   m_precond_method   = "AMG";
   m_precond_tol      = 0.;
   m_precond_max_iter = 1;
   m_precond_verbose  = false;
   parseMethodAndParams(pp_precond, m_precond_method, m_precond_tol,
                        m_precond_max_iter, m_precond_verbose);

   ParmParse pp_precond_precond( ((string)pp_precond.prefix() + ".precond").c_str());

   // Defaults; over-ridden if present in ParmParse object
   m_precond_precond_method   = "";
   m_precond_precond_tol      = 0.;
   m_precond_precond_max_iter = 0;
   m_precond_precond_verbose  = false;
   parseMethodAndParams(pp_precond_precond, m_precond_precond_method, m_precond_precond_tol,
                        m_precond_precond_max_iter, m_precond_precond_verbose);

   if ( m_method == "BiCGStab" ) {
      m_Chombo_solver = new BiCGStabSolver< LevelData<FArrayBox> >;
   }
   else {
      m_Chombo_solver = new GMRESSolver< LevelData<FArrayBox> >;
   }

   m_Chombo_solver->define(this, true);
}



void
FieldSolver::parseMethodAndParams( const ParmParse& a_pp,
                                   string&          a_method,
                                   double&          a_tol,
                                   int&             a_max_iter,
                                   bool&            a_verbose ) const
{
   if (a_pp.contains("method")) {
      a_pp.get("method", a_method);
   }

   if (a_pp.contains("tol")) {
      a_pp.get("tol", a_tol);
   }

   if (a_pp.contains("max_iter")) {
      a_pp.get("max_iter", a_max_iter);
   }

   if (a_pp.contains("verbose")) {
      a_pp.get("verbose", a_verbose);
   }
}



void
FieldSolver::computePotential( LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_charge_density )
{
   LevelData<FArrayBox> rhs;
   rhs.define(a_charge_density);

   // Add any inhomogeneous boundary values to the right-hand side
   subtractBcDivergence(rhs);

   solve( rhs, a_phi );
}



void
FieldSolver::solve( const LevelData<FArrayBox>& a_rhs,
                    LevelData<FArrayBox>&       a_solution )
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

   setPreconditionerConvergenceParams(m_precond_tol, m_precond_max_iter,
                                      m_precond_precond_tol, m_precond_precond_max_iter);

   setToZero(a_solution);
   m_Chombo_solver->solve(a_solution, a_rhs);

   if (procID() == 0) {
      if ( m_method == "BiCGStab" ) {
         int exit_status = ((BiCGStabSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_exitStatus;
         if ( exit_status == 1 ) {
            cout << "      BiCGStab converged successfully" << endl;
         }
         else {
            cout << "      BiCGStab solver returned " << exit_status << endl;
         }
      }
      else if ( m_method == "GMRES" ) {
         int exit_status = ((GMRESSolver< LevelData<FArrayBox> >*)m_Chombo_solver)->m_exitStatus;
         if ( exit_status == 1 ) {
            cout << "      GMRES converged successfully" << endl;
         }
         else {
            cout << "      GMRES solver returned " << exit_status << endl;
         }
      }
   }
}



void
FieldSolver::computePoloidalFieldWithBCs( const LevelData<FArrayBox>& a_phi,
                                          LevelData<FluxBox>&         a_field,
                                          const bool                  a_homogeneousBCs ) const
{
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FluxBox> mapped_field(m_geometry.grids(), 2, a_field.ghostVect());
   computeMappedPoloidalFieldWithBCs(a_phi, mapped_field, a_homogeneousBCs);

   unmapPoloidalGradient(mapped_field, a_field);
}



void
FieldSolver::computeMappedPoloidalFieldWithBCs( const LevelData<FArrayBox>& a_phi,
                                                LevelData<FluxBox>&         a_field,
                                                const bool                  a_homogeneousBCs ) const
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
FieldSolver::computePoloidalField( const LevelData<FArrayBox>& a_phi,
                                   LevelData<FArrayBox>&       a_field ) const
{
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FArrayBox> mapped_field(m_geometry.grids(), 2, IntVect::Unit);
   computeMappedPoloidalField(a_phi, mapped_field);

   unmapPoloidalGradient(mapped_field, a_field);
}



void
FieldSolver::computePoloidalField( const LevelData<FArrayBox>& a_phi,
                                   LevelData<FluxBox>&         a_field ) const
{
   CH_assert(a_field.nComp() == 2);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FluxBox> mapped_field(m_geometry.grids(), 2, IntVect::Unit);
   computeMappedPoloidalField(a_phi, mapped_field);

   unmapPoloidalGradient(mapped_field, a_field);
}



void
FieldSolver::computeMappedPoloidalField( const LevelData<FArrayBox>& a_phi,
                                         LevelData<FArrayBox>&       a_field ) const
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

   // Extrapolate the potential to two ghost cells at physical boundaries

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = mag_coord_sys->boundaries();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = mag_coord_sys->whichBlock(grids[dit]);
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
      const ProblemDomain& domain = block_coord_sys.domain();
      for (int dir=0; dir<SpaceDim; ++dir) {
         if ( !domain.isPeriodic(dir) ) {
            IntVect grow_vec = 2*IntVect::Unit;
            grow_vec[dir] = 0;
            Box interior_box = grow(grids[dit],grow_vec);
            Box domain_box = grow(domain.domainBox(),grow_vec);

            extrapBoundaryGhostsForCC(phi_cell[dit], interior_box, domain_box, dir, boundaries[block_number]);
         }
      }
   }

   int order = m_second_order? 2: 4;

   computeMappedPoloidalFieldWithGhosts(phi_cell, a_field, order);
}



void
FieldSolver::computeMappedPoloidalField( const LevelData<FArrayBox>& a_phi,
                                         LevelData<FluxBox>&         a_field ) const
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

   // Fourth-order extrapolate the potential to two cells at physical boundaries

   int tmp_order = 4;

   const MagCoordSys* mag_coord_sys = m_geometry.getCoordSys();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = mag_coord_sys->boundaries();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = mag_coord_sys->whichBlock(grids[dit]);
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(block_number);
      const ProblemDomain& domain = block_coord_sys.domain();
      for (int dir=0; dir<SpaceDim; ++dir) {
         if ( !domain.isPeriodic(dir) ) {
            IntVect grow_vec = phi_cell.ghostVect();
            grow_vec[dir] = 0;
            Box interior_box = grow(grids[dit],grow_vec);
            Box domain_box = grow(domain.domainBox(),grow_vec);

            extrapBoundaryGhostsForFC(phi_cell[dit], interior_box, domain_box, dir, tmp_order, boundaries[block_number]);
         }
      }
   }

   int order = m_second_order? 2: 4;

   computeMappedPoloidalFieldWithGhosts(phi_cell, a_field, order);
}



void
FieldSolver::computeField( const LevelData<FArrayBox>& a_phi,
                           LevelData<FArrayBox>&       a_field ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FArrayBox> mapped_field(a_field.disjointBoxLayout(), 3, a_field.ghostVect());
   computeMappedField( a_phi, mapped_field );

   unmapGradient(mapped_field, a_field);
}



void
FieldSolver::computeField( const LevelData<FArrayBox>& a_phi,
                           LevelData<FluxBox>&         a_field ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FluxBox> mapped_field(a_field.disjointBoxLayout(), 3, a_field.ghostVect());
   computeMappedField( a_phi, mapped_field );

   unmapGradient(mapped_field, a_field);
}



void
FieldSolver::computeMappedField( const LevelData<FArrayBox>& a_phi,
                                 LevelData<FArrayBox>&       a_field ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FArrayBox> poloidal_field(a_field.disjointBoxLayout(), 2, IntVect::Unit);
   computeMappedPoloidalField(a_phi, poloidal_field);

   injectPoloidalVector(poloidal_field, a_field);
}



void
FieldSolver::computeMappedField( const LevelData<FArrayBox>& a_phi,
                                 LevelData<FluxBox>&         a_field ) const
{
   CH_assert(a_field.nComp() == 3);
   CH_assert(a_field.ghostVect() == IntVect::Unit);

   LevelData<FluxBox> poloidal_field(a_field.disjointBoxLayout(), 2, IntVect::Unit);
   computeMappedPoloidalField(a_phi, poloidal_field);

   injectPoloidalVector(poloidal_field, a_field);
}



void
FieldSolver::computeFluxDivergence( const LevelData<FArrayBox>& a_in,
                                    LevelData<FArrayBox>&       a_out,
                                    const bool                  a_homogeneous_bcs )
{
   const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

   LevelData<FArrayBox> phi(grids, 1, 3*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi[dit].copy(a_in[dit]);
   }
   fillInternalGhosts(phi);

   LevelData<FluxBox> flux(grids, 2, IntVect::Unit);
   computePoloidalFieldWithBCs(phi, flux, a_homogeneous_bcs);

   // Multiply the field by the unmapped, face-centered GKP coefficients
   multiplyUnmappedCoefficients(flux);

   m_geometry.fillTransversePhysicalGhosts(flux);

   m_geometry.applyAxisymmetricCorrection(flux);

   // Convert to face-averaged
   if (!m_second_order) fourthOrderAverage(flux);

   m_geometry.computeMappedGridDivergence(flux, a_out, !m_second_order);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_out[dit] /= m_volume[dit];
   }
}



void
FieldSolver::computeBcDivergence( const PotentialBC&    a_bc,
                                  LevelData<FArrayBox>& a_out ) const
{
   const DisjointBoxLayout& grids = a_out.disjointBoxLayout();
   LevelData<FArrayBox> phi(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phi[dit].setVal(0.);
   }

   LevelData<FluxBox> mapped_field(grids, 2, IntVect::Unit);
   computeMappedPoloidalFieldWithBCs( phi, mapped_field, false );
   
   LevelData<FluxBox> flux(grids, 2, IntVect::Unit);
   unmapPoloidalGradient(mapped_field, flux);

   // Multiply the field by the unmapped, face-centered coefficients
   multiplyUnmappedCoefficients(flux);

   m_geometry.applyAxisymmetricCorrection(flux);

   // Convert to face-averaged if we're fourth-order
   if (!m_second_order) fourthOrderAverage(flux);

   m_geometry.computeMappedGridDivergence(flux, a_out, !m_second_order);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_out[dit] /= m_volume[dit];
   }
}



void
FieldSolver::setBc( const PotentialBC& a_bc )
{
   m_codim1_stencils.clear();
   m_codim2_stencils.clear();
   constructBoundaryStencils(!m_second_order, a_bc, m_codim1_stencils, m_codim2_stencils );
}



void
FieldSolver::addBcDivergence( LevelData<FArrayBox>& a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit] += m_bc_divergence[dit];
   }
}



void
FieldSolver::subtractBcDivergence( LevelData<FArrayBox>& a_data ) const
{
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      a_data[dit] -= m_bc_divergence[dit];
   }
}



void
FieldSolver::computeRadialFSAverage( const LevelData<FluxBox>& a_in,
                                     double&                   a_lo_value,
                                     double&                   a_hi_value,
                                     LevelData<FArrayBox>&     a_out ) const
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
                
            if ( block_number < 2 ) {
                    
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
   unmapPoloidalGradient(mapped_flux_even, flux_even);

   LevelData<FluxBox> flux_odd(grids, SpaceDim, IntVect::Unit);
   unmapPoloidalGradient(mapped_flux_odd, flux_odd);
    
   // Multiply the flux by the unmapped, face-centered coefficients
   multiplyUnmappedCoefficients(flux_even);
   multiplyUnmappedCoefficients(flux_odd);
   
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
      if (block_number < 2) {
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
FieldSolver::fillInternalGhosts( LevelData<FArrayBox>& a_phi ) const
{
   if (m_mblex_potential_Ptr) {
      if ( m_geometry.extrablockExchange() ) {
         m_geometry.exchangeExtraBlockGhosts(a_phi);
      }
      else {
         m_mblex_potential_Ptr->interpGhosts(a_phi);
      }
   }

   a_phi.exchange();
}



void
FieldSolver::unmapGradient( const LevelData<FArrayBox>& a_mapped_gradient,
                            LevelData<FArrayBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 3);
   CH_assert(a_gradient.nComp() == 3);
   const DisjointBoxLayout & grids = m_geometry.grids();

   LevelData<FArrayBox> mapped_poloidal_gradient(grids, 2, a_mapped_gradient.ghostVect());
   projectPoloidalVector(a_mapped_gradient, mapped_poloidal_gradient);

   LevelData<FArrayBox> poloidal_gradient(grids, 2, a_gradient.ghostVect());
   unmapPoloidalGradient(mapped_poloidal_gradient, poloidal_gradient);

   injectPoloidalVector(poloidal_gradient, a_gradient);
}



void
FieldSolver::unmapGradient( const LevelData<FluxBox>& a_mapped_gradient,
                            LevelData<FluxBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 3);
   CH_assert(a_gradient.nComp() == 3);
   const DisjointBoxLayout & grids = m_geometry.grids();

   LevelData<FluxBox> mapped_poloidal_gradient(grids, 2, a_mapped_gradient.ghostVect());
   projectPoloidalVector(a_mapped_gradient, mapped_poloidal_gradient);

   LevelData<FluxBox> poloidal_gradient(grids, 2, a_gradient.ghostVect());
   unmapPoloidalGradient(mapped_poloidal_gradient, poloidal_gradient);

   injectPoloidalVector(poloidal_gradient, a_gradient);
}



void
FieldSolver::unmapPoloidalGradient( const LevelData<FArrayBox>& a_mapped_gradient,
                                    LevelData<FArrayBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 2);
   CH_assert(a_gradient.nComp() == 2);
   const DisjointBoxLayout & grids = m_geometry.grids();

   // Multiply by NJInverse
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);

      Box box = a_mapped_gradient[dit].box();
      FArrayBox NJInverse(box, 4);
      block_coord_sys.getPointwiseNJInverse(NJInverse);

      FORT_MULT_NJINVERSE(CHF_BOX(box),
                          CHF_CONST_FRA(a_mapped_gradient[dit]),
                          CHF_CONST_FRA(NJInverse),
                          CHF_FRA(a_gradient[dit]));
   }
}



void
FieldSolver::unmapPoloidalGradient( const LevelData<FluxBox>& a_mapped_gradient,
                                    LevelData<FluxBox>&       a_gradient ) const
{
   CH_assert(a_mapped_gradient.nComp() == 2);
   CH_assert(a_gradient.nComp() == 2);
   const DisjointBoxLayout & grids = m_geometry.grids();

   // Multiply by NJInverse
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);

      FluxBox NJInverse(a_mapped_gradient[dit].box(), 4);
      block_coord_sys.getPointwiseNJInverse(NJInverse);

      for (int dir=0; dir<2; dir++) {
         FORT_MULT_NJINVERSE(CHF_BOX(a_mapped_gradient[dit][dir].box()),
                             CHF_CONST_FRA(a_mapped_gradient[dit][dir]),
                             CHF_CONST_FRA(NJInverse[dir]),
                             CHF_FRA(a_gradient[dit][dir]));
      }
   }
}



void
FieldSolver::extractNormalComponent( const LevelData<FluxBox>& a_in,
                                     const int                 a_dir,
                                     LevelData<FluxBox>&       a_out ) const
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
FieldSolver::preCond( LevelData<FArrayBox>&       a_cor,
                      const LevelData<FArrayBox>& a_residual )
{
   setToZero(a_cor);
   solvePreconditioner(a_residual, a_cor);
}



void
FieldSolver::applyOp( LevelData<FArrayBox>& a_out,
                      const LevelData<FArrayBox>& a_in,
                      bool a_homogeneous )
{
   // We've already accounted for inhomogeneous boundary values, so if
   // a_homogeneous is ever passed in as false we need to know why.
   CH_assert( a_homogeneous );

   computeFluxDivergence(a_in, a_out, a_homogeneous);
}



void
FieldSolver::residual(  LevelData<FArrayBox>&       a_lhs,
                        const LevelData<FArrayBox>& a_phi,
                        const LevelData<FArrayBox>& a_rhs,
                        bool a_homogeneous)
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
FieldSolver::create( LevelData<FArrayBox>&       a_lhs,
                     const LevelData<FArrayBox>& a_rhs )
{
   a_lhs.define(a_rhs.disjointBoxLayout(), a_rhs.nComp(), a_rhs.ghostVect());
}



void
FieldSolver::assign( LevelData<FArrayBox>&       a_lhs,
                     const LevelData<FArrayBox>& a_rhs )
{
   for (DataIterator dit(a_lhs.dataIterator()); dit.ok(); ++dit) {
      a_lhs[dit].copy(a_rhs[dit]);
   }
}



Real
FieldSolver::dotProduct( const LevelData<FArrayBox>& a_1,
                         const LevelData<FArrayBox>& a_2 )
{
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
FieldSolver::incr( LevelData<FArrayBox>&       a_lhs,
                   const LevelData<FArrayBox>& a_x,
                   Real                        a_scale )
{
   for (DataIterator dit(a_lhs.dataIterator()); dit.ok(); ++dit) {
      a_lhs[dit].plus(a_x[dit], a_scale);
   }
}



void
FieldSolver::axby( LevelData<FArrayBox>&       a_lhs,
                   const LevelData<FArrayBox>& a_x,
                   const LevelData<FArrayBox>& a_y,
                   Real                        a_a,
                   Real                        a_b )
{
   for (DataIterator dit(a_lhs.dataIterator()); dit.ok(); ++dit) {
      FArrayBox& this_fab = a_lhs[dit];
      this_fab.copy(a_x[dit]);
      this_fab *= a_a;
      this_fab.plus(a_y[dit], a_b);
   }
}



void
FieldSolver::scale( LevelData<FArrayBox>& a_lhs,
                    const Real&           a_scale )
{
   for (DataIterator dit(a_lhs.dataIterator()); dit.ok(); ++dit) {
      a_lhs[dit] *= a_scale;
   }
}



Real
FieldSolver::norm( const LevelData<FArrayBox>& a_rhs,
                   int                         a_ord )
{
   return CH_XD::norm(a_rhs, a_rhs.interval(), a_ord);
}



void
FieldSolver::setToZero( LevelData<FArrayBox>& a_lhs )
{
   for (DataIterator dit(a_lhs.dataIterator()); dit.ok(); ++dit) {
      a_lhs[dit].setVal(0.);
   }
}



// <---- End LinOp virtuals



void
FieldSolver::computeMappedPoloidalFieldWithGhosts( const LevelData<FArrayBox>& a_phi,
                                                   LevelData<FArrayBox>&       a_field,
                                                   const int                   a_order ) const
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

      for (int dir=0; dir<2; dir++) {
         FORT_CELL_CENTERED_FIELD_COMPONENT(CHF_BOX(box),
                                            CHF_CONST_INT(dir),
                                            CHF_CONST_FRA1(a_phi[dit],0),
                                            CHF_CONST_REALVECT(dx),
                                            CHF_CONST_INT(tmp_order),
                                            CHF_FRA1(a_field[dit],dir));
      }
   }

   a_field.exchange();

   // If fourth-order, recompute the field at valid cell centers to fourth-order

   if ( a_order == 4 ) {

      for (DataIterator dit(grids); dit.ok(); ++dit) {
         const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
         RealVect dx = block_coord_sys.dx();

         for (int dir=0; dir<2; ++dir) {
            FORT_CELL_CENTERED_FIELD_COMPONENT(CHF_BOX(grids[dit]),
                                               CHF_CONST_INT(dir),
                                               CHF_CONST_FRA1(a_phi[dit],0),
                                               CHF_CONST_REALVECT(dx),
                                               CHF_CONST_INT(a_order),
                                               CHF_FRA1(a_field[dit],dir));
         }
      }
   }
}



void
FieldSolver::computeMappedPoloidalFieldWithGhosts( const LevelData<FArrayBox>& a_phi,
                                                   LevelData<FluxBox>&         a_field,
                                                   const int                   a_order ) const
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
      for (int dir=0; dir<2; dir++) {
         IntVect grow_vect = 2*IntVect::Unit;
         grow_vect[dir] = 0;
         Box box = grow(grids[dit],grow_vect);

         FORT_FACE_INTERPOLATE(CHF_CONST_INT(dir),
                               CHF_BOX(surroundingNodes(box,dir)),
                               CHF_CONST_INT(a_order),
                               CHF_CONST_FRA1(a_phi[dit],0),
                               CHF_FRA1(phi_face[dit][dir],0));
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

      for (int dir=0; dir<2; dir++) {
         Box box_dir = surroundingNodes(box,dir);
         FORT_FACE_CENTERED_FIELD_COMPONENT(CHF_BOX(box_dir),
                                            CHF_CONST_INT(dir),
                                            CHF_CONST_FRA1(a_phi[dit],0),
                                            CHF_CONST_REALVECT(dx),
                                            CHF_CONST_INT(tmp_order),
                                            CHF_FRA1(a_field[dit][dir],dir));

         for (int tdir=0; tdir<2; ++tdir) {
            if (tdir != dir) {
               FORT_CELL_CENTERED_FIELD_COMPONENT(CHF_BOX(box_dir),
                                                  CHF_CONST_INT(tdir),
                                                  CHF_CONST_FRA1(phi_face[dit][dir],0),
                                                  CHF_CONST_REALVECT(dx),
                                                  CHF_CONST_INT(tmp_order),
                                                  CHF_FRA1(a_field[dit][dir],tdir));
            }
         }
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

         for (int dir=0; dir<2; dir++) {
            Box box_dir = surroundingNodes(grids[dit],dir);
            FORT_FACE_CENTERED_FIELD_COMPONENT(CHF_BOX(box_dir),
                                               CHF_CONST_INT(dir),
                                               CHF_CONST_FRA1(a_phi[dit],0),
                                               CHF_CONST_REALVECT(dx),
                                               CHF_CONST_INT(a_order),
                                               CHF_FRA1(a_field[dit][dir],dir));

            for (int tdir=0; tdir<2; ++tdir) {
               if (tdir != dir) {
                  FORT_CELL_CENTERED_FIELD_COMPONENT(CHF_BOX(box_dir),
                                                     CHF_CONST_INT(tdir),
                                                     CHF_CONST_FRA1(phi_face[dit][dir],0),
                                                     CHF_CONST_REALVECT(dx),
                                                     CHF_CONST_INT(a_order),
                                                     CHF_FRA1(a_field[dit][dir],tdir));
               }
            }
         }
      }
   }

   a_field.exchange();
}



void
FieldSolver::constructBoundaryStencils( const bool                        a_fourth_order,
                                        const PotentialBC&                a_bc,
                                        Vector< Vector<CoDim1Stencil> >&  a_codim1_stencils,
                                        Vector< Vector<CoDim2Stencil> >&  a_codim2_stencils ) const
{
   CH_assert(a_codim1_stencils.size() == 0);
   CH_assert(a_codim2_stencils.size() == 0);

   int order = a_fourth_order? 4: 2;

   const MultiBlockCoordSys* coord_sys_ptr = m_geometry.coordSysPtr();

   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = coord_sys_ptr->boundaries();

   int num_blocks = coord_sys_ptr->numBlocks();
   a_codim1_stencils.resize(num_blocks);
   a_codim2_stencils.resize(num_blocks);

   const Vector<Box>& mapping_blocks = coord_sys_ptr->mappingBlocks();

   for (int block_number=0; block_number<num_blocks; ++block_number) {
      const NewCoordSys* block_coord_sys = coord_sys_ptr->getCoordSys(block_number);
      Box domain_box = mapping_blocks[block_number];
      RealVect dx = block_coord_sys->dx();

      const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];

      Vector<CoDim1Stencil>& codim1_stencil = a_codim1_stencils[block_number];      
      Vector<CoDim2Stencil>& codim2_stencil = a_codim2_stencils[block_number];      

      int num_codim1_neighbors = 0;
      int num_codim2_neighbors = 0;

      for (int dir=0; dir<SpaceDim; ++dir) {
         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();

            if (this_block_boundaries[dir + side*SpaceDim].isDomainBoundary()) {

               codim1_stencil.resize(++num_codim1_neighbors);

               Box box = bdryBox(domain_box, dir, side, 1);
               FArrayBox bv(box,1);

               RefCountedPtr<GridFunction> bc_func = a_bc.getBCFunction(block_number, dir, side );
               if (bc_func) {
                  bc_func->assign(bv, m_geometry, box, 0., false);
               }
               else {
                  for (BoxIterator bit(box); bit.ok(); ++bit) {
                     bv(bit(),0) = a_bc.getBCValue(block_number, dir, side);
                  }
               }

               int bc_type = a_bc.getBCType(block_number, dir, side);

               codim1_stencil[num_codim1_neighbors-1].define(bv, bc_type, dx, dir, side, order);

               Box codim1_box = codim1_stencil[num_codim1_neighbors-1].box();

               for (int tdir=0; tdir < SpaceDim; ++tdir) {
                  if ( tdir != dir ) {
                     for (SideIterator sit2; sit2.ok(); ++sit2) {
                        Side::LoHiSide side2 = sit2();

                        bool transverse_boundary = this_block_boundaries[tdir + side2*SpaceDim].isDomainBoundary();

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



void
FieldSolver::accumPhysicalGhosts( const Vector< Vector<CoDim1Stencil> >&  a_codim1_stencils,
                                  const Vector< Vector<CoDim2Stencil> >&  a_codim2_stencils,
                                  const bool                              a_extrapolate_from_interior,
                                  const bool                              a_include_bvs,
                                  LevelData<FArrayBox>&                   a_data ) const
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
FieldSolver::extrapBoundaryGhostsForCC(FArrayBox&                              a_data,
                                       const Box&                              a_interiorbox,
                                       const Box&                              a_domain_box,
                                       const int                               a_dir,
                                       const Tuple<BlockBoundary, 2*SpaceDim>& a_block_boundaries) const
{
   // This function fourth-order extrapolates to fill two layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));

   for (int side=-1; side<2; side+=2) {

      if ( a_block_boundaries[a_dir + (side+1)*SpaceDim/2].isDomainBoundary() ) {

         Box dstbox;
         bool isbc = false;
      
         switch(side) 
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(a_dir) == a_domain_box.smallEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellLo(a_interiorbox, a_dir, 2);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(a_dir) == a_domain_box.bigEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellHi(a_interiorbox, a_dir, 2);
               }
               break;
            }

         if (isbc) {
            CH_assert(a_interiorbox.size(a_dir)>=5);
            CH_assert(bx.contains(dstbox));
            FORT_EXTRAP_FOR_CC_OPS(CHF_CONST_INT(a_dir),
                                   CHF_CONST_INT(side),
                                   CHF_BOX(dstbox),
                                   CHF_BOX(a_interiorbox),
                                   CHF_FRA(a_data));
         }
      }
   }
} 



void
FieldSolver::extrapBoundaryGhostsForFC(FArrayBox&                              a_data,
                                       const Box&                              a_interiorbox,
                                       const Box&                              a_domain_box,
                                       const int                               a_dir,
                                       const int                               a_order,
                                       const Tuple<BlockBoundary, 2*SpaceDim>& a_block_boundaries) const
{
   // This function fourth-order extrapolates to fill two layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));

   for (int side=-1; side<2; side+=2) {

      if ( a_block_boundaries[a_dir + (side+1)*SpaceDim/2].isDomainBoundary() ) {

         Box dstbox;
         bool isbc = false;
      
         switch(side) 
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(a_dir) == a_domain_box.smallEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellLo(a_interiorbox, a_dir, 2);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(a_dir) == a_domain_box.bigEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellHi(a_interiorbox, a_dir, 2);
               }
               break;
            }

         if (isbc) {
            CH_assert( (a_order==4 && a_interiorbox.size(a_dir)>=5)
                    || (a_order==2 && a_interiorbox.size(a_dir)>=3));
            CH_assert(bx.contains(dstbox));
            FORT_EXTRAP_FOR_FC_OPS(CHF_CONST_INT(a_dir),
                                   CHF_CONST_INT(side),
                                   CHF_CONST_INT(a_order),
                                   CHF_BOX(dstbox),
                                   CHF_BOX(a_interiorbox),
                                   CHF_FRA(a_data));
         }
      }
   }
} 



void
FieldSolver::injectPoloidalVector( const LevelData<FArrayBox>& a_poloidal_vector,
                                   LevelData<FArrayBox>&       a_vector ) const
{
   CH_assert(a_poloidal_vector.nComp() == 2);
   CH_assert(a_vector.nComp() == 3);

   const DisjointBoxLayout& grids = m_geometry.grids();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FArrayBox& this_poloidal_vector = a_poloidal_vector[dit];
      FArrayBox& this_vector = a_vector[dit];

      this_vector.copy(this_poloidal_vector, 0, 0, 1);
      this_vector.setVal(0., 1);
      this_vector.copy(this_poloidal_vector, 1, 2, 1);
   }   
}



void
FieldSolver::injectPoloidalVector( const LevelData<FluxBox>& a_poloidal_vector,
                                   LevelData<FluxBox>&       a_vector ) const
{
   CH_assert(a_poloidal_vector.nComp() == 2);
   CH_assert(a_vector.nComp() == 3);

   const DisjointBoxLayout& grids = m_geometry.grids();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FluxBox& this_poloidal_vector = a_poloidal_vector[dit];
      FluxBox& this_vector = a_vector[dit];
      for (int dir=0; dir<2; ++dir) {
         const FArrayBox& this_poloidal_vector_dir = this_poloidal_vector[dir];
         FArrayBox& this_vector_dir = this_vector[dir];

         this_vector_dir.copy(this_poloidal_vector_dir, 0, 0, 1);
         this_vector_dir.setVal(0., 1);
         this_vector_dir.copy(this_poloidal_vector_dir, 1, 2, 1);
      }
   }   
}



void
FieldSolver::projectPoloidalVector( const LevelData<FArrayBox>& a_vector,
                                    LevelData<FArrayBox>&       a_poloidal_vector ) const
{
   CH_assert(a_vector.nComp() == 3);
   CH_assert(a_poloidal_vector.nComp() == 2);

   const DisjointBoxLayout& grids = m_geometry.grids();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FArrayBox& this_vector = a_vector[dit];
      FArrayBox& this_poloidal_vector = a_poloidal_vector[dit];

      this_poloidal_vector.copy(this_vector, 0, 0, 1);
      this_poloidal_vector.copy(this_vector, 2, 1, 1);
   }   
}



void
FieldSolver::projectPoloidalVector( const LevelData<FluxBox>& a_vector,
                                    LevelData<FluxBox>&       a_poloidal_vector ) const
{
   CH_assert(a_vector.nComp() == 3);
   CH_assert(a_poloidal_vector.nComp() == 2);

   const DisjointBoxLayout& grids = m_geometry.grids();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const FluxBox& this_vector = a_vector[dit];
      FluxBox& this_poloidal_vector = a_poloidal_vector[dit];
      for (int dir=0; dir<2; ++dir) {
         const FArrayBox& this_vector_dir = this_vector[dir];
         FArrayBox& this_poloidal_vector_dir = this_poloidal_vector[dir];

         this_poloidal_vector_dir.copy(this_vector_dir, 0, 0, 1);
         this_poloidal_vector_dir.copy(this_vector_dir, 2, 1, 1);
      }
   }   
}



double
FieldSolver::L2Norm( const LevelData<FArrayBox>& a_data )
{
   return sqrt(dotProduct(a_data, a_data));
}



double
FieldSolver::globalMax(const double a_data) const
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
