#include "RosenbluthPotentials.H"
#include "Directions.H"

#include "KineticSpecies.H"
#include "RosenbluthPotentialsF_F.H"

#include "NamespaceHeader.H" 

RosenbluthPotentials::RosenbluthPotentials(LevelData<FArrayBox>& a_phi_one,
                      LevelData<FArrayBox>& a_phi_two,
                      const LevelData<FArrayBox>& a_rho,
                      const PhaseGeom& a_phase_geom,
                      const Real a_mass,
                      const Real a_pcg_tol,
                      const Real a_pcg_maxiter,
                      const int a_mult_num,
                      const int a_verbocity) 


   : m_verbosity(a_verbocity),
     m_phase_geom(a_phase_geom),
     m_mass(a_mass),
     m_pcg_tol(a_pcg_tol),
     m_pcg_maxiter(a_pcg_maxiter),
     m_mult_num(a_mult_num)
{

   const DisjointBoxLayout& grids( a_rho.getBoxes() );
   const int n_comp( a_rho.nComp() );
   LevelData<FArrayBox> rhs_withBC( grids, n_comp, IntVect::Zero );
   
   imposeBC1(rhs_withBC, a_rho);
   solve(a_phi_one,rhs_withBC);

   imposeBC2(rhs_withBC, a_phi_one, a_rho);
   solve(a_phi_two,rhs_withBC);

}

RosenbluthPotentials::~RosenbluthPotentials()
{
}


void RosenbluthPotentials::imposeBC1(LevelData<FArrayBox>& a_rhsBC,
                                 const LevelData<FArrayBox>& a_rho) const

{

   //Compute multipole coefficients 
   const CFG::MagGeom& mag_geom( m_phase_geom.magGeom() );
   CFG::LevelData<CFG::FArrayBox> multipole_coeff( mag_geom.grids(), m_mult_num, CFG::IntVect::Zero );
   computeMultipoleCoeff(multipole_coeff, a_rho); 
   LevelData<FArrayBox> injected_mult_coeff;
   m_phase_geom.injectConfigurationToPhase( multipole_coeff, injected_mult_coeff);

   // Get velocity coordinate system parameters
   const VEL::VelCoordSys& vel_coords = m_phase_geom.velSpaceCoordSys();
   const VEL::ProblemDomain& vel_domain = vel_coords.domain();
   const VEL::Box& domain_box = vel_domain.domainBox();
   const VEL::RealVect& vel_dx = vel_coords.dx();
   const int num_vpar_cells = domain_box.size(0);
   const int num_mu_cells = domain_box.size(1);

   const LevelData<FArrayBox>& injected_B = m_phase_geom.getBFieldMagnitude();

   //Impose multipole boundary conditions to the rhs (i.e., b vector) of Ax=b 
   for (DataIterator dit( a_rhsBC.dataIterator() ); dit.ok(); ++dit) {
      a_rhsBC[dit].copy( a_rho[dit] );
      FArrayBox& this_rhsBC = a_rhsBC[dit];
      FArrayBox& this_mult_coeff = injected_mult_coeff[dit];
      const FArrayBox& this_b = injected_B[dit];

      FORT_IMPOSE_MULTIPOLE_BC1(CHF_FRA1(this_rhsBC,0),
                               CHF_CONST_FRA(this_mult_coeff),
                               CHF_CONST_FRA1(this_b,0),
                               CHF_BOX(this_rhsBC.box()),
                               CHF_CONST_REALVECT(vel_dx),
                               CHF_CONST_INT(num_vpar_cells),
                               CHF_CONST_INT(num_mu_cells),
                               CHF_CONST_REAL(m_mass));
   }
}

void RosenbluthPotentials::imposeBC2(LevelData<FArrayBox>& a_rhsBC,
                                 const LevelData<FArrayBox>& a_phi_one,
                                 const LevelData<FArrayBox>& a_rho) const

{
   // Get velocity coordinate system parameters and Bfield
   const VEL::VelCoordSys& vel_coords = m_phase_geom.velSpaceCoordSys();
   const VEL::ProblemDomain& vel_domain = vel_coords.domain();
   const VEL::Box& domain_box = vel_domain.domainBox();
   const VEL::RealVect& vel_dx = vel_coords.dx();
   const int num_vpar_cells = domain_box.size(0);
   const int num_mu_cells = domain_box.size(1);
   const LevelData<FArrayBox>& injected_B = m_phase_geom.getBFieldMagnitude();

   //Separate phi_one into inner (depends on a_rho profile, 
   //and is confined inside Vr_domain_min) and outer parts (analytic)
   LevelData<FArrayBox> phi_inner( a_phi_one.getBoxes(), 1, IntVect::Zero );

   for (DataIterator dit( a_rhsBC.dataIterator() ); dit.ok(); ++dit) {
      FArrayBox& this_phi_inner = phi_inner[dit];
      const FArrayBox& this_phi = a_phi_one[dit];
      const FArrayBox& this_b = injected_B[dit];

      FORT_SEPARATE_PHI1_INNER(CHF_FRA1(this_phi_inner,0),
                              CHF_CONST_FRA1(this_phi,0),
                              CHF_CONST_FRA1(this_b,0),
                              CHF_BOX(this_phi.box()),
                              CHF_CONST_REALVECT(vel_dx),
                              CHF_CONST_INT(num_vpar_cells),
                              CHF_CONST_INT(num_mu_cells),
                              CHF_CONST_REAL(m_mass));
   }

   //Calculate multipole coefficients for the inner part of phi_one
   const CFG::MagGeom& mag_geom( m_phase_geom.magGeom() );
   CFG::LevelData<CFG::FArrayBox> mult_coeff_phi1( mag_geom.grids(), m_mult_num, CFG::IntVect::Zero );
   computeMultipoleCoeff(mult_coeff_phi1, phi_inner); 
   LevelData<FArrayBox> inj_mcoeff_phi1;
   m_phase_geom.injectConfigurationToPhase( mult_coeff_phi1, inj_mcoeff_phi1);

   //Compute multipole coefficients for rho (used for analytic expression for BC driven by phi_one_outer)
   CFG::LevelData<CFG::FArrayBox> mult_coeff_rho( mag_geom.grids(), m_mult_num, CFG::IntVect::Zero );
   computeMultipoleCoeff(mult_coeff_rho, a_rho); 
   LevelData<FArrayBox> inj_mcoeff_rho; 
   m_phase_geom.injectConfigurationToPhase( mult_coeff_rho, inj_mcoeff_rho);

   //Impose multipole boundary conditions to the rhs (i.e., b vector) of Ax=b 
   for (DataIterator dit( a_rhsBC.dataIterator() ); dit.ok(); ++dit) {
      a_rhsBC[dit].copy( a_phi_one[dit] );
      FArrayBox& this_rhsBC = a_rhsBC[dit];
      FArrayBox& this_mcoeff_rho = inj_mcoeff_rho[dit];
      FArrayBox& this_mcoeff_phi1 = inj_mcoeff_phi1[dit];
      const FArrayBox& this_b = injected_B[dit];

      FORT_IMPOSE_MULTIPOLE_BC2(CHF_FRA1(this_rhsBC,0),
                               CHF_CONST_FRA(this_mcoeff_rho),
                               CHF_CONST_FRA(this_mcoeff_phi1),
                               CHF_CONST_FRA1(this_b,0),
                               CHF_BOX(this_rhsBC.box()),
                               CHF_CONST_REALVECT(vel_dx),
                               CHF_CONST_INT(num_vpar_cells),
                               CHF_CONST_INT(num_mu_cells),
                               CHF_CONST_REAL(m_mass));
   }
}

void RosenbluthPotentials::computeMultipoleCoeff(CFG::LevelData<CFG::FArrayBox>& a_mult_coeff,
                                                 const LevelData<FArrayBox>& a_rho) const

{

   //Create the untility species object for calculating moments
   string name;
   name.assign("utility_species");

   RefCountedPtr<KineticSpecies> utility_spec 
      = RefCountedPtr<KineticSpecies>(
              new KineticSpecies( name, m_mass, 1, m_phase_geom ) );

   utility_spec->distributionFunction().define(a_rho.getBoxes(), m_mult_num, IntVect::Zero);
   LevelData<FArrayBox>& mult_kernels( utility_spec->distributionFunction() );


   // Get velocity coordinate system parameters and magnetic field
   const VEL::VelCoordSys& vel_coords = m_phase_geom.velSpaceCoordSys();
   const VEL::RealVect& vel_dx = vel_coords.dx();
   const LevelData<FArrayBox>& injected_B = m_phase_geom.getBFieldMagnitude();

   //Compute kernels for multipole calculations
   for (DataIterator dit( mult_kernels.dataIterator() ); dit.ok(); ++dit) {
      FArrayBox& this_kernel = mult_kernels[dit];
      const FArrayBox& this_rho = a_rho[dit];
      const FArrayBox& this_b = injected_B[dit];

      FORT_COMPUTE_MULT_KERNELS(CHF_FRA(this_kernel),
                               CHF_CONST_FRA1(this_rho,0),
                               CHF_CONST_FRA1(this_b,0),
                               CHF_BOX(this_kernel.box()),
                               CHF_CONST_REALVECT(vel_dx),
                               CHF_CONST_REAL(m_mass));
   }

   utility_spec->chargeDensity( a_mult_coeff );   
}

void RosenbluthPotentials::solve(LevelData<FArrayBox>& a_solution,
                                 const LevelData<FArrayBox>& a_rhs) const
{
   // Get coordinate system parameters 
   const PhaseGrid& phase_grid = m_phase_geom.phaseGrid();
   const DisjointBoxLayout& dbl = phase_grid.disjointBoxLayout();

   const VEL::VelCoordSys& vel_geom = m_phase_geom.velSpaceCoordSys();
   
   const LevelData<FArrayBox>& injected_B = m_phase_geom.getBFieldMagnitude();

   int num_config_boxes = phase_grid.numConfigBoxes();
   
   int max_iterations = 0;
   double max_final_norm = 0.;

   for (int k=0; k<num_config_boxes; ++k) {

      const MPI_Comm& config_box_comm = phase_grid.configBoxComm(k);
      const List<VEL::Box>& velocity_slice = phase_grid.velocitySlice(k);

      // Make the Hypre solver and vectors

      HYPRE_StructMatrix A;
      HYPRE_StructVector x;
      HYPRE_StructVector b;
      HYPRE_StructSolver solver;
      createHypreData(config_box_comm, velocity_slice, A, x, b, solver);

      // Loop over configuration space

      const CFG::Box& config_box = phase_grid.configBox(k);
      CFG::BoxIterator bit(config_box);
      for (bit.begin(); bit.ok(); ++bit) {
         CFG::IntVect cfg_iv = bit();

         List<Box> phase_slices;
         getSliceBoxes(cfg_iv, velocity_slice, phase_slices);

         //Get local B 
         double B_local;

         bool box_found = false;
         for (ListIterator<Box> it(phase_slices); it.ok(); ++it) {
            for (DataIterator dit(dbl); dit.ok(); ++dit) {
               Box overlap = dbl[dit] & it();
               if (overlap.ok()) {

                 const FArrayBox& this_B = injected_B[dit];   
                 const Box& Bbox = this_B.box();
                 int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
                 int mu_index = Bbox.smallEnd(MU_DIR);
                 IntVect ivB;
                 ivB[VPARALLEL_DIR] = vp_index;
                 ivB[MU_DIR] = mu_index;
                 ivB[RADIAL_DIR] = cfg_iv[0];
                 ivB[POLOIDAL_DIR] = cfg_iv[1];

                 B_local = this_B(ivB,0);

                 box_found = true;
                 break;
               }

              if (box_found) break;
            }   
         }

         // Fill the matrix
         constructMatrix(vel_geom.domain(), vel_geom.dx(), velocity_slice, B_local, A);

         // Load the Hypre right-hand side and set the initial guess to zero

         for (ListIterator<Box> it(phase_slices); it.ok(); ++it) {
            for (DataIterator dit(dbl); dit.ok(); ++dit) {
               Box overlap = dbl[dit] & it();

               if (overlap.ok()) {
                  FArrayBox this_data(overlap,1);
                  this_data.copy(a_rhs[dit], overlap);

                  VEL::IntVect vlo = m_phase_geom.vel_restrict(overlap.smallEnd());
                  VEL::IntVect vhi = m_phase_geom.vel_restrict(overlap.bigEnd());
                  HYPRE_StructVectorSetBoxValues(b, vlo.dataPtr(), vhi.dataPtr(), this_data.dataPtr());

                  this_data.copy(a_solution[dit], overlap);
                  HYPRE_StructVectorSetBoxValues(x, vlo.dataPtr(), vhi.dataPtr(), this_data.dataPtr());
               }

            }
         }
         HYPRE_StructVectorAssemble(b);
         HYPRE_StructVectorAssemble(x);

         // Solve the linear system         

         HYPRE_StructPCGSetup(solver, A, b, x);
         HYPRE_StructPCGSolve(solver, A, b, x);

         // Get the number of iterations used and the final residual

         int num_iterations;
         HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
         
         if (num_iterations > max_iterations) max_iterations = num_iterations;

         double final_norm;
         HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_norm);

         if (final_norm > max_final_norm) max_final_norm = final_norm;

         //         cout << "      RosenbluthPotentials PCG residual at " << cfg_iv << " = " << final_norm << " after " << num_iterations << " iterations" << endl;

         // Get the data from Hypre solution vector

         for (ListIterator<Box> it(phase_slices); it.ok(); ++it) {
            for (DataIterator dit(dbl); dit.ok(); ++dit) {
               Box overlap = dbl[dit] & it();

               if (overlap.ok()) {
                  FArrayBox this_data(overlap,1);
                  VEL::IntVect vlo = m_phase_geom.vel_restrict(overlap.smallEnd());
                  VEL::IntVect vhi = m_phase_geom.vel_restrict(overlap.bigEnd());
                  HYPRE_StructVectorGetBoxValues(x, vlo.dataPtr(), vhi.dataPtr(), this_data.dataPtr());

                  a_solution[dit].copy(this_data, overlap);
               }
            }
         }
      }

      destroyHypreData(A, x, b, solver);
   }

   if (m_verbosity) {
#ifdef CH_MPI
      double local_max_final_norm = max_final_norm;
      int local_max_iterations = max_iterations;
      
      MPI_Allreduce(&local_max_final_norm, &max_final_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&local_max_iterations, &max_iterations, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

      if (procID()==0) {
         cout << "      RosenbluthPotentials PCG max residual = " << max_final_norm << " after " << max_iterations << " iterations (max)" << endl;
      }
   }

}

void RosenbluthPotentials::createHypreData(const MPI_Comm&       a_comm,
                                   const List<VEL::Box>& a_boxes,
                                   HYPRE_StructMatrix&   a_matrix,
                                   HYPRE_StructVector&   a_x,
                                   HYPRE_StructVector&   a_b,
                                   HYPRE_StructSolver&   a_solver) const
{
   /* Create the grid */

   HYPRE_StructGrid grid;
   HYPRE_StructGridCreate(a_comm, VEL_DIM, &grid);

   for (ListIterator<VEL::Box> it(a_boxes); it.ok(); ++it) {
      const VEL::Box& box = it();

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      HYPRE_StructGridSetExtents(grid, lower.dataPtr(), upper.dataPtr());
   }

   HYPRE_StructGridAssemble(grid);

   /* Define the discretization stencil for the system matrix A */

   /* Define the geometry of the stencil. Each represents a
      relative offset (in the index space). */

   int offsets[][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};
   int stencil_size = 5;

   /* Create an empty 2D stencil object */
   HYPRE_StructStencil stencil;
   HYPRE_StructStencilCreate(VEL_DIM, stencil_size, &stencil);

   /* Assign each of the stencil entries */
   int entry;
   for (entry = 0; entry < stencil_size; entry++)
      HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
   
   /* Create an empty matrix object */
   HYPRE_StructMatrixCreate(a_comm, grid, stencil, &a_matrix);
   HYPRE_StructMatrixInitialize(a_matrix);


   /* Set up Struct Vectors for b and x. */

   /* Create an empty vector object */
   HYPRE_StructVectorCreate(a_comm, grid, &a_b);
   HYPRE_StructVectorCreate(a_comm, grid, &a_x);

   HYPRE_StructVectorInitialize(a_b);
   HYPRE_StructVectorInitialize(a_x);

   /* Set up the solver (See the Hypre Reference Manual for
      descriptions of all of the options.) */

   // Make the PCG solver

   HYPRE_StructPCGCreate(a_comm, &a_solver);

   // Set the PCG solver parameters
   HYPRE_StructPCGSetTol(a_solver, m_pcg_tol); /* convergence tolerance */
   HYPRE_StructPCGSetMaxIter(a_solver, m_pcg_maxiter);
   HYPRE_StructPCGSetLogging(a_solver, 1);
   HYPRE_StructPCGSetPrintLevel(a_solver, 0); /* amount of info. printed */

/*
   // Make the preconditioner

   HYPRE_StructSolver precond;
   HYPRE_StructPFMGCreate(a_comm, &precond);
   HYPRE_StructPFMGSetMaxIter(precond, 1);
   HYPRE_StructPFMGSetTol(precond, 0.0);
   HYPRE_StructPFMGSetZeroGuess(precond);
   HYPRE_StructPFMGSetRAPType(precond, 0);
   HYPRE_StructPFMGSetRelaxType(precond, 1);
   HYPRE_StructPFMGSetNumPreRelax(precond, 1);
   HYPRE_StructPFMGSetNumPostRelax(precond, 1);
   HYPRE_StructPFMGSetPrintLevel(precond, 0);
   HYPRE_StructPFMGSetLogging(precond, 0);
   HYPRE_StructPCGSetPrecond(a_solver,
                             HYPRE_StructPFMGSolve,
                             HYPRE_StructPFMGSetup,
                             precond);

*/
   HYPRE_StructStencilDestroy(stencil);
   HYPRE_StructGridDestroy(grid);
}


void RosenbluthPotentials::constructMatrix( const VEL::ProblemDomain&  a_domain,
                                    const VEL::RealVect&       a_dx,
                                    const List<VEL::Box>&      a_boxes,
                                    const double&              a_Bfield,
                                    HYPRE_StructMatrix&        a_matrix) const
{
   /*
     Fills a_matrix with the coeffients corresponding to the standard
     second-order discretization of the negative Laplacian with homogeneous
     Dirichlet boundary conditions.
   */

   VEL::RealVect dx_inv2;
   for (int n=0; n<VEL_DIM; ++n) {
      dx_inv2[n] = 1. / (a_dx[n] * a_dx[n]);
   } 

   int nentries = 5;
   int* stencil_indices = new int[nentries];
   for (int j = 0; j < nentries; j++) {
      stencil_indices[j] = j;
   }

   const VEL::Box& domain_box = a_domain.domainBox();
   VEL::IntVect domain_lo = domain_box.smallEnd();
   VEL::IntVect domain_hi = domain_box.bigEnd();

   for (ListIterator<VEL::Box> it(a_boxes); it.ok(); ++it) {
      const VEL::Box& box = it();

      IntVect lower(box.loVect());
      IntVect upper(box.hiVect());

      int n = box.numPts();
      int nvalues = nentries*n;

      double* values = new double[nvalues];

      int i = 0;
      VEL::BoxIterator bit(box);
      for (bit.begin(); bit.ok(); ++bit) {
         VEL::IntVect iv = bit();

         double mu_coeff = 4.0 * (m_mass/a_Bfield) * (iv[1]+0.5) * a_dx[1];
         double mu_coeff_HiFace = 4.0 * (m_mass/a_Bfield) * (iv[1]+1) * a_dx[1];
         double mu_coeff_LoFace = 4.0 * (m_mass/a_Bfield) * (iv[1]) * a_dx[1];

         // Second-order centered difference Laplacian 
         values[i] = -2.0 * (dx_inv2[0] + mu_coeff * dx_inv2[1]);
         values[i+1] = values[i+2] = dx_inv2[0];
         values[i+3] = mu_coeff_LoFace * dx_inv2[1];
         values[i+4] = mu_coeff_HiFace * dx_inv2[1];

         // Homogeneous Dirichlet boundary condition modifications
         if (iv[0] == domain_lo[0]) {
             values[i+1] = 0.;         
         }
         if (iv[0] == domain_hi[0]) {
            values[i+2] = 0.;
         }
         if (iv[1] == domain_lo[1]) {
            values[i+3] = 0.;
         }
         if (iv[1] == domain_hi[1]) {
            values[i+4] = 0.;
         }

         i += nentries;
      }

      HYPRE_StructMatrixSetBoxValues(a_matrix, lower.dataPtr(), upper.dataPtr(),
                                     nentries, stencil_indices, values);

      delete [] values;
   }

   HYPRE_StructMatrixAssemble(a_matrix);

   delete [] stencil_indices;
}

void RosenbluthPotentials::destroyHypreData(HYPRE_StructMatrix& a_matrix,
                                    HYPRE_StructVector& a_x,
                                    HYPRE_StructVector& a_b,
                                    HYPRE_StructSolver& a_solver) const
{
   HYPRE_StructPCGDestroy(a_solver);
   HYPRE_StructVectorDestroy(a_b);
   HYPRE_StructVectorDestroy(a_x);
   HYPRE_StructMatrixDestroy(a_matrix);
}



void RosenbluthPotentials::getSliceBoxes(const CFG::IntVect     a_cfg_iv,
                                 const List<VEL::Box>&  a_vel_boxes,
                                 List<Box>&             a_phase_boxes) const
{
   IntVect lo, hi;
   for (int n=0; n<CFG_DIM; ++n) {
      lo[n] = hi[n] = a_cfg_iv[n];
   }

   for (ListIterator<VEL::Box> it(a_vel_boxes); it.ok(); ++it) {
      const VEL::Box& vel_box = it();

      VEL::IntVect vel_box_lo = vel_box.smallEnd();
      VEL::IntVect vel_box_hi = vel_box.bigEnd();
      for (int n=CFG_DIM; n<SpaceDim; ++n) {
         lo[n] = vel_box_lo[n-CFG_DIM];
         hi[n] = vel_box_hi[n-CFG_DIM];
      }

      a_phase_boxes.add(Box(lo,hi));
   }
}
                                 
#include "NamespaceFooter.H"
