#include "MBTridiagonalSolver.H"
#include "Directions.H"
#include "LogRectEllipticOpBC.H"
#include "SNCoreEllipticOpBC.H"

#include "NamespaceHeader.H"

extern "C" {
   void dgtsv_(const int&, const int &, double *, double *, double *, double *, const int &, int &);
   void dgesv_(const int&, const int &, double *, const int &, int *, double *, const int &, int &);
}



MBTridiagonalSolver::MBTridiagonalSolver( const MultiBlockLevelGeom&  a_geom,
                                          const int                   a_discretization_order )
   : MBSolver(a_geom, a_discretization_order),
     m_flux_surface((MagGeom&)a_geom)

{
   IntVect stencil_box_lo(IntVect::Zero);
   IntVect stencil_box_hi;
   if (m_discretization_order == 4) {
      stencil_box_hi = 4*IntVect::Unit;
   }
   else {
      stencil_box_hi = 2*IntVect::Unit;
   }
   Box stencil_box(stencil_box_lo, stencil_box_hi);
          
   m_A_stencil_values.define(stencil_box, 1);

   m_A_diagonal_offset = (stencil_box.numPts() - 1) / 2;

   m_A_radial.define(m_geometry.gridsFull(), 3, IntVect::Zero);
}
      


void
MBTridiagonalSolver::constructMatrixGeneral( LevelData<FArrayBox>&  a_alpha_coefficient,
                                             LevelData<FluxBox>&    a_tensor_coefficient,
                                             LevelData<FArrayBox>&  a_beta_coefficient,
                                             const EllipticOpBC&    a_bc )
{
   bool fourthOrder = (m_discretization_order == 4);

   m_periodic_or_neumann = isCoreRadialPeriodicOrNeumannBC(a_bc);

   constructTridiagonalMatrix(a_alpha_coefficient, a_tensor_coefficient, a_beta_coefficient, a_bc,
                              m_A_stencil_values, m_A_radial, m_A_diagonal_offset, fourthOrder,
                              m_rhs_from_bc);
}



void
MBTridiagonalSolver::constructMatrixBlockGeneral( const int              a_block_row,
                                                  const int              a_block_column,
                                                  LevelData<FArrayBox>&  a_alpha_coefficient,
                                                  LevelData<FluxBox>&    a_tensor_coefficient,
                                                  LevelData<FArrayBox>&  a_beta_coefficient,
                                                  const EllipticOpBC&    a_bc )
{
   MayDay::Error("MBTridiagonalSolver::constructMatrixBlockGeneral() is not implemented");
}



void
MBTridiagonalSolver::multiplyMatrix( const LevelData<FArrayBox>&  a_in,
                                     LevelData<FArrayBox>&        a_out ) const
{
   MayDay::Error("MBTridiagonalSolver::multiplyMatrix() not implemented" );
}



void
MBTridiagonalSolver::solve( const LevelData<FArrayBox>&  a_rhs,
                            LevelData<FArrayBox>&        a_solution,
                            bool                         a_homogeneous_bcs )
{
  CH_TIME("MBTridiagonalSolver::solve");

   // This function assumes that boundary values have already been 
   // included in the right-hand side, so the a_homogeneous_bcs argument
   // plays no role.  We therefore require it to be set to false and
   // perform an assertion to ensure that no attempt is made to try to
   // use it.
   CH_assert(!a_homogeneous_bcs);

   for (DataIterator dit(a_solution.dataIterator()); dit.ok(); ++dit) {
      a_solution[dit].copy(a_rhs[dit]);
   }

   solveFluxSurfaceAverage(a_solution);
}



void
MBTridiagonalSolver::dumpMatrix(const string&  a_file_name) const
{
   MayDay::Error("MBTridiagonalSolver::dumpMatrix() not implemented");
}



void
MBTridiagonalSolver::constructTridiagonalMatrix( LevelData<FArrayBox>&               a_alpha_coefficient, 
                                                 LevelData<FluxBox>&                 a_tensor_coefficient,
                                                 LevelData<FArrayBox>&               a_beta_coefficient,
                                                 const EllipticOpBC&                 a_bc,
                                                 FArrayBox&                          a_stencil_values,
                                                 LevelData<FArrayBox>&               a_radial,
                                                 const int                           a_diagonal_offset,
                                                 const bool                          a_fourth_order,
                                                 LevelData<FArrayBox>&               a_rhs_from_bc ) const
{
   // Make sure the coefficient ghost cell values are set.
   // The ghost cell values at physical boundaries aren't used.
   a_tensor_coefficient.exchange();

   int stencil_size = a_stencil_values.box().numPts();
    
   double * stencil_values = a_stencil_values.dataPtr(0);

   const DisjointBoxLayout & grids = a_tensor_coefficient.disjointBoxLayout();

   setZero(a_rhs_from_bc);

   LevelData<FArrayBox> structured_values(grids, stencil_size, IntVect::Zero);
   FArrayBox tmp_stencil_values(a_stencil_values.box(), a_stencil_values.nComp());

   Vector< Vector<CoDim1Stencil> > codim1_stencils;
   Vector< Vector<CoDim2Stencil> > codim2_stencils;

   constructBoundaryStencils(a_fourth_order, a_bc, codim1_stencils, codim2_stencils );

   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_coord_sys_ptr->boundaries();

   const Vector<Box>& mapping_blocks = m_coord_sys_ptr->mappingBlocks();

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     const Box & box = grids[dit];
     int block_number = m_coord_sys_ptr->whichBlock(box);
     const NewCoordSys& block_coord_sys = *(m_coord_sys_ptr->getCoordSys(block_number));
     RealVect dx = block_coord_sys.dx();

     const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];
     bool lo_block_interface[SpaceDim];
     bool hi_block_interface[SpaceDim];
     if (m_mblex_potential_Ptr) {
        for (int dir=0; dir<SpaceDim; ++dir) {
           lo_block_interface[dir] = this_block_boundaries[dir].isInterface();
           hi_block_interface[dir] = this_block_boundaries[dir + SpaceDim].isInterface();
        }
     }
     else {
        for (int dir=0; dir<SpaceDim; ++dir) {
           lo_block_interface[dir] = false;
           hi_block_interface[dir] = false;
        }
     }

     const Box& domainBox = mapping_blocks[block_number];

     FluxBox& this_coef = a_tensor_coefficient[dit];

     const FArrayBox* alpha = NULL;
     if ( a_alpha_coefficient.isDefined() ) {
        alpha = &a_alpha_coefficient[dit];
     }

     const FArrayBox* beta = NULL;
     if ( a_beta_coefficient.isDefined() ) {
        beta = &a_beta_coefficient[dit];
     }

     // Set the intra-block matrix entries

     BoxIterator bit(box);
     for (bit.begin();bit.ok();++bit) {
       IntVect iv = bit();

       a_stencil_values.setVal(0.);

       // Compute the stencil couplings to cell iv without regard to boundary conditions, skipping
       // the contributions due to cell faces on block interfaces, which are handled later by
       // addUnstructuredMatrixEntries().

       for (int dir=0; dir<SpaceDim; ++dir) {
          bool at_lo_interface = (lo_block_interface[dir] && iv[dir] == domainBox.smallEnd(dir));
          bool at_hi_interface = (hi_block_interface[dir] && iv[dir] == domainBox.bigEnd(dir));
          bool force_codim2_condense = lo_block_interface[dir] || hi_block_interface[dir];

          for (SideIterator sit; sit.ok(); ++sit) {
             Side::LoHiSide side = sit();

             bool update_rhs_from_bc_only = (at_lo_interface && side == Side::LoHiSide::Lo)
                || (at_hi_interface && side == Side::LoHiSide::Hi);

             for (int dir2=0; dir2<SpaceDim; ++dir2) {

                tmp_stencil_values.setVal(0.);

                accumStencilMatrixEntries(iv, dir, side, dir2, this_coef, dx,
                                          a_fourth_order, tmp_stencil_values);

                modifyStencilForBCs( codim1_stencils[block_number], codim2_stencils[block_number],
                                     iv, tmp_stencil_values, a_rhs_from_bc[dit],
                                     update_rhs_from_bc_only, force_codim2_condense );

                a_stencil_values += tmp_stencil_values;
             }
          }
       }

       if ( alpha ) {
          a_stencil_values *= alpha->operator()(iv);
          a_rhs_from_bc[dit](iv,0) *= alpha->operator()(iv);
       }

       if ( beta ) {
          stencil_values[a_diagonal_offset] += beta->operator()(iv);
       }

       for (int jj=0; jj<stencil_size; jj++) {
          structured_values[dit](iv,jj) = stencil_values[jj];
       }
     }
   }

   // Average poloidally
   FluxSurface flux_surface_op((MagGeom&)m_geometry);

   LevelData<FArrayBox> fs_average(flux_surface_op.grids(), stencil_size, IntVect::Zero);
   flux_surface_op.average(structured_values, fs_average);
   flux_surface_op.spread(fs_average, structured_values);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     const Box & box = grids[dit];

     a_radial[dit].setVal(0.);

     double * radial_ldiag = a_radial[dit].dataPtr(0);
     double * radial_diag  = a_radial[dit].dataPtr(1);
     double * radial_udiag = a_radial[dit].dataPtr(2);

     int diagonal_index = 0;

     BoxIterator bit(box);
     for (bit.begin();bit.ok();++bit) {
       IntVect iv = bit();

       for (int jj=0; jj<stencil_size; jj++) {
          stencil_values[jj] = structured_values[dit](iv,jj);
       }

       for (int i=0; i<=6*(3*SpaceDim-5); i+=3) {
          radial_ldiag[diagonal_index] += stencil_values[i];
          radial_diag[diagonal_index]  += stencil_values[i+1];
          radial_udiag[diagonal_index] += stencil_values[i+2];
       }

       diagonal_index++;

     }
   }

#if 0
   testMatrixConstruct( a_alpha_coefficient, 
                        a_tensor_coefficient,
                        a_beta_coefficient,
                        block_boundaries,
                        codim1_stencils,
                        codim2_stencils,
                        a_stencil_values,
                        a_fourth_order,
                        a_rhs_from_bc );
#endif
}



void
MBTridiagonalSolver::solveFluxSurfaceAverage( LevelData<FArrayBox>& a_data ) const
{
   const ProblemDomain& domain = m_geometry.gridsFull().physDomain();
   const Box& domain_box = domain.domainBox();
   const bool radially_periodic = domain.isPeriodic(RADIAL_DIR);
   const int n = domain_box.size(RADIAL_DIR);

   int radial_lo = domain_box.smallEnd(RADIAL_DIR);

   Box domain_box_fs_intersection = adjCellLo(domain_box, POLOIDAL_DIR, -1);
#if CFG_DIM==3
   domain_box_fs_intersection = adjCellLo(domain_box_fs_intersection, TOROIDAL_DIR, -1);   
#endif   
   
   const DisjointBoxLayout& fs_grids = m_flux_surface.grids();
   LevelData<FArrayBox> Pr(fs_grids, 3, IntVect::Zero);
   m_flux_surface.average(m_A_radial, Pr);

   double * temp  = new double[4*n];
   for (int i=0; i<4*n; ++i) temp[i] = 0.;

   double * ldiag = &(temp[0]);
   double * diag  = &(temp[n]);
   double * udiag = &(temp[2*n]);
   double * bx    = &(temp[3*n]);

   double *send_buffer = new double[n];

   DataIterator dit = fs_grids.dataIterator();
   for (int comp=0; comp<Pr.nComp(); ++comp) {

      for (int i=0; i<n; ++i) send_buffer[i] = 0.;

      for (dit.begin(); dit.ok(); ++dit) {
         Box overlap = fs_grids[dit] & domain_box_fs_intersection;
         if ( overlap.ok() ) {
            FArrayBox& this_fab = Pr[dit];
            for (BoxIterator bit(overlap); bit.ok(); ++bit) {
               IntVect iv = bit();
               send_buffer[iv[RADIAL_DIR]-radial_lo] = this_fab(iv,comp);
            }
         }
      }

#ifdef CH_MPI
      MPI_Allreduce(send_buffer, &temp[comp*n], n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      for (int i=0; i<n; ++i) {
        temp[comp*n+i] = send_buffer[i];
      }
#endif
   }

   for (int i=0; i<n; ++i) send_buffer[i] = 0.;

   for (dit.begin(); dit.ok(); ++dit) {
      Box overlap = fs_grids[dit] & domain_box_fs_intersection;
      if ( overlap.ok() ) {
         FArrayBox& this_fab = a_data[dit];
         for (BoxIterator bit(overlap); bit.ok(); ++bit) {
            IntVect iv = bit();
            send_buffer[iv[RADIAL_DIR]-radial_lo] = this_fab(iv);
         }
      }
   }

#ifdef CH_MPI
   MPI_Allreduce(send_buffer, bx, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
   for (int i=0; i<n; ++i) {
     bx[i] = send_buffer[i];
   }
#endif

   if ( m_periodic_or_neumann ) {

      // With periodic or Neumann radial boundary conditions, the resulting tridiagonal
      // matrix is singular, which is problematic for the direct solver.  We
      // therefore solve the system as a quadratic minimization problem
      // with a linear constraint, using a Lagrange multiplier.  The constraint
      // is that the integral of the solution over the radial domain is zero,
      // so we need the flux surface areas to compute the integral.

      const LevelData<FArrayBox>& fs_areas = m_flux_surface.areas();

      for (dit.begin(); dit.ok(); ++dit) {
         Box overlap = fs_grids[dit] & domain_box_fs_intersection;
         if ( overlap.ok() ) {
            const FArrayBox& this_fab = fs_areas[dit];
            for (BoxIterator bit(overlap); bit.ok(); ++bit) {
               IntVect iv = bit();
               send_buffer[iv[RADIAL_DIR]-radial_lo] = this_fab(iv);
            }
         }
      }

      double * weights = new double[n];

#ifdef CH_MPI
      MPI_Allreduce(send_buffer, weights, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      for (int i=0; i<n; ++i) {
         weights[i] = send_buffer[i];
      }
#endif

      double total_weight = 0.;
      for (int i=0; i<n; ++i) total_weight += weights[i];
      for (int i=0; i<n; ++i) weights[i] /= total_weight;

      solveTridiagonalNeumann(n, diag, udiag, ldiag, bx, weights);

      delete [] weights;
   }
   else {
      solveTridiagonal(radially_periodic, n, diag, udiag, ldiag, bx);
   }

   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_fab = a_data[dit];
      BoxIterator bit(fs_grids[dit]);
      for (bit.begin(); bit.ok(); ++bit) {
         IntVect iv = bit();
         this_fab(iv) = bx[iv[RADIAL_DIR]-radial_lo];
      }
   }

   delete [] send_buffer;
   delete [] temp;
}



void
MBTridiagonalSolver::solveTridiagonal( const bool    a_periodic,
                                       const int            a_n,
                                       double *          a_diag,
                                       double *         a_udiag,
                                       double *         a_ldiag,
                                       double *    a_fs_average ) const
{
  int i, info, nrhs, ldb;

  if (!a_periodic) {

    nrhs = 1;
    ldb = a_n;

    // Solve tridiagonal
    dgtsv_(a_n, nrhs, a_ldiag+1, a_diag, a_udiag, a_fs_average, ldb, info);

    if (info != 0) {
       cout << "dgtsv failed, returning " << info << ", n = " << a_n << endl;
    }
  }
  else {

    int nm1 = a_n-1;
    ldb = nm1;

    double a11 = a_diag[0];
    double a21 = a_ldiag[1];
    double am1 = a_udiag[nm1];
    double a12 = a_udiag[0];
    double a1m = a_ldiag[0];
    double b1 = a_fs_average[0];

    nrhs = 2;
    double * rhs = new double[nrhs*nm1];
    double * y = rhs;
    double * z = rhs + nm1;
    for (i=0; i<nm1; i++) {
      y[i] = a_fs_average[i+1];
    }
    z[0] = -a21;
    for (i=1; i<nm1-1; i++) {
      z[i] = 0.;
    }
    z[nm1-1] = -am1;

    // Solve tridiagonals
    dgtsv_(nm1, nrhs, a_ldiag+2, a_diag+1, a_udiag+1, rhs, ldb, info);

    if (info != 0) {
       cout << "dgtsv failed, returning " << info << ", n = " << a_n << endl;
    }

    double denom = a11 + a12*z[0] + a1m*z[nm1-1];
    double beta;

    if (denom != 0.) {
      beta = (b1 - a12*y[0] - a1m*y[nm1-1]) / (a11 + a12*z[0] + a1m*z[nm1-1]);
      }
    else {
      cout << "Zero denominator in tridiagonal solve" << endl;
      exit(1);
    }

    a_fs_average[0] = beta;
    for (i=1; i<a_n; i++) {
      a_fs_average[i] = y[i-1] + beta * z[i-1];
    }

    delete [] rhs;
  }
}



void
MBTridiagonalSolver::solveTridiagonalNeumann( const int            a_n,
                                              double *          a_diag,
                                              double *         a_udiag,
                                              double *         a_ldiag,
                                              double *    a_fs_average,
                                              const double * a_weights ) const
{
  int n_augmented = a_n + 1;
  int num_mat_elems = n_augmented * n_augmented;

  double * A = new double[num_mat_elems];
  for (int k=0; k<num_mat_elems; ++k) A[k] = 0.;

  // Load the diagonal
  for (int i=0, k=0; i<a_n; ++i, k+=n_augmented+1) {
     A[k] = a_diag[i];
  }

  // Load the lower diagonal
  for (int i=1, k=1; i<a_n; ++i, k+=n_augmented+1) {
     A[k] = a_ldiag[i];
  }

  // Load the upper diagonal
  for (int i=0, k=n_augmented; i<a_n-1; ++i, k+=n_augmented+1) {
     A[k] = a_udiag[i];
  }

  // Load the right-hand border
  for (int k=n_augmented*a_n, i=0; k<num_mat_elems-1; ++k, ++i) {
     A[k] = a_weights[i];
  }

  for (int k=a_n, i=0; k<n_augmented*a_n; k+=n_augmented, ++i) {
     A[k] = a_weights[i];
  }
  // Load the lower border

  // Create the augmented rhs
  int nrhs = 1;
  double * b = new double[n_augmented];
  for (int k=0; k<a_n; ++k) b[k] = a_fs_average[k];
  b[a_n] = 0.;

  // Remove the average to make well-posed
  double b_av = 0.;
  for (int k=0; k<a_n; ++k) {
     b_av += a_weights[k]*b[k];
  }
  for (int k=0; k<a_n; ++k) b[k] -= b_av;

  int * ipvt = new int[n_augmented];
  int info;

  dgesv_( n_augmented, nrhs, A, n_augmented, ipvt, b, n_augmented, info);

  if (info != 0) {
     cout << "dgesv failed, returning " << info << endl;
  }

  for (int i=0; i<a_n; ++i) a_fs_average[i] = b[i];

  //     cout << "      Lagrange multiplier = " << b[a_n] << endl;

  delete [] ipvt;
  delete [] b;
  delete [] A;
}



bool
MBTridiagonalSolver::isCoreRadialPeriodicOrNeumannBC( const EllipticOpBC& a_bc ) const
{
   bool flag = false;

   if ( typeid(a_bc) == typeid(LogRectEllipticOpBC) ) {
      const MagGeom& geom = (MagGeom&)m_geometry;
      const ProblemDomain& domain = (geom.getBlockCoordSys(0)).domain();

      if ( domain.isPeriodic(RADIAL_DIR) ||
           (a_bc.getBCType(0, RADIAL_DIR, 0) == EllipticOpBC::NEUMANN &&
            a_bc.getBCType(0, RADIAL_DIR, 1) == EllipticOpBC::NEUMANN) ) {
         flag = true;
      }
   }

   if ( typeid(a_bc) == typeid(SNCoreEllipticOpBC) ) {
      if ( a_bc.getBCType(0, RADIAL_DIR, 0) == EllipticOpBC::NEUMANN &&
           a_bc.getBCType(0, RADIAL_DIR, 1) == EllipticOpBC::NEUMANN ) {
         flag = true;
      }
   }

   return flag;
}


#include "NamespaceFooter.H"
