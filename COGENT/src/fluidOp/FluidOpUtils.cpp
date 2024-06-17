#include "FluidOpUtils.H"
#include "Directions.H"

#include "FourthOrderUtil.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"

#include "FluidOpF_F.H"

#include "NamespaceHeader.H" 



FluidOpUtils::FluidOpUtils(const ParmParse&  a_pp,
                           const MagGeom&    a_geometry,
                           const double      a_larmor,
                           const int         a_verbosity )
   : m_verbosity(a_verbosity),
     m_geometry(a_geometry),
     m_larmor(a_larmor),
     m_flux_surface(a_geometry),
     m_relaxation_coefficient(0.),
     m_hyperviscosity_op(NULL)

{
   parseParameters( a_pp );
   if (m_verbosity>0) {
      printParameters();
   }

   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   
   m_volume.define(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(m_volume);
   
}


FluidOpUtils::~FluidOpUtils()
{
   if (m_hyperviscosity_op) delete m_hyperviscosity_op;
}


void FluidOpUtils::initializeHyperViscosityOp(const ParmParse& a_pp,
                                              EllipticOpBC&    a_bc,
                                              bool             a_use_hyperviscosity_bcs)
{
   
   bool include_parallel_hyperviscosity(false);
   bool include_perpendicular_hyperviscosity(false);

   m_use_hyperviscosity_bcs = a_use_hyperviscosity_bcs;
   m_mapped_hyperviscosity = false;
   m_hyperviscosity_order = 4;
   
   GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
   std::string grid_function_name;
   
   a_pp.query("mapped_hyperviscosity", m_mapped_hyperviscosity);
   a_pp.query("hyperviscosity_order", m_hyperviscosity_order);
   
   // Get parallel coefficient
   RefCountedPtr<GridFunction> D_par_func;
   if (a_pp.contains("D_par")) {
      include_parallel_hyperviscosity = true;
      a_pp.get("D_par", grid_function_name );
      D_par_func = grid_library->find( grid_function_name );
   }

   // Get perpendicular coefficient
   RefCountedPtr<GridFunction> D_perp_func;
   if (a_pp.contains("D_perp")) {
      include_perpendicular_hyperviscosity = true;
      a_pp.get("D_perp", grid_function_name );
      D_perp_func = grid_library->find( grid_function_name );
   }
   
   if (!(include_parallel_hyperviscosity || include_perpendicular_hyperviscosity)) {
      MayDay::Error("FluidOpUtils::either parallel or perpendicular hyperviscosity coefficient must be specified");
   }
   
   
   // Get distjoint box layout
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   // Get geometry order
   int order = (m_geometry.secondOrder()) ? 2 : 4;

   // Get parallel coefficient profile function
   LevelData<FluxBox> D_par_fc(grids, 1, IntVect::Zero);
   if (include_parallel_hyperviscosity) {
     LevelData<FArrayBox> D_par_cc(grids, 1, 2*IntVect::Unit);
     D_par_func->assign( D_par_cc, m_geometry, 0.0);
     m_geometry.fillInternalGhosts( D_par_cc );
     SpaceUtils::interpToFaces(D_par_fc, D_par_cc, order);
   }

   // Get perpendicular coefficient profile function
   LevelData<FluxBox> D_perp_fc(grids, 1, IntVect::Zero);
   if (include_perpendicular_hyperviscosity) {
     LevelData<FArrayBox> D_perp_cc(grids, 1, 2*IntVect::Unit);
     D_perp_func->assign( D_perp_cc, m_geometry, 0.0);
     m_geometry.fillInternalGhosts( D_perp_cc );
     SpaceUtils::interpToFaces(D_perp_fc, D_perp_cc, order);
   }
   
   // Get parallel elliptic coefficients
   const LevelData<FluxBox>& par_coeff = m_geometry.getEllipticOpParCoeff();
   const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();

   // Get perpendicular elliptic coefficients
   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped  = m_geometry.getEllipticOpPerpCoeffMapped();
   
   // Compute diffusion coefficient
   LevelData<FluxBox> D_tensor(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   LevelData<FluxBox> D_tensor_mapped(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const RealVect& dx= coord_sys.getMappedCellSize();

      RealVect mapped_fac(dx);
      mapped_fac *= dx;
      
      D_tensor[dit].setVal(0.);
      D_tensor_mapped[dit].setVal(0.);

      // Add parallel contribution
      if (include_parallel_hyperviscosity) {
         FluxBox par_coeff_tensor(par_coeff[dit].box(),SpaceDim*SpaceDim);
         FluxBox par_coeff_tensor_mapped(par_coeff_mapped[dit].box(),SpaceDim*SpaceDim);
         if (!m_mapped_hyperviscosity) {
            par_coeff_tensor.copy(par_coeff[dit]);
            par_coeff_tensor_mapped.copy(par_coeff_mapped[dit]);
         }
         else {
            par_coeff_tensor.setVal(0.);
            par_coeff_tensor_mapped.setVal(0.);
         }
         for (int dir = 0; dir < SpaceDim; dir++) {
            for (int n = 0; n < SpaceDim*SpaceDim; n++) {
               if (!m_mapped_hyperviscosity) {
                  par_coeff_tensor[dir].mult(D_par_fc[dit][dir],0,n,1);
                  par_coeff_tensor_mapped[dir].mult(D_par_fc[dit][dir],0,n,1);
               }
               else {
                  if ((SpaceDim == 2 && n == 3) || (SpaceDim == 3 && n == 4)) {
                     par_coeff_tensor_mapped[dir].copy(D_par_fc[dit][dir],0,n,1);
                     par_coeff_tensor_mapped[dir].mult(mapped_fac[dir],n,1);
                  }
               }
            }
            D_tensor[dit][dir].plus(par_coeff_tensor[dir]);
            D_tensor_mapped[dit][dir].plus(par_coeff_tensor_mapped[dir]);
         }
      }
      
      // Add perpendicular contribution
      if (include_perpendicular_hyperviscosity) {
         FluxBox perp_coeff_tensor(perp_coeff[dit].box(),SpaceDim*SpaceDim);
         FluxBox perp_coeff_tensor_mapped(perp_coeff_mapped[dit].box(),SpaceDim*SpaceDim);
         if (!m_mapped_hyperviscosity) {
            perp_coeff_tensor.copy(perp_coeff[dit]);
            perp_coeff_tensor_mapped.copy(perp_coeff_mapped[dit]);
         }
         else {
            perp_coeff_tensor.setVal(0.);
            perp_coeff_tensor_mapped.setVal(0.);
         }
         for (int dir = 0; dir < SpaceDim; dir++) {
            for (int n = 0; n < SpaceDim*SpaceDim; n++) {
               if (!m_mapped_hyperviscosity) {
                  perp_coeff_tensor[dir].mult(D_perp_fc[dit][dir],0,n,1);
                  perp_coeff_tensor_mapped[dir].mult(D_perp_fc[dit][dir],0,n,1);
               }
               else {
                  if (n == 0 || n == 8) {
                     perp_coeff_tensor_mapped[dir].copy(D_perp_fc[dit][dir],0,n,1);
                     perp_coeff_tensor_mapped[dir].mult(mapped_fac[dir],n,1);
                  }
               }
            }
            D_tensor[dit][dir].plus(perp_coeff_tensor[dir]);
            D_tensor_mapped[dit][dir].plus(perp_coeff_tensor_mapped[dir]);
         }
      }
   }
   
   m_hyperviscosity_op = new Diffusion(a_pp, m_geometry);
   m_hyperviscosity_op->setOperatorCoefficients(D_tensor, D_tensor_mapped, a_bc);
   
   //Check that low_pollution option is used for mapped_hyperviscosity
   if (!(m_hyperviscosity_op->lowPollution()) && m_mapped_hyperviscosity) {
      MayDay::Error("FluidOpUtils::mapped hyperviscosity requires the low_pollution option");
   }
}

void FluidOpUtils::addHyperViscosity(LevelData<FArrayBox>&        a_rhs,
                                     const LevelData<FArrayBox>&  a_soln,
                                     const bool                   a_is_mapped)
{
   
   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   
   LevelData<FArrayBox> hypervisc(grids, 1, IntVect::Zero);
   
   // Compute 2nd order HV
   m_hyperviscosity_op->computeFluxDivergence( a_soln, hypervisc, false, !m_use_hyperviscosity_bcs);
   if (m_mapped_hyperviscosity) m_geometry.multJonValid(hypervisc);

   // Compute 4th order HV
   if (m_hyperviscosity_order == 4) {
      LevelData<FArrayBox> tmp(grids, 1, IntVect::Zero);
      for (DataIterator dit(grids); dit.ok(); ++dit) {
          tmp[dit].copy(hypervisc[dit]);
      }
      m_hyperviscosity_op->computeFluxDivergence( tmp, hypervisc, false, !m_use_hyperviscosity_bcs);
      if (m_mapped_hyperviscosity) m_geometry.multJonValid(hypervisc);
   }
   
   if (a_is_mapped) {
      m_geometry.multJonValid(hypervisc);
   }
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
       a_rhs[dit] -= hypervisc[dit];
   }

}


void FluidOpUtils::initializeBoundaryBuffers(const ParmParse& a_pp,
                                             const bool       a_include_full_relaxation)
{
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
   std::string grid_function_name;
   
   // Get spatial profile for boundary buffer
   RefCountedPtr<GridFunction> boundary_buffer_profile_func;
   if (a_pp.contains("boundary_buffer_profile")) {
       a_pp.get("boundary_buffer_profile", grid_function_name );
       boundary_buffer_profile_func = grid_library->find( grid_function_name );
   }
   else {
      MayDay::Error("FluidOpUtils::boundary_buffer_profile must be specified");
   }
   
   
   m_boundary_buffer_profile.define(grids, 1, IntVect::Zero);
   boundary_buffer_profile_func->assign(m_boundary_buffer_profile,
                                        m_geometry,
                                        0.0);
   
   
   // Get reference function and relaxation coefficient

   if (a_include_full_relaxation) {

      RefCountedPtr<GridFunction> ref_soln_func;

      if (a_pp.contains("reference_function")) {
         a_pp.get("reference_function", grid_function_name );
         ref_soln_func = grid_library->find( grid_function_name );
      }
      else {
         MayDay::Error("FluidOpUtils::reference_function must be specified");
      }
      
      if (a_pp.contains("relaxation_coefficient")) {
         a_pp.get("relaxation_coefficient", m_relaxation_coefficient );
      }
      else {
         MayDay::Error("FluidOpUtils::relaxation_coefficient must be specified");
      }

      m_ref_soln.define(grids, 1, IntVect::Zero);
      ref_soln_func->assign( m_ref_soln, m_geometry, 0.0);
   }
}


void FluidOpUtils::suppressNonZonalCompInBoundaryBuffers(LevelData<FArrayBox>& a_data)
{
   /*
    This function suppresses non zonal component in
    the boundary buffers
    */
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
      
   LevelData<FArrayBox> data_fs(grids, 1, IntVect::Zero);
   m_flux_surface.averageAndSpread(a_data, data_fs);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     FArrayBox& this_data = a_data[dit];
     FArrayBox& this_data_fs = data_fs[dit];
     
     FArrayBox this_data_osc(grids[dit],1);
     this_data_osc.copy(this_data);
     this_data_osc.minus(this_data_fs);
     
     FArrayBox this_data_osc_pp(grids[dit],1);
     this_data_osc_pp.copy(this_data_osc);
     this_data_osc_pp.mult(m_boundary_buffer_profile[dit]);
     this_data_osc_pp.negate();
     this_data_osc_pp.plus(this_data_osc);
     
     this_data.copy(this_data_fs);
     this_data.plus(this_data_osc_pp);
   }
}

void FluidOpUtils::addRelaxationInBoundaryBuffers(LevelData<FArrayBox>&       a_rhs,
                                                  const LevelData<FArrayBox>& a_soln_phys,
                                                  const bool                  a_is_mapped)
{
   /*
    Adds relaxation to reference solution in the boundary buffers
    */
   
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   LevelData<FArrayBox> relaxation(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      relaxation[dit].copy(a_soln_phys[dit]);
      relaxation[dit].minus(m_ref_soln[dit]);
      relaxation[dit].mult(m_boundary_buffer_profile[dit]);
      relaxation[dit].mult(m_relaxation_coefficient);
   }

   if (a_is_mapped) {
      m_geometry.multJonValid(relaxation);
   }
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_rhs[dit] -= relaxation[dit];
   }
}

void FluidOpUtils::computeExBAdvection(LevelData<FArrayBox>& a_ExB_advection,
                                       const LevelData<FArrayBox>& a_soln,
                                       const LevelData<FluxBox>& a_E_field,
                                       const std::string& a_advection_scheme,
                                       const bool a_extrap_to_ghosts)
{
  /*
    Computes physical div(ExB/B^2 * a_soln);
    Has some 4th order elements, but not fully 4th order 
   */
  
   // Get geometry parameters
   const DisjointBoxLayout& grids = m_geometry.grids();
   bool fourth_order = (m_geometry.secondOrder()) ? false : true;
   int order = (m_geometry.secondOrder()) ? 2 : 4;
   
   // Compute the ExB drift on face centers
   LevelData<FluxBox> ExB_drift(grids, 3, IntVect::Unit);
   m_geometry.computeEXBDrift(a_E_field, ExB_drift);

   IntVect ghosts;
   if (a_soln.ghostVect() >= 2*IntVect::Unit) {
     ghosts = a_soln.ghostVect();
   }
   else {
     if (m_geometry.secondOrder()) {
       ghosts = 2*IntVect::Unit;
     }
     else {
       ghosts = 3*IntVect::Unit;
     }
   }
   LevelData<FArrayBox> soln_with_ghosts(grids, 1, ghosts);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     soln_with_ghosts[dit].copy(a_soln[dit]);
   }

   if (a_extrap_to_ghosts) {
     if (m_geometry.secondOrder()) {
       // Use 4th order here, because we need at least two ghost cells
       // for most of our advection schemes (e.g., uw3)
       m_geometry.extrapolateToPhysicalGhosts(soln_with_ghosts, true);
     }
     else {
       // For 4th order we need information in corner ghosts.
       // Presently use second-order extrapolation into all ghosts
       m_geometry.fillInternalGhosts(soln_with_ghosts);
       m_geometry.extrapolateAtPhysicalBoundaries(soln_with_ghosts, 4, 3);
     }
   }
   
   LevelData<FluxBox> soln_face(grids, 1, IntVect::Unit);
   if (a_advection_scheme == "c2") {
      SpaceUtils::interpToFaces(soln_face, soln_with_ghosts, order);
   }
   
   else {
     LevelData<FluxBox> ExB_drift_mapped(grids, SpaceDim, IntVect::Unit);
     for (DataIterator dit(ExB_drift_mapped.dataIterator()); dit.ok(); ++dit) {
       if (SpaceDim == 3) {
         ExB_drift_mapped[dit].copy(ExB_drift[dit]);
       }
       if (SpaceDim == 2) {
         ExB_drift_mapped[dit].copy(ExB_drift[dit],0,0,1);
         ExB_drift_mapped[dit].copy(ExB_drift[dit],2,1,1);
       }
     }
     m_geometry.multNTransposePointwise(ExB_drift_mapped);
   
     SpaceUtils::upWindToFaces(soln_face,
			       soln_with_ghosts,
			       ExB_drift_mapped,
			       a_advection_scheme,
			       order);
   }

   LevelData<FluxBox> flux_ExB(grids, SpaceDim, IntVect::Unit);
   for (DataIterator dit(flux_ExB.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         for (int nComp=0; nComp<SpaceDim; nComp++) {
            flux_ExB[dit][dir].copy(soln_face[dit][dir],0,nComp);
         }
         flux_ExB[dit][dir].mult(ExB_drift[dit][dir],0,0);
         if (SpaceDim == 2) {
            flux_ExB[dit][dir].mult(ExB_drift[dit][dir],2,1);
         }
         if (SpaceDim == 3) {
            flux_ExB[dit][dir].mult(ExB_drift[dit][dir],1,1);
            flux_ExB[dit][dir].mult(ExB_drift[dit][dir],2,2);
         }
         flux_ExB[dit][dir] *= m_larmor;
      }
   }
   
   m_geometry.applyAxisymmetricCorrection(flux_ExB);
   
   m_geometry.computeMappedGridDivergence(flux_ExB, a_ExB_advection, fourth_order);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_ExB_advection[dit] /= m_volume[dit];
   }
}

void FluidOpUtils::computeExBGradData(LevelData<FArrayBox>&       a_ExB_grad_data,
                                      const LevelData<FArrayBox>& a_data,
                                      const LevelData<FArrayBox>& a_E_field,
                                      const std::string&          a_advection_scheme,
                                      const bool                  a_extrap_to_ghosts)
{
  /*
    Computes physical ExB/B^2 * \nabla(a_data);
    Has some 4th order elements, but not fully 4th order
   */
  
   // Get geometry parameters
   const DisjointBoxLayout& grids = m_geometry.grids();
   bool fourth_order = (m_geometry.secondOrder()) ? false : true;
   int order = (m_geometry.secondOrder()) ? 2 : 4;
   
   // Compute the ExB drift on face centers
   LevelData<FArrayBox> ExB_drift(grids, 3, IntVect::Zero);
   m_geometry.computeEXBDrift(a_E_field, ExB_drift);

   IntVect ghosts;
   if (a_data.ghostVect() >= 2*IntVect::Unit) {
     ghosts = a_data.ghostVect();
   }
   else {
     ghosts = 2*IntVect::Unit;
   }
   LevelData<FArrayBox> data_with_ghosts(grids, 1, ghosts);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     data_with_ghosts[dit].copy(a_data[dit]);
   }

   if (a_extrap_to_ghosts) {
     // Use 4th order here, because we need at least two ghost cells
     // for most of our advection schemes (e.g., uw3)
     m_geometry.extrapolateToPhysicalGhosts(data_with_ghosts, true);
   }
   
   LevelData<FArrayBox> ExB_drift_tmp(grids, SpaceDim, IntVect::Zero);
   for (DataIterator dit(ExB_drift_tmp.dataIterator()); dit.ok(); ++dit) {
      if (SpaceDim == 3) {
         ExB_drift_tmp[dit].copy(ExB_drift[dit]);
      }
      if (SpaceDim == 2) {
         ExB_drift_tmp[dit].copy(ExB_drift[dit],0,0,1);
         ExB_drift_tmp[dit].copy(ExB_drift[dit],2,1,1);
      }
   }
   LevelData<FArrayBox> ExB_drift_mapped(grids, SpaceDim, IntVect::Zero);
   m_geometry.multiplyNTranspose(ExB_drift_mapped,
                                 ExB_drift_tmp);
   
   LevelData<FArrayBox> grad_data_mapped(grids, SpaceDim, IntVect::Zero);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      const MagBlockCoordSys& coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const RealVect& dx= coord_sys.getMappedCellSize();
      
      SpaceUtils::upWindMappedGradient(grad_data_mapped[dit],
                                       data_with_ghosts[dit],
                                       ExB_drift_mapped[dit],
                                       dx,
                                       a_advection_scheme);
   }

   LevelData<FArrayBox> grad_data(grids, SpaceDim, IntVect::Zero);
   if (SpaceDim == 3) {
      m_geometry.unmapGradient(grad_data_mapped, grad_data );
   }
   if (SpaceDim == 2) {
      m_geometry.unmapPoloidalGradient(grad_data_mapped, grad_data );
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_ExB_grad_data[dit].setVal(0.);
      FArrayBox tmp(grids[dit],1);
      
      for (int n=0; n<SpaceDim; n++) {
         tmp.copy(ExB_drift_tmp[dit],n,0,1);
         tmp.mult(grad_data[dit],n,0,1);
         a_ExB_grad_data[dit].plus(tmp);
      }
      
      a_ExB_grad_data[dit].mult(m_larmor);
   }
}

void
FluidOpUtils::subtractFSaverage(LevelData<FArrayBox>&  a_data)
{
   
   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   LevelData<FArrayBox> data_fs_averaged(grids, a_data.nComp(), a_data.ghostVect());
   m_flux_surface.averageAndSpread(a_data, data_fs_averaged);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
     a_data[dit] -= data_fs_averaged[dit];
   }
}


Real
FluidOpUtils::computeAdvectionDt(const LevelData<FluxBox>&  a_faceVel_mapped,
                                 const std::string&         a_flux_type,
                                 bool                       a_time_step_diagnostics)
{
   
   // Get velocities normal to cell faces and compute area-weighted normal
   // velocities -- 2nd order should be good enough.
   // Note: For axisymmetric 2-spatial-D configuration geometries, the area-weighted
   // velocities returned by the following function call contain a factor
   // of 2piR times major radius R; a consequence of the fact that this
   // factor is included in all of the metric factors.  However, the
   // physical cell volume being divided out below also contains
   // a 2piR factor (having been derived from the metrics), and therefore
   // this extra factor has no net effect (to second order).
   // (Of course this is irrrelevant in 3D)

   const MagCoordSys& coord_sys = *(m_geometry.getCoordSys());
   const DisjointBoxLayout& grids( m_geometry.grids() );

   LevelData<FluxBox> faceNormalVel(grids, 1);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys( m_geometry.getBlockCoordSys( grids[dit] ) );
      RealVect face_area = block_coord_sys.getMappedFaceArea();

      const FluxBox& thisFaceVel = a_faceVel_mapped[dit];
      FluxBox& thisNormalVel = faceNormalVel[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& thisNormalVel_dir = thisNormalVel[dir];
         thisNormalVel_dir.copy(thisFaceVel[dir],dir,0,1);
         thisNormalVel_dir *= face_area[dir];
      }
   }

   // now average back to cell-centers and divide by cell volumes.
   // instead of averaging face->cell we pick
   // the max absolute value on the two opposing faces.

   LevelData<FArrayBox> cellVolumes(grids, 1, IntVect::Zero);
   m_geometry.getCellVolumes(cellVolumes);

   struct {
      double val;
      int rank;
   } pair_in, pair_out;

   Real maxVelLoc = -1.;
   int maxblockLoc(-1);
   IntVect maxindLoc;
   int maxDirLoc(-1);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coord_sys.whichBlock(grids[dit]);

      FArrayBox cellVel(grids[dit], 1);
      FArrayBox cellVelDir(grids[dit],SpaceDim);
      cellVel.setVal(0.);
      for (int dir=0; dir<SpaceDim; dir++) {
         // average face velocity on the two faces in this direction
         int faceComp = 0;
         int cellComp = dir;
         EdgeToCell(faceNormalVel[dit], faceComp,
                    cellVelDir, cellComp,
                    dir);
         cellVelDir.abs(dir,1);
         cellVel.plus(cellVelDir, dir, 0, 1);
      }

      // Divide by the cell physical volume, which is J times the computational cell volume
      cellVel.divide(cellVolumes[dit], 0, 0, 1);

      // now compute maxVelLoc on this box
      // note that this is essentially an estimate of max(vel/dx)
      Real thisMax = cellVel.norm(0,0,1);
      if (thisMax > maxVelLoc) {
         maxVelLoc = thisMax;
         if (a_time_step_diagnostics) {
            maxblockLoc = block_number;
            maxindLoc = cellVel.maxIndex();

            // Figure out which direction made the biggest contribution
            double max_cv = 0.;
            for (int dir=0; dir<SpaceDim; ++dir) {
               if (cellVelDir(maxindLoc,dir) > max_cv) {
                  max_cv = cellVelDir(maxindLoc,dir);
                  maxDirLoc = dir;
               }
            }
         }
      }
   }

   RealVect X;
   if (a_time_step_diagnostics) {
      
     if (maxblockLoc < 0) MayDay::Error( "FluidOpUtils::computeAdvectionDt: time step diagnostic calculation failed, check if the number of processes used is consistent with the domain decomposition " );
     
      RealVect dx = coord_sys.getCoordSys(maxblockLoc)->dx();
      RealVect xi = dx*maxindLoc;
      xi += 0.5*dx;
      X = coord_sys.getCoordSys(maxblockLoc)->realCoord(xi);
   }

   Real maxVel = maxVelLoc;
   IntVect maxind = maxindLoc;
   int maxDir = maxDirLoc;
#ifdef CH_MPI
   if (a_time_step_diagnostics) {
      pair_in.val = maxVel;
      pair_in.rank = procID();
      MPI_Allreduce(&pair_in, &pair_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      maxVel = pair_out.val;
      MPI_Bcast(X.dataPtr(), SpaceDim, MPI_DOUBLE, pair_out.rank, MPI_COMM_WORLD);
      MPI_Bcast(maxind.dataPtr(), SpaceDim, MPI_INT, pair_out.rank, MPI_COMM_WORLD);
      MPI_Bcast(&maxDir, 1, MPI_INT, pair_out.rank, MPI_COMM_WORLD);
   }
   else {
      MPI_Allreduce(&maxVelLoc, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   }
#endif

   Real smallVal = 1.0e-15;
   Real newDt = 0;
   Real dt_dim_factor = sqrt(SpaceDim);
   Real stability_bound;
   
   if (maxVel > smallVal) {
      if (a_flux_type.compare("c2")==0) {
         stability_bound = 2.8;
      }
      else if (a_flux_type.compare("uw1")==0) {
         stability_bound = 2.7852;
      }
      else if (a_flux_type.compare("uw3")==0) {
         stability_bound = 1.7453;
      }
      else if (a_flux_type.compare("uw5")==0) {
         stability_bound = 1.7320;
      }
      else if (a_flux_type.compare("weno5")==0) {
         stability_bound = 1.7320;
      }
      else if (a_flux_type.compare("bweno")==0) {
         stability_bound = 1.7453;
      }
      else {
         MayDay::Error("FluidOpUtils::computeAdvectionDt: unknown flux type");
      }
   
      newDt = stability_bound / (maxVel * dt_dim_factor);
   }

   if (a_time_step_diagnostics && procID()==0) {
#if CFG_DIM==2
      cout << "Vlasov operator time step ("<< newDt <<") was limited by the velocity at (R,Z) = " << X << " and mapped coordinate = " << maxind << endl;
#endif
#if CFG_DIM==3
      cout << "Vlasov operator time step ("<< newDt <<") was limited by the velocity at (X,Y,Z) = " << X << " and mapped coordinate = " << maxind << endl;
#endif
      cout << "The ";
      switch(maxDir)
         {
         case RADIAL_DIR:
            cout << "radial";
            break;
         case POLOIDAL_DIR:
            cout << "poloidal";
            break;
#if CFG_DIM==3
         case TOROIDAL_DIR:
            cout << "toroidal";
            break;
#endif
         }
      cout << " direction makes the largest contribution to the stable dt reciprocal at that point" << endl;
   }

   return newDt;
}


void FluidOpUtils::parseParameters( const ParmParse& a_pp )
{

}


void FluidOpUtils::printParameters()
{

}

#include "NamespaceFooter.H"
