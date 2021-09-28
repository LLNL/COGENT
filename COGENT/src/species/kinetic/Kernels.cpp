#include "Kernels.H"
#include "Directions.H"
#include "PhaseBlockCoordSys.H"
#include "FourthOrderUtil.H"
#include "EdgeToCell.H"
#include "CellToEdge.H"
#include "ConstFact.H"

#include "KernelsF_F.H"

#include "NamespaceHeader.H"


template <> void
Kernel<FArrayBox>::computeVelCfgComp(LevelData<FArrayBox>&      a_vel_cfg_comp,
                                     const PhaseGeom&           a_phase_geom,
                                     const LevelData<FluxBox>&  a_field,
                                     const int                  a_option     ) const
{
   /*
    Computes cell-centered physical components of GK velocity in configuration directions
    */

   const DisjointBoxLayout & grids = a_vel_cfg_comp.getBoxes();
   LevelData<FluxBox> pointwiseFaceVel(grids, SpaceDim, IntVect::Unit);
   a_phase_geom.computeGKVelocities(a_field, pointwiseFaceVel, false, a_option);

   for (DataIterator dit(a_vel_cfg_comp.dataIterator()); dit.ok(); ++dit) {
     a_vel_cfg_comp[dit].setVal(0.);
     for (int dir=0; dir<CFG_DIM; dir++) {
          FORT_COMPUTE_VEL_CELL(CHF_FRA(a_vel_cfg_comp[dit]),
                                CHF_FRA(pointwiseFaceVel[dit][dir]),
                                CHF_CONST_INT(dir),
                                CHF_BOX(a_vel_cfg_comp[dit].box()));
     }
   }
}


template <> void
Kernel<FluxBox>::computeVelCfgComp(LevelData<FluxBox>&        a_vel_cfg_comp,
                                   const PhaseGeom&           a_phase_geom,
                                   const LevelData<FluxBox>&  a_field,
                                   const int                  a_option     ) const
{
   /*
    Computes face-centered physical components of GK velocity in configuration directions
    */

   const DisjointBoxLayout & grids = a_vel_cfg_comp.getBoxes();
   LevelData<FluxBox> pointwiseFaceVel(grids, SpaceDim, IntVect::Unit);
   a_phase_geom.computeGKVelocities(a_field, pointwiseFaceVel, false, a_option);

   for (DataIterator dit(a_vel_cfg_comp.dataIterator()); dit.ok(); ++dit) {
     a_vel_cfg_comp[dit].setVal(0.);
     for (int dir=0; dir<CFG_DIM; dir++) {
        a_vel_cfg_comp[dit][dir].copy(pointwiseFaceVel[dit][dir]);
     }
   }
}


template <> void
Kernel<FArrayBox>::computeVelCfgCompNormals(LevelData<FArrayBox>&      a_vel_cfg_comp,
                                            const PhaseGeom&           a_phase_geom,
                                            const LevelData<FluxBox>&  a_field,
                                            const int                  a_option     ) const
{
   /*
    Computes cell-centered normal components of GK velocity
    (i.e., projection of physical velocity onto normal directions)
    in configuration directions
    */

   const DisjointBoxLayout & grids = a_vel_cfg_comp.getBoxes();
   LevelData<FluxBox> pointwiseFaceVel(grids, SpaceDim, IntVect::Unit);
   a_phase_geom.computeGKVelocities(a_field, pointwiseFaceVel, false, a_option);
   a_phase_geom.multNTransposePointwise( pointwiseFaceVel );
   
   const CFG::MagGeom& mag_geom( a_phase_geom.magGeom() );
   const CFG::DisjointBoxLayout& grids_cfg = mag_geom.grids();
   CFG::LevelData<CFG::FluxBox> face_areas(grids_cfg, 1, CFG::IntVect::Unit);
   mag_geom.getPointwiseFaceAreas(face_areas);

   LevelData<FluxBox> face_areas_inj;
   a_phase_geom.injectConfigurationToPhase(face_areas, face_areas_inj);
   
   for (DataIterator dit(a_vel_cfg_comp.dataIterator()); dit.ok(); ++dit) {
     for (int dir=0; dir<CFG_DIM; dir++) {
          
         const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
         block_coord_sys.divideInjectedData(face_areas_inj[dit][dir],
                                            pointwiseFaceVel[dit][dir]);

         RealVect dx = block_coord_sys.dx();
         Real mapped_volume = block_coord_sys.getMappedCellVolume();
         Real mapped_volume_cfg = mapped_volume / dx[VPARALLEL_DIR] / dx[MU_DIR];
         Real mapped_face_area_cfg = mapped_volume_cfg / dx[dir];
         pointwiseFaceVel[dit][dir].mult(mapped_face_area_cfg);
        
         FORT_COMPUTE_VEL_CELL_NORMALS(CHF_FRA1(a_vel_cfg_comp[dit],dir),
                                       CHF_FRA(pointwiseFaceVel[dit][dir]),
                                       CHF_CONST_INT(dir),
                                       CHF_BOX(a_vel_cfg_comp[dit].box()));
     }
   }
}


template <> void
Kernel<FluxBox>::computeVelCfgCompNormals(LevelData<FluxBox>&        a_vel_cfg_comp,
                                          const PhaseGeom&           a_phase_geom,
                                          const LevelData<FluxBox>&  a_field,
                                          const int                  a_option     ) const
{
   /*
    Computes face-centered normal components of GK velocity
    (i.e., projection of physical velocity onto normal directions)
    in configuration directions
    */

   const DisjointBoxLayout & grids = a_vel_cfg_comp.getBoxes();
   LevelData<FluxBox> pointwiseFaceVel(grids, SpaceDim, IntVect::Unit);
   a_phase_geom.computeGKVelocities(a_field, pointwiseFaceVel, false, a_option);
   a_phase_geom.multNTransposePointwise( pointwiseFaceVel );
   
   const CFG::MagGeom& mag_geom( a_phase_geom.magGeom() );
   const CFG::DisjointBoxLayout& grids_cfg = mag_geom.grids();
   CFG::LevelData<CFG::FluxBox> face_areas(grids_cfg, 1, CFG::IntVect::Unit);
   mag_geom.getPointwiseFaceAreas(face_areas);

   LevelData<FluxBox> face_areas_inj;
   a_phase_geom.injectConfigurationToPhase(face_areas, face_areas_inj);
   
   for (DataIterator dit(a_vel_cfg_comp.dataIterator()); dit.ok(); ++dit) {
     for (int dir=0; dir<CFG_DIM; dir++) {
          
         const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
         block_coord_sys.divideInjectedData(face_areas_inj[dit][dir],
                                            pointwiseFaceVel[dit][dir]);

         RealVect dx = block_coord_sys.dx();
         Real mapped_volume = block_coord_sys.getMappedCellVolume();
         Real mapped_volume_cfg = mapped_volume / dx[VPARALLEL_DIR] / dx[MU_DIR];
         Real mapped_face_area_cfg = mapped_volume_cfg / dx[dir];
         pointwiseFaceVel[dit][dir].mult(mapped_face_area_cfg);

         a_vel_cfg_comp[dit][dir].copy(pointwiseFaceVel[dit][dir]);
     }
   }
}


template <> void
DensityKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                const KineticSpecies& a_kinetic_species ) const
{
   /*
     Nothing to do.  a_result has already been initialized by the distribution function,
     and the moment kernel is unity.
   */
}


template <> void
DensityKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                              const KineticSpecies& a_kinetic_species ) const
{
   /*
     Nothing to do.  a_result has already been initialized by the distribution function,
     and the moment kernel is unity.
   */
}


template <> void
ChargeDensityKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                      const KineticSpecies& a_kinetic_species ) const
{
   /*
     Nothing to do.  a_result has already been initialized by the distribution function,
     and the moment kernel is unity.  Multiplication by the charge (the scale factor
     for this Kernel class) is performed by the MomentOp after integrating over velocity.
   */
}


template <> void
ChargeDensityKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                    const KineticSpecies& a_kinetic_species ) const
{
   /*
     Nothing to do.  a_result has already been initialized by the distribution function,
     and the moment kernel is unity.  Multiplication by the charge (the scale factor
     for this Kernel class) is performed by the MomentOp after integrating over velocity.
   */
}


template <> void
MassDensityKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                    const KineticSpecies& a_kinetic_species ) const
{
   /*
     Nothing to do.  a_result has already been initialized by the distribution function,
     and the moment kernel is unity.  Multiplication by the mass (the scale factor
     for this Kernel class) is performed by the MomentOp after integrating over velocity.
   */
}


template <> void
MassDensityKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                  const KineticSpecies& a_kinetic_species ) const
{
   /*
     Nothing to do.  a_result has already been initialized by the distribution function,
     and the moment kernel is unity.  Multiplication by the mass (the scale factor
     for this Kernel class) is performed by the MomentOp after integrating over velocity.
   */
}


template <> void
MomentumDensityKernel<FArrayBox>::eval( LevelData<FArrayBox>&  a_result,
                                        const KineticSpecies&  a_kinetic_species ) const
{
  const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
  const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
  const DisjointBoxLayout & grids = dfn.getBoxes();

  //Compute cell-centered velocity configurational components                                                                                                                                                     
  LevelData<FArrayBox> cellVelCFG(grids, CFG_DIM, IntVect::Unit);
  computeVelCfgComp( cellVelCFG, phase_geom, m_field, PhaseGeom::DIAGNOSTICS);

  DataIterator dit = dfn.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
     for (int dir=0; dir<CFG_DIM; ++dir) {
        a_result[dit].copy(dfn[dit],0,dir,1);
     }
     a_result[dit].mult(cellVelCFG[dit], 0, 0, CFG_DIM);
  }
}


template <> void
MomentumDensityKernel<FluxBox>::eval( LevelData<FluxBox>&    a_result,
                                      const KineticSpecies&  a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const DisjointBoxLayout & grids = dfn.getBoxes();

   LevelData<FluxBox> dfn_face(grids, 1, dfn.ghostVect());
   CellToEdge(dfn, dfn_face);

   // Compute face-centered velocity components in configuration directions
   
   LevelData<FluxBox> vel_cfg(grids, CFG_DIM, IntVect::Unit);
   computeVelCfgComp( vel_cfg, phase_geom, m_field, PhaseGeom::DIAGNOSTICS);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      FluxBox& this_result = a_result[dit];
      FluxBox& this_dfn = dfn_face[dit];
      FluxBox& this_vel = vel_cfg[dit];

      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result_dir = this_result[dir];
         FArrayBox& this_dfn_dir = this_dfn[dir];

         this_result_dir.copy(this_vel[dir], 0, 0, CFG_DIM);
         for (int dir2=0; dir2<CFG_DIM; ++dir2) {
            this_result_dir.mult(this_dfn_dir, 0, dir2, 1);
         }
      }
   }
}


template <> void
FourthMomentKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                     const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& B_injected = geometry.getBFieldMagnitude();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     //const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);

      FArrayBox& this_result = a_result[dit];
      const FArrayBox& this_B = B_injected[dit];

      const Box& Bbox = this_B.box();
      int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
      int mu_index = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

      BoxIterator bit(this_result.box());
      for (bit.begin(); bit.ok(); ++bit) {
         IntVect iv = bit();
         IntVect ivB = iv;
         ivB[VPARALLEL_DIR] = vp_index;
         ivB[MU_DIR] = mu_index;
         double v_parallel = velocityRealCoords(iv,0);
         double mu = velocityRealCoords(iv,1);
         double v_perp2 =  mu * this_B(ivB) ;
         double v2 = (mass * v_parallel * v_parallel + v_perp2);
         double v4 = v2 * v2;

         // Multiply the a_result by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) *= v4/4.0;
         }
      }
   }
}


template <> void
FourthMomentKernel<FluxBox>::eval( LevelData<FluxBox>&    a_result,
                                   const KineticSpecies&  a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FluxBox>& B_injected_face = geometry.getBFieldMagnitudeFace();
   
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) { 
      
         FArrayBox& this_result = a_result[dit][dir];
         const FArrayBox& this_B = B_injected_face[dit][dir];

         const Box& Bbox = this_B.box();
         int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
         int mu_index = Bbox.smallEnd(MU_DIR);

         // Get the physical velocity coordinates for this part of phase space
         const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

         BoxIterator bit(velocityRealCoords.box());
         for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            IntVect ivB = iv;
            ivB[VPARALLEL_DIR] = vp_index;
            ivB[MU_DIR] = mu_index;
            double v_parallel = velocityRealCoords(iv,0);
            double mu = velocityRealCoords(iv,1);
            double v_perp2 =  mu * this_B(ivB) ;
            double v2 = (mass * v_parallel * v_parallel + v_perp2);
            double v4 = v2 * v2;
            
            // Multiply the a_result by the velocity square.
            for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
               this_result(iv,n_comp) *= v4/4.0;
            }
         }
      }
   }
}


template <> void
PressureKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                 const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& B_injected = geometry.getBFieldMagnitude();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();

   LevelData<FArrayBox> vpar_inj;
   geometry.injectConfigurationToPhase( m_v_parallel_shift, vpar_inj );

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      //PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);

      FArrayBox& this_result = a_result[dit];
      const FArrayBox& this_B = B_injected[dit];
      const FArrayBox& this_vparshift = vpar_inj[dit]; 

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];


      FORT_COMPUTE_PRESSURE(CHF_FRA(this_result),
                            CHF_CONST_FRA(velocityRealCoords),
                            CHF_CONST_FRA1(this_vparshift,0),
                            CHF_CONST_FRA1(this_B,0),
                            CHF_CONST_REAL(mass),
                            CHF_BOX(this_result.box()));
   }
}


template <> void
PressureKernel<FluxBox>::eval( LevelData<FluxBox>&    a_result,
                               const KineticSpecies&  a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FluxBox>& B_injected_face = geometry.getBFieldMagnitudeFace();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();

   LevelData<FArrayBox> vpar_inj;
   geometry.injectConfigurationToPhase( m_v_parallel_shift, vpar_inj );

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];
         const FArrayBox& this_B = B_injected_face[dit][dir];
         const FArrayBox& this_vparshift = vpar_inj[dit]; 

         // Get the physical velocity coordinates for this part of phase space
         const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

         FORT_COMPUTE_PRESSURE(CHF_FRA(this_result),
                               CHF_CONST_FRA(velocityRealCoords),
                               CHF_CONST_FRA1(this_vparshift,0),
                               CHF_CONST_FRA1(this_B,0),
                               CHF_CONST_REAL(mass),
                               CHF_BOX(this_result.box()));
      }
   }
}


template <> void
ParallelPressureKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                         const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& B_injected = geometry.getBFieldMagnitude();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();

   LevelData<FArrayBox> vpar_inj;
   geometry.injectConfigurationToPhase( m_v_parallel_shift, vpar_inj );

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_result = a_result[dit];
      const FArrayBox& this_B = B_injected[dit];
      const FArrayBox& this_vparshift = vpar_inj[dit];

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

      FORT_COMPUTE_PARALLEL_PRESSURE(CHF_FRA(this_result),
                                     CHF_CONST_FRA(velocityRealCoords),
                                     CHF_CONST_FRA1(this_vparshift,0),
                                     CHF_CONST_FRA1(this_B,0),
                                     CHF_CONST_REAL(mass),
                                     CHF_BOX(this_result.box()));
   }
}


template <> void
ParallelPressureKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                       const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FluxBox>& B_injected_face = geometry.getBFieldMagnitudeFace();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();

   LevelData<FArrayBox> vpar_inj;
   geometry.injectConfigurationToPhase( m_v_parallel_shift, vpar_inj );

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];
         const FArrayBox& this_B = B_injected_face[dit][dir];
         const FArrayBox& this_vparshift = vpar_inj[dit];

         // Get the physical velocity coordinates for this part of phase space
         const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

         FORT_COMPUTE_PARALLEL_PRESSURE(CHF_FRA(this_result),
                                        CHF_CONST_FRA(velocityRealCoords),
                                        CHF_CONST_FRA1(this_vparshift,0),
                                        CHF_CONST_FRA1(this_B,0),
                                        CHF_CONST_REAL(mass),
                                        CHF_BOX(this_result.box()));
      }
   }
}


template <> void
PerpPressureKernel<FArrayBox>::eval(LevelData<FArrayBox>& a_result,
                                    const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& B_injected = geometry.getBFieldMagnitude();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_result = a_result[dit];
      const FArrayBox& this_B = B_injected[dit];

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

      FORT_COMPUTE_PERP_PRESSURE(CHF_FRA(this_result),
                                 CHF_CONST_FRA(velocityRealCoords),
                                 CHF_CONST_FRA1(this_B,0),
                                 CHF_CONST_REAL(mass),
                                 CHF_BOX(this_result.box()));
   }
}


template <> void
PerpPressureKernel<FluxBox>::eval(LevelData<FluxBox>& a_result,
                                  const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FluxBox>& B_injected_face = geometry.getBFieldMagnitudeFace();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];
         const FArrayBox& this_B = B_injected_face[dit][dir];

         // Get the physical velocity coordinates for this part of phase space
         const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

         FORT_COMPUTE_PERP_PRESSURE(CHF_FRA(this_result),
                                    CHF_CONST_FRA(velocityRealCoords),
                                    CHF_CONST_FRA1(this_B,0),
                                    CHF_CONST_REAL(mass),
                                    CHF_BOX(this_result.box()));
      }
   }
}


template <> void
ParallelHeatFluxKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                         const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   
   const LevelData<FArrayBox>& B_injected = geometry.getBFieldMagnitude();
   
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();
   
   LevelData<FArrayBox> vpar_inj;
   geometry.injectConfigurationToPhase( m_v_parallel_shift, vpar_inj );
   
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      //PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);
      
      FArrayBox& this_result = a_result[dit];
      const FArrayBox& this_B = B_injected[dit];
      const FArrayBox& this_vparshift = vpar_inj[dit];
      
      const Box& Bbox = this_B.box();
      int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
      int mu_index = Bbox.smallEnd(MU_DIR);
      
      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];
      
      BoxIterator bit(this_result.box());
      for (bit.begin(); bit.ok(); ++bit) {
         IntVect iv = bit();
         IntVect ivB = iv;
         ivB[VPARALLEL_DIR] = vp_index;
         ivB[MU_DIR] = mu_index;
         double v_parallel = velocityRealCoords(iv,0);
         double mu = velocityRealCoords(iv,1);
         double v_perp2 =  mu * this_B(ivB) ;
         double v2 = (mass * pow(v_parallel-this_vparshift(ivB,0),2) + v_perp2);
         
         // Multiply the a_result by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) *= (v_parallel-this_vparshift(ivB,0)) * v2 / 2.0;
         }
      }
   }
}


template <> void
ParallelHeatFluxKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                       const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FluxBox>& B_injected_face = geometry.getBFieldMagnitudeFace();
   
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();
   
   LevelData<FArrayBox> vpar_inj;
   geometry.injectConfigurationToPhase( m_v_parallel_shift, vpar_inj );
   
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];
         const FArrayBox& this_B = B_injected_face[dit][dir];
         const FArrayBox& this_vparshift = vpar_inj[dit];
      
         const Box& Bbox = this_B.box();
         int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
         int mu_index = Bbox.smallEnd(MU_DIR);
      
         // Get the physical velocity coordinates for this part of phase space
         const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];
      
         BoxIterator bit(velocityRealCoords.box());
         for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            IntVect ivB = iv;
            ivB[VPARALLEL_DIR] = vp_index;
            ivB[MU_DIR] = mu_index;
            double v_parallel = velocityRealCoords(iv,0);
            double mu = velocityRealCoords(iv,1);
            double v_perp2 =  mu * this_B(ivB) ;
            double v2 = (mass * pow(v_parallel-this_vparshift(ivB,0),2) + v_perp2);
         
            // Multiply the a_result by the velocity square.
            for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
               this_result(iv,n_comp) *= (v_parallel-this_vparshift(ivB,0)) * v2 / 2.0;
            }
         }
      }
   }
}


template <> void
EnergyKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                               const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom&             geometry    = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>&  B_injected  = geometry.getBFieldMagnitude();
   const LevelData<FArrayBox>&  dfn         = a_kinetic_species.distributionFunction();
   const Real mass = a_kinetic_species.mass();
   const Real charge = a_kinetic_species.charge();

   LevelData<FArrayBox> phi_injected;
   geometry.injectConfigurationToPhase(m_phi, phi_injected);
  
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

      //PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);
      FArrayBox&                this_result     = a_result[dit];
      const FArrayBox&          this_B          = B_injected[dit];
      const FArrayBox&          this_phi        = phi_injected[dit];

      const Box& Bbox = this_B.box();
      int   vp_index  = Bbox.smallEnd(VPARALLEL_DIR);
      int   mu_index  = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

      BoxIterator bit(this_result.box());
      for (bit.begin(); bit.ok(); ++bit) {

         IntVect iv   = bit();
         IntVect ivB  = iv;

         ivB[VPARALLEL_DIR] = vp_index;
         ivB[MU_DIR]        = mu_index;

         double v_parallel  = velocityRealCoords(iv,0);
         double mu          = velocityRealCoords(iv,1);
         double mv_perp2     = mu * this_B(ivB) ;
         double mv2          = mass*v_parallel*v_parallel + mv_perp2;
         double energy       = mv2/2.0 + charge*this_phi(ivB,0);

         // Multiply the a_result by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) *= energy;
         }
      }
   }
}


template <> void
EnergyKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                             const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FluxBox>& B_injected_face = geometry.getBFieldMagnitudeFace();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const Real mass = a_kinetic_species.mass();
   const Real charge = a_kinetic_species.charge();
  
   const CFG::DisjointBoxLayout& grids_cfg = m_phi.getBoxes();
   CFG::LevelData<CFG::FArrayBox> phi_tmp(grids_cfg, 1, 2*CFG::IntVect::Unit);

   CFG::DataIterator dit_cfg = grids_cfg.dataIterator();
   for (dit_cfg.begin(); dit_cfg.ok(); ++dit_cfg) {
     phi_tmp[dit_cfg].copy(m_phi[dit_cfg]);
   }
  
   const CFG::MagGeom& mag_geom( geometry.magGeom() );
   mag_geom.extrapolateToPhysicalGhosts(phi_tmp, !mag_geom.secondOrder());
  
   CFG::IntVect ghost_cfg = (geometry.secondOrder()) ? CFG::IntVect::Zero : CFG::IntVect::Unit;
   CFG::LevelData<CFG::FluxBox> phi_face(grids_cfg, 1, ghost_cfg);
   mag_geom.convertCellToFace(phi_face, phi_tmp);
  
   LevelData<FluxBox> phi_face_injected;
   geometry.injectConfigurationToPhase(phi_face, phi_tmp, phi_face_injected);
  
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   LevelData<FluxBox> realCoordsOnFaces(grids, SpaceDim, IntVect::Zero);
   geometry.getFaceCenteredRealCoords(realCoordsOnFaces);
   
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];
         const FArrayBox& this_B = B_injected_face[dit][dir];
         const FArrayBox& this_phi = phi_face_injected[dit][dir];
         const FArrayBox& this_real_coords = realCoordsOnFaces[dit][dir];

         const Box& Bbox = this_B.box();
         int   vp_index  = Bbox.smallEnd(VPARALLEL_DIR);
         int   mu_index  = Bbox.smallEnd(MU_DIR);


         BoxIterator bit(this_result.box());
         for (bit.begin(); bit.ok(); ++bit) {

            IntVect iv   = bit();
            IntVect ivB  = iv;

            ivB[VPARALLEL_DIR] = vp_index;
            ivB[MU_DIR]        = mu_index;

            double v_parallel  = this_real_coords(iv,VPARALLEL_DIR);
            double mu          = this_real_coords(iv,MU_DIR);
            double mv_perp2     = mu * this_B(ivB) ;
            double mv2          = mass*v_parallel*v_parallel + mv_perp2;
            double energy       = mv2/2.0 + charge*this_phi(ivB,0);

            // Multiply the a_result by the velocity square.
            for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
               this_result(iv,n_comp) *= energy;
            }
         }
      }
   }
}

template <> void
KineticEnergyKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                     const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom&             geometry    = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>&  B_injected  = geometry.getBFieldMagnitude();
   const LevelData<FArrayBox>&  dfn         = a_kinetic_species.distributionFunction();
   const Real mass = a_kinetic_species.mass();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

      //PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);
      FArrayBox&                this_result     = a_result[dit];
      const FArrayBox&          this_B          = B_injected[dit];

      const Box& Bbox = this_B.box();
      int   vp_index  = Bbox.smallEnd(VPARALLEL_DIR);
      int   mu_index  = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

      BoxIterator bit(this_result.box());
      for (bit.begin(); bit.ok(); ++bit) {

         IntVect iv   = bit();
         IntVect ivB  = iv;

         ivB[VPARALLEL_DIR] = vp_index;
         ivB[MU_DIR]        = mu_index;

         double v_parallel  = velocityRealCoords(iv,0);
         double mu          = velocityRealCoords(iv,1);
         double mv_perp2     = mu * this_B(ivB) ;
         double mv2          = mass*v_parallel*v_parallel + mv_perp2;

         // Multiply the a_result by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) *= mv2 / 2.0;
         }
      }
   }
}


template <> void
KineticEnergyKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                   const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FluxBox>& B_injected_face = geometry.getBFieldMagnitudeFace();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const Real mass = a_kinetic_species.mass();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];
         const FArrayBox& this_B = B_injected_face[dit][dir];

         const Box& Bbox = this_B.box();
         int   vp_index  = Bbox.smallEnd(VPARALLEL_DIR);
         int   mu_index  = Bbox.smallEnd(MU_DIR);

         // Get the physical velocity coordinates for this part of phase space
         const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

         BoxIterator bit(this_result.box());
         for (bit.begin(); bit.ok(); ++bit) {

            IntVect iv   = bit();
            IntVect ivB  = iv;

            ivB[VPARALLEL_DIR] = vp_index;
            ivB[MU_DIR]        = mu_index;

            double v_parallel  = velocityRealCoords(iv,0);
            double mu          = velocityRealCoords(iv,1);
            double mv_perp2     = mu * this_B(ivB) ;
            double mv2          = mass*v_parallel*v_parallel + mv_perp2;

            // Multiply the a_result by the velocity square.
            for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
               this_result(iv,n_comp) *= mv2 / 2.0;
            }
         }
      }
   }
}

template <> void
PerpEnergyKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                   const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& B_injected = geometry.getBFieldMagnitude();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      //PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);

      FArrayBox& this_result = a_result[dit];
      const FArrayBox& this_B = B_injected[dit];

      const Box& Bbox = this_B.box();
      int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
      int mu_index = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

      BoxIterator bit(this_result.box());
      for (bit.begin(); bit.ok(); ++bit) {
         IntVect iv = bit();
         IntVect ivB = iv;
         ivB[VPARALLEL_DIR] = vp_index;
         ivB[MU_DIR] = mu_index;
         double mu = velocityRealCoords(iv,1);
         double mv_perp2 =  mu * this_B(ivB) ;

         // Multiply the a_result by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) *= mv_perp2/2.0;
         }
      }
   }
}


template <> void
PerpEnergyKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                 const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FluxBox>& B_injected_face = geometry.getBFieldMagnitudeFace();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];
         const FArrayBox& this_B = B_injected_face[dit][dir];

         const Box& Bbox = this_B.box();
         int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
         int mu_index = Bbox.smallEnd(MU_DIR);

         // Get the physical velocity coordinates for this part of phase space
         const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

         BoxIterator bit(this_result.box());
         for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            IntVect ivB = iv;
            ivB[VPARALLEL_DIR] = vp_index;
            ivB[MU_DIR] = mu_index;
            double mu = velocityRealCoords(iv,1);
            double mv_perp2 =  mu * this_B(ivB) ;

            // Multiply the a_result by the velocity square.
            for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
               this_result(iv,n_comp) *= mv_perp2/2.0;
            }
         }
      }
   }
}


template <> void
ParallelEnergyKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                       const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom&             geometry    = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>&  B_injected  = geometry.getBFieldMagnitude();
   const LevelData<FArrayBox>&  dfn         = a_kinetic_species.distributionFunction();
   const Real mass = a_kinetic_species.mass();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

      FArrayBox&                this_result     = a_result[dit];
      const FArrayBox&          this_B          = B_injected[dit];

      const Box& Bbox = this_B.box();
      int   vp_index  = Bbox.smallEnd(VPARALLEL_DIR);
      int   mu_index  = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];
     
      BoxIterator bit(this_result.box());
      for (bit.begin(); bit.ok(); ++bit) {

         IntVect iv   = bit();
         IntVect ivB  = iv;

         ivB[VPARALLEL_DIR] = vp_index;
         ivB[MU_DIR]        = mu_index;

         double v_parallel  = velocityRealCoords(iv,0);
         double mv2          = mass*v_parallel*v_parallel;

         // Multiply the a_result by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) *= mv2 / 2.0;
         }
      }
   }
}


template <> void
ParallelEnergyKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                     const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom&             geometry    = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FluxBox>&  B_injected_face  = geometry.getBFieldMagnitudeFace();

   const LevelData<FArrayBox>&  dfn         = a_kinetic_species.distributionFunction();
   const Real mass = a_kinetic_species.mass();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];
         const FArrayBox& this_B = B_injected_face[dit][dir];

         const Box& Bbox = this_B.box();
         int   vp_index  = Bbox.smallEnd(VPARALLEL_DIR);
         int   mu_index  = Bbox.smallEnd(MU_DIR);

         // Get the physical velocity coordinates for this part of phase space
         const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

         BoxIterator bit(this_result.box());
         for (bit.begin(); bit.ok(); ++bit) {

            IntVect iv   = bit();
            IntVect ivB  = iv;
            
            ivB[VPARALLEL_DIR] = vp_index;
            ivB[MU_DIR]        = mu_index;

            double v_parallel  = velocityRealCoords(iv,0);
            double mv2          = mass*v_parallel*v_parallel;

            // Multiply the a_result by the velocity square.
            for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
               this_result(iv,n_comp) *= mv2 / 2.0;
            }
         }
      }
   }
}


template <> void
ParallelMomKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                    const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
 
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_result = a_result[dit];

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

      FORT_COMPUTE_PAR_MOM(CHF_FRA(this_result),
			   CHF_CONST_FRA1(velocityRealCoords,0),
			   CHF_BOX(this_result.box()));
   }
}


template <> void
ParallelMomKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                  const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
 
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];

         // Get the physical velocity coordinates for this part of phase space
         const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];

         FORT_COMPUTE_PAR_MOM(CHF_FRA(this_result),
                              CHF_CONST_FRA1(velocityRealCoords,0),
                              CHF_BOX(this_result.box()));
      }
   }
}


template <> void
ParticleFluxKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                     const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const DisjointBoxLayout & grids = dfn.getBoxes();

   //Compute CFG componenents of cell-centered normal velocities
   LevelData<FArrayBox> cellVelCFG(grids, CFG_DIM, IntVect::Zero);
   computeVelCfgCompNormals( cellVelCFG, phase_geom, m_field, m_velocity_option );

   //Compute radial particle flux
   DataIterator dit = dfn.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_result = a_result[dit];
      for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
        this_result.mult(cellVelCFG[dit], RADIAL_DIR, n_comp, 1);
      }
   }
}


template <> void
ParticleFluxKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                   const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const DisjointBoxLayout & grids = dfn.getBoxes();

   //Compute CFG componenents of face-centered normal velocities
   LevelData<FluxBox> vel_cfg(grids, CFG_DIM, IntVect::Zero);
   computeVelCfgCompNormals( vel_cfg, phase_geom, m_field, m_velocity_option );

   //Compute radial particle flux
   DataIterator dit = dfn.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
           this_result.mult(vel_cfg[dit][dir], RADIAL_DIR, n_comp, 1);
         }
      }
   }
}


template <> void
HeatFluxKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                 const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& B_injected = phase_geom.getBFieldMagnitude();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();
   double charge = a_kinetic_species.charge();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   
   //Compute CFG componenents of cell-centered normal velocities
   LevelData<FArrayBox> cellVelCFG(grids, CFG_DIM, IntVect::Zero);
   computeVelCfgCompNormals( cellVelCFG, phase_geom, m_field, m_velocity_option );
   
   //Compute radial heat flux
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_result = a_result[dit];

      const FArrayBox& this_B = B_injected[dit];
      const FArrayBox& this_m_phi = m_phi[dit]; 

      const Box& Bbox = this_B.box();
      int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
      int mu_index = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = phase_geom.getVelocityRealCoords()[dit];

      BoxIterator bit(this_result.box());
      for (bit.begin(); bit.ok(); ++bit) {
         IntVect iv = bit();
         IntVect ivB = iv;
         ivB[VPARALLEL_DIR] = vp_index;
         ivB[MU_DIR] = mu_index;
         double v_parallel = velocityRealCoords(iv,0);
         double mu = velocityRealCoords(iv,1);
         double v_perp2 =  mu * this_B(ivB) ;
         double v2 = (mass * v_parallel * v_parallel + v_perp2);

         // Multiply the distribution function by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) *= (v2 / 2.0 + charge*this_m_phi(ivB,0));
         }
      }
      for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
         this_result.mult(cellVelCFG[dit], RADIAL_DIR, n_comp, 1);
      }
   }
}


template <> void
HeatFluxKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                               const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FluxBox>& B_injected_face = phase_geom.getBFieldMagnitudeFace();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();
   double charge = a_kinetic_species.charge();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   
   //Compute CFG componenents of face-centered normal velocities
   LevelData<FluxBox> vel_cfg(grids, CFG_DIM, IntVect::Zero);
   computeVelCfgCompNormals( vel_cfg, phase_geom, m_field, m_velocity_option );
   
   //Compute radial heat flux
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];

         const FArrayBox& this_B = B_injected_face[dit][dir];
         const FArrayBox& this_m_phi = m_phi[dit]; 

         const Box& Bbox = this_B.box();
         int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
         int mu_index = Bbox.smallEnd(MU_DIR);

         // Get the physical velocity coordinates for this part of phase space
         const FArrayBox& velocityRealCoords = phase_geom.getVelocityRealCoords()[dit];

         BoxIterator bit(velocityRealCoords.box());
         for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            IntVect ivB = iv;
            ivB[VPARALLEL_DIR] = vp_index;
            ivB[MU_DIR] = mu_index;
            double v_parallel = velocityRealCoords(iv,0);
            double mu = velocityRealCoords(iv,1);
            double v_perp2 =  mu * this_B(ivB) ;
            double v2 = (mass * v_parallel * v_parallel + v_perp2);

            // Multiply the distribution function by the velocity square.
            for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
               this_result(iv,n_comp) *= (v2 / 2.0 + charge*this_m_phi(ivB,0));
            }
         }

         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result.mult(vel_cfg[dit][dir], RADIAL_DIR, n_comp, 1);
         }
      }
   }
}


template <> void
GuidingCenterPoloidalMomKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                                 const KineticSpecies& a_kinetic_species ) const
{

   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const DisjointBoxLayout & grids = dfn.getBoxes();
   
   //Compute cell-centered velocity configuration components
   LevelData<FArrayBox> cellVelCFG(grids, CFG_DIM, IntVect::Zero);
   computeVelCfgComp( cellVelCFG, phase_geom, m_field, PhaseGeom::DIAGNOSTICS );
   
   //Compute poloidal velocity
   LevelData<FArrayBox> cellVel_theta(grids, 1, IntVect::Zero);
   phase_geom.computePoloidalProjection( cellVel_theta, cellVelCFG );
   
   //Compute radial particle flux
   DataIterator dit = dfn.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_result = a_result[dit];
      FArrayBox& this_vtheta = cellVel_theta[dit];
      
      this_result.mult(this_vtheta, 0, 0, 1);
   }

}


template <> void
GuidingCenterPoloidalMomKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                               const KineticSpecies& a_kinetic_species ) const
{

   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const DisjointBoxLayout & grids = dfn.getBoxes();
   
   //Compute cell-centered velocity configuration components
   LevelData<FluxBox> cellVelCFG(grids, CFG_DIM, IntVect::Zero);
   computeVelCfgComp( cellVelCFG, phase_geom, m_field, PhaseGeom::DIAGNOSTICS );
   
   //Compute poloidal velocity
   LevelData<FluxBox> faceVel_theta(grids, 1, IntVect::Zero);
   phase_geom.computePoloidalProjection( faceVel_theta, cellVelCFG );
   
   //Compute radial particle flux
   DataIterator dit = dfn.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];
         FArrayBox& this_vtheta = faceVel_theta[dit][dir];
      
         this_result.mult(this_vtheta, 0, 0, 1);
      }
   }
}


template <> void
MagnetizationKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                      const KineticSpecies& a_kinetic_species ) const
{
   
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_result = a_result[dit];
      
      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];
      
      BoxIterator bit(this_result.box());
      for (bit.begin(); bit.ok(); ++bit) {
         IntVect iv = bit();
         double mu = velocityRealCoords(iv,1);
         
         // Multiply the a_result by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) *= -mu;
         }
      }
   }
   
}


template <> void
MagnetizationKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                    const KineticSpecies& a_kinetic_species ) const
{
   
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<CFG_DIM; ++dir) {
         FArrayBox& this_result = a_result[dit][dir];
      
         // Get the physical velocity coordinates for this part of phase space
         const FArrayBox& velocityRealCoords = geometry.getVelocityRealCoords()[dit];
      
         BoxIterator bit(velocityRealCoords.box());
         for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            double mu = velocityRealCoords(iv,1);
         
            // Multiply the a_result by the velocity square.
            for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
               this_result(iv,n_comp) *= -mu;
            }
         }
      }
   }
   
}


template <> void
DeltaFKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                               const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();

   LevelData<FArrayBox> dfn_tmp(grids, dfn.nComp(), dfn.ghostVect()+IntVect::Unit);

   for (dit.begin(); dit.ok(); ++dit) {
      dfn_tmp[dit].copy( dfn[dit] );
   }
   dfn_tmp.exchange();

   //Remove BstarParallel factor from dfn_tmp 
   //NB: this division is not 4th order accurate
   phase_geom.divideBStarParallel( dfn_tmp );

   //Convert dfn_tmp to cell-centeres
   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& coord_sys = phase_geom.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = coord_sys.domain();

      fourthOrderAverageCell(dfn_tmp[dit], domain, grids[dit], -1);
   }
   dfn_tmp.exchange();

   //Compute deltaF (i.e., dfn - Maxwellian distribution with the 3 first moments)
   const LevelData<FArrayBox>& B_injected = phase_geom.getBFieldMagnitude();
   double mass = a_kinetic_species.mass();
   double pi = Constants::PI;
   double factor = sqrt(pi * pow((2.0/mass),3));


   LevelData<FArrayBox> n_inj;
   phase_geom.injectConfigurationToPhase( m_density, n_inj );
   LevelData<FArrayBox> T_inj;
   phase_geom.injectConfigurationToPhase( m_temperature, T_inj );
   LevelData<FArrayBox> vpar_inj;
   phase_geom.injectConfigurationToPhase( m_v_parallel_shift, vpar_inj );

   for (dit.begin(); dit.ok(); ++dit) {

      FArrayBox& this_result = a_result[dit];
      FArrayBox& this_dfn_tmp = dfn_tmp[dit];
      const FArrayBox& this_B = B_injected[dit];

      const FArrayBox& this_n = n_inj[dit];
      const FArrayBox& this_T = T_inj[dit]; 
      const FArrayBox& this_vparshift = vpar_inj[dit]; 
 
      const Box& Bbox = this_B.box();
      int vpB = Bbox.smallEnd(VPARALLEL_DIR);
      int muB = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = phase_geom.getVelocityRealCoords()[dit];

      BoxIterator bit(this_result.box());
      for (bit.begin(); bit.ok(); ++bit) {
         IntVect iv = bit();
         IntVect ivB = iv;
         ivB[VPARALLEL_DIR] = vpB;
         ivB[MU_DIR] = muB;
         double v_parallel = velocityRealCoords(iv,0);
         double mu = velocityRealCoords(iv,1);
         double eparnorm = 0.5 * mass * pow(v_parallel-this_vparshift(ivB,0),2) / this_T(ivB,0);
         double munorm   = 0.5 * this_B(ivB) * mu / this_T(ivB,0);
         double val      = exp( -( eparnorm + munorm ) );
         val    = val * this_n(ivB,0) / ( factor * sqrt(this_T(ivB,0)) * this_T(ivB,0) );


         // Multiply the distribution function by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) = this_dfn_tmp(iv,n_comp) - val ;
         }
      }
   }

   //Convert a_result to cell-averages 
   for (dit.begin(); dit.ok(); ++dit) {
      dfn_tmp[dit].copy( a_result[dit] );
   }
   dfn_tmp.exchange();

   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& coord_sys = phase_geom.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = coord_sys.domain();

      fourthOrderAverageCell(dfn_tmp[dit], domain, grids[dit]);
   }
   dfn_tmp.exchange();

   for (dit.begin(); dit.ok(); ++dit) {
      a_result[dit].copy(dfn_tmp[dit]);
   }
   
   //Restore BstarParallel factor
   phase_geom.multBStarParallel( a_result );
}


template <> void
DeltaFKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                             const KineticSpecies& a_kinetic_species ) const
{
   // Implementation of this function requires BStarParallel on faces
   MayDay::Error("DeltaFKernel<FluxBox>::eval() has not been implemented");
}


template <> void
MaxwellianKernel<FArrayBox>::eval( LevelData<FArrayBox>& a_result,
                                   const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();


   //Compute Maxwellian distribution with the 3 first moments
   const LevelData<FArrayBox>& B_injected = phase_geom.getBFieldMagnitude();
   double mass = a_kinetic_species.mass();
   //   double pi = Constants::PI;
   //   double factor = sqrt(pi * pow((2.0/mass),3));


   LevelData<FArrayBox> n_inj;
   phase_geom.injectConfigurationToPhase( m_density, n_inj );
   LevelData<FArrayBox> T_inj;
   phase_geom.injectConfigurationToPhase( m_temperature, T_inj );
   LevelData<FArrayBox> vpar_inj;
   phase_geom.injectConfigurationToPhase( m_v_parallel_shift, vpar_inj );

   for (dit.begin(); dit.ok(); ++dit) {

      FArrayBox& this_result = a_result[dit];
      const FArrayBox& this_B = B_injected[dit];

      const FArrayBox& this_n = n_inj[dit];
      const FArrayBox& this_T = T_inj[dit]; 
      const FArrayBox& this_vparshift = vpar_inj[dit]; 
 
      //      const Box& Bbox = this_B.box();
      //      int vpB = Bbox.smallEnd(VPARALLEL_DIR);
      //      int muB = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      const FArrayBox& velocityRealCoords = phase_geom.getVelocityRealCoords()[dit];

      FORT_COMPUTE_MAXWELLIAN(CHF_FRA(this_result),
                              CHF_CONST_FRA(velocityRealCoords),
                              CHF_CONST_FRA1(this_vparshift,0),
                              CHF_CONST_FRA1(this_n,0),
                              CHF_CONST_FRA1(this_T,0),
                              CHF_CONST_FRA1(this_B,0),
                              CHF_CONST_REAL(mass),
                              CHF_BOX(this_result.box()));
   }

   if (!phase_geom.secondOrder()) {

     //Convert a_result to cell-averages 
     LevelData<FArrayBox> dfn_tmp(grids, dfn.nComp(), dfn.ghostVect()+IntVect::Unit);
     for (dit.begin(); dit.ok(); ++dit) {
       dfn_tmp[dit].copy( a_result[dit] );
     }
     phase_geom.fillInternalGhosts(dfn_tmp);

     for (dit.begin(); dit.ok(); ++dit) {
       const PhaseBlockCoordSys& coord_sys = phase_geom.getBlockCoordSys(grids[dit]);
       const ProblemDomain& domain = coord_sys.domain();

       fourthOrderAverageCell(dfn_tmp[dit], domain, grids[dit]);
     }

     phase_geom.fillInternalGhosts(dfn_tmp);

     for (dit.begin(); dit.ok(); ++dit) {
       a_result[dit].copy(dfn_tmp[dit]);
     }

   }
   
   //Multiply BstarParallel factor
   phase_geom.multBStarParallel( a_result );
}


template <> void
MaxwellianKernel<FluxBox>::eval( LevelData<FluxBox>& a_result,
                                 const KineticSpecies& a_kinetic_species ) const
{
   // Implementation of this function requires BStarParallel on faces
   MayDay::Error("MaxwellianKernel<FluxBox>::eval() has not been implemented");
}

#include "NamespaceFooter.H"
