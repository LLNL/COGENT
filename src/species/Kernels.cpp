enum SingleNullBlockType {LCORE,RCORE,LCSOL,RCSOL,LSOL,RSOL,LPF,RPF,MCORE,MCSOL,NUM_SINGLE_NULL_BLOCKS};

#include "Kernels.H"
#include "Directions.H"
#include "PhaseBlockCoordSys.H"
#include "FourthOrderUtil.H"
#include "EdgeToCell.H"
#include "ConstFact.H"

#include "KernelsF_F.H"

#include "NamespaceHeader.H"



void
DensityKernel::eval( LevelData<FArrayBox>& a_result,
                     const KineticSpecies& a_kinetic_species ) const
{
   /*
     Nothing to do.  a_result has already been initialized by the distribution function,
     and the moment kernel is unity.
   */
}



void
ChargeDensityKernel::eval( LevelData<FArrayBox>& a_result,
                           const KineticSpecies& a_kinetic_species ) const
{
   /*
     Nothing to do.  a_result has already been initialized by the distribution function,
     and the moment kernel is unity.  Multiplication by the charge (the scale factor
     for this Kernel class) is performed by the MomentOp after integrating over velocity.
   */
}



void
MassDensityKernel::eval( LevelData<FArrayBox>& a_result,
                         const KineticSpecies& a_kinetic_species ) const
{
   /*
     Nothing to do.  a_result has already been initialized by the distribution function,
     and the moment kernel is unity.  Multiplication by the mass (the scale factor
     for this Kernel class) is performed by the MomentOp after integrating over velocity.
   */
}



void
MomentumDensityKernel::eval( LevelData<FArrayBox>& a_result,
                             const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const int n_kernel_comp = nComponents();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);

      FArrayBox& this_result = a_result[dit];

      // Get the physical velocity coordinates for this part of phase space
      FArrayBox velocityRealCoords(this_result.box(), n_kernel_comp);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);

      // Multiply the a_result by the velocity.  Multiplication by the
      // mass is performed later after integrating over velocity.
      for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
         this_result.mult(velocityRealCoords, 0, n_comp*n_kernel_comp, n_kernel_comp);
      }
   }
}

void

FourthMomentKernel::eval( LevelData<FArrayBox>& a_result,
                      const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& B_injected = geometry.getBFieldMagnitude();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);

      FArrayBox& this_result = a_result[dit];
      const FArrayBox& this_B = B_injected[dit];

      const Box& Bbox = this_B.box();
      int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
      int mu_index = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      FArrayBox velocityRealCoords(this_result.box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);

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

PressureKernel::PressureKernel(const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift)
   : m_v_parallel_shift(v_parallel_shift)
{
}


void
PressureKernel::eval( LevelData<FArrayBox>& a_result,
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
      const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);

      FArrayBox& this_result = a_result[dit];
      const FArrayBox& this_B = B_injected[dit];
      const FArrayBox& this_vparshift = vpar_inj[dit]; 

      const Box& Bbox = this_B.box();
      int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
      int mu_index = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      FArrayBox velocityRealCoords(this_result.box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);

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
            this_result(iv,n_comp) *= v2 / 3.0;
         }
      }
   }
}


void
PerpEnergyKernel::eval( LevelData<FArrayBox>& a_result,
                      const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& B_injected = geometry.getBFieldMagnitude();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);

      FArrayBox& this_result = a_result[dit];
      const FArrayBox& this_B = B_injected[dit];

      const Box& Bbox = this_B.box();
      int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
      int mu_index = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      FArrayBox velocityRealCoords(this_result.box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);

      BoxIterator bit(velocityRealCoords.box());
      for (bit.begin(); bit.ok(); ++bit) {
         IntVect iv = bit();
         IntVect ivB = iv;
         ivB[VPARALLEL_DIR] = vp_index;
         ivB[MU_DIR] = mu_index;
         double mu = velocityRealCoords(iv,1);
         double v_perp2 =  mu * this_B(ivB) ;

         // Multiply the a_result by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) *= v_perp2/2.0;
         }
      }
   }
}

void
ParallelMomKernel::eval( LevelData<FArrayBox>& a_result,
                      const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& B_injected = geometry.getBFieldMagnitude();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
 
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);

      FArrayBox& this_result = a_result[dit];
      const FArrayBox& this_B = B_injected[dit];

      const Box& Bbox = this_B.box();
      int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
      int mu_index = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      FArrayBox velocityRealCoords(this_result.box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);

      BoxIterator bit(velocityRealCoords.box());
      for (bit.begin(); bit.ok(); ++bit) {
         IntVect iv = bit();
         IntVect ivB = iv;
         ivB[VPARALLEL_DIR] = vp_index;
         ivB[MU_DIR] = mu_index;
         double v_parallel = velocityRealCoords(iv,0);

         // Multiply the a_result by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) *= v_parallel;
         }
      }
   }
}


ParticleFluxKernel::ParticleFluxKernel (const LevelData<FluxBox>& field)
   : m_field(field)
{
}


void
ParticleFluxKernel::eval( LevelData<FArrayBox>& a_result,
                      const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const DisjointBoxLayout & grids = dfn.getBoxes();

   //Compute cell-centered psi and theta velocity projections
   LevelData<FArrayBox> CellVel(grids, SpaceDim, IntVect::Zero);
   computeVelCell( CellVel, phase_geom, m_field );

   //Compute radial particle flux
   DataIterator dit = dfn.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_result = a_result[dit];
      FArrayBox& this_vel = CellVel[dit];

      this_result.mult(this_vel, 0, 0, 1);
   }
}

HeatFluxKernel::HeatFluxKernel (const LevelData<FluxBox>& field,
                                const LevelData<FArrayBox>& phi )
   : m_field(field),
     m_phi(phi)
{
}


void
HeatFluxKernel::eval( LevelData<FArrayBox>& a_result,
                      const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();

   const LevelData<FArrayBox>& B_injected = phase_geom.getBFieldMagnitude();

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   double mass = a_kinetic_species.mass();
   double charge = a_kinetic_species.charge();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();

   //Compute cell-centered psi and theta velocity projections
   LevelData<FArrayBox> CellVel(grids, SpaceDim, IntVect::Zero);
   computeVelCell( CellVel, phase_geom, m_field );

   //Compute radial heat flux
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(grids[dit]);

      FArrayBox& this_result = a_result[dit];
      FArrayBox& this_vel = CellVel[dit];

      const FArrayBox& this_B = B_injected[dit];
      const FArrayBox& this_m_phi = m_phi[dit]; 

      const Box& Bbox = this_B.box();
      int vp_index = Bbox.smallEnd(VPARALLEL_DIR);
      int mu_index = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      FArrayBox velocityRealCoords(this_result.box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);

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
            this_result(iv,n_comp) *= (v2 / 2.0 + charge*this_m_phi(ivB,0)) ;
         }
      }

     this_result.mult(this_vel, 0, 0, 1);

   }
}

PoloidalMomKernel::PoloidalMomKernel (const LevelData<FluxBox>& field, const double larmor) 
   : m_field(field),
     m_larmor(larmor)
{
}


void
PoloidalMomKernel::eval( LevelData<FArrayBox>& a_result,
                      const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const CFG::MagGeom & mag_geom = phase_geom.magGeom();
   const MultiBlockCoordSys& multi_block_coord_sys( *(phase_geom.coordSysPtr()) );

   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   double charge = a_kinetic_species.charge();

   //Create dfn_tmp with two layers of ghost cells filled by extrapolation of the interior data
   LevelData<FArrayBox> dfn_tmp(grids, dfn.nComp(), dfn.ghostVect()+2*IntVect::Unit);
   for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      dfn_tmp[dit].copy(dfn[dit]);
   }
   dfn_tmp.exchange();

   for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(grids[dit]);
      const ProblemDomain& block_domain = block_coord_sys.domain();

      fourthOrderCellExtrapAtDomainBdry(dfn_tmp[dit], block_domain, grids[dit]);
   }
   dfn_tmp.exchange();

   //Get magnetic field parameters in the configuration space
   const CFG::LevelData<CFG::FArrayBox>& cfg_BFieldCell = mag_geom.getCCBField();
   const CFG::LevelData<CFG::FArrayBox>& cfg_BFieldMagCell = mag_geom.getCCBFieldMag();
   const CFG::LevelData<CFG::FArrayBox>& cfg_curlBFieldDirCell = mag_geom.getCCCurlBFieldDir();

   //Create the 2-component (R, Z) transerve version of curlb 
   CFG::LevelData<CFG::FArrayBox> cfg_curlbTwoComp(mag_geom.grids(), 2, CFG::IntVect::Zero);
   CFG::DataIterator dit = cfg_curlbTwoComp.dataIterator();
   for ( dit.begin(); dit.ok(); ++dit ) {
      CFG::FArrayBox& this_curlbTwoComp = cfg_curlbTwoComp[dit];
      const CFG::FArrayBox& this_curlbThreeComp = cfg_curlBFieldDirCell[dit];
      CFG::BoxIterator bit(this_curlbTwoComp.box());
      for (bit.begin(); bit.ok(); ++bit) {
         CFG::IntVect iv = bit();
         this_curlbTwoComp(iv,0) = this_curlbThreeComp(iv,0);
         this_curlbTwoComp(iv,1) = this_curlbThreeComp(iv,2);
      }   
   }

   //Inject the magnetic field parameters into the phase space
   LevelData<FArrayBox> curlbTwoComp;
   phase_geom.injectConfigurationToPhase( cfg_curlbTwoComp, curlbTwoComp);

   LevelData<FArrayBox> BFieldCell;
   phase_geom.injectConfigurationToPhase( cfg_BFieldCell, BFieldCell);

   LevelData<FArrayBox> BFieldMagCell;
   phase_geom.injectConfigurationToPhase( cfg_BFieldMagCell, BFieldMagCell);

   //Compute cell-centered psi and theta velocity projections
   LevelData<FArrayBox> CellVel(grids, SpaceDim, IntVect::Zero);
   computeVelCell( CellVel, phase_geom, m_field );

   //Compute poloidal flow velocity 
   const VEL::VelCoordSys& vel_coords = phase_geom.velSpaceCoordSys();
   const VEL::RealVect& vel_dx = vel_coords.dx();

   LevelData<FArrayBox> gradFFactor(grids, 1, IntVect::Zero);
   
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {

      const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(grids[dit]);
      const RealVect& coord_dx =  block_coord_sys.dx();

      FArrayBox& this_result = a_result[dit];
      FArrayBox& this_vel = CellVel[dit];
      const FArrayBox& this_dfn = dfn_tmp[dit];
      const FArrayBox& this_Bmagn = BFieldMagCell[dit];
      const FArrayBox& this_Bvect = BFieldCell[dit];
      
      FArrayBox& this_curlb = curlbTwoComp[dit];
      block_coord_sys.computePsiThetaProjections(this_curlb);

      FArrayBox& this_GradFFactor = gradFFactor[dit];
      block_coord_sys.computeGradFFactor(this_GradFFactor);

      const int block_type( multi_block_coord_sys.whichBlock( grids[dit] ) );
      
      //Multiply GradFFactor by -1.0 for the PF region, because the mapped coordinate
      //increases toward separatrix in the PF blocks (This multiplication is 
      // required to to make the calculation of gradF consistent between blocks)
      if (block_type == RPF || block_type == LPF) {this_GradFFactor.mult(-1.0);}

      FORT_GET_POL_MOMENTUM(CHF_BOX(this_result.box()),
                            CHF_CONST_REALVECT(vel_dx),
                            CHF_CONST_REALVECT(coord_dx),
                            CHF_FRA1(this_result,0),
                            CHF_CONST_FRA1(this_dfn,0),
                            CHF_FRA(this_vel),
                            CHF_CONST_FRA1(this_Bmagn,0),
                            CHF_CONST_FRA(this_Bvect),
                            CHF_CONST_FRA(this_curlb),
                            CHF_CONST_FRA1(this_GradFFactor,0),
                            CHF_CONST_REAL(m_larmor),
                            CHF_CONST_REAL(charge));
    }
}



DeltaFKernel::DeltaFKernel (const CFG::LevelData<CFG::FArrayBox>& density,
                            const CFG::LevelData<CFG::FArrayBox>& temperature,
                            const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift)
   : m_density(density),
     m_temperature(temperature),
     m_v_parallel_shift(v_parallel_shift)

{
}


void
DeltaFKernel::eval( LevelData<FArrayBox>& a_result,
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

      const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(grids[dit]);

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
      FArrayBox velocityRealCoords(this_result.box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);

      BoxIterator bit(velocityRealCoords.box());
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



void
Kernel::computeVelCell(LevelData<FArrayBox>& a_velCell,
                       const PhaseGeom& a_phase_geom,
                       const LevelData<FluxBox>& a_field ) const
{
   const DisjointBoxLayout & grids = a_velCell.getBoxes();
   LevelData<FluxBox> FaceVel(grids, SpaceDim, IntVect::Unit);
   a_phase_geom.updateVelocities(a_field, FaceVel);

   for (DataIterator dit(a_velCell.dataIterator()); dit.ok(); ++dit) {
     const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
     FluxBox& this_FaceVel = FaceVel[dit];
     block_coord_sys.computePsiThetaProjections(this_FaceVel);
 
     FArrayBox& this_velCell = a_velCell[dit];
     for (int dir=0; dir<SpaceDim; dir++) {
          FArrayBox& this_FaceVelDir = FaceVel[dit][dir];
          FORT_COMPUTE_VEL_CELL(CHF_FRA(this_velCell),
                                CHF_FRA(this_FaceVelDir),                                 
                                CHF_CONST_INT(dir),
                                CHF_BOX(this_velCell.box()));
     }
   }
}

#include "NamespaceFooter.H"
