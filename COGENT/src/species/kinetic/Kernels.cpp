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


MomentumDensityKernel::MomentumDensityKernel(const LevelData<FluxBox>& field)
  : m_field(field)
{
}

void
MomentumDensityKernel::eval( LevelData<FArrayBox>& a_result,
                             const KineticSpecies& a_kinetic_species ) const
{
  const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
  const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
  const DisjointBoxLayout & grids = dfn.getBoxes();

  //Compute cell-centered velocity configurational components                                                                                                                                                     
  LevelData<FArrayBox> cellVelCFG(grids, CFG_DIM, IntVect::Unit);
  computeCellVelCfgComp( cellVelCFG, phase_geom, m_field);

  DataIterator dit = dfn.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    a_result[dit].copy(dfn[dit],0,0,1);
    a_result[dit].copy(dfn[dit],0,1,1);
    a_result[dit].mult(cellVelCFG[dit], 0, 0, CFG_DIM);
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

      // Get the physical velocity coordinates for this part of phase space
      FArrayBox velocityRealCoords(this_result.box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);

      FORT_COMPUTE_PRESSURE(CHF_FRA(this_result),
                            CHF_CONST_FRA(velocityRealCoords),
			    CHF_CONST_FRA1(this_vparshift,0),
			    CHF_CONST_FRA1(this_B,0),
			    CHF_CONST_REAL(mass),
                            CHF_BOX(this_result.box()));
   }
}

ParallelHeatFluxKernel::ParallelHeatFluxKernel(const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift)
: m_v_parallel_shift(v_parallel_shift)
{
}


void
ParallelHeatFluxKernel::eval( LevelData<FArrayBox>& a_result,
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
            this_result(iv,n_comp) *= (v_parallel-this_vparshift(ivB,0)) * v2 / 2.0;
         }
      }
   }
}


void
EnergyKernel::eval( LevelData<FArrayBox>& a_result,
                    const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom&             geometry    = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>&  B_injected  = geometry.getBFieldMagnitude();
   const LevelData<FArrayBox>&  dfn         = a_kinetic_species.distributionFunction();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

      const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);
      FArrayBox&                this_result     = a_result[dit];
      const FArrayBox&          this_B          = B_injected[dit];

      const Box& Bbox = this_B.box();
      int   vp_index  = Bbox.smallEnd(VPARALLEL_DIR);
      int   mu_index  = Bbox.smallEnd(MU_DIR);

      // Get the physical velocity coordinates for this part of phase space
      FArrayBox velocityRealCoords(this_result.box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);

      BoxIterator bit(velocityRealCoords.box());
      for (bit.begin(); bit.ok(); ++bit) {

         IntVect iv   = bit();
         IntVect ivB  = iv;

         ivB[VPARALLEL_DIR] = vp_index;
         ivB[MU_DIR]        = mu_index;

         double v_parallel  = velocityRealCoords(iv,0);
         double mu          = velocityRealCoords(iv,1);
         double v_perp2     = mu * this_B(ivB) ;
         double v2          = v_parallel*v_parallel + v_perp2;

         // Multiply the a_result by the velocity square.
         for (int n_comp=0; n_comp<a_result.nComp(); ++n_comp) {
            this_result(iv,n_comp) *= v2 / 2.0;
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
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
 
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);

      FArrayBox& this_result = a_result[dit];

      // Get the physical velocity coordinates for this part of phase space
      FArrayBox velocityRealCoords(this_result.box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);

      FORT_COMPUTE_PAR_MOM(CHF_FRA(this_result),
			   CHF_CONST_FRA1(velocityRealCoords,0),
			   CHF_BOX(this_result.box()));
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

   //Compute cell-centered velocity configuration components
   LevelData<FArrayBox> cellVelCFG(grids, CFG_DIM, IntVect::Zero);
   computeCellVelCfgComp( cellVelCFG, phase_geom, m_field );

   //Compute radial velocity
   LevelData<FArrayBox> cellVel_r(grids, 1, IntVect::Zero);
   phase_geom.computeRadialProjection( cellVel_r, cellVelCFG );
   
   //Compute radial particle flux
   DataIterator dit = dfn.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_result = a_result[dit];
      FArrayBox& this_vr = cellVel_r[dit];

      this_result.mult(this_vr, 0, 0, 1);
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

   //Compute cell-centered velocity configuration components
   LevelData<FArrayBox> cellVelCFG(grids, CFG_DIM, IntVect::Zero);
   computeCellVelCfgComp( cellVelCFG, phase_geom, m_field );
   
   //Compute radial velocity
   LevelData<FArrayBox> cellVel_r(grids, 1, IntVect::Zero);
   phase_geom.computeRadialProjection( cellVel_r, cellVelCFG );
   
   //Compute radial heat flux
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(grids[dit]);

      FArrayBox& this_result = a_result[dit];
      FArrayBox& this_vr = cellVel_r[dit];

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
            this_result(iv,n_comp) *= (v2 / 2.0 + charge*this_m_phi(ivB,0));
         }
      }

      this_result.mult(this_vr, 0, 0, 1);
   }
}


PerpCurrentDensityKernel::PerpCurrentDensityKernel(const LevelData<FluxBox>& field)
: m_field(field)
{
}


void
PerpCurrentDensityKernel::eval( LevelData<FArrayBox>& a_result,
                                const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const DisjointBoxLayout & grids = dfn.getBoxes();
   
   //Compute cell-centered velocity configurational components
   LevelData<FArrayBox> cellVelCFG(grids, CFG_DIM, IntVect::Unit);
   computeCellVelCfgComp( cellVelCFG, phase_geom, m_field, true );
   
   DataIterator dit = dfn.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     a_result[dit].copy(dfn[dit],0,0,1);
     a_result[dit].copy(dfn[dit],0,1,1);
     a_result[dit].mult(cellVelCFG[dit], 0, 0, CFG_DIM);
   }
}


GuidingCenterPoloidalMomKernel::GuidingCenterPoloidalMomKernel(const LevelData<FluxBox>& field)
   : m_field(field)
{
}


void
GuidingCenterPoloidalMomKernel::eval( LevelData<FArrayBox>& a_result,
                                      const KineticSpecies& a_kinetic_species ) const
{

   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   const DisjointBoxLayout & grids = dfn.getBoxes();
   
   //Compute cell-centered velocity configuration components
   LevelData<FArrayBox> cellVelCFG(grids, CFG_DIM, IntVect::Zero);
   computeCellVelCfgComp( cellVelCFG, phase_geom, m_field );
   
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


void
MagnetizationKernel::eval( LevelData<FArrayBox>& a_result,
                           const KineticSpecies& a_kinetic_species ) const
{
   
   const PhaseGeom& geometry = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();
   
   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[dit]);
      
      FArrayBox& this_result = a_result[dit];
      
      // Get the physical velocity coordinates for this part of phase space
      FArrayBox velocityRealCoords(this_result.box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);
      
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


MaxwellianKernel::MaxwellianKernel (const CFG::LevelData<CFG::FArrayBox>& density,
                                    const CFG::LevelData<CFG::FArrayBox>& temperature,
                                    const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift)
   : m_density(density),
     m_temperature(temperature),
     m_v_parallel_shift(v_parallel_shift)

{
}


void
MaxwellianKernel::eval( LevelData<FArrayBox>& a_result,
                        const KineticSpecies& a_kinetic_species ) const
{
   const PhaseGeom& phase_geom = a_kinetic_species.phaseSpaceGeometry();
   const LevelData<FArrayBox>& dfn = a_kinetic_species.distributionFunction();

   const DisjointBoxLayout& grids = dfn.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();


   //Compute Maxwellian distribution with the 3 first moments
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


void
Kernel::computeCellVelCfgComp(LevelData<FArrayBox>& a_velCellCfgComp,
                              const PhaseGeom& a_phase_geom,
                              const LevelData<FluxBox>& a_field,
                              const bool a_mag_drifts_only     ) const
{
   const DisjointBoxLayout & grids = a_velCellCfgComp.getBoxes();
   LevelData<FluxBox> pointwiseFaceVel(grids, SpaceDim, IntVect::Unit);
   a_phase_geom.computeGKVelocities(a_field, pointwiseFaceVel, a_mag_drifts_only);

   for (DataIterator dit(a_velCellCfgComp.dataIterator()); dit.ok(); ++dit) {
     a_velCellCfgComp[dit].setVal(0.);
     for (int dir=0; dir<CFG_DIM; dir++) {
          FORT_COMPUTE_VEL_CELL(CHF_FRA(a_velCellCfgComp[dit]),
                                CHF_FRA(pointwiseFaceVel[dit][dir]),
                                CHF_CONST_INT(dir),
                                CHF_BOX(a_velCellCfgComp[dit].box()));
     }
   }
}

#include "NamespaceFooter.H"
