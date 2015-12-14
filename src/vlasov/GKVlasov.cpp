#include "GKVlasov.H"
#include "Directions.H"

#include "FourthOrderUtil.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "MomentOp.H"

#include "mappedLimiterF_F.H"
#include "altFaceAverages.H"
#include "mappedAdvectionFlux.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "SingleNullCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

using std::max;

const char* GKVlasov::pp_name = {"gkvlasov"};

Real GKVlasov::s_stability_bound[NUM_FLUX] = {2.06,2.7852,1.7453,1.7320,1.7320,1.7453};

GKVlasov::GKVlasov( ParmParse& a_pp,
                    const Real a_larmor_number )
  : m_larmor_number(a_larmor_number),
    m_saved_dt(-1.0),
    m_face_avg_type(INVALID),
    m_dt_dim_factor(1.0)
{
   if (a_pp.contains("limiter")) {
      if ( procID()==0 ) MayDay::Warning("GKVlasov: Use of input flag 'limiter' deprecated");
   }

   if (a_pp.contains("verbose") && procID()==0) {
      a_pp.get("verbose", m_verbose);
   }
   else {
      m_verbose = false;
   }

   if (a_pp.contains("time_step_diagnostics")) {
      a_pp.get("time_step_diagnostics", m_time_step_diagnostics);
   }
   else {
      m_time_step_diagnostics = false;
   }

   if (a_pp.contains("face_avg_type")) {
      std::string dummy;
      a_pp.get("face_avg_type", dummy);
      if (dummy.compare("uw1")==0) {
         m_face_avg_type = UW1;
      }
      else if (dummy.compare("uw3")==0) {
         m_face_avg_type = UW3;
      }
      else if (dummy.compare("uw5")==0) {
         m_face_avg_type = UW5;
      }
      else if (dummy.compare("weno5")==0) {
         m_face_avg_type = WENO5;
      }
      else if (dummy.compare("bweno")==0) {
         m_face_avg_type = BWENO;
      }
   }
   else {
      if ( procID()==0 ) MayDay::Warning("Flux type unspecified; using default (BWENO)");
      m_face_avg_type = BWENO;
   }
   // In theory, advection is in (PDIM-1) dimensions, so we could relax this
   // a little, but for now, let's be conservative.
   m_dt_dim_factor = (m_face_avg_type>PPM) ? sqrt(PDIM) : 1.0;
}



GKVlasov::~GKVlasov()
{
}


void
GKVlasov::evalRHS( KineticSpecies&           a_rhs_species,
                   const KineticSpecies&     a_soln_species,
                   const LevelData<FluxBox>& a_Efield,
                   const Real                a_time )
{
   /*
     Evaluates the (negated) phase space divergence:
        rhs = - divergence_R ( R_dot soln ) - divergence_v ( v_dot soln )
   */
   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& dbl( soln_dfn.getBoxes() );

   LevelData<FluxBox> velocity( dbl, SpaceDim, IntVect::Unit );
   a_soln_species.computeVelocity( velocity, a_Efield );

   const PhaseGeom& geometry( a_rhs_species.phaseSpaceGeometry() );
   LevelData<FluxBox> flux( dbl, SpaceDim, IntVect::Zero );
   computeFlux( soln_dfn, velocity, flux, geometry );

   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   const bool OMIT_NT(false);
   geometry.mappedGridDivergence( rhs_dfn, flux, OMIT_NT );

   // Divide by cell volume and negate
   for (DataIterator dit( rhs_dfn.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys&
         block_coord_sys( geometry.getBlockCoordSys( dbl[dit] ) );
      double fac( -1.0 / block_coord_sys.getMappedCellVolume() );
      rhs_dfn[dit].mult( fac );
   }
}



void
GKVlasov::evalRHS( KineticSpecies&                 a_rhs_species,
                   double&                         a_lo_value,
                   double&                         a_hi_value,
                   CFG::LevelData<CFG::FArrayBox>& a_radial_flux_divergence_average,
                   const KineticSpecies&           a_soln_species,
                   const LevelData<FluxBox>&       a_Efield,
                   const Real                      a_time )
{
   /*
     Evaluates the (negated) phase space divergence:
     rhs = - divergence_R ( R_dot soln ) - divergence_v ( v_dot soln )
   */
   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& dbl( soln_dfn.getBoxes() );
    
   LevelData<FluxBox> velocity( dbl, SpaceDim, IntVect::Unit );
   a_soln_species.computeVelocity( velocity, a_Efield );
    
   const PhaseGeom& geometry( a_rhs_species.phaseSpaceGeometry() );
   LevelData<FluxBox> flux( dbl, SpaceDim, IntVect::Zero );
   computeFlux( soln_dfn, velocity, flux, geometry );
    
   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   const bool OMIT_NT(false);
   geometry.mappedGridDivergence( rhs_dfn, flux, OMIT_NT );
    
   // Divide by cell volume and negate
   for (DataIterator dit( rhs_dfn.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys&
         block_coord_sys( geometry.getBlockCoordSys( dbl[dit] ) );
      double fac( -1.0 / block_coord_sys.getMappedCellVolume() );
      rhs_dfn[dit].mult( fac );
   }
    
   computeRadialFluxDivergence(geometry, flux, a_soln_species.mass(),
                               a_soln_species.charge(), a_lo_value, a_hi_value,
                               a_radial_flux_divergence_average );
}


void
GKVlasov::computeRadialFluxDivergence(const PhaseGeom&                a_geometry,
                                      LevelData<FluxBox>&             a_flux,
                                      double                          a_mass,
                                      double                          a_charge,
                                      double&                         a_lo_value,
                                      double&                         a_hi_value,
                                      CFG::LevelData<CFG::FArrayBox>& a_radial_flux_divergence_average) const
{
   const DisjointBoxLayout& grids( a_flux.getBoxes() );
   const MultiBlockCoordSys* coords = a_geometry.coordSysPtr();
    
   LevelData<FluxBox> flux_even(grids, SpaceDim, IntVect::Zero);
   LevelData<FluxBox> flux_odd(grids, SpaceDim, IntVect::Zero);
    
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      flux_even[dit].copy(a_flux[dit]);
      flux_odd[dit].copy(a_flux[dit]);
   }
    
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
         FArrayBox& this_flux_even_dir = flux_even[dit][dir];
         FArrayBox& this_flux_odd_dir = flux_odd[dit][dir];
         int block_number = coords->whichBlock(grids[dit]);
            
         if (dir == RADIAL_DIR) {
            if ( block_number < 2 ) {
               Box box( this_flux_even_dir.box() );
               BoxIterator bit(box);
                    
               for (bit.begin(); bit.ok(); ++bit) {
                  IntVect iv = bit();
                  for (int n=0; n<SpaceDim; n=n+1) {
                     if (iv[0]%2 == 0 ) {
                        this_flux_even_dir(iv,n) = 0.0;
                     }
                     else {
                        this_flux_odd_dir(iv,n) = 0.0;
                     }
                  }
               }
            }
            else {
               this_flux_even_dir.setVal(0.);
               this_flux_odd_dir.setVal(0.);
            }
         }
            
         else {
            this_flux_even_dir.setVal(0.);
            this_flux_odd_dir.setVal(0.);
         }
      }
   }
    
   // Create fake species objects so that we can use the MomentOp mechanism
   // to integrate their distribution functions (set to the phase space flux
   // divergence) over velocity space
   KineticSpecies temp_even("temp_even", a_mass, a_charge, a_geometry);
   LevelData<FArrayBox>& phase_divergence_even = temp_even.distributionFunction();
   phase_divergence_even.define(grids, 1, IntVect::Zero);
   KineticSpecies temp_odd("temp_odd", a_mass, a_charge, a_geometry);
   LevelData<FArrayBox>& phase_divergence_odd = temp_odd.distributionFunction();
   phase_divergence_odd.define(grids, 1, IntVect::Zero);
    
   a_geometry.mappedGridDivergence( phase_divergence_even, flux_even, false );
   a_geometry.mappedGridDivergence( phase_divergence_odd, flux_odd, false );
    
   LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
   a_geometry.getCellVolumes(volume);
    
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phase_divergence_even[dit] /= volume[dit];
      phase_divergence_odd[dit] /= volume[dit];
   }

   MomentOp& moment_op = MomentOp::instance();

   CFG::LevelData<CFG::FArrayBox> config_divergence_even( a_geometry.magGeom().grids(), 1, CFG::IntVect::Zero);
   moment_op.compute( config_divergence_even, temp_even, ChargeDensityKernel() );

   CFG::LevelData<CFG::FArrayBox> config_divergence_odd( a_geometry.magGeom().grids(), 1, CFG::IntVect::Zero);
   moment_op.compute( config_divergence_odd, temp_odd, ChargeDensityKernel() );

   const CFG::MagGeom& config_geom = a_geometry.magGeom();
   const CFG::MagCoordSys* config_coord_sys = config_geom.getCoordSys();
   const CFG::DisjointBoxLayout& config_grids = config_geom.grids();

   CFG::FluxSurface fs(config_geom);

   CFG::LevelData<CFG::FArrayBox> fs_average_even(config_grids, 1, CFG::IntVect::Zero);
   fs.averageAndSpread(config_divergence_even, fs_average_even);

   CFG::LevelData<CFG::FArrayBox> fs_average_odd(config_grids, 1, CFG::IntVect::Zero);
   fs.averageAndSpread(config_divergence_odd, fs_average_odd);
    
   for (CFG::DataIterator dit(config_grids); dit.ok(); ++dit) {
      CFG::Box box( config_grids[dit] );
      CFG::BoxIterator bit(box);
      for (bit.begin(); bit.ok(); ++bit) {
         CFG::IntVect iv = bit();
         if (iv[0]%2 == 0 ) {
            a_radial_flux_divergence_average[dit](iv,0) = -fs_average_odd[dit](iv,0);
         }
         else {
            a_radial_flux_divergence_average[dit](iv,0) = -fs_average_even[dit](iv,0);
         }
      }
   }
    
   const PhaseBlockCoordSys& block0_coord_sys = (const PhaseBlockCoordSys&)(*(coords->getCoordSys(0)));
   const Box& block0_domain_box = block0_coord_sys.domain().domainBox();

   int lower_core_radial_index = block0_domain_box.smallEnd(RADIAL_DIR);
   int upper_core_radial_index = block0_domain_box.bigEnd(RADIAL_DIR);

   a_lo_value = a_hi_value = -DBL_MAX;

   for (CFG::DataIterator dit(config_grids); dit.ok(); ++dit) {
      int block_number = config_coord_sys->whichBlock(config_grids[dit]);

      if ( block_number < 2 ) {
         const CFG::Box& box = config_grids[dit];

         int radial_index = box.smallEnd(RADIAL_DIR);
         if ( radial_index == lower_core_radial_index ) {
            if (radial_index%2 == 0) {
               a_lo_value = -fs_average_odd[dit](box.smallEnd());
            }
            else {
               a_lo_value = -fs_average_even[dit](box.smallEnd());
            }
         }
         radial_index = box.bigEnd(RADIAL_DIR);
         if ( radial_index == upper_core_radial_index ) {
            if (radial_index%2 == 0) {
               a_hi_value = fs_average_even[dit](box.bigEnd());
            }
            else {
               a_hi_value = fs_average_odd[dit](box.bigEnd());
            }
         }
      }
   }

   a_lo_value = globalMax(a_lo_value);

   double fac = ( typeid(*(config_geom.getCoordSys())) != typeid(CFG::SingleNullCoordSys) ) ? 1.0 : 2.0 ;
   a_hi_value = fac * globalMax(a_hi_value); //NB: factor of two is needed to compensate for the flux averaging at block boundaries

}

void
GKVlasov::computeFlux( const LevelData<FArrayBox>& a_dist_fn,
                       const LevelData<FluxBox>&   a_velocity,
                       LevelData<FluxBox>&         a_flux,
                       const PhaseGeom&            a_phase_geom )
{

   /*
      Compute the phase space flux given the input phase space advection
      velocity.   This is where the hyperbolic stuff connects.
   */

   // Construct appropriately accurate face-averages of phi and advVel
   LevelData<FluxBox> faceDist(a_dist_fn.getBoxes(), a_dist_fn.nComp(), a_dist_fn.ghostVect() );

   // If we're limiting the face-centered values, do it here
   if (m_face_avg_type==PPM) {
      computeFaceAverages( faceDist, a_dist_fn, a_phase_geom.secondOrder() );
      faceDist.exchange();
      applyMappedLimiter( faceDist, a_dist_fn, a_velocity, a_phase_geom );
   }
   else if (m_face_avg_type==UW1) {
      uw1FaceAverages( faceDist, a_dist_fn, a_velocity, a_phase_geom );
   }
   else if (m_face_avg_type==UW3) {
      uw3FaceAverages( faceDist, a_dist_fn, a_velocity, a_phase_geom );
   }
   else if (m_face_avg_type==UW5) {
      uw5FaceAverages( faceDist, a_dist_fn, a_velocity, a_phase_geom );
   }
   else if (m_face_avg_type==WENO5) {
      weno5FaceAverages( faceDist, a_dist_fn, a_velocity, a_phase_geom );
   }
   else if (m_face_avg_type==BWENO) {
      bwenoFaceAverages( faceDist, a_dist_fn, a_velocity, a_phase_geom );
   }

   if ( a_phase_geom.secondOrder() ) {

      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxes( a_flux, faceDist, a_velocity, false );
   }
   else {

      // Fill transverse ghosts at the block boundaries in preparation for
      // computing the flux using the fourth-order product formula
      a_phase_geom.fillTransverseGhosts(faceDist, false);

      // Compute computational-space fluxes; in mappedAdvectionFlux.cpp
      computeCompFaceFluxes( a_flux, faceDist, a_velocity, true );
   }

   a_flux.exchange();
}



Real
GKVlasov::computeDt( const LevelData<FluxBox>&    a_Efield,
                     const KineticSpeciesPtrVect& a_species_vect )
{
   Real dt(BASEFAB_REAL_SETVAL);

   for (int s(0); s<a_species_vect.size(); s++) {

      const KineticSpecies& species( *(a_species_vect[s]) );
      const LevelData<FArrayBox>& dfn( species.distributionFunction() );

      LevelData<FluxBox> velocity( dfn.getBoxes(), SpaceDim, IntVect::Unit );
      species.computeMappedVelocity( velocity, a_Efield );

      const Real UNIT_CFL(1.0);
      const PhaseGeom& geometry( species.phaseSpaceGeometry() );
      Real speciesDt( computeMappedDtSpecies( velocity, geometry, UNIT_CFL ) );
      CH_assert(speciesDt >= 0);

      dt = Min( dt, speciesDt );
   }

   return dt;
}

Real
GKVlasov::computeTimeScale( const LevelData<FluxBox>&    a_Efield,
                            const KineticSpeciesPtrVect& a_species_vect )
{
   Real dt(BASEFAB_REAL_SETVAL);

   for (int s(0); s<a_species_vect.size(); s++) {

      const KineticSpecies& species( *(a_species_vect[s]) );
      const LevelData<FArrayBox>& dfn( species.distributionFunction() );

      LevelData<FluxBox> velocity( dfn.getBoxes(), SpaceDim, IntVect::Unit );
      species.computeMappedVelocity( velocity, a_Efield );

      const Real UNIT_CFL(1.0);
      const PhaseGeom& geometry( species.phaseSpaceGeometry() );
      Real speciesDt( computeMappedTimeScaleSpecies( velocity, geometry) );
      CH_assert(speciesDt >= 0);

      dt = Min( dt, speciesDt );
   }

   return dt;
}


void
GKVlasov::initialize( KineticSpeciesPtrVect& soln,
                      const Real      time )
{
   MayDay::Error( "not implemented!" );
}



inline void
GKVlasov::computeFaceAverages( LevelData<FluxBox>&         a_face_data,
                               const LevelData<FArrayBox>& a_cell_data,
                               const bool                  a_second_order  ) const
{
   if ( a_second_order) {
      CellToEdge( a_cell_data, a_face_data );
   }
   else {
      fourthOrderCellToFace( a_face_data, a_cell_data );
   }
}


#undef PLOT_STABLEDT

Real
GKVlasov::computeMappedDtSpecies( const LevelData<FluxBox>& a_faceVel,
                                  const PhaseGeom&          a_geom,
                                  Real                      a_cfl )
{
   const DisjointBoxLayout grids = a_faceVel.getBoxes();
   CH_assert(grids == a_geom.gridsFull());

   MultiBlockCoordSys* coords = a_geom.coordSysPtr();

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

   LevelData<FluxBox> faceNormalVel(grids, 1);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coords->whichBlock(grids[dit]);
      const PhaseBlockCoordSys* block_coord_sys
         = static_cast<const PhaseBlockCoordSys*>(coords->getCoordSys(block_number));
      RealVect face_area = block_coord_sys->getMappedFaceArea();

      const FluxBox& thisFaceVel = a_faceVel[dit];
      FluxBox& thisNormalVel = faceNormalVel[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& thisNormalVel_dir = thisNormalVel[dir];
         thisNormalVel_dir.copy(thisFaceVel[dir],dir,0,1);
         thisNormalVel_dir *= face_area[dir];
      }
   }

#ifdef PLOT_STABLEDT
   LevelData<FArrayBox> stableDt(grids, 1);
#endif

   // now average back to cell-centers and divide by cell volumes.
   // instead of averaging face->cell we pick
   // the max absolute value on the two opposing faces.
   LevelData<FArrayBox> cellVolumes(grids, 1, IntVect::Zero);
   a_geom.getCellVolumes(cellVolumes);

   struct {
      double val;
      int rank;
   } pair_in, pair_out;

   Real maxVelLoc = 0.;
   int maxblockLoc(-1);
   IntVect maxindLoc;
   int maxDirLoc(-1);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coords->whichBlock(grids[dit]);

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
         if (m_time_step_diagnostics) {
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

#ifdef PLOT_STABLEDT
      stableDt[dit].copy(cellVel);
#endif
   }

   RealVect X;
   if (m_time_step_diagnostics) {
      RealVect dx = coords->getCoordSys(maxblockLoc)->dx();
      RealVect xi = dx*maxindLoc;
      xi += 0.5*dx;
      X = coords->getCoordSys(maxblockLoc)->realCoord(xi);
   }

   Real maxVel = maxVelLoc;
   IntVect maxind = maxindLoc;
   int maxDir = maxDirLoc;
#ifdef CH_MPI
   if (m_time_step_diagnostics) {
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

   if (maxVel > smallVal)
      {
         // 2.06 factor per our JCP paper, equation (80), however approximate
         newDt = a_cfl * s_stability_bound[m_face_avg_type] / (maxVel * m_dt_dim_factor);
         //       newDt = a_cfl * 2.06 / maxVel;
      }

   if (m_time_step_diagnostics && procID()==0) {
#if PDIM==4
      cout << "Vlasov operator time step was limited by the velocity at (R,Z,vparallel,mu) = " << X << " and mapped coordinate = " << maxind << endl;
#endif
#if PDIM==5
      cout << "Vlasov operator time step was limited by the velocity at (R,phi,Z,vparallel,mu) = " << X << " and mapped coordinate = " << maxind << endl;
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
#if PDIM==5
         case TOROIDAL_DIR:
            cout << "toroidal";
            break;
#endif
         case VPARALLEL_DIR:
            cout << "vparallel";
            break;
         case MU_DIR:
            cout << "mu";
            break;
         }
      cout << " direction makes the largest contribution to the stable dt reciprocal at that point" << endl;
   }

#ifdef PLOT_STABLEDT
   a_geom.plotAtVelocityIndex( "stabledt", a_geom.vel_restrict(maxind), stableDt );
   exit(1);
#endif

   return newDt;
}

Real
GKVlasov::computeMappedTimeScaleSpecies( const LevelData<FluxBox>& a_faceVel,
                                         const PhaseGeom&          a_geom )
{
   const DisjointBoxLayout grids = a_faceVel.getBoxes();
   CH_assert(grids == a_geom.gridsFull());

   MultiBlockCoordSys* coords = a_geom.coordSysPtr();

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

   LevelData<FluxBox> faceNormalVel(grids, 1);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coords->whichBlock(grids[dit]);
      const PhaseBlockCoordSys* block_coord_sys
         = static_cast<const PhaseBlockCoordSys*>(coords->getCoordSys(block_number));
      RealVect face_area = block_coord_sys->getMappedFaceArea();

      const FluxBox& thisFaceVel = a_faceVel[dit];
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
   a_geom.getCellVolumes(cellVolumes);

   struct {
      double val;
      int rank;
   } pair_in, pair_out;

   Real maxVelLoc = 0.;
   int maxblockLoc(-1);
   IntVect maxindLoc;
   int maxDirLoc(-1);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      int block_number = coords->whichBlock(grids[dit]);

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
         if (m_time_step_diagnostics) {
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
   if (m_time_step_diagnostics) {
      RealVect dx = coords->getCoordSys(maxblockLoc)->dx();
      RealVect xi = dx*maxindLoc;
      xi += 0.5*dx;
      X = coords->getCoordSys(maxblockLoc)->realCoord(xi);
   }

   Real maxVel = maxVelLoc;
   IntVect maxind = maxindLoc;
   int maxDir = maxDirLoc;
#ifdef CH_MPI
   if (m_time_step_diagnostics) {
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
   if (maxVel > smallVal) newDt = 1.0 / maxVel;

   return newDt;
}



void
GKVlasov::applyMappedLimiter( LevelData<FluxBox>&         a_facePhi,
                              const LevelData<FArrayBox>& a_cellPhi,
                              const LevelData<FluxBox>&   a_faceVel,
                              const PhaseGeom&            a_geom )
{
   int nComp = a_facePhi.nComp();

   // this specifies the number of ghost faces we will need in
   // each transverse and normal direction
   int transverseGrow = 2;
   int normalGrow = 0;

   // use C value from Colella and Sekora
   Real limiterConst =  1.25;

   // may need to do exchange on facePhi
   a_facePhi.exchange();

   const DisjointBoxLayout& grids = a_facePhi.getBoxes();

   // in order to do upwinding, need normal velocities in
   // computational space. computing Fourth-order face
   // averages just in order to upwinding is over-working the
   // issue, but it's the most convenient thing to do here.
   LevelData<FluxBox> normalVel(grids, 1, a_faceVel.ghostVect());

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = a_geom.getBlockCoordSys(grids[dit]);
      RealVect face_area = block_coord_sys.getMappedFaceArea();

      const FluxBox& thisFaceVel = a_faceVel[dit];
      FluxBox& thisNormalVel = normalVel[dit];
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& thisNormalVel_dir = thisNormalVel[dir];
         thisNormalVel_dir.copy(thisFaceVel[dir],dir,0,1);
         thisNormalVel_dir *= face_area[dir];
      }
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = a_geom.getBlockCoordSys(grids[dit]);
      const RealVect& dx = block_coord_sys.dx();

      const Box& gridBox = grids[dit];
      Box lapBox(gridBox);
      // need two laplacian ghost cells for this limiter
      int lapBoxGrow = max(normalGrow+2, transverseGrow);
      lapBox.grow(lapBoxGrow);
      FArrayBox ccLaplacians(lapBox, nComp);
      FluxBox& thisFacePhi = a_facePhi[dit];
      const FArrayBox& thisCellPhi = a_cellPhi[dit];
      const FluxBox& thisNormalVel = normalVel[dit];
      // check to be sure that we have enough phi to compute
      // all of these Laplacians
      {
         Box LapPhiBox(lapBox);
         LapPhiBox.grow(1);
         CH_assert(thisCellPhi.box().contains(LapPhiBox));
      }
      
      for (int dir=0; dir<SpaceDim; dir++)
         {
            // box of valid edges for this grid
            Box faceBox(gridBox);
            // for 4th order, we need extra faces in the transverse
            // directions. (for mapped-grid case, we actually need
            // _2_ transverse faces

            // need extra faces in the tangential directions in order
            // to handle 4th-order multiplies
            faceBox.grow(transverseGrow);
            // however, need different number of extra faces in normal dir
            faceBox.grow(dir,normalGrow - transverseGrow);

            faceBox.surroundingNodes(dir);
            FArrayBox& thisFacePhiDir = thisFacePhi[dir];
            const FArrayBox& thisNormalVelDir = thisNormalVel[dir];
            {
               // context for computing limited face values
               Box grownFaceBox(faceBox);
               // need an extra face's worth of the FC laplacians
               grownFaceBox.grow(dir,1);

               FArrayBox centeredLaplacian(grownFaceBox, nComp);

               // compute centered Laplacian
               centeredLaplacian.setVal(0.0);
               FORT_CENTEREDLAPLACIAN(CHF_FRA(centeredLaplacian),
                                      CHF_CONST_FRA(thisFacePhiDir),
                                      CHF_CONST_FRA(thisCellPhi),
                                      CHF_BOX(grownFaceBox),
                                      CHF_CONST_INT(dir),
                                      CHF_CONST_REAL(dx[dir]));

               // compute cell-centered Laplacians
               ccLaplacians.setVal(0.0);
               FORT_CCLAPLACIAN(CHF_FRA(ccLaplacians),
                                CHF_CONST_FRA(thisCellPhi),
                                CHF_BOX(lapBox),
                                CHF_CONST_INT(dir),
                                CHF_CONST_REAL(dx[dir]));

               // now compute limited face value
               FORT_LIMITFACEVALUES(CHF_FRA(thisFacePhiDir),
                                    CHF_CONST_FRA(thisCellPhi),
                                    CHF_CONST_FRA(centeredLaplacian),
                                    CHF_CONST_FRA(ccLaplacians),
                                    CHF_BOX(grownFaceBox),
                                    CHF_CONST_INT(dir),
                                    CHF_CONST_REAL(dx[dir]),
                                    CHF_CONST_REAL(limiterConst));
               // end context for computing limited face values
            } // (centeredLaplacian goes out of scope)

            // now compute parabolic interpolants
            // storage for cell-centered D^2a
            Box lapBoxDir(gridBox);
            // need transverse faces for 4th order
            lapBoxDir.grow(transverseGrow);
            // need this to be grown by one in normal dir
            lapBoxDir.grow(dir,normalGrow-transverseGrow+1);

            FArrayBox D2a(lapBoxDir,nComp);
            FArrayBox D2aLim(lapBoxDir,nComp);

            // initialize D2a to be -a_6/(h^2)
            // first compute a_6...
            FORT_COMPUTEA6(CHF_FRA(D2a),
                           CHF_CONST_FRA(thisCellPhi),
                           CHF_CONST_FRA(thisFacePhiDir),
                           CHF_BOX(lapBoxDir),
                           CHF_CONST_INT(dir));

            // then multiply by -2/(h^2)
            // (dfm 5/5/09 -- was missing the factor of 2)
            Real mult = -2.0/(dx[dir]*dx[dir]);
            D2a *= mult;

            // now limit D2a w/r/t other 2nd-derivatives
            FORT_LIMITD2A(CHF_FRA(D2aLim),
                          CHF_CONST_FRA(D2a),
                          CHF_CONST_FRA(ccLaplacians),
                          CHF_BOX(lapBoxDir),
                          CHF_CONST_INT(dir),
                          CHF_CONST_REAL(limiterConst));

            // storage for left and right interpolants
            //   note that while leftPhi and rightPhi are associated
            //   with face indicies, it is more natural to compute the
            //   PPM limited face values on a per-cell basis.  Thus,
            //   each of these FArrayBoxes is one cell longer than it
            //   need to be
            Box growBox( faceBox );
            growBox.grow(dir,1);
            FArrayBox leftPhi(growBox, nComp);
            FArrayBox rightPhi(growBox, nComp);

            // We operate on the cells of the domain plus one ghost on each
            // end in the current direction
            Box cellBox( growBox );
            cellBox.growHi(dir,-1);
            FORT_LEFTRIGHTSTATES(CHF_FRA(leftPhi),
                                 CHF_FRA(rightPhi),
                                 CHF_CONST_FRA(thisFacePhiDir),
                                 CHF_CONST_FRA(thisCellPhi),
                                 CHF_CONST_FRA(D2aLim),
                                 CHF_CONST_FRA(D2a),
                                 CHF_BOX(cellBox),
                                 CHF_CONST_INT(dir));

            // need to do something about boundary values here

            // now pick upwind state
            // note that normalVel only has a single component, so we
            // use 0 as the component argument for the CHF_FRA1 macro
            // as a debugging check, setThisFacePhiDir to a bogus value
            // first
            thisFacePhiDir.setVal(1000000000);
            Box overlap( faceBox & thisNormalVelDir.box() );
            FORT_SELECTUPWIND(CHF_FRA(thisFacePhiDir),
                              CHF_CONST_FRA(leftPhi),
                              CHF_CONST_FRA(rightPhi),
                              CHF_CONST_FRA1(thisNormalVelDir,0),
                              CHF_BOX(overlap));

         } // end loop over directions

   } // end loop over grids
}


double
GKVlasov::globalMax( const double a_data ) const
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
