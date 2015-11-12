#include "Advect.H"
#include "AdvectF_F.H"
#include "CONSTANTS.H"

#include "FourthOrderUtil.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"

#include "MagBlockCoordSys.H"
#include "altFaceAverages.H"
#include "mappedAdvectionFlux.H"
#include "newMappedGridIO.H"
#include "BlockRegister.H"
#include "inspect.H"

#include "NamespaceHeader.H"

using std::max;

char* Advect::pp_name = {"Advect"};

Advect::Advect( ParmParse pp )
  : m_limit_face_values(true),
    m_saved_dt(-1.0),
    m_flux_type(INVALID)
{
   if (pp.contains("verbose") && procID()==0) {
      pp.get("verbose", m_verbose);
   }
   else {
      m_verbose = false;
   }

   Vector<Real> velocity_array;
   if (pp.contains("velocity"))  {
     pp.queryarr("velocity", velocity_array, 0, 2);
     m_velocity = RealVect(velocity_array);
   }
   else {
     MayDay::Error("Velocity vector unspecified");
   }

   if (pp.contains("flux_type")) {
      std::string dummy;
      pp.get("flux_type", dummy);
      if (dummy.compare("uw5")==0) {
         m_flux_type = UW5;
      }
      else if (dummy.compare("weno5")==0) {
         m_flux_type = WENO5;
      }
      else if (dummy.compare("bweno")==0) {
         m_flux_type = BWENO;
      }
   }
   else {
      MayDay::Warning("Flux type unspecified; using default (PPM)");
      m_flux_type = PPM;
   }
   
}



Advect::~Advect()
{
}



void Advect::evalRHS( const MagGeom&              a_geom,
                      LevelData<FArrayBox>&       a_rhs,
                      const LevelData<FArrayBox>& a_soln,
                      const Real                  a_time )
{
  CH_assert(a_soln.ghostVect() == 4*IntVect::Unit);

   /*
     Evaluates the (negated) divergence:

     rhs = - divergence ( velocity * soln )

   */

   Real localDt = -1.0;

   /*
     Compute the advection velocity .
   */

   // Compute the advection velocity.
   const DisjointBoxLayout & dbl = a_soln.getBoxes();
   LevelData<FluxBox> velocity(dbl, SpaceDim, IntVect::Unit);
   updateVelocity( a_geom, velocity, false );

   // since we've got the advection velocities here, this is a
   // good place to compute dt
   // use a cfl of 1 here, and we can apply the real cfl elsewhere.
   Real fakeCFL = 1.0;
   LevelData<FluxBox> mapped_velocity(dbl, SpaceDim, IntVect::Unit);
   updateVelocity( a_geom, mapped_velocity, true );
   Real speciesDt = computeMappedDtSpecies( mapped_velocity, a_geom, fakeCFL );

   CH_assert(speciesDt > 0);

   if (localDt < 0) {
      localDt = speciesDt;
   }
   else if (speciesDt < localDt) {
      localDt = speciesDt;
   }

   //   BCInterface& bc = rhs_species.boundaryConditions();
   //   bc.ghostCellBC(dfn, geom, 0, time);

   //   LevelData<FluxBox> flux( dbl, SpaceDim, IntVect::Zero );
   LevelData<FluxBox> flux( dbl, SpaceDim, IntVect::Unit );
   computeFlux( a_soln, velocity, flux, a_geom );

   // set boundary conditions on flux, if appropriate
   //bc.fluxBC(flux, velocity, dfn, geom, 0, time);

   // Compute the flux divergence
   a_geom.computeMappedGridDivergence( flux, a_rhs, true );

   // Divide by cell volume and move the average flux divergence to the right-hand side
   DataIterator it = dbl.dataIterator();
   for (it.begin(); it.ok(); ++it) {
      const MagBlockCoordSys& block_coord_sys( a_geom.getBlockCoordSys( dbl[it] ) );

      double fac = -1. / block_coord_sys.getMappedCellVolume();

      a_rhs[it].mult(fac);
   }

   // set dt_save
   m_saved_dt = localDt;
}


void Advect::computeFlux( const LevelData<FArrayBox>& dist_fn,
                          const LevelData<FluxBox>&   velocity,
                          LevelData<FluxBox>&         flux,
                          const MagGeom&  a_geom)
{
  CH_assert(dist_fn.ghostVect() == 4*IntVect::Unit);

   /*
      Compute the flux given the input advection
      velocity.   This is where the hyperbolic stuff connects.
   */

   const DisjointBoxLayout& grids = dist_fn.getBoxes();
   int ncomp = dist_fn.nComp();

   CH_assert(velocity.getBoxes() == grids);
   CH_assert(flux.getBoxes() == grids);

   // Compute computational-grid cell averages
   LevelData<FArrayBox> dist_cg(grids, ncomp, dist_fn.ghostVect());

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      dist_cg[dit].copy(dist_fn[dit]);
   }

   // construct appropriately accurate face-averages of phi and advVel
   LevelData<FluxBox> faceDist(grids, dist_cg.nComp(), dist_cg.ghostVect() );

   // if we're limiting the face-centered values, do it here
   if (m_limit_face_values) {
      if (m_flux_type==PPM) {
         computeFaceAverages( faceDist, dist_cg );
         faceDist.exchange();
         MultiBlockCoordSys* coords = a_geom.coordSysPtr();
         applyMappedLimiter( faceDist, dist_cg, velocity, *coords );
      }
      else if (m_flux_type==UW5) {
         uw5FaceAverages( faceDist, dist_fn, velocity, a_geom );
      }
      else if (m_flux_type==WENO5) {
         weno5FaceAverages( faceDist, dist_fn, velocity, a_geom );
      }
      else if (m_flux_type==BWENO) {
         bwenoFaceAverages( faceDist, dist_fn, velocity, a_geom );
      }

#if 0

   // The face distribution values are now only good on two layers of
   // ghost cells due to the fourth-order interpolation just performed.

   // The alternative limiter needs another layer of velocity ghost cells.
   // Can we avoid this copy somehow?

//   LevelData<FluxBox> grown_velocity(velocity.disjointBoxLayout(),
//                                     velocity.nComp(),
//                                     velocity.ghostVect() + IntVect::Unit);

//   for (dit.begin(); dit.ok(); ++dit) {
//      for ( int dir=0; dir<SpaceDim; ++dir ) {
//         grown_velocity[dit][dir].setVal(0.);
//         grown_velocity[dit][dir].copy(velocity[dit][dir]);
//      }
//   }
//   grown_velocity.exchange();

#endif
   }
   else {
      computeFaceAverages( faceDist, dist_fn );
      faceDist.exchange();
   }
   
   // compute computational-space fluxes -
   // The fourth-order product formula is used here, so faceDist and
   // velocity need to have data in one more ghost cell than flux.
   // flux has zero, velocity has data in at least one and
   // faceDist has data in two.
   computeCompFaceFluxes( flux, faceDist, velocity, true );
      
   // Extrapolate to transverse cell faces at boundary
   for (dit.begin(); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys( a_geom.getBlockCoordSys(grids[dit]) );
      const ProblemDomain& domain( block_coord_sys.domain() );
      
      FluxBox& this_F( flux[dit] );
      for (int dir(0); dir<SpaceDim; dir++) {
         Box box_dir( surroundingNodes( grids[dit], dir) );
         secondOrderTransExtrapAtDomainBdry( this_F[dir], dir, box_dir, domain );
      }
   }
   
   flux.exchange();
}



Real Advect::computeDt( const MagGeom&              a_geom,
                        const LevelData<FArrayBox>& a_soln )
{
  Real dt = 0;

  // if we've already computed a dt in an earlier flux computation:
  if (m_saved_dt > 0)
    {
      dt = m_saved_dt;
    }
  else
    {
      // at the intial time, before we've done any work
      // this code is more or less copied from the
      // evalRHS function
      dt = -1;

      /*
        Compute the advection velocity.
      */
      if (m_verbose) cout << "   Computing velocities" << endl;

      LevelData<FluxBox> mapped_velocity( a_soln.getBoxes(), SpaceDim, IntVect::Unit );
      updateVelocity( a_geom, mapped_velocity, true );

      // since we've got the advection velocities here, this is a
      // good place to compute dt
      // use a cfl of 1 here, and we can apply the real cfl elsewhere.
      Real fakeCFL = 1.0;
      Real speciesDt = computeMappedDtSpecies( mapped_velocity,
                                               a_geom,
                                               fakeCFL );

      CH_assert(speciesDt >= 0);

      if (dt < 0)
         {
            dt = speciesDt;
         }
      else if (speciesDt < dt)
         {
            dt = speciesDt;
         }

    } // end if we have to compute dt from scratch here

  return dt;
}



void Advect::initialize( LevelData<FArrayBox>& soln,
                         const Real      time )
{
   MayDay::Error( "not implemented!" );
}



inline void
Advect::computeFaceAverages(LevelData<FluxBox>&         a_face_data,
                            const LevelData<FArrayBox>& a_cell_data) const
{
  bool useFourthOrder = true;
  bool useSecondOrder = !useFourthOrder;

  if (useFourthOrder)
    {
      fourthOrderCellToFace( a_face_data, a_cell_data );
    }
  else if (useSecondOrder)
    {
      CellToEdge( a_cell_data, a_face_data );
    }
  else
    {
      MayDay::Error("Bad Space Order in LevelAdvectOperator");
    }
}


#undef PLOT_STABLEDT

Real
Advect::computeMappedDtSpecies( const LevelData<FluxBox>& a_faceVel,
                                const MagGeom&            a_geom,
                                Real                      a_cfl )
{
  const DisjointBoxLayout grids = a_faceVel.getBoxes();
  CH_assert(grids == a_geom.grids());

  const MagCoordSys* coords( a_geom.getCoordSys() );

  // compute area-weighted normal velocities -- 2nd order should
  // be good enough

  LevelData<FluxBox> faceNormalVel(grids, 1);
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
     int block_number = coords->whichBlock(grids[dit]);
     const MagBlockCoordSys* block_coord_sys( coords->getCoordSys(block_number) );
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

  LevelData<FArrayBox> cellVolumes( grids, 1, IntVect::Zero );
  a_geom.getCellVolumes( cellVolumes );
  
  Real maxVelLoc = 0.0;
  for (dit.begin(); dit.ok(); ++dit) {
     FArrayBox cellVel(grids[dit], 1);
     cellVel.setVal(0.0);
     for (int dir=0; dir<SpaceDim; ++dir) {
        // average face velocity on the two faces in this direction
        int faceComp = 0;
//        int cellComp = 0;
        int cellComp = dir;
//        FArrayBox cellVelDir(grids[dit], 1);
        FArrayBox cellVelDir( grids[dit], SpaceDim );
        EdgeToCell( faceNormalVel[dit], faceComp,
                    cellVelDir, cellComp,
                    dir );
        cellVelDir.abs(dir,1);
        cellVel.plus(cellVelDir, dir, 0, 1);
//        cellVelDir.abs();
//        cellVel += cellVelDir;
     }
     
     // Divide by the cell physical volume, which is J times the computational cell volume
     cellVel.divide( cellVolumes[dit], 0, 0, 1 );
     
     Real thisMax = cellVel.norm(0,0,1);
     if (thisMax > maxVelLoc) {
        maxVelLoc = thisMax;
     }
     
#ifdef PLOT_STABLEDT
      stableDt[dit].copy(cellVel);
      stableDt[dit].invert(a_cfl*2.06);
#endif
  }

  Real maxVel = maxVelLoc;
#ifdef CH_MPI
  MPI_Allreduce(&maxVelLoc, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  Real smallVal = 1.0e-15;
  Real newDt = 0;

  if (maxVel > smallVal) {
    // 2.06 factor per our JCP paper, equation (80), however approximate
    newDt = a_cfl * 2.06 / maxVel;
  }

#ifdef PLOT_STABLEDT
  plotData( "stabledt", a_geom, stableDt );
  exit(1);
#endif

  return newDt;
}



void
Advect::applyMappedLimiter( LevelData<FluxBox>&         a_facePhi,
                            const LevelData<FArrayBox>& a_cellPhi,
                            const LevelData<FluxBox>&   a_faceVel,
                            const MultiBlockCoordSys&   a_coords)
{
}


void
Advect::updateVelocity(const MagGeom&              a_geom,
                       LevelData<FluxBox>&         a_velocity,
                       const bool                  a_mapped)
{
   CH_assert(a_velocity.ghostVect() == IntVect::Unit);

   MagCoordSys* coords = static_cast<MagCoordSys*>(a_geom.coordSysPtr());
   const DisjointBoxLayout& grids = a_velocity.disjointBoxLayout();

   LevelData<FluxBox> N_face;
   if (a_mapped) {
      N_face.define(a_velocity.disjointBoxLayout(), SpaceDim*SpaceDim, IntVect::Unit);
      a_geom.getPointwiseN(N_face);
   }

   const LevelData<FluxBox>& BFieldFace = a_geom.getFCBField();

   DataIterator dit = a_velocity.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = a_geom.getBlockCoordSys(grids[dit]);

      FluxBox& this_velocity = a_velocity[dit];

      RealVect dx = block_coord_sys.dx();

      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_velocity_dir = this_velocity[dir];
         Box vel_box_dir(this_velocity_dir.box());

#if 1
         FORT_VEL_TEST(CHF_CONST_INT(dir),
                       CHF_BOX(vel_box_dir),
                       CHF_REALVECT(m_velocity),
                       CHF_FRA(this_velocity_dir));
#else

         this_velocity_dir.copy(BFieldFace[dit][dir],0,0,1);
         this_velocity_dir.copy(BFieldFace[dit][dir],2,1,1);

#endif

         if (a_mapped) {

            // Multiply the velocity by NTranspose to convert the velocities
            // components to the computational coordinate frame
            FORT_MULT_CFG_NT(
                             CHF_BOX(vel_box_dir),
                             CHF_CONST_FRA(N_face[dit][dir]),
                             CHF_FRA(this_velocity_dir)
                             );
         }
      }
   }

   a_geom.fillTransverseGhosts(a_velocity);

   // We now have the velocity on valid cell faces plus one ghost cell layer
   // Convert face-centered values to face averages.
   fourthOrderAverage(a_velocity);

   a_geom.fillTransverseGhosts(a_velocity);
}




void
Advect::mappedGridDivergenceWithoutNT( const MagCoordSys* coord_sys,
                                       LevelData<FArrayBox>& a_divF,
                                       const LevelData<FluxBox>& a_F ) const
{
   // Create a temporary to store the normal component of the input
   // flux on cell faces.  We only need enough ghost cells to compute the divergence,
   // even though the input flux may have more.
   CH_assert(a_F.ghostVect() >= a_divF.ghostVect());
   CH_assert(a_F.nComp() == SpaceDim);
   LevelData<FluxBox> FluxNormal(a_F.getBoxes(),
                                 a_F.nComp()/SpaceDim,
                                 a_divF.ghostVect());

   const DisjointBoxLayout& grids = a_divF.disjointBoxLayout();

   // Now compute the divergence in the usual way
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = coord_sys->whichBlock(grids[dit]);
      const MagBlockCoordSys* block_coord_sys( coord_sys->getCoordSys(block_number) );
      RealVect dx = block_coord_sys->dx();

      FArrayBox& thisDiv = a_divF[dit];
      const FluxBox& thisFlux = a_F[dit];
      FluxBox& thisFluxNormal = FluxNormal[dit];
      // First, set divF to 0
      thisDiv.setVal(0.0);
      // since we're computing volume*div, and the fluxes
      // are multiplied by the appropriate areas, don't need
      // a dx here.
      Real fakeDx = 1.0;

      // now loop over directions and increment with directional
      // derivative
      for (int dir=0; dir<SpaceDim; dir++) {
         const Box& thisBox = thisDiv.box();
         FArrayBox& thisFluxNormal_dir = thisFluxNormal[dir];

         double face_area = 1.;
         for (int tdir=0; tdir<SpaceDim; ++tdir) {
            if (tdir != dir) face_area *= dx[tdir];
         }

         thisFluxNormal_dir.copy(thisFlux[dir],dir,0,1);
         thisFluxNormal_dir *= face_area;

         // use fortran from divergence operator here
         FORT_CFGDIVERGENCE(CHF_CONST_FRA(thisFluxNormal_dir),
                             CHF_FRA(thisDiv),
                             CHF_BOX(thisBox),
                             CHF_CONST_REAL(fakeDx),
                             CHF_INT(dir));
      }
   }
}



void
Advect::plotData( const string               a_file_name,
                  const MagGeom& a_geom,
                  const LevelData<FArrayBox>&  a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   Box domain_box = grids.physDomain().domainBox();

   MultiBlockCoordSys* mag_coord_sys = (MultiBlockCoordSys*)a_geom.coordSysPtr();
   WriteMappedUGHDF5(a_file_name, grids, a_data, *mag_coord_sys, domain_box);

}

#include "NamespaceFooter.H"
