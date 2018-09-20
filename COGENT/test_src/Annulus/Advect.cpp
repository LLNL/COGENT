#include "Advect.H"
#include "AdvectF_F.H"

#include "FourthOrderUtil.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "CornerCopier.H"

#include "MagCoordSys.H"
#include "MagBlockCoordSys.H"
#include "altFaceAverages.H"
#include "mappedAdvectionFlux.H"
#include "newMappedGridIO.H"

#include "InOutFlowBC.H"

#include "NamespaceHeader.H"

using std::max;

char* Advect::pp_name = {"Advect"};

Advect::Advect( ParmParse pp)
  : m_limit_face_values(false),
    m_saved_dt(-1.0)
{
   if (pp.contains("verbose") && procID()==0) {
      pp.get("verbose", m_verbose);
   }
   else {
      m_verbose = false;
   }
}



Advect::~Advect()
{
}



inline
void fixPeriodicCorners( LevelData<FArrayBox>& a_soln,
                         const MultiBlockCoordSys& a_coord_sys )
{
   // to handle copying of periodic corners, since these are not "valid"
   // cells that are exchanged using a normal exchange() call
   for (int block(0); block<a_coord_sys.numBlocks(); ++block) {
      const MagBlockCoordSys*
         block_coord_sys(
            static_cast<const MagBlockCoordSys*>(
               a_coord_sys.getCoordSys( block ) ) );
      const ProblemDomain& domain( block_coord_sys->domain() );
      if (domain.isPeriodic()) {
         CornerCopier ccopier( a_soln.disjointBoxLayout(),
                               a_soln.disjointBoxLayout(),
                               domain,
                               a_soln.ghostVect() );
         a_soln.exchange( ccopier );
      }
   }
}


void Advect::evalRHS( const MultiBlockLevelGeom&  geom,
                      LevelData<FArrayBox>&       rhs,
                      const LevelData<FArrayBox>& soln,
                      const Real                  time )
{
   /*
     Evaluates the (negated) divergence:

     rhs = - divergence ( velocity * soln )

   */

   if (m_verbose) cout << "   Accumulating RHS" << endl;

   Real localDt = -1.0;

   const MultiBlockCoordSys* coords = geom.coordSysPtr();

   /*
     Compute the advection velocity .
   */
   if (m_verbose) cout << "   Computing velocities" << endl;

   // Compute the advection velocity.
   const DisjointBoxLayout & dbl = soln.getBoxes();
   LevelData<FluxBox> velocity(dbl, SpaceDim, IntVect::Unit);
   updateVelocity(geom, velocity);

   if (m_verbose) cout << "   Done computing velocities" << endl;

   // since we've got the advection velocities here, this is a
   // good place to compute dt
   // use a cfl of 1 here, and we can apply the real cfl elsewhere.
   Real fakeCFL = 1.0;
   Real speciesDt = computeMappedDtSpecies(velocity, geom, fakeCFL);

   CH_assert(speciesDt > 0);

   if (localDt < 0)
      {
         localDt = speciesDt;
      }
   else if (speciesDt < localDt)
      {
         localDt = speciesDt;
      }

   if (m_verbose) cout << "   Applying BCs" << endl;

   // Copy const soln to a temporary so we can set ghost cells
   LevelData<FArrayBox> tmp;
   tmp.define( soln );
   InOutFlowBC iobc;
   iobc.setFaceVelocity( &velocity );
   iobc.ghostCellBC( tmp, geom, 0, time );
   iobc.clearFaceVelocity();
   fixPeriodicCorners( tmp, *coords );
   tmp.exchange();

   if (m_verbose) cout << "   Computing fluxes" << endl;

   LevelData<FluxBox> flux( dbl, SpaceDim, IntVect::Zero );
   computeFlux( tmp, velocity, flux, geom );

   // set boundary conditions on flux, if appropriate
   //bc.fluxBC(flux, velocity, dfn, geom, 0, time);

   if (m_verbose) cout << "   Computing divergence" << endl;

   // Compute the flux divergence
   mappedGridDivergenceWithoutNT(coords, rhs, flux);

   // Divide by cell volume and move the average flux divergence to the right-hand side

   DataIterator it = dbl.dataIterator();
   for (it.begin(); it.ok(); ++it) {
      int block_number = coords->whichBlock(dbl[it]);
      const MagBlockCoordSys* block_coord_sys
         = static_cast<const MagBlockCoordSys*>(coords->getCoordSys(block_number));

      double fac = -1. / block_coord_sys->getMappedCellVolume();

      rhs[it].mult(fac);
   }

   // set dt_save
   m_saved_dt = localDt;

   if (m_verbose) cout << "   Done accumulating RHS" << endl;
}


void Advect::computeFlux( const LevelData<FArrayBox>& dist_fn,
                          const LevelData<FluxBox>&   velocity,
                          LevelData<FluxBox>&         flux,
                          const MultiBlockLevelGeom&  a_geom)
{
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
   //   a_geom.divideJ(dist_cg);  // already did this before evalRHS

   MultiBlockCoordSys* coords = a_geom.coordSysPtr();

   // construct appropriately accurate face-averages of phi and advVel
   LevelData<FluxBox> faceDist(grids, dist_cg.nComp(), dist_cg.ghostVect() );
#if 0
   computeFaceAverages( faceDist, dist_cg );
   faceDist.exchange();

   // if we're limiting the face-centered values, do it here
   if (m_limit_face_values)
      {
         applyMappedLimiter(faceDist, dist_cg, velocity, *coords);
      }
#else

   // The alternative limiter needs another layer of velocity ghost cells.
   // Can we avoid this copy somehow?

   LevelData<FluxBox> grown_velocity(velocity.disjointBoxLayout(),
                                     velocity.nComp(),
                                     velocity.ghostVect() + IntVect::Unit);

   for (dit.begin(); dit.ok(); ++dit) {
      for ( int dir=0; dir<SpaceDim; ++dir ) {
         grown_velocity[dit][dir].setVal(0.);
         grown_velocity[dit][dir].copy(velocity[dit][dir]);
      }
   }
   grown_velocity.exchange();

   altFaceAverages( faceDist, dist_fn, grown_velocity, *coords );

#endif

   // compute computational-space fluxes -
   computeCompFaceFluxes( flux, faceDist, velocity, true );

   flux.exchange();
}



Real Advect::computeDt( const MultiBlockLevelGeom& geom,
                        LevelData<FArrayBox>&      soln )
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

      LevelData<FluxBox> velocity(soln.getBoxes(), SpaceDim, IntVect::Unit);
      updateVelocity(geom, velocity);

      // since we've got the advection velocities here, this is a
      // good place to compute dt
      // use a cfl of 1 here, and we can apply the real cfl elsewhere.
      Real fakeCFL = 1.0;
      Real speciesDt = computeMappedDtSpecies(velocity,
                                              geom,
                                              fakeCFL);

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
  bool useSecondOrder = false;

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



Real
Advect::computeMappedDtSpecies(const LevelData<FluxBox>&  a_faceVel,
                               const MultiBlockLevelGeom& a_geom,
                               Real                       a_cfl)
{
  const DisjointBoxLayout grids = a_faceVel.getBoxes();
  CH_assert(grids == a_geom.grids());

  DataIterator dit = grids.dataIterator();

  // compute area-weighted normal velocities -- 2nd order should
  // be good enough
  bool doFourthOrder = false;
  LevelData<FluxBox> faceNormalVel(grids, 1);

  //  const MultiBlockCoordSys* coords = a_geom.coordSysPtr();
  MagCoordSys* coords = static_cast<MagCoordSys*>(a_geom.coordSysPtr());

  for (dit.begin(); dit.ok(); ++dit) {
     int block_number = coords->whichBlock(grids[dit]);
     const MagBlockCoordSys* block_coord_sys
        = static_cast<const MagBlockCoordSys*>(coords->getCoordSys(block_number));
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
  coords->getCellVolumes(cellVolumes);

  Real maxVelLoc = 0;
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox cellVel(grids[dit], 1);
      for (int dir=0; dir<SpaceDim; dir++)
        {
          int faceComp = 0;
          int cellComp = 0;
          // take max face velocity on the two faces
          // use absolute value
          faceNormalVel[dit][dir].abs();
          EdgeToCellMax(faceNormalVel[dit], faceComp,
                        cellVel, cellComp,
                        dir);
          cellVel.divide(cellVolumes[dit],
                         0,0, 1);

          // now compute maxVelLoc on this box
          // note that this is essentially an estimate of max(vel/dx)
          Real thisMax = cellVel.norm(0,0,1);
          if (thisMax > maxVelLoc) maxVelLoc = thisMax;
        }

    }

  Real maxVel = maxVelLoc;
#ifdef CH_MPI
  MPI_Allreduce(&maxVelLoc, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  Real smallVal = 1.0e-15;
  Real newDt = 0;

  if (maxVel > smallVal)
    {
      newDt = a_cfl/maxVel;
    }

  return newDt;
}



void
Advect::applyMappedLimiter( LevelData<FluxBox>&         a_facePhi,
                            const LevelData<FArrayBox>& a_cellPhi,
                            const LevelData<FluxBox>&   a_faceVel,
                            const MultiBlockCoordSys&   a_coords)
{
#if 0
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
  DataIterator dit = grids.dataIterator();

  // in order to do upwinding, need normal velocities in
  // computational space. computing Fourth-order face
  // averages just in order to upwinding is over-working the
  // issue, but it's the most convenient thing to do here.
  LevelData<FluxBox> normalVel(grids, 1, a_faceVel.ghostVect());

  for (dit.begin(); dit.ok(); ++dit) {
     int block_number = a_coords.whichBlock(grids[dit]);
     const PhaseBlockCoordSys* block_coord_sys
        = static_cast<const PhaseBlockCoordSys*>(a_coords.getCoordSys(block_number));
     RealVect face_area = block_coord_sys->getMappedFaceArea();

     const FluxBox& thisFaceVel = a_faceVel[dit];
     FluxBox& thisNormalVel = normalVel[dit];
     for (int dir=0; dir<SpaceDim; ++dir) {
        FArrayBox& thisNormalVel_dir = thisNormalVel[dir];
        thisNormalVel_dir.copy(thisFaceVel[dir],dir,0,1);
        thisNormalVel_dir *= face_area[dir];
     }
  }

  for (dit.begin(); dit.ok(); ++dit) {
     int block_number = a_coords.whichBlock(grids[dit]);
     const PhaseBlockCoordSys* block_coord_sys
        = static_cast<const PhaseBlockCoordSys*>(a_coords.getCoordSys(block_number));
     const RealVect& dx = block_coord_sys->dx();

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
          FORT_SELECTUPWIND(CHF_FRA(thisFacePhiDir),
                            CHF_CONST_FRA(leftPhi),
                            CHF_CONST_FRA(rightPhi),
                            CHF_CONST_FRA1(thisNormalVelDir,0),
                            CHF_BOX(faceBox));

        } // end loop over directions

    } // end loop over grids
#endif
}



void
Advect::updateVelocity(const MultiBlockLevelGeom&  geom,
                       LevelData<FluxBox>&         velocity)
{
   MagCoordSys* coords = static_cast<MagCoordSys*>(geom.coordSysPtr());
   const DisjointBoxLayout& grids = velocity.disjointBoxLayout();

   LevelData<FluxBox> N_face(velocity.disjointBoxLayout(), SpaceDim*SpaceDim, velocity.ghostVect());
   coords->getPointwiseN(N_face);

   DataIterator dit = velocity.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = coords->whichBlock(grids[dit]);
      const MagBlockCoordSys* block_coord_sys
         = static_cast<const MagBlockCoordSys*>(coords->getCoordSys(block_number));

      const RealVect& dx = block_coord_sys->dx();

      FluxBox& this_velocity = velocity[dit];

      for (int dir=0; dir<SpaceDim; ++dir) {
         FORT_VEL_TEST(CHF_CONST_INT(dir),
                       CHF_BOX(this_velocity[dir].box()),
                       CHF_CONST_REALVECT(dx),
                       CHF_FRA(this_velocity[dir]));
      }

      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_velocity_dir = this_velocity[dir];
         Box vel_box_dir(grids[dit]);
         vel_box_dir.surroundingNodes(dir);

         // Multiply the velocity by NTranspose to convert the velocities
         // components to the physical frame
         FORT_MULT_CFG_NT(
                          CHF_BOX(vel_box_dir),
                          CHF_CONST_FRA(N_face[dit][dir]),
                          CHF_FRA(this_velocity_dir)
                          );

         secondOrderTransExtrap(this_velocity_dir, dir, vel_box_dir);
      }
   }

   //   plotData("velocity", geom, velocity);
   //   exit(1);
}


void
Advect::mappedGridDivergenceWithoutNT( const MultiBlockCoordSys* coord_sys,
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
      const MagBlockCoordSys* block_coord_sys
         = static_cast<const MagBlockCoordSys*>(coord_sys->getCoordSys(block_number));
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
                  const MultiBlockLevelGeom& a_geom,
                  const LevelData<FluxBox>&  a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   LevelData<FArrayBox> data_cell(grids, a_data.nComp(), a_data.ghostVect());

   DataIterator dit2 = grids.dataIterator();
   for (dit2.begin(); dit2.ok(); ++dit2) {
      const FluxBox& this_data = a_data[dit2];
      FArrayBox& this_data_cell = data_cell[dit2];
      this_data_cell.setVal(0.);
      for (int dir=0; dir<SpaceDim; dir++) {
         const FArrayBox& this_data_dir = this_data[dir];
         FArrayBox shift_fab(this_data_dir.box(), this_data_dir.nComp());
         shift_fab.copy(this_data_dir);

         shift_fab.shiftHalf(dir,-1);
         this_data_cell.plus(shift_fab);
         shift_fab.shiftHalf(dir,2);
         this_data_cell.plus(shift_fab);
      }
      this_data_cell.mult(0.5/SpaceDim);
   }

   // bit of a hack to get a reasonable domainBox
   // this is silly, but quick and easy
   IntVectSet gridIVS;
   LayoutIterator lit = grids.layoutIterator();
   for (lit.begin(); lit.ok(); ++lit) {
      gridIVS |= grids[lit];
   }
   Box domain_box = gridIVS.minBox();

   const MagCoordSys* mag_coord_sys = static_cast<const MagCoordSys*>(a_geom.coordSysPtr());
   WriteMappedUGHDF5(a_file_name, grids, data_cell, *mag_coord_sys, domain_box);
}


#include "NamespaceFooter.H"
