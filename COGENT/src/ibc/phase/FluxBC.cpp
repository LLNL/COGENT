#include "FluxBC.H"
#include "Directions.H"
#include "ConstFact.H"
#include "FluxBCF_F.H"
#include "inspect.H"

#include "NamespaceHeader.H"

const string FluxBC::pp_name = "flux_bc";

FluxBC::FluxBC( const BoundaryBoxLayoutPtrVect&  a_bdry_layouts,
                const ParmParse&                 a_pp)

: m_bdry_layouts(a_bdry_layouts),
  m_total_fluxes(3, 0.)

{

   parseParameters( a_pp );
   
   normalizeFluxes( m_total_fluxes );
}

void FluxBC::applyBC(KineticSpecies& a_species,
                     const LevelData<FluxBox>& a_velocity_mapped,
                     const CFG::LevelData<CFG::FArrayBox>& a_phi,
                     const Real& a_time) const
{
   /*
    Fill ghost cell data to enforce given fluxes particle,
    parallel momentum, and energy
    */

   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   
   LevelData<FArrayBox>& soln( a_species.distributionFunction() );
   
   if (!m_BstarParallel.isDefined()) {
      m_BstarParallel.define(soln.getBoxes(),1,IntVect::Unit);
   }
   geometry.getBStarParallel(m_BstarParallel);
   
   LevelData<FArrayBox> phi_injected;
   geometry.injectConfigurationToPhase(a_phi, phi_injected);
   
   // Get (n, Vpar, T) values and (Gamma, M, Q) total fluxes
   // at the physical boundary
   Vector<Real> outflow_moments(NUM_MOMS, 0.);
   
   // Loop over block boundaries included in the physical boundary
   for (int i_bnd=0; i_bnd<m_bdry_layouts.size(); i_bnd++) {
      
      //Get (n, nVpar, nT) moments and (Gamma, M, Q) fluxes
      Vector<Real> this_moments(NUM_MOMS, 0.);
      
      computeOutflowMoments(this_moments,
                            soln,
                            a_velocity_mapped,
                            phi_injected,
                            a_species,
                            m_bdry_layouts[i_bnd]);
      
      for (int n=0; n<NUM_MOMS; n++) {
         outflow_moments[n] += this_moments[n];
      }
   }

   // Obtain averages of dens, vpar, and temperature
   for (int n=0; n<PART_FLUX; n++) {
      outflow_moments[n] /= outflow_moments[AREA];
      if (n != DENSITY) {
         outflow_moments[n] /= outflow_moments[DENSITY];
      }
   }
   
   // Compute the flux (3x3) matrix and the background Maxwellian contribution.
   // Matrix elements correspond to Gamma, M, Q fluxes driven by delta_(n,vpar,T) perturbations
   Vector<Real> flux_matrix(9, 0.);
   Vector<Real> inlfow_bkgr_fluxes(3, 0.);
   
   for (int i_bnd=0; i_bnd<m_bdry_layouts.size(); i_bnd++) {
      
      Vector<Real> this_flux_matrix(9, 0.);
      Vector<Real> this_inlfow_bkgr_fluxes(3, 0.);
      
      computeInflowFluxMatrix(this_flux_matrix,
                              this_inlfow_bkgr_fluxes,
                              outflow_moments,
                              a_velocity_mapped,
                              phi_injected,
                              a_species,
                              m_bdry_layouts[i_bnd]);
      
      for (int n=0; n<9; n++) {
         flux_matrix[n] += this_flux_matrix[n];
      }

      for (int n=0; n<3; n++) {
         inlfow_bkgr_fluxes[n] += this_inlfow_bkgr_fluxes[n];
      }
   }

   // Compute delta_(n,vpar,T) perturbations
   Vector<Real> pert_moments(3, 0.);
   Vector<Real> rhs(3, 0.);
   for (int n=0; n<3; n++) {
      rhs[n] = m_total_fluxes[n] - inlfow_bkgr_fluxes[n] + outflow_moments[n+3];
   }
   
   solveTrilinearSystem(pert_moments, flux_matrix, rhs);
  
   // Fill ghost data
   fillGhosts(a_species, a_velocity_mapped, outflow_moments, pert_moments);
}

void
FluxBC::fillGhosts(KineticSpecies&           a_species,
                   const LevelData<FluxBox>& a_mapped_velocity,
                   const Vector<Real>&       a_outflow_moments,
                   const Vector<Real>&       a_pert_moments) const
{
   /*
    Fill ghosts. For outflow -- extrapolate dfn. For inflow, fill ghosts
    to enforce the perturbed Maxwellian at the boundary face.
    */
   
   const PhaseGeom& geometry = a_species.phaseSpaceGeometry();
   LevelData<FArrayBox>& soln( a_species.distributionFunction() );
   const DisjointBoxLayout& interior_grids = soln.getBoxes();
   const Real mass = a_species.mass();
   IntVect ghosts = soln.ghostVect();
   
   const LevelData<FArrayBox>& injected_B = geometry.getBFieldMagnitude();
   
   // Loop over block boundaries comprising in the physical boundary
    for (int i_bnd=0; i_bnd<m_bdry_layouts.size(); i_bnd++) {
       
       const int& dir( m_bdry_layouts[i_bnd]->dir() );
       const Side::LoHiSide& side( m_bdry_layouts[i_bnd]->side() );
       const int ISIDE(side);
       
       // Loop over (boundary) boxes
       const DisjointBoxLayout& bdry_grids( m_bdry_layouts[i_bnd]->disjointBoxLayout() );
       for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

          CH_assert(bdry_grids[dit].size(dir) == 1);
          CH_assert(ghosts[dir] > 1);
          
          Box boundaryBox(bdry_grids[dit]);
          boundaryBox.growDir(dir, side, ghosts[dir] - 1 );
          
          const DataIndex& interior_dit( m_bdry_layouts[i_bnd]->dataIndex(dit) );
          
          const PhaseBlockCoordSys& phase_coords = geometry.getBlockCoordSys(interior_grids[interior_dit]);

          FArrayBox velocityRealCoords(boundaryBox, VEL_DIM);
          phase_coords.getVelocityRealCoords(velocityRealCoords);
          
          const FArrayBox& this_face_vel = a_mapped_velocity[interior_dit][dir];
          FORT_SET_FLUX_INFLOW_BC(CHF_FRA1(soln[interior_dit],0),
                                  CHF_BOX(boundaryBox),
                                  CHF_CONST_VR(a_outflow_moments),
                                  CHF_CONST_VR(a_pert_moments),
                                  CHF_CONST_FRA(velocityRealCoords),
                                  CHF_CONST_REAL(mass),
                                  CHF_CONST_FRA1(this_face_vel,dir),
                                  CHF_CONST_FRA1(injected_B[interior_dit],0),
                                  CHF_CONST_FRA1(m_BstarParallel[interior_dit],0),
                                  CHF_CONST_INT(dir),
                                  CHF_CONST_INT(ISIDE));
       
       }
    }
}


void
FluxBC::computeOutflowMoments(Vector<Real>&                 a_moments,
                              const LevelData<FArrayBox>&   a_dfn,
                              const LevelData<FluxBox>&     a_mapped_velocity,
                              const LevelData<FArrayBox>&   a_phi_inj,
                              const KineticSpecies&         a_species,
                              const BoundaryBoxLayoutPtr&   a_bdry_layout) const
{
   /*
    Compute FS-averaged (but arbitrary weighted) density, parallel momentum (nVpar),
    and  pressure (nT) of the last valid cell distribution (a_dfn = BstarPar * f)
    and outflow fluxes of the particle, parallel momentum, and energy.
    */
   
   // Get geometry parameters
   const PhaseGeom& geometry = a_species.phaseSpaceGeometry();
   //const DisjointBoxLayout& interior_grids = geometry.gridsFull();
   const DisjointBoxLayout& interior_grids = a_dfn.getBoxes();
   const VEL::VelCoordSys& vel_coords = geometry.velSpaceCoordSys();
   const Real detJ = vel_coords.pointwiseJ( VEL::RealVect::Zero );
   const VEL::RealVect& vel_dx( vel_coords.dx() );
   const Real vel_area = detJ * vel_dx.product();
   const Real mass = a_species.mass();
   const Real charge = a_species.charge();
   
   const LevelData<FArrayBox>& injected_B = geometry.getBFieldMagnitude();
   
   const int& dir( a_bdry_layout->dir() );
   const Side::LoHiSide& side( a_bdry_layout->side() );
   const int ISIDE(side);
   
   // Create storage for local (this rank) and global (_all) data
   double *moments_loc = new double[a_moments.size()];
   double *moments_all = new double[a_moments.size()];
   
   for (int i=0; i<a_moments.size(); ++i) {
       moments_loc[i]=0.0;
       moments_all[i]=0.0;
   }
   
   // Loop over (boundary) boxes
   const DisjointBoxLayout& bdry_grids( a_bdry_layout->disjointBoxLayout() );
   for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
      const Box& boundaryBox = bdry_grids[dit];
      
      const DataIndex& interior_dit( a_bdry_layout->dataIndex(dit) );
      const PhaseBlockCoordSys& phase_coords = geometry.getBlockCoordSys(interior_grids[interior_dit]);
      const CFG::MagBlockCoordSys& cfg_coords = geometry.getMagBlockCoordSys(interior_grids[interior_dit]);
      const CFG::RealVect& cfg_dx = cfg_coords.dx();
      
      CH_assert(boundaryBox.size(dir) == 1);
      Real cfg_area_mapped(1.0);
      for (int tdir = 0; tdir<CFG_DIM; tdir++) {
         if (tdir != dir) {
            cfg_area_mapped *= cfg_dx[tdir];
         }
      }
      
      FArrayBox velocityRealCoords(boundaryBox, VEL_DIM);
      phase_coords.getVelocityRealCoords(velocityRealCoords);
      
      const FluxBox& this_mapped_vel = a_mapped_velocity[interior_dit];
            
      FArrayBox kernels( boundaryBox, a_moments.size() );
      kernels.setVal(0.);
      FORT_COMPUTE_OUTFLOW_MOMENTS(CHF_FRA(kernels),
                                   CHF_BOX(boundaryBox),
                                   CHF_CONST_FRA1(a_dfn[interior_dit],0),
                                   CHF_CONST_FRA1(this_mapped_vel[dir],dir),
                                   CHF_CONST_FRA1(a_phi_inj[interior_dit],0),
                                   CHF_CONST_FRA1(injected_B[interior_dit],0),
                                   CHF_CONST_FRA1(m_BstarParallel[interior_dit],0),
                                   CHF_CONST_FRA(velocityRealCoords),
                                   CHF_CONST_REAL(mass),
                                   CHF_CONST_REAL(charge),
                                   CHF_CONST_REAL(vel_area),
                                   CHF_CONST_REAL(cfg_area_mapped),
                                   CHF_CONST_INT(dir),
                                   CHF_CONST_INT(ISIDE));
      
      // Sum over patches (boxes) that belong to the same processor (rank)
      for (int n = 0; n<a_moments.size(); n++) {
         moments_loc[n] += kernels.sum(n);
      }
   }

   // Sum over all ranks.
   // We could also use mpi_comm local to this boundary date (if needed for speed-up)
   int size = a_moments.size();
   MPI_Allreduce(moments_loc, moments_all, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   for (int i=DENSITY; i<NUM_MOMS; ++i) {
      a_moments[i] = moments_all[i];
   }
   
   delete [] moments_loc;
   delete [] moments_all;
}

void
FluxBC::computeInflowFluxMatrix(Vector<Real>&               a_M,
                                Vector<Real>&               a_background_fluxes,
                                const Vector<Real>&         a_outflow_moms,
                                const LevelData<FluxBox>&   a_mapped_velocity,
                                const LevelData<FArrayBox>& a_phi_inj,
                                const KineticSpecies&       a_species,
                                const BoundaryBoxLayoutPtr& a_bdry_layout) const
{
   /*
    Computes the (3x3) flux matrix. That is Gamma, M, Q fluxes
    driven by delta_(n,vpar,T) perturbations
    */
   
   // Get geometry parameters
   const PhaseGeom& geometry = a_species.phaseSpaceGeometry();
   const DisjointBoxLayout& interior_grids = geometry.gridsFull();
   const VEL::VelCoordSys& vel_coords = geometry.velSpaceCoordSys();
   const Real detJ = vel_coords.pointwiseJ( VEL::RealVect::Zero );
   const VEL::RealVect& vel_dx( vel_coords.dx() );
   const Real vel_area = detJ * vel_dx.product();
   const Real mass = a_species.mass();
   const Real charge = a_species.charge();
   
   const LevelData<FArrayBox>& injected_B = geometry.getBFieldMagnitude();
   
   const int& dir( a_bdry_layout->dir() );
   const Side::LoHiSide& side( a_bdry_layout->side() );
   const int ISIDE(side);
   
   // Create storage for local (this rank) and global (_all) data
   int size = a_M.size() + a_background_fluxes.size();
   double *elements_loc = new double[size];
   double *elements_all = new double[size];
   
   for (int i=0; i<size; ++i) {
       elements_loc[i]=0.0;
       elements_all[i]=0.0;
   }
   
   // Loop over (boundary) boxes
   const DisjointBoxLayout& bdry_grids( a_bdry_layout->disjointBoxLayout() );
   for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {
      const Box& boundaryBox = bdry_grids[dit];
      
      const DataIndex& interior_dit( a_bdry_layout->dataIndex(dit) );
      const PhaseBlockCoordSys& phase_coords = geometry.getBlockCoordSys(interior_grids[interior_dit]);
      const CFG::MagBlockCoordSys& cfg_coords = geometry.getMagBlockCoordSys(interior_grids[interior_dit]);
      const CFG::RealVect& cfg_dx = cfg_coords.dx();
      
      CH_assert(boundaryBox.size(dir) == 1);

      // Get mapped face area
      Real cfg_area_mapped(1.0);
      for (int tdir = 0; tdir<CFG_DIM; tdir++) {
         if (tdir != dir) {
            cfg_area_mapped *= cfg_dx[tdir];
         }
      }
      
      FArrayBox velocityRealCoords(boundaryBox, VEL_DIM);
      phase_coords.getVelocityRealCoords(velocityRealCoords);
      
      const FluxBox& this_mapped_vel = a_mapped_velocity[interior_dit];
            
      FArrayBox kernels( boundaryBox, size );
      kernels.setVal(0.0);
      FORT_COMPUTE_FLUX_MATRIX(CHF_FRA(kernels),
                               CHF_BOX(boundaryBox),
                               CHF_CONST_VR(a_outflow_moms),
                               CHF_CONST_FRA1(this_mapped_vel[dir],dir),
                               CHF_CONST_FRA1(a_phi_inj[interior_dit],0),
                               CHF_CONST_FRA1(injected_B[interior_dit],0),
                               CHF_CONST_FRA1(m_BstarParallel[interior_dit],0),
                               CHF_CONST_FRA(velocityRealCoords),
                               CHF_CONST_REAL(mass),
                               CHF_CONST_REAL(charge),
                               CHF_CONST_REAL(vel_area),
                               CHF_CONST_REAL(cfg_area_mapped),
                               CHF_CONST_INT(dir),
                               CHF_CONST_INT(ISIDE));
      
      // Sum over patches (boxes) that belong to the same processor (rank)
      //const Real factor = vel_area / mass;
      for (int n = 0; n<size; n++) {
         elements_loc[n] += kernels.sum(n);
      }
   }

   // Sum over all ranks.
   // We could also use mpi_comm local to this boundary date (if needed for speed-up)
   MPI_Allreduce(elements_loc, elements_all, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   for (int i=0; i<a_M.size(); ++i) {
      a_M[i] = elements_all[i];
   }

   for (int i=0; i<3; ++i) {
      a_background_fluxes[i] = elements_all[a_M.size()+i];
   }
   
   delete [] elements_loc;
   delete [] elements_all;
}

void
FluxBC::solveTrilinearSystem(Vector<Real>&         a_x,
                             const Vector<Real>&   a_M,
                             const Vector<Real>&   a_b) const
{
   /*
    Solves linear linear system
    a0*x + a1*y + a2*z == b0,
    a3*x + a4*y + a5*z == b1,
    a6*x + a7*y + a8*z == b2
    */

   Real det = a_M[2]*a_M[4]*a_M[6] - a_M[1]*a_M[5]*a_M[6]
            - a_M[2]*a_M[3]*a_M[7] + a_M[0]*a_M[5]*a_M[7]
            + a_M[1]*a_M[3]*a_M[8] - a_M[0]*a_M[4]*a_M[8];
   
   a_x[0] = -a_M[5]*a_M[7]*a_b[0] + a_M[4]*a_M[8]*a_b[0]
          + a_M[2]*a_M[7]*a_b[1] - a_M[1]*a_M[8]*a_b[1]
          - a_M[2]*a_M[4]*a_b[2] + a_M[1]*a_M[5]*a_b[2];
   
   a_x[1] = a_M[5]*a_M[6]*a_b[0] - a_M[3]*a_M[8]*a_b[0]
          - a_M[2]*a_M[6]*a_b[1] + a_M[0]*a_M[8]*a_b[1]
          + a_M[2]*a_M[3]*a_b[2] - a_M[0]*a_M[5]*a_b[2];

   a_x[2] = -a_M[4]*a_M[6]*a_b[0] + a_M[3]*a_M[7]*a_b[0]
          + a_M[1]*a_M[6]*a_b[1] - a_M[0]*a_M[7]*a_b[1]
          - a_M[1]*a_M[3]*a_b[2] + a_M[0]*a_M[4]*a_b[2];
   
   for (int n=0; n<3; ++n) {
      a_x[n] /= -det;
   }
}

void
FluxBC::normalizeFluxes (Vector<Real>& a_fluxes) const
{
   /*
    Convert fluxes specified in Watt=Joule/second
    to COGENT units
    */
   
   //Get normalization parameters (units)
   Real temp, length, dens;
   ParmParse ppunits( "units" );
   ppunits.get("temperature",temp);     //[eV]
   ppunits.get("length",length);        //[m]
   ppunits.get("number_density",dens);  //[m^{-3}]
 
   Real charge = Constants::ELEMENTARY_CHARGE;
   Real tempJoules = charge * temp;

   Real masskg = Constants::MASS_OF_PROTON;
   Real velocity = sqrt( tempJoules / masskg );

   //Renormalize fluxes
   for (int i=0; i<a_fluxes.size(); ++i) {
     a_fluxes[i] /= dens * tempJoules * pow(length,2) * velocity;
   }
}

void
FluxBC::parseParameters(const ParmParse& a_pp)
{
   a_pp.query("particle_flux", m_total_fluxes[0]);
   a_pp.query("parallel_momentum_flux", m_total_fluxes[1]);
   a_pp.query("heat_flux", m_total_fluxes[2]);
}

#include "NamespaceFooter.H"

