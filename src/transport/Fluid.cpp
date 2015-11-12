#include <math.h>

#include "FourthOrderUtil.H"
#include "CONSTANTS.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "CartesianCoordSys.H"
#include "MagCoordSys.H"
#include "MagBlockCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "PhaseGeom.H"
#include "PhaseGeomF_F.H"
#include "PhaseBlockCoordSys.H"

#include "Fluid.H"
#include "TransportF_F.H"

#include "NamespaceHeader.H" //Should be the last one

Fluid::Fluid( const string& a_species_name, ParmParse& a_pptpm, const int a_verbosity )
    : D_m(0),
      m_verbosity(0),
      m_model_only(false),
      m_first_step(true)
{
   m_verbosity = a_verbosity;
   m_pptpm = a_pptpm;
   m_species_name = a_species_name;

   if (m_pptpm.contains("verbose") && procID()==0) {
     m_pptpm.get("verbose", m_verbosity);
   }
   if (m_pptpm.contains("D_matrix")){
     m_pptpm.getarr("D_matrix",D_m,0,4);
   }
   if (m_pptpm.contains("model_only")) {
     m_pptpm.get("model_only",m_model_only);
   }

}

Fluid::~Fluid()
{
}

void Fluid::evalTpmRHS( KineticSpeciesPtrVect&       rhs,
                        const KineticSpeciesPtrVect& soln,
                        const int                    species,
                        const Real                   time )
{

  /*
    Evaluates the fluid radial flux across magnetic surfaces: df/dt_anomalous=-div(flux_transport).
    Currently, this computes second order diffusion of density with constant Diffusion coefficient
    [i.e. partial n/ partial t = -div( -D*grad_r(n)) ] by computing the following eq. for fB:
    df/dt = -div( -D*grad_r(f) ) where all spatial derivatives are taken at constant velocity.
    [i.e. df(r, mu(r,v))/dr|v=const = df/dr|mu=const + (dmu/dr|v=const)*df/dmu|r=const ]
  */
      // parse the pptpm database for initial condition items for this species
      ParseParameters();

      // print parameters at the first time step
      if ((m_verbosity) && (m_first_step)) {printParameters(soln);}

      // get vlasov RHS for the current species
      KineticSpecies& rhs_species( *(rhs[species]) );
      LevelData<FArrayBox> & rhs_dfn = rhs_species.distributionFunction();

      // get solution distribution function (f*Bstarpar) for the current species
      const KineticSpecies& soln_species( *(soln[species]) );
      const LevelData<FArrayBox> & soln_fB = soln_species.distributionFunction();

      // get coordinate system parameters
      const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
      const MultiBlockCoordSys* coords = phase_geom.coordSysPtr();
      const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
      const ProblemDomain& phase_domain = phase_geom.domain();
      const Box& domain_box = phase_domain.domainBox();
      int num_r_cells = domain_box.size(0);
      int num_mu_cells = domain_box.size(3);

      // copy const soln_fB to a temporary 
      const IntVect ghostVect(2*IntVect::Unit);
      const DisjointBoxLayout& dbl = soln_fB.getBoxes();
      LevelData<FArrayBox> fB(dbl, 1, ghostVect); 
      DataIterator sdit = fB.dataIterator();
      for (sdit.begin(); sdit.ok(); ++sdit) {
        fB[sdit].copy(soln_fB[sdit]);
      }
      
      fB.exchange();

      // create cell-centered dfB/dr and dfB/dmu
      LevelData<FArrayBox> dfB_dr_cc(dbl, 1, IntVect::Unit);
      LevelData<FArrayBox> dfB_dmu_cc(dbl, 1, IntVect::Unit);
      DataIterator dit0 = dfB_dmu_cc.dataIterator();
      for (dit0.begin(); dit0.ok(); ++dit0)
      {
        // get phase space dx on this patch
        int block_number = coords->whichBlock(dbl[dit0]);
        const PhaseBlockCoordSys* block_coord_sys = static_cast<const PhaseBlockCoordSys*>(coords->getCoordSys(block_number));
        const RealVect& coord_dx =  block_coord_sys->dx();

        const FArrayBox& fB_on_patch = fB[dit0];
        FArrayBox& dfB_dmu_on_patch = dfB_dmu_cc[dit0];
        FArrayBox& dfB_dr_on_patch  = dfB_dr_cc[dit0];

        FORT_DFB_CELL_CENTER(CHF_BOX(dfB_dmu_on_patch.box()),
                             CHF_CONST_REALVECT(coord_dx),
                             CHF_CONST_INT(num_r_cells),
                             CHF_CONST_INT(num_mu_cells),
                             CHF_CONST_FRA1(fB_on_patch,0),
                             CHF_FRA1(dfB_dmu_on_patch,0),
                             CHF_FRA1(dfB_dr_on_patch,0));
      }
      dfB_dr_cc.exchange();
      dfB_dmu_cc.exchange();

      // create transport flux box and iterate over patches
      LevelData<FluxBox> fluxA(dbl, 1, IntVect::Zero);
      DataIterator dit = fluxA.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
        // get phase space dx on this patch
        const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[dit]);
        const RealVect& phase_dx =  block_coord_sys.dx();

        // get face/cell centered metrics h_r, h_theta, and h_phi on this patch
        Box phase_box(dbl[dit]);
        FluxBox metrics_faces(phase_box, 3);   // to put 3 real space metric componenents on each CFG DIM face
        FArrayBox metrics_cells(phase_box, 3); // to put 3 real space metric componenents in each CFG DIM cell
        getFaceCenteredMetrics(metrics_faces, phase_geom, dbl[dit]);
        getCellCenteredMetrics(metrics_cells, phase_geom, dbl[dit]);

        // create face-centered fluxes on this patch
        const FArrayBox& fB_on_patch = fB[dit];
        const FArrayBox& b_on_patch = inj_B[dit];
        const FArrayBox& dfB_dmu_on_patch = dfB_dmu_cc[dit];
        const FArrayBox& dfB_dr_on_patch  = dfB_dr_cc[dit];

        for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisFluxA = fluxA[dit][dir];

          FORT_EVAL_FLUID_FLUX(CHF_CONST_INT(dir),
                               CHF_BOX(thisFluxA.box()),
                               CHF_CONST_REALVECT(phase_dx),
                               CHF_FRA1(thisFluxA,0),
                               CHF_CONST_FRA(metrics_faces[dir]),
                               CHF_CONST_FRA(metrics_cells),
                               CHF_CONST_FRA1(fB_on_patch,0),
                               CHF_CONST_FRA1(b_on_patch,0),
                               CHF_CONST_FRA1(dfB_dmu_on_patch,0),
                               CHF_CONST_FRA1(dfB_dr_on_patch,0),
                               CHF_CONST_INT(num_r_cells),
                               CHF_CONST_INT(num_mu_cells),
                               CHF_CONST_REAL(D_m[0]));
        }
      }

      // calculate div(flux)
      phase_geom.averageAtBlockBoundaries(fluxA);
      LevelData<FArrayBox> rhs_transport;
      rhs_transport.define(rhs_dfn);
      phase_geom.mappedGridDivergenceFromFluxNormals(rhs_transport, fluxA);
      DataIterator rdit = rhs_transport.dataIterator();

      // add (or overwrite) transport RHS to Vlasov RHS
      for (rdit.begin(); rdit.ok(); ++rdit)
      {
        if(m_model_only){
        rhs_dfn[rdit].copy(rhs_transport[rdit]);}
        else{
        rhs_dfn[rdit].plus(rhs_transport[rdit]);}
      }

  m_first_step = false;  

}

void Fluid::getFaceCenteredMetrics(FluxBox&            metrics_faces,
                                     const PhaseGeom&    a_phase_geom,
                                     const Box&          a_dbl)
{
   // get face centered N components on this patch to formulate face centered metrics h_r and h_theta
   const CFG::MagBlockCoordSys& mag_block_coord_sys = a_phase_geom.getMagBlockCoordSys(a_dbl);
   const CFG::RealVect& real_dx = mag_block_coord_sys.dx();
   Box phase_box(a_dbl);
   CFG::Box box_config;
   a_phase_geom.projectPhaseToConfiguration(phase_box, box_config);
   box_config.grow(2);
   CFG::FluxBox N_cfg_face(box_config, CFG_DIM*CFG_DIM);  //4 total components of N on each CFG DIM face
   mag_block_coord_sys.getPointwiseN(N_cfg_face);
   FluxBox N_face;
   a_phase_geom.injectConfigurationToPhase(N_cfg_face, N_face);

   // get face-centered 2*pi*Rmaj on this patch (note that 2piRmaj=h_toroidalangle)
   CFG::FluxBox TwoPiR_cfg_face(box_config, 1);
   for (int dir=0; dir<CFG_DIM; ++dir)
   {
     CFG::RealVect offset = real_dx;
     offset *= 0.5;
     offset[dir] = 0.0;
     CFG::FArrayBox& this_fab = TwoPiR_cfg_face[dir];
     CFG::BoxIterator bit(this_fab.box());
     for (bit.begin(); bit.ok(); ++bit)
     {
       CFG::IntVect iv = bit();
       CFG::RealVect mapped_loc = iv*real_dx + offset;
       double TwoPiRmaj = 2. * Pi * mag_block_coord_sys.majorRadius(mapped_loc);
       this_fab(iv) = TwoPiRmaj;
     }
   }
   FluxBox TwoPiR_face;
   a_phase_geom.injectConfigurationToPhase(TwoPiR_cfg_face,TwoPiR_face);

   for (int dir=0; dir<CFG_DIM; dir++)
   {
      FArrayBox& metrics_on_patch = metrics_faces[dir];

      FORT_METRICS_FACES( CHF_BOX(metrics_on_patch.box()),
                         CHF_CONST_INT(dir),
                         CHF_CONST_FRA1(TwoPiR_face[dir],0),
                         CHF_CONST_FRA(N_face[dir]),
                         CHF_FRA(metrics_on_patch));
   }

}

void Fluid::getCellCenteredMetrics(FArrayBox&          metrics_cells,
                                     const PhaseGeom&    a_phase_geom,
                                     const Box&          a_dbl)
{
   // get cell centered N components on this patch to formulate cell centered metrics h_r and h_theta
   const CFG::MagBlockCoordSys& mag_block_coord_sys = a_phase_geom.getMagBlockCoordSys(a_dbl);
   const CFG::RealVect& real_dx = mag_block_coord_sys.dx();
   Box phase_box(a_dbl);
   CFG::Box box_config;
   a_phase_geom.projectPhaseToConfiguration(phase_box, box_config);
   box_config.grow(2);
   CFG::FArrayBox N_cfg_cent(box_config, CFG_DIM*CFG_DIM);  //4 total components of N on each cell center
   mag_block_coord_sys.getPointwiseN(N_cfg_cent);
   FArrayBox N_cent;
   a_phase_geom.injectConfigurationToPhase(N_cfg_cent, N_cent);

   // get cell centered 2*pi*Rmaj on this patch (note that 2piRmaj=h_toroidalangle)
   CFG::FArrayBox TwoPiR_cfg_cent(box_config, 1);
   CFG::BoxIterator bit(TwoPiR_cfg_cent.box());
   CFG::RealVect offset = real_dx;
   offset *= 0.5;
   for (bit.begin(); bit.ok(); ++bit)
   {
     CFG::IntVect iv = bit();
     CFG::RealVect mapped_loc = iv*real_dx + offset;
     double TwoPiRmaj = 2. * Pi * mag_block_coord_sys.majorRadius(mapped_loc);
     TwoPiR_cfg_cent(iv) = TwoPiRmaj;
   }
   FArrayBox TwoPiR_cent;
   a_phase_geom.injectConfigurationToPhase(TwoPiR_cfg_cent,TwoPiR_cent);

   FORT_METRICS_CELLS( CHF_BOX(metrics_cells.box()),
                       CHF_CONST_FRA1(TwoPiR_cent,0),
                       CHF_CONST_FRA(N_cent),
                       CHF_FRA(metrics_cells));

}

void Fluid::ParseParameters()
{
   m_pptpm.queryarr("D_matrix",D_m,0,4);
   CH_assert( D_m[0] >= 0);
   m_pptpm.query("model_only", m_model_only);
}

inline void Fluid::printParameters(const KineticSpeciesPtrVect& soln)
{
   // get stuff to calculate stability parameters on time step and number of mu cells
   const KineticSpecies& soln_species( *(soln[0]) );
   const LevelData<FArrayBox>& soln_dfn = soln_species.distributionFunction();
   const DisjointBoxLayout& dbl = soln_dfn.getBoxes();
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
   const ProblemDomain& phase_domain = phase_geom.domain();
   const Box& domain_box = phase_domain.domainBox();
   int num_r_cells = domain_box.size(0);
   int num_mu_cells = domain_box.size(3);

   // get dr and dlnB/dr, and hr
   const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
   const DisjointBoxLayout& grids = inj_B.getBoxes();
   LevelData<FArrayBox> hr(grids, 1, IntVect::Zero); // r metrics
   DataIterator bdit= inj_B.dataIterator();
   LevelData<FArrayBox> dlnB_dr;
   dlnB_dr.define(inj_B);
   Real dr;
   for (bdit.begin(); bdit.ok(); ++bdit)
   {
     const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[bdit]);
     const RealVect& phase_dx =  block_coord_sys.dx();
     dr = phase_dx[0];

     // get face centered metrics h_r, h_theta, and h_phi on this patch
     Box phase_box(dbl[bdit]);
     FArrayBox metrics_cells(phase_box, 3); // to put 3 real space metric componenents in each CFG DIM cell
     getCellCenteredMetrics(metrics_cells, phase_geom, dbl[bdit]);
     hr[bdit].copy(metrics_cells,0,0,1);

     const FArrayBox& b_on_patch = inj_B[bdit];
     FArrayBox& dlnB_dr_on_patch = dlnB_dr[bdit];

     FORT_DLOGB_DR( CHF_BOX(b_on_patch.box()),
                    CHF_CONST_REAL(dr),
                    CHF_CONST_FRA1(b_on_patch,0),
                    CHF_FRA1(dlnB_dr_on_patch,0));
   }
   hr.exchange();
   dlnB_dr.exchange();

   // get the minimum and maximum value of dlnB/dr and min value of hr
   Real local_minimum(100);     // 100 is just a starting point
   Real local_minimum_hr(100);  // 100 is just a starting point
   for (bdit.begin(); bdit.ok(); ++bdit)
   {
      Box box( grids[bdit] );
      Real box_min( dlnB_dr[bdit].min( box ) );
      Real box_min_hr( hr[bdit].min (box) );
      local_minimum = Min( local_minimum, box_min );
      local_minimum_hr = Min( local_minimum_hr, box_min_hr );
   }
   Real min_dlnB_dr( local_minimum );
   Real min_hr( local_minimum_hr );
#ifdef CH_MPI
   MPI_Allreduce( &local_minimum, &min_dlnB_dr, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
   MPI_Allreduce( &local_minimum_hr, &min_hr, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif

// cout << "min =" << min_dlnB_dr << endl;
// cout << "min =" << min_hr << endl;

   const Real dt_stable = 0.5*dr*dr*min_hr*min_hr/D_m[0]/(1-min_dlnB_dr*num_mu_cells*dr)/(1-min_dlnB_dr*num_mu_cells*dr);

   if (procID() == 0 ) {
     cout << "Fluid Transport Model being defined for species " << m_species_name <<endl;
     cout << "Diff Coef.  =  " << D_m[0] << ", model only = " <<  m_model_only << endl;
     cout << "Stable time step < " << dt_stable << endl;
   }

}


#include "NamespaceFooter.H"
