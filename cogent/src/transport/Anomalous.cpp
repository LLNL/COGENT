#include <math.h>
#include "CONSTANTS.H"

#include "Directions.H"
#include "PhaseGeom.H"
#include "PhaseGeomF_F.H"
#include "PhaseBlockCoordSys.H"
#include "Anomalous.H"
#include "TransportF_F.H"

#include "MomentOp.H"
#include "Kernels.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "MagBlockCoordSys.H"
#include "FourthOrderUtil.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H" // has to be the last one

Anomalous::Anomalous( const string& a_species_name, ParmParse& a_pptpm, const int a_verbosity )
    : D_fluid(4,0),
      D_kinet(4,0),
      DN0(0),
      model_only(false),
      const_coeff(true),
      m_simple_diffusion(false),
      moment_op( MomentOp::instance() ),
      m_first_step(true),
      m_arbitrary_grid(true),
      verbosity(0)
{

   verbosity = a_verbosity;
   pptpm = a_pptpm;
   species_name = a_species_name;

}

Anomalous::~Anomalous()
{
}

void Anomalous::evalTpmRHS( KineticSpecies&               rhs_species,
                            const KineticSpeciesPtrVect&  soln,
                            const int                     species,
                            const Real                    time )
{
  /*
    Evaluates the anomalous flux normal to magnetic surfaces as
    described in Ch. 5 of Justin Angus' PhD Thesis:
    df/dt_anomalous=-div(flux_anom)
    flux_anom = flux_psi*psi_hat
    flux_psi = U_psi*f - D_psi*nabla_psi(f)
    nabla_psi =  nabla_psi|v=const = nabla_psi|mu=const - mu/B*nabla_psi(B)*d/dmu
    U_psi = Dn0*nabla_psi(ln(n)) + DT0*nabla_psi(ln(T))
          + (vp^2/Vth^2-3/2)*(Dn2*nabla_psi(ln(n))+DT2*nabla_psi(ln(T)))
    D_psi = D0
  */
      // parse the pptpm database for initial condition items for this species
      if ( m_first_step ) {ParseParameters();}

      // print parameters at the first time step
      const KineticSpecies& soln_species( *(soln[species]) );
      if ((verbosity) && (m_first_step)) {printParameters(soln_species);}

      // get vlasov RHS for the current species
      //KineticSpecies& rhs_species( *(rhs[species]) );
      LevelData<FArrayBox>& rhs_dfn = rhs_species.distributionFunction();

      // get solution distribution function (f*Bstarpar) for the current species
      //const KineticSpecies& soln_species( *(soln[species]) );
      const LevelData<FArrayBox>& soln_fB = soln_species.distributionFunction();
      double mass = soln_species.mass();
      const DisjointBoxLayout& dbl = soln_fB.getBoxes();

      // get coordinate system parameters
      const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
      const CFG::MagGeom& mag_geom = phase_geom.magGeom();
      const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
      const ProblemDomain& phase_domain = phase_geom.domain();
      const Box& domain_box = phase_domain.domainBox();
      int num_r_cells = domain_box.size(0);
      int num_mu_cells = domain_box.size(3);

      // copy const soln_fB (along wiht the boundary values) to a temporary
      const IntVect ghostVect(4*IntVect::Unit);
      LevelData<FArrayBox> fB(dbl, 1, ghostVect);
      DataIterator sdit = fB.dataIterator();
      for (sdit.begin(); sdit.ok(); ++sdit) {
        fB[sdit].copy(soln_fB[sdit]);
      }
      fB.exchange();

      // get density, mean parallel velocity, temperature, and fourth moment coefficient
      const CFG::IntVect cfg_ghostVect(4*CFG::IntVect::Unit);
      CFG::LevelData<CFG::FArrayBox> dens_cfg(mag_geom.grids(), 1, cfg_ghostVect);
      CFG::LevelData<CFG::FArrayBox> Upar_cfg(mag_geom.grids(), 1, cfg_ghostVect);
      CFG::LevelData<CFG::FArrayBox> temp_cfg(mag_geom.grids(), 1, cfg_ghostVect);
      CFG::LevelData<CFG::FArrayBox> four_cfg(mag_geom.grids(), 1, cfg_ghostVect);
      CFG::LevelData<CFG::FArrayBox> perp_cfg(mag_geom.grids(), 1, cfg_ghostVect);
      moment_op.compute(dens_cfg, soln_species, fB, DensityKernel());
      moment_op.compute(Upar_cfg, soln_species, fB, ParallelMomKernel());
      moment_op.compute(four_cfg, soln_species, fB, FourthMomentKernel());
      moment_op.compute(perp_cfg, soln_species, fB, PerpEnergyKernel());
      CFG::DataIterator cfg_dit = dens_cfg.dataIterator();
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
        Upar_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
      }
      moment_op.compute(temp_cfg, soln_species, fB, PressureKernel(Upar_cfg));

      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
        four_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by pressure
        perp_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by pressure
      }
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
        temp_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
      }
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
        four_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by temperature
      }
      LevelData<FArrayBox> density;
      LevelData<FArrayBox> Upar;
      LevelData<FArrayBox> temperature;
      LevelData<FArrayBox> fourth_coef;
      LevelData<FArrayBox> perp_coef;
      phase_geom.injectConfigurationToPhase(dens_cfg, density);
      phase_geom.injectConfigurationToPhase(Upar_cfg, Upar);
      phase_geom.injectConfigurationToPhase(temp_cfg, temperature);
      phase_geom.injectConfigurationToPhase(four_cfg, fourth_coef);
      phase_geom.injectConfigurationToPhase(perp_cfg, perp_coef);

   
      //get NJinverse and B-field data for dealigned grid calculations
      CFG::LevelData<CFG::FluxBox> pointwiseNJinv_cfg(mag_geom.grids(), CFG_DIM*CFG_DIM, cfg_ghostVect);
      mag_geom.getPointwiseNJinverse(pointwiseNJinv_cfg);
      LevelData<FluxBox> inj_pointwiseNJinv;
      phase_geom.injectConfigurationToPhase(pointwiseNJinv_cfg, inj_pointwiseNJinv);
      const CFG::LevelData<CFG::FluxBox>& bunit_cfg = mag_geom.getFCBFieldDir();
      LevelData<FluxBox> inj_bunit;
      phase_geom.injectConfigurationToPhase(bunit_cfg, inj_bunit);
   

      if (m_first_step) {
         
         // get face centered h_r, h_theta, and h_phi on each CFG_DIM face
         //
         m_lame_faces.define(dbl, 3, ghostVect);
         getFaceCenteredLameCoefficients(m_lame_faces, phase_geom, dbl);
         

         // set spatial depedence of transport coefficients
         // 
         CFG::LevelData<CFG::FArrayBox> D_cfg( mag_geom.grids(), 5, cfg_ghostVect );
         if (const_coeff) {
            for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
               for (int iD=0; iD<4; iD++) {
                  D_cfg[cfg_dit].setVal(D_kinet[iD],iD);
               }
               D_cfg[cfg_dit].setVal(DN0,4);
            } 
         } 
         else {
            setTransportCoeffSpatialDependence( D_cfg, mag_geom, time);
         }
         CFG::LevelData<CFG::FluxBox> D_cfg_faces(mag_geom.grids(), D_cfg.nComp(), cfg_ghostVect);
         fourthOrderCellToFaceCenters(D_cfg_faces, D_cfg);
         phase_geom.injectConfigurationToPhase( D_cfg_faces, D_kinet_faces);
         
      
         // compute the preconditioner coeffficient for the implicit solver
         //
         m_precond_D.define(mag_geom.grids(), CFG_DIM * CFG_DIM, CFG::IntVect::Unit);
         computePrecondCoefficient(m_precond_D, mag_geom, D_cfg_faces);
  
      }


      // create cell-centered dfB/dmu
      LevelData<FArrayBox> dfBdmu_cc(dbl, 1, 2*IntVect::Unit);
      DataIterator dit0 = dfBdmu_cc.dataIterator();
      for (dit0.begin(); dit0.ok(); ++dit0)
      {  
        // get phase space dx on this patch
        const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[dit0]);
        const RealVect& phase_dx =  block_coord_sys.dx();

        const FArrayBox& fB_on_patch = fB[dit0];
        FArrayBox& dfBdmu_on_patch = dfBdmu_cc[dit0];

        FORT_DFBDMU_CELL_CENTER( CHF_BOX(dfBdmu_on_patch.box()),
                                 CHF_CONST_REALVECT(phase_dx),
                                 CHF_CONST_INT(num_mu_cells),
                                 CHF_CONST_FRA1(fB_on_patch,0),
                                 CHF_FRA1(dfBdmu_on_patch,0) );
      }
      dfBdmu_cc.exchange(); 

      // calculate face-averaged flux
      LevelData<FluxBox> fluxNorm(dbl, 1, IntVect::Unit);
      LevelData<FluxBox> flux(dbl, SpaceDim, IntVect::Unit);
      const int simpleDiff = m_simple_diffusion? 1: 0;
      for (DataIterator dit( flux.dataIterator() ); dit.ok(); ++dit)
      {
        // get phase space dx on this patch
        const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[dit]);
        const RealVect& phase_dx =  block_coord_sys.dx();

        const FArrayBox& fB_on_patch = fB[dit];
        const FArrayBox& dfBdmu_on_patch = dfBdmu_cc[dit];
        const FArrayBox& B_on_patch  = inj_B[dit];
        const FArrayBox& N_on_patch  = density[dit];
        const FArrayBox& U_on_patch  = Upar[dit];
        const FArrayBox& T_on_patch  = temperature[dit];
        const FArrayBox& C_on_patch  = fourth_coef[dit];
        const FArrayBox& P_on_patch  = perp_coef[dit];
         
        for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisfluxNorm = fluxNorm[dit][dir];
          FArrayBox& thisflux = flux[dit][dir];
          const FArrayBox& lame_on_patch = m_lame_faces[dit][dir];
          const FArrayBox& bunit_on_patch  = inj_bunit[dit][dir];
          const FArrayBox& NJinv_on_patch = inj_pointwiseNJinv[dit][dir];
          const FArrayBox& D_kinet_on_patch = D_kinet_faces[dit][dir];

          if(!m_arbitrary_grid) {
            FORT_EVAL_ANOM_FLUX(CHF_CONST_INT(dir),
                          CHF_BOX(thisfluxNorm.box()),
                          CHF_CONST_REALVECT(phase_dx),
                          CHF_CONST_FRA(lame_on_patch),
                          CHF_CONST_FRA(D_kinet_on_patch),
                          CHF_CONST_REAL(mass),
                          CHF_CONST_INT(num_r_cells),
                          CHF_CONST_INT(simpleDiff),
                          CHF_CONST_FRA1(fB_on_patch,0),
                          CHF_CONST_FRA1(dfBdmu_on_patch,0),
                          CHF_CONST_FRA1(B_on_patch,0),
                          CHF_CONST_FRA1(N_on_patch,0),
                          CHF_CONST_FRA1(U_on_patch,0),
                          CHF_CONST_FRA1(T_on_patch,0),
                          CHF_CONST_FRA1(C_on_patch,0),
                          CHF_CONST_FRA1(P_on_patch,0),
                          CHF_FRA1(thisfluxNorm,0) );
          } else {
            //cout << "arb grid check 0 " << endl;
            FORT_EVAL_ANOM_FLUX_NOT_ALIGNED(CHF_CONST_INT(dir),
                          CHF_BOX(thisflux.box()),
                          CHF_CONST_REALVECT(phase_dx),
                          CHF_CONST_FRA(NJinv_on_patch),
                          CHF_CONST_FRA(bunit_on_patch),
                          CHF_CONST_FRA(D_kinet_on_patch),
                          CHF_CONST_REAL(mass),
                          CHF_CONST_INT(num_r_cells),
                          CHF_CONST_INT(simpleDiff),
                          CHF_CONST_FRA1(fB_on_patch,0),
                          CHF_CONST_FRA1(dfBdmu_on_patch,0),
                          CHF_CONST_FRA1(B_on_patch,0),
                          CHF_CONST_FRA1(N_on_patch,0),
                          CHF_CONST_FRA1(U_on_patch,0),
                          CHF_CONST_FRA1(T_on_patch,0),
                          CHF_CONST_FRA1(C_on_patch,0),
                          CHF_CONST_FRA1(P_on_patch,0),
                          CHF_FRA(thisflux) );
          }
        }
      }

      // calculate div(flux)
      LevelData<FArrayBox> rhs_transport;
      rhs_transport.define(rhs_dfn);
      if (!m_arbitrary_grid) {
         phase_geom.averageAtBlockBoundaries(fluxNorm);
         phase_geom.mappedGridDivergenceFromFluxNormals(rhs_transport, fluxNorm);
      }
      else {
         //cout << "arb grid check 1 " << endl;
         phase_geom.applyAxisymmetricCorrection(flux);
         phase_geom.averageAtBlockBoundaries(flux);
         const bool OMIT_NT(false);
         phase_geom.mappedGridDivergence( rhs_transport, flux, OMIT_NT );
         for (DataIterator dit( rhs_transport.dataIterator() ); dit.ok(); ++dit) {
            const PhaseBlockCoordSys& block_coord_sys( phase_geom.getBlockCoordSys( dbl[dit] ) );
            double fac( 1.0 / block_coord_sys.getMappedCellVolume() );
            rhs_transport[dit].mult( fac );
         }

      }
   
      // add (or overwrite) transport RHS to Vlasov RHS
      DataIterator rdit = rhs_transport.dataIterator();
      for (rdit.begin(); rdit.ok(); ++rdit)
      {
        if(model_only) {
           rhs_dfn[rdit].copy(rhs_transport[rdit]);
        }
        else  {
           //cout << "model plus " << endl;
           rhs_dfn[rdit].plus(rhs_transport[rdit]);
        }
      }

      m_first_step = false;
}

void Anomalous::getFaceCenteredLameCoefficients( LevelData<FluxBox>&       a_lame_faces,
                                                 const PhaseGeom&          a_phase_geom,
                                                 const DisjointBoxLayout&  a_dbl )
{
   // FluxBox to be filled needs 3 components on each face (hr, htheta, and hphi)
   CH_assert(a_lame_faces.nComp() == 3);
   DataIterator dit = a_lame_faces.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {

   // get face centered N and dX/dq components on this patch
   const CFG::MagBlockCoordSys& mag_block_coord_sys = a_phase_geom.getMagBlockCoordSys(a_dbl[dit]);
   const CFG::RealVect& real_dx = mag_block_coord_sys.dx();
   Box phase_box(a_dbl[dit]);
   CFG::Box box_config;
   a_phase_geom.projectPhaseToConfiguration(phase_box, box_config);
   box_config.grow(4);
   CFG::FluxBox N_cfg_face(box_config, CFG_DIM*CFG_DIM);  //4 components of N on each CFG_DIM face
   CFG::FluxBox dXdq_cfg_face(box_config, CFG_DIM*CFG_DIM);  //4 components on each CFG_DIM face
   mag_block_coord_sys.getPointwiseN(N_cfg_face);
   mag_block_coord_sys.getFaceCentereddXdxi(dXdq_cfg_face);
   FluxBox N_face, dXdq_face;
   a_phase_geom.injectConfigurationToPhase(N_cfg_face, N_face);
   a_phase_geom.injectConfigurationToPhase(dXdq_cfg_face, dXdq_face);
      
   for (int dir=0; dir<CFG_DIM; dir++)
   {
      FArrayBox& lame_on_patch = a_lame_faces[dit][dir];

      FORT_LAME_COEFFICIENTS(CHF_BOX(lame_on_patch.box()),
                             CHF_CONST_FRA(N_face[dir]),
                             CHF_CONST_FRA(dXdq_face[dir]),
                             CHF_FRA(lame_on_patch));
   }
   }

}

void Anomalous::getCellCenteredLameCoefficients( FArrayBox&        a_lame_cells,
                                                 const PhaseGeom&  a_phase_geom,
                                                 const Box&        a_dbl )
{
   // get cell centered N and dX/dq components
   const CFG::MagBlockCoordSys& mag_block_coord_sys = a_phase_geom.getMagBlockCoordSys(a_dbl);
   const CFG::RealVect& real_dx = mag_block_coord_sys.dx();
   Box phase_box(a_dbl);
   CFG::Box box_config;
   a_phase_geom.projectPhaseToConfiguration(phase_box, box_config);
   box_config.grow(2);
   CFG::FArrayBox N_cfg_cent(box_config, CFG_DIM*CFG_DIM);  //4 components of N on each cell center
   CFG::FArrayBox dXdq_cfg_cent(box_config, CFG_DIM*CFG_DIM);  //4 components on each cell center
   mag_block_coord_sys.getPointwiseN(N_cfg_cent);
   mag_block_coord_sys.getCellCentereddXdxi(dXdq_cfg_cent);
   FArrayBox N_cent, dXdq_cent;
   a_phase_geom.injectConfigurationToPhase(N_cfg_cent, N_cent);
   a_phase_geom.injectConfigurationToPhase(dXdq_cfg_cent, dXdq_cent);

   FORT_LAME_COEFFICIENTS( CHF_BOX(a_lame_cells.box()),
                           CHF_CONST_FRA(N_cent),
                           CHF_CONST_FRA(dXdq_cent),
                           CHF_FRA(a_lame_cells));

}

void Anomalous::setTransportCoeffSpatialDependence( CFG::LevelData<CFG::FArrayBox>&      a_D_kinet_cfg,
                                                    const CFG::MagGeom&  a_mag_geom,
                                                    const Real           a_time )
{
   CH_assert(a_D_kinet_cfg.nComp() == 5);
   CFG::DataIterator cfg_dit = a_D_kinet_cfg.dataIterator();
   const CFG::IntVect cfg_ghostVect(4*CFG::IntVect::Unit);
   CFG::LevelData<CFG::FArrayBox> D0_cfg, D1_cfg, D2_cfg, D3_cfg, DN_cfg;
   CFG::LevelData<CFG::FArrayBox> D0_fluid, D1_fluid, D2_fluid, D3_fluid;
   
   D0_cfg.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D1_cfg.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D2_cfg.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D3_cfg.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   DN_cfg.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   //
   D0_fluid.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D1_fluid.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D2_fluid.define( a_mag_geom.grids(), 1, cfg_ghostVect );
   D3_fluid.define( a_mag_geom.grids(), 1, cfg_ghostVect );


   // assign shape functions to CFG vectors
   //
   if( pptpm.contains("shape_function_D0") ) {
      m_shape_function_D0->assign( D0_cfg, a_mag_geom, a_time, false);
      m_shape_function_D0->assign( D0_fluid, a_mag_geom, a_time, false);
   } 
   else {
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D0_cfg[cfg_dit].setVal( 1.0 );
         D0_fluid[cfg_dit].setVal( 1.0 );
      }
   }
   if( pptpm.contains("shape_function_D1") ) {
      m_shape_function_D1->assign( D1_cfg, a_mag_geom, a_time, false);
      m_shape_function_D0->assign( D1_fluid, a_mag_geom, a_time, false);
   } 
   else {
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D1_cfg[cfg_dit].setVal( 1.0 );
         D1_fluid[cfg_dit].setVal( 1.0 );
      }
   }
   if( pptpm.contains("shape_function_D2") ) {
      m_shape_function_D2->assign( D2_cfg, a_mag_geom, a_time, false);
      m_shape_function_D0->assign( D2_fluid, a_mag_geom, a_time, false);
   } 
   else {
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D2_cfg[cfg_dit].setVal( 1.0 );
         D2_fluid[cfg_dit].setVal( 1.0 );
      }
   }
   if( pptpm.contains("shape_function_D3") ) {
      m_shape_function_D3->assign( D3_cfg, a_mag_geom, a_time, false);
      m_shape_function_D0->assign( D3_fluid, a_mag_geom, a_time, false);
   } 
   else {
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D3_cfg[cfg_dit].setVal( 1.0 );
         D3_fluid[cfg_dit].setVal( 1.0 );
      }
   }
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      DN_cfg[cfg_dit].setVal( 0.0 );
   }


   // multiply by coefficients set in input file
   //
   if (pptpm.contains("D_kinet")) { // spatial functions refer to kinetic transport coeffs
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D0_cfg[cfg_dit].mult( D_kinet[0] );
         D1_cfg[cfg_dit].mult( D_kinet[1] );
         D2_cfg[cfg_dit].mult( D_kinet[2] );
         D3_cfg[cfg_dit].mult( D_kinet[3] );
      }
   }


   // if spatial functions are for fluid coefficients, perform arithmetic
   //
   if (pptpm.contains("D_fluid")) { // spatial functions refer to fluid transport coeffs
      //CFG::LevelData<CFG::FArrayBox> D0_fluid, D1_fluid, D2_fluid, D3_fluid;
      //D0_fluid.define(D0_cfg);
      //D1_fluid.define(D1_cfg);
      //D2_fluid.define(D2_cfg);
      //D3_fluid.define(D3_cfg);
     
      for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
         D0_fluid[cfg_dit].mult( D_fluid[0] );
         D1_fluid[cfg_dit].mult( D_fluid[1] );
         D2_fluid[cfg_dit].mult( D_fluid[2] );
         D3_fluid[cfg_dit].mult( D_fluid[3] );
         if (D_fluid[0]==0) {
            D0_cfg[cfg_dit].copy( D3_fluid[cfg_dit].mult( 2.0/3.0 ) );
            D1_cfg[cfg_dit].copy( D1_fluid[cfg_dit].mult( -1.0 ) );
            D2_cfg[cfg_dit].copy( D2_fluid[cfg_dit].mult( -2.0/3.0 ) );
            D3_cfg[cfg_dit].copy( D1_fluid[cfg_dit] );
            DN_cfg[cfg_dit].copy( D0_fluid[cfg_dit] ); // this is free parameter used to insure D0 is not zero
         }
         else {
            D0_cfg[cfg_dit].copy( D0_fluid[cfg_dit] );               // DN0
            D1_cfg[cfg_dit].copy( D1_fluid[cfg_dit].mult( -1.0 ) );  // DT0
            D2_cfg[cfg_dit].copy( D2_fluid[cfg_dit].mult( -2.0/3.0 ) );
            D2_cfg[cfg_dit].plus( D0_fluid[cfg_dit] );            // DN2
            D3_cfg[cfg_dit].copy( D3_fluid[cfg_dit].mult( -2.0/3.0 ) );
            D3_cfg[cfg_dit].plus( D0_fluid[cfg_dit] );            
            D3_cfg[cfg_dit].plus( D1_fluid[cfg_dit] );            // DT2
         }
      }
   }
     

   // insert coefficients into appropriate components of Fluxbox
   //
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
      a_D_kinet_cfg[cfg_dit].copy(D0_cfg[cfg_dit],0,0,1);
      a_D_kinet_cfg[cfg_dit].copy(D1_cfg[cfg_dit],0,1,1);
      a_D_kinet_cfg[cfg_dit].copy(D2_cfg[cfg_dit],0,2,1);
      a_D_kinet_cfg[cfg_dit].copy(D3_cfg[cfg_dit],0,3,1);
      a_D_kinet_cfg[cfg_dit].copy(DN_cfg[cfg_dit],0,4,1);
   }

}


void Anomalous::computePrecondCoefficient( CFG::LevelData<CFG::FluxBox>&         a_D,
                                          const CFG::MagGeom&                    a_mag_geom,
                                          const CFG::LevelData<CFG::FluxBox>&    a_D_faces )
{
   
   CFG::LevelData<CFG::FluxBox> NJinverse(a_mag_geom.grids(), CFG_DIM*CFG_DIM, CFG::IntVect::Unit);
   a_mag_geom.getPointwiseNJinverse(NJinverse);
   
   CFG::LevelData<CFG::FluxBox> N(a_mag_geom.grids(), CFG_DIM*CFG_DIM, CFG::IntVect::Unit);
   a_mag_geom.getPointwiseN(N);
   
   const CFG::LevelData<CFG::FluxBox>& bunit = a_mag_geom.getFCBFieldDir();
         
   CFG::LevelData<CFG::FluxBox> a_D0(a_mag_geom.grids(), 1, CFG::IntVect::Unit );

   for (CFG::DataIterator dit(a_D.dataIterator()); dit.ok(); ++dit) {
      a_D0[dit].copy(a_D_faces[dit],0,0,1);      

      for (int dir=0; dir<CFG_DIM; ++dir) {

         const CFG::FArrayBox& this_D0 = a_D0[dit][dir];
         CFG::FArrayBox& this_coeff = a_D[dit][dir];
         CFG::FArrayBox& this_N = N[dit][dir];
         CFG::FArrayBox& this_NJinverse = NJinverse[dit][dir];
         const CFG::FArrayBox& this_b = bunit[dit][dir];
         const CFG::Box& box = this_coeff.box();

         //Get coefficient in the cylindrical coordinate frame
         CFG::FArrayBox phys_coeff(box, CFG_DIM * CFG_DIM);
         CFG::BoxIterator bit(box);
         for (bit.begin(); bit.ok(); ++bit) {
            CFG::IntVect iv = bit();
            
            double bNmag = sqrt(pow(this_b(iv,0),2) + pow(this_b(iv,2),2));
            CFG::RealVect bN;
            bN[RADIAL_DIR] = this_b(iv,2) / bNmag;
            bN[POLOIDAL_DIR] = -this_b(iv,0) / bNmag;
            
            phys_coeff(iv,0) = this_D0(iv,0) * bN[RADIAL_DIR]   * bN[RADIAL_DIR];
            phys_coeff(iv,1) = this_D0(iv,0) * bN[POLOIDAL_DIR] * bN[RADIAL_DIR];
            phys_coeff(iv,2) = this_D0(iv,0) * bN[RADIAL_DIR]   * bN[POLOIDAL_DIR];
            phys_coeff(iv,3) = this_D0(iv,0) * bN[POLOIDAL_DIR] * bN[POLOIDAL_DIR];
         
         }
         
         CFG::FArrayBox tmp(box, 1);
         // Multiply the coefficient matrix times the NJInverse matrix
         CFG::FArrayBox phys_coeff_times_NJinv(box, CFG_DIM * CFG_DIM);
         phys_coeff_times_NJinv.setVal(0.);
         for (int row=0; row<CFG_DIM; ++row) {
            for (int col=0; col<CFG_DIM; ++col) {
               for (int m=0; m<CFG_DIM; ++m) {
                  tmp.copy(phys_coeff,CFG_DIM*row+m,0,1);
                  tmp.mult(this_NJinverse,CFG_DIM*m+col,0,1);
                  phys_coeff_times_NJinv.plus(tmp,0,CFG_DIM*row+col,1);
               }
            }
         }
         
         //Premultiply by the NTranspose matrix
         this_coeff.setVal(0.);
         for (int row=0; row<CFG_DIM; ++row) {
            for (int col=0; col<CFG_DIM; ++col) {
               for (int m=0; m<CFG_DIM; ++m) {
                  tmp.copy(this_N,CFG_DIM*m+row,0,1);
                  tmp.mult(phys_coeff_times_NJinv,CFG_DIM*m+col,0,1);
                  this_coeff.plus(tmp,0,CFG_DIM*row+col,1);
               }
            }
         }
         
      }
   
   }
   
}


void Anomalous::ParseParameters()

{

   pptpm.query("model_only",model_only);
   pptpm.query("simple_diffusion",m_simple_diffusion);
   pptpm.query("verbosity",verbosity);
   CH_assert( verbosity == 0 || verbosity == 1);
   pptpm.query("arbitrary_grid",m_arbitrary_grid);

   
   if (pptpm.contains("const_coeff") ) {
     pptpm.get("const_coeff", const_coeff);

     if (!const_coeff) {
 
       CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
       std::string grid_function_name;
       //pptpm.query( "shape_function", grid_function_name );
       //m_shape_function = grid_library->find( grid_function_name );
      
       if( pptpm.query( "shape_function_D0", grid_function_name ) ) {
         m_shape_function_D0 = grid_library->find( grid_function_name );
       }
       if( pptpm.query( "shape_function_D1", grid_function_name ) ) {
         m_shape_function_D1 = grid_library->find( grid_function_name );
       }
       if ( pptpm.query( "shape_function_D2", grid_function_name ) ) {
         m_shape_function_D2 = grid_library->find( grid_function_name );
       }
       if( pptpm.query( "shape_function_D3", grid_function_name ) ) {
         m_shape_function_D3 = grid_library->find( grid_function_name );
       }
     
     }

   }

   // check to make sure coefficients properly set in input file
   if ( procID()==0 && !pptpm.contains("D_fluid") && !pptpm.contains("D_kinet") ) {
     cout << "NEITHER FLUID NOR KINETIC TRANSPORT COEFFICIENTS ARE SET FOR SPECIES " << species_name << endl;
     MayDay::Error( "TRANSPORT MODEL BEING USED, BUT COEFFICIENTS NOT PROPERLY DEFINED" );
   }
   if ( procID()==0 && pptpm.contains("D_fluid") && pptpm.contains("D_kinet") ) {
     cout << "BOTH FLUID AND KINETIC TRANSPORT COEFFICIENTS ARE SET FOR SPECIES " << species_name << endl;
     MayDay::Error( "TRANSPORT MODEL BEING USED, BUT COEFFICIENTS NOT PROPERLY DEFINED" );
   }
   if (pptpm.contains("D_fluid")) {
     pptpm.getarr("D_fluid",D_fluid,0,4);
     // convert fluid matrix components to kinetic coefficients
     if (D_fluid[0]==0) {
       D_kinet[0] = 2.0/3.0*D_fluid[3];
       D_kinet[1] = -D_fluid[1];
       D_kinet[2] = -2.0/3.0*D_fluid[2];
       D_kinet[3] = D_fluid[1];
       DN0 = D_kinet[0]; // this is free parameter used to insure D0 is not zero
     }
     else {
       D_kinet[0] = D_fluid[0];                                // D0
       D_kinet[1] = -D_fluid[1];                               // DT0
       D_kinet[2] = -2.0/3.0*D_fluid[2]+D_fluid[0];            // DN2
       D_kinet[3] = -2.0/3.0*D_fluid[3]+D_fluid[0]+D_fluid[1]; // DT2
     }
   }
   if (pptpm.contains("D_kinet")) {
     pptpm.getarr("D_kinet",D_kinet,0,4);
     // convert kinetic matrix components to fluid coefficients
     D_fluid[0] = D_kinet[0];                                  // DNN
     D_fluid[1] = -D_kinet[1];                                 // DNT
     D_fluid[2] = 3.0/2.0*(D_kinet[0]-D_kinet[2]);             // DTN
     D_fluid[3] = 3.0/2.0*(D_kinet[0]-D_kinet[1]-D_kinet[3]);  // DTT
   }

   // make sure that Anomalous Diff. coeff. D0 is positive
   if ( procID()==0 && (D_kinet[0] <= 0) ) {
     cout << "DIFFUSION COEFFICIENT D0 FOR SPECIES " << species_name << " IS <= ZERO " << endl;
     if (pptpm.contains("D_fluid")) {
       MayDay::Error( "FLUID COEFFICIENT D11 = D_fluid[0] MUST BE >= ZERO " );
     }
     else {
       MayDay::Error( "KINETIC DIFFUSION COEFFICIENT D0 = D_kinet[0] MUST BE > ZERO " );
     }
   }

}

/*
void Anomalous::computeStableDt( Real& stableDt, const KineticSpecies& soln ) {

   stableDt = computeDt(soln);
   // JRA prevent -nan return that occurs before first step and causes code not to run
   if(stableDt != stableDt) stableDt = 123456.7;

}
*/

inline void Anomalous::printParameters( const KineticSpecies& soln )
{

   Real dt_stable = computeStableDt(soln); 

   std::string tpt_only_str("false");
   std::string tpt_arbitraryGrid_str("true");
   std::string tpt_simpleDiffusion_str("false");
   if (model_only){ tpt_only_str = "true";}
   if (!m_arbitrary_grid){ tpt_arbitraryGrid_str = "false";}
   if (m_simple_diffusion){ tpt_simpleDiffusion_str = "true";}

   if (procID() == 0 ) {
     cout << "--------------------------------------------------------------------- " << endl;
     cout << "Anomalous Transport Model being defined for species " << species_name <<endl;
     cout << "tpt_only = " <<  tpt_only_str << endl;
     cout << "Arbitrary grid calculation  = " <<  tpt_arbitraryGrid_str << endl;
     cout << "simple diffusion model  = " <<  tpt_simpleDiffusion_str << endl;
     cout << "Fluid Tranpsort Matrix:  " << " D11 = " << D_fluid[0] << ", D12 = " << D_fluid[1] << endl;
     cout << "                         " << " D21 = " << D_fluid[2] << ", D22 = " << D_fluid[3] << endl;
     cout << "Kinetic Diffusion Coef:  " << " D0  = " << D_kinet[0] << endl;
     cout << "Kinetic Advection Coefs: " << " DN0 = " << DN0        << ", DT0  = " << D_kinet[1] << endl;
     cout << "                         " << " DN2 = " << D_kinet[2] << ", DT2  = " << D_kinet[3] << endl;
     //cout << "Stable Diffusive time step           = " << dt_stable_diff << endl;
     //cout << "Stable Advective-Diffusive time step = " << dt_stable_adv_diff << endl;
     cout << "Stable time step < " << dt_stable << endl;
     cout << "--------------------------------------------------------------------- " << endl;
   }
}

Real Anomalous::computeStableDt( const KineticSpecies& soln_species )
{
   // get stuff to calculate stability parameters on time step and number of mu cells
   //const KineticSpecies& soln_species( *(soln[0]) );
   const LevelData<FArrayBox>& soln_fB = soln_species.distributionFunction();
   const DisjointBoxLayout& dbl = soln_fB.getBoxes();
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
   const CFG::MultiBlockLevelGeom & mag_geom = phase_geom.magGeom();
   const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
   double mass = soln_species.mass();
   const ProblemDomain& phase_domain = phase_geom.domain();
   const Box& domain_box = phase_domain.domainBox();
   int num_r_cells = domain_box.size(0);

   // copy const soln_fB to a temporary 
   const IntVect ghostVect(4*IntVect::Unit);
   LevelData<FArrayBox> fB(dbl, 1, ghostVect); 
   DataIterator sdit = fB.dataIterator();
   for (sdit.begin(); sdit.ok(); ++sdit) {
     fB[sdit].copy(soln_fB[sdit]);
   }
   
   // get density, mean parallel velocity, temperature, and fourth moment coefficient
   const CFG::IntVect cfg_ghostVect(4*CFG::IntVect::Unit);
   CFG::LevelData<CFG::FArrayBox> dens_cfg(mag_geom.grids(), 1, cfg_ghostVect);
   CFG::LevelData<CFG::FArrayBox> Upar_cfg(mag_geom.grids(), 1, cfg_ghostVect);
   CFG::LevelData<CFG::FArrayBox> temp_cfg(mag_geom.grids(), 1, cfg_ghostVect);
   CFG::LevelData<CFG::FArrayBox> four_cfg(mag_geom.grids(), 1, cfg_ghostVect);
   CFG::LevelData<CFG::FArrayBox> perp_cfg(mag_geom.grids(), 1, cfg_ghostVect);

   moment_op.compute(dens_cfg, soln_species, fB, DensityKernel());
   moment_op.compute(Upar_cfg, soln_species, fB, ParallelMomKernel());
   moment_op.compute(four_cfg, soln_species, fB, FourthMomentKernel());
   moment_op.compute(perp_cfg, soln_species, fB, PerpEnergyKernel());
   CFG::DataIterator cfg_dit = dens_cfg.dataIterator();
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
     Upar_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
   }
   moment_op.compute(temp_cfg, soln_species, fB, PressureKernel(Upar_cfg));

   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
     four_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by pressure
     perp_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by pressure
   }
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
     temp_cfg[cfg_dit].divide(dens_cfg[cfg_dit]);
   }
   for (cfg_dit.begin(); cfg_dit.ok(); ++cfg_dit) {
     four_cfg[cfg_dit].divide(temp_cfg[cfg_dit]); // dividing by temperature
   }
   LevelData<FArrayBox> density;
   LevelData<FArrayBox> Upar;
   LevelData<FArrayBox> temperature;
   LevelData<FArrayBox> fourth_coef;
   LevelData<FArrayBox> perp_coef;
   phase_geom.injectConfigurationToPhase(dens_cfg, density);
   phase_geom.injectConfigurationToPhase(Upar_cfg, Upar);
   phase_geom.injectConfigurationToPhase(temp_cfg, temperature);
   phase_geom.injectConfigurationToPhase(four_cfg, fourth_coef);
   phase_geom.injectConfigurationToPhase(perp_cfg, perp_coef);

   // calculate cell-centered beta=U_psi^2/D_psi
   LevelData<FArrayBox> beta(dbl, 1, IntVect::Zero);
   DataIterator bdit = beta.dataIterator();
   for (bdit.begin(); bdit.ok(); ++bdit)
   {
     // get phase space dx on this patch
     const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[bdit]);
     const RealVect& phase_dx =  block_coord_sys.dx();

     const FArrayBox& B_on_patch  = inj_B[bdit];
     const FArrayBox& N_on_patch  = density[bdit];
     const FArrayBox& U_on_patch  = Upar[bdit];
     const FArrayBox& T_on_patch  = temperature[bdit];
     const FArrayBox& C_on_patch  = fourth_coef[bdit];
     const FArrayBox& P_on_patch  = perp_coef[bdit];

     FORT_EVAL_BETA( CHF_BOX(beta[bdit].box()),
                     CHF_CONST_REALVECT(phase_dx),
                     CHF_CONST_VR(D_kinet),
                     CHF_CONST_REAL(mass),
                     CHF_CONST_INT(num_r_cells),
                     CHF_CONST_FRA1(B_on_patch,0),
                     CHF_CONST_FRA1(N_on_patch,0),
                     CHF_CONST_FRA1(U_on_patch,0),
                     CHF_CONST_FRA1(T_on_patch,0),
                     CHF_CONST_FRA1(C_on_patch,0),
                     CHF_CONST_FRA1(P_on_patch,0),
                     CHF_FRA1(beta[bdit],0) );
   }

   // get r lame coefficient
   LevelData<FArrayBox> hr(dbl, 1, IntVect::Zero);
   DataIterator dit= hr.dataIterator();
   Real dr;
   for (dit.begin(); dit.ok(); ++dit)
   {
     const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[dit]);
     const RealVect& phase_dx =  block_coord_sys.dx();
     dr = phase_dx[0];

     // get cell-centered h_r, h_theta, and h_phi on this patch
     Box phase_box(dbl[dit]);
     FArrayBox lame_cells(phase_box, 3);
     getCellCenteredLameCoefficients(lame_cells, phase_geom, dbl[dit]);
     hr[dit].copy(lame_cells,0,0);

   }
   hr.exchange();

   // convert fluid matrix components to Anomalous coefficients
   // const Real D2 = 2.0/3.0*D_fluid[2]-D_fluid[0];
   //const Real D0 = 2.0/3.0*D_fluid[3]-D_fluid[1] - 7.0/2.0*D2;
   //const Real D1 = 2.0/3.0*D_fluid[3]-D_fluid[1] + D_fluid[0] - four/three*D_fluid[2];
   //const Real D3 = 3.0/2.0*D2 - D_fluid[1];
   //const Real D4 = -2.0*D2;
   //const Real D5 = -4.0/3.0*( 2.0/3.0*D_fluid[3]-D_fluid[1] );

   // get the min value of hr and Dpsi and max values of Dpsi and Upsi
   Real local_min_hr(1000);  // 1000 is just a starting point
   Real local_max_beta(-1000);

   for (dit.begin(); dit.ok(); ++dit)
   {
     Box box( dbl[dit] );
     Real box_min_hr( hr[dit].min (box) );
     Real box_max_beta( beta[dit].max (box) );
     local_min_hr = Min( local_min_hr, box_min_hr );
     local_max_beta = Max( local_max_beta, box_max_beta );
   }

   Real min_hr( local_min_hr );
   Real max_beta( local_max_beta );
#ifdef CH_MPI
   MPI_Allreduce( &local_min_hr, &min_hr, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
   MPI_Allreduce( &local_max_beta, &max_beta, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif

   // calculate min time step set by purely diffusive parts
   // dt*sum_i( D_i/dx_i^2 ) << 1/2
   const Real alpha_r = D_kinet[0]/(min_hr*min_hr*dr*dr);
   const Real dt_stable_diff = 0.5/alpha_r;

   // calculate min time step set by advective-diffusive
   // dt*sum_i( U_i^2/D_i ) << 2
   //const Real beta_r = max_Upsi*max_Upsi/min_Dpsi;
   const Real dt_stable_adv_diff = 2.0/max_beta;
   
   // take mininum stable time step
   Real stableDt = Min(dt_stable_diff, dt_stable_adv_diff);

   // JRA prevent -nan return that occurs before first step and causes code not to run   
   if(stableDt != stableDt) stableDt = 123456.7;

   return stableDt;

}

#include "NamespaceFooter.H"
