#include "FokkerPlanck.H"
#include "RosenbluthPotentials.H"


#include "FourthOrderUtil.H"
#include "Directions.H"
#include "EdgeToCell.H"
#include "PhaseGeom.H"
#include "KineticFunctionLibrary.H"
#include "inspect.H"
#include "ConstFact.H"

#include "FokkerPlanckF_F.H"

#include "NamespaceHeader.H" 


FokkerPlanck::FokkerPlanck( ParmParse& a_ppcls, const int a_verbosity )
   : m_verbosity(a_verbosity),
     m_ppcls(a_ppcls),
     m_cls_freq(-1.0),
     m_pcg_tol(1.0e-5),
     m_pcg_maxiter(100),
     m_mult_num(1),
     m_nD(5),
     m_update_freq(-1),
     m_it_counter(0),
     m_first_step(true),
     m_subtract_background(false),
     m_debug(false)
{
   m_verbosity = true;
   parseParameters( a_ppcls );
   if (m_verbosity>0) {
      printParameters();
   }

   m_fixed_cls_freq = (m_cls_freq<0.0) ? false : true ;

}


FokkerPlanck::~FokkerPlanck()
{
}

void FokkerPlanck::preTimeStep(const KineticSpeciesPtrVect& a_soln, 
                               const int  a_species,
                               const Real a_time)
{

   // Get solution distribution function (J*Bstar_par*dfn) for the current species
   const KineticSpecies& soln_species( *(a_soln[a_species]) );
   const LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );

   const DisjointBoxLayout& grids( soln_dfn.getBoxes() );
   const int n_comp( soln_dfn.nComp() );

   //Get test-particle (tp) and field-particle (fp) masses
   const double mass_tp = soln_species.mass();
   const double mass_fp = soln_species.mass();

   // Get coordinate system parameters 
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();

   //Compute normalization
   if (m_first_step) {computeClsNorm(m_cls_norm, soln_species.mass(), soln_species.charge());}
 
   //Update phi?
   bool update_phi(false);
   if (m_update_freq < 0) {update_phi = true;}
   else if (m_it_counter % m_update_freq == 0 ) { update_phi=true; }

   if (!m_subtract_background) {

     LevelData<FArrayBox> dfn( grids, n_comp, IntVect::Zero );
     for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      dfn[dit].copy( soln_dfn[dit] );
     }
     //Covert delta_dfn to cell centers
     convertToCellCenters(phase_geom, dfn);

     //Define m_phi and set it to zero at the first time step
     if (m_first_step) {
       m_phi.define( grids, 2, IntVect::Zero );
       m_D.define( grids, m_nD, 2*IntVect::Unit );
       for (DataIterator dit(soln_dfn.dataIterator()); dit.ok(); ++dit) {
         m_phi[dit].setVal(0.0);
         m_D[dit].setVal(0.0);
       }
     }
     
     if ( update_phi || m_first_step ) {
       evalRosenbluthPotentials(m_phi, phase_geom, dfn, mass_tp); 
       evalCoefficients(m_D,m_phi,phase_geom,mass_tp,mass_fp);
     }

   } else {

     // Create reference (J*Bstar_par*dfn_bckg) distribution
     KineticSpeciesPtr ref_species( soln_species.clone( IntVect::Zero, false ) );
     m_ref_func->assign( *ref_species, a_time );
     LevelData<FArrayBox>& ref_dfn( ref_species->distributionFunction() );

     // Compute the difference from the reference (background) solution
     LevelData<FArrayBox> delta_dfn( grids, n_comp, IntVect::Zero );
     for (DataIterator dit(soln_dfn.dataIterator()); dit.ok(); ++dit) {
        delta_dfn[dit].copy( soln_dfn[dit] );
        delta_dfn[dit].minus( ref_dfn[dit] );
     }

     //Covert delta_dfn to cell centers
     convertToCellCenters(phase_geom, delta_dfn);

     //Define m_phi and m_phi_F0, and compute m_phi_F0 
     if (m_first_step) {
       m_phi.define( grids, 2, IntVect::Zero );
       m_phi_F0.define( grids, 2, IntVect::Zero );
       m_D.define( grids, m_nD, 2*IntVect::Unit );
       m_D_F0.define( grids, m_nD, 2*IntVect::Unit );
       for (DataIterator dit(soln_dfn.dataIterator()); dit.ok(); ++dit) {
         m_phi[dit].setVal(0.0);
         m_phi_F0[dit].setVal(0.0);
         m_D[dit].setVal(0.0);
         m_D_F0[dit].setVal(0.0);
       }
       evalRosenbluthPotentials(m_phi_F0, phase_geom, ref_dfn, mass_tp);
       evalCoefficients(m_D_F0,m_phi_F0,phase_geom,mass_tp,mass_fp);
     } 

     if ( update_phi ) {
       evalRosenbluthPotentials(m_phi, phase_geom, delta_dfn, mass_tp); 
       evalCoefficients(m_D,m_phi,phase_geom,mass_tp,mass_fp);
     }
   }

   //Update iteration counter
   m_it_counter+=1;
   m_first_step = false;
}

void FokkerPlanck::evalClsRHS( KineticSpeciesPtrVect& a_rhs,
                        const KineticSpeciesPtrVect&  a_soln,
                        const int                     a_species,
                        const Real                    a_time,
                        const int                     a_flag )
{
  if (m_debug) evalClsRHS_LowOrder(a_rhs,a_soln,a_species,a_time,a_flag);
  else         evalClsRHS_Main    (a_rhs,a_soln,a_species,a_time,a_flag);
  return;
}

// NB: a_soln is on the computational grid and has 4 ghost cells (passed here as Nans)
// a_rhs has probably zero ghost cells (from accertion in computedivergence)  (double check)
void FokkerPlanck::evalClsRHS_Main( KineticSpeciesPtrVect& a_rhs,
                        const KineticSpeciesPtrVect&  a_soln,
                        const int                     a_species,
                        const Real                    a_time,
                        const int                     a_flag )
{

   // Get solution distribution function (J*Bstar_par*dfn) for the current species
   const KineticSpecies& soln_species( *(a_soln[a_species]) );
   const LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );

   const DisjointBoxLayout& grids( soln_dfn.getBoxes() );
   const int n_comp( soln_dfn.nComp() );

   // Get coordinate system parameters 
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();

   // Create collisional flux 
   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Zero);

   if (!m_subtract_background) {

     LevelData<FArrayBox> dfn( grids, n_comp, IntVect::Zero );
     for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      dfn[dit].copy( soln_dfn[dit] );
     }
     //Covert delta_dfn to cell centers
     convertToCellCenters(phase_geom, dfn);

     computeFlux(flux, phase_geom, m_D, dfn);

   } else {

     // Create reference (J*Bstar_par*dfn_bckg) distribution
     KineticSpeciesPtr ref_species( soln_species.clone( IntVect::Zero, false ) );
     m_ref_func->assign( *ref_species, a_time );
     LevelData<FArrayBox>& ref_dfn( ref_species->distributionFunction() );

     // Compute the difference from the reference (background) solution
     LevelData<FArrayBox> delta_dfn( grids, n_comp, IntVect::Zero );
     for (DataIterator dit(soln_dfn.dataIterator()); dit.ok(); ++dit) {
        delta_dfn[dit].copy( soln_dfn[dit] );
        delta_dfn[dit].minus( ref_dfn[dit] );
     }

     //Covert delta_dfn to cell centers
     convertToCellCenters(phase_geom, delta_dfn);

     //Compute C[F1,F0]+C[F0,F1]+C[F1,F1]
     LevelData<FluxBox> flux_tmp(grids, SpaceDim, IntVect::Zero);

     //Compute C[F1,F0]
     computeFlux(flux_tmp, phase_geom, m_D_F0, delta_dfn);
     for (DataIterator dit(soln_dfn.dataIterator()); dit.ok(); ++dit) {
        flux[dit].copy( flux_tmp[dit] );
     }

     computeFlux(flux_tmp, phase_geom, m_D, ref_dfn);
     for (DataIterator dit(soln_dfn.dataIterator()); dit.ok(); ++dit) {
       for (int dir=0; dir<SpaceDim; dir++) {
         flux[dit][dir].plus( flux_tmp[dit][dir] );
       }
     }

     //Compute C[F1,F0] + C[F0,F1] + C[F1,F1]
     computeFlux(flux_tmp, phase_geom, m_D, delta_dfn);
     for (DataIterator dit(soln_dfn.dataIterator()); dit.ok(); ++dit) {
       for (int dir=0; dir<SpaceDim; dir++) {
         flux[dit][dir].plus( flux_tmp[dit][dir] );
       }
     }
   }

   //Add collisional contribution to rhs
   KineticSpecies& rhs_species( *(a_rhs[a_species]) );
   LevelData<FArrayBox>& rhs_dfn( rhs_species.distributionFunction() );
  
   LevelData<FArrayBox> rhs_cls(grids, n_comp, IntVect::Zero);
   phase_geom.mappedGridDivergence(rhs_cls, flux, true);
  
   for (DataIterator dit( rhs_cls.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(grids[dit]);
      double fac = 1. / block_coord_sys.getMappedCellVolume();
      rhs_cls[dit].mult( fac );
      if (m_fixed_cls_freq) {rhs_cls[dit].mult( m_cls_freq );}
      else {rhs_cls[dit].mult( m_cls_norm );}     
      rhs_dfn[dit].plus( rhs_cls[dit] ); 
    }
}

void FokkerPlanck::evalClsRHS_LowOrder( KineticSpeciesPtrVect& a_rhs,
                        const KineticSpeciesPtrVect&  a_soln,
                        const int                     a_species,
                        const Real                    a_time,
                        const int                     a_flag )
{ 
   // Get solution distribution function (J*Bstar_par*dfn) for the current species
   const KineticSpecies& soln_species( *(a_soln[a_species]) );
   const LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );

   const DisjointBoxLayout& grids( soln_dfn.getBoxes() );
   const int n_comp( soln_dfn.nComp() );

   // Get coordinate system parameters 
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();

   // Create collisional flux 
   LevelData<FArrayBox> flux(grids, 2, IntVect::Unit);

   // Create reference (J*Bstar_par*dfn_bckg) distribution
   KineticSpeciesPtr ref_species( soln_species.clone( IntVect::Zero, false ) );
   m_ref_func->assign( *ref_species, a_time );
   LevelData<FArrayBox>& ref_dfn( ref_species->distributionFunction() );

   if (!m_subtract_background) {

     LevelData<FArrayBox> dfn( grids, n_comp, IntVect::Zero );
     for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      dfn[dit].copy( soln_dfn[dit] );
     }

     computeFluxCellCentered(flux, phase_geom, m_D, dfn, ref_dfn);

   } else {

     // Compute the difference from the reference (background) solution
     LevelData<FArrayBox> delta_dfn( grids, n_comp, IntVect::Zero );
     for (DataIterator dit(soln_dfn.dataIterator()); dit.ok(); ++dit) {
        delta_dfn[dit].copy( soln_dfn[dit] );
        delta_dfn[dit].minus( ref_dfn[dit] );
     }

     //Compute C[F1,F0]+C[F0,F1]+C[F1,F1]
     LevelData<FArrayBox> flux_tmp(grids, 2, IntVect::Unit);

     //Compute C[F1,F0]
     computeFluxCellCentered(flux_tmp, phase_geom, m_D_F0, delta_dfn, ref_dfn);
     for (DataIterator dit(soln_dfn.dataIterator()); dit.ok(); ++dit) {
        flux[dit].copy( flux_tmp[dit] );
     }

     //Compute C[F1,F0] + C[F0,F1]
     computeFluxCellCentered(flux_tmp, phase_geom, m_D, ref_dfn, ref_dfn);
     for (DataIterator dit(soln_dfn.dataIterator()); dit.ok(); ++dit) {
       flux[dit].plus( flux_tmp[dit] );
     }

     //Compute C[F1,F0] + C[F0,F1] + C[F1,F1]
     computeFluxCellCentered(flux_tmp, phase_geom, m_D, delta_dfn, ref_dfn);
     for (DataIterator dit(soln_dfn.dataIterator()); dit.ok(); ++dit) {
       flux[dit].plus( flux_tmp[dit] );
     }
   }
  
   //Add collisional contribution to rhs
   KineticSpecies& rhs_species( *(a_rhs[a_species]) );
   LevelData<FArrayBox>& rhs_dfn( rhs_species.distributionFunction() );
  
   LevelData<FArrayBox> rhs_cls(grids, n_comp, IntVect::Zero);
   computeDivergence(rhs_cls,phase_geom,flux);
  
   for (DataIterator dit( rhs_cls.dataIterator() ); dit.ok(); ++dit) {
      if (m_fixed_cls_freq) {rhs_cls[dit].mult( m_cls_freq );}
      else {rhs_cls[dit].mult( m_cls_norm );}     
      rhs_dfn[dit].plus( rhs_cls[dit] ); 
   }
}

void FokkerPlanck::assemblePrecondMatrix( void *a_P,
                                          const KineticSpeciesPtrVect& a_soln,
                                          int a_species,
                                          const Mapping& a_mapping)
{
  BandedMatrix *Pmat = (BandedMatrix*) a_P;

  const KineticSpecies&       soln_species(*(a_soln[a_species]));
  const LevelData<FArrayBox>& soln_dfn    (soln_species.distributionFunction());
  const DisjointBoxLayout&    grids       (soln_dfn.disjointBoxLayout());
  const PhaseGeom&            phase_geom  (soln_species.phaseSpaceGeometry());
  const int                   n_comp      (soln_dfn.nComp());
  const VEL::VelCoordSys&     vel_coords  (phase_geom.velSpaceCoordSys());
  const VEL::RealVect&        vel_dx      (vel_coords.dx());
  const LevelData<FArrayBox>& pMapping    (a_mapping.getPointMapping()); 

  Real  dv = 1.0/vel_dx[0], dmu = 1.0/vel_dx[1], dv_sq = dv*dv, dmu_sq = dmu*dmu;
  DataIterator dit = grids.dataIterator();
  int count = 0;
  for (dit.begin(); dit.ok(); ++dit) {
    const Box& grid = grids[dit];
    const FArrayBox& pMap = pMapping[dit];

    /* grid size */
    IntVect bigEnd   = grid.bigEnd(),
            smallEnd = grid.smallEnd();
    IntVect gridSize(bigEnd); gridSize -= smallEnd; gridSize += 1;

    BoxIterator bit(grid);
    for (bit.begin(); bit.ok(); ++bit) {
      /* this point */
      IntVect ic = bit();
      /* neighboring points */
      IntVect ie(ic),
              iw(ic),
              in(ic),
              is(ic),
              ine(ic),
              inw(ic),
              ise(ic),
              isw(ic);
      /* north-south is along mu; east-west is along v|| */
      ie[SpaceDim-2]++;                          /* east  */
      iw[SpaceDim-2]--;                          /* west  */
      in[SpaceDim-1]++;                          /* north */
      is[SpaceDim-1]--;                          /* south */
      ine[SpaceDim-2]++; ine[SpaceDim-1]++;   /* northeast */
      inw[SpaceDim-2]--; inw[SpaceDim-1]++;   /* northwest */
      ise[SpaceDim-2]++; ise[SpaceDim-1]--;   /* southeast */
      isw[SpaceDim-2]--; isw[SpaceDim-1]--;   /* southwest */
      /* global row/column numbers */
      int pc, pe, pw, pn, ps, pne, pnw, pse, psw, pcl;
      pc  = (int) pMap.get(ic ,0);
      pn  = (int) pMap.get(in ,0);
      ps  = (int) pMap.get(is ,0);
      pe  = (int) pMap.get(ie ,0);
      pw  = (int) pMap.get(iw ,0);
      pne = (int) pMap.get(ine,0);
      pnw = (int) pMap.get(inw,0);
      pse = (int) pMap.get(ise,0);
      psw = (int) pMap.get(isw,0);
      pcl = (int) pMap.get(ic ,1);

      /* coefficients */
      Real *D_c,
           *D_e,
           *D_w,
           *D_n,
           *D_s,
           *D_ne,
           *D_nw,
           *D_se,
           *D_sw,
           ac, an, as, ae, aw, ane, anw, ase, asw;

      D_c  = new Real[m_nD];
      D_e  = new Real[m_nD];
      D_w  = new Real[m_nD];
      D_n  = new Real[m_nD];
      D_s  = new Real[m_nD];
      D_ne = new Real[m_nD];
      D_nw = new Real[m_nD];
      D_se = new Real[m_nD];
      D_sw = new Real[m_nD];

      m_D[dit].getVal(&D_c [0],ic );
      m_D[dit].getVal(&D_e [0],ie );
      m_D[dit].getVal(&D_w [0],iw );
      m_D[dit].getVal(&D_n [0],in );
      m_D[dit].getVal(&D_s [0],is );
      m_D[dit].getVal(&D_ne[0],ine);
      m_D[dit].getVal(&D_nw[0],inw);
      m_D[dit].getVal(&D_se[0],ise);
      m_D[dit].getVal(&D_sw[0],isw);

      if (m_subtract_background) {
        Real *D0_c,
             *D0_e,
             *D0_w,
             *D0_n,
             *D0_s,
             *D0_ne,
             *D0_nw,
             *D0_se,
             *D0_sw;

        D0_c  = new Real[m_nD];
        D0_e  = new Real[m_nD];
        D0_w  = new Real[m_nD];
        D0_n  = new Real[m_nD];
        D0_s  = new Real[m_nD];
        D0_ne = new Real[m_nD];
        D0_nw = new Real[m_nD];
        D0_se = new Real[m_nD];
        D0_sw = new Real[m_nD];

        m_D_F0[dit].getVal(&D0_c [0],ic );
        m_D_F0[dit].getVal(&D0_e [0],ie );
        m_D_F0[dit].getVal(&D0_w [0],iw );
        m_D_F0[dit].getVal(&D0_n [0],in );
        m_D_F0[dit].getVal(&D0_s [0],is );
        m_D_F0[dit].getVal(&D0_ne[0],ine);
        m_D_F0[dit].getVal(&D0_nw[0],inw);
        m_D_F0[dit].getVal(&D0_se[0],ise);
        m_D_F0[dit].getVal(&D0_sw[0],isw);

        for (int v=0; v<m_nD; v++) {
          D_c [v] += D0_c [v];
          D_n [v] += D0_n [v];
          D_s [v] += D0_s [v];
          D_e [v] += D0_e [v];
          D_w [v] += D0_w [v];
          D_ne[v] += D0_ne[v];
          D_nw[v] += D0_nw[v];
          D_se[v] += D0_se[v];
          D_sw[v] += D0_sw[v];
        }

        delete[] D0_c;
        delete[] D0_e;
        delete[] D0_w;
        delete[] D0_n;
        delete[] D0_s;
        delete[] D0_ne;
        delete[] D0_nw;
        delete[] D0_se;
        delete[] D0_sw;
      }

      /*
       * D[0] -> D_v
       * D[1] -> D_mu
       * D[2] -> D_v_v
       * D[3] -> D_v_mu
       * D[4] -> D_mu_mu
       */

      ac = ( - D_c[0]*dv - 2*D_c[1]*dmu - D_c[2]*dv_sq - D_e[2]*dv_sq 
             - 4*D_c[4]*dmu_sq - 4*D_n[4]*dmu_sq - 4*D_c[3]*dv*dmu );
      an = ( 2*D_n[1]*dmu + 2*D_n[3]*dv*dmu + 4*D_n[4]*dmu_sq );
      as = ( 4*D_c[4]*dmu_sq + 2*D_c[3]*dv*dmu );
      ae = ( D_e[0]*dv + D_e[2]*dv_sq + 2*D_e[3]*dv*dmu );
      aw = ( D_c[2]*dv_sq + 2*D_c[3]*dv*dmu );
      ane =  0.0;
      anw =  -2*D_n[3]*dv*dmu;
      ase =  -2*D_e[3]*dv*dmu;
      asw =  0.0;

      delete[] D_c;
      delete[] D_e;
      delete[] D_w;
      delete[] D_n;
      delete[] D_s;
      delete[] D_ne;
      delete[] D_nw;
      delete[] D_se;
      delete[] D_sw;

      if (m_fixed_cls_freq) {
        ac  *= m_cls_freq;
        ae  *= m_cls_freq;
        aw  *= m_cls_freq;
        an  *= m_cls_freq;
        as  *= m_cls_freq;
        ane *= m_cls_freq;
        anw *= m_cls_freq;
        ase *= m_cls_freq;
        asw *= m_cls_freq;
      } else {
        ac  *= m_cls_norm;
        ae  *= m_cls_norm;
        aw  *= m_cls_norm;
        an  *= m_cls_norm;
        as  *= m_cls_norm;
        ane *= m_cls_norm;
        anw *= m_cls_norm;
        ase *= m_cls_norm;
        asw *= m_cls_norm;
      }

      int  ncols = 9, bs = n_comp*n_comp, ix = 0;
      int  *icols = (int*)  calloc (ncols,sizeof(int));
      Real *data  = (Real*) calloc (ncols*bs,sizeof(Real));

      /* center element */
      icols[ix] = pc; 
      for (int v=0; v<bs; v++) data[ix*bs+v] = ac;
      ix++;
      /* east element */
      if (pe >= 0) {
        icols[ix] = pe;
        for (int v=0; v<bs; v++) data[ix*bs+v] = ae;
        ix++;
      }
      /* west element */
      if (pw >= 0) {
        icols[ix] = pw;
        for (int v=0; v<bs; v++) data[ix*bs+v] = aw;
        ix++;
      }
      /* north element */
      if (pn >= 0) {
        icols[ix] = pn;
        for (int v=0; v<bs; v++) data[ix*bs+v] = an;
        ix++;
      }
      /* south element */
      if (ps >= 0) {
        icols[ix] = ps;
        for (int v=0; v<bs; v++) data[ix*bs+v] = as;
        ix++;
      }
      /* north east element */
      if (pne >= 0) {
        icols[ix] = pne;
        for (int v=0; v<bs; v++) data[ix*bs+v] = ane;
        ix++;
      }
      /* north west element */
      if (pnw >= 0) {
        icols[ix] = pnw;
        for (int v=0; v<bs; v++) data[ix*bs+v] = anw;
        ix++;
      }
      /* south east element */
      if (pse >= 0) {
        icols[ix] = pse;
        for (int v=0; v<bs; v++) data[ix*bs+v] = ase;
        ix++;
      }
      /* south west element */
      if (psw >= 0) {
        icols[ix] = psw;
        for (int v=0; v<bs; v++) data[ix*bs+v] = asw;
        ix++;
      }

      Pmat->setRowValues(pcl,pc,ix,icols,data);
      free(data);
      free(icols);
    }
    count++;
  }

  return;
}

void FokkerPlanck::evalRosenbluthPotentials( LevelData<FArrayBox>& a_phi,
                                             const PhaseGeom& a_phase_geom,
                                             const LevelData<FArrayBox>& a_dfn,
                                             const double a_mass ) const 

{

   const DisjointBoxLayout& grids( a_phi.getBoxes() );

   //Create temporary phi_one and phi_two 
   LevelData<FArrayBox> phi_one(a_dfn.disjointBoxLayout(),1,IntVect::Zero);
   LevelData<FArrayBox> phi_two(a_dfn.disjointBoxLayout(),1,IntVect::Zero);
   for (DataIterator dit( a_phi.dataIterator() ); dit.ok(); ++dit) {
      phi_one[dit].copy(a_phi[dit], grids[dit], 0 , grids[dit], 0, 1);
      phi_two[dit].copy(a_phi[dit], grids[dit], 1 , grids[dit], 0, 1);
   }


   //Calculate Rosenbluth Potentials 
   RosenbluthPotentials RosenbluthPotentials(phi_one, phi_two, a_dfn, a_phase_geom, 
                                             a_mass, m_pcg_tol, m_pcg_maxiter, m_mult_num, m_verbosity);

   for (DataIterator dit( a_phi.dataIterator() ); dit.ok(); ++dit) {
      a_phi[dit].copy(phi_one[dit], grids[dit], 0 , grids[dit], 0, 1);
      a_phi[dit].copy(phi_two[dit], grids[dit], 0 , grids[dit], 1, 1);
   }
}

void FokkerPlanck::evalCoefficients( LevelData<FArrayBox>& a_D,
                                     const LevelData<FArrayBox>& a_phi,
                                     const PhaseGeom& a_phase_geom,
                                     const double a_mass_tp,
                                     const double a_mass_fp) const
{
   // Get velocity coordinate system parameters
   const VEL::VelCoordSys&    vel_coords = a_phase_geom.velSpaceCoordSys();
   const VEL::ProblemDomain&  vel_domain = vel_coords.domain();
   const VEL::Box&            domain_box = vel_domain.domainBox();
   const VEL::RealVect&       vel_dx     = vel_coords.dx();

   const int num_vpar_cells = domain_box.size(0);
   const int num_mu_cells   = domain_box.size(1);

   //Create temporary phi with two extra layers of ghost cells (filled with zeros)
   LevelData<FArrayBox> phi_tmp(a_phi.disjointBoxLayout(),
                                a_phi.nComp(),
                                a_phi.ghostVect()+2*IntVect::Unit);

   const DisjointBoxLayout& grids( a_phi.getBoxes() );
   for (DataIterator dit( a_phi.dataIterator() ); dit.ok(); ++dit) {
      phi_tmp[dit].setVal(0.0);
      phi_tmp[dit].copy(a_phi[dit],grids[dit]);
   }
   phi_tmp.exchange();
   FillGhostCells(a_phase_geom, phi_tmp);
   
   const LevelData<FArrayBox>& injected_B = a_phase_geom.getBFieldMagnitude();
   for (DataIterator dit( a_D.dataIterator() ); dit.ok(); ++dit) {
      
      // Check if the have more than 4 cells on a box in MU_DIR 
      // used for calcuation of second derivatives in mu at mu=0 bnd   
      const Box& grid_box = grids[dit];
      if (grid_box.size(3)<5) { 
       MayDay::Error(" Box size in MU_DIR must be greater that 4 for FP operator calculations");
      }  

      FArrayBox& this_D = a_D[dit];
      const FArrayBox& this_phi_tmp = phi_tmp[dit];
      const FArrayBox& this_b = injected_B[dit];
      
      FORT_EVALUATE_COEFFICIENTS( CHF_FRA(this_D),
                                  CHF_CONST_FRA(this_phi_tmp),
                                  CHF_CONST_FRA1(this_b,0),
                                  CHF_BOX(a_phi[dit].box()),
                                  CHF_CONST_REALVECT(vel_dx),
                                  CHF_CONST_REAL(a_mass_tp),
                                  CHF_CONST_REAL(a_mass_fp),
                                  CHF_CONST_INT(num_vpar_cells),
                                  CHF_CONST_INT(num_mu_cells));
   }
   /* Divide by JB* */
   a_phase_geom.divideJonValid(a_D);
   a_phase_geom.divideBStarParallel(a_D);
   FillGhostCells(a_phase_geom,a_D);
}

void FokkerPlanck::computeFlux( LevelData<FluxBox>& a_flux,
                                const PhaseGeom& a_phase_geom,
                                const LevelData<FArrayBox>& a_D,
                                const LevelData<FArrayBox>& a_dfn ) const
{

   // Get velocity coordinate system parameters
   const VEL::VelCoordSys& vel_coords = a_phase_geom.velSpaceCoordSys();
   const VEL::ProblemDomain& vel_domain = vel_coords.domain();
   const VEL::Box& domain_box = vel_domain.domainBox();
   const VEL::RealVect& vel_dx = vel_coords.dx();
   const int num_vpar_cells = domain_box.size(0);
   const int num_mu_cells = domain_box.size(1);


   //Create temporary dfn and phi with two extra layers of ghost cells (filled with zeros)
   LevelData<FArrayBox> dfn_tmp(a_dfn.disjointBoxLayout(),
                                a_dfn.nComp(),
                                a_dfn.ghostVect()+2*IntVect::Unit);

   const DisjointBoxLayout& grids( a_dfn.getBoxes() );
   for (DataIterator dit( a_dfn.dataIterator() ); dit.ok(); ++dit) {
      dfn_tmp[dit].setVal(0.0);
      dfn_tmp[dit].copy(a_dfn[dit],grids[dit]);
   }
   dfn_tmp.exchange();

   //Compute cell-cenetered collision fluxes (0 comp - vpar_dir, 1 comp - mu_dir)
   LevelData<FArrayBox> flux_cell(grids, 2, IntVect::Zero);
   const LevelData<FArrayBox>& injected_B = a_phase_geom.getBFieldMagnitude();

   for (DataIterator dit( a_dfn.dataIterator() ); dit.ok(); ++dit) {
      
      // Check if the have more than 4 cells on a box in MU_DIR 
      // used for calcuation of second derivatives in mu at mu=0 bnd   
      const Box& grid_box = grids[dit];
      if (grid_box.size(3)<5) { 
       MayDay::Error(" Box size in MU_DIR must be greater that 4 for FP operator calculations");
      }  

      FArrayBox& this_flux_cell = flux_cell[dit];
      FArrayBox& this_dfn_tmp   = dfn_tmp[dit];
      const FArrayBox& this_D   = a_D[dit];
      
      FORT_EVALUATE_FLUX_CELL(CHF_FRA(this_flux_cell),
                              CHF_FRA1(this_dfn_tmp,0),
                              CHF_CONST_FRA(this_D),
                              CHF_BOX(this_flux_cell.box()),
                              CHF_CONST_REALVECT(vel_dx),
                              CHF_CONST_INT(num_vpar_cells),
                              CHF_CONST_INT(num_mu_cells));
   }

   //Calculate collision flux cell-average
   convertToCellAverage(a_phase_geom, flux_cell);

   //Calculate collision fluxes on cell faces 
   LevelData<FluxBox> flux_face(grids, 2, IntVect::Zero);
   convertToCellFaces(flux_face, flux_cell, a_phase_geom);

   //Create final (combined r , theta, mu, and v_par) flux
   for (DataIterator dit(a_flux.dataIterator() ); dit.ok(); ++dit) {
      const FArrayBox& this_flux_cell = flux_cell[dit];

      FArrayBox tmp_flux_cell(grow(this_flux_cell.box(),1),2);
      tmp_flux_cell.setVal(0.);
      tmp_flux_cell.copy(this_flux_cell);

      for (int dir=0; dir<SpaceDim; dir++) {
          FArrayBox& this_flux = a_flux[dit][dir];
          FArrayBox& this_flux_face = flux_face[dit][dir];
          FORT_EVALUATE_FLUX_COMBINED(CHF_FRA(this_flux),
                                           CHF_CONST_INT(dir),
                                           CHF_BOX(this_flux.box()),
                                           CHF_CONST_FRA(this_flux_face),
                                           CHF_CONST_FRA(tmp_flux_cell),
                                           CHF_CONST_INT(num_vpar_cells),
                                           CHF_CONST_INT(num_mu_cells));
      }
   }
}

void FokkerPlanck::computeFluxCellCentered( LevelData<FArrayBox>& a_flux,
                                            const PhaseGeom& a_phase_geom,
                                            const LevelData<FArrayBox>& a_D,
                                            const LevelData<FArrayBox>& a_dfn,
                                            const LevelData<FArrayBox>& a_dfn_ref) const
{
   // Get velocity coordinate system parameters
   const VEL::VelCoordSys&    vel_coords  = a_phase_geom.velSpaceCoordSys();
   const VEL::ProblemDomain&  vel_domain  = vel_coords.domain();
   const VEL::Box&            domain_box  = vel_domain.domainBox();
   const VEL::RealVect&       vel_dx      = vel_coords.dx();

   const int num_vpar_cells = domain_box.size(0);
   const int num_mu_cells   = domain_box.size(1);

   //Create temporary dfn with two extra layers of ghost cells (filled with zeros)
   LevelData<FArrayBox> dfn_tmp(a_dfn.disjointBoxLayout(),
                                a_dfn.nComp(),
                                2*IntVect::Unit);

   const DisjointBoxLayout& grids( a_dfn.getBoxes() );
   for (DataIterator dit( a_dfn.dataIterator() ); dit.ok(); ++dit) {
      dfn_tmp[dit].setVal(0.0);
      dfn_tmp[dit].copy(a_dfn_ref[dit],grids[dit]);
   }
   FillGhostCells(a_phase_geom,dfn_tmp);
   for (DataIterator dit( a_dfn.dataIterator() ); dit.ok(); ++dit) {
      dfn_tmp[dit].copy(a_dfn[dit],grids[dit]);
   }
   dfn_tmp.exchange();

   //Compute cell-cenetered collision fluxes (0 comp - vpar_dir, 1 comp - mu_dir)
   for (DataIterator dit( a_dfn.dataIterator() ); dit.ok(); ++dit) {
      
      FArrayBox& this_flux_cell     = a_flux[dit];
      const FArrayBox& this_dfn_tmp = dfn_tmp[dit];
      const FArrayBox& this_D       = a_D[dit];
      
      FORT_EVALUATE_FLUX_CELL_LOW_ORDER(CHF_FRA(this_flux_cell),
                                        CHF_CONST_FRA1(this_dfn_tmp,0),
                                        CHF_CONST_FRA(this_D),
                                        CHF_BOX(this_flux_cell.box()),
                                        CHF_CONST_REALVECT(vel_dx),
                                        CHF_CONST_INT(num_vpar_cells),
                                        CHF_CONST_INT(num_mu_cells));
   }
   a_flux.exchange();
}

void FokkerPlanck::computeDivergence( LevelData<FArrayBox>& a_rhs,
                                      const PhaseGeom& a_phase_geom,
                                      const LevelData<FArrayBox>& a_flux) const
{
   // Get velocity coordinate system parameters
   const VEL::VelCoordSys&    vel_coords  = a_phase_geom.velSpaceCoordSys();
   const VEL::ProblemDomain&  vel_domain  = vel_coords.domain();
   const VEL::Box&            domain_box  = vel_domain.domainBox();
   const VEL::RealVect&       vel_dx      = vel_coords.dx();

   const int num_vpar_cells = domain_box.size(0);
   const int num_mu_cells   = domain_box.size(1);

   //Compute cell-cenetered divergence 
   for (DataIterator dit( a_rhs.dataIterator() ); dit.ok(); ++dit) {
      
      const FArrayBox& this_flux = a_flux[dit];
      FArrayBox& this_rhs  = a_rhs[dit];
      
      FORT_EVALUATE_FOKKERPLANCK_RHS( CHF_FRA1(this_rhs,0),
                                      CHF_CONST_FRA(this_flux),
                                      CHF_BOX(this_rhs.box()),
                                      CHF_CONST_REALVECT(vel_dx),
                                      CHF_CONST_INT(num_vpar_cells),
                                      CHF_CONST_INT(num_mu_cells));
   }
}

void FokkerPlanck::convertToCellAverage(const PhaseGeom&  a_phase_geom,
                                      LevelData<FArrayBox>&      a_dfn) const
{
   // Convert from cell centered values to cell averages.  We need to use
   // a temporary with an extra ghost cell for the interpolation

   LevelData<FArrayBox> dfn_tmp(a_dfn.disjointBoxLayout(),
                                a_dfn.nComp(),
                                a_dfn.ghostVect()+IntVect::Unit);

   const DisjointBoxLayout& grids = a_dfn.disjointBoxLayout();

   DataIterator dit(grids.dataIterator());
   for (dit.begin(); dit.ok(); ++dit) {
      dfn_tmp[dit].copy(a_dfn[dit]);
   }
   dfn_tmp.exchange();

   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = coord_sys.domain();

      fourthOrderAverageCell(dfn_tmp[dit], domain, grids[dit]);
   }
   dfn_tmp.exchange();

   for (dit.begin(); dit.ok(); ++dit) {
      a_dfn[dit].copy(dfn_tmp[dit]);
   }
}

void FokkerPlanck::convertToCellCenters(const PhaseGeom&  a_phase_geom,
                                      LevelData<FArrayBox>&      a_dfn) const
{
   // Convert from cell centered values to cell averages.  We need to use
   // a temporary with an extra ghost cell for the interpolation

   LevelData<FArrayBox> dfn_tmp(a_dfn.disjointBoxLayout(),
                                a_dfn.nComp(),
                                a_dfn.ghostVect()+IntVect::Unit);

   const DisjointBoxLayout& grids = a_dfn.disjointBoxLayout();

   DataIterator dit(grids.dataIterator());
   for (dit.begin(); dit.ok(); ++dit) {
      dfn_tmp[dit].copy(a_dfn[dit]);
   }
   dfn_tmp.exchange();

   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = coord_sys.domain();

      fourthOrderAverageCell(dfn_tmp[dit], domain, grids[dit], -1);
   }
   dfn_tmp.exchange();

   for (dit.begin(); dit.ok(); ++dit) {
      a_dfn[dit].copy(dfn_tmp[dit]);
   }
}

void FokkerPlanck::convertToCellFaces(LevelData<FluxBox>& a_faceData,
                                 const LevelData<FArrayBox>& a_cellData,
                                 const PhaseGeom&  a_phase_geom) const
{
   // Convert from cell centered values to cell averages.  We need to use
   // a temporary with two extra ghost cells for the interpolation. The ghost cells at the domain 
   // boundaries are filled with zeros for our purposes. Since, the boundary fluxes in the velocity directions 
   // will be overwrited later, and the fluxes in the configuration domansions do not contribute to the RHS. 

   LevelData<FArrayBox> dfn_tmp(a_cellData.disjointBoxLayout(),
                                a_cellData.nComp(),
                                a_cellData.ghostVect()+2*IntVect::Unit);

   const DisjointBoxLayout& grids = a_cellData.disjointBoxLayout();

   DataIterator dit(grids.dataIterator());
   for (dit.begin(); dit.ok(); ++dit) {
      dfn_tmp[dit].setVal(0.0);
      dfn_tmp[dit].copy(a_cellData[dit]);
   }
   dfn_tmp.exchange();   

// FillGhostCells(a_phase_geom, dfn_tmp);

   fourthOrderCellToFace(a_faceData, dfn_tmp);
}

void FokkerPlanck::FillGhostCells(const PhaseGeom&  a_phase_geom,
                                  LevelData<FArrayBox>& a_dfn) const
{

// NB: fourthOrderCellExtrapAtDomainBdry does not fill "corners"
 
   const DisjointBoxLayout& grids = a_dfn.disjointBoxLayout();

   DataIterator dit(grids.dataIterator());
   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = coord_sys.domain();

      fourthOrderCellExtrapAtDomainBdry(a_dfn[dit], domain, grids[dit]);
   }
   a_dfn.exchange();
}

void FokkerPlanck::computeClsNorm(Real&                       a_cls_norm,
                                  const double                a_mass,
                                  const double                a_charge ) const

{
    //Get normalization parameters (units)
    double N, T, L, M;
    ParmParse ppunits( "units" );
    ppunits.get("number_density",N);  //[m^{-3}]
    ppunits.get("length",L);          //[m]
    ppunits.get("temperature",T);     //[eV]
    ppunits.get("mass",M);            //[Amu]

    //Fundamental constants  
    double pi = Constants::PI;
    double ech = Constants::ELEMENTARY_CHARGE;
    double eps0 = Constants::VACUUM_PERMITTIVITY;
    double Mkg = M * Constants::MASS_OF_PROTON;

    //Derived normalizations
    double V = sqrt (ech * T / Mkg);  //[m/s]
    double F = N / pow(V, 3) / pi;    //[s^3/m^6]
    double tau = L / V;               //[s]
    
    //Compute collisional normalization 
    double Coulomb_Lg = 23 - log( sqrt(2.0) * pow(a_charge,3) * sqrt(N)/1000.0 / pow(T, 3.0/2.0) ); 
    a_cls_norm = F * tau * Coulomb_Lg * pow(a_charge * ech, 4) / pow( a_mass * Mkg * eps0, 2);

}


inline
void FokkerPlanck::parseParameters( ParmParse& a_ppcls )
{
   a_ppcls.query( "cls_freq", m_cls_freq );
   a_ppcls.query( "subtract_background", m_subtract_background );
   a_ppcls.query( "update_frequency", m_update_freq);
   a_ppcls.query( "verbose", m_verbosity);
   a_ppcls.query( "debug",m_debug);

   a_ppcls.query( "convergence_tolerance", m_pcg_tol );
   a_ppcls.query( "max_interation_number", m_pcg_maxiter);
   a_ppcls.query( "multipole_number", m_mult_num);

   KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   std::string function_name;
   a_ppcls.query( "ref_function", function_name );
   m_ref_func = library->find( function_name );

}


inline
void FokkerPlanck::printParameters()
{
   if (procID()==0) {
      std::cout << "FokkerPlanck collisions parameters:" << std::endl;
      std::cout << "  cls_freq  =  " << m_cls_freq
                << ", subtract_background = " << m_subtract_background;
      std::cout << "  Reference Function:" << std::endl;
      m_ref_func->printParameters();
   }
}


Real FokkerPlanck::computeDt(const KineticSpeciesPtrVect& soln)
{
   return DBL_MAX;
}

#include "NamespaceFooter.H"
