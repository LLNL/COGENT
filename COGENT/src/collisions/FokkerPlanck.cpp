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
#include "FPLimitersF_F.H"

#include "GKVector.H"
#include "GKOps.H"

#include "NamespaceHeader.H" 


FokkerPlanck::FokkerPlanck( const std::string& a_ppcls_str, const int a_verbosity )
   : m_verbosity((const bool)a_verbosity),
     m_time_implicit(true),
     m_limiters(false),
     m_cls_freq(-1.0),
     m_pcg_tol(1.0e-5),
     m_pcg_maxiter(100),
     m_mult_num(1),
     m_nD(5),
     m_update_freq(-1),
     m_it_counter(0),
     m_nbands(9),
     m_my_pc_idx(-1),
     m_rosenbluth_skip_stage(false),
     m_subtract_background(false),
     m_ref_updated(false),
     m_compute_maxwellian(false),
     m_conserve_energy(false),
     m_fp_energy_cons_epsilon(1.0e-16)
{
   m_verbosity = true;
   m_opt_string = a_ppcls_str;
   ParmParse m_ppcls(a_ppcls_str.c_str());
   parseParameters( m_ppcls );

   m_fixed_cls_freq = (m_cls_freq<0.0) ? false : true ;

   if (m_subtract_background) {
     /* Need to do this till a more robust energy-conservation fix
      * is implemented that will work for subtract-background option
      * as well */
     m_conserve_energy = false;
   }

   if (m_verbosity == true) {
      printParameters();
   }
}

FokkerPlanck::~FokkerPlanck()
{
}

static inline int sign(Real x)
{
  return(x < 0 ? -1 : 1);
}

void FokkerPlanck::defineBlockPC( std::vector<Preconditioner<GKVector,GKOps>*>& a_pc,
                                  std::vector<DOFList>&                         a_dof_list,
                                  const GKVector&                               a_X,
                                  GKOps&                                        a_ops,
                                  const std::string&                            a_out_string,
                                  const std::string&                            a_opt_string,
                                  bool                                          a_im,
                                  const KineticSpecies&                         a_species,
                                  const GlobalDOFKineticSpecies&                a_global_dofs,
                                  const int                                     a_species_index )
{
  if (a_im && m_time_implicit) {
  
    CH_assert(a_pc.size() == a_dof_list.size());

    if (!procID()) {
      std::cout << "  Kinetic Species " << a_species_index
                << " - FP collision term: "
                << " creating " << _BASIC_GK_PC_ << " preconditioner "
                << " (index = " << a_pc.size() << ").\n";
    }

    Preconditioner<GKVector, GKOps> *pc;
    pc = new BasicGKPreconditioner<GKVector,GKOps>;
    DOFList dof_list(0);

    const LevelData<FArrayBox>& soln_dfn    (a_species.distributionFunction());
    const DisjointBoxLayout&    grids       (soln_dfn.disjointBoxLayout());
    const int                   n_comp      (soln_dfn.nComp());
    const LevelData<FArrayBox>& pMapping    (a_global_dofs.data()); 
  
    for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const Box& grid = grids[dit];
      const FArrayBox& pMap = pMapping[dit];
      for (BoxIterator bit(grid); bit.ok(); ++bit) {
        IntVect ic = bit();
        for (int n(0); n < n_comp; n++) {
          dof_list.push_back((int) pMap.get(ic ,n) - a_global_dofs.mpiOffset());
        }
      }
    }

    m_my_pc_idx = a_pc.size();

    dynamic_cast<BasicGKPreconditioner<GKVector,GKOps>*>
      (pc)->define(a_X, a_ops, m_opt_string, m_opt_string, a_im, dof_list);
    a_pc.push_back(pc);
    a_dof_list.push_back(dof_list);

  }

  return;
}

void FokkerPlanck::updateBlockPC( std::vector<Preconditioner<GKVector,GKOps>*>& a_pc,
                                  const KineticSpecies&                         a_species,
                                  const GlobalDOFKineticSpecies&                a_global_dofs,
                                  const Real                                    a_time,
                                  const Real                                    a_shift,
                                  const bool                                    a_im,
                                  const int                                     a_species_index )
{
  if (a_im && m_time_implicit) {
    CH_assert(m_my_pc_idx >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx);

    if (!procID()) {
      std::cout << "    ==> Updating " <<_BASIC_GK_PC_ << " preconditioner "
                << " for FP collision term of kinetic species " << a_species_index << ".\n";
    }

    BasicGKPreconditioner<GKVector,GKOps> *pc = dynamic_cast<BasicGKPreconditioner<GKVector,GKOps>*>
                                                  (a_pc[m_my_pc_idx]);
    CH_assert(pc != NULL);

    BandedMatrix *pc_mat((BandedMatrix*)pc->getBandedMatrix());
    pc_mat->setToIdentityMatrix();
    assemblePrecondMatrix((void*)pc_mat, a_species, a_global_dofs, a_shift);
    pc_mat->finalAssembly();
  }
  return;
}

void FokkerPlanck::assemblePrecondMatrix( void                            *a_P,
                                          const KineticSpecies&           a_species,
                                          const GlobalDOFKineticSpecies&  a_global_dofs,
                                          const Real                      a_shift )
{
  /*
   * To ensure that the entries of the preconditioning matrix are correct, do the 
   * following:-
   * 1) Comment out convertToCellCenters() and convertToCellAverages() in evalClsRHS().
   * 2) Set efac = 1.0 in eval_rhs() in FokkerPlanckF.ChF
   * 3) In FokkerPlanckF.ChF and FPLimitersF.ChF, uncomment the low-order discretizations
   *    and comment out the high-order ones.
   * 4) Comment out all fillGhostCellsDfn*() function calls.
   * 5) Ensure that the option "rosenbluth_skip_stage" for this collision operator is
   *    set to "true" in the input file.
   * Then, the preconditioning matrix assembled here will be the exact Jacobian, and 
   * if solved exactly (say, with large number of Gauss-Seidel iterations), will cause
   * GMRES to converge in 1 iteration.
  */

  if (m_time_implicit) {

    BandedMatrix *Pmat = (BandedMatrix*) a_P;
  
    const LevelData<FArrayBox>& soln_dfn    (a_species.distributionFunction());
    const DisjointBoxLayout&    grids       (soln_dfn.disjointBoxLayout());
    const PhaseGeom&            phase_geom  (a_species.phaseSpaceGeometry());
    const int                   n_comp      (soln_dfn.nComp());
    const VEL::VelCoordSys&     vel_coords  (phase_geom.velSpaceCoordSys());
    const VEL::RealVect&        vel_dx      (vel_coords.dx());
    //    const VEL::ProblemDomain&   vel_domain  = vel_coords.domain();
    //    const VEL::Box&             domain_box  = vel_domain.domainBox();
    const LevelData<FArrayBox>& pMapping    (a_global_dofs.data()); 
  
    //    const int nvpar = domain_box.size(0);
    //    const int nmu   = domain_box.size(1);
  
    Real  dv = 1.0/vel_dx[0], dmu = 1.0/vel_dx[1], dv_sq = dv*dv, dmu_sq = dmu*dmu;
  
    DataIterator dit = grids.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
      const Box& grid = grids[dit];
      const FArrayBox& pMap = pMapping[dit];
  
      /* grid size */
      IntVect bigEnd   = grid.bigEnd(),
              smallEnd = grid.smallEnd();
      IntVect gridSize(bigEnd); gridSize -= smallEnd; gridSize += 1;
  
      FArrayBox *Dv, *Dmu, *Dvv, *Dmumu, *Dvmu;
      Dv    = m_coeffs.D_v[dit].dataPtr();
      Dmu   = m_coeffs.D_mu[dit].dataPtr();
      Dvv   = m_coeffs.D_vv[dit].dataPtr();
      Dmumu = m_coeffs.D_mumu[dit].dataPtr();
      Dvmu  = m_coeffs.D_vmu[dit].dataPtr();
  
      FArrayBox *Dv_F0, *Dmu_F0, *Dvv_F0, *Dmumu_F0, *Dvmu_F0;
      if (m_subtract_background) {
        Dv_F0    = m_coeffs_F0.D_v[dit].dataPtr();
        Dmu_F0   = m_coeffs_F0.D_mu[dit].dataPtr();
        Dvv_F0   = m_coeffs_F0.D_vv[dit].dataPtr();
        Dmumu_F0 = m_coeffs_F0.D_mumu[dit].dataPtr();
        Dvmu_F0  = m_coeffs_F0.D_vmu[dit].dataPtr();
      } else {
        Dv_F0 = Dmu_F0 = Dvv_F0 = Dmumu_F0 = Dvmu_F0 = NULL;
      }
  
      BoxIterator bit(grid);
      for (bit.begin(); bit.ok(); ++bit) {
        /* this point */
        IntVect ic = bit();
        /* neighboring points */
        IntVect ie(ic);
        IntVect iw(ic);
        IntVect in(ic);
        IntVect is(ic);
        IntVect ine(ic);
        IntVect inw(ic);
        IntVect ise(ic);
        IntVect isw(ic);
        /* north-south is along mu; east-west is along v|| */
        ie[_VPAR_DIM_]++;                     /* east  */
        iw[_VPAR_DIM_]--;                     /* west  */
        in[_MU_DIM_]++;                       /* north */
        is[_MU_DIM_]--;                       /* south */
        ine[_VPAR_DIM_]++; ine[_MU_DIM_]++;   /* northeast */
        inw[_VPAR_DIM_]--; inw[_MU_DIM_]++;   /* northwest */
        ise[_VPAR_DIM_]++; ise[_MU_DIM_]--;   /* southeast */
        isw[_VPAR_DIM_]--; isw[_MU_DIM_]--;   /* southwest */
  
        Real DvL, DvR, DmuL, DmuR, DvvL, DvvR, DmumuL, DmumuR;
        Real DvmuWS, DvmuES, DvmuWN, DvmuEN;
  
        DvL  = Dv->get(ic,0);
        DvR  = Dv->get(ie,0);
        DmuL = Dmu->get(ic,0);
        DmuR = Dmu->get(in,0);
  
        DvvL   = Dvv->get(ic,0);
        DvvR   = Dvv->get(ie,0);
        DmumuL = Dmumu->get(ic,0);
        DmumuR = Dmumu->get(in,0);
  
        DvmuWS = Dvmu->get(ic ,0);
        DvmuES = Dvmu->get(ie ,0);
        DvmuWN = Dvmu->get(in ,0);
        DvmuEN = Dvmu->get(ine,0);
  
        Real DvL_F0, DvR_F0, DmuL_F0, DmuR_F0, DvvL_F0, DvvR_F0, DmumuL_F0, DmumuR_F0;
        Real DvmuWS_F0, DvmuES_F0, DvmuWN_F0, DvmuEN_F0;
  
        if (m_subtract_background) {
          DvL_F0  = Dv_F0->get(ic,0);
          DvR_F0  = Dv_F0->get(ie,0);
          DmuL_F0 = Dmu_F0->get(ic,0);
          DmuR_F0 = Dmu_F0->get(in,0);
  
          DvvL_F0   = Dvv_F0->get(ic,0);
          DvvR_F0   = Dvv_F0->get(ie,0);
          DmumuL_F0 = Dmumu_F0->get(ic,0);
          DmumuR_F0 = Dmumu_F0->get(in,0);
  
          DvmuWS_F0 = Dvmu_F0->get(ic ,0);
          DvmuES_F0 = Dvmu_F0->get(ie ,0);
          DvmuWN_F0 = Dvmu_F0->get(in ,0);
          DvmuEN_F0 = Dvmu_F0->get(ine,0);
        } else {
          DvL_F0  = 0.0;
          DvR_F0  = 0.0;
          DmuL_F0 = 0.0;
          DmuR_F0 = 0.0;
  
          DvvL_F0   = 0.0;
          DvvR_F0   = 0.0;
          DmumuL_F0 = 0.0;
          DmumuR_F0 = 0.0;
  
          DvmuWS_F0 = 0.0;
          DvmuES_F0 = 0.0;
          DvmuWN_F0 = 0.0;
          DvmuEN_F0 = 0.0;
        }
  
        for (int n(0); n < n_comp; n++) {
          /* global row/column numbers */
          int pc, pe, pw, pn, ps, pne, pnw, pse, psw;//, pee, pww, pnn, pss;
          pc  = (int) pMap.get(ic ,n);
          pn  = (int) pMap.get(in ,n);
          ps  = (int) pMap.get(is ,n);
          pe  = (int) pMap.get(ie ,n);
          pw  = (int) pMap.get(iw ,n);
          pne = (int) pMap.get(ine,n);
          pnw = (int) pMap.get(inw,n);
          pse = (int) pMap.get(ise,n);
          psw = (int) pMap.get(isw,n);
  
          /* coefficients */
          Real nu = (m_fixed_cls_freq ? m_cls_freq : m_cls_norm);
  
          Real ac = 0;
          Real an = 0;
          Real as = 0;
          Real ae = 0;
          Real aw = 0;
          Real ane = 0;
          Real anw = 0;
          Real ase = 0;
          Real asw = 0;
          
          /* Advection terms along vpar */
          if ((DvL < 0.0) && (DvR < 0.0)) {
            ac += nu * DvR * dv;
            aw -= nu * DvL * dv;
          } else if ((DvL >= 0.0) && (DvR < 0.0)) {
            ac += nu * (DvR - DvL) * dv;
          } else if ((DvL < 0.0) && (DvR >= 0.0)) {
            ae += nu * DvR * dv;
            aw -= nu * DvL * dv;
          } else {
            ae += nu * DvR * dv;
            ac -= nu * DvL * dv;
          }
          if (m_subtract_background) {
            if ((DvL_F0 < 0.0) && (DvR_F0 < 0.0)) {
              ac += nu * DvR_F0 * dv;
              aw -= nu * DvL_F0 * dv;
            } else if ((DvL_F0 >= 0.0) && (DvR_F0 < 0.0)) {
              ac += nu * (DvR_F0 - DvL_F0) * dv;
            } else if ((DvL_F0 < 0.0) && (DvR_F0 >= 0.0)) {
              ae += nu * DvR_F0 * dv;
              aw -= nu * DvL_F0 * dv;
            } else {
              ae += nu * DvR_F0 * dv;
              ac -= nu * DvL_F0 * dv;
            }
          }
          
          /* Advection terms along mu */
          if ((DmuL < 0.0) && (DmuR < 0.0)) {
            ac += 2.0 * nu * DmuR * dmu;
            as -= 2.0 * nu * DmuL * dmu;
          } else if ((DmuL >= 0.0) && (DmuR < 0.0)) {
            ac += 2.0 * nu * (DmuR - DmuL) * dmu;
          } else if ((DmuL < 0.0) && (DmuR >= 0.0)) {
            an += 2.0 * nu * DmuR * dmu;
            as -= 2.0 * nu * DmuL * dmu;
          } else {
            an += 2.0 * nu * DmuR * dmu;
            ac -= 2.0 * nu * DmuL * dmu;
          }
          if (m_subtract_background) {
            if ((DmuL_F0 < 0.0) && (DmuR_F0 < 0.0)) {
              ac += 2.0 * nu * DmuR_F0 * dmu;
              as -= 2.0 * nu * DmuL_F0 * dmu;
            } else if ((DmuL_F0 >= 0.0) && (DmuR_F0 < 0.0)) {
              ac += 2.0 * nu * (DmuR_F0 - DmuL_F0) * dmu;
            } else if ((DmuL_F0 < 0.0) && (DmuR_F0 >= 0.0)) {
              an += 2.0 * nu * DmuR_F0 * dmu;
              as -= 2.0 * nu * DmuL_F0 * dmu;
            } else {
              an += 2.0 * nu * DmuR_F0 * dmu;
              ac -= 2.0 * nu * DmuL_F0 * dmu;
            }
          }
  
          /* Laplacian term along vpar */
          ae += nu * DvvR * dv_sq;
          aw += nu * DvvL * dv_sq;
          ac -= nu * (DvvR + DvvL) * dv_sq;
          if (m_subtract_background) {
            ae += nu * DvvR_F0 * dv_sq;
            aw += nu * DvvL_F0 * dv_sq;
            ac -= nu * (DvvR_F0 + DvvL_F0) * dv_sq;
          }
          
          /* Laplacian term along mu */
          an += 4.0 * nu * DmumuR * dmu_sq;
          as += 4.0 * nu * DmumuL * dmu_sq;
          ac -= 4.0 * nu * (DmumuR + DmumuL) * dmu_sq;
          if (m_subtract_background) {
            an += 4.0 * nu * DmumuR_F0 * dmu_sq;
            as += 4.0 * nu * DmumuL_F0 * dmu_sq;
            ac -= 4.0 * nu * (DmumuR_F0 + DmumuL_F0) * dmu_sq;
          }
  
          /* Cross term 1: d/dv(Dvmu df/dmu) */
          ac  += 0.5 * nu * dv * dmu * (DvmuES - DvmuEN - DvmuWS + DvmuWN);
          ae  += 0.5 * nu * dv * dmu * (DvmuES - DvmuEN);
          aw  += 0.5 * nu * dv * dmu * (-DvmuWS + DvmuWN);
          an  += 0.5 * nu * dv * dmu * (DvmuEN - DvmuWN);
          as  += 0.5 * nu * dv * dmu * (-DvmuES + DvmuWS);
          asw += 0.5 * nu * dv * dmu * (DvmuWS);
          ase += 0.5 * nu * dv * dmu * (-DvmuES);
          anw += 0.5 * nu * dv * dmu * (-DvmuWN);
          ane += 0.5 * nu * dv * dmu * (DvmuEN);
          if (m_subtract_background) {
            ac  += 0.5 * nu * dv * dmu * (DvmuES_F0 - DvmuEN_F0 - DvmuWS_F0 + DvmuWN_F0);
            ae  += 0.5 * nu * dv * dmu * (DvmuES_F0 - DvmuEN_F0);
            aw  += 0.5 * nu * dv * dmu * (-DvmuWS_F0 + DvmuWN_F0);
            an  += 0.5 * nu * dv * dmu * (DvmuEN_F0 - DvmuWN_F0);
            as  += 0.5 * nu * dv * dmu * (-DvmuES_F0 + DvmuWS_F0);
            asw += 0.5 * nu * dv * dmu * (DvmuWS_F0);
            ase += 0.5 * nu * dv * dmu * (-DvmuES_F0);
            anw += 0.5 * nu * dv * dmu * (-DvmuWN_F0);
            ane += 0.5 * nu * dv * dmu * (DvmuEN_F0);
          }
  
          /* Cross term 2: d/dmu(Dvmu df/dv) */
          ac  += 0.5 * nu * dv * dmu * (DvmuWN - DvmuEN - DvmuWS + DvmuES);
          ae  += 0.5 * nu * dv * dmu * (DvmuEN - DvmuES);
          aw  += 0.5 * nu * dv * dmu * (-DvmuWN + DvmuWS);
          an  += 0.5 * nu * dv * dmu * (DvmuWN - DvmuEN);
          as  += 0.5 * nu * dv * dmu * (-DvmuWS + DvmuES);
          asw += 0.5 * nu * dv * dmu * (DvmuWS);
          ase += 0.5 * nu * dv * dmu * (-DvmuES);
          anw += 0.5 * nu * dv * dmu * (-DvmuWN);
          ane += 0.5 * nu * dv * dmu * (DvmuEN);
          if (m_subtract_background) {
            ac  += 0.5 * nu * dv * dmu * (DvmuWN_F0 - DvmuEN_F0 - DvmuWS_F0 + DvmuES_F0);
            ae  += 0.5 * nu * dv * dmu * (DvmuEN_F0 - DvmuES_F0);
            aw  += 0.5 * nu * dv * dmu * (-DvmuWN_F0 + DvmuWS_F0);
            an  += 0.5 * nu * dv * dmu * (DvmuWN_F0 - DvmuEN_F0);
            as  += 0.5 * nu * dv * dmu * (-DvmuWS_F0 + DvmuES_F0);
            asw += 0.5 * nu * dv * dmu * (DvmuWS_F0);
            ase += 0.5 * nu * dv * dmu * (-DvmuES_F0);
            anw += 0.5 * nu * dv * dmu * (-DvmuWN_F0);
            ane += 0.5 * nu * dv * dmu * (DvmuEN_F0);
          }
          
          int  ncols = m_nbands, ix = 0;
          int  *icols = (int*)  calloc (ncols,sizeof(int));
          Real *data  = (Real*) calloc (ncols,sizeof(Real));
  
          /* center element */
          icols[ix] = pc; 
          data[ix] = a_shift - ac;
          ix++;
  
          /* east element */
          if (pe >= 0) {
            icols[ix] = pe;
            data[ix] = -ae;
            ix++;
          }
          /* west element */
          if (pw >= 0) {
            icols[ix] = pw;
            data[ix] = -aw;
            ix++;
          }
          /* north element */
          if (pn >= 0) {
            icols[ix] = pn;
            data[ix] = -an;
            ix++;
          }
          /* south element */
          if (ps >= 0) {
            icols[ix] = ps;
            data[ix] = -as;
            ix++;
          }
          /* north east element */
          if (pne >= 0) {
            icols[ix] = pne;
            data[ix] = -ane;
            ix++;
          }
          /* north west element */
          if (pnw >= 0) {
            icols[ix] = pnw;
            data[ix] = -anw;
            ix++;
          }
          /* south east element */
          if (pse >= 0) {
            icols[ix] = pse;
            data[ix] = -ase;
            ix++;
          }
          /* south west element */
          if (psw >= 0) {
            icols[ix] = psw;
            data[ix] = -asw;
            ix++;
          }
  
          CH_assert(ix <= m_nbands);
          CH_assert(ix <= Pmat->getNBands());
          Pmat->setRowValues(pc,ix,icols,data);
          free(data);
          free(icols);
        }
      }
    }
  }
  return;
}

void FokkerPlanck::preTimeStep( const KineticSpeciesPtrVect& a_soln_mapped, 
                                const int a_species, 
                                const Real a_time,
                                const KineticSpeciesPtrVect& a_soln_physical )
{ 
  static bool first_call = true;
  computeReferenceSolution  ( a_soln_physical, a_species, a_time);

  if (first_call) {
    const KineticSpecies& soln_species(*(a_soln_mapped[a_species]));
    const PhaseGeom&      phase_geom(soln_species.phaseSpaceGeometry());
    const CFG::MagGeom&   mag_geom(phase_geom.magGeom());

    m_fp_kinetic_energy.define(mag_geom.gridsFull(), 2, CFG::IntVect::Zero);
    for (CFG::DataIterator dit(m_fp_kinetic_energy.dataIterator()); dit.ok(); ++dit) {
      m_fp_kinetic_energy[dit].setVal(0.0);
    }
    m_energy_cons.define(mag_geom.gridsFull(), 1, CFG::IntVect::Zero);
    for (CFG::DataIterator dit(m_energy_cons.dataIterator()); dit.ok(); ++dit) {
      m_energy_cons[dit].setVal(1.0);
    }
  }
  
  const KineticSpecies&       soln_species(*(a_soln_mapped[a_species]));
  const LevelData<FArrayBox>& soln_dfn(soln_species.distributionFunction());
  const DisjointBoxLayout&    grids(soln_dfn.getBoxes());
  const PhaseGeom&            phase_geom(soln_species.phaseSpaceGeometry());
  const int                   n_comp(soln_dfn.nComp());
    
  LevelData<FArrayBox> dfn(grids,n_comp,IntVect::Zero);
  for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
    dfn[dit].copy(soln_dfn[dit]);
  }
  convertToCellCenters(phase_geom,dfn);

  const RefCountedPtr<PhaseGeom>& phase_geom_ptr(soln_species.phaseSpaceGeometryPtr());

  computePotentialsAndCoeffs( soln_species, dfn, grids, phase_geom_ptr); 

  m_it_counter+=1;
  first_call = false;
  return;
}

void FokkerPlanck::evalRosenbluthPotentials( LevelData<FArrayBox>&            a_phi,
                                             const RefCountedPtr<PhaseGeom>&  a_phase_geom,
                                             const LevelData<FArrayBox>&      a_dfn,
                                             const double                     a_mass ) const 
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
  a_ppcls.query( "time_implicit", m_time_implicit);
  a_ppcls.query( "cls_freq", m_cls_freq );
  a_ppcls.query( "update_frequency", m_update_freq);
  a_ppcls.query( "verbose", m_verbosity);
  a_ppcls.query( "rosenbluth_skip_stage",m_rosenbluth_skip_stage);
  a_ppcls.query( "convergence_tolerance", m_pcg_tol );
  a_ppcls.query( "max_interation_number", m_pcg_maxiter);
  a_ppcls.query( "multipole_number", m_mult_num);

  a_ppcls.query( "conserve_energy",m_conserve_energy);
  a_ppcls.query( "energy_conserve_eps",m_fp_energy_cons_epsilon);
  a_ppcls.query( "use_limiters", m_limiters);

  a_ppcls.query( "subtract_background", m_subtract_background );
  if (m_subtract_background) {
    a_ppcls.query( "compute_maxwellian",m_compute_maxwellian);
  }
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
    std::cout << "  cls_freq            = " << m_cls_freq << std::endl
              << "  subtract_background = " << m_subtract_background    << std::endl
              << "  compute_maxwellian  = " << m_compute_maxwellian     << std::endl
              << "  conserve_energy     = " << m_conserve_energy        << std::endl
              << "  energy_conserve_eps = " << m_fp_energy_cons_epsilon << std::endl
              << "  implicit in time    = " << m_time_implicit          << std::endl
              << "  use_limiters        = " << m_limiters               << std::endl;
    if ((m_subtract_background) && !m_compute_maxwellian) {
      std::cout << "  Reference Function:" << std::endl;
      m_ref_func->printParameters();
    }
  }
  return;
}

void FokkerPlanck::convertToCellCenters(const PhaseGeom&      a_phase_geom,
                                        LevelData<FArrayBox>& a_x) const
{
  const DisjointBoxLayout& grids = a_x.disjointBoxLayout();
  const int n_comp = a_x.nComp();

  LevelData<FArrayBox> y(grids, n_comp, a_x.ghostVect()+IntVect::Unit);

  DataIterator dit(grids.dataIterator());

  for(dit.begin(); dit.ok(); ++dit) {
    y[dit].copy(a_x[dit]);
  }
  y.exchange();

  for (dit.begin(); dit.ok(); ++dit) {
    const PhaseBlockCoordSys& coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
    const ProblemDomain& domain = coord_sys.domain();

    fourthOrderAverageCell(y[dit], domain, grids[dit], -1);
  }
  y.exchange();

  for(dit.begin(); dit.ok(); ++dit) {
    a_x[dit].copy(y[dit]);
  }

  return;
}

void FokkerPlanck::convertToCellAverages( const PhaseGeom&      a_phase_geom,
                                          LevelData<FArrayBox>& a_x) const
{
  const DisjointBoxLayout& grids = a_x.disjointBoxLayout();
  const int n_comp = a_x.nComp();

  LevelData<FArrayBox> y(grids, n_comp, a_x.ghostVect()+IntVect::Unit);

  DataIterator dit(grids.dataIterator());

  for(dit.begin(); dit.ok(); ++dit) {
    y[dit].copy(a_x[dit]);
  }
  y.exchange();

  for (dit.begin(); dit.ok(); ++dit) {
    const PhaseBlockCoordSys& coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
    const ProblemDomain& domain = coord_sys.domain();

    fourthOrderAverageCell(y[dit], domain, grids[dit], 1);
  }
  y.exchange();

  for(dit.begin(); dit.ok(); ++dit) {
    a_x[dit].copy(y[dit]);
  }

  return;
}

void FokkerPlanck::fillGhostCellsPhi( const PhaseGeom& a_phase_geom, LevelData<FArrayBox>& a_x) const
{
  const DisjointBoxLayout& grids = a_x.disjointBoxLayout();
  const IntVect& ghosts(a_x.ghostVect());
  const int nvars(a_x.nComp());
  DataIterator dit(grids.dataIterator());

  for (dit.begin(); dit.ok(); ++dit) {
    const PhaseBlockCoordSys& coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
    const ProblemDomain& domain = coord_sys.domain();

    fourthOrderCellExtrapAtDomainBdry(a_x[dit],domain,grids[dit]);
  }

  a_x.exchange();

  /* fill velocity space corners */
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox& x(a_x[dit]);
    const Box& grid(x.box());

    IntVect bigEnd = grid.bigEnd(), smallEnd = grid.smallEnd();
    for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
      for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
#if CH_SPACEDIM==5
      for (int m = smallEnd[2]; m <bigEnd[2]; m++) {
#endif

        IntVect ic, i1, i2, i3;
        ic[0] = i1[0] = i2[0] = i3[0] = i;
        ic[1] = i1[1] = i2[1] = i3[1] = j;
#if CH_SPACEDIM==5
        ic[2] = i1[2] = i2[2] = i3[2] = m;
#endif

        for (int k=smallEnd[_VPAR_DIM_]+ghosts[_VPAR_DIM_]-1; k >= smallEnd[_VPAR_DIM_]; k--) {
          for (int l=smallEnd[_MU_DIM_]+ghosts[_MU_DIM_]-1; l >= smallEnd[_MU_DIM_]; l--) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k+1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l+1;
            i3[_VPAR_DIM_]  = k+1;
            i3[_MU_DIM_]    = l+1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val1 + val2 - val3;
              x.set(ic,n,val);
            }
          }
        }
      
        for (int k=bigEnd[_VPAR_DIM_]-ghosts[_VPAR_DIM_]+1; k <= bigEnd[_VPAR_DIM_]; k++) {
          for (int l=smallEnd[_MU_DIM_]+ghosts[_MU_DIM_]-1; l >= smallEnd[_MU_DIM_]; l--) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k-1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l+1;
            i3[_VPAR_DIM_]  = k-1;
            i3[_MU_DIM_]    = l+1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val1 + val2 - val3;
              x.set(ic,n,val);
            }
          }
        }
      
        for (int k=smallEnd[_VPAR_DIM_]+ghosts[_VPAR_DIM_]-1; k >= smallEnd[_VPAR_DIM_]; k--) {
          for (int l=bigEnd[_MU_DIM_]-ghosts[_MU_DIM_]+1; l <= bigEnd[_MU_DIM_]; l++) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k+1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l-1;
            i3[_VPAR_DIM_]  = k+1;
            i3[_MU_DIM_]    = l-1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val1 + val2 - val3;
              x.set(ic,n,val);
            }
          }
        }
      
        for (int k=bigEnd[_VPAR_DIM_]-ghosts[_VPAR_DIM_]+1; k <= bigEnd[_VPAR_DIM_]; k++) {
          for (int l=bigEnd[_MU_DIM_]-ghosts[_MU_DIM_]+1; l <= bigEnd[_MU_DIM_]; l++) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k-1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l-1;
            i3[_VPAR_DIM_]  = k-1;
            i3[_MU_DIM_]    = l-1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val1 + val2 - val3;
              x.set(ic,n,val);
            }
          }
        }
      
#if CH_SPACEDIM==5
      }
#endif
      }
    }
  }
  a_x.exchange();
  return;
}

void FokkerPlanck::fillGhostCellsDfn( const PhaseGeom& a_phase_geom, LevelData<FArrayBox>& a_x) const
{
  const DisjointBoxLayout& grids = a_x.disjointBoxLayout();
  //  const IntVect& ghosts(a_x.ghostVect());
  //  const int nvars(a_x.nComp());
  DataIterator dit(grids.dataIterator());

  for (dit.begin(); dit.ok(); ++dit) {
    const PhaseBlockCoordSys& coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
    const ProblemDomain& domain = coord_sys.domain();

    fourthOrderCellExtrapAtDomainBdry(a_x[dit],domain,grids[dit]);
  }

  a_x.exchange();
  return;
}

void FokkerPlanck::fillGhostCellsDfnVPar( const PhaseGeom& a_phase_geom, LevelData<FArrayBox>& a_x) const
{
  const DisjointBoxLayout& grids = a_x.disjointBoxLayout();
  const IntVect& ghosts(a_x.ghostVect());
  const int nvars(a_x.nComp());
  DataIterator dit(grids.dataIterator());

  for (dit.begin(); dit.ok(); ++dit) {
    const PhaseBlockCoordSys& coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
    const ProblemDomain& domain = coord_sys.domain();

    fourthOrderCellExtrapAtDomainBdry(a_x[dit],domain,grids[dit]);
  }

  a_x.exchange();

  /* fill velocity space corners */
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox& x(a_x[dit]);
    const Box& grid(x.box());

    IntVect bigEnd = grid.bigEnd(), smallEnd = grid.smallEnd();
    for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
      for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
#if CH_SPACEDIM==5
      for (int m = smallEnd[2]; m <= bigEnd[2]; m++) {
#endif

        IntVect ic, i1, i2, i3;
        ic[0] = i1[0] = i2[0] = i3[0] = i;
        ic[1] = i1[1] = i2[1] = i3[1] = j;
#if CH_SPACEDIM==5
        ic[2] = i1[2] = i2[2] = i3[2] = m;
#endif

        for (int k=smallEnd[_VPAR_DIM_]+ghosts[_VPAR_DIM_]-1; k >= smallEnd[_VPAR_DIM_]; k--) {
          for (int l=smallEnd[_MU_DIM_]+ghosts[_MU_DIM_]-1; l >= smallEnd[_MU_DIM_]; l--) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k+1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l+1;
            i3[_VPAR_DIM_]  = k+1;
            i3[_MU_DIM_]    = l+1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val1 + val3 - val2;
              x.set(ic,n,val);
            }
          }
        }
      
        for (int k=bigEnd[_VPAR_DIM_]-ghosts[_VPAR_DIM_]+1; k <= bigEnd[_VPAR_DIM_]; k++) {
          for (int l=smallEnd[_MU_DIM_]+ghosts[_MU_DIM_]-1; l >= smallEnd[_MU_DIM_]; l--) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k-1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l+1;
            i3[_VPAR_DIM_]  = k-1;
            i3[_MU_DIM_]    = l+1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val1 + val3 - val2;
              x.set(ic,n,val);
            }
          }
        }
      
        for (int k=smallEnd[_VPAR_DIM_]+ghosts[_VPAR_DIM_]-1; k >= smallEnd[_VPAR_DIM_]; k--) {
          for (int l=bigEnd[_MU_DIM_]-ghosts[_MU_DIM_]+1; l <= bigEnd[_MU_DIM_]; l++) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k+1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l-1;
            i3[_VPAR_DIM_]  = k+1;
            i3[_MU_DIM_]    = l-1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val1 + val3 - val2;
              x.set(ic,n,val);
            }
          }
        }
      
        for (int k=bigEnd[_VPAR_DIM_]-ghosts[_VPAR_DIM_]+1; k <= bigEnd[_VPAR_DIM_]; k++) {
          for (int l=bigEnd[_MU_DIM_]-ghosts[_MU_DIM_]+1; l <= bigEnd[_MU_DIM_]; l++) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k-1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l-1;
            i3[_VPAR_DIM_]  = k-1;
            i3[_MU_DIM_]    = l-1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val1 + val3 - val2;
              x.set(ic,n,val);
            }
          }
        }
     
#if CH_SPACEDIM==5
      }
#endif
      }
    }
  }
  a_x.exchange();
  return;
}

void FokkerPlanck::fillGhostCellsDfnMu( const PhaseGeom& a_phase_geom, LevelData<FArrayBox>& a_x) const
{
  const DisjointBoxLayout& grids = a_x.disjointBoxLayout();
  const IntVect& ghosts(a_x.ghostVect());
  const int nvars(a_x.nComp());
  DataIterator dit(grids.dataIterator());

  for (dit.begin(); dit.ok(); ++dit) {
    const PhaseBlockCoordSys& coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
    const ProblemDomain& domain = coord_sys.domain();

    fourthOrderCellExtrapAtDomainBdry(a_x[dit],domain,grids[dit]);
  }

  a_x.exchange();

  /* fill velocity space corners */
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox& x(a_x[dit]);
    const Box& grid(x.box());

    IntVect bigEnd = grid.bigEnd(), smallEnd = grid.smallEnd();
    for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
      for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
#if CH_SPACEDIM==5
      for (int m = smallEnd[2]; m <= bigEnd[2]; m++) {
#endif

        IntVect ic, i1, i2, i3;
        ic[0] = i1[0] = i2[0] = i3[0] = i;
        ic[1] = i1[1] = i2[1] = i3[1] = j;
#if CH_SPACEDIM==5
        ic[2] = i1[2] = i2[2] = i3[2] = m;
#endif

        for (int k=smallEnd[_VPAR_DIM_]+ghosts[_VPAR_DIM_]-1; k >= smallEnd[_VPAR_DIM_]; k--) {
          for (int l=smallEnd[_MU_DIM_]+ghosts[_MU_DIM_]-1; l >= smallEnd[_MU_DIM_]; l--) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k+1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l+1;
            i3[_VPAR_DIM_]  = k+1;
            i3[_MU_DIM_]    = l+1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val2 + val3 - val1;
              x.set(ic,n,val);
            }
          }
        }
      
        for (int k=bigEnd[_VPAR_DIM_]-ghosts[_VPAR_DIM_]+1; k <= bigEnd[_VPAR_DIM_]; k++) {
          for (int l=smallEnd[_MU_DIM_]+ghosts[_MU_DIM_]-1; l >= smallEnd[_MU_DIM_]; l--) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k-1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l+1;
            i3[_VPAR_DIM_]  = k-1;
            i3[_MU_DIM_]    = l+1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val2 + val3 - val1;
              x.set(ic,n,val);
            }
          }
        }
      
        for (int k=smallEnd[_VPAR_DIM_]+ghosts[_VPAR_DIM_]-1; k >= smallEnd[_VPAR_DIM_]; k--) {
          for (int l=bigEnd[_MU_DIM_]-ghosts[_MU_DIM_]+1; l <= bigEnd[_MU_DIM_]; l++) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k+1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l-1;
            i3[_VPAR_DIM_]  = k+1;
            i3[_MU_DIM_]    = l-1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val2 + val3 - val1;
              x.set(ic,n,val);
            }
          }
        }
      
        for (int k=bigEnd[_VPAR_DIM_]-ghosts[_VPAR_DIM_]+1; k <= bigEnd[_VPAR_DIM_]; k++) {
          for (int l=bigEnd[_MU_DIM_]-ghosts[_MU_DIM_]+1; l <= bigEnd[_MU_DIM_]; l++) {
            ic[_VPAR_DIM_]  = k;
            ic[_MU_DIM_]    = l;
            i1[_VPAR_DIM_]  = k-1;
            i1[_MU_DIM_]    = l;
            i2[_VPAR_DIM_]  = k;
            i2[_MU_DIM_]    = l-1;
            i3[_VPAR_DIM_]  = k-1;
            i3[_MU_DIM_]    = l-1;

            for (int n = 0; n < nvars; n++) {
              Real val, val1, val2, val3;
              val1  = x.get(i1,n);
              val2  = x.get(i2,n);
              val3  = x.get(i3,n);
              val   = val2 + val3 - val1;
              x.set(ic,n,val);
            }
          }
        }
     
#if CH_SPACEDIM==5
      }
#endif
      }
    }
  }
  a_x.exchange();
  return;
}

void FokkerPlanck::computeReferenceSolution(const KineticSpeciesPtrVect& a_soln,
                                            const int a_species,
                                            const Real a_time)
{
  static bool first_call = true;
  const KineticSpecies& soln_species( *(a_soln[a_species]) );
  const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();

  if (first_call) {
    m_F0.define( phase_geom.gridsFull(), 1, IntVect::Zero );
  }

  if (m_compute_maxwellian) {

    /* Compute the Maxwellian and use that as the reference distribution F0 */
    const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
  
    CFG::LevelData<CFG::FArrayBox> density( mag_geom.grids(), 1, CFG::IntVect::Zero );
    soln_species.numberDensity( density );
  
    CFG::LevelData<CFG::FArrayBox> ParallelMom( mag_geom.grids(), 1, CFG::IntVect::Zero );
    soln_species.ParallelMomentum( ParallelMom );
  
    for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      ParallelMom[dit].divide(density[dit]);
    }
    CFG::LevelData<CFG::FArrayBox> pressure( mag_geom.grids(), 1, CFG::IntVect::Zero );
    soln_species.pressureMoment(pressure, ParallelMom);
  
    for (CFG::DataIterator dit(density.dataIterator()); dit.ok(); ++dit) {
      pressure[dit].divide(density[dit]);
    }
  
    MaxwellianKernel<FArrayBox> maxwellian(density,pressure,ParallelMom);
    maxwellian.eval(m_F0,soln_species);
    phase_geom.multJonValid(m_F0);
    convertToCellCenters(phase_geom, m_F0);

    m_ref_updated = true;

  } else {

    /* use the specified reference distribution function */
    if (first_call) {

      KineticSpeciesPtr ref_species( soln_species.clone( IntVect::Unit, false ) );
      m_ref_func->assign( *ref_species, a_time );
      LevelData<FArrayBox>& ref_dfn( ref_species->distributionFunction() );
      for (DataIterator dit(m_F0.dataIterator()); dit.ok(); ++dit) {
        m_F0[dit].copy(ref_dfn[dit]);
      }
      convertToCellCenters(phase_geom, m_F0);
      m_ref_updated = true;

    } else {

      m_ref_updated = false;

    }

  }

  first_call = false;
  return;
}

void FokkerPlanck::evalCoefficients(  FokkerPlanckCoeffs&         a_coeffs,
                                      const LevelData<FArrayBox>& a_phi,
                                      const PhaseGeom&            a_phase_geom,
                                      const Real                  a_mass_tp,
                                      const Real                  a_mass_fp ) const
{
  const VEL::VelCoordSys&   vel_coords  = a_phase_geom.velSpaceCoordSys();
  const VEL::ProblemDomain& vel_domain  = vel_coords.domain();
  const VEL::Box&           domain_box  = vel_domain.domainBox();
  const VEL::RealVect&      vel_dx      = vel_coords.dx();

  const int num_vpar_cells  = domain_box.size(0);
  const int num_mu_cells    = domain_box.size(1);

  LevelData<FArrayBox> phi_tmp( a_phi.disjointBoxLayout(),
                                a_phi.nComp(),
                                a_phi.ghostVect()+4*IntVect::Unit );

  const DisjointBoxLayout& grids(a_phi.getBoxes());
  for (DataIterator dit(a_phi.dataIterator()); dit.ok(); ++dit) {
    phi_tmp[dit].setVal(0.0);
    phi_tmp[dit].copy(a_phi[dit],grids[dit]);
  }
  /* divide by (J*B_star_parallel) because the distribution function from
   * which phi is computed is multiplied by this factor */
  a_phase_geom.divideJonValid(phi_tmp);
  a_phase_geom.divideBStarParallel(phi_tmp);
  fillGhostCellsPhi(a_phase_geom, phi_tmp);

  const LevelData<FArrayBox>& injected_B = a_phase_geom.getBFieldMagnitude();

  /* Compute D_v */
  for (DataIterator dit(a_coeffs.D_v.dataIterator()); dit.ok(); ++dit) {

    const FArrayBox&  this_phi  = phi_tmp[dit];
    const FArrayBox&  this_b    = injected_B[dit];
    FArrayBox&        Dv        = a_coeffs.D_v[dit].data();

    FORT_EVAL_DV(  CHF_FRA1(Dv,0),
                   CHF_CONST_FRA(this_phi),
                   CHF_CONST_FRA1(this_b,0),
                   CHF_BOX(Dv.box()),
                   CHF_CONST_REALVECT(vel_dx),
                   CHF_CONST_REAL(a_mass_tp),
                   CHF_CONST_REAL(a_mass_fp),
                   CHF_CONST_INT(num_vpar_cells),
                   CHF_CONST_INT(num_mu_cells));
  }

  /* Compute D_mu */
  for (DataIterator dit(a_coeffs.D_mu.dataIterator()); dit.ok(); ++dit) {

    const FArrayBox&  this_phi  = phi_tmp[dit];
    const FArrayBox&  this_b    = injected_B[dit];
    FArrayBox&        Dmu       = a_coeffs.D_mu[dit].data();

    FORT_EVAL_DMU( CHF_FRA1(Dmu,0),
                   CHF_CONST_FRA(this_phi),
                   CHF_CONST_FRA1(this_b,0),
                   CHF_BOX(Dmu.box()),
                   CHF_CONST_REALVECT(vel_dx),
                   CHF_CONST_REAL(a_mass_tp),
                   CHF_CONST_REAL(a_mass_fp),
                   CHF_CONST_INT(num_vpar_cells),
                   CHF_CONST_INT(num_mu_cells));
  }

  /* Compute D_vv */
  for (DataIterator dit(a_coeffs.D_vv.dataIterator()); dit.ok(); ++dit) {

    const FArrayBox&  this_phi  = phi_tmp[dit];
    const FArrayBox&  this_b    = injected_B[dit];
    FArrayBox&        Dvv       = a_coeffs.D_vv[dit].data();

    FORT_EVAL_DVV( CHF_FRA1(Dvv,0),
                   CHF_CONST_FRA(this_phi),
                   CHF_CONST_FRA1(this_b,0),
                   CHF_BOX(Dvv.box()),
                   CHF_CONST_REALVECT(vel_dx),
                   CHF_CONST_REAL(a_mass_tp),
                   CHF_CONST_REAL(a_mass_fp),
                   CHF_CONST_INT(num_vpar_cells),
                   CHF_CONST_INT(num_mu_cells));
  }

  /* Compute D_mumu */
  for (DataIterator dit(a_coeffs.D_mumu.dataIterator()); dit.ok(); ++dit) {

    const FArrayBox&  this_phi  = phi_tmp[dit];
    const FArrayBox&  this_b    = injected_B[dit];
    FArrayBox&        Dmumu     = a_coeffs.D_mumu[dit].data();

    FORT_EVAL_DMUMU( CHF_FRA1(Dmumu,0),
                     CHF_CONST_FRA(this_phi),
                     CHF_CONST_FRA1(this_b,0),
                     CHF_BOX(Dmumu.box()),
                     CHF_CONST_REALVECT(vel_dx),
                     CHF_CONST_REAL(a_mass_tp),
                     CHF_CONST_REAL(a_mass_fp),
                     CHF_CONST_INT(num_vpar_cells),
                     CHF_CONST_INT(num_mu_cells));
  }

  /* Compute D_vmu */
  for (DataIterator dit(a_coeffs.D_vmu.dataIterator()); dit.ok(); ++dit) {

    const FArrayBox&  this_phi  = phi_tmp[dit];
    const FArrayBox&  this_b    = injected_B[dit];
    FArrayBox&        Dvmu      = a_coeffs.D_vmu[dit].data();

    Box tempBox(Dvmu.box());
    tempBox.growHi(0,-1);
    tempBox.growHi(1,-1);

    FORT_EVAL_DVMU( CHF_FRA1(Dvmu,0),
                    CHF_CONST_FRA(this_phi),
                    CHF_CONST_FRA1(this_b,0),
                    CHF_BOX(tempBox),
                    CHF_CONST_REALVECT(vel_dx),
                    CHF_CONST_REAL(a_mass_tp),
                    CHF_CONST_REAL(a_mass_fp),
                    CHF_CONST_INT(num_vpar_cells),
                    CHF_CONST_INT(num_mu_cells));
  }

  return;
}

void FokkerPlanck::computePotentialsAndCoeffs(const KineticSpecies&            a_species, 
                                              const LevelData<FArrayBox>&      a_dfn,
                                              const DisjointBoxLayout&         a_grids,
                                              const RefCountedPtr<PhaseGeom>&  a_phase_geom )
{
  static bool first_call = true;

  //Get test-particle (tp) and field-particle (fp) masses
  const double mass_tp = a_species.mass();
  const double mass_fp = a_species.mass();

  if (first_call) {
    //Compute normalization
    computeClsNorm(m_cls_norm, a_species.mass(), a_species.charge());

    // Define the data structures for potentials
    m_phi.define(a_grids, 2, IntVect::Zero);
    for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
      m_phi[dit].setVal(0.0);
    }

    // Define the data structures for coefficients
    m_coeffs.D_v.define   (a_grids, 1, IntVect::Zero);
    m_coeffs.D_mu.define  (a_grids, 1, IntVect::Zero);
    m_coeffs.D_vv.define  (a_grids, 1, IntVect::Zero);
    m_coeffs.D_mumu.define(a_grids, 1, IntVect::Zero);
    m_coeffs.D_vmu.define (a_grids, 1, IntVect::Zero);
    for (DataIterator dit(a_grids.dataIterator()); dit.ok(); ++dit) {
      m_coeffs.D_v[dit].reDir(_VPAR_DIM_);
      m_coeffs.D_mu[dit].reDir(_MU_DIM_);
      m_coeffs.D_vv[dit].reDir(_VPAR_DIM_);
      m_coeffs.D_mumu[dit].reDir(_MU_DIM_);
      m_coeffs.D_v[dit].setVal(0.0);
      m_coeffs.D_mu[dit].setVal(0.0);
      m_coeffs.D_vv[dit].setVal(0.0);
      m_coeffs.D_mumu[dit].setVal(0.0);
      m_coeffs.D_vmu[dit].setVal(0.0);
    }
    
    if (m_subtract_background) {
      m_phi_F0.define(a_grids, 2, IntVect::Zero);
      for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
        m_phi_F0[dit].setVal(0.0);
      }
    
      m_coeffs_F0.D_v.define   (a_grids, 1, IntVect::Zero);
      m_coeffs_F0.D_mu.define  (a_grids, 1, IntVect::Zero);
      m_coeffs_F0.D_vv.define  (a_grids, 1, IntVect::Zero);
      m_coeffs_F0.D_mumu.define(a_grids, 1, IntVect::Zero);
      m_coeffs_F0.D_vmu.define (a_grids, 1, IntVect::Zero);
      for (DataIterator dit(a_grids.dataIterator()); dit.ok(); ++dit) {
        m_coeffs_F0.D_v[dit].reDir(_VPAR_DIM_);
        m_coeffs_F0.D_mu[dit].reDir(_MU_DIM_);
        m_coeffs_F0.D_vv[dit].reDir(_VPAR_DIM_);
        m_coeffs_F0.D_mumu[dit].reDir(_MU_DIM_);
        m_coeffs_F0.D_v[dit].setVal(0.0);
        m_coeffs_F0.D_mu[dit].setVal(0.0);
        m_coeffs_F0.D_vv[dit].setVal(0.0);
        m_coeffs_F0.D_mumu[dit].setVal(0.0);
        m_coeffs_F0.D_vmu[dit].setVal(0.0);
      }
    }

  }

  if ((m_update_freq < 0) || (m_it_counter % m_update_freq == 0) || first_call) {

    if (m_subtract_background) {

      LevelData<FArrayBox> delta_f(a_grids, a_dfn.nComp(), IntVect::Zero);
      for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
        delta_f[dit].copy(a_dfn[dit]);
        delta_f[dit].minus(m_F0[dit]);
      }

      if (m_ref_updated || first_call) {
        if (!procID() && m_verbosity) {
          std::cout<<"  Note: Reference distribution function updated (by recomputing Maxwellian).\n";
          std::cout<<"  Note: Recomputing Rosenbluth potentials for reference distribution function.\n";
        }
        evalRosenbluthPotentials(m_phi_F0, a_phase_geom, m_F0, mass_tp);
        evalCoefficients(m_coeffs_F0, m_phi_F0, *a_phase_geom, mass_tp, mass_fp);
        m_ref_updated = false;
      }

      evalRosenbluthPotentials(m_phi, a_phase_geom, delta_f, mass_tp);
      evalCoefficients(m_coeffs, m_phi, *a_phase_geom, mass_tp, mass_fp);

    } else {

      evalRosenbluthPotentials(m_phi, a_phase_geom, a_dfn, mass_tp);
      evalCoefficients(m_coeffs, m_phi, *a_phase_geom, mass_tp, mass_fp);

    }

  }

  first_call = false;
  return;
}

void FokkerPlanck::computeFokkerPlanckFlux( LevelData<FaceBox>&         a_flux_vpar,
                                            LevelData<FaceBox>&         a_flux_mu,
                                            const FokkerPlanckCoeffs&   a_coeffs,
                                            const LevelData<FArrayBox>& a_dfn,
                                            const PhaseGeom&            a_phase_geom )
{
  const VEL::VelCoordSys&   vel_coords  = a_phase_geom.velSpaceCoordSys();
  const VEL::ProblemDomain& vel_domain  = vel_coords.domain();
  const VEL::Box&           domain_box  = vel_domain.domainBox();
  const VEL::RealVect&      vel_dx      = vel_coords.dx();

  const int num_vpar_cells  = domain_box.size(0);
  const int num_mu_cells    = domain_box.size(1);

  LevelData<FArrayBox> dfn_tmp( a_dfn.disjointBoxLayout(),
                                a_dfn.nComp(),
                                a_dfn.ghostVect()+4*IntVect::Unit );
  const DisjointBoxLayout& grids(a_dfn.disjointBoxLayout());

  for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
    dfn_tmp[dit].setVal(0.0);
    dfn_tmp[dit].copy(a_dfn[dit],grids[dit]);
  }
  dfn_tmp.exchange();
  fillGhostCellsDfn(a_phase_geom, dfn_tmp);

  LevelData<FaceBox> dfn_face_vpar, dfn_face_mu;
  dfn_face_vpar.define(grids,1,IntVect::Zero);
  dfn_face_mu.define  (grids,1,IntVect::Zero);
  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    dfn_face_vpar[dit].reDir(_VPAR_DIM_);
    dfn_face_mu[dit].reDir(_MU_DIM_);
    dfn_face_vpar[dit].setVal(0.0);
    dfn_face_mu[dit].setVal(0.0);
  }

  int flag = (m_limiters ? 1 : 0);

  /* reconstruct dfn at the faces */
  for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {

    const FArrayBox&  dfn       = dfn_tmp[dit];
    const FArrayBox&  Dv        = a_coeffs.D_v[dit].data();
    const FArrayBox&  Dmu       = a_coeffs.D_mu[dit].data();
    FArrayBox&        dfn_vpar  = dfn_face_vpar[dit].data();
    FArrayBox&        dfn_mu    = dfn_face_mu[dit].data();

    FORT_RECONSTRUCT_DFN_FACE_VPAR( CHF_FRA1(dfn_vpar,0),
                                    CHF_CONST_FRA1(Dv,0),
                                    CHF_CONST_FRA1(dfn,0),
                                    CHF_BOX(dfn_vpar.box()),
                                    CHF_CONST_INT(flag),
                                    CHF_CONST_INT(num_vpar_cells));

    FORT_RECONSTRUCT_DFN_FACE_MU  ( CHF_FRA1(dfn_mu,0),
                                    CHF_CONST_FRA1(Dmu,0),
                                    CHF_CONST_FRA1(dfn,0),
                                    CHF_BOX(dfn_mu.box()),
                                    CHF_CONST_INT(flag),
                                    CHF_CONST_INT(num_mu_cells));
  }

  for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
    dfn_tmp[dit].setVal(0.0);
    dfn_tmp[dit].copy(a_dfn[dit],grids[dit]);
  }
  dfn_tmp.exchange();
  fillGhostCellsDfnMu(a_phase_geom, dfn_tmp);

  /* Compute flux along vpar */
  for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {

    const FArrayBox&  dfn       = dfn_tmp[dit];
    const FArrayBox&  Dv        = a_coeffs.D_v[dit].data();
    const FArrayBox&  Dvv       = a_coeffs.D_vv[dit].data();
    const FArrayBox&  Dvmu      = a_coeffs.D_vmu[dit].data();
    const FArrayBox&  dfn_face  = dfn_face_vpar[dit].data();
    FArrayBox&        flux      = a_flux_vpar[dit].data();

    FORT_FLUX_VPAR(  CHF_FRA(flux),
                     CHF_CONST_FRA1(Dv,0),
                     CHF_CONST_FRA1(Dvv,0),
                     CHF_CONST_FRA1(Dvmu,0),
                     CHF_CONST_FRA(dfn),
                     CHF_CONST_FRA(dfn_face),
                     CHF_BOX(flux.box()),
                     CHF_CONST_REALVECT(vel_dx),
                     CHF_CONST_INT(num_vpar_cells),
                     CHF_CONST_INT(num_mu_cells),
                     CHF_CONST_INT(flag));
  }

  for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
    dfn_tmp[dit].setVal(0.0);
    dfn_tmp[dit].copy(a_dfn[dit],grids[dit]);
  }
  dfn_tmp.exchange();
  fillGhostCellsDfnVPar(a_phase_geom, dfn_tmp);

  /* Compute flux along mu */
  for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {

    const FArrayBox&  dfn       = dfn_tmp[dit];
    const FArrayBox&  Dmu       = a_coeffs.D_mu[dit].data();
    const FArrayBox&  Dmumu     = a_coeffs.D_mumu[dit].data();
    const FArrayBox&  Dvmu      = a_coeffs.D_vmu[dit].data();
    const FArrayBox&  dfn_face  = dfn_face_mu[dit].data();
    FArrayBox&        flux      = a_flux_mu[dit].data();

    FORT_FLUX_MU(  CHF_FRA(flux),
                   CHF_CONST_FRA1(Dmu,0),
                   CHF_CONST_FRA1(Dmumu,0),
                   CHF_CONST_FRA1(Dvmu,0),
                   CHF_CONST_FRA(dfn),
                   CHF_CONST_FRA(dfn_face),
                   CHF_BOX(flux.box()),
                   CHF_CONST_REALVECT(vel_dx),
                   CHF_CONST_INT(num_vpar_cells),
                   CHF_CONST_INT(num_mu_cells),
                   CHF_CONST_INT(flag));
  }

  return;
}

void FokkerPlanck::computeEnergyConservationFactor( const KineticSpecies&       a_species,
                                                    const LevelData<FaceBox>&   a_flux_vpar,
                                                    const LevelData<FaceBox>&   a_flux_mu,
                                                    const DisjointBoxLayout&    a_grids,
                                                    const PhaseGeom&            a_phase_geom )
{
  const VEL::VelCoordSys&   vel_coords  = a_phase_geom.velSpaceCoordSys();
  const VEL::RealVect&      vel_dx      = vel_coords.dx();

  if (m_conserve_energy) {
    
    LevelData<FArrayBox> untamed_rhs(a_grids,2,IntVect::Zero);
    for (DataIterator dit(untamed_rhs.dataIterator()); dit.ok(); ++dit) {
      FArrayBox&        rhs       = untamed_rhs[dit];
      const FArrayBox&  flux_vpar = a_flux_vpar[dit].data();
      const FArrayBox&  flux_mu   = a_flux_mu[dit].data();

      
      FORT_EVAL_RHS_COMPONENTS( CHF_FRA(rhs),
                                CHF_CONST_FRA(flux_vpar),
                                CHF_CONST_FRA(flux_mu),
                                CHF_BOX(rhs.box()),
                                CHF_CONST_REALVECT(vel_dx) );
    }
    
    a_species.energyMoment(m_fp_kinetic_energy, untamed_rhs);
    
    for (CFG::DataIterator dit(m_energy_cons.dataIterator()); dit.ok(); ++dit) {
      m_fp_kinetic_energy[dit].abs();
      m_fp_kinetic_energy[dit].plus(m_fp_energy_cons_epsilon);
      m_energy_cons[dit].copy(m_fp_kinetic_energy[dit],0,0);
      m_energy_cons[dit].divide(m_fp_kinetic_energy[dit],1,0);
    }

  } else {

    for (CFG::DataIterator dit(m_energy_cons.dataIterator()); dit.ok(); ++dit) {
      m_energy_cons[dit].setVal(1.0);
    }

  }

}

void FokkerPlanck::computeFokkerPlanckRHS(  LevelData<FArrayBox>&       a_rhs,
                                            const LevelData<FaceBox>&   a_flux_vpar,
                                            const LevelData<FaceBox>&   a_flux_mu,
                                            const PhaseGeom&            a_phase_geom )
{
  const VEL::VelCoordSys&   vel_coords  = a_phase_geom.velSpaceCoordSys();
  const VEL::ProblemDomain& vel_domain  = vel_coords.domain();
  const VEL::Box&           domain_box  = vel_domain.domainBox();
  const VEL::RealVect&      vel_dx      = vel_coords.dx();

  const int num_vpar_cells  = domain_box.size(0);
  const int num_mu_cells    = domain_box.size(1);

  LevelData<FArrayBox> energy_cons_factor;
  a_phase_geom.injectConfigurationToPhase(m_energy_cons, energy_cons_factor);

  for (DataIterator dit(a_rhs.dataIterator()); dit.ok(); ++dit) {

    FArrayBox&        rhs       = a_rhs[dit];
    const FArrayBox&  flux_vpar = a_flux_vpar[dit].data();
    const FArrayBox&  flux_mu   = a_flux_mu[dit].data();
    const FArrayBox&  energy    = energy_cons_factor[dit];
    const Real        nu        = ( m_fixed_cls_freq ? m_cls_freq : m_cls_norm );

    FORT_EVAL_RHS(  CHF_FRA1(rhs,0),
                    CHF_CONST_FRA(flux_vpar),
                    CHF_CONST_FRA(flux_mu),
                    CHF_CONST_FRA1(energy,0),
                    CHF_BOX(rhs.box()),
                    CHF_CONST_REALVECT(vel_dx),
                    CHF_CONST_INT(num_vpar_cells),
                    CHF_CONST_INT(num_mu_cells),
                    CHF_CONST_REAL(nu) );
  }

  return;
}

void FokkerPlanck::evalClsRHS(  KineticSpeciesPtrVect&        a_rhs,
                                const KineticSpeciesPtrVect&  a_soln,
                                const int                     a_species,
                                const int                     a_species_bkgr,
                                const Real                    a_time )
{
  const KineticSpecies&           soln_species(*(a_soln[a_species]));
  const LevelData<FArrayBox>&     soln_dfn(soln_species.distributionFunction());
  const DisjointBoxLayout&        grids(soln_dfn.getBoxes());
  const PhaseGeom&                phase_geom(soln_species.phaseSpaceGeometry());
  const RefCountedPtr<PhaseGeom>& phase_geom_ptr(soln_species.phaseSpaceGeometryPtr());
  const int                       n_comp(soln_dfn.nComp());

  LevelData<FArrayBox> dfn(grids,n_comp,IntVect::Zero);
  for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
    dfn[dit].copy(soln_dfn[dit]);
  }
  convertToCellCenters(phase_geom,dfn);

  if (!m_rosenbluth_skip_stage) {
    computePotentialsAndCoeffs(soln_species, dfn, grids, phase_geom_ptr);
  }

  LevelData<FaceBox> fp_flux_vpar, fp_flux_mu;
  fp_flux_vpar.define (grids,2,IntVect::Zero);
  fp_flux_mu.define   (grids,2,IntVect::Zero);
  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    fp_flux_vpar[dit].reDir(_VPAR_DIM_);
    fp_flux_mu[dit].reDir(_MU_DIM_);
    fp_flux_vpar[dit].setVal(0.0);
    fp_flux_mu[dit].setVal(0.0);
  }

  if (m_subtract_background) {

    LevelData<FArrayBox> delta_f(grids, dfn.nComp(), IntVect::Zero);
    for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      delta_f[dit].copy(dfn[dit]);
      delta_f[dit].minus(m_F0[dit]);
    }

    LevelData<FaceBox> tmp_flux_vpar, tmp_flux_mu;
    tmp_flux_vpar.define (grids,2,IntVect::Zero);
    tmp_flux_mu.define   (grids,2,IntVect::Zero);
    for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      tmp_flux_vpar[dit].reDir(_VPAR_DIM_);
      tmp_flux_mu[dit].reDir(_MU_DIM_);
      tmp_flux_vpar[dit].setVal(0.0);
      tmp_flux_mu[dit].setVal(0.0);
    }

    computeFokkerPlanckFlux(tmp_flux_vpar, tmp_flux_mu, m_coeffs_F0, delta_f, phase_geom);
    for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      fp_flux_vpar[dit].plus(tmp_flux_vpar[dit]);
      fp_flux_mu[dit].plus(tmp_flux_mu[dit]);
    }

    computeFokkerPlanckFlux(tmp_flux_vpar, tmp_flux_mu, m_coeffs, m_F0, phase_geom);
    for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      fp_flux_vpar[dit].plus(tmp_flux_vpar[dit]);
      fp_flux_mu[dit].plus(tmp_flux_mu[dit]);
    }

    computeFokkerPlanckFlux(tmp_flux_vpar, tmp_flux_mu, m_coeffs, delta_f, phase_geom);
    for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
      fp_flux_vpar[dit].plus(tmp_flux_vpar[dit]);
      fp_flux_mu[dit].plus(tmp_flux_mu[dit]);
    }

  } else {
    
    computeFokkerPlanckFlux(fp_flux_vpar, fp_flux_mu, m_coeffs, dfn, phase_geom);

  }

  /* compute energy conservation fix */
  computeEnergyConservationFactor(soln_species, fp_flux_vpar, fp_flux_mu, grids, phase_geom);

  /* reusing dfn to store the RHS now - save memory */
  for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
    dfn[dit].setVal(0.0);
  }
  computeFokkerPlanckRHS(dfn, fp_flux_vpar, fp_flux_mu, phase_geom);
  convertToCellAverages(phase_geom,dfn);

  KineticSpecies& rhs_species(*(a_rhs[a_species]));
  LevelData<FArrayBox>& rhs(rhs_species.distributionFunction());
  for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
    rhs[dit].plus(dfn[dit]);
  }

  return;
}

Real FokkerPlanck::computeTimeScale(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  const KineticSpecies&       soln_species(*(a_soln[a_species]));
  const LevelData<FArrayBox>& dfn(soln_species.distributionFunction());
  const PhaseGeom&            phase_geom(soln_species.phaseSpaceGeometry());
  const VEL::VelCoordSys&     vel_coords(phase_geom.velSpaceCoordSys());
  const VEL::RealVect&        vel_dx(vel_coords.dx());

  Real dvpar = vel_dx[0], dmu = vel_dx[1];
  Real tau = DBL_MAX;

  for (DataIterator dit(dfn.dataIterator()); dit.ok(); ++dit) {
    const FArrayBox& Dv     = m_coeffs.D_v[dit].data();
    const FArrayBox& Dmu    = m_coeffs.D_mu[dit].data();
    const FArrayBox& Dvv    = m_coeffs.D_vv[dit].data();
    const FArrayBox& Dmumu  = m_coeffs.D_mumu[dit].data();
    const FArrayBox& Dvmu   = m_coeffs.D_vmu[dit].data();

    Real Dv_max     = std::max(abs(Dv.max()),abs(Dv.min()));
    Real Dmu_max    = std::max(abs(Dmu.max()),abs(Dmu.min()));
    Real Dvv_max    = std::max(abs(Dvv.max()),abs(Dvv.min()));
    Real Dmumu_max  = std::max(abs(Dmumu.max()),abs(Dmumu.min()));
    Real Dvmu_max   = std::max(abs(Dvmu.max()),abs(Dvmu.min()));

    if (m_subtract_background) {
      const FArrayBox& Dv     = m_coeffs_F0.D_v[dit].data();
      const FArrayBox& Dmu    = m_coeffs_F0.D_mu[dit].data();
      const FArrayBox& Dvv    = m_coeffs_F0.D_vv[dit].data();
      const FArrayBox& Dmumu  = m_coeffs_F0.D_mumu[dit].data();
      const FArrayBox& Dvmu   = m_coeffs_F0.D_vmu[dit].data();

      Dv_max    += std::max(abs(Dv.max()),abs(Dv.min()));
      Dmu_max   += std::max(abs(Dmu.max()),abs(Dmu.min()));
      Dvv_max   += std::max(abs(Dvv.max()),abs(Dvv.min()));
      Dmumu_max += std::max(abs(Dmumu.max()),abs(Dmumu.min()));
      Dvmu_max  += std::max(abs(Dvmu.max()),abs(Dvmu.min()));
    }

    if (dvpar/Dv_max < tau)           tau = dvpar/Dv_max;
    if (0.5*dmu/Dmu_max < tau)        tau = 0.5*dmu/Dmu_max;
    if (dvpar*dvpar/Dvv_max < tau)    tau = dvpar*dvpar/Dvv_max;
    if (0.25*dmu*dmu/Dmumu_max < tau) tau = 0.25*dmu*dmu/Dmumu_max;
    if (0.5*dvpar*dmu/Dvmu_max < tau) tau = 0.5*dvpar*dmu/Dvmu_max;
  }

    const Real nu = ( m_fixed_cls_freq ? m_cls_freq : m_cls_norm );
    tau /= nu;

  Real global_tau = 0.0;
#ifdef CH_MPI
  MPI_Allreduce(&tau,&global_tau,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#else
  global_tau = tau;
#endif

  return global_tau;
}

Real FokkerPlanck::computeDtExplicitTI(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  return computeTimeScale(a_soln, a_species);
}

Real FokkerPlanck::computeDtImExTI(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  if (m_time_implicit) {
    return DBL_MAX;
  } else {
    return computeTimeScale(a_soln, a_species);
  }
}

#include "NamespaceFooter.H"
