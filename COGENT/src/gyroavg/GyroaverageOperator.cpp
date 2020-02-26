#include "GyroaverageOperator.H"
#include "Directions.H"
#include "PhaseGeom.H"
#include "SpaceUtils.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "SpaceUtils.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

static bool closerPoint(  const pair<CFG::IntVect, Real>& a_i,
                          const pair<CFG::IntVect, Real>& a_j )
{
  return (a_i.second < a_j.second);
}

static bool isColinear( const CFG::IntVect&               a_point,
                        const std::vector<CFG::IntVect>&  a_set )
{
  if (a_set.size() < 2) { 
    return false;
  } else {
    bool retval = false;
    for (int d=0; d<CFG_DIM; d++) {
      int count = 0;
#if CFG_DIM==2
      for (int i=0; i<a_set.size(); i++) {
        if (a_set[i][d] == a_point[d]) count++;
      }
#else
      int d1 = d+1; if (d1 >= CFG_DIM) d1 -= CFG_DIM;
      int d2 = d+2; if (d2 >= CFG_DIM) d2 -= CFG_DIM;
      for (int i=0; i<a_set.size(); i++) {
        if ((a_set[i][d1] == a_point[d1]) && (a_set[i][d2] == a_point[d2])) count++;
      }
#endif
      if (count >= 2) retval = true;
    }
    return retval;
  }
}

#if CFG_DIM==3
static bool isCoplanar( const CFG::IntVect&               a_point,
                        const std::vector<CFG::IntVect>&  a_set )
{
  if (a_set.size() < 4) { 
    return false;
  } else {
    bool retval = false;
    for (int d=0; d<CFG_DIM; d++) {
      int count = 0;
      for (int i=0; i<a_set.size(); i++) {
        if (a_set[i][d] == a_point[d]) count++;
      }
      if (count >= 4) retval = true;
    }
    return retval;
  }
}
#endif

void GyroaverageOperator::parseParams(const ParmParse& a_pp)
{
  a_pp.query("npts_interp", m_npts_interp);
  a_pp.query("npts_ghost", m_gpt);
}

void GyroaverageOperator::printParams() const
{
  if (!procID()) {
    std::cout << "Gyro-averaging operator inputs:\n";
    std::cout << "  species = " << m_name 
              << " (Z=" << m_charge
              << ", m=" << m_mass << ")\n";
    std::cout << "  npts_interp = " << m_npts_interp << "\n";
    std::cout << "  npts_ghost = " << m_gpt << "\n";
    std::cout << "\n";
  }
}

void GyroaverageOperator::define( const PhaseGeom&          a_phase_geom,
                                  const std::string&        a_name,
                                  Real                      a_mass,
                                  Real                      a_charge )
{
  CH_assert(!isDefined());
  m_name = a_name;
  m_mass = a_mass;
  m_charge = a_charge;
  m_phase_geom = &a_phase_geom;

  /* default values */
  m_npts_interp = 5;
  m_gpt = 5;
  m_tolerance = 1e-12;

  /* read inputs */
  ParmParse pp("gyroavg");
  parseParams(pp);

  /* sanity check */
  CH_assert(m_npts_interp > 0);

  /* angular spacing of gyroaveraging points */
  Real delta_alpha = 2*(4.0*atan(1.0)) / m_npts_interp;

  /* print useful things to screen */
  printParams();

  /* start time */
#ifdef CH_MPI
  double wc_define_start = MPI_Wtime();
#endif

  /* get dmu */
  const VEL::VelCoordSys& vel_coords = a_phase_geom.velSpaceCoordSys();
  const VEL::RealVect& vel_dx = vel_coords.dx();
  const Real dmu = vel_dx[_MU_];

  /* Larmor number */
  Real larmor = a_phase_geom.larmorNumber();

  /* get magnetic geometry */
  const CFG::MagGeom& mag_geom(a_phase_geom.magGeom());

  /* get magnetic field magnitude data */
  const CFG::LevelData<CFG::FArrayBox>& cfg_b_field_magn(mag_geom.getCCBFieldMag());
  const CFG::LevelData<CFG::FArrayBox>& cfg_b_field(mag_geom.getCCBField());

  /* create CFG leveldatas for real coordinates */
  const CFG::MultiBlockCoordSys& mag_coord_sys( *(mag_geom.coordSysPtr()) );
  const CFG::DisjointBoxLayout& cfg_grids(mag_geom.grids());
  CFG::LevelData<CFG::FArrayBox> cfg_rcoords(cfg_grids, CFG_DIM, m_gpt*CFG::IntVect::Unit);

  for (CFG::DataIterator dit(cfg_grids.dataIterator()); dit.ok(); ++dit) {

    const CFG::Box& bx(cfg_grids[dit]);
    const int block_number( mag_coord_sys.whichBlock(bx) );
    const CFG::MagBlockCoordSys& block_coord_sys = static_cast<const CFG::MagBlockCoordSys&>
                                                    (*(mag_coord_sys.getCoordSys(block_number)));
    /* get the physical coordinates */
    block_coord_sys.getCellCenteredRealCoords( cfg_rcoords[dit] );
  }

  /* Inject the configuration space variables to phase space variables */
  LevelData<FArrayBox> inj_b_field, inj_b_field_magn, inj_rcoords, inj_mcoords;
  a_phase_geom.injectConfigurationToPhase(cfg_b_field, inj_b_field);
  a_phase_geom.injectConfigurationToPhase(cfg_b_field_magn, inj_b_field_magn);
  a_phase_geom.injectConfigurationToPhaseNoExchange(cfg_rcoords, inj_rcoords);

  /* define the gyroaveraging grids - phase space grids collapsed in vpar direction */
  a_phase_geom.getConfigurationPlusMuGrids( m_grids, -1 );

  /* define and allocate the gyroaveraging operator */
  m_gyroavg_op.define(m_grids, 1);

  int count_all_points = 0, count_bad_points = 0;
  
  for (DataIterator dit(m_grids.dataIterator()); dit.ok(); ++dit) {

    /* create configuration space boxes for physical and mapped coordinates */
    CFG::Box rcoords_bx;
    a_phase_geom.projectPhaseToConfiguration(inj_rcoords[dit].box(), rcoords_bx);
    /* create configuration space physical and mapped coordinates*/
    CFG::FArrayBox  rcoords(rcoords_bx, inj_rcoords[dit].nComp());
    a_phase_geom.projectPhaseToConfigurationLocal(inj_rcoords[dit], rcoords);

    /* configuration space magnetic field */
    CFG::Box mag_bx;
    a_phase_geom.projectPhaseToConfiguration(m_grids[dit], mag_bx);
    CFG::FArrayBox  b_field_magn(mag_bx, inj_b_field_magn[dit].nComp()),
                    b_field(mag_bx, inj_b_field[dit].nComp());
    a_phase_geom.projectPhaseToConfigurationLocal(inj_b_field_magn[dit], b_field_magn);
    a_phase_geom.projectPhaseToConfigurationLocal(inj_b_field[dit], b_field);

    /* gyroavg operator */
    BaseFab<Stencil>& op_fab = m_gyroavg_op[dit];

    /* now create the gyroaverage operator at each grid point */
    BoxIterator bit(m_grids[dit]);
    for (bit.begin(); bit.ok(); ++bit) {

      /* index of this grid point */
      const IntVect iv(bit());
      /* configuration space index of this grid point */
      CFG::IntVect i_cfg = a_phase_geom.config_restrict(iv);

      /* physical coordinates of this grid point (x,y) */
      CFG::RealVect x;
      for (int d=0; d<CFG_DIM; d++) x[d] = rcoords(i_cfg, d);

      /* magnetic field and its magnitude at this grid point */
      Real b_magn = b_field_magn(i_cfg, 0);
      Real bx = b_field(i_cfg,0),
           by = b_field(i_cfg,1),
           bz = b_field(i_cfg,2);

      /* initialize the gyroavg operator at this grid point */
      Stencil& op = op_fab(iv, 0);
      op.clear();

      /* mu */
      Real mu = iv[MU_DIR] * dmu;

      /* calculate gyroradius for this mu */
      Real rho = computeGyroradius(larmor, a_mass, a_charge, b_magn, mu);

      if (rho < m_tolerance) {

        op.push_back(StencilPoint(i_cfg, 1.0));

      } else {

#if CFG_DIM==2
        /* calculate major and minor radii of the elliptical gyro-orbit */
        Real rho_a, rho_b;
        if (b_magn == 0.0) {
          rho_a = rho_b = 0.0;
        } else {
          rho_a = rho;
          rho_b = rho * (by/b_magn);
        }
#else
        /* normalized magnetic field vector */
        CFG::RealVect nb;
        nb[0] = bx / sqrt(bx*bx+by*by+bz*bz);
        nb[1] = by / sqrt(bx*bx+by*by+bz*bz);
        nb[2] = bz / sqrt(bx*bx+by*by+bz*bz);
#endif
  
        /* create interpolation points */
        for (int n = 0; n < m_npts_interp; n++) {
          
          /* compute angle */
          Real alpha = n * delta_alpha;
  
          /* compute coordinates of interpolation point */
          CFG::RealVect x_pts;
#if CFG_DIM==2
          {
            x_pts[0] = x[0] + rho_a * cos(alpha);
            x_pts[1] = x[1] + rho_b * sin(alpha);
          }
#else
          {
            /* find one normal vector to the magnetic field */
            CFG::RealVect v1;
            if (sqrt(bx*bx+bz*bz) > m_tolerance) {
              v1[0] = nb[2];
              v1[1] = 0.0;
              v1[2] = -nb[0];
            } else {
              v1[0] = 0.0;
              v1[1] = nb[2];
              v1[2] = -nb[1];
            }

            /* find the other normal vector to the magnetic field */
            CFG::RealVect v2;
            v2[0] = nb[1]*v1[2] - nb[2]*v1[1];
            v2[1] = nb[2]*v1[0] - nb[0]*v1[2];
            v2[2] = nb[0]*v1[1] - nb[1]*v1[0];

            /* compute coordinates of interpolation point */
            x_pts[0] = x[0] + rho * (v1[0]*cos(alpha) + v2[0]*sin(alpha));
            x_pts[1] = x[1] + rho * (v1[1]*cos(alpha) + v2[1]*sin(alpha));
            x_pts[2] = x[2] + rho * (v1[2]*cos(alpha) + v2[2]*sin(alpha));
          }
#endif
  
          /* weight of this point */
          Real pt_weight = 1.0 / ((Real)m_npts_interp);
  
          /* let's find the grid point closest to the interpolation point */
          CFG::IntVect i_pts(i_cfg);
          Real delta = distanceFunction(rcoords, i_pts, x_pts);
          while(1) {
            bool did_i_move = false;
            for (int d=0; d<CFG_DIM; d++) {
              for (int sign=-1; sign<=1; sign+=2) {
                CFG::IntVect i_new = i_pts + sign*CFG::BASISV(d);
                if (rcoords_bx.contains(i_new)) {
                  Real delta_new = distanceFunction(rcoords, i_new, x_pts);
                  if (delta_new < delta) {
                    delta = delta_new;
                    i_pts = i_new;
                    did_i_move = true;
                  }
                }
              }
            }
            if (!did_i_move) break;
          }
  
          if (distanceFunction(rcoords, i_pts, x_pts) > distanceFunction(rcoords, i_cfg, x_pts)) {
            std::cout << "Rank " << procID() << ": Error in GyroaverageOperator::define() - "
                      << "probe moved away from interp point on gyro-orbit!\n";
            MayDay::Error("Error encountered; see message above");
          }
  
          /* Now, let's find the other nearby points for the interpolation */
          std::vector<int> sign(CFG_DIM,0);
          for (int d=0; d<CFG_DIM; d++) {

            CFG::IntVect j1 = i_pts + CFG::BASISV(d);
            Real distance_1 = (rcoords_bx.contains(j1) ? distanceFunction(rcoords, j1, x_pts) : DBL_MAX);

            CFG::IntVect j2 = i_pts - CFG::BASISV(d);
            Real distance_2 = (rcoords_bx.contains(j2) ? distanceFunction(rcoords, j2, x_pts) : DBL_MAX);

            CH_assert(rcoords_bx.contains(j1) || rcoords_bx.contains(j2));
            sign[d] = (distance_1 < distance_2 ? 1 : -1);
          }

          std::vector<CFG::IntVect> index_set(0);
#if CFG_DIM==2
          {
            index_set.push_back(i_pts);
            index_set.push_back(i_pts + sign[0]*CFG::BASISV(0));
            index_set.push_back(i_pts + sign[1]*CFG::BASISV(1));
            index_set.push_back(i_pts + sign[0]*CFG::BASISV(0) + sign[1]*CFG::BASISV(1));
          }
#else
          {
            index_set.push_back(i_pts);
            index_set.push_back(i_pts + sign[0]*CFG::BASISV(0));
            index_set.push_back(i_pts + sign[1]*CFG::BASISV(1));
            index_set.push_back(i_pts + sign[0]*CFG::BASISV(0) + sign[1]*CFG::BASISV(1));
            index_set.push_back(i_pts + sign[2]*CFG::BASISV(2));
            index_set.push_back(i_pts + sign[0]*CFG::BASISV(0) + sign[2]*CFG::BASISV(2));
            index_set.push_back(i_pts + sign[1]*CFG::BASISV(1) + sign[2]*CFG::BASISV(2));
            index_set.push_back(i_pts + sign[0]*CFG::BASISV(0) + sign[1]*CFG::BASISV(1) + sign[2]*CFG::BASISV(2));
          }
#endif
  
          /* get the stencil to interpolate a variable at this point from grid data */
          Stencil op_pt;
          op_pt.clear();
          bool status;
          getInterpStencil(op_pt, rcoords, index_set, x_pts, status);

          if (status == false) {
            count_bad_points++;
          }
          count_all_points++;
  
          /* scale interp weights by the weight of this point 
           * and push to the overall averaging operator */
          for (int k = 0; k < op_pt.size(); k++) {
            op_pt[k].second *= pt_weight;
            op.push_back(op_pt[k]);
          }
  
        }

      }

    }

  }

  int count_all_points_global = 0, count_bad_points_global = 0;
#ifdef CH_MPI
  MPI_Allreduce(&count_all_points, &count_all_points_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_bad_points, &count_bad_points_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  count_all_points_global = count_all_points;
  count_bad_points_global = count_bad_points;
#endif

  /* stop time */
#ifdef CH_MPI
  double wc_define_stop = MPI_Wtime();

  if (!procID()) {
    std::cout << "Gyro-averaging operator for species " 
              << m_name
              << ":\n"
              << "  " << count_bad_points_global
              << " of " << count_all_points_global
              << " were outside local domains with ghost cells.\n"
              << "  setup time = "
              << wc_define_stop-wc_define_start
              << " (seconds) \n\n";
  }
#endif

  m_is_defined = true;
  return;
}

void GyroaverageOperator::define(const GyroaverageOperator& a_gyroavg)
{
  CH_assert(!isDefined());
  m_name = a_gyroavg.m_name;
  m_mass = a_gyroavg.m_mass;
  m_charge = a_gyroavg.m_charge;
  m_npts_interp = a_gyroavg.m_npts_interp;
  m_gpt = a_gyroavg.m_gpt;
  m_phase_geom = a_gyroavg.m_phase_geom;

  const LevelData< BaseFab<Stencil> >& copy_from = a_gyroavg.op();
  const DisjointBoxLayout& grids = copy_from.disjointBoxLayout();
  const IntVect& gv = copy_from.ghostVect();
  int ncomp = copy_from.nComp();

  m_grids.define(grids);
  m_gyroavg_op.define(m_grids, ncomp, gv);
  copy_from.copyTo(m_gyroavg_op);

  m_is_defined = true;
}

void GyroaverageOperator::getInterpStencil( Stencil&                          a_stencil,
                                            const CFG::FArrayBox&             a_coords,
                                            const std::vector<CFG::IntVect>&  a_idx_set,
                                            const CFG::RealVect&              a_x,
                                            bool&                             a_status)
{
  a_stencil.clear();

  Real  ub = 1.0+m_tolerance,
        lb = -m_tolerance;

#if CFG_DIM==2

  UnstructuredInterpUtils::bilinearInterpStencil( a_stencil,
                                                  a_coords,
                                                  a_idx_set,
                                                  a_x );
  /* checking the bilinear interp coefficients */
  CH_assert(a_stencil.size() == 4);
  Real  c0 = a_stencil[0].second,
        c1 = a_stencil[1].second,
        c2 = a_stencil[2].second,
        c3 = a_stencil[3].second;
  if (    (c0 > ub) || (c1 > ub) || (c2 > ub) || (c3 > ub)
      ||  (c0 < lb) || (c1 < lb) || (c2 < lb) || (c3 < lb) ) {
    a_status = false;
  } else {
    a_status = true;
  }

#else 

  UnstructuredInterpUtils::trilinearInterpStencil(  a_stencil,
                                                    a_coords,
                                                    a_idx_set,
                                                    a_x );
  /* checking the bilinear interp coefficients */
  CH_assert(a_stencil.size() == 8);
  Real  c0 = a_stencil[0].second,
        c1 = a_stencil[1].second,
        c2 = a_stencil[2].second,
        c3 = a_stencil[3].second,
        c4 = a_stencil[4].second,
        c5 = a_stencil[5].second,
        c6 = a_stencil[6].second,
        c7 = a_stencil[7].second;
  if (    (c0 > ub) || (c1 > ub) || (c2 > ub) || (c3 > ub)
      ||  (c4 > ub) || (c5 > ub) || (c6 > ub) || (c7 > ub)
      ||  (c0 < lb) || (c1 < lb) || (c2 < lb) || (c3 < lb) 
      ||  (c4 < lb) || (c5 < lb) || (c6 < lb) || (c7 < lb) ) {
    a_status = false;
  } else {
    a_status = true;
  }

#endif

  /* done */
  return;
}

void GyroaverageOperator::applyOp(FArrayBox&               a_var_bar,
                                  const FArrayBox&         a_var,
                                  const BaseFab<Stencil>&  a_op) const
{
  CH_assert(isDefined());
  CH_assert(a_var_bar.nComp() == a_var.nComp());
  
  const Box& var_bar_bx = a_var_bar.box();
  const Box& op_bx = a_op.box();
  const Box& var_bx = a_var.box();
  CH_assert(var_bar_bx.contains(op_bx));

  const int ncomp = a_var.nComp();

  BoxIterator bit(var_bar_bx);
  for (bit.begin(); bit.ok(); ++bit) {

    const IntVect iv(bit());

    IntVect stencil_idx(iv);
    stencil_idx[VPARALLEL_DIR] = op_bx.smallEnd(VPARALLEL_DIR);
    const Stencil& op = a_op(stencil_idx, 0);

    VEL::IntVect vel_idx = m_phase_geom->vel_restrict(iv);
    for (int d = VPARALLEL_DIR; d <= MU_DIR; d++) {
      int imin = var_bx.smallEnd(d),
          imax = var_bx.bigEnd(d),
          i    = vel_idx[d-VPARALLEL_DIR];

      if (imin == imax) {
        vel_idx[d-VPARALLEL_DIR] = imin;
      } else {
        if ((i-imin)*(i-imax) > 0) {
          MayDay::Abort("Error in GyroaverageOperator::applyOp(): invalid velocity space index!");
        }
      }
    }

    std::vector<Real> var_avg;
    evalStencil(var_avg, a_var, op, vel_idx);

    for (int v=0; v<ncomp; v++) {
      a_var_bar(iv,v) = var_avg[v];
    }
  }

  return;
}

void GyroaverageOperator::applyOp(LevelData<FArrayBox>&                 a_var_bar,
                                  const CFG::LevelData<CFG::FArrayBox>& a_var ) const
{
  CH_assert(isDefined());
  CH_assert(a_var.isDefined());

  /* if the gyroaveraged var is already defined, clear it */
  if (a_var_bar.isDefined()) a_var_bar.clear();
  /* and now define it */
  a_var_bar.define(m_grids, a_var.nComp());

  /* inject var into phase space */
  LevelData<FArrayBox> injected_var;
  m_phase_geom->injectConfigurationToPhase(a_var, injected_var);

  /* create a version with ghost points in the configuration dimensions */
  IntVect ghost_vec(IntVect::Zero); for (int d=0; d<CFG_DIM; d++) ghost_vec[d] = m_gpt;
  LevelData<FArrayBox> var_wghosts;
  var_wghosts.define(injected_var.disjointBoxLayout(), a_var.nComp(), ghost_vec);
  SpaceUtils::copyAndFillGhostCellsSimple(var_wghosts, injected_var);

  /* now apply the gyroaveraging */
  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    applyOp(a_var_bar[dit], var_wghosts[dit], m_gyroavg_op[dit]);
  }

  return;
}

void GyroaverageOperator::applyOp(LevelData<FArrayBox>&       a_var_bar,
                                  const LevelData<FArrayBox>& a_var ) const
{
  CH_assert(isDefined());
  CH_assert(a_var.isDefined());

  /* if the gyroaveraged var is already defined, clear it */
  if (a_var_bar.isDefined()) a_var_bar.clear();
  /* and now define it */
  a_var_bar.define(a_var.disjointBoxLayout(), a_var.nComp());

  /* create a version with ghost points in the configuration dimensions */
  IntVect ghost_vec(IntVect::Zero); for (int d=0; d<CFG_DIM; d++) ghost_vec[d] = m_gpt;
  LevelData<FArrayBox> var_wghosts;
  var_wghosts.define(a_var.disjointBoxLayout(), a_var.nComp(), ghost_vec);
  SpaceUtils::copyAndFillGhostCellsSimple(var_wghosts, a_var);

  /* now apply the gyroaveraging */
  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    applyOp(a_var_bar[dit], var_wghosts[dit], m_gyroavg_op[dit]);
  }

  return;
}

void GyroaverageOperator::gyroaveragedEField( LevelData<FArrayBox>&                 a_E_field,
                                              const CFG::LevelData<CFG::FArrayBox>& a_phi,
                                              const int                             a_order ) const
{
  CH_assert(isDefined());
  CH_assert(a_phi.nComp() == 1);

  /* even for CFG_DIM=2, efield has 3 components.
   * The second component (toroidal) is 0 */
  int n_efield_comp = 3;

  /* compute gyroaveraged potential */
  LevelData<FArrayBox> phi_bar;
  applyOp(phi_bar, a_phi);
  CH_assert(phi_bar.nComp() == 1);

  /* define the Efield */
  if (a_E_field.isDefined()) a_E_field.clear();
  IntVect gvec = IntVect::Unit;
  a_E_field.define(m_grids, n_efield_comp, gvec);

  /* compute the gradient in configuration space */
  m_phase_geom->computeGradientCfgSpace(phi_bar, a_E_field, a_order);

  for (DataIterator dit(a_E_field.dataIterator()); dit.ok(); ++dit) {
    /* negate it since E-field is -grad */
    a_E_field[dit].negate();

#if CFG_DIM==2
    /* if CFG_DIM = 2, then PhaseGeom::computeGradientCfgSpace() would have
     * populated the first two components with the radial and poloidal gradients.
     * So, push the 2nd component to 3rd (since the 3rd component is the 
     * poloidal one, and set 2nd component to zero */
    a_E_field[dit].copy(a_E_field[dit], CFG_DIM-1, CFG_DIM);
    a_E_field[dit].setVal(0.0, CFG_DIM-1);
#endif
  }

  return;
}

void GyroaverageOperator::gyroaveragedEField( LevelData<FluxBox>&                   a_E_field,
                                              const CFG::LevelData<CFG::FArrayBox>& a_phi,
                                              const int                             a_order ) const
{
  CH_assert(isDefined());
  CH_assert(a_phi.nComp() == 1);

  /* even for CFG_DIM=2, efield has 3 components.
   * The second component (toroidal) is 0 */
  int n_efield_comp = 3;

  /* compute gyroaveraged potential */
  LevelData<FArrayBox> phi_bar;
  applyOp(phi_bar, a_phi);
  CH_assert(phi_bar.nComp() == 1);

  /* define the Efield */
  if (a_E_field.isDefined()) a_E_field.clear();
  IntVect gvec = IntVect::Unit;
  a_E_field.define(m_grids, n_efield_comp, gvec);

  /* compute the gradient in configuration space */
  m_phase_geom->computeGradientCfgSpace(phi_bar, a_E_field, a_order);

  for (DataIterator dit(a_E_field.dataIterator()); dit.ok(); ++dit) {
    for (int dir=0; dir<SpaceDim; dir++) {

      /* negate it since E-field is -grad */
      a_E_field[dit][dir].negate();

#if CFG_DIM==2
      /* if CFG_DIM = 2, then PhaseGeom::computeGradientCfgSpace() would have
       * populated the first two components with the radial and poloidal gradients.
       * So, push the 2nd component to 3rd (since the 3rd component is the 
       * poloidal one, and set 2nd component to zero */
      a_E_field[dit][dir].copy(a_E_field[dit][dir], CFG_DIM-1, CFG_DIM);
      a_E_field[dit][dir].setVal(0.0, CFG_DIM-1);
#endif
    }
  }


  return;
}

void GyroaverageOperator::evalStencil(  std::vector<Real>&    a_var_avg,
                                        const FArrayBox&      a_var,
                                        const Stencil&        a_op,
                                        const VEL::IntVect&   a_vel_idx ) const
{
  a_var_avg = std::vector<Real>(a_var.nComp(),0.0);
  for (int i=0; i<a_op.size(); i++) {
    IntVect idx = m_phase_geom->tensorProduct(a_op[i]._point_, a_vel_idx);
    for (int v=0; v<a_var.nComp(); v++) {
      a_var_avg[v] += a_op[i]._weight_ * a_var(idx, v);
    }
  }
  return;
}

#include "NamespaceFooter.H"
