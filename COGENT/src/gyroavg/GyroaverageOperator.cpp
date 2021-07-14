#include "GyroaverageOperator.H"
#include "Directions.H"
#include "PhaseGeom.H"
#include "SpaceUtils.H"

#include "NamespaceHeader.H"

#if 0
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
#endif

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
  a_pp.query("interp_order", m_interp_order);
  a_pp.query("npts_interp", m_npts_interp);
  a_pp.query("npts_ghost", m_gpt);
  a_pp.query("reduce_interp_order", m_reduce_interp_order);
  a_pp.query("interp_mode_mapped", m_interp_mode_mapped);

  a_pp.query("debug_mode", m_debug_mode);
  if (m_debug_mode) {
    a_pp.query("debug_mode_fake_gyroradius", m_debug_fakegyroradius);
    a_pp.query("debug_mode_fake_gyroradius_dr", m_debug_fakegyroradius_dr);
  }
}

void GyroaverageOperator::printParams() const
{
  if (!procID()) {
    std::cout << "Gyro-averaging operator inputs:\n";
    std::cout << "  species = " << m_name 
              << " (Z=" << m_charge
              << ", m=" << m_mass << ")\n";
    std::cout << "  npts_interp = " << m_npts_interp << "\n";
    std::cout << "  interp_order = " << m_interp_order << "\n";
    std::cout << "  npts_ghost = " << m_gpt << "\n";
    std::cout << "  interp mode = " 
              << (m_interp_mode_mapped ? "mapped space" : "physical space")
              << "\n";
    std::cout << "  reduce order at degenerate cells = " 
              << (m_reduce_interp_order ? "yes" : "no")
              << "\n";
    if (m_debug_mode) {
      std::cout << "  **debug mode: " << "true" << "\n";
      if (m_debug_fakegyroradius) {
        std::cout << "  **debug mode fake gyroradius: " << "true" << "\n";
        std::cout << "  **debug mode fake gyroradius_dr: " 
                  << m_debug_fakegyroradius_dr << "\n";
      }
    }
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

  m_max_muidx = -1;
  int max_muidx_lcl = -1;

  m_max_gyroradius = -DBL_MAX;
  Real max_gyroradius_lcl = -DBL_MAX;

  /* default values */
  m_npts_interp = 5;
  m_reduce_interp_order = false;
  m_interp_mode_mapped = true;
  m_interp_order = 3;
  m_gpt = 5;
  m_tolerance = 1e-10;

  m_debug_mode = false;
  m_debug_fakegyroradius = false;
  m_debug_fakegyroradius_dr = 0.0;

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

  /* get magnetic field and its magnitude */
  const LevelData<FArrayBox>& b_field_magn(a_phase_geom.getBFieldMagnitude());
  const LevelData<FArrayBox>& b_field(a_phase_geom.getBField());

  /* define the gyroaveraging grids - phase space grids collapsed in vpar direction */
  a_phase_geom.getConfigurationPlusMuGrids( m_grids, -1 );

  if (m_debug_mode) {
    /* define and allocate the gyroradius array */
    m_gyroradius.define(m_grids, 1);
  }

  /* define and allocate the gyroaveraging operator */
  m_gyroavg_op.define(m_grids, 1);

  int count_all_points = 0, count_bad_points = 0;
  Real max_badness = 0.0;
  int stencil_size_naive = 0, stencil_size_compressed = 0;
  
  for (DataIterator dit(m_grids.dataIterator()); dit.ok(); ++dit) {

    const PhaseBlockCoordSys& phase_coord_sys = a_phase_geom.getBlockCoordSys(m_grids[dit]);

    CFG::Box coords_bx;
    CFG::FArrayBox rcoords, mcoords, icoords;

    {
      Box bx_phase = m_grids[dit];
      bx_phase.grow(m_gpt);
      a_phase_geom.projectPhaseToConfiguration(bx_phase, coords_bx);

      rcoords.define(coords_bx, PDIM);
      mcoords.define(coords_bx, PDIM);
      icoords.define(coords_bx, PDIM);

      FArrayBox rcoords_phase(bx_phase, PDIM);
      FArrayBox mcoords_phase(bx_phase, PDIM);
      phase_coord_sys.getCellCenteredRealCoords(rcoords_phase);
      phase_coord_sys.getCellCenteredMappedCoords(mcoords_phase);

      a_phase_geom.projectPhaseToConfigurationLocal(rcoords_phase, rcoords);
      a_phase_geom.projectPhaseToConfigurationLocal(mcoords_phase, mcoords);

      if (m_interp_mode_mapped) {
        icoords.copy(mcoords);
      } else {
        icoords.copy(rcoords);
      }
    }

    /* gyroavg operator */
    BaseFab<Stencil>& op_fab = m_gyroavg_op[dit];

    IntVect idx_min;
    {
      Box bbox = b_field_magn[dit].box();
      idx_min = bbox.smallEnd();
    }

    /* now create the gyroaverage operator at each grid point */
    BoxIterator bit(m_grids[dit]);
    for (bit.begin(); bit.ok(); ++bit) {

      /* index of this grid point */
      const IntVect iv(bit());

      IntVect iv_b(bit()); 
      iv_b[VPARALLEL_DIR] = idx_min[VPARALLEL_DIR]; 
      iv_b[MU_DIR] = idx_min[MU_DIR];

      /* configuration space index of this grid point */
      CFG::IntVect i_cfg = a_phase_geom.config_restrict(iv);

      /* physical coordinates of this grid point (x,y) */
      CFG::RealVect x;
      for (int d=0; d<CFG_DIM; d++) x[d] = rcoords(i_cfg, d);

      /* magnetic field magnitude at this grid point */
      Real b_magn = b_field_magn[dit](iv_b, 0);

      /* initialize the gyroavg operator at this grid point */
      Stencil& op = op_fab(iv, 0);
      op.clear();

      /* mu */
      if (iv[MU_DIR] > max_muidx_lcl) max_muidx_lcl = iv[MU_DIR];
      Real mu = (iv[MU_DIR]+0.5) * dmu;

      /* calculate gyroradius for this mu */
      Real rho = 0;
      if (m_debug_mode && m_debug_fakegyroradius) {
        rho = m_debug_fakegyroradius_dr * iv[MU_DIR];
      } else {
        rho = computeGyroradius(larmor, a_mass, a_charge, b_magn, mu);
      }
      if (rho > max_gyroradius_lcl) max_gyroradius_lcl = rho;
      if (m_debug_mode) m_gyroradius[dit](iv,0) = rho;

      if (rho < m_tolerance) {

        op.push_back(StencilPoint(i_cfg, 1.0));

      } else {

#if CFG_DIM==2
        /* calculate major and minor radii of the elliptical gyro-orbit */
        Real by = b_field[dit](iv_b,POLOIDAL_DIR);
        Real rho_a, rho_b;
        if (b_magn == 0.0) {
          rho_a = rho_b = 0.0;
        } else {
          rho_a = rho;
          rho_b = rho * (by/b_magn);
        }
#else
        /* find normal vector to magnetic field */
        CFG::RealVect nb, nb_abs;
        for (int dir=0; dir < CFG_DIM; dir++) {
          nb[dir] = b_field[dit](iv_b,dir)/b_magn;
          nb_abs[dir] = (nb[dir] < 0 ? -nb[dir] : nb[dir]);
        }

        /* find one normal vector to the magnetic field */
        CFG::RealVect v1(CFG::RealVect::Zero);
        if (      (nb_abs[POLOIDAL_DIR] < m_tolerance)
              &&  (nb_abs[TOROIDAL_DIR] < m_tolerance) ) {
          v1[RADIAL_DIR] = nb[POLOIDAL_DIR];
          v1[POLOIDAL_DIR] = -nb[RADIAL_DIR];
        } else {
          v1[TOROIDAL_DIR] = nb[POLOIDAL_DIR];
          v1[POLOIDAL_DIR] = -nb[TOROIDAL_DIR];
        }

        /* normalize v1 */
        Real v1_magn = sqrt(v1.dotProduct(v1));
        v1 /= v1_magn;

        /* find the other normal vector to the magnetic field */
        CFG::RealVect v2(CFG::RealVect::Zero);
        v2[0] = nb[1]*v1[2] - nb[2]*v1[1];
        v2[1] = nb[2]*v1[0] - nb[0]*v1[2];
        v2[2] = nb[0]*v1[1] - nb[1]*v1[0];

        /* normalize v2 */
        Real v2_magn = sqrt(v2.dotProduct(v2));
        v2 /= v2_magn;
#endif
  
        /* create interpolation points */
        for (int n = 0; n < m_npts_interp; n++) {
          
          /* compute angle */
          Real alpha = n * delta_alpha;
  
          /* compute coordinates of interpolation point */
          CFG::RealVect x_pts;
#if CFG_DIM==2
          x_pts[0] = x[0] + rho_a * cos(alpha);
          x_pts[1] = x[1] + rho_b * sin(alpha);
#else
          for (int d=0; d<CFG_DIM; d++) {
            x_pts[d] = x[d] + rho * (v1[d]*cos(alpha) + v2[d]*sin(alpha));
          }
#endif

          CFG::RealVect xi_pts;
          if (m_interp_mode_mapped) {
            RealVect x_pts_phase = a_phase_geom.config_inject(x_pts);
            RealVect xi_pts_phase = phase_coord_sys.mappedCoord(x_pts_phase);
            xi_pts = a_phase_geom.config_restrict(xi_pts_phase);
          } else {
            xi_pts = x_pts;
          }
  
          /* weight of this point */
          Real pt_weight = 1.0 / ((Real)m_npts_interp);
  
          /* let's find the grid point closest to the interpolation point */
          CFG::IntVect i_pts(i_cfg);
          Real delta = distanceFunction(icoords, i_pts, xi_pts);
          while(1) {
            bool did_i_move = false;
            for (int d=0; d<CFG_DIM; d++) {
              for (int sign=-1; sign<=1; sign+=2) {
                CFG::IntVect i_new = i_pts + sign*CFG::BASISV(d);
                if (coords_bx.contains(i_new)) {
                  Real delta_new = distanceFunction(icoords, i_new, xi_pts);
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
  
          if (  distanceFunction(icoords, i_pts, xi_pts) 
                > distanceFunction(icoords, i_cfg, xi_pts)) {
            std::cout << "Rank " << procID() << ": Error in GyroaverageOperator::define() - "
                      << "probe moved away from interp point on gyro-orbit!\n";
            MayDay::Error("Error encountered; see message above");
          }
  
          /* Now, let's find the other nearby points for the interpolation */
          std::vector<int> sign(CFG_DIM,0);
          for (int d=0; d<CFG_DIM; d++) {

            CFG::IntVect j1 = i_pts + CFG::BASISV(d);
            Real distance_1 = (coords_bx.contains(j1) ? distanceFunction(icoords, j1, xi_pts) : DBL_MAX);

            CFG::IntVect j2 = i_pts - CFG::BASISV(d);
            Real distance_2 = (coords_bx.contains(j2) ? distanceFunction(icoords, j2, xi_pts) : DBL_MAX);

            CH_assert(coords_bx.contains(j1) || coords_bx.contains(j2));
            sign[d] = (distance_1 < distance_2 ? 1 : -1);
          }

          Stencil op_pt;
          std::vector<CFG::IntVect> index_set(0);

          int min_interp_order;
          if (m_reduce_interp_order) {
            min_interp_order = 2;
          } else {
            min_interp_order = m_interp_order;
          }
          bool status;
          Real badness;

          for (int try_interp_order = m_interp_order; try_interp_order >= min_interp_order; try_interp_order--) {

            CFG::IntVect i0 = i_pts;
            
            index_set.clear();
            if ((try_interp_order == 3) || (try_interp_order == 4)) {
              for (int d=0; d<CFG_DIM; d++) {
                if (coords_bx.contains(i0 - sign[d]*CFG::BASISV(d))) {
                  i0 = i0 - sign[d]*CFG::BASISV(d);
                }
              }
              for (int d=0; d<CFG_DIM; d++) {
                while (!coords_bx.contains
                        (i0 + (try_interp_order-1)*sign[d]*CFG::BASISV(d)) 
                      ) {
                  if (coords_bx.contains(i0 - sign[d]*CFG::BASISV(d))) {
                    i0 = i0 - sign[d]*CFG::BASISV(d);
                  } else {
                    break;
                  }
                }
              }
            }
#if CFG_DIM==2
            for (int d1 = 0; d1 < try_interp_order; d1++) {
              for (int d2 = 0; d2 < try_interp_order; d2++) {
                index_set.push_back(  i0 
                                      + d1*sign[1]*CFG::BASISV(1) 
                                      + d2*sign[0]*CFG::BASISV(0));
              }
            }
#else
            for (int d1 = 0; d1 < try_interp_order; d1++) {
              for (int d2 = 0; d2 < try_interp_order; d2++) {
                for (int d3 = 0; d3 < try_interp_order; d3++) {
                  index_set.push_back(  i0 
                                        + d1*sign[2]*CFG::BASISV(2) 
                                        + d2*sign[1]*CFG::BASISV(1) 
                                        + d3*sign[0]*CFG::BASISV(0) );
                }
              }
            }
#endif
            for (int n = 0; n < index_set.size(); n++) {
              CH_assert(coords_bx.contains(index_set[n]));
            }
    
            /* get the stencil to interpolate a variable at this point 
             * from grid data */
            op_pt.clear();
            getInterpStencil( op_pt, 
                              icoords, 
                              index_set, 
                              xi_pts, 
                              try_interp_order, 
                              status, 
                              badness );
            if ((badness < 1e-14) && (status == true)) break;
          }

          if (status == false) {
            count_bad_points++;
          }
          count_all_points++;
          if (badness >= max_badness) {
            max_badness = badness;
          }
  
          /* scale interp weights by the weight of this point 
           * and push to the overall averaging operator */
          for (int k = 0; k < op_pt.size(); k++) {
            op_pt[k].second *= pt_weight;
            op.push_back(op_pt[k]);
          }
  
        }

      }

      stencil_size_naive += op.size();
      StencilUtils::compress(op);
      stencil_size_compressed += op.size();

    }

  }

  int count_all_points_global = 0, count_bad_points_global = 0;
  int stencil_size_naive_global = 0, stencil_size_compressed_global = 0;
  Real max_badness_global = 0;
#ifdef CH_MPI
  MPI_Allreduce(&count_all_points, &count_all_points_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_bad_points, &count_bad_points_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&stencil_size_naive, &stencil_size_naive_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&stencil_size_compressed, &stencil_size_compressed_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&max_badness, &max_badness_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  count_all_points_global = count_all_points;
  count_bad_points_global = count_bad_points;
  stencil_size_naive_global = stencil_size_naive;
  stencil_size_compressed_global = stencil_size_compressed;
  max_badness_global = max_badness;
#endif

  /* stop time */
#ifdef CH_MPI
  double wc_define_stop = MPI_Wtime();
#endif

#ifdef CH_MPI
  MPI_Allreduce(&max_muidx_lcl, &m_max_muidx, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&max_gyroradius_lcl, &m_max_gyroradius, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  m_max_muidx = max_muidx_lcl;
  m_max_gyroradius = max_gyroradius_lcl;
#endif

#ifdef CH_MPI
  if (!procID()) {
#endif
    printf("Gyro-averaging operator for species %s:\n", m_name.c_str());
    printf("  maximum gyroradius: %1.6e\n", m_max_gyroradius);
    printf("  %d of %d (%2.0f%%) were too far away from local domains with ghost cells.\n",
            count_bad_points_global, count_all_points_global,
            ((double)count_bad_points_global) / ((double)count_all_points_global) * 100);
    printf("  stencil size compression: %2.0f%%\n",
           (((double) stencil_size_compressed_global)/((double) stencil_size_naive_global) * 100) );
    if (max_badness_global > m_tolerance) {
      printf("  maximum badness: %1.4e\n", max_badness_global );
    }
#ifdef CH_MPI
    printf("  setup time: %lf (seconds).\n", wc_define_stop-wc_define_start );
#endif
    printf("\n");
#ifdef CH_MPI
  }
#endif

  CH_assert(m_max_muidx > 0);
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
                                            const int                         a_interp_order,
                                            bool&                             a_status,
                                            Real&                             a_badness )
{
  a_stencil.clear();

  Real  ub = 2.0,
        lb = -1.0;

#if CFG_DIM==2
  if (a_interp_order == 2) {
    UnstructuredInterpUtils::bilinearInterpStencil( a_stencil,
                                                    a_coords,
                                                    a_idx_set,
                                                    a_x,
                                                    m_tolerance );
    CH_assert(a_stencil.size() == 4);
  } else if (a_interp_order == 3) {
    UnstructuredInterpUtils::bicubicInterpStencil(  a_stencil,
                                                    a_coords,
                                                    a_idx_set,
                                                    a_x,
                                                    m_tolerance );
    CH_assert(a_stencil.size() == 9);
  } else if (a_interp_order == 4) {
    UnstructuredInterpUtils::biquarticInterpStencil(a_stencil,
                                                    a_coords,
                                                    a_idx_set,
                                                    a_x,
                                                    m_tolerance );
    CH_assert(a_stencil.size() == 16);
  } else {
    std::cout << "Invalid interpolation order for gyroaveraging!"
              << a_interp_order << "\n";
    MayDay::Error("See error message");
  }
#else 
  if (a_interp_order == 2) {
    UnstructuredInterpUtils::trilinearInterpStencil(  a_stencil,
                                                      a_coords,
                                                      a_idx_set,
                                                      a_x,
                                                      m_tolerance );
    CH_assert(a_stencil.size() == 8);
  } else if (a_interp_order == 3) {
    UnstructuredInterpUtils::tricubicInterpStencil( a_stencil,
                                                    a_coords,
                                                    a_idx_set,
                                                    a_x,
                                                    m_tolerance );
    CH_assert(a_stencil.size() == 27);
  } else {
    std::cout << "Invalid interpolation order for gyroaveraging!"
              << a_interp_order << "\n";
    MayDay::Error("See error message");
  }
#endif

  Real sum = 0;
  for (int v = 0; v < a_stencil.size(); v++) sum += a_stencil[v].second;
  a_badness = (1.0-sum);
  a_badness = (a_badness < 0 ? -a_badness : a_badness);

  a_status = true;
  for (int v = 0; v < a_stencil.size(); v++) {
    if ((a_stencil[v].second > ub) || (a_stencil[v].second < lb)) {
      a_status = false;
    }
  }

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

/* Here, if a_var_bar is not already defined, it will be defined as a
 * LevelData<FArrayBox> in phase space with the v_parallel dimension flat */
void GyroaverageOperator::applyOp(LevelData<FArrayBox>&                 a_var_bar,
                                  const CFG::LevelData<CFG::FArrayBox>& a_var ) const
{
  CH_assert(isDefined());
  CH_assert(a_var.isDefined());

  /* if a_var_bar is not defined, define it */
  if (!a_var_bar.isDefined()) {
    a_var_bar.define(m_grids, a_var.nComp());
  }

  /* inject var into phase space */
  LevelData<FArrayBox> injected_var;
  m_phase_geom->injectConfigurationToPhase(a_var, injected_var);

  /* create a version with ghost points in the configuration dimensions */
  IntVect ghost_vec(IntVect::Zero); 
  for (int d=0; d<CFG_DIM; d++) {
    ghost_vec[d] = m_gpt;
  }
  LevelData<FArrayBox> var_wghosts;
  var_wghosts.define( injected_var.disjointBoxLayout(), 
                      a_var.nComp(), 
                      ghost_vec );
  SpaceUtils::copyAndFillGhostCellsSimple(var_wghosts, injected_var);

  /* now apply the gyroaveraging */
  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    applyOp(a_var_bar[dit], var_wghosts[dit], m_gyroavg_op[dit]);
  }

  if (m_debug_mode) {
    /* useful to debug: plot the input var and the gyroaveraged
     * var at various mu to see if the operator is working as 
     * expected. We'll remove these lines once we think the 
     * gyroaverage operator is perfect. */
    m_phase_geom->plotConfigurationData("var", a_var, 0.0);
    for (int mu_idx = 0; mu_idx <= m_max_muidx; mu_idx++) {
      int vpar_idx = 0;
      VEL::IntVect velidx(vpar_idx, mu_idx);
  
      CFG::LevelData<CFG::FArrayBox> var_bar_cfg;
      m_phase_geom->getConfigurationData(velidx, a_var_bar, var_bar_cfg);
  
      char idx_str[50];
      sprintf(idx_str, "var_bar_%02d", mu_idx);
      m_phase_geom->plotConfigurationData(  std::string(idx_str), 
                                            var_bar_cfg, 
                                            0.0 );
    }
  }
  if (m_debug_mode) {
    /* useful to debug: plot the gyroradius information.
     * We'll remove these lines once we think the 
     * gyroaverage operator is perfect. */
    for (int mu_idx = 0; mu_idx <= m_max_muidx; mu_idx++) {
      int vpar_idx = 0;
      VEL::IntVect velidx(vpar_idx, mu_idx);
  
      CFG::LevelData<CFG::FArrayBox> rho_cfg;
      m_phase_geom->getConfigurationData(velidx, m_gyroradius, rho_cfg);
  
      char idx_str[50];
      sprintf(idx_str, "rho_%02d", mu_idx);
      m_phase_geom->plotConfigurationData(  std::string(idx_str), 
                                            rho_cfg, 
                                            0.0 );
    }
  }

  return;
}

/* Here, if a_var_bar is not already defined, it will be defined as a
 * LevelData<FArrayBox> in phase space with same layout as a_var */
void GyroaverageOperator::applyOp(LevelData<FArrayBox>&       a_var_bar,
                                  const LevelData<FArrayBox>& a_var ) const
{
  CH_assert(isDefined());
  CH_assert(a_var.isDefined());

  /* if a_var_bar is not defined, define it */
  if (!a_var_bar.isDefined()) {
    a_var_bar.define(a_var.disjointBoxLayout(), a_var.nComp());
  }

  /* create a version with ghost points in the configuration dimensions */
  IntVect ghost_vec(IntVect::Zero); 
  for (int d=0; d<CFG_DIM; d++) {
    ghost_vec[d] = m_gpt;
  }
  LevelData<FArrayBox> var_wghosts;
  var_wghosts.define(a_var.disjointBoxLayout(), a_var.nComp(), ghost_vec);
  SpaceUtils::copyAndFillGhostCellsSimple(var_wghosts, a_var);

  /* now apply the gyroaveraging */
  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    applyOp(a_var_bar[dit], var_wghosts[dit], m_gyroavg_op[dit]);
  }

  if (m_debug_mode) {
    /* useful to debug: plot the input var and the gyroaveraged
     * var at various mu to see if the operator is working as 
     * expected. We'll remove these lines once we think the 
     * gyroaverage operator is perfect. */
    for (int mu_idx = 0; mu_idx <= m_max_muidx; mu_idx++) {
      int vpar_idx = 0;
      VEL::IntVect velidx(vpar_idx, mu_idx);
  
      CFG::LevelData<CFG::FArrayBox> var_cfg;
      m_phase_geom->getConfigurationData(velidx, a_var, var_cfg);
  
      char idx_str[50];
      sprintf(idx_str, "var_%02d", mu_idx);
      m_phase_geom->plotConfigurationData(  std::string(idx_str), 
                                            var_cfg, 
                                            0.0 );
    }
    for (int mu_idx = 0; mu_idx <= m_max_muidx; mu_idx++) {
      int vpar_idx = 0;
      VEL::IntVect velidx(vpar_idx, mu_idx);
  
      CFG::LevelData<CFG::FArrayBox> var_bar_cfg;
      m_phase_geom->getConfigurationData(velidx, a_var_bar, var_bar_cfg);
  
      char idx_str[50];
      sprintf(idx_str, "var_bar_%02d", mu_idx);
      m_phase_geom->plotConfigurationData(  std::string(idx_str), 
                                            var_bar_cfg, 
                                            0.0 );
    }
  }
  if (m_debug_mode) {
    /* useful to debug: plot the gyroradius information.
     * We'll remove these lines once we think the 
     * gyroaverage operator is perfect. */
    for (int mu_idx = 0; mu_idx <= m_max_muidx; mu_idx++) {
      int vpar_idx = 0;
      VEL::IntVect velidx(vpar_idx, mu_idx);
  
      CFG::LevelData<CFG::FArrayBox> rho_cfg;
      m_phase_geom->getConfigurationData(velidx, m_gyroradius, rho_cfg);
  
      char idx_str[50];
      sprintf(idx_str, "rho_%02d", mu_idx);
      m_phase_geom->plotConfigurationData(  std::string(idx_str), 
                                            rho_cfg, 
                                            0.0 );
    }
  }

  return;
}

/* Here, if a_var_tilde is not already defined, it will be defined as a
 * LevelData<FArrayBox> in phase space with the v_parallel dimension flat */
void GyroaverageOperator::applyOpTwice( LevelData<FArrayBox>&       a_var_tilde,
                                        const LevelData<FArrayBox>& a_var ) const
{
  CH_assert(isDefined());
  CH_assert(a_var.isDefined());

  LevelData<FArrayBox> var_bar(m_grids, a_var.nComp());
  applyOp(var_bar, a_var);

  if (!a_var_tilde.isDefined()) {
    a_var_tilde.define(m_grids, a_var.nComp());
  }
  applyOp(a_var_tilde, var_bar);

  if (m_debug_mode) {
    /* useful to debug: plot the input var and the gyroaveraged
     * var at various mu to see if the operator is working as 
     * expected. We'll remove these lines once we think the 
     * gyroaverage operator is perfect. */
    for (int mu_idx = 0; mu_idx <= m_max_muidx; mu_idx++) {
      int vpar_idx = 0;
      VEL::IntVect velidx(vpar_idx, mu_idx);
  
      CFG::LevelData<CFG::FArrayBox> var_cfg;
      m_phase_geom->getConfigurationData(velidx, a_var, var_cfg);
  
      char idx_str[50];
      sprintf(idx_str, "var_%02d", mu_idx);
      m_phase_geom->plotConfigurationData(  std::string(idx_str), 
                                            var_cfg, 
                                            0.0 );
    }
    for (int mu_idx = 0; mu_idx <= m_max_muidx; mu_idx++) {
      int vpar_idx = 0;
      VEL::IntVect velidx(vpar_idx, mu_idx);
  
      CFG::LevelData<CFG::FArrayBox> var_tilde_cfg;
      m_phase_geom->getConfigurationData(velidx, a_var_tilde, var_tilde_cfg);
  
      char idx_str[50];
      sprintf(idx_str, "var_tilde_%02d", mu_idx);
      m_phase_geom->plotConfigurationData(  std::string(idx_str), 
                                            var_tilde_cfg, 
                                            0.0 );
    }
  }
  if (m_debug_mode) {
    /* useful to debug: plot the gyroradius information.
     * We'll remove these lines once we think the 
     * gyroaverage operator is perfect. */
    for (int mu_idx = 0; mu_idx <= m_max_muidx; mu_idx++) {
      int vpar_idx = 0;
      VEL::IntVect velidx(vpar_idx, mu_idx);
  
      CFG::LevelData<CFG::FArrayBox> rho_cfg;
      m_phase_geom->getConfigurationData(velidx, m_gyroradius, rho_cfg);
  
      char idx_str[50];
      sprintf(idx_str, "rho_%02d", mu_idx);
      m_phase_geom->plotConfigurationData(  std::string(idx_str), 
                                            rho_cfg, 
                                            0.0 );
    }
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
 
  LevelData<FArrayBox> phi_phase;
  m_phase_geom->injectAndExpandConfigurationToPhase(a_phi, phi_phase);

  /* compute gyroaveraged potential */
  LevelData<FArrayBox> phi_bar;
  applyOp(phi_bar, phi_phase);
  CH_assert(phi_bar.nComp() == 1);

  /* define a temporary Efield */
  IntVect gvec = IntVect::Unit;
  LevelData<FluxBox> E_field_wghosts;
  E_field_wghosts.define(phi_bar.disjointBoxLayout(), SpaceDim, gvec);

  /* compute the gradient in configuration space */
  m_phase_geom->computeGradient(phi_bar, E_field_wghosts, a_order);

  /* define the Efield */
  if (a_E_field.isDefined()) a_E_field.clear();
  a_E_field.define(phi_bar.disjointBoxLayout(), n_efield_comp, IntVect::Zero);
  for (DataIterator dit(a_E_field.dataIterator()); dit.ok(); ++dit) {
    for (int dir=0; dir<SpaceDim; dir++) {

      a_E_field[dit][dir].setVal(0.0);
#if CFG_DIM==2
      a_E_field[dit][dir].copy(E_field_wghosts[dit][dir], RADIAL_DIR, 0);
      a_E_field[dit][dir].copy(E_field_wghosts[dit][dir], POLOIDAL_DIR, 2);
#else
      a_E_field[dit][dir].copy( E_field_wghosts[dit][dir], 
                                RADIAL_DIR, 
                                RADIAL_DIR, 
                                CFG_DIM );
#endif

      /* negate it since E-field is -grad */
      a_E_field[dit][dir].negate();
    }
  }

  if (m_debug_mode) {
    /* useful to debug: plot the input phi and the gyroaveraged
     * E field at various mu to see if the operator is working as 
     * expected. We'll remove these lines once we think the 
     * gyroaverage operator is perfect. */
    m_phase_geom->plotConfigurationData("phi", a_phi, 0.0);
    for (int mu_idx = 0; mu_idx <= m_max_muidx; mu_idx++) {
      int vpar_idx = 0;
      VEL::IntVect velidx(vpar_idx, mu_idx);
  
      CFG::LevelData<CFG::FluxBox> Efield_cfg;
      m_phase_geom->getConfigurationData(velidx, a_E_field, Efield_cfg);
  
      char idx_str[50];
      sprintf(idx_str, "Efield_%02d_", mu_idx);
      m_phase_geom->plotConfigurationData(  std::string(idx_str), 
                                            Efield_cfg, 
                                            0.0 );
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
