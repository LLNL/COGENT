#include "GyroaverageOperator.H"
#include "Directions.H"
#include "PhaseGeom.H"

#include "NamespaceHeader.H"

static void fillGhostCells( FArrayBox& a_phi,
                            const Box& a_bx_int,
                            const int  a_dir )
{
  const Box& phi_bx = a_phi.box();

  for (int side=-1; side<2; side+=2) {

    int ii;
    Box bdrybox;

    if (side == -1) {
      int n_gpt = a_bx_int.smallEnd(a_dir) - phi_bx.smallEnd(a_dir);
      ii = a_bx_int.smallEnd(a_dir);
      bdrybox = adjCellLo(a_bx_int, a_dir, n_gpt);
    } else {
      int n_gpt = phi_bx.bigEnd(a_dir) - a_bx_int.bigEnd(a_dir);
      ii = a_bx_int.bigEnd(a_dir);
      bdrybox = adjCellHi(a_bx_int, a_dir, n_gpt);
    }

    BoxIterator bit(bdrybox);
    for (bit.begin(); bit.ok(); ++bit) {

      IntVect iv(bit());
      int j = (Real) (side < 0 ? ii-iv[a_dir] : iv[a_dir]-ii );

      IntVect i_int_0(iv); i_int_0[a_dir] = ii;
      IntVect i_int_1(iv); i_int_1[a_dir] = ii - side;
      IntVect i_int_2(iv); i_int_2[a_dir] = ii - 2*side;
      IntVect i_int_3(iv); i_int_3[a_dir] = ii - 3*side;

      Real c0 = ((1.0+j)*(2.0+j)*(3.0+j))/6.0;
      Real c1 = -(j*(2.0+j)*(3.0+j))/2.0;
      Real c2 = (j*(1.0+j)*(3.0+j))/2.0;
      Real c3 = -(j*(1.0+j)*(2.0+j))/6.0;

      for (int n=0; n<a_phi.nComp(); n++) {
//        a_phi(iv,n) =   a_phi(i_int_0,n);
//        a_phi(iv,n) =   a_phi(i_int_0,n) + j*(a_phi(i_int_0,n)-a_phi(i_int_1,n));
        a_phi(iv,n) =   c0 * a_phi(i_int_0,n)
                      + c1 * a_phi(i_int_1,n)
                      + c2 * a_phi(i_int_2,n)
                      + c3 * a_phi(i_int_3,n);
      }

    }
  }

  return;
}

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

  /* read inputs */
  ParmParse pp("gyroavg");
  parseParams(pp);

  /* sanity check */
  CH_assert(m_npts_interp > 0);

  /* angular spacing of gyroaveraging points */
  Real delta_alpha = 2*(4.0*atan(1.0)) / m_npts_interp;

  /* print useful things to screen */
  printParams();

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

  /* create CFG leveldatas for mapped and real coordinates */
  const CFG::MultiBlockCoordSys& mag_coord_sys( *(mag_geom.coordSysPtr()) );
  const CFG::DisjointBoxLayout& cfg_grids(mag_geom.grids());
  CFG::LevelData<CFG::FArrayBox> cfg_rcoords(cfg_grids, CFG_DIM, m_gpt*CFG::IntVect::Unit);
  CFG::LevelData<CFG::FArrayBox> cfg_mcoords(cfg_grids, CFG_DIM, m_gpt*CFG::IntVect::Unit);

  for (CFG::DataIterator dit(cfg_grids.dataIterator()); dit.ok(); ++dit) {

    const CFG::Box& bx(cfg_grids[dit]);
    const int block_number( mag_coord_sys.whichBlock(bx) );
    const CFG::MagBlockCoordSys& block_coord_sys = static_cast<const CFG::MagBlockCoordSys&>
                                                    (*(mag_coord_sys.getCoordSys(block_number)));

    /* get the physical coordinates */
    block_coord_sys.getCellCenteredRealCoords( cfg_rcoords[dit] );
    /* get the mapped coordinates */
    block_coord_sys.getCellCenteredMappedCoords( cfg_mcoords[dit] );
  }

  /* Inject the configuration space variables to phase space variables */
  LevelData<FArrayBox> inj_b_field, inj_b_field_magn, inj_rcoords, inj_mcoords;
  a_phase_geom.injectConfigurationToPhase(cfg_b_field, inj_b_field);
  a_phase_geom.injectConfigurationToPhase(cfg_b_field_magn, inj_b_field_magn);
  a_phase_geom.injectConfigurationToPhase(cfg_rcoords, inj_rcoords);
  a_phase_geom.injectConfigurationToPhase(cfg_mcoords, inj_mcoords);

  /* define the gyroaveraging grids - phase space grids collapsed in vpar direction */
  a_phase_geom.getConfigurationPlusMuGrids( m_grids, -1 );

  /* define and allocate the gyroaveraging operator */
  m_gyroavg_op.define(m_grids, 1);
  
  for (DataIterator dit(m_grids.dataIterator()); dit.ok(); ++dit) {

    /* create configuration space boxes for physical and mapped coordinates */
    CFG::Box rcoords_bx, mcoords_bx;
    a_phase_geom.projectPhaseToConfiguration(inj_rcoords[dit].box(), rcoords_bx);
    a_phase_geom.projectPhaseToConfiguration(inj_mcoords[dit].box(), mcoords_bx);
    /* create configuration space physical and mapped coordinates*/
    CFG::FArrayBox  rcoords(rcoords_bx, inj_rcoords[dit].nComp()), 
                    mcoords(mcoords_bx, inj_mcoords[dit].nComp());
    a_phase_geom.projectPhaseToConfigurationLocal(inj_rcoords[dit], rcoords);
    a_phase_geom.projectPhaseToConfigurationLocal(inj_mcoords[dit], mcoords);

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

      /* mapped coordinates of this grid point (x,y) */
      CFG::RealVect xi;
      for (int d=0; d<CFG_DIM; d++) xi[d] = mcoords(i_cfg, d);

      /* min and max indices of mapped coordinates box */
      CFG::IntVect  i_max(mcoords_bx.bigEnd()),
                    i_min(mcoords_bx.smallEnd());

#if CFG_DIM == 2
      /* magnetic field and its magnitude at this grid point */
      Real b_magn = b_field_magn(i_cfg, 0);
      Real by = b_field(i_cfg,1);

      /* initialize the gyroavg operator at this grid point */
      Stencil& op = op_fab(iv, 0);
      op.clear();

      /* mu */
      Real mu = iv[MU_DIR] * dmu;

      /* calculate gyroradius for this mu */
      Real rho = computeGyroradius(larmor, a_mass, a_charge, b_magn, mu);

      /* calculate major and minor radii of the elliptical gyro-orbit */
      Real rho_a, rho_b;
      if (b_magn == 0.0) {
        rho_a = rho_b = 0.0;
      } else {
        rho_a = rho;
        rho_b = rho * (by/b_magn);
      }

      /* create interpolation points */
      for (int n = 0; n < m_npts_interp; n++) {
        
        /* compute angle */
        Real alpha = n * delta_alpha;

        /* compute coordinates of interpolation points */
        CFG::RealVect x_pts;
        x_pts[0] = x[0] + rho_a * cos(alpha);
        x_pts[1] = x[1] + rho_b * sin(alpha);

        /* weight of this point */
        Real pt_weight = 1.0 / ((Real)m_npts_interp);

        /* find this point in the grid - i.e. locate the (i,j) index of the 
         * left-bottom grid point of the cell that contains this point*/

        /* find its mapped coordinates */
        const PhaseBlockCoordSys& block_coord_sys(a_phase_geom.getBlockCoordSys(m_grids[dit]));
        RealVect x_pts_inj = a_phase_geom.config_inject(x_pts);
        RealVect xi_pts_inj = block_coord_sys.mappedCoord(x_pts_inj);
        CFG::RealVect xi_pts = a_phase_geom.config_restrict(xi_pts_inj);
        
        /* initialize to current grid point */
        CFG::IntVect i_pts(i_cfg);

        /* along each configuration space dimension, find the grid point index */
        for (int d=0; d < CFG_DIM; d++) {
          if (xi_pts[d] > xi[d]) {
            while (mcoords(i_pts+CFG::BASISV(d),d) < xi_pts[d]) {
              i_pts[d]++;
              if (i_pts[d] == i_max[d]) {
                std::cout << "Rank " << procID() << ", dim = " << d
                          << ", Error in GyroaverageOperator::define(): "
                          << "search for interpolation point reached end of box. "
                          << "Consider increasing number of ghost points for gyroaverage operator.\n";
                MayDay::Error("");
              }
            }
          } else {
            while (mcoords(i_pts,d) > xi_pts[d]) {
              i_pts[d]--;
              if (i_pts[d] < i_min[d]) {
                std::cout << "Rank " << procID() << ", dim = " << d
                          << ", Error in GyroaverageOperator::define(): "
                          << "search for interpolation point reached end of box. "
                          << "Consider increasing number of ghost points for gyroaverage operator.\n";
                MayDay::Error("");
              }
            }
          }
          CH_assert( (mcoords(i_pts,d) - xi_pts[d]) * (mcoords(i_pts+CFG::BASISV(d),d) - xi_pts[d]) <= 0 );

          /* check if this point is within reasonable distance of the grid point */
          if (abs(i_pts[d]-i_cfg[d]) > m_gpt) {
            MayDay::Error("Error in GyroaverageOperator::define() - interpolation point too far away!");
          }
        }

        /* get the stencil to interpolate a variable at this point from grid data */
        Stencil op_pt;
        op_pt.clear();
        getInterpStencil(op_pt, rcoords, i_pts, x_pts);

        /* scale interp weights by the weight of this point 
         * and push to the overall averaging operator */
        for (int k = 0; k < op_pt.size(); k++) {
          op_pt[k].second *= pt_weight;
          op.push_back(op_pt[k]);
        }

      }

#else
      MayDay::Error("Gyroaveraging operator not yet implemented for CFG_DIM=3.");
#endif

    }

  }

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

void GyroaverageOperator::applyOp(FArrayBox&               a_phi_bar,
                                  const FArrayBox&         a_phi,
                                  const BaseFab<Stencil>&  a_op) const
{
  CH_assert(isDefined());
  CH_assert(a_phi_bar.nComp() == a_phi.nComp());
  
  const Box& phi_bar_bx = a_phi_bar.box();
  const Box& op_bx = a_op.box();
  const Box& phi_bx = a_phi.box();
  CH_assert(op_bx == phi_bar_bx);

  const int ncomp = a_phi.nComp();

  /* project a_phi to a configuration space variable */
  CFG::Box cfg_bx;
  m_phase_geom->projectPhaseToConfiguration(phi_bx, cfg_bx);
  CFG::FArrayBox phi(cfg_bx, a_phi.nComp());
  m_phase_geom->projectPhaseToConfigurationLocal(a_phi, phi);

  BoxIterator bit(phi_bar_bx);
  for (bit.begin(); bit.ok(); ++bit) {

    const IntVect iv(bit());
    const Stencil& op = a_op(iv, 0);

    std::vector<Real> phi_avg;
    evalStencil(phi_avg, phi, op);

    for (int v=0; v<ncomp; v++) {
      a_phi_bar(iv,v) = phi_avg[v];
    }
  }

  return;
}

void GyroaverageOperator::applyOp(LevelData<FArrayBox>&                 a_phi_bar,
                                  const CFG::LevelData<CFG::FArrayBox>& a_phi ) const
{
  CH_assert(isDefined());
  CH_assert(a_phi.isDefined());

  if (a_phi_bar.isDefined()) {
    a_phi_bar.clear();
  }
  /* define it */
  a_phi_bar.define(m_grids, a_phi.nComp());

  /* inject phi into phase space */
  LevelData<FArrayBox> injected_phi;
  m_phase_geom->injectConfigurationToPhase(a_phi, injected_phi);

  /* create a version with ghost points in the configuration dimensions */
  IntVect ghost_vec(IntVect::Zero); for (int d=0; d<CFG_DIM; d++) ghost_vec[d] = m_gpt;
  LevelData<FArrayBox> phi_wghosts;
  phi_wghosts.define(injected_phi.disjointBoxLayout(), a_phi.nComp(), ghost_vec);

  /* copy phi to the version with ghost points and then fill the ghost points */
  /* interior */
  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    phi_wghosts[dit].setVal(0.0);
    phi_wghosts[dit].copy(injected_phi[dit], injected_phi[dit].box());
  }
  /* codim-1 boundaries */
  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    for (int dir=0; dir<CFG_DIM; dir++) {
      fillGhostCells(phi_wghosts[dit], injected_phi[dit].box(), dir);
    }
  }
  phi_wghosts.exchange();
  /* higher codim boundaries */
  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    Box bx_int(injected_phi[dit].box());
    const Box& bx_phi(phi_wghosts[dit].box());
    for (int dir=0; dir<CFG_DIM-1; dir++) {
      bx_int.growLo(dir, (bx_int.smallEnd(dir)-bx_phi.smallEnd(dir)) );
      bx_int.growHi(dir, (bx_phi.bigEnd(dir)-bx_int.bigEnd(dir)) );
      fillGhostCells(phi_wghosts[dit], bx_int, dir+1);
    }
  }
  phi_wghosts.exchange();
  /* done - hopefully! */

  /* now apply the gyroaveraging */
  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    applyOp(a_phi_bar[dit], phi_wghosts[dit], m_gyroavg_op[dit]);
  }

  return;
}

void GyroaverageOperator::getInterpStencil( Stencil&              a_stencil,
                                            const CFG::FArrayBox& a_coords,
                                            const CFG::IntVect&   a_i0,
                                            const CFG::RealVect&  a_x )
{
  a_stencil.clear();

#if CFG_DIM==2

  /* 
   * Bilinear interpolation using 4 arbitrary points in 2D space 
  */

  /* set the (i,j) grid coordinates of the 4 points */
  CFG::IntVect  i0(a_i0),
                i1(a_i0[0]+1,a_i0[1]),
                i2(a_i0[0],a_i0[1]+1),
                i3(a_i0[0]+1,a_i0[1]+1);

  /* extract the (x,y) coords of the interpolation point */
  Real  xc = a_x[0],
        yc = a_x[1];

  /* extract the (x,y) coords of the four points */
  Real  x0 = a_coords(i0,0),
        y0 = a_coords(i0,1);
  Real  x1 = a_coords(i1,0),
        y1 = a_coords(i1,1);
  Real  x2 = a_coords(i2,0),
        y2 = a_coords(i2,1);
  Real  x3 = a_coords(i3,0),
        y3 = a_coords(i3,1);

  /* bilinear interp coefficients - computed using Mathematica */
  Real c0, c1, c2, c3;
  c0 = (x1*(-(xc*(y2 - y3)*(y1 - yc)) + x3*(y1 - y3)*(y2 - yc) - x2*(y1 - y2)*(y3 - yc)) - 
       x2*x3*(y2 - y3)*(y1 - yc) + x2*xc*(y1 - y3)*(y2 - yc) - x3*xc*(y1 - y2)*(y3 - yc))/
          (x0*(x3*(y1 - y2)*(y0 - y3) - x2*(y0 - y2)*(y1 - y3) + x1*(y0 - y1)*(y2 - y3)) + 
               x1*x2*(y1 - y2)*(y0 - y3) - x1*x3*(y0 - y2)*(y1 - y3) + x2*x3*(y0 - y1)*(y2 - y3));
  c1 = (x0*(xc*(y2 - y3)*(y0 - yc) - x3*(y0 - y3)*(y2 - yc) + x2*(y0 - y2)*(y3 - yc)) + 
       x2*x3*(y2 - y3)*(y0 - yc) - x2*xc*(y0 - y3)*(y2 - yc) + x3*xc*(y0 - y2)*(y3 - yc))/
          (x0*(x3*(y1 - y2)*(y0 - y3) - x2*(y0 - y2)*(y1 - y3) + x1*(y0 - y1)*(y2 - y3)) + 
               x1*x2*(y1 - y2)*(y0 - y3) - x1*x3*(y0 - y2)*(y1 - y3) + x2*x3*(y0 - y1)*(y2 - y3));
  c2 = (x0*(-(xc*(y1 - y3)*(y0 - yc)) + x3*(y0 - y3)*(y1 - yc) - x1*(y0 - y1)*(y3 - yc)) - 
       x1*x3*(y1 - y3)*(y0 - yc) + x1*xc*(y0 - y3)*(y1 - yc) - x3*xc*(y0 - y1)*(y3 - yc))/
          (x0*(x3*(y1 - y2)*(y0 - y3) - x2*(y0 - y2)*(y1 - y3) + x1*(y0 - y1)*(y2 - y3)) + 
               x1*x2*(y1 - y2)*(y0 - y3) - x1*x3*(y0 - y2)*(y1 - y3) + x2*x3*(y0 - y1)*(y2 - y3));
  c3 = (x0*(xc*(y1 - y2)*(y0 - yc) - x2*(y0 - y2)*(y1 - yc) + x1*(y0 - y1)*(y2 - yc)) + 
       x1*x2*(y1 - y2)*(y0 - yc) - x1*xc*(y0 - y2)*(y1 - yc) + x2*xc*(y0 - y1)*(y2 - yc))/
          (x0*(x3*(y1 - y2)*(y0 - y3) - x2*(y0 - y2)*(y1 - y3) + x1*(y0 - y1)*(y2 - y3)) + 
               x1*x2*(y1 - y2)*(y0 - y3) - x1*x3*(y0 - y2)*(y1 - y3) + x2*x3*(y0 - y1)*(y2 - y3));

  /* push these coefficients and coordinates to the stencil */
  a_stencil.push_back(StencilPoint(i0,c0));
  a_stencil.push_back(StencilPoint(i1,c1));
  a_stencil.push_back(StencilPoint(i2,c2));
  a_stencil.push_back(StencilPoint(i3,c3));

#else
  MayDay::Error("Gyroaveraging operator not yet implemented for CFG_DIM=3.");
#endif

  /* done */
  return;
}

#include "NamespaceFooter.H"
