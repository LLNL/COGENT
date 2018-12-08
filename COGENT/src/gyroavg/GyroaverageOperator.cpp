#include "GyroaverageOperator.H"
#include "PhaseGeom.H"

#include "NamespaceHeader.H"

static void fillGhostCells( CFG::FArrayBox& a_phi,
                            const CFG::Box& a_bx_int,
                            const int       a_dir )
{
  const CFG::Box& phi_bx = a_phi.box();

  for (int side=-1; side<2; side+=2) {

    int ii;
    CFG::Box bdrybox;

    if (side == -1) {
      int n_gpt = a_bx_int.smallEnd(a_dir) - phi_bx.smallEnd(a_dir);
      ii = a_bx_int.smallEnd(a_dir);
      bdrybox = adjCellLo(a_bx_int, a_dir, n_gpt);
    } else {
      int n_gpt = phi_bx.bigEnd(a_dir) - a_bx_int.bigEnd(a_dir);
      ii = a_bx_int.bigEnd(a_dir);
      bdrybox = adjCellHi(a_bx_int, a_dir, n_gpt);
    }

    CFG::BoxIterator bit(bdrybox);
    for (bit.begin(); bit.ok(); ++bit) {

      CFG::IntVect iv(bit());
      int j = (Real) (side < 0 ? ii-iv[a_dir] : iv[a_dir]-ii );

      CFG::IntVect i_int_0(iv); i_int_0[a_dir] = ii;
      CFG::IntVect i_int_1(iv); i_int_1[a_dir] = ii - side;
      CFG::IntVect i_int_2(iv); i_int_2[a_dir] = ii - 2*side;
      CFG::IntVect i_int_3(iv); i_int_3[a_dir] = ii - 3*side;

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

  /* get the velocity coordinate system */
  const VEL::VelCoordSys&   vel_coords    = a_phase_geom.velSpaceCoordSys();
  const VEL::ProblemDomain& vel_domain    = vel_coords.domain();
  const VEL::RealVect&      vel_dx        = vel_coords.dx();
  const VEL::Box&           vel_domain_bx = vel_domain.domainBox();

  /* get number of grids points along mu and dmu */
  m_n_mu = vel_domain_bx.size(_MU_);
  const Real dmu = vel_dx[_MU_];

  /* create a vector of mu */
  std::vector<Real> mu_vec(m_n_mu, 0.0);
  for (int i=0; i<m_n_mu; i++) {
    mu_vec[i] = i * dmu;
  }

  /* Larmor number */
  Real larmor = a_phase_geom.larmorNumber();

  const CFG::MagGeom& mag_geom(a_phase_geom.magGeom());
  const CFG::MultiBlockCoordSys& mag_coord_sys( *(mag_geom.coordSysPtr()) );
  const CFG::DisjointBoxLayout& grids(mag_geom.grids());

  /* get magnetic field magnitude data */
  const CFG::LevelData<CFG::FArrayBox>& b_field_magn_level(mag_geom.getCCBFieldMag());
  const CFG::LevelData<CFG::FArrayBox>& b_field_level(mag_geom.getCCBField());

  /* define and allocate the gyroaveraging operator */
  m_gyroavg_op.define(grids, m_n_mu);
  
  for (CFG::DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {

    const CFG::Box& bx(grids[dit]);
    const int block_number( mag_coord_sys.whichBlock(bx) );
    const CFG::MagBlockCoordSys& block_coord_sys 
            = static_cast<const CFG::MagBlockCoordSys&>
              (*(mag_coord_sys.getCoordSys(block_number)));

    /* create a box with many ghost points */
    CFG::Box bx_wghosts = grow(bx, m_gpt);
    CFG::IntVect i_max(bx_wghosts.bigEnd()),
                 i_min(bx_wghosts.smallEnd());

    /* get the physical coordindates */
    CFG::FArrayBox rcoords(bx_wghosts, CFG_DIM);
    block_coord_sys.getCellCenteredRealCoords( rcoords );

    /* get the mapped coordindates */
    CFG::FArrayBox mcoords(bx_wghosts, CFG_DIM);
    block_coord_sys.getCellCenteredMappedCoords( mcoords );

    const CFG::FArrayBox& b_field_magn_fab = b_field_magn_level[dit];
    const CFG::FArrayBox& b_field_fab = b_field_level[dit];
    CFG::BaseFab<Stencil>& op_fab = m_gyroavg_op[dit];

    CFG::BoxIterator bit(bx);
    for (bit.begin(); bit.ok(); ++bit) {

      const CFG::IntVect iv(bit());

      /* physical coordinates of this grid point (x,y) */
      CFG::RealVect x;
      for (int d=0; d<CFG_DIM; d++) x[d] = rcoords(iv, d);

      /* mapped coordinates of this grid point (x,y) */
      CFG::RealVect xi;
      for (int d=0; d<CFG_DIM; d++) xi[d] = mcoords(iv, d);

#if CFG_DIM == 2
      /* magnetic field and its magnitude at this grid point */
      Real b_magn = b_field_magn_fab.get(iv, 0);
      Real by = b_field_fab.get(iv,1);

      const CFG::IntVect i(iv);
      for (int l=0; l<m_n_mu; l++) {

        Stencil& op = op_fab(iv, l);
        op.clear();

        /* get mu */
        Real mu = mu_vec[l];

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
          CFG::RealVect xi_pts = block_coord_sys.mappedCoord(x_pts);
          CFG::IntVect i_pts(i);
          for (int dim=0; dim < CFG_DIM; dim++) {
            if (xi_pts[dim] > xi[dim]) {
              while (mcoords(i_pts+CFG::BASISV(dim),dim) < xi_pts[dim]) {
                i_pts[dim]++;
                if (i_pts[dim] == i_max[dim]) {
                  std::cout << "Rank " << procID() << ", dim = " << dim
                            << ", Error in GyroaverageOperator::define(): "
                            << "search for interpolation point reached end of box. "
                            << "Consider increasing number of ghost points for gyroaverage operator.\n";
                  MayDay::Error("");
                }
              }
            } else {
              while (mcoords(i_pts,dim) > xi_pts[dim]) {
                i_pts[dim]--;
                if (i_pts[dim] < i_min[dim]) {
                  std::cout << "Rank " << procID() << ", dim = " << dim
                            << ", Error in GyroaverageOperator::define(): "
                            << "search for interpolation point reached end of box. "
                            << "Consider increasing number of ghost points for gyroaverage operator.\n";
                  MayDay::Error("");
                }
              }
            }
            CH_assert( (mcoords(i_pts,dim) - xi_pts[dim]) * (mcoords(i_pts+CFG::BASISV(dim),dim) - xi_pts[dim]) <= 0 );

            /* check if this point is within reasonable distance of the grid point */
            if (abs(i_pts[dim]-i[dim]) > m_gpt) {
              MayDay::Error("Error in GyroaverageOperator::define() - interpolation point too far away!");
            }
          }

          /* get the stencil to interpolate a variable at this point from grid data */
          Stencil op_pt;
          op_pt.clear();
          getInterpStencil( op_pt, rcoords, i_pts, x_pts );

          /* scale interp weights by the weight of this point 
           * and push to the overall averaging operator */
          for (int k = 0; k < op_pt.size(); k++) {
            op_pt[k].second *= pt_weight;
            op.push_back(op_pt[k]);
          }

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
  m_n_mu = a_gyroavg.m_n_mu;
  m_phase_geom = a_gyroavg.m_phase_geom;

  const CFG::LevelData< CFG::BaseFab<Stencil> >& copy_from = a_gyroavg.op();
  const CFG::DisjointBoxLayout& grids = copy_from.disjointBoxLayout();
  const CFG::IntVect& gv = copy_from.ghostVect();
  int ncomp = copy_from.nComp();

  m_gyroavg_op.define(grids, ncomp, gv);
  copy_from.copyTo(m_gyroavg_op);

  m_is_defined = true;
}

void GyroaverageOperator::applyOp(CFG::FArrayBox&               a_phi_bar,
                                  const CFG::FArrayBox&         a_phi,
                                  const CFG::BaseFab<Stencil>&  a_op) const
{
  CH_assert(isDefined());
  CH_assert(a_phi_bar.nComp() == m_n_mu*a_phi.nComp());
  
  const CFG::Box& phi_bx = a_phi.box();
  const CFG::Box& phi_bar_bx = a_phi_bar.box();
  const CFG::Box& op_bx = a_op.box();
  CH_assert(phi_bx.contains(phi_bar_bx));
  CH_assert(op_bx.contains(phi_bar_bx));

  const int ncomp = a_phi.nComp();

  CFG::BoxIterator bit(phi_bar_bx);
  for (bit.begin(); bit.ok(); ++bit) {

    const CFG::IntVect iv(bit());

    for (int l=0; l<m_n_mu; l++) {

      const Stencil& op = a_op(iv, l);
      std::vector<Real> phi_avg;
      evalStencil(phi_avg, a_phi, op);

      for (int v=0; v<ncomp; v++) {
        a_phi_bar(iv,v*m_n_mu+l) = phi_avg[v];
      }
    }
  }

  return;
}

void GyroaverageOperator::applyOp(CFG::LevelData<CFG::FArrayBox>&               a_phi_bar,
                                  const CFG::LevelData<CFG::FArrayBox>&         a_phi ) const
{
  CH_assert(isDefined());
  CH_assert(a_phi.isDefined());

  const CFG::DisjointBoxLayout& phi_grids = a_phi.disjointBoxLayout();
  const CFG::DisjointBoxLayout& op_grids  = m_gyroavg_op.disjointBoxLayout();
  CH_assert(phi_grids == op_grids);

  if (a_phi_bar.isDefined()) {
    const CFG::DisjointBoxLayout& phi_bar_grids = a_phi_bar.disjointBoxLayout();
    CH_assert(phi_grids == phi_bar_grids);
    CH_assert(a_phi_bar.ghostVect() == CFG::IntVect::Zero);
  } else {
    a_phi_bar.define(phi_grids, m_n_mu*a_phi.nComp());
  }

  CFG::LevelData<CFG::FArrayBox> phi_wghosts;
  phi_wghosts.define(phi_grids, a_phi.nComp(), m_gpt*CFG::IntVect::Unit);

  for (CFG::DataIterator dit(phi_grids); dit.ok(); ++dit) {

    phi_wghosts[dit].setVal(0.0);
    phi_wghosts[dit].copy(a_phi[dit], phi_grids[dit]);

  }
  
  for (CFG::DataIterator dit(phi_grids); dit.ok(); ++dit) {

    CFG::Box bx_int(phi_grids[dit]);
    for (int dir=0; dir<CFG_DIM; dir++) {
      fillGhostCells(phi_wghosts[dit], bx_int, dir);
    }

  }

  phi_wghosts.exchange();

  for (CFG::DataIterator dit(phi_grids); dit.ok(); ++dit) {

    CFG::Box bx_int(phi_grids[dit]);
    const CFG::Box& bx_phi(phi_wghosts[dit].box());

    for (int dir=0; dir<CFG_DIM-1; dir++) {
      bx_int.growLo(dir, (bx_int.smallEnd(dir)-bx_phi.smallEnd(dir)) );
      bx_int.growHi(dir, (bx_phi.bigEnd(dir)-bx_int.bigEnd(dir)) );
      fillGhostCells(phi_wghosts[dit], bx_int, dir+1);
    }

  }

  phi_wghosts.exchange();

  m_phase_geom->plotConfigurationData("phi_wg", phi_wghosts, 0.0);

  for (CFG::DataIterator dit(phi_grids); dit.ok(); ++dit) {
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
