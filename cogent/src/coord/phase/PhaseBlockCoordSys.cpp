#include "PhaseBlockCoordSys.H"
#include "PhaseBlockCoordSysF_F.H"
#include "FourthOrderUtil.H"

#include "Directions.H"

#undef CH_SPACEDIM
#include "Injection.H.transdim"
#include "Slicing.H.transdim"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

using namespace CH_MultiDim;



PhaseBlockCoordSys::PhaseBlockCoordSys( const CFG::MagBlockCoordSys& a_configuration_space_coord_sys,
                                        const VEL::VelCoordSys&      a_velocity_space_coord_sys,
                                        const ProblemDomain&         a_domain )
   : m_domain(a_domain),
     m_configuration_coords(a_configuration_space_coord_sys),
     m_velocity_coords(a_velocity_space_coord_sys)
{
   //  Construct the phase space mesh size from the configuration and velocity space values.
   //  Compute the mapped cell volumes and face areas.
   const CFG::RealVect & configuration_dx = m_configuration_coords.dx();
   const VEL::RealVect & velocity_dx = m_velocity_coords.dx();

   m_mapped_cell_volume = 1.;
   for (int i=0; i<CFG_DIM; i++) {
      m_dx[i] = configuration_dx[i];
      m_mapped_cell_volume *= m_dx[i];
   }
   for (int i=VPARALLEL_DIR; i<=MU_DIR; i++) {
      m_dx[i] = velocity_dx[i-VPARALLEL_DIR];
      m_mapped_cell_volume *= m_dx[i];
   }

   for (int dir=0; dir<SpaceDim; ++dir) {
      m_mapped_face_area[dir] = 1.;
      for (int tdir=0; tdir<SpaceDim; ++tdir) {
         if (tdir != dir) m_mapped_face_area[dir] *= m_dx[tdir];
      }
   }
}



PhaseBlockCoordSys::~PhaseBlockCoordSys()
{
}



void
PhaseBlockCoordSys::getN( FluxBox&   a_N,
                          const Box& a_box ) const
{
#if 0
   bool found_box = false;

   DataIterator dit = m_grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

      if ((*m_configuration_metrics)[dit].box().contains(configurationFlatten(m_grids[dit],a_box))
          && (*m_velocity_metrics)[dit].box().contains(velocityFlatten(m_grids[dit],a_box)) ) {

         a_N.setVal(0.);

         for (int dir=0; dir<SpaceDim; dir++) {
            FArrayBox& thisNDir = a_N[dir];

            if (dir < CFG_DIM) {
               const FArrayBox& faceMetricDir = (*m_configuration_metrics)[dit][dir];
               int ncomp = faceMetricDir.nComp();

               for (int i=0; i<ncomp; ++i) {
                  for (int j=0; j<ncomp; ++j) {
                     thisNDir.copy(faceMetricDir, a_box, this->getNcomponent(i,j),
                                   a_box, m_configuration_coords.getNcomponent(i,j), 1);
                  }
               }

               // Apply area factors to get phase space face integrals
               const FArrayBox& thisvol = (*m_velocity_volumes)[dit];
               FORT_MULT_VEL(CHF_BOX(a_box),
                             CHF_FRA1(thisvol,0),
                             CHF_FRA(thisNDir));
            }
            else {
               const FArrayBox& faceMetricDir = (*m_velocity_metrics)[dit][dir];
               int ncomp = faceMetricDir.nComp();

               for (int i=0; i<ncomp; ++i) {
                  for (int j=0; j<ncomp; ++j) {
                     thisNDir.copy(faceMetricDir, a_box, this->getNcomponent(i+CFG_DIM,j+CFG_DIM),
                                   a_box, m_velocity_coords.getNcomponent(i,j), 1);
                  }
               }

               // Apply area factors to get phase space face integrals
               const FArrayBox& thisvol = (*m_configuration_volumes)[dit];
               FORT_MULT_CFG(CHF_BOX(a_box),
                             CHF_FRA1(thisvol,0),
                             CHF_FRA(thisNDir));
            }
         } // end loop over face directions

         found_box = true;
      }
   }

   if (!found_box) {
      MayDay::Error("PhaseBlockCoordSys::getN(): Could not fill input array");
   }
#else
   // MultiblockLevelGeom wants to call this function to compute cell volumes, but it doesn't
   // matter what we return since we are overriding the cell volume calculation anyway.
   a_N.setVal(0., a_box);
#endif
}



void
PhaseBlockCoordSys::cellVol( FArrayBox&     a_vol,
                             const FluxBox& a_N,
                             const Box&     a_box ) const
{
#if 0
   bool found_box = false;

   DataIterator dit = m_grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

      if ((*m_configuration_volumes)[dit].box().contains(configurationFlatten(m_grids[dit],a_box))
          && (*m_velocity_volumes)[dit].box().contains(velocityFlatten(m_grids[dit],a_box)) ) {

         a_vol.setVal(1.);

         FORT_MULT_CFG(CHF_BOX(a_box),
                       CHF_FRA1((*m_configuration_volumes)[dit],0),
                       CHF_FRA(a_vol));
         FORT_MULT_VEL(CHF_BOX(a_box),
                       CHF_FRA1((*m_velocity_volumes)[dit],0),
                       CHF_FRA(a_vol));
         found_box = true;
      }
   }

   if (!found_box) {
      MayDay::Error("PhaseBlockCoordSys::cellVol(): Could not fill input array");
   }
#else
   // MultiblockLevelGeom calls this function to compute cell volumes, but never
   // uses them
   a_vol.setVal(0.);
#endif
}



void
PhaseBlockCoordSys::getAvgJ( FArrayBox& a_J,
                             const Box& a_box ) const
{
  CH_assert(a_J.box().contains(a_box));

  IntVect lo = a_box.smallEnd();
  IntVect hi = a_box.bigEnd();

  CFG::Box cfg_box(config_restrict(lo),config_restrict(hi));
  CFG::FArrayBox cfg_J(cfg_box,1);
  m_configuration_coords.getAvgJ(cfg_J, cfg_box);

  VEL::Box vel_box(vel_restrict(lo),vel_restrict(hi));
  VEL::FArrayBox vel_J(vel_box,1);
  m_velocity_coords.getAvgJ(vel_J, vel_box);

  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit) {
    IntVect iv = bit();
    CFG::IntVect cfg_iv = config_restrict(iv);
    VEL::IntVect vel_iv = vel_restrict(iv);

    a_J(iv) = cfg_J(cfg_iv) * vel_J(vel_iv);
  }
}


void
PhaseBlockCoordSys::pointwiseJ( FArrayBox& a_J,
                                const Box& a_box ) const
{
  CH_assert(a_J.box().contains(a_box));

  IntVect lo(a_box.smallEnd());
  IntVect hi(a_box.bigEnd());

  CFG::Box cfg_box(config_restrict(lo),config_restrict(hi));
  CFG::FArrayBox cfg_J(cfg_box,1);
  CFG::FArrayBox cfg_Xi(cfg_box,CFG_DIM);
  m_configuration_coords.getCellCenteredMappedCoords(cfg_Xi);
  m_configuration_coords.pointwiseJ(cfg_J, cfg_Xi, cfg_box);

  VEL::Box vel_box(vel_restrict(lo),vel_restrict(hi));
  VEL::FArrayBox vel_J(vel_box,1);
  VEL::FArrayBox vel_Xi(vel_box,VEL_DIM);
  const RealVect offset(0.5 * RealVect::Unit);
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
     IntVect iv(bit());
     VEL::RealVect mappedLoc(vel_restrict(m_dx * iv + offset));
     VEL::IntVect vel_iv(vel_restrict(iv));
     for (int dir(0); dir<VEL_DIM; ++dir) {
        vel_Xi(vel_iv,dir) = mappedLoc[dir];
     }
  }
  m_velocity_coords.pointwiseJ(vel_J, vel_Xi, vel_box);

  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
     IntVect iv(bit());
     CFG::IntVect cfg_iv(config_restrict(iv));
     VEL::IntVect vel_iv(vel_restrict(iv));
     a_J(iv) = cfg_J(cfg_iv) * vel_J(vel_iv);
  }
}



void
PhaseBlockCoordSys::projectOntoFluxSurface( FArrayBox& a_data ) const
{
   CH_assert(a_data.nComp() == SpaceDim);

   const Box& box = a_data.box();
   CFG::Box cfg_box(config_restrict(box.smallEnd()),config_restrict(box.bigEnd()));
   CFG::FArrayBox cfg_unit_tangent(cfg_box,CFG_DIM);
   m_configuration_coords.computeFluxSurfaceUnitTangent(cfg_unit_tangent);

   FArrayBox unit_tangent;
   injectConfigurationToPhase(cfg_unit_tangent, unit_tangent);

   FORT_PHASE_BLOCK_PROJECT(CHF_BOX(a_data.box()),
                            CHF_CONST_FRA(unit_tangent),
                            CHF_FRA(a_data));
}



void
PhaseBlockCoordSys::projectOntoFluxSurface( FluxBox& a_data ) const
{
   CH_assert(a_data.nComp() == SpaceDim);

   const Box& box = a_data.box();
   CFG::Box cfg_box(config_restrict(box.smallEnd()),config_restrict(box.bigEnd()));
   CFG::FluxBox cfg_unit_tangent(cfg_box,CFG_DIM);
   m_configuration_coords.computeFluxSurfaceUnitTangent(cfg_unit_tangent);

   FluxBox unit_tangent;
   injectConfigurationToPhase(cfg_unit_tangent, unit_tangent);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_data = a_data[dir];
      FORT_PHASE_BLOCK_PROJECT(CHF_BOX(this_data.box()),
                               CHF_CONST_FRA(unit_tangent[dir]),
                               CHF_FRA(this_data));
   }
}



void
PhaseBlockCoordSys::projectOntoFluxSurfaceNormal( FArrayBox& a_data ) const
{
   CH_assert(a_data.nComp() == SpaceDim);

   const Box& box = a_data.box();
   CFG::Box cfg_box(config_restrict(box.smallEnd()),config_restrict(box.bigEnd()));
   CFG::FArrayBox cfg_unit_normal(cfg_box,CFG_DIM);
   m_configuration_coords.computeFluxSurfaceUnitNormal(cfg_unit_normal);

   FArrayBox unit_normal;
   injectConfigurationToPhase(cfg_unit_normal, unit_normal);

   FORT_PHASE_BLOCK_PROJECT(CHF_BOX(a_data.box()),
                            CHF_CONST_FRA(unit_normal),
                            CHF_FRA(a_data));
}



void
PhaseBlockCoordSys::projectOntoFluxSurfaceNormal( FluxBox& a_data ) const
{
   CH_assert(a_data.nComp() == SpaceDim);

   const Box& box = a_data.box();
   CFG::Box cfg_box(config_restrict(box.smallEnd()),config_restrict(box.bigEnd()));
   CFG::FluxBox cfg_unit_normal(cfg_box,CFG_DIM);
   m_configuration_coords.computeFluxSurfaceUnitNormal(cfg_unit_normal);

   FluxBox unit_normal;
   injectConfigurationToPhase(cfg_unit_normal, unit_normal);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_data = a_data[dir];
      FORT_PHASE_BLOCK_PROJECT(CHF_BOX(this_data.box()),
                               CHF_CONST_FRA(unit_normal[dir]),
                               CHF_FRA(this_data));
   }
}


void
PhaseBlockCoordSys::computePsiThetaProjections( FArrayBox& a_data ) const
{
   const Box& box = a_data.box();
   CFG::Box cfg_box(config_restrict(box.smallEnd()),config_restrict(box.bigEnd()));

   CFG::FArrayBox cfg_unit_normal(cfg_box,CFG_DIM);
   m_configuration_coords.computeFluxSurfaceUnitNormal(cfg_unit_normal);
   FArrayBox unit_normal;
   injectConfigurationToPhase(cfg_unit_normal, unit_normal);

   CFG::FArrayBox cfg_unit_tangent(cfg_box,CFG_DIM);
   m_configuration_coords.computeFluxSurfaceUnitTangent(cfg_unit_tangent);
   FArrayBox unit_tangent;
   injectConfigurationToPhase(cfg_unit_tangent, unit_tangent);

   FORT_PHASE_BLOCK_PSITHETA_PROJECTIONS(CHF_BOX(a_data.box()),
                                         CHF_CONST_FRA(unit_normal),
                                         CHF_CONST_FRA(unit_tangent),
                                         CHF_FRA(a_data));
}



void
PhaseBlockCoordSys::computePsiThetaProjections( FluxBox& a_data ) const
{
   const Box& box = a_data.box();
   CFG::Box cfg_box(config_restrict(box.smallEnd()),config_restrict(box.bigEnd()));

   CFG::FluxBox cfg_unit_normal(cfg_box,CFG_DIM);
   m_configuration_coords.computeFluxSurfaceUnitNormal(cfg_unit_normal);
   FluxBox unit_normal;
   injectConfigurationToPhase(cfg_unit_normal, unit_normal);

   CFG::FluxBox cfg_unit_tangent(cfg_box,CFG_DIM);
   m_configuration_coords.computeFluxSurfaceUnitTangent(cfg_unit_tangent);
   FluxBox unit_tangent;
   injectConfigurationToPhase(cfg_unit_tangent, unit_tangent);

   for (int dir=0; dir<SpaceDim; ++dir) {
      FArrayBox& this_data = a_data[dir];
      FORT_PHASE_BLOCK_PSITHETA_PROJECTIONS(CHF_BOX(this_data.box()),
                                            CHF_CONST_FRA(unit_normal[dir]),
                                            CHF_CONST_FRA(unit_tangent[dir]),
                                            CHF_FRA(this_data));
   }
}


CFG::RealVect
PhaseBlockCoordSys::config_restrict( const RealVect& a_x ) const
{
   CFG::RealVect restrict_x;

   for (int i=0; i<CFG_DIM; i++) {
      restrict_x[i] = a_x[i];
   }

   return restrict_x;
}



CFG::IntVect
PhaseBlockCoordSys::config_restrict( const IntVect& a_x ) const
{
   CFG::IntVect restrict_x;

   for (int i=0; i<CFG_DIM; i++) {
      restrict_x[i] = a_x[i];
   }

   return restrict_x;
}



VEL::IntVect
PhaseBlockCoordSys::vel_restrict( const IntVect& a_x ) const
{
   VEL::IntVect restrict_x;

   for (int i=VPARALLEL_DIR; i<=MU_DIR; i++) {
      restrict_x[i-VPARALLEL_DIR] = a_x[i];
   }

   return restrict_x;
}



VEL::RealVect
PhaseBlockCoordSys::vel_restrict( const RealVect& a_x ) const
{
   VEL::RealVect restrict_x;

   for (int i=VPARALLEL_DIR; i<=MU_DIR; i++) {
      restrict_x[i-VPARALLEL_DIR] = a_x[i];
   }

   return restrict_x;
}



RealVect
PhaseBlockCoordSys::realCoord( const RealVect& a_Xi ) const
{
   RealVect x;
   {
      CFG::RealVect Xi_restricted( config_restrict( a_Xi ) );
      CFG::RealVect x_restricted( m_configuration_coords.realCoord( Xi_restricted ) );
      for (int i=0; i<CFG_DIM; i++) x[i] = x_restricted[i];
   }
   {
      VEL::RealVect Xi_restricted( vel_restrict( a_Xi ) );
      VEL::RealVect x_restricted( m_velocity_coords.realCoord( Xi_restricted ) );
      for (int i=VPARALLEL_DIR; i<=MU_DIR; i++) x[i] = x_restricted[i-VPARALLEL_DIR];
   }

   return x;
}



RealVect
PhaseBlockCoordSys::mappedCoord( const RealVect& a_x ) const
{
   RealVect Xi;
   {
      CFG::RealVect x_restricted( config_restrict( a_x ) );
      CFG::RealVect Xi_restricted( m_configuration_coords.mappedCoord( x_restricted ) );
      for (int i=0; i<CFG_DIM; i++) Xi[i] = Xi_restricted[i];
   }
   {
      VEL::RealVect x_restricted( vel_restrict( a_x ) );
      VEL::RealVect Xi_restricted( m_velocity_coords.mappedCoord( x_restricted ) );
      for (int i=VPARALLEL_DIR; i<=MU_DIR; i++) Xi[i] = Xi_restricted[i-VPARALLEL_DIR];
   }

   return Xi;
}


Real
PhaseBlockCoordSys::dXdXi( const RealVect& a_X,
                           int             a_dirX,
                           int             a_dirXi ) const
{
   Real dXdxi;

   if (a_dirX < CFG_DIM) {
      if (a_dirXi < CFG_DIM) {
         dXdxi = m_configuration_coords.dXdXi(config_restrict(a_X), a_dirX, a_dirXi);
      }
      else {
         dXdxi = 0.;
      }
   }
   else {
      if (a_dirXi < CFG_DIM) {
         dXdxi = 0.;
      }
      else {
         dXdxi = m_velocity_coords.dXdXi(vel_restrict(a_X), a_dirX-CFG_DIM, a_dirXi-CFG_DIM);
      }
   }

   return dXdxi;
}



void
PhaseBlockCoordSys::getVelocityRealCoords( FArrayBox& a_velocityRealCoords ) const
{
   CH_assert(a_velocityRealCoords.nComp() == VEL_DIM);

   // start with assumption of a cell-centred fab
   RealVect offset = 0.5*RealVect::Unit;
   offset *= m_dx;
   for (int dir=0; dir<SpaceDim; dir++) {
      if (a_velocityRealCoords.box().type(dir) == IndexType::NODE) {
         offset[dir] = 0.0;
      }
   }

   BoxIterator bit(a_velocityRealCoords.box());
   for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect mappedLoc = m_dx*iv + offset;
      RealVect realLoc = realCoord(mappedLoc);
      for (int dir=0; dir<VEL_DIM; ++dir) {
         a_velocityRealCoords(iv,dir) = realLoc[CFG_DIM + dir];
      }
   }
}


void
PhaseBlockCoordSys::getCellCenteredRealCoords( FArrayBox& a_x ) const
{
   const Box& box( a_x.box() );
   FArrayBox xi_array( box, a_x.nComp() );
   getCellCenteredMappedCoords( xi_array );

   for (BoxIterator bit( box ); bit.ok(); ++bit) {
      IntVect iv( bit() );
      RealVect xi;
      for (int dir(0); dir<PDIM; ++dir) {
         xi[dir] = xi_array( iv, dir );
      }
      RealVect real_loc( realCoord( xi ) );
      for (int dir(0); dir<PDIM; ++dir) {
         a_x( iv, dir ) = real_loc[dir];
      }
   }
}


void
PhaseBlockCoordSys::getCellCenteredMappedCoords( FArrayBox& a_xi ) const
{
   CH_assert(a_xi.nComp()==PDIM);
   const Box& box( a_xi.box() );
   for (int dir(0); dir<PDIM; dir++) {
      CH_assert( a_xi.box().type(dir)!=IndexType::NODE );
   }

   RealVect shift( 0.5 * RealVect::Unit );
   for (BoxIterator bit( box ); bit.ok(); ++bit) {
      IntVect iv( bit() );
      RealVect mapped_loc( m_dx * (iv + shift) );
      for (int dir(0); dir<PDIM; ++dir) {
         a_xi( iv, dir ) = mapped_loc[dir];
      }
   }
}


void
PhaseBlockCoordSys::getFaceCenteredMappedCoords( FluxBox& a_xi ) const
{
   CH_assert(a_xi.nComp()==PDIM);

   for (int dir=0; dir<SpaceDim; ++dir) {
     FArrayBox& this_xi = a_xi[dir];
     RealVect shift( 0.5 * RealVect::Unit );
     shift[dir] = 0;
     for (BoxIterator bit( a_xi[dir].box() ); bit.ok(); ++bit) {
       IntVect iv( bit() );
       RealVect mapped_loc( m_dx * (iv + shift) );
       for (int dir2(0); dir2<PDIM; ++dir2) {
         this_xi( iv, dir2 ) = mapped_loc[dir2];
       }
     }
   }
}


/*
  Transdimensional utilities
*/


void
PhaseBlockCoordSys::injectConfigurationToPhase( const CFG::FArrayBox& a_src,
                                                FArrayBox&            a_dst ) const
{
   const Box& domainBox = m_domain.domainBox();
   const CP1::SliceSpec slice_vp(VPARALLEL_DIR,domainBox.smallEnd(VPARALLEL_DIR));
   const SliceSpec slice_mu(MU_DIR,domainBox.smallEnd(MU_DIR));

   CP1::Box CP1_injected_box;
   injectBox(CP1_injected_box, a_src.box(), slice_vp);
   Box injected_box;
   injectBox(injected_box, CP1_injected_box, slice_mu);

   a_dst.define(injected_box, a_src.nComp());

   CP1::BaseFab<Real> CP1_injected_fab;
   injectBaseFab(CP1_injected_fab, a_src, slice_vp);
   BaseFab<Real> injected_fab;
   injectBaseFab(injected_fab, CP1_injected_fab, slice_mu);
   a_dst.copy(injected_fab);
}



void
PhaseBlockCoordSys::injectConfigurationToPhase( const CFG::FArrayBox& a_src,
                                                const VEL::IntVect&   a_vpmu,
                                                FArrayBox&            a_dst ) const
{
   const CP1::SliceSpec slice_vp(VPARALLEL_DIR,a_vpmu[VPARALLEL_DIR-CFG_DIM]);
   const SliceSpec slice_mu(MU_DIR,a_vpmu[MU_DIR-CFG_DIM]);

   CP1::Box CP1_injected_box;
   injectBox(CP1_injected_box, a_src.box(), slice_vp);
   Box injected_box;
   injectBox(injected_box, CP1_injected_box, slice_mu);

   a_dst.define(injected_box, a_src.nComp());

   CP1::BaseFab<Real> CP1_injected_fab;
   injectBaseFab(CP1_injected_fab, a_src, slice_vp);
   BaseFab<Real> injected_fab;
   injectBaseFab(injected_fab, CP1_injected_fab, slice_mu);
   a_dst.copy(injected_fab);
}



void
PhaseBlockCoordSys::injectConfigurationToPhase( const CFG::FluxBox& a_src,
                                                FluxBox&            a_dst ) const
{
   const Box& domainBox = m_domain.domainBox();
   const CP1::SliceSpec slice_vp(VPARALLEL_DIR,domainBox.smallEnd(VPARALLEL_DIR));
   const SliceSpec slice_mu(MU_DIR,domainBox.smallEnd(MU_DIR));

   CP1::Box CP1_injected_box;
   injectBox(CP1_injected_box, a_src.box(), slice_vp);
   Box injected_box;
   injectBox(injected_box, CP1_injected_box, slice_mu);

   a_dst.define(injected_box, a_src.nComp());

   for (int dir=0; dir<CFG_DIM; ++dir) {
      CP1::BaseFab<Real> CP1_injected_fab;
      injectBaseFab(CP1_injected_fab, a_src[dir], slice_vp);
      CP1_injected_fab.shiftHalf(dir,-1);
      BaseFab<Real> injected_fab;
      injectBaseFab(injected_fab, CP1_injected_fab, slice_mu);
      injected_fab.shiftHalf(dir,-1);  // BaseFab injection seems to result in cell centering; fix that
      a_dst[dir].copy(injected_fab);
   }
}



void
PhaseBlockCoordSys::injectConfigurationToPhase( const CFG::FluxBox& a_src,
                                                const VEL::IntVect& a_vpmu,
                                                FluxBox&            a_dst ) const
{
   const CP1::SliceSpec slice_vp(VPARALLEL_DIR,a_vpmu[VPARALLEL_DIR-CFG_DIM]);
   const SliceSpec slice_mu(MU_DIR,a_vpmu[MU_DIR-CFG_DIM]);

   CP1::Box CP1_injected_box;
   injectBox(CP1_injected_box, a_src.box(), slice_vp);
   Box injected_box;
   injectBox(injected_box, CP1_injected_box, slice_mu);

   a_dst.define(injected_box, a_src.nComp());

   for (int dir=0; dir<CFG_DIM; ++dir) {
      CP1::BaseFab<Real> CP1_injected_fab;
      injectBaseFab(CP1_injected_fab, a_src[dir], slice_vp);
      CP1_injected_fab.shiftHalf(dir,-1);
      BaseFab<Real> injected_fab;
      injectBaseFab(injected_fab, CP1_injected_fab, slice_mu);
      injected_fab.shiftHalf(dir,-1);  // BaseFab injection seems to result in cell centering; fix that
      a_dst[dir].copy(injected_fab);
   }
}



void
PhaseBlockCoordSys::projectPhaseToConfiguration( const Box& a_src,
                                                 CFG::Box&  a_dst ) const
{
   const Box& domainBox = m_domain.domainBox();
   const CP1::SliceSpec slice_vp(VPARALLEL_DIR,domainBox.smallEnd(VPARALLEL_DIR));
   const SliceSpec slice_mu(MU_DIR,domainBox.smallEnd(MU_DIR));

   CP1::Box CP1_sliced_box;
   sliceBox(CP1_sliced_box, a_src, slice_mu);
   sliceBox(a_dst, CP1_sliced_box, slice_vp);
}



Box
PhaseBlockCoordSys::configurationFlatten( const Box& a_disjointBox,
                                          const Box& a_src_box ) const
{
   Box dst_box(a_src_box);

   for (int dir=CFG_DIM; dir<SpaceDim; ++dir) {
      int small = a_disjointBox.smallEnd(dir);
      dst_box.setSmall(dir, small);
      dst_box.setBig(dir, small);
   }

   return dst_box;
}



Box
PhaseBlockCoordSys::velocityFlatten( const Box& a_disjointBox,
                                     const Box& a_src_box ) const
{
   Box dst_box(a_src_box);

   for (int dir=0; dir<CFG_DIM; ++dir) {
      int small = a_disjointBox.smallEnd(dir);
      dst_box.setSmall(dir, small);
      dst_box.setBig(dir, small);
   }

   return dst_box;
}



#include "NamespaceFooter.H"




