#include "PhaseGeom.H"
#include "PhaseGeomF_F.H"
#include "Directions.H"
#include "FourthOrderUtil.H"
#include "newMappedGridIO.H"
#include "CONSTANTS.H"
#include "inspect.H"

#include "EdgeToCell.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "MagCoordSys.H"
#include "MagBlockCoordSys.H"
#include "MillerBlockCoordSys.H"
#include "newMappedGridIO.H"
#include "inspect.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM VEL_DIM
#include "newMappedGridIO.H"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM CP1_DIM
#include "AMRIO.H"
#include "RealVect.H"
#include "IntVectSet.H"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#include "Injection.H.transdim"
#include "Slicing.H.transdim"
#include "ReductionOps.H.multidim"
#include "SpreadingCopier.H.multidim"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif
#define CH_SPACEDIM PDIM

#include "SPACE.H"

#undef FIX_RADIAL_BOUNDARY_FLUX

#include "NamespaceHeader.H"

using namespace CH_MultiDim;


PhaseGeom::PhaseGeom( ParmParse&                      a_parm_parse,
                      const PhaseCoordSys*            a_coord_sys,
                      const PhaseGrid&                a_grids,
                      CFG::MagGeom&                   a_mag_geom,
                      const VEL::VelCoordSys&         a_vel_coords,
                      int                             a_ghosts,
                      double                          a_larmor_number )
   : MultiBlockLevelGeom(a_coord_sys, a_grids.disjointBoxLayout(), a_ghosts, Interval(CFG_DIM,SpaceDim-1), Vector<int>(VEL_DIM,0)),
     m_ghostVect(4*IntVect::Unit),
     m_mag_geom(a_mag_geom),
     m_vel_coords(a_vel_coords),
     m_phase_coords(*a_coord_sys),
     m_phase_grid(a_grids),
     m_domain(a_grids.domain()),
     m_speciesDefined(false),
     m_larmor_number(a_larmor_number)
{
   ParmParse psm( "phase_space_mapping" );
   if (psm.contains("velocity_type")) {
      psm.get("velocity_type", m_velocity_type);
   }
   else {
      m_velocity_type = "gyrokinetic";
   }

   if (m_velocity_type == "free_stream") {
     if (psm.contains("free_stream_components")) {
       Vector<Real> read_coords(SpaceDim);
       psm.getarr( "free_stream_components", read_coords, 0, SpaceDim );
       m_freestream_components = RealVect(read_coords);
     }
     else {
       for (int dir=0; dir<CFG_DIM; ++dir) {
         m_freestream_components[dir] = 1.;
       }
       for (int dir=CFG_DIM; dir<SpaceDim; ++dir) {
         m_freestream_components[dir] = 0.;
       }
     }
   }

   if (psm.contains("no_drifts")) {
      psm.get("no_drifts", m_no_drifts);
   }
   else {
      m_no_drifts = false;
   }

   if (psm.contains("no_parallel_streaming")) {
      psm.get("no_parallel_streaming", m_no_parallel_streaming);
   }
   else {
      m_no_parallel_streaming = false;
   }

   if (psm.contains("second_order")) {
      psm.get("second_order", m_second_order);
   }
   else {
      m_second_order = false;
   }

   define();
}



PhaseGeom::PhaseGeom( const PhaseGeom& a_phase_geom,
                      double           a_mass,
                      double           a_charge )
   : MultiBlockLevelGeom(a_phase_geom.coordSysPtr(),a_phase_geom.gridsFull(),a_phase_geom.ghosts(),Interval(CFG_DIM,SpaceDim-1),Vector<int>(VEL_DIM,0)),
     m_ghostVect(a_phase_geom.m_ghostVect),
     m_mag_geom(a_phase_geom.m_mag_geom),
     m_vel_coords(a_phase_geom.m_vel_coords),
     m_phase_coords(a_phase_geom.phaseCoordSys()),
     m_phase_grid(a_phase_geom.m_phase_grid),
     m_domain(a_phase_geom.m_domain),
     m_configuration_metrics(a_phase_geom.m_configuration_metrics),
     m_configuration_tangrad_metrics(a_phase_geom.m_configuration_tangrad_metrics),
     m_configuration_volumes(a_phase_geom.m_configuration_volumes),
     m_configuration_J(a_phase_geom.m_configuration_J),
     m_configuration_face_areas(a_phase_geom.m_configuration_face_areas),
     m_velocity_metrics(a_phase_geom.m_velocity_metrics),
     m_velocity_tangrad_metrics(a_phase_geom.m_velocity_tangrad_metrics),
     m_velocity_volumes(a_phase_geom.m_velocity_volumes),
     m_velocity_J(a_phase_geom.m_velocity_J),
     m_velocity_face_areas(a_phase_geom.m_velocity_face_areas),
     m_BCell(a_phase_geom.m_BCell),
     m_gradBCell(a_phase_geom.m_gradBCell),
     m_curlbCell(a_phase_geom.m_curlbCell),
     m_BMagCell(a_phase_geom.m_BMagCell),
     m_bdotcurlbCell(a_phase_geom.m_bdotcurlbCell),
     m_BFace(a_phase_geom.m_BFace),
     m_gradBFace(a_phase_geom.m_gradBFace),
     m_curlbFace(a_phase_geom.m_curlbFace),
     m_BMagFace(a_phase_geom.m_BMagFace),
     m_bdotcurlbFace(a_phase_geom.m_bdotcurlbFace),
     m_mblexPtr(a_phase_geom.m_mblexPtr),
     m_exchange_transverse_block_register(a_phase_geom.m_exchange_transverse_block_register),
     m_velocity_type(a_phase_geom.m_velocity_type),
     m_no_drifts(a_phase_geom.m_no_drifts),
     m_no_parallel_streaming(a_phase_geom.m_no_parallel_streaming),
     m_speciesDefined(false),
     m_freestream_components(a_phase_geom.m_freestream_components),
     m_second_order(a_phase_geom.m_second_order),
     m_larmor_number(a_phase_geom.m_larmor_number)
{
   defineSpeciesState(a_mass, a_charge);
}



void
PhaseGeom::define()
{
   /*
    *  Inject the configuration and velocity space metrics, as well as the magnetic field
    *  quantities at cell and face centers, so that we have them on the phase space domain
    *  decomposition.  Since these quantities never change but are used frequently, we
    *  precompute and store them here.
    */

   CFG::IntVect cfg_ghostVect(config_restrict(m_ghostVect));
   VEL::IntVect vel_ghostVect(vel_restrict(m_ghostVect));

   /*
    *  Get the configuration space data and inject into phase space
    */

   const CFG::DisjointBoxLayout& cfg_grids = m_mag_geom.gridsFull();

   int cfg_NumN = CFG_DIM * CFG_DIM;
   CFG::LevelData<CFG::FluxBox> cfg_N(cfg_grids, cfg_NumN, cfg_ghostVect);
   CFG::LevelData<CFG::FluxBox> cfg_tanGradN(cfg_grids, cfg_NumN, cfg_ghostVect);
   m_mag_geom.getMetricTerms(cfg_N, cfg_tanGradN);

   m_configuration_metrics = new LevelData<FluxBox>;
   m_configuration_tangrad_metrics = new LevelData<FluxBox>;
   injectConfigurationToPhase(cfg_N, *m_configuration_metrics);
   injectConfigurationToPhase(cfg_tanGradN, *m_configuration_tangrad_metrics);

   CFG::LevelData<CFG::FArrayBox> cfg_volume(cfg_grids, 1, cfg_ghostVect);
   m_mag_geom.getCellVolumes(cfg_volume);
   m_configuration_volumes = new LevelData<FArrayBox>;
   injectConfigurationToPhase(cfg_volume, *m_configuration_volumes);

   // Include an extra ghost cell layer since functions using this data will
   // sometimes be using it with the fourth order product/quotient formula
   CFG::LevelData<CFG::FArrayBox> cfg_J(cfg_grids, 1,
                                        cfg_ghostVect + CFG::IntVect::Unit);
   m_mag_geom.getJ(cfg_J);
   m_configuration_J = new LevelData<FArrayBox>;
   injectConfigurationToPhase(cfg_J, *m_configuration_J);

   CFG::LevelData<CFG::FluxBox> cfg_face_areas(cfg_grids, 1, CFG::IntVect::Zero);
   m_mag_geom.getFaceAreas(cfg_face_areas);
   m_configuration_face_areas = new LevelData<FluxBox>;
   injectConfigurationToPhase(cfg_face_areas, *m_configuration_face_areas);

   /*
    *  Get the velocity space data and inject into phase space
    */

   const VEL::DisjointBoxLayout& vel_grids = m_vel_coords.grids();
   int vel_NumN = VEL_DIM * VEL_DIM;

   VEL::LevelData<VEL::FluxBox> vel_N(vel_grids, vel_NumN, vel_ghostVect);
   VEL::LevelData<VEL::FluxBox> vel_tanGradN(vel_grids, vel_NumN, vel_ghostVect);
   m_vel_coords.getMetricTerms(vel_N, vel_tanGradN);

   m_velocity_metrics = new LevelData<FluxBox>;
   m_velocity_tangrad_metrics = new LevelData<FluxBox>;
   injectVelocityToPhase(vel_N, *m_velocity_metrics);
   injectVelocityToPhase(vel_tanGradN, *m_velocity_tangrad_metrics);

   VEL::LevelData<VEL::FArrayBox> vel_volume(vel_grids, 1, vel_ghostVect);
   m_vel_coords.getCellVolumes(vel_volume);

   m_velocity_volumes = new LevelData<FArrayBox>;
   injectVelocityToPhase(vel_volume, *m_velocity_volumes);

   // Include an extra ghost cell layer since functions using this data will
   // sometimes be using it with the fourth order product/quotient formula
   VEL::LevelData<VEL::FArrayBox> vel_J(vel_grids, 1,
                                        vel_ghostVect + VEL::IntVect::Unit);
   m_vel_coords.getJ(vel_J);
   m_velocity_J = new LevelData<FArrayBox>;
   injectVelocityToPhase(vel_J, *m_velocity_J);

   VEL::LevelData<VEL::FluxBox> vel_face_areas(vel_grids, 1, VEL::IntVect::Zero);
   m_vel_coords.getFaceAreas(vel_face_areas);
   m_velocity_face_areas = new LevelData<FluxBox>;
   injectVelocityToPhase(vel_face_areas, *m_velocity_face_areas);

   /*
    *  Get the cell- and face-centered field data and inject into phase space
    */

   const CFG::LevelData<CFG::FArrayBox>& cfg_BFieldCell = m_mag_geom.getCCBField();
   CH_assert(cfg_BFieldCell.ghostVect() >= cfg_ghostVect);

   const CFG::LevelData<CFG::FArrayBox>& cfg_BFieldMagCell = m_mag_geom.getCCBFieldMag();
   CH_assert(cfg_BFieldMagCell.ghostVect() >= cfg_ghostVect);

   const CFG::LevelData<CFG::FArrayBox>& cfg_BFieldDirCell = m_mag_geom.getCCBFieldDir();
   CH_assert(cfg_BFieldDirCell.ghostVect() >= cfg_ghostVect);

   const CFG::LevelData<CFG::FArrayBox>& cfg_gradBFieldMagCell = m_mag_geom.getCCGradBFieldMag();
   CH_assert(cfg_gradBFieldMagCell.ghostVect() >= cfg_ghostVect);

   const CFG::LevelData<CFG::FArrayBox>& cfg_curlBFieldDirCell = m_mag_geom.getCCCurlBFieldDir();
   CH_assert(cfg_curlBFieldDirCell.ghostVect() >= cfg_ghostVect);

   const CFG::LevelData<CFG::FArrayBox>& cfg_BFieldDirdotcurlBFieldDirCell = m_mag_geom.getCCBFieldDirdotCurlBFieldDir();
   CH_assert(cfg_BFieldDirdotcurlBFieldDirCell.ghostVect() >= cfg_ghostVect);

   m_BCell = new LevelData<FArrayBox>;
   m_gradBCell = new LevelData<FArrayBox>;
   m_curlbCell = new LevelData<FArrayBox>;
   m_BMagCell = new LevelData<FArrayBox>;
   m_bdotcurlbCell = new LevelData<FArrayBox>;
   injectConfigurationToPhase(cfg_BFieldCell, *m_BCell);
   injectConfigurationToPhase(cfg_gradBFieldMagCell, *m_gradBCell);
   injectConfigurationToPhase(cfg_curlBFieldDirCell, *m_curlbCell);
   injectConfigurationToPhase(cfg_BFieldMagCell, *m_BMagCell);
   injectConfigurationToPhase(cfg_BFieldDirdotcurlBFieldDirCell, *m_bdotcurlbCell);

   const CFG::LevelData<CFG::FluxBox>& cfg_BFieldFace = m_mag_geom.getFCBField();
   CH_assert(cfg_BFieldFace.ghostVect() >= cfg_ghostVect);

   const CFG::LevelData<CFG::FluxBox>& cfg_BFieldMagFace = m_mag_geom.getFCBFieldMag();
   CH_assert(cfg_BFieldMagFace.ghostVect() >= cfg_ghostVect);

   const CFG::LevelData<CFG::FluxBox>& cfg_BFieldDirFace = m_mag_geom.getFCBFieldDir();
   CH_assert(cfg_BFieldDirFace.ghostVect() >= cfg_ghostVect);

   const CFG::LevelData<CFG::FluxBox>& cfg_gradBFieldMagFace = m_mag_geom.getFCGradBFieldMag();
   CH_assert(cfg_gradBFieldMagFace.ghostVect() >= cfg_ghostVect);

   const CFG::LevelData<CFG::FluxBox>& cfg_curlBFieldDirFace = m_mag_geom.getFCCurlBFieldDir();
   CH_assert(cfg_curlBFieldDirFace.ghostVect() >= cfg_ghostVect);

   const CFG::LevelData<CFG::FluxBox>& cfg_BFieldDirdotcurlBFieldDirFace = m_mag_geom.getFCBFieldDirdotCurlBFieldDir();
   CH_assert(cfg_BFieldDirdotcurlBFieldDirFace.ghostVect() >= cfg_ghostVect);
   
   m_BFace = new LevelData<FluxBox>;
   m_gradBFace = new LevelData<FluxBox>;
   m_curlbFace = new LevelData<FluxBox>;
   m_BMagFace = new LevelData<FluxBox>;
   m_bdotcurlbFace = new LevelData<FluxBox>;
   injectConfigurationToPhase(cfg_BFieldFace, cfg_BFieldCell, *m_BFace);
   injectConfigurationToPhase(cfg_gradBFieldMagFace, cfg_gradBFieldMagCell, *m_gradBFace);
   injectConfigurationToPhase(cfg_curlBFieldDirFace, cfg_curlBFieldDirCell, *m_curlbFace);
   injectConfigurationToPhase(cfg_BFieldMagFace, cfg_BFieldMagCell, *m_BMagFace);
   injectConfigurationToPhase(cfg_BFieldDirdotcurlBFieldDirFace, cfg_BFieldDirdotcurlBFieldDirCell, *m_bdotcurlbFace);

   if ( m_coordSysPtr->numBlocks() > 1 ) {

      m_mblexPtr = new MultiBlockLevelExchangeAverage();

      int nghost = 0;
      for (int dim=0; dim<SpaceDim; ++dim) {
         if (m_ghostVect[dim] > nghost) nghost = m_ghostVect[dim];
      }
      int numGhost(nghost);
      int spaceOrder(4);
      m_mblexPtr->define( this, numGhost, spaceOrder );

      RefCountedPtr<MultiBlockCoordSys> coordSysRCP((MultiBlockCoordSys*)m_coordSysPtr);
      coordSysRCP.neverDelete();  // Works around some problem with RefCountedPtr

      m_exchange_transverse_block_register = new BlockRegister(coordSysRCP, m_gridsFull, 1);
   }
   else {
      m_mblexPtr = NULL;
      m_exchange_transverse_block_register = NULL;
   }
}


void
PhaseGeom::defineSpeciesState( double a_mass,
                               double a_charge_state )
{
   m_mass = a_mass;
   m_charge_state = a_charge_state;

   // Since BStar never changes, we precompute and store it here.
   m_BStar.define(m_gridsFull, 3, m_ghostVect);
   m_BStarParallel.define(m_gridsFull, 1, m_ghostVect);

   computeBStar(m_no_drifts,
                m_larmor_number,
                m_mass,
                m_charge_state,
                *m_BCell,
                *m_BMagCell,
                *m_curlbCell,
                *m_bdotcurlbCell,
                m_BStar,
                m_BStarParallel);

   m_speciesDefined = true;
}



PhaseGeom::~PhaseGeom()
{
   if (m_speciesDefined) {
      delete m_bdotcurlbFace;
      delete m_BMagFace;
      delete m_curlbFace;
      delete m_gradBFace;
      delete m_BFace;
      delete m_bdotcurlbCell;
      delete m_BMagCell;
      delete m_curlbCell;
      delete m_gradBCell;
      delete m_BCell;
      delete m_velocity_J;
      delete m_velocity_volumes;
      delete m_velocity_tangrad_metrics;
      delete m_velocity_metrics;
      delete m_configuration_J;
      delete m_configuration_volumes;
      delete m_configuration_tangrad_metrics;
      delete m_configuration_metrics;
   }

   if (m_exchange_transverse_block_register) delete m_exchange_transverse_block_register;
   if (m_mblexPtr) delete m_mblexPtr;
}



void
PhaseGeom::updateVelocities( const LevelData<FluxBox>& a_Efield,
                             LevelData<FluxBox>&       a_velocity,
                             const bool                a_apply_axisymmetric_correction ) const
{
   // This function expects a_Efield to contain the face-centered field including
   // one layer of transverse cell faces on all block boundaries.
   CH_assert(a_velocity.ghostVect() == IntVect::Unit);
   CH_assert(config_restrict(a_Efield.ghostVect()) == CFG::IntVect::Unit);
   CH_assert(vel_restrict(a_Efield.ghostVect()) == VEL::IntVect::Zero);

   if (m_velocity_type == "gyrokinetic" || m_velocity_type == "ExB") {
      computeGKVelocities(a_Efield, a_velocity);
   }
   else {
      computeTestVelocities(a_velocity);
   }

   //Corrects for the extra R factor added to the metric coefficients (e.g., cellVoll)
   //for the purposes of 4h order. The model ("analytic") geometry does
   //not have that extra factor.
   if ( a_apply_axisymmetric_correction) {
      applyAxisymmetricCorrection(a_velocity);
   }

   if ( !m_second_order ) {

      // Convert face-centered values to face averages on valid cell faces.
      // The values on ghost cell faces remain second-order.
      fourthOrderAverage(a_velocity);
   }

   a_velocity.exchange();
}



void
PhaseGeom::computeGKVelocities( const LevelData<FluxBox>& a_Efield,
                                LevelData<FluxBox>&       a_velocity ) const
{
   const CFG::IntVect& velocity_cfg_ghosts = config_restrict(a_velocity.ghostVect());
   CH_assert(config_restrict(a_Efield.ghostVect()) >= velocity_cfg_ghosts);
   CH_assert(config_restrict(m_BFace->ghostVect()) >= velocity_cfg_ghosts);
   CH_assert(config_restrict(m_gradBFace->ghostVect()) >= velocity_cfg_ghosts);
   CH_assert(config_restrict(m_curlbFace->ghostVect()) >= velocity_cfg_ghosts);

   // Set options to turn certain terms on or off
   const int include_drifts = m_no_drifts? 0: 1;
   const int include_gradb = (m_velocity_type == "ExB")? 0: 1;
   const int include_par_streaming = (m_no_parallel_streaming || m_velocity_type == "ExB")? 0: 1;
   //use field alignment presently turns off parallel streaming and EXB drift on radial faces (assumes no Epol!!!)
   const int use_field_alignment = (m_mag_geom.modelGeom())? 1: 0;

   const DisjointBoxLayout& grids = a_velocity.disjointBoxLayout();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      RealVect dx = block_coord_sys.dx();

      const FluxBox& this_Efield = a_Efield[dit];
      const FluxBox& this_B = (*m_BFace)[dit];
      const FluxBox& this_gradB = (*m_gradBFace)[dit];
      const FluxBox& this_curlb = (*m_curlbFace)[dit];
      FluxBox& this_velocity = a_velocity[dit];

      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_velocity_dir = this_velocity[dir];

         FORT_COMPUTE_GK_VELOCITY(
                                  CHF_CONST_INT(dir),
                                  CHF_BOX(this_velocity_dir.box()),
                                  CHF_CONST_REALVECT(dx),
                                  CHF_CONST_REAL(m_charge_state),
                                  CHF_CONST_REAL(m_mass),
                                  CHF_CONST_REAL(m_larmor_number),
                                  CHF_CONST_INT(include_drifts),
                                  CHF_CONST_INT(include_par_streaming),
                                  CHF_CONST_INT(include_gradb),
				  CHF_CONST_INT(use_field_alignment),
                                  CHF_CONST_FRA(this_Efield[dir]),
                                  CHF_CONST_FRA(this_B[dir]),
                                  CHF_CONST_FRA(this_gradB[dir]),
                                  CHF_CONST_FRA(this_curlb[dir]),
                                  CHF_FRA(this_velocity_dir)
                                  );
      }
   }
}



void
PhaseGeom::computeTestVelocities( LevelData<FluxBox>& a_velocity ) const
{
   const DisjointBoxLayout& grids = a_velocity.disjointBoxLayout();

   if (m_velocity_type == "bfield") {
      for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         const FluxBox& this_B = (*m_BFace)[dit];
         FluxBox& this_velocity = a_velocity[dit];
         for (int dir=0; dir<SpaceDim; dir++) {
            FArrayBox& this_velocity_dir = this_velocity[dir];

            FORT_COMPUTE_BFIELD_VELOCITY(
                                         CHF_BOX(this_velocity_dir.box()),
                                         CHF_CONST_FRA(this_B[dir]),
                                         CHF_FRA(this_velocity_dir)
                                         );
         }
      }
   }
   else if (m_velocity_type == "annular_poloidal_velocity") {

      DataIterator dit = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
         const CFG::MillerBlockCoordSys* mag_block_coord_sys
           = dynamic_cast<const CFG::MillerBlockCoordSys*>(&getMagBlockCoordSys(grids[dit]));
         const PhaseBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);

         if (mag_block_coord_sys != NULL) {

           double rmin = mag_block_coord_sys->getInnerFluxLabel();
           Real Rmajor0 = mag_block_coord_sys->centralMajorRadius();
           Real rbar = mag_block_coord_sys->getAvMinorRad();
           int const_minorrad = mag_block_coord_sys->getConstMinorrad();
           RealVect dx = block_coord_sys.dx();

           for (int dir=0; dir<SpaceDim; dir++) {
             FArrayBox& thisVel = a_velocity[dit][dir];
             FORT_ANNULUS_POLVEL_TEST(CHF_CONST_INT(dir),
                                      CHF_BOX(thisVel.box()),
                                      CHF_CONST_REALVECT(dx),
                                      CHF_CONST_REAL(rmin),
                                      CHF_CONST_REAL(rbar),
                                      CHF_CONST_REAL(Rmajor0),
                                      CHF_FRA(thisVel),
                                      CHF_INT(const_minorrad));
           }
         }
         else {
           MayDay::Error("PhaseGeom::computeTestVelocities(): annular_poloidal_velocity is only implemented for Miller geometry");
         }
      }
   }
   else if (m_velocity_type == "annular_radial_velocity") {

      DataIterator dit = a_velocity.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
         const CFG::MillerBlockCoordSys* mag_block_coord_sys
           = dynamic_cast<const CFG::MillerBlockCoordSys*>(&getMagBlockCoordSys(grids[dit]));
         const PhaseBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);

         if (mag_block_coord_sys != NULL) {

           double rmin = mag_block_coord_sys->getInnerFluxLabel();
           Real Rmajor0 = mag_block_coord_sys->centralMajorRadius();
           Real rbar = mag_block_coord_sys->getAvMinorRad();
           int const_minorrad = mag_block_coord_sys->getConstMinorrad();
           RealVect dx = block_coord_sys.dx();

           for (int dir=0; dir<SpaceDim; dir++) {
             FArrayBox& thisVel = a_velocity[dit][dir];
             FORT_ANNULUS_RADVEL_TEST(CHF_CONST_INT(dir),
                                      CHF_BOX(thisVel.box()),
                                      CHF_CONST_REALVECT(dx),
                                      CHF_CONST_REAL(rmin),
                                      CHF_CONST_REAL(rbar),
                                      CHF_CONST_REAL(Rmajor0),
                                      CHF_FRA(thisVel),
                                      CHF_INT(const_minorrad));
           }
         }
         else {
           MayDay::Error("PhaseGeom::computeTestVelocities(): annular_radial_velocity is only implemented for Miller geometry");
         }
      }
   }
   else if (m_velocity_type == "annular_radpol_velocity") {

      DataIterator dit = a_velocity.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
         const CFG::MillerBlockCoordSys* mag_block_coord_sys
           = dynamic_cast<const CFG::MillerBlockCoordSys*>(&getMagBlockCoordSys(grids[dit]));
         const PhaseBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);

         if(mag_block_coord_sys != NULL) {

           double rmin = mag_block_coord_sys->getInnerFluxLabel();
           Real Rmajor0 = mag_block_coord_sys->centralMajorRadius();
           Real rbar = mag_block_coord_sys->getAvMinorRad();
           int const_minorrad = mag_block_coord_sys->getConstMinorrad();
           RealVect dx = block_coord_sys.dx();

           for (int dir=0; dir<SpaceDim; dir++) {
             FArrayBox& thisVel = a_velocity[dit][dir];
             FORT_ANNULUS_RADPOLVEL_TEST(CHF_CONST_INT(dir),
                                         CHF_BOX(thisVel.box()),
                                         CHF_CONST_REALVECT(dx),
                                         CHF_CONST_REAL(rmin),
                                         CHF_CONST_REAL(rbar),
                                         CHF_CONST_REAL(Rmajor0),
                                         CHF_FRA(thisVel),
                                         CHF_INT(const_minorrad));
           }
         }
         else {
           MayDay::Error("PhaseGeom::computeTestVelocities(): annular_radpol_velocity is only implemented for Miller geometry");
         }
      }
   }
   else if (m_velocity_type == "major_radial_velocity") {

      DataIterator dit = a_velocity.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
//       const CFG::MagBlockCoordSys* mag_block_coord_sys
//         = dynamic_cast<const CFG::MillerBlockCoordSys*>(&getMagBlockCoordSys(grids[dit]));
         const PhaseBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);

//       if(mag_block_coord_sys != NULL) {

//         double rmin = mag_block_coord_sys->getInnerFluxLabel();
//         Real Rmajor0 = mag_block_coord_sys->centralMajorRadius();
//         Real rbar = mag_block_coord_sys->getAvMinorRad();
//         int const_minorrad = mag_block_coord_sys->getConstMinorrad();
           RealVect dx = block_coord_sys.dx();

           for (int dir=0; dir<SpaceDim; dir++) {
             FArrayBox& thisVel = a_velocity[dit][dir];
             FORT_MAJOR_RADIAL_VEL_TEST(CHF_CONST_INT(dir),
                                        CHF_BOX(thisVel.box()),
//                             CHF_CONST_REALVECT(dx),
//                             CHF_CONST_REAL(rmin),
//                             CHF_CONST_REAL(rbar),
//                             CHF_CONST_REAL(Rmajor0),
                                        CHF_FRA(thisVel) );
//                             CHF_INT(const_minorrad));
           }
//       }
//       else {
//         MayDay::Error("PhaseGeom::computeTestVelocities(): major_radial_velocity is only implemented for Miller geometry");
//       }
      }
   }
   else if (m_velocity_type == "free_stream") {

     DataIterator dit = a_velocity.dataIterator();
     for (dit.begin(); dit.ok(); ++dit) {
       const CFG::MagBlockCoordSys& mag_block_coord_sys = getMagBlockCoordSys(grids[dit]);

       FluxBox& this_vel = a_velocity[dit];
       Box box(this_vel.box());
       CFG::Box box_config;
       projectPhaseToConfiguration(box, box_config);
       CFG::FluxBox R_cfg(box_config, 1);

       mag_block_coord_sys.getPointwiseMajorRadius(R_cfg);
       FluxBox R;
       injectConfigurationToPhase(R_cfg, R);

       int axisymmetric = mag_block_coord_sys.isAxisymmetric()? 1: 0;

       for (int dir=0; dir<SpaceDim; dir++) {
         FORT_FREE_STREAM_VEL_TEST(CHF_BOX(this_vel[dir].box()),
                                   CHF_CONST_REALVECT(m_freestream_components),
                                   CHF_CONST_FRA1(R[dir],0),
                                   CHF_CONST_INT(axisymmetric),
                                   CHF_FRA(a_velocity[dit][dir]) );
       }
     }
   }
   else {
     MayDay::Error("PhaseBlockCoordSys::computeTestVelocities: Unknown velocity type");
   }
}



void
PhaseGeom::updateMappedVelocities( const LevelData<FluxBox>& a_Efield,
                                   LevelData<FluxBox>&       a_velocity ) const
{
   // N.B.: this is only second-order

   updateVelocities(a_Efield, a_velocity, false);

   multNTransposePointwise( a_velocity );
}



void
PhaseGeom::applyAxisymmetricCorrection( LevelData<FluxBox>& a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();

   for (dit.begin(); dit.ok(); ++dit) {
     const CFG::MagBlockCoordSys& mag_block_coord_sys = getMagBlockCoordSys(grids[dit]);

     if ( mag_block_coord_sys.isAxisymmetric() && !(m_mag_geom.modelGeom()) ) {

       FluxBox& this_data = a_data[dit];
       Box box(this_data.box());
       CFG::Box box_config;
       projectPhaseToConfiguration(box, box_config);
       CFG::FluxBox R_cfg(box_config, 1);

       mag_block_coord_sys.getPointwiseMajorRadius(R_cfg);
       FluxBox R;
       injectConfigurationToPhase(R_cfg, R);

       for (int dir=0; dir<CFG_DIM; ++dir) {
         FORT_MULT_RADIAL_COMPONENT(CHF_BOX(surroundingNodes(box,dir)),
                                    CHF_CONST_FRA1(R[dir],0),
                                    CHF_FRA1(this_data[dir],RADIAL_DIR));
       }
     }
   }
}



void
PhaseGeom::multNTransposePointwise( LevelData<FluxBox>& a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();

   for (dit.begin(); dit.ok(); ++dit) {
     const CFG::MagBlockCoordSys& mag_block_coord_sys = getMagBlockCoordSys(grids[dit]);

     FluxBox& this_data = a_data[dit];
     Box box(this_data.box());
     CFG::Box box_config;
     projectPhaseToConfiguration(box, box_config);
     CFG::FluxBox N_cfg_face(box_config, CFG_DIM*CFG_DIM);

     mag_block_coord_sys.getPointwiseN(N_cfg_face);
     FluxBox N_face;
     injectConfigurationToPhase(N_cfg_face, N_face);

     for (int dir=0; dir<CFG_DIM; ++dir) {
       FORT_MULT_CFG_NT(CHF_BOX(surroundingNodes(box,dir)),
                        CHF_CONST_FRA(N_face[dir]),
                        CHF_FRA(this_data[dir]));
     }

     CFG::FArrayBox TwoPirRmaj_cfg_cell(box_config, 1);
     mag_block_coord_sys.getTwoPirRmaj(TwoPirRmaj_cfg_cell);
     FArrayBox TwoPirRmaj_cell;
     injectConfigurationToPhase(TwoPirRmaj_cfg_cell, TwoPirRmaj_cell);

     for (int dir=CFG_DIM; dir<SpaceDim; ++dir) {
       FORT_MULT_VEL_NT(CHF_BOX(surroundingNodes(box,dir)),
                        CHF_CONST_FRA1(TwoPirRmaj_cell,0),
                        CHF_FRA(this_data[dir]));
     }
   }
}



void
PhaseGeom::computeBStar(bool                        a_no_drifts,
                        double                      a_larmor_number,
                        double                      a_mass,
                        double                      a_charge_state,
                        const LevelData<FArrayBox>& a_B,
                        const LevelData<FArrayBox>& a_BMag,
                        const LevelData<FArrayBox>& a_curlb,
                        const LevelData<FArrayBox>& a_bdotcurlb,
                        LevelData<FArrayBox>&       a_BStar,
                        LevelData<FArrayBox>&       a_BStarParallel) const
{
   double prefactor = a_no_drifts? 0.: a_larmor_number * a_mass / a_charge_state;

   const DisjointBoxLayout& grids = a_BStar.disjointBoxLayout();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      const RealVect& dx = block_coord_sys.dx();

      FORT_COMPUTE_BSTAR(CHF_BOX(a_BStar[dit].box()),
                         CHF_CONST_REALVECT(dx),
                         CHF_CONST_REAL(prefactor),
                         CHF_CONST_FRA(a_B[dit]),
                         CHF_CONST_FRA1(a_BMag[dit],0),
                         CHF_CONST_FRA(a_curlb[dit]),
                         CHF_CONST_FRA1(a_bdotcurlb[dit],0),
                         CHF_FRA(a_BStar[dit]),
                         CHF_FRA1(a_BStarParallel[dit],0));
   }
   a_BStar.exchange();
   a_BStarParallel.exchange();
}



void
PhaseGeom::computeMetricTermProductAverage( LevelData<FluxBox>&       a_Product,
                                            const LevelData<FluxBox>& a_F,
                                            bool                      a_fourthOrder ) const
{
   CH_assert( !a_fourthOrder || a_F.ghostVect() >= IntVect::Unit );

   // Determine whether only face normals are being computed
   bool face_normals_only;
   if (a_F.nComp() == SpaceDim) {
      face_normals_only = true;
   }
   else if (a_F.nComp() == m_configuration_metrics->nComp()
                         + m_velocity_metrics->nComp()) {
      face_normals_only = false;
   }
   else {
      MayDay::Error("PhaseBlockCoordSys::computeMetricTermProductAverage: Improper number of components in a_F");
   }

   const DisjointBoxLayout& grids = a_F.getBoxes();

   // temp storage
   IntVect grownGhost2(a_Product.ghostVect());
   grownGhost2 += 2*IntVect::Unit;

   int nGradFTerms_cfg = CFG_DIM*CFG_DIM*(CFG_DIM-1);
   int nGradFTerms_vel = VEL_DIM*VEL_DIM*(VEL_DIM-1);
   int nGradFTerms = Max(nGradFTerms_cfg,nGradFTerms_vel);

   LevelData<FluxBox> tanGradF(grids,
                               nGradFTerms,
                               grownGhost2);

   LevelData<FluxBox> dotTanGrads(a_Product.getBoxes(),
                                  a_Product.nComp(),
                                  grownGhost2);

   if (a_fourthOrder) {
      // compute tangential gradients of F
      computeTangentialGradSpecial(tanGradF, a_F);
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      const CFG::MagBlockCoordSys& mag_block_coord_sys = getMagBlockCoordSys(grids[dit]);

      const FluxBox& thisF = a_F[dit];
      FluxBox& thisFlux = a_Product[dit];

      // precomputed face-centered metrics
      const FluxBox& thisFaceMetric_cfg = (*m_configuration_metrics)[dit];
      const FluxBox& thisFaceMetric_vel = (*m_velocity_metrics)[dit];

      // precomputed face-centered metric transverse gradients
      const FluxBox& thisGradN_cfg = (*m_configuration_tangrad_metrics)[dit];
      const FluxBox& thisGradN_vel = (*m_velocity_tangrad_metrics)[dit];

      FluxBox& thisGradF = tanGradF[dit];
      FluxBox& thisDotGrad = dotTanGrads[dit];

      for (int dir=0; dir<SpaceDim; dir++) {
         bool configuration_dir = dir < CFG_DIM;
         int rel_dir = configuration_dir? dir: dir-CFG_DIM;  // dir relative to cfg or vel space

         FArrayBox& thisFluxDir = thisFlux[dir];
         const FArrayBox& thisFDir = thisF[dir];

         // get the appropriate metrics for this direction
         const FArrayBox * faceMetrics;
         if (configuration_dir) {
            faceMetrics = &thisFaceMetric_cfg[dir];
         }
         else {
            faceMetrics = &thisFaceMetric_vel[dir];
         }
         const FArrayBox& thisFaceMetricDir = *faceMetrics;

         // first part is easy -- dot product of F and metrics
         int metricStartComp = -1;
         int fStartComp = configuration_dir? 0: CFG_DIM;
         int ncomp = -1;

         // two possibilities here
         if (face_normals_only) {
            // one vector component of F
            // in this case, we only use the face-normal terms here if
            // there are only SpaceDim components to F.
            // note that since we're multiplying by N^T instead of N.
            // we're doing (0,dir) instead of (dir,0).
            // because of the way that the metric terms are indexed,
            // metricStartComp will multiply by F1,
            // metricStartComp+1 will multiply by F2, etc...
            if (configuration_dir) {
               metricStartComp = mag_block_coord_sys.getNcomponent(0,dir);
               ncomp = CFG_DIM;
            }
            else {
               metricStartComp = m_vel_coords.getNcomponent(0,rel_dir);
               ncomp = VEL_DIM;
            }
         }
         // otherwise, do all of the components
         else {
            metricStartComp = 0;
            if (configuration_dir) {
               ncomp = m_configuration_metrics->nComp();
            }
            else {
               ncomp = m_velocity_metrics->nComp();
            }
         }

         if (configuration_dir) {
            FORT_PSCSPOINTDOTPRODCFG(CHF_FRA1(thisFluxDir, 0),
                                     CHF_CONST_FRA(thisFDir),
                                     CHF_CONST_INT(fStartComp),
                                     CHF_CONST_FRA(thisFaceMetricDir),
                                     CHF_CONST_INT(metricStartComp),
                                     CHF_CONST_INT(ncomp),
                                     CHF_BOX(thisFluxDir.box()) );
         }
         else {
            FORT_PSCSPOINTDOTPRODVEL(CHF_FRA1(thisFluxDir, 0),
                                     CHF_CONST_FRA(thisFDir),
                                     CHF_CONST_INT(fStartComp),
                                     CHF_CONST_FRA(thisFaceMetricDir),
                                     CHF_CONST_INT(metricStartComp),
                                     CHF_CONST_INT(ncomp),
                                     CHF_BOX(thisFluxDir.box()) );
         }

         if (a_fourthOrder) {
            // now compute correction term
            // first do pointwise dot products of gradF and gradN

            // get the appropriate metric transverse gradients for this direction
            const FArrayBox* gradNdir;
            if (configuration_dir) {
               gradNdir = &thisGradN_cfg[dir];
            }
            else {
               gradNdir = &thisGradN_vel[dir];
            }
            const FArrayBox& thisGradNDir = *gradNdir;

            FArrayBox& thisGradFDir = thisGradF[dir];
            FArrayBox& thisDotGradDir = thisDotGrad[dir];
            int fGradStartComp = 0;
            // this little bit of wierdness is due to the fact that
            // the tanGrad direction is not allowed to be the same as the
            // face direction (i.e. dir can't equal nGradDir)
            int nGradDir = 0;
            if (rel_dir == 0) nGradDir = 1;
            if (configuration_dir) {
               int nGradStartComp = mag_block_coord_sys.tanGradComp(rel_dir, nGradDir, metricStartComp);
               int nGradComp = ncomp*(CFG_DIM-1);

               FORT_PSCSPOINTWISEREDUCEDDOTPRODCFG(CHF_FRA1(thisDotGradDir, 0),
                                                   CHF_CONST_FRA(thisGradNDir),
                                                   CHF_CONST_INT(nGradStartComp),
                                                   CHF_CONST_FRA(thisGradFDir),
                                                   CHF_CONST_INT(fGradStartComp),
                                                   CHF_CONST_INT(nGradComp),
                                                   CHF_BOX(thisDotGradDir.box()) );
            }
            else {
               int nGradStartComp = m_vel_coords.tanGradComp(rel_dir, nGradDir, metricStartComp);
               int nGradComp = ncomp*(VEL_DIM-1);

               FORT_PSCSPOINTWISEREDUCEDDOTPRODVEL(CHF_FRA1(thisDotGradDir, 0),
                                                   CHF_CONST_FRA(thisGradNDir),
                                                   CHF_CONST_INT(nGradStartComp),
                                                   CHF_CONST_FRA(thisGradFDir),
                                                   CHF_CONST_INT(fGradStartComp),
                                                   CHF_CONST_INT(nGradComp),
                                                   CHF_BOX(thisDotGradDir.box()) );
            }

            // since we computed undivided differences for the tanGrad's,
            //dx should also be one
            Real mult = 1.0/12.0;

            thisFluxDir.plus(thisDotGradDir, mult);
         } // end if we're doing the 4th-order corrections

         // Finish multiplication by face areas.  Like the base class from which it
         // is derived, this function is to return a face integral (rather than
         // a face average).  However, the configuration and velocity space metric
         // factors don't know about phase space, and therefore the area factors
         // included in them need to be augmented to phase space areas.  That
         // happens here.
         if (configuration_dir) {
            const FArrayBox& thisvol = (*m_velocity_volumes)[dit];
            FORT_MULT_VEL(CHF_BOX(thisFluxDir.box()),
                          CHF_FRA1(thisvol,0),
                          CHF_FRA(thisFluxDir));
         }
         else {
            const FArrayBox& thisvol = (*m_configuration_volumes)[dit];
            FORT_MULT_CFG(CHF_BOX(thisFluxDir.box()),
                          CHF_FRA1(thisvol,0),
                          CHF_FRA(thisFluxDir));
         }
      } // end loop over face directions

   } // end loop over grids
}



void
PhaseGeom::computeTangentialGradSpecial( LevelData<FluxBox>&       a_gradPhi,
                                         const LevelData<FluxBox>& a_phiFace ) const
{
   // want to compute undivided differences here
   Real fakeDx = 1.0;
   int nPhiComp = a_phiFace.nComp();
   CH_assert(a_gradPhi.nComp() >= (CFG_DIM-1)*nPhiComp && a_gradPhi.nComp() >= (VEL_DIM-1)*nPhiComp);

   const DisjointBoxLayout& grids = a_gradPhi.disjointBoxLayout();

   for (DataIterator dit(a_gradPhi.dataIterator()); dit.ok(); ++dit) {
      const CFG::MagBlockCoordSys& mag_block_coord_sys = getMagBlockCoordSys(grids[dit]);

      FluxBox& thisGradPhi = a_gradPhi[dit];
      const FluxBox& thisPhi = a_phiFace[dit];

      for (int faceDir=0; faceDir<SpaceDim; faceDir++) {
         bool configuration_dir = faceDir < CFG_DIM;
         int gradDirlo, gradDirhi;
         if (configuration_dir) {
            gradDirlo = 0;
            gradDirhi = CFG_DIM;
         }
         else {
            gradDirlo = CFG_DIM;
            gradDirhi = SpaceDim;
         }

         FArrayBox& thisGradPhiDir = thisGradPhi[faceDir];
         const FArrayBox& thisPhiDir = thisPhi[faceDir];

         // now loop over tangential directions
         for (int gradDir=gradDirlo; gradDir<gradDirhi; gradDir++) {
            if (gradDir != faceDir) {
               for (int n=0; n<nPhiComp; n++) {
                  int gradComp;
                  if (configuration_dir) {
                     gradComp = mag_block_coord_sys.tanGradComp(faceDir, gradDir, n);
                  }
                  else {
                     gradComp = m_vel_coords.tanGradComp(faceDir-CFG_DIM, gradDir-CFG_DIM, n);
                  }

                  Box gradBox(thisPhiDir.box());
                  gradBox.grow(gradDir, -1);
                  gradBox &= thisGradPhiDir.box();
                  FORT_PSCSTANGRADFACE(CHF_FRA1(thisGradPhiDir, gradComp),
                                       CHF_CONST_FRA1(thisPhiDir, n),
                                       CHF_BOX(gradBox),
                                       CHF_INT(gradDir),
                                       CHF_CONST_REAL(fakeDx));
               } // end loop over components
            } // end if gradDir is a tangential dir
         } // end loop over grad directions
      } // end loop over face directions
   } // end loop over boxes
}



void
PhaseGeom::mappedGridDivergence(LevelData<FArrayBox>&     a_divF,
                                const LevelData<FluxBox>& a_F,
                                const bool                a_omit_NT) const
{
   // We only need enough ghost cells to compute the divergence,
   // even though the input flux may have more.
   CH_assert(a_F.ghostVect() >= a_divF.ghostVect());
   LevelData<FluxBox> FluxNormal(a_F.getBoxes(),
                                           a_F.nComp()/SpaceDim,
                                           a_divF.ghostVect());

   const DisjointBoxLayout& grids = a_divF.disjointBoxLayout();
   DataIterator dit = grids.dataIterator();

   if ( a_omit_NT ) {
     // Get the normal component of the input flux and multiply it by
     // the face area

     for (dit.begin(); dit.ok(); ++dit) {
       const PhaseBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
       RealVect dx = block_coord_sys.dx();
       FluxBox& thisFluxNormal = FluxNormal[dit];
       const FluxBox& thisFlux = a_F[dit];

       for (int dir=0; dir<SpaceDim; dir++) {
         FArrayBox& thisFluxNormal_dir = thisFluxNormal[dir];

         double face_area = 1.;
         for (int tdir=0; tdir<SpaceDim; ++tdir) {
           if (tdir != dir) face_area *= dx[tdir];
         }

         thisFluxNormal_dir.copy(thisFlux[dir],dir,0,1);
         thisFluxNormal_dir *= face_area;
       }
     }
   }
   else {
      // Compute the normal component of NTranspose times the input flux

      bool fourth_order = !m_second_order && a_F.ghostVect() >= IntVect::Unit;

      computeMetricTermProductAverage(FluxNormal, a_F, fourth_order);
   }

   // Enforce conservation
   averageAtBlockBoundaries(FluxNormal);

#ifdef FIX_RADIAL_BOUNDARY_FLUX
   if ( !m_saved_flux.isDefined() ) {
      m_saved_flux.define(FluxNormal);
   }
   else {
      setRadialBoundaryFlux(FluxNormal);
   }
#endif

   // Now compute the divergence in the usual way
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& thisDiv = a_divF[dit];
      FluxBox& thisFluxNormal = FluxNormal[dit];
      // First, set divF to 0
      thisDiv.setVal(0.0);
      // since we're computing volume*div, and the fluxes
      // are multiplied by the appropriate areas, don't need
      // a dx here.
      Real fakeDx = 1.0;

      // now loop over directions and increment with directional derivative
      for (int dir=0; dir<SpaceDim; dir++) {
         FORT_PSCSDIVERGENCE(CHF_CONST_FRA(thisFluxNormal[dir]),
                             CHF_FRA(thisDiv),
                             CHF_BOX(thisDiv.box()),
                             CHF_CONST_REAL(fakeDx),
                             CHF_INT(dir));
      }
   }
}


void
PhaseGeom::mappedGridDivergenceFromFluxNormals( LevelData<FArrayBox>&      a_divF,
                                              const LevelData<FluxBox>&  a_F ) const
{
   CH_assert(a_F.nComp() == a_divF.nComp());
   const DisjointBoxLayout& grids = a_divF.disjointBoxLayout();
   DataIterator dit = a_divF.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
     const PhaseBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
     RealVect dx = block_coord_sys.dx();
     FArrayBox& thisDiv = a_divF[dit];
     thisDiv.setVal(0.0);

     for (int dir=0; dir<SpaceDim; dir++)
     {
       const FArrayBox& thisFlux = a_F[dit][dir];
       const Box& thisBox = thisDiv.box();

       // use fortran to calculate mapped cell average div(a_f)
       // assuming a_f are face-averaged values
       FORT_FLUXNORMALDIVERGENCE(CHF_CONST_FRA(thisFlux),
                                 CHF_FRA(thisDiv),
                                 CHF_BOX(thisBox),
                                 CHF_CONST_REALVECT(dx),
                                 CHF_INT(dir));
     }
   }
}



void
PhaseGeom::fillTransverseGhosts( LevelData<FluxBox>& a_data,
                                 const bool          a_do_block_exchange ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_phase_coords.boundaries();

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      const ProblemDomain& block_domain = block_coord_sys.domain();
      const int block_number = m_phase_coords.whichBlock(grids[dit]);
      const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];

      for (int dir=0; dir<SpaceDim; ++dir) {
         Box interior = surroundingNodes(grids[dit],dir);

         if ( a_do_block_exchange ) {
            for (int tdir=0; tdir<CFG_DIM; ++tdir) {  // block interfaces are only in configuration space
               if ( tdir != dir ) {  // transverse directions only

                  // If the low or high boundaries in this transverse direction are block
                  // interfaces, then grow the interior box to trick the extrapolation utility
                  // into ignoring those sides
                  if ( this_block_boundaries[tdir].isInterface() ) {
                     interior.growLo(tdir,1);
                  }
                  if ( this_block_boundaries[tdir + SpaceDim].isInterface() ) {
                     interior.growHi(tdir,1);
                  }
               }
            }
         }

         secondOrderTransExtrapAtDomainBdry(a_data[dit][dir],
                                            dir,
                                            interior,
                                            block_domain);
      }
   }

   a_data.exchange();

   if ( a_do_block_exchange ) {
      exchangeTransverseAtBlockBoundaries(a_data);
   }
}



void
PhaseGeom::averageAtBlockBoundaries(LevelData<FluxBox>& a_data) const
{
   if ( m_coordSysPtr->numBlocks() > 1 ) {

      const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

      RefCountedPtr<MultiBlockCoordSys> coordSysRCP((MultiBlockCoordSys*)m_coordSysPtr);
      coordSysRCP.neverDelete();

      BlockRegister blockRegister(coordSysRCP, grids, 0);

      DataIterator dit = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
         for (int idir = 0; idir < SpaceDim; idir++) {
            for (SideIterator sit; sit.ok(); ++sit) {
               Side::LoHiSide side = sit();
               if (blockRegister.hasInterface(dit(), idir, side)) {
                  FArrayBox flux_comp(a_data[dit][idir].box(), a_data.nComp());
                  flux_comp.copy(a_data[dit][idir]);
                  blockRegister.storeFlux(flux_comp, dit(), idir, side);
               }
            }
         }
      }
      blockRegister.close();

      const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries = m_coordSysPtr->boundaries();

      for (dit.begin(); dit.ok(); ++dit) {
         const Box& baseBox = grids[dit];
         int block_number = m_coordSysPtr->whichBlock(grids[dit]);
         int faceID = 0;
         for (SideIterator sit; sit.ok(); ++sit) {
            Side::LoHiSide side = sit();
            for (int idir = 0; idir < SpaceDim; idir++) {
               if (blockRegister.hasInterface(dit(), idir, side)) {
                  // maybe better if this is done inside BlockRegister
                  const BlockBoundary& bb = boundaries[block_number][faceID];
                  int reorientFace = bb.reorientFace(idir);
                  Box faceBox = adjCellBox(baseBox, idir, side, 1);
                  // if Lo, then shift +1; if Hi, then shift -1
                  faceBox.shiftHalf(idir, -sign(side));
                  Side::LoHiSide sideOther = flip(side);
                  // Need to define these FABs.
                  FArrayBox fluxThisFab(faceBox, a_data.nComp());
                  FArrayBox fluxOtherFab(faceBox, a_data.nComp());
                  blockRegister.getFlux(fluxThisFab, dit(),
                                        idir, side, side);
                  fluxThisFab.mult(reorientFace * 0.5);
                  blockRegister.getFlux(fluxOtherFab, dit(),
                                        idir, side, sideOther);
                  fluxOtherFab.mult(0.5);
                  fluxThisFab += fluxOtherFab;
                  a_data[dit][idir].copy(fluxThisFab);
               }
               faceID++;
            } // iterate over dimensions
         } // iterate over sides
      }
   }
}


void
PhaseGeom::exchangeTransverseAtBlockBoundaries(LevelData<FluxBox>& a_data) const
{
   if ( m_coordSysPtr->numBlocks() > 1 ) {

      CH_assert(a_data.ghostVect() >= IntVect::Unit);

      const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

      for ( int tdir=0; tdir<SpaceDim; ++tdir ) {   // Transverse directions

         m_exchange_transverse_block_register->setToZero( a_data.nComp() );

         for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
            for (int idir = 0; idir < CFG_DIM; idir++) {
               if ( idir != tdir ) {
                  for (SideIterator sit; sit.ok(); ++sit) {
                     Side::LoHiSide side = sit();
                     if (m_exchange_transverse_block_register->hasInterface(dit(), idir, side)) {
                        FArrayBox flux_comp(a_data[dit][tdir].box(), a_data.nComp());
                        flux_comp.copy(a_data[dit][tdir]);
                        flux_comp.shiftHalf(tdir,1);
                        flux_comp.shiftHalf(idir, sign(side));
                        m_exchange_transverse_block_register->storeAux(flux_comp, dit(), idir, side);
                     }
                  }
               }
            }
         }
         m_exchange_transverse_block_register->close();

         for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
            for (int idir = 0; idir < CFG_DIM; idir++) {
               if ( idir != tdir ) {
                  for (SideIterator sit; sit.ok(); ++sit) {
                     Side::LoHiSide side = sit();
                     if (m_exchange_transverse_block_register->hasInterface(dit(), idir, side)) {
                        Box faceBox = adjCellBox(grids[dit], idir, side, 1);
                        faceBox.growHi(tdir,1);
                        faceBox.shiftHalf(idir, -sign(side));
                        FArrayBox fluxOtherFab(faceBox, a_data.nComp());
                        m_exchange_transverse_block_register->getAux(fluxOtherFab, dit(),
                                                                     idir, side, side);
                        fluxOtherFab.shiftHalf(idir, sign(side));
                        fluxOtherFab.shiftHalf(tdir,-1);
                        a_data[dit][tdir].copy(fluxOtherFab);
                     }
                  } // iterate over sides
               }
            } // iterate over dimensions
         }
      }

      a_data.exchange();
   }
}


void
PhaseGeom::getBStarParallel( LevelData<FArrayBox>& a_BStarParallel ) const
{
   CH_assert(a_BStarParallel.isDefined());

   for (DataIterator dit(m_BStarParallel.dataIterator()); dit.ok(); ++dit) {
      a_BStarParallel[dit].copy(m_BStarParallel[dit]);
   }
}

void
PhaseGeom::multBStarParallel( LevelData<FArrayBox>& a_dfn ) const
{
   const DisjointBoxLayout& grids = a_dfn.disjointBoxLayout();

   if ( grids.compatible( m_gridsFull ) ) {

      DataIterator dit = m_BStarParallel.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
        for ( int n=0; n<a_dfn.nComp(); ++n ) {
           a_dfn[dit].mult(m_BStarParallel[dit], 0, n, 1);
        }
      }
   }
   else {

      Copier copier;
      copier.ghostDefine(m_gridsFull, grids, m_domain, m_ghostVect);

      // LevelData to copy to.  Since the copier just defined only
      // copies into valid cells, don't bother making any ghost cells
      LevelData<FArrayBox> tmp( grids, 1, IntVect::Zero );

      m_BStarParallel.copyTo(tmp, copier);

      DataIterator dit = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
        for ( int n=0; n<a_dfn.nComp(); ++n ) {
          a_dfn[dit].mult(tmp[dit], 0, n, 1);
        }
      }
   }
}


void PhaseGeom::multBStarParallel( LevelData<FArrayBox>& dfn,
                                   const BoundaryBoxLayout& bdry_layout ) const
{
   const DisjointBoxLayout& bdry_grids( bdry_layout.disjointBoxLayout() );
   DataIterator dit( bdry_grids.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit) {
      const DataIndex bdin( bdry_layout.dataIndex( dit ) );
      FArrayBox tmp_B( dfn[dit].box(), dfn[dit].nComp() );
      tmp_B.copy( m_BStarParallel[bdin], dfn[dit].box() );
      dfn[dit] *= tmp_B;
   }
}


void
PhaseGeom::divideBStarParallel( LevelData<FArrayBox>& a_dfn ) const
{
   const DisjointBoxLayout& grids = a_dfn.disjointBoxLayout();

   if ( grids.compatible( m_gridsFull ) ) {

      DataIterator dit = m_BStarParallel.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
        for ( int n=0; n<a_dfn.nComp(); ++n ) {
           a_dfn[dit].divide(m_BStarParallel[dit], 0, n, 1);
        }
      }
   }
   else {

      Copier copier;
      copier.ghostDefine(m_gridsFull, grids, m_domain, m_ghostVect);

      // LevelData to copy to.  Since the copier just defined only
      // copies into valid cells, don't bother making any ghost cells
      LevelData<FArrayBox> tmp( grids, 1, IntVect::Zero );

      m_BStarParallel.copyTo(tmp, copier);

      DataIterator dit = grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit) {
        for ( int n=0; n<a_dfn.nComp(); ++n ) {
          a_dfn[dit].divide(tmp[dit], 0, n, 1);
        }
      }
   }
}



void
PhaseGeom::getJ( const LevelData<FArrayBox>& a_configuration_J,
                 const LevelData<FArrayBox>& a_velocity_J,
                 LevelData<FArrayBox>&       a_J ) const
{
   const DisjointBoxLayout& grids = a_J.disjointBoxLayout();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_J = a_J[dit];
      Box this_box = this_J.box();

      CH_assert(a_configuration_J[dit].box().contains(configurationFlatten(grids[dit],this_box)));
      CH_assert(a_velocity_J[dit].box().contains(velocityFlatten(grids[dit],this_box)));

      this_J.setVal(1.);

      FORT_MULT_CFG(CHF_BOX(this_box),
                    CHF_FRA1(a_configuration_J[dit],0),
                    CHF_FRA(this_J));
      FORT_MULT_VEL(CHF_BOX(this_box),
                    CHF_FRA1(a_velocity_J[dit],0),
                    CHF_FRA(this_J));

   }

   a_J.exchange();
}


void
PhaseGeom::multJonValid( LevelData<FArrayBox>& a_dfn ) const
{
   // Multiplies the argument by J on valid cells

   const DisjointBoxLayout& grids = a_dfn.disjointBoxLayout();

   if ( m_second_order ) {

      LevelData<FArrayBox> J(grids, 1, IntVect::Zero);
      getJ(*m_configuration_J, *m_velocity_J, J);

      for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
         a_dfn[dit] *= J[dit];
      }
   }
   else {

      IntVect grown_vect = IntVect::Unit;
      LevelData<FArrayBox> J(grids, 1, grown_vect);
      getJ(*m_configuration_J, *m_velocity_J, J);

      LevelData<FArrayBox> grown_dfn(grids, 1, grown_vect);

      for (DataIterator dit(grown_dfn.dataIterator()); dit.ok(); ++dit) {
         grown_dfn[dit].copy(a_dfn[dit]);
      }
      grown_dfn.exchange();

      fourthOrderCellProd(a_dfn, J, grown_dfn);

      // Compute the fourth-order product
      for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
         const PhaseBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
         fourthOrderCellProd(a_dfn[dit], grown_dfn[dit], J[dit], grids[dit], coord_sys.domain(), true);
      }
   }
}


void PhaseGeom::multPointwiseJ( LevelData<FArrayBox>& a_u,
                                const BoundaryBoxLayout& a_bdry_layout ) const
{
   const DisjointBoxLayout& bdry_grids( a_bdry_layout.disjointBoxLayout() );
   for (DataIterator dit( bdry_grids.dataIterator() ); dit.ok(); ++dit) {
      const Box interior_box( a_bdry_layout.interiorBox( dit ) );
      const PhaseBlockCoordSys& coord_sys( getBlockCoordSys( interior_box ) );
      const Box& fill_box( a_u[dit].box() );
      FArrayBox J( fill_box, 1 );
      coord_sys.pointwiseJ( J, fill_box );
      a_u[dit] *= J;
   }
}


void
PhaseGeom::divideJonValid( LevelData<FArrayBox>& a_dfn ) const
{
   // Divides the argument by J on valid cells

   const DisjointBoxLayout& grids = a_dfn.disjointBoxLayout();

   if ( m_second_order ) {

      LevelData<FArrayBox> J(grids, 1, IntVect::Zero);
      getJ(*m_configuration_J, *m_velocity_J, J);

      for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
         a_dfn[dit] /= J[dit];
      }
   }
   else {

      IntVect grown_vect = IntVect::Unit;
      LevelData<FArrayBox> J(grids, 1, grown_vect);
      getJ(*m_configuration_J, *m_velocity_J, J);

      LevelData<FArrayBox> grown_dfn(grids, 1, grown_vect);

      for (DataIterator dit(grown_dfn.dataIterator()); dit.ok(); ++dit) {
         grown_dfn[dit].copy(a_dfn[dit]);
      }
      grown_dfn.exchange();

      // Compute fourth-order quotient
      for (DataIterator dit(a_dfn.dataIterator()); dit.ok(); ++dit) {
         const PhaseBlockCoordSys& coord_sys = getBlockCoordSys(grids[dit]);
         cellFGToCellF(a_dfn[dit], grown_dfn[dit], J[dit], grids[dit], coord_sys.domain(), true);
      }
   }
}


void
PhaseGeom::getCellVolumes( LevelData<FArrayBox>& a_volume ) const
{
   const DisjointBoxLayout& grids = a_volume.disjointBoxLayout();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_volume = a_volume[dit];
      const Box& box(this_volume.box());

      CH_assert((*m_configuration_volumes)[dit].box().contains(configurationFlatten(grids[dit],box)));
      CH_assert((*m_velocity_volumes)[dit].box().contains(velocityFlatten(grids[dit],box)));

      this_volume.setVal(1.);

      FORT_MULT_CFG(CHF_BOX(box),
                    CHF_FRA1((*m_configuration_volumes)[dit],0),
                    CHF_FRA(this_volume));
      FORT_MULT_VEL(CHF_BOX(box),
                    CHF_FRA1((*m_velocity_volumes)[dit],0),
                    CHF_FRA(this_volume));
   }
   a_volume.exchange();
}



void
PhaseGeom::multFaceAreas( LevelData<FluxBox>& a_data ) const
{
   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FluxBox& this_data = a_data[dit];
      Box this_box = this_data.box();

      // Note: cfg_volumes includes a 2piR factor for axisymmetric systems

      const FArrayBox& cfg_volumes = (*m_configuration_volumes)[dit];
      const FArrayBox& vel_volumes = (*m_velocity_volumes)[dit];

      const FluxBox& cfg_face_areas = (*m_configuration_face_areas)[dit];
      CH_assert(cfg_face_areas.box().contains(configurationFlatten(grids[dit],this_box)));

      const FluxBox& vel_face_areas = (*m_velocity_face_areas)[dit];
      CH_assert(vel_face_areas.box().contains(velocityFlatten(grids[dit],this_box)));

      for (int dir=0; dir<CFG_DIM; ++dir) {
        FArrayBox& this_data_dir = this_data[dir];
        const Box& this_box_dir = this_data_dir.box();

        FORT_MULT_CFG(CHF_BOX(this_box_dir),
                      CHF_FRA1(cfg_face_areas[dir],0),
                      CHF_FRA(this_data_dir));
        FORT_MULT_VEL(CHF_BOX(this_box_dir),
                      CHF_FRA1(vel_volumes,0),
                      CHF_FRA(this_data_dir));
      }

      for (int dir=CFG_DIM; dir<SpaceDim; ++dir) {
        FArrayBox& this_data_dir = this_data[dir];
        const Box& this_box_dir = this_data_dir.box();

        FORT_MULT_VEL(CHF_BOX(this_box_dir),
                      CHF_FRA1(vel_face_areas[dir],0),
                      CHF_FRA(this_data_dir));
        FORT_MULT_CFG(CHF_BOX(this_box_dir),
                      CHF_FRA1(cfg_volumes,0),
                      CHF_FRA(this_data_dir));
      }
   }
}


void PhaseGeom::zeroBoundaryFlux( const int           a_dir,
                                  LevelData<FluxBox>& a_flux ) const
{
   const DisjointBoxLayout grids = a_flux.disjointBoxLayout();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_phase_coords.boundaries();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      const Box& domainBox = block_coord_sys.domain().domainBox();
      const int block_number = m_phase_coords.whichBlock(grids[dit]);
      const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];

      FArrayBox& this_flux = a_flux[dit][a_dir];

      if ( this_block_boundaries[a_dir].isDomainBoundary() ) {
         Box zero_box = bdryLo(surroundingNodes(domainBox,a_dir), a_dir, 1);
         zero_box &= this_flux.box();

         if (zero_box.ok()) {
            this_flux.setVal(0., zero_box, 0, this_flux.nComp()-1);
         }
      }

      if ( this_block_boundaries[a_dir + SpaceDim].isDomainBoundary() ) {
         Box zero_box = bdryHi(surroundingNodes(domainBox,a_dir), a_dir, 1);
         zero_box &= this_flux.box();

         if (zero_box.ok()) {
            this_flux.setVal(0., zero_box, 0, this_flux.nComp()-1);
         }
      }
   }
}



void PhaseGeom::setRadialBoundaryFlux( LevelData<FluxBox>& a_flux ) const
{
   const DisjointBoxLayout grids = a_flux.disjointBoxLayout();
   const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& block_boundaries = m_phase_coords.boundaries();

   int dir = RADIAL_DIR;

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const PhaseBlockCoordSys& block_coord_sys = getBlockCoordSys(grids[dit]);
      const Box& domainBox = block_coord_sys.domain().domainBox();
      const int block_number = m_phase_coords.whichBlock(grids[dit]);
      const Tuple<BlockBoundary, 2*SpaceDim>& this_block_boundaries = block_boundaries[block_number];

      FArrayBox& this_flux = a_flux[dit][dir];
      const FArrayBox& this_saved_flux = m_saved_flux[dit][dir];

      if ( this_block_boundaries[dir].isDomainBoundary() ) {
         Box copy_box = bdryLo(surroundingNodes(domainBox,dir), dir, 1);
         copy_box &= this_flux.box();

         if (copy_box.ok()) {
            this_flux.copy(this_saved_flux, copy_box);
         }
      }

      if ( this_block_boundaries[dir + SpaceDim].isDomainBoundary() ) {
         Box copy_box = bdryHi(surroundingNodes(domainBox,dir), dir, 1);
         copy_box &= this_flux.box();

         if (copy_box.ok()) {
            this_flux.copy(this_saved_flux, copy_box);
         }
      }
   }
}


void
PhaseGeom::fillInternalGhosts( LevelData<FArrayBox>& a_data ) const
{
   CH_assert(a_data.ghostVect() >= m_ghostVect);

   a_data.exchange();

   if (m_mblexPtr && !(m_mag_geom.extrablockExchange())) {
      m_mblexPtr->interpGhosts(a_data);
   }
   
   else if (m_mag_geom.extrablockExchange()) {
      exchangeExtraBlockGhosts(a_data);
      a_data.exchange();
   }
}

void
PhaseGeom::exchangeExtraBlockGhosts( LevelData<FArrayBox>& a_data ) const
{
   /*
    Fills extrablock ghost cells using a pure exchange from valid data
    rather than interpolation (for experimental purposes only).
    */
   
   if ( m_coordSysPtr->numBlocks() > 1 ) {
      
      RefCountedPtr<MultiBlockCoordSys> coordSysRCP((MultiBlockCoordSys*)m_coordSysPtr);
      coordSysRCP.neverDelete();
      
      const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
      
      int ncomp = a_data.nComp();
      
      IntVect ghost_vect = a_data.ghostVect();
      int max_ghost = -1;
      for (int dir=0; dir<SpaceDim; ++dir) {
         if (ghost_vect[dir] > max_ghost) max_ghost = ghost_vect[dir];
      }
      
      BlockRegister blockRegister(coordSysRCP, grids, max_ghost);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir = 0; dir < SpaceDim; dir++) {
            IntVect grow_vect = ghost_vect;
            grow_vect[dir] = 0;
            int nghost = ghost_vect[dir];
            for (SideIterator sit; sit.ok(); ++sit) {
               Side::LoHiSide side = sit();
               if (blockRegister.hasInterface(dit(), dir, side)) {
                  Box fill_box = adjCellBox(grids[dit], dir, side, 1);
                  fill_box.grow(grow_vect);
                  FArrayBox temp(fill_box, ncomp*nghost);
                  for (int n=0; n<nghost; ++n) {
                     temp.shift(dir, -sign(side));
                     temp.copy(a_data[dit],0,n*ncomp,ncomp);
                  }
                  temp.shiftHalf(dir, (2*nghost-1)*sign(side));
                  blockRegister.storeAux(temp, dit(), dir, side);
               }
            }
         }
      }
      blockRegister.close();
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir = 0; dir < SpaceDim; dir++) {
            IntVect grow_vect = ghost_vect;
            grow_vect[dir] = 0;
            int nghost = ghost_vect[dir];
            for (SideIterator sit; sit.ok(); ++sit) {
               Side::LoHiSide side = sit();
               if (blockRegister.hasInterface(dit(), dir, side)) {
                  Box fill_box = adjCellBox(grids[dit], dir, side, -1);
                  fill_box.grow(grow_vect);
                  FArrayBox temp(fill_box, ncomp*nghost);
                  temp.shiftHalf(dir, sign(side));
                  blockRegister.getAux(temp, dit(),
                                       dir, side, side);
                  temp.shiftHalf(dir, -sign(side));
                  for (int n=0; n<nghost; ++n) {
                     temp.shift(dir, sign(side));
                     a_data[dit].copy(temp,n*ncomp,0,ncomp);
                  }
               }
            } // iterate over dimensions
         } // iterate over sides
      }
   }
}




/*
  Transdimensional utilities
*/


void
PhaseGeom::injectConfigurationToPhase( const CFG::LevelData<CFG::FArrayBox>& a_src,
                                       LevelData<FArrayBox>&                 a_dst ) const
{
   if (a_dst.isDefined()) {
      MayDay::Error("PhaseBlockCoordSys: a_dst LevelData passed to injectConfigurationToPhase is already defined");
   }

   const Box& domainBox = m_domain.domainBox();

   // locate these at the low end of the vp- and mu- domains
   const CP1::SliceSpec slice_vp(VPARALLEL_DIR,domainBox.smallEnd(VPARALLEL_DIR));
   const SliceSpec slice_mu(MU_DIR,domainBox.smallEnd(MU_DIR));

   IntVect ghostVect = config_inject(a_src.ghostVect());

   // flatten phase space grid in the mu direction
   DisjointBoxLayout mu_flattened_dbl;
   adjCellLo(mu_flattened_dbl, m_gridsFull, MU_DIR, -1);

   // then flatten it in the vparallel coordinate
   DisjointBoxLayout vpmu_flattened_dbl;
   adjCellLo(vpmu_flattened_dbl, mu_flattened_dbl, VPARALLEL_DIR, -1);

   // create CP1 injection of CFG src data
   CP1::LevelData<CP1::FArrayBox> CP1_temp;
   injectLevelData(CP1_temp, a_src, slice_vp);

   // create injection of CP1 data
   LevelData<FArrayBox> vpmu_flattened_dst;
   injectLevelData(vpmu_flattened_dst, CP1_temp, slice_mu);

   // define destination
   a_dst.define(vpmu_flattened_dbl, a_src.nComp(), ghostVect);

   Vector<int> spreadingDirs;
   spreadingDirs.push_back(VPARALLEL_DIR);
   spreadingDirs.push_back(MU_DIR);

   const SpreadingOp spreadOp(spreadingDirs);

   for (int block=0; block<m_coordSysPtr->numBlocks(); ++block) {

     const ProblemDomain& thisMappingDomain
       = ((const PhaseBlockCoordSys&)(*m_coordSysPtr->getCoordSys(block))).domain();

     SpreadingCopier spreadCopier(vpmu_flattened_dst.getBoxes(),
                                  vpmu_flattened_dbl,
                                  thisMappingDomain,
                                  ghostVect,
                                  spreadingDirs);

     // spread injected data in the vp and mu directions
     vpmu_flattened_dst.copyTo(vpmu_flattened_dst.interval(),
                               a_dst, a_dst.interval(),
                               spreadCopier,
                               spreadOp);
   }

   a_dst.exchange();
}



void
PhaseGeom::injectConfigurationToPhase( const CFG::LevelData<CFG::FluxBox>& a_src,
                                       LevelData<FluxBox>&                 a_dst ) const
{
   if (a_dst.isDefined()) {
      MayDay::Error("PhaseBlockCoordSys: a_dst LevelData passed to injectConfigurationToPhase is already defined");
   }

   const Box& domainBox = m_domain.domainBox();
   const CP1::SliceSpec slice_vp(VPARALLEL_DIR,domainBox.smallEnd(VPARALLEL_DIR));
   const SliceSpec slice_mu(MU_DIR,domainBox.smallEnd(MU_DIR));

   IntVect ghostVect = config_inject(a_src.ghostVect());

   // flatten phase space grid in the mu direction
   DisjointBoxLayout mu_flattened_dbl;
   adjCellLo(mu_flattened_dbl, m_gridsFull, MU_DIR, -1);

   // then flatten it in the vparallel coordinate
   DisjointBoxLayout vpmu_flattened_dbl;
   adjCellLo(vpmu_flattened_dbl, mu_flattened_dbl, VPARALLEL_DIR, -1);

   // create CP1 injection of CFG src data
   CP1::LevelData<CP1::FluxBox> CP1_temp;
   injectLevelFlux(CP1_temp, a_src, slice_vp);

   // create injection of CP1 data
   LevelData<FluxBox> vpmu_flattened_dst;
   injectLevelFlux(vpmu_flattened_dst, CP1_temp, slice_mu);

   // Compute cell-centered data to set values on new directions

   CFG::LevelData<CFG::FArrayBox> src_cell(a_src.getBoxes(), a_src.nComp(), a_src.ghostVect());

   CFG::DataIterator cdit = src_cell.dataIterator();
   for (cdit.begin(); cdit.ok(); ++cdit) {
      const CFG::FluxBox& this_src = a_src[cdit];
      CFG::FArrayBox& this_src_cell = src_cell[cdit];
      this_src_cell.setVal(0.);
      for (int dir=0; dir<CFG_DIM; dir++) {
         const CFG::FArrayBox& this_src_dir = this_src[dir];
         CFG::FArrayBox shift_fab(this_src_dir.box(), this_src_dir.nComp());
         shift_fab.copy(this_src_dir);

         shift_fab.shiftHalf(dir,-1);
         this_src_cell.plus(shift_fab);
         shift_fab.shiftHalf(dir,2);
         this_src_cell.plus(shift_fab);
      }
      this_src_cell.mult(0.5/CFG_DIM);
   }

   DataIterator dit = vpmu_flattened_dst.dataIterator();
   for (dit.begin(),cdit.begin(); dit.ok()&&cdit.ok(); ++dit,++cdit) {

      CP1::BaseFab<Real> CP1_injected_fab;
      injectBaseFab(CP1_injected_fab, src_cell[cdit], slice_vp);
      BaseFab<Real> injected_fab;
      injectBaseFab(injected_fab, CP1_injected_fab, slice_mu);

      injected_fab.shiftHalf(VPARALLEL_DIR,-1);
      vpmu_flattened_dst[dit][VPARALLEL_DIR].copy(injected_fab);
      injected_fab.shiftHalf(VPARALLEL_DIR,1);

      injected_fab.shiftHalf(MU_DIR,-1);
      vpmu_flattened_dst[dit][MU_DIR].copy(injected_fab);
      injected_fab.shiftHalf(MU_DIR,1);
   }

   // define destination
   a_dst.define(vpmu_flattened_dbl, a_src.nComp(), ghostVect);

   Vector<int> spreadingDirs;
   spreadingDirs.push_back(VPARALLEL_DIR);
   spreadingDirs.push_back(MU_DIR);

   const FaceSpreadingOp spreadOp(spreadingDirs);

   for (int block=0; block<m_coordSysPtr->numBlocks(); ++block) {

     const ProblemDomain& thisMappingDomain
       = ((const PhaseBlockCoordSys&)(*m_coordSysPtr->getCoordSys(block))).domain();

     SpreadingCopier spreadCopier(vpmu_flattened_dst.getBoxes(),
                                  vpmu_flattened_dbl,
                                  thisMappingDomain,
                                  ghostVect,
                                  spreadingDirs);

     // spread injected data in the vp and mu directions
     vpmu_flattened_dst.copyTo(vpmu_flattened_dst.interval(),
                               a_dst, a_dst.interval(),
                               spreadCopier,
                               spreadOp);
   }

   a_dst.exchange();
}



void
PhaseGeom::injectConfigurationToPhase( const CFG::LevelData<CFG::FluxBox>&   a_src,
                                       const CFG::LevelData<CFG::FArrayBox>& a_src_cell,
                                       LevelData<FluxBox>&                   a_dst ) const
{
   const Box& domainBox = m_domain.domainBox();
   const CP1::SliceSpec slice_vp(VPARALLEL_DIR,domainBox.smallEnd(VPARALLEL_DIR));
   const SliceSpec slice_mu(MU_DIR,domainBox.smallEnd(MU_DIR));

   IntVect ghostVect = config_inject(a_src.ghostVect());

   // flatten phase space grid in the mu direction
   DisjointBoxLayout mu_flattened_dbl;
   adjCellLo(mu_flattened_dbl, m_gridsFull, MU_DIR, -1);

   // then flatten it in the vparallel coordinate
   DisjointBoxLayout vpmu_flattened_dbl;
   adjCellLo(vpmu_flattened_dbl, mu_flattened_dbl, VPARALLEL_DIR, -1);

   // create CP1 injection of CFG src data
   CP1::LevelData<CP1::FluxBox> CP1_temp;
   injectLevelFlux(CP1_temp, a_src, slice_vp);

   // create injection of CP1 data
   LevelData<FluxBox> vpmu_flattened_dst;
   injectLevelFlux(vpmu_flattened_dst, CP1_temp, slice_mu);

   DataIterator dit = vpmu_flattened_dst.dataIterator();
   CFG::DataIterator cdit = a_src_cell.dataIterator();
   for (dit.begin(),cdit.begin(); dit.ok()&&cdit.ok(); ++dit,++cdit) {

      CP1::BaseFab<Real> CP1_injected_fab;
      injectBaseFab(CP1_injected_fab, a_src_cell[cdit], slice_vp);
      BaseFab<Real> injected_fab;
      injectBaseFab(injected_fab, CP1_injected_fab, slice_mu);

      injected_fab.shiftHalf(VPARALLEL_DIR,-1);
      vpmu_flattened_dst[dit][VPARALLEL_DIR].copy(injected_fab);
      injected_fab.shiftHalf(VPARALLEL_DIR,1);

      injected_fab.shiftHalf(MU_DIR,-1);
      vpmu_flattened_dst[dit][MU_DIR].copy(injected_fab);
      injected_fab.shiftHalf(MU_DIR,1);
   }

   // define destination
   a_dst.define(vpmu_flattened_dbl, a_src.nComp(), ghostVect);

   Vector<int> spreadingDirs;
   spreadingDirs.push_back(VPARALLEL_DIR);
   spreadingDirs.push_back(MU_DIR);

   const FaceSpreadingOp spreadOp(spreadingDirs);

   for (int block=0; block<m_coordSysPtr->numBlocks(); ++block) {

     const ProblemDomain& thisMappingDomain
       = ((const PhaseBlockCoordSys&)(*m_coordSysPtr->getCoordSys(block))).domain();

     SpreadingCopier spreadCopier(vpmu_flattened_dst.getBoxes(),
                                  vpmu_flattened_dbl,
                                  thisMappingDomain,
                                  ghostVect,
                                  spreadingDirs);

     // spread injected data in the vp and mu directions
     vpmu_flattened_dst.copyTo(vpmu_flattened_dst.interval(),
                               a_dst, a_dst.interval(),
                               spreadCopier,
                               spreadOp);
   }

   a_dst.exchange();
}



void
PhaseGeom::injectConfigurationToPhase( const CFG::FArrayBox& a_src,
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
PhaseGeom::injectConfigurationToPhase( const CFG::FluxBox& a_src,
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
PhaseGeom::injectVelocityToPhase( const VEL::LevelData<VEL::FArrayBox>& a_src,
                                  LevelData<FArrayBox>&                 a_dst ) const
{
   if (a_dst.isDefined()) {
      MayDay::Error("PhaseBlockCoordSys: a_dst LevelData passed to injectVelocityToPhase is already defined");
   }

   const Box& domainBox = m_domain.domainBox();

   // locate these at the low end of the x- and y- domains
   const D3::SliceSpec slice_x(RADIAL_DIR,domainBox.smallEnd(RADIAL_DIR));
   const SliceSpec slice_y(POLOIDAL_DIR,domainBox.smallEnd(POLOIDAL_DIR));

   IntVect ghostVect = vel_inject(a_src.ghostVect());

   // flatten phase space grid in the x direction
   DisjointBoxLayout x_flattened_dbl;
   adjCellLo(x_flattened_dbl, m_gridsFull, RADIAL_DIR, -1);

   // then flatten it in the y direction
   DisjointBoxLayout xy_flattened_dbl;
   adjCellLo(xy_flattened_dbl, x_flattened_dbl, POLOIDAL_DIR, -1);

   // create x injection of src data
   D3::LevelData<D3::FArrayBox> x_injected_src;
   injectLevelData(x_injected_src, a_src, slice_x);

   // create xy injection of src data
   LevelData<FArrayBox> xy_injected_src;
   injectLevelData(xy_injected_src, x_injected_src, slice_y);

   Vector<int> spreadingDirs;
   spreadingDirs.push_back(RADIAL_DIR);
   spreadingDirs.push_back(POLOIDAL_DIR);
   SpreadingCopier spreadCopier(xy_injected_src.getBoxes(),
                                xy_flattened_dbl,
                                m_domain,
                                ghostVect,
                                spreadingDirs);
   const SpreadingOp spreadOp(spreadingDirs);

   // define the destination using the flattened dbl
   a_dst.define(xy_flattened_dbl, a_src.nComp(), ghostVect);

   // spread injected data in the x and y directions
   xy_injected_src.copyTo(xy_injected_src.interval(),
                          a_dst, a_dst.interval(),
                          spreadCopier,
                          spreadOp);

   a_dst.exchange();
}



void
PhaseGeom::injectVelocityToPhase( const VEL::LevelData<VEL::FluxBox>& a_src,
                                  LevelData<FluxBox>&                 a_dst ) const
{
   if (a_dst.isDefined()) {
      MayDay::Error("PhaseBlockCoordSys: a_dst LevelData passed to injectVelocityToPhase is already defined");
   }

#if CFG_DIM==2
   const Box domainBox= m_domain.domainBox();

   // locate these at the low end of the x- and y- domains
   const D3::SliceSpec slice_x(RADIAL_DIR,domainBox.smallEnd(RADIAL_DIR));
   const SliceSpec slice_y(POLOIDAL_DIR,domainBox.smallEnd(POLOIDAL_DIR));

   IntVect ghostVect = vel_inject(a_src.ghostVect());

   // flatten phase space grid in the x direction
   DisjointBoxLayout x_flattened_dbl;
   adjCellLo(x_flattened_dbl, m_gridsFull, RADIAL_DIR, -1);

   // then flatten it in the y direction
   DisjointBoxLayout xy_flattened_dbl;
   adjCellLo(xy_flattened_dbl, x_flattened_dbl, POLOIDAL_DIR, -1);

   // define the destination using the flattened dbl
   a_dst.define(xy_flattened_dbl, a_src.nComp(), ghostVect);

   // create x injection of src data
   D3::LevelData<D3::FluxBox> x_injected_src;
   injectLevelFlux(x_injected_src, a_src, slice_x);

   // create xy injection of src data
   LevelData<FluxBox> xy_injected_src;
   injectLevelFlux(xy_injected_src, x_injected_src, slice_y);

   Vector<int> spreadingDirs;
   spreadingDirs.push_back(RADIAL_DIR);
   spreadingDirs.push_back(POLOIDAL_DIR);
   SpreadingCopier spreadCopier(xy_injected_src.getBoxes(),
                                xy_flattened_dbl,
                                m_domain,
                                ghostVect,
                                spreadingDirs);
   const FaceSpreadingOp spreadOp(spreadingDirs);

   // spread injected data in the x and y directions
   xy_injected_src.copyTo(xy_injected_src.interval(),
                          a_dst, a_dst.interval(),
                          spreadCopier,
                          spreadOp);

   a_dst.exchange();
#else
      MayDay::Error("PhaseBlockCoordSys:injectVelocityToPhase only implemented for 2D configuration space");
#endif
}


void
PhaseGeom::projectPhaseToConfiguration( const LevelData<FArrayBox>&     a_src,
                                        CFG::LevelData<CFG::FArrayBox>& a_dst ) const
{
   const Box& domainBox = m_domain.domainBox();

   // Slice in the mu direction at the low mu coordinate

   SliceSpec slice_mu(MU_DIR,domainBox.smallEnd(MU_DIR));
   CP1::LevelData<CP1::FArrayBox> temp1;
   sliceLevelDataLocalOnly(temp1, a_src, slice_mu);

   // Slice in the v_parallel direction at the low v_parallel coordinate

   CP1::SliceSpec slice_vp(VPARALLEL_DIR,domainBox.smallEnd(VPARALLEL_DIR));
   CFG::LevelData<CFG::FArrayBox> temp2;
   sliceLevelDataLocalOnly(temp2, temp1, slice_vp);

   const CFG::DisjointBoxLayout& src_dbl = temp2.disjointBoxLayout();
   const CFG::DisjointBoxLayout& dst_dbl = a_dst.disjointBoxLayout();

   const CFG::ProblemDomain& config_domain = dst_dbl.physDomain();

   CFG::Copier copier;
   copier.ghostDefine(src_dbl,
                      dst_dbl,
                      config_domain,
                      temp2.ghostVect(),
                      a_dst.ghostVect());

   temp2.copyTo(a_dst, copier);
}


void
PhaseGeom::projectPhaseToConfiguration( const LevelData<FluxBox>&     a_src,
                                       CFG::LevelData<CFG::FluxBox>& a_dst ) const
{
    const Box& domainBox = m_domain.domainBox();
    
   // Slice in the mu direction at the low mu coordinate
    
   SliceSpec slice_mu(MU_DIR,domainBox.smallEnd(MU_DIR));
   CP1::LevelData<CP1::FluxBox> temp1;
   sliceLevelFluxLocalOnly(temp1, a_src, slice_mu);
    
   // Slice in the v_parallel direction at the low v_parallel coordinate
    
   CP1::SliceSpec slice_vp(VPARALLEL_DIR,domainBox.smallEnd(VPARALLEL_DIR));
   CFG::LevelData<CFG::FluxBox> temp2;
   sliceLevelFluxLocalOnly(temp2, temp1, slice_vp);

   const CFG::DisjointBoxLayout& src_dbl = temp2.disjointBoxLayout();
   const CFG::DisjointBoxLayout& dst_dbl = a_dst.disjointBoxLayout();

   const CFG::ProblemDomain& config_domain = dst_dbl.physDomain();

   CFG::Copier copier;
   copier.ghostDefine(src_dbl,
                      dst_dbl,
                      config_domain,
                      temp2.ghostVect(),
                      a_dst.ghostVect());

   temp2.copyTo(a_dst, copier);
}


void
PhaseGeom::projectPhaseToConfiguration( const Box& a_src,
                                        CFG::Box&  a_dst ) const
{
  a_dst.define(config_restrict(a_src.smallEnd()),
               config_restrict(a_src.bigEnd()),
               config_restrict(a_src.type()));
}



void
PhaseGeom::projectPhaseToVelocity( const Box& a_src,
                                   VEL::Box&  a_dst ) const
{
  a_dst.define(vel_restrict(a_src.smallEnd()),
               vel_restrict(a_src.bigEnd()),
               vel_restrict(a_src.type()));
}



CFG::IntVect
PhaseGeom::config_restrict( const IntVect& a_x ) const
{
   CFG::IntVect restrict_x;

   for (int i=0; i<CFG_DIM; i++) {
      restrict_x[i] = a_x[i];
   }

   return restrict_x;
}



IntVect
PhaseGeom::config_inject( const CFG::IntVect& a_v ) const
{
   IntVect injected_v;

   for (int dir=0; dir<CFG_DIM; dir++) {
      injected_v[dir] = a_v[dir];
   }
   for (int dir=CFG_DIM; dir<SpaceDim; dir++) {
      injected_v[dir] = 0;
   }

   return injected_v;
}



VEL::IntVect
PhaseGeom::vel_restrict( const IntVect& a_x ) const
{
   VEL::IntVect restrict_x;

   for (int i=VPARALLEL_DIR; i<=MU_DIR; i++) {
      restrict_x[i-VPARALLEL_DIR] = a_x[i];
   }

   return restrict_x;
}



IntVect
PhaseGeom::vel_inject( const VEL::IntVect& a_v ) const
{
   IntVect injected_v;

   for (int dir=0; dir<CFG_DIM; dir++) {
      injected_v[dir] = 0;
   }
   for (int dir=CFG_DIM; dir<SpaceDim; dir++) {
      injected_v[dir] = a_v[dir-CFG_DIM];
   }

   return injected_v;
}



IntVect
PhaseGeom::tensorProduct( const CFG::IntVect& cfg_iv,
                          const VEL::IntVect& vel_iv ) const
{
  IntVect iv;

  for (int dir=0; dir<CFG_DIM; ++dir) {
    iv[dir] = cfg_iv[dir];
  }
  for (int dir=CFG_DIM; dir<SpaceDim; ++dir) {
    iv[dir] = vel_iv[dir-CFG_DIM];
  }

  return iv;
}



void
PhaseGeom::getConfigurationData( const VEL::IntVect &          vspace_index,
                                 const LevelData<FluxBox>&     src_phase,
                                 CFG::LevelData<CFG::FluxBox>& dst_config ) const
{
   // Slice in the mu direction at the current mu coordinate

   SliceSpec slice_mu(MU_DIR,vspace_index[1]);
   CP1::LevelData<CP1::FluxBox> temp;
   sliceLevelFlux(temp, src_phase, slice_mu);

   // Slice in the v_parallel direction at the current v_parallel coordinate

   CP1::SliceSpec slice_vp(VPARALLEL_DIR,vspace_index[0]);
   sliceLevelFlux(dst_config, temp, slice_vp);
}



void
PhaseGeom::getConfigurationData( const VEL::IntVect &            a_vspace_index,
                                 const LevelData<FArrayBox>&     a_src,
                                 CFG::LevelData<CFG::FArrayBox>& a_dst ) const
{
   // Slice in the mu direction at the current mu coordinate

   SliceSpec slice_mu(MU_DIR,a_vspace_index[1]);
   CP1::LevelData<CP1::FArrayBox> temp;
   sliceLevelData(temp, a_src, slice_mu);

   // Slice in the v_parallel direction at the current v_parallel coordinate

   CP1::SliceSpec slice_vp(VPARALLEL_DIR,a_vspace_index[0]);
   sliceLevelData(a_dst, temp, slice_vp);
}



void
PhaseGeom::setConfigurationData( const VEL::IntVect &                  a_vspace_index,
                                 const CFG::LevelData<CFG::FArrayBox>& a_src,
                                 LevelData<FArrayBox>&                 a_dst ) const
{
   // Inject in the v_parallel direction at the current v_parallel coordinate

   CP1::SliceSpec slice_vp(VPARALLEL_DIR,a_vspace_index[0]);
   CP1::LevelData<CP1::FArrayBox> temp;
   injectLevelData(temp, a_src, slice_vp);

   // Inject in the mu direction at the current mu coordinate

   SliceSpec slice_mu(MU_DIR,a_vspace_index[1]);
   injectLevelData(a_dst, temp, slice_mu);
}



void
PhaseGeom::getVelocityData( const CFG::IntVect &          a_cspace_index,
                            const LevelData<FluxBox>&     a_src,
                            VEL::LevelData<VEL::FluxBox>& a_dst ) const
{
   // Slice in the x direction at the current x coordinate

   SliceSpec slice_x(RADIAL_DIR,a_cspace_index[0]);
   CP1::LevelData<CP1::FluxBox> temp;
   sliceLevelFlux(temp, a_src, slice_x);

   // Slice in the y direction at the current y coordinate

   CP1::SliceSpec slice_y(POLOIDAL_DIR,a_cspace_index[1]);
   sliceLevelFlux(a_dst, temp, slice_y);
}


void
PhaseGeom::getVelocityData( const CFG::IntVect &            a_cspace_index,
                            const LevelData<FArrayBox>&     a_src,
                            VEL::LevelData<VEL::FArrayBox>& a_dst ) const
{
   // Slice in the y direction at the current y coordinate

   SliceSpec slice_y(POLOIDAL_DIR,a_cspace_index[1]);
   CP1::LevelData<CP1::FArrayBox> temp;
   sliceLevelData(temp, a_src, slice_y);

   // Slice in the x direction at the current x coordinate

   CP1::SliceSpec slice_x(RADIAL_DIR,a_cspace_index[0]);
   sliceLevelData(a_dst, temp, slice_x);
}



void
PhaseGeom::setVelocityData( const CFG::IntVect &                  a_cspace_index,
                            const VEL::LevelData<VEL::FArrayBox>& a_src,
                            LevelData<FArrayBox>&                 a_dst ) const
{
   // Inject in the x direction at the current x coordinate

   CP1::SliceSpec slice_x(RADIAL_DIR,a_cspace_index[0]);
   CP1::LevelData<CP1::FArrayBox> temp;
   injectLevelData(temp, a_src, slice_x);

   // Inject in the y direction at the current y coordinate

   SliceSpec slice_y(POLOIDAL_DIR,a_cspace_index[1]);
   injectLevelData(a_dst, temp, slice_y);
}

/*
  Plotting (MOVE THIS ELSEWHERE)
*/

void
PhaseGeom::plotAtVelocityIndex( const string               a_file_name,
                                const VEL::IntVect&        a_vspace_index,
                                const LevelData<FluxBox>&  a_data,
                                const double&              a_time ) const
{
   CH_assert( validVelocityIndex(a_vspace_index) );

   // Since the passed data is a LevelData<FluxBox>, we average the face
   // values to cell centers
   LevelData<FArrayBox> tmp(a_data.getBoxes(), SpaceDim*a_data.nComp(), a_data.ghostVect());
   EdgeToCell(a_data, tmp);

   // Evaluate data at the specified velocity space index.  We would have done
   // this before averaging faces to cell centers, but we don't have an
   // EdgeToCell function that works only in configuration space.
   CFG::LevelData<CFG::FArrayBox> tmp_cfg;
   getConfigurationData(a_vspace_index, tmp, tmp_cfg);

   const CFG::DisjointBoxLayout& grids = tmp_cfg.disjointBoxLayout();

   // Next, average the directions
   CFG::LevelData<CFG::FArrayBox> data_at_vpt(grids, a_data.nComp(), tmp_cfg.ghostVect());

   CFG::DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      CFG::FArrayBox& this_tmp_cfg = tmp_cfg[dit];
      CFG::FArrayBox& this_data = data_at_vpt[dit];

      this_data.setVal(0.);
      for (int comp=0; comp<data_at_vpt.nComp(); comp++) {
        for (int dir=0; dir<SpaceDim; dir++) {
          this_data.plus(this_tmp_cfg,comp*SpaceDim + dir,comp,1);
        }
      }
      this_data.mult(1./SpaceDim);
   }

   CFG::Box domain_box = configurationDomainBox(grids);
   domain_box.grow(4);
   const CFG::MultiBlockCoordSys* mag_coord_sys = m_mag_geom.coordSysPtr();
   WriteMappedUGHDF5(a_file_name, grids, data_at_vpt, *mag_coord_sys, domain_box, a_time);
}



void
PhaseGeom::plotAtVelocityIndex( const string                 a_file_name,
                                const VEL::IntVect&          a_vspace_index,
                                const LevelData<FArrayBox>&  a_data,
                                const double&                a_time ) const
{
   CH_assert( validVelocityIndex(a_vspace_index) );

   // Evaluate data at the specified velocity space index.
   CFG::LevelData<CFG::FArrayBox> data_at_vpt;
   getConfigurationData(a_vspace_index, a_data, data_at_vpt);

   const CFG::DisjointBoxLayout& grids = data_at_vpt.disjointBoxLayout();

   CFG::Box domain_box = configurationDomainBox(grids);
   domain_box.grow(4);
   const CFG::MultiBlockCoordSys* mag_coord_sys = m_mag_geom.coordSysPtr();
   WriteMappedUGHDF5(a_file_name, grids, data_at_vpt, *mag_coord_sys, domain_box, a_time);
}



void
PhaseGeom::plotAtConfigurationIndex( const string               a_file_name,
                                     const CFG::IntVect&        cspace_index,
                                     const LevelData<FluxBox>&  a_data,
                                     const double&              a_time ) const
{
   CH_assert( validConfigurationIndex(cspace_index) );

   // Since the passed data is a LevelData<FluxBox>, we average the face
   // values to cell centers
   LevelData<FArrayBox> tmp(a_data.getBoxes(), SpaceDim*a_data.nComp(), a_data.ghostVect());
   EdgeToCell(a_data, tmp);

   // Evaluate data at the specified configuration space index.  We would have done
  // this before averaging faces to cell centers, but we don't have an
   // EdgeToCell function that works only in velocity space.
   VEL::LevelData<VEL::FArrayBox> tmp_vel;
   getVelocityData(cspace_index, tmp, tmp_vel);

   // Next, average the directions
   VEL::LevelData<VEL::FArrayBox> data_at_cpt(tmp_vel.getBoxes(), a_data.nComp(), tmp_vel.ghostVect());

   const VEL::DisjointBoxLayout& grids = data_at_cpt.disjointBoxLayout();

   VEL::DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      VEL::FArrayBox& this_tmp_vel = tmp_vel[dit];
      VEL::FArrayBox& this_data = data_at_cpt[dit];

      this_data.setVal(0.);
      for (int comp=0; comp<data_at_cpt.nComp(); comp++) {
         for (int dir=CFG_DIM; dir<SpaceDim; dir++) {
            this_data.plus(this_tmp_vel,comp*SpaceDim + dir,comp,1);
         }
      }
      this_data.mult(1./VEL_DIM);
   }

   VEL::Box domain_box = velocityDomainBox(grids);
   WriteMappedUGHDF5(a_file_name, grids, data_at_cpt, m_vel_coords, domain_box, a_time);
}



void
PhaseGeom::plotAtConfigurationIndex( const string                 a_file_name,
                                     const CFG::IntVect&          cspace_index,
                                     const LevelData<FArrayBox>&  a_data,
                                     const double&                a_time ) const
{
   CH_assert( validConfigurationIndex(cspace_index) );

   // Evaluate data at the specified configuration space index.
   VEL::LevelData<VEL::FArrayBox> data_at_cpt;
   getVelocityData(cspace_index, a_data, data_at_cpt);

   const VEL::DisjointBoxLayout& grids = data_at_cpt.disjointBoxLayout();

   VEL::Box domain_box = velocityDomainBox(grids);
   WriteMappedUGHDF5(a_file_name, grids, data_at_cpt, m_vel_coords, domain_box, a_time);
}



void
PhaseGeom::plotConfigurationData( const string                           a_file_name,
                                  const CFG::LevelData<CFG::FArrayBox>&  a_data,
                                  const double&                          a_time ) const
{
   const CFG::DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   CFG::Box domain_box = configurationDomainBox(grids);
   const CFG::MultiBlockCoordSys* mag_coord_sys = m_mag_geom.coordSysPtr();
   WriteMappedUGHDF5(a_file_name, grids, a_data, *mag_coord_sys, domain_box, a_time);
}



void
PhaseGeom::plotAtMuIndex( const string                 a_file_name,
                          const int                    a_mu_index,
                          const LevelData<FArrayBox>&  a_data,
                          const double&                a_time ) const
{
   CH_assert( m_domain.domainBox().smallEnd(MU_DIR) <= a_mu_index && a_mu_index <= m_domain.domainBox().bigEnd(MU_DIR) );

   SliceSpec slice_mu(MU_DIR, a_mu_index);
   CP1::LevelData<CP1::FArrayBox> data_at_mu_index;
   sliceLevelData(data_at_mu_index, a_data, slice_mu);

   const CP1::DisjointBoxLayout& grids = data_at_mu_index.disjointBoxLayout();

   CP1::IntVectSet gridIVS;
   CP1::LayoutIterator lit = grids.layoutIterator();
   for (lit.begin(); lit.ok(); ++lit) {
      gridIVS |= grids[lit];
   }
   CP1::Box domain_box = gridIVS.minBox();

   Vector<CP1::DisjointBoxLayout> gridsv;
   gridsv.push_back(grids);

   Vector<CP1::LevelData<CP1::FArrayBox>* > datav;
   datav.push_back(&data_at_mu_index);

   Vector<string> namev;
   namev.push_back("component_0");

   const PhaseBlockCoordSys& block_zero_coord_sys = (const PhaseBlockCoordSys&)(*m_coordSysPtr->getCoordSys(0));

   double dvp = block_zero_coord_sys.dx()[VPARALLEL_DIR];
   double dt = 1.;

   Vector<int> ratiov;
   ratiov.push_back(1);
   int numlevels = 1;

   string full_file_name = a_file_name + "3d.hdf5";

   CP1::WriteAMRHierarchyHDF5(full_file_name,
                              gridsv,
                              datav,
                              namev,
                              domain_box,
                              dvp,
                              dt,
                              a_time,
                              ratiov,
                              numlevels);
}



void
PhaseGeom::plotConfigurationData( const string                         a_file_name,
                                  const CFG::LevelData<CFG::FluxBox>&  a_data,
                                  const double&                        a_time ) const
{
   const CFG::DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   CFG::LevelData<CFG::FArrayBox> data_cell(grids, a_data.nComp(), a_data.ghostVect());

   CFG::DataIterator cdit = grids.dataIterator();
   for (cdit.begin(); cdit.ok(); ++cdit) {
      const CFG::FluxBox& this_data = a_data[cdit];
      CFG::FArrayBox& this_data_cell = data_cell[cdit];
      this_data_cell.setVal(0.);
      for (int dir=0; dir<CFG_DIM; dir++) {
         const CFG::FArrayBox& this_data_dir = this_data[dir];
         CFG::FArrayBox shift_fab(this_data_dir.box(), this_data_dir.nComp());
         shift_fab.copy(this_data_dir);

         shift_fab.shiftHalf(dir,-1);
         this_data_cell.plus(shift_fab);
         shift_fab.shiftHalf(dir,2);
         this_data_cell.plus(shift_fab);
      }
      this_data_cell.mult(0.5/CFG_DIM);
   }

   CFG::Box domain_box = configurationDomainBox(grids);
   const CFG::MultiBlockCoordSys* mag_coord_sys = m_mag_geom.coordSysPtr();
   WriteMappedUGHDF5(a_file_name, grids, data_cell, *mag_coord_sys, domain_box, a_time);
}


void
PhaseGeom::plotVelocityData( const string                           a_file_name,
                             const VEL::LevelData<VEL::FArrayBox>&  a_data,
                             const double&                          a_time ) const
{
   const VEL::DisjointBoxLayout& grids = a_data.disjointBoxLayout();

   VEL::Box domain_box = velocityDomainBox(grids);
   WriteMappedUGHDF5(a_file_name, grids, a_data, m_vel_coords, domain_box, a_time);
}



void
PhaseGeom::plotVParPoloidalData( const string                a_file_name,
                                 const int                   a_radial_index,
                                 const int                   a_toroidal_index,
                                 const int                   a_mu_index,
                                 const LevelData<FArrayBox>& a_data,
                                 const double&               a_time ) const
{
  // Slice in the mu direction at the specified mu coordinate

  SliceSpec slice_mu(MU_DIR, a_mu_index);
  CP1::LevelData<CP1::FArrayBox> mu_sliced_data;
  sliceLevelData(mu_sliced_data, a_data, slice_mu);

  CP1::Box mu_sliced_domainbox;
  sliceBox(mu_sliced_domainbox, m_domain.domainBox(), slice_mu);

#if CFG_DIM==3

  // Slice in the toroidal direction at the specified toroidal coordinate

  CP1::SliceSpec slice_z(z_DIR, a_toroidal_index);
  D3::LevelData<D3::FArrayBox> muz_sliced_data;
  sliceLevelData(muz_sliced_data, mu_sliced_data, slice_z);

  D3::Box muz_sliced_domainbox;
  sliceBox(muz_sliced_domainbox, mu_sliced_domainbox, slice_z);

#endif

  // Slice in the radial direction at the specified radial coordinate

  D3::SliceSpec slice_r(RADIAL_DIR, a_radial_index);
  D2::LevelData<D2::FArrayBox> sliced_data;
  sliceLevelData(sliced_data, mu_sliced_data, slice_r);

  D2::Box sliced_domainbox;
  sliceBox(sliced_domainbox, mu_sliced_domainbox, slice_r);

  const Vector<D2::DisjointBoxLayout> vectGrids(1, sliced_data.getBoxes());
  const Vector<D2::LevelData<D2::FArrayBox>* > vectData(1, const_cast<D2::LevelData<D2::FArrayBox>* >(&sliced_data) );
  Vector<int> refRatio(1,1);
  int numLevels = 1;

  // create names for the variables and placeholder values for dt and time,
  Real dt = 1.0;
  Real time = 1.0;
  Real dx = 1.0;
  int nComp = sliced_data.nComp();
  Vector<string> compNames(nComp);
  for (int n=0; n<nComp; n++)
     {
        char labelChSt[80];
        sprintf(labelChSt, "component_%d", n);
        string label(labelChSt);
        compNames[n] = label;
     }

  WriteAMRHierarchyHDF5(a_file_name,
                        vectGrids,
                        vectData,
                        compNames,
                        sliced_domainbox,
                        dx,
                        dt,
                        time,
                        refRatio,
                        numLevels);
}



void
PhaseGeom::plotDivergence( const string&             a_file_name,
                           const LevelData<FluxBox>& a_flux,
                           const double              a_time ) const
{
   CH_assert(a_flux.ghostVect() >= IntVect::Unit);
   const DisjointBoxLayout& grids = a_flux.disjointBoxLayout();

   LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      flux[dit].copy(a_flux[dit]);
   }

   applyAxisymmetricCorrection( flux );

   if (m_second_order) fourthOrderAverage(flux);

   LevelData<FArrayBox> divergence(grids, 1, IntVect::Zero);
   bool OMIT_NT = false;
   mappedGridDivergence( divergence, flux, OMIT_NT );

   LevelData<FArrayBox> cell_volume(grids, 1, IntVect::Zero);
   getCellVolumes(cell_volume);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      divergence[dit] /= cell_volume[dit];
   }

   plotData( a_file_name, divergence, a_time );
}



void
PhaseGeom::plotData( const string&                a_file_name,
                     const LevelData<FArrayBox>&  a_data,
                     const double&                a_time ) const
{
   WriteMappedUGHDF5(a_file_name, m_gridsFull, a_data, m_phase_coords, m_domain.domainBox(), a_time);
}



CFG::Box
PhaseGeom::configurationDomainBox(const CFG::DisjointBoxLayout& grids) const
{
   // bit of a hack to get a reasonable domainBox
   // this is silly, but quick and easy
   CFG::IntVectSet gridIVS;
   CFG::LayoutIterator lit = grids.layoutIterator();
   for (lit.begin(); lit.ok(); ++lit) {
      gridIVS |= grids[lit];
   }
   return gridIVS.minBox();
}



VEL::Box
PhaseGeom::velocityDomainBox(const VEL::DisjointBoxLayout& grids) const
{
   // bit of a hack to get a reasonable domainBox
   // this is silly, but quick and easy
   VEL::IntVectSet gridIVS;
   VEL::LayoutIterator lit = grids.layoutIterator();
   for (lit.begin(); lit.ok(); ++lit) {
      gridIVS |= grids[lit];
   }
   return gridIVS.minBox();
}



bool
PhaseGeom::validConfigurationIndex(const CFG::IntVect& index) const
{
   const CFG::MultiBlockCoordSys* mag_coord_sys = m_mag_geom.coordSysPtr();
   return mag_coord_sys->whichBlock(index) != -1;
}



bool
PhaseGeom::validVelocityIndex(const VEL::IntVect& index) const
{
   return m_vel_coords.domain().domainBox().contains(index);
}



Box
PhaseGeom::configurationFlatten( const Box& a_disjointBox,
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
PhaseGeom::velocityFlatten( const Box& a_disjointBox,
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


const PhaseBlockCoordSys&
PhaseGeom::getBlockCoordSys(const Box& a_box) const
{
  int block_number = m_coordSysPtr->whichBlock(a_box);
  return (const PhaseBlockCoordSys&)(*m_coordSysPtr->getCoordSys(block_number));
}



const CFG::MagBlockCoordSys&
PhaseGeom::getMagBlockCoordSys(const Box& a_box) const
{
  int block_number = m_coordSysPtr->whichBlock(a_box);
  const CFG::MagCoordSys* mag_coord_sys = (const CFG::MagCoordSys*)m_mag_geom.coordSysPtr();
  return *mag_coord_sys->getCoordSys(block_number);
}



#include "NamespaceFooter.H"
