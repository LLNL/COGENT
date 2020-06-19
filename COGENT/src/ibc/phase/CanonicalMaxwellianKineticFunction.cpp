#include "CanonicalMaxwellianKineticFunction.H"

#include <iostream>
#include <typeinfo>
#include <array>

#include "DataIterator.H"
#include "Directions.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "FourthOrderUtil.H"
#include "LevelData.H"
#include "MaxwellianKineticFunctionF_F.H"
#include "MayDay.H"
#include "LogRectPhaseCoordSys.H"
#include "SingleNullPhaseCoordSys.H"
#include "SNCorePhaseCoordSys.H"
#include "MultiBlockCoordSys.H"
#include "PhaseBlockCoordSys.H"
#include "Vector.H"

#include "KineticFunctionUtils.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "ToroidalBlockCoordSys.H"
#include "GridFunctionLibrary.H"
#include "Constant.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"



CanonicalMaxwellianKineticFunction::CanonicalMaxwellianKineticFunction( ParmParse& a_pp,
                                                      const int& a_verbosity )
   : m_verbosity(a_verbosity),
     m_profile_option(1),
     m_zero_larmor(false)

{
   parseParameters( a_pp );
}


template<class Func> bool zeroCheck( const LevelData<FArrayBox>& a_data,
                                     Func a_func )
{
   Real min( BASEFAB_REAL_SETVAL );
   const DisjointBoxLayout& grids( a_data.disjointBoxLayout() );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      min = std::min( a_data[dit].min( grids[dit] ), min );
   }
   return a_func( min, 0.0 );
}

inline
bool isPositiveDefinite( const LevelData<FArrayBox>& a_data )
{
   return zeroCheck( a_data, std::greater<Real>() );
}


inline
bool isNonNegative( const LevelData<FArrayBox>& a_data )
{
   return zeroCheck( a_data, std::greater_equal<Real>() );
}


void CanonicalMaxwellianKineticFunction::assign( KineticSpecies& a_species,
                                                 const Real& a_time ) const
{
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   checkGeometryValidity( geometry );
   
   LevelData<FArrayBox>& dfn( a_species.distributionFunction() );
   const DisjointBoxLayout& grids( dfn.disjointBoxLayout() );

   const LevelData<FArrayBox>& injected_B( geometry.getBFieldMagnitude() );
   const LevelData<FArrayBox>& real_coords( geometry.getCellCenteredRealCoords() );

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      setPointValues( dfn[dit],
                      dfn[dit].box(),
                      geometry.getBlockCoordSys( grids[dit] ),
                      geometry.getMagBlockCoordSys( grids[dit] ),
                      injected_B[dit],
                      real_coords[dit],
                      a_species.mass(),
                      a_species.charge(),
                      geometry.getLarmorNumber());
   }

   geometry.multBStarParallel( dfn );
   if ( !(geometry.secondOrder()) )  {
      KineticFunctionUtils::convertToCellAverage( geometry, dfn );
   }
   geometry.multJonValid( dfn );
   dfn.exchange();
}


void CanonicalMaxwellianKineticFunction::assign( KineticSpecies& a_species,
                                                 const BoundaryBoxLayout& a_bdry_layout,
                                                 const Real& a_time ) const
{
   CH_TIMERS("CanonicalMaxwellianKineticFunction::assign");
   CH_TIMER("setPointValues", t_set_point_values);
   CH_TIMER("multBStarParallel", t_mult_bstar_parallel);

   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   checkGeometryValidity( geometry );

   LevelData<FArrayBox>& dfn( a_species.distributionFunction() );
   const DisjointBoxLayout& grids( dfn.disjointBoxLayout() );

   const LevelData<FArrayBox>& injected_B( geometry.getBFieldMagnitude() );
   const LevelData<FArrayBox>& real_coords( geometry.getCellCenteredRealCoords() );

   // NB: This is a cheat - there's one too many cells at the (dir,side) face
   // of the boundary box, but it doesn't matter because one-sided difference
   // will be used at that face to construct the cell average.  We do want the
   // extra cell in all other directions.
   LevelData<FArrayBox> dfn_tmp( grids, dfn.nComp(), IntVect::Unit );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {

      const Box& interior_box( a_bdry_layout.interiorBox( dit ) );
      Box fill_box( dfn_tmp[dit].box() );
      fill_box.growDir( a_bdry_layout.dir(), a_bdry_layout.side(), -1 );

      const PhaseBlockCoordSys& phase_coord_sys( geometry.getBlockCoordSys( interior_box ) );
      const CFG::MagBlockCoordSys& mag_coord_sys( geometry.getMagBlockCoordSys( interior_box ) );
      
      const DataIndex& internal_dit( a_bdry_layout.dataIndex( dit ) );
      CH_START(t_set_point_values);
      setPointValues( dfn_tmp[dit],
                      fill_box,
                      phase_coord_sys,
                      mag_coord_sys,
                      injected_B[internal_dit],
                      real_coords[internal_dit],
                      a_species.mass(),
                      a_species.charge(),
                      geometry.getLarmorNumber());
      CH_STOP(t_set_point_values);
   }

   CH_START(t_mult_bstar_parallel);
   geometry.multBStarParallel( dfn_tmp, a_bdry_layout );
   CH_STOP(t_mult_bstar_parallel);
   
   if ( !(geometry.secondOrder()) )  {
      for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
         Box domain_box( dfn_tmp[dit].box() );
         domain_box.growDir( a_bdry_layout.dir(), a_bdry_layout.side(), -1 );
         ProblemDomain domain( domain_box );
         fourthOrderAverageCell( dfn_tmp[dit], domain, grids[dit] );
      }
   }
   
   dfn_tmp.copyTo( dfn );
   dfn.exchange();
}


inline
void CanonicalMaxwellianKineticFunction::setPointValues(FArrayBox&                    a_dfn,
                                                        const Box&                    a_box,
                                                        const PhaseBlockCoordSys&     a_phase_coord_sys,
                                                        const CFG::MagBlockCoordSys&  a_mag_coord_sys,
                                                        const FArrayBox&              a_B,
                                                        const FArrayBox&              a_cc_coords,
                                                        const Real&                   a_mass,
                                                        const Real&                   a_charge,
                                                        const Real&                   a_larmor_number ) const
{
   CH_TIMERS("CanonicalMaxwellianKineticFunction::setPointValues");
   CH_TIMER("getCellCenteredRealCoords", t_get_real_coords);

#if CFG_DIM == 3
   if (typeid(a_mag_coord_sys) != typeid(CFG::ToroidalBlockCoordSys) &&
       typeid(a_mag_coord_sys) != typeid(CFG::SingleNullBlockCoordSys)) {
      const std::string msg( "CanonicalMaxwellianKineticFunction: only toroidal geometry implementation exists ");
      MayDay::Error( msg.c_str() );
   }
#endif
   
   IntVect lo = a_box.smallEnd();
   IntVect hi = a_box.bigEnd();
   
   CFG::Box cfg_box(a_phase_coord_sys.config_restrict(lo),
                    a_phase_coord_sys.config_restrict(hi));
   
   CFG::FArrayBox cc_cfg_coords( cfg_box, CFG_DIM );
   a_mag_coord_sys.getCellCenteredRealCoords( cc_cfg_coords );
   
   CFG::FArrayBox magnetic_flux( cfg_box, 1 );
   a_mag_coord_sys.getMagneticFlux( cc_cfg_coords, magnetic_flux );
   
   FArrayBox magnetic_flux_inj;
   a_phase_coord_sys.injectConfigurationToPhase(magnetic_flux, magnetic_flux_inj);

   CH_START(t_get_real_coords);
   FArrayBox cc_phase_coords( a_box, SpaceDim );
   CH_assert((a_cc_coords.box()).contains(a_box));
   cc_phase_coords.copy(a_cc_coords);
   CH_STOP(t_get_real_coords);
   
   Real R0;          //Major radius coordinate of the magnetic axis
   Real RBtor;       //R*Btor constant factor
   Real psi_mid;     //Magnetic flux (unnormalized) at the midpoint
   Real dpsidr_mid;  //dpsi/dr at the midpoint
   
   FArrayBox toroidal_coords_inj;

#if CFG_DIM == 3
   if (typeid(a_mag_coord_sys) == typeid(CFG::ToroidalBlockCoordSys)) {

      CFG::FArrayBox cfg_toroidal_coords( cfg_box, CFG_DIM );
      ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).getToroidalCoords( cfg_toroidal_coords );
      
      a_phase_coord_sys.injectConfigurationToPhase(cfg_toroidal_coords, toroidal_coords_inj);

      if (m_r_mid < 0) {
         m_r_mid = ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).getAvMinorRad();
      }
      
      RBtor = ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).getRBtoroidal();
      R0 = ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).centralMajorRadius();
      psi_mid = ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).getMagneticFlux(m_r_mid);
      dpsidr_mid = ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).dpsidr(m_r_mid);
      
      //Presently not used; but might be useful for some other models
      Real q_mid = ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).getSafetyFactor(m_r_mid/R0);
   }
#endif
   
   if (typeid(a_mag_coord_sys) == typeid(CFG::SingleNullBlockCoordSys)) {

      CFG::FArrayBox cfg_toroidal_coords( cfg_box, CFG_DIM );
      ((const CFG::SingleNullBlockCoordSys&)a_mag_coord_sys).getToroidalCoords( cfg_toroidal_coords );
      
      a_phase_coord_sys.injectConfigurationToPhase(cfg_toroidal_coords, toroidal_coords_inj);
      
      if (m_r_mid < 0) {
         const std::string msg( "CanonicalMaxwellian:: midpoint_minor_radius must be specified for SN geometry ");
         MayDay::Error( msg.c_str() );
      }
      
      RBtor = ((const CFG::SingleNullBlockCoordSys&)a_mag_coord_sys).getRBtoroidal();
      CFG::RealVect axis = ((const CFG::SingleNullBlockCoordSys&)a_mag_coord_sys).getMagAxis();
      R0 = axis[RADIAL_DIR];
      
#if CFG_DIM == 3
      CFG::RealVect X_mid(axis[RADIAL_DIR]+m_r_mid,0,axis[POLOIDAL_DIR]);
#else
      CFG::RealVect X_mid(axis[RADIAL_DIR]+m_r_mid, axis[POLOIDAL_DIR]);
#endif
      array<double,3> B_mid = ((const CFG::SingleNullBlockCoordSys&)a_mag_coord_sys).computeBField(X_mid);
#if CFG_DIM == 3
      Real Bpol = sqrt(pow(B_mid[RADIAL_DIR],2)+pow(B_mid[POLOIDAL_DIR],2));
#else
      Real Bpol = sqrt(pow(B_mid[RADIAL_DIR],2)+pow(B_mid[2],2));
#endif
      psi_mid = ((const CFG::SingleNullBlockCoordSys&)a_mag_coord_sys).getMagneticFlux(X_mid);
      dpsidr_mid = X_mid[RADIAL_DIR] * Bpol;
   }
   
   Real larmor_number = (m_zero_larmor) ? 0.0 : a_larmor_number;

   //NB:For convenience the input for kappa is always in 1/R0 units,
   //thus need to renormalize to be consitent the paticular coordinate choice
   //which is done inside Fortran
   
   FORT_SET_CANONICAL_MAXWELL(CHF_FRA1(a_dfn,0),
                              CHF_BOX(a_box),
                              CHF_CONST_FRA(cc_phase_coords),
                              CHF_CONST_FRA(toroidal_coords_inj),
                              CHF_CONST_FRA1(a_B,0),
                              CHF_CONST_FRA1(magnetic_flux_inj,0),
                              CHF_CONST_REAL(RBtor),
                              CHF_CONST_REAL(R0),
                              CHF_CONST_REAL(dpsidr_mid),
                              CHF_CONST_REAL(psi_mid),
                              CHF_CONST_REAL(m_r_mid),
                              CHF_CONST_VR(m_density_parm),
                              CHF_CONST_VR(m_temperature_parm),
                              CHF_CONST_VR(m_perturbation_parm),
                              CHF_CONST_REAL(a_mass),
                              CHF_CONST_REAL(a_charge),
                              CHF_CONST_REAL(larmor_number),
                              CHF_CONST_INT(m_profile_option));
}

inline
void CanonicalMaxwellianKineticFunction::checkGeometryValidity( const PhaseGeom& a_geometry ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   bool unknown_geom( typeid(coord_sys) != typeid(LogRectPhaseCoordSys) );
   unknown_geom &= (typeid(coord_sys) != typeid(SingleNullPhaseCoordSys));
   unknown_geom &= (typeid(coord_sys) != typeid(SNCorePhaseCoordSys));
   
   if ( unknown_geom ) {
      const std::string msg( "attempt to use unknown geometry. ");
      MayDay::Error( msg.c_str() );
   }
}

inline
void CanonicalMaxwellianKineticFunction::parseParameters( ParmParse& a_pp )
{

   a_pp.query( "functional_option", m_profile_option);
   a_pp.query( "zero_larmor", m_zero_larmor); 
  
   if (a_pp.contains("midpoint_minor_radius")) {
      a_pp.get("midpoint_minor_radius", m_r_mid);
   }
   else {
      m_r_mid = -1.0;
   }
      
   
   int num_parm;
   if (m_profile_option == 1) {
      num_parm = 4;
   }
   else if (m_profile_option == 2) {
      num_parm = 3;
   }
   else {
      const std::string msg( "only profile option 1 or 2 are supported ");
      MayDay::Error( msg.c_str() );
   }

   m_temperature_parm.resize(num_parm);
   m_density_parm.resize(num_parm);

   a_pp.get( "density_midpoint_value", m_density_parm[0] );
   a_pp.get( "density_kappa",          m_density_parm[1] );
   a_pp.get( "density_radial_width",   m_density_parm[2] );
   
   CH_assert( m_density_parm[2]>=0 );
   
   a_pp.get( "temperature_midpoint_value", m_temperature_parm[0] );
   a_pp.get( "temperature_kappa",          m_temperature_parm[1] );
   a_pp.get( "temperature_radial_width",   m_temperature_parm[2] );
   
   CH_assert( m_temperature_parm[2]>=0 );

   if (m_profile_option == 1) {
      a_pp.get( "density_flattop_fac",     m_density_parm[3] );
      a_pp.get( "temperature_flattop_fac", m_temperature_parm[3] );
   }
   
   m_perturbation_parm.resize(4);
   m_perturbation_parm.assign( 0.0 );
   a_pp.getarr( "perturbation", m_perturbation_parm, 0, 4 );
   
   CH_assert( m_perturbation_parm[1]>=0 );
   
   if (m_verbosity) {
      printParameters();
   }
}

inline
void CanonicalMaxwellianKineticFunction::printParameters() const
{
   if (procID()==0) {
      std::cout << "Canonical Maxwellian kinetic function parameters:" << std::endl;
      std::cout << "  density_midpoint_value: "       << m_density_parm[0]       << std::endl;
      std::cout << "  density_kappa: "                << m_density_parm[1]       << std::endl;
      std::cout << "  density_radial_width: "         << m_density_parm[2]       << std::endl;
      if (m_profile_option == 1) {
         std::cout << "  density_flattop_fac: "       << m_density_parm[3]       << std::endl;
      }
      std::cout << "  temperature_midpoint_value: "   << m_temperature_parm[0]   << std::endl;
      std::cout << "  temperature_kappa: "            << m_temperature_parm[1]   << std::endl;
      std::cout << "  temperature_radial_width: "     << m_temperature_parm[2]   << std::endl;
      if (m_profile_option == 1) {
         std::cout << "  temeprature_flattop_fac: "   << m_temperature_parm[3]   << std::endl;
      }
     
      std::cout << "  perturbation_parameters: "      << m_perturbation_parm     << std::endl;
      std::cout << std::endl;
   }
}

#include "NamespaceFooter.H"
