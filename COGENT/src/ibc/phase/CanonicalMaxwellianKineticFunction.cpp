#include "CanonicalMaxwellianKineticFunction.H"

#include <iostream>
#include <typeinfo>

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
     m_density_val(1.0),
     m_density_kappa(0.0),
     m_density_width(0.0),
     m_temperature_val(1.0),
     m_temperature_kappa(0.0),
     m_temperature_width(0.0)
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

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      setPointValues( dfn[dit],
                      dfn[dit].box(),
                      geometry.getBlockCoordSys( grids[dit] ),
                      geometry.getMagBlockCoordSys( grids[dit] ),
                      injected_B[dit],
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
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   checkGeometryValidity( geometry );

   LevelData<FArrayBox>& dfn( a_species.distributionFunction() );
   const DisjointBoxLayout& grids( dfn.disjointBoxLayout() );

   const LevelData<FArrayBox>& injected_B( geometry.getBFieldMagnitude() );

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
      setPointValues( dfn_tmp[dit],
                      fill_box,
                      phase_coord_sys,
                      mag_coord_sys,
                      injected_B[internal_dit],
                      a_species.mass(),
                      a_species.charge(),
                      geometry.getLarmorNumber());
   }

   geometry.multBStarParallel( dfn_tmp, a_bdry_layout );
   
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
                                                        const Real&                   a_mass,
                                                        const Real&                   a_charge,
                                                        const Real&                   a_larmor_number ) const
{

#if CFG_DIM == 3
   if (typeid(a_mag_coord_sys) != typeid(CFG::ToroidalBlockCoordSys)) {
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

   FArrayBox cc_phase_coords( a_box, SpaceDim );
   a_phase_coord_sys.getCellCenteredRealCoords( cc_phase_coords );

#if CFG_DIM == 3
   
   CFG::FArrayBox cfg_toroidal_coords( cfg_box, CFG_DIM );
   ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).getToroidalCoords( cfg_toroidal_coords );

   FArrayBox toroidal_coords_inj;
   a_phase_coord_sys.injectConfigurationToPhase(cfg_toroidal_coords, toroidal_coords_inj);

   Real RBtor = ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).getRBtoroidal();
   Real R0 = ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).centralMajorRadius();
   Real r_p = ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).getAvMinorRad();
   Real psi_p = ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).getMagneticFlux(r_p);
   Real dpsidr_p = ((const CFG::ToroidalBlockCoordSys&)a_mag_coord_sys).dpsidr(r_p);

   //For convenience the input for kappa is in 1/R0 units, 
   //thus need to renormalize to be consitent with 1/psi units
   Real n_kappa_psi = m_density_kappa/(dpsidr_p * R0);
   Real T_kappa_psi = m_temperature_kappa/(dpsidr_p * R0);
   
   FORT_SET_CANONICAL_MAXWELL(CHF_FRA(a_dfn),
                              CHF_BOX(a_box),
                              CHF_CONST_FRA(cc_phase_coords),
			      CHF_CONST_FRA(toroidal_coords_inj),
                              CHF_CONST_FRA1(a_B,0),
                              CHF_CONST_FRA1(magnetic_flux_inj,0),
                              CHF_CONST_REAL(RBtor),
                              CHF_CONST_REAL(psi_p),
                              CHF_CONST_REAL(dpsidr_p),
                              CHF_CONST_REAL(m_density_val),
                              CHF_CONST_REAL(n_kappa_psi),
                              CHF_CONST_REAL(m_density_width),
                              CHF_CONST_REAL(m_temperature_val),
                              CHF_CONST_REAL(T_kappa_psi),
                              CHF_CONST_REAL(m_temperature_width),
      			      CHF_CONST_REALVECT(m_mode),
                              CHF_CONST_REAL(a_mass),
                              CHF_CONST_REAL(a_charge),
                              CHF_CONST_REAL(a_larmor_number));
#endif
}

inline
void CanonicalMaxwellianKineticFunction::checkGeometryValidity( const PhaseGeom& a_geometry ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   bool unknown_geom( typeid(coord_sys) != typeid(LogRectPhaseCoordSys) );
   unknown_geom &= (typeid(coord_sys) != typeid(SingleNullPhaseCoordSys));
   unknown_geom &= (typeid(coord_sys) != typeid(SNCorePhaseCoordSys));
   
   if ( unknown_geom ) {
      const std::string msg( "CanonicalMaxwellianKineticFunction: Attempt to use unknown geometry. ");
      MayDay::Error( msg.c_str() );
   }
}

inline
void CanonicalMaxwellianKineticFunction::parseParameters( ParmParse& a_pp )
{
   a_pp.get( "density_midpoint_value", m_density_val );
   a_pp.get( "density_radial_width",   m_density_width );
   a_pp.get( "density_kappa",          m_density_kappa );

   CH_assert( m_density_width>=0 );
   
   a_pp.get( "temperature_midpoint_value", m_temperature_val );
   a_pp.get( "temperature_radial_width",   m_temperature_width );
   a_pp.get( "temperature_kappa",          m_temperature_kappa );

   CH_assert( m_temperature_width>=0 );
   
   Vector<Real> temp( CFG_DIM );
   temp.assign( 0.0 );
   a_pp.getarr( "mode", temp, 0, CFG_DIM );
   m_mode = CFG::RealVect( temp );

   if (m_verbosity) {
      printParameters();
   }
}

inline
void CanonicalMaxwellianKineticFunction::printParameters() const
{
   if (procID()==0) {
      std::cout << "Canonical Maxwellian kinetic function parameters:" << std::endl;
      std::cout << "  density_midpoint_value: "       << m_density_val           << std::endl;
      std::cout << "  density_radial_width: "         << m_density_width         << std::endl;
      std::cout << "  density_kappa: "                << m_density_kappa         << std::endl;
      std::cout << "  temperature_midpoint_value: "   << m_temperature_val       << std::endl;
      std::cout << "  temperature_radial_width: "     << m_temperature_width     << std::endl;
      std::cout << "  temperature_kappa: "            << m_temperature_kappa     << std::endl;
      std::cout << "  mode_coefficients: "            << m_mode                  << std::endl;
      std::cout << std::endl;
   }
}

#include "NamespaceFooter.H"
