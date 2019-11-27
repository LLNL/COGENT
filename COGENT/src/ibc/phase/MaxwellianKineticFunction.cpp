#include "MaxwellianKineticFunction.H"

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
#include "GridFunctionLibrary.H"
#include "Constant.H"
#include "inspect.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#include "Slicing.H.transdim"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif
#define CH_SPACEDIM PDIM


#include "NamespaceHeader.H"

using namespace CH_MultiDim;


MaxwellianKineticFunction::MaxwellianKineticFunction( ParmParse& a_pp,
                                                      const int& a_verbosity )
   : m_verbosity(a_verbosity),
     m_enforce_input_density_profile(false)
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


void MaxwellianKineticFunction::assign( KineticSpecies& a_species,
                                        const Real& a_time ) const
{
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   checkGeometryValidity( geometry );

   LevelData<FArrayBox>& dfn( a_species.distributionFunction() );
   const DisjointBoxLayout& grids( dfn.disjointBoxLayout() );

   CFG::IntVect nghosts( CFG::IntVect::Zero );
   for (int d(0); d<CFG_DIM; d++) {
      nghosts[d] = (dfn.ghostVect())[d];
   }
   
   LevelData<FArrayBox> injected_density;
   initializeField( injected_density, *m_ic_density, geometry, nghosts, a_time );
   CH_assert( isNonNegative( injected_density ) );

   LevelData<FArrayBox> injected_temperature;
   initializeField( injected_temperature, *m_ic_temperature, geometry, nghosts, a_time );
   CH_assert( isPositiveDefinite( injected_temperature ) );

   LevelData<FArrayBox> injected_vparallel;
   initializeField( injected_vparallel, *m_ic_vparallel, geometry, nghosts, a_time );
   //CH_assert( isPositiveDefinite( injected_vparallel ) );
 
   const LevelData<FArrayBox>& injected_B( geometry.getBFieldMagnitude() );
   const LevelData<FArrayBox>& real_coords( geometry.getCellCenteredRealCoords() );

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      setPointValues( dfn[dit],
                      dfn[dit].box(),
                      geometry.getBlockCoordSys( grids[dit] ),
                      injected_density[dit],
                      injected_temperature[dit],
                      injected_vparallel[dit],
                      injected_B[dit],
                      real_coords[dit],
                      a_species.mass() );
   }

   geometry.multBStarParallel( dfn );
   if ( !(geometry.secondOrder()) )  {
      KineticFunctionUtils::convertToCellAverage( geometry, dfn );
   }
   
   // This slightly changes dfn to exactly match
   // the density profile in the input; this works only if mult and
   // divide J on valid are exact inversion of each other (which
   // may not be the case for a high-order). Presently only applied in
   // the valid cells, perhaps extend to ghost region later.
   if (m_enforce_input_density_profile)
   {
      const CFG::MultiBlockLevelGeom& mag_geometry( geometry.magGeom() );
      const CFG::DisjointBoxLayout& mag_grids( mag_geometry.grids() );

      CFG::LevelData<CFG::FArrayBox> density_moment(mag_grids,1,CFG::IntVect::Zero);
      a_species.numberDensity(density_moment);

      LevelData<FArrayBox> density_moment_inj;
      geometry.injectConfigurationToPhase( density_moment, density_moment_inj );

      for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
         FORT_ENFORCE_INPUT_DENS_PROF(CHF_FRA(dfn[dit]),
				      CHF_BOX(grids[dit]),
                                      CHF_CONST_FRA1(density_moment_inj[dit],0),
                                      CHF_CONST_FRA1(injected_density[dit],0));
      }
   }
   
   geometry.multJonValid( dfn );
   dfn.exchange();
}


void MaxwellianKineticFunction::assign( KineticSpecies& a_species,
                                        const BoundaryBoxLayout& a_bdry_layout,
                                        const Real& a_time ) const
{
   CH_TIMERS("MaxwellianKineticFunction::assign");
   CH_TIMER("setPointValues", t_set_point_values);
   CH_TIMER("initializeDensity", t_initialize_density);
   CH_TIMER("initializeTemperature", t_initialize_temperature);
   CH_TIMER("initializeVparallel", t_initialize_vparallel);
   CH_TIMER("multBStarParallel", t_mult_bstar_parallel);

   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   checkGeometryValidity( geometry );

   LevelData<FArrayBox>& dfn( a_species.distributionFunction() );
   const DisjointBoxLayout& grids( dfn.disjointBoxLayout() );

   const LevelData<FArrayBox>& injected_B( geometry.getBFieldMagnitude() );
   const LevelData<FArrayBox>& real_coords( geometry.getCellCenteredRealCoords() );
   const LevelData<FArrayBox>& normalized_flux( geometry.getNormalizedMagneticFluxCell() );

   // NB: This is a cheat - there's one too many cells at the (dir,side) face
   // of the boundary box, but it doesn't matter because one-sided difference
   // will be used at that face to construct the cell average.  We do want the
   // extra cell in all other directions.
   LevelData<FArrayBox> dfn_tmp( grids, dfn.nComp(), IntVect::Unit );
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      const DataIndex& internal_dit( a_bdry_layout.dataIndex( dit ) );

      const Box& interior_box( a_bdry_layout.interiorBox( dit ) );
      Box fill_box( dfn_tmp[dit].box() );
      fill_box.growDir( a_bdry_layout.dir(), a_bdry_layout.side(), -1 );

      CH_START(t_initialize_density);
      FArrayBox density( fill_box, 1 );
      initializeField( density, *m_ic_density, geometry, real_coords[internal_dit], normalized_flux[internal_dit], interior_box, a_time );
      CH_STOP(t_initialize_density);

      CH_START(t_initialize_temperature);
      FArrayBox temperature( fill_box, 1 );
      initializeField( temperature, *m_ic_temperature, geometry, real_coords[internal_dit], normalized_flux[internal_dit], interior_box, a_time );
      CH_STOP(t_initialize_temperature);

      CH_START(t_initialize_vparallel);
      FArrayBox vparallel( fill_box, 1 );
      initializeField( vparallel, *m_ic_vparallel, geometry, real_coords[internal_dit], normalized_flux[internal_dit], interior_box, a_time );
      CH_STOP(t_initialize_vparallel);

      const PhaseBlockCoordSys& coord_sys( geometry.getBlockCoordSys( interior_box ) );
      CH_START(t_set_point_values);
      setPointValues( dfn_tmp[dit],
                      fill_box,
                      coord_sys,
                      density,
                      temperature,
                      vparallel,
                      injected_B[internal_dit],
                      real_coords[internal_dit],
                      a_species.mass() );
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
RefCountedPtr<CFG::GridFunction> createFieldFunction( const std::string& a_name,
                                                      const ParmParse& a_pp )
{
   CFG::GridFunctionLibrary* library = CFG::GridFunctionLibrary::getInstance();
   const std::string subprefix( "." + a_name );
   const std::string prefix( a_pp.prefix() + subprefix );
   ParmParse fpp( prefix.c_str() );
   std::string function_name;
   fpp.query( "function", function_name );
   return library->find( function_name );
}


inline
void MaxwellianKineticFunction::parseParameters( ParmParse& a_pp )
{
   m_ic_temperature = createFieldFunction( "temperature", a_pp );
   m_ic_density = createFieldFunction( "density", a_pp );
   
   if (a_pp.contains("vparallel.function")) {
       m_ic_vparallel = createFieldFunction( "vparallel", a_pp );
   }
   else {
       m_ic_vparallel = RefCountedPtr<CFG::GridFunction>( new CFG::Constant( 0.0, m_verbosity ) );
   }
   
   a_pp.query("enforce_input_density_profile", m_enforce_input_density_profile);
   
   if (m_verbosity) {
      printParameters();
   }
}


inline
void MaxwellianKineticFunction::initializeField( LevelData<FArrayBox>& a_field,
                                                 const CFG::GridFunction& a_ic,
                                                 const PhaseGeom& a_phase_geometry,
                                                 const CFG::IntVect& a_nghosts, 
                                                 const Real& a_time ) const
{
   const CFG::MultiBlockLevelGeom& mag_geometry( a_phase_geometry.magGeom() );
   const CFG::DisjointBoxLayout& grids( mag_geometry.grids() );
   const int SCALAR(1);
   const bool cell_average( false );

   CFG::LevelData<CFG::FArrayBox> cfg_field( grids, SCALAR, a_nghosts );
   a_ic.assign( cfg_field, mag_geometry, a_time, cell_average );
   a_phase_geometry.injectConfigurationToPhase( cfg_field, a_field );
}


inline
void MaxwellianKineticFunction::initializeField( FArrayBox& a_field,
                                                 const CFG::GridFunction& a_ic,
                                                 const PhaseGeom& a_geometry,
                                                 const FArrayBox& a_real_coords,
                                                 const FArrayBox& a_normalized_flux,
                                                 const Box& a_interior_box,
                                                 const Real& a_time ) const
{
   CH_TIMERS("MaxwellianKineticFunction::initializeField");
   CH_TIMER("initializeFieldAssign", t_assign);

   const CFG::MultiBlockLevelGeom& mag_geometry( a_geometry.magGeom() );

   const Box& box( a_field.box() );
   CFG::Box box_cfg;
   a_geometry.projectPhaseToConfiguration( box, box_cfg );

   CFG::Box interior_box_cfg;
   a_geometry.projectPhaseToConfiguration( a_interior_box, interior_box_cfg );

   CFG::FArrayBox real_coords_cfg(box_cfg,CFG_DIM);

   // Slice in the mu direction at the low mu coordinate
   SliceSpec slice_mu(MU_DIR,box.smallEnd(MU_DIR));
   CP1::FArrayBox temp1;
   sliceBaseFab((CP1::BaseFab<Real>&)temp1, (BaseFab<Real>&)a_real_coords, slice_mu);

   // Slice in the v_parallel direction at the low v_parallel coordinate
   CP1::SliceSpec slice_vp(VPARALLEL_DIR,box.smallEnd(VPARALLEL_DIR));
   CFG::FArrayBox temp2;
   sliceBaseFab((CFG::BaseFab<Real>&)temp2, (CP1::BaseFab<Real>&)temp1, slice_vp);

   real_coords_cfg.copy(temp2);
   
   const Box& flux_box = a_normalized_flux.box();

   CFG::FArrayBox normalized_flux_cfg(box_cfg,1);

   SliceSpec slice_mu2(MU_DIR, flux_box.smallEnd(MU_DIR));
   CP1::FArrayBox temp3;
   sliceBaseFab((CP1::BaseFab<Real>&)temp3, (BaseFab<Real>&)a_normalized_flux, slice_mu2);

   // Slice in the v_parallel direction at the low v_parallel coordinate
   CP1::SliceSpec slice_vp2(VPARALLEL_DIR, flux_box.smallEnd(VPARALLEL_DIR));
   CFG::FArrayBox temp4;
   sliceBaseFab((CFG::BaseFab<Real>&)temp4, (CP1::BaseFab<Real>&)temp3, slice_vp2);

   normalized_flux_cfg.copy(temp4);

   CFG::FArrayBox cfg_field( box_cfg, 1 );
   const bool cell_average( false );
   CH_START(t_assign);
   
   const CFG::MultiBlockCoordSys& coord_sys( *(mag_geometry.coordSysPtr()) );
   const int block_number( coord_sys.whichBlock( interior_box_cfg ) );
   
   a_ic.assign( cfg_field, mag_geometry, real_coords_cfg, normalized_flux_cfg, block_number, a_time, cell_average );
   CH_STOP(t_assign);

   a_geometry.injectConfigurationToPhase( cfg_field, a_field );
}


inline
void MaxwellianKineticFunction::checkGeometryValidity( const PhaseGeom& a_geometry ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   bool unknown_geom( typeid(coord_sys) != typeid(LogRectPhaseCoordSys) );
   unknown_geom &= (typeid(coord_sys) != typeid(SingleNullPhaseCoordSys));
   unknown_geom &= (typeid(coord_sys) != typeid(SNCorePhaseCoordSys));
   
   if ( unknown_geom ) {
      const std::string msg( "MaxwellianKineticFunction: Attempt to use unknown geometry. ");
      MayDay::Error( msg.c_str() );
   }
}

inline
void MaxwellianKineticFunction::setPointValues(
   FArrayBox&                a_dfn,
   const Box&                a_box,
   const PhaseBlockCoordSys& a_coord_sys,
   const FArrayBox&          a_density,
   const FArrayBox&          a_temperature,
   const FArrayBox&          a_vparallel,
   const FArrayBox&          a_B,
   const FArrayBox&          a_real_coords,
   const Real&               a_mass ) const
{
   CH_TIME("MaxwellianKineticFunction::setPointValues");

   FORT_SET_MAXWELL4D( CHF_FRA(a_dfn),
                       CHF_BOX(a_box),
                       CHF_CONST_FRA(a_real_coords),
                       CHF_CONST_FRA1(a_density,0),
                       CHF_CONST_FRA1(a_temperature,0),
                       CHF_CONST_FRA1(a_vparallel,0),
                       CHF_CONST_FRA1(a_B,0),
                       CHF_CONST_REAL(a_mass) );
}


inline
void MaxwellianKineticFunction::printParameters() const
{
   if (procID()==0) {
      std::cout << "Maxwellian kinetic function parameters:" << std::endl;
      std::cout << "-- Density --" << std::endl;
      m_ic_density->printParameters();
      std::cout << "-- Temperature --" << std::endl;
      m_ic_temperature->printParameters();
      std::cout << "-- Vparallel --" << std::endl;
      m_ic_vparallel->printParameters();
      std::cout << std::endl;
   }
}


#include "NamespaceFooter.H"
