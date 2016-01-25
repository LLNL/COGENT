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
#include "MillerPhaseCoordSys.H"
#include "SlabPhaseCoordSys.H"
#include "MultiBlockCoordSys.H"
#include "PhaseBlockCoordSys.H"
#include "SNCorePhaseCoordSys.H"
#include "SingleNullPhaseCoordSys.H"
#include "Vector.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunctionLibrary.H"
#include "Constant.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"



MaxwellianKineticFunction::MaxwellianKineticFunction( ParmParse& a_pp,
                                                      const int& a_verbosity )
   : m_verbosity(a_verbosity)
{
   parseParameters( a_pp );
}


void MaxwellianKineticFunction::convertToCellAverage(
   const MultiBlockCoordSys&  a_coord_sys,
   LevelData<FArrayBox>&      a_dfn ) const
{
   LevelData<FArrayBox> dfn_tmp(a_dfn.disjointBoxLayout(),
                                a_dfn.nComp(),
                                a_dfn.ghostVect()+IntVect::Unit);

   const DisjointBoxLayout& grids( a_dfn.disjointBoxLayout() );

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      dfn_tmp[dit].copy( a_dfn[dit] );
   }
   dfn_tmp.exchange();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const int block_number( a_coord_sys.whichBlock( grids[dit] ) );
      const PhaseBlockCoordSys* coord_sys
         = dynamic_cast<const PhaseBlockCoordSys*>( a_coord_sys.getCoordSys( block_number ) );

      fourthOrderAverageCell( dfn_tmp[dit], coord_sys->domain(), grids[dit] );
   }
   dfn_tmp.exchange();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_dfn[dit].copy( dfn_tmp[dit] );
   }
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

   LevelData<FArrayBox> injected_density;
   initializeField( injected_density, *m_ic_density, geometry, a_time );
   CH_assert( isNonNegative( injected_density ) );

   LevelData<FArrayBox> injected_temperature;
   initializeField( injected_temperature, *m_ic_temperature, geometry, a_time );
   CH_assert( isPositiveDefinite( injected_temperature ) );

   LevelData<FArrayBox> injected_vparallel;
   initializeField( injected_vparallel, *m_ic_vparallel, geometry, a_time );
   //CH_assert( isPositiveDefinite( injected_vparallel ) );
 
   const LevelData<FArrayBox>& injected_B( geometry.getBFieldMagnitude() );

   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
      setPointValues( dfn[dit],
                      dfn[dit].box(),
                      geometry.getBlockCoordSys( grids[dit] ),
                      injected_density[dit],
                      injected_temperature[dit],
                      injected_vparallel[dit],
                      injected_B[dit],
                      a_species.mass() );
   }

   geometry.multBStarParallel( dfn );
   if ( !(geometry.secondOrder()) )  {
      convertToCellAverage( *geometry.coordSysPtr(), dfn );
   }
   geometry.multJonValid( dfn );
   dfn.exchange();
}


void MaxwellianKineticFunction::assign( KineticSpecies& a_species,
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

      FArrayBox density( fill_box, 1 );
      initializeField( density, *m_ic_density, geometry, interior_box, a_time );

      FArrayBox temperature( fill_box, 1 );
      initializeField( temperature, *m_ic_temperature, geometry, interior_box, a_time );

      FArrayBox vparallel( fill_box, 1 );
      initializeField( vparallel, *m_ic_vparallel, geometry, interior_box, a_time );

      const PhaseBlockCoordSys& coord_sys( geometry.getBlockCoordSys( interior_box ) );
      const DataIndex& internal_dit( a_bdry_layout.dataIndex( dit ) );
      setPointValues( dfn_tmp[dit],
                      fill_box,
                      coord_sys,
                      density,
                      temperature,
                      vparallel,
                      injected_B[internal_dit],
                      a_species.mass() );
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
    
   if (m_verbosity) {
      printParameters();
   }
}


inline
void MaxwellianKineticFunction::initializeField( LevelData<FArrayBox>& a_field,
                                                 const CFG::GridFunction& a_ic,
                                                 const PhaseGeom& a_phase_geometry,
                                                 const Real& a_time ) const
{
   const CFG::MultiBlockLevelGeom& mag_geometry( a_phase_geometry.magGeom() );
   const CFG::DisjointBoxLayout& grids( mag_geometry.grids() );
   const int SCALAR(1);
   const CFG::IntVect nghosts(CFG::IntVect::Unit);
   const bool cell_average( false );

   CFG::LevelData<CFG::FArrayBox> cfg_field( grids, SCALAR, nghosts );
   a_ic.assign( cfg_field, mag_geometry, a_time, cell_average );
   a_phase_geometry.injectConfigurationToPhase( cfg_field, a_field );
}


inline
void MaxwellianKineticFunction::initializeField( FArrayBox& a_field,
                                                 const CFG::GridFunction& a_ic,
                                                 const PhaseGeom& a_geometry,
                                                 const Box& a_interior_box,
                                                 const Real& a_time ) const
{
   const CFG::MultiBlockLevelGeom& mag_geometry( a_geometry.magGeom() );

   const Box& box( a_field.box() );
   CFG::Box box_cfg;
   a_geometry.projectPhaseToConfiguration( box, box_cfg );

   CFG::Box interior_box_cfg;
   a_geometry.projectPhaseToConfiguration( a_interior_box, interior_box_cfg );

   CFG::FArrayBox cfg_field( box_cfg, 1 );
   const bool cell_average( false );
   a_ic.assign( cfg_field, mag_geometry, interior_box_cfg, a_time, cell_average );

   a_geometry.injectConfigurationToPhase( cfg_field, a_field );
}


inline
void MaxwellianKineticFunction::checkGeometryValidity( const PhaseGeom& a_geometry ) const
{
   const MultiBlockCoordSys& an_coord_sys( *(a_geometry.coordSysPtr()) );
   bool not_annular (typeid(an_coord_sys) != typeid(MillerPhaseCoordSys));
   not_annular &= (typeid(an_coord_sys) != typeid(SlabPhaseCoordSys));

   const MultiBlockCoordSys& sn_coord_sys( *(a_geometry.coordSysPtr()) );
   bool not_single_null( typeid(sn_coord_sys) != typeid(SingleNullPhaseCoordSys) );
   bool not_sncore( typeid(sn_coord_sys) != typeid(SNCorePhaseCoordSys) );

   if ( not_annular && not_single_null && not_sncore ) {
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
   const Real&               a_mass ) const
{
   FArrayBox cell_center_coords( a_box, PDIM );
   a_coord_sys.getCellCenteredRealCoords( cell_center_coords );
   FORT_SET_MAXWELL4D( CHF_FRA(a_dfn),
                       CHF_BOX(a_box),
                       CHF_CONST_FRA(cell_center_coords),
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
