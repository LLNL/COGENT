#include "SingleNullPhaseBC.H"

#include "ConstFact.H"
#include "Directions.H"
#include "FourthOrderUtil.H"
#include "KineticFunctionLibrary.H"
#include "LevelData.H"
#include "PhaseGeom.H"
#include "SingleNullPhaseCoordSys.H"
#include "SPMD.H"
#include "PhaseBCUtils.H"

#include <sstream>

#include "flipGrids.H"
#include "PoloidalBCF_F.H"

#include "EdgeToCell.H"

#include "FourthOrderBC.H.multidim"
#include "CodimBC.H.multidim"
#include "BCUtils.H.multidim"

#include "SingleNullPhaseBCF_F.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "SingleNullBlockCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM


#include "NamespaceHeader.H"

/////// BEGIN INLINE METHODS ///////////////////////////////////////////////////

inline
KineticFunction& SingleNullPhaseBC::radialInflowFunc( const Side::LoHiSide& a_side, const int& a_block_type )
{
   if ( a_block_type==MCORE ||
        a_block_type==LCORE ||
        a_block_type==RCORE ){
      return *(m_inflow_function[RADIAL_CORE]);
   }
   else if ( a_block_type==LPF ||
             a_block_type==RPF ){
      return *(m_inflow_function[RADIAL_PF]);
   }
   else if (a_block_type==LSOL  ||
            a_block_type==MCSOL ||
            a_block_type==LCSOL ||
            a_block_type==RCSOL ||
            a_block_type==RSOL  ){
      return *(m_inflow_function[RADIAL_SOL]);
   }
   else {
      MayDay::Error("SingleNullPhaseBC::radialInflowData(): What on earth did you just do?  There should be no way to get this message!");
   }

   // This return will never be reached, but the compiler wants to see something returned
   return *(m_inflow_function[RADIAL_CORE]);
}


inline
KineticFunction& SingleNullPhaseBC::poloidalInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[POLOIDAL_OUTER_DIV]);
   }
   return *(m_inflow_function[POLOIDAL_INNER_DIV]);
}


inline
KineticFunction& SingleNullPhaseBC::vParallelInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[VPAR_LOWER]);
   }
   return *(m_inflow_function[VPAR_UPPER]);
}


inline
KineticFunction& SingleNullPhaseBC::muInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[MU_LOWER]);
   }
   return *(m_inflow_function[MU_UPPER]);
}


inline
KineticFunction& SingleNullPhaseBC::inflowFunc( const int& a_dir,
                                                const Side::LoHiSide& a_side )
{
   KineticFunction* inflow_func;
   if (a_dir==RADIAL_DIR) {
      MayDay::Error( "SingleNullPhaseBC::inflowFunc(): Does not handle RADIAL inflow functions!" );
   }
   else if (a_dir==POLOIDAL_DIR) {
      inflow_func = &poloidalInflowFunc( a_side );
   }
   else if (a_dir==VPARALLEL_DIR) {
      inflow_func = &vParallelInflowFunc( a_side );
   }
   else if (a_dir==MU_DIR) {
      inflow_func = &muInflowFunc( a_side );
   }
   else {
      MayDay::Error( "SingleNullPhaseBC::inflowFunc(): Toroidal BCs not implemented!" );
   }
   return *(inflow_func);
}

/////// END INLINE METHODS /////////////////////////////////////////////////////


SingleNullPhaseBC::SingleNullPhaseBC( const std::string& a_name,
                                      ParmParse& a_pp,
                                      const int& a_verbosity )
   : m_name(a_name),
     m_verbosity(a_verbosity)
{
   m_inflow_function.resize( NUM_INFLOW );

   m_bdry_name.resize( NUM_INFLOW );
   m_bdry_name[RADIAL_CORE] = "radial_core";
   m_bdry_name[RADIAL_SOL] = "radial_sol";
   m_bdry_name[RADIAL_PF] = "radial_pf";
   m_bdry_name[POLOIDAL_INNER_DIV] = "poloidal_inner_div";
   m_bdry_name[POLOIDAL_OUTER_DIV] = "poloidal_outer_div";
   m_bdry_name[VPAR_LOWER] = "vpar_lower";
   m_bdry_name[VPAR_UPPER] = "vpar_upper";
   m_bdry_name[MU_LOWER] = "mu_lower";
   m_bdry_name[MU_UPPER] = "mu_upper";

   parseParameters( a_pp );
}


SingleNullPhaseBC::~SingleNullPhaseBC()
{
#if 0
   for (int i(0); i<m_inflow_function.size(); i++) {
      if (!(m_inflow_function[i])) {
         delete m_inflow_function[i];
      }
   }
#endif
}


inline
void SingleNullPhaseBC::fillInflowData( KineticSpeciesPtrVect& a_bdry_data,
                                        const BoundaryBoxLayoutPtrVect& a_bdry_layout,
                                        const Real& a_time )
{
   for (int i(0); i<a_bdry_layout.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(a_bdry_layout[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      KineticSpecies& bdry_data( *(a_bdry_data[i]) );

      if (dir==RADIAL_DIR) {
         const PhaseGeom& geometry( bdry_data.phaseSpaceGeometry() );
         const SingleNullPhaseCoordSys& coord_sys(
            dynamic_cast<const SingleNullPhaseCoordSys&>( geometry.phaseCoordSys()) );
         const DisjointBoxLayout& grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
            const Box& interior_box( bdry_layout.interiorBox( dit ) );
            const int block( coord_sys.whichBlock( interior_box ) );
            const int block_type( coord_sys.blockType( block ) );
            KineticFunction& inflow_func( radialInflowFunc( side, block_type ) );
            inflow_func.assign( bdry_data, bdry_layout, a_time );
         }
      }
      else {
         KineticFunction& inflow_func( inflowFunc( dir, side ) );
         inflow_func.assign( bdry_data, bdry_layout, a_time );
      }
   }
}



void SingleNullPhaseBC::apply( KineticSpecies& a_species_comp,
                               const CFG::LevelData<CFG::FArrayBox>& a_phi,
                               const LevelData<FluxBox>& a_velocity,
                               const Real& a_time )
{
   const PhaseGeom& geometry( a_species_comp.phaseSpaceGeometry() );
   const SingleNullPhaseCoordSys& coord_sys(
      dynamic_cast<const SingleNullPhaseCoordSys&>( geometry.phaseCoordSys()) );

   LevelData<FArrayBox>& u( a_species_comp.distributionFunction() );
   const DisjointBoxLayout& grids( u.disjointBoxLayout() );
   const IntVect& ghost_vect( u.ghostVect() );

   BoundaryBoxLayoutPtrVect all_bdry_layouts;
   PhaseBCUtils::defineBoundaryBoxLayouts( all_bdry_layouts,
                                           grids,
                                           coord_sys,
                                           ghost_vect );

   KineticSpeciesPtrVect all_bdry_data;
   PhaseBCUtils::defineInflowDataStorage( all_bdry_data,
                                          all_bdry_layouts,
                                          a_species_comp );

   fillInflowData( all_bdry_data, all_bdry_layouts, a_time );

   PhaseBCUtils::setInflowOutflowBC( u,
                                     all_bdry_layouts,
                                     all_bdry_data,
                                     coord_sys,
                                     a_velocity );
   
   // interpolate all other codim boundaries
   CodimBC::setCodimBoundaryValues( u, coord_sys );
}


void SingleNullPhaseBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "SingleNullPhaseBC =============================" << std::endl;
      std::cout << "- variable: "  << m_name << "-------------" << std::endl;
      for (int i(0); i<m_inflow_function.size(); i++) {
         std::cout << "  " << m_bdry_name[i] << ": " << std::endl;
         m_inflow_function[i]->printParameters();
      }
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << "===============================================" << std::endl;
   }
}


inline
void SingleNullPhaseBC::parseParameters( ParmParse& a_pp )
{
   KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   for (int i(0); i<m_inflow_function.size(); i++) {
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );
      std::string function_name;
      fpp.query( "function", function_name );
      m_inflow_function[i] = library->find( function_name );
   }

   if (m_verbosity) {
      printParameters();
   }
}



#if 0
void SingleNullPhaseBC::ghostCellBC( LevelData<FArrayBox>&      a_soln,
                                     const MultiBlockLevelGeom& a_geometry,
                                     const Real                 a_dx,
                                     const Real                 a_time,
                                     const BCAuxData&           a_aux_data  )
{
   const LevelData<FluxBox>& face_vel( a_aux_data.getFluxBoxData( "Velocity" ) );

   BoundaryBoxData inflow_data( a_soln, a_geometry );
   fillInflowData( inflow_data, a_soln, a_geometry );
   applyPoloidalReflections( inflow_data, a_soln, a_geometry, a_aux_data );

   const MultiBlockCoordSys* coord_sys( a_geometry.coordSysPtr() );
   const DisjointBoxLayout& grids( a_soln.disjointBoxLayout() );
   for (DataIterator dit( grids ); dit.ok(); ++dit) {

      for (int dir(RADIAL_DIR); dir<=MU_DIR; dir++) {
         for (SideIterator si; si.ok(); ++si) {
            Side::LoHiSide side( si() );

            const Box& this_box(grids[dit]);
            if (BCUtils::isPhysicalBoundary( *coord_sys, this_box, dir, side )) {

               const Box& domain_box( BCUtils::blockDomain( this_box, *coord_sys ) );
               const IntVect& ghost_vect( a_soln.ghostVect() );
               Box boundary_box( BCUtils::getGhostBox( domain_box,
                                                       this_box,
                                                       dir,
                                                       side,
                                                       ghost_vect[dir]) );

               if (!boundary_box.isEmpty()) {
#if (PDIM>4)
                  if (dir!=TOROIDAL_DIR) {
#endif
                     FArrayBox& this_soln( a_soln[dit] );
                     FArrayBox& this_inflow_data( inflow_data.data( dit, dir, side ) );
                     this_soln.copy( this_inflow_data, fill_box );
#if 1
                     const FluxBox& this_face_vel( face_vel[dit] );
                     FourthOrderBC::setInflowOutflowBC( this_soln,
                                                        boundary_box,
                                                        this_inflow_data,
                                                        this_face_vel,
                                                        dir,
                                                        side );
#endif
#if (PDIM>4)
                  }
                  else {
                     MayDay::Error( "SingleNullPhaseBC: Toroidal BCs not implemented!" );
                  }
#endif
               }
            }
         }
      }

      CodimBC::setCodimBoundaryValues( a_soln, *coord_sys );
   }
}


void SingleNullPhaseBC::createReflectedData( LevelData<FArrayBox>&       a_rflct_data,
                                             AuxDataMap&                 a_rflct_data_map,
                                             const LevelData<FArrayBox>& a_soln,
                                             const MultiBlockCoordSys&   a_coord_sys,
                                             const Side::LoHiSide&       a_side ) const
{
   const int REFLECT_DIR(VPARALLEL_DIR);
   const int REFLECT_COORD(0);
   const ReflectBoxMap rbm( REFLECT_DIR, REFLECT_COORD );
   const DisjointBoxLayout& grids( a_soln.disjointBoxLayout() );
   const int DIR(POLOIDAL_DIR);
   DisjointBoxLayout rflct_grids;
   buildAuxGhostBoxLayout( rflct_grids, grids, a_coord_sys, rbm, a_soln.ghostVect(), DIR, a_side );
   a_rflct_data_map.define( grids, rflct_grids, rbm, DIR, a_side );
   a_rflct_data.define( rflct_grids, a_soln.nComp() );
   a_soln.copyTo( a_rflct_data );
}


void SingleNullPhaseBC::applyPoloidalReflections(
   BoundaryBoxData&            a_inflow_data,
   const LevelData<FArrayBox>& a_soln,
   const MultiBlockLevelGeom&  a_geometry,
   const BCAuxData&            a_aux_data ) const
{
   const LevelData<FluxBox>& face_vel( a_aux_data.getFluxBoxData( "Velocity" ) );
   const LevelData<FArrayBox>& qphi_injected( a_aux_data.getFArrayBoxData( "qPhi" ) );
   const MultiBlockCoordSys* coord_sys( a_geometry.coordSysPtr() );

   for (SideIterator si; si.ok(); ++si) {
      const Side::LoHiSide side(si());

      LevelData<FArrayBox>& this_inflow_data( a_inflow_data.data( POLOIDAL_DIR, side ) );
      for (DataIterator idit( this_inflow_data.dataIterator() ); idit.ok(); ++idit) {

         const Box& this_box( this_inflow_data[idit].box() );
         if (BCUtils::isPhysicalBoundary( *coord_sys, this_box, POLOIDAL_DIR, side )) {

            LevelData<FArrayBox> rflct_data;
            AuxDataMap rflct_data_map;
            createReflectedData( rflct_data, rflct_data_map, a_soln, *coord_sys, side );

            DataIndex sdit( a_inflow_data.mapToSrcDataIndex( idit, POLOIDAL_DIR, side ) );
            DataIndex rdit( rflct_data_map.auxDataIteratorToDataIndex( sdit ) );

            const FluxBox& this_face_vel( face_vel[sdit] );
            Box face_box( computeFaceBox( this_box, this_face_vel.box(), POLOIDAL_DIR, side ) );
            FArrayBox face_velocity( face_box, 1 );
            computeNormalVelocity( face_velocity, this_face_vel, POLOIDAL_DIR );

            const int SIDE(side);
            FORT_SET_POLOIDAL_DIVERTER_BC( CHF_FRA(this_inflow_data[idit]),
                                           CHF_BOX(this_box),
                                           CHF_CONST_FRA(rflct_data[rdit]),
                                           CHF_CONST_FRA(face_velocity),
                                           CHF_CONST_FRA1(qphi_injected[sdit],1),
                                           CHF_CONST_INT(SIDE) );
         }
      }
   }
}
#endif

#include "NamespaceFooter.H"



