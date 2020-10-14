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
#include "LogicalSheathBC.H"
#include "FluxBC.H"

#include <sstream>

#include "FlipGrids.H"
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
   if ( a_block_type==CFG::SingleNullBlockCoordSys::MCORE ||
        a_block_type==CFG::SingleNullBlockCoordSys::LCORE ||
        a_block_type==CFG::SingleNullBlockCoordSys::RCORE ){
      return *(m_inflow_function[RADIAL_CORE]);
   }
   else if ( a_block_type==CFG::SingleNullBlockCoordSys::LPF ||
             a_block_type==CFG::SingleNullBlockCoordSys::RPF ){
      return *(m_inflow_function[RADIAL_PF]);
   }
   else if (a_block_type==CFG::SingleNullBlockCoordSys::LSOL  ||
            a_block_type==CFG::SingleNullBlockCoordSys::MCSOL ||
            a_block_type==CFG::SingleNullBlockCoordSys::LCSOL ||
            a_block_type==CFG::SingleNullBlockCoordSys::RCSOL ||
            a_block_type==CFG::SingleNullBlockCoordSys::RSOL  ){
      return *(m_inflow_function[RADIAL_SOL]);
   }
   else {
      MayDay::Error("SingleNullPhaseBC::radialInflowData(): What on earth did you just do?  There should be no way to get this message!");
   }

   // This return will never be reached, but the compiler wants to see something returned
   return *(m_inflow_function[RADIAL_CORE]);
}


#if CFG_DIM==3
inline
KineticFunction& SingleNullPhaseBC::toroidalInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[TOROIDAL_LOWER]);
   }
   return *(m_inflow_function[TOROIDAL_UPPER]);
}
#endif


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
#if CFG_DIM==3
   else if (a_dir==TOROIDAL_DIR) {
      inflow_func = &toroidalInflowFunc( a_side );
   }
#endif
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
      MayDay::Error( "SingleNullPhaseBC::inflowFunc(): BCs in a_dir direction not implemented!" );
   }
   return *(inflow_func);
}


inline
std::string SingleNullPhaseBC::radialBcType( const Side::LoHiSide& a_side, const int& a_block_type )
{
   if ( a_block_type==CFG::SingleNullBlockCoordSys::MCORE ||
        a_block_type==CFG::SingleNullBlockCoordSys::LCORE ||
        a_block_type==CFG::SingleNullBlockCoordSys::RCORE ){
      return m_bc_type[RADIAL_CORE];
   }
   else if ( a_block_type==CFG::SingleNullBlockCoordSys::LPF ||
            a_block_type==CFG::SingleNullBlockCoordSys::RPF ){
      return m_bc_type[RADIAL_PF];
   }
   else if (a_block_type==CFG::SingleNullBlockCoordSys::LSOL  ||
            a_block_type==CFG::SingleNullBlockCoordSys::MCSOL ||
            a_block_type==CFG::SingleNullBlockCoordSys::LCSOL ||
            a_block_type==CFG::SingleNullBlockCoordSys::RCSOL ||
            a_block_type==CFG::SingleNullBlockCoordSys::RSOL  ){
      return m_bc_type[RADIAL_SOL];
   }
   else {
      MayDay::Error("SingleNullPhaseBC::radialInflowData(): What on earth did you just do?  There should be no way to get this message!");
   }
   
   // This return will never be reached, but the compiler wants to see something returned
   return m_bc_type[RADIAL_CORE];
}


#if CFG_DIM==3
inline
std::string SingleNullPhaseBC::toroidalBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[TOROIDAL_LOWER];
   }
   return m_bc_type[TOROIDAL_UPPER];
}
#endif


inline
std::string SingleNullPhaseBC::poloidalBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[POLOIDAL_INNER_DIV];
   }
   return m_bc_type[POLOIDAL_OUTER_DIV];
}

inline
std::string SingleNullPhaseBC::vParallelBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[VPAR_LOWER];
   }
   return m_bc_type[VPAR_UPPER];
}

inline
std::string SingleNullPhaseBC::muBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[MU_LOWER];
   }
   return m_bc_type[MU_UPPER];
}

inline
std::string SingleNullPhaseBC::getBcType( const int& a_dir,
                                     const Side::LoHiSide& a_side )
{
   std::string bc_type;
   if (a_dir==RADIAL_DIR) {
      MayDay::Error( "SingleNullPhaseBC::inflowFunc(): Does not handle RADIAL inflow functions!" );
   }
#if CFG_DIM==3
   else if (a_dir==TOROIDAL_DIR) {
      bc_type = toroidalBcType( a_side );
   }
#endif
   else if (a_dir==POLOIDAL_DIR) {
      bc_type = poloidalBcType( a_side );
   }
   else if (a_dir==VPARALLEL_DIR) {
      bc_type = vParallelBcType( a_side );
   }
   else if (a_dir==MU_DIR) {
      bc_type = muBcType( a_side );
   }
   else {
      MayDay::Error( "SingleNullPhaseBC::inflowFunc(): BCs in direction a_dir not implemented!" );
   }
   return bc_type;
}



/////// END INLINE METHODS /////////////////////////////////////////////////////


SingleNullPhaseBC::SingleNullPhaseBC( const std::string& a_name,
                                      ParmParse& a_pp,
                                      const int& a_verbosity )
   : m_name(a_name),
     m_verbosity(a_verbosity),
     m_logical_sheath(false),
     m_flux_bc(false),
     m_all_bdry_defined(false),
     m_pp(a_pp)
{
   m_inflow_function.resize( NUM_INFLOW );
   m_bc_type.resize( NUM_INFLOW );
   m_bdry_name.resize( NUM_INFLOW );

   m_bdry_name.resize( NUM_INFLOW );
   m_bdry_name[RADIAL_CORE] = "radial_core";
   m_bdry_name[RADIAL_SOL] = "radial_sol";
   m_bdry_name[RADIAL_PF] = "radial_pf";
#if CFG_DIM==3
   m_bdry_name[TOROIDAL_LOWER] = "toroidal_lower";
   m_bdry_name[TOROIDAL_UPPER] = "toroidal_upper";
#endif
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
void SingleNullPhaseBC::fillInflowData( const Real& a_time )
{
   CH_TIMERS("SingleNullPhaseBC::fillInflowData");
   CH_TIMER("inflow_func.assign", t_inflow_func_assign);

   for (int i(0); i<m_all_bdry_layouts.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(m_all_bdry_layouts[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      KineticSpecies& bdry_data( *(m_all_bdry_data[i]) );

      if (dir==RADIAL_DIR) {
         const PhaseGeom& geometry( bdry_data.phaseSpaceGeometry() );
         const SingleNullPhaseCoordSys& coord_sys(
         dynamic_cast<const SingleNullPhaseCoordSys&>( geometry.phaseCoordSys()) );
         const DisjointBoxLayout& grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
            const Box& interior_box( bdry_layout.interiorBox( dit ) );
            const int block( coord_sys.whichBlock( interior_box ) );
            const int block_type( coord_sys.poloidalBlock( block ) );
            KineticFunction& inflow_func( radialInflowFunc( side, block_type ) );
            CH_START(t_inflow_func_assign);
            inflow_func.assign( bdry_data, bdry_layout, a_time );
            CH_STOP(t_inflow_func_assign);
         }
      }
      else {
         KineticFunction& inflow_func( inflowFunc( dir, side ) );
         CH_START(t_inflow_func_assign);
         inflow_func.assign( bdry_data, bdry_layout, a_time );
         CH_STOP(t_inflow_func_assign);
      }
   }
}

void SingleNullPhaseBC::apply( KineticSpeciesPtrVect&  a_species,
                               const int&  a_species_index,
                               const CFG::LevelData<CFG::FArrayBox>& a_phi,
                               const LevelData<FluxBox>& a_velocity,
                               const Real& a_time )
{
   CH_TIMERS("SingleNullPhaseBC::apply");
   CH_TIMER("defineBoundaryBoxLayouts", t_define_boundary_box_layouts);
   CH_TIMER("defineInflowDataStorage", t_define_inflow_data_storage);
   CH_TIMER("fillInflowData", t_fill_inflow_data);
   CH_TIMER("setInflowOutflowBC", t_set_inflow_outflow_BC);
   CH_TIMER("setCodimBoundaryValues", t_set_codim_boundary_values);

   KineticSpecies& species_physical( *(a_species[a_species_index]) );
   
   const PhaseGeom& geometry( species_physical.phaseSpaceGeometry() );
   const SingleNullPhaseCoordSys& coord_sys(
      dynamic_cast<const SingleNullPhaseCoordSys&>( geometry.phaseCoordSys()) );

   LevelData<FArrayBox>& Bf( species_physical.distributionFunction() );
   const DisjointBoxLayout& grids( Bf.disjointBoxLayout() );
   const IntVect& ghost_vect( Bf.ghostVect() );

   if (m_all_bdry_defined == false) {
      CH_START(t_define_boundary_box_layouts);
      PhaseBCUtils::defineBoundaryBoxLayouts(m_all_bdry_layouts,
                                             grids,
                                             coord_sys,
                                             ghost_vect );
      CH_STOP(t_define_boundary_box_layouts);
      
      CH_START(t_define_inflow_data_storage);
      PhaseBCUtils::defineInflowDataStorage(m_all_bdry_data,
                                            m_all_bdry_layouts,
                                            species_physical );
      CH_STOP(t_define_inflow_data_storage);
      
      CH_START(t_fill_inflow_data);
      Vector<std::string> all_bc_type;
      fillInflowData( a_time );
      setAllBcType( m_all_bdry_layouts );
      m_all_bdry_defined = true;
      CH_STOP(t_fill_inflow_data);
   }
   
   CH_START(t_set_inflow_outflow_BC);
   PhaseBCUtils::setInflowOutflowBC(Bf,
                                    m_all_bdry_layouts,
                                    m_all_bdry_data,
                                    m_all_bc_type,
                                    coord_sys,
                                    a_velocity );
   CH_STOP(t_set_inflow_outflow_BC);

   if (m_logical_sheath) {
       applySheathBC(a_species, a_species_index, a_phi, a_velocity);
   }
   
   if (m_flux_bc) {
      applyFluxBC(species_physical, a_phi, a_velocity, a_time);
   }
   
   CH_START(t_set_codim_boundary_values);
   // interpolate all other codim boundaries
   CodimBC::setCodimBoundaryValues( Bf, coord_sys );
   CH_STOP(t_set_codim_boundary_values);
}

inline
void SingleNullPhaseBC::setAllBcType( const BoundaryBoxLayoutPtrVect&  a_bdry_layout )
{
   m_all_bc_type.resize(a_bdry_layout.size());
   for (int i(0); i<a_bdry_layout.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(a_bdry_layout[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      if (dir == RADIAL_DIR) {
         KineticSpecies& bdry_data( *(m_all_bdry_data[i]) );
         const PhaseGeom& geometry( bdry_data.phaseSpaceGeometry() );
         const SingleNullPhaseCoordSys& coord_sys(dynamic_cast<const SingleNullPhaseCoordSys&>( geometry.phaseCoordSys()) );
         const DisjointBoxLayout& grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
            const Box& interior_box( bdry_layout.interiorBox( dit ) );
            const int block( coord_sys.whichBlock( interior_box ) );
            const int block_type( coord_sys.poloidalBlock( block ) );
            m_all_bc_type[i] = radialBcType(side, block_type);
         }
      }
      else {
         m_all_bc_type[i] = getBcType(dir, side);
      }
   }
}


void SingleNullPhaseBC::applySheathBC( KineticSpeciesPtrVect&  a_species,
                                       const int&  a_species_index,
                                       const CFG::LevelData<CFG::FArrayBox>& a_phi,
                                       const LevelData<FluxBox>& a_velocity)
{
   // Get coordinate system parameters
   KineticSpecies& species_physical( *(a_species[a_species_index]) );
   const PhaseGeom& geometry( species_physical.phaseSpaceGeometry() );
   const MultiBlockCoordSys& coord_sys( *(geometry.coordSysPtr()) );
   const PhaseGrid& phase_grid = geometry.phaseGrid();
   const DisjointBoxLayout& dbl = phase_grid.disjointBoxLayout();
   
   //For sheath BC calculations we only need one-cell-wide
   //boundary storage, so define it here
   BoundaryBoxLayoutPtrVect all_bdry_layouts;
   PhaseBCUtils::defineBoundaryBoxLayouts(all_bdry_layouts,
                                          dbl,
                                          coord_sys,
                                          IntVect::Unit );

   // Loop over boundaries
   for (int i(0); i<all_bdry_layouts.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[i]) );
      const int& dir( bdry_layout.dir() );
      if (dir == POLOIDAL_DIR) {
         LogicalSheathBC sheathBC(all_bdry_layouts[i], m_pp);
         sheathBC.applyBC(a_species, a_species_index, a_velocity, a_phi);
      }
   }
}


void SingleNullPhaseBC::applyFluxBC(KineticSpecies& a_species,
                                    const CFG::LevelData<CFG::FArrayBox>& a_phi,
                                    const LevelData<FluxBox>& a_velocity,
                                    const Real& a_time)
{
   // Get coordinate system parameters
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   const MultiBlockCoordSys& coord_sys( *(geometry.coordSysPtr()) );
   const PhaseGrid& phase_grid = geometry.phaseGrid();
   const DisjointBoxLayout& dbl = phase_grid.disjointBoxLayout();
   
   //For flux BC calculations we only need one-cell-wide
   //boundary storage, so define it here
   BoundaryBoxLayoutPtrVect all_bdry_layouts;
   PhaseBCUtils::defineBoundaryBoxLayouts(all_bdry_layouts,
                                          dbl,
                                          coord_sys,
                                          IntVect::Unit );

   
   //Get boundary layout for the inner core boundary
   BoundaryBoxLayoutPtrVect flux_bdry_layouts;
   for (int i(0); i<all_bdry_layouts.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      bool include_bdry(false);
      if (dir == RADIAL_DIR && side == Side::Lo) {
          KineticSpecies& bdry_data( *(m_all_bdry_data[i]) );
          const PhaseGeom& geometry( bdry_data.phaseSpaceGeometry() );
          const SingleNullPhaseCoordSys& coord_sys(dynamic_cast<const SingleNullPhaseCoordSys&>( geometry.phaseCoordSys()) );
          const DisjointBoxLayout& grids( bdry_layout.disjointBoxLayout() );
          for (LayoutIterator lit( grids.layoutIterator() ); lit.ok(); ++lit) {
             // This is one-cell-wide boundary box
             Box box_tmp = grids[lit];
             // Shift boudary box into the physical domain (in the radial dir)
             box_tmp.shift(RADIAL_DIR, 1);
             const int block( coord_sys.whichBlock( box_tmp ) );
             const int block_type( coord_sys.poloidalBlock( block ) );
             
             if (block_type==CFG::SingleNullBlockCoordSys::MCORE  ||
                 block_type==CFG::SingleNullBlockCoordSys::LCORE  ||
                 block_type==CFG::SingleNullBlockCoordSys::RCORE) {
                
                include_bdry = true;
             }
          }
         
          if (include_bdry) {
            flux_bdry_layouts.push_back(all_bdry_layouts[i]);
          }
      }
   }
   
   //Apply flux BC to the inner radial boundary
   std::string prefix( m_pp.prefix() );
   prefix += "." + m_bdry_name[RADIAL_CORE];
   ParmParse pp( prefix.c_str() );
   FluxBC fluxBC(flux_bdry_layouts, pp);
   fluxBC.applyBC(a_species, a_velocity, a_phi, a_time);
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
      std::cout << "  logical_sheath  =  " << m_logical_sheath << std::endl;
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
      
      fpp.query( "type", m_bc_type[i] );
      
      if (m_bc_type[i] == "flux") {
         m_flux_bc = true;
      }
   }

   a_pp.query("logical_sheath",m_logical_sheath);

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



