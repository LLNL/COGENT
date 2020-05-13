#include "LogRectPhaseBC.H"
#include "ConstFact.H"
#include "Directions.H"
#include "FourthOrderUtil.H"
#include "KineticFunctionLibrary.H"
#include "LevelData.H"
#include "PhaseGeom.H"
#include "SPMD.H"
#include "LogicalSheathBC.H"
#include "FluxBC.H"
#include "PhaseBCUtils.H"

#include <sstream>

#include "EdgeToCell.H"

#include "inspect.H"

#include "FourthOrderBC.H.multidim"
#include "CodimBC.H.multidim"
#include "BCUtils.H.multidim"


#include "NamespaceHeader.H"

/////// BEGIN INLINE METHODS ///////////////////////////////////////////////////

inline
KineticFunction& LogRectPhaseBC::radialInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[RADIAL_LOWER]);
   }
   return *(m_inflow_function[RADIAL_UPPER]);
}

inline
KineticFunction& LogRectPhaseBC::poloidalInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[POLOIDAL_LOWER]);
   }
   return *(m_inflow_function[POLOIDAL_UPPER]);
}

#if CFG_DIM==3
inline
KineticFunction& LogRectPhaseBC::toroidalInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[TOROIDAL_LOWER]);
   }
   return *(m_inflow_function[TOROIDAL_UPPER]);
}
#endif

inline
KineticFunction& LogRectPhaseBC::vParallelInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[VPAR_LOWER]);
   }
   return *(m_inflow_function[VPAR_UPPER]);
}


inline
KineticFunction& LogRectPhaseBC::muInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[MU_LOWER]);
   }
   return *(m_inflow_function[MU_UPPER]);
}


inline
KineticFunction& LogRectPhaseBC::inflowFunc( const int& a_dir,
                                              const Side::LoHiSide& a_side )
{
   KineticFunction* inflow_func;
   if (a_dir==RADIAL_DIR) {
      inflow_func = &radialInflowFunc( a_side );
   }
   else if (a_dir==POLOIDAL_DIR) {
      inflow_func = &poloidalInflowFunc( a_side );
   }
#if CFG_DIM==3
   else if (a_dir==TOROIDAL_DIR) {
      inflow_func = &toroidalInflowFunc( a_side );
   }
#endif
   else if (a_dir==VPARALLEL_DIR) {
      inflow_func = &vParallelInflowFunc( a_side );
   }
   else if (a_dir==MU_DIR) {
      inflow_func = &muInflowFunc( a_side );
   }
   else {
      MayDay::Error( "LogRectPhaseBC: BCs not implemented!" );
   }
   return *(inflow_func);
}

inline
std::string LogRectPhaseBC::radialBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[RADIAL_LOWER];
   }
   return m_bc_type[RADIAL_UPPER];
}

inline
std::string LogRectPhaseBC::poloidalBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[POLOIDAL_LOWER];
   }
   return m_bc_type[POLOIDAL_UPPER];
}

#if CFG_DIM==3
inline
std::string LogRectPhaseBC::toroidalBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[TOROIDAL_LOWER];
   }
   return m_bc_type[TOROIDAL_UPPER];
}
#endif

inline
std::string LogRectPhaseBC::vParallelBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[VPAR_LOWER];
   }
   return m_bc_type[VPAR_UPPER];
}

inline
std::string LogRectPhaseBC::muBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[MU_LOWER];
   }
   return m_bc_type[MU_UPPER];
}

inline
std::string LogRectPhaseBC::getBcType( const int& a_dir,
                                       const Side::LoHiSide& a_side )
{
   std::string bc_type;
   if (a_dir==RADIAL_DIR) {
      bc_type = radialBcType( a_side );
   }
   else if (a_dir==POLOIDAL_DIR) {
      bc_type = poloidalBcType( a_side );
   }
#if CFG_DIM==3
   else if (a_dir==TOROIDAL_DIR) {
      bc_type = toroidalBcType( a_side );
   }
#endif
   else if (a_dir==VPARALLEL_DIR) {
      bc_type = vParallelBcType( a_side );
   }
   else if (a_dir==MU_DIR) {
      bc_type = muBcType( a_side );
   }
   else {
      MayDay::Error( "LogRectPhaseBC: BCs not implemented!" );
   }
   return bc_type;
}

/////// END INLINE METHODS /////////////////////////////////////////////////////


LogRectPhaseBC::LogRectPhaseBC( const std::string& a_name,
                                ParmParse& a_pp,
                                const int& a_verbosity )
   : m_name(a_name),
     m_verbosity(a_verbosity),
     m_pp(a_pp),
     m_all_bdry_defined(false),
     m_logical_sheath_bc(false),
     m_flux_bc(false)

{
   m_inflow_function.resize( NUM_INFLOW );
   m_bc_type.resize( NUM_INFLOW );
   m_bdry_name.resize( NUM_INFLOW );
   
   m_bdry_name[RADIAL_LOWER] = "radial_lower";
   m_bdry_name[RADIAL_UPPER] = "radial_upper";
   m_bdry_name[POLOIDAL_LOWER] = "poloidal_lower";
   m_bdry_name[POLOIDAL_UPPER] = "poloidal_upper";
#if CFG_DIM==3
   m_bdry_name[TOROIDAL_LOWER] = "toroidal_lower";
   m_bdry_name[TOROIDAL_UPPER] = "toroidal_upper";
#endif
   m_bdry_name[VPAR_LOWER] = "vpar_lower";
   m_bdry_name[VPAR_UPPER] = "vpar_upper";
   m_bdry_name[MU_LOWER] = "mu_lower";
   m_bdry_name[MU_UPPER] = "mu_upper";

   parseParameters( a_pp );
}


LogRectPhaseBC::~LogRectPhaseBC()
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
void LogRectPhaseBC::fillInflowData(KineticSpecies& a_species,
                                    const CFG::LevelData<CFG::FArrayBox>& a_phi,
                                    const LevelData<FluxBox>& a_velocity,
                                    const Real& a_time )
{
   CH_TIMERS("LogRectPhaseBC::fillInflowData");
   CH_TIMER("inflow_func.assign", t_inflow_func_assign);

   for (int i(0); i<m_all_bdry_layouts.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(m_all_bdry_layouts[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      KineticSpecies& bdry_data( *(m_all_bdry_data[i]) );
      KineticFunction& inflow_func( inflowFunc( dir, side ) );
      CH_START(t_inflow_func_assign);
      inflow_func.assign( bdry_data, bdry_layout, a_time );
      CH_STOP(t_inflow_func_assign);
   }
}

inline
Box computeFillBox( const Box&                a_interior_box,
                    const MultiBlockCoordSys& a_coord_sys,
                    const int&                a_dir,
                    const Side::LoHiSide&     a_side,
                    const int&                a_n_ghost )
{
   const Box& domain_box( BCUtils::blockDomain( a_interior_box, a_coord_sys ) );
   return Box( BCUtils::getGhostBox( domain_box,
                                     a_interior_box,
                                     a_dir,
                                     a_side,
                                     a_n_ghost ) );
}

void LogRectPhaseBC::apply(KineticSpecies& a_species_comp,
                           const CFG::LevelData<CFG::FArrayBox>& a_phi,
                           const LevelData<FluxBox>& a_velocity,
                           const Real& a_time )
{
   CH_TIMERS("LogRectPhaseBC::apply");
   CH_TIMER("defineBoundaryBoxLayouts", t_define_boundary_box_layouts);
   CH_TIMER("defineInflowDataStorage", t_define_inflow_data_storage);
   CH_TIMER("fillInflowData", t_fill_inflow_data);
   CH_TIMER("setInflowOutflowBC", t_set_inflow_outflow_BC);
   CH_TIMER("setCodimBoundaryValues", t_set_codim_boundary_values);

   const PhaseGeom& geometry( a_species_comp.phaseSpaceGeometry() );
   const MultiBlockCoordSys& coord_sys( *(geometry.coordSysPtr()) );
   
   LevelData<FArrayBox>& Bf( a_species_comp.distributionFunction() );
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
                                            a_species_comp );
      CH_STOP(t_define_inflow_data_storage);

      CH_START(t_fill_inflow_data);
      Vector<std::string> all_bc_type;
      fillInflowData( a_species_comp, a_phi, a_velocity, a_time );
      setAllBcType( m_all_bdry_layouts );
      m_all_bdry_defined = true;
      CH_STOP(t_fill_inflow_data);
   }
      
   CH_START(t_set_inflow_outflow_BC);
   PhaseBCUtils::setInflowOutflowBC( Bf,
                                     m_all_bdry_layouts,
                                     m_all_bdry_data,
                                     m_all_bc_type,
                                     coord_sys,
                                     a_velocity );
   CH_STOP(t_set_inflow_outflow_BC);
   
   // interpolate all other codim boundaries
   CH_START(t_set_codim_boundary_values);
   CodimBC::setCodimBoundaryValues( Bf, coord_sys );
   CH_STOP(t_set_codim_boundary_values);
   
   if (m_logical_sheath_bc) {
      applySheathBC(a_species_comp, a_phi, a_velocity);
   }
      
   if (m_flux_bc) {
      applyFluxBC(a_species_comp, a_phi, a_velocity, a_time);
   }
}

inline
void LogRectPhaseBC::setAllBcType( const BoundaryBoxLayoutPtrVect&  a_bdry_layout )
{
   m_all_bc_type.resize(a_bdry_layout.size());
   for (int i(0); i<a_bdry_layout.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(a_bdry_layout[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      m_all_bc_type[i] = getBcType(dir, side);
   }
}

void LogRectPhaseBC::applySheathBC( KineticSpecies& a_species,
                                    const CFG::LevelData<CFG::FArrayBox>& a_phi,
                                    const LevelData<FluxBox>& a_velocity)
{
   
   // Get coordinate system parameters
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
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
      const Side::LoHiSide& side( bdry_layout.side() );
      
      if (getBcType(dir, side) == "logical_sheath") {

         std::string prefix( m_pp.prefix() );
         prefix += "." + m_bdry_name[i];
         ParmParse pp( prefix.c_str() );

         LogicalSheathBC sheathBC(all_bdry_layouts[i], pp);
         sheathBC.applyBC(a_species, a_velocity, a_phi);
      }
   }
}

void LogRectPhaseBC::applyFluxBC( KineticSpecies& a_species,
                                  const CFG::LevelData<CFG::FArrayBox>& a_phi,
                                  const LevelData<FluxBox>& a_velocity,
                                  const Real& a_time)
{
   // Get coordinate system parameters
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
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
   BoundaryBoxLayoutPtrVect flux_bdry_layouts;
   for (int i(0); i<all_bdry_layouts.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(all_bdry_layouts[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );

      if (getBcType(dir, side) == "flux") {
         
         flux_bdry_layouts.push_back(all_bdry_layouts[i]);
         
         std::string prefix( m_pp.prefix() );
         prefix += "." + m_bdry_name[i];
         ParmParse pp( prefix.c_str() );
         
         FluxBC fluxBC(flux_bdry_layouts, pp);
         fluxBC.applyBC(a_species, a_velocity, a_phi, a_time);
      }
   }
}

void LogRectPhaseBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "LogRectPhaseBC ================================" << std::endl;
      std::cout << "- variable: "  << m_name << "-------------" << std::endl;
      for (int i(0); i<m_inflow_function.size(); i++) {
         std::cout << "  " << m_bdry_name[i] << ": " << std::endl;
         if (m_inflow_function[i].isNull()) std::cout << "no BC specified " << std::endl;
         else m_inflow_function[i]->printParameters();
      }
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << "===============================================" << std::endl;
   }
}


inline
void LogRectPhaseBC::parseParameters( ParmParse& a_pp )
{
   KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   for (int i(0); i<m_inflow_function.size(); i++) {
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );

      //Parsing inflow function
      std::string function_name;
      fpp.query( "function", function_name );
      
      if (function_name.length() == 0) {
         m_inflow_function[i] = RefCountedPtr<KineticFunction>(NULL);
      }
      else {
         m_inflow_function[i] = library->find( function_name );
      }
   
   
      //Parsing BC type
      fpp.query( "type", m_bc_type[i] );
      
      if (m_bc_type[i] == "flux") {
         m_flux_bc = true;
      }

      if (m_bc_type[i] == "logical_sheath") {
         m_logical_sheath_bc = true;
      }
   }
   
   if (m_verbosity) {
      printParameters();
   }
}


#include "NamespaceFooter.H"

