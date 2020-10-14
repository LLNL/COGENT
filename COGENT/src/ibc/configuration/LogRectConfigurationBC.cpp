#include "LogRectConfigurationBC.H"

#include "GridFunctionLibrary.H"
#include "LogRectCoordSys.H"

//#include "FlipGrids.H"
//#include "PoloidalBCF_F.H"

#include "CodimBC.H"
#include "FluidBCUtils.H"


#include "NamespaceHeader.H"

/////// BEGIN INLINE METHODS ///////////////////////////////////////////////////

inline
GridFunction& LogRectConfigurationBC::radialInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[RADIAL_LOWER]);
   }
   return *(m_inflow_function[RADIAL_UPPER]);
}

inline
GridFunction& LogRectConfigurationBC::poloidalInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[POLOIDAL_LOWER]);
   }
   return *(m_inflow_function[POLOIDAL_UPPER]);
}

#if CFG_DIM==3
inline
GridFunction& LogRectConfigurationBC::toroidalInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[TOROIDAL_LOWER]);
   }
   return *(m_inflow_function[TOROIDAL_UPPER]);
}
#endif


inline
GridFunction& LogRectConfigurationBC::inflowFunc( const int&            a_dir,
                                                  const Side::LoHiSide& a_side )
{
   GridFunction* inflow_func;
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
   else {
      MayDay::Error( "LogRectConfigurationBC: BCs not implemented!" );
   }
   return *(inflow_func);
}


inline
std::string LogRectConfigurationBC::radialBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[RADIAL_LOWER];
   }
   return m_bc_type[RADIAL_UPPER];
}

inline
std::string LogRectConfigurationBC::poloidalBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[POLOIDAL_LOWER];
   }
   return m_bc_type[POLOIDAL_UPPER];
}

#if CFG_DIM==3
inline
std::string LogRectConfigurationBC::toroidalBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[TOROIDAL_LOWER];
   }
   return m_bc_type[TOROIDAL_UPPER];
}
#endif

inline
std::string LogRectConfigurationBC::getBcType( const int&             a_dir,
                                               const Side::LoHiSide&  a_side )
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
   else {
      MayDay::Error( "LogRectPhaseBC: BCs not implemented!" );
   }
   return bc_type;
}

/////// END INLINE METHODS /////////////////////////////////////////////////////


LogRectConfigurationBC::LogRectConfigurationBC( const std::string&  a_species_name,
                                                const std::string&  a_variable_name,
                                                const int&          a_verbosity )
   : m_species_name(a_species_name),
     m_variable_name(a_variable_name),
     m_verbosity(a_verbosity)
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

   parseParameters();

}


LogRectConfigurationBC::~LogRectConfigurationBC()
{
}


inline
void LogRectConfigurationBC::fillInflowData( const MagGeom&  a_geometry,
                                             const Real      a_time )
{
   for (int i(0); i<m_all_bdry_layouts.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(m_all_bdry_layouts[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      GridFunction& inflow_func( inflowFunc( dir, side ) );
      LevelData<FArrayBox>& bdry_data( *(m_all_bdry_data[i]) );
      inflow_func.assign( bdry_data, a_geometry, bdry_layout, a_time );
   }
}

void LogRectConfigurationBC::apply( FluidSpecies&  a_species_comp,
                                    const Real&    a_time )
{
   CH_TIME("LogRectConfigurationBC::apply()");
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const LogRectCoordSys& coord_sys( dynamic_cast<const LogRectCoordSys&>( *geometry.getCoordSys()) );

   LevelData<FArrayBox>& u( a_species_comp.cell_var( m_variable_name ) );
   const DisjointBoxLayout& grids( u.disjointBoxLayout() );
   const IntVect& ghost_vect( u.ghostVect() );

   if(m_all_bdry_defined==false) {
      FluidBCUtils::defineBoundaryBoxLayouts( m_all_bdry_layouts,
                                              grids,
                                              coord_sys,
                                              ghost_vect );
   
      FluidBCUtils::defineInflowDataStorage( m_all_bdry_data,
                                             m_all_bdry_layouts,
                                             m_variable_name,
                                             a_species_comp );
   
      fillInflowData( geometry, a_time );
      setAllBcType( m_all_bdry_layouts );    
      m_all_bdry_defined = true;
   }
   setAllBcType( m_all_bdry_layouts );    

   const LevelData<FluxBox>& velocity= a_species_comp.velocity();
   FluidBCUtils::setInflowOutflowBC( u,
                                     m_all_bdry_layouts,
                                     m_all_bdry_data,
                                     m_all_bc_type,
                                     geometry,
                                     velocity );
   
   // interpolate all other codim boundaries
   CodimBC::setCodimBoundaryValues( u, coord_sys );
  
}

void LogRectConfigurationBC::applyFluxBC( const FluidSpecies&  a_species_comp,
                                          LevelData<FluxBox>&  a_dst,
                                    const Real                 a_time )
{
   CH_TIME("LogRectConfigurationBC::applyFluxBC()");
      
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const LogRectCoordSys& coord_sys( dynamic_cast<const LogRectCoordSys&>( *geometry.getCoordSys()) );
  
   const DisjointBoxLayout& grids( a_dst.disjointBoxLayout() );
   //const IntVect& ghost_vect( IntVect::Unit ); // JRA, should this be a_dst.ghostVect() ?
   const IntVect& ghost_vect( a_dst.ghostVect() ); // JRA, should this be a_dst.ghostVect() ?
   
   if(m_all_bdry_layouts_flux_defined==false) {
      FluidBCUtils::defineBoundaryBoxLayouts( m_all_bdry_layouts_flux,
                                              grids,
                                              coord_sys,
                                              ghost_vect );
      m_all_bdry_layouts_flux_defined = true;
   }
   setAllBcType( m_all_bdry_layouts_flux );    

   FluidBCUtils::setFluxBC( a_dst,
                            m_all_bdry_layouts_flux,
                            m_all_bc_type,
                            geometry ) ;
   
}

void LogRectConfigurationBC::applyEdgeBC( const FluidSpecies&      a_species_comp,
                                          LevelData<EdgeDataBox>&  a_dst,
                                    const Real                     a_time )
{
   CH_TIME("LogRectConfigurationBC::applyEdgeBC()");
   
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const LogRectCoordSys& coord_sys( dynamic_cast<const LogRectCoordSys&>( *geometry.getCoordSys()) );
   
   const DisjointBoxLayout& grids( a_dst.disjointBoxLayout() );
   //const IntVect& ghost_vect( IntVect::Unit );
   const IntVect& dst_ghost_vect( a_dst.ghostVect() );
   if(dst_ghost_vect > IntVect::Zero) {   

      if(m_all_bdry_layouts_edge_defined==false) {
         FluidBCUtils::defineBoundaryBoxLayouts( m_all_bdry_layouts_edge,
                                                 grids,
                                                 coord_sys,
                                                 dst_ghost_vect );
         //setAllBcType( m_all_bdry_layouts_edge );    
         m_all_bdry_layouts_edge_defined = true;
      }
      setAllBcType( m_all_bdry_layouts_edge );    

      FluidBCUtils::setEdgeBC( a_dst,
                               m_all_bdry_layouts_edge,
                               m_all_bc_type );

   }

}

void LogRectConfigurationBC::setFluxBC( const FluidSpecies&  a_species_comp,
                                        LevelData<FluxBox>&  a_dst,
                                  const LevelData<FluxBox>&  a_src,
                                  const Real                 a_time )
{
   CH_TIME("LogRectConfigurationBC::setFluxBC()");
      
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const LogRectCoordSys& coord_sys( dynamic_cast<const LogRectCoordSys&>( *geometry.getCoordSys()) );
  
   const DisjointBoxLayout& grids( a_dst.disjointBoxLayout() );
   //const IntVect& ghost_vect( IntVect::Unit ); // JRA, should this be a_dst.ghostVect() ?
   const IntVect& ghost_vect( a_dst.ghostVect() ); // JRA, should this be a_dst.ghostVect() ?
   
   if(m_all_bdry_layouts_flux_defined==false) {
      FluidBCUtils::defineBoundaryBoxLayouts( m_all_bdry_layouts_flux,
                                              grids,
                                              coord_sys,
                                              ghost_vect );
      m_all_bdry_layouts_flux_defined = true;
   }
   setAllBcType( m_all_bdry_layouts_flux );    

   FluidBCUtils::setFluxBC( a_dst,
                            m_all_bdry_layouts_flux,
                            a_src ) ;
   
}

void LogRectConfigurationBC::setEdgeBC( const FluidSpecies&      a_species_comp,
                                        LevelData<EdgeDataBox>&  a_dst,
                                  const LevelData<EdgeDataBox>&  a_src,
                                  const Real                     a_time )
{
   CH_TIME("LogRectConfigurationBC::setEdgeBC()");
   
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const LogRectCoordSys& coord_sys( dynamic_cast<const LogRectCoordSys&>( *geometry.getCoordSys()) );
   
   const DisjointBoxLayout& grids( a_dst.disjointBoxLayout() );
   //const IntVect& ghost_vect( IntVect::Unit );
   const IntVect& ghost_vect( a_dst.ghostVect() );
   
   if(m_all_bdry_layouts_edge_defined==false) {
      FluidBCUtils::defineBoundaryBoxLayouts( m_all_bdry_layouts_edge,
                                              grids,
                                              coord_sys,
                                              ghost_vect );
      //setAllBcType( m_all_bdry_layouts_edge );    
      m_all_bdry_layouts_edge_defined = true;
   }
   setAllBcType( m_all_bdry_layouts_edge );    

   FluidBCUtils::setEdgeBC( a_dst,
                            m_all_bdry_layouts_edge,
                            a_src );
   
}

void LogRectConfigurationBC::applyBC( const FluidSpecies&    a_species_comp,
                                      LevelData<FArrayBox>&  a_dst,
                                const Real                   a_time )
{
   CH_TIME("LogRectConfigurationBC::applyBC()");

   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const LogRectCoordSys& coord_sys( dynamic_cast<const LogRectCoordSys&>( *geometry.getCoordSys()) );

   const DisjointBoxLayout& grids( a_dst.disjointBoxLayout() );
   const IntVect& ghost_vect( a_dst.ghostVect() );

   if(m_all_bdry_defined==false) {
      FluidBCUtils::defineBoundaryBoxLayouts( m_all_bdry_layouts,
                                              grids,
                                              coord_sys,
                                              ghost_vect );
   
      setAllBcType( m_all_bdry_layouts );    
      m_all_bdry_defined = true;
   }
   setAllBcType( m_all_bdry_layouts );    
   
   FluidBCUtils::setBC( a_dst,
                        m_all_bdry_layouts,
                        m_all_bc_type,
                        coord_sys );

}

void LogRectConfigurationBC::setBC( const FluidSpecies&    a_species_comp,
                                    LevelData<FArrayBox>&  a_dst,
                              const LevelData<FArrayBox>&  a_src,
                              const Real                   a_time )
{
   CH_TIME("LogRectConfigurationBC::setBC()");
   
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   
   cout << "LogRectConfigurationBC::setBC() doesn't do anything yet!" << endl;
   
}


inline
void LogRectConfigurationBC::setAllBcType( const BoundaryBoxLayoutPtrVect&  a_bdry_layout )
{
   m_all_bc_type.resize(a_bdry_layout.size());
   for (int i(0); i<a_bdry_layout.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(a_bdry_layout[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      m_all_bc_type[i] = getBcType(dir, side);
   }
}

void LogRectConfigurationBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "LogRectConfigurationBC =============================" << std::endl;
      std::cout << "- variable: "  << m_variable_name << "-------------" << std::endl;
      for (int i(0); i<m_inflow_function.size(); i++) {
         std::cout << "  " << m_bdry_name[i] << ": " << std::endl;
         if ( m_inflow_function[i] ) m_inflow_function[i]->printParameters();
      }
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << "===============================================" << std::endl;
   }
}


inline
void LogRectConfigurationBC::parseParameters()
{
   string prefix = "BC." + m_species_name + "." + m_variable_name;
   ParmParse pp(prefix.c_str());

   GridFunctionLibrary* library = GridFunctionLibrary::getInstance();
   for (int i(0); i<m_inflow_function.size(); i++) {
      std::string prefix( pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );
      std::string function_name;
      if ( fpp.contains("function") ) {
         fpp.query( "function", function_name );
         m_inflow_function[i] = library->find( function_name );
      }
      fpp.query( "type", m_bc_type[i] );
   
   }

   if (m_verbosity) {
      printParameters();
   }
}



#include "NamespaceFooter.H"



