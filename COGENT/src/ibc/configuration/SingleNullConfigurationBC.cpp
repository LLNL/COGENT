#include "SingleNullConfigurationBC.H"

#include "GridFunctionLibrary.H"
#include "SingleNullCoordSys.H"

//#include "FlipGrids.H"
//#include "PoloidalBCF_F.H"

#include "CodimBC.H"
#include "FluidBCUtils.H"


#include "NamespaceHeader.H"

/////// BEGIN INLINE METHODS ///////////////////////////////////////////////////

inline
GridFunction& SingleNullConfigurationBC::radialInflowFunc( const Side::LoHiSide& a_side, const int& a_block_type )
{
   if ( a_block_type==SingleNullBlockCoordSys::MCORE ||
        a_block_type==SingleNullBlockCoordSys::LCORE ||
        a_block_type==SingleNullBlockCoordSys::RCORE ){
      return *(m_inflow_function[RADIAL_CORE]);
   }
   else if ( a_block_type==SingleNullBlockCoordSys::LPF ||
             a_block_type==SingleNullBlockCoordSys::RPF ){
      return *(m_inflow_function[RADIAL_PF]);
   }
   else if (a_block_type==SingleNullBlockCoordSys::LSOL  ||
            a_block_type==SingleNullBlockCoordSys::MCSOL ||
            a_block_type==SingleNullBlockCoordSys::LCSOL ||
            a_block_type==SingleNullBlockCoordSys::RCSOL ||
            a_block_type==SingleNullBlockCoordSys::RSOL  ){
      return *(m_inflow_function[RADIAL_SOL]);
   }
   else {
      MayDay::Error("SingleNullConfigurationBC::radialInflowData(): What on earth did you just do?  There should be no way to get this message!");
   }

   // This return will never be reached, but the compiler wants to see something returned
   return *(m_inflow_function[RADIAL_CORE]);
}


inline
GridFunction& SingleNullConfigurationBC::poloidalInflowFunc( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return *(m_inflow_function[POLOIDAL_OUTER_DIV]);
   }
   return *(m_inflow_function[POLOIDAL_INNER_DIV]);
}


inline
GridFunction& SingleNullConfigurationBC::inflowFunc( const int& a_dir,
                                             const Side::LoHiSide& a_side )
{
   GridFunction* inflow_func;
   if (a_dir==RADIAL_DIR) {
      MayDay::Error( "SingleNullConfigurationBC::inflowFunc(): Does not handle RADIAL inflow functions!" );
   }
   else if (a_dir==POLOIDAL_DIR) {
      inflow_func = &poloidalInflowFunc( a_side );
   }
   else {
      MayDay::Error( "SingleNullConfigurationBC::inflowFunc(): BCs not implemented!" );
   }
   return *(inflow_func);
}

inline
std::string SingleNullConfigurationBC::radialBcType( const Side::LoHiSide& a_side, const int& a_block_type )
{
   if ( a_block_type==SingleNullBlockCoordSys::MCORE ||
        a_block_type==SingleNullBlockCoordSys::LCORE ||
        a_block_type==SingleNullBlockCoordSys::RCORE ){
      return m_bc_type[RADIAL_CORE];
   }
   else if ( a_block_type==SingleNullBlockCoordSys::LPF ||
             a_block_type==SingleNullBlockCoordSys::RPF ){
      return m_bc_type[RADIAL_PF];
   }
   else if (a_block_type==SingleNullBlockCoordSys::LSOL  ||
            a_block_type==SingleNullBlockCoordSys::MCSOL ||
            a_block_type==SingleNullBlockCoordSys::LCSOL ||
            a_block_type==SingleNullBlockCoordSys::RCSOL ||
            a_block_type==SingleNullBlockCoordSys::RSOL  ){
      return m_bc_type[RADIAL_SOL];
   }
   else {
      MayDay::Error("SingleNullConfigurationBC::radialInflowData(): What on earth did you just do?  There should be no way to get this message!");
   }
   
   // This return will never be reached, but the compiler wants to see something returned
   return m_bc_type[RADIAL_CORE];
}

inline
std::string SingleNullConfigurationBC::poloidalBcType( const Side::LoHiSide& a_side )
{
   if (a_side==Side::Lo) {
      return m_bc_type[POLOIDAL_INNER_DIV];
   }
   return m_bc_type[POLOIDAL_OUTER_DIV];
}

inline
std::string SingleNullConfigurationBC::getBcType( const int& a_dir,
                                         const Side::LoHiSide& a_side )
{
   std::string bc_type;
   if (a_dir==RADIAL_DIR) {
      MayDay::Error( "SingleNullConfigurationBC::inflowFunc(): Does not handle RADIAL inflow functions!" );
   }
   else if (a_dir==POLOIDAL_DIR) {
      bc_type = poloidalBcType( a_side );
   }
   else {
      MayDay::Error( "SingleNullConfigurationBC::inflowFunc(): Toroidal BCs not implemented!" );
   }
   return bc_type;
}



/////// END INLINE METHODS /////////////////////////////////////////////////////


SingleNullConfigurationBC::SingleNullConfigurationBC( const std::string& a_species_name,
                                                      const std::string& a_variable_name,
                                                      const int& a_verbosity )
   : m_species_name(a_species_name),
     m_variable_name(a_variable_name),
     m_verbosity(a_verbosity),
     m_logical_sheath(false)
{
   m_inflow_function.resize( NUM_INFLOW );
   m_bc_type.resize( NUM_INFLOW );
   m_bdry_name.resize( NUM_INFLOW );
   
   m_bdry_name[RADIAL_CORE] = "radial_core";
   m_bdry_name[RADIAL_SOL] = "radial_sol";
   m_bdry_name[RADIAL_PF] = "radial_pf";
   m_bdry_name[POLOIDAL_INNER_DIV] = "poloidal_inner_div";
   m_bdry_name[POLOIDAL_OUTER_DIV] = "poloidal_outer_div";

   parseParameters();
}


SingleNullConfigurationBC::~SingleNullConfigurationBC()
{
}


inline
void SingleNullConfigurationBC::fillInflowData( const MagGeom&  a_geometry,
                                                const Real      a_time )
{
   m_all_bc_type.resize(m_all_bdry_layouts.size());
   for (int i(0); i<m_all_bdry_layouts.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(m_all_bdry_layouts[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      LevelData<FArrayBox>& bdry_data( *(m_all_bdry_data[i]) );

      if (dir==RADIAL_DIR) {
         const SingleNullCoordSys& coord_sys(
            dynamic_cast<const SingleNullCoordSys&>( *a_geometry.getCoordSys()) );
         const DisjointBoxLayout& grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
            const Box& interior_box( bdry_layout.interiorBox( dit ) );
            const int block( coord_sys.whichBlock( interior_box ) );
            const int block_type( coord_sys.blockType( block ) );
            GridFunction& inflow_func( radialInflowFunc( side, block_type ) );
            inflow_func.assign( bdry_data, a_geometry, bdry_layout, a_time );
            m_all_bc_type[i] = radialBcType(side, block_type);
         }
      }
      else {
         GridFunction& inflow_func( inflowFunc( dir, side ) );
         inflow_func.assign( bdry_data, a_geometry, bdry_layout, a_time );
         m_all_bc_type[i] = getBcType(dir, side);
      }
   }
}

void SingleNullConfigurationBC::apply( FluidSpecies&  a_species_comp,
                                       const Real&    a_time )
{
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const SingleNullCoordSys& coord_sys(
      dynamic_cast<const SingleNullCoordSys&>( *geometry.getCoordSys()) );

   LevelData<FArrayBox>& u( a_species_comp.cell_var(a_species_comp.cell_var_component(m_variable_name)) );
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
      m_all_bdry_defined = true;
   }

   const LevelData<FluxBox>& velocity= a_species_comp.velocity();
   FluidBCUtils::setInflowOutflowBC( u,
                                     m_all_bdry_layouts,
                                     m_all_bdry_data,
                                     m_all_bc_type,
                                     geometry,
                                     velocity );

#if 0
   if (m_logical_sheath) {
      applyLogicalSheathBC(a_species_comp, all_bdry_layouts, velocity, a_phi, a_time);
   }
#endif
   
   // interpolate all other codim boundaries
   CodimBC::setCodimBoundaryValues( u, coord_sys );
}

void SingleNullConfigurationBC::applyFluxBC(const FluidSpecies&  a_species_comp,
                                         LevelData<FluxBox>& a_dst,
                                         const LevelData<FluxBox>& a_src,
                                         const Real&    a_time )
{
}

void SingleNullConfigurationBC::applyEdgeBC(const FluidSpecies&  a_species_comp,
                                         LevelData<EdgeDataBox>& a_dst,
                                         const LevelData<EdgeDataBox>& a_src,
                                         const Real&    a_time )
{
}

void SingleNullConfigurationBC::applyBC( const FluidSpecies&          a_species_comp,
                                               LevelData<FArrayBox>&  a_dst,
                                         const Real                   a_time )
{
}


void SingleNullConfigurationBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "SingleNullConfigurationBC =============================" << std::endl;
      std::cout << "- variable: "  << m_variable_name << "-------------" << std::endl;
      for (int i(0); i<m_inflow_function.size(); i++) {
         std::cout << "  " << m_bdry_name[i] << ": " << std::endl;
         if ( m_inflow_function[i] ) m_inflow_function[i]->printParameters();
      }
      std::cout << "  logical_sheath  =  " << m_logical_sheath << std::endl;
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << "===============================================" << std::endl;
   }
}


inline
void SingleNullConfigurationBC::parseParameters()
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

   pp.query("logical_sheath",m_logical_sheath);

   if (m_verbosity) {
      printParameters();
   }
}

#include "NamespaceFooter.H"



