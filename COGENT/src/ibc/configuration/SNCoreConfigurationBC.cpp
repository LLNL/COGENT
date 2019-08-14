#include "SNCoreConfigurationBC.H"

#include "GridFunctionLibrary.H"
#include "SNCoreCoordSys.H"

//#include "FlipGrids.H"
//#include "PoloidalBCF_F.H"

#include "CodimBC.H"
#include "FluidBCUtils.H"


#include "NamespaceHeader.H"

/////// BEGIN INLINE METHODS ///////////////////////////////////////////////////

inline
GridFunction& SNCoreConfigurationBC::radialInflowFunc( const Side::LoHiSide& a_side, const int& a_block_type )
{
   if (a_block_type==SNCoreBlockCoordSys::LCORE||SNCoreBlockCoordSys::RCORE) {
      if (a_side == 0) {
         return *(m_inflow_function[RADIAL_INNER]);
      }
      else if (a_side == 1) {
         return *(m_inflow_function[RADIAL_OUTER]);
      }
   }
   else {
      MayDay::Error("SNCoreConfigurationBC::radialInflowData(): What on earth did you just do?  There should be no way to get this message!");
   }
   
   // This return will never be reached, but the compiler wants to see something returned
   return *(m_inflow_function[RADIAL_INNER]);
}


inline
GridFunction& SNCoreConfigurationBC::inflowFunc( const int& a_dir,
                                           const Side::LoHiSide& a_side )
{
   GridFunction* inflow_func;
   if (a_dir==RADIAL_DIR) {
      MayDay::Error( "SNCoreConfigurationBC::inflowFunc(): Does not handle RADIAL inflow functions!" );
   }
   else {
      MayDay::Error( "SNCoreConfigurationBC::inflowFunc(): Toroidal BCs not implemented!" );
   }
   return *(inflow_func);
}


inline
std::string SNCoreConfigurationBC::radialBcType( const Side::LoHiSide& a_side, const int& a_block_type )
{
   if (a_block_type==SNCoreBlockCoordSys::LCORE||SNCoreBlockCoordSys::RCORE) {
      if (a_side == 0) {
         return m_bc_type[RADIAL_INNER];
      }
      else if (a_side == 1) {
         return m_bc_type[RADIAL_OUTER];
      }
   }
   else {
      MayDay::Error("SNCoreConfigurationBC::radialInflowData(): What on earth did you just do?  There should be no way to get this message!");
   }
   
   // This return will never be reached, but the compiler wants to see something returned
   return m_bc_type[RADIAL_INNER];
}

inline
std::string SNCoreConfigurationBC::getBcType( const int& a_dir,
                                      const Side::LoHiSide& a_side )
{
   std::string bc_type;
   if (a_dir==RADIAL_DIR) {
      MayDay::Error( "SNCoreConfigurationBC::inflowFunc(): Does not handle RADIAL inflow functions!" );
   }
   else {
      MayDay::Error( "SNCoreConfigurationBC::inflowFunc(): Toroidal BCs not implemented!" );
   }
   return bc_type;
}

/////// END INLINE METHODS /////////////////////////////////////////////////////


SNCoreConfigurationBC::SNCoreConfigurationBC( const std::string& a_species_name,
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
   
   m_bdry_name[RADIAL_INNER] = "radial_inner";
   m_bdry_name[RADIAL_OUTER] = "radial_outer";

   parseParameters();
}


SNCoreConfigurationBC::~SNCoreConfigurationBC()
{
}


inline
void SNCoreConfigurationBC::fillInflowData( const  MagGeom&  a_geometry,
                                            const  Real      a_time )
{
   m_all_bc_type.resize( m_all_bdry_layouts.size() );
   for (int i(0); i<m_all_bdry_layouts.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(m_all_bdry_layouts[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      LevelData<FArrayBox>& bdry_data( *(m_all_bdry_data[i]) );

      if (dir==RADIAL_DIR) {
         const SNCoreCoordSys& coord_sys(
            dynamic_cast<const SNCoreCoordSys&>( *a_geometry.getCoordSys()) );
         const DisjointBoxLayout& grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
            const Box& interior_box( bdry_layout.interiorBox( dit ) );
            const int block( coord_sys.whichBlock( interior_box ) );
            const int block_type( coord_sys.blockType( block ) );
            GridFunction& inflow_func( radialInflowFunc( side, block_type ) );
            inflow_func.assign( bdry_data, a_geometry, bdry_layout, a_time );
            m_all_bc_type[i] = radialBcType( side, block_type );
         }
      }
   }
}

void SNCoreConfigurationBC::apply( FluidSpecies&  a_species_comp,
                                   const Real&    a_time )
{
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const SNCoreCoordSys& coord_sys(
      dynamic_cast<const SNCoreCoordSys&>( *geometry.getCoordSys()) );

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

   // interpolate all other codim boundaries
   CodimBC::setCodimBoundaryValues( u, coord_sys );
}

void SNCoreConfigurationBC::applyFluxBC(const FluidSpecies&  a_species_comp,
                                         LevelData<FluxBox>& a_dst,
                                         const LevelData<FluxBox>& a_src,
                                         const Real&    a_time )
{
}

void SNCoreConfigurationBC::applyEdgeBC(const FluidSpecies&  a_species_comp,
                                         LevelData<EdgeDataBox>& a_dst,
                                         const LevelData<EdgeDataBox>& a_src,
                                         const Real&    a_time )
{
}

void SNCoreConfigurationBC::applyBC( const FluidSpecies&          a_species_comp,
                                           LevelData<FArrayBox>&  a_dst,
                                     const Real                   a_time )
{
}





void SNCoreConfigurationBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "SNCoreConfigurationBC =============================" << std::endl;
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
void SNCoreConfigurationBC::parseParameters()
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
      fpp.query( "type", m_bc_type[i]);
   }

   pp.query("logical_sheath",m_logical_sheath);

   if (m_verbosity) {
      printParameters();
   }
}



#include "NamespaceFooter.H"



