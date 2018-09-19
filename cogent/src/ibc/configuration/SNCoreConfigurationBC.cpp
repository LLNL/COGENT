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
   if (a_block_type==L_CORE||R_CORE) {
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
   if (a_block_type==L_CORE||R_CORE) {
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


SNCoreConfigurationBC::SNCoreConfigurationBC( const std::string& a_name,
                                      ParmParse& a_pp,
                                      const int& a_verbosity )
   : m_name(a_name),
     m_verbosity(a_verbosity),
     m_logical_sheath(false)
{
   m_inflow_function.resize( NUM_INFLOW );
   m_bc_type.resize( NUM_INFLOW );
   m_bdry_name.resize( NUM_INFLOW );
   
   m_bdry_name[RADIAL_INNER] = "radial_inner";
   m_bdry_name[RADIAL_OUTER] = "radial_outer";

   parseParameters( a_pp );
}


SNCoreConfigurationBC::~SNCoreConfigurationBC()
{
}


inline
void SNCoreConfigurationBC::fillInflowData( FluidSpeciesPtrVect& a_bdry_data,
                                    Vector<std::string>& a_bc_type,
                                    const BoundaryBoxLayoutPtrVect& a_bdry_layout,
                                    const Real& a_time )
{
   a_bc_type.resize(a_bdry_layout.size());
   for (int i(0); i<a_bdry_layout.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(a_bdry_layout[i]) );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      FluidSpecies& bdry_data( static_cast<FluidSpecies&>(*(a_bdry_data[i])) );
      const MagGeom& geometry( bdry_data.configurationSpaceGeometry() );

      if (dir==RADIAL_DIR) {
         const SNCoreCoordSys& coord_sys(
            dynamic_cast<const SNCoreCoordSys&>( *geometry.getCoordSys()) );
         const DisjointBoxLayout& grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
            const Box& interior_box( bdry_layout.interiorBox( dit ) );
            const int block( coord_sys.whichBlock( interior_box ) );
            const int block_type( coord_sys.blockType( block ) );
            GridFunction& inflow_func( radialInflowFunc( side, block_type ) );
            inflow_func.assign( bdry_data.cell_data(), geometry, bdry_layout, a_time );
            a_bc_type[i] = radialBcType( side, block_type );
         }
      }
   }
}

void SNCoreConfigurationBC::apply( FluidSpecies& a_species_comp,
                               const Real& a_time )
{
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const SNCoreCoordSys& coord_sys(
      dynamic_cast<const SNCoreCoordSys&>( *geometry.getCoordSys()) );

   LevelData<FArrayBox>& u( a_species_comp.cell_data() );
   const DisjointBoxLayout& grids( u.disjointBoxLayout() );
   const IntVect& ghost_vect( u.ghostVect() );

   const LevelData<FluxBox>& velocity= a_species_comp.velocity();

   BoundaryBoxLayoutPtrVect all_bdry_layouts;
   FluidBCUtils::defineBoundaryBoxLayouts( all_bdry_layouts,
                                           grids,
                                           coord_sys,
                                           ghost_vect );

   FluidSpeciesPtrVect all_bdry_data;
   FluidBCUtils::defineInflowDataStorage( all_bdry_data,
                                          all_bdry_layouts,
                                          a_species_comp );
   Vector<std::string> all_bc_type;
   fillInflowData( all_bdry_data, all_bc_type, all_bdry_layouts, a_time );

  
   FluidBCUtils::setInflowOutflowBC( u,
                                     all_bdry_layouts,
                                     all_bdry_data,
                                     all_bc_type,
                                     coord_sys,
                                     velocity );

   // interpolate all other codim boundaries
   CodimBC::setCodimBoundaryValues( u, coord_sys );
}



void SNCoreConfigurationBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "SNCoreConfigurationBC =============================" << std::endl;
      std::cout << "- variable: "  << m_name << "-------------" << std::endl;
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
void SNCoreConfigurationBC::parseParameters( ParmParse& a_pp )
{
   GridFunctionLibrary* library = GridFunctionLibrary::getInstance();
   for (int i(0); i<m_inflow_function.size(); i++) {
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );
      std::string function_name;
      if ( fpp.contains("function") ) {
         fpp.query( "function", function_name );
         m_inflow_function[i] = library->find( function_name );
      }
      fpp.query( "type", m_bc_type[i]);
   }

   a_pp.query("logical_sheath",m_logical_sheath);

   if (m_verbosity) {
      printParameters();
   }
}



#include "NamespaceFooter.H"



