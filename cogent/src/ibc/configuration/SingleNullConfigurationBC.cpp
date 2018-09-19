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
   if ( a_block_type==MCORE ||
       a_block_type==LCORE ||
       a_block_type==RCORE ){
      return m_bc_type[RADIAL_CORE];
   }
   else if ( a_block_type==LPF ||
            a_block_type==RPF ){
      return m_bc_type[RADIAL_PF];
   }
   else if (a_block_type==LSOL  ||
            a_block_type==MCSOL ||
            a_block_type==LCSOL ||
            a_block_type==RCSOL ||
            a_block_type==RSOL  ){
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


SingleNullConfigurationBC::SingleNullConfigurationBC( const std::string& a_name,
                                      ParmParse& a_pp,
                                      const int& a_verbosity )
   : m_name(a_name),
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

   parseParameters( a_pp );
}


SingleNullConfigurationBC::~SingleNullConfigurationBC()
{
}


inline
void SingleNullConfigurationBC::fillInflowData( FluidSpeciesPtrVect& a_bdry_data,
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
         const SingleNullCoordSys& coord_sys(
            dynamic_cast<const SingleNullCoordSys&>( *geometry.getCoordSys()) );
         const DisjointBoxLayout& grids( bdry_layout.disjointBoxLayout() );
         for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
            const Box& interior_box( bdry_layout.interiorBox( dit ) );
            const int block( coord_sys.whichBlock( interior_box ) );
            const int block_type( coord_sys.blockType( block ) );
            GridFunction& inflow_func( radialInflowFunc( side, block_type ) );
            inflow_func.assign( bdry_data.cell_data(), geometry, bdry_layout, a_time );
            a_bc_type[i] = radialBcType(side, block_type);
         }
      }
      else {
         GridFunction& inflow_func( inflowFunc( dir, side ) );
         inflow_func.assign( bdry_data.cell_data(), geometry, bdry_layout, a_time );
         a_bc_type[i] = getBcType(dir, side);
      }
   }
}

void SingleNullConfigurationBC::apply( FluidSpecies& a_species_comp,
                               const Real& a_time )
{
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const SingleNullCoordSys& coord_sys(
      dynamic_cast<const SingleNullCoordSys&>( *geometry.getCoordSys()) );

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

#if 0
   if (m_logical_sheath) {
      applyLogicalSheathBC(a_species_comp, all_bdry_layouts, velocity, a_phi, a_time);
   }
#endif
   
   // interpolate all other codim boundaries
   CodimBC::setCodimBoundaryValues( u, coord_sys );
}



void SingleNullConfigurationBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "SingleNullConfigurationBC =============================" << std::endl;
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
void SingleNullConfigurationBC::parseParameters( ParmParse& a_pp )
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
      fpp.query( "type", m_bc_type[i] );
   }

   a_pp.query("logical_sheath",m_logical_sheath);

   if (m_verbosity) {
      printParameters();
   }
}

#include "NamespaceFooter.H"



