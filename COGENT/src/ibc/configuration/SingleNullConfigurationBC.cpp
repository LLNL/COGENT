#include "SingleNullConfigurationBC.H"
#include "SingleNullCoordSys.H"
#include "NamespaceHeader.H"

SingleNullConfigurationBC::SingleNullConfigurationBC( const std::string& a_species_name,
                                                      const std::string& a_variable_name,
                                                      const int&         a_poloidal_blocks_per_sector,
                                                      const int&         a_verbosity )
   : FluidVarBC(a_species_name,
                a_variable_name,
                a_verbosity,
                NUM_BOUNDARIES),
     m_poloidal_blocks_per_sector(a_poloidal_blocks_per_sector),
     m_verbosity(a_verbosity)
{

   setNames();
   
   string pp_prefix = "BC." + m_species_name + "." + m_variable_name;
   ParmParse pp(pp_prefix.c_str());
   parseParameters(pp);
   
   if (m_verbosity) {
      printParameters();
   }
   
}


SingleNullConfigurationBC::~SingleNullConfigurationBC()
{
}

void
SingleNullConfigurationBC::setNames()
{
   m_bdry_name[RADIAL_CORE] = "radial_core";
   m_bdry_name[RADIAL_SOL] = "radial_sol";
   m_bdry_name[RADIAL_PF] = "radial_pf";
   m_bdry_name[POLOIDAL_INNER_DIV] = "poloidal_inner_div";
   m_bdry_name[POLOIDAL_OUTER_DIV] = "poloidal_outer_div";
}

std::string
SingleNullConfigurationBC::getBCType(const int  a_block_number,
                                     const int  a_dir,
                                     const int  a_side )
{
   CH_assert(a_side == 0 || a_side == 1);
   std::string bc_type;

   switch( a_block_number % m_poloidal_blocks_per_sector )
      {
      case SingleNullBlockCoordSys::MCORE:
      case SingleNullBlockCoordSys::LCORE:
      case SingleNullBlockCoordSys::RCORE:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            bc_type = m_bc_type[RADIAL_CORE];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            bc_type = m_bc_type[TOROIDAL_CORE];
         }
#endif
         else {
            MayDay::Error("SingleNullConfigurationBC::getBCType(): Invalid argument");
         }
         break;
      case SingleNullBlockCoordSys::MCSOL:
      case SingleNullBlockCoordSys::LCSOL:
      case SingleNullBlockCoordSys::RCSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            bc_type = m_bc_type[RADIAL_SOL];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            bc_type = m_bc_type[TOROIDAL_SOL];
         }
#endif
         else {
            MayDay::Error("SingleNullConfigurationBC::getBCType(): Invalid argument");
         }
         break;
      case SingleNullBlockCoordSys::LSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            bc_type = m_bc_type[RADIAL_SOL];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            bc_type = m_bc_type[TOROIDAL_SOL];
         }
#endif
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            bc_type = m_bc_type[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullConfigurationBC::getBCType(): Invalid argument");
         }
         break;
      case SingleNullBlockCoordSys::RSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            bc_type = m_bc_type[RADIAL_SOL];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            bc_type = m_bc_type[TOROIDAL_SOL];
         }
#endif
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            bc_type = m_bc_type[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullConfigurationBC::getBCType(): Invalid argument");
         }
         break;
      case SingleNullBlockCoordSys::LPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            bc_type = m_bc_type[RADIAL_PF];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            bc_type = m_bc_type[TOROIDAL_PF];
         }
#endif
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            bc_type = m_bc_type[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullConfigurationBC::getBCType(): Invalid argument");
         }
         break;
      case SingleNullBlockCoordSys::RPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            bc_type = m_bc_type[RADIAL_PF];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            bc_type = m_bc_type[TOROIDAL_PF];
         }
#endif
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            bc_type = m_bc_type[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullConfigurationBC::getBCType(): Invalid argument");
         }
         break;
      default:
         MayDay::Error("SingleNullConfigurationBC::getBCType(): Unrecognized block number");
      }

   return bc_type;
}

RefCountedPtr<GridFunction>
SingleNullConfigurationBC::getBCFunction(const int  a_block_number,
                                         const int  a_dir,
                                         const int  a_side )
{
   CH_assert(a_side == 0 || a_side == 1);
   RefCountedPtr<GridFunction> function;

   switch( a_block_number % m_poloidal_blocks_per_sector )
      {
      case SingleNullBlockCoordSys::MCORE:
      case SingleNullBlockCoordSys::LCORE:
      case SingleNullBlockCoordSys::RCORE:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            function = m_bc_function[RADIAL_CORE];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            function = m_bc_function[TOROIDAL_CORE];
         }
#endif
         else {
            MayDay::Error("SingleNullConfigurationBC::getBCFunction(): Invalid argument");
         }
         break;
      case SingleNullBlockCoordSys::MCSOL:
      case SingleNullBlockCoordSys::LCSOL:
      case SingleNullBlockCoordSys::RCSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            function = m_bc_function[RADIAL_SOL];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            function = m_bc_function[TOROIDAL_SOL];
         }
#endif
         break;
      case SingleNullBlockCoordSys::LSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            function = m_bc_function[RADIAL_SOL];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            function = m_bc_function[TOROIDAL_SOL];
         }
#endif
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            function = m_bc_function[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullConfigurationBC::getBCFunction(): Invalid argument");
         }
         break;
      case SingleNullBlockCoordSys::RSOL:
         if (a_dir == RADIAL_DIR && a_side == 1) {
            function = m_bc_function[RADIAL_SOL];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            function = m_bc_function[TOROIDAL_SOL];
         }
#endif
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            function = m_bc_function[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullConfigurationBC::getBCFunction(): Invalid argument");
         }
         break;
      case SingleNullBlockCoordSys::LPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            function = m_bc_function[RADIAL_PF];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            function = m_bc_function[TOROIDAL_PF];
         }
#endif
         else if (a_dir == POLOIDAL_DIR && a_side == 1) {
            function = m_bc_function[POLOIDAL_INNER_DIV];
         }
         else {
            MayDay::Error("SingleNullConfigurationBC::getBCFunction(): Invalid argument");
         }
         break;
      case SingleNullBlockCoordSys::RPF:
         if (a_dir == RADIAL_DIR && a_side == 0) {
            function = m_bc_function[RADIAL_PF];
         }
#if CFG_DIM==3
         else if (a_dir == TOROIDAL_DIR) {
            function = m_bc_function[TOROIDAL_PF];
         }
#endif
         else if (a_dir == POLOIDAL_DIR && a_side == 0) {
            function = m_bc_function[POLOIDAL_OUTER_DIV];
         }
         else {
            MayDay::Error("SingleNullConfigurationBC::getBCFunction(): Invalid argument");
         }
         break;
      default:
         MayDay::Error("SingleNullConfigurationBC::getBCFunction(): Unrecognized block number");
      }

   return function;
}


void SingleNullConfigurationBC::parseParameters(ParmParse& a_pp )
{

   GridFunctionLibrary* library = GridFunctionLibrary::getInstance();
   for (int i(0); i<m_bc_function.size(); i++) {
      std::string prefix( a_pp.prefix() );
      prefix += "." + m_bdry_name[i];
      ParmParse fpp( prefix.c_str() );
      std::string function_name;
      if ( fpp.contains("function") ) {
         fpp.query( "function", function_name );
         m_bc_function[i] = library->find( function_name );
      }
      fpp.query( "type", m_bc_type[i] );
      
      if (m_bc_type[i] == "recycling") {
         m_recycling_bc = true;
      }
   }

   a_pp.query( "recycling_coefficient", m_recycling_coefficient );

}

void SingleNullConfigurationBC::printParameters() const
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "SingleNullConfigurationBC =============================" << std::endl;
      std::cout << "- variable: "  << m_variable_name << "-------------" << std::endl;
      for (int i(0); i<m_bc_function.size(); i++) {
         std::cout << "  " << m_bdry_name[i] << ": " << std::endl;
         if ( m_bc_function[i] ) m_bc_function[i]->printParameters();
      }
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << "===============================================" << std::endl;
   }
}

#include "NamespaceFooter.H"



