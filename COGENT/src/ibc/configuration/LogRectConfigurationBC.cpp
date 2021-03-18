#include "LogRectConfigurationBC.H"
#include "NamespaceHeader.H"

LogRectConfigurationBC::LogRectConfigurationBC( const std::string&  a_species_name,
                                                const std::string&  a_variable_name,
                                                const int&          a_verbosity )
   : FluidVarBC(a_species_name,
                a_variable_name,
                a_verbosity,
                NUM_BOUNDARIES),
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


LogRectConfigurationBC::~LogRectConfigurationBC()
{
}


void
LogRectConfigurationBC::setNames()
{
   m_bdry_name[RADIAL_LOWER] = "radial_lower";
   m_bdry_name[RADIAL_UPPER] = "radial_upper";
   m_bdry_name[POLOIDAL_LOWER] = "poloidal_lower";
   m_bdry_name[POLOIDAL_UPPER] = "poloidal_upper";
#if CFG_DIM==3
   m_bdry_name[TOROIDAL_LOWER] = "toroidal_lower";
   m_bdry_name[TOROIDAL_UPPER] = "toroidal_upper";
#endif
}

std::string
LogRectConfigurationBC::getBCType(const int  a_block_number,
                                  const int  a_dir,
                                  const int  a_side )
{
   std::string bc_type;

   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         bc_type = m_bc_type[RADIAL_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_type = m_bc_type[RADIAL_UPPER];
      }
      else {
         MayDay::Error("LogRectConfigurationBC::getBCType(): Invalid side argument");
      }
   }
   else if ( a_dir == POLOIDAL_DIR ) {
      if ( a_side == 0 ) {
         bc_type = m_bc_type[POLOIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_type = m_bc_type[POLOIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectConfigurationBC::getBCType(): Invalid side argument");
      }
   }
#if CFG_DIM==3
   else if ( a_dir == TOROIDAL_DIR ) {
      if ( a_side == 0 ) {
         bc_type = m_bc_type[TOROIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         bc_type = m_bc_type[TOROIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectConfigurationBC::getBCType(): Invalid side argument");
      }
   }
#endif
   else {
      MayDay::Error("LogRectConfigurationBC::getBCType(): Invalid direction argument");
   }

   return bc_type;
}

RefCountedPtr<GridFunction>
LogRectConfigurationBC::getBCFunction(const int  a_block_number,
                                      const int  a_dir,
                                      const int  a_side )
{
   RefCountedPtr<GridFunction> function;

   if ( a_dir == RADIAL_DIR ) {
      if ( a_side == 0 ) {
         function = m_bc_function[RADIAL_LOWER];
      }
      else if ( a_side == 1 ) {
         function = m_bc_function[RADIAL_UPPER];
      }
      else {
         MayDay::Error("LogRectConfigurationBC::getBCFunction(): Invalid side argument");
      }
   }
   else if ( a_dir == POLOIDAL_DIR ) {
      if ( a_side == 0 ) {
         function = m_bc_function[POLOIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         function = m_bc_function[POLOIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectConfigurationBC::getBCFunction(): Invalid side argument");
      }
   }
#if CFG_DIM==3
   else if ( a_dir == TOROIDAL_DIR ) {
      if ( a_side == 0 ) {
         function = m_bc_function[TOROIDAL_LOWER];
      }
      else if ( a_side == 1 ) {
         function = m_bc_function[TOROIDAL_UPPER];
      }
      else {
         MayDay::Error("LogRectConfigurationBC::getBCFunction(): Invalid side argument");
      }
   }
#endif
   else {
      MayDay::Error("LogRectConfigurationBC::getBCFunction(): Invalid direction argument");
   }

   return function;
}

void LogRectConfigurationBC::parseParameters(ParmParse& a_pp )
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

void LogRectConfigurationBC::printParameters() const
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



