#include "SNCoreConfigurationBC.H"
#include "SNCoreCoordSys.H"
#include "NamespaceHeader.H"

SNCoreConfigurationBC::SNCoreConfigurationBC( const std::string& a_species_name,
                                              const std::string& a_variable_name,
                                              const int& a_verbosity )
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


SNCoreConfigurationBC::~SNCoreConfigurationBC()
{
}

void
SNCoreConfigurationBC::setNames()
{
   m_bdry_name[RADIAL_INNER] = "radial_inner";
   m_bdry_name[RADIAL_OUTER] = "radial_outer";
}

std::string
SNCoreConfigurationBC::getBCType(const int  a_block_number,
                                 const int  a_dir,
                                 const int  a_side )
{
   CH_assert(a_side == 0 || a_side == 1);
   std::string bc_type;

   switch( a_block_number )
      {
      case SNCoreBlockCoordSys::MCORE:
      case SNCoreBlockCoordSys::LCORE:
      case SNCoreBlockCoordSys::RCORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               bc_type = m_bc_type[RADIAL_INNER];
            }
            else if ( a_side == 1 ) {
               bc_type = m_bc_type[RADIAL_OUTER];
            }
            else {
               MayDay::Error("SNCoreConfigurationBC::getBCType(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCoreConfigurationBC::getBCType(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCoreConfigurationBC::getBCType(): Unrecognized block number");
      }

   return bc_type;
}

RefCountedPtr<GridFunction>
SNCoreConfigurationBC::getBCFunction(const int  a_block_number,
                                     const int  a_dir,
                                     const int  a_side )
{
   CH_assert(a_side == 0 || a_side == 1);
   RefCountedPtr<GridFunction> function;

   switch( a_block_number )
      {
      case SNCoreBlockCoordSys::MCORE:
      case SNCoreBlockCoordSys::LCORE:
      case SNCoreBlockCoordSys::RCORE:
         if ( a_dir == RADIAL_DIR ) {
            if ( a_side == 0 ) {
               function = m_bc_function[RADIAL_INNER];
            }
            else if ( a_side == 1 ) {
               function = m_bc_function[RADIAL_OUTER];
            }
            else {
               MayDay::Error("SNCoreConfigurationBC::getBCFunction(): Invalid side argument");
            }
         }
         else {
            MayDay::Error("SNCoreConfigurationBC::getBCFunction(): Invalid direction argument");
         }
         break;
      default:
         MayDay::Error("SNCoreConfigurationBC::getBCFunction(): Unrecognized block number");
      }

   return function;
}

void SNCoreConfigurationBC::parseParameters(ParmParse& a_pp )
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

void SNCoreConfigurationBC::printParameters() const
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



