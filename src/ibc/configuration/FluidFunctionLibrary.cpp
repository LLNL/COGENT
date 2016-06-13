#include "FluidFunctionLibrary.H"

#include "ParmParse.H"
#include "FluidFunctionFactory.H"

#include <sstream>

#include "NamespaceHeader.H"

FluidFunctionLibrary*
FluidFunctionLibrary::s_fluid_function_library_instance( NULL );


FluidFunctionLibrary* FluidFunctionLibrary::getInstance()
{
   if (!s_fluid_function_library_instance) {
      s_fluid_function_library_instance = new FluidFunctionLibrary();
      s_fluid_function_library_instance->parseFromInput();
   }
   return s_fluid_function_library_instance;
}


const RefCountedPtr<FluidFunction>
FluidFunctionLibrary::find( const std::string& a_name )
{
   Library::const_iterator it = m_library.find( a_name );
   if ( it == m_library.end() ) {
      stringstream msg( "FluidFunction ", ios_base::out|ios_base::ate );
      msg << a_name << " not found in library!";
      MayDay::Error( msg.str().c_str() );
   }
   return it->second;
}


FluidFunctionLibrary::~FluidFunctionLibrary() {
   if (!s_fluid_function_library_instance) {
      delete s_fluid_function_library_instance;
   }
}


void
FluidFunctionLibrary::parseFromInput()
{
   ParmParse gfpp( "fluid_function_library" );

   gfpp.query( "verbosity", m_verbosity );

   int n(0);
   gfpp.query( "number", n );
   if (n==0) {
      MayDay::Warning( "FluidFunctionLibrary: No fluid functions found in input file" );
      return;
   }

   std::vector<std::string> name( n );
   gfpp.queryarr( "list", name, 0, n );
   for (int i(0); i<n; i++) {
      std::string prefix( gfpp.prefix() );
      prefix += "." + name[i];
      ParmParse fpp( prefix.c_str() );
      FluidFunctionFactory factory;
      registerFluidFunction( name[i], factory.create( fpp, 0 ) );
   }

   if (m_verbosity) {
      printLibrary();
   }
}


void
FluidFunctionLibrary::printLibrary()
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "Fluid Function Library =========================" << std::endl;
      std::cout << std::endl;

      Library::const_iterator it = m_library.begin();
      for (; it!=m_library.end(); it++) {
         std::cout << "\""<< it->first << "\":" <<std::endl;
         (it->second)->printParameters();
      }
      std::cout << "===============================================" << std::endl;
   }
}


void
FluidFunctionLibrary::registerFluidFunction(
   const std::string& a_name,
   RefCountedPtr<FluidFunction> a_function )
{
   Library::const_iterator it = m_library.find( a_name );
   if ( it == m_library.end() ) {
      m_library.insert( Library::value_type( a_name, a_function ) );
   }
   else {
      stringstream msg( "FluidFunction ", ios_base::out|ios_base::ate );
      msg << a_name << " already exists in library!";
      MayDay::Error( msg.str().c_str() );
   }
}

#include "NamespaceFooter.H"
