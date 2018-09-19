#include "GridFunctionLibrary.H"

#include "ParmParse.H"
#include "GridFunctionFactory.H"

#include <sstream>

#include "NamespaceHeader.H"

GridFunctionLibrary*
GridFunctionLibrary::s_grid_function_library_instance( NULL );


GridFunctionLibrary* GridFunctionLibrary::getInstance()
{
   if (!s_grid_function_library_instance) {
      s_grid_function_library_instance = new GridFunctionLibrary();
      s_grid_function_library_instance->parseFromInput();
   }
   return s_grid_function_library_instance;
}


const RefCountedPtr<GridFunction>
GridFunctionLibrary::find( const std::string& a_name )
{
   Library::const_iterator it = m_library.find( a_name );
   if ( it == m_library.end() ) {
      stringstream msg( "GridFunction ", ios_base::out|ios_base::ate );
      msg << a_name << " not found in library!";
      MayDay::Error( msg.str().c_str() );
   }
   return it->second;
}


GridFunctionLibrary::~GridFunctionLibrary() {
   if (!s_grid_function_library_instance) {
      delete s_grid_function_library_instance;
   }
}


void
GridFunctionLibrary::parseFromInput()
{
   ParmParse gfpp( "grid_function_library" );

   gfpp.query( "verbosity", m_verbosity );

   int n(0);
   gfpp.query( "number", n );
   if (n==0) {
      MayDay::Warning( "GridFunctionLibrary: No grid functions found in input file" );
      return;
   }

   std::vector<std::string> name( n );
   gfpp.queryarr( "list", name, 0, n );
   for (int i(0); i<n; i++) {
      std::string prefix( gfpp.prefix() );
      prefix += "." + name[i];
      ParmParse fpp( prefix.c_str() );
      GridFunctionFactory factory;
      registerGridFunction( name[i], factory.create( fpp, 0 ) );
   }

   if (m_verbosity) {
      printLibrary();
   }
}


void
GridFunctionLibrary::printLibrary()
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "Grid Function Library =========================" << std::endl;
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
GridFunctionLibrary::registerGridFunction(
   const std::string& a_name,
   RefCountedPtr<GridFunction> a_function )
{
   Library::const_iterator it = m_library.find( a_name );
   if ( it == m_library.end() ) {
      m_library.insert( Library::value_type( a_name, a_function ) );
   }
   else {
      stringstream msg( "GridFunction ", ios_base::out|ios_base::ate );
      msg << a_name << " already exists in library!";
      MayDay::Error( msg.str().c_str() );
   }
}

#include "NamespaceFooter.H"
