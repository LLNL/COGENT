#include "KineticFunctionLibrary.H"

#include "ParmParse.H"
#include "KineticFunctionFactory.H"

#include <sstream>

#include "NamespaceHeader.H"

KineticFunctionLibrary*
KineticFunctionLibrary::s_kinetic_function_library_instance( NULL );


KineticFunctionLibrary* KineticFunctionLibrary::getInstance()
{
   if (!s_kinetic_function_library_instance) {
      s_kinetic_function_library_instance = new KineticFunctionLibrary();
      s_kinetic_function_library_instance->parseFromInput();
   }
   return s_kinetic_function_library_instance;
}


const RefCountedPtr<KineticFunction>
KineticFunctionLibrary::find( const std::string& a_name )
{
   Library::const_iterator it = m_library.find( a_name );
   if ( it == m_library.end() ) {
      stringstream msg( "KineticFunction ", ios_base::out|ios_base::ate );
      msg << a_name << " not found in library!";
      MayDay::Error( msg.str().c_str() );
   }
   return it->second;
}


KineticFunctionLibrary::~KineticFunctionLibrary() {
   if (!s_kinetic_function_library_instance) {
      delete s_kinetic_function_library_instance;
   }
}


void
KineticFunctionLibrary::parseFromInput()
{
   ParmParse gfpp( "kinetic_function_library" );

   gfpp.query( "verbosity", m_verbosity );

   int n(0);
   gfpp.query( "number", n );
   if (n==0) {
      MayDay::Warning( "KineticFunctionLibrary: No kinetic functions found in input file" );
      return;
   }

   std::vector<std::string> name( n );
   gfpp.queryarr( "list", name, 0, n );
   for (int i(0); i<n; i++) {
      std::string prefix( gfpp.prefix() );
      prefix += "." + name[i];
      ParmParse fpp( prefix.c_str() );
      KineticFunctionFactory factory;
      registerKineticFunction( name[i], factory.create( fpp, 0 ) );
   }

   if (m_verbosity) {
      printLibrary();
   }
}


void
KineticFunctionLibrary::printLibrary()
{
   if (procID()==0) {
      std::cout << std::endl;
      std::cout << "Kinetic Function Library =========================" << std::endl;
      std::cout << std::endl;

      Library::const_iterator it = m_library.begin();
      for (; it!=m_library.end(); it++) {
         std::cout << "\""<< it->first << "\":" <<std::endl;
         (it->second)->printParameters();
      }
      std::cout << "==================================================" << std::endl;
   }
}


void
KineticFunctionLibrary::registerKineticFunction(
   const std::string& a_name,
   RefCountedPtr<KineticFunction> a_function )
{
   Library::const_iterator it = m_library.find( a_name );
   if ( it == m_library.end() ) {
      m_library.insert( Library::value_type( a_name, a_function ) );
   }
   else {
      stringstream msg( "KineticFunction ", ios_base::out|ios_base::ate );
      msg << a_name << " already exists in library!";
      MayDay::Error( msg.str().c_str() );
   }
}

#include "NamespaceFooter.H"
