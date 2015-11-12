#ifndef _KINETIC_FUNCTION_LIBRARY_H_
#define _KINETIC_FUNCTION_LIBRARY_H_

#include <string>
#include <map>

#include "KineticFunction.H"
#include "RefCountedPtr.H"

#include "NamespaceHeader.H"

/**
 * Generic grid function library.
*/
class KineticFunctionLibrary
{
public:

   static KineticFunctionLibrary* getInstance();

   const RefCountedPtr<KineticFunction> find( const std::string& name );

private:

   KineticFunctionLibrary() : m_verbosity(0) {;}

   ~KineticFunctionLibrary();

   void parseFromInput();

   void printLibrary();

   void registerKineticFunction( const std::string& name,
                                 RefCountedPtr<KineticFunction> function );

   static KineticFunctionLibrary* s_kinetic_function_library_instance;

   typedef std::map<std::string,RefCountedPtr<KineticFunction> > Library;

   std::map<std::string,RefCountedPtr<KineticFunction> > m_library;

   int m_verbosity;
};

#include "NamespaceFooter.H"

#endif
