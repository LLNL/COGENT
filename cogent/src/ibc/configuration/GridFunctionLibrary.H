#ifndef _GRID_FUNCTION_LIBRARY_H_
#define _GRID_FUNCTION_LIBRARY_H_

#include <string>
#include <map>

#include "GridFunction.H"
#include "RefCountedPtr.H"

#include "NamespaceHeader.H"

/**
 * Generic grid function library.
*/
class GridFunctionLibrary
{
public:

   static GridFunctionLibrary* getInstance();

   const RefCountedPtr<GridFunction> find( const std::string& name );

private:

   GridFunctionLibrary() : m_verbosity(0) {;}

   ~GridFunctionLibrary();

   void parseFromInput();

   void printLibrary();

   void registerGridFunction( const std::string& name,
                              RefCountedPtr<GridFunction> function );

   static GridFunctionLibrary* s_grid_function_library_instance;

   typedef std::map<std::string,RefCountedPtr<GridFunction> > Library;

   std::map<std::string,RefCountedPtr<GridFunction> > m_library;

   int m_verbosity;
};

#include "NamespaceFooter.H"

#endif
