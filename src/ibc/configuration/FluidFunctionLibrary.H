#ifndef _FLUID_FUNCTION_LIBRARY_H_
#define _FLUID_FUNCTION_LIBRARY_H_

#include <string>
#include <map>

#include "FluidFunction.H"
#include "RefCountedPtr.H"

#include "NamespaceHeader.H"

/**
 * Generic fluid function library.
*/
class FluidFunctionLibrary
{
public:

   static FluidFunctionLibrary* getInstance();

   const RefCountedPtr<FluidFunction> find( const std::string& name );

private:

   FluidFunctionLibrary() : m_verbosity(0) {;}

   ~FluidFunctionLibrary();

   void parseFromInput();

   void printLibrary();

   void registerFluidFunction( const std::string& name,
                              RefCountedPtr<FluidFunction> function );

   static FluidFunctionLibrary* s_fluid_function_library_instance;

   typedef std::map<std::string,RefCountedPtr<FluidFunction> > Library;

   std::map<std::string,RefCountedPtr<FluidFunction> > m_library;

   int m_verbosity;
};

#include "NamespaceFooter.H"

#endif
