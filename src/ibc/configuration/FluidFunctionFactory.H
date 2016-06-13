#ifndef _Fluid_Function_FACTORY_
#define _Fluid_Function_FACTORY_

#include "ParmParse.H"
#include "MultiBlockLevelGeom.H"
#include "FluidFunction.H"

#include "ConstantFluidFunction.H"

#include "NamespaceHeader.H"

/**
 * Configuration space function factory class.
 * Allows calling routine to create a function of abstract class
 * FluidFunction with no knowledge of the specific type.  In this way,
 * the calling code is insulted from modifications due to the addition of
 * new functions.
*/
class FluidFunctionFactory
{
   public:

      /// Constructor.
      /**
       */
      FluidFunctionFactory() {;}

      /// Destructor.
      /**
       */
      ~FluidFunctionFactory() {;}

      /// Create new FluidFunction object
      /**
       * Returns a dynamically allocated new FluidFunction object of the
       * type indicated in the passed ParmParse database.  The calling
       * routine is responsible for de-allocating the memory.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new FluidFunction object.
       * @param[in] verbosity Flag to turn on additional levels of output
       *
       * The following represents a sample input entry:
       *
       * \verbatim
       * type = "localized"
       */
      inline RefCountedPtr<FluidFunction> create( ParmParse& pp,
                                                 const int& verbosity = 0 );
};


inline
RefCountedPtr<FluidFunction> FluidFunctionFactory::create(
   ParmParse& a_pp,
   const int& a_verbosity )
{
   std::string model( "unknown" );
   a_pp.get( "type", model );

   std::string name( "foo" );

   RefCountedPtr<FluidFunction> result(NULL);
   if (model == "constant") {
      result = RefCountedPtr<FluidFunction>( new ConstantFluidFunction( a_pp, a_verbosity ) );
   }
   else if (model == "zero") {
      result = RefCountedPtr<FluidFunction>( new ConstantFluidFunction( 0.0, a_verbosity ) );
   }
   else {
      MayDay::Error( "FluidFunctionFactory: Unimplemented fluid function type selected" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
