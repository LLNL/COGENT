#ifndef _KINETIC_SPECIES_FACTORY_
#define _KINETIC_SPECIES_FACTORY_

#include "ParmParse.H"
#include "RefCountedPtr.H"

#include "MultiBlockLevelGeom.H"
#include "KineticFunction.H"

#include "ConstantKineticFunction.H"
#include "LocalizedKineticFunction.H"
#include "MaxwellianKineticFunction.H"
#include "CanonicalMaxwellianKineticFunction.H"
#include "TensorProductSineKineticFunction.H"
#include "ArbitraryKineticFunction.H"

#include "NamespaceHeader.H"

/**
 * Distribution function factory class.
 * Allows calling routine to create a function of abstract class
 * KineticFunction with no knowledge of the specific type.  In this way,
 * the calling code is insulted from modifications due to the addition of
 * new functions.
*/
class KineticFunctionFactory
{
   public:

      /// Constructor.
      /**
       */
      KineticFunctionFactory() {;}

      /// Destructor.
      /**
       */
      ~KineticFunctionFactory() {;}

      /// Create new KineticFunction object
      /**
       * Returns a dynamically allocated new KineticFunction object of the
       * type indicated in the passed ParmParse database.  The calling
       * routine is responsible for de-allocating the memory.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new GridFunction object.
       * @param[in] verbosity Flag to turn on additional levels of output
       *
       * The following represents a sample input entry:
       *
       * \verbatim
       * type = "maxwellian"
       */
   RefCountedPtr<KineticFunction> create( ParmParse& a_pp,
                                          const int& a_verbosity = 0 );
};


inline
RefCountedPtr<KineticFunction> KineticFunctionFactory::create(
   ParmParse& a_pp,
   const int& a_verbosity )
{
   std::string model( "unknown" );
   a_pp.get( "type", model );

   RefCountedPtr<KineticFunction> result(NULL);
   if (model == "localized") {
      result = RefCountedPtr<KineticFunction>(
         new LocalizedKineticFunction( a_pp, a_verbosity ) );
   }
   else if (model == "maxwellian") {
      result = RefCountedPtr<KineticFunction>(
         new MaxwellianKineticFunction( a_pp, a_verbosity ) );
   }
   else if (model == "canonical_maxwellian") {
      result = RefCountedPtr<KineticFunction>(
         new CanonicalMaxwellianKineticFunction( a_pp, a_verbosity ) );
   }
   else if (model == "constant") {
      result = RefCountedPtr<KineticFunction>(
         new ConstantKineticFunction( a_pp, a_verbosity ) );
   }
   else if (model == "zero") {
      result = RefCountedPtr<KineticFunction>(
         new ConstantKineticFunction( 0.0, a_verbosity ) );
   }
   else if (model == "tensor_product_sine") {
      result = RefCountedPtr<KineticFunction>(
         new TensorProductSineKineticFunction( a_pp, a_verbosity ) );
   }
   else if (model == "arbitrary") {
      result = RefCountedPtr<KineticFunction>(
         new ArbitraryKineticFunction( a_pp, a_verbosity ) );
   }
   else {
      MayDay::Error( "KineticFunctionFactory: Unimplemented kinetic species function selected" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
