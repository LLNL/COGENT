#ifndef _Grid_Function_FACTORY_
#define _Grid_Function_FACTORY_

#include "ParmParse.H"
#include "MultiBlockLevelGeom.H"
#include "GridFunction.H"

#include "Cosine.H"
//#include "CosTensorProd.H"
#include "Tanh.H"
#include "SingleNullDistr.H"
#include "Localized.H"
#include "LocalizedPinch.H"
#include "Constant.H"
#include "DataArray.H"
#include "RZdata.H"
#include "OneDimData.H"
#include "Arbitrary.H"

#include "NamespaceHeader.H"

/**
 * Configuration space function factory class.
 * Allows calling routine to create a function of abstract class
 * GridFunction with no knowledge of the specific type.  In this way,
 * the calling code is insulted from modifications due to the addition of
 * new functions.
*/
class GridFunctionFactory
{
   public:

      /// Constructor.
      /**
       */
      GridFunctionFactory() {;}

      /// Destructor.
      /**
       */
      ~GridFunctionFactory() {;}

      /// Create new GridFunction object
      /**
       * Returns a dynamically allocated new GridFunction object of the
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
       * type = "localized"
       */
      inline RefCountedPtr<GridFunction> create( ParmParse& pp,
                                                 const int& verbosity = 0 );
};


inline
RefCountedPtr<GridFunction> GridFunctionFactory::create(
   ParmParse& a_pp,
   const int& a_verbosity )
{
   std::string model( "unknown" );
   a_pp.get( "type", model );

   std::string name( "foo" );

   RefCountedPtr<GridFunction> result(NULL);
   if (model == "cosine") {
      result = RefCountedPtr<GridFunction>( new Cosine( a_pp, a_verbosity ) );
   }
   else if (model == "localized") {
      result = RefCountedPtr<GridFunction>( new Localized( a_pp, a_verbosity ) );
   }
   else if (model == "localizedPinch") {
      result = RefCountedPtr<GridFunction>( new LocalizedPinch( a_pp, a_verbosity ) );
   }
   else if (model == "tanh") {
      result = RefCountedPtr<GridFunction>( new Tanh( a_pp, a_verbosity ) );
   }
   else if (model == "SingleNullDistr") {
      result = RefCountedPtr<GridFunction>( new SingleNullDistr( a_pp, a_verbosity ) );
   }
   else if (model == "constant") {
      result = RefCountedPtr<GridFunction>( new Constant( a_pp, a_verbosity ) );
   }
//   else if (model == "costensorprod") {
//      result = RefCountedPtr<GridFunction>( new CosTensorProd( a_pp, a_verbosity ) );
//   }
   else if (model == "zero") {
      result = RefCountedPtr<GridFunction>( new Constant( 0.0, a_verbosity ) );
   }
   else if (model == "data_array") {
      result = RefCountedPtr<GridFunction>( new DataArray( a_verbosity ) );
   }
   else if (model == "RZdata") {
      result = RefCountedPtr<GridFunction>( new RZdata( a_pp, a_verbosity ) );
   }
   else if (model == "OneDimData") {
     result = RefCountedPtr<GridFunction>( new OneDimData( a_pp, a_verbosity ) );
   }
   else if (model == "arbitrary") {
      result = RefCountedPtr<GridFunction>( new Arbitrary( a_pp, a_verbosity ) );
   }
   else {
      MayDay::Error( "GridFunctionFactory: Unimplemented grid function type selected" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
