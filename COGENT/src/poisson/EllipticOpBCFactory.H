#ifndef _ELLIPTICOPBCFACTORY_H_
#define _ELLIPTICOPBCFACTORY_H_


#include "ParmParse.H"
#include "EllipticOpBC.H"
#include "LogRectEllipticOpBC.H"
#include "SingleNullEllipticOpBC.H"
#include "SNCoreEllipticOpBC.H"

#include <string>

#include "NamespaceHeader.H"

class EllipticOpBCFactory
{
   public:

      /// Constructor.
      /**
       */
      EllipticOpBCFactory() {;}

      /// Destructor.
      /**
       */
      ~EllipticOpBCFactory() {;}

      /// Create new EllipticOpBC object
      /**
       * Returns a dynamically allocated new EllipticOpBC object of the
       * type indicated in the passed ParmParse database.  The calling
       * routine is responsible for de-allocating the memory.
       *
       * @param[in] name String identifier of this EllipticOpBC.
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new EllipticOpBC object.
       * @param[in] verbosity Flag to turn on additional levels of output
       *
       * The following represents a sample input entry:
       *
       * \verbatim
       * type = "localized"
       */
      inline EllipticOpBC* create( const std::string& name,
                                   ParmParse& pp,
                                   const std::string& coord_sys_type,
                                   const int& verbosity = 0 );
};

inline
EllipticOpBC* EllipticOpBCFactory::create(
   const std::string& a_name,
   ParmParse& a_pp,
   const std::string& a_coord_sys_type,
   const int& a_verbosity )
{
   EllipticOpBC* result(NULL);
   if (a_coord_sys_type == "LogicallyRectangular" ) {
      result = new LogRectEllipticOpBC( a_name, a_pp, a_verbosity );
   }
   else if (a_coord_sys_type == "SingleNull" ) {
      result = new SingleNullEllipticOpBC( a_name, a_pp, a_verbosity );
   }
   else if (a_coord_sys_type == "SNCore" ) {
      result = new SNCoreEllipticOpBC( a_name, a_pp, a_verbosity );
   }
   else {
      MayDay::Error( "EllipticOpBCFactory: Unimplemented boundary condition for given coordinate system type" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
