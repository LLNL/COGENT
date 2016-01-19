#ifndef _POTENTIALBCFACTORY_H_
#define _POTENTIALBCFACTORY_H_


#include "ParmParse.H"
#include "PotentialBC.H"
#include "AnnulusPotentialBC.H"
#include "SlabPotentialBC.H"
#include "RectangularTorusPotentialBC.H"
#include "SingleNullPotentialBC.H"
#include "SNCorePotentialBC.H"

#include <string>

#include "NamespaceHeader.H"

class PotentialBCFactory
{
   public:

      /// Constructor.
      /**
       */
      PotentialBCFactory() {;}

      /// Destructor.
      /**
       */
      ~PotentialBCFactory() {;}

      /// Create new PotentialBC object
      /**
       * Returns a dynamically allocated new PotentialBC object of the
       * type indicated in the passed ParmParse database.  The calling
       * routine is responsible for de-allocating the memory.
       *
       * @param[in] name String identifier of this PotentialBC.
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new PotentialBC object.
       * @param[in] verbosity Flag to turn on additional levels of output
       *
       * The following represents a sample input entry:
       *
       * \verbatim
       * type = "localized"
       */
      inline PotentialBC* create( const std::string& name,
                                  ParmParse& pp,
                                  const std::string& coord_sys_type,
                                  const int& verbosity = 0 );
};

inline
PotentialBC* PotentialBCFactory::create(
   const std::string& a_name,
   ParmParse& a_pp,
   const std::string& a_coord_sys_type,
   const int& a_verbosity )
{
   PotentialBC* result(NULL);
   if (a_coord_sys_type == "Annulus" ) {
      result = new AnnulusPotentialBC( a_name, a_pp, a_verbosity );
   }
   else if (a_coord_sys_type == "Slab" ) {
      result = new SlabPotentialBC( a_name, a_pp, a_verbosity );
   }
   else if (a_coord_sys_type == "RectangularTorus" ) {
      result = new RectangularTorusPotentialBC( a_name, a_pp, a_verbosity );
   }
   else if (a_coord_sys_type == "SingleNull" ) {
      result = new SingleNullPotentialBC( a_name, a_pp, a_verbosity );
   }
   else if (a_coord_sys_type == "SNCore" ) {
      result = new SNCorePotentialBC( a_name, a_pp, a_verbosity );
   }
   else {
      MayDay::Error( "PotentialBCFactory: Unimplemented boundary condition for given coordinate system type" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
