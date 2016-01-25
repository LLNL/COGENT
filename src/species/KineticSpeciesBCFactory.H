#ifndef _KINETIC_SPECIES_BC_FACTORY_H_
#define _KINETIC_SPECIES_BC_FACTORY_H_

#define CH_SPACEDIM PDIM

#include "ParmParse.H"
#include "KineticSpecies.H"
#include "KineticSpeciesBC.H"
#include "AnnulusPhaseBC.H"
#include "SlabPhaseBC.H"
#include "SingleNullPhaseBC.H"
#include "SNCorePhaseBC.H"

#include <string>

#include "NamespaceHeader.H"

class KineticSpeciesBCFactory
{
   public:

      /// Constructor.
      /**
       */
      KineticSpeciesBCFactory() {;}

      /// Destructor.
      /**
       */
      ~KineticSpeciesBCFactory() {;}

      /// Create new InitialCondition object
      /**
       * Returns a dynamically allocated new InitialCondition object of the
       * type indicated in the passed ParmParse database.  The calling
       * routine is responsible for de-allocating the memory.
       *
       * @param[in] name String identifier of this initial condition.
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new InitialCondition object.
       * @param[in] verbosity Flag to turn on additional levels of output
       *
       * The following represents a sample input entry:
       *
       * \verbatim
       * type = "localized"
       */
      inline KineticSpeciesBC* create( const std::string& name,
                                       ParmParse& pp,
                                       const std::string& coord_sys_type,
                                       const int& verbosity = 0 );
};

inline
KineticSpeciesBC* KineticSpeciesBCFactory::create(
   const std::string& a_name,
   ParmParse& a_pp,
   const std::string& a_coord_sys_type,
   const int& a_verbosity )
{
   KineticSpeciesBC* result(NULL);
   if (a_coord_sys_type == "Annulus" ) {
      result = new AnnulusPhaseBC( a_name, a_pp, a_verbosity );
   }
   else if (a_coord_sys_type == "SingleNull" ) {
      result = new SingleNullPhaseBC( a_name, a_pp, a_verbosity );
   }
   else if (a_coord_sys_type == "SNCore" ) {
      result = new SNCorePhaseBC( a_name, a_pp, a_verbosity );
   }
   else if (a_coord_sys_type == "Slab" ) {
      result = new SlabPhaseBC( a_name, a_pp, a_verbosity );
   }
   else {
      MayDay::Error( "KineticSpeciesBCFactory: Unimplemented boundary condition for given coordinate system type" );
   }
   return result;
}

#include "NamespaceFooter.H"

#endif
