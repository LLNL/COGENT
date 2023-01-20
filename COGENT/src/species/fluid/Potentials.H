#ifndef _POTENTIALS_H_
#define _POTENTIALS_H_

#include "CFGVars.H"
#include "NamespaceHeader.H"


/// This class contains electrostatic and magnetic vector potentials
/**
 * Potentials interface
 */
class Potentials
   : public CFGVars
{
   public:

      /// Constructor.
      /**
       */
      Potentials(const string&       pp_prefix,
                 const std::string&  name,
                 const MagGeom&      geometry,
                 const IntVect&      ghost_vect,
                 const bool          evolve_A_parallel);

      /// Destructor.
      /**
       */
      ~Potentials() {;}
   
      virtual
      RefCountedPtr<CFGVars> clone( const IntVect& ghost_vect = IntVect::Zero,
                                    const bool copy_data = true ) const;

      bool evolvingApar() const {return m_evolve_A_parallel;}
   
      bool m_evolve_A_parallel;


};

#include "NamespaceFooter.H"

#endif
