#ifndef _TIME_FUNCTION_H_
#define _TIME_FUNCTION_H_

#include "REAL.H"

#include "NamespaceHeader.H"

/**
 * Generic time (0D) function interface.
*/
class TimeFunction
{
   public:

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new TimeFunction object.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
      TimeFunction( const int& verbosity );
   
      /// Destructor.
      /**
       */
      virtual ~TimeFunction() {;}

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const = 0;
   
      /// get a_value = fun(a_time)
      /**
       */
      virtual void getValue( Real&  a_value,
                       const Real&  a_time ) const = 0;
   
   protected:
   
      bool m_verbosity;
  
};

#include "NamespaceFooter.H"

#endif
