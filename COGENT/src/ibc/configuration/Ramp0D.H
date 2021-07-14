#ifndef _RAMP0D_H_
#define _RAMP0D_H_

#include <cmath>

#include "TimeFunction.H"

#include "NamespaceHeader.H"

/**
 * Generic 0D constant function class.
 * Assigns Real value to a constant.
 *
 * Required input keys and data types:
 * \b value
 * real constant value to set
 *
 * The following represents a sample input entry for this function choice:
 *
 * \verbatim
 * type = "constant"
 * value = 1.0
*/
class Ramp0D : public TimeFunction
{
   public:

      /// Constructor using specified values.
      /**
       * Constructs the object using a passed-in parameter.
       *
       */
      Ramp0D( const Real&  a_rise_time,
              const Real&  a_amplitude, 
              const int&   verbosity )
         : TimeFunction(verbosity), 
           m_rise_time(a_rise_time),
           m_amplitude(a_amplitude) {;}

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       */
      Ramp0D( ParmParse& pp, const int& verbosity )
         : TimeFunction(verbosity), 
           m_rise_time(0.0),
           m_amplitude(0.0)
      {
         parseParameters( pp );
      }

      /// Destructor.
      /**
       */
      virtual ~Ramp0D() {;}

      /// Print object parameters.
      /**
       */
      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Ramp0D function parameters:" << std::endl;
            std::cout << " rise time = " << m_rise_time  << std::endl;
            std::cout << " amplitude = " << m_amplitude  << std::endl << std::endl;
         }
      }
      
      /// get a_value = fun(a_time)
      /**
       */
      virtual void getValue( Real&  a_value,
                       const Real&  a_time ) const
      {
         a_value = a_time/m_rise_time*m_amplitude;
         a_value = std::min(a_value,m_amplitude);
      }

   private:

      // prohibit copying
      Ramp0D( const Ramp0D& );
      Ramp0D& operator=( const Ramp0D& );

      /// Parse the input database for parameters.
      /**
       */
      void parseParameters( ParmParse& pp )
      {
         pp.get( "rise_time", m_rise_time );
         pp.get( "amplitude", m_amplitude );
         if (m_verbosity) {
            printParameters();
         }
      }

      Real m_rise_time;
      Real m_amplitude;
};

#include "NamespaceFooter.H"

#endif

