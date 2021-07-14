#ifndef _SINSQ0D_H_
#define _SINSQ0D_H_

#include <cmath>

#include "TimeFunction.H"

#include "NamespaceHeader.H"

/**
 * Sin^2 0D constant function class.
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
class SinSq0D : public TimeFunction
{
   public:

      /// Constructor using specified values.
      /**
       * Constructs the object using a passed-in parameter.
       *
       */
      SinSq0D( const Real&  a_rise_time,
               const Real&  a_amplitude, 
               const int&   verbosity )
         : TimeFunction(verbosity), 
           m_rise_time(a_rise_time),
           m_exponent(2.0),
           m_amplitude(a_amplitude) {;}

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       */
      SinSq0D( ParmParse& pp, const int& verbosity )
         : TimeFunction(verbosity), 
           m_rise_time(0.0),
           m_exponent(2.0),
           m_amplitude(0.0)
      {
         parseParameters( pp );
      }

      /// Destructor.
      /**
       */
      virtual ~SinSq0D() {;}

      /// Print object parameters.
      /**
       */
      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "Sinsq0D function parameters:" << std::endl;
            std::cout << " rise time = " << m_rise_time  << std::endl;
            std::cout << " exponent  = " << m_exponent  << std::endl << std::endl;
            std::cout << " amplitude = " << m_amplitude  << std::endl << std::endl;
         }
      }
      
      /// get a_value = fun(a_time)
      /**
       */
      virtual void getValue( Real&  a_value,
                       const Real&  a_time ) const
      {
         a_value = sin(3.14159/2.0*a_time/m_rise_time);
         //a_value = a_value*a_value*m_amplitude;
         a_value = pow(a_value,m_exponent)*m_amplitude;
      }

   private:

      // prohibit copying
      SinSq0D( const SinSq0D& );
      SinSq0D& operator=( const SinSq0D& );

      /// Parse the input database for parameters.
      /**
       */
      void parseParameters( ParmParse& pp )
      {
         pp.get( "rise_time", m_rise_time );
         m_exponent = 2.0;
         pp.query( "exponent", m_exponent );
         pp.get( "amplitude", m_amplitude );
         if (m_verbosity) {
            printParameters();
         }
      }

      Real m_rise_time;
      Real m_exponent;
      Real m_amplitude;
};

#include "NamespaceFooter.H"

#endif

