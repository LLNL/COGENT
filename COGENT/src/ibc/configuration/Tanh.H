#ifndef _TANH_H_
#define _TANH_H_

#include "GridFunction.H"

#include "NamespaceHeader.H"

/**
 * Tanh function class.
 * Assigns hyperbolic tangent in flux coordnates across an annular
 * geometry:
 *
 *   f = a \tanh( (\psi-\psi_m) / w ) + b
 *
 * where
 *
 *   a = (f_i-f_o) / (\tanh(\phi_i-\phi_m)-\tanh(\phi_o-\phi_m)),
 *   b = f_o - a \tanh(\phi_o-\phi_m),
 *
 * f_i and f_o are the values at the inner (\psi_i) and outer (\psi_o)
 * flux surfaces, and \psi_m is the flux coordinate at which the inflection
 * occurs.
 *
 * This function has only been implemented for annular geometries.
 *
 * Required input keys and data types:
 * \b inner_radial_value
 * real value at inner flux surface
 *
 * \b outer_radial_value
 * real value at outer flux surface
 *
 * \b radial_midpoint
 * real value, which specify the inflection flux surface by 
 * \psi_m = \psi_i + radial_midpoit*(\psi_o -\psi_i)
 *
 * \b radial_width
 * real width scale in normalized flux coordinates
 *
 * The following represents a sample input entry for this function choice.
 *
 * \verbatim
 * type = "tanh"
 * inner_radial_value = 0.0
 * outer_radial_value = 1.0
 * radial_midpoint = 0.5
 * radial_width = 0.25
 *
*/
class Tanh : public GridFunction
{
   public:

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new GridFunction object.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
      Tanh( ParmParse& a_pp, const int& a_verbosity );

       /// Destructor.
      /**
       */
      virtual ~Tanh() {;}

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      Tanh( const Tanh& );
      Tanh& operator=( const Tanh& );

      /// Set pointwise values at cell centers.
      /**
       */
      virtual void setPointwise(FArrayBox&                  dfn,
                                const MultiBlockLevelGeom&  geometry,
                                const FArrayBox&            real_coords,
                                const FArrayBox&            normalized_flux,
                                const int                   block_number ) const;
   
      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& pp );
   
      /// Check geometry validity.
      /**
       */
      virtual void checkGeometryValidity( const MultiBlockLevelGeom& geometry ) const;

      Real m_inner_radial_value;
      Real m_outer_radial_value;
      Real m_radial_midpoint;
      Real m_radial_width;
};

#include "NamespaceFooter.H"

#endif
