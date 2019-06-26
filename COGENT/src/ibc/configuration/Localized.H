#ifndef _LOCALIZED_H_
#define _LOCALIZED_H_

#include "GridFunction.H"
#include "NamespaceHeader.H"

/**
 * Localized function class.
 * Assigns Gaussian functions of all space variables in the form
 *
 *   f = A \exp( \sum_{i=1}^D ((x_i-c_i)/w_i)^2 ) + floor
 *
 * where A is the peak scale, floor is the background floor, 
 * the vector c is the center of the pulse, and the vector w
 * give the width scales for each direction.
 *
 * Implemented for both annular and single null geometries.
 *
 * Required input keys and data types:
 * \b amplitude
 * real scale of the peak of the Gaussian function
 *
 * \b location
 * real coordinates of peak of Gaussian in physical coordinates (R,Z)
 *
 * \b width
 * real vector of 1/e widths of Gaussian in each direction
 *
 * \b floor
 * constant real floor
 *
 * The following represents a sample input entry for this function choice.
 *
 * \verbatim
 * type = "localized"
 * amplitude = 1.0
 * location = 0.5 1.0
 * width = 0.025 0.1
 * floor = 1.0
 *
*/
class Localized : public GridFunction
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
      Localized( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      virtual ~Localized() {;}

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      Localized( const Localized& );
      Localized& operator=( const Localized& );

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

      Real        m_amplitude;
      RealVect    m_location;
      RealVect    m_width;
      Real        m_floor;
};

#include "NamespaceFooter.H"

#endif
