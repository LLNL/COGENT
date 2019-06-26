#ifndef _ARBITRARY_H_
#define _ARBITRARY_H_

#include "GridFunction.H"
#include "ParsingCore.H"

#include "NamespaceHeader.H"

/**
 * Arbitrary function class.
 * Assigns arbitrary function in configuration space:
 * 
 * coordinate_type="mapped" uses the mapped r and theta coordinates as "x" and "y"  (default)
 * coordinate_type="physical" uses the physical R and Z coordinates as "x" and "y"
 * coordinate_type="flux" uses the normalized flux as "x" and mapped theta coordinate as "y"
 * coordinate_type="outer_midplane" uses the shifter outer radial coordinate (R-Rsep) as "x" and mapped theta coordinate as "y"
 *
 * The following represents a sample input entry for this function choice.
 *
 * \verbatim
 * type = "arbitrary"
 * coordinate_type="flux"
 * function="sin(x)+sin(y)"
 *
*/
class Arbitrary : public GridFunction
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
      Arbitrary( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      virtual ~Arbitrary() {delete m_pscore;}

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      Arbitrary( const Arbitrary& );
      Arbitrary& operator=( const Arbitrary& );

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

      ParsingCore *m_pscore;
      ParsingCore *m_pscore2;

      std::string m_function;
      std::string m_function2;
      std::string m_coord_type;

};

#include "NamespaceFooter.H"

#endif
