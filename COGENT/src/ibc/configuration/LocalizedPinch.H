#ifndef _LOCALIZEDPINCH_H_
#define _LOCALIZEDPINCH_H_

#include "GridFunction.H"
#include "NamespaceHeader.H"

/**
 * LocalizedPinch function class.
 *
 * Assigns specific localized radial function for magnetic field (B)
 * Corresponding to a pinch equilibrium in cylindrical coordinates:
 *
 *   B is solution of d(r^2*)B^2/2 + P))/dr = 2*r*P, 
 *   where P = P00 + P01*exp(-c*r^2) + P10*exp(-a*r^2-b*r^4) + P11*exp(-d*r^2-b*r^4) 
 *
 *   *** NOTE the solution requires all coefficients in exponetials to be nonzero ***
 *
 * The following represents a sample input entry for this function choice.
 *
 * \verbatim
 * type = "localizedPinch"
 * PreCoeffs = P00 P01 P10 P11
 * expCoeffs = a b c d
 *
*/
class LocalizedPinch : public GridFunction
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
      LocalizedPinch( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      virtual ~LocalizedPinch() {;}

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      LocalizedPinch( const LocalizedPinch& );
      LocalizedPinch& operator=( const LocalizedPinch& );

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

      Vector<Real>  m_PreCoeffs;
      Vector<Real>  m_expCoeffs;
      std::string   m_varType;

};

#include "NamespaceFooter.H"

#endif
