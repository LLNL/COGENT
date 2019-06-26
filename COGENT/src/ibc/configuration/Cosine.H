#ifndef _COSINE_H_
#define _COSINE_H_

#include "GridFunction.H"
#include "NamespaceHeader.H"

/**
 * Cosine function class.
 * Assigns cosine tensor product function in configuration space:
 *
 *   f = C + A \cdot \prod_{i=1}^D \cos(2\pi k_i x_i / L_i + phi_i)
 *
 * where C is a constant, A is the amplitude, the components k_i of the mode
 * vector are the wavenumbers in the corresponding directions, and the
 * components \phi_i of the phase vector are the phase shifts in the
 * corresponding directions.
 *
 * This function has only been implemented for annular geometries.
 *
 * Required input keys and data types:
 * \b constant
 * real constant offset C
 *
 * \b amplitude
 * real amplitude A
 *
 * \b mode
 * real vector of wavenumbers k
 *
 * \b phase
 * real vector fo phase shifts \phi
 *
 * The following represents a sample input entry for this function choice.
 *
 * \verbatim
 * type = "cosine"
 * constant = 1.0
 * amplitude = 1.0
 * mode = 1.0 1.0
 * phase = 0.0 0.0
 *
*/
class Cosine : public GridFunction
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
      Cosine( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      virtual ~Cosine() {;}

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      Cosine( const Cosine& );
      Cosine& operator=( const Cosine& );

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

      Real     m_constant;
      Real     m_amplitude;
      RealVect m_mode;
      RealVect m_phase;
};

#include "NamespaceFooter.H"

#endif
