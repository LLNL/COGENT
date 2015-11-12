#ifndef _COSINE_H_
#define _COSINE_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "GridFunction.H"
#include "LevelData.H"
#include "ParmParse.H"
#include "MagBlockCoordSys.H"
#include "MultiBlockLevelGeom.H"
#include "RealVect.H"

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

      /// Set the grid function.
      /**
       * Sets the grid function on the LevelData for the given geometry
       * at the given time.  It is assumed that derived classes will check
       * that the specification is suitable for the given geometry.  Cell
       * average values are set unless the optional parameter is false, in
       * which case, pointwise values are set.
       *
       * @param[out] data LevelData to set.
       * @param[in] geometry Corresponding level geometry.
       * @param[in] time Time value at which to evaluate the function.
       * @param[in] as_cell_name Boolean switch for cell average values or
       *                         point values.
       */
      virtual void assign( LevelData<FArrayBox>& data,
                           const MultiBlockLevelGeom& geometry,
                           const Real& time,
                           const bool& cell_averages = true ) const;

      /// Set the grid function.
      /**
       * Sets the grid function on the box for the given geometry
       * at the given time.  It is assumed that derived classes will check
       * that the specification is suitable for the given geometry.  Cell
       * average values are set unless the optional parameter is false, in
       * which case, pointwise values are set.
       *
       * @param[out] data LevelData to set.
       * @param[in] geometry Corresponding level geometry.
       * @param[in] box Box from which to get coordinate system.
       * @param[in] time Time value at which to evaluate the function.
       * @param[in] as_cell_name Boolean switch for cell average values or
       *                         point values.
       */
      virtual void assign( FArrayBox& data,
                           const MultiBlockLevelGeom& geometry,
                           const Box& box,
                           const Real& time,
                           const bool& cell_averages = true ) const;

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      Cosine( const Cosine& );
      Cosine& operator=( const Cosine& );

      /// Set cell-average values.
      /**
       */
      inline void setCellAverages( FArrayBox&              a_dfn,
                                   const MagBlockCoordSys& a_coord_sys ) const;

      /// Set pointwise values at cell centers.
      /**
       */
      inline void setPointwise( FArrayBox&              a_dfn,
                                const MagBlockCoordSys& a_coord_sys ) const;

      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& a_pp );

      /// Check that the geometry is annular.
      /**
       */
      inline void checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const;

      int m_verbosity;

      Real     m_constant;
      Real     m_amplitude;
      RealVect m_mode;
      RealVect m_phase;
};

#include "NamespaceFooter.H"

#endif
