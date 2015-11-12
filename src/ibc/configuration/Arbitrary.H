#ifndef _ARBITRARY_H_
#define _ARBITRARY_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "GridFunction.H"
#include "LevelData.H"
#include "ParmParse.H"
#include "MagBlockCoordSys.H"
#include "MultiBlockLevelGeom.H"
#include "RealVect.H"

#include "ParsingCore.H"

#include "NamespaceHeader.H"

/**
 * Arbitrary function class.
 * Assigns arbitrary function in configuration space:
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

      ParsingCore *m_pscore;

   private:

      // prohibit copying
      Arbitrary( const Arbitrary& );
      Arbitrary& operator=( const Arbitrary& );

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
      std::string m_function;

};

#include "NamespaceFooter.H"

#endif
