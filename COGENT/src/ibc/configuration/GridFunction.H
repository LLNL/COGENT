#ifndef _GRID_FUNCTION_H_
#define _GRID_FUNCTION_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "MultiBlockLevelGeom.H"
#include "MagBlockCoordSys.H"
#include "BoundaryBoxLayout.H.multidim"

#include "NamespaceHeader.H"

/**
 * Generic grid function interface.
*/
class GridFunction
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
      GridFunction( const int& verbosity );
   
      /// Destructor.
      /**
       */
      virtual ~GridFunction() {;}

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
      void assign(LevelData<FArrayBox>& data,
                  const MultiBlockLevelGeom& geometry,
                  const Real& time,
                  const bool& as_cell_averages = true ) const;

      void assign(LevelData<FluxBox>& data,
                  const MultiBlockLevelGeom& geometry,
                  const Real& time,
                  const bool& as_face_averages = true ) const;

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
      void assign(FArrayBox& dfn,
                  const MultiBlockLevelGeom& geometry,
                  const FArrayBox& real_coords,
                  const FArrayBox& normalized_flux,
                  const int  block_number,
                  const Real& time,
                  const bool& cell_averages = true ) const;

      void assign(FluxBox& dfn,
                  const MultiBlockLevelGeom& geometry,
                  const FluxBox& real_coords,
                  const FluxBox& normalized_flux,
                  const int block_number,
                  const Real& time,
                  const bool& face_averages = true ) const;

      /// Set the function on the boundary boxes.
      /**
       * Sets the cell-averaged function value for the grid function
       * at the given time on the boxes specified in the BoundaryBoxLayout.
       *
       * NB: The cell-averaged values include the J and BStarParallel factors,
       * i.e., the internal LevelData contains the average of fJB^*_{\parallel}
       *
       * @param[inout] data  Grid function to set.
       * @param[in] bdry_layout BoundaryBoxLayout defining boxes to fill.
       * @param[in] time Time value at which to evaluate the function.
       */
      void assign(LevelData<FArrayBox>& data,
                  const MultiBlockLevelGeom& geometry,
                  const BoundaryBoxLayout& bdry_layout,
                  const Real& time ) const;

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const = 0;
   
   protected:
   
      /// Check that the geometry is consistent with a particular grid function realization.
      /**
       */
      virtual void checkGeometryValidity( const MultiBlockLevelGeom& geometry ) const {;}
   
      /// Set pointwise values at cell centers.
      /**
       */
      virtual void setPointwise(FArrayBox&                  dfn,
                                const MultiBlockLevelGeom&  geometry,
                                const FArrayBox&            real_coords,
                                const FArrayBox&            normalized_flux,
                                const int                   block_number ) const = 0;
      /// Set cell-average values.
      /**
       */
      void setCellAverages(FArrayBox&                  data,
                           const MultiBlockLevelGeom&  geometry,
                           const FArrayBox&            real_coords,
                           const FArrayBox&            normalized_flux,
                           const int                   block_number) const;
      /// Set face-average values.
      /**
       */
      void setFaceAverages(FluxBox&                    data,
                           const MultiBlockLevelGeom&  geometry,
                           const FluxBox&              real_coords,
                           const FluxBox&              normalized_flux,
                           const int                   block_number) const;
   

      /// Returns block coordinate system
      const MagBlockCoordSys& getCoordSys(const MultiBlockLevelGeom& geometry,
                                          const Box&                 box ) const;


      /// Returns block coordinate system
      const MagBlockCoordSys& getCoordSys(const MultiBlockLevelGeom& geometry,
                                          const int                  block_number ) const;
   
      bool m_verbosity;
  
};

#include "NamespaceFooter.H"

#endif
