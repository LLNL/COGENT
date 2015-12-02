#ifndef _SINGLENULLDISRT_H_
#define _SINGLENULLDISRT_H_

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
 * SingleNullDistr function class.
 * Assigns flux-function-dependent distribution f(psi) for a single-null geometry 
 *
 * if subtype = "Lorentz"
 *
 *   f(psi) = core_value/((Z(psi)-Z_core)^2/width^2 + 1.0) for core and sol regions
 *
 *   f(psi) = sep_value + (sep_value - pf_value)  
 *                       * (Z(psi) - Z_xpt)/(Z_xpt - Z_pf) for a private flux region
 *       
 *
 * where core_value, sep_value, and pf_value (input paramteres) are values of f at the
 * core boundary, separatrix and outer PF flux surfaces. For the core and sol blocks, 
 * Z (psi) corresponds to the physical vertical coordinate at the top of a given
 * flux surface, psi. For the pf block case, Z(psi) is the vertical coordinate
 * at the mutual boundary of LPF and RPF blocks (i.e., the poloidal median of the 
 * private-flux block). Z_core corresponds to Z(psi) for the inner core boundary, 
 * Z_xpt is the vertical physical coordinate of the x-point, Z_pf corresponds to 
 * Z(psi) for the outer PF boundary.The "width" parameter is determined self-consistently
 * to satisfy the condition that f at the separatrix (given by the core/sol specification) 
 * is consistent with sep_value. 
 *
 * Implemented only for a single null geometry.
 *
 * Required input keys and data types:
 * \b core_value
 * value of f at the (inner) core boundary
 *
 * \b sep_value
 *  value of f at the separatrix
 *  
 * \b pf_vale
 * value of f at the flux-surface corresponding to the bottom of the PF region
 *
 * The following represents a sample input entry for this function choice.
 *
 * \verbatim
 * type = "SingleNullDistr"
 * subtype = "Lorentz"
 * core_value = 1.0
 * sep_value = 0.5
 * pf_value = 0.0
 *
 * if subtype = "Tanh"
 * See the definitions and verbatim in Tanh.H, 
 * NB: midpoint_fraction in this Tanh is related to radial_midpoint in Tanh.H as
 * radial_midpoint = x_inner + (x_outer-x_inner) * midpoint_fraction.
 * If boltzmann_equilibirum is set to true (default = false), then assigns distribution 
 * corresponding to the boltzmann_equilibirum for a potential for a density given by "Tanh".
 * 
 * if subtype = "Localized"
 * See the definitions and verbatim in Localized.H, 
 * NB: midpoint_fraction in this Localized is related to the location variable in Localized.H as
 * location = x_inner + (x_outer-x_inner) * midpoint_fraction.
 *
*/
class SingleNullDistr : public GridFunction
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
      SingleNullDistr( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      virtual ~SingleNullDistr() {;}

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
      virtual void assign( LevelData<FArrayBox>& data,
                           const MultiBlockLevelGeom& geometry,
                           const BoundaryBoxLayout& bdry_layout,
                           const Real& time ) const;

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      SingleNullDistr( const SingleNullDistr& );
      SingleNullDistr& operator=( const SingleNullDistr& );

      /// Set cell-average values.
      /**
       */
      inline void setCellAverages( FArrayBox&              data,
                                   const MultiBlockCoordSys& coords,
                                   const MagBlockCoordSys& block_coord_sys,
                                   const int               block_number ) const;

      /// Set pointwise values at cell centers.
      /**
       */
      inline void setPointwise( FArrayBox&              data,
                                const MultiBlockCoordSys& coords,
                                const MagBlockCoordSys& block_coord_sys,
                                const int               block_number ) const;

      void getGeomParameters(const MultiBlockCoordSys& a_coord_sys) const;

      /// Set radial distribution for the CORE and SOL regions.
      /**
       */
      inline
      Real radialCoreSol( const Real                a_Z,
                          const MultiBlockCoordSys& coords)  const;

      /// Set Radial distribution for the PF region.
      /**
       */
      inline
      Real radialPF( const Real                a_Z,
                     const MultiBlockCoordSys& coords)  const;


      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& pp );

      /// Check that the geometry is annular or single null.
      /**
       */
      inline void checkGeometryValidity( const MultiBlockLevelGeom& geometry ) const;


      /// Computes and returns tanh function (as in TanhF.ChF)
      /**
      * 
      */
      double formTanh(double x, double x1, double x2) const;

      int m_verbosity;

      Real        m_core_value;
      Real        m_sep_value;
      Real        m_pf_value;
      
      std::string m_subtype;

      Real m_inner_radial_value;
      Real m_outer_radial_value;
      Real m_midpoint_fraction;
      Real m_radial_width;
      Real m_amplitude;
      Real m_floor; 
      bool m_boltzmann_equilibrium;


};

#include "NamespaceFooter.H"

#endif
