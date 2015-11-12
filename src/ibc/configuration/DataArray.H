#ifndef _DATAARRAY_H_
#define _DATAARRAY_H_

#include "MultiBlockLevelGeom.H"

#include "GridFunction.H"

#include "NamespaceHeader.H"

/**
 * Class implementing a time-dependent function on a physical boundary
 * section defined by an FArrayBox of data
 * Simply fills an FArrayBox with previously stored data
 *
 * \verbatim
*/
class DataArray : public GridFunction
{
   public:

      /// Default constructor.
      /**
       * Default constructor
       *
       * @param[in] verbosity Flag to turn on additional levels of output
       */
      DataArray(const bool& verbosity) : m_verbosity(verbosity) {};
     
      /// Destructor.
      /**
       */
      virtual ~DataArray() {;}

      /// Sets the object data.
      /**
       * Sets the object data using a passed-in FArrayBox.  Any previously
       * set data is destroyed.
       *
       * @param[in] data           FArrayBox containing data from which to set
       * @param[in] cell_averages  Whether data is cell-averaged (assumed cell-centerd if false)
       */
      void setData( const FArrayBox& data,
                    const bool&      cell_averages = true );

      void setData( const DataArray& data );

      /// Returns the box over which the object data is defined
      /**
       * Returns the box over which the object data is defined
       *
       * @param[out] box           Box of the FArrayBox data
       */
      Box getDefinedBox() const {return m_data.box();}

      /// Set the grid function.
      /**
       * Sets the grid function on the LevelData for the given geometry
       * at the given time.  It is assumed that derived classes will check
       * that the specification is suitable for the given geometry.  Cell
       * average values are set unless the optional parameter is false, in
       * which case, pointwise values are set.
       *
       * NB: No check is made to ensure that sufficient data is available
       * to fill the target.
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
       * NB: No check is made to ensure that sufficient data is available
       * to fill the target.
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
      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "DataArray parameters:" << std::endl;
         }
      }

   private:

      // prohibit copying
      DataArray( const DataArray& );
      DataArray& operator=( const DataArray& );

      FArrayBox m_data;
      bool m_cell_averages;
      int m_verbosity;
};

#include "NamespaceFooter.H"

#endif

