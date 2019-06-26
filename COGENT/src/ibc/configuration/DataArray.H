#ifndef _DATAARRAY_H_
#define _DATAARRAY_H_

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
      DataArray(const bool& verbosity) : GridFunction(verbosity) {};
     
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

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      DataArray( const DataArray& );
      DataArray& operator=( const DataArray& );

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
      virtual void checkGeometryValidity( const MultiBlockLevelGeom& geometry ) const {;}
   
      FArrayBox m_data;
      bool m_cell_averages;
};

#include "NamespaceFooter.H"

#endif

