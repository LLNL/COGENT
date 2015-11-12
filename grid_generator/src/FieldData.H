#ifndef _FIELDDATA_H_
#define _FIELDDATA_H_

//#define USE_DCT_FIELD

#include "BSplineInterp.H"

#include "NamespaceHeader.H"

class FieldData
{
  public:

   FieldData(const string& geometry_file_name,
             const string& field_coefficients_file_name,
             const string& block_name);

    ~FieldData();

   void getPhysicalCoordinates(const FArrayBox& xi,
                               FArrayBox&       physical_coordinates) const;

   void getFieldUnitVectorFromMappingFile(const FArrayBox& physical_coordinates,
                                          FArrayBox&       field_unit_vector) const;

   void getFieldUnitVectorFromDCT(const FArrayBox& physical_coordinates,
                                  FArrayBox&       field_unit_vector) const;

   void convertToMappedDerivatives( const FArrayBox& a_mapped_coordinates,
                                    FArrayBox&       a_field_unit_vector ) const;

   void writePhysicalCoordinates(const FArrayBox& physical_coordinates) const;

   void writeVectors(const FArrayBox& physical_coordinates,
                     const FArrayBox& vector_data) const;

   RealVect getMappedCoordinate( const RealVect& X ) const;

  private:

   void readRZMapping(const string& file_name,
                      const int     block_number);

   void readFieldFromMappingFile(const string& file_name,
                                 const int     block_number);

   void readFieldFromDCTFile(const string& file_name);

   void precomputeCoarseMappingData(const int n);

   RealVect mappedCoordNewton( const RealVect& X,
                               const RealVect& xi_initial ) const;

   int n_radial_points, n_poloidal_points;
   int n_radial_extend, n_poloidal_extend;

   BSplineInterp* m_RZ_interpolator;
   BSplineInterp* m_field_RZ_interpolator;
   BSplineInterp* m_field_interpolator;

   FArrayBox m_physical_coords_coarse;
   FArrayBox m_mapped_coords_coarse;

   FArrayBox m_psi_coefs;
   double m_Rmin, m_Rmax, m_Zmin, m_Zmax;

   int m_block_number;
};

#include "NamespaceFooter.H"

#endif
