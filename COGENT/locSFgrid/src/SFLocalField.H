#ifndef _SFLOCALFIELD_H_
#define _SFLOCALFIELD_H_

#include "FArrayBox.H"
#include "RealVect.H"
#include <string>

#include "NamespaceHeader.H"

class SFLocalField
{
  public:

   SFLocalField(const double& Rmaj      , 
                const double& zmax      , 
                const double& zm0       ,
                const string& block_name);

    ~SFLocalField();

   void getPhysicalCoordinates(const FArrayBox& xi                  ,
                               FArrayBox      & physical_coordinates) const;

   void getFieldUnitVector(const FArrayBox& physical_coordinates,
                           FArrayBox&       field_unit_vector) const;

   void writePhysicalCoordinates(const FArrayBox& physical_coordinates) const;

   void writeVectors(const FArrayBox& physical_coordinates,
                     const FArrayBox& vector_data) const;

   RealVect getMappedCoordinate(const RealVect& X ) const;


  private:

   void setIdealSFCoefs();

   int m_block_number;

   double m_zmax  ;

   double m_zm0   ;

   double m_Rmaj  ;

   double m_phimax;

   double m_beta  ;

   double third   ;
  
  // rotate to vertical orientation
  RealVect rotateToVert(const RealVect& X ) const;  

  // rotate from vertial orientation
  RealVect rotateFromVert(const RealVect& X) const; 

  double newtonSolveCubic(const double a, 
                          const double b, 
                          const double c, 
                          const double d, 
                          double guess  ) const;

   /* the coefficients here for the snowflake expansion */
   double m_l1;
   double m_l2;
   double m_q1;
   double m_q2;
   double m_q3;
   double m_c1;
   double m_c2;
   double m_c3;
   double m_c4;

  /* coefficients for rotations for the 12 blocks */
  double m_costh [12];
  double m_sinth [12];

};

#include "NamespaceFooter.H"

#endif
