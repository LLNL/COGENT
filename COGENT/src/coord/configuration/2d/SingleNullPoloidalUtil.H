#ifndef _SINGLENULL_POLOIDAL_UTIL_
#define _SINGLENULL_POLOIDAL_UTIL_

#include "ParmParse.H"
#include "FArrayBox.H"
#include "RealVect.H"
#include "Interp.H"


#include "NamespaceHeader.H"

class SingleNullPoloidalUtil
{
public:

   SingleNullPoloidalUtil( const ParmParse&  pp,
                           const bool        spectral_field,
                           double            RB_toroidal );

   ~SingleNullPoloidalUtil() {;}

   void setRZInterp( const ParmParse&  pp,
                     const string&     interp_method,
                     const FArrayBox&  interp_node_coords,
                     const FArrayBox&  RZ_data );

   void setFieldInterp( const ParmParse&  pp,
                        const string&     interp_method,
                        const FArrayBox&  interp_node_coords,
                        const FArrayBox&  field_data );

   void computeFieldDataSpectral( const FArrayBox&  RZ,
                                  const int         dir,
                                  FArrayBox&        BField,
                                  FArrayBox&        BFieldMag,
                                  FArrayBox&        BFieldDir,
                                  FArrayBox&        gradBFieldMag,
                                  FArrayBox&        curlBFieldDir,
                                  FArrayBox&        BFieldDirdotcurlBFieldDir,
                                  const bool        derived_data_only ) const;

   array<double,3> computeBField(const RealVect& X) const;
   
   array<double,3> computeBUnit(const RealVect& X) const;
   
   void getNodalFieldData( const FArrayBox&  RZ,
                           FArrayBox&        A,
                           FArrayBox&        b,
                           FArrayBox&        Bmag) const;

   void getMagneticFluxFromDCT( const FArrayBox& physical_coordinates,
                                FArrayBox&       magnetic_flux ) const;
   
   double getMagneticFluxFromDCT( const RealVect& physical_coordinate ) const;

   RealVect mappedCoordNewton( const RealVect&  X,
                               const RealVect&  xi_initial,
                               const IntVect&   iv_initial ) const;

   void computeFieldDataFromMappingFile( const FArrayBox&  Xi,
                                         const int         dir,
                                         FArrayBox&        BField,
                                         FArrayBox&        BFieldMag,
                                         FArrayBox&        BFieldDir,
                                         FArrayBox&        gradBFieldMag,
                                         FArrayBox&        curlBFieldDir,
                                         FArrayBox&        BFieldDirdotcurlBFieldDir,
                                         const bool        derived_data_only ) const;

   double interpolate(int component, int mode, const RealVect& xi) const
       {return m_RZ_interp->interpolate(component, mode, xi);}

   void interpolate(const int src_component, const int dst_component, const int mode, const FArrayBox& xi,
                    const Box& box, FArrayBox& dXdxi) const
       {m_RZ_interp->interpolate(src_component, dst_component, mode, xi, box, dXdxi);}

   IntVect getMappingDims() const {return IntVect(m_RZ_interp->dim(0), m_RZ_interp->dim(1));}

   double getRmin() const {return m_Rmin;}

   double getRmax() const {return m_Rmax;}

   double getZmin() const {return m_Zmin;}

   double getZmax() const {return m_Zmax;}

   RealVect getXpoint() const {return m_Xpoint;}

   RealVect getMagAxis() const {return m_magAxis;}

   enum field_vars {RBR_VAR, RBZ_VAR, RB_VAR, UNIT_BR_VAR, UNIT_BPHI_VAR, UNIT_BZ_VAR, NUM_FIELD_VARS};
   enum rz_vars {R_VAR, Z_VAR};

private:

   FArrayBox m_psi_coefs;
   double m_Rmin, m_Rmax, m_Zmin, m_Zmax;
   int m_psi_expansion_order;
   double m_RB_toroidal;
   RealVect m_Xpoint;
   RealVect m_magAxis;
   Interp* m_RZ_interp;
   Interp* m_field_interp;
};

#include "NamespaceFooter.H"

#endif
