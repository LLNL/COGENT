#ifndef _INSULATORCONDUCTORBC_H_
#define _INSULATORCONDUCTORBC_H_

#include "TimeFunction.H"
#include "ParsingCore.H"
#include "LevelData.H"
#include "MagGeom.H"
#include "EdgeDataBox.H"

#include "NamespaceHeader.H"


/// Insulator/Conductor BC class
/**
 * A utility class for performing operation with insulator/conductor BCs.
 */
class InsulatorConductorBC
{
  public:

   /// Constructor
   //
   InsulatorConductorBC( const string&          a_species_name,
                         const string&          a_variable_name,
                         const Vector<string>&  a_bdry_name );

   /// Destructor
   //
   ~InsulatorConductorBC() {;}

   void defineInsulatorConductorBinary( const LevelData<FArrayBox>&      a_data,
                                        const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
                                        const MagGeom&                   a_geometry ) const;
  
   void defineInsulatorConductorBinary( const LevelData<FluxBox>&        a_data,
                                        const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
                                        const MagGeom&                   a_geometry ) const;
   
   void defineInsulatorConductorBinary( const LevelData<EdgeDataBox>&    a_data,
                                        const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
                                        const MagGeom&                   a_geometry ) const;
   
   void applyBC( LevelData<FArrayBox>&      a_data,
           const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
           const Vector<string>&            a_all_bc_type,
           const MagGeom&                   a_geometry,   
           const Real                       a_time ) const;
   
   void applyFluxBC( LevelData<FluxBox>&    a_data,
           const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
           const Vector<string>&            a_all_bc_type,
           const MagGeom&                   a_geometry,   
           const Real                       a_time ) const;
   
   void setInsulatorBC( LevelData<FArrayBox>&      a_dst,
                  const LevelData<FArrayBox>&      a_src,
                  const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
                  const Vector<string>&            a_all_bc_type,
                  const MagGeom&                   a_geometry,   
                  const Real                       a_time ) const;
   
   void applyEdgeBC( LevelData<EdgeDataBox>&    a_data,
               const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
               const Vector<string>&            a_all_bc_type,
               const MagGeom&                   a_geometry ) const;
   
private:

   void parseParameters(const Vector<string>&  a_bdry_name);
  
   std::string m_species_name;
   std::string m_variable_name;
 
   std::vector<int> m_bdry_is_conductor;
   std::vector<int> m_bdry_is_insulator_conductor;
   RealVect m_Xmin_insulator;
   RealVect m_Xmax_insulator;
   mutable LevelData<FArrayBox> m_ICbinary;      // insulator conductor binary
   mutable LevelData<FluxBox> m_ICbinary_fc;     // insulator conductor binary on faces
   mutable LevelData<EdgeDataBox> m_ICbinary_ec; // insulator conductor binary on edges
   RefCountedPtr<TimeFunction> m_timeFunction;
      
};

#include "NamespaceFooter.H"

#endif
