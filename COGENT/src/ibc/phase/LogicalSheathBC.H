#ifndef _LOGICALSHEATHBC_H_
#define _LOGICALSHEATHBC_H_

#include "BoundaryBoxLayout.H.multidim"
#include "KineticSpeciesBC.H"

#ifdef CH_MULTIDIM

#include "NamespaceHeader.H"


/// Logical BC class
/**
 * A utility class for performing operation with logical sheath BCs.
 */
class LogicalSheathBC
{
  public:

   /// Constructor
   /** Constructor
    * @param[in] coords  the boundary box layout
    */
   LogicalSheathBC( const BoundaryBoxLayoutPtr& bdry_layout,
                    const ParmParse&            pp );

   /// Destructor
   /** Destructor
    */
   ~LogicalSheathBC() {;}

   
   void computeSheathBC(LevelData<FArrayBox>& phi_bc,
                        const KineticSpeciesPtrVect& species,
                        const int& species_index) const;
   
   void fillBoundaryData(LevelData<FArrayBox>& bdry_data,
                         const LevelData<FArrayBox>& dfn) const;
   
   void computeBoundaryIonCurrent(LevelData<FArrayBox>& bdry_ion_current,
                                  const KineticSpeciesPtrVect& species) const;
   
   void applyBC(KineticSpeciesPtrVect& species,
                const int& species_index,
                const LevelData<FluxBox>& velocity,
                const CFG::LevelData<CFG::FArrayBox>& phi ) const;

   /// Sum the src array over velocity space
   /** Sum the src array over velocity space
    *
    * @param[in]  src  Data to be summed over velocity space
    * @param[out] dst  velocity-space sums
    */
   void sum( const LevelData<FArrayBox>& src,
             LevelData<FArrayBox>&       dst ) const;
   
   static const string pp_name;
   
private:

   void parseParameters(const ParmParse& pp);
   
   DisjointBoxLayout m_grids_inj;
   const BoundaryBoxLayoutPtr m_bdry_layout;
   bool m_compute_potential_BC;
   std::string m_sheath_bc_type;
   
};

#include "NamespaceFooter.H"

#endif /* CH_MULTIDIM */

#endif
