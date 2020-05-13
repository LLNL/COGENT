#ifndef _FLUXBC_H_
#define _FLUXBC_H_

#include "BoundaryBoxLayout.H.multidim"
#include "KineticSpeciesBC.H"

#ifdef CH_MULTIDIM

#include "NamespaceHeader.H"


/// Logical BC class
/**
 * A utility class for performing operation with logical sheath BCs.
 */
class FluxBC
{
  public:

   /// Constructor
   /** Constructor
    * @param[in] coords  the boundary box layout
    */
   FluxBC( const BoundaryBoxLayoutPtrVect&    bdry_layout,
           const ParmParse&                   pp );

   /// Destructor
   /** Destructor
    */
   ~FluxBC() {;}

      
   void applyBC(KineticSpecies& species,
                const LevelData<FluxBox>& velocity,
                const CFG::LevelData<CFG::FArrayBox>& phi,
                const Real& time) const;
   
   static const string pp_name;
   
private:

   void fillGhosts(KineticSpecies&           species,
                   const LevelData<FluxBox>& mapped_velocity,
                   const Vector<Real>&       outflow_moments,
                   const Vector<Real>&       pert_moments) const;
   
   void computeOutflowMoments(Vector<Real>&                 moments,
                              const LevelData<FArrayBox>&   dfn,
                              const LevelData<FluxBox>&     mapped_velocity,
                              const LevelData<FArrayBox>&   phi_inj,
                              const KineticSpecies&         species,
                              const BoundaryBoxLayoutPtr&   bdry_layout) const;
   
   void computeInflowFluxMatrix(Vector<Real>&               a_M,
                                Vector<Real>&               a_background_fluxes,
                                const Vector<Real>&         a_outflow_moms,
                                const LevelData<FluxBox>&   a_mapped_velocity,
                                const LevelData<FArrayBox>& a_phi_inj,
                                const KineticSpecies&       a_species,
                                const BoundaryBoxLayoutPtr& a_bdry_layout) const;
   
   void solveTrilinearSystem(Vector<Real>&         x,
                             const Vector<Real>&   M,
                             const Vector<Real>&   b) const;
   
   void normalizeFluxes (Vector<Real>& fluxes) const;
   
   void parseParameters(const ParmParse& pp);
   

   const BoundaryBoxLayoutPtrVect m_bdry_layouts;
   Vector<Real> m_total_fluxes;
   mutable LevelData<FArrayBox> m_BstarParallel;
   enum MomType {DENSITY,VPAR,TEMPERATURE,PART_FLUX,MOM_FLUX,ENERGY_FLUX,AREA,NUM_MOMS};
};

#include "NamespaceFooter.H"

#endif /* CH_MULTIDIM */

#endif
