#ifndef _EFIELD_H_
#define _EFIELD_H_

#include "CFGVar.H"
#include "GKPoisson.H"
#include "BoltzmannElectron.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "GKTimeIntegration.H"
#include "KineticSpecies.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H"
namespace PS = PS_NAMESPACE;


/// EField interface.
/**
 * EField variable representation.
 */
class EField
   : public CFGVar
{
public:
   
   /// Constructor.
   /**
       */
   EField( const string&       pp_prefix,
           const std::string&  name,
           const MagGeom&      geometry,
           const IntVect&      ghost_vect )
      : CFGVar(pp_prefix, name, geometry),
        m_poisson(NULL),
        m_defined(false)
   {
      cell_data().define(geometry.gridsFull(), 3, IntVect::Unit);
      face_data().define(geometry.gridsFull(), 3, IntVect::Unit);
   }

   /// Copy Constructor.
   /**
    * @param[in] foo EField from which to initialize.
    */
   EField( const EField& foo );

   /// Destructor.
   /**
    */
   ~EField() {;}

   virtual void define( const double                      larmor,
                        const double                      debye,
                        const PS::KineticSpeciesPtrVect&  kinetic_species,
                        BoltzmannElectron*                boltzmann_electron,
                        const bool                        fixed_efield,
                        const bool                        support_divfree_phase_vel,
                        const int                         cur_step );

   void fillInternalGhosts( LevelData<FArrayBox>& phi ) const
   {
      CH_assert(m_poisson);
      m_poisson->fillInternalGhosts(phi);
   }

   void computeField( const LevelData<FArrayBox>& phi,
                      LevelData<FArrayBox>&       field_cell ) const
   {
      CH_assert(m_poisson);
      m_poisson->computeField(phi, field_cell );
   }

   void computeField( const LevelData<FArrayBox>& phi,
                      LevelData<FluxBox>&         field_face ) const
   {
      CH_assert(m_poisson);
      m_poisson->computeField(phi, field_face );
   }

   virtual void computeEField( const PS::GKState&                state,
                               const PS::KineticSpeciesPtrVect&  kinetic_species,
                               const FluidSpeciesPtrVect&        fluid_species,
                               const PS::ScalarPtrVect&          scalars,
                               LevelData<FArrayBox>&             phi,
                               PotentialBC&                      bc,
                               const bool                        update_potential,
                               const bool                        initial_time );

   void setCoreBC( const double  core_inner_bv,
                   const double  core_outer_bv,
                   PotentialBC&  bc ) const;

   void computePoloidalField( const LevelData<FArrayBox>& phi,                                                                                                                        LevelData<FluxBox>&         field ) const {m_poisson->computePoloidalField(phi, field);}

   bool supportsDivFreePhaseVel() const {return m_phi_node.isDefined();}

   void interpToNodes(const LevelData<FArrayBox>& phi)
     {
        if ( m_phi_node.isDefined() ) {
           m_poisson->interpToNodes(phi, m_phi_node);
        }
     }

   bool secondOrder() const {return m_poisson->secondOrder();}

   LevelData<FArrayBox>& getPhiNode() {return m_phi_node;}

   const LevelData<FArrayBox>& getPhiNode() const {return m_phi_node;}

   void computeQuasiNeutralElectronDensity( LevelData<FArrayBox>&        quasineutral_density,
                                            LevelData<FArrayBox>&        potential,
                                            const PotentialBC&           bc, 
                                            const LevelData<FArrayBox>&  ion_density) const;

   virtual void updateImplicitPotential( LevelData<FArrayBox>&             phi,
                                         const PS::KineticSpeciesPtrVect&  kinetic_species,
                                         const Vector<Real>&               scalar_data,
                                         LevelData<FArrayBox>&             divJperp,
                                         PotentialBC&                      bc,
                                         const Real                        dt );

   virtual
   RefCountedPtr<CFGVar> clone( const bool copy_data = true ) const;

protected:

   void computeIonChargeDensity( LevelData<FArrayBox>&             ion_charge_density,
                                 const PS::KineticSpeciesPtrVect&  species ) const;

   void computeIonMassDensity( LevelData<FArrayBox>&             mass_density,
                               const PS::KineticSpeciesPtrVect&  species ) const;

   void computeTotalChargeDensity( LevelData<FArrayBox>&             charge_density,
                                   const PS::KineticSpeciesPtrVect&  species ) const;

   void computeIonParallelCurrentDensity( LevelData<FArrayBox>&             ion_current_density,
                                          const PS::KineticSpeciesPtrVect&  species ) const;

   LevelData<FArrayBox> m_phi_node;

   BoltzmannElectron* m_boltzmann_electron;
   bool m_fixed_efield;
   GKPoisson* m_poisson;
   bool m_defined;
};


#include "NamespaceFooter.H"

#endif
