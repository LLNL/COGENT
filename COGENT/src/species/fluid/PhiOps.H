#ifndef _PHIOPS_H_
#define _PHIOPS_H_

#include "CFGVars.H"
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


/// PhiOps interface.
/**
 * PhiOps interface
 */
class PhiOps
{
public:
   
   /// Constructor.
   /**
       */
   PhiOps( const string&    a_pp_prefix,
           const MagGeom&   a_geometry );

   /// Copy Constructor.
   /**
    * @param[in] foo PhiOps from which to initialize.
    */
   PhiOps( const PhiOps& foo );

   /// Destructor.
   /**
    */
   virtual ~PhiOps()
   {
      if (m_poisson) delete m_poisson;
   }

   virtual void define( const double                      larmor,
                        const double                      debye,
                        const PS::KineticSpeciesPtrVect&  kinetic_species,
                        BoltzmannElectron*                boltzmann_electron,
                        const bool                        fixed_efield,
                        const bool                        support_divfree_phase_vel,
                        const int                         potential_ghosts,
                        const int                         cur_step );

   void defineEField( EMFields&  EM_fields ) const;

   void definePotential( const LevelData<FArrayBox>&  phi,
                         const bool                   include_phi_node,
                         EMFields&                    EM_fields ) const;

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

   void setCoreBC( const double  core_inner_bv,
                   const double  core_outer_bv,
                   EllipticOpBC& bc ) const;

   void computePoloidalField( const LevelData<FArrayBox>& phi,                                                                                                           LevelData<FluxBox>&         field ) const
               {m_poisson->computePoloidalField(phi, field);}

   void interpToNodes(const LevelData<FArrayBox>& phi,
                      LevelData<FArrayBox>&       phi_node)
     {
        if ( phi_node.isDefined() ) {
           m_poisson->interpToNodes(phi, phi_node);
        }
     }

   //   bool secondOrder() const {return m_poisson->secondOrder();}

   void computeQuasiNeutralElectronDensity( LevelData<FArrayBox>&        quasineutral_density,
                                            LevelData<FArrayBox>&        potential,
                                            EllipticOpBC&                bc, 
                                            const LevelData<FArrayBox>&  ion_density) const;

   virtual void updateEField( EMFields&                         EM_fields,
                              const bool                        update_potential,
                              const PS::KineticSpeciesPtrVect&  kinetic_species,
                              const FluidSpeciesPtrVect&        fluid_species,
                              const PS::ScalarPtrVect&          scalars,
                              EllipticOpBC&                     bc,
                              const bool                        inject_field,
                              const bool                        initial_time );

protected:

   void computeIonChargeDensity( LevelData<FArrayBox>&             ion_charge_density,
                                 const PS::KineticSpeciesPtrVect&  species) const;

   void computeIonChargeDensity( LevelData<FArrayBox>&               a_ion_charge_density,
                                 const PS::KineticSpeciesPtrVect&    a_kinetic_species_phys,
                                 const FluidSpeciesPtrVect&          a_fluid_species_comp) const;

   void computePolarizationMassDensity(LevelData<FArrayBox>&             mass_density,
                                       const PS::KineticSpeciesPtrVect&  kinetic_species,
                                       const FluidSpeciesPtrVect&        fluid_species) const;

   void computeTotalChargeDensity( LevelData<FArrayBox>&             charge_density,
                                   const PS::KineticSpeciesPtrVect&  kinetic_species,
                                   const FluidSpeciesPtrVect&        fluid_species) const;

   void computeIonParallelCurrentDensity( LevelData<FArrayBox>&             ion_current_density,
                                          const PS::KineticSpeciesPtrVect&  species ) const;

   void inject( PS::LevelData<PS::FluxBox>&       injected_data,
                const PS::KineticSpeciesPtrVect&  kinetic_species,
                const LevelData<FluxBox>&         data_face,
                const LevelData<FArrayBox>&       data_cell );

   const MagGeom& m_geometry;

   BoltzmannElectron* m_boltzmann_electron;
   bool m_fixed_efield;
   GKPoisson* m_poisson;
   bool m_FLR_charge_density;
   bool m_defined;
   bool m_apply_harm_filtering;
   int m_harm_filtering_dir;
   bool m_include_pol_dens_i;
   bool m_include_pol_dens_e;
   bool m_consistent_upper_bc_only;
   bool m_consistent_lower_bc_only;
   bool m_fixed_krho2;
   Real m_krho2_value;
   
private:
   
   void applyBoussinesqCorrections(LevelData<FArrayBox>& total_charge_density,
                                   const LevelData<FArrayBox>& ion_mass_density) const;
   
   void parseParameters( ParmParse& pp);
   
};


#include "NamespaceFooter.H"

#endif
