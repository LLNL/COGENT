#ifndef _PHIAMPEREOPS_H_
#define _PHIAMPEREOPS_H_

#include "PhiOps.H"

#include "NamespaceHeader.H"

/// PhiAmpereOps interface.
/**
 * PhiAmpereOps variable representation.
 */
class PhiAmpereOps
   : public PhiOps
{
public:
   
   /// Constructor.
   /**
       */
   PhiAmpereOps( const string&      pp_prefix,
                 const MagGeom&     geometry );

   /// Copy Constructor.
   /**
    * @param[in] foo PhiAmpereOps from which to initialize.
    */
   PhiAmpereOps( const PhiAmpereOps& foo );

   /// Destructor.
   /**
    */
   virtual ~PhiAmpereOps() {;}

   virtual void define( const double                      larmor,
                        const double                      debye,
                        const PS::KineticSpeciesPtrVect&  kinetic_species,
                        BoltzmannElectron*                boltzmann_electron,
                        const bool                        fixed_efield,
                        const bool                        support_divfree_phase_vel,
                        const int                         phi_ghosts,
                        const int                         cur_step );

   void setErAverage( const LevelData<FArrayBox>& Er_average ) {m_Er_average_cell.define(Er_average);}

   void setErAverage( const LevelData<FluxBox>& Er_average ) {m_Er_average_face.define(Er_average);}

   void setETilde( const LevelData<FArrayBox>& E_tilde ) {m_E_tilde_cell.define(E_tilde);}

   void setETilde( const LevelData<FluxBox>& E_tilde ) {m_E_tilde_face.define(E_tilde);}

   const LevelData<FArrayBox>& getErAverageCell() const {return m_Er_average_cell;}

   LevelData<FArrayBox>& getErAverageCell() {return m_Er_average_cell;}

   const LevelData<FluxBox>& getErAverageFace() const {return m_Er_average_face;}

   LevelData<FluxBox>& getErAverageFace() {return m_Er_average_face;}

   const LevelData<FArrayBox>& getETildeCell() const {return m_E_tilde_cell;}

   const LevelData<FluxBox>& getETildeFace() const {return m_E_tilde_face;}

   const LevelData<FArrayBox>& getRadialGKPDivergenceAverage() const {return m_radial_gkp_divergence_average;}

   double getRadialGKPDivergenceAverageLo() const {return m_radial_gkp_divergence_average_lo;}

   double getRadialGKPDivergenceAverageHi() const {return m_radial_gkp_divergence_average_hi;}

   void updateAverageAndPerturbation( const Vector<Real>&               scalar_data,
                                      const LevelData<FArrayBox>&       Er,
                                      const PS::KineticSpeciesPtrVect&  species_phys,
                                      const Real                        time );

   virtual void updateEField( EMFields&                         EM_fields,
                              const bool                        update_potential,
                              const PS::KineticSpeciesPtrVect&  kinetic_species,
                              const FluidSpeciesPtrVect&        fluid_species,
                              const PS::ScalarPtrVect&          scalars,
                              EllipticOpBC&                     bc,
                              const bool                        inject_Efield,
                              const bool                        initial_time );

protected:

   void updateAveragedEfield( LevelData<FluxBox>&          Er_average_face,
                              LevelData<FArrayBox>&        Er_average_cell,
                              const LevelData<FArrayBox>&  Er,
                              const double                 Er_lo,
                              const double                 Er_hi ) const;

   void updateEfieldPoloidalVariation( const LevelData<FluxBox>&    E_tilde_mapped_face,
                                       const LevelData<FArrayBox>&  E_tilde_mapped_cell,
                                       LevelData<FluxBox>&          E_tilde_phys_face,
                                       LevelData<FArrayBox>&        E_tilde_phys_cell,
                                       LevelData<FArrayBox>&        phi_tilde_fs_average,
                                       double&                      phi_tilde_fs_average_lo,
                                       double&                      phi_tilde_fs_average_hi,
                                       const LevelData<FArrayBox>&  radial_gkp_divergence_average,
                                       const double                 radial_gkp_divergence_average_lo,
                                       const double                 radial_gkp_divergence_average_hi ) const;

   void computePhiTilde( const PS::KineticSpeciesPtrVect&  kinetic_species,
                         const BoltzmannElectron&          ne,
                         LevelData<FArrayBox>&             phi_tilde ) const;

   LevelData<FArrayBox> m_Er_average_cell;
   LevelData<FluxBox> m_Er_average_face;
   LevelData<FArrayBox> m_E_tilde_cell;
   LevelData<FluxBox> m_E_tilde_face;

   LevelData<FArrayBox> m_radial_gkp_divergence_average;
   double m_radial_gkp_divergence_average_lo;
   double m_radial_gkp_divergence_average_hi;

   bool m_cold_electrons;
};


#include "NamespaceFooter.H"

#endif
