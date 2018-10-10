#ifndef _EFIELDSELFCONSISTENTBC_H_
#define _EFIELDSELFCONSISTENTBC_H_

#include "EField.H"

#include "NamespaceHeader.H"

/// EFieldSelfConsistentBC interface.
/**
 * EFieldSelfConsistentBC variable representation.
 */
class EFieldSelfConsistentBC
   : public EField
{
public:
   
   /// Constructor.
   /**
       */
   EFieldSelfConsistentBC( const string&      pp_prefix,
                           const std::string& name,
                           const MagGeom&     geometry,
                           const IntVect&     ghost_vect );

   /// Copy Constructor.
   /**
    * @param[in] foo EFieldSelfConsistentBC from which to initialize.
    */
   EFieldSelfConsistentBC( const EFieldSelfConsistentBC& foo );

   /// Destructor.
   /**
    */
   ~EFieldSelfConsistentBC() {;}

   virtual void computeEField( const PS::GKState&                state,
                               const PS::KineticSpeciesPtrVect&  kinetic_species,
                               const FluidSpeciesPtrVect&        fluid_species,
                               const PS::ScalarPtrVect&          scalars,
                               LevelData<FArrayBox>&             phi,
                               EllipticOpBC&                     bc,
                               const bool                        update_potential,
                               const bool                        initial_time );

   virtual void updateImplicitPotential( LevelData<FArrayBox>&             phi,
                                         const PS::KineticSpeciesPtrVect&  kinetic_species,
                                         const Vector<Real>&               scalar_data,
                                         LevelData<FArrayBox>&             divJperp,
                                         EllipticOpBC&                     bc,
                                         const Real                        dt );

   virtual
   RefCountedPtr<CFGVars> clone( const bool copy_data = true ) const;

   double getRadialGKPDivergenceAverageLo() const {return m_radial_gkp_divergence_average_lo;}

   double getRadialGKPDivergenceAverageHi() const {return m_radial_gkp_divergence_average_hi;}

protected:

   LevelData<FArrayBox> m_radial_gkp_divergence_average;
   double m_radial_gkp_divergence_average_lo;
   double m_radial_gkp_divergence_average_hi;
};


#include "NamespaceFooter.H"

#endif
