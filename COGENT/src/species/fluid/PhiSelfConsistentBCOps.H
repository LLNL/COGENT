#ifndef _PHISELFCONSISTENTBCOPS_H_
#define _PHISELFCONSISTENTBCOPS_H_

#include "PhiOps.H"

#include "NamespaceHeader.H"

/// PhiSelfConsistentBCOps interface.
/**
 * PhiSelfConsistentBCOps variable representation.
 */
class PhiSelfConsistentBCOps
   : public PhiOps
{
public:
   
   /// Constructor.
   /**
       */
   PhiSelfConsistentBCOps( const string&      pp_prefix,
                           const MagGeom&     geometry );

   /// Copy Constructor.
   /**
    * @param[in] foo PhiSelfConsistentBCOps from which to initialize.
    */
   PhiSelfConsistentBCOps( const PhiSelfConsistentBCOps& foo );

   /// Destructor.
   /**
    */
   virtual ~PhiSelfConsistentBCOps() {;}

   double getRadialGKPDivergenceAverageLo() const {return m_radial_gkp_divergence_average_lo;}

   double getRadialGKPDivergenceAverageHi() const {return m_radial_gkp_divergence_average_hi;}

   virtual void updateEField( EMFields&                         EM_fields,
                              const bool                        update_potential,
                              const PS::KineticSpeciesPtrVect&  kinetic_species,
                              const FluidSpeciesPtrVect&        fluid_species,
                              const PS::ScalarPtrVect&          scalars,
                              EllipticOpBC&                     bc,
                              const bool                        inject_Efield,
                              const bool                        initial_time );
protected:

   LevelData<FArrayBox> m_radial_gkp_divergence_average;
   double m_radial_gkp_divergence_average_lo;
   double m_radial_gkp_divergence_average_hi;

private:

};


#include "NamespaceFooter.H"

#endif
