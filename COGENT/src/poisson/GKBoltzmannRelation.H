#ifndef _GK_BOLTZMANN_RELATION_H_
#define _GK_BOLTZMANN_RELATION_H_

#include "MagBlockCoordSys.H"

#include "Dimensions.H"

//#undef CH_SPACEDIM
//#define CH_SPACEDIM PDIM
#include "BoltzmannElectron.H"
//#undef CH_SPACEDIM
//#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H"
//namespace PS = PS_NAMESPACE;

/**
 */
class GKBoltzmannRelation
{
public:

   /// Constructor.
   /**
    */
   GKBoltzmannRelation() {;}

   /// Destructor.
   /**
    */
   ~GKBoltzmannRelation() {;}

   /// Updates number density in BoltzmannElectron species.
   /**
    * Given the pointwise, cell-centered potential phi, the number
    * density in BoltzmannElectron species is updated using
    * n_e = exp(phi / Te)
    * @param[in] phi the pointwise, cell-centered potential field.
    * @param[out] boltz_elect a Boltzmann electron species.
    */
   void updateDensity( const LevelData<FArrayBox>& phi,
                       BoltzmannElectron& boltz_elect ) const;

   /// Computes the derivative of the Boltzmann relation w.r.t. potential.
   /**
    * Given the Boltzmann species, the derivative of the Boltzmann
    * relation with respect to potential is computed using the current
    * Boltzmann species electron denisty.
    *
    * @param[in] boltz_elect a Boltzmann electron species.
    * @param[out] derivative derivative of the Boltzmann relation.
    */
   void phiDerivative( const BoltzmannElectron&    boltz_elect,
                       LevelData<FArrayBox>&           derivative ) const;

};

#include "NamespaceFooter.H"

#endif
