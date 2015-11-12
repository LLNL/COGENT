#include "GKBoltzmannRelation.H"
#include "GKBoltzmannRelationF_F.H"

#include "NamespaceHeader.H"



inline Real
pointwiseVal( Real denavg, Real potential, Real temperature )
{
   return denavg * exp( potential / temperature );
}



inline Real
pointwiseDeriv( Real denavg, Real potential, Real temperature )
{
   return pointwiseVal( denavg, potential, temperature ) / temperature;
}



void
GKBoltzmannRelation::updateDensity( const LevelData<FArrayBox>& a_phi,
                                    BoltzmannElectron&      a_boltz_elect ) const
{
   LevelData<FArrayBox>& ne = a_boltz_elect.numberDensity();
   LevelData<FArrayBox>& Te = a_boltz_elect.temperature();

   const DisjointBoxLayout& grids = ne.getBoxes();

   DataIterator dit = ne.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FORT_COMPUTE_ELECTRON_DENSITY(CHF_BOX(grids[dit]),
                                    CHF_CONST_FRA1(a_phi[dit],0),
                                    CHF_CONST_FRA1(Te[dit],0),
                                    CHF_FRA1(ne[dit],0));
   }
}



void
GKBoltzmannRelation::phiDerivative( const BoltzmannElectron& a_boltz_elect,
                                    LevelData<FArrayBox>&        a_derivative ) const
{
   const LevelData<FArrayBox>& ne = a_boltz_elect.numberDensity();
   const LevelData<FArrayBox>& Te = a_boltz_elect.temperature();

   DataIterator dit = a_derivative.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& this_derivative = a_derivative[dit];
      this_derivative.copy(ne[dit]);
      this_derivative.divide(Te[dit]);
   }
}



#include "NamespaceFooter.H"
