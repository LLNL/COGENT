#ifndef MILLERCOORDSYS_H
#define MILLERCOORDSYS_H

#include "MagBlockCoordSys.H"
#include "MagCoordSys.H"
#include "AnnulusPotentialBC.H"

#include "NamespaceHeader.H"

class MillerCoordSys : public MagCoordSys
{
public:

   /**
      \name Constructors, destructor and defines
   */
   /*@{*/
   MillerCoordSys( ParmParse&               a_pp,
                   const std::vector<int>&  a_numcells,
                   const std::vector<bool>& a_is_periodic,
                   const std::vector<int>&  a_decomposition );

   ///
   /** Destructor.
    */
   virtual ~MillerCoordSys();

   /*@}*/

   /// block mapping conversion function
   /** given xi_k in block k, returns which block (k_valid)
       contains that location in its valid region and its
       location in that block's computational coordinate space.
   */
   virtual void blockRemapping(RealVect&       xi_valid,
                               int&            n_valid,
                               const RealVect& xiSrc,
                               int             nSrc) const;

   virtual void defineBoundaries();

   virtual IntVect getDecomposition(int block) const {return m_decomposition;}

   virtual const PotentialBC& getDivergenceCleaningBC() const {return m_divergence_cleaning_bc;}

   bool isPeriodic( const int dir ) const;

protected:

   virtual void postProcessMetricData( LevelData<FArrayBox>& data ) const;

   virtual void postProcessMetricData( LevelData<FluxBox>& data ) const;

   void spreadRadially( LevelData<FArrayBox>& data ) const;

   void spreadRadially( LevelData<FluxBox>& data ) const;

private:

   ProblemDomain m_block_zero_domain;
   bool m_spread_radially;
   IntVect m_decomposition;
   AnnulusPotentialBC m_divergence_cleaning_bc;

   void operator = (const MillerCoordSys&){;}
};



#include "NamespaceFooter.H"

#endif
