#ifndef SLABCOORDSYS_H
#define SLABCOORDSYS_H

#include "MagBlockCoordSys.H"
#include "MagCoordSys.H"
//#include "AnnulusPotentialBC.H"
#include "SlabPotentialBC.H"

#include "NamespaceHeader.H"

class SlabCoordSys : public MagCoordSys
{
public:

   /**
      \name Constructors, destructor and defines
   */
   /*@{*/
   SlabCoordSys( ParmParse&               a_pp,
                   const std::vector<int>&  a_numcells,
                   const std::vector<bool>& a_is_periodic,
                   const std::vector<int>&  a_decomposition );

   ///
   /** Destructor.
    */
   virtual ~SlabCoordSys();

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

//protected:

//   virtual void postProcessMetricData( LevelData<FArrayBox>& data ) const;

//   virtual void postProcessMetricData( LevelData<FluxBox>& data ) const;

private:

   ProblemDomain m_block_zero_domain;
   IntVect m_decomposition;
//   AnnulusPotentialBC m_divergence_cleaning_bc;
   SlabPotentialBC m_divergence_cleaning_bc;

   void operator = (const SlabCoordSys&){;}
};



#include "NamespaceFooter.H"

#endif
