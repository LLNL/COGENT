#ifndef RECTANGULARTORUSCOORDSYS_H
#define RECTANGULARTORUSCOORDSYS_H

#include "MagBlockCoordSys.H"
#include "MagCoordSys.H"
#include "RectangularTorusPotentialBC.H"

#include "RectangularTorusDirections.H"

#include "NamespaceHeader.H"

class RectangularTorusCoordSys : public MagCoordSys
{
public:

   /**
      \name Constructors, destructor and defines
   */
   /*@{*/
   RectangularTorusCoordSys( ParmParse&               a_pp,
                   const std::vector<int>&  a_numcells,
                   const std::vector<bool>& a_is_periodic,
                   const std::vector<int>&  a_decomposition );

   ///
   /** Destructor.
    */
   virtual ~RectangularTorusCoordSys();

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


private:

   ProblemDomain m_block_zero_domain;
   IntVect m_decomposition;
   RectangularTorusPotentialBC m_divergence_cleaning_bc;

   void operator = (const RectangularTorusCoordSys&){;}
};



#include "NamespaceFooter.H"

#endif
