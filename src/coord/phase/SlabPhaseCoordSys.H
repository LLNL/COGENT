#ifndef _SLABPHASECOORDSYS_H_
#define _SLABPHASECOORDSYS_H_

#include "PhaseCoordSys.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "SlabCoordSys.H"
#include "newMappedGridIO.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

class SlabPhaseCoordSys : public PhaseCoordSys
{
public:

   /**
      \name Constructors, destructor and defines
   */
   /*@{*/

   ///
   /** Default constructor.
    */
   SlabPhaseCoordSys();

   SlabPhaseCoordSys( ParmParse&                   pp,
                        const CFG::SlabCoordSys&   mag_coords,
                        const VEL::VelCoordSys&      vel_coords,
                        const Vector<ProblemDomain>& domains );

   ///
   /** Destructor.
    */
   virtual ~SlabPhaseCoordSys() {;}

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

   virtual IntVect getDecomposition(int block) const {return m_decomposition;}

   virtual void defineBoundaries();

   virtual std::string type() const { return "Slab"; }

   virtual bool isType( const std::string& a_type ) const
   {
      return ( a_type.compare( this->type() )==0 );
   }

private:

   void getDecompositionParams( ParmParse& pp );

   IntVect m_decomposition;
      bool m_is_periodic[SpaceDim];
};


#include "NamespaceFooter.H"

#endif
