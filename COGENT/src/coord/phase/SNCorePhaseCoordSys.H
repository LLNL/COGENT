#ifndef _SNCOREPHASECOORDSYS_H_
#define _SNCOREPHASECOORDSYS_H_

#include "PhaseCoordSys.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "SNCoreCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"


class SNCorePhaseCoordSys : public PhaseCoordSys
{
public:

   /**
      \name Constructors, destructor and defines
   */
   /*@{*/

   ///
   /** Default constructor.
    */
   SNCorePhaseCoordSys();

   SNCorePhaseCoordSys( ParmParse&                   pp,
                        const CFG::SNCoreCoordSys&   mag_coords,
                        const VEL::VelCoordSys&      vel_coords,
                        const Vector<ProblemDomain>& domains );

   ///
   /** Destructor.
    */
   virtual ~SNCorePhaseCoordSys() {;}

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

   /// general block mapping conversion function
   /** Given a_xiSrc in mapped space of block a_nSrc,
       and block index a_nDst,
       returns coordinates in mapped space of block a_nDst.
   */
   virtual RealVect blockRemappingGeneral(int             nDst,
                                          const RealVect& xiSrc,
                                          int             nSrc) const;

   virtual IntVect getDecomposition(int block) const;

   virtual void defineBoundaries();
   void defineBoundaries2();
   void defineBoundaries3();

   int blockType(int block) const {return ((CFG::SNCoreBlockCoordSys*)m_config_coords.getCoordSys(block))->blockType();}

   virtual std::string type() const { return "SNCore"; }

   virtual bool isType( const std::string& a_type ) const
   {
      return ( a_type.compare( this->type() )==0 );
   }

private:

   void getDecompositionParams( ParmParse& pp );

   const CFG::SNCoreCoordSys& m_config_coords;

   IntVect m_decomp_mcore;
   IntVect m_decomp_lcore;
   IntVect m_decomp_rcore;
};


#include "NamespaceFooter.H"

#endif
