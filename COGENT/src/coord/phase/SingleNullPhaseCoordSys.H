#ifndef _SINGLENULLPHASECOORDSYS_H_
#define _SINGLENULLPHASECOORDSYS_H_

#include "PhaseCoordSys.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "SingleNullCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"


class SingleNullPhaseCoordSys : public PhaseCoordSys
{
public:

   /**
      \name Constructors, destructor and defines
   */
   /*@{*/

   ///
   /** Default constructor.
    */
   SingleNullPhaseCoordSys();

   SingleNullPhaseCoordSys( ParmParse&                     pp,
                            const CFG::SingleNullCoordSys& mag_coords,
                            const VEL::VelCoordSys&        vel_coords,
                            const Vector<ProblemDomain>&   domains );

   ///
   /** Destructor.
    */
   virtual ~SingleNullPhaseCoordSys() {;}

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
   virtual RealVect blockRemappingGeneral(int a_nDst,
                                          const RealVect& a_xiSrc,
                                          int a_nSrc) const;

   virtual IntVect getDecomposition(int block) const;

   virtual void defineBoundaries();

   int blockType(int block) const {return ((CFG::SingleNullBlockCoordSys*)m_config_coords.getCoordSys(block))->blockType();}

   virtual std::string type() const { return "SingleNull"; }

   virtual bool isType( const std::string& a_type ) const
   {
      return ( a_type.compare( this->type() )==0 );
   }

private:

   void defineBoundaries8();

   void defineBoundaries10();

   void getDecompositionParams( ParmParse& pp );

   const CFG::SingleNullCoordSys& m_config_coords;

   IntVect m_decomp_mcore;
   IntVect m_decomp_lcore;
   IntVect m_decomp_rcore;
   IntVect m_decomp_lpf;
   IntVect m_decomp_rpf;
   IntVect m_decomp_mcsol;
   IntVect m_decomp_lcsol;
   IntVect m_decomp_rcsol;
   IntVect m_decomp_lsol;
   IntVect m_decomp_rsol;
};


#include "NamespaceFooter.H"

#endif
