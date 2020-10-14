#ifndef _PHASECOORDSYS_H_
#define _PHASECOORDSYS_H_

#include "MultiBlockCoordSys.H"
#include "PhaseBlockCoordSys.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "MagCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

#include <string>

class PhaseCoordSys : public MultiBlockCoordSys
{
public:

   /**
      \name Constructors, destructor and defines
   */
   /*@{*/

   ///
   /** Default constructor.
    */
   PhaseCoordSys();

   PhaseCoordSys( const RefCountedPtr<CFG::MagCoordSys>&  mag_coords,
                  const RefCountedPtr<VEL::VelCoordSys>&  vel_coords,
                  const Vector<ProblemDomain>&            domains );

   ///
   /** Destructor.
    */
   virtual ~PhaseCoordSys();

   /*@}*/

   /// block mapping conversion function
   /** given xi_k in block k, returns which block (k_valid)
       contains that location in its valid region and its
       location in that block's computational coordinate space.
   */
   virtual void blockRemapping(RealVect&       xi_valid,
                               int&            n_valid,
                               const RealVect& xiSrc,
                               int             nSrc) const = 0;

   /// Return the mapping which corresponds to the block with the given index.
   virtual const PhaseBlockCoordSys* getCoordSys(int a_n) const
           {return (const PhaseBlockCoordSys*)m_coordSysVect[a_n];}

   /// Return the mapping which corresponds to the given box.
   /// \param a_box The Box for which the coordinate system is retrieved. This
   ///              Box must be contained completely within a block or this
   ///              method will fail.
   virtual const PhaseBlockCoordSys* getCoordSys(const Box& a_box) const
           { return (const PhaseBlockCoordSys*)m_coordSysVect[whichBlock(a_box)]; }

   virtual IntVect getDecomposition(int block) const = 0;

   virtual void defineBoundaries() = 0;

   int numBlocks() const {return m_mappingBlocks.size();}

   virtual std::string type() const = 0;

   virtual bool isType( const std::string& name ) const = 0;

   virtual bool containsPhysicalBoundary( int                    block_number,
                                          int                    dir,
                                          const Side::LoHiSide&  side ) const;
protected:

   const RefCountedPtr<CFG::MagCoordSys> m_mag_coords;
   const RefCountedPtr<VEL::VelCoordSys> m_vel_coords;

};


#include "NamespaceFooter.H"

#endif
