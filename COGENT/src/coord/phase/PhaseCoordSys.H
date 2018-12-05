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

   virtual IntVect getDecomposition(int block) const = 0;

   virtual void defineBoundaries() = 0;

   int numBlocks() const {return m_mappingBlocks.size();}

   virtual std::string type() const = 0;

   virtual bool isType( const std::string& name ) const = 0;

protected:

   const RefCountedPtr<CFG::MagCoordSys> m_mag_coords;
   const RefCountedPtr<VEL::VelCoordSys> m_vel_coords;

};


#include "NamespaceFooter.H"

#endif
