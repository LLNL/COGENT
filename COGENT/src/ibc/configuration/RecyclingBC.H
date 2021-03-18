#ifndef _RECYCLINGBC_H_
#define _RECYCLINGBC_H_

#include "BoundaryBoxLayout.H.multidim"
#include "MagGeom.H"

#include "NamespaceHeader.H"


/// Logical BC class
/**
 * A utility class for performing operation with logical sheath BCs.
 */
class RecyclingBC
{
  public:

   /// Constructor
   /** Constructor
    * @param[in] coords  the boundary box layout
    */
   RecyclingBC(const BoundaryBoxLayoutPtr& bdry_layout,
               const ParmParse&            pp );

   /// Destructor
   /** Destructor
    */
   ~RecyclingBC() {;}

      
   void applyBC(LevelData<FluxBox>&         outflux,
                const LevelData<FluxBox>&   influx,
                const Real                  R,
                const Real                  time ) const;

   
   void applyBC(EllipticOpBC&             elliptic_bc,
                const LevelData<FluxBox>& face_dist,
                const Real                R,
                const Real                time ) const;
   
   void applyBC(LevelData<FArrayBox>&       a_dst,
                const LevelData<FArrayBox>& a_src,
                const Real                  a_R,
                const Real                  a_time ) const;
   
   static const string pp_name;
   
private:

   void parseParameters(const ParmParse& pp);
   
   const BoundaryBoxLayoutPtr m_bdry_layout;
};

#include "NamespaceFooter.H"

#endif
