#ifndef SNCORECOORDSYS_H
#define SNCORECOORDSYS_H

#include "MagCoordSys.H"
#include "SNCoreBlockCoordSys.H"
#include "SNCorePotentialBC.H"
#include "SingleNullCoordSys.H"

#include "NamespaceHeader.H"


class SNCoreCoordSys : public MagCoordSys
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/
  SNCoreCoordSys( ParmParse& pp_grid,
                  ParmParse& pp_geom );

  SNCoreCoordSys( ParmParse&                pp_geom,
                  const SingleNullCoordSys& single_null_coord_sys,
                  const DisjointBoxLayout&  single_null_dbl,
                  DisjointBoxLayout&        dbl );
  ///
  /** Destructor.
   */
  virtual ~SNCoreCoordSys();

  /*@}*/

  /// block mapping conversion function
  /** given xi_k in block k, returns which block (k_valid)
      contains that location in its valid region and its
      location in that block's computational coordinate space.
  */
   virtual void blockRemapping(RealVect& xi_valid,
                               int& n_valid,
                               const RealVect& xiSrc,
                               int nSrc) const;

   int findBlock(const RealVect& X) const;

   virtual void defineBoundaries();

   virtual IntVect getDecomposition(int block) const;

   virtual const PotentialBC& getDivergenceCleaningBC() const {return m_divergence_cleaning_bc;}

   virtual void definePotentialBC( PotentialBC& a_bc ) const;

protected:

   void define( ParmParse& a_pp_geom );

   void defineBoundaries2();

   void defineBoundaries3();

   void blockRemapping2(RealVect& xi_valid,
                        int& n_valid,
                        const RealVect& xiSrc,
                        int nSrc) const;
   void blockRemapping3(RealVect& xi_valid,
                        int& n_valid,
                        const RealVect& xiSrc,
                        int nSrc) const;

   void setXPointNeighborhood();

   void readGridParams( ParmParse& a_pp );

   void printInvalidPointDiagnostics( const int                  nSrc,
                                      const RealVect&            xiSrc,
                                      const RealVect&            X,
                                      const RealVect&            xi_valid,
                                      const int                  n_valid,
                                      const SNCoreBlockCoordSys* src_cs,
                                      const SNCoreBlockCoordSys* valid_cs ) const;

   IntVect lo_mapped_index(int block) const;

   IntVect hi_mapped_index(int block) const;

   bool m_original_two_blocks;

   int m_numcells_core_radial;
   int m_numcells_mcore_poloidal;
   int m_numcells_lcore_poloidal;
   int m_numcells_rcore_poloidal;
   IntVect m_decomp_mcore;
   IntVect m_decomp_lcore;
   IntVect m_decomp_rcore;

   RealVect m_Xpoint;
   double m_xpoint_radius;
   double m_ray_angle[NUM_SNCORE_BLOCKS+1];
   RealVect m_ray[NUM_SNCORE_BLOCKS+1];
   int m_sector_block[NUM_SNCORE_BLOCKS];

   SNCorePotentialBC m_divergence_cleaning_bc;
};



#include "NamespaceFooter.H"

#endif
