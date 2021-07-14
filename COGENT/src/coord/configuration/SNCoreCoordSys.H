#ifndef SNCORECOORDSYS_H
#define SNCORECOORDSYS_H

#include "MagCoordSys.H"
#include "SNCoreBlockCoordSys.H"
#include "SNCoreBlockCoordSysModel.H"
#include "SNCoreEllipticOpBC.H"
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

   int findBlockNearXpt( const POL::RealVect& X_pol ) const;

   virtual void defineBoundaries();

   virtual IntVect getDecomposition(int block) const;

   virtual const EllipticOpBC& getDivergenceCleaningBC() const {return *m_divergence_cleaning_bc;}

   virtual void defineEllipticOpBC( EllipticOpBC& a_bc ) const;

   int blockType(int a_block) const {return ((SNCoreBlockCoordSys*)m_coordSysVect[a_block])->poloidalBlock();}

   int numPoloidalBlocks() const {return m_num_poloidal_blocks;}

   int numToroidalSectors() const {return m_num_toroidal_sectors;}

   int toroidalBlockNumber(const int block_number) const {return block_number / m_num_poloidal_blocks;}

   // Object ParmParse name.
   static const std::string pp_name;

   virtual std::string type() const { return "SNCore"; }

#if CFG_DIM ==3
   /// toroidal block mapping conversion function
   /** Given a_xiSrc in mapped space of block a_nSrc,
    and a_side, returns coordinate in mapped space of the
    destination block (determined from a_side argument)
    and the  second-order interpolation stencil.
    */
   
   virtual void toroidalBlockRemapping( IntVect&               ivDst,
                                        int&                   nDst,
                                        Vector<Real>&          interpStecil,
                                        Vector<int>&           interpStecilOffsets,
                                        const RealVect&        xiSrc,
                                        const int              nSrc,
                                        const Side::LoHiSide&  side ) const;
   
#endif

   virtual Vector<RealVect> displacements(const Vector<RealVect>&   dstCoords,
                                          const Vector<int>&        dstBlocks,
                                          const RealVect&           srcCoords,
                                          int                       srcBlock) const;

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
   
   bool m_model_geometry;

   int m_numcells_core_radial;
   int m_numcells_mcore_poloidal;
   int m_numcells_lcore_poloidal;
   int m_numcells_rcore_poloidal;
   int m_numcells_toroidal;
#if CFG_DIM==3
   double m_toroidal_width;
#endif
   IntVect m_decomp_mcore;
   IntVect m_decomp_lcore;
   IntVect m_decomp_rcore;

   int m_num_poloidal_blocks;
   int m_num_toroidal_sectors;

   POL::RealVect m_Xpoint;
   double m_xpoint_radius;
   double m_ray_angle[SNCoreBlockCoordSys::NUM_SNCORE_BLOCKS+1];
   POL::RealVect m_ray[SNCoreBlockCoordSys::NUM_SNCORE_BLOCKS+1];
   int m_sector_block[SNCoreBlockCoordSys::NUM_SNCORE_BLOCKS];

   SNCoreEllipticOpBC* m_divergence_cleaning_bc;

   Vector<MagBlockCoordSys *> m_coord_vec;
};



#include "NamespaceFooter.H"

#endif
