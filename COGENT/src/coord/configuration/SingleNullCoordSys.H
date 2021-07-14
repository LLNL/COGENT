#ifndef SINGLENULLCOORDSYS_H
#define SINGLENULLCOORDSYS_H

#include "MagCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "SingleNullEllipticOpBC.H"
#include "SingleNullBlockCoordSys.H"

#include "NamespaceHeader.H"

#define POLOIDAL_BLOCK_SEP 2
#define TOROIDAL_BLOCK_SEP 8

class SingleNullCoordSys : public MagCoordSys
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/
  SingleNullCoordSys( ParmParse& pp_grid,
                      ParmParse& pp_geom );

  ///
  /** Destructor.
   */
  virtual ~SingleNullCoordSys();

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

   int findBlockNearXpt(const POL::RealVect& X_pol ) const;

   virtual void defineBoundaries();

   virtual IntVect getDecomposition(int block) const;

   virtual const EllipticOpBC& getDivergenceCleaningBC() const {return *m_divergence_cleaning_bc;}

   virtual void defineEllipticOpBC( EllipticOpBC& bc ) const;

   virtual void setBoundaryBoxes();

   int blockType(int a_block) const {return ((SingleNullBlockCoordSys*)m_coordSysVect[a_block])->poloidalBlock();}

   int numPoloidalBlocks() const {return m_num_poloidal_blocks;}

   int numToroidalSectors() const {return m_num_toroidal_sectors;}

   int numCellsCoreRadial() const {return m_numcells_core_radial;}

   int numCellsMcorePoloidal() const {return m_numcells_mcore_poloidal;}

   int numCellsLcorePoloidal() const {return m_numcells_lcore_poloidal;}

   int numCellsRcorePoloidal() const {return m_numcells_rcore_poloidal;}
   
   int numCellsToroidal() const {return m_numcells_toroidal;}
   
   bool isModelGeom() const {return m_model_geometry;}

     //ABSOLETE: Extrapolate E-field into the open-field-line region
   void extrapolateEfield( LevelData<FluxBox>&   Er_average_face,
                           LevelData<FArrayBox>& Er_average_cell ) const;
   POL::RealVect m_Xpoint;

   void defineStencilsUe( const DisjointBoxLayout&  grids,
                         const int                 nghost );

   
   const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& getStencilsUe() const {return m_stencil_ue;}
   
   const LayoutData< IntVectSet >& getGhostCellsUe() const {return m_ghostCells_ue;}
   
   // Object ParmParse name.
   static const std::string pp_name;
   
   virtual std::string type() const { return "SingleNull"; }

   virtual bool containsPhysicalBoundary( int                    block_number,
                                          int                    dir,
                                          const Side::LoHiSide&  side ) const;

#if CFG_DIM==3
   virtual void toroidalBlockRemapping( IntVect&               ivDst,
                                        int&                   nDst,
                                        Vector<Real>&          interpStecil,
                                        Vector<int>&           interpStecilOffsets,
                                        const RealVect&        xiSrc,
                                        const int              nSrc,
                                        const Side::LoHiSide&  side) const;

   void getInterpolationCoefficients(Vector<Real>&    coeff,
                                     Vector<int>&     offsets,
                                     const RealVect&  xi0_dst,
                                     const IntVect&   iv0_dst,
                                     const RealVect&  dx_dst,
                                     const int        order) const;

   void getInterpolationCoefficients(Vector<Real>&                  a_coeff,
                                     Vector<int>&                   a_offsets,
                                     const RealVect&                a_X,
                                     const IntVect&                 a_iv0_dst,
                                     const SingleNullBlockCoordSys& a_coords) const;
   
   Real getInterpMappedPoloidalCoord(const RealVect& X,
                                     const RealVect& X_dst_pol_cent,
                                     const RealVect& X_dst_pol_hi,
                                     const RealVect& X_dst_pol_lo,
                                     const RealVect& xi_dst_cent,
                                     const RealVect& dx) const;

   void applyToroidalPeriodicity(RealVect& X,
				 const Side::LoHiSide&  side) const;

   Vector<Box> getToroidalBoundaryBoxes( const int              block_number,
                                         const Side::LoHiSide&  side,
                                         const int              shear_interp_order ) const;
#endif

   virtual Vector<RealVect> displacements(const Vector<RealVect>&   dstCoords,
                                          const Vector<int>&        dstBlocks,
                                          const RealVect&           srcCoords,
                                          int                       srcBlock) const;

   bool sameMagFluxSurf(const int& pol_block_src,
			int&       pol_block) const;


   int poloidalBlockNumber(const int block_number) const {return block_number % m_num_poloidal_blocks;}

   int toroidalBlockNumber(const int block_number) const {return block_number / m_num_poloidal_blocks;}

   int blockNumber(const int poloidal_block_number, const int toroidal_block_number) const 
              {return m_num_poloidal_blocks * toroidal_block_number + poloidal_block_number;}

#if CFG_DIM==3
   double toroidalWidth() const {return m_toroidal_width;}
#endif
   
   
   bool isCORE(const int block_number) const;
   bool isSOL(const int block_number) const;
   bool isPF(const int block_number) const;
   
protected:

#if CFG_DIM==3
   bool fixPoloidalBlock( int&      num_block,
                          IntVect&  iv ) const;
#endif
   
   void defineBoundaries8();

   void defineBoundaries10();

   void blockRemapping8(RealVect& xi_valid,
                        int& n_valid,
                        const RealVect& xiSrc,
                        int nSrc) const;

   void blockRemapping10(RealVect& xi_valid,
                        int& n_valid,
                        const RealVect& xiSrc,
                        int nSrc) const;

   void setXPointNeighborhood();
   
   void defineStencilsUe8( const DisjointBoxLayout&  grids,
                          const int                 nghost );


   void readParams( ParmParse& a_pp_grid,
                    ParmParse& a_pp_geom );

   void printInvalidPointDiagnostics( const int                      nSrc,
                                      const RealVect&                xiSrc,
                                      const RealVect&                X,
                                      const RealVect&                xi_valid,
                                      const int                      n_valid,
                                      const SingleNullBlockCoordSys* src_cs,
                                      const SingleNullBlockCoordSys* valid_cs ) const;
    
   

   void checkGridConsistency() const;

   IntVect lo_mapped_index(int block_number) const;

   IntVect hi_mapped_index(int block_number) const;

   string m_geom_type;

   bool m_model_geometry;

   bool m_original_eight_blocks;

   int m_numcells_core_radial;
   int m_numcells_mcore_poloidal;
   int m_numcells_lcore_poloidal;
   int m_numcells_rcore_poloidal;
   int m_numcells_sol_radial;
   int m_numcells_lpf_poloidal;
   int m_numcells_rpf_poloidal;
   int m_numcells_pf_radial;
   int m_numcells_toroidal;

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

   int m_num_poloidal_blocks;
   int m_num_toroidal_sectors;

   double m_xpoint_radius;
   double m_toroidal_width;
   int    m_shear_interp_order;
        
   double m_ray_angle[SingleNullBlockCoordSys::NUM_POLOIDAL_BLOCKS+1];
   POL::RealVect m_ray[SingleNullBlockCoordSys::NUM_POLOIDAL_BLOCKS];
   int m_sector_block[SingleNullBlockCoordSys::NUM_POLOIDAL_BLOCKS];

   SingleNullEllipticOpBC* m_divergence_cleaning_bc;
   
   LayoutData< RefCountedPtr< IVSFAB<MBStencil> > > m_stencil_ue;
   LayoutData< IntVectSet > m_ghostCells_ue;

   Vector<MagBlockCoordSys *> m_coord_vec;

   std::string m_interpolation_scheme;
};



#include "NamespaceFooter.H"

#endif
