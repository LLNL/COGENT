#ifndef SINGLENULLCOORDSYS_H
#define SINGLENULLCOORDSYS_H

#include "MagCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "SingleNullEllipticOpBC.H"

#include "NamespaceHeader.H"

#define BLOCK_SEPARATION 2

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

   int findBlockNearXpt(const RealVect& X ) const;

   virtual void defineBoundaries();

   virtual IntVect getDecomposition(int block) const;

   virtual const EllipticOpBC& getDivergenceCleaningBC() const {return m_divergence_cleaning_bc;}

   virtual void defineEllipticOpBC( EllipticOpBC& bc ) const;

   int blockType(int a_block) const {return ((SingleNullBlockCoordSys*)m_coordSysVect[a_block])->blockType();}

   int numCellsCoreRadial() const {return m_numcells_core_radial;}

   int numCellsMcorePoloidal() const {return m_numcells_mcore_poloidal;}

   int numCellsLcorePoloidal() const {return m_numcells_lcore_poloidal;}

   int numCellsRcorePoloidal() const {return m_numcells_rcore_poloidal;}
   
   bool isModelGeom() const {return m_model_geometry;}
   
     //ABSOLETE: Extrapolate E-field into the open-field-line region
   void extrapolateEfield( LevelData<FluxBox>&   Er_average_face,
                           LevelData<FArrayBox>& Er_average_cell ) const;
   RealVect m_Xpoint;

   void defineStencilsUe( const DisjointBoxLayout&  grids,
                         const int                 nghost );

   
   const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& getStencilsUe() const {return m_stencil_ue;}
   
   const LayoutData< IntVectSet >& getGhostCellsUe() const {return m_ghostCells_ue;}
   
   // Object ParmParse name.
   static const std::string pp_name;
   
   virtual std::string type() const { return "SingleNull"; }

#if CFG_DIM ==3
   /// toroidal block mapping conversion function
   /** Given a_xiSrc in mapped space of block a_nSrc,
    and a_side, returns coordinate in mapped space of the
    destination block (determined from a_side argument)
    and the  second-order interpolation stencil.
    */
   
   virtual void toroidalBlockRemapping(IntVect& a_ivDst,
                                       Vector<Real>& a_interpStecil,
                                       const RealVect& a_xiSrc,
                                       const int a_nSrc,
                                       const Side::LoHiSide& a_side) const;
   
#endif

protected:

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


   void readGridParams( ParmParse& a_pp );

   void printInvalidPointDiagnostics( const int                      nSrc,
                                      const RealVect&                xiSrc,
                                      const RealVect&                X,
                                      const RealVect&                xi_valid,
                                      const int                      n_valid,
                                      const SingleNullBlockCoordSys* src_cs,
                                      const SingleNullBlockCoordSys* valid_cs ) const;
    
   

   void checkGridConsistency() const;

   IntVect lo_mapped_index(int block) const;

   IntVect hi_mapped_index(int block) const;

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

   double m_xpoint_radius;
        
   double m_ray_angle[NUM_SINGLE_NULL_BLOCKS+1];
   RealVect m_ray[NUM_SINGLE_NULL_BLOCKS];
   int m_sector_block[NUM_SINGLE_NULL_BLOCKS];

   SingleNullEllipticOpBC m_divergence_cleaning_bc;
   
   LayoutData< RefCountedPtr< IVSFAB<MBStencil> > > m_stencil_ue;
   LayoutData< IntVectSet > m_ghostCells_ue;

   Vector<MagBlockCoordSys *> m_coord_vec;
};



#include "NamespaceFooter.H"

#endif
