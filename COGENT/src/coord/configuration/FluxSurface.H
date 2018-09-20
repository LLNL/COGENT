#ifndef _FLUXSURFACE_H_
#define _FLUXSURFACE_H_

enum BlockType {CORE, PF, SOL};

#include "MagGeom.H"

#ifdef CH_MULTIDIM

#include "NamespaceHeader.H"


/// Flux surface operation class
/**
 * A utility class for performing operations on the flux surfaces
 * of a magnetic geometry.
 */
class FluxSurface
{
public:

   /// Constructor
   /** Constructor
    * @param[in] coords  The magnetic geometry
    */
   FluxSurface( const MagGeom& geom,
                const bool     shell_average = true );

   /// Destructor
   /** Destructor
    */
   ~FluxSurface() {;}

   /// Compute the flux surface areas
   /** Compute the flux surface areas
    */
   void computeAreas();

   /// Integrate the src array over flux surfaces
   /** Integrate the src array over flux surfaces
    *
    * @param[in]  src  Data to be integrated over flux surfaces
    * @param[out] dst  Flux surface integrals
    */
   void integrate( const LevelData<FArrayBox>& src,
                   LevelData<FArrayBox>&       dst ) const;

   /// Average the src array over flux surfaces
   /** Average the src array over flux surfaces
    *
    * @param[in]  src  Data to be averaged over flux surfaces
    * @param[out] dst  Flux surface averages
    */
   void average( const LevelData<FArrayBox>& src,
                 LevelData<FArrayBox>&       dst ) const;

   /// Sum the src array over flux surfaces
   /** Sum the src array over flux surfaces
    *
    * @param[in]  src  Data to be summed over flux surfaces
    * @param[out] dst  Flux surface sums
    */
   void sum( const LevelData<FArrayBox>& src,
             LevelData<FArrayBox>&       dst ) const;

   /// Spread the src array over flux surfaces
   /** Spread the src array over flux surfaces
    *
    * @param[in]  src  Data to be spread over flux surfaces
    * @param[out] dst  Spread data
    */
   void spread( const LevelData<FArrayBox>& src,
                LevelData<FArrayBox>&       dst ) const;

   /// Average and spread the src array over flux surfaces (includes single-null topology)
   /** Average and spread the src array over flux surfaces (includes single-null topology)
    *
    * @param[in]  src  Data to be spread over flux surfaces 
    * @param[out] dst  Spread data
    */
   void averageAndSpread( const LevelData<FArrayBox>& src,
                          LevelData<FArrayBox>&       dst ) const;


   /// Add flux surface data
   /** Add flux surface data
    *
    * @param[in]     src  Flux surface data
    * @param[in,out] dst  Data to be added to
    */
   void add( const LevelData<FArrayBox>& flux_surface_data,
             LevelData<FArrayBox>&       data              ) const;

   /// Subtract flux surface data
   /** Subtract flux surface data
    *
    * @param[in]     src  Flux surface data
    * @param[in,out] dst  Data to be subtracted from
    */
   void subtract( const LevelData<FArrayBox>& flux_surface_data,
                  LevelData<FArrayBox>&       data              ) const;

   /// Multiply by flux surface data
   /** Multiply by flux surface data
    *
    * @param[in]     src  Flux surface data
    * @param[in,out] dst  Data to be multiplied
    */
   void multiply( const LevelData<FArrayBox>& flux_surface_data,
                  LevelData<FArrayBox>&       data              ) const;

   /// Divide by flux surface data
   /** Divide by flux surface data
    *
    * @param[in]     src  Flux surface data
    * @param[in,out] dst  Data to be divided
    */
   void divide( const LevelData<FArrayBox>& flux_surface_data,
                LevelData<FArrayBox>&       data              ) const;

   /// Set values in the corresponding blocks to zero
   /** Set values in the corresponding blocks to zero
    *
    * @param[in]     block_type  Block type
    * @param[in]     src         Data to be cropped
    * @param[in,out] dst         Cropped data
    */
   void zeroBlockData( const int                   block_type,
                       const LevelData<FArrayBox>& src,
                       LevelData<FArrayBox>&       dst        ) const;

   /// Multiplies by gradPsi
    void multGradPsi(LevelData<FArrayBox>& a_src ) const;

    /// Const reference to the flux surface grids
   /** Const reference to the flux surface grids
    *
    * @param[out]  grids  Flux surface grids
    */
   const DisjointBoxLayout& grids() const {return m_grids;}

   /// Const reference to the flux surface areas
   /** Const reference to the flux surface areas
    *
    * @param[out]  areas  Flux surface areas
    */
   const LevelData<FArrayBox>& areas() const {return m_areas;}

protected:

   const MagGeom* m_magnetic_geometry;

   bool m_single_null;
   bool m_shell_average;

   DisjointBoxLayout m_grids;

   LevelData<FArrayBox> m_areas;
   LevelData<FArrayBox> m_areas_core;
   LevelData<FArrayBox> m_areas_pf;

   LevelData<FArrayBox> m_volume;
};

#include "NamespaceFooter.H"

#endif /* CH_MULTIDIM */

#endif   /* _FLUXSURFACE_H_ */
