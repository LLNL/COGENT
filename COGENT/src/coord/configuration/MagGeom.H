#ifndef MAGGEOM_H
#define MAGGEOM_H

#include "MultiBlockLevelGeom.H"
#include "MagCoordSys.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "BlockRegister.H"
#include "EllipticOpBC.H"
#include "EdgeDataBox.H"
#include "NodeFArrayBox.H"
#include "CoDimCopyManager.H"

#include "BoundaryBoxLayout.H.multidim"
#include "SpaceUtils.H"

#include "NamespaceHeader.H"

class MagGeom : public MultiBlockLevelGeom
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

   MagGeom(ParmParse&                         pp,
           const RefCountedPtr<MagCoordSys>&  coord_sys,
           const DisjointBoxLayout&           grids,
           int                                ghosts);

  ///
  /** Destructor.
   */
  virtual ~MagGeom();

  /*@}*/

   void getMetricTerms( LevelData<FluxBox>& N,
                        LevelData<FluxBox>& tanGradN ) const;

   void getMetricTerms( LevelData<FluxBox>& N) const;
   
   void setMetricTerms( const DisjointBoxLayout&, 
                        const IntVect&           ) const;
   
   void getCellVolumes( LevelData<FArrayBox>& volume ) const;
   
   void setCellVolumes( const DisjointBoxLayout&, 
                        const IntVect&           ) const;

   void getJ( LevelData<FArrayBox>& J ) const;
   
   void setJ( const DisjointBoxLayout&, 
              const IntVect&           ) const;

   void getJ( LevelData<FluxBox>& J ) const;
   
   void setJonFaces( const DisjointBoxLayout&, 
                     const IntVect&           ) const;
   
   void getJ( LevelData<EdgeDataBox>& J ) const;
   
   void setJonEdges( const DisjointBoxLayout&, 
                     const IntVect&           ) const;
   
   void getJ( LevelData<NodeFArrayBox>& J ) const;
   
   void setJonNodes( const DisjointBoxLayout&, 
                     const IntVect&           ) const;
   
   void getN( LevelData<FArrayBox>& N ) const;
   
   void setN( const DisjointBoxLayout&, 
              const IntVect&           ) const;
   
   void getdXdxi( LevelData<FArrayBox>& dXdxi ) const;
   
   void setdXdxiOnCells( const DisjointBoxLayout&, 
                         const IntVect&           ) const;
   
   void getdXdxi( LevelData<FluxBox>& dXdxi ) const;
   
   void setdXdxiOnFaces( const DisjointBoxLayout&, 
                         const IntVect&           ) const;
   
   void getdXdxi( LevelData<EdgeDataBox>& dXdxi ) const;
   
   void setdXdxiOnEdges( const DisjointBoxLayout&, 
                         const IntVect&           ) const;
   
   void getXphys( LevelData<FArrayBox>& Xphys ) const;
   
   void setXphys( const DisjointBoxLayout&, 
                  const IntVect&           ) const;

   void getPointwisePoloidalJ( LevelData<FArrayBox>& a_data ) const;

   void getPointwiseN( LevelData<FluxBox>& N ) const;

   void getPointwiseNJInverse( LevelData<FluxBox>& NJinverse ) const;
   
   void setPointwiseNJInverseOnFaces( const DisjointBoxLayout&, 
                                      const IntVect&           ) const;
   
   void getPointwiseNJInverse( LevelData<FArrayBox>& NJinverse ) const;
   
   void setPointwiseNJInverseOnCells( const DisjointBoxLayout&, 
                                      const IntVect&           ) const;

   void getPointwiseFaceAreas( LevelData<FluxBox>& face_areas ) const;

   void setPointwiseFaceAreas( const DisjointBoxLayout&,
                               const IntVect&           ) const;

   void getPointwiseFaceAreas( LevelData<FArrayBox>& face_areas_cc ) const;

   void setPointwiseFaceAreasOnCells( const DisjointBoxLayout&,
                                     const IntVect&           ) const;
   
   void multJonValid( LevelData<FArrayBox>& dfn ) const;
   
   void multJonFaces( LevelData<FluxBox>& dfn ) const;
   
   void multJonEdges( LevelData<EdgeDataBox>& dfn ) const;

   void divideJonValid( LevelData<FArrayBox>& dfn ) const;
   
   void divideJonFaces( LevelData<FluxBox>& dfn ) const;
   
   void divideJonEdges( LevelData<EdgeDataBox>& dfn ) const;

   void applyOnAxisCurlCorrection( LevelData<EdgeDataBox>& a_Curl,
                             const LevelData<FArrayBox>&   a_val) const;

   void multPointwiseJ( LevelData<FArrayBox>& dfn,
                        const BoundaryBoxLayout& bdry_layout ) const;

   // Computes physical face areas (no 2piR factors for axisymmetric)
   void getFaceAreas( LevelData<FluxBox>& areas ) const;

   void computeMetricTermProductAverage( LevelData<FluxBox>&       product,
                                         const LevelData<FluxBox>& data,
                                         const bool                fourthOrder ) const;
   
   void computedxidXProductAverage( LevelData<FluxBox>&  product,
                              const LevelData<FluxBox>&  data,
                              const bool                 fourthOrder ) const;

   void computeMappedGridDivergence( const LevelData<FluxBox>& F,
                                     LevelData<FArrayBox>&     divergence,
                                     const bool                fourthOrder ) const;
   
   void mappedGridDivergenceFromFluxNorms( const LevelData<FluxBox>&    F,
                                                 LevelData<FArrayBox>&  divergence) const;

   void mappedGridDivergenceFromFluxNorms( const FluxBox&    F,
                                                 FArrayBox&  divergence ) const;
   
   void physicalDivergenceFromEdgeTans( LevelData<NodeFArrayBox>&  divF, 
                                  const LevelData<EdgeDataBox>&    covar_F ) const;

   void physicalDivergenceAtEdges( LevelData<EdgeDataBox>&    a_divF_ce,
                             const LevelData<NodeFArrayBox>&  a_F_nc,
                             const LevelData<FArrayBox>&      a_F_cc ) const;

   void ParaGradAtNodesFromEdges( LevelData<NodeFArrayBox>&  a_gradF,
                            const LevelData<EdgeDataBox>&    a_F ) const;
   
   void PerpGradAtCellsFromEdges( LevelData<FArrayBox>&    a_gradF,
                            const LevelData<EdgeDataBox>&  a_F ) const;

   void mappedGridCurlFromEdgeTans( const LevelData<EdgeDataBox>&  E,
                                    const int                      face_dir,
                                          LevelData<FArrayBox>&    face_curl) const;
   
   void mappedGridCurlofVirtComp( const LevelData<FArrayBox>&    F,
                                        LevelData<EdgeDataBox>&  edge_curl) const;
   
   void mappedGridCurlofVirtComp( const LevelData<FArrayBox>&  F,
                                        LevelData<FArrayBox>&  cell_curl) const;
   
   void convertPhysToCovar( LevelData<FArrayBox>&  a_F,
                      const int                    a_inverse) const;
   
   void convertPhysToCovar( LevelData<FluxBox>&  a_F,
                      const int                  a_inverse) const;
   
   void convertPhysToContravar( LevelData<FArrayBox>&  a_F,
                          const int                    a_inverse) const;
   
   void convertPhysToContravar( LevelData<FluxBox>&  a_F,
                          const int                  a_inverse) const;
   
   void convertContravarToCovar( LevelData<EdgeDataBox>&  a_F,
                           const int                      a_inverse) const;
   
   void computeJaStrainTensorPhys( LevelData<FluxBox>&    a_JaW_phys, 
                             const LevelData<FArrayBox>&  a_velocity_phys ) const;
   
   void computeJaPhysStrainFlux( LevelData<FluxBox>&    a_JaW0_Flux_phys, 
                                 LevelData<FluxBox>&    a_JaW1_Flux_phys,
                           const LevelData<FArrayBox>&  a_velocity_phys ) const;
   
   // Project the vector from phys directions to unit normal
   // coord sys with 0-component normal to the a_dir direction 
   // (or inverse)
   //
   void projectOntoFluxSurfaceDir( FArrayBox&  a_data,
                             const DataIndex&  a_dit,
                             const int         a_dir,
                             const int         a_inverse ) const;

   void averageAtBlockBoundaries(LevelData<FluxBox>& data) const;

   void exchangeTransverseAtBlockBoundaries(LevelData<FluxBox>& data) const;

   void applyAxisymmetricCorrection( LevelData<FluxBox>& data ) const;

   void checkMultiblockMappingConsistency() const;

   double maxBlockBoundaryDifference(LevelData<FluxBox>& data) const;

   double maxMappingError() const;

   void getMagneticFlux(LevelData<FArrayBox>& psi) const;

   void getMagneticFlux(LevelData<FluxBox>& psi) const;
   
   const LevelData<FArrayBox>& getCellCenteredRealCoords() const {return m_cell_centered_real_coords;}

   const LevelData<FluxBox>& getFaceCenteredRealCoords() const {return m_face_centered_real_coords;}
   
   const LevelData<EdgeDataBox>& getEdgeCenteredRealCoords() const {return m_edge_centered_real_coords;}

   const LevelData<FArrayBox>& getNormalizedMagneticFluxCell() const {return m_normalized_magnetic_flux_cell;}

   const LevelData<FluxBox>& getNormalizedMagneticFluxFace() const {return m_normalized_magnetic_flux_face;}

   const LevelData<FArrayBox>& getCCBField() const {return m_BField_cc;}

   const LevelData<FArrayBox>& getCCBFieldMag() const {return m_BFieldMag_cc;}

   const LevelData<FArrayBox>& getCCBFieldDir() const {return m_BFieldDir_cc;}

   const LevelData<FArrayBox>& getCCGradBFieldMag() const {return m_gradBFieldMag_cc;}

   const LevelData<FArrayBox>& getCCCurlBFieldDir() const {return m_curlBFieldDir_cc;}

   const LevelData<FArrayBox>& getCCBFieldDirdotCurlBFieldDir() const {return m_BFieldDirdotcurlBFieldDir_cc;}

   const LevelData<FluxBox>& getFCBField() const {return m_BField_fc;}

   const LevelData<FluxBox>& getFCBFieldMag() const {return m_BFieldMag_fc;}

   const LevelData<FluxBox>& getFCBFieldDir() const {return m_BFieldDir_fc;}

   const LevelData<FluxBox>& getFCGradBFieldMag() const {return m_gradBFieldMag_fc;}

   const LevelData<FluxBox>& getFCCurlBFieldDir() const {return m_curlBFieldDir_fc;}

   const LevelData<FluxBox>& getFCBFieldDirdotCurlBFieldDir() const {return m_BFieldDirdotcurlBFieldDir_fc;}
   
   const LevelData<FluxBox>& getEllipticOpPerpCoeff() const {return m_perp_coeff;}
   
   const LevelData<FluxBox>& getEllipticOpPerpCoeffMapped() const {return m_perp_coeff_mapped;}

   const LevelData<FluxBox>& getEllipticOpParCoeff() const {return m_par_coeff;}
   
   const LevelData<FluxBox>& getEllipticOpParCoeffMapped() const {return m_par_coeff_mapped;}
   
   void getEllipticOpRadCoeff(LevelData<FluxBox>& rad_coeff) const;
   
   void getEllipticOpRadCoeffMapped(LevelData<FluxBox>& rad_coeff_mapped) const;

   void plotMagneticFieldData(const double& time) const;
   
   void plotFieldAlignment(const double& time) const;

   void plotFieldDivergence( const LevelData<FluxBox>& field,
                             const double             time ) const;

   void plotCellData( const string&               file_name,  
                      const LevelData<FArrayBox>& data,
                      const double&               time ) const;

   void plotFaceData( const string&             file_name,  
                      const LevelData<FluxBox>& data,
                      const double&             time ) const;

   void writeGeometryData( const DisjointBoxLayout& grids, const double& time ) const;

   bool fullyPeriodic() const {return m_fully_periodic;}
    
   bool extrablockExchange() const {return m_extrablock_exchange;}

   const RefCountedPtr<MagCoordSys>& getCoordSys() const {return m_coord_sys;}

   void fillInternalGhosts( LevelData<FArrayBox>& data ) const;
   
   void exchangeExtraBlockGhosts( LevelData<FArrayBox>& data ) const;

   void exchangeExtraBlockGhosts( LevelData<FArrayBox>& data,
                                  const int codim,
                                  const IntVect& boundary_dirs,
                                  const IntVect& exchange_dirs ) const;

   void getExtrablockExchangeLayout( const int                 dir,
                                     const int                 codim,
                                     const DisjointBoxLayout&  grids,
                                     const IntVect&            ghost_vect,
                                     const IntVect&            exchange_dirs,
                                     BoxLayout&                ghosts_layout,
                                     LayoutData<DataIndex>&    index_map,
                                     LayoutData<IntVect>&      shift_map ) const;
  
   void fillCorners( LevelData<FArrayBox>&  a_data,
                     const IntVect&         a_nghosts,
                     const int              a_max_codim ) const;

   void fillTransversePhysicalGhosts( LevelData<FArrayBox>& data ) const;

   void fillTransversePhysicalGhosts( LevelData<FluxBox>& data ) const;

   void fillTransverseGhosts( LevelData<FluxBox>& data, const bool do_block_exchange ) const;

   void defineEllipticOpBC( EllipticOpBC& bc ) const {m_coord_sys->defineEllipticOpBC(bc);}

   const MagBlockCoordSys& getBlockCoordSys(const Box& box) const;

   const MagBlockCoordSys& getBlockCoordSys(const int block_number) const;

   const string getParmParsePrefix() const {return m_pp_prefix;}
    
   void getIntegralsForVelocity( const LevelData<FArrayBox>& Efield_cell,
                                 const LevelData<FArrayBox>& phi_node,
                                 const bool                  fourth_order,
                                 LevelData<FArrayBox>&       nodal_integrals,
                                 LevelData<FluxBox>&         face_integrals,
                                 LevelData<FArrayBox>&       volume_integrals ) const;

   void computeBxEIntegrals( const LevelData<FArrayBox>& phi_node,
                             bool                        fourth_order,
                             LevelData<FluxBox>&         Exb_integral ) const;

   void computeEXBDrift( const LevelData<FluxBox>& E_field,
                         LevelData<FluxBox>&       ExB_drift ) const;
   
   //Extracts a radial electric field at the top cut, and then interpolates
   //it on the grid (taking into account the the magnetic flux expansion)
   void interpolateErFromMagFS(LevelData<FluxBox>&   Er_face,
                               LevelData<FArrayBox>& Er_cell) const;
    
    
   //Returns the magnetic surface lalel (Z coordinate along the top cut)
   double getMagFS( const RealVect& xi ) const;

   // Convert a 3-dimensional gradient in mapped coordinates to the physical frame
   void unmapGradient( const LevelData<FArrayBox>& mapped_gradient,
                       LevelData<FArrayBox>&       gradient ) const;

   // Convert a 3-dimensional gradient in mapped coordinates to the physical frame
   void unmapGradient( const LevelData<FluxBox>& mapped_gradient,
                       LevelData<FluxBox>&       gradient ) const;

   // Convert a 3D gradient in mapped coordinates to the physical the frame
   void unmap3DGradient( const LevelData<FArrayBox>& mapped_gradient,
                               LevelData<FArrayBox>&       gradient ) const;

   // Convert a 3D gradient in mapped coordinates to the physical the frame
   void unmap3DGradient( const LevelData<FluxBox>& mapped_gradient,
                               LevelData<FluxBox>&       gradient ) const;

   // Convert a poloidal gradient in mapped coordinates to the physical the frame
   void unmapPoloidalGradient( const LevelData<FArrayBox>& mapped_gradient,
                               LevelData<FArrayBox>&       gradient ) const;

   // Convert a poloidal gradient in mapped coordinates to the physical the frame
   void unmapPoloidalGradient( const LevelData<FluxBox>& mapped_gradient,
                               LevelData<FluxBox>&       gradient ) const;

   void injectPoloidalVector( const LevelData<FArrayBox>& poloidal_vector,
                              LevelData<FArrayBox>&       vector ) const;

   void injectPoloidalVector( const LevelData<FluxBox>& poloidal_vector,
                              LevelData<FluxBox>&       vector ) const;

   void projectPoloidalVector( const LevelData<FArrayBox>& vector,
                               LevelData<FArrayBox>&       poloidal_vector ) const;

   void projectPoloidalVector( const LevelData<FluxBox>& vector,
                               LevelData<FluxBox>&       poloidal_vector ) const;

   //Project vector onto the parallel (to the magneitc field) direction
   void projectOntoParallel( LevelData<FArrayBox>& vector ) const;
   
   void projectOntoParallel( LevelData<FluxBox>& vector ) const;
   
   //Compute vector projection onto the direction of phi x bpol
   void computeRadialProjection( LevelData<FArrayBox>& radComp,
                                 const LevelData<FArrayBox>& vector) const;

   void computeRadialProjection( LevelData<FluxBox>& radComp,
                                 const LevelData<FluxBox>& vector) const;

   //Compute vector projection onto the direction of bpol
   void computePoloidalProjection( LevelData<FArrayBox>& polComp,
                                   const LevelData<FArrayBox>& vector) const;
   
   void computePoloidalProjection( LevelData<FluxBox>& polComp,
                                   const LevelData<FluxBox>& vector) const;
   
   void computeMappedPoloidalGradientWithGhosts( const LevelData<FArrayBox>& phi,
                                                 LevelData<FArrayBox>&       field,
                                                 const int                   order ) const;

   void computeMappedPoloidalGradientWithGhosts( const LevelData<FArrayBox>& phi,
                                                LevelData<FluxBox>&          field,
                                                const int                    order ) const;
   
   void computeMappedGradient( const LevelData<FArrayBox>&  phi,
                                     LevelData<FluxBox>&    field,
                               const int                    order ) const;

   bool shearedMBGeom() const {return m_sheared_mb_geom;}
  
   bool fieldAlignedMapping() const {return m_coord_sys->getCoordSys(0)->isFieldAlignedMapping();}
  
   bool secondOrder() const {return m_second_order;}

   bool mixedBoundaries() const;
   
#if CFG_DIM == 3

   int shearedInterpOrder() const {return m_sheared_interp_order;}
   int shearedGhosts() const {return m_sheared_ghosts;}
   int getMultiblockDir() const {return m_mb_dir;}
   const LevelData<FArrayBox>& getShearedRemappedIndex() const {return m_sheared_remapped_index;}
   const LevelData<FArrayBox>& getShearedInterpStencil() const {return m_sheared_interp_stencil;}
   const LevelData<FArrayBox>& getShearedInterpStencilOffsets() const {return m_sheared_interp_stencil_offsets;}
   void interpolateFromShearedGhosts(LevelData<FArrayBox>& a_data) const;
   
#endif

   Vector<Box> getBoundaryBoxes( const int block_number, const int dir, const Side::LoHiSide& side ) const
      {
         return m_coord_sys->getBoundaryBoxes(block_number, dir, side);
      }

   bool providesFlux() const {return m_coord_sys->providesFlux();}
   
private:

   void setRealCoords();

   void setNormalizedMagneticFlux();

   void computeFieldData( LevelData<FArrayBox>& BField,
                          LevelData<FArrayBox>& BFieldMag,
                          LevelData<FArrayBox>& BFieldDir,
                          LevelData<FArrayBox>& gradBFieldMag,
                          LevelData<FArrayBox>& curlBFieldDir,
                          LevelData<FArrayBox>& BFieldDirdotcurlBFieldDir ) const;

   void computeFieldData( LevelData<FluxBox>& BField,
                          LevelData<FluxBox>& BFieldMag,
                          LevelData<FluxBox>& BFieldDir,
                          LevelData<FluxBox>& gradBFieldMag,
                          LevelData<FluxBox>& curlBFieldDir,
                          LevelData<FluxBox>& BFieldDirdotcurlBFieldDir ) const;

   void cellCenter( const LevelData<FluxBox>& fc_data,
                    LevelData<FArrayBox>&     cc_data ) const;

   // Below is the set of routines that address the dealignment between the grid and the magnetic flux surfaces.
   // The magnietc flux surfaces are labeled by the Z coordinate along the top cut (boundaries between LCORE/RCORE and LCSOL/RCSOL).
   // Methods are presently of second order.
    
   // The computeMagFluxMappingCell/Face functions retrun mapping from the arbitrary grid point (face-centered of cell-centered)
   // comp0 = Z, magneitc surface label (i.e., Z at the top cut )
   // comp1 = irad, cell-centered index (center of the poloidal face) at the top cut, where Z[irad-1]<Z<Z[irad]
   // comp2 = interp_coeff, the interpolation coefficient defined as
   // for Z[0]<Z<Z[nrad-1], F[Z] = F[irad]-(F[irad]-F[irad-1])*interp_coeff, interp_coeff = (Z[irad]-Z)/(Z[irad]-Z[irad-1])
   // for Z>Z[nrad-1], comp1 = nrad, F[Z]=F[nrad-1]+(F[nrad-1]-F[nrad-2])*interp_coeff, interp_coeff = (Z-Z[nrad-1])/(Z[nrad-1]-Z[nrad-2])
   // for Z<Z[0], comp1 = 0, F[Z]=F[0]-(F[1]-F[0])*interp_coeff, interp_coeff = (Z[0]-Z)/(Z[1]-Z[0])
   // nrad - number of the radial cells
   void computeMagFluxMappingCell( LevelData<FArrayBox>& magFS_mapping_cell) const;
   void computeMagFluxMappingFace(LevelData<FluxBox>& magFS_mapping_face) const;
    
   //Creates a flux-function by extracting data along the top cut, and then deposits it onto the grid (second order)
   void interpolateFromMagFS( const LevelData<FluxBox>& data,
                              LevelData<FluxBox>& interp) const;
    
   void interpolateFromMagFS( const LevelData<FluxBox>& data_face,
                              const LevelData<FArrayBox>& data_cell,
                              LevelData<FArrayBox>& interp) const;
    
    
   /// Returns the magnetic flux function along the top cut
   double magFluxZ( const double Z ) const;

   /// Returns the magnetic flux for the model flux function
   double modelMagFlux( const RealVect& xi) const;
   
   
   /// Computes perpendicular [(I - b * bT)] and parallel [b * bT] matrices that are used in eliptic operators
   void computeEllipticOpCoefficients(LevelData<FluxBox>& a_perp_coeff,
                                      LevelData<FluxBox>& a_par_coeff ) const;

   /// Computes mapped perpendicular [NT*(I - b * bT)*NJinv] and parallel [NT*(b * bT)*NJinv] matrices that are used in eliptic operators
   void computeEllipticOpCoefficientsMapped(LevelData<FluxBox>& a_perp_coeff,
                                            LevelData<FluxBox>& a_par_coeff ) const;
  
   /// Computes mapped face-averaged coefficients by subgrid integrations
   void computeSubGridEllipticOpCoefficientsMapped(FluxBox& perp_coeff,
                                                   FluxBox& par_coeff,
                                                   const MagBlockCoordSys& coord_sys) const;
   
   /// Plots geometry and B-field data
   void plotInitializationData(const ParmParse& a_pp,
                               const DisjointBoxLayout& grids,
                               const double& time) const;

   
   
   LevelData<FArrayBox> m_magFS_mapping_cell;
   LevelData<FluxBox> m_magFS_mapping_face;
   
   LevelData<FluxBox> m_perp_coeff;
   LevelData<FluxBox> m_par_coeff;
   LevelData<FluxBox> m_perp_coeff_mapped;
   LevelData<FluxBox> m_par_coeff_mapped;
   
   mutable LevelData<FluxBox> m_rad_coeff;
   mutable LevelData<FluxBox> m_rad_coeff_mapped;
    
   mutable LevelData<FArrayBox> m_cell_volume;
   mutable LevelData<FluxBox> m_face_areas;
   mutable LevelData<FArrayBox> m_face_areas_cc;

   mutable LevelData<FArrayBox> m_J;
   mutable LevelData<FluxBox> m_JonFaces;
   mutable LevelData<EdgeDataBox> m_JonEdges;
   mutable LevelData<NodeFArrayBox> m_JonNodes;
   mutable LevelData<FArrayBox> m_N;
   mutable LevelData<FArrayBox> m_dXdxi;
   mutable LevelData<FluxBox> m_dXdxi_onFaces;
   mutable LevelData<EdgeDataBox> m_dXdxi_onEdges;
   mutable LevelData<FArrayBox> m_Xphys;

   mutable LevelData<FluxBox> m_metrics;
   mutable LevelData<FluxBox> m_N_face_centered;
   mutable LevelData<FluxBox> m_NJinverse_face_centered;
   mutable LevelData<FArrayBox> m_NJinverse_cell_centered;
   mutable LevelData<FArrayBox> m_poloidal_J;

   const RefCountedPtr<MagCoordSys> m_coord_sys;

   mutable LevelData<FluxBox> m_BFieldCorrection;

   MultiBlockLevelExchangeAverage* m_mblexPtr;
   BlockRegister* m_exchange_transverse_block_register;
   
   LevelData<FArrayBox> m_BField_cc;
   LevelData<FArrayBox> m_BFieldMag_cc;
   LevelData<FArrayBox> m_BFieldDir_cc;
   LevelData<FArrayBox> m_gradBFieldMag_cc;
   LevelData<FArrayBox> m_curlBFieldDir_cc;
   LevelData<FArrayBox> m_BFieldDirdotcurlBFieldDir_cc;

   LevelData<FArrayBox> m_cell_centered_real_coords;
   LevelData<FluxBox> m_face_centered_real_coords;
   LevelData<EdgeDataBox> m_edge_centered_real_coords;

   LevelData<FArrayBox> m_normalized_magnetic_flux_cell;
   LevelData<FluxBox> m_normalized_magnetic_flux_face;
   
   LevelData<FluxBox> m_BField_fc;
   LevelData<FluxBox> m_BFieldMag_fc;
   LevelData<FluxBox> m_BFieldDir_fc;
   LevelData<FluxBox> m_gradBFieldMag_fc;
   LevelData<FluxBox> m_curlBFieldDir_fc;
   LevelData<FluxBox> m_BFieldDirdotcurlBFieldDir_fc;
   
   mutable LevelData<FArrayBox> m_nodal_integrals;
   mutable LevelData<FluxBox> m_face_Binverse_integrals;
   mutable LevelData<FluxBox> m_face_b_integrals;
   mutable LevelData<FluxBox> m_face_bXgradB_integrals;
   mutable LevelData<FArrayBox> m_volume_B_integrals;
   mutable LevelData<FArrayBox> m_volume_BdotGradB_integrals;
   mutable bool m_velocity_data_cached;

   bool m_correct_field;
   bool m_fully_periodic;
   bool m_extrablock_exchange;
   bool m_sheared_mb_geom;
   bool m_dealignment_corrections;
   bool m_model_flux_function; 
   bool m_second_order;
   
   //Copier m_exchangeCopier2nd;

   const string m_pp_prefix;
  
#if CFG_DIM == 3
   
   void initializeShearedMBGeom(const DisjointBoxLayout& grids);
   void getShearedGhostBoxLayout();
   void getMagShearInterpCoeff(LevelData<FArrayBox>& remapped_iv,
                               LevelData<FArrayBox>& interp_stencil,
                               LevelData<FArrayBox>& interp_stencil_offsets,
                               const DisjointBoxLayout& grids,
                               const Side::LoHiSide& side) const;
   
   int m_sheared_interp_order;
   int m_sheared_ghosts;
   int m_mb_dir;
   
   BoxLayout m_shearedGhostBLLoEnd;
   BoxLayout m_shearedGhostBLHiEnd;
   
   LevelData<FArrayBox> m_sheared_remapped_index;
   LevelData<FArrayBox> m_sheared_interp_stencil;
   LevelData<FArrayBox> m_sheared_interp_stencil_offsets;
   
#endif
   
   // Cached data supporting extrablock exchanges
   struct EBE_Data
   {
      BoxLayout             ghostLayout;
      LayoutData<DataIndex> index_map;
      LayoutData<IntVect>   shift_map;
      IntVect               ghost_vect;
      IntVect               exchange_dirs;
      Copier                copier;
   };
   mutable list<EBE_Data*> m_ebe_data_cache[CH_SPACEDIM][CH_SPACEDIM];
   mutable list<CoDimCopyManager<FArrayBox>*> m_ebe_copy_manager_cache;
};



#include "NamespaceFooter.H"

#endif
