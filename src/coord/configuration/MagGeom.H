#ifndef MAGGEOM_H
#define MAGGEOM_H

#include "MultiBlockLevelGeom.H"
#include "MagCoordSys.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "BlockRegister.H"
#include "PotentialBC.H"

#include "BoundaryBoxLayout.H.multidim"

#include "NamespaceHeader.H"

class MagGeom : public MultiBlockLevelGeom
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

   MagGeom(ParmParse&                 pp,
           const MagCoordSys*         coordSysPtr,
           const DisjointBoxLayout&   grids,
           int                        ghosts);

  ///
  /** Destructor.
   */
  virtual ~MagGeom();

  /*@}*/

   void getMetricTerms( LevelData<FluxBox>& N,
                        LevelData<FluxBox>& tanGradN ) const;

   void getCellVolumes( LevelData<FArrayBox>& volume ) const;

   void getJ( LevelData<FArrayBox>& J ) const;

   void getPointwiseN( LevelData<FluxBox>& N ) const;

   void getPointwiseNJinverse( LevelData<FluxBox>& NJinverse ) const;

   void multJonValid( LevelData<FArrayBox>& dfn ) const;

   void multPointwiseJ( LevelData<FArrayBox>& dfn,
                        const BoundaryBoxLayout& bdry_layout ) const;

   // Computes physical face areas (no 2piR factors for axisymmetric)
   void getFaceAreas( LevelData<FluxBox>& areas ) const;

   void computeMetricTermProductAverage( LevelData<FluxBox>&       product,
                                         const LevelData<FluxBox>& data,
                                         const bool                fourthOrder ) const;

   void computeMappedGridDivergence( const LevelData<FluxBox>& F,
                                     LevelData<FArrayBox>&     divergence,
                                     const bool                fourthOrder ) const;

   void averageAtBlockBoundaries(LevelData<FluxBox>& data) const;

   void exchangeTransverseAtBlockBoundaries(LevelData<FluxBox>& data) const;

   void applyAxisymmetricCorrection( LevelData<FluxBox>& data ) const;

   void checkMultiblockMappingConsistency() const;

   double maxBlockBoundaryDifference(LevelData<FluxBox>& data) const;

   void getMagneticFlux(LevelData<FArrayBox>& psi) const;

   void getMagneticFlux(LevelData<FluxBox>& psi) const;

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

   void plotFieldAlignment(const double& time) const;

   void plotFieldDivergence( const LevelData<FluxBox>& field,
                             const double              time ) const;

   void plotCellData( const string&               file_name,  
                      const LevelData<FArrayBox>& data,
                      const double&               time ) const;

   void plotFaceData( const string&             file_name,  
                      const LevelData<FluxBox>& data,
                      const double&             time ) const;

   void writeGeometryData( const DisjointBoxLayout& grids, const double& time ) const;

   bool fullyPeriodic() const {return m_fully_periodic;}
    
   bool extrablockExchange() const {return m_extrablock_exchange;}

   bool modelGeom() const {return m_model_geometry;}

   const MagCoordSys* getCoordSys() const {return m_coord_sys;}

   void fillInternalGhosts( LevelData<FArrayBox>& data ) const;
   
   void exchangeExtraBlockGhosts( LevelData<FArrayBox>& data ) const;

   void fillTransversePhysicalGhosts( LevelData<FArrayBox>& data ) const;

   void fillTransversePhysicalGhosts( LevelData<FluxBox>& data ) const;

   void fillTransverseGhosts( LevelData<FluxBox>& data, const bool do_block_exchange ) const;

   void definePotentialBC( PotentialBC& bc ) const {m_coord_sys->definePotentialBC(bc);}

   const MagBlockCoordSys& getBlockCoordSys(const Box& box) const;

   const MagBlockCoordSys& getBlockCoordSys(const int block_number) const;

   const string getParmParsePrefix() const {return m_pp_prefix;}
    
   //Extracts a radial electric field at the top cut, and then interpolates
   //it on the grid (taking into account the the magnetic flux expansion)
   void interpolateErFromMagFS(LevelData<FluxBox>&   Er_face,
                               LevelData<FArrayBox>& Er_cell) const;
    
    
   //Returns the magnetic surface lalel (Z coordinate along the top cut)
   double getMagFS( const RealVect& xi ) const;


private:

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
    
    
   //Returns the magnetic flux function
   double magFlux( const RealVect& xi ) const;
    
   //Returns the magnetic flux function along the top cut
   double magFluxZ( const double Z ) const;
    
   LevelData<FArrayBox> m_magFS_mapping_cell;
   LevelData<FluxBox> m_magFS_mapping_face;
    
   mutable LevelData<FArrayBox> m_cell_volume;

   mutable LevelData<FArrayBox> m_J;

   mutable LevelData<FluxBox> m_N_face_centered;
   mutable LevelData<FluxBox> m_NJinverse_face_centered;

   const MagCoordSys* m_coord_sys;

   mutable LevelData<FluxBox> m_BFieldCorrection;

   MultiBlockLevelExchangeAverage* m_mblexPtr;
   BlockRegister* m_exchange_transverse_block_register;

   LevelData<FArrayBox> m_BField_cc;
   LevelData<FArrayBox> m_BFieldMag_cc;
   LevelData<FArrayBox> m_BFieldDir_cc;
   LevelData<FArrayBox> m_gradBFieldMag_cc;
   LevelData<FArrayBox> m_curlBFieldDir_cc;
   LevelData<FArrayBox> m_BFieldDirdotcurlBFieldDir_cc;

   LevelData<FluxBox> m_BField_fc;
   LevelData<FluxBox> m_BFieldMag_fc;
   LevelData<FluxBox> m_BFieldDir_fc;
   LevelData<FluxBox> m_gradBFieldMag_fc;
   LevelData<FluxBox> m_curlBFieldDir_fc;
   LevelData<FluxBox> m_BFieldDirdotcurlBFieldDir_fc;

   bool m_correct_field;

   bool m_fully_periodic;
   
   bool m_extrablock_exchange;
 
   bool m_model_geometry;

   const string m_pp_prefix;
};



#include "NamespaceFooter.H"

#endif
