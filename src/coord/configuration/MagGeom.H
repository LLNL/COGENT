#ifndef MAGGEOM_H
#define MAGGEOM_H

#include "MultiBlockLevelGeom.H"
#include "MagCoordSys.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "BlockRegister.H"

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

   void getNJInverse( LevelData<FluxBox>& NJInverse ) const;

   void getPointwiseN( LevelData<FluxBox>& a_N ) const;

   void getPointwiseNJinverse( LevelData<FluxBox>& a_NJinverse ) const;

   void multPointwiseJ( LevelData<FArrayBox>& dfn,
                        const BoundaryBoxLayout& bdry_layout ) const;

   // Computes physical face areas (no 2piR factors for axisymmetric)
   void getFaceAreas( LevelData<FluxBox>& a_areas ) const;

   void computeMetricTermProductAverage( LevelData<FluxBox>&       a_product,
                                         const LevelData<FluxBox>& a_data,
                                         const bool                a_fourthOrder ) const;

   void computeMappedGridDivergence( const LevelData<FluxBox>& F,
                                     LevelData<FArrayBox>&     divergence,
                                     const bool                fourthOrder ) const;

   void averageAtBlockBoundaries(LevelData<FluxBox>& a_data) const;

   void exchangeTransverseAtBlockBoundaries(LevelData<FluxBox>& a_data) const;

   void applyAxisymmetricCorrection( LevelData<FluxBox>& data ) const;

   void checkMultiblockMappingConsistency() const;

   double maxBlockBoundaryDifference(LevelData<FluxBox>& a_data) const;

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

   void plotFieldAlignment() const;

   void plotCellData( const string&               file_name,  
                      const LevelData<FArrayBox>& data ) const;

   void plotFaceData( const string&             file_name,  
                      const LevelData<FluxBox>& data ) const;

   void writeGeometryData( const DisjointBoxLayout& grids ) const;

   bool fullyPeriodic() const {return m_fully_periodic;}

   const MagCoordSys* getCoordSys() const {return m_coord_sys;}

   void fillInternalGhosts( LevelData<FArrayBox>& a_data ) const;

   void fillTransverseGhosts( LevelData<FluxBox>& a_data, const bool do_block_exchange ) const;

   void definePotentialBC( PotentialBC& bc ) const {m_coord_sys->definePotentialBC(bc);}

   const MagBlockCoordSys& getBlockCoordSys(const Box& box) const;

   const MagBlockCoordSys& getBlockCoordSys(const int block_number) const;

   const string getParmParsePrefix() const {return m_pp_prefix;}

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

   const string m_pp_prefix;
};



#include "NamespaceFooter.H"

#endif
