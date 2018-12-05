#ifndef PHASEGEOM_H
#define PHASEGEOM_H

#include "Dimensions.H"
#include "BoundaryBoxLayout.H.multidim"
#include "MultiBlockLevelGeom.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "PhaseCoordSys.H"
#include "PhaseBlockCoordSys.H"
#include "PhaseGrid.H"
#include "BlockRegister.H"


#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "MagGeom.H"
#include "MultiBlockLevelGeom.H"
#include "MultiBlockLevelExchangeAverage.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "GyroaverageOperator.H"

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;
namespace VEL = VEL_NAMESPACE;
namespace CP1 = CP1_NAMESPACE;

class PhaseGeom : public MultiBlockLevelGeom
{
public:

   /**
      Constructors, destructor and defines
   */
   PhaseGeom( ParmParse&                              parm_parse,
              const RefCountedPtr<PhaseCoordSys>&     coord_sys,
              const RefCountedPtr<PhaseGrid>&         grids,
              const RefCountedPtr<CFG::MagGeom>&      mag_geom,
              const RefCountedPtr<VEL::VelCoordSys>&  vel_coords,
              int                                     ghosts,
              double                                  larmor_number );

   // Species-specific constructor
   PhaseGeom( const PhaseGeom&                        phase_geom,
              const RefCountedPtr<PhaseCoordSys>&     coord_sys,
              const RefCountedPtr<VEL::VelCoordSys>&  vel_coords,
              double                                  mass,
              double                                  charge );
   /**
      Destructor.
    */
   virtual ~PhaseGeom();

   void define();

   void defineSpeciesState( double mass, double charge_state );

   /// Computes the phase space velocities R_dot and v_dot
   /**
    * Computes the phase space velocities
    *
    *    R_dot = v_parallel * (BStar/BStarPar) + (larmor_number / charge_state) * (b/BStarPar) X G
    *
    * and
    *
    *    v_dot = - (1/mass) * (BStar/BStarPar) dot G
    *
    * where
    *
    *    BStar = B + prefactor * v_parallel * curl_b
    *
    *    BStarPar = BStar dot b
    *
    *    G = - charge_state * EField + 0.5 * mu * gradB
    *
    * @param[in]  phi   Potential
    * @param[out] velocity Phase space velocity
    */
   void updateVelocities( const LevelData<FluxBox>& Efield,
                          LevelData<FluxBox>&       velocity,
                          const int                 option,
                          const bool                apply_axisymmetric_correction = false) const;

   void computeGKVelocities( const LevelData<FluxBox>& Efield,
                             LevelData<FluxBox>&       velocity_grown,
                             const int                 option ) const;

   void computeTestVelocities( LevelData<FluxBox>& velocity ) const;

   void updateMappedVelocities( const LevelData<FluxBox>& Efield,
                                LevelData<FluxBox>&       velocity ) const;

   //Updates normal velocities, velocity options correspond to
   //0 --  full GK velocity
   //1 --  drifts only (includes both gradB/curlb and ExB)
   //2 --  no 0-order parallel terms (that cancel each other for an FS-averaged Maxwellian)
   void updateVelocityNormals( const CFG::LevelData<CFG::FArrayBox>& Efield_cell,
                               const CFG::LevelData<CFG::FArrayBox>& phi_node,
                               const bool                fourth_order_Efield,
                               LevelData<FluxBox>&       velocity,
                               const int                 velocity_option = 0 ) const;

   void computeXStarOmega( const CFG::LevelData<CFG::FArrayBox>& Efield_cell,
                           const CFG::LevelData<CFG::FArrayBox>& phi_node,
                           const bool                fourth_order_Efield,
                           LevelData<FluxBox>&       XStarOmega,
                           LevelData<FluxBox>&       UhatNormal,
                           const int                 velocity_option = 0 ) const;

   /// Computes the metric 3-vector BStar
   /**
    * Computes the metric 3-vector
    *
    *    BStar = B + (larmor_number*m_mass/charge_state) * v_parallel * curl_b
    *
    * and its parallel component
    *
    *    BStarParallel = b dot BStar
    *
    */
   void computeBStar(bool                        no_drifts,
                     double                      larmor_number,
                     double                      mass,
                     double                      charge_state,
                     const LevelData<FArrayBox>& B,
                     const LevelData<FArrayBox>& BMag,
                     const LevelData<FArrayBox>& curlb,
                     const LevelData<FArrayBox>& bdotcurlb,
                     LevelData<FArrayBox>&       BStar,
                     LevelData<FArrayBox>&       BStarParallel) const;

   /// Computes product of the phase space metric factors and the input flux,
   /// to second- or fourth-order
   /**
    * Computes the product of the phase space metric factors and the input
    * flux to second- or fourth-order.
    *
    * @param[out] Product product of N^T*F
    * @param[in]  F             flux
    * @param[in]  fourthOrder   product is fourth-order if true, second-order otherwise
    */
   void computeMetricTermProductAverage( LevelData<FluxBox>&       Product,
                                         const LevelData<FluxBox>& F,
                                         const bool                fourthOrder ) const;

   void computeTangentialGradSpecial( LevelData<FluxBox>&       gradPhi,
                                      const LevelData<FluxBox>& phiFace ) const;

   void mappedGridDivergence(LevelData<FArrayBox>&     divF,
                             const LevelData<FluxBox>& F,
                             const bool                omit_NT) const;

   void zeroBoundaryFlux( const int dir, LevelData<FluxBox>& a_flux ) const;

   /// Compute the phase space flux divergence(F) assuming that the input
   /// flux F already contains metric and Jacobian values. Difference between
   /// this function and mappgedGridDivergnceWithoutNT is 1) this function
   /// accepts fluxes with any number of components (as long as it is the same
   /// as the number of components in divF) and 2) computes div(F) rather than
   /// div(F)*mapped_volume. This avoids dividing by mapped_volume after.
   /**
    * @param[out] divF mapped phase space divergence
    * @param[in]  F    flux (assumed same number of component as divF)
    */
   void mappedGridDivergenceFromFluxNormals( LevelData<FArrayBox>&       divF,
                                             const LevelData<FluxBox>&   F ) const;

   void mappedGridDivergenceFromIntegratedFluxNormals( LevelData<FArrayBox>&       divF,
                                                       const LevelData<FluxBox>&   F ) const;

   void fillTransverseGhosts( LevelData<FluxBox>& a_data, const bool do_block_exchange ) const;

   void averageAtBlockBoundaries(LevelData<FluxBox>& data) const;

   void exchangeTransverseAtBlockBoundaries(LevelData<FluxBox>& a_data) const;

   void applyAxisymmetricCorrection( LevelData<FluxBox>& data ) const;

   void multNTransposePointwise( LevelData<FluxBox>& data ) const;

   void getBStarParallel( LevelData<FArrayBox>& BStarParallel ) const;

   void multBStarParallel( LevelData<FArrayBox>& dfn ) const;

   void multBStarParallel( LevelData<FArrayBox>& dfn,
                           const BoundaryBoxLayout& bdry_layout ) const;

   void divideBStarParallel( LevelData<FArrayBox>& dfn ) const;

   /// Fills the argument with the cell-averaged mapping Jacobian
   /**
    * Fills the argument with the cell-averaged mapping Jacobian J
    *
    * @param[out] J  cell-averaged Jacobian
    */
   void getJ( const LevelData<FArrayBox>& configuration_J,
              const LevelData<FArrayBox>& velocityJ,
              LevelData<FArrayBox>&       J ) const;

   /// Multiplies the argument by the cell-averaged mapping Jacobian
   /**
    * Multiplies the argument by the cell-averaged mapping Jacobian J
    *
    * @param[in/out] data  data to be multiplied by cell-averaged Jacobian
    */
   void multJonValid( LevelData<FArrayBox>& dfn ) const;

   void multPointwiseJ( LevelData<FArrayBox>& dfn,
                        const BoundaryBoxLayout& bdry_layout ) const;

   /// Divides the argument by the cell-averaged mapping Jacobian
   /**
    * Divides the argument by the cell-averaged mapping Jacobian J
    *
    * @param[in/out] data  data to be divided by cell-averaged Jacobian
    */
   void divideJonValid( LevelData<FArrayBox>& dfn ) const;

   void getCellVolumes( LevelData<FArrayBox>& volume ) const;

   void fillInternalGhosts( LevelData<FArrayBox>& a_data ) const;
   
   void exchangeExtraBlockGhosts( LevelData<FArrayBox>& a_data ) const;

   const LevelData<FArrayBox>& getBFieldMagnitude() const {return *m_BMagCell;}

   const CFG::MagGeom& magGeom() const {return *m_mag_geom;}

   /// Returns \f$v_{\parallel}-\mu\f$ space mapping.
   /**
    */
   const RefCountedPtr<PhaseCoordSys>& phaseCoordSysPtr() const {return m_phase_coords;}

   const PhaseCoordSys& phaseCoordSys() const {return *m_phase_coords;}

   const RefCountedPtr<VEL::VelCoordSys>& velSpaceCoordSysPtr() const { return m_vel_coords; }

   const VEL::VelCoordSys& velSpaceCoordSys() const { return *m_vel_coords; }

   const PhaseBlockCoordSys& getBlockCoordSys(const Box& box) const;

   const CFG::MagBlockCoordSys& getMagBlockCoordSys(const Box& box) const;

   bool secondOrder() const {return m_second_order;}

   bool divFreeVelocity() const {return m_divfree_velocity;}

   // Multiply by physical face areas (the components in the velocity directions
   // also contain a 2piR factor for axisymmetric problems)
   void multFaceAreas( LevelData<FluxBox>& data ) const;

   void injectConfigurationToPhase( const CFG::LevelData<CFG::FArrayBox>& src,
                                    LevelData<FArrayBox>&                 dst ) const;

   void injectConfigurationToPhase( const CFG::LevelData<CFG::FluxBox>& src,
                                    LevelData<FluxBox>&                 dst ) const;

   void injectConfigurationToPhase( const CFG::LevelData<CFG::FluxBox>&   src,
                                    const CFG::LevelData<CFG::FArrayBox>& src_cell,
                                    LevelData<FluxBox>&                   dst ) const;

   void injectConfigurationToPhase( const CFG::FArrayBox& src,
                                    FArrayBox&            dst ) const;

   void injectConfigurationToPhase( const CFG::FluxBox& src,
                                    FluxBox&            dst ) const;

   void injectVelocityToPhase( const VEL::LevelData<VEL::FArrayBox>& src,
                               LevelData<FArrayBox>&                 dst ) const;

   void injectVelocityToPhase( const VEL::LevelData<VEL::FluxBox>& src,
                               LevelData<FluxBox>&                 dst ) const;

   void projectPhaseToConfiguration( const LevelData<FArrayBox>&     src,
                                     CFG::LevelData<CFG::FArrayBox>& dst ) const;

   void projectPhaseToConfiguration( const LevelData<FluxBox>&     src,
                                     CFG::LevelData<CFG::FluxBox>& dst ) const;
 
   void projectPhaseToConfiguration( const FArrayBox&     a_src,
                                     CFG::FArrayBox&      a_dst ) const;

   void projectPhaseToConfiguration( const Box& src,
                                     CFG::Box&  dst ) const;

   void projectPhaseToVelocity( const Box& src,
                                VEL::Box&  dst ) const;

   CFG::IntVect config_restrict( const IntVect& x ) const;

   IntVect config_inject( const CFG::IntVect& v ) const;

   VEL::IntVect vel_restrict( const IntVect& x ) const;

   IntVect vel_inject( const VEL::IntVect& v ) const;

   IntVect tensorProduct( const CFG::IntVect& cfg_iv,
                          const VEL::IntVect& vel_iv) const;

   /// Get the configuration space LevelData<FluxBox> contained in the
   /// phase space LevelData<FluxBox> at the velocity space coordinate
   /// vspace_coord.
   /**
    * Get the configuration space LevelData<FluxBox> contained in the
    * phase space LevelData<FluxBox> at the velocity space coordinate
    * vspace_index.
    *
    * @param[in]  vspace_index  velocity space index
    * @param[in]  src_phase     source phase space data
    * @param[out] dst_config    destination configuration space data
    */
   void getConfigurationData( const VEL::IntVect& vspace_index,
                              const LevelData<FluxBox>& src_phase,
                              CFG::LevelData<CFG::FluxBox>& dst_config ) const;
   void getConfigurationData( const VEL::IntVect& vspace_index,
                              const LevelData<FArrayBox>& src_phase,
                              CFG::LevelData<CFG::FArrayBox>& dst_config ) const;

   /// Set the configuration space LevelData<FluxBox> contained in the
   /// phase space LevelData<FluxBox> at the velocity space coordinate
   /// vspace_coord.
   /**
    * Set the configuration space LevelData<FluxBox> contained in the
    * phase space LevelData<FluxBox> at the velocity space coordinate
    * vspace_index.
    *
    * @param[in]  vspace_index  velocity space index
    * @param[in]  src_config    source configuration space data
    * @param[out] dst_phase     destination phase space data
    */
   void setConfigurationData( const VEL::IntVect& vspace_index,
                              const CFG::LevelData<CFG::FArrayBox>& src_config,
                              LevelData<FArrayBox>& dst_phase ) const;

   /// Get the velocity space LevelData<FluxBox> contained in the
   /// phase space LevelData<FluxBox> at the configuration space coordinate
   /// cspace_coord.
   /**
    * Get the velocity space LevelData<FluxBox> contained in the
    * phase space LevelData<FluxBox> at the configuration space coordinate
    * cspace_index.
    *
    * @param[in]  cspace_index  configuration space index
    * @param[in]  src_phase     source phase space data
    * @param[out] dst_vel       destination velocity space data
    */
   void getVelocityData( const CFG::IntVect& cspace_index,
                         const LevelData<FluxBox>& src_phase,
                         VEL::LevelData<VEL::FluxBox>& dst_vel ) const;
   void getVelocityData( const CFG::IntVect& cspace_index,
                         const LevelData<FArrayBox>& src_phase,
                         VEL::LevelData<VEL::FArrayBox>& dst_vel ) const;

   /// Set the velocity space LevelData<FluxBox> contained in the
   /// phase space LevelData<FluxBox> at the configuration space coordinate
   /// cspace_coord.
   /**
    * Set the velocity space LevelData<FluxBox> contained in the
    * phase space LevelData<FluxBox> at the configuration space coordinate
    * cspace_index.
    *
    * @param[in]  cspace_index  configuration space index
    * @param[in]  src_vel       source phase space data
    * @param[out] dst_phase     destination velocity space data
    */
   void setVelocityData( const CFG::IntVect& cspace_index,
                         const VEL::LevelData<VEL::FArrayBox>& src_vel,
                         LevelData<FArrayBox>& dst_phase ) const;

   IntVect ghostVect() const {return m_ghostVect;}

   /*
    *  Plotting routines
    */

   void plotAtVelocityIndex(const string               file_name,
                            const VEL::IntVect&        vspace_index,
                            const LevelData<FluxBox>&  data,
                            const double&              time) const;

   void plotAtVelocityIndex(const string                 file_name,
                            const VEL::IntVect&          vspace_index,
                            const LevelData<FArrayBox>&  data,
                            const double&                time) const;

   void plotAtConfigurationIndex(const string               file_name,
                                 const CFG::IntVect&        cspace_index,
                                 const LevelData<FluxBox>&  data,
                                 const double&              time) const;

   void plotAtConfigurationIndex(const string                 file_name,
                                 const CFG::IntVect&          cspace_index,
                                 const LevelData<FArrayBox>&  data,
                                 const double&                time) const;

   void plotAtMuIndex( const string                file_name,
                       const int                   mu_index,
                       const LevelData<FArrayBox>& data,
                       const double&               time) const;

   void plotConfigurationData(const string                          file_name,
                              const CFG::LevelData<CFG::FArrayBox>& data,
                              const double&                         time) const;

   void plotConfigurationData(const string                        file_name,
                              const CFG::LevelData<CFG::FluxBox>& data,
                              const double&                       time) const;

   void plotVelocityData(const string                          file_name,
                         const VEL::LevelData<VEL::FArrayBox>& data,
                         const double&                         time) const;

   void plotVelocityAtXYIndex(const string                 file_name,
                              const CFG::IntVect&          cspace_index,
                              const LevelData<FArrayBox>&  data,
                              const double&                time) const;

   void plotVParPoloidalData(const string                file_name,
                             const int                   radial_index,
                             const int                   toroidal_index,
                             const int                   mu_index,
                             const LevelData<FArrayBox>& data,
                             const double&               time) const;

   void plotDivergence( const string&             file_name,
                        const LevelData<FluxBox>& flux,
                        const double              time ) const;

   void plotData( const string&               file_name,
                  const LevelData<FArrayBox>& data,
                  const double&               time) const;

   CFG::Box configurationDomainBox(const CFG::DisjointBoxLayout& grids) const;

   VEL::Box velocityDomainBox(const VEL::DisjointBoxLayout& grids) const;

   bool validConfigurationIndex(const CFG::IntVect& index) const;

   bool validVelocityIndex(const VEL::IntVect& index) const;

   const ProblemDomain& domain() const {return m_domain;}

   const PhaseGrid& phaseGrid() const {return *m_phase_grid;}

   //Extract configuraion component of a 4D vector (e.g., GKVelocity)
   void getConfigurationComponents( LevelData<FArrayBox>& configComp,
                                    const LevelData<FArrayBox>& vector) const;
   
   //Compute vector projection onto the direction of phi x bpol
   void computeRadialProjection( LevelData<FArrayBox>& radComp,
                                 const LevelData<FArrayBox>& vector) const;

   //Compute vector projection onto the direction of bpol
   void computePoloidalProjection( LevelData<FArrayBox>& poloidalComp,
                                   const LevelData<FArrayBox>& vector) const;

   enum velocity_option {FULL_VELOCITY, PARALLEL_STREAMING_VELOCITY, DRIFT_VELOCITY, EXB_DRIFT_VELOCITY,
                         MAGNETIC_DRIFT_VELOCITY, NO_ZERO_ORDER_PARALLEL_VELOCITY};
   
   double getLarmorNumber() const {return m_larmor_number;}

   inline Real larmorNumber() const { return m_larmor_number; }

   inline void gyroaverageOp(const GyroaverageOperator* const a_op) const
   {
     m_gyroavg_op = a_op;
   }

   inline const GyroaverageOperator* const gyroaverageOp() const
   {
     return m_gyroavg_op;
   }

   inline void gyroaverage(CFG::LevelData<CFG::FArrayBox>& a_phi_bar,
                           const CFG::LevelData<CFG::FArrayBox>& a_phi) const
   {
     if (m_gyroavg_op == NULL) {
       MayDay::Error("Error in PhaseGeom::gyroaverage() - m_gyroavg_op is NULL.");
     }
     m_gyroavg_op->applyOp(a_phi_bar, a_phi);
   }

protected:

   Box configurationFlatten( const Box& disjointBox,
                             const Box& src_box ) const;

   Box velocityFlatten( const Box& disjointBox,
                        const Box& src_box ) const;

   IntVect m_ghostVect;

   const RefCountedPtr<CFG::MagGeom>     m_mag_geom;
   const RefCountedPtr<VEL::VelCoordSys> m_vel_coords;

private:

   const RefCountedPtr<PhaseCoordSys> m_phase_coords;

   const RefCountedPtr<PhaseGrid> m_phase_grid;

   ProblemDomain m_domain;

   LevelData<FluxBox> *   m_configuration_metrics;
   LevelData<FluxBox> *   m_configuration_tangrad_metrics;
   LevelData<FArrayBox> * m_configuration_volumes;
   LevelData<FArrayBox> * m_configuration_J;
   LevelData<FluxBox> *   m_configuration_face_areas;

   LevelData<FluxBox> *   m_velocity_metrics;
   LevelData<FluxBox> *   m_velocity_tangrad_metrics;
   LevelData<FArrayBox> * m_velocity_volumes;
   LevelData<FArrayBox> * m_velocity_J;
   LevelData<FluxBox> *   m_velocity_face_areas;

   LevelData<FArrayBox> * m_BCell;
   LevelData<FluxBox> * m_BFace;
   LevelData<FArrayBox> * m_gradBCell;
   LevelData<FluxBox> * m_gradBFace;
   LevelData<FArrayBox> * m_curlbCell;
   LevelData<FluxBox> * m_curlbFace;
   LevelData<FArrayBox> * m_BMagCell;
   LevelData<FluxBox> * m_BMagFace;
   LevelData<FArrayBox> * m_bdotcurlbCell;
   LevelData<FluxBox> * m_bdotcurlbFace;

   LevelData<FArrayBox> m_BStar;
   LevelData<FArrayBox> m_BStarParallel;
   LevelData<FArrayBox> m_BStarParallel_cell_averaged;

   MultiBlockLevelExchangeAverage* m_mblexPtr;
   BlockRegister* m_exchange_transverse_block_register;

   string m_velocity_type;
   bool m_divfree_velocity;
   bool m_no_drifts;
   bool m_no_parallel_streaming;
   bool m_speciesDefined;
   RealVect m_freestream_components;

   bool m_second_order;   // experimental

   double m_larmor_number;
   double m_mass;
   double m_charge_state;
   
#if CFG_DIM == 3
   
   void getShearedGhostBoxLayout();
   void interpolateFromShearedGhosts(LevelData<FArrayBox>& a_data) const;
   
   BoxLayout m_shearedGhostBLLoEnd;
   BoxLayout m_shearedGhostBLHiEnd;

#endif   
   LevelData<FArrayBox> * m_sheared_remapped_index;
   LevelData<FArrayBox> * m_sheared_interp_stencil;
   
   mutable GyroaverageOperator const *m_gyroavg_op;
};

typedef Vector<RefCountedPtr<PhaseGeom> > PhaseGeomPtrVect;

typedef RefCountedPtr<PhaseGeom>  PhaseGeomPtr;

#include "NamespaceFooter.H"

#endif
