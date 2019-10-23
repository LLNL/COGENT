#ifndef _SINGLENULLBLOCKCOORDSYS_H_
#define _SINGLENULLBLOCKCOORDSYS_H_

#undef PLOT_INVERSE_POINTS

#include "MagBlockCoordSys.H"
#include "Dimensions.H"

#if CFG_DIM==3
#include "BSplineInterp3D.H"
#include "LevelData.H"
#endif

#define POL_NAMESPACE D2
#define POL_DIM 2

#undef CH_SPACEDIM
#define CH_SPACEDIM POL_DIM
#include "SingleNullPoloidalUtil.H"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H"
namespace POL = POL_NAMESPACE;


class SingleNullBlockCoordSys
   : public MagBlockCoordSys
{
   public:

      /// Constructor
      /**
       * Constructor with initialization.
       *
       * @param[in] pp               the ParmParse database.
       * @param[in] domain           the configuration space domain.
       * @param[in] dx               mapped grid cell size
       * @param[in] poloidal_block   poloidal block
       * @param[in] toroidal_sector  toroidal sector
       */
      SingleNullBlockCoordSys( ParmParse&            pp,
                               const ProblemDomain&  domain,
                               const RealVect&       dx,
                               const int             poloidal_block,
                               const int             toroidal_sector,
                               const int             toroidal_block_separation );

      /// Destructor
      /**
       */
      virtual ~SingleNullBlockCoordSys();

      virtual void printInit() const;

      virtual bool blockNameIsValid( const string& a_block_name ) const;

      void readFiles( ParmParse& pp );

      virtual void definePoints( const ParmParse&     pp,
                                 int&                 block_poloidal,
                                 const int            block_full_poloidal,
                                 const POL::IntVect&  mapping_block_size,
                                 const int&           n_poloidal_extend,
                                 double &             dtheta,
                                 double *             theta_pts ) const;

      /// Returns real coordinate.
      /**
       * Returns location \f$x\f$ in real-world coordinates of
       * mapped-coordinate location \f$\xi\f$.
       *
       * @param[in] xi mapped grid coordinate.
       */
      virtual RealVect realCoord( const RealVect& xi ) const;

      /// Returns computational coordinate.
      /**
       * Returns location \f$\xi\f$ in mapped coordinate space of real-world
       * location \f$x\f$.
       *
       * @param[in] x physical space coordinate.
       */
      virtual RealVect mappedCoord( const RealVect& x ) const;

      /// Returns the derivatives of the physical coordinates with respect to
      /// the computational coordinate.
      /**
       * Returns the derivatives of the physical coordinates with respect to
       * the computational coordinates at location \f$xi\f$.
       *
       * @param[in] xi mapped coordinate.
       */
      virtual Real  dXdXi(const RealVect& Xi, int dirX, int dirXi) const;

      /// Fills the destComp component of a_dxdXi with the derivative of
      /// x w/ respect to Xi in the dirX direction
      /**
       * Fills the destComp component of a_dxdXi with the derivative of
       * x w/ respect to Xi in the dirX direction
       *
       * @param[in] Xi mapped coordinates
       * @param[in] destComp destination component
       * @param[in] dirX physical coordinate to be differentiated
       * @param[in] dirXi computational coordinate
       * @param[in] box box within which derivatives are computed
       * @param[out] dxdXi output derivatives
       */
      void dXdXi(FArrayBox&       dxdXi,
                 const FArrayBox& Xi,
                 int              destComp,
                 int              dirX,
                 int              dirXi,
                 const Box&       box) const;

#if CFG_DIM==3
      void setInterp( const ParmParse&       pp,
                      const POL::FArrayBox&  interp_node_coords_pol,
                      const POL::FArrayBox&  RZ_data );

      RealVect traceField( const RealVect&  X,
                           const double&    toroidal_increment ) const;

      RealVect mappedCoordNewton( const RealVect&  X,
                                  const RealVect&  xi_initial,
                                  const IntVect&   iv_initial ) const;
   
      void getPoloidalDisjointBoxLayout(DisjointBoxLayout&   a_grids,
                                        Vector<Box>&         a_boxes,
                                        int&                 a_decomp_num,
                                        const Box&           a_box) const;

      void assembleDecomposedData(FArrayBox&                    a_interp_node_data,
                                  FArrayBox&                    a_interp_node_coords,
                                  const LevelData<FArrayBox>&   a_interp_node_data_decomp,
                                  const LevelData<FArrayBox>&   a_interp_node_coords_decomp,
                                  const Vector<Box>&            a_grid_boxes,
                                  const int&                    a_decomp_num ) const;
#endif

      bool containsPoloidalPoint( const POL::RealVect& X, const bool verbose ) const;
   
      virtual array<double,3> computeBField(const RealVect& a_X) const;

      virtual void computeFieldData( const int  dir,
                                     FArrayBox& BField,
                                     FArrayBox& BFieldMag,
                                     FArrayBox& BFieldDir,
                                     FArrayBox& gradBFieldMag,
                                     FArrayBox& curlBFieldDir,
                                     FArrayBox& BFieldDirdotcurlBFieldDir,
                                     const bool derived_data_only = false ) const;

      //Computes magnetic flux Psi
      virtual void getMagneticFlux( const FArrayBox& a_physical_coordinates,
                                    FArrayBox&       a_magnetic_flux ) const;

      virtual double getMagneticFlux( const RealVect& a_physical_coordinate ) const;
   
      virtual void getNodalFieldData(FArrayBox& points, FArrayBox& A, FArrayBox& b, FArrayBox& Bmag) const;

      ///Return toroidal coordinates
      void getToroidalCoords(FArrayBox& coords) const;
  
      //Computes normalized magnetic flux
      virtual void getNormMagneticFlux( const FArrayBox& a_physical_coordinates,
                                        FArrayBox&       a_magnetic_flux ) const;
   
      virtual double getNormMagneticFlux( const RealVect& a_physical_coordinate ) const;
   
      //Returns the physical radial coordinate (R) for a given value of normalized psi
      virtual double getOuterMidplaneCoord(const double psiNorm) const;
   
      //Returns the radial coordinate of the separatrix at the outer midplane
      virtual double getOuterRsep() const;
   
      virtual double getRBtoroidal() const {return m_RB_toroidal;}

#if CFG_DIM==3
      double getToroidalMappedOffset() const {return m_toroidal_mapped_offset;}
#endif

      static const std::string pp_name;

      virtual int poloidalBlock() const {return m_poloidal_block;}

      virtual void getCellCenterRealCoords();

      bool truncated() const {return m_poloidally_truncated;}

      bool isValid(const RealVect& xi, const bool omit_toroidal ) const;
   
      bool isFieldAlignedMapping() const {return m_field_aligned_mapping;}

      virtual string geometryType() const {return "SingleNull";}
   
      virtual RealVect getMagAxis() const {return injectToroidal(m_poloidal_util->getMagAxis());}
   
      virtual RealVect getXpoint() const {return injectToroidal(m_poloidal_util->getXpoint());}

      virtual string blockName( int poloidal_block ) const;

      enum SingleNullBlockType {LCORE,RCORE,LCSOL,RCSOL,LSOL,RSOL,LPF,RPF,MCORE,MCSOL,NUM_POLOIDAL_BLOCKS};

      enum PhysCoordType {CYLINDRICAL, CARTESIAN};

      void convertCylindricalToCartesian( const RealVect&   cartesian_vect,
                                          array<double,3>&  vect ) const;

      void convertCylindricalToCartesian( const FArrayBox&  cartesian_coord,
                                          FArrayBox&        coord ) const;

      void convertCartesianToToroidal(RealVect& coord) const;

      // Removes the toroidal component of the input if there is one, otherwise just copies
      virtual POL::IntVect restrictToPoloidal( const IntVect& v ) const;

      // Removes the toroidal component of the input if there is one, otherwise just copies
      virtual POL::RealVect restrictToPoloidal( const RealVect& v ) const;

      // Returns a box with the toroidal dimension removed if there is one, otherwise just copies
      virtual POL::Box restrictToPoloidal( const Box& box ) const;

      // Restricts the input physical coordinate to the poloidal plane
      virtual POL::RealVect restrictPhysCoordToPoloidal( const RealVect& a_X ) const;

      // Restricts the input physical coordinates to the poloidal plane and returns
      // the result in the second argument
      virtual void restrictPhysCoordsToPoloidal( const FArrayBox&  coords,
                                                 POL::FArrayBox&   coords_pol ) const;

      virtual void restrictPhysCoordsToPoloidal( const FArrayBox&  coords,
                                                 const int         toroidal_index,
                                                 POL::FArrayBox&   coords_pol ) const;

   // Restricts the input mapped coordinates to the poloidal plane and returns
      // the result in the second argument
      virtual void restrictMappedCoordsToPoloidal( const FArrayBox&  coords,
                                                   POL::FArrayBox&   coords_pol ) const;

      // Spreads the input in the toroidal direction if there is one; otherwise copies the input
      // to the output
      virtual RealVect injectToroidal( const POL::RealVect  a_v ) const;

      void injectToroidal( const POL::FArrayBox&  fab_pol,
                           const int              toroidal_index,
                           FArrayBox&             fab ) const;

   protected:

      int m_poloidal_block;
      int m_toroidal_sector;

#if CFG_DIM==2
      POL::FArrayBox m_realPoloidalCoords;
      POL::IntVect m_rpc_coarsen_ratio;
#endif
#if CFG_DIM==3
      FArrayBox m_realCoords;
      IntVect m_rc_coarsen_ratio;
#endif
   
      double m_RB_toroidal;

      bool m_spectral_field;

      POL::SingleNullPoloidalUtil* m_poloidal_util;

      bool m_field_aligned_mapping;
   
   private:
   
      bool m_poloidally_truncated;

      PhysCoordType m_phys_coord_type;
   
#ifdef PLOT_INVERSE_POINTS
      FILE *m_ipt_file[SingleNullBlockCoordSys::NUM_POLOIDAL_BLOCKS];
#endif

#if CFG_DIM==3
      Interp3D* m_interp;

      double m_toroidal_mapped_offset;
#endif

      int m_toroidal_ghosts;   // Number of mapped ghost cells in the toroidal direction
      int m_toroidal_mapping_refinement; //Refinemnt of mapping in toroidal dir relative to grid
   
      // Field tracing parameters; field lines are traced until they hit a prescribed poloidal plane.
      int m_field_trace_step_num;     //Number of steps in toroidal dir used to ray trace field between poloidal planes
};

#include "NamespaceFooter.H"

#endif
