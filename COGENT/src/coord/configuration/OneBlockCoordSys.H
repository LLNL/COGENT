#ifndef _ONEBLOCKCOORDSYS_H_
#define _ONEBLOCKCOORDSYS_H_

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

/**
 * arbitrary one block geometry class
 *
 * Nomenclature: indices are still RADIAL_DIR <--> R,  TOROIDAL_DIR <--> Y,  POLOIDAL_DIR <--> Z,
 * to maintain contact with toroidal simulations
*/
class OneBlockCoordSys
   : public MagBlockCoordSys
{
   public:

   // Constructor
   OneBlockCoordSys(ParmParse&            parm_parse,
                    const ProblemDomain&  domain,
                    const int             numcells_toroidal,
                    const int             toroidal_sector,
                    const int             toroidal_block_separation);

   // Destructor
   virtual ~OneBlockCoordSys();

   void readFiles( ParmParse& pp );
      
   virtual RealVect realCoord( const RealVect& Xi ) const;
   
   virtual RealVect mappedCoord( const RealVect& X ) const;
   
   void getCellCenterRealCoords();
   
   bool isValid(const RealVect& Xi, const bool omit_toroidal) const;									 
   virtual Real dXdXi(const RealVect& Xi, int dirX, int dirXi) const;
 
   virtual void dXdXi( FArrayBox&  dxdXi,
                 const FArrayBox&  Xi,
                       int         destComp,
                       int         dirX,
                       int         dirXi,
                 const Box&        box ) const;

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
   
   //
   //  below are functions related to the magnetic field 
   //
   
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
  
   virtual bool isFieldAlignedMapping() const {return m_field_aligned_mapping;}
   
  
   static const std::string pp_name;
  
   virtual string geometryType() const {return "OneBlock";}

   double rmin() const {return m_rmin;}
   double rmax() const {return m_rmax;}
   
   double thetamin() const {return m_thetamin;}
   double thetamax() const {return m_thetamax;}
   
   virtual RealVect getMagAxis() const {return injectToroidal(m_poloidal_util->getMagAxis());}

   virtual RealVect getXpoint() const {return injectToroidal(m_poloidal_util->getXpoint());}
   
   enum PhysCoordType {CYLINDRICAL, CARTESIAN};

   void convertCylindricalToCartesian( const RealVect&   cartesian_vect,
                                       array<double,3>&  vect ) const;

   void convertCylindricalToCartesian( const FArrayBox&  cartesian_coord,
                                       FArrayBox&        coord ) const;

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

      double m_rmin, m_rmax;
      double m_thetamin, m_thetamax;
#if CFG_DIM==3
      double m_toroidal_width;
#endif
   
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
      bool m_mapping_includes_RBdata;
   
   private:
   
      PhysCoordType m_phys_coord_type;
   
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
