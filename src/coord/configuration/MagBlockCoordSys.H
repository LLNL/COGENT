#ifndef _MAGBLOCKCOORDSYS_H_
#define _MAGBLOCKCOORDSYS_H_

#include "ParmParse.H"
#include "DisjointBoxLayout.H"
#include "NewFourthOrderCoordSys.H"

#include "NamespaceHeader.H"

/// Magnetic Geometry abstract base class
/**
 * Abstract base class for flux-surface-aligned, magnetic geometry-based
 * coordinate system.
*/
class MagBlockCoordSys
   : public NewFourthOrderCoordSys
{
   public:

      /// Constructors.
      /**
       * Default (required for an abstract base class)
       *
       */
      MagBlockCoordSys() {;}

      /**
       * Constructor with initialization.  Input data is read from
       * the parm_parse object, the cellSpacing is computed and
       * the define() member is called.
       *
       * @param[in] parm_parse the ParmParse database.
       */
      MagBlockCoordSys( ParmParse& parm_parse );

      /// Destructor
      /**
       */
      virtual ~MagBlockCoordSys() {;}

      /**
       * Defines the object.
       *
       * @param[in] domain_box the mapped index domain box
       * @param[in] cellSpacing the computational space mesh size.
       */
      virtual void define(const ProblemDomain& domain,
                          const RealVect& cellSpacing);

      /// Returns real coordinate.
      /**
       * Returns location \f$x\f$ in real-world coordinates of
       * mapped-coordinate location \f$\xi\f$.
       *
       * @param[in] xi mapped grid coordinate.
       */
      virtual RealVect realCoord( const RealVect& xi ) const = 0;

      /// Returns computational coordinate.
      /**
       * Returns location \f$\xi\f$ in mapped coordinate space of real-world
       * location \f$x\f$.
       *
       * @param[in] x physical space coordinate.
       */
      virtual RealVect mappedCoord( const RealVect& x ) const = 0;

      /// Returns the derivatives of the physical coordinates with respect to
      /// the computational coordinate.
      /**
       * Returns the derivatives of the physical coordinates with respect to
       * the computational coordinates at location \f$xi\f$.
       *
       * @param[in] xi mapped coordinate.
       */
      virtual Real dXdXi(const RealVect& Xi, int dirX, int dirXi) const = 0;

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
      virtual void dXdXi(FArrayBox& dxdXi,
                         const FArrayBox& Xi,
                         int destComp,
                         int dirX,
                         int dirXi,
                         const Box& box) const = 0;

      virtual void cellVol(FArrayBox& vol, const FluxBox& N, const Box& box) const;

      virtual void getCellCenteredMappedCoords(FArrayBox& xi) const;

      virtual void getFaceCenteredMappedCoords(const int dir, FArrayBox& xi) const;

      virtual void pointwiseJ(FArrayBox& a_J, const FArrayBox& a_Xi, const Box& a_box) const;

#if 1
void
      virtual incrementFaceMetricWithEdgeTerm(FArrayBox& a_N,
                                              int a_faceDir,
                                              int a_edgeDir,
                                              const Box& a_box,
                                              const Side::LoHiSide& a_side) const;
#endif

#if 1
      virtual void computeTransverseFaceMetric(FArrayBox& a_faceMetrics,
                                               const Box& a_box,
                                               int a_faceDir,
                                               int a_dDir) const;
#endif

      /// Returns major radius
      /**
       * Returns the major radius at mapped-coordinate location \f$\xi\f$.
       *
       * @param[in] xi mapped grid coordinate.
       */
      Real majorRadius( const RealVect& xi ) const {return realCoord(xi)[0];}

      /// Returns major radius adjusted for constant minor radius option
      /**
       * Returns the major radius at mapped-coordinate location \f$\xi\f$ adjusted
       * for the constant minor radius option (in a derived class).
       *
       * @param[in] xi mapped grid coordinate.
       */
      virtual double majorRadiusConstMinorRadAdj(const RealVect a_xi) const {return majorRadius(a_xi);}

      /// Returns the mapped cell volume
      /**
       * Returns the mapped cell volume
       *
       */
      Real getMappedCellVolume() const {return m_mapped_cell_volume;}

      /// Returns the areas of mapped cell faces
      /**
       * Returns the areas of mapped cell faces
       *
       */
      RealVect getMappedFaceArea() const {return m_mapped_face_area;}

      double lowerMappedCoordinate(int dir) const {
         return m_dx[dir] * m_domain.domainBox().smallEnd(dir);
      };

      double upperMappedCoordinate(int dir) const {
         return m_dx[dir] * (m_domain.domainBox().bigEnd(dir)+1);
      };

      void getPointwiseMajorRadius( FluxBox& a_R ) const;

      /**
       * Fill argument with the derivatives dXdxi evaluated at cell centers
       *
       * @param[out] dXdxi  dXdxi
      */
      void getCellCentereddXdxi(FArrayBox& dXdxi) const;

      /**
       * Fill argument with the derivatives dXdxi evaluated at face centers
       *
       * @param[out] dXdxi   dXdxi
      */
      void getFaceCentereddXdxi(FluxBox& dXdxi) const;

      /**
       * Fill argument with the metric quantity N evaluated at cell centers
       *
       * @param[out] N  Metric N
      */
      void getPointwiseN(FArrayBox& N) const;

      /**
       * Fill argument with the metric quantity N evaluated at face centers
       *
       * @param[out] N   Metric N
      */
      void getPointwiseN(FluxBox& N) const;

      /**
       * Fill argument with the metric quantity NJInverse evaluated at cell centers
       *
       * @param[out] NJInverse  Product of N and JInverse
      */
      void getPointwiseNJInverse(FArrayBox& NJInverse) const;

      /**
       * Fill argument with the metric quantity NJInverse evaluated at face centers
       *
       * @param[out] NJInverse  Product of N and JInverse
      */
      void getPointwiseNJInverse(FluxBox& NJInverse) const;

      /**
       * Fill argument with 2 pi times the minor radius times the major radius
       *
       * @param[out] TwoPirRmaj  TwoPi * r * Rmaj
      */
      void getTwoPirRmaj(FArrayBox& TwoPirRmaj) const;

      /**
       * Returns the value of the innermost flux label
       */
      double getInnerFluxLabel() const {return m_inner_radial_boundary;}

      /// Fills the input FArrayBox with the real coordinates at cell centers
      /**
       */
      void getCellCenteredRealCoords( FArrayBox& cell_centered_coords ) const;

      /// Fills the input FArrayBox with the real coordinates at dir face centers
      /**
       */
      virtual void getFaceCenteredRealCoords( const int dir, FArrayBox& a_x ) const;

      /// Fills the input FArrayBox with the real coordinates at cell nodes
      /**
       */
      void getNodeCenteredRealCoords( FArrayBox&  a_X ) const;

      /// Fills the input FArrayBox with the flux coordinates at cell centers
      /**
       */
      void getCellCenteredFluxCoords( FArrayBox& cell_centered_coords ) const;

      /// Computes the flux surface unit tangent vector
      /**
       * @param[out] data  flux surface unit tangent vector
       */
      void computeFluxSurfaceUnitTangent( FArrayBox& data ) const;
      void computeFluxSurfaceUnitTangent( FluxBox& data ) const;

      /// Computes the flux surface outward unit normal
      /**
       * @param[out] data  flux surface outward unit normal
       */
      void computeFluxSurfaceUnitNormal( FArrayBox& data ) const;
      void computeFluxSurfaceUnitNormal( FluxBox& data ) const;

      /// Computes the unit vector in the direction of the
      /// mapped grid radial coordinate
      /**
       * @param[out] data  radial unit vector
       */
      void computeUnitRadial( FArrayBox& data ) const;
      void computeUnitRadial( FluxBox& data ) const;

      /// Project the vector argument onto the flux surface
      /**
       * @param[in] data    vector to be projected
       * @param[out] data   projected vector
       */
      void projectOntoFluxSurface( FArrayBox& data ) const;
      void projectOntoFluxSurface( FluxBox& data ) const;

      /// Project the vector argument onto the flux surface unit normal
      /**
       * @param[in] data    vector to be projected
       * @param[out] data   projected vector
       */
      void projectOntoFluxSurfaceNormal( FArrayBox& data ) const;
      void projectOntoFluxSurfaceNormal( FluxBox& data ) const;

      /// Project the vector argument onto the unit vector in the
      /// direction of the mapped radial coordinate
      /**
       * @param[in] data    vector to be projected
       * @param[out] data   projected vector
       */
      void projectOntoRadial( FArrayBox& data ) const;
      void projectOntoRadial( FluxBox& data ) const;

      /// Compute projections of the configurational components of vector 
      /// argument onto the psi(normal) and theta (tangental) directions
      /**
       * @param[in] data    vector to be projected
       * @param[out] data   projected vector
       */
      void computePsiThetaProjections( FArrayBox& data ) const;
      void computePsiThetaProjections( FluxBox& data ) const;


      void computeGradFFactor( FArrayBox& a_data ) const;
      void computeGradFFactor( FluxBox& a_data ) const;

      /**
       * Returns the value of the outermost flux label
       */
      double getOuterFluxLabel() const {return m_outer_radial_boundary;}

      /// Returns mapped-space problem domain
      const ProblemDomain& domain() const {return m_domain;}

      virtual double getRBtoroidal() const = 0;

      virtual void computeFieldData( const int  dir,
                                     FArrayBox& BField,
                                     FArrayBox& BFieldMag,
                                     FArrayBox& BFieldDirl,
                                     FArrayBox& gradBFieldMag,
                                     FArrayBox& curlBFieldDir,
                                     FArrayBox& BFieldDirdotcurlBFieldDir,
                                     const bool derived_data_only = false ) const = 0;

      virtual void getMagneticFlux( const FArrayBox& physical_coordinates,
                                    FArrayBox&       magnetic_flux ) const = 0;

      virtual double getMagneticFlux( const RealVect& a_physical_coordinate ) const = 0;

      bool isAxisymmetric() const {return m_axisymmetric;}

      virtual string geometryType () const = 0 ;

   protected:

      ProblemDomain m_domain;
      Real m_mapped_cell_volume;
      RealVect m_mapped_face_area;
      bool m_axisymmetric;
      bool m_verbose;
      double m_inner_radial_boundary;
      double m_outer_radial_boundary;
};



#include "NamespaceFooter.H"

#endif
