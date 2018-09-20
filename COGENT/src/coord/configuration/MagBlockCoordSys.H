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

      virtual void getNodeCenteredMappedCoords(FArrayBox& xi) const;

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
   
      RealVect getMappedCellSize() const {return m_dx;}

      void getPointwiseMajorRadius( const FArrayBox& xi,
                                    FArrayBox&       R ) const;

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
       * Fill argument with the metric quantity NJInverse evaluated at a given mapped coordinate
       *
       * @param[out] NJInverse  Product of N and JInverse
       */

       Vector<Real> getPointwiseNJInverse(const RealVect& Xi) const;

      /**
       * Fill argument with the metric quantity NT evaluated at a given mapped coordinate
       *
       * @param[out] NT transpose of the 3D N (used in divergence calculation)
       */

      Vector<Real> getPointwiseNT(const RealVect& Xi) const;
   
  
      /**
       * Fill argument with 2 pi times the minor radius times the major radius
       *
       * @param[out] TwoPirRmaj  TwoPi * r * Rmaj
      */
      void getTwoPirRmaj(FArrayBox& TwoPirRmaj) const;

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


      /// Compute projections of the configurational components of vector 
      /// argument onto the psi(normal) and theta (tangental) directions
      /**
       * @param[in] data    vector to be projected
       * @param[out] data   projected vector
       */
      void computePsiThetaProjections( FArrayBox& data ) const;
      void computePsiThetaProjections( FluxBox& data ) const;

      /// Returns mapped-space problem domain
      const ProblemDomain& domain() const {return m_domain;}

      virtual double getRBtoroidal() const = 0;

      /// Returns B field at a given physical coordinate
      virtual Vector<Real> computeBField(const RealVect& a_X) const = 0;
   
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

      virtual void getNodalFieldData(FArrayBox& points,
                                     FArrayBox& A,
                                     FArrayBox& b,
                                     FArrayBox& Bmag) const = 0;

      bool isAxisymmetric() const {return m_axisymmetric;}
      bool isFieldAligned() const {return m_field_aligned;}
      bool isPointwiseMetrics() const {return m_pointwise_metrics;}
      bool isSubGridGeom() const {return m_subgrid_geometry;}
      IntVect getSubGridData() const {return m_N_mr;}
   
      virtual string geometryType () const = 0 ;
   
   private:
   
      /// Inline method to compute difference of products needed for cofactors.
      /**
       * comp1*comp2 minus comp3*comp4 of a_tmp, put into destcomp component of cofarr
      */
      void TRM1x2less3x4(FArrayBox& a_out,
			 const FArrayBox& a_in,
			 const Box& box, int destcomp,
			 int i1, int i2, int i3, int i4) const;

      Real TRM1x2less3x4(const Vector<Real>&  a_in,
			 int i1, int i2, int i3, int i4) const;

   
      /// Method to compute terms in dXdXi, in dXdXiArr
      void compute_dXdXi ( FArrayBox& dXdXiArr,
			   const FArrayBox& xi,
			   const Box& box  ) const;

      Vector<Real> compute_dXdXi (const RealVect& Xi) const;
      
      /// Method to compute cofactor terms of dXdXi, in cofarr
      void compute_dXdXiCofac ( FArrayBox& cofarr,
				const FArrayBox& xi,
				const Box& box  ) const;

      Vector<Real> compute_dXdXiCofac ( const RealVect& Xi) const;

      /// Computes the physical cell volume by a subgrid method
      /**
       * Computes the physical cell volume by refining the grid and
       * then summing over refined (subGrid) elements
       *
      */
      void cellVolSubGrid(FArrayBox& vol, const FArrayBox& Xi) const;

   protected:

      ProblemDomain m_domain;
      Real m_mapped_cell_volume;
      RealVect m_mapped_face_area;
      bool m_axisymmetric;
      bool m_field_aligned;
      bool m_pointwise_metrics;
      bool m_verbose;
      bool m_subgrid_geometry;
      IntVect m_N_mr;
};



#include "NamespaceFooter.H"

#endif
