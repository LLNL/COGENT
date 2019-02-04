#ifndef _TOROIDALBLOCKCOORDSYS_H_
#define _TOROIDALBLOCKCOORDSYS_H_

#if CFG_DIM==3

#include "MagBlockCoordSys.H"

#include "NamespaceHeader.H"

/// Toroidal magnetic geometry class
/**
 * Implements an analytic toroidal model with concentric circular flux surfaces
 *
 * only 3D case is supported:
 * physical coordinates X <-> x[0], Y <-> x[1], Z <-> x[2]
 * mapped coordinates minor radius (xi[0]), toroidal_angle (xi[1]), poloidal angle (xi[2])
 *
 * If field-algined option is true:
 * the toroidal angle measures distanse along field lines (a.k.a. field-algined coordinate)
 * the poloidal angle labels the field line in the poloidal plane going through the middle of a block
 *
 *
 * separation_shift = block_separation * dx[1] * block_number
*/

class ToroidalBlockCoordSys
   : public MagBlockCoordSys
{
   public:

      /// Constructor
      /**
       * Constructor with initialization.  Input data is read from
       * the parm_parse object, the cellSpacing is computed and
       * the define() member is called.
       *
       * @param[in] parm_parse the ParmParse database.
       * @param[in] domain the configuration space domain.
       * @param[in] block block id.
       * @param[in] num_blocks total number or blocks.
       * @param[in] number of separating cells (between blocks).
       */
   
      ToroidalBlockCoordSys( ParmParse&               parm_parse,
                             const ProblemDomain&     domain,
                             const int                block,
                             const int                num_blocks,
                             const int                block_separation);


      /// Destructor
      /**
       */
      virtual ~ToroidalBlockCoordSys() {;}

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
      virtual void dXdXi(FArrayBox& dxdXi,
                         const FArrayBox& Xi,
                         int destComp,
                         int dirX,
                         int dirXi,
                         const Box& box) const;
   
      virtual Vector<Real> computeBField(const RealVect& a_X) const;

      virtual void computeFieldData( const int  dir,
                                     FArrayBox& BField,
                                     FArrayBox& BFieldMag,
                                     FArrayBox& BFieldDir,
                                     FArrayBox& gradBFieldMag,
                                     FArrayBox& curlBFieldDir,
                                     FArrayBox& BFieldDirdotcurlBFieldDir,
                                     const bool derived_data_only = false ) const;

      virtual void getMagneticFlux( const FArrayBox& physical_coordinates,
                                    FArrayBox&       magnetic_flux ) const;

      virtual double getMagneticFlux( const RealVect& a_physical_coordinate ) const;
   
      Real getMagneticFlux( const Real r ) const;
   
      virtual void getNodalFieldData(FArrayBox& points, FArrayBox& A, FArrayBox& b, FArrayBox& Bmag) const;

      ///Return toroidal coordinates
      void getToroidalCoords(FArrayBox& coords) const;

      /// Return geometry parameters
      /**
       * Return geometry parameters
       */
      virtual double getRBtoroidal() const {return m_Btor_scale;}
   
      /// Returns major radius at center of toroidal cross section (a.k.a. R0)
      /**
       * Returns the major radius at the middle of the toroidal cross section (R0)
       *
       * @param[in] xi mapped grid coordinate.
       */
      Real centralMajorRadius( ) const { return m_R0; }
   
      /// Returns average minor radius
      /**
       * Returns average minor radius
       */
      Real getAvMinorRad () const { return 0.5*(m_rmax + m_rmin); }

      //Returns magnetic safety factor
      /**
       * Returns magnetic safety factor as a function of rnorm=r/R0
       *
       * @param[in] rnorm = r / R0 normalized radius.
       */
      Real getSafetyFactor( const Real rnorm ) const;

      //Returns safety factor derivative
      /**
       * Returns safety factor derivative dq(r_norm)/dr_norm
       *
       * @param[in] rnorm = r / R0 normalized radius.
       */
      Real getSafetyFactorDerivative( const Real rnorm ) const;
   
      // Object ParmParse name.
      static const std::string pp_name;

      virtual string geometryType() const {return "Toroidal";}

      double getToroidalWedgeFraction() const {return m_phi_frac;}

      //Convert (physical) Cartesian coordinates (x, y, z) of
      //a vector (vect) to the toroidal components (r, phi, theta) 
      void convertCartesianToToroidal(RealVect& vect,
				      const RealVect& coord) const;

      //convert Cartesian (x,y,z) coordinates to toroidal (r, phi, theta) 
      void convertCartesianToToroidal(RealVect& coord) const;
			
      //convert cylindrical components (R, Phi, Z) of a vector 
      //(vect) to the Cartesian (physical) components (x, y, z)
      void convertCylindricalToCartesian(RealVect& vect,
					 const RealVect& coord) const;
  
      void convertCylindricalToCartesian(FArrayBox& vect,
					 const FArrayBox& coord) const;

      /// Compute the radial derivative of the flux function
      /**
      * Public method to compute the radial derivative of the flux function.
      */
      Real dpsidr(const Real r) const;

      /// Compute the second radial derivative (d2psi/dr^2) of the flux function
      /**
       * Public method to compute the second radial derivative of the flux function.
       */
      Real d2psidr2(const Real r) const;


    private:

      /// Parse parameters.
      /**
       * Private method to obtain control parameters from 
       * "gksystem.magnetic_geometry_mapping.slab" section
       *
       * \param[in] ppcls slab geometry input database.
       */
      void parseParameters( ParmParse& pp );
   
      /// Print parameters.
      /**
       * Private method to print parameters.
       */
      void printParameters();


      /// Compute toroidal angle from phys coords.
      /**
       * Private method to compute toroidal angle.
       */
      double computePhi(const Real x,
                        const Real y ) const;


      /// Compute poloidal angle mapped coord (xi2) from phys coords.
      /**
       * Private method to compute poloidal angle.
       */
      double computeTheta(const Real z,
                          const Real r,
                          const Real rxy,
                          const Real phi) const;

      /// Compute poloidal angle corresponding to the given toroidal angle.
      /**
       * Private method to compute poloidal angle along field lines.
       */
      double computeFieldAlignedTheta(const Real r,
                                      const Real phi,
                                      const Real xi2) const;
  
      /// Compute derivatives of the poloidal angle wrt to the mapped coordinates
      RealVect getThetaDerivatives(double r,
                                   double phi,
                                   double xi2) const;

      /// Compute vector object to store magnetic flux data at face and cell centers
      void initializeMagneticFluxFunction(const Real& rmin,
					  const Real& rmax) const;

  
      double m_rmin, m_rmax, m_phi_frac;
      double m_R0, m_a;
  
      double m_Btor_scale;
      Vector<Real> m_q;
      mutable Vector<Real> m_magnetic_flux_interp;
   
      int m_block;
      int m_num_blocks;
      double m_separation;
      double m_eps;
      mutable bool m_is_flux_defined;
  
      RealVect m_mapped_block_size;

      double m_mb_center;
};

#include "NamespaceFooter.H"

#endif

#endif
