#ifndef _MILLERBLOCKCOORDSYS_H_
#define _MILLERBLOCKCOORDSYS_H_

#include "MagBlockCoordSys.H"

#include "NamespaceHeader.H"

/// Miller et al. magnetic geometry class
/**
 * Implements the analytic, magnetic equilibrium geometry described in
 *
 * R. L. Miller et al., "Noncircular, finite aspect ratio, local
 * equilibrium model", Phys. Plasmas, Vol. 5, No. 4, April 1998.
*/
class MillerBlockCoordSys
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
       * @param[in] ghostVect vector of ghost cells stored for mapping.
       */
      MillerBlockCoordSys( ParmParse&               parm_parse,
                           const ProblemDomain&     domain);

      /// Destructor
      /**
       */
      virtual ~MillerBlockCoordSys() {;}

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

      Real dX2dXi2( const RealVect& Xi, int dirX, int dirXi0, int dirXi1 ) const;

      virtual array<double,3> computeBField(const RealVect& a_X) const;
   
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

      virtual void getNodalFieldData(FArrayBox& points,
                                     FArrayBox& A,
                                     FArrayBox& b,
                                     FArrayBox& Bmag) const;

      virtual void getNodeCenteredMappedCoords(FArrayBox& a_xi) const;

      virtual void getCellCenteredMappedCoords(FArrayBox& xi) const;

      virtual void getFaceCenteredMappedCoords(const int dir, FArrayBox& xi) const;

      virtual double majorRadiusConstMinorRadAdj(const RealVect a_xi) const;

      /// Return geometry parameters
      /**
       * Return geometry parameters
       */

      /**
       * Returns the value of the innermost flux label
       */

      double getInnerFluxLabel() const {return m_inner_radial_boundary;}
      double getOuterFluxLabel() const {return m_outer_radial_boundary;}
   
      double getBeta() const {return  m_beta;}
      double getKappa() const {return  m_kappa;}
      double getDpsidr() const {return  m_dpsidr;}
      double getR0() const {return  m_origin[0];}
      double getDrR0() const {return  m_drR0;}
      double getSKappa() const {return  m_s_kappa;}
      double getSDelta() const {return  m_s_delta;}
      double getBtorScale() const {return  m_Btor_scale;}

      void getRBpoloidal( const int dir, FArrayBox& RBpol, FArrayBox& dRBpoldt ) const;

      virtual double getRBtoroidal() const {return m_Btor_scale;}

      /// Returns major radius at center of toroidal cross section (a.k.a. R0)
      /**
       * Returns the major radius at the middle of the toroidal cross section (R0)
       *
       * @param[in] xi mapped grid coordinate.
       */
      Real centralMajorRadius( ) const { return m_origin[0]; }

      /// Returns average minor radius
      /**
       * Returns average minor radius
       */
      Real getAvMinorRad () const { return m_average_radius; }

      /**
       * Returns an integer flag that tells whether a constant-minor-radius approximation
       * is used
       */
      int getConstMinorrad() const { return m_const_minorrad; }

      // Object ParmParse name.
      static const std::string pp_name;

      virtual bool isFieldAlignedMapping() const {return false;}

      virtual string geometryType() const {return "Miller";}

   protected:

      double m_inner_radial_boundary;
      double m_outer_radial_boundary;
      double m_beta;
      double m_kappa;
      double m_dpsidr;
      double m_drR0;
      double m_s_kappa;
      double m_s_delta;
      double m_Btor_scale;
      double m_average_radius;
      int m_const_minorrad;
      Vector<double> m_origin;
};

#include "NamespaceFooter.H"

#endif
