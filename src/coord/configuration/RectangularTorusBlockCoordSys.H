#ifndef _RECTANGULARTORUSBLOCKCOORDSYS_H_
#define _RECTANGULARTORUSBLOCKCOORDSYS_H_

#include "MagBlockCoordSys.H"

#include "NamespaceHeader.H"

/// Miller et al. magnetic geometry class
/**
 * Implements the analytic, magnetic equilibrium geometry described in
 *
 * R. L. Miller et al., "Noncircular, finite aspect ratio, local
 * equilibrium model", Phys. Plasmas, Vol. 5, No. 4, April 1998.
*/
class RectangularTorusBlockCoordSys
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
      RectangularTorusBlockCoordSys( ParmParse&               parm_parse,
                           const ProblemDomain&     domain);

      /// Destructor
      /**
       */
      virtual ~RectangularTorusBlockCoordSys() {;}

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

      virtual void computeFieldData( const int  dir,
                                     FArrayBox& BField,
                                     FArrayBox& BFieldMag,
                                     FArrayBox& BFieldDir,
                                     FArrayBox& gradBFieldMag,
                                     FArrayBox& curlBFieldDir,
                                     FArrayBox& BFieldDirdotcurlBFieldDir,
                                     const bool derived_data_only = false ) const;

      virtual void getCellCenteredMappedCoords(FArrayBox& xi) const;

      virtual void getFaceCenteredMappedCoords(const int dir, FArrayBox& xi) const;

      //virtual double majorRadiusConstMinorRadAdj(const RealVect a_xi) const;

      /// Return geometry parameters
      /**
       * Return geometry parameters
       */
      double getR0() const {return  m_origin[0];}
      double getBtorScale() const {return  m_Btor_scale;}

      void getRBvertical( const int dir, FArrayBox& RBvert) const;

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

      int getInfMajorrad() const { return m_inf_majorrad; }
//      virtual static const string& getParmParseName() {return pp_name;}

      // Object ParmParse name.
//      static char* pp_name;
      static const std::string pp_name;

      virtual string geometryType() const {return "RectangularTorus";}

      double getLowerZLabel() const {return m_lower_z_boundary;}
      double getUpperZLabel() const {return m_upper_z_boundary;}

#if CFG_DIM==3
      double getMinTorAngle() const {return m_min_phi_boundary;}
      double getMaxTorAngle() const {return m_max_phi_boundary;}
#endif

   protected:

      double m_Btor_scale;
      double m_bz_inner;
      double m_bz_outer;
      double m_average_radius;
      Vector<double> m_origin;
      double m_lower_z_boundary;
      double m_upper_z_boundary;
      double m_min_phi_boundary;
      double m_max_phi_boundary;
      int m_inf_majorrad;
};

#include "NamespaceFooter.H"

#endif
