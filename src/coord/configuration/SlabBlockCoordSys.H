#ifndef _SLABBLOCKCOORDSYS_H_
#define _SLABBLOCKCOORDSYS_H_

#include "MagBlockCoordSys.H"

#include "NamespaceHeader.H"

/// Slab et al. magnetic geometry class
/**
 * Implements an analytic slab model, rectangular mesh, field in Z direction
 * with optional gradient in x
 *
 * Nomenclature: indices are still RADIAL_DIR (x), POLOIDAL_DIR (y), TOROIDAL_DIR (z)
 * to maintain contact with toroidal simulations
*/
class SlabBlockCoordSys
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
      SlabBlockCoordSys( ParmParse&               parm_parse,
                           const ProblemDomain&     domain);

      /// Destructor
      /**
       */
      virtual ~SlabBlockCoordSys() {;}

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

      virtual void getCellCenteredMappedCoords(FArrayBox& xi) const;

      virtual void getFaceCenteredMappedCoords(const int dir, FArrayBox& xi) const;

      // virtual double majorRadiusConstMinorRadAdj(const RealVect a_xi) const;

      /// Return geometry parameters
      /**
       * Return geometry parameters
       */

      double getBzInner( ) const {return m_BzInner;}

      double getBzOuter( ) const {return m_BzOuter;}

     /// Returns By (poloidal) field at inner value of x.  ByOuter then set.
     /**
      * Returns By field at inner value of x.  ByOuter then set by BzInner and BzOuter
      */
      double getByInner( ) const {return m_ByInner;}

      double getRBtoroidal() const {return 0.0;}

//      virtual static const string& getParmParseName() {return pp_name;}

      // Object ParmParse name.
//      static char* pp_name;
      static const std::string pp_name;

      virtual string geometryType() const {return "Slab";}

   protected:

      double m_ximax;
      double m_xmax;
      double m_ymax;
#if CFG_DIM==3
      double m_zmax;
#endif
      RealVect m_xmaxvec; 
      double m_BzInner;
      double m_BzOuter;
      double m_ByInner;
      double m_twopi;

  /*      double m_beta;
      double m_kappa;
      double m_dpsidr;
      double m_drR0;
      double m_s_kappa;
      double m_s_delta;
      double m_Btor_scale;
      double m_average_radius;
      int m_const_minorrad;
      Vector<double> m_origin;
  */
};

#include "NamespaceFooter.H"

#endif
