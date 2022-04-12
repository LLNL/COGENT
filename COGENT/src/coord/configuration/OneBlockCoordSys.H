#ifndef _ONEBLOCKCOORDSYS_H_
#define _ONEBLOCKCOORDSYS_H_

#include "MagBlockCoordSys.H"
#include "RadialUtils.H"

#include "Interp.H"

#include "NamespaceHeader.H"

/// New arbitrary one block geometry class
/**
 * Implements a(n?) one block geometry model via tables
 *
 * Nomenclature: indices are still RADIAL_DIR <--> R, POLOIDAL_DIR <--> Z, TOROIDAL_DIR <--> Y
 * to maintain contact with toroidal simulations
*/
class OneBlockCoordSys
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
      OneBlockCoordSys( ParmParse&      parm_parse,
                         const ProblemDomain&  domain );

      /// Destructor
      /**
       */
      virtual ~OneBlockCoordSys();


      void init( ParmParse& pp );
      void readFiles( ParmParse& pp );
      virtual RealVect mappedCoordNewton( const RealVect& x, const RealVect& xi_initial, const IntVect& iv_initial ) const;
      void getCellCenterRealCoords();
      bool isValid(const RealVect& xi) const;

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
      virtual void dXdXi( FArrayBox&  dxdXi,
                    const FArrayBox&  Xi,
                          int         destComp,
                          int         dirX,
                          int         dirXi,
                    const Box&        box ) const;
   
      virtual array<double,3> computeBField(const RealVect& a_X) const;

      virtual void computeFieldData( const int  dir,
                                     FArrayBox& BField,
                                     FArrayBox& BFieldMag,
                                     FArrayBox& BFieldDir,
                                     FArrayBox& gradBFieldMag,
                                     FArrayBox& curlBFieldDir,
                                     FArrayBox& BFieldDirdotcurlBFieldDir,
                                     const bool derived_data_only = false ) const;
									 
  /// Pointer to class methods
  void (OneBlockCoordSys::*calcB)(double R, double Z, double* p_B) const;
		
  /// Magnetic field types: Const, Toroidal, Bennett, BRData (custom profile for B(r))
  void calcBConst(double, double, double*) const;

  virtual void getMagneticFlux( const FArrayBox& physical_coordinates,
                                    FArrayBox&       magnetic_flux ) const;

  virtual double getMagneticFlux( const RealVect& a_physical_coordinate ) const;

  virtual void getNodalFieldData(FArrayBox& points, FArrayBox& A, FArrayBox& b, FArrayBox& Bmag) const;

  // Object ParmParse name.
  static const std::string pp_name;
  
  virtual bool isFieldAlignedMapping() const {return false;}
  
  virtual string geometryType() const {return "OneBlock";}
  
  double getRBtoroidal() const;

  double rmin() const {return m_rmin;}
  double rmax() const {return m_rmax;}

   protected:


   double m_rmin, m_rmax;
   double m_thetamin, m_thetamax;
#if CFG_DIM==3
   double m_phimin, m_phimax;
#endif
   double m_Btor_scale;	               /// Magnetic constant for toroidal model: B_0*r_0
   double m_Btor_0;	               /// Magnetic field amplitude for constant and bennett models
   double m_r0;	                       /// Radial size of Bennett model
   HermiteSpline m_B_spline;   /// Hermite spline representation for magnetic field B(r)

   FArrayBox m_realCoords;
   int m_rc_coarsen_ratio;

   IntVect m_mapping_block_lo;
   IntVect m_mapping_block_hi;
   IntVect m_mapping_block_size;

   Interp* m_RZ_interp;
	  
};

#include "NamespaceFooter.H"

#endif
