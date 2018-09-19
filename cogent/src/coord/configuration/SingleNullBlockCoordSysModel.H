#ifndef _SINGLENULLBLOCKCOORDSYSMODEL_H_
#define _SINGLENULLBLOCKCOORDSYSMODEL_H_

#include "MagBlockCoordSys.H"
#include "Interp.H"

#include "NamespaceHeader.H"

class SingleNullBlockCoordSysModel
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
      SingleNullBlockCoordSysModel( ParmParse&           parm_parse,
                               const ProblemDomain& domain,
                               const RealVect&      dx,
                               const int            block_type);

      /// Destructor
      /**
       */
      virtual ~SingleNullBlockCoordSysModel();

      void init( ParmParse& pp );

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
   
   
      virtual Vector<Real> computeBField(const RealVect& Xi) const;

      virtual void computeFieldData( const int  dir,
                                     FArrayBox& BField,
                                     FArrayBox& BFieldMag,
                                     FArrayBox& BFieldDir,
                                     FArrayBox& gradBFieldMag,
                                     FArrayBox& curlBFieldDir,
                                     FArrayBox& BFieldDirdotcurlBFieldDir,
                                     const bool derived_data_only = false ) const;

      //Computes magnetic flux Psi
      void getMagneticFlux( const FArrayBox& a_physical_coordinates,
                            FArrayBox&       a_magnetic_flux ) const;

      double getMagneticFlux( const RealVect& a_physical_coordinate ) const;
   
      virtual void getNodalFieldData(FArrayBox& points, FArrayBox& A, FArrayBox& b, FArrayBox& Bmag) const;

      //Computes normalized magnetic flux
      void getNormMagneticFlux( const FArrayBox& a_physical_coordinates,
                                FArrayBox&       a_magnetic_flux ) const;
   
      double getNormMagneticFlux( const RealVect& a_physical_coordinate ) const;
   
      //Returns the physical radial coordinate (R) for a given value of normalized psi
      double getOuterMidplaneCoord(const double psiNorm) const;
   
      //Returns the radial coordinate of the separatrix at the outer midplane
      double getOuterRsep() const;
   
      virtual double getRBtoroidal() const {return m_RB_toroidal;}

      virtual const string getParmParseName() const {return pp_name;}

      int blockType() const {return m_block_type;}

      // Object ParmParse name.
      static const std::string pp_name;

      void getCellCenterRealCoords();

      bool truncated() const {return m_poloidally_truncated;}

      bool isValid(const RealVect& xi) const;

      virtual string geometryType() const {return "SingleNull";}
   
      RealVect getMagAxis() const {return m_magAxis;}
      RealVect getXpoint() const {return m_Xpoint;}

      RealVect realCoordPointwise( const RealVect& a_xi ) const;
      double getCoreWidth() const {return m_core_width;}
      double getSolWidth() const {return m_sol_width;}
      double getPFWidth() const {return m_pf_width;}
      double getPFArcLength () const {return m_div_leg_length;}
      double getCoreArcLength () const {return m_pol_ref_length;}

   protected:

      int m_block_type;

   private:

      void computeFieldDataSpectral( const int  dir,
                                     FArrayBox& BField,
                                     FArrayBox& BFieldMag,
                                     FArrayBox& BFieldDir,
                                     FArrayBox& gradBFieldMag,
                                     FArrayBox& curlBFieldDir,
                                     FArrayBox& BFieldDirdotcurlBFieldDir,
                                     const bool derived_data_only = false ) const;

      void computeFieldDataFromMappingFile( const int  dir,
                                            FArrayBox& BField,
                                            FArrayBox& BFieldMag,
                                            FArrayBox& BFieldDir,
                                            FArrayBox& gradBFieldMag,
                                            FArrayBox& curlBFieldDir,
                                            FArrayBox& BFieldDirdotcurlBFieldDir,
                                            const bool derived_data_only = false ) const;

      void getMagneticFluxFromDCT( const FArrayBox& a_physical_coordinates,
                                   FArrayBox&       a_magnetic_flux ) const;

      double getMagneticFluxFromDCT( const RealVect& a_physical_coordinate ) const;

  
      //In case we are doing a problem involving x-point vicinity only
      bool m_poloidally_truncated;

      double m_RB_toroidal;

      IntVect m_mapping_block_lo;
      IntVect m_mapping_block_hi;
      IntVect m_mapping_block_size;

      Interp* m_RZ_interp;

      FArrayBox m_psi_coefs;
      double m_Rmin, m_Rmax, m_Zmin, m_Zmax;   
      FArrayBox m_RBR_coefs;
      FArrayBox m_RBZ_coefs;
      Interp* m_field_interp;
      int m_psi_expansion_order;

      bool m_spectral_field;
   
      RealVect m_Xpoint;
      RealVect m_magAxis;


      //Returns (shifted) physical coordinates for the intersection of the gridlines given by psi(R,Z)=psiVal, theta(R,Z)=thetaVal
      //NB: the function returns solution only in the RCORE and RCSOL blocks, symmetry has to used to get solution in other blocks
      RealVect gridLinesIntersection( const double psiVal, const double thetaVal, const int side ) const;
   
      //Returns (shifted) Z coordinate at the lower CORE block boundary corresponding to psiVal
      double getZlo(double psiVal ) const;
   
      //Returns (shifted) X coordinate at the RCSOL- RSOL block boundary corresponding to psiVal
      double getXlo(double psiVal ) const;
   
      //Returns (shifted) Z coordinate at the upper CORE block boundary corresponding to psiVal
      double getZhi(double psiVal ) const;
   
      //Returns (shifted) physical coordinates on the reference flux surface (m_pol_ref_psi) correpsponding to the input value of the arc length.
      //"ArcLength" is measured by the polar angle going in the counterclockwise direction with the origin at the magnetic axis and zero at the X-point.
      RealVect arcLengthInverse(const double length) const;
   
      //Magnetic-flux countour-function (defined only for Z>m_Zx + m_Z0)
      double psi(const RealVect& a_x) const;
   
      //Poloidal coordinate countour-function (defined only for Z>m_Zx + m_Z0)
      double theta(const RealVect& a_x) const;
   
      //Returns value of psi at R=m_R0 (defined only for Z>m_Zx + m_Z0)
      double psiAtR0(const double a_Z) const;
   
      //Computes normalization, wuch that the input BpScale corresponds to Bpol (T) at the separatrix top.
      void normalizeBpScale();
   
      //Return partial Z and X derivatives f psi and theta (defined only for Z>m_Zx + m_Z0)
      double psiZ(const RealVect& a_x) const;
      double psiX(const RealVect& a_x) const;
      double thetaZ(const RealVect& a_x) const;
      double thetaX(const RealVect& a_x) const;
   
      double m_BpScale; //in the input that corresponds to Bpol(T) at the separatrix top (later it is internally normalized by normalizeBpScale())
   
      // Geometry parameters
      double m_a;       //shape coeffiecient of a flux surface
      double m_b;       //shape coeffiecient of a flux surface
      double m_c;       //shape coeffiecient of a flux surface
   
      double m_Zx;      //Z-location of the X-point (unshifted)
      double m_Zc;      //Z-location of the O-point (unshifted)
      double m_Zsep_hi; //Z-location of the separatrix top (unshifted)
   
      double m_R0;      //radial shift
      double m_Z0;      //vertical shift
   
      double m_core_width; //vertical size of the core block at the top of a tokamak (i.e., CORE top block boundary)
      double m_sol_width;  //vertical size of the sol block at the top of a tokamak
      double m_pf_width;   //vertical size of the pf block when projected at the top of a tokamak (we have up-down symmetry w/r/t the X-point)
      double m_div_leg_length; //poloidal length of divertor legs (measured by the poloidal angle; origin at magnetic axis --makes use of the vertical symmetry about Zx)
   
      double m_pol_ref_z_lo;      //z-coordinate at the CORE bottom block boundary that designate the flux surface, on which the poloidal coordinate is define as the arc length)
      double m_pol_ref_length;    //poloidal length of the reference flux surface from bottom to top (half of the full closed-surface length) -- here, measured as PI
      double m_pol_ref_length_mp; //poloidal length of the reference flux surface from bottom to midplane -- here, measured as PI/2
      double m_pol_ref_psi;       //value of psi at the reference flux surface
   
      double m_tol;
      int m_max_iterations;

};

#include "NamespaceFooter.H"

#endif
