#ifndef _SNCOREBLOCKCOORDSYS_H_
#define _SNCOREBLOCKCOORDSYS_H_

enum SNCoreBlockType {L_CORE,R_CORE,M_CORE,NUM_SNCORE_BLOCKS};

#undef PLOT_INVERSE_POINTS

#include "MagBlockCoordSys.H"
#include "Interp.H"

#include "NamespaceHeader.H"

class SNCoreBlockCoordSys
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
      SNCoreBlockCoordSys( ParmParse&           parm_parse,
                           const ProblemDomain& domain,
                           const RealVect&      dx,
                           const int            block_type );

      /// Destructor
      /**
       */
      virtual ~SNCoreBlockCoordSys();

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

      /// Returns computational coordinate computed by a Newton iteration
      /**
       * Returns location \f$\xi\f$ in mapped coordinate space of real-world
       * location \f$x\f$.
       *
       * @param[in] x physical space coordinate.
       * @param[in] xi_initial  initial guess for Newton iteration
       */
      virtual RealVect mappedCoordNewton( const RealVect& x, const RealVect& xi_initial, const IntVect& iv_initial ) const;
   
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

     virtual Vector<Real> computeBField(const RealVect& a_X) const;
   
     virtual void computeFieldData( const int  dir,
                                     FArrayBox& BField,
                                     FArrayBox& BFieldMag,
                                     FArrayBox& BFieldDir,
                                     FArrayBox& gradBFieldMag,
                                     FArrayBox& curlBFieldDir,
                                     FArrayBox& BFieldDirdotcurlBFieldDir,
                                     const bool derived_data_only = false ) const;

      //Computes magnetic flux Psi
      virtual void getMagneticFlux( const FArrayBox& physical_coordinates,
                                    FArrayBox&       magnetic_flux ) const;

      virtual double getMagneticFlux( const RealVect& a_physical_coordinate ) const;
   
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

      bool isValid(const RealVect& xi) const;

      virtual string geometryType() const {return "SNCore";}
   
      RealVect getMagAxis() const {return m_magAxis;}

      RealVect getXpoint() const {return m_Xpoint;}
   
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

      void getMagneticFluxFromDCT( const FArrayBox& physical_coordinates,
                                   FArrayBox&       magnetic_flux ) const;

      double getMagneticFluxFromDCT( const RealVect& physical_coordinate ) const;

      FArrayBox m_realCoords;
      int m_rc_coarsen_ratio;

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

#ifdef PLOT_INVERSE_POINTS
      FILE *m_ipt_file[NUM_SNCORE_BLOCKS];
#endif
   
};

#include "NamespaceFooter.H"

#endif
