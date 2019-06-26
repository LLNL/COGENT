#ifndef _SLABBLOCKCOORDSYS_H_
#define _SLABBLOCKCOORDSYS_H_

#include "MagBlockCoordSys.H"

#include "NamespaceHeader.H"

/// Slab magnetic geometry class
/**
 * Implements an analytic slab model, rectangular mesh, field in y and z directions
 * with optional gradient in x
 *
 * 2D nomenclature: indices are RADIAL_DIR (x), POLOIDAL_DIR (z)
 * 3D nomenclature: indices are RADIAL_DIR (x), TOROIDAL_DIR (y), POLOIDAL_DIR (z)
 * 
 * If field-algined option is true the physical coordinates are given in 3D by
 *
 * x[0] = xi[0] * m_phys_block_size[0]/m_mapped_block_size[0]
 * x[1] = (xi[1] - separation_shift) * m_phys_block_size[1]/m_mapped_block_size[1]
 * x[2] = xi[2] * m_phys_block_size[2]/m_mapped_block_size[2] + 
 *      + Bz/By * (xi[1] - xi[1]_mid) * m_phys_block_size[1]/m_mapped_block_size[1]
 *
 * Bz/By = Bz/By[x0] = Bz/By [xi[0] * m_phys_block_size[0]/m_mapped_block_size[0]],
 * 
 * separation_shift = block_separation * dx[mb_dir] * block
 * 
 * 0, 1, and 2 directions correspond to x,y,and z
 * 
 * In other words the mapped coordinates xi[0] and xi[2] designate the location
 * of the field-line at the mid-plane (in xi[1] direction), and xi[1] is a 
 * field-line coordinate measure by the "toroidal distance"
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
       * @param[in] block block id.
       * @param[in] num_blocks total number or blocks.
       * @param[in] number of separating cells (between blocks).
       */
   
      SlabBlockCoordSys( ParmParse&               parm_parse,
                         const ProblemDomain&     domain,
                         const int                block,
                         const int                num_blocks,
                         const int                block_separation);


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

      virtual void getNodalFieldData(FArrayBox& points, FArrayBox& A, FArrayBox& b, FArrayBox& Bmag) const;

      /// Return geometry parameters
      /**
       * Return geometry parameters
       */
      //Need to include this function, as required by MagBlockCoordSys
      double getRBtoroidal() const {return 0.0;}

      //Returns magnetic field pitch
      double getPitch( const RealVect& physical_coordinate ) const;

      //Returns magnetic field pitch "shear" in mapped coordinates
      double getPitchShear( const RealVect& physical_coordinate ) const;

      // Object ParmParse name.
      static const std::string pp_name;

      virtual string geometryType() const {return "Slab";}

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

      double m_xmax;
      double m_zmax;
#if CFG_DIM==3
      double m_ymax;
#endif
   
      double m_ByInner;
      double m_ByOuter;
      double m_BzInner;
      double m_BzOuter;

      int m_block;
      int m_num_blocks;
      int m_mb_dir;
      double m_separation;
   
      RealVect m_phys_block_size;
      RealVect m_mapped_block_size;

      double m_mb_center;
};

#include "NamespaceFooter.H"

#endif
