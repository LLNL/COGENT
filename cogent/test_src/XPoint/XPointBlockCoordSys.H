#ifndef _XPOINTBLOCKCOORDSYS_H_
#define _XPOINTBLOCKCOORDSYS_H_

#undef PLOT_INVERSE_POINTS

#include "MagBlockCoordSys.H"

#include "NamespaceHeader.H"

class XPointBlockCoordSys
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
      XPointBlockCoordSys( ParmParse&           parm_parse,
                               const ProblemDomain&     domain,
                               const int                block_type);

      /// Destructor
      /**
       */
      virtual ~XPointBlockCoordSys();

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

      virtual Real dX2dXi2( const RealVect& Xi, int dirX, int dirXi0, int dirXi1 ) const;

      virtual void getRBpoloidal( const int dir, FArrayBox& RBpol, FArrayBox& dRBpoldt ) const;

      virtual double getRBtoroidal() const {return 1.;}

      virtual const string getParmParseName() const {return pp_name;}

      int blockType() const {return m_block_type;}

      IntVect domainSmallEnd() const {return m_domain.domainBox().smallEnd();}

      Real centralMajorRadius( ) const { return m_origin[0]; }

      virtual void computeFieldData( const int  dir,
                                     FArrayBox& BField,
                                     FArrayBox& BFieldMag,
                                     FArrayBox& BFieldDirl,
                                     FArrayBox& gradBFieldMag,
                                     FArrayBox& curlBFieldDir,
                                     FArrayBox& BFieldDirdotcurlBFieldDir ) const {;}

       // Object ParmParse name.
      static const char* pp_name;

   protected:

      int m_block_type;
      Vector<double> m_origin;

      double m_alpha;
#ifdef PLOT_INVERSE_POINTS
      FILE *m_ipt_file[8];
#endif
};

#include "NamespaceFooter.H"

#endif
