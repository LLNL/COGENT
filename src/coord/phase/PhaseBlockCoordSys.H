#ifndef _PHASE_BLOCK_COORD_SYS_H_
#define _PHASE_BLOCK_COORD_SYS_H_

#include "Dimensions.H"
#include "NewFourthOrderCoordSys.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "MagBlockCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM VEL_DIM
#include "VelCoordSys.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#include "LevelData.H.multidim"  // Why do I need to include this???
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;
namespace VEL = VEL_NAMESPACE;
namespace CP1 = CP1_NAMESPACE;

/// Phase block coordinate system class
/**
 * Coordinate system for phase space block.
*/
class PhaseBlockCoordSys : public NewFourthOrderCoordSys
{
   public:
      /// Constructor
      /**
       * Constructor
       *
       * @param[in] config_space_coord_sys     configuration space coord system.
       * @param[in] velocity_space_coord_sys   velocity space coord system.
       * @param[in] domain                     phase space domain.
       */
      PhaseBlockCoordSys( const CFG::MagBlockCoordSys& config_space_coord_sys,
                          const VEL::VelCoordSys&      velocity_space_coord_sys,
                          const ProblemDomain&         domain );

      /// Destructor
      /**
       */
      virtual ~PhaseBlockCoordSys();

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
      virtual RealVect mappedCoord(const RealVect& x) const;

      /// Returns the derivative of the spatial coordinate X with respect to the
      /// mapped coordinate Xi
      /**
       * Returns the derivative of the spatial coordinate \f$X\f$ with respect to
       * the mapped coordinate \f$\Xi\f$.
       *
       * @param[in] X physical space coordinate.
       * @param[in] dirX component of X to be differentiated.
       * @param[in] dirXi component of Xi with respect to which derivative is taken.
       */
      virtual Real dXdXi(const RealVect& X, int dirX, int dirXi) const;

      virtual void getN(FluxBox& N, const Box& box) const;

      virtual void getAvgJ(FArrayBox& data, const Box& box) const;

      virtual void pointwiseJ(FArrayBox& data, const Box& box) const;

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

      /// Project the vector argument onto the unit vector in the
      /// direction of the mapped radial coordinate
      /**
       * @param[in] data    vector to be projected
       * @param[out] data   projected vector
       */
      void projectOntoRadial( FArrayBox& data ) const;
      void projectOntoRadial( FluxBox& data ) const;

      /// Compute projections of the configurational components of vector 
      /// argument onto the psi(normal) and theta (tangental) directions
      /**
       * @param[in] data    vector to be projected
       * @param[out] data   projected vector
       */
      void computePsiThetaProjections( FArrayBox& data ) const;
      void computePsiThetaProjections( FluxBox& data ) const;

      void computeGradFFactor( FArrayBox& a_data ) const;
      void computeGradFFactor( FluxBox& a_data ) const;

      virtual void cellVol(FArrayBox& vol, const FluxBox& N, const Box& box) const;

      /// Returns configuration space mapping.
      /**
       */
      const CFG::MagBlockCoordSys& configSpaceCoordSys() const { return m_configuration_coords; }

      /// Returns \f$v_{\parallel}-\mu\f$ space mapping.
      /**
       */
      const VEL::VelCoordSys& velSpaceCoordSys() const { return m_velocity_coords; }

      /// Fills the input FArrayBox with the real (physical) velocities
      /**
       */
      void getVelocityRealCoords(FArrayBox& velocityRealCoords) const;

      /// Fills the input FArrayBox with the real coordinates at cell centers
      /**
       */
      void getCellCenteredRealCoords( FArrayBox& cell_centered_coords ) const;

      /// Fills the input FArrayBox with the mapped coordinates at cell centers
      /**
       */
      void getCellCenteredMappedCoords( FArrayBox& cell_centered_coords ) const;

      /// Fills the input FluxBox with the mapped coordinates at cell face centers
      /**
       */
      void getFaceCenteredMappedCoords( FluxBox& face_centered_coords ) const;

      /// Fills the input FArrayBox with the flux coordinates at cell centers
      /**
       */
      void getCellCenteredFluxCoords( FArrayBox& cell_centered_coords ) const;

      /// Returns species definition status
      /**
       * Returns whether or not the object has been defined for a particular species
       *
       */
   //      bool speciesDefined() { return m_speciesDefined; }

      /// Returns the volume of mapped cells
      /**
       * Returns the volume of mapped cells
       *
       */
      double getMappedCellVolume() const {return m_mapped_cell_volume;}

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

      // Transdimensional utilities:

      void injectConfigurationToPhase( const CFG::FArrayBox& src,
                                       FArrayBox&            dst ) const;
      void injectConfigurationToPhase( const CFG::FArrayBox& src,
                                       const VEL::IntVect&   vpmu,
                                       FArrayBox&            dst ) const;
      void injectConfigurationToPhase( const CFG::FluxBox& src,
                                       FluxBox&            dst ) const;
      void injectConfigurationToPhase( const CFG::FluxBox& src,
                                       const VEL::IntVect&  vpmu,
                                       FluxBox&            dst ) const;
      void projectPhaseToConfiguration( const Box& src,
                                        CFG::Box& dst ) const;

      Box configurationFlatten(const Box& disjointBox, const Box& src_box) const;

      Box velocityFlatten(const Box& disjointBox, const Box& src_box) const;

      CFG::RealVect config_restrict( const RealVect& x ) const;

      CFG::IntVect config_restrict( const IntVect& x ) const;

      VEL::RealVect vel_restrict( const RealVect& x ) const;

      VEL::IntVect vel_restrict( const IntVect& x ) const;

      const ProblemDomain& domain() const {return m_domain;}

   private:

      ProblemDomain m_domain;

      const CFG::MagBlockCoordSys& m_configuration_coords;
      const VEL::VelCoordSys&      m_velocity_coords;

      double   m_mapped_cell_volume;
      RealVect m_mapped_face_area;
};

#include "NamespaceFooter.H"

#endif
