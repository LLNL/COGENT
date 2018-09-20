#include "REAL.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "MultiBlockLevelGeom.H"
#include "BCInterface.H"

#include "NamespaceHeader.H"

/**
 * Boundary condition interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set ghost cell boundary conditions on solution LevelData.
*/
class InOutFlowBC : BCInterface
{
   public:

      InOutFlowBC();

      /// Destructor.
      /**
       */
      virtual ~InOutFlowBC() {;}

      /// Set inflow data
      /**
       *  Virtual function for setting boundary data.  Not used; empty
       *  implementation.
       *
       *  @param[in,out] soln a prototype of the data to which the bc will be applied
       *  @param[in] geometry the geometry.
       *  @param[in] dx the computational mesh size.
       *  @param[in] time the time at which the initial conditions are to
       * be set
       */
      virtual void setBoundaryData( const LevelData<FArrayBox>& soln,
                                    const MultiBlockLevelGeom&  geometry,
                                    const Real                  dx,
                                    const Real                  time ) {;}

      /// Set ghost cell boundary conditions.
      /**
       *  Pure virtual function that will be resposible for setting ghost
       *  cell boundary conditions on the soln data.
       *
       *  @param[in,out] soln the data for which ghost cells will be set.
       *  @param[in] geometry the geometry.
       *  @param[in] dx the computational mesh size.
       *  @param[in] time the time at which the boundary conditions are to
       * be set
       */
      virtual void ghostCellBC( LevelData<FArrayBox>&      soln,
                                const MultiBlockLevelGeom& geometry,
                                const Real                 dx,
                                const Real                 time );


      // set flux boundary conditions
      /**
       * virtual function which sets flux boundary conditions at the boundary
       * This can be useful when a boundary condition is simple to prescribe
       * on the fluxes, but difficult to express in terms of ghost cells.
       *
       *  @param[in,out] flux the face-centered data on which the BC will be set.
       *  @param[in] velocity the velocity at the boundary
       *  @param[in] the cell-centered solution
       *  @param[in] geometry the geometry.
       *  @param[in] dx the computational mesh size.
       *  @param[in] time the time at which the boundary conditions are to
       * be set
       */
      virtual void fluxBC( LevelData<FluxBox>&          flux,
                           const LevelData<FluxBox>&    velocity,
                           const LevelData<FArrayBox>&  soln,
                           const MultiBlockLevelGeom&   geometry,
                           const Real                   dx,
                           const Real                   time ) {;}

      void setFaceVelocity( LevelData<FluxBox>* a_face_vel )
      {
         CH_assert(a_face_vel!=NULL);
         m_face_vel = a_face_vel;
      }

      void clearFaceVelocity() { m_face_vel = NULL; }

   protected:

      Box getGhostBox( const Box& a_domain_box,
                       const Box& a_box,
                       const int& a_dir,
                       const Side::LoHiSide& a_side,
                       const int& a_num_ghosts ) const;

      void computeNormalVelocity( FArrayBox&      a_normalFaceVel,
                                  const FluxBox&  a_faceVel,
                                  const RealVect& a_faceArea,
                                  const int&      a_dir);

      void setRadialBoundaryValues( LevelData<FArrayBox>&     a_soln,
                                    const ProblemDomain&      a_domain,
                                    const LevelData<FluxBox>& a_face_vel,
                                    const RealVect&           a_face_area );

      LevelData<FluxBox>* m_face_vel;
};

#include "NamespaceFooter.H"
