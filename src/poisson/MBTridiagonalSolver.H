#ifndef _MBTRIDIAGONALSOLVER_H_
#define _MBTRIDIAGONALSOLVER_H_

#include "MBSolver.H"
#include "FluxSurface.H"

#include "NamespaceHeader.H"



class MBTridiagonalSolver
   : public MBSolver
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] geom the geometry
    *  @param[in] discretization_ordr discretization order
    */
   MBTridiagonalSolver( const MultiBlockLevelGeom&   geom,
                        const int                    discretization_order);
      
   /// Destructor.
   /**
    */
   virtual ~MBTridiagonalSolver() {};

   virtual void constructMatrixGeneral( LevelData<FArrayBox>& alpha_coefficient,
                                        LevelData<FluxBox>&   tensor_coefficient,
                                        LevelData<FArrayBox>& beta_coefficient,
                                        const PotentialBC&    bc );

   virtual void multiplyMatrix( const LevelData<FArrayBox>& in,
                                LevelData<FArrayBox>&       out ) const;

   virtual void solve( const LevelData<FArrayBox>&  rhs,
                       LevelData<FArrayBox>&        solution,
                       bool                         homogeneous_bcs );

   virtual void dumpMatrix(const string& file_name) const;

protected:

   void constructTridiagonalMatrix( LevelData<FArrayBox>&    alpha_coefficient, 
                                    LevelData<FluxBox>&      tensor_coefficient,
                                    LevelData<FArrayBox>&    beta_coefficient,
                                    const PotentialBC&       bc,
                                    FArrayBox&               stencil_values,
                                    LevelData<FArrayBox>&    radial,
                                    const int                diagonal_offset,
                                    const bool               fourthOrder,
                                    const bool               dropOrder,
                                    LevelData<FArrayBox>&    rhs_from_bc ) const;

   void solveFluxSurfaceAverage( LevelData<FArrayBox>& data ) const;

   void solveTridiagonal( const bool     periodic,
                          const int      n,
                          double *       diag,
                          double *       udiag,
                          double *       ldiag,
                          double *       fs_average ) const;

   void solveTridiagonalNeumann( const int      n,
                                 double *       diag,
                                 double *       udiag,
                                 double *       ldiag,
                                 double *       fs_average,
                                 const double * weights ) const;

   bool isCoreRadialPeriodicOrNeumannBC( const PotentialBC& bc ) const;

   bool m_periodic_or_neumann;

   FluxSurface m_flux_surface;

   LevelData<FArrayBox> m_A_radial;

   FArrayBox m_A_stencil_values;

   int m_A_diagonal_offset;
};



#include "NamespaceFooter.H"

#endif
