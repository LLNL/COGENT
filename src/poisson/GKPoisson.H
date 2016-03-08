#ifndef _GKPOISSON_H_
#define _GKPOISSON_H_

//#if CFG_DIM==2
//#define PS_NAMESPACE D4
//#define PDIM 4
//#endif
//#if CFG_DIM==3
//#define PS_NAMESPACE D5
//#define PDIM 5
//#endif

#include "FieldSolver.H"
#include "MultiBlockLevelExchangeAverage.H"
#ifdef with_petsc
#include "MBPETScSolver.H"
#else
#include "MBHypreSolver.H"
#endif

#include "NamespaceHeader.H"
//namespace PS = PS_NAMESPACE;

/// Gyrokintic Poisson operator class.
/**
 * This operator class solves the mapped gyrokinetic Poisson equation:
 *
 * \f[
 *  {\nabla} \cdot \left ( N^T D N/J  {\nabla} \Phi \right ) = \rho,
 * \f]
 * where
 * \f{eqnarray*}
 * {\bf D} & \equiv & (De)^2  \mathbf{I} + \frac{(La)^2}{ B^2 }
 *  \left ( \sum_i Z_i m_i \bar{n}_i \right ) \left ( \mathbf{I} -
 *  \mathbf{b}\mathbf{b}^T \right ),\\  
 * \rho & \equiv &  n_e - \sum_i Z_i \bar{n}_i.
 * \f}
*/
class GKPoisson
   : public FieldSolver
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] geom the geometry data.
    *  @param[in] larmor_number the Larmor number for the operator.
    *  @param[in] debye_number the Debye number for the operator.
    */
   GKPoisson( const ParmParse&   pp,
              const MagGeom&     geom,
              const Real         larmor_number,
              const Real         debye_number );
      
   /// Destructor.
   /**
    */
   virtual ~GKPoisson();

   /// Construct the GKP operator coefficients
   /**
    * Computes the GKP operator coefficients.  The polarization density
    * term requires the mass-weighted sum of the ion charge densities,
    * which is obtained from the species vector
    *
    * @param[in] species   Array of kinetic species
    */

   virtual void setOperatorCoefficients( const LevelData<FArrayBox>& ni,
                                         const PotentialBC&          bc );

   void computeCoefficients( const LevelData<FArrayBox>& ni,
                             LevelData<FluxBox>& mapped_coefficients,
                             LevelData<FluxBox>& unmapped_coefficients );

   virtual void setPreconditionerConvergenceParams( const double tol,
                                                    const int    max_iter,
                                                    const double precond_tol,
                                                    const int    precond_max_iter );

   virtual void solvePreconditioner( const LevelData<FArrayBox>& in,
                                     LevelData<FArrayBox>&       out );

   virtual void multiplyUnmappedCoefficients( LevelData<FluxBox>& data ) const;

   /// ParmParse object name
   /**
    * ParmParse object name.  Used to locate input file parameters of interest to
    * this class.
    */
   static const char* pp_name;

protected:

   void getMinMax(LevelData<FArrayBox>& density, double& min, double& max) const;

   void getMinMax(LevelData<FluxBox>& density, double& min, double& max) const;

   void fillDensityGhosts(LevelData<FArrayBox>& a_density) const;

#ifdef with_petsc
   MBPETScSolver* m_preconditioner;
#else
   MBHypreSolver* m_preconditioner;
#endif

   Real m_larmor_number2;
   Real m_debye_number2;

   MultiBlockLevelExchangeAverage* m_mblex_density_Ptr; 
   int m_density_interp_ghosts;

   LevelData<FluxBox> m_mapped_coefficients;
   LevelData<FluxBox> m_unmapped_coefficients;
};



#include "NamespaceFooter.H"

#endif
