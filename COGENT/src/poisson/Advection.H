#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include "EllipticOp.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "MBHypreSolver.H"

#include "MagFluxAlignedMBLEXCenter.H"

#include "GridFunction.H"
#include "GridFunctionLibrary.H"

#include "NamespaceHeader.H"

/// Advection operator class.
/**
 * This operator class solves the mapped gyrokinetic Poisson equation:
 *
 * \f[
 *  {\nabla} \cdot \left ( N^T D N/J  {\nabla} \Phi \right ) = \rho,
 * \f]
*/
class Advection
   : public EllipticOp
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] geom the geometry data.
    */
   Advection( const ParmParse&  pp,
              const MagGeom&    geom );
      
   /// Destructor.
   /**
    */
   virtual ~Advection();

   /// Construct the operator coefficients
   /**
    * Computes the operator coefficients.
    *  @param[in] pointwise face-centered unmapped_coefficients.
    *  @param[in] pointwise face-centered mapped_coefficients.
    *  @param[in] boundary conditions
    *
    */

   virtual void setOperatorCoefficients( const LevelData<FluxBox>&  velocity_phys);

   virtual void updatePreconditioner(const EllipticOpBC&    bc );

   virtual void solvePreconditioner( const LevelData<FArrayBox>&  in,
                                     LevelData<FArrayBox>&        out );

   virtual void applyOp(LevelData<FArrayBox>& out,
                        const LevelData<FArrayBox>& in,
                        bool homogeneous = false);

   virtual void computeBcDivergence( LevelData<FArrayBox>& out );
   
   virtual void multiplyCoefficients(LevelData<FluxBox>& data,
                                     const bool mapped_coeff ) const {};

   void applyPCOp(const LevelData<FArrayBox>& in_vec,
                  LevelData<FArrayBox>&       matvec);

   const string getPrecondMethod() const {return m_preconditioner->getMethod();}
   const Real   getPrecondTol() const {return m_preconditioner->getTol();}
   const Real   getPrecondMaxIter() const {return m_preconditioner->getMaxIter();}
   const bool   getPrecondVerbose() const {return m_preconditioner->getVerbose();}
   
   
   /// parseParameters function
   /**
    * parseParameters function
    */
   void parseParameters( const ParmParse& pp );

   /// printParameters function
   /**
    * printParameters function
    */
   void printParameters();

   /// ParmParse object name
   /**
    * ParmParse object name.  Used to locate input file parameters of interest to
    * this class.
    */
   static const char* pp_name;
   
protected:

   MBHypreSolver* m_preconditioner;

   LevelData<FluxBox> m_velocity;
   LevelData<FArrayBox> m_volume_reciprocal;
   
   MultiBlockLevelExchangeCenter* m_mblx_ptr;

   bool m_verbosity;
   bool m_precond_verbose;
   bool m_precond_precond_verbose;
   
   Real m_linear_shift;
   Real m_diffusion_coeff;
};



#include "NamespaceFooter.H"

#endif
