#ifndef _POISSON_H_
#define _POISSON_H_

#include "FieldSolver.H"
#ifdef USE_PETSC
#include "MBPETScSolver.H"
#else
#include "MBHypreSolver.H"
#endif

#include "NamespaceHeader.H"

/// Poisson solver class.
/**
 * This operator class solves the mapped Poisson equation:
 *
 * \f[
 *  {\nabla} \cdot \left ( N^T N/J  {\nabla} \Phi \right ) = \rho,
 * \f]
*/
class Poisson
   : public FieldSolver
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] geom the geometry
    */
   Poisson( const ParmParse& pp,
            const MagGeom&   geom );
      
   /// Destructor.
   /**
    */
   virtual ~Poisson() {;}

   /// Construct the operator coefficients
   /**
    * Computes the operator coefficients.
    */
   void setOperatorCoefficients( const PotentialBC& bc );

   void computeCoefficients( LevelData<FluxBox>& coefficients );

   virtual void solvePreconditioner( const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out );

   // Multiply by the unmapped coefficients
   virtual void multiplyUnmappedCoefficients( LevelData<FluxBox>& a_data ) const {;}

   /// ParmParse object name
   /**
    * ParmParse object name.  Used to locate input file parameters of interest to
    * this class.
    */
   static char* pp_name;

protected:

#ifdef USE_PETSC
   MBPETScSolver m_solver;
#else
   MBHypreSolver m_solver;
#endif

   LevelData<FluxBox> m_mapped_coefficients;
};



#include "NamespaceFooter.H"

#endif
