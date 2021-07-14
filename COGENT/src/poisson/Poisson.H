#ifndef _POISSON_H_
#define _POISSON_H_

#include "EllipticOp.H"
#include "MBHypreSolver.H"

#include "MagFluxAlignedMBLEXCenter.H"

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
   : public EllipticOp
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
   virtual ~Poisson();

   /// Construct the operator coefficients
   /**
    * Computes the operator coefficients.
    */
   void setOperatorCoefficients( const EllipticOpBC& bc );

   void computeCoefficients( LevelData<FluxBox>& coefficients );

   virtual void solvePreconditioner( const LevelData<FArrayBox>& in,
                                     LevelData<FArrayBox>&       out );

   // Multiply by the unmapped coefficients
   virtual void multiplyCoefficients( LevelData<FluxBox>& data,
                                      const bool mapped_coeff) const {;}

   /// ParmParse object name
   /**
    * ParmParse object name.  Used to locate input file parameters of interest to
    * this class.
    */
   static const std::string pp_name;

protected:

   MBHypreSolver* m_preconditioner;

   LevelData<FluxBox> m_mapped_coefficients;

   MagFluxAlignedMBLEXCenter* m_mblx_ptr;

   bool m_precond_verbose;
   bool m_precond_precond_verbose;
};



#include "NamespaceFooter.H"

#endif
