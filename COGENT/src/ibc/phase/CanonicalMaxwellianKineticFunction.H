#ifndef _CANONICAL_MAXWELLIAN_KINETIC_FUNCTION_H_
#define _CANONICAL_MAXWELLIAN_KINETIC_FUNCTION_H_

#include "BoundaryBoxLayout.H.multidim"
#include "REAL.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "LevelData.H"
#include "PhaseGeom.H"
#include "ParmParse.H"
#include "RealVect.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunction.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include <string>

#include "NamespaceHeader.H"

/**
 * Canonical Maxwellian distribution function initial condition class.
 * Initializes canonical Maxwellian distribution as a function of
 * the energy invariant: 
 *
 * eps = (m v_par^2 + mu * B)
 *
 * and the "toroidlal invariant":
 * 
 * psi_inv = psi + (larmor_number * mass/charge)*RBtor/B*vpar 
 *
 * The canonical Maxwellian funciton with a (m,n) perturbation is defined as
 *
 * dfn = fac * n(psi_inv) * exp(-0.5*(m*vpar^2 + mu*B)/T(psi_inv)) * (1 + perturbation)
 * 
 * where 
 *
 * perturbation = mode_coeff[0]*cos(mode_coeff[1]*phi + mode_coeff[2]*theta) 
 *
 * fac = 1/dsqrt(PI*(two/mass)**3) / T(r_bar)^(3/2)
 *
 * One-dimensional density and temperature functions are specified as
 *
 * A(psi_inv) = A_val * Exp[-A_kappa * A_width * Tanh[(psi_inv-psi_mid)/A_width]]
 *
 * For convenience (e.g., comparison with CBC) kappa is specified in the input
 * in units of 1/R0. NB: A_kappa = A_kappa_input/(dpsidr_p * R0) 
 *
 * The following represents a sample input entry for this Function choice.
 *
 * \verbatim
 * type = "canonical_maxwellian"
 * density_val = 1.0
 * density_width = 0.2
 * density_kappa = 1.0
 * temperature_val = 1.0
 * temperature_width = 0.2
 * temperature_kappa = 1.0
 * mode_coefficients = 0.01 2 3
 *
*/
class CanonicalMaxwellianKineticFunction : public KineticFunction
{
   public:

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new GridFunction object.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
     CanonicalMaxwellianKineticFunction( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
     virtual ~CanonicalMaxwellianKineticFunction() {;}

      /// Set the function on the interior.
      /**
       * Sets the cell-averaged function value for the kinetic species
       * at the given time.  It is assumed that derived classes will check
       * that the specification is suitable for the given geometry embedded in
       * the KineticSpecies object.
       *
       * NB: The cell-averaged values include the J and BStarParallel factors,
       * i.e., the internal LevelData contains the average of fJB^*_{\parallel}
       *
       * @param[inout] species Kinetic species to intiialize.
       * @param[in] time Time value at which to evaluate the function.
       */
      virtual void assign( KineticSpecies& species, const Real& time ) const;

      /// Set the function on the boundary boxes.
      /**
       * Sets the cell-averaged function value for the kinetic species
       * at the given time on the boxes specified in the BoundaryBoxLayout.
       * It is assumed that derived classes will check that the specification is
       * suitable for the given geometry embedded in the KineticSpecies object.
       *
       * NB: The cell-averaged values include the J and BStarParallel factors,
       * i.e., the internal LevelData contains the average of fJB^*_{\parallel}
       *
       * @param[inout] species Kinetic species to intiialize.
       * @param[in] bdry_layout BoundaryBoxLayout defining boxes to fill.
       * @param[in] time Time value at which to evaluate the function.
       */
      virtual void assign( KineticSpecies& species,
                           const BoundaryBoxLayout& bdry_layout,
                           const Real& time ) const;

      /// Print object parameters.
      /**
       */
      inline void printParameters() const;

   private:

      // prohibit copying
      CanonicalMaxwellianKineticFunction( const CanonicalMaxwellianKineticFunction& );
      CanonicalMaxwellianKineticFunction& operator=( const CanonicalMaxwellianKineticFunction& );

      /// Set pointwise values at cell centers.
      /**
       */
      inline void setPointValues(FArrayBox&                    a_dfn,
                                 const Box&                    a_box,
                                 const PhaseBlockCoordSys&     a_phase_coord_sys,
                                 const CFG::MagBlockCoordSys&  a_mag_coord_sys,
                                 const FArrayBox&              a_B,
                                 const Real&                   a_mass,
                                 const Real&                   a_charge,
                                 const Real&                   a_larmor_number ) const;

      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& a_pp );

      /// Check that the geometry is annular or single null.
      /**
       */
      inline void checkGeometryValidity( const PhaseGeom& a_geometry ) const;

      int m_verbosity;
   
      int m_profile_option;
   
      Vector<Real> m_density_parm;
      Vector<Real> m_temperature_parm;
      Vector<Real> m_perturbation_parm;
};

#include "NamespaceFooter.H"

#endif
