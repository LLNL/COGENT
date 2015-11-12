#ifndef  _FIXEDBCKGR_H_
#define  _FIXEDBCKGR_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "REAL.H"
#include "NTRInterface.H"
#include "ParmParse.H"
#include <sstream>

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunction.H"
#include "GridFunctionLibrary.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM



#include "NamespaceHeader.H"

/**
 * FixedBckgr neutrals operator class.
*/
class FixedBckgr
   : public NTRInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  FixedBckgr( ParmParse&  a_ppntr,
        const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~FixedBckgr();

  /// Accumulates the RHS of the collision operator.
  /**
   * @param[out] rhs data holder for rhs (without ghost cells).
   * @param[in] soln current solution with which to evaluate RHS.
   * @param[in] species index of species to which collision is applied.
   * @param[in] time current time.
   */
  virtual void evalNtrRHS( KineticSpeciesPtrVect&       rhs,
                           const KineticSpeciesPtrVect& soln,
                           const int                    species,
                           const Real                   time );

  /// Compute a stable time step.
  /**
   * Computes and returns an estimate of the maximum stable time step.
   */
  Real computeDt(const KineticSpeciesPtrVect& soln);

   
private:

   int m_verbosity;

   Real m_ntr_freq;
   bool m_fixed_ntr_freq;
   bool m_first_step;

   RefCountedPtr<KineticFunction> m_ref_func;
   RefCountedPtr<CFG::GridFunction> m_ref_temp;
   RefCountedPtr<CFG::GridFunction> m_ref_dens;

   LevelData<FArrayBox> m_temperature;
   LevelData<FArrayBox> m_sc_ntr_freq;

   
   /// Compute the test-particle collision RHS.
   /**
    * Computes the test-particle part of the collision operator.
    */
   void testPartCollRHS( LevelData<FArrayBox>& rhs_coll,
                         const LevelData<FArrayBox>& delta_F ) const;
   

  ///Computes the self-consistent ntr_freq from the density and temperature profiles
  /**
   * Computes the dimensionless normalization ntr_norm * n/T^(3/2), wher n and T are the 
   * Maxwellian density and temeprature profiles specified by grid functions, and
   * 
   * ntr_norm[s] = (4*Pi*Z^4*Coulomb_Lg/2^(3/2)) * N[1/m^3]*10^(-6) * a_charge^4 *ech[CGS]^4 
   *                * L[m]*10^2 / (a_mass^(1/2) * T[ev]^2 * (1.6*10^(-12))^2)
   *
   * or, equvivalently in SI units
   *
   * ntr_norm[s] = (1/(4*Pi*eps0^2)*Z^4*Coulomb_Lg/2^(3/2)) * N[1/m^3] * a_charge^4 * ech[SI]^2 
   *                * L[m] / (a_mass^(1/2) * T[eV]^2) 
   * 
   */

  void computeSelfConsistFreq(LevelData<FArrayBox>&       a_ntr_freq,
                              const LevelData<FArrayBox>& a_density,
                              const double                a_mass,
                              const double                a_charge      ) const;


   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "NTR.species" section
    * of ParmParse database, where "species" is replaced by the species name
    *
    * \param[in] ppntr NTR.species input database.
    */
   void parseParameters( ParmParse& pp );
   
   /// Print parameters.
   /**
    * Private method to print parameters.
    */
   void printParameters();
   
};


#include "NamespaceFooter.H"

#endif
