#ifndef  _MY_KROOK_H_
#define  _MY_KROOK_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "REAL.H"
#include "CLSInterface.H"
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
 * Krook collisions operator class.
*/
class MyKrook
   : public CLSInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  MyKrook(const std::string&  a_ppcls_str,
         const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~MyKrook();

  /// Accumulates the RHS of the collision operator.
  /**
   * @param[out] rhs data holder for rhs (without ghost cells).
   * @param[in] soln current solution with which to evaluate RHS.
   * @param[in] species index of species to which collision is applied.
   * @param[in] time current time.
   */
  virtual void evalClsRHS( KineticSpeciesPtrVect&       rhs,
                           const KineticSpeciesPtrVect& soln,
                           const int                    species,
                           const Real                   time );

  inline
  virtual void evalClsRHSExplicit( KineticSpeciesPtrVect&        a_rhs,
                                   const KineticSpeciesPtrVect&  a_soln,
                                   const int                     a_species,
                                   const Real                    a_time )
  {
    if (!m_time_implicit) evalClsRHS(a_rhs, a_soln, a_species, a_time);
  }

  inline
  virtual void evalClsRHSImplicit( KineticSpeciesPtrVect&        a_rhs,
                                   const KineticSpeciesPtrVect&  a_soln,
                                   const int                     a_species,
                                   const Real                    a_time )
  {
    if (m_time_implicit) evalClsRHS(a_rhs, a_soln, a_species, a_time);
  }

  /// Compute a stable time step.
  /**
   * Computes and returns an estimate of the maximum stable time step.
   */
   Real computeDt(const KineticSpeciesPtrVect& soln);
   
   Real collisionFrequency() {return m_cls_freq;}
   
   void addReferenceDfn( KineticSpecies& result,
                         const Real      time, 
                         const Real      scale );

private:

   int m_verbosity;

   bool m_time_implicit;

   bool m_first_step;
   
   Real m_cls_freq;

   RefCountedPtr<KineticFunction> m_ref_func;

   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "CLS.species" section
    * of ParmParse database, where "species" is replaced by the species name
    *
    * \param[in] ppcls CLS.species input database.
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
