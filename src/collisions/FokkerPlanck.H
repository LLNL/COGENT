#ifndef  _FOKKERPLANCK_H_
#define  _FOKKERPLANCK_H_

#include "BandedMatrix.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "REAL.H"
#include "CLSInterface.H"
#include "ParmParse.H"
#include <sstream>

#include "NamespaceHeader.H"


/**
 * FokkerPlanck collisions operator class.
*/
class FokkerPlanck
   : public CLSInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  FokkerPlanck( ParmParse&  a_ppcls,
        const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~FokkerPlanck();

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
                           const Real                   time,
                           const int                    flag = 1);

  /// Compute a stable time step.
  /**
   * Computes and returns an estimate of the maximum stable time step.
   */
  Real computeDt(const KineticSpeciesPtrVect& soln);

  /// Report that collision operator is nonlinear
  /**
   * Reports that this collision operator is nonlinear
   */
  inline bool isLinear() { return false; }

  bool setupPrecondMatrix(void *a_P,int a_N,int a_bs)
    {
      BandedMatrix *P = (BandedMatrix*) a_P;
      int nbands = 9; /* central differences of derivatives */
      P->define(a_N,nbands,a_bs);
      return(P->isDefined());
    }

  void assemblePrecondMatrix(void*,const KineticSpeciesPtrVect&,int);

private:

   int m_verbosity;
   ParmParse m_ppcls;

   bool m_fixed_cls_freq;
   Real m_cls_freq;
   Real m_cls_norm;
   Real m_pcg_tol;
   Real m_pcg_maxiter;
   int  m_mult_num;
   int  m_nD;

   int  m_update_freq;
   int  m_it_counter;

   bool m_first_step;
   bool m_subtract_background;

   LevelData<FArrayBox> m_phi;
   LevelData<FArrayBox> m_phi_F0;
   LevelData<FArrayBox> m_D;
   LevelData<FArrayBox> m_D_F0;

   RefCountedPtr<KineticFunction> m_ref_func;
   
   /// Computes the Rosenbluth potentials

   void evalRosenbluthPotentials( LevelData<FArrayBox>& phi,
                                  const PhaseGeom& phase_geom,
                                  const LevelData<FArrayBox>& dfn,
                                  const double mass ) const;

  /// Computes the coefficients D_v, D_mu, D_v_v, D_v_mu, D_mu_mu
  //  from the Rosenbluth potentials
  
  void evalCoefficients( LevelData<FArrayBox>& D,
                         const LevelData<FArrayBox>& phi,
                         const PhaseGeom& phase_geom,
                         const double mass_tp,
                         const double mass_fp) const;
  

  /// Computes the collision flux 

   void computeFlux( LevelData<FluxBox>& flux,
                     const PhaseGeom& phase_geom,
                     const LevelData<FArrayBox>& D,
                     const LevelData<FArrayBox>& dfn) const;


  ///Computes the normalization factor 
  /**
   * Computes the normalization factor 
   * 
   * cls_norm = dfn_norm[s^3/m^6] * tau_norm[s] * 
   *            Coulomb_Lg * a_charge^4 * ech[SI]^4 / (a_mass * M[kg] * eps0 )^2
   *
   * dfn_norm = (1/Pi) * N[1/m^3] / (T[J]/M[kg])^(3/2)
   *
   * tau_norm = L[m] / Sqrt(T[J]/M[kg])
   * 
   */
  void computeClsNorm(Real&                       a_cls_norm,
                      const double                a_mass,
                      const double                a_charge      ) const;


  /// Fourth-order conversion utilities

  void convertToCellAverage(const PhaseGeom&  coord_sys,
                            LevelData<FArrayBox>&   dfn) const;

  void convertToCellCenters(const PhaseGeom&  coord_sys,
                            LevelData<FArrayBox>&  dfn) const;
  
  void convertToCellFaces(LevelData<FluxBox>& faceData,
                          const LevelData<FArrayBox>& cellData,
                          const PhaseGeom& coord_sys) const;

  void FillGhostCells(const PhaseGeom&  coord_sys,
                            LevelData<FArrayBox>&   dfn) const;

   
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

#define _IndexMapping_(p,i,N,ndims,offset) \
  { \
    p = i[ndims-1]; \
    for (int c=ndims-2; c>=0; c--) p = (p*N[c]) + i[c];\
    p += offset; \
  }

#include "NamespaceFooter.H"

#endif
