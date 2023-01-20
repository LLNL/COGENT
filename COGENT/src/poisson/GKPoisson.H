#ifndef _GKPOISSON_H_
#define _GKPOISSON_H_

#include "EllipticOp.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "MBHypreSolver.H"

#include "MagFluxAlignedMBLEXCenter.H"

#include "GridFunction.H"
#include "GridFunctionLibrary.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "KineticSpecies.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H"
namespace PS = PS_NAMESPACE;
namespace VEL = VEL_NAMESPACE;

typedef enum { gkp_ct_no_flr, gkp_ct_flr1, gkp_ct_flr2 } LHSTermType;


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
   : public EllipticOp
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
              const Real         debye_number,
              const bool         allocate_precond );

   GKPoisson( const ParmParse&,
              const MagGeom&,
              const PS::KineticSpeciesPtrVect&,
              const Real,
              const Real,
              const bool );

   GKPoisson( const ParmParse&   pp,
              const MagGeom&     geom,
              const Real         larmor_number,
              const Real         debye_number,
              const bool         allocate_precond,
              const bool         second_order,
              const bool         low_pollution );
      
   /// Destructor.
   /**
    */
   virtual ~GKPoisson();

   virtual MBHypreSolver* allocatePreconditioner( const MagGeom&                  geom,
                                                  const int                       discretization_order,
                                                  MultiBlockLevelExchangeCenter*  mblx_ptr );

   /// Construct the GKP operator coefficients
   /**
    * Computes the GKP operator coefficients.  The polarization density
    * term requires the mass-weighted sum of the ion charge densities,
    * which is obtained from the species vector
    *
    * @param[in] species   Array of kinetic species
    */

   virtual void setOperatorCoefficients( const LevelData<FArrayBox>& ni,
                                         EllipticOpBC&               bc,
                                         const bool                  update_preconditioner );

   virtual void setOperatorCoefficients( const PS::KineticSpeciesPtrVect& kin_species,
                                         const LevelData<FArrayBox>&      ni,
                                         EllipticOpBC&                    bc,
                                         const bool                       update_preconditioner );

   virtual void setOperatorCoefficients( const LevelData<FArrayBox>& ion_mass_density,
                                         EllipticOpBC&               bc,
                                         const bool                  update_preconditioner,
                                         double&                     lo_value,
                                         double&                     hi_value,
                                         LevelData<FArrayBox>&       radial_gkp_divergence_average);
   
   virtual void computeCoefficients( const LevelData<FArrayBox>&,
                                     const EllipticOpBC&  bc,
                                     LevelData<FluxBox>&,
                                     LevelData<FluxBox>& );

   void computeCoefficientsFLR_1( const LevelData<FArrayBox>&,
                                  LevelData<FluxBox>&,
                                  LevelData<FluxBox>& );

   void computeCoefficientsFLR_2( const LevelData<FArrayBox>&,
                                  LevelData<FluxBox>&,
                                  LevelData<FluxBox>& );

   void updateBoundaryDataAddingCoefs( EllipticOpBC&  bc,
                                       const bool     homogeneous_bcs = false)
   {
      updateBoundaryData(m_unmapped_coefficients, bc, homogeneous_bcs);
   }
   
   virtual void updatePreconditioner( MBSolver*              preconditioner,
                                      LevelData<FArrayBox>&  volume_reciprocal,
                                      LevelData<FluxBox>&    mapped_coefficients,
                                      const EllipticOpBC&    bc );

   virtual void solvePreconditioner( const LevelData<FArrayBox>& in,
                                     LevelData<FArrayBox>&       out );
   
   virtual void applyPCOp(const LevelData<FArrayBox>& in_vec,
                          LevelData<FArrayBox>& matvec);

   virtual void multiplyCoefficients( LevelData<FluxBox>& data,
                                      const bool mapped_coeff ) const;

   virtual void applyOp(  LevelData<FArrayBox>& out,
                          const LevelData<FArrayBox>& in,
                          bool homogeneous = false);

   void setFLRObjects( const PS::KineticSpeciesPtrVect&);

   bool isBoussinesq() const {return m_boussinesq;}
   
   inline void setModel(const std::string& model) {m_model = model;};

   /// parseParameters function
   /**
    * parseParameters function
    */
   void parseParameters( const ParmParse&   a_ppntr );

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
   
   Real m_debye_number2;
   Real m_alpha;
   
protected:

   void getMinMax(LevelData<FArrayBox>& density, double& min, double& max) const;

   void getMinMax(LevelData<FluxBox>& density, double& min, double& max) const;

   void fillDensityGhosts(LevelData<FArrayBox>& a_density) const;
   
   void computeParallelConductivity(const LevelData<FluxBox>& Te,
                                    LevelData<FluxBox>&       parallel_conductivity ) const;
   
   MBHypreSolver* m_preconditioner;

   Real m_larmor_number2;

   int m_density_interp_ghosts;
   int m_discretization_order;
   
   LevelData<FluxBox> m_mapped_coefficients;
   LevelData<FluxBox> m_unmapped_coefficients;

   LevelData<FluxBox> m_mapped_coeffs_FLR_1;
   LevelData<FluxBox> m_mapped_coeffs_FLR_2;
   LevelData<FluxBox> m_unmapped_coeffs_FLR_1;
   LevelData<FluxBox> m_unmapped_coeffs_FLR_2;

   std::vector<PS::LevelData<PS::FArrayBox>* > m_FLR_integrand_factor;

   RefCountedPtr<GridFunction> m_FLR_switch_func;
   LevelData<FArrayBox> m_FLR_bc_factor;
   int m_FLR_bc_buffer;
   bool m_plot_FLR_switch;
   
   MultiBlockLevelExchangeCenter* m_mblx_ptr;

   bool m_include_FLR_effects;
   PS::KineticSpeciesPtrVect m_species_vec;
   const PS::PhaseGeom* m_phase_geom;
   const PS::DisjointBoxLayout* m_phase_grids;
   const PS::MomentOp& m_moment_op;

   void computeFLRIntegrandFactor();

   string m_precond_method;
   double m_precond_tol;
   int m_precond_max_iter;
   bool m_precond_verbose;

   string m_precond_precond_method;
   double m_precond_precond_tol;
   int m_precond_precond_max_iter;
   bool m_precond_precond_verbose;

private:
   void init( const ParmParse& pp,
              const Real       larmor_number,
              const Real       debye_number,
              const bool       initialize_precond );
   
   int m_verbosity;
   bool m_boussinesq;
   std::string m_model;

   LHSTermType m_LHS_term_type;

   inline void setLHSTermType(const LHSTermType a_arg) 
   { 
      m_LHS_term_type = a_arg; 
   }
   inline void unsetLHSTermType() { m_LHS_term_type = gkp_ct_no_flr; }
   inline LHSTermType lhsTermType() const { return m_LHS_term_type; }

};



#include "NamespaceFooter.H"

#endif
