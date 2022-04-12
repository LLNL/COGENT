#ifndef _GKVORTICITY_H_
#define _GKVORTICITY_H_

#include "EllipticOp.H"
#include "EllipticOpBCFactory.H"
#include "GKPoisson.H"
#include "FluxSurface.H"
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

/// Gyrokintic Vorticity operator class.
/**
 * This operator class solves the mapped gyrokinetic Vorticity equation:
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
class GKVorticity
   : public GKPoisson
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] geom the geometry data.
    *  @param[in] larmor_number the Larmor number for the operator.
    *  @param[in] debye_number the Debye number for the operator.
    */
   
   GKVorticity(const ParmParse&   pp,
               const ParmParse&   pp_base,
               const MagGeom&     geom,
               const Real         larmor_number,
               const Real         debye_number,
               const bool         second_order,
               const bool         low_pollution,
               const bool         include_pol_den_correction,
               const std::string& model = "Vorticity" );
      
   /// Destructor.
   /**
    */
   virtual ~GKVorticity();

   virtual MBHypreSolver* allocatePreconditioner( const MagGeom&                  geom,
                                                  const int                       discretization_order,
                                                  MultiBlockLevelExchangeCenter*  mblx_ptr );
   
   virtual void solvePreconditioner( const LevelData<FArrayBox>& in,
                                     LevelData<FArrayBox>&       out );

   virtual void applyOp(LevelData<FArrayBox>& out,
                        const LevelData<FArrayBox>& in,
                        bool homogeneous = false);

   virtual void computeBcDivergence( LevelData<FArrayBox>& out );
   
   virtual void updatePreconditioner( MBSolver*              preconditioner,
                                      LevelData<FArrayBox>&  volume_reciprocal,
                                      LevelData<FluxBox>&    mapped_coefficients,
                                      const EllipticOpBC&    potential_bc );


   void setVorticityOperatorCoefficients( const LevelData<FArrayBox>&  ion_mass_density,
                                          const LevelData<FArrayBox>&  ion_charge_density,
                                          const LevelData<FArrayBox>&  electron_temeprature,
                                          EllipticOpBC&                bc,
                                          const bool                   update_preconditioner );

   void computeVorticityCoefficients(const LevelData<FArrayBox>& ion_mass_density,
                                     const LevelData<FArrayBox>& ion_charge_density,
                                     const LevelData<FArrayBox>& electron_temeprature,
                                     const EllipticOpBC&         bc,
                                     LevelData<FluxBox>&         mapped_coefficients,
                                     LevelData<FluxBox>&         unmapped_coefficients );

   void verifyMatrix( const MBSolver* matrix );


   void applyHighOrderCorrectionOp(LevelData<FArrayBox>&        out,
                                   const LevelData<FArrayBox>&  in,
                                   bool                         homogeneous );

   void computeFluxDivergence(const LevelData<FArrayBox>& in,
                              LevelData<FArrayBox>&       out,
                              const bool                  homogeneous_bcs,
                              const bool                  extrap_to_ghosts = false );

   void  computeFluxDivergenceWithCoeff(LevelData<FArrayBox>&       out,
                                        const LevelData<FArrayBox>& in,
                                        const LevelData<FluxBox>&   coeff_mapped,
                                        const bool                  homogeneous_bcs,
                                        const bool                  subtract_fs_par_div,
                                        const bool                  extrap_to_ghosts = false);
   
   /// Substract flux-averaged parallel divergence; used for 3D MB vorticity
   void subtractFSAverParDiv( LevelData<FArrayBox>&      divergence,
                              const  LevelData<FluxBox>& flux) const;

   
   /// Set implicit dt value
   inline void setImplicitDt(const Real a_value) {m_dt_implicit = a_value;};
   
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
   
   void computeParallelConductivity(const LevelData<FluxBox>& Te,
                                    LevelData<FluxBox>&       parallel_conductivity ) const;
   
   void computeDiffusionCoeff(LevelData<FluxBox>& diffusion_coeff_unmapped,
                              LevelData<FluxBox>& diffusion_coeff_mapped) const;
   
   void convertCellToFace(LevelData<FluxBox>&          face_data,
                          const LevelData<FArrayBox>&  cell_data ) const;

   RefCountedPtr<GridFunction> m_charge_exchange_func;
   RefCountedPtr<GridFunction> m_pol_diffusion_func;
   RefCountedPtr<GridFunction> m_perp_diffusion_func;
   RefCountedPtr<GridFunction> m_parallel_conductivity_func;

   LevelData<FluxBox> m_parallel_cond_face;
   LevelData<FluxBox> m_charge_exchange_face;
   LevelData<FluxBox> m_pol_diffusion_face;
   LevelData<FluxBox> m_perp_diffusion_face;

   bool m_include_pol_den_correction;
   bool m_include_diffusion;
   bool m_include_high_order_corr;
   bool m_symmetrized_preconditioner;

private:
   void init( const ParmParse&  pp,
              const ParmParse&  pp_base );
   
   LevelData<FluxBox> m_diffusion_coeff_mapped;
   LevelData<FluxBox> m_diffusion_coeff_unmapped;
   LevelData<FluxBox> m_M_mapped;
   LevelData<FluxBox> m_M_unmapped;
   LevelData<FluxBox> m_N2_mapped;
   LevelData<FluxBox> m_N2_unmapped;
   LevelData<FArrayBox> m_ne_over_Te;

   int m_verbosity;
   std::string m_model;

   Real m_dt_implicit;
   Real m_parallel_conductivity_limit;
   Real m_precond_scale_fac;
   bool m_subtract_fs_par_div;
   bool m_use_vorticity_bcs;
  
   FluxSurface m_flux_surface;
      
   EllipticOpBC* m_vorticity_bc;
   RefCountedPtr<EllipticOpBC> m_potential_bc;
};



#include "NamespaceFooter.H"

#endif
