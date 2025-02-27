#ifndef _GKVORTICITY_H_
#define _GKVORTICITY_H_

#include "GKPoisson.H"
#include "FluxSurface.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "MagFluxAlignedMBLEXCenter.H"
#include "MBHypreSolver.H"

#include "GridFunction.H"
#include "GridFunctionLibrary.H"

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
               const bool         include_diffusion,
               const std::string& model = "Vorticity",
               const bool         allocate_precond = true);
      
   /// Destructor.
   /**
    */
   virtual ~GKVorticity();

   virtual MBHypreSolver* allocatePreconditioner( const MagGeom&                  geom,
                                                  const int                       discretization_order,
                                                  const int                       nvar,
                                                  MultiBlockLevelExchangeCenter*  mblx_ptr );
   
   virtual void solvePreconditioner( const LevelData<FArrayBox>& in,
                                     LevelData<FArrayBox>&       out );
   
   virtual void updatePreconditioner( MBSolver*              preconditioner,
                                      LevelData<FArrayBox>&  volume_reciprocal,
                                      LevelData<FluxBox>&    mapped_coefficients,
                                      const EllipticOpBC&    potential_bc );

   virtual void setVorticityOperatorCoefficients( const LevelData<FArrayBox>&  ion_mass_density,
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
   
   void setDiffusionOperatorCoeff(const LevelData<FluxBox>& pol_diffusion_face,
                                  const LevelData<FluxBox>& perp_diffusion_face);
 
   void updateParallelConductivity(const LevelData<FluxBox>& Te,
                                   LevelData<FluxBox>&       parallel_conductivity ) const;

   void updateParallelConductivity(const LevelData<FArrayBox>& Te,
                                   LevelData<FArrayBox>&       parallel_conductivity ) const;
   
   void computeParallelConductivity(const LevelData<FluxBox>& Te,
                                    LevelData<FluxBox>&       parallel_conductivity ) const;

   void computeParallelConductivity(const LevelData<FArrayBox>& Te,
                                    LevelData<FArrayBox>&       parallel_conductivity ) const;

   void setVorticityBc(EllipticOpBC& bc);
   
   void computeFluxDivergence(const LevelData<FArrayBox>& in,
                              LevelData<FArrayBox>&       out,
                              const bool                  homogeneous_bcs,
                              const bool                  extrap_to_ghosts = false );
   
   /// Substract flux-averaged parallel divergence; used for 3D MB vorticity
   void subtractFSAverParDiv( LevelData<FArrayBox>&      divergence,
                              const  LevelData<FluxBox>& flux) const;

   
   /// Set implicit dt value
   inline void setImplicitDt(const Real a_value) {m_dt_implicit = a_value;};
   
   /// Set beta value
   inline void setBeta(const Real a_value) {m_beta = a_value;};
   
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
   
   
   void convertCellToFace(LevelData<FluxBox>&          face_data,
                          const LevelData<FArrayBox>&  cell_data ) const;

   RefCountedPtr<GridFunction> m_parallel_conductivity_func;
   
   FluxSurface m_flux_surface;
   RefCountedPtr<EllipticOpBC> m_vorticity_bc;

   LevelData<FluxBox> m_M_mapped;
   LevelData<FluxBox> m_M_unmapped;
   LevelData<FluxBox> m_N2_mapped;
   LevelData<FluxBox> m_N2_unmapped;
   LevelData<FArrayBox> m_ne_over_Te;

   LevelData<FArrayBox> m_parallel_cond_cell;
   LevelData<FluxBox> m_parallel_cond_face;
   LevelData<FluxBox> m_charge_exchange_face;
   
   LevelData<FluxBox> m_mass_density_face;
   LevelData<FluxBox> m_charge_density_face;
   LevelData<FluxBox> m_electron_temperature_face;
   
   LevelData<FluxBox> m_diffusion_coeff_mapped;
   LevelData<FluxBox> m_diffusion_coeff_unmapped;

   bool m_symmetrized_preconditioner;

   bool m_subtract_fs_par_div;
   bool m_use_vorticity_bcs;

   bool m_include_high_order_corr;
   bool m_include_diffusion;
   bool m_include_pol_den_correction;
   
   std::string m_model;
   
   Real m_beta;
   
   int m_verbosity;

private:
   
   
   void init( const ParmParse&  pp,
              const ParmParse&  pp_base,
              const bool        allocate_precond);
      

   RefCountedPtr<GridFunction> m_charge_exchange_func;

   
   Real m_dt_implicit;
   Real m_parallel_conductivity_limit;
   Real m_precond_scale_fac;
   
};



#include "NamespaceFooter.H"

#endif
