#ifndef _GKVORTICITY_H_
#define _GKVORTICITY_H_

#include "EllipticOp.H"
#include "GKPoisson.H"
#include "FluxSurface.H"
#include "MultiBlockLevelExchangeAverage.H"
#ifdef with_petsc
#include "MBPETScSolver.H"
#else
#include "MBHypreSolver.H"
#endif

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
               const bool         include_pol_den_correction = false );
   GKVorticity(const ParmParse&   pp,
               const ParmParse&   pp_base,
               const MagGeom&     geom,
               const Real         larmor_number,
               const Real         debye_number,
               const bool         second_order,
               const bool         low_pollution,
               const bool         include_pol_den_correction = false );
      
   /// Destructor.
   /**
    */
   virtual ~GKVorticity();

#ifdef with_petsc
   virtual MBPETScSolver* allocatePreconditioner( const MagGeom&                  geom,
                                                  const int                       discretization_order,
                                                  MultiBlockLevelExchangeCenter*  mblx_ptr );
   }
#else
   virtual MBHypreSolver* allocatePreconditioner( const MagGeom&                  geom,
                                                  const int                       discretization_order,
                                                  MultiBlockLevelExchangeCenter*  mblx_ptr );
#endif

   void setVorticityOperatorCoefficients( const LevelData<FArrayBox>&  ion_mass_density,
                                          const LevelData<FArrayBox>&  ion_charge_density,
                                          const EllipticOpBC&          bc,
                                          const bool                   update_preconditioner );

   /// Construct the GKP operator coefficients
   /**
    * Computes the GKP operator coefficients.  The polarization density
    * term requires the mass-weighted sum of the ion charge densities,
    * which is obtained from the species vector
    *
    * @param[in] species   Array of kinetic species
    */

   
   virtual void computeCoefficients(const LevelData<FArrayBox>& ni,
                                    LevelData<FluxBox>& mapped_coefficients,
                                    LevelData<FluxBox>& unmapped_coefficients );

   virtual void updatePreconditioner( MBSolver*              preconditioner,
                                      LevelData<FArrayBox>&  volume_reciprocal,
                                      LevelData<FluxBox>&    mapped_coefficients,
                                      const EllipticOpBC&    bc );

   void verifyMatrix( const MBSolver* matrix );

   virtual void solvePreconditioner( const LevelData<FArrayBox>& in,
                                     LevelData<FArrayBox>&       out );

   virtual void applyOp(LevelData<FArrayBox>& out,
                        const LevelData<FArrayBox>& in,
                        bool homogeneous = false);

   void applyPolDenCorrectionOp(LevelData<FArrayBox>&        out,
                                const LevelData<FArrayBox>&  in,
                                bool                         homogeneous );

   virtual void computeBcDivergence( LevelData<FArrayBox>& out );

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
   
   void setPolDenCorrectionCoeff(const LevelData<FArrayBox>& ion_charge_density );

   /// Substract flux-averaged parallel divergence; used for 3D MB vorticity
   void subtractFSAverParDiv( LevelData<FArrayBox>&      divergence,
                              const  LevelData<FluxBox>& flux) const;

   inline void setModel(const std::string& model) {m_model = model;};
   
   /// Return perpendicular coefficeints (those of the perpendicular Laplacian operator, M)
   const LevelData<FluxBox>& getPerpCoeffMapped() const {return m_perp_coeff_mapped;}
   
   /// Return parallel coefficeints (those of the parallel current operator)
   const LevelData<FluxBox>& getParCoeffMapped() const {return m_par_coeff_mapped;}

   /// Return face-cenetered electron temperature
   const LevelData<FluxBox>& getFCElectronTemp() const {return m_electron_temperature_face;}
   
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

   RefCountedPtr<GridFunction> m_electron_temperature;
   RefCountedPtr<GridFunction> m_charge_exchange_coeff;
   RefCountedPtr<GridFunction> m_parallel_conductivity;

   LevelData<FArrayBox> m_electron_temperature_cell;
   LevelData<FluxBox> m_electron_temperature_face;
   LevelData<FluxBox> m_parallel_cond_face;
   LevelData<FluxBox> m_charge_exchange_coeff_face;

   bool m_include_pol_den_correction;

private:
   void init( const ParmParse&  pp,
              const ParmParse&  pp_base );
   
   LevelData<FluxBox> m_perp_coeff_mapped;
   LevelData<FluxBox> m_par_coeff_mapped;
   LevelData<FluxBox> m_perp_coeff_unmapped;
   LevelData<FluxBox> m_par_coeff_unmapped;
   LevelData<FArrayBox> m_ne_over_Te;
   int m_verbosity;
   std::string m_model;

   Real m_dt_implicit;
   bool m_subtract_fs_par_div;
   FluxSurface m_flux_surface;
};



#include "NamespaceFooter.H"

#endif
