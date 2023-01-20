#ifndef _GKEMVORTICITY_H_
#define _GKEMVORTICITY_H_

#include "GKVorticity.H"

#include "NamespaceHeader.H"
namespace PS = PS_NAMESPACE;
namespace VEL = VEL_NAMESPACE;

/// EM  Vorticity operator class.
/**
 * This operator class solves the EM Vorticity equation
 *
*/
class GKEMVorticity
   : public GKVorticity
{
public:

   /// Constructor with initialization.
   /**
    *  @param[in] pp the input database.
    *  @param[in] geom the geometry data.
    *  @param[in] larmor_number the Larmor number for the operator.
    *  @param[in] debye_number the Debye number for the operator.
    */
   
   GKEMVorticity(const ParmParse&   pp,
                 const ParmParse&   pp_base,
                 const MagGeom&     geom,
                 const Real         larmor_number,
                 const Real         ion_skin_depth,
                 const Real         electron_skin_depth,
                 const bool         second_order,
                 const bool         low_pollution,
                 const bool         include_pol_den_correction,
                 const bool         include_diffusion,
                 const std::string& pc_model = "EMVorticity");
      
   /// Destructor.
   /**
    */
   virtual ~GKEMVorticity();
   
   virtual void updatePreconditioner( MBSolver*              preconditioner,
                                      LevelData<FArrayBox>&  volume_reciprocal,
                                      LevelData<FluxBox>&    mapped_coefficients,
                                      const EllipticOpBC&    bc );
   
   
   void updateEMPreconditioner(MBSolver*              preconditioner,
                               const EllipticOpBC&    phi_bc,
                               const EllipticOpBC&    Apar_bc);
   
   void updateESPreconditioner(MBSolver*              preconditioner,
                               MBSolver*              pc_Apar,
                               const EllipticOpBC&    phi_bc,
                               const EllipticOpBC&    Apar_bc);

   void updateAlfvenWavePreconditioner(MBSolver*                     preconditioner,
                                       const LevelData<FArrayBox>&   ion_mass_density,
                                       const EllipticOpBC&           phi_bc,
                                       const EllipticOpBC&           Apar_bc);
   
   void computeDivbFlux(const LevelData<FArrayBox>&  in,
                        LevelData<FArrayBox>&        out,
                        const bool                   fourth_order) const;

   void computeAltParGrad( const LevelData<FArrayBox>&  in,
                           LevelData<FArrayBox>&        out ) const;

   virtual void solvePreconditioner( const LevelData<FArrayBox>& in,
                                     LevelData<FArrayBox>&       out );
   
   virtual void applyOp(LevelData<FArrayBox>& out,
                        const LevelData<FArrayBox>& in,
                        bool homogeneous = false);

   virtual void computeBcDivergence( LevelData<FArrayBox>& out );

   virtual void setVorticityOperatorCoefficients(const LevelData<FArrayBox>&  ion_mass_density,
                                                 const LevelData<FArrayBox>&  ion_charge_density,
                                                 const LevelData<FArrayBox>&  electron_temeprature,
                                                 EllipticOpBC&                bc,
                                                 const bool                   update_preconditioner );

   void setVorticityOperatorCoefficients(const LevelData<FArrayBox>&  ion_mass_density,
                                         const LevelData<FArrayBox>&  ion_charge_density,
                                         const LevelData<FArrayBox>&  electron_temeprature,
                                         EllipticOpBC&                bc_phi,
                                         EllipticOpBC&                bc_Apar,
                                         const bool                   update_preconditioner );
   
   void computeEllipticOpCoefficients(const LevelData<FArrayBox>& ion_mass_density,
                                      const LevelData<FArrayBox>& ion_charge_density,
                                      const LevelData<FArrayBox>& electron_temperature,
                                      const EllipticOpBC&         bc,
                                      const std::string&          model,
                                      LevelData<FluxBox>&         mapped_coefficients,
                                      LevelData<FluxBox>&         unmapped_coefficients );

   void applyNegativePerpendicularLaplacian( const LevelData<FArrayBox>&  a_in,
                                             LevelData<FArrayBox>&        a_out ) const
   {
      m_neg_Laplacian_perp->multiplyMatrix(a_in, a_out);
   }
   
   void applyDivbOperator( const LevelData<FArrayBox>&  a_in,
                           LevelData<FArrayBox>&        a_out ) const
   {
      m_divb_solver->multiplyMatrix(a_in, a_out);
   }
   
   void setJparallelBc(EllipticOpBC& bc);
   
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
   
   void  updateAparRHS(LevelData<FArrayBox>&       r,
                       const LevelData<FArrayBox>& z_phi);
   
   MBHypreSolver* m_pc_Apar;
   MBHypreSolver* m_neg_Laplacian_perp;
   MBHypreSolver* m_divb_solver;
   
   Real m_ion_skin_depth_sq;
   Real m_electron_skin_depth_sq;
      
   LevelData<FArrayBox> m_charge_density_cell;
   LevelData<FArrayBox> m_electron_temperature_cell;
   LevelData<FArrayBox> m_eta;
   LevelData<FArrayBox> m_divb;
   
   RefCountedPtr<EllipticOpBC> m_jpar_bc;
   RefCountedPtr<EllipticOpBC> m_Apar_bc;
   bool m_use_jpar_bcs;
   
private:
   
   LevelData<FluxBox> m_elliptic_Apar_mapped;
   LevelData<FluxBox> m_elliptic_Apar_unmapped;
   LevelData<FluxBox> m_elliptic_phi_mapped;
   LevelData<FluxBox> m_elliptic_phi_unmapped;
   
   std::string m_pc_model;
   bool m_pc_updated;
};



#include "NamespaceFooter.H"

#endif
