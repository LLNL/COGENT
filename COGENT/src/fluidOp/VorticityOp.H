#ifndef  _VORTICITYOP_H_
#define  _VORTICITYOP_H_

#include "FluidOpInterface.H"
#include "GKFluidOp.H"
#include "GKPoisson.H"
#include "CellVar.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "GKVlasov.H"
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include <string>


#include "NamespaceHeader.H"

/**
 * VorticityOp operator class.
*/
class VorticityOp
   : public FluidOpInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  VorticityOp( const std::string& a_pp_str,
               const MagGeom&     a_geometry,
               const double       a_larmor,
               const int          a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~VorticityOp();

   /// Evaluates the RHS.
   /**
    *  Pure virtual function that will be resposible for evaluating
    *  field RHS on the given field.
    *
    *  @param[in,out] rhs           -  contribution to be added to rhs.
    *  @param[in] fields            -  current solution for fields.
    *  @param[in] fluids            -  current solution for fluids.
    *  @param[in] kinetic_specties  -  current solution for kinetic species.
    *  @param[in] fieldVecComp      -  component of the field vector to which operator is applied.
    *  @param[in] time              -  the time at which the field RHS is to be evaluated
    */
   virtual void accumulateExplicitRHS( FluidSpeciesPtrVect&               rhs,
                                       const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                       const FluidSpeciesPtrVect&         fluid_species,
                                       const PS::ScalarPtrVect&           scalars,
                                       const EField&                      E_field,
                                       const int                          fluid_vec_comp,
                                       const Real                         time );

   virtual void accumulateImplicitRHS( FluidSpeciesPtrVect&               rhs,
                                       const PS::KineticSpeciesPtrVect&   kinetic_species,
                                       const FluidSpeciesPtrVect&         fluid_species,
                                       const PS::ScalarPtrVect&           scalars,
                                       const EField&                      E_field,
                                       const int                          fluid_vec_comp,
                                       const Real                         time );

   virtual void preSolutionOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species,
                                   const FluidSpeciesPtrVect&         fluid_species,
                                   const PS::ScalarPtrVect&           scalars,
                                   const EField&                      E_field,
                                   const double                       time );

   virtual void preOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species,
                           const FluidSpeciesPtrVect&         fluid_species,
                           const PS::ScalarPtrVect&           scalars,
                           const EField&                      E_field,
                           const double                       time );

   virtual void defineBlockPC( std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>&,
                               std::vector<PS::DOFList>&,
                               const PS::GKVector&,
                               PS::GKOps&,
                               const std::string&,
                               const std::string&,
                               bool,
                               const FluidSpecies&,
                               const PS::GlobalDOFFluidSpecies&,
                               int );

   virtual void updateBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>&,
                                const PS::KineticSpeciesPtrVect&,
                                const FluidSpeciesPtrVect&,
                                const Real,
                                const Real,
                                const bool,
                                const int );

   virtual void evalSolutionOp( FluidSpeciesPtrVect&               rhs,
                                const PS::KineticSpeciesPtrVect&   kinetic_species_comp,
                                const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                const FluidSpeciesPtrVect&         fluid_species_comp,
                                const FluidSpeciesPtrVect&         fluid_species_phys,
                                const PS::ScalarPtrVect&           scalars,
                                const int                          component,
                                const Real                         time );

   virtual void solveSolutionPC( FluidSpeciesPtrVect&              a_fluid_species_solution,
                                 const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                 const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                 const int                         a_component );

   virtual void updatePCImEx( const PS::KineticSpeciesPtrVect& kinetic_species,
                              double                           time,
                              double                           shift );

   virtual void solvePCImEx( FluidSpeciesPtrVect&              fluid_species_solution,
                             const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                             const FluidSpeciesPtrVect&        fluid_species_rhs,
                             const int                         component );

   virtual bool trivialSolutionOp() const {return false;}

   virtual CFGVars* createStateVariable( const string&       a_pp_prefix,
                                         const std::string&  a_name,
                                         const MagGeom&      a_geometry,
                                         const IntVect&      a_ghost_vect ) const
   {
      return new CellVar(a_pp_prefix, a_name, a_geometry, 1, a_ghost_vect);
   }

   virtual void fillGhostCells( FluidSpecies&  species_phys,
                                const double   time );

   // Compute the divergence of the perpedicular ion current due solely to magnetic drifts
   void computeDivPerpIonMagCurrentDensity( LevelData<FArrayBox>&             div_Jperp_mag,
                                            const EField&                     E_field,
                                            const PS::KineticSpeciesPtrVect&  species,
                                            const Real&                       time );

   void computeDivPerpIonCurrentDensity( LevelData<FArrayBox>&             div_Jperp,
                                         const EField&                     E_field,
                                         const PS::KineticSpeciesPtrVect&  species,
                                         const LevelData<FArrayBox>&       ion_charge_density,
                                         const Real&                       time,
                                         const bool                        compute_mag_contrib = true );

   static const char* pp_name;

private:

   //Computes ion mass density
   void computeIonMassDensity( LevelData<FArrayBox>&             mass_density,
                               const PS::KineticSpeciesPtrVect&  species ) const;

   //Computes ion charge density
   void computeIonChargeDensity( LevelData<FArrayBox>&            ion_charge_density,
                                const PS::KineticSpeciesPtrVect&  species ) const;

   void setCoreBC( const double   core_inner_bv,
                   const double   core_outer_bv,
                   EllipticOpBC&  bc ) const ;

   void setZero( LevelData<FArrayBox>& data ) const;
   
   int m_verbosity;
   const MagGeom& m_geometry;
   double m_larmor;
   bool m_par_cond_op_coefs_defined;

   GKPoisson*    m_parallel_current_divergence_op;
   GKPoisson*    m_par_cond_op;
   GKPoisson*    m_imex_pc_op;
   GKPoisson*    m_gyropoisson_op;

   PS::GKVlasov*  m_vlasov;

   std::string m_opt_string;
   int m_my_pc_idx_e, m_my_pc_idx_i;

   EllipticOpBC* m_parallel_current_divergence_op_bcs;
   EllipticOpBC* m_par_cond_op_bcs;
   EllipticOpBC* m_imex_pc_op_bcs;
   EllipticOpBC* m_gyropoisson_op_bcs;

   LevelData<FArrayBox> m_volume;
   LevelData<FArrayBox> m_negativeDivJpar;
   LevelData<FArrayBox> m_divJperp;
   LevelData<FArrayBox> m_divJperp_mag;
   LevelData<FluxBox> m_face_bXE_integrals;
   LevelData<FArrayBox> m_ion_charge_density;
   LevelData<FArrayBox> m_ion_mass_density;

   bool m_sigma_div_e_coefs_set;
};


#include "NamespaceFooter.H"

#endif
