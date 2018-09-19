#ifndef  _VORTICITYOP_H_
#define  _VORTICITYOP_H_

#include "FluidOpInterface.H"
#include "GKFluidOp.H"
#include "GKPoisson.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "GKVlasov.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

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
  VorticityOp( ParmParse&      a_pp,
               const MagGeom&  a_geometry,
               const double    a_larmor,
               const int       a_verbosity );

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
                                       const bool                         recompute_kinetic_terms,
                                       const Real                         time );

   virtual void evalSolutionOp( FluidSpeciesPtrVect&               rhs,
                                const PS::KineticSpeciesPtrVect&   kinetic_species,
                                const FluidSpeciesPtrVect&         fluid_species,
                                const PS::ScalarPtrVect&           scalars,
                                const int                          component,
                                const Real                         time );

   virtual void solveSolutionPC( FluidSpeciesPtrVect&              a_fluid_species_solution,
                                 const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                 const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                 const int                         a_component );

   virtual void updatePCImEx( const PS::KineticSpeciesPtrVect& kinetic_species,
                              double                           mshift );

   virtual void solvePCImEx( FluidSpeciesPtrVect&              fluid_species_solution,
                             const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                             const FluidSpeciesPtrVect&        fluid_species_rhs,
                             const int                         component );

   virtual bool trivialSolutionOp() const {return false;}

   void computeDivPerpIonCurrentDensity( LevelData<FArrayBox>&             div_Jperp,
                                         const EField&                     E_field,
                                         const PS::KineticSpeciesPtrVect&  species,
                                         const Real&                       time ) const;

private:

   //Computes ion mass density
   void computeIonMassDensity( LevelData<FArrayBox>&             mass_density,
                               const PS::KineticSpeciesPtrVect&  species ) const;

   //Computes ion charge density
   void computeIonChargeDensity( LevelData<FArrayBox>&            ion_charge_density,
                                const PS::KineticSpeciesPtrVect&  species ) const;

   void setCoreBC( const double  core_inner_bv,
                   const double  core_outer_bv,
                   PotentialBC&  bc ) const ;

   void setZero( LevelData<FArrayBox>& data ) const;
   
   int m_verbosity;
   const MagGeom& m_geometry;
   bool m_par_cond_op_coefs_defined;

   GKPoisson*    m_parallel_current_divergence_op;
   GKPoisson*    m_par_cond_op;
   GKPoisson*    m_imex_pc_op;
   GKPoisson*    m_gyropoisson_op;

   PS::GKVlasov*  m_vlasov;

   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "field_model.field_name" section
    * of ParmParse database, where "field_name" is replaced by the field name (e.g., field_model.VectPotential)
    *
    * \param[in] ppcls field_model.field input database.
    */
   void parseParameters( ParmParse& pp );
   
   /// Print parameters.
   /**
    * Private method to print parameters.
    */
   void printParameters();

   PotentialBC* m_parallel_current_divergence_op_bcs;
   PotentialBC* m_par_cond_op_bcs;
   PotentialBC* m_imex_pc_op_bcs;
   PotentialBC* m_gyropoisson_op_bcs;

   LevelData<FArrayBox> m_negativeDivJpar;
   LevelData<FArrayBox> m_divJperp;

   bool m_sigma_div_e_coefs_set;
};


#include "NamespaceFooter.H"

#endif
