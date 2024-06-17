#ifndef  _EMVORTICITYOP_H_
#define  _EMVORTICITYOP_H_

#include "VorticityOp.H"
#include "Diffusion.H"
#include "Advection.H"
#include "GKEMVorticity.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include <string>
#include "NamespaceHeader.H"

/**
 * VorticityOp operator class.
*/
class EMVorticityOp
   : public VorticityOp
{
public:

   
  /// Null constructor (used for state variable creation)
  EMVorticityOp( const MagGeom&  geometry )
      :VorticityOp(geometry),
       m_geometry(geometry),
       m_A_parallel_bcs(NULL),
       m_imex_pc_op(NULL),
       m_Apar_soln_op(NULL),
       m_vorticity_diffusion_op(NULL),
       m_nlp_Apar_op(NULL),
       m_div_jpar_op(NULL) {};
   
  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  EMVorticityOp(const ParmParse&   pp,
                const MagGeom&     geometry,
                const double       larmor,
                const double       ion_skin_depth,
                const double       electron_skin_depth,
                const int          verbosity );

  /// Destructor.
  /**
   */
  virtual ~EMVorticityOp();

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
                                       const EMFields&                    EM_fields,
                                       const int                          fluid_vec_comp,
                                       const Real                         time );

   virtual void accumulateImplicitRHS( FluidSpeciesPtrVect&               rhs,
                                       const PS::KineticSpeciesPtrVect&   kinetic_species,
                                       const FluidSpeciesPtrVect&         fluid_species,
                                       const PS::ScalarPtrVect&           scalars,
                                       const EMFields&                    EM_fields,
                                       const int                          fluid_vec_comp,
                                       const Real                         time );

   virtual void preSolutionOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species,
                                   const FluidSpeciesPtrVect&         fluid_species,
                                   const PS::ScalarPtrVect&           scalars,
                                   const EMFields&                    EM_fields,
                                   const double                       time );

   virtual void preOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species,
                           const FluidSpeciesPtrVect&         fluid_species,
                           const PS::ScalarPtrVect&           scalars,
                           const EMFields&                    EM_fields,
                           const double                       time );


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

   virtual void updatePCImEx( const FluidSpeciesPtrVect&       fluid_species,
                              const PS::KineticSpeciesPtrVect& kinetic_species,
                              const double                     time,
                              const int                        step,
                              const int                        stage,
                              const double                     shift,
                              const int                        component,
                              const std::string&               op_name );

   virtual void solvePCImEx( FluidSpeciesPtrVect&              fluid_species_solution,
                             const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                             const FluidSpeciesPtrVect&        fluid_species_rhs,
                             const std::string&                op_name,
                             const int                         component );
   
   virtual void postTimeEval( FluidSpecies&  species_comp,
                              FluidSpecies&  species_phys,
                              const Real     dt,
                              const Real     time,
                              const int      stage );

   void setAparDerivatives( const PS::KineticSpeciesPtrVect&  kinetic_species_phys,
                            const LevelData<FArrayBox>&       Apar,
                            EMFields&                         EM_fields );

   virtual bool trivialSolutionOp() const {return false;}


   virtual CFGVars* createStateVariable( const string&       a_pp_prefix,
                                         const std::string&  a_name,
                                         const MagGeom&      a_geometry,
                                         const IntVect&      a_ghost_vect ) const
   {
      return new Potentials(a_pp_prefix, a_name, a_geometry, a_ghost_vect, true);
   }

   virtual void fillGhostCells( FluidSpecies&  species_phys,
                                const double   time );

   virtual bool isInitializationConstrained(const FluidSpecies& fluid_phys,
					    const int           step);
   
   virtual void applyInitializationConstraints(FluidSpeciesPtrVect&               fluid_comp,
                                               FluidSpeciesPtrVect&               fluid_phys,
                                               const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                               const EMFields&                    EM_fields,
                                               const int                          component,
                                               const double                       time );

   static const char* pp_name;
   
protected:
   
   /// Computes RHS of the EMVorticity model
   void accumulateRHS(LevelData<FArrayBox>&              rhs_phi,
                      LevelData<FArrayBox>&              rhs_Apar,
                      const LevelData<FArrayBox>&        soln_phi,
                      const LevelData<FArrayBox>&        soln_Apar,
                      const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                      const EMFields&                    EM_fields,
                      const Real                         time );
   
   void accumulatePhiRHS(LevelData<FArrayBox>&              rhs_phi,
                         const LevelData<FArrayBox>&        soln_phi,
                         const LevelData<FArrayBox>&        soln_Apar,
                         const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                         const EMFields&                    EM_fields,
                         const Real                         time );
   
   void accumulateAparRHS(LevelData<FArrayBox>&              rhs_Apar,
                          const LevelData<FArrayBox>&        soln_phi,
                          const LevelData<FArrayBox>&        soln_Apar,
                          const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                          const EMFields&                    EM_fields,
                          const Real                         time );
   
   EllipticOpBC*  m_A_parallel_bcs;
      
   Real m_ion_skin_depth;
   Real m_electron_skin_depth;
   
   const MagGeom& m_geometry;
   
   LevelData<FArrayBox> m_volume;
   
private:
   
   /// Parse parameters.
   void parseParameters( const ParmParse& pp );
   
   /// Print parameters.
   void printParameters();
   
   
   void addDivJpar(LevelData<FArrayBox>&        rhs,
                   const LevelData<FArrayBox>&  Apar,
                   const LevelData<FArrayBox>&  jpar);
      
   void verifyMatrix(const double   shift );
      
   LevelData<FArrayBox> m_par_cond;
   
   GKVorticity*    m_vorticity_diffusion_op;
   
   GKEMVorticity* m_Apar_soln_op;
   GKEMVorticity* m_imex_pc_op;

   Diffusion*   m_nlp_Apar_op;
   Advection*   m_div_jpar_op;
   
   RefCountedPtr<EllipticOpBC>  m_jpar_bcs;
   
   bool m_is_time_implicit;
   
};


#include "NamespaceFooter.H"

#endif
