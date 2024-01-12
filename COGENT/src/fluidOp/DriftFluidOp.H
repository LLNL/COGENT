#ifndef  _DRIFTFLUIDOP_H_
#define  _DRIFTFLUIDOP_H_

#include "FluidOpInterface.H"
#include "FluidVarBC.H"
#include "Diffusion.H"
#include "FluidOpUtils.H"
#include "SpaceUtils.H.multidim"

#include "NamespaceHeader.H"

/**
 * Drif- fluid equations operator class.
 *
 *  Presently solves
 *    Physics model: d(n)/dt = -1/c * ExB/B^2 * nabla(n0)
 *    COGENT equation: d(J*n)/dt = -J * Larmor * ExB/B^2 * nabla(n0)
 *
 *
*/
class DriftFluidOp
   : public FluidOpInterface
{
public:

  /// Null constructor (used for state variable creation) 
  DriftFluidOp( const MagGeom&  geometry )
   : m_geometry(geometry),
     m_fluid_variable_bc(NULL),
     m_fluid_op_utils(NULL)  {};

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  DriftFluidOp(const string&   pp_str,
               const string&   species_name,
               const MagGeom&  geometry,
               const Real      larmor,
               const int       verbosity );

  /// Destructor.
  /**
   */
  virtual ~DriftFluidOp();

   /// Evaluates the RHS.
   /**
    *  Pure virtual function that will be resposible for evaluating
    *  field RHS on the given field.
    *
    *  @param[in,out] rhs           -  contribution to be added to rhs.
    *  @param[in] fields            -  current solution for fields.
    *  @param[in] fluids            -  current solution for fluids.
    *  @param[in] kinetic_specties  -  current solution for kinetic species.
    *  @param[in] component         -  component of the field vector being computed
    *  @param[in] time              -  the time at which the field RHS is to be evaluated
    */
   virtual void accumulateExplicitRHS(FluidSpeciesPtrVect&               rhs,
                                      const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                      const FluidSpeciesPtrVect&         fluid_species,
                                      const PS::ScalarPtrVect&           scalars,
                                      const EMFields&                    EM_fields,
                                      const int                          component,
                                      const Real                         time);
   
   virtual void accumulateImplicitRHS(FluidSpeciesPtrVect&               rhs,
                                      const PS::KineticSpeciesPtrVect&   kinetic_species,
                                      const FluidSpeciesPtrVect&         fluid_species,
                                      const PS::ScalarPtrVect&           scalars,
                                      const EMFields&                    EM_fields,
                                      const int                          component,
                                      const Real                         time);
   
   virtual void defineBlockPC( std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>&,
                               std::vector<PS::DOFList>&,
                               const PS::ODEVector&,
                               void*,
                               const std::string&,
                               const std::string&,
                               bool,
                               const FluidSpecies&,
                               const PS::GlobalDOFFluidSpecies&,
                               const int,
                               const int );

   virtual void updateBlockPC(  std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>&,
                                const PS::KineticSpeciesPtrVect&,
                                const FluidSpeciesPtrVect&,
                                const Real,
                                const int,
                                const int,
                                const Real,
                                const bool,
                                const int );

   virtual void updatePCImEx( const FluidSpeciesPtrVect&        fluid_species,
                              const PS::KineticSpeciesPtrVect&  kinetic_species,
                              const double                      time,
                              const int                         step,
                              const int                         stage,
                              const double                      shift,
                              const int                         component,
                              const std::string& );

   virtual void solvePCImEx( FluidSpeciesPtrVect&              fluid_species_solution,
                             const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                             const FluidSpeciesPtrVect&        fluid_species_rhs,
                             const std::string&,
                             const int                         component );

   virtual bool trivialSolutionOp() const {return true;}

   virtual CFGVars* createStateVariable( const string&       a_pp_prefix,
                                         const std::string&  a_name,
                                         const MagGeom&      a_geometry,
                                         const IntVect&      a_ghost_vect ) const
   {
      return new FluidSpecies(a_pp_prefix, a_name, a_geometry, a_ghost_vect);
   }

   virtual void fillGhostCells( FluidSpecies&  species_phys,
                                const double   time );

   /// Compute a stable time step.
   /**
    * Computes and returns an estimate of the maximum stable time step.
    */
   Real computeDtExplicitTI( const FluidSpeciesPtrVect&  fluid_comp );
   Real computeDtImExTI( const FluidSpeciesPtrVect&  fluid_comp );
   Real computeCourantTimeStep( const FluidSpeciesPtrVect&  fluid_comp );
  
   virtual void postTimeEval( FluidSpecies&  species_comp,
                              FluidSpecies&  species_phys,
                              const Real     dt,
                              const Real     time,
                              const int      stage );
   
   static const char* pp_name;

   Vector<int> m_harmonic_filtering;
   
private:
   
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
   
   
   virtual void accumulateRHS(FluidSpeciesPtrVect&               rhs,
                               const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                               const FluidSpeciesPtrVect&         fluid_species,
                               const PS::ScalarPtrVect&           scalars,
                               const EMFields&                    EM_fields,
                               const int                          component,
                               const Real                         time);

   void addAdvection(LevelData<FArrayBox>&         rhs,
                     const LevelData<FArrayBox>&   soln,
                     const EMFields&               EM_fields,
                     const Real                    time);
  
   void subtractExBGradN(LevelData<FArrayBox>&        rhs,
                         const LevelData<FArrayBox>&  soln,
                         const EMFields&              EM_fields,
                         const Real                   time);
   
   void subtractExBGradN0(LevelData<FArrayBox>&  rhs,
                          const EMFields&        EM_fields,
                          const Real             time);
   
   RefCountedPtr<GridFunction> m_n0_func;
   RefCountedPtr<GridFunction> m_D_par_func;
   LevelData<FArrayBox> m_n0;
   LevelData<FArrayBox> m_n0_gradient;
   LevelData<FluxBox> m_D_par_fc;
   
   LevelData<FluxBox> m_E_field;
   
   const MagGeom& m_geometry;
   
   Real m_larmor;
   Real m_max_dt;
   
   int m_verbosity;
   int m_my_pc_idx;
   
   bool m_first_call;
   bool m_is_time_implicit;
   bool m_is_linearized;
   bool m_include_hyperviscosity;
   bool m_use_hyperviscosity_bcs;
   bool m_include_boundary_relaxation;
   bool m_time_step_diagnostics;

   std::string m_species_name;
   std::string m_opt_string;
   std::string m_advection_scheme;
    
   RefCountedPtr<EllipticOpBC>   m_hyperviscosity_op_bcs;

   FluidVarBC* m_fluid_variable_bc;
   
   FluidOpUtils* m_fluid_op_utils;
     
};


#include "NamespaceFooter.H"

#endif
