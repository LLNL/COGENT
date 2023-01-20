#ifndef  _BURGERSOP_H_
#define  _BURGERSOP_H_

#include "FluidOpInterface.H"
#include "FluidVarBC.H"
#include "Diffusion.H"
#include "SpaceUtils.H.multidim"

#include "NamespaceHeader.H"

/**
 * multi-dim burgers equation operator class.
 *
 * df/dt + d/dx(cx*f^2/2) + d/dy(cy*f^2/2) + d/dz(cz*f^2/2) = D*nabla^2f
 *
*/
class BurgersOp
   : public FluidOpInterface
{
public:

  /// Null constructor (used for state variable creation) 
  BurgersOp( const MagGeom&  geometry )
     : m_geometry(geometry) {};

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  BurgersOp( const string&   pp_str,
             const string&   species_name,
             const MagGeom&  geometry,
             const int       verbosity );

  /// Destructor.
  /**
   */
  virtual ~BurgersOp();

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
   Real computeDt( const FluidSpeciesPtrVect&  fluid_species );
   
   static const char* pp_name;

private:

   int m_verbosity;
   const MagGeom& m_geometry;
   int m_is_time_implicit;
   std::string m_advScheme;
   bool m_advection, m_useLaxSplitting;
   Vector<Real> m_advCoef;

   std::string m_opt_string;
   int m_my_pc_idx;
   
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
   
   Diffusion*        m_diffusion_op;
   EllipticOpBC*     m_bc;
   FluidVarBC*       m_fluid_variable_bc;
  
   RefCountedPtr<GridFunction> m_D_rad;
   RefCountedPtr<GridFunction> m_D_par;
   RefCountedPtr<GridFunction> m_D_perp;
   RefCountedPtr<GridFunction> m_electron_dens;
   RefCountedPtr<GridFunction> m_electron_temp;

   virtual void accumulateRHS(FluidSpeciesPtrVect&               rhs,
                              const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                              const FluidSpeciesPtrVect&         fluid_species,
                              const PS::ScalarPtrVect&           scalars,
                              const EMFields&                    EM_fields,
                              const int                          component,
                              const Real                         time);
   
   void computeDiffusionCoefficients(LevelData<FluxBox>& D_tensor,
                                     LevelData<FluxBox>& D_tensor_mapped,
                                     const Real          time) const;

};


#include "NamespaceFooter.H"

#endif
