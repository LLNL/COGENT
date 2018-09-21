#ifndef  _ONEFIELDOP_H_
#define  _ONEFIELDOP_H_

#include "FluidOpInterface.H"
#include "Diffusion.H"

#include "NamespaceHeader.H"

/**
 * Single fluid operator class.
*/
class OneFieldOp
   : public FluidOpInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  OneFieldOp( ParmParse&      pp,
              const MagGeom&  geometry,
              const int       verbosity );

  /// Destructor.
  /**
   */
  virtual ~OneFieldOp();

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
                                      const EField&                      E_field,
                                      const int                          component,
                                      const Real                         time);
   
   virtual void accumulateImplicitRHS(FluidSpeciesPtrVect&               rhs,
                                      const PS::KineticSpeciesPtrVect&   kinetic_species,
                                      const FluidSpeciesPtrVect&         fluid_species,
                                      const PS::ScalarPtrVect&           scalars,
                                      const EField&                      E_field,
                                      const int                          component,
                                      const bool                         recompute_kinetic_terms,
                                      const Real                         time);
   
   virtual void updatePCImEx( const PS::KineticSpeciesPtrVect&  kinetic_species,
                              const double                      mshift );

   virtual void solvePCImEx( FluidSpeciesPtrVect&              fluid_species_solution,
                             const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                             const FluidSpeciesPtrVect&        fluid_species_rhs,
                             const int                         component );

   virtual bool trivialSolutionOp() const {return true;}

   /// Compute a stable time step.
   /**
    * Computes and returns an estimate of the maximum stable time step.
    */
   Real computeDt( const FluidSpeciesPtrVect&  fluid_species );
   

private:

  
   int m_verbosity;
   const MagGeom& m_geometry;
   int m_is_time_implicit;
   
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
   
   Diffusion* m_diffusion_op;
   EllipticOpBC* m_bc;
   
   double m_D;
   RefCountedPtr<GridFunction> m_D_shape;
   RefCountedPtr<GridFunction> m_source;

   virtual void accumulateRHS(FluidSpeciesPtrVect&               rhs,
                              const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                              const FluidSpeciesPtrVect&         fluid_species,
                              const PS::ScalarPtrVect&           scalars,
                              const EField&                      E_field,
                              const int                          component,
                              const Real                         time);
};


#include "NamespaceFooter.H"

#endif
