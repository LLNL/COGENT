#ifndef _GKFLUIDOP_H_
#define _GKFLUIDOP_H_

#include "REAL.H"
#include "FluidOpInterface.H"
#include "FluidSpecies.H"
#include "ParmParse.H"
#include "EMFields.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "Scalar.H"
#include "GlobalDOF.H"
#include "GKUnit.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include <string>
#include <vector>

//#include "SPACE.H"
//#include "FluidOpF_F.H"

#include "NamespaceHeader.H"
namespace PS = PS_NAMESPACE; 

/**
 * Gyrokinetic fluid operator class.
*/
class GKFluidOp
{
   public:

      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      GKFluidOp( const MagGeom&      geometry,
                 const PS::GKUnits&  units,
                 const int           verbose );

      /// Destructor.
      /**
       */
      virtual ~GKFluidOp();

      /// Accumulates the RHS of the collision operator.
      /**
       * @param[out] rhs               - data holder for rhs.
       * @param[in]  fluids            - current solution for fluids.
       * @param[in]  kinetic_specties  - current solution for kinetic species.
       * @param[in]  time current time.
       */
      virtual void accumulateRHS( FluidSpeciesPtrVect&               rhs,
                                  const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                  const FluidSpeciesPtrVect&         fluid_species,
                                  const PS::ScalarPtrVect&           scalars,
                                  const EMFields&                    EM_fields,
                                  const bool                         implicit,
                                  const Real                         time);

      void preSolutionOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                              const FluidSpeciesPtrVect&         fluid_species,
                              const PS::ScalarPtrVect&           scalars,
                              const EMFields&                    EM_fields,
                              const double                       time );

      void preOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                      const FluidSpeciesPtrVect&         fluid_species,
                      const PS::ScalarPtrVect&           scalars,
                      const EMFields&                    EM_fields,
                      const double                       time );
      
      void postTimeEval( FluidSpeciesPtrVect&  a_species_comp,
                         FluidSpeciesPtrVect&  a_species_phys,
                         const double          a_dt,
                         const double          a_time );

      void enforcePositivity( FluidSpeciesPtrVect&  a_species_comp );

      void postStageEval( FluidSpeciesPtrVect&  a_species_comp,
                          FluidSpeciesPtrVect&  a_species_phys,
                          const int             a_stage,
                          const double          a_dt,
                          const double          a_time );
      
      void getMemberVarForPlotting( LevelData<FArrayBox>&  a_Var,
                              const CFGVars&               a_fluid_vars,
                              const string                 a_var_name ) const;

      virtual void convertToPhysical( FluidSpeciesPtrVect&  a_fluid_species_phys,
                                const FluidSpeciesPtrVect&  a_fluid_species_comp ) const;

      virtual void evalSolutionOp( FluidSpeciesPtrVect&               a_rhs,
                                   const PS::KineticSpeciesPtrVect&   a_kinetic_species_comp,
                                   const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                   const FluidSpeciesPtrVect&         a_fluid_species_comp,
                                   const FluidSpeciesPtrVect&         a_fluid_species_phys,
                                   const PS::ScalarPtrVect&           a_scalars,
                                   const Real                         a_time );

      void defineMultiPhysicsPC(std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>&,
                                std::vector<PS::DOFList>&,
                                const FluidSpeciesPtrVect&,
                                const PS::GlobalDOFFluidSpeciesPtrVect&,
                                const PS::ODEVector&,
                                void*,
                                const std::string&,
                                const std::string&,
                                bool,
                                const int );

      void updateMultiPhysicsPC(std::vector<PS::Preconditioner<PS::ODEVector,PS::AppCtxt>*>&,
                                const PS::KineticSpeciesPtrVect&,
                                const FluidSpeciesPtrVect&,
                                const Real,
                                const int,
                                const int,
                                const Real,
                                const bool );

      virtual void solveSolutionPC( FluidSpeciesPtrVect&              a_fluid_species_solution,
                                    const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                    const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                    int                               a_idx = -1 );

      virtual void solvePCImEx( FluidSpeciesPtrVect&              a_fluid_species_solution,
                                const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                const std::string&                a_name,
                                int                               a_idx = -1 );

      bool trivialSolutionOp( const FluidSpeciesPtrVect& fluid_species );

      void initialize( FluidSpeciesPtrVect&  fluid_species,
                 const Real                  time );
      
      bool isInitializationConstrained(const FluidSpeciesPtrVect&  fluid_phys,
				       const int step);
   
      void applyInitializationConstraints(FluidSpeciesPtrVect&               fluid_comp,
                                          FluidSpeciesPtrVect&               fluid_phys,
                                          const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                          const EMFields&                    EM_fields,
                                          const double                       time );
   
      void fillGhostCells( FluidSpeciesPtrVect&  fluid_species,
                           const double          time );

      /// Compute a stable time step.
      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDtExplicitTI( const FluidSpeciesPtrVect&  fluids_comp );
      Real computeDtImExTI( const FluidSpeciesPtrVect&  fluids_comp );
   

      /// Compute time scale.
      /**
       * Computes and returns an estimate of the time scale of the collision term.
       */
      Real computeTimeScale( const FluidSpeciesPtrVect&  fluids);

      /// returns the fluid model associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      FluidOpInterface& fluidModel( const std::string& name ) const;

      /// returns the fluid model name associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      std::string fluidModelName( const std::string& name );

      void setStepConstKinCoeff(  const bool,
                                  const FluidSpeciesPtrVect& );

   private:

      // prevent copying
      GKFluidOp( const GKFluidOp& );
      const GKFluidOp& operator=( const GKFluidOp& );

      inline
      int getFluidModelIndex( std::map<std::string,int>::iterator& it )
      {
         return it->second;
      }
   
      bool m_verbose;
      mutable std::map<std::string,int> m_species_map;
      std::vector<FluidOpInterface*> m_fluid_model;
      std::map<std::string,int> m_fluid_model_name;
};

#include "NamespaceFooter.H"

#endif
