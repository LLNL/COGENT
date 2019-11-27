#ifndef _GKFLUIDOP_H_
#define _GKFLUIDOP_H_

#include "REAL.H"
#include "FluidOpInterface.H"
#include "FluidSpecies.H"
#include "ParmParse.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "Scalar.H"
#include "GlobalDOF.H"
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
      GKFluidOp( const MagGeom&  geometry,
                 const double    larmor,
                 const int       verbose );

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
                                  const EField&                      E_field,
                                  const bool                         implicit,
                                  const Real                         time);

      void preSolutionOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                              const FluidSpeciesPtrVect&         fluid_species,
                              const PS::ScalarPtrVect&           scalars,
                              const EField&                      E_field,
                              const double                       time );

      void preOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                      const FluidSpeciesPtrVect&         fluid_species,
                      const PS::ScalarPtrVect&           scalars,
                      const EField&                      E_field,
                      const double                       time );

      virtual void convertToPhysical( FluidSpeciesPtrVect&  a_fluid_species_phys,
                                const FluidSpeciesPtrVect&  a_fluid_species_comp ) const;

      virtual void evalSolutionOp( FluidSpeciesPtrVect&               a_rhs,
                                   const PS::KineticSpeciesPtrVect&   a_kinetic_species_comp,
                                   const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                   const FluidSpeciesPtrVect&         a_fluid_species_comp,
                                   const FluidSpeciesPtrVect&         a_fluid_species_phys,
                                   const PS::ScalarPtrVect&           a_scalars,
                                   const Real                         a_time );

      void defineMultiPhysicsPC(std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>&,
                                std::vector<PS::DOFList>&,
                                const FluidSpeciesPtrVect&,
                                const PS::GlobalDOFFluidSpeciesPtrVect&,
                                const PS::GKVector&,
                                PS::GKOps&,
                                const std::string&,
                                const std::string&,
                                bool);

      void updateMultiPhysicsPC(std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>&,
                                const PS::KineticSpeciesPtrVect&,
                                const FluidSpeciesPtrVect&,
                                const Real,
                                const Real,
                                const bool );

      virtual void solveSolutionPC( FluidSpeciesPtrVect&              a_fluid_species_solution,
                                    const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                    const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                    int                               a_idx = -1 );

      void updatePC( const PS::KineticSpeciesPtrVect&  kinetic_species,
                     const FluidSpeciesPtrVect&        fluid_species,
                     const double                      time,
                     const double                      shift,
                     int                               a_idx = -1 );

      virtual void solvePCImEx( FluidSpeciesPtrVect&              a_fluid_species_solution,
                                const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                int                               a_idx = -1 );

      bool trivialSolutionOp( const FluidSpeciesPtrVect& fluid_species );

      void initialize( FluidSpeciesPtrVect&  fluid_species,
                       const Real            time );

      void fillGhostCells( FluidSpeciesPtrVect&  fluid_species,
                           const double          time );

      /// Compute a stable time step.
      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDtExplicitTI( const FluidSpeciesPtrVect&  fluids);
      Real computeDtImExTI( const FluidSpeciesPtrVect&  fluids);
   

      /// Compute time scale.
      /**
       * Computes and returns an estimate of the time scale of the collision term.
       */
      Real computeTimeScale( const FluidSpeciesPtrVect&  fluids);

      /// returns the fluid model associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      FluidOpInterface& fluidModel( const std::string& name );

      /// returns the fluid model name associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      std::string fluidModelName( const std::string& name );

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
      std::map<std::string,int> m_species_map;
      std::vector<FluidOpInterface*> m_fluid_model;
      std::map<std::string,int> m_fluid_model_name;
};

#include "NamespaceFooter.H"

#endif
