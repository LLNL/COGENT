#ifndef  _FLUIDOPINTERFACE_H_
#define  _FLUIDOPINTERFACE_H_

#include "FluidSpecies.H"
#include "EField.H"
#include "inspect.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "KineticSpecies.H"
#include "Scalar.H"
#include "GlobalDOF.H"
#include "MultiPhysicsPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "NamespaceHeader.H"
class GKVector;
class GKOps;
#include "NamespaceFooter.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H"
namespace PS = PS_NAMESPACE;

/**
 * FluidOp interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set operators on fluid species.
*/
class FluidOpInterface
{
   public:

      /// Destructor.
      /**
       */
      virtual ~FluidOpInterface() {;}

      /// Evaluates the RHS.
      /**
       *  Pure virtual function that will be resposible for evaluating
       *  fluid RHS on the given species.
       *
       *  @param[in,out] rhs           -  contribution to be added to rhs.
       *  @param[in] kinetic_specties  -  current solution for kinetic species.
       *  @param[in] fluid_species     -  current solution for fluid species.
       *  @param[in] component         -  fluid component to be computed
       *  @param[in] time              -  the time at which the fluid RHS is to be evaluated.
       */
      virtual void accumulateExplicitRHS( FluidSpeciesPtrVect&               rhs,
                                          const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                          const FluidSpeciesPtrVect&         fluid_species,
                                          const PS::ScalarPtrVect&           scalars,
                                          const EField&                      E_field,
                                          const int                          component,
                                          const Real                         time) = 0;

      virtual void accumulateImplicitRHS( FluidSpeciesPtrVect&               rhs,
                                          const PS::KineticSpeciesPtrVect&   kinetic_species,
                                          const FluidSpeciesPtrVect&         fluid_species,
                                          const PS::ScalarPtrVect&           scalars,
                                          const EField&                      E_field,
                                          const int                          component,
                                          const Real                         time) = 0;

      virtual void preSolutionOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species,
                                      const FluidSpeciesPtrVect&         fluid_species,
                                      const PS::ScalarPtrVect&           scalars,
                                      const EField&                      E_field,
                                      const double                       time ) {;}

      virtual void preOpEval( const PS::KineticSpeciesPtrVect&   kinetic_species,
                              const FluidSpeciesPtrVect&         fluid_species,
                              const PS::ScalarPtrVect&           scalars,
                              const EField&                      E_field,
                              const double                       time ) {;}

      virtual void defineBlockPC( std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>&,
                                  std::vector<PS::DOFList>&,
                                  const PS::GKVector&,
                                  PS::GKOps&,
                                  const std::string&,
                                  const std::string&,
                                  bool,
                                  const FluidSpecies&,
                                  const PS::GlobalDOFFluidSpecies&,
                                  int ) { }

      virtual void updateBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>&,
                                   const PS::KineticSpeciesPtrVect&,
                                   const FluidSpeciesPtrVect&,
                                   const Real,
                                   const Real,
                                   const bool,
                                   const int ) { }

      virtual void evalSolutionOp( FluidSpeciesPtrVect&               rhs,
                                   const PS::KineticSpeciesPtrVect&   kinetic_species_comp,
                                   const PS::KineticSpeciesPtrVect&   kinetic_species_phys,
                                   const FluidSpeciesPtrVect&         fluid_species_comp,
                                   const FluidSpeciesPtrVect&         fluid_species_phys,
                                   const PS::ScalarPtrVect&           scalars,
                                   const int                          component,
                                   const Real                         time )
      {
         // The default solution operator is just the identity
         rhs[component]->copy(*(fluid_species_comp[component]));
      }

      virtual void updatePCImEx( const FluidSpeciesPtrVect&       fluid_species,
                                 const PS::KineticSpeciesPtrVect& kinetic_species,
                                 const double                     time,
                                 const double                     shift,
                                 const int                        component) {}

      virtual void solveSolutionPC( FluidSpeciesPtrVect&              fluid_species_solution,
                                    const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                                    const FluidSpeciesPtrVect&        fluid_species_rhs,
                                    const int                         component )
      {
         // The default solution operator is just the identity
         fluid_species_solution[component]->copy(*(fluid_species_rhs[component]));
      }

      virtual void updatePC( const PS::KineticSpeciesPtrVect&  kinetic_species,
                             const EllipticOpBC&               bc ) {};

      virtual void solvePCImEx( FluidSpeciesPtrVect&              fluid_species_solution,
                                const PS::KineticSpeciesPtrVect&  kinetic_species_rhs,
                                const FluidSpeciesPtrVect&        fluid_species_rhs,
                                const int                         component )
      {
         // The default solution operator is just the identity
         fluid_species_solution[component]->copy(*(fluid_species_rhs[component]));
      }

      virtual CFGVars* createStateVariable( const string&       pp_prefix,
                                            const std::string&  name,
                                            const MagGeom&      geometry,
                                            const IntVect&      ghost_vect ) const = 0;

      virtual void initialize( CFGVars&    a_species,
                               const Real  a_time )
      {
         CH_TIME("FluidOpInterface::initialize()");
         for (int n=0; n<a_species.num_cell_vars(); ++n) {
            
            const int num_comps = a_species.cell_var(n).nComp();
            if (num_comps>1) {
               string this_cell_var_comp_name;
               for (int m=0; m<num_comps; m++) {
                  if(m==0) this_cell_var_comp_name = a_species.cell_var_name(n)+"_0";
                  if(m==1) this_cell_var_comp_name = a_species.cell_var_name(n)+"_1";
                  if(m==2) this_cell_var_comp_name = a_species.cell_var_name(n)+"_2";
                  const GridFunction& fsic( fluidSpeciesIC( this_cell_var_comp_name ) );
                  LevelData<FArrayBox> this_cell_var_comp( a_species.cell_var(n).getBoxes(), 1, 
                                                           a_species.cell_var(n).ghostVect() );
                  fsic.assign( this_cell_var_comp, a_species.configurationSpaceGeometry(), a_time );
                  LevelData<FArrayBox>& cell_var_ref( a_species.cell_var(n) );
                  for (DataIterator dit( this_cell_var_comp.dataIterator() ); dit.ok(); ++dit) {
                     cell_var_ref[dit].copy(this_cell_var_comp[dit],0,m,1);
                  }
               }
            } 
            else {
               const GridFunction& fsic( fluidSpeciesIC( a_species.cell_var_name(n) ) );
               fsic.assign( a_species.cell_var(n), a_species.configurationSpaceGeometry(), a_time );
            }
            
         }
         for (int n=0; n<a_species.num_face_vars(); ++n) {
            const GridFunction& fsic( fluidSpeciesIC( a_species.face_var_name(n) ) );
            fsic.assign( a_species.face_var(n), a_species.configurationSpaceGeometry(), a_time );
         }
         if ( typeid(a_species) == typeid(FluidSpecies&) ) {
            (static_cast<FluidSpecies&>(a_species)).convertFromPhysical();
         }
      }

      virtual void fillGhostCells( FluidSpecies&  a_species_phys,
                                   const double   a_time ) = 0;

      virtual Real computeDtExplicitTI(const FluidSpeciesPtrVect& fluids)
      {
         return DBL_MAX;
      }
   
      virtual Real computeDtImExTI(const FluidSpeciesPtrVect& fluids)
      {
         return DBL_MAX;
      }
   
      virtual Real TimeScale(const FluidSpeciesPtrVect& fluids)
      {
         return DBL_MAX;
      }

      virtual bool trivialSolutionOp() const = 0;

      void parseInitialConditions( const string&  a_species_name,
                                   const string&  a_variable_name )
      {
         GridFunctionLibrary* library = GridFunctionLibrary::getInstance();

         std::string ic_prefix;
         if ( a_species_name.length() == 0 ) {
            ic_prefix = "IC." + a_variable_name;
         }
         else {
            ic_prefix = "IC." + a_species_name + "." + a_variable_name;
         }
         ParmParse ppic( ic_prefix.c_str() );
         std::string function_name;
         ppic.query( "function", function_name );
         RefCountedPtr<GridFunction> ic( library->find( function_name ) );
         m_fluid_ics.insert( FluidSpeciesICMap::value_type( a_variable_name, ic ) );
      }

      const GridFunction& fluidSpeciesIC( const std::string& a_name ) const
      {
         FluidSpeciesICMap::const_iterator it = m_fluid_ics.find( a_name );
         if ( it == m_fluid_ics.end() ) {
            const std::string msg( "GKSystemIC: Initial condition " + a_name + " not found!" );
            MayDay::Error( msg.c_str() );
         }
         return *(it->second);
      }

      typedef std::map<std::string,RefCountedPtr<GridFunction> > FluidSpeciesICMap;
      FluidSpeciesICMap   m_fluid_ics;
};

#include "NamespaceFooter.H"
#endif
