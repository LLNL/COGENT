#include "GKFluidOp.H"
#include "OneFieldOp.H"
#include "BurgersOp.H"
#include "IdealMhdOp.H"
#include "FullFluidOp.H"
#include "TwoFieldNeutralsOp.H"
#include "NullFluidOp.H"
#include "VorticityOp.H"
#include "AmpereErAverageOp.H"

#include <float.h>
#include <sstream>

#include "NamespaceHeader.H"

GKFluidOp::GKFluidOp( const MagGeom&  a_geometry,
                      const double    a_larmor,
                      const int       a_verbose )
   : m_verbose(a_verbose)
{
   bool more_fluid_species(true);
   int count(0);
   while (more_fluid_species) {

      // look for data specifying another kinetic species
      std::stringstream s;
      s << "fluid_species." << count+1;
      ParmParse ppspecies( s.str().c_str() );

      std::string species_name("Invalid");
      if (ppspecies.contains("name")) {
         ppspecies.get("name", species_name);
         
         std::string op_type;
         FluidOpInterface* model(NULL);
         
         if (ppspecies.contains( "operator_type" )) {
            ppspecies.get( "operator_type", op_type );
            const std::string prefix( "fluid_model." + species_name );
            
            if (op_type == "NullFluidOp") {
               model = new NullFluidOp();
            }
            else if (op_type == "VorticityOp") {
               model = new VorticityOp( prefix, a_geometry, a_larmor, m_verbose );
            }
            else if (op_type == "AmpereErAverageOp") {
               model = new AmpereErAverageOp();
            }
            else if (op_type == "OneFieldOp") {
               model = new OneFieldOp( prefix, species_name, a_geometry, m_verbose );
            }
            else if (op_type == "BurgersOp") {
               model = new BurgersOp( prefix, species_name, a_geometry, m_verbose );
            }
            else if (op_type == "IdealMhdOp") {
               model = new IdealMhdOp( prefix, species_name, a_geometry, m_verbose );
            }
            else if (op_type == "FullFluidOp") {
               model = new FullFluidOp( prefix, species_name, a_geometry, m_verbose );
            }
            else if (op_type == "TwoFieldNeutralsOp") {
               model = new TwoFieldNeutralsOp( prefix, species_name, a_geometry, m_verbose );
            }
            else {
               MayDay::Error("Unknown fluid operator type specified for a fluid species");
            }
         }
         else {
            MayDay::Error("No fluid operator type specified for a fluid species");
         }
         
         m_fluid_model.push_back( model );
         typedef std::map<std::string,int>::value_type valType;
         m_species_map.insert( valType( species_name, count ) );
         m_fluid_model_name.insert( valType( op_type, count ) );
         if (!procID()) {
           cout << "Fluid model for " << count << "\t" << species_name << ":\t";
           cout << op_type << "\n" ;
         }
         count++;
      }
      else {
         more_fluid_species = false;
      }
   }
}


GKFluidOp::~GKFluidOp()
{
   for (int i(0); i<m_fluid_model.size(); i++ ) {
      delete m_fluid_model[i];
   }
}


FluidOpInterface& GKFluidOp::fluidModel( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_species_map.find( a_name ) );
   CH_assert(it!=m_species_map.end());
   const int index((*it).second);
   return *(m_fluid_model[index]);
}


std::string GKFluidOp::fluidModelName( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_species_map.find( a_name ) );
   CH_assert(it!=m_species_map.end());
   const int index((*it).second);

   for (std::map<std::string,int>::iterator it=m_fluid_model_name.begin(); it!=m_fluid_model_name.end(); it++) {
     if (it->second == index) return it->first;
   }
   return("NullFluidOp");
}


void GKFluidOp::accumulateRHS( FluidSpeciesPtrVect&               a_rhs,
                               const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                               const FluidSpeciesPtrVect&         a_fluid_species,
                               const PS::ScalarPtrVect&           a_scalars,
                               const EField&                      a_E_field,
                               const bool                         a_implicit,
                               const Real                         a_time)
                              
{
   for (int species(0); species<a_rhs.size(); species++) {
      FluidSpecies& rhs_species( static_cast<FluidSpecies&>(*(a_rhs[species])) );
      const std::string species_name( rhs_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      if ( a_implicit ) {
         fluidOp.accumulateImplicitRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, species, a_time );
      }
      else {
         fluidOp.accumulateExplicitRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, species, a_time );
      }
   }
}


void GKFluidOp::preSolutionOpEval( const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                   const FluidSpeciesPtrVect&         a_fluid_species,
                                   const PS::ScalarPtrVect&           a_scalars,
                                   const EField&                      a_E_field,
                                   const double                       a_time )
{
   for (int species(0); species<a_fluid_species.size(); species++) {
      FluidSpecies& rhs_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
      const std::string species_name( rhs_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.preSolutionOpEval( a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, a_time );
   }
}


void GKFluidOp::preOpEval( const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                           const FluidSpeciesPtrVect&         a_fluid_species,
                           const PS::ScalarPtrVect&           a_scalars,
                           const EField&                      a_E_field,
                           const double                       a_time )
{
   for (int species(0); species<a_fluid_species.size(); species++) {
      FluidSpecies& rhs_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
      const std::string species_name( rhs_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.preOpEval( a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, a_time );
   }
}


void GKFluidOp::convertToPhysical( FluidSpeciesPtrVect&  a_fluid_species_phys,
                             const FluidSpeciesPtrVect&  a_fluid_species_comp ) const
{
   CH_TIME("GKFluidOp::convertToPhysical");
   const int num_species = a_fluid_species_phys.size();
   CH_assert( num_species == a_fluid_species_comp.size());

   for (int s=0; s<num_species; ++s) {
      a_fluid_species_phys[s]->copy(*a_fluid_species_comp[s]);
      a_fluid_species_phys[s]->convertToPhysical();
   }

}


void GKFluidOp::evalSolutionOp( FluidSpeciesPtrVect&               a_rhs,
                                const PS::KineticSpeciesPtrVect&   a_kinetic_species_comp,
                                const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                const FluidSpeciesPtrVect&         a_fluid_species_comp,
                                const FluidSpeciesPtrVect&         a_fluid_species_phys,
                                const PS::ScalarPtrVect&           a_scalars,
                                const Real                         a_time )
{
   for (int species(0); species<a_rhs.size(); species++) {
      FluidSpecies& rhs_species( static_cast<FluidSpecies&>(*(a_rhs[species])) );
      const std::string species_name( rhs_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.evalSolutionOp( a_rhs, 
                              a_kinetic_species_comp, 
                              a_kinetic_species_phys, 
                              a_fluid_species_comp, 
                              a_fluid_species_phys, 
                              a_scalars, 
                              species, 
                              a_time );
   }
}


void GKFluidOp::solveSolutionPC( FluidSpeciesPtrVect&              a_fluid_species_solution,
                                 const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                 const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                                 int                               a_idx )
{
   if (a_idx < 0) {
      for (int species(0); species<a_fluid_species_solution.size(); species++) {
         FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species_solution[species])) );
         const std::string species_name( fluid_species.name() );
         FluidOpInterface& fluidOp( fluidModel( species_name ) );
         fluidOp.solveSolutionPC( a_fluid_species_solution, a_kinetic_species_rhs, a_fluid_species_rhs, species );
      }
   } else {
      FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species_solution[a_idx])) );
      const std::string species_name( fluid_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.solveSolutionPC( a_fluid_species_solution, a_kinetic_species_rhs, a_fluid_species_rhs, a_idx );
   }
}


void GKFluidOp::updatePC( const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                          const FluidSpeciesPtrVect&        a_fluid_species,
                          const double                      a_time,
                          const double                      a_shift,
                          int                               a_idx )
{
   if (a_idx < 0) {
      for (int species(0); species<a_fluid_species.size(); species++) {
         FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
         const std::string species_name( fluid_species.name() );
         FluidOpInterface& fluidOp( fluidModel( species_name ) );
         fluidOp.updatePCImEx( a_fluid_species, a_kinetic_species, a_time, a_shift, a_idx );
      }
   } else {
      FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species[a_idx])) );
      const std::string species_name( fluid_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.updatePCImEx( a_fluid_species, a_kinetic_species, a_time, a_shift, a_idx );
   }
}


void GKFluidOp::solvePCImEx( FluidSpeciesPtrVect&              a_fluid_species_solution,
                             const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                             const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                             int                               a_idx )
{
   if (a_idx < 0 ) {
      for (int species(0); species<a_fluid_species_solution.size(); species++) {
         FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species_solution[species])) );
         const std::string species_name( fluid_species.name() );
         FluidOpInterface& fluidOp( fluidModel( species_name ) );
         fluidOp.solvePCImEx( a_fluid_species_solution, a_kinetic_species_rhs, a_fluid_species_rhs, species );
      }
   } else {
      FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species_solution[a_idx])) );
      const std::string species_name( fluid_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.solvePCImEx( a_fluid_species_solution, a_kinetic_species_rhs, a_fluid_species_rhs, a_idx );
   }
}

void GKFluidOp::defineMultiPhysicsPC( std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
                                      std::vector<PS::DOFList>&                                 a_dof_list,
                                      const FluidSpeciesPtrVect&                                a_fluid_species,
                                      const PS::GlobalDOFFluidSpeciesPtrVect&                   a_global_dofs,
                                      const PS::GKVector&                                       a_soln_vec,
                                      PS::GKOps&                                                a_gkops,
                                      const std::string&                                        a_out_string,
                                      const std::string&                                        a_opt_string,
                                      bool                                                      a_im )
{
   for (int species(0); species<a_fluid_species.size(); species++) {
      FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
      const PS::GlobalDOFFluidSpecies& gdofs_species(*(a_global_dofs[species]));
      const std::string species_name( fluid_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );

      fluidOp.defineBlockPC(a_pc, a_dof_list, a_soln_vec, a_gkops, a_out_string, a_opt_string, a_im,
                            fluid_species, gdofs_species, species );
   }
   return;
}

void GKFluidOp::updateMultiPhysicsPC( std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
                                      const PS::KineticSpeciesPtrVect&                          a_kin_species_phys,
                                      const FluidSpeciesPtrVect&                                a_fluid_species,
                                      const Real                                                a_time,
                                      const Real                                                a_shift,
                                      const bool                                                a_im )
{
   for (int species(0); species<a_fluid_species.size(); species++) {
      FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
      const std::string species_name( fluid_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );

      fluidOp.updateBlockPC(a_pc, a_kin_species_phys, a_fluid_species, a_time, a_shift, a_im, species );
   }
   return;
}

bool GKFluidOp::trivialSolutionOp( const FluidSpeciesPtrVect& a_fluid_species )
{
   bool trivial = true;

   for (int species(0); species<a_fluid_species.size(); species++) {
      const FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
      const std::string species_name( fluid_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      trivial &= fluidOp.trivialSolutionOp();
   }

   return trivial;
}



void GKFluidOp::initialize( FluidSpeciesPtrVect&  a_fluid_species,
                            const Real            a_time )
{
   CH_TIME("GKFluidOp::initialize()");
   for (int species(0); species<a_fluid_species.size(); species++) {
      //      FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
      CFGVars& fluid_species = *(a_fluid_species[species]);
      const std::string species_name( fluid_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.initialize(fluid_species, a_time);
   }
}



void GKFluidOp::fillGhostCells( FluidSpeciesPtrVect&  a_fluid_species,
                                const double          a_time )
{
   for (int species(0); species<a_fluid_species.size(); species++) {
      FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
      const std::string species_name( fluid_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.fillGhostCells(fluid_species, a_time);
   }
}



Real GKFluidOp::computeDtExplicitTI(const FluidSpeciesPtrVect& fluids)
{
  Real dt(DBL_MAX);

  std::map<std::string,int>::iterator it;
  for (it=m_fluid_model_name.begin(); it!=m_fluid_model_name.end(); ++it) {
    Real tmp = m_fluid_model[it->second]->computeDtExplicitTI(fluids);
    dt = (tmp < dt ? tmp : dt);
  }
  return dt;
}

Real GKFluidOp::computeDtImExTI(const FluidSpeciesPtrVect& fluids)
{
  Real dt(DBL_MAX);

  std::map<std::string,int>::iterator it;
  for (it=m_fluid_model_name.begin(); it!=m_fluid_model_name.end(); ++it) {
    Real tmp = m_fluid_model[it->second]->computeDtImExTI(fluids);
    dt = (tmp < dt ? tmp : dt);
  }
  return dt;
}

Real GKFluidOp::computeTimeScale( const FluidSpeciesPtrVect& fluids )
{
  std::map<std::string,int>::iterator it;
  Real scale = DBL_MAX;
  int count = 0;
  for (it=m_fluid_model_name.begin(); it!=m_fluid_model_name.end(); ++it) {
    Real tmp = m_fluid_model[it->second]->TimeScale(fluids);
    scale = (tmp < scale ? tmp : scale);
    count++;
  }
  return (count ? scale : -1);
}


#include "NamespaceFooter.H"
