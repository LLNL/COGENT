#include "GKFluidOp.H"
#include "OneFieldOp.H"
#include "NullFluidOp.H"
#include "VorticityOp.H"

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
            ParmParse ppfluid( prefix.c_str() );
            
            if (op_type == "NullFluidOp") {
               model = new NullFluidOp();
            }
            else if (op_type == "VorticityOp") {
               model = new VorticityOp( ppfluid, a_geometry, a_larmor, m_verbose );
            }
            else if (op_type == "AmpereErAverageOp") {
               // We don't need a fluid operator in this case because the
               // Vlasov operator is doing all the work
               model = new NullFluidOp();
            }
            else if (op_type == "OneFieldOp") {
               model = new OneFieldOp( ppfluid, a_geometry, m_verbose );
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
                               const bool                         a_recompute_kinetic_terms,
                               const Real                         a_time)
                              
{
   for (int species(0); species<a_rhs.size(); species++) {
      FluidSpecies& rhs_species( static_cast<FluidSpecies&>(*(a_rhs[species])) );
      const std::string species_name( rhs_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      if ( a_implicit ) {
         fluidOp.accumulateImplicitRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, species, a_recompute_kinetic_terms, a_time );
      }
      else {
         fluidOp.accumulateExplicitRHS( a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, species, a_time );
      }
   }
}


void GKFluidOp::convertToPhysical( FluidSpeciesPtrVect&        a_fluid_species_phys,
                                   const FluidSpeciesPtrVect&  a_fluid_species_comp,
                                   const IntVect&              a_ghosts ) const
{
   int num_species = a_fluid_species_comp.size();

   a_fluid_species_phys.resize(num_species);
   
   for (int species(0); species<num_species; species++) {
      a_fluid_species_phys[species] = a_fluid_species_comp[species]->convertToPhysical(a_ghosts);
   }
}


void GKFluidOp::evalSolutionOp( FluidSpeciesPtrVect&               a_rhs,
                                const PS::KineticSpeciesPtrVect&   a_kinetic_species,
                                const FluidSpeciesPtrVect&         a_fluid_species,
                                const PS::ScalarPtrVect&           a_scalars,
                                const Real                         a_time )
{
   for (int species(0); species<a_rhs.size(); species++) {
      FluidSpecies& rhs_species( static_cast<FluidSpecies&>(*(a_rhs[species])) );
      const std::string species_name( rhs_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.evalSolutionOp( a_rhs, a_kinetic_species, a_fluid_species, a_scalars, species, a_time );
   }
}


void GKFluidOp::solveSolutionPC( FluidSpeciesPtrVect&              a_fluid_species_solution,
                                 const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                                 const FluidSpeciesPtrVect&        a_fluid_species_rhs )
{
   for (int species(0); species<a_fluid_species_solution.size(); species++) {
      FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species_solution[species])) );
      const std::string species_name( fluid_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.solveSolutionPC( a_fluid_species_solution, a_kinetic_species_rhs, a_fluid_species_rhs, species );
   }
}


void GKFluidOp::updatePC( const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                          const FluidSpeciesPtrVect&        a_fluid_species,
                          const double                      a_mshift )
{
   for (int species(0); species<a_fluid_species.size(); species++) {
      FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species[species])) );
      const std::string species_name( fluid_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.updatePCImEx( a_kinetic_species, a_mshift );
   }
}


void GKFluidOp::solvePCImEx( FluidSpeciesPtrVect&              a_fluid_species_solution,
                             const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                             const FluidSpeciesPtrVect&        a_fluid_species_rhs )
{
   for (int species(0); species<a_fluid_species_solution.size(); species++) {
      FluidSpecies& fluid_species( static_cast<FluidSpecies&>(*(a_fluid_species_solution[species])) );
      const std::string species_name( fluid_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.solvePCImEx( a_fluid_species_solution, a_kinetic_species_rhs, a_fluid_species_rhs, species );
   }
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



Real GKFluidOp::computeDt(const FluidSpeciesPtrVect& fluids)
{
  Real dt(DBL_MAX);

  std::map<std::string,int>::iterator it;
  for (it=m_fluid_model_name.begin(); it!=m_fluid_model_name.end(); ++it) {
    Real tmp = m_fluid_model[it->second]->computeDt(fluids);
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
