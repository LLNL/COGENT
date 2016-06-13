#include "GKFluidOp.H"

#include <float.h>
#include <sstream>

#include "NamespaceHeader.H"

GKFluidOp::GKFluidOp( const int a_verbose )
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
         
         std::string model_type(_FLUID_NONE_);
         FluidOpInterface* model(NULL);
         
         if (ppspecies.contains( "model" )) {
            ppspecies.get( "model", model_type );
            const std::string prefix( "fluid_model." + species_name );
            ParmParse ppfluid( prefix.c_str() );
            
            if (model_type == _FLUID_NONE_) {
               model = new NullFluidOp();
            }

            //else if (model_type == _ONE_FIELD_) {
            //   model = new OneField( ppfluid, m_verbose );
            //}

            else {
               if (procID()==0) { 
                 cout << "Unknown collision model specified for " << species_name;
                 cout << ". Using none.\n";
               }
               model_type = _FLUID_NONE_;
               model = new NullFluidOp();
            }
         }
         else {
            if (procID()==0) cout << "No fluid model specified for " << species_name << ".\n";
            model = new NullFluidOp();
         }
         
         m_fluid_model.push_back( model );
         typedef std::map<std::string,int>::value_type valType;
         m_species_map.insert( valType( species_name, count ) );
         m_fluid_model_name.insert( valType( model_type, count ) );
         if (!procID()) {
           cout << "Fluid model for " << count << "\t" << species_name << ":\t";
           cout << model_type << "\n" ;
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
   return(_FLUID_NONE_);
}


void GKFluidOp::accumulateRHS( FluidSpeciesPtrVect&               a_rhs,
                               const FieldPtrVect&                a_fields,
                               const FluidSpeciesPtrVect&         a_fluids,
                               const PS::KineticSpeciesPtrVect&   a_kinetic_species,
                               const LevelData<FluxBox>&          a_E_field,
                               const Real                         a_time)
                              
{
   for (int species(0); species<a_rhs.size(); species++) {
      FluidSpecies& rhs_species( *(a_rhs[species]) );
      const std::string species_name( rhs_species.name() );
      FluidOpInterface& fluidOp( fluidModel( species_name ) );
      fluidOp.evalFluidRHS( a_rhs, a_fields, a_fluids, a_kinetic_species, a_E_field, species, a_time );
   }
}


Real GKFluidOp::computeDt( const FieldPtrVect&        fields,
                           const FluidSpeciesPtrVect& fluids)
{
  Real dt(DBL_MAX);

  std::map<std::string,int>::iterator it;
  for (it=m_fluid_model_name.begin(); it!=m_fluid_model_name.end(); ++it) {
    Real tmp = m_fluid_model[it->second]->computeDt(fields, fluids);
    dt = (tmp < dt ? tmp : dt);
  }
  return dt;
}

Real GKFluidOp::computeTimeScale( const FieldPtrVect&        fields,
                                  const FluidSpeciesPtrVect& fluids )
{
  std::map<std::string,int>::iterator it;
  Real scale = DBL_MAX;
  int count = 0;
  for (it=m_fluid_model_name.begin(); it!=m_fluid_model_name.end(); ++it) {
    Real tmp = m_fluid_model[it->second]->TimeScale(fields, fluids);
    scale = (tmp < scale ? tmp : scale);
    count++;
  }
  return (count ? scale : -1);
}


#include "NamespaceFooter.H"
