#include "GKFieldOp.H"
#include "Vorticity.H"

#include <float.h>
#include <sstream>

#include "NamespaceHeader.H"

GKFieldOp::GKFieldOp( const int a_verbose )
   : m_verbose(a_verbose)
{
   bool more_fields(true);
   int count(0);
   while (more_fields) {

      // look for data specifying another kinetic species
      std::stringstream s;
      s << "field." << count+1;
      ParmParse ppfield( s.str().c_str() );

      std::string field_name("Invalid");
      if (ppfield.contains("name")) {
         ppfield.get("name", field_name);
         
         std::string model_type(_FIELD_NONE_);
         FieldOpInterface* model(NULL);
         
         if (ppfield.contains( "model" )) {
            ppfield.get( "model", model_type );
            const std::string prefix( "field_model." + field_name );
            ParmParse ppfield_model( prefix.c_str() );
            
            if (model_type == _FIELD_NONE_) {
               model = new NullFieldOp();
            }

            else if (model_type == _VORTICITY_) {
               model = new Vorticity( ppfield_model, m_verbose );
            }

            else {
               if (procID()==0) { 
                 cout << "Unknown field model specified for " << field_name;
                 cout << ". Using none.\n";
               }
               model_type = _FIELD_NONE_;
               model = new NullFieldOp();
            }
         }
         else {
            if (procID()==0) cout << "No field model specified for " << field_name << ".\n";
            model = new NullFieldOp();
         }
         
         m_field_model.push_back( model );
         typedef std::map<std::string,int>::value_type valType;
         m_fieldVecComp_map.insert( valType( field_name, count ) );
         m_field_model_name.insert( valType( model_type, count ) );
         if (!procID()) {
           cout << "Field model for " << count << "\t" << field_name << ":\t";
           cout << model_type << "\n" ;
         }
         count++;
      }
      else {
         more_fields = false;
      }
   }
}


GKFieldOp::~GKFieldOp()
{
   for (int i(0); i<m_field_model.size(); i++ ) {
      delete m_field_model[i];
   }
}


FieldOpInterface& GKFieldOp::fieldModel( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_fieldVecComp_map.find( a_name ) );
   CH_assert(it!=m_fieldVecComp_map.end());
   const int index((*it).second);
   return *(m_field_model[index]);
}


std::string GKFieldOp::fieldModelName( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_fieldVecComp_map.find( a_name ) );
   CH_assert(it!=m_fieldVecComp_map.end());
   const int index((*it).second);

   for (std::map<std::string,int>::iterator it=m_field_model_name.begin(); it!=m_field_model_name.end(); it++) {
     if (it->second == index) return it->first;
   }
   return(_FIELD_NONE_);
}


void GKFieldOp::accumulateRHS( FieldPtrVect&                      a_rhs,
                               const FieldPtrVect&                a_fields,
                               const FluidSpeciesPtrVect&         a_fluids,
                               const PS::KineticSpeciesPtrVect&   a_kinetic_species,
                               const LevelData<FluxBox>&          a_E_field,
                               const Real                         a_time)
                              
{
   for (int component(0); component<a_rhs.size(); component++) {
      Field& rhs_component( *(a_rhs[component]) );
      const std::string component_name( rhs_component.name() );
      FieldOpInterface& fieldOp( fieldModel( component_name ) );
      fieldOp.evalFieldRHS( a_rhs, a_fields, a_fluids, a_kinetic_species, a_E_field, component, a_time );
   }
}


Real GKFieldOp::computeDt( const FieldPtrVect&        fields,
                           const FluidSpeciesPtrVect& fluids)
{
  Real dt(DBL_MAX);

  std::map<std::string,int>::iterator it;
  for (it=m_field_model_name.begin(); it!=m_field_model_name.end(); ++it) {
    Real tmp = m_field_model[it->second]->computeDt(fields, fluids);
    dt = (tmp < dt ? tmp : dt);
  }
  return dt;
}

Real GKFieldOp::computeTimeScale( const FieldPtrVect&        fields,
                                  const FluidSpeciesPtrVect& fluids)
{
  std::map<std::string,int>::iterator it;
  Real scale = DBL_MAX;
  int count = 0;
  for (it=m_field_model_name.begin(); it!=m_field_model_name.end(); ++it) {
    Real tmp = m_field_model[it->second]->TimeScale(fields, fluids);
    scale = (tmp < scale ? tmp : scale);
    count++;
  }
  return (count ? scale : -1);
}


#include "NamespaceFooter.H"
