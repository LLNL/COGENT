#include "GKScalarOp.H"
#include "NullScalarOp.H"
#include "SelfConsistentBCOp.H"

#include "NamespaceHeader.H"

GKScalarOp::GKScalarOp(const PhaseGeom& a_geometry,
                       const GKUnits&   a_units,
                       const int        a_verbose )
   : m_verbose(a_verbose)
{
   bool more_scalars(true);
   int count(0);
   while (more_scalars) {

      // look for data specifying another scalar species
      std::stringstream s;
      s << "scalar." << count+1;
      ParmParse ppscalar( s.str().c_str() );

      std::string scalar_name("Invalid");
      if (ppscalar.contains("name")) {
         ppscalar.get("name", scalar_name);
         
         std::string op_type;
         ScalarOpInterface* op(NULL);
         
         if (ppscalar.contains( "operator_type" )) {
            ppscalar.get( "operator_type", op_type );
            const std::string prefix( "scalar_model." + scalar_name );
            
            if (op_type == "NullScalarOp") {
               op = new NullScalarOp();
            }
            else if (op_type == "SelfConsistentBCOp") {
               const double larmor( a_units.larmorNumber() );
               op = new SelfConsistentBCOp(prefix, a_geometry, larmor, m_verbose);
            }
            else {
               MayDay::Error("GKScalarOp::GKScalarOp(): Unrecognized scalar operator");
            }
         }
         else {
            MayDay::Error("GKScalarOp::GKScalarOp: Unspecified scalar operator");
         }
         
         m_scalar_op.push_back( op );
         typedef std::map<std::string,int>::value_type valType;
         m_scalarVecComp_map.insert( valType( scalar_name, count ) );
         m_scalar_op_name.insert( valType( op_type, count ) );
         if (!procID()) {
           cout << "Scalar operator for " << count << "\t" << scalar_name << ":\t";
           cout << op_type << "\n" ;
         }
         count++;
      }
      else {
         more_scalars = false;
      }
   }
}


GKScalarOp::~GKScalarOp()
{
   for (int i(0); i<m_scalar_op.size(); i++ ) {
      delete m_scalar_op[i];
   }
}


ScalarOpInterface& GKScalarOp::scalarOp( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_scalarVecComp_map.find( a_name ) );
   CH_assert(it!=m_scalarVecComp_map.end());
   const int index((*it).second);
   return *(m_scalar_op[index]);
}


std::string GKScalarOp::scalarOpName( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_scalarVecComp_map.find( a_name ) );
   CH_assert(it!=m_scalarVecComp_map.end());
   const int index((*it).second);

   for (std::map<std::string,int>::iterator it=m_scalar_op_name.begin(); it!=m_scalar_op_name.end(); it++) {
     if (it->second == index) return it->first;
   }
   return("NullScalarOp");
}


void GKScalarOp::accumulateRHS( GKRHSData&                         a_rhs,
                                const KineticSpeciesPtrVect&       a_kinetic_species,
                                const CFG::FluidSpeciesPtrVect&    a_fluid_species,
                                const ScalarPtrVect&               a_scalars,
                                const CFG::EField&                 a_E_field,
                                const bool                         a_implicit,
                                const bool                         a_recompute_kinetic_terms,
                                const Real                         a_time)
                              
{
   for (int component(0); component<a_scalars.size(); component++) {
      Scalar& rhs_component( *(a_scalars[component]) );
      const std::string component_name( rhs_component.name() );
      ScalarOpInterface& scalarOperator( scalarOp( component_name ) );
      if ( a_implicit ) {
         scalarOperator.accumulateImplicitRHS( a_rhs, a_kinetic_species, a_fluid_species, a_scalars,
                                               a_E_field, component, a_recompute_kinetic_terms, a_time );
      }
      else {
         scalarOperator.accumulateExplicitRHS( a_rhs, a_kinetic_species, a_fluid_species, a_scalars,
                                               a_E_field, component, a_time );
      }
   }
}


Real GKScalarOp::computeDtExplicitTI(const ScalarPtrVect& a_scalars)
{
  Real dt(DBL_MAX);

  std::map<std::string,int>::iterator it;
  for (it=m_scalar_op_name.begin(); it!=m_scalar_op_name.end(); ++it) {
    Real tmp = m_scalar_op[it->second]->computeDtExplicitTI(a_scalars);
    dt = (tmp < dt ? tmp : dt);
  }
  return dt;
}

Real GKScalarOp::computeDtImExTI(const ScalarPtrVect& a_scalars)
{
  Real dt(DBL_MAX);

  std::map<std::string,int>::iterator it;
  for (it=m_scalar_op_name.begin(); it!=m_scalar_op_name.end(); ++it) {
    Real tmp = m_scalar_op[it->second]->computeDtImExTI(a_scalars);
    dt = (tmp < dt ? tmp : dt);
  }
  return dt;
}

Real GKScalarOp::computeTimeScale(const ScalarPtrVect& a_scalars)
{
  std::map<std::string,int>::iterator it;
  Real scale = DBL_MAX;
  int count = 0;
  for (it=m_scalar_op_name.begin(); it!=m_scalar_op_name.end(); ++it) {
    Real tmp = m_scalar_op[it->second]->TimeScale(a_scalars);
    scale = (tmp < scale ? tmp : scale);
    count++;
  }
  return (count ? scale : -1);
}


#include "NamespaceFooter.H"
