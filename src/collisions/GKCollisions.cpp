#include "GKCollisions.H"

#include <float.h>
#include <sstream>

#include "NamespaceHeader.H"

GKCollisions::GKCollisions( const int a_verbose )
   : m_verbose(a_verbose)
{
   bool more_kinetic_species(true);
   int count(0);
   while (more_kinetic_species) {

      // look for data specifying another kinetic species
      std::stringstream s;
      s << "kinetic_species." << count+1;
      ParmParse ppspecies( s.str().c_str() );

      std::string species_name("Invalid");
      if (ppspecies.contains("name")) {
         ppspecies.get("name", species_name);
         
         std::string cls_type(_CLS_NONE_);
         CLSInterface* cls(NULL);
         
         if (ppspecies.contains( "cls" )) {
            ppspecies.get( "cls", cls_type );
            const std::string prefix( "CLS." + species_name );
            ParmParse ppcls( prefix.c_str() );
            
            if (cls_type == _CLS_KROOK_) {
               cls = new Krook( ppcls, m_verbose );
            }
            else if (cls_type == _CLS_MYKROOK_) {
               cls = new MyKrook( ppcls, m_verbose );
            }
            else if (cls_type == _CLS_LORENTZ_) {
               cls = new Lorentz( ppcls, m_verbose );
            }
            else if (cls_type == _CLS_LINEARIZED_) {
               cls = new Linearized( ppcls, m_verbose );
            }
            else if (cls_type == _CLS_FOKKERPLANCK_) {
               cls = new FokkerPlanck( ppcls, m_verbose );
            }
            else if (cls_type == _CLS_CONSDRAGDIFF_) {
               cls = new ConsDragDiff( species_name, ppcls, m_verbose );
            }
            else if (cls_type == _CLS_NONE_) {
               cls = new NullCLS();
            }
            else {
               if (procID()==0) { 
                 cout << "Unknown collision model specified for " << species_name;
                 cout << ". Using none.\n";
               }
               cls_type = _CLS_NONE_;
               cls = new NullCLS();
            }
         }
         else {
            if (procID()==0) cout << "No collision model specified for " << species_name << ".\n";
            cls = new NullCLS();
         }
         
         m_collision_model.push_back( cls );
         typedef std::map<std::string,int>::value_type valType;
         m_species_map.insert( valType( species_name, count ) );
         m_collision_model_name.insert( valType( cls_type, count ) );
         if (!procID()) {
           cout << "Collision model for " << count << "\t" << species_name << ":\t";
           cout << cls_type << "\n" ;
         }
         count++;
      }
      else {
         more_kinetic_species = false;
      }
   }
}


GKCollisions::~GKCollisions()
{
   for (int i(0); i<m_collision_model.size(); i++ ) {
      delete m_collision_model[i];
   }
}


CLSInterface& GKCollisions::collisionModel( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_species_map.find( a_name ) );
   CH_assert(it!=m_species_map.end());
   const int index((*it).second);
   return *(m_collision_model[index]);
}


std::string GKCollisions::collisionModelName( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_species_map.find( a_name ) );
   CH_assert(it!=m_species_map.end());
   const int index((*it).second);

   for (std::map<std::string,int>::iterator it=m_collision_model_name.begin(); it!=m_collision_model_name.end(); it++) {
     if (it->second == index) return it->first;
   }
   return(_CLS_NONE_);
}


void GKCollisions::accumulateRHS( KineticSpeciesPtrVect&       a_rhs,
                                  const KineticSpeciesPtrVect& a_soln,
                                  const Real                   a_time,
                                  const int                    a_flag )
{
   for (int species(0); species<a_rhs.size(); species++) {
      KineticSpecies& rhs_species( *(a_rhs[species]) );
      const std::string species_name( rhs_species.name() );
      CLSInterface& CLS( collisionModel( species_name ) );
      CLS.evalClsRHS( a_rhs, a_soln, species, a_time, a_flag );
   }
}


Real GKCollisions::computeDt( const KineticSpeciesPtrVect& soln )
{
  Real dt(DBL_MAX);
//  int count(0);
  std::map<std::string,int>::iterator it;
  for (it=m_collision_model_name.begin(); it!=m_collision_model_name.end(); ++it) {
    Real tmp = m_collision_model[it->second]->computeDt(soln);
    dt = (tmp < dt ? tmp : dt);
//    count++;
  }
//  return (count ? dt : -1);
  return dt;
}

Real GKCollisions::computeTimeScale( const KineticSpeciesPtrVect& soln )
{
  std::map<std::string,int>::iterator it;
  Real scale = DBL_MAX;
  int count = 0;
  for (it=m_collision_model_name.begin(); it!=m_collision_model_name.end(); ++it) {
    Real tmp = m_collision_model[it->second]->TimeScale(soln);
    scale = (tmp < scale ? tmp : scale);
    count++;
  }
  return (count ? scale : -1);
}

bool GKCollisions::isLinear()
{
  bool linear_flag = true;
  std::map<std::string,int>::iterator it;
  for (it=m_collision_model_name.begin(); it!=m_collision_model_name.end(); ++it) {
    int index = getCollisionModelIndex( it );
    linear_flag = (linear_flag && m_collision_model[index]->isLinear());
  }
  return linear_flag;
}

//bool GKCollisions::setupPrecondMatrix(void *a_P,int a_N,int a_bs)
bool GKCollisions::setupPrecondMatrix( void *a_P, int a_N )
{
   int count = 0;
   bool flag = true;
   if (!procID()) {
      cout << "Setting up preconditioner matrix for collision model(s)... ";
   }
   std::map<std::string,int>::iterator it;
   for (it=m_collision_model_name.begin(); it!=m_collision_model_name.end(); ++it) {
      int index = getCollisionModelIndex( it );
//     bool done = m_collision_model[it->second]->setupPrecondMatrix(a_P,a_N,a_bs);
      bool done = m_collision_model[index]->setupPrecondMatrix( a_P, a_N );
      flag = (flag && done );
      count++;
  }
  CH_assert(count<=1);
  if (!procID()) {
    cout << "Done.\n";
  }
  return flag;
}

void GKCollisions::assemblePrecondMatrix( void *a_P,
                                          const KineticSpeciesPtrVect& a_soln,
                                          const Mapping& a_mapping)
{
  CH_assert(a_soln.size() == 1);
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies& soln_species(*(a_soln[species]));
    const std::string species_name(soln_species.name());
    CLSInterface& CLS(collisionModel(species_name));
    CLS.assemblePrecondMatrix(a_P,a_soln,species,a_mapping);
  }
}

void GKCollisions::preTimeStep( const KineticSpeciesPtrVect& a_soln, const Real a_time )
{
   for (int species(0); species<a_soln.size(); species++) {
      KineticSpecies& soln_species( *(a_soln[species]) );
      const std::string species_name( soln_species.name() );
      CLSInterface& CLS( collisionModel( species_name ) );
      CLS.preTimeStep( a_soln, species, a_time );
   }
}

void GKCollisions::postTimeStage( const KineticSpeciesPtrVect& a_soln, const Real a_time, const int a_stage )
{
   for (int species(0); species<a_soln.size(); species++) {
      KineticSpecies& soln_species( *(a_soln[species]) );
      const std::string species_name( soln_species.name() );
      CLSInterface& CLS( collisionModel( species_name ) );
      CLS.postTimeStage( a_soln, species, a_time, a_stage );
   }
}

#include "NamespaceFooter.H"
