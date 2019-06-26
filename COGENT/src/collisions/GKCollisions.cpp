#include "GKCollisions.H"

#include <float.h>
#include <sstream>
#include <algorithm>

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
            
            if (cls_type == _CLS_KROOK_) {
               cls = new Krook( prefix, m_verbose );
            }
            else if (cls_type == _CLS_MYKROOK_) {
               cls = new MyKrook( prefix, m_verbose );
            }
            else if (cls_type == _CLS_LORENTZ_) {
               cls = new Lorentz( prefix, m_verbose );
            }
            else if (cls_type == _CLS_LINEARIZED_) {
               cls = new Linearized( prefix, m_verbose );
            }
            else if (cls_type == _CLS_FOKKERPLANCK_) {
               cls = new FokkerPlanck( prefix, m_verbose );
            }
            else if (cls_type == _CLS_CONSDRAGDIFF_) {
               cls = new ConsDragDiff( species_name, prefix, m_verbose );
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
                                  const bool                   a_implicit,
                                  const Real                   a_time )
{
   for (int species(0); species<a_rhs.size(); species++) {
      KineticSpecies& rhs_species( *(a_rhs[species]) );
      const std::string species_name( rhs_species.name() );
      CLSInterface& CLS( collisionModel( species_name ) );
      if (a_implicit) {
        CLS.evalClsRHSImplicit( a_rhs, a_soln, species, a_time );
      } else {
        CLS.evalClsRHSExplicit( a_rhs, a_soln, species, a_time );
      }
   }
}

Real GKCollisions::computeDtExplicitTI( const KineticSpeciesPtrVect& soln )
{
  Real dt(DBL_MAX);
  int count(0);
  std::map<std::string,int>::iterator it;
  for (it=m_collision_model_name.begin(); it!=m_collision_model_name.end(); ++it) {
    Real tmp = m_collision_model[it->second]->computeDtExplicitTI(soln,it->second);
    dt = (tmp < dt ? tmp : dt);
    count++;
  }
  return (count ? dt : -1);
}

Real GKCollisions::computeDtImExTI( const KineticSpeciesPtrVect& soln )
{
  Real dt(DBL_MAX);
  int count(0);
  std::map<std::string,int>::iterator it;
  for (it=m_collision_model_name.begin(); it!=m_collision_model_name.end(); ++it) {
    Real tmp = m_collision_model[it->second]->computeDtImExTI(soln,it->second);
    dt = (tmp < dt ? tmp : dt);
    count++;
  }
  return (count ? dt : -1);
}

Real GKCollisions::computeTimeScale( const KineticSpeciesPtrVect& soln )
{
  std::map<std::string,int>::iterator it;
  Real scale = DBL_MAX;
  int count = 0;
  for (it=m_collision_model_name.begin(); it!=m_collision_model_name.end(); ++it) {
    Real tmp = m_collision_model[it->second]->computeTimeScale(soln,it->second);
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

int GKCollisions::precondMatrixBands()
{
  int max_bands = 0;
  std::map<std::string,int>::iterator it;
  for (it=m_collision_model_name.begin(); it!=m_collision_model_name.end(); ++it) {
    int index = getCollisionModelIndex( it );
    max_bands = std::max(max_bands, m_collision_model[index]->precondMatrixBands());
  }
  return max_bands;
}

void GKCollisions::assemblePrecondMatrix( void *a_P,
                                          const KineticSpeciesPtrVect& a_soln,
                                          const GlobalDOFKineticSpeciesPtrVect& a_gdofs,
                                          const Real a_shift)
{
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies&           soln_species(*(a_soln[species]));
    GlobalDOFKineticSpecies&  gdofs_species(*(a_gdofs[species]));
    const std::string         species_name(soln_species.name());
    CLSInterface&             CLS(collisionModel(species_name));

    CLS.assemblePrecondMatrix(a_P,soln_species,gdofs_species,a_shift);
  }
}

void GKCollisions::defineMultiPhysicsPC(std::vector<Preconditioner<GKVector,GKOps>*>& a_pc,
                                        std::vector<DOFList>&                         a_dof_list,
                                        const KineticSpeciesPtrVect&                  a_soln,
                                        const GlobalDOFKineticSpeciesPtrVect&         a_gdofs,
                                        const GKVector&                               a_x,
                                        GKOps&                                        a_ops,
                                        const std::string&                            a_out_string,
                                        const std::string&                            a_opt_string,
                                        bool                                          a_im )
{
  for (int species(0); species<a_soln.size(); species++) {
    const KineticSpecies&           soln_species(*(a_soln[species]));
    const GlobalDOFKineticSpecies&  gdofs_species(*(a_gdofs[species]));
    const std::string               species_name(soln_species.name());
    CLSInterface&                   CLS(collisionModel(species_name));

    CLS.defineBlockPC(a_pc, a_dof_list, a_x, a_ops, a_out_string, a_opt_string, a_im, 
                      soln_species, gdofs_species, species);
  }
}

void GKCollisions::updateMultiPhysicsPC(std::vector<Preconditioner<GKVector,GKOps>*>& a_pc,
                                        const KineticSpeciesPtrVect&                  a_soln,
                                        const GlobalDOFKineticSpeciesPtrVect&         a_gdofs,
                                        const Real                                    a_time,
                                        const Real                                    a_shift,
                                        const bool                                    a_im )
{
  for (int species(0); species<a_soln.size(); species++) {
    const KineticSpecies&           soln_species(*(a_soln[species]));
    const GlobalDOFKineticSpecies&  gdofs_species(*(a_gdofs[species]));
    const std::string               species_name(soln_species.name());
    CLSInterface&                   CLS(collisionModel(species_name));

    CLS.updateBlockPC(a_pc, soln_species, gdofs_species, a_time, a_shift, a_im, species);
  }
}

void GKCollisions::preTimeStep( const KineticSpeciesPtrVect& a_soln, 
                                const Real a_time,
                                const KineticSpeciesPtrVect& a_soln_physical )

{
   for (int species(0); species<a_soln.size(); species++) {
      KineticSpecies& soln_species( *(a_soln[species]) );
      const std::string species_name( soln_species.name() );
      CLSInterface& CLS( collisionModel( species_name ) );
      CLS.preTimeStep( a_soln, species, a_time, a_soln_physical );
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
