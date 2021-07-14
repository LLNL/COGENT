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
   std::vector<std::string> species_names_vec;
   
   //create species map
   while (more_kinetic_species) {
      std::stringstream s;
      s << "kinetic_species." << count+1;
      ParmParse ppspecies( s.str().c_str() );

      std::string this_species_name("Invalid");
      if (ppspecies.contains("name")) {
         ppspecies.get("name", this_species_name);
         species_names_vec.push_back(this_species_name);
         typedef std::map<std::string,int>::value_type valType;
         m_species_map.insert( valType( this_species_name, count) );
         count++;
      }
      else {
         more_kinetic_species = false;
      }
   }
   m_num_species=count;
   
   //get all collisional models for all species
   count=0;
   while (count<m_num_species) {
      std::stringstream ss;
      // look for all collisional operator for all the kinetic species
      ss << "kinetic_species." << count+1;
      std::string species_name("Invalid");
      ParmParse ppspecies( ss.str().c_str() );
      ppspecies.get("name", species_name);
      std::string cls_type(_CLS_NONE_);
      CLSInterface* cls(NULL);

      for (int count_bkgr(0); count_bkgr<m_num_species; count_bkgr++){
         std::stringstream cls_s;
         if (count == count_bkgr) {
            //input verbatim to sepcify like-species collision model
            //kinetic_species.1.cls = "Krook"
            cls_s << "cls";
         }
         else {
            //input verbatim to specify unlike-species collision model
            //kinetic_species.1.cls_3 = "Krook"
            cls_s << "cls_" << count_bkgr + 1;
         }

         std::string cls_species(cls_s.str());
         std::string species_name_bkgr(species_names_vec[count_bkgr]);
         
         // parse input for collision models; use NullCLS if unspecified
         if (ppspecies.contains( cls_species )) {
            ppspecies.get( cls_species.c_str(), cls_type );
            
            std::string prefix;
            //prefix verbatim to specify cls model parameters for like-collisions
            if (count == count_bkgr) {
               prefix = "CLS."+ species_name;
            }
             //prefix verbatim to specify cls model parameters for unlike-collisions
            else {
               prefix = "CLS."+ species_name+"."+species_name_bkgr;
            }

            if (cls_type == _CLS_KROOK_) {
               cls = new Krook( prefix, m_verbose );
            }
            else if (cls_type == _CLS_LORENTZ_) {
	       cls = new Lorentz( prefix, m_verbose );
            }
            else if (cls_type == _CLS_LINEARIZED_) {
               cls = new Linearized( prefix, m_verbose );
            }
            else if (cls_type == _CLS_LINEARIZEDUNLIKE_) {
               cls = new LinearizedUnlike( prefix, m_verbose );
            }
            else if (cls_type == _CLS_FOKKERPLANCK_) {
               cls = new FokkerPlanck( prefix, m_verbose);
            }
            else if (cls_type == _CLS_CONSDRAGDIFF_) {
               cls = new ConsDragDiff( species_name, prefix, m_verbose );
            }
            else if (cls_type == _CLS_NONE_) {
               cls = new NullCLS();
            }
            else {
               if (procID()==0) {
                  cout << "Unknown collision model specified for " << species_name <<" with "<<species_name_bkgr;
                  cout << ". Using none.\n";
               }
               cls_type = _CLS_NONE_;
               cls = new NullCLS();
            }
         }
         else {
            if (procID()==0) cout << "No collision model specified for " << species_name <<" with "<<species_name_bkgr<< ".\n";
            cls_type = _CLS_NONE_;
            cls = new NullCLS();
         }
         
         m_collision_model.push_back( cls );
         m_collision_model_name.push_back(cls_type);

         if (!procID()) {
            cout << "Collision model for " << species_name <<" with " << species_name_bkgr<< "\t"  << ":\t";
            cout << cls_type << "\n" ;
         }
      }
      
      count++;
   }
}


GKCollisions::~GKCollisions()
{
   for (int i(0); i<m_collision_model.size(); i++ ) {
      delete m_collision_model[i];
   }
}


CLSInterface& GKCollisions::collisionModel( const std::string& a_name,const std::string& a_name_bkgr )

{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_species_map.find( a_name ) );
   CH_assert(it!=m_species_map.end());
   const int index_a((*it).second);

   const mapIterator it_bkgr( m_species_map.find( a_name_bkgr ) );
   CH_assert(it_bkgr!=m_species_map.end());
   const int index_b((*it_bkgr).second);
   return *(m_collision_model[index_a*m_num_species+index_b]);
}


std::string GKCollisions::collisionModelName( const std::string& a_name, const std::string& a_name_bkgr )
{
  typedef std::map<std::string,int>::iterator mapIterator;
  const mapIterator it( m_species_map.find( a_name ) );
  CH_assert(it!=m_species_map.end());
  const int index_a((*it).second);

  const mapIterator it_bkgr( m_species_map.find( a_name_bkgr ) );
  CH_assert(it_bkgr!=m_species_map.end());
  const int index_b((*it_bkgr).second);
   return(m_collision_model_name[index_a*m_num_species+index_b]);
}


void GKCollisions::accumulateRHS( KineticSpeciesPtrVect&       a_rhs,
                                  const KineticSpeciesPtrVect& a_soln,
                                  const bool                   a_implicit,
                                  const Real                   a_time )
{
   for (int species(0); species<a_rhs.size(); species++) {
      KineticSpecies& rhs_species( *(a_rhs[species]) );
      const std::string species_name( rhs_species.name() );
      for(int species_bkgr(0); species_bkgr<a_rhs.size(); species_bkgr++){
        KineticSpecies& rhs_species_bkgr( *(a_rhs[species_bkgr]) );
        const std::string species_name_bkgr( rhs_species_bkgr.name() );
        CLSInterface& CLS( collisionModel( species_name,species_name_bkgr ) );

        if (a_implicit) {
          CLS.evalClsRHSImplicit( a_rhs, a_soln, species, species_bkgr, a_time );
        } else {
          CLS.evalClsRHSExplicit( a_rhs, a_soln, species, species_bkgr, a_time );
        }
     }
   }
}



Real GKCollisions::computeDtExplicitTI( const KineticSpeciesPtrVect& soln )
{
  Real dt(DBL_MAX);
  int count(0);
  for (int species(0);species<m_num_species;species++) {
    for (int species_bkgr(0);species_bkgr<m_num_species;species_bkgr++) {
      Real tmp = m_collision_model[species*m_num_species+species_bkgr]->computeDtExplicitTI(soln,species);
      dt = (tmp < dt ? tmp : dt);
      count++;
    }
  }
  return (count ? dt : -1);
}

Real GKCollisions::computeDtImExTI( const KineticSpeciesPtrVect& soln )
{
  Real dt(DBL_MAX);
  int count(0);
  for (int species(0);species<m_num_species;species++) {
    for (int species_bkgr(0);species_bkgr<m_num_species;species_bkgr++) {
      Real tmp = m_collision_model[species*m_num_species+species_bkgr]->computeDtImExTI(soln,species);
      dt = (tmp < dt ? tmp : dt);
      count++;
    }
  }
  return (count ? dt : -1);
}
Real GKCollisions::computeTimeScale( const KineticSpeciesPtrVect& soln )
{
  Real scale = DBL_MAX;
  int count = 0;
  for (int species(0);species<m_num_species;species++) {
    for (int species_bkgr(0);species_bkgr<m_num_species;species_bkgr++) {
      Real tmp = m_collision_model[species*m_num_species+species_bkgr]->computeTimeScale(soln,species);
      scale = (tmp < scale ? tmp : scale);
      count++;
    }
  }
  return (count ? scale : -1);
}

bool GKCollisions::isLinear()
{
  bool linear_flag = true;
  for (int species(0);species<m_num_species;species++) {
    for (int species_bkgr(0);species_bkgr<m_num_species;species_bkgr++) {
      linear_flag = (linear_flag && m_collision_model[species*m_num_species+species_bkgr]->isLinear());
     }
  }
  return linear_flag;
}


int GKCollisions::precondMatrixBands()
{
  int max_bands = 0;
  for (int species(0);species<m_num_species;species++) {
    for (int species_bkgr(0);species_bkgr<m_num_species;species_bkgr++) {
      max_bands = std::max(max_bands, m_collision_model[species*m_num_species+species_bkgr]->precondMatrixBands());
    }
  }
  return max_bands;
}


void GKCollisions::assemblePrecondMatrix( void *a_P,
                                          const KineticSpeciesPtrVect& a_soln,
                                          const GlobalDOFKineticSpeciesPtrVect& a_gdofs,
                                          const int                                     a_step,
                                          const int                                     a_stage,
                                          const Real a_shift)
{
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies&           soln_species(*(a_soln[species]));
    GlobalDOFKineticSpecies&  gdofs_species(*(a_gdofs[species]));
    const std::string         species_name(soln_species.name());
    for (int species_bkgr(0); species_bkgr<a_soln.size(); species_bkgr++) {
      KineticSpecies&           soln_species_bkgr(*(a_soln[species_bkgr]));
      const std::string         species_name_bkgr(soln_species_bkgr.name());

      CLSInterface&             CLS(collisionModel(species_name,species_name_bkgr));

      CLS.assemblePrecondMatrix(  a_P,
                                  soln_species,
                                  gdofs_species,
                                  a_step,
                                  a_stage,
                                  a_shift );
    }
  }
}

void GKCollisions::defineMultiPhysicsPC(std::vector<Preconditioner<ODEVector,AppCtxt>*>& a_pc,
                                        std::vector<DOFList>&                             a_dof_list,
                                        const KineticSpeciesPtrVect&                      a_soln,
                                        const GlobalDOFKineticSpeciesPtrVect&             a_gdofs,
                                        const ODEVector&                                  a_x,
                                        void*                                             a_ops,
                                        const std::string&                                a_out_string,
                                        const std::string&                                a_opt_string,
                                        bool                                              a_im )
{
  for (int species(0); species<a_soln.size(); species++) {
    const KineticSpecies&           soln_species(*(a_soln[species]));
    const GlobalDOFKineticSpecies&  gdofs_species(*(a_gdofs[species]));
    const std::string         species_name(soln_species.name());
    for (int species_bkgr(0); species_bkgr<a_soln.size(); species_bkgr++) {
      const KineticSpecies&     soln_species_bkgr(*(a_soln[species_bkgr]));
      const std::string         species_name_bkgr(soln_species_bkgr.name());

      CLSInterface&             CLS(collisionModel(species_name,species_name_bkgr));

      CLS.defineBlockPC(a_pc, a_dof_list, a_x, a_ops, a_out_string, a_opt_string, a_im,
                      soln_species, gdofs_species, species);
     }
   }
}

void GKCollisions::updateMultiPhysicsPC(std::vector<Preconditioner<ODEVector,AppCtxt>*>& a_pc,
                                        const KineticSpeciesPtrVect&                      a_soln,
                                        const GlobalDOFKineticSpeciesPtrVect&             a_gdofs,
                                        const Real                                        a_time,
                                        const int                                         a_step,
                                        const int                                         a_stage,
                                        const Real                                        a_shift,
                                        const bool                                        a_im )
{
  for (int species(0); species<a_soln.size(); species++) {
    const KineticSpecies&           soln_species(*(a_soln[species]));
    const GlobalDOFKineticSpecies&  gdofs_species(*(a_gdofs[species]));
    const std::string         species_name(soln_species.name());
    for (int species_bkgr(0); species_bkgr<a_soln.size(); species_bkgr++) {
      const KineticSpecies&     soln_species_bkgr(*(a_soln[species_bkgr]));
      const std::string         species_name_bkgr(soln_species_bkgr.name());

      CLSInterface&             CLS(collisionModel(species_name,species_name_bkgr));

      CLS.updateBlockPC(  a_pc, 
                          soln_species, 
                          gdofs_species, 
                          a_time, 
                          a_step,
                          a_stage,
                          a_shift, 
                          a_im, 
                          species);
     }
   }
}

void GKCollisions::preTimeStep( const KineticSpeciesPtrVect& a_soln,
                                const Real a_time,
                                const KineticSpeciesPtrVect& a_soln_physical )

{
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies&           soln_species(*(a_soln[species]));
    const std::string         species_name(soln_species.name());
    for (int species_bkgr(0); species_bkgr<a_soln.size(); species_bkgr++) {
    KineticSpecies&           soln_species_bkgr(*(a_soln[species_bkgr]));
    const std::string         species_name_bkgr(soln_species_bkgr.name());

    CLSInterface&             CLS(collisionModel(species_name,species_name_bkgr));
    CLS.preTimeStep( a_soln, species, a_time, a_soln_physical );
   }
 }
}

void GKCollisions::postTimeStage( const KineticSpeciesPtrVect& a_soln, const Real a_time, const int a_stage )
{
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies&           soln_species(*(a_soln[species]));
    const std::string         species_name(soln_species.name());
    for (int species_bkgr(0); species_bkgr<a_soln.size(); species_bkgr++) {
    KineticSpecies&           soln_species_bkgr(*(a_soln[species_bkgr]));
    const std::string         species_name_bkgr(soln_species_bkgr.name());

    CLSInterface&             CLS(collisionModel(species_name,species_name_bkgr));
    CLS.postTimeStage( a_soln, species, a_time, a_stage );
   }
  }
}

#include "NamespaceFooter.H"
