#include "SelfConsistentBCOp.H"

#include "NamespaceHeader.H"

const char* SelfConsistentBCOp::pp_name = {"SelfConsistentBC_op"};

SelfConsistentBCOp::SelfConsistentBCOp(const ParmParse&    a_pp,
                                       const PhaseGeom&    a_geometry,
                                       const double        a_larmor,
                                       const int           a_verbosity )
   : m_vlasov(NULL),
     m_verbosity(a_verbosity),
     m_geometry(a_geometry),
     m_larmor(a_larmor),
     m_is_vlasov_implicit(false)
{
}

SelfConsistentBCOp::~SelfConsistentBCOp()
{
   delete m_vlasov;
}

void SelfConsistentBCOp::accumulateExplicitRHS(GKRHSData&                         a_rhs,
                                               const KineticSpeciesPtrVect&       a_kinetic_species,
                                               const CFG::FluidSpeciesPtrVect&    a_fluid_species,
                                               const ScalarPtrVect&               a_scalars,
                                               const CFG::EMFields&               a_EM_fields,
                                               const CFG::PhiOps&                 a_phi_ops,
                                               const int                          a_component,
                                               const Real                         a_time )
{
  //Presently, this us this development only if vlasov is implicit.
  //If vlasov is explicit, BC updated as part of the call to applyVlasovOperator
  
  if (m_vlasov == NULL) {
    int n_species = a_kinetic_species.size();
    std::vector<std::string> name_list(n_species);
    for (int n=0; n<n_species; n++) {
      const KineticSpecies& species = *(a_kinetic_species[n]);
      name_list[n] = species.name();
    }
    
    ParmParse pp_vlasov(GKVlasov::pp_name);
    m_vlasov = new GKVlasovAmpere(pp_vlasov, m_larmor, name_list, true);
    std::vector<std::string> implicit_vlasov_species(0);
    if (m_vlasov) m_vlasov->implicitSpecies(implicit_vlasov_species);
    m_is_vlasov_implicit = (implicit_vlasov_species.size() > 0);
  }
  
  if (m_is_vlasov_implicit) {
    
    if (!m_rhs.isDefined()) {
      m_rhs.define(a_kinetic_species, a_fluid_species, a_scalars);
    }
    
    m_vlasov->accumulateRHS(m_rhs, a_kinetic_species, a_EM_fields, a_phi_ops, false, a_time);
    
    ScalarPtrVect& rhs_scalar = a_rhs.dataScalar();
    Vector<Real>& rhs_scalar_data = rhs_scalar[a_rhs.getScalarComponent("Er_boundary")]->data();
 
    ScalarPtrVect& local_rhs_scalar = m_rhs.dataScalar();
    Vector<Real>& local_rhs_scalar_data = local_rhs_scalar[a_rhs.getScalarComponent("Er_boundary")]->data();
    
    rhs_scalar_data[0] = local_rhs_scalar_data[0];
    rhs_scalar_data[1] = local_rhs_scalar_data[1];
  }
}

void SelfConsistentBCOp::accumulateImplicitRHS(GKRHSData&                         a_rhs,
                                               const KineticSpeciesPtrVect&       a_kinetic_species,
                                               const CFG::FluidSpeciesPtrVect&    a_fluid_species,
                                               const ScalarPtrVect&               a_scalars,
                                               const CFG::EMFields&               a_EM_fields,
                                               const int                          a_component,
                                               const bool                         a_recompute_kinetic_terms,
                                               const Real                         a_time )
{
  //Presently, only explicit integration of consistent BC is implemented. Thus do nothing here.
}
#include "NamespaceFooter.H"

