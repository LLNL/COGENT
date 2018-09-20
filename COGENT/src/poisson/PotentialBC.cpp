#include "PotentialBC.H"

#include "NamespaceHeader.H"


PotentialBC::PotentialBC(const int a_num_boundaries)
{
   m_bc_function.resize(a_num_boundaries);
   m_bdry_name.resize(a_num_boundaries);
   m_bc_type.resize(a_num_boundaries);
   m_bc_value.resize(a_num_boundaries);

   for (int i=0; i<a_num_boundaries; ++i) {
      m_bc_type[i] = UNDEFINED;
      m_bc_value[i] = 0.;
      CH_assert(m_bc_function[i].isNull());
   }
}


PotentialBC::~PotentialBC()
{
   m_bc_function.resize(0);
   m_bdry_name.resize(0);
   m_bc_type.resize(0);
   m_bc_value.resize(0);
}


bool
PotentialBC::hasNeumannCondition() const
{
   int has = false;

   for (int i=0; i<m_bc_type.size(); ++i) {
      if (m_bc_type[i] == NEUMANN) {
         has = true;
         break;
      }
   }

   return has;
}




#include "NamespaceFooter.H"
