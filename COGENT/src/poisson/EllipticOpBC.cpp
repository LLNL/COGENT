#include "EllipticOpBC.H"

#include "NamespaceHeader.H"


EllipticOpBC::EllipticOpBC(const int a_num_boundaries,
                           const int a_num_blocks)
   : m_num_boundaries(a_num_boundaries),
     m_num_blocks(a_num_blocks)
{
   m_bc_function.resize(a_num_boundaries);
   m_bdry_name.resize(a_num_boundaries);
   m_bc_type.resize(a_num_boundaries);
   m_bc_value.resize(a_num_boundaries);
   m_bc_subtype.resize(a_num_boundaries);

   for (int i=0; i<a_num_boundaries; ++i) {
      CH_assert(m_bc_function[i].isNull());
      m_bdry_name[i] = "";
      m_bc_type[i] = UNDEFINED;
      m_bc_value[i] = 0.;
      m_bc_subtype[i] = "";
   }
   
   m_bc_block_data.resize(2*a_num_blocks*SpaceDim);
   for (int i=0; i<m_bc_block_data.size(); ++i) {
      CH_assert(m_bc_block_data[i].isNull());
   }

   m_bc_neu_nat_factor.resize(2*a_num_blocks*SpaceDim);
   for (int i=0; i<m_bc_neu_nat_factor.size(); ++i) {
      CH_assert(m_bc_neu_nat_factor[i].isNull());
   }
}


EllipticOpBC::~EllipticOpBC()
{
   m_bc_function.resize(0);
   m_bdry_name.resize(0);
   m_bc_type.resize(0);
   m_bc_value.resize(0);
   m_bc_block_data.resize(0);
   m_bc_neu_nat_factor.resize(0);
   m_bc_subtype.resize(0);
}


bool
EllipticOpBC::hasNeumannCondition() const
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


bool
EllipticOpBC::hasNaturalCondition() const
{
   int has = false;

   for (int i=0; i<m_bc_type.size(); ++i) {
      if (m_bc_type[i] == NATURAL) {
         has = true;
         break;
      }
   }

   return has;
}


bool
EllipticOpBC::hasCoupledBoundary() const
{
   int has = false;

   for (int i=0; i<m_bc_subtype.size(); ++i) {
      if (m_bc_subtype[i] == "coupled") {
         has = true;
         break;
      }
   }

   return has;
}


void
EllipticOpBC::copyBaseData( const EllipticOpBC& a_data )
{
   for (int i=0; i<m_bc_function.size(); ++i) {
      m_bc_function[i] = a_data.m_bc_function[i];
   }
   for (int i=0; i<m_bc_block_data.size(); ++i) {
      m_bc_block_data[i] = a_data.m_bc_block_data[i];
   }
   for (int i=0; i<m_bc_neu_nat_factor.size(); ++i) {
      m_bc_neu_nat_factor[i] = a_data.m_bc_neu_nat_factor[i];
   }
   for (int i=0; i<m_bc_type.size(); ++i) {
      m_bc_type[i] = a_data.m_bc_type[i];
   }
   for (int i=0; i<m_bc_subtype.size(); ++i) {
      m_bc_subtype[i] = a_data.m_bc_subtype[i];
   }
   for (int i=0; i<m_bc_value.size(); ++i) {
      m_bc_value[i] = a_data.m_bc_value[i];
   }
}


void
EllipticOpBC::setExtrapolatedType()
{
   for (int i=0; i<m_bc_type.size(); ++i) {
      m_bc_type[i] = EXTRAPOLATED;
   }
}


#include "NamespaceFooter.H"
