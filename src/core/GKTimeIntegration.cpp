#include "GKTimeIntegration.H"
#include "inspect.H"

#include "GKUtils.H.multidim"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GKUtils.H.multidim"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

////////////////////////////////////////////////////////////////////

inline
void copyData( KineticSpeciesPtrVect& a_kinetic_species_dst,
               const KineticSpeciesPtrVect& a_kinetic_species_src,
               CFG::FluidSpeciesPtrVect& a_fluid_species_dst,
               const CFG::FluidSpeciesPtrVect& a_fluid_species_src,
               CFG::FieldPtrVect& a_field_dst,
               const CFG::FieldPtrVect& a_field_src )
{
   a_kinetic_species_dst.resize( a_kinetic_species_src.size() );
   for (int s(0); s<a_kinetic_species_dst.size(); s++) {
      a_kinetic_species_dst[s]->copy( *(a_kinetic_species_src[s]) );
   }
   a_fluid_species_dst.resize( a_fluid_species_src.size() );
   for (int s(0); s<a_fluid_species_dst.size(); s++) {
      a_fluid_species_dst[s]->copy( *(a_fluid_species_src[s]) );
   }
   a_field_dst.resize( a_field_src.size() );
   for (int s(0); s<a_field_dst.size(); s++) {
      a_field_dst[s]->copy( *(a_field_src[s]) );
   }
}
   
void GKState::copy( const GKState& a_state )
{
   copyData( m_kinetic_species, a_state.dataKinetic(),
             m_fluid_species, a_state.dataFluid(),
             m_fields, a_state.dataField() );
}

void GKState::copy( const GKRHSData& a_rhs )
{
   copyData( m_kinetic_species, a_rhs.dataKinetic(),
             m_fluid_species, a_rhs.dataFluid(),
             m_fields, a_rhs.dataField() );
}

void GKRHSData::copy( const GKRHSData& a_rhs )
{
   copyData( m_kinetic_species, a_rhs.dataKinetic(),
             m_fluid_species, a_rhs.dataFluid(),
             m_fields, a_rhs.dataField() );
}

void GKRHSData::copy( const GKState& a_state )
{
   copyData( m_kinetic_species, a_state.dataKinetic(),
             m_fluid_species, a_state.dataFluid(),
             m_fields, a_state.dataField() );
}

////////////////////////////////////////////////////////////////////

inline
int computeVectorSize( const KineticSpeciesPtrVect& a_kinetic_species,
                       const CFG::FluidSpeciesPtrVect& a_fluid_species,
                       const CFG::FieldPtrVect& a_fields )
{
   int total_size(0);
   for (int s(0); s<a_kinetic_species.size(); s++) {
      total_size += a_kinetic_species[s]->size();
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      total_size += a_fluid_species[s]->size();
   }
   for (int s(0); s<a_fields.size(); s++) {
      total_size += a_fields[s]->size();
   }
   return total_size;
}

int GKState::getVectorSize()
{
   CH_assert( isDefined() );
   return computeVectorSize( m_kinetic_species, m_fluid_species, m_fields );
}

int GKRHSData::getVectorSize()
{
   CH_assert( isDefined() );
   return computeVectorSize( m_kinetic_species, m_fluid_species, m_fields );
}

////////////////////////////////////////////////////////////////////

inline
void copyToArray( Real* a_dst,
                  const KineticSpeciesPtrVect& a_kinetic_species,
                  const CFG::FluidSpeciesPtrVect& a_fluid_species,
                  const CFG::FieldPtrVect& a_fields )
{
   int offset(0);
   for (int s(0); s<a_kinetic_species.size(); s++) {
      offset += GKUtils::copyFromLevelData(
         a_dst + offset,
         a_kinetic_species[s]->distributionFunction() );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      offset += CFG::GKUtils::copyFromLevelData( a_dst + offset,
                                                 a_fluid_species[s]->data() );
   }
   for (int s(0); s<a_fields.size(); s++) {
      offset += CFG::GKUtils::copyFromLevelData( a_dst + offset,
                                                 a_fields[s]->data() );
   }
}

void GKState::copyTo( Real* a_vector ) const
{
   CH_assert( isDefined() );
   copyToArray( a_vector, m_kinetic_species, m_fluid_species, m_fields );
}

void GKRHSData::copyTo( Real* a_vector ) const 
{
   CH_assert( isDefined() );
   copyToArray( a_vector, m_kinetic_species, m_fluid_species, m_fields );
}

////////////////////////////////////////////////////////////////////

inline
void copyFromArray( KineticSpeciesPtrVect& a_kinetic_species,
                    CFG::FluidSpeciesPtrVect& a_fluid_species,
                    CFG::FieldPtrVect& a_fields,
                    Real* a_src )
{
   int offset(0);
   for (int s(0); s<a_kinetic_species.size(); s++) {
      offset += GKUtils::copyToLevelData(
         a_kinetic_species[s]->distributionFunction(),
         a_src + offset );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      offset += CFG::GKUtils::copyToLevelData( a_fluid_species[s]->data(),
                                               a_src + offset );
   }
   for (int s(0); s<a_fields.size(); s++) {
      offset += CFG::GKUtils::copyToLevelData( a_fields[s]->data(),
                                               a_src + offset );
   }
}

void GKState::copyFrom( Real* a_vector )
{
   CH_assert( isDefined() );
   copyFromArray( m_kinetic_species, m_fluid_species, m_fields, a_vector );
}

void GKRHSData::copyFrom( Real* a_vector )
{
   CH_assert( isDefined() );
   copyFromArray( m_kinetic_species, m_fluid_species, m_fields, a_vector );
}

////////////////////////////////////////////////////////////////////

inline
void addFromArray( KineticSpeciesPtrVect& a_kinetic_species,
                   CFG::FluidSpeciesPtrVect& a_fluid_species,
                   CFG::FieldPtrVect& a_fields,
                   Real* a_src,
                   const Real& a_factor )
{
   int offset(0);
   for (int s(0); s<a_kinetic_species.size(); s++) {
      offset += GKUtils::addToLevelData(
         a_kinetic_species[s]->distributionFunction(),
         a_src + offset,
         a_factor );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      offset += CFG::GKUtils::addToLevelData( a_fluid_species[s]->data(),
                                              a_src + offset,
                                              a_factor );
   }
   for (int s(0); s<a_fields.size(); s++) {
      offset += CFG::GKUtils::addToLevelData( a_fields[s]->data(),
                                              a_src + offset,
                                              a_factor );
   }
}

void GKState::addFrom( Real* a_vector, const Real& a_factor )
{
   CH_assert( isDefined() );
   addFromArray( m_kinetic_species, m_fluid_species, m_fields, a_vector, a_factor );
}

////////////////////////////////////////////////////////////////////

inline
void scaleData( KineticSpeciesPtrVect& a_kinetic_species,
                CFG::FluidSpeciesPtrVect& a_fluid_species,
                CFG::FieldPtrVect& a_fields,
                const Real& a_factor )
{
   for (int s(0); s<a_kinetic_species.size(); s++) {
      GKUtils::scaleLevelData( a_kinetic_species[s]->distributionFunction(), a_factor );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      CFG::GKUtils::scaleLevelData( a_fluid_species[s]->data(), a_factor );
   }
   for (int s(0); s<a_fields.size(); s++) {
      CFG::GKUtils::scaleLevelData( a_fields[s]->data(), a_factor );
   }
}

void GKState::scale( const Real& a_factor )
{
   CH_assert( isDefined() );
   scaleData( m_kinetic_species, m_fluid_species, m_fields, a_factor );
}

void GKRHSData::scale( const Real& a_factor )
{
   CH_assert( isDefined() );
   scaleData( m_kinetic_species, m_fluid_species, m_fields, a_factor );
}

////////////////////////////////////////////////////////////////////

inline
void incrementData( KineticSpeciesPtrVect& a_kinetic_species,
                    const KineticSpeciesPtrVect& a_dkinetic_species,
                    CFG::FluidSpeciesPtrVect& a_fluid_species,
                    const CFG::FluidSpeciesPtrVect& a_dfluid_species,
                    CFG::FieldPtrVect& a_fields,
                    const CFG::FieldPtrVect& a_dfields,
                    const Real& a_factor )
{
   CH_assert( a_kinetic_species.size()==a_dkinetic_species.size() );
   CH_assert( a_fluid_species.size()==a_dfluid_species.size() );
   CH_assert( a_fields.size()==a_dfields.size() );
   for (int s(0); s<a_kinetic_species.size(); s++) {
//      GKUtils::incrementLevelData( *(a_kinetic_species[s]),
//                                   *(a_dkinetic_species[s]),
//                                   a_factor );
      a_kinetic_species[s]->addData( *(a_dkinetic_species[s]), a_factor );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
//      CFG::GKUtils::incrementLevelData( *(a_fluid_species[s]),
//                                        *(a_dfluid_species[s]),
//                                        a_factor );
      a_fluid_species[s]->addData( *(a_dfluid_species[s]), a_factor );
   }
   for (int s(0); s<a_fields.size(); s++) {
//      CFG::GKUtils::incrementLevelData( *(a_field[s]),
//                                        *(a_dfield[s]),
//                                        a_factor );
      a_fields[s]->addData( *(a_dfields[s]), a_factor );
   }
}

void GKState::increment( const GKState& a_delta,
                         const Real& a_factor,
                         const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_kinetic_species, a_delta.dataKinetic(),
                  m_fluid_species, a_delta.dataFluid(),
                  m_fields, a_delta.dataField(),
                  a_factor );
}

void GKState::increment( const GKRHSData& a_delta,
                         const Real& a_factor,
                         const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_kinetic_species, a_delta.dataKinetic(),
                  m_fluid_species, a_delta.dataFluid(),
                  m_fields, a_delta.dataField(),
                  a_factor );
}

void GKRHSData::increment( const GKRHSData& a_delta,
                           const Real& a_factor,
                           const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_kinetic_species, a_delta.dataKinetic(),
                  m_fluid_species, a_delta.dataFluid(),
                  m_fields, a_delta.dataField(),
                  a_factor );
}

void GKRHSData::increment( const GKState& a_delta,
                           const Real& a_factor,
                           const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_kinetic_species, a_delta.dataKinetic(),
                  m_fluid_species, a_delta.dataFluid(),
                  m_fields, a_delta.dataField(),
                  a_factor );
}

////////////////////////////////////////////////////////////////////

Real GKRHSData::dotProduct( const GKRHSData& a_vector )
{
   CH_assert( isDefined() );
   double sum_local(0);

   const KineticSpeciesPtrVect& b_kinetic_species( a_vector.dataKinetic() );
   CH_assert( m_kinetic_species.size()==b_kinetic_species.size() );
   for (int s(0); s<m_kinetic_species.size(); s++) {
      sum_local += GKUtils::innerProductLevelData(
         m_kinetic_species[s]->distributionFunction(),
         b_kinetic_species[s]->distributionFunction() );
   }
   
   const  CFG::FluidSpeciesPtrVect& b_fluid_species( a_vector.dataFluid() );
   CH_assert( m_fluid_species.size()==b_fluid_species.size() );
   for (int s(0); s<m_fluid_species.size(); s++) {
      sum_local += CFG::GKUtils::innerProductLevelData(
         m_fluid_species[s]->data(),
         b_fluid_species[s]->data() );
   }
   
   const  CFG::FieldPtrVect& b_fields( a_vector.dataField() );
   CH_assert( m_fields.size()==b_fields.size() );
   for (int s(0); s<m_fields.size(); s++) {
      sum_local += CFG::GKUtils::innerProductLevelData(
         m_fields[s]->data(),
         b_fields[s]->data() );
   }
   
   double sum(0.0);

#ifdef CH_MPI
   MPI_Allreduce( &sum_local, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
   sum = sum_local;
#endif

   return sum;
}
   
////////////////////////////////////////////////////////////////////
   
inline
Real normLp( const KineticSpeciesPtrVect& a_kinetic_species,
             const CFG::FluidSpeciesPtrVect& a_fluid_species,
             const CFG::FieldPtrVect& a_fields,
             const int a_p )
{
   double accum(0);
   for (int s(0); s<a_kinetic_species.size(); s++) {
      accum += GKUtils::accumulateLevelData( a_kinetic_species[s]->distributionFunction(), a_p );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      accum += CFG::GKUtils::accumulateLevelData( a_fluid_species[s]->data(), a_p );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      accum += CFG::GKUtils::accumulateLevelData( a_fields[s]->data(), a_p );
   }
   return GKUtils::reduce( accum, a_p );
}

   
Real GKState::computeNorm( const int a_p ) const
{
   CH_assert(isDefined()); 
   CH_assert((a_p <= 2) && (a_p >= 0));
   return normLp( m_kinetic_species, m_fluid_species, m_fields, a_p );
}

Real GKRHSData::computeNorm( const int a_p ) const
{
   CH_assert(isDefined());
   CH_assert((a_p <= 2) && (a_p >= 0));
   return normLp( m_kinetic_species, m_fluid_species, m_fields, a_p );
}

#include "NamespaceFooter.H"
