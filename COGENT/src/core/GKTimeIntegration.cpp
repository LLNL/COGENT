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
int computeVectorSize( const KineticSpeciesPtrVect&     a_kinetic_species,
                       const CFG::FluidSpeciesPtrVect&  a_fluid_species,
                       const ScalarPtrVect&             a_scalars)
{
   int total_size(0);
   for (int s(0); s<a_kinetic_species.size(); s++) {
      total_size += a_kinetic_species[s]->size();
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      total_size += a_fluid_species[s]->size();
   }
   for (int s(0); s<a_scalars.size(); s++) {
      total_size += a_scalars[s]->size();
   }
   return total_size;
}

int GKState::getVectorSize() const
{
   CH_assert( isDefined() );
   return computeVectorSize( m_kinetic_species, m_fluid_species, m_scalars );
}

int GKRHSData::getVectorSize() const
{
   CH_assert( isDefined() );
   return computeVectorSize( m_kinetic_species, m_fluid_species, m_scalars );
}

////////////////////////////////////////////////////////////////////

inline
void copyToArray( Real* a_dst,
                  const KineticSpeciesPtrVect&     a_kinetic_species,
                  const CFG::FluidSpeciesPtrVect&  a_fluid_species,
                  const ScalarPtrVect&             a_scalars )
{
   int offset(0);
   for (int s(0); s<a_kinetic_species.size(); s++) {
      offset += GKUtils::copyFromLevelData(
         a_dst + offset,
         a_kinetic_species[s]->distributionFunction() );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      //      offset += CFG::GKUtils::copyFromLevelData( a_dst + offset,
      //                                                 a_fluid_species[s]->data() );
      offset += CFG::GKUtils::copyFromLevelData( a_dst + offset,
                                                 a_fluid_species[s]->cell_data() );
   }
   for (int s(0); s<a_scalars.size(); s++) {
      Vector<Real>& data = a_scalars[s]->data();
      for (int i=0; i<data.size(); ++i) {
         a_dst[offset++] = data[i];
      }
   }
}

void GKState::copyTo( Real* a_vector ) const
{
   CH_assert( isDefined() );
   copyToArray( a_vector, m_kinetic_species, m_fluid_species, m_scalars );
}

void GKRHSData::copyTo( Real* a_vector ) const 
{
   CH_assert( isDefined() );
   copyToArray( a_vector, m_kinetic_species, m_fluid_species, m_scalars );
}

////////////////////////////////////////////////////////////////////

inline
void copyFromArray( KineticSpeciesPtrVect&     a_kinetic_species,
                    CFG::FluidSpeciesPtrVect&  a_fluid_species,
                    ScalarPtrVect&             a_scalars,
                    const Real* a_src )
{
   int offset(0);
   for (int s(0); s<a_kinetic_species.size(); s++) {
      offset += GKUtils::copyToLevelData(
         a_kinetic_species[s]->distributionFunction(),
         a_src + offset );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      //      offset += CFG::GKUtils::copyToLevelData( a_fluid_species[s]->data(),
      //                                               a_src + offset );
      offset += CFG::GKUtils::copyToLevelData( a_fluid_species[s]->cell_data(),
                                               a_src + offset );
   }
   for (int s(0); s<a_scalars.size(); s++) {
      Vector<Real>& data = a_scalars[s]->data();
      for (int i=0; i<data.size(); ++i) {
         data[i] = a_src[offset++];
      }
   }
}

void GKState::copyFrom( const Real* a_vector )
{
   CH_assert( isDefined() );
   copyFromArray( m_kinetic_species, m_fluid_species, m_scalars, a_vector );
}

void GKRHSData::copyFrom( const Real* a_vector )
{
   CH_assert( isDefined() );
   copyFromArray( m_kinetic_species, m_fluid_species, m_scalars, a_vector );
}

////////////////////////////////////////////////////////////////////

inline
void addFromArray( KineticSpeciesPtrVect&     a_kinetic_species,
                   CFG::FluidSpeciesPtrVect&  a_fluid_species,
                   ScalarPtrVect&             a_scalars,
                   Real*                      a_src,
                   const Real&                a_factor )
{
   int offset(0);
   for (int s(0); s<a_kinetic_species.size(); s++) {
      offset += GKUtils::addToLevelData(
         a_kinetic_species[s]->distributionFunction(),
         a_src + offset,
         a_factor );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      //      offset += CFG::GKUtils::addToLevelData( a_fluid_species[s]->data(),
      //                                              a_src + offset,
      //                                              a_factor );
      offset += CFG::GKUtils::addToLevelData( a_fluid_species[s]->cell_data(),
                                              a_src + offset,
                                              a_factor );
   }
   for (int s(0); s<a_scalars.size(); s++) {
      Vector<Real>& data = a_scalars[s]->data();
      for (int i=0; i<data.size(); ++i) {
         data[i] += a_src[offset++];
      }
   }
}

void GKState::addFrom( Real* a_vector, const Real& a_factor )
{
   CH_assert( isDefined() );
   addFromArray( m_kinetic_species, m_fluid_species, m_scalars, a_vector, a_factor );
}

////////////////////////////////////////////////////////////////////

inline
void scaleData( KineticSpeciesPtrVect&     a_kinetic_species,
                CFG::FluidSpeciesPtrVect&  a_fluid_species,
                ScalarPtrVect&             a_scalars,
                const Real&                a_factor )
{
   for (int s(0); s<a_kinetic_species.size(); s++) {
      GKUtils::scaleLevelData( a_kinetic_species[s]->distributionFunction(), a_factor );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      //      CFG::GKUtils::scaleLevelData( a_fluid_species[s]->data(), a_factor );
      CFG::GKUtils::scaleLevelData( a_fluid_species[s]->cell_data(), a_factor );
   }
   for (int s(0); s<a_scalars.size(); s++) {
      Vector<Real>& data = a_scalars[s]->data();
      for (int i=0; i<data.size(); ++i) {
         data[i] *= a_factor;
      }
   }
}

void GKState::scale( const Real& a_factor )
{
   CH_assert( isDefined() );
   scaleData( m_kinetic_species, m_fluid_species, m_scalars, a_factor );
}

void GKRHSData::scale( const Real& a_factor )
{
   CH_assert( isDefined() );
   scaleData( m_kinetic_species, m_fluid_species, m_scalars, a_factor );
}

////////////////////////////////////////////////////////////////////

inline
void incrementData( KineticSpeciesPtrVect&            a_kinetic_species,
                    const KineticSpeciesPtrVect&      a_dkinetic_species,
                    CFG::FluidSpeciesPtrVect&         a_fluid_species,
                    const CFG::FluidSpeciesPtrVect&   a_dfluid_species,
                    ScalarPtrVect&                    a_scalars,
                    const ScalarPtrVect&              a_dscalars,
                    const Real&                       a_factor )
{
   CH_assert( a_kinetic_species.size()==a_dkinetic_species.size() );
   CH_assert( a_fluid_species.size()==a_dfluid_species.size() );
   CH_assert( a_scalars.size()==a_dscalars.size() );
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
   for (int s(0); s<a_scalars.size(); s++) {
      const Vector<Real>& ddata = a_dscalars[s]->data();
      Vector<Real>& data = a_scalars[s]->data();
      for (int i=0; i<data.size(); ++i) {
         data[i] += a_factor * ddata[i];
      }
   }
}

void GKState::increment( const GKState& a_delta,
                         const Real& a_factor,
                         const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_kinetic_species, a_delta.dataKinetic(),
                  m_fluid_species, a_delta.dataFluid(),
                  m_scalars, a_delta.dataScalar(),
                  a_factor );
}

void GKState::increment( const GKRHSData& a_delta,
                         const Real& a_factor,
                         const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_kinetic_species, a_delta.dataKinetic(),
                  m_fluid_species, a_delta.dataFluid(),
                  m_scalars, a_delta.dataScalar(),
                  a_factor );
}

void GKRHSData::increment( const GKRHSData& a_delta,
                           const Real& a_factor,
                           const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_kinetic_species, a_delta.dataKinetic(),
                  m_fluid_species, a_delta.dataFluid(),
                  m_scalars, a_delta.dataScalar(),
                  a_factor );
}

void GKRHSData::increment( const GKState& a_delta,
                           const Real& a_factor,
                           const bool a_update_flux_register )
{
   CH_assert( isDefined() );
   incrementData( m_kinetic_species, a_delta.dataKinetic(),
                  m_fluid_species, a_delta.dataFluid(),
                  m_scalars, a_delta.dataScalar(),
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
                                                       //         m_fluid_species[s]->data(),
                                                       //         b_fluid_species[s]->data() );
         m_fluid_species[s]->cell_data(),
         b_fluid_species[s]->cell_data() );
   }
   
   const ScalarPtrVect& b_scalars( a_vector.dataScalar() );
   CH_assert( m_scalars.size()==b_scalars.size() );
   for (int s(0); s<m_scalars.size(); s++) {
      Vector<Real>& mdata = m_scalars[s]->data();
      Vector<Real>& bdata = b_scalars[s]->data();
      for (int i=0; i<mdata.size(); ++i) {
         sum_local += mdata[i] * bdata[i];
      }
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
Real normLp( const KineticSpeciesPtrVect&      a_kinetic_species,
             const CFG::FluidSpeciesPtrVect&   a_fluid_species,
             const ScalarPtrVect&              a_scalars,
             const int                         a_p )
{
   double accum(0);
   for (int s(0); s<a_kinetic_species.size(); s++) {
      accum += GKUtils::accumulateLevelData( a_kinetic_species[s]->distributionFunction(), a_p );
   }
   for (int s(0); s<a_fluid_species.size(); s++) {
      //      accum += CFG::GKUtils::accumulateLevelData( a_fluid_species[s]->data(), a_p );
      accum += CFG::GKUtils::accumulateLevelData( a_fluid_species[s]->cell_data(), a_p );
   }
   for (int s(0); s<a_scalars.size(); s++) {
      const Vector<Real>& data = a_scalars[s]->data();
      if (a_p == 0) {
         for (int i=0; i<data.size(); ++i) {
            if (fabs(data[i]) > accum) accum = fabs(data[i]);
         }
      }
      else {
         for (int i=0; i<data.size(); ++i) {
            accum += pow(data[i],a_p);
         }
      }
   }
   return GKUtils::reduce( accum, a_p );
}

   
Real GKState::computeNorm( const int a_p ) const
{
   CH_assert(isDefined()); 
   CH_assert((a_p <= 2) && (a_p >= 0));
   return normLp( m_kinetic_species, m_fluid_species, m_scalars, a_p );
}

Real GKRHSData::computeNorm( const int a_p ) const
{
   CH_assert(isDefined());
   CH_assert((a_p <= 2) && (a_p >= 0));
   return normLp( m_kinetic_species, m_fluid_species, m_scalars, a_p );
}

////////////////////////////////////////////////////////////////////

static inline 
const Real& getVal(
                    int                             a_idx,
                    GlobalDOF                       *a_global_dof,
                    const KineticSpeciesPtrVect&    a_kinetic_species,
                    const CFG::FluidSpeciesPtrVect& a_fluid_species
                  )
{
  if (a_global_dof->savedAddress(a_idx)) {

    const GKAddress& address = a_global_dof->getAddress(a_idx);
    if (address.family() == "kinetic") {
      
      int i_kinetic_species = address.number();
      const GlobalDOFKineticSpeciesPtrVect& gdof = a_global_dof->dataKinetic();
      return( gdof[i_kinetic_species]->getVal(*a_kinetic_species[i_kinetic_species], address) );

    } else if (address.family() == "fluid" ) {

      int i_fluid_species = address.number();
      const GlobalDOFFluidSpeciesPtrVect& gdof = a_global_dof->dataFluid();
      ///      return( gdof[i_fluid_species]->getVal(*a_fluid_species[i_fluid_species], address) );
      return( gdof[i_fluid_species]->getVal(static_cast<const CFG::FluidSpecies&>(*a_fluid_species[i_fluid_species]), address) );
    } else {
      std::cout << "Error in GKState/GKRHSData::operator[]: Invalid value for "
                << "\"family\" in saved address.\n";
      exit(1);
    }

  } else {

    /* a_idx is local, so add MPI offset */
    int gidx = a_idx + a_global_dof->mpiOffset();

    int i_kinetic_species = a_global_dof->isIndexKineticSpecies (gidx,a_kinetic_species);
    int i_fluid_species   = a_global_dof->isIndexFluidSpecies   (gidx,a_fluid_species  );
    CH_assert(    ( (i_kinetic_species != -1) && (i_fluid_species == -1) )
              ||  ( (i_kinetic_species == -1) && (i_fluid_species != -1) )
              ||  ( (i_kinetic_species == -1) && (i_fluid_species == -1) ) );

    GKAddress new_address;

    if (i_kinetic_species != -1) {

      new_address.family("kinetic");
      new_address.number(i_kinetic_species);

      const GlobalDOFKineticSpeciesPtrVect& gdof = a_global_dof->dataKinetic();
      const Real& val = gdof[i_kinetic_species]->getVal(  gidx, *a_kinetic_species[i_kinetic_species],
                                                          new_address );
      a_global_dof->setAddress(a_idx, new_address);
      return(val);

    } else if (i_fluid_species != -1) {

      new_address.family("fluid");
      new_address.number(i_fluid_species);

      const GlobalDOFFluidSpeciesPtrVect& gdof = a_global_dof->dataFluid();
      //      const Real& val = gdof[i_fluid_species]->getVal(  gidx, *a_fluid_species[i_fluid_species],
      //                                                        new_address );
      const Real& val = gdof[i_fluid_species]->getVal(  gidx, static_cast<const CFG::FluidSpecies&>(*a_fluid_species[i_fluid_species]),
                                                        new_address );
      a_global_dof->setAddress(a_idx, new_address);
      return(val);
    } else {

      /* something went wrong! */
      std::cout << "Error in GKState/GKRHSData::operator[]: Index " << a_idx << " not found!\n";
      exit(1);
    }

  }
}

const Real& GKState::operator[] (int a_idx) const
{
  CH_assert(a_idx >= 0);
  CH_assert(a_idx < m_vector_size);

  return( getVal(a_idx, m_global_dof, m_kinetic_species, m_fluid_species) );
}

const Real& GKRHSData::operator[] (int a_idx) const
{
  CH_assert(a_idx >= 0);
  CH_assert(a_idx < m_vector_size);

  return( getVal(a_idx, m_global_dof, m_kinetic_species, m_fluid_species) );
}

static inline 
Real& getVal(
              int                       a_idx,
              GlobalDOF                 *a_global_dof,
              KineticSpeciesPtrVect&    a_kinetic_species,
              CFG::FluidSpeciesPtrVect& a_fluid_species
           )
{
  if (a_global_dof->savedAddress(a_idx)) {

    const GKAddress& address = a_global_dof->getAddress(a_idx);
    if (address.family() == "kinetic") {
      
      int i_kinetic_species = address.number();
      const GlobalDOFKineticSpeciesPtrVect& gdof = a_global_dof->dataKinetic();
      return( gdof[i_kinetic_species]->getVal(*a_kinetic_species[i_kinetic_species], address) );

    } else if (address.family() == "fluid" ) {

      int i_fluid_species = address.number();
      const GlobalDOFFluidSpeciesPtrVect& gdof = a_global_dof->dataFluid();
      //      return( gdof[i_fluid_species]->getVal(*a_fluid_species[i_fluid_species], address) );
      return( gdof[i_fluid_species]->getVal(static_cast<CFG::FluidSpecies&>(*a_fluid_species[i_fluid_species]), address) );
    } else {
      std::cout << "Error in GKState/GKRHSData::operator[]: Invalid value for "
                << "\"family\" in saved address.\n";
      exit(1);
    }

  } else {

    /* a_idx is local, so add MPI offset */
    int gidx = a_idx + a_global_dof->mpiOffset();

    int i_kinetic_species = a_global_dof->isIndexKineticSpecies (gidx,a_kinetic_species);
    int i_fluid_species   = a_global_dof->isIndexFluidSpecies   (gidx,a_fluid_species  );
    CH_assert(    ( (i_kinetic_species != -1) && (i_fluid_species == -1) )
              ||  ( (i_kinetic_species == -1) && (i_fluid_species != -1) )
              ||  ( (i_kinetic_species == -1) && (i_fluid_species == -1) ) );

    GKAddress new_address;

    if (i_kinetic_species != -1) {

      new_address.family("kinetic");
      new_address.number(i_kinetic_species);

      const GlobalDOFKineticSpeciesPtrVect& gdof = a_global_dof->dataKinetic();
      Real& val = gdof[i_kinetic_species]->getVal(  gidx, *a_kinetic_species[i_kinetic_species],
                                                    new_address );
      a_global_dof->setAddress(a_idx, new_address);
      return(val);

    } else if (i_fluid_species != -1) {

      new_address.family("fluid");
      new_address.number(i_fluid_species);

      const GlobalDOFFluidSpeciesPtrVect& gdof = a_global_dof->dataFluid();
      //      Real& val = gdof[i_fluid_species]->getVal(  gidx, *a_fluid_species[i_fluid_species],
      //                                                  new_address );
      Real& val = gdof[i_fluid_species]->getVal(  gidx, static_cast<CFG::FluidSpecies&>(*a_fluid_species[i_fluid_species]),
                                                  new_address );
      a_global_dof->setAddress(a_idx, new_address);
      return(val);
    } else {

      /* something went wrong! */
      std::cout << "Error in GKState/GKRHSData::operator[]: Index " << a_idx << " not found!\n";
      exit(1);
    }

  }
}

Real& GKState::operator[] (int a_idx)
{
  CH_assert(a_idx >= 0);
  CH_assert(a_idx < m_vector_size);

  return( getVal(a_idx, m_global_dof, m_kinetic_species, m_fluid_species) );
}

Real& GKRHSData::operator[] (int a_idx)
{
  CH_assert(a_idx >= 0);
  CH_assert(a_idx < m_vector_size);

  return( getVal(a_idx, m_global_dof, m_kinetic_species, m_fluid_species) );
}

#include "NamespaceFooter.H"
