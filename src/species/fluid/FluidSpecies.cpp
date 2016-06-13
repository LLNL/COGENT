#include "FluidSpecies.H"

#include "NamespaceHeader.H"

FluidSpecies::FluidSpecies( const string&  a_name,
                            const Real     a_mass,
                            const Real     a_charge,
                            const string&  a_model,
                            const MagGeom& a_geometry )
  : m_name( a_name ),
    m_model( a_model ),
    m_mass( a_mass ),
    m_charge( a_charge ),
    m_geometry( a_geometry )
{
}


FluidSpecies::FluidSpecies( const FluidSpecies& a_foo )
  : m_name( a_foo.m_name ),
    m_model( a_foo.m_model ),
    m_mass( a_foo.m_mass ),
    m_charge( a_foo.m_charge ),
    m_geometry( a_foo.m_geometry )
{
   m_state.define( a_foo.m_state );
}

void FluidSpecies::numberDensity( LevelData<FArrayBox>& a_rho ) const
{
   massDensity( a_rho );
   for (DataIterator dit(a_rho.dataIterator()); dit.ok(); ++dit ) {
      a_rho[dit].divide( m_mass );
   }
}


void FluidSpecies::massDensity( LevelData<FArrayBox>& a_rho ) const
{
   DisjointBoxLayout dbl( m_state.disjointBoxLayout() );
   CH_assert( dbl.compatible( a_rho.disjointBoxLayout() ) );
   Interval src_interval( MASS_DENSITY, MASS_DENSITY );
   Interval dst_interval( 0 , 0 );
   m_state.copyTo( src_interval, a_rho, dst_interval );
}


void FluidSpecies::chargeDensity( LevelData<FArrayBox>& a_rho ) const
{
   massDensity( a_rho );
   Real ratio( m_charge / m_mass );
   for (DataIterator dit( a_rho.dataIterator()); dit.ok(); ++dit ) {
      a_rho[dit].mult( ratio );
   }
}


void FluidSpecies::pressure( LevelData<FArrayBox>& a_pressure ) const
{
}

void FluidSpecies::temperature( LevelData<FArrayBox>& a_temperature ) const
{
}

void FluidSpecies::velocity( LevelData<FArrayBox>& a_velocity ) const
{
}


bool FluidSpecies::isSpecies( const string& a_name ) const
{
   if (m_name == a_name) return true;
   return false;
}


const FluidSpecies& FluidSpecies::operator=( const FluidSpecies& a_rhs )
{
   if (&a_rhs != this)
   {
      CH_assert( (m_state.disjointBoxLayout()).compatible( a_rhs.m_state.disjointBoxLayout() ) );
      m_name = a_rhs.m_name;
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;
      m_state.define( a_rhs.m_state );
   }
   return *this;
}


void FluidSpecies::copy( const FluidSpecies& a_rhs )
{
   if (&a_rhs != this)
   {
      CH_assert( (m_state.disjointBoxLayout()).compatible( a_rhs.m_state.disjointBoxLayout() ) );
      m_name = a_rhs.m_name;
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;

      for (DataIterator dit( m_state.dataIterator() ); dit.ok(); ++dit) {
         m_state[dit].copy( a_rhs.m_state[dit] );
      }
   }
}


void FluidSpecies::zeroData()
{
   DataIterator dit( m_state.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit) {
      m_state[dit].setVal( 0.0 );
   }
}


void FluidSpecies::addData( const FluidSpecies& a_rhs,
                            const Real a_factor )
{
   CH_assert( (m_state.disjointBoxLayout()).compatible( a_rhs.m_state.disjointBoxLayout() ) );
   for (DataIterator dit( m_state.dataIterator() ); dit.ok(); ++dit) {
      m_state[dit].plus( a_rhs.m_state[dit], a_factor );
   }
}


bool FluidSpecies::conformsTo( const FluidSpecies& a_rhs,
                               const bool a_include_ghost_cells ) const
{
   const LevelData<FArrayBox>& thisData( m_state );
   const LevelData<FArrayBox>& rhsData( a_rhs.m_state );
   
   const DisjointBoxLayout& thisBoxes( thisData.disjointBoxLayout() );
   const DisjointBoxLayout& rhsBoxes( rhsData.disjointBoxLayout() );
   
   bool status( true );
   status &= thisBoxes.compatible( rhsBoxes );
   status &= ( thisData.nComp() == rhsData.nComp() );
   
   if ( a_include_ghost_cells ) {
      status &= ( thisData.ghostVect() == rhsData.ghostVect() );
   }
   
   return status;
}


RefCountedPtr<FluidSpecies>
FluidSpecies::clone( const IntVect a_ghost_vect, const bool a_copy_soln_data ) const
{
   RefCountedPtr<FluidSpecies> result
      = RefCountedPtr<FluidSpecies>(
         new FluidSpecies( m_name, m_mass, m_charge, m_model, m_geometry ) );

   result->m_state.define( m_state.disjointBoxLayout(),
                           m_state.nComp(),
                           a_ghost_vect );

   if (a_copy_soln_data) {
      LevelData<FArrayBox>& result_state( result->m_state );
      for (DataIterator dit( result_state.dataIterator() ); dit.ok(); ++dit) {
         result_state[dit].copy( m_state[dit] );
      }
   }

   return result;
}


#include "NamespaceFooter.H"
