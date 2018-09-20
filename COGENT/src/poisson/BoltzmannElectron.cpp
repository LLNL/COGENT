#include "BoltzmannElectron.H"

#include "NamespaceHeader.H"

const string BoltzmannElectron::s_NAME( "Boltzmann Electron" );



BoltzmannElectron::BoltzmannElectron(
         const Real                      a_mass,
         const Real                      a_charge,
         const MultiBlockLevelGeom&      a_geom,
         const LevelData<FArrayBox>&     a_temperature )
   : m_mass( a_mass ),
     m_charge( a_charge ),
     m_geom( a_geom )
{
   m_density.define( a_geom.grids(), 1, IntVect::Zero );

   CH_assert( a_temperature.ghostVect() >= IntVect::Unit );
   m_temperature.define( a_temperature.disjointBoxLayout(), 1, a_temperature.ghostVect() );
   DataIterator dit = m_temperature.dataIterator();
   for ( dit.begin(); dit.ok(); ++dit ) {
      m_temperature[dit].copy(a_temperature[dit]);
   }
   m_temperature.exchange();
}



BoltzmannElectron::BoltzmannElectron( const BoltzmannElectron& a_foo )
   : m_mass( a_foo.m_mass ),
     m_charge( a_foo.m_charge ),
     m_geom( a_foo.m_geom )
{
   m_temperature.define( a_foo.m_temperature );
   m_density.define( a_foo.m_density );
}



void BoltzmannElectron::chargeDensity( LevelData<FArrayBox>& a_rho ) const
{
   a_rho.define( m_density );
}



bool BoltzmannElectron::isSpecies( const string& a_name ) const
{
   if (a_name == s_NAME) return true;
   return false;
}



const BoltzmannElectron& BoltzmannElectron::operator=(
   const BoltzmannElectron& a_rhs )
{
   if (&a_rhs != this)
   {
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;
      m_temperature.define( a_rhs.m_temperature );
      m_density.define( a_rhs.m_density );
   }
   return *this;
}



void BoltzmannElectron::copy( const BoltzmannElectron& a_rhs )
{
   if (&a_rhs != this)
   {
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;

      DataIterator dit1( m_temperature.dataIterator() );
      for (dit1.begin(); dit1.ok(); ++dit1) {
         m_temperature[dit1].copy( a_rhs.m_temperature[dit1] );
      }

      DataIterator dit2( m_density.dataIterator() );
      for (dit2.begin(); dit2.ok(); ++dit2) {
         m_density[dit2].copy( a_rhs.m_density[dit2] );
      }
   }
}



void BoltzmannElectron::zeroData()
{
   DataIterator dit( m_density.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit) {
      m_density[dit].setVal( 0.0 );
   }
}



void BoltzmannElectron::addData( const BoltzmannElectron& a_rhs,
                                 const Real          a_factor )
{
   try {
      DataIterator dit( m_density.dataIterator() );
      for (dit.begin(); dit.ok(); ++dit) {
         m_density[dit].plus( a_rhs.m_density[dit], a_factor );
      }
   }
   catch (std::bad_cast) {
      MayDay::Error( "Invalid BoltzmannElectron passed to BoltzmannElectron::addData!" );
   }
}



bool BoltzmannElectron::conformsTo( const BoltzmannElectron& a_rhs,
                                    const bool a_include_ghost_cells ) const
{
   try {
      const LevelData<FArrayBox>& thisData = m_density;
      const LevelData<FArrayBox>& rhsData = a_rhs.m_density;

      const DisjointBoxLayout& thisBoxes = thisData.disjointBoxLayout();
      const DisjointBoxLayout& rhsBoxes = rhsData.disjointBoxLayout();

      bool status( true );
      status &= thisBoxes.compatible( rhsBoxes );
      status &= ( thisData.nComp() == rhsData.nComp() );

      if ( a_include_ghost_cells) {
         status &= ( thisData.ghostVect() == rhsData.ghostVect() );
      }

      return status;
   }
   catch (std::bad_cast) {
      MayDay::Error( "Invalid BoltzmannElectron passed to BoltzmannElectron::conformsTo!" );
   }
   return false;
}


#if 0
RefCountedPtr<BoltzmannElectron>
BoltzmannElectron::clone( const IntVect ghostVect,
                          const bool copy_soln_data ) const
{
   RefCountedPtr<BoltzmannElectron> result
      = RefCountedPtr<BoltzmannElectron>(
                                         new BoltzmannElectron( m_mass, m_charge, m_geom, m_temperature) );

   IntVect cfg_ghostVect;
   for (int dir=0; dir<CFG_DIM; ++dir) {
      cfg_ghostVect[dir] = ghostVect[dir];
   }
   result->m_density.define( m_density.disjointBoxLayout(),
                             m_density.nComp(),
                             cfg_ghostVect );

   if (copy_soln_data) {

      LevelData<FArrayBox>& result_charge_density = result->m_density;
      DataIterator dit1( result_charge_density.dataIterator() );
      for (dit1.begin(); dit1.ok(); ++dit1) {
         result_charge_density[dit1].copy( m_density[dit1] );
      }

      LevelData<FArrayBox>& result_temperature = result->m_temperature;
      DataIterator dit2( result_temperature.dataIterator() );
      for (dit2.begin(); dit2.ok(); ++dit2) {
         result_temperature[dit2].copy( m_temperature[dit2] );
      }
   }

   return result;
}
#endif



#include "NamespaceFooter.H"
