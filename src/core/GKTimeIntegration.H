#ifndef _GK_TimeIntegration_H_
#define _GK_TimeIntegration_H_

#include "REAL.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "PhaseGeom.H"

#include "NamespaceHeader.H"

class GKRHSData;

class GKState
{
public:

   /// Default constructor
   /**
      Object requires define() to be called before all other functions.
   */
   GKState()
      : m_is_defined( false ),
        m_ghost_vect( IntVect::Zero )
   {
   }
   
   /// Constructor
   /**
      Object requires define() to be called before all other functions.
   */
   inline
   GKState( const IntVect& a_ghost_vect )
      : m_is_defined( false ),
        m_ghost_vect( a_ghost_vect )
   {
   }
   
   inline
   GKState( const GKState& a_state )
      : m_is_defined( a_state.m_is_defined ),
        m_ghost_vect( a_state.m_ghost_vect )
   {
      if (&a_state!=this) {
         define( a_state.m_species_mapped, a_state.m_geometry );
      }
   }
   
   /// Destructor
   virtual ~GKState() {;}
  
   inline
   void define( const KineticSpeciesPtrVect& a_species,
                const RefCountedPtr<PhaseGeom>& a_geometry )
   {
      const int n( a_species.size() );
      m_species_mapped.resize( n );
      for (int s(0); s<n; s++) {
         m_species_mapped[s] = a_species[s]->clone( m_ghost_vect );
      }
      m_geometry = a_geometry;
      CH_assert(!m_geometry.isNull());
      m_is_defined = true;
   }
   
   
   inline
   void define( const GKState& a_state )
   {
      m_ghost_vect = a_state.m_ghost_vect;
      define( a_state.m_species_mapped, a_state.m_geometry );
   }
   
   
   inline
   void copy( const GKState& a_state )
   {
      m_species_mapped.resize( a_state.m_species_mapped.size() );
      for (int s(0); s<m_species_mapped.size(); s++) {
         m_species_mapped[s]->copy( *(a_state.m_species_mapped[s]) );
      }
   }
   
   
   inline
   void zero()
   {
      for (int s(0); s<m_species_mapped.size(); s++) {
         m_species_mapped[s]->zeroData();
      }
   }

   
   void increment( const GKState& a_state,
                   Real a_factor = 1.0,
                   bool a_update_flux_register = false );

   void increment( const GKRHSData& a_rhs,
                   Real a_factor = 1.0,
                   bool a_update_flux_register = false );


   void computeNorms( Real& a_norm_1, Real& a_norm_2, Real& a_norm_inf);

   
   inline
   bool isDefined() const
   {
      return m_is_defined;
   }

   
   inline
   KineticSpeciesPtrVect& data()
   {
      CH_assert( isDefined() );
      return m_species_mapped;
   }
   
   inline
   const KineticSpeciesPtrVect& data() const
   {
      CH_assert( isDefined() );
      return m_species_mapped;
   }

   inline
   RefCountedPtr<PhaseGeom>& geometry()
   {
      CH_assert( isDefined() );
      return m_geometry;
   }
   
   inline
   const RefCountedPtr<PhaseGeom>& geometry() const
   {
      CH_assert( isDefined() );
      return m_geometry;
   }

private:

   // Disallowed for all the usual reasons
   void operator=(const GKState&);
//   GKState(const GKState&);

   bool                     m_is_defined;
   IntVect                  m_ghost_vect;
   KineticSpeciesPtrVect    m_species_mapped;
   RefCountedPtr<PhaseGeom> m_geometry;
};


class GKRHSData
{
public:

   
   inline
   GKRHSData()
      : m_is_defined( false )
   {
   }

   
   inline
   void define( const KineticSpeciesPtrVect& a_species )
   {
      const int n( a_species.size() );
      m_species_mapped.resize( n );
      for (int s(0); s<n; s++) {
         m_species_mapped[s] = a_species[s]->clone( IntVect::Zero );
      }
      m_is_defined = true;
   }

   
   inline
   void define( const GKState& a_state )
   {
      define( a_state.data() );
   }

   
   inline
   void copy( const GKRHSData& a_rhs )
   {
      m_species_mapped.resize( a_rhs.m_species_mapped.size() );
      for (int s(0); s<m_species_mapped.size(); s++) {
         m_species_mapped[s]->copy( *(a_rhs.m_species_mapped[s]) );
      }
   }

   
   inline
   void copy( const GKState& a_rhs )
   {
      const KineticSpeciesPtrVect& rhs_species_mapped( a_rhs.data() ); 
      m_species_mapped.resize( rhs_species_mapped.size() );
      for (int s(0); s<m_species_mapped.size(); s++) {
         m_species_mapped[s]->copy( *(rhs_species_mapped[s]) );
      }
   }

   
   inline
   void zero()
   {
      for (int s(0); s<m_species_mapped.size(); s++) {
         m_species_mapped[s]->zeroData();
      }
   }

   
   void increment( const GKRHSData& a_increment,
                   Real a_factor = 1.0,
                   bool a_update_flux_register = false );


   void computeNorms( Real& a_norm_1, Real& a_norm_2, Real& a_norm_inf);

   
   inline
   bool isDefined() const
   {
      return m_is_defined;
   }

   
   inline
   KineticSpeciesPtrVect& data()
   {
      CH_assert( isDefined() );
      return m_species_mapped;
   }
   
   inline
   const KineticSpeciesPtrVect& data() const
   {
      CH_assert( isDefined() );
      return m_species_mapped;
   }
   
private:

   // Disallowed for all the usual reasons
   void operator=(const GKRHSData&);
   GKRHSData(const GKRHSData&);

   bool                  m_is_defined;
   KineticSpeciesPtrVect m_species_mapped;
};


#include "NamespaceFooter.H"

#endif

