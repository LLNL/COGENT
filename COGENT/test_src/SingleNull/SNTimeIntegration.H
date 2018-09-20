#ifndef _SN_TimeIntegration_H_
#define _SN_TimeIntegration_H_

#include "REAL.H"
#include "RefCountedPtr.H"
#include "FluidSpecies.H"
#include "MagGeom.H"
#include "SNGlobalDOF.H"

#include "NamespaceHeader.H"

class SNRHSData;

class SNState
{
public:

   /// Default constructor
   /**
      Object requires define() to be called before all other functions.
   */
   SNState()
      : m_is_defined( false ),
        m_ghost_vect( IntVect::Zero )
   {
   }
   
   /// Constructor
   /**
      Object requires define() to be called before all other functions.
   */
   inline
   SNState( const IntVect& a_ghost_vect )
      : m_is_defined( false ),
        m_ghost_vect( a_ghost_vect )
   {
   }
   
   inline
   SNState( const SNState& a_state )
      : m_is_defined( a_state.m_is_defined ),
        m_ghost_vect( a_state.m_ghost_vect )
   {
      if (&a_state!=this) {
         define( a_state.m_species, a_state.m_geometry );
      }
   }
   
   /// Destructor
   virtual ~SNState() {;}
  
   inline
   void define( const FluidSpeciesPtrVect& a_species,
                const RefCountedPtr<MagGeom>& a_geometry )
   {
      m_species.resize( a_species.size() );
      for (int s(0); s<a_species.size(); s++) {
         m_species[s] = a_species[s]->clone( m_ghost_vect );
      }

      m_geometry = a_geometry;
      CH_assert(!m_geometry.isNull());
      m_is_defined = true;
   }
   
   
   inline
   void define( const SNState& a_state )
   {
      m_ghost_vect = a_state.m_ghost_vect;
      define( a_state.m_species, a_state.m_geometry );
   }
   
   
   inline
   void define( const SNState& a_state, const IntVect& a_ghost_vect )
   {
      m_ghost_vect = a_ghost_vect;
      define( a_state.m_species, a_state.m_geometry );
   }
   
   
   void copy( const SNState& a_state );

   void copy( const SNRHSData& a_rhs );
   
   inline
   void setGlobalDOF( SNGlobalDOF *a_global_dof)
   {
     m_global_dof = a_global_dof;
   }

   inline
   SNGlobalDOF* getGlobalDOF() const
   {
     return m_global_dof;
   }
   

   
   inline
   void zero()
   {
      for (int s(0); s<m_species.size(); s++) {
         m_species[s]->zeroData();
      }
   }

   int getVectorSize() const;
   
   void increment( const SNState& a_state,
                   const Real& a_factor = 1.0,
                   const bool a_update_flux_register = false );

   void increment( const SNRHSData& a_rhs,
                   const Real& a_factor = 1.0,
                   const bool a_update_flux_register = false );

   Real computeNorm ( const int a_p ) const; 

   void scale( const Real& a_factor );

   
   inline
   bool isDefined() const
   {
      return m_is_defined;
   }

   
   inline
   FluidSpeciesPtrVect& data()
   {
      CH_assert( isDefined() );
      return m_species;
   }

   inline
   const FluidSpeciesPtrVect& data() const
   {
      CH_assert( isDefined() );
      return m_species;
   }

   inline
   RefCountedPtr<MagGeom>& geometry()
   {
      CH_assert( isDefined() );
      return m_geometry;
   }
   
   inline
   const RefCountedPtr<MagGeom>& geometry() const
   {
      CH_assert( isDefined() );
      return m_geometry;
   }

   void copyTo( Real* a_Y ) const;

   void copyFrom( const Real* a_Y );

   void addFrom( Real *a_Y, const Real& a_a=1.0 );

private:
   
   // Disallowed for all the usual reasons
   void operator=(const SNState&);
   
   bool                   m_is_defined;
   IntVect                m_ghost_vect;
   FluidSpeciesPtrVect    m_species;
   RefCountedPtr<MagGeom> m_geometry;
   SNGlobalDOF                *m_global_dof;
};


class SNRHSData
{
public:

   
   inline
   SNRHSData()
      : m_is_defined( false )
   {
   }

   
   inline
   void define( const FluidSpeciesPtrVect& a_species )
   {
      m_species.resize( a_species.size() );
      for (int s(0); s<a_species.size(); s++) {
         m_species[s] = a_species[s]->clone( IntVect::Zero );
      }
      m_is_defined = true;
   }

   
   inline
   void define( const SNState& a_state )
   {
      define( a_state.data() );
   }

   
   inline
   void define( const SNRHSData& a_rhs )
   {
      define( a_rhs.data() );
   }

   void scale(const Real& a_factor);
   
   void copy( const SNRHSData& a_rhs );

   void copy( const SNState& a_rhs );
   
   inline
   void zero()
   {
      for (int s(0); s<m_species.size(); s++) {
         m_species[s]->zeroData();
      }
   }
   
   void increment( const SNRHSData& a_increment,
                   const Real& a_factor = 1.0,
                   const bool a_update_flux_register = false );


   void increment( const SNState& a_increment,
                   const Real& a_factor = 1.0,
                   const bool a_update_flux_register = false );


   Real computeNorm ( const int a_p ) const;

   Real dotProduct (const SNRHSData& a_Y);
   
   inline
   bool isDefined() const
   {
      return m_is_defined;
   }

   
   inline
   FluidSpeciesPtrVect& data()
   {
      CH_assert( isDefined() );
      return m_species;
   }
   
   inline
   const FluidSpeciesPtrVect& data() const
   {
      CH_assert( isDefined() );
      return m_species;
   }

   void copyTo( Real* a_Y ) const;

   void copyFrom( const Real* a_Y );

   int getVectorSize();
   
private:

   // Disallowed for all the usual reasons
   void operator=(const SNRHSData&);
   SNRHSData(const SNRHSData&);

   bool                m_is_defined;
   FluidSpeciesPtrVect m_species;
};


#include "NamespaceFooter.H"

#endif

