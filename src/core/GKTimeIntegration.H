#ifndef _GK_TimeIntegration_H_
#define _GK_TimeIntegration_H_

#include "REAL.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "PhaseGeom.H"
#include "GlobalDOF.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluidSpecies.H"
#include "Field.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

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
         define( a_state.m_kinetic_species,
                 a_state.m_fluid_species,
                 a_state.m_fields,
                 a_state.m_geometry );
      }
   }
   
   /// Destructor
   virtual ~GKState() {;}
  
   inline
   void define( const KineticSpeciesPtrVect&    a_kinetic_species,
                const CFG::FluidSpeciesPtrVect& a_fluid_species,
                const CFG::FieldPtrVect&        a_fields,
                const RefCountedPtr<PhaseGeom>& a_geometry )
   {
      m_kinetic_species.resize( a_kinetic_species.size() );
      for (int s(0); s<a_kinetic_species.size(); s++) {
         m_kinetic_species[s] = a_kinetic_species[s]->clone( m_ghost_vect );
      }

      const CFG::IntVect ghost_vect( configSpaceGhostVector() );
      m_fluid_species.resize( a_fluid_species.size() );
      for (int s(0); s<a_fluid_species.size(); s++) {
         m_fluid_species[s] = a_fluid_species[s]->clone( ghost_vect );
      }
      
      m_fields.resize( a_fields.size() );
      for (int s(0); s<a_fields.size(); s++) {
         m_fields[s] = a_fields[s]->clone( ghost_vect );
      }

      m_geometry = a_geometry;
      CH_assert(!m_geometry.isNull());
      m_is_defined = true;

      int nv = getVectorSize();
      m_global_dof.define( nv,m_kinetic_species,m_fluid_species,m_fields );
   }
   
   
   inline
   void define( const GKState& a_state )
   {
      m_ghost_vect = a_state.m_ghost_vect;
      define( a_state.m_kinetic_species,
              a_state.m_fluid_species,
              a_state.m_fields,
              a_state.m_geometry );
   }
   
   
   inline
   void define( const GKState& a_state, const IntVect& a_ghost_vect )
   {
      m_ghost_vect = a_ghost_vect;
      define( a_state.m_kinetic_species,
              a_state.m_fluid_species,
              a_state.m_fields,
              a_state.m_geometry );
   }
   
   
   void copy( const GKState& a_state );

   void copy( const GKRHSData& a_rhs );
   
   
   inline
   void zero()
   {
      for (int s(0); s<m_kinetic_species.size(); s++) {
         m_kinetic_species[s]->zeroData();
      }
      for (int s(0); s<m_fluid_species.size(); s++) {
         m_fluid_species[s]->zeroData();
      }
      for (int s(0); s<m_kinetic_species.size(); s++) {
         m_fields[s]->zeroData();
      }
   }

   int getVectorSize();
   
   
   void increment( const GKState& a_state,
                   const Real& a_factor = 1.0,
                   const bool a_update_flux_register = false );

   void increment( const GKRHSData& a_rhs,
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
   KineticSpeciesPtrVect& dataKinetic()
   {
      CH_assert( isDefined() );
      return m_kinetic_species;
   }

   inline 
   const GlobalDOF& globalDOF() const
   {
     CH_assert(isDefined());
     return m_global_dof;
   }
   
   inline
   const KineticSpeciesPtrVect& dataKinetic() const
   {
      CH_assert( isDefined() );
      return m_kinetic_species;
   }

   inline
   CFG::FluidSpeciesPtrVect& dataFluid()
   {
      CH_assert( isDefined() );
      return m_fluid_species;
   }
   
   inline
   const CFG::FluidSpeciesPtrVect& dataFluid() const
   {
      CH_assert( isDefined() );
      return m_fluid_species;
   }

   inline
   CFG::FieldPtrVect& dataField()
   {
      CH_assert( isDefined() );
      return m_fields;
   }
   
   inline
   const CFG::FieldPtrVect& dataField() const
   {
      CH_assert( isDefined() );
      return m_fields;
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

   void copyTo( Real* a_Y ) const;

   void copyFrom( Real* a_Y );

   void addFrom( Real *a_Y, const Real& a_a=1.0 );

private:

   CFG::IntVect configSpaceGhostVector()
   {
      CFG::IntVect ghost_vect;
      for (int d(0); d<CFG_DIM; d++) {
         ghost_vect[d] = m_ghost_vect[d];
      }
      return ghost_vect;
   }
   
   // Disallowed for all the usual reasons
   void operator=(const GKState&);
//   GKState(const GKState&);

   bool                     m_is_defined;
   IntVect                  m_ghost_vect;
   KineticSpeciesPtrVect    m_kinetic_species;
   CFG::FluidSpeciesPtrVect m_fluid_species;
   CFG::FieldPtrVect        m_fields;
   RefCountedPtr<PhaseGeom> m_geometry;
   GlobalDOF                m_global_dof;
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
   void define( const KineticSpeciesPtrVect&    a_kinetic_species,
                const CFG::FluidSpeciesPtrVect& a_fluid_species,
                const CFG::FieldPtrVect&        a_fields )
   {
      m_kinetic_species.resize( a_kinetic_species.size() );
      for (int s(0); s<a_kinetic_species.size(); s++) {
         m_kinetic_species[s] = a_kinetic_species[s]->clone( IntVect::Zero );
      }
      m_fluid_species.resize( a_fluid_species.size() );
      for (int s(0); s<a_fluid_species.size(); s++) {
         m_fluid_species[s] = a_fluid_species[s]->clone( CFG::IntVect::Zero );
      }
      m_fields.resize( a_fields.size() );
      for (int s(0); s<a_fields.size(); s++) {
         m_fields[s] = a_fields[s]->clone( CFG::IntVect::Zero );
      }
      m_is_defined = true;

      int nv = getVectorSize();
      m_global_dof.define( nv,m_kinetic_species,m_fluid_species,m_fields );
   }

   
   inline
   void define( const GKState& a_state )
   {
      define( a_state.dataKinetic(), a_state.dataFluid(), a_state.dataField() );
   }

   
   inline
   void define( const GKRHSData& a_rhs )
   {
      define( a_rhs.dataKinetic(), a_rhs.dataFluid(), a_rhs.dataField() );
   }

   void scale(const Real& a_factor);

   
   void copy( const GKRHSData& a_rhs );

   void copy( const GKState& a_rhs );

   
   inline
   void zero()
   {
      for (int s(0); s<m_kinetic_species.size(); s++) {
         m_kinetic_species[s]->zeroData();
      }
      for (int s(0); s<m_fluid_species.size(); s++) {
         m_fluid_species[s]->zeroData();
      }
      for (int s(0); s<m_fields.size(); s++) {
         m_fields[s]->zeroData();
      }
   }

   
   void increment( const GKRHSData& a_increment,
                   const Real& a_factor = 1.0,
                   const bool a_update_flux_register = false );


   void increment( const GKState& a_increment,
                   const Real& a_factor = 1.0,
                   const bool a_update_flux_register = false );


   Real computeNorm ( const int a_p ) const;

   Real dotProduct (const GKRHSData& a_Y);
   
   inline
   bool isDefined() const
   {
      return m_is_defined;
   }

   
   inline
   KineticSpeciesPtrVect& dataKinetic()
   {
      CH_assert( isDefined() );
      return m_kinetic_species;
   }
   
   inline 
   const GlobalDOF& globalDOF() const
   {
     CH_assert(isDefined());
     return m_global_dof;
   }
   
   inline
   const KineticSpeciesPtrVect& dataKinetic() const
   {
      CH_assert( isDefined() );
      return m_kinetic_species;
   }

   inline
   CFG::FluidSpeciesPtrVect& dataFluid()
   {
      CH_assert( isDefined() );
      return m_fluid_species;
   }
   
   inline
   const CFG::FluidSpeciesPtrVect& dataFluid() const
   {
      CH_assert( isDefined() );
      return m_fluid_species;
   }

   inline
   CFG::FieldPtrVect& dataField()
   {
      CH_assert( isDefined() );
      return m_fields;
   }
   
   inline
   const CFG::FieldPtrVect& dataField() const
   {
      CH_assert( isDefined() );
      return m_fields;
   }

   void copyTo( Real* a_Y ) const;

   void copyFrom( Real* a_Y );

   int getVectorSize();
   
private:

   // Disallowed for all the usual reasons
   void operator=(const GKRHSData&);
   GKRHSData(const GKRHSData&);

   bool                     m_is_defined;
   KineticSpeciesPtrVect    m_kinetic_species;
   CFG::FluidSpeciesPtrVect m_fluid_species;
   CFG::FieldPtrVect        m_fields;
   GlobalDOF                m_global_dof;
};


#include "NamespaceFooter.H"

#endif

