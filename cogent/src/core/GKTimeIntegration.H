#ifndef _GK_TimeIntegration_H_
#define _GK_TimeIntegration_H_

#include "REAL.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "PhaseGeom.H"
#include "GlobalDOF.H"
#include "Scalar.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluidSpecies.H"
#include "CFGVar.H"
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
//        m_is_defined_ref( false ),
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
//        m_is_defined_ref( false ),
        m_ghost_vect( a_ghost_vect ),
        m_global_dof ( NULL )
   {
   }
   
   inline
   GKState( const GKState& a_state )
      : m_is_defined( a_state.m_is_defined ),
//        m_is_defined_ref( a_state.m_is_defined_ref ),
        m_ghost_vect( a_state.m_ghost_vect )
   {
      if (&a_state!=this) {
         define( a_state.m_kinetic_species,
                 a_state.m_fluid_species,
                 a_state.m_scalars,
                 a_state.m_geometry );
//         if (a_state.isDefined_ref()){
//             define_ref( a_state.m_species_mapped_ref, a_state.m_geometry );
//         }
         m_global_dof = a_state.m_global_dof;
      }
   }
   
   /// Destructor
   virtual ~GKState() {;}
  
   inline
   void define( const KineticSpeciesPtrVect&    a_kinetic_species,
                const CFG::FluidSpeciesPtrVect& a_fluid_species,
                const ScalarPtrVect&            a_scalars,
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
      
      m_scalars.resize( a_scalars.size() );
      for (int s(0); s<a_scalars.size(); s++) {
         m_scalars[s] = a_scalars[s]->clone();
      }

      m_geometry = a_geometry;
      CH_assert(!m_geometry.isNull());
      m_is_defined = true;
      m_vector_size = getVectorSize();
   }
   
//   inline
//   void define_ref( const KineticSpeciesPtrVect& a_species,
//                const RefCountedPtr<PhaseGeom>& a_geometry )
//   {
//      const int n( a_species.size() );
//      m_species_mapped_ref.resize( n );
//      for (int s(0); s<n; s++) {
//         m_species_mapped_ref[s] = a_species[s]->clone( m_ghost_vect );
//      }
//      m_geometry = a_geometry;
//      CH_assert(!m_geometry.isNull());
//      m_is_defined_ref = true;
//   }  

   
   inline
   void define( const GKState& a_state )
   {
      m_ghost_vect = a_state.m_ghost_vect;
      define( a_state.m_kinetic_species,
              a_state.m_fluid_species,
              a_state.m_scalars,
              a_state.m_geometry );
//      if(a_state.isDefined_ref() ){
//          define_ref( a_state.m_species_mapped_ref, a_state.m_geometry );
//      }
      m_global_dof = a_state.m_global_dof;
   }
   
//   inline
//   void define_ref( const GKState& a_state )
//   {
//      m_ghost_vect = a_state.m_ghost_vect;
//      define_ref( a_state.m_species_mapped_ref, a_state.m_geometry );
//   }
  
   
   inline
   void define( const GKState& a_state, const IntVect& a_ghost_vect )
   {
      m_ghost_vect = a_ghost_vect;
      define( a_state.m_kinetic_species,
              a_state.m_fluid_species,
              a_state.m_scalars,
              a_state.m_geometry );
      m_global_dof = a_state.m_global_dof;
   }
   
   
   inline
   void setGlobalDOF( GlobalDOF *a_global_dof)
   {
     m_global_dof = a_global_dof;
   }

   inline
   GlobalDOF* getGlobalDOF() const
   {
     return m_global_dof;
   }
   
   
   inline
   void zero()
   {
      for (int s(0); s<m_kinetic_species.size(); s++) {
         m_kinetic_species[s]->zeroData();
      }
      for (int s(0); s<m_fluid_species.size(); s++) {
         m_fluid_species[s]->zeroData();
      }
      for (int s(0); s<m_scalars.size(); s++) {
         m_scalars[s]->zeroData();
      }
   }

   int getVectorSize() const;
   
   
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

//   inline
//   bool isDefined_ref() const
//   {
//      return m_is_defined_ref;
//   }
   
   inline
   KineticSpeciesPtrVect& dataKinetic()
   {
      CH_assert( isDefined() );
      return m_kinetic_species;
   }

   inline
   const KineticSpeciesPtrVect& dataKinetic() const
   {
      CH_assert( isDefined() );
      return m_kinetic_species;
   }

   inline
   int getKineticComponent( const string& name ) const
   {
      int component = -1;
      for (int i=0; i<m_kinetic_species.size(); ++i) {
         if (m_kinetic_species[i]->name() == name) {
            component = i;
            break;
         }
      }
      CH_assert(component >= 0);
      return component;
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
   int getFluidComponent( const string& name ) const
   {
      int component = -1;
      for (int i=0; i<m_fluid_species.size(); ++i) {
         if (m_fluid_species[i]->name() == name) {
            component = i;
            break;
         }
      }
      CH_assert(component >= 0);
      return component;
   }

   inline
   ScalarPtrVect& dataScalar()
   {
      CH_assert( isDefined() );
      return m_scalars;
   }
   
   inline
   const ScalarPtrVect& dataScalar() const
   {
      CH_assert( isDefined() );
      return m_scalars;
   }

   inline
   int getScalarComponent( const string& name ) const
   {
      int component = -1;
      for (int i=0; i<m_scalars.size(); ++i) {
         if (m_scalars[i]->name() == name) {
            component = i;
            break;
         }
      }
      CH_assert(component >= 0);
      return component;
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

   void copyFrom( const Real* a_Y );

   void addFrom( Real *a_Y, const Real& a_a=1.0 );

   const Real& operator[] (int idx) const;
   
   Real& operator[] (int idx);
   
   CFG::IntVect configSpaceGhostVector()
   {
      CFG::IntVect ghost_vect;
      for (int d(0); d<CFG_DIM; d++) {
         ghost_vect[d] = m_ghost_vect[d];
      }
      return ghost_vect;
   }
   
private:

   // Disallowed for all the usual reasons
   void operator=(const GKState&);
//   GKState(const GKState&);

   bool                     m_is_defined;
//   bool                     m_is_defined_ref;
   int                      m_vector_size;
   IntVect                  m_ghost_vect;
   KineticSpeciesPtrVect    m_kinetic_species;
   CFG::FluidSpeciesPtrVect m_fluid_species;
   ScalarPtrVect            m_scalars;
   RefCountedPtr<PhaseGeom> m_geometry;
   GlobalDOF                *m_global_dof;
};


class GKRHSData
{
public:

   
   inline
   GKRHSData()
      : m_is_defined( false ),
        m_global_dof ( NULL )
   {
   }

   
   inline
   void define( const KineticSpeciesPtrVect&    a_kinetic_species,
                const CFG::FluidSpeciesPtrVect& a_fluid_species,
                const ScalarPtrVect&            a_scalars )
   {
      m_kinetic_species.resize( a_kinetic_species.size() );
      for (int s(0); s<a_kinetic_species.size(); s++) {
         m_kinetic_species[s] = a_kinetic_species[s]->clone( IntVect::Zero );
      }
      m_fluid_species.resize( a_fluid_species.size() );
      for (int s(0); s<a_fluid_species.size(); s++) {
         m_fluid_species[s] = a_fluid_species[s]->clone( CFG::IntVect::Zero );
      }
      m_scalars.resize( a_scalars.size() );
      for (int s(0); s<a_scalars.size(); s++) {
         m_scalars[s] = a_scalars[s]->clone();
      }
      m_is_defined = true;
      m_vector_size = getVectorSize();
   }

   
   inline
   void define( const GKState& a_state )
   {
      define( a_state.dataKinetic(), a_state.dataFluid(), a_state.dataScalar() );
      m_global_dof = a_state.getGlobalDOF();
   }

   
   inline
   void define( const GKRHSData& a_rhs )
   {
      define( a_rhs.dataKinetic(), a_rhs.dataFluid(), a_rhs.dataScalar() );
      m_global_dof = a_rhs.getGlobalDOF();
   }

   void scale(const Real& a_factor);

   
   inline
   void setGlobalDOF( GlobalDOF *a_global_dof)
   {
     m_global_dof = a_global_dof;
   }

   inline
   GlobalDOF* getGlobalDOF() const
   {
     return m_global_dof;
   }
   
   inline
   void zero()
   {
      for (int s(0); s<m_kinetic_species.size(); s++) {
         m_kinetic_species[s]->zeroData();
      }
      for (int s(0); s<m_fluid_species.size(); s++) {
         m_fluid_species[s]->zeroData();
      }
      for (int s(0); s<m_scalars.size(); s++) {
         m_scalars[s]->zeroData();
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
   const KineticSpeciesPtrVect& dataKinetic() const
   {
      CH_assert( isDefined() );
      return m_kinetic_species;
   }

   inline
   int getKineticComponent( const string& name ) const
   {
      int component = -1;
      for (int i=0; i<m_kinetic_species.size(); ++i) {
         if (m_kinetic_species[i]->name() == name) {
            component = i;
            break;
         }
      }
      CH_assert(component >= 0);
      return component;
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
   int getFluidComponent( const string& name ) const
   {
      int component = -1;
      for (int i=0; i<m_fluid_species.size(); ++i) {
         if (m_fluid_species[i]->name() == name) {
            component = i;
            break;
         }
      }
      CH_assert(component >= 0);
      return component;
   }

   inline
   ScalarPtrVect& dataScalar()
   {
      CH_assert( isDefined() );
      return m_scalars;
   }
   
   inline
   const ScalarPtrVect& dataScalar() const
   {
      CH_assert( isDefined() );
      return m_scalars;
   }

   inline
   int getScalarComponent( const string& name ) const
   {
      int component = -1;
      for (int i=0; i<m_scalars.size(); ++i) {
         if (m_scalars[i]->name() == name) {
            component = i;
            break;
         }
      }
      CH_assert(component >= 0);
      return component;
   }

   void copyTo( Real* a_Y ) const;

   void copyFrom( const Real* a_Y );

   int getVectorSize() const;

   const Real& operator[] (int idx) const;
   
   Real& operator[] (int idx);
   
private:

   // Disallowed for all the usual reasons
   void operator=(const GKRHSData&);
   GKRHSData(const GKRHSData&);

   bool                     m_is_defined;
   int                      m_vector_size;
   KineticSpeciesPtrVect    m_kinetic_species;
   CFG::FluidSpeciesPtrVect m_fluid_species;
   ScalarPtrVect            m_scalars;
   GlobalDOF                *m_global_dof;
};


#include "NamespaceFooter.H"

#endif

