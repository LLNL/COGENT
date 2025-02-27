#ifndef  _SINGLENULLPHASEBC_H_
#define  _SINGLENULLPHASEBC_H_

#include "PhaseBlockCoordSys.H"

#include "BoundaryBoxLayout.H.multidim"
#include "REAL.H"
#include "RefCountedPtr.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "ParmParse.H"
#include "KineticFunction.H"
#include "KineticSpeciesBC.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM VEL_DIM
#include "CartesianCS.H"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

/**
 * Boundary condition interface class for phase space quantities (distribution functions)
 *
*/
class SingleNullPhaseBC
   : public KineticSpeciesBC
{
   public:

      // Constructor
      SingleNullPhaseBC( const std::string& name,
                         ParmParse& pp,
                         const int& verbosity = 0 );

      /// Destructor.
      /**
       */
      virtual ~SingleNullPhaseBC();

      /// Set the boundary conditions.
      /**
       * Sets the boundary condition on the species for the given geometry
       * at the given time.  Assumes distribution function is in computational
       * coordinates.
       *
       * @param[out] species  Physical kinetic species vector
       * @param[out] species_index  Species to set (index)
       * @param[in] phi           Potential field.
       * @param[in] velocity      Phase space velocity field.
       * @param[in] time          Time value at which to evaluate the function.
       */
      virtual void apply( KineticSpeciesPtrVect& species,
                          const int& species_index,
                          const CFG::LevelData<CFG::FArrayBox>& phi,
                          const LevelData<FluxBox>& velocity,
                          const Real& time );

      /// Check association with variable
      /**
       * Returns true if this boundary condition corresponds to the variable
       * identified by the name parameter.  This assumes that all derived
       * types store the name of the variable to which they apply.
       *
       * @param[in] name String identifier of the name.
       */
      virtual bool isForVariable( const std::string& a_name ) const
      {
         return (m_name==a_name);
      }

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      SingleNullPhaseBC( const SingleNullPhaseBC& );
      SingleNullPhaseBC& operator=( const SingleNullPhaseBC& );

      void fillInflowData( const Real& a_time );
   
//      void createReflectedData( LevelData<FArrayBox>&       a_rflct_data,
//                              AuxDataMap&                 a_rflct_data_map,
//                              const LevelData<FArrayBox>& a_soln,
//                              const MultiBlockCoordSys&   a_coord_sys,
//                              const Side::LoHiSide&       a_side ) const;

//      void applyPoloidalReflections( BoundaryBoxData&            a_inflow_data,
//                                   const LevelData<FArrayBox>& a_soln,
//                                   const MultiBlockLevelGeom&  a_geometry,
//                                   const BCAuxData&            a_aux_data ) const;

      /// Parse the input database for parameters.
      /**
       */
   
      void applySheathBC(KineticSpeciesPtrVect& species,
                         const int& species_index,
                         const CFG::LevelData<CFG::FArrayBox>& a_phi,
                         const LevelData<FluxBox>& a_velocity);
   
      void applyFluxBC( KineticSpecies& species,
                       const CFG::LevelData<CFG::FArrayBox>& phi,
                       const LevelData<FluxBox>& velocity,
                       const Real& time);
   
      inline void parseParameters( ParmParse& a_pp );

      inline KineticFunction& radialInflowFunc( const Side::LoHiSide& a_side,
                                                 const int& a_block_type );

#if CFG_DIM==3
      inline KineticFunction& toroidalInflowFunc( const Side::LoHiSide& a_side );
#endif

      inline KineticFunction& poloidalInflowFunc( const Side::LoHiSide& a_side );

      inline KineticFunction& vParallelInflowFunc( const Side::LoHiSide& a_side );

      inline KineticFunction& muInflowFunc( const Side::LoHiSide& a_side );

      inline KineticFunction& inflowFunc( const int& a_dir,
                                           const Side::LoHiSide& a_side );

      void setAllBcType( const BoundaryBoxLayoutPtrVect&  a_bdry_layout );
   
      inline std::string radialBcType( const Side::LoHiSide& a_side,
                                       const int& a_block_type );
   
#if CFG_DIM==3
      inline std::string toroidalBcType( const Side::LoHiSide& a_side );
#endif
   
      inline std::string poloidalBcType( const Side::LoHiSide& a_side );
   
      inline std::string vParallelBcType( const Side::LoHiSide& a_side );
   
      inline std::string muBcType( const Side::LoHiSide& a_side );
   
      inline std::string getBcType( const int& a_dir,
                                    const Side::LoHiSide& a_side );
   

      const std::string m_name;
      int m_verbosity;
      bool m_logical_sheath;
      bool m_flux_bc;

      enum {INVALID=-1,
            RADIAL_CORE, RADIAL_SOL, RADIAL_PF,
#if CFG_DIM==3
            TOROIDAL_LOWER, TOROIDAL_UPPER,
#endif
            POLOIDAL_INNER_DIV, POLOIDAL_OUTER_DIV,
            VPAR_LOWER, VPAR_UPPER,
            MU_LOWER, MU_UPPER,
            NUM_INFLOW};

      Vector<std::string> m_all_bc_type;
      BoundaryBoxLayoutPtrVect m_all_bdry_layouts;
      KineticSpeciesPtrVect m_all_bdry_data;
      Vector<RefCountedPtr<KineticFunction> > m_inflow_function;
      Vector<std::string> m_bdry_name;
      Vector<std::string> m_bc_type;
      bool m_all_bdry_defined;
      ParmParse m_pp;
};

#include "NamespaceFooter.H"

#endif

