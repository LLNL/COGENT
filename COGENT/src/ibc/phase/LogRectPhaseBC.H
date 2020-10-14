#ifndef  _LOGRECTPHASEBC_H_
#define  _LOGRECTPHASEBC_H_

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

#include "NamespaceHeader.H"

/**
 * Boundary condition interface class for phase space quantities (distribution functions)
 *
*/
class LogRectPhaseBC
   : public KineticSpeciesBC
{
   public:

      // Constructor
      LogRectPhaseBC( const std::string& name,
                      ParmParse& pp,
                      const int& verbosity = 0 );

      /// Destructor.
      /**
       */
      virtual ~LogRectPhaseBC();

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
      LogRectPhaseBC( const LogRectPhaseBC& );
      LogRectPhaseBC& operator=( const LogRectPhaseBC& );

      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& a_pp );

      KineticFunction& radialInflowFunc( const Side::LoHiSide& a_side );

      KineticFunction& poloidalInflowFunc( const Side::LoHiSide& a_side );

#if CFG_DIM==3
      KineticFunction& toroidalInflowFunc( const Side::LoHiSide& a_side );
#endif

      KineticFunction& vParallelInflowFunc( const Side::LoHiSide& a_side );

      KineticFunction& muInflowFunc( const Side::LoHiSide& a_side );

      KineticFunction& inflowFunc( const int& a_dir,
                                    const Side::LoHiSide& a_side );

      void fillInflowData(KineticSpecies& species_comp,
                          const CFG::LevelData<CFG::FArrayBox>& phi,
                          const LevelData<FluxBox>& velocity,
                          const Real& time);
   
      void setAllBcType( const BoundaryBoxLayoutPtrVect&  a_bdry_layout );
   
      void applySheathBC( KineticSpeciesPtrVect& species,
                          const int& species_index,
                          const CFG::LevelData<CFG::FArrayBox>& phi,
                          const LevelData<FluxBox>& velocity);
   
      void applyFluxBC( KineticSpecies& species,
                       const CFG::LevelData<CFG::FArrayBox>& phi,
                       const LevelData<FluxBox>& velocity,
                       const Real& time);

      inline std::string radialBcType( const Side::LoHiSide& a_side );
      inline std::string poloidalBcType( const Side::LoHiSide& a_side );

#if CFG_DIM==3
      inline std::string toroidalBcType( const Side::LoHiSide& a_side );
#endif

      inline std::string vParallelBcType( const Side::LoHiSide& a_side );
      inline std::string muBcType( const Side::LoHiSide& a_side );
   
      inline std::string getBcType( const int& a_dir,
                                    const Side::LoHiSide& a_side );

      inline std::string getBcName( const int& a_dir,
                                    const Side::LoHiSide& a_side );

      const std::string m_name;
      int m_verbosity;

      enum {INVALID=-1,
            RADIAL_LOWER, RADIAL_UPPER,
            POLOIDAL_LOWER, POLOIDAL_UPPER,
#if CFG_DIM==3
            TOROIDAL_LOWER, TOROIDAL_UPPER,
#endif
            VPAR_LOWER, VPAR_UPPER,
            MU_LOWER, MU_UPPER,
            NUM_INFLOW};

      Vector<std::string> m_all_bc_type;
      BoundaryBoxLayoutPtrVect m_all_bdry_layouts;
      KineticSpeciesPtrVect m_all_bdry_data;
      Vector<RefCountedPtr<KineticFunction> > m_inflow_function;
      Vector<std::string> m_bc_type;
      Vector<std::string> m_bdry_name;
      bool m_all_bdry_defined;
      bool m_logical_sheath_bc;
      bool m_flux_bc;
      ParmParse m_pp;
   
 };

#include "NamespaceFooter.H"

#endif
