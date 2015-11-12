#ifndef  _ANNULUSPHASEBC_H_
#define  _ANNULUSPHASEBC_H_

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
class AnnulusPhaseBC
   : public KineticSpeciesBC
{
   public:

      // Constructor
      AnnulusPhaseBC( const std::string& name,
                      ParmParse& pp,
                      const int& verbosity = 0 );

      /// Destructor.
      /**
       */
      virtual ~AnnulusPhaseBC();

      /// Set the boundary conditions.
      /**
       * Sets the boundary condition on the species for the given geometry
       * at the given time.  Assumes distribution function is in computational
       * coordinates.
       *
       * @param[out] species_phys Kinetic species to set.
       * @param[in] phi           Potential field.
       * @param[in] velocity      Phase space velocity field.
       * @param[in] time          Time value at which to evaluate the function.
       */
      virtual void apply( KineticSpecies& species_comp,
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
      AnnulusPhaseBC( const AnnulusPhaseBC& );
      AnnulusPhaseBC& operator=( const AnnulusPhaseBC& );

      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& a_pp );

      KineticFunction& radialInflowFunc( const Side::LoHiSide& a_side );

      KineticFunction& vParallelInflowFunc( const Side::LoHiSide& a_side );

      KineticFunction& muInflowFunc( const Side::LoHiSide& a_side );

      KineticFunction& inflowFunc( const int& a_dir,
                                    const Side::LoHiSide& a_side );

      void fillInflowData( KineticSpeciesPtrVect& a_bdry_data,
                           const BoundaryBoxLayoutPtrVect& a_bdry_layout,
                           const Real& a_time );

      const std::string m_name;
      int m_verbosity;

      enum {INVALID=-1,
            RADIAL_INNER, RADIAL_OUTER,
            VPAR_LOWER, VPAR_UPPER,
            MU_LOWER, MU_UPPER,
            NUM_INFLOW};

      Vector<RefCountedPtr<KineticFunction> > m_inflow_function;
      Vector<std::string> m_bdry_name;
};

#include "NamespaceFooter.H"

#endif
