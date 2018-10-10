#ifndef  _SNCORECONFIGURATIONBC_H_
#define  _SNCORECONFIGURATIONBC_H_


#include "GridFunction.H"
#include "FluidVarBC.H"

#include "NamespaceHeader.H"

/**
 * Boundary condition interface class for fluid quantities in single null geometries
 *
*/
class SNCoreConfigurationBC
   : public FluidVarBC
{
   public:

      // Constructor
      SNCoreConfigurationBC( const std::string&  species_name,
                             const std::string&  variable_name,
                             const int&          verbosity = 0 );

      /// Destructor.
      /**
       */
      virtual ~SNCoreConfigurationBC();

      /// Set the boundary conditions.
      /**
       * Sets the boundary condition on the species for the given geometry
       * at the given time.  Assumes distribution function is in computational
       * coordinates.
       *
       * @param[out] species_phys Fluid species to set.
       * @param[in] time          Time value at which to evaluate the function.
       */
      virtual void apply( FluidSpecies& species_comp,
                          const Real&   time );

      /// Check association with variable
      /**
       * Returns true if this boundary condition corresponds to the variable
       * identified by the variable_name parameter.  This assumes that all derived
       * types store the name of the variable to which they apply.
       *
       * @param[in] variable_name String identifier of the variable name.
       */
      virtual bool isForVariable( const std::string& a_variable_name ) const
      {
         return (m_variable_name==a_variable_name);
      }

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

   private:

      // prohibit copying
      SNCoreConfigurationBC( const SNCoreConfigurationBC& );
      SNCoreConfigurationBC& operator=( const SNCoreConfigurationBC& );

      void fillInflowData( FluidSpeciesPtrVect& a_bdry_data,
                           Vector<std::string>& a_bc_type,
                           const BoundaryBoxLayoutPtrVect& a_bdry_layout,
                           const Real& a_time );
   

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
   
   
      inline void parseParameters();

      inline GridFunction& radialInflowFunc( const Side::LoHiSide& a_side,
                                              const int& a_block_type );

      inline GridFunction& inflowFunc( const int& a_dir,
                                        const Side::LoHiSide& a_side );
   
      inline std::string radialBcType( const Side::LoHiSide& a_side,
                                       const int& a_block_type );
   
      inline std::string getBcType( const int& a_dir,
                                    const Side::LoHiSide& a_side );


      const std::string m_species_name;
      const std::string m_variable_name;
      int m_verbosity;
      bool m_logical_sheath; 

      enum {INVALID=-1,
            RADIAL_INNER,
            RADIAL_OUTER,
            NUM_INFLOW};

      Vector<RefCountedPtr<GridFunction> > m_inflow_function;
      Vector<std::string> m_bdry_name;
      Vector<std::string> m_bc_type;
};

#include "NamespaceFooter.H"

#endif

