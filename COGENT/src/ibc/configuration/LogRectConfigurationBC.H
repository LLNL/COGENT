#ifndef  _LOGRECTCONFIGURATIONBC_H_
#define  _LOGRECTCONFIGURATIONBC_H_

#include "TimeFunction.H"
#include "GridFunction.H"
#include "ParsingCore.H"
#include "FluidVarBC.H"
#include "EdgeDataBox.H"

#include "NamespaceHeader.H"

/**
 * Boundary condition interface class for fluid quantities in single null geometries
 *
*/
class LogRectConfigurationBC
   : public FluidVarBC
{
   public:

      // Constructor
      LogRectConfigurationBC( const std::string&  species_name,
                              const std::string&  variable_name,
                              const int&          verbosity = 0 );

      /// Destructor.
      /**
       */
      virtual ~LogRectConfigurationBC();

      /// Set the boundary conditions.
      /**
       * Sets the boundary condition on the species for the given geometry
       * at the given time.  Assumes distribution function is in computational
       * coordinates.
       *
       * @param[out] species_phys Fluid species to set.
       * @param[in] time          Time value at which to evaluate the function.
       */
      virtual void apply( FluidSpecies&  species_comp,
                          const Real&    time );
       
      virtual void applyBC( const FluidSpecies&    species_comp,
                            LevelData<FArrayBox>&  a_dst,
                      const Real                   time );
      
      virtual void applyFluxBC( const FluidSpecies&  a_species_comp,
                                LevelData<FluxBox>&  a_dst,
                          const Real                 a_time );
      
      virtual void applyEdgeBC( const FluidSpecies&      a_species_comp,
                                LevelData<EdgeDataBox>&  a_dst,
                          const Real                     a_time );
      
      virtual void setBC( const FluidSpecies&    a_species_comp,
                          LevelData<FArrayBox>&  a_dst,
                    const LevelData<FArrayBox>&  a_src,
                    const Real                   a_time );
      
      virtual void setFluxBC( const FluidSpecies&  a_species_comp,
                              LevelData<FluxBox>&  a_dst,
                        const LevelData<FluxBox>&  a_src,
                        const Real                 a_time );
      
      virtual void setEdgeBC( const FluidSpecies&      a_species_comp,
                              LevelData<EdgeDataBox>&  a_dst,
                        const LevelData<EdgeDataBox>&  a_src,
                        const Real                     a_time );
      

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
      LogRectConfigurationBC( const LogRectConfigurationBC& );
      LogRectConfigurationBC& operator=( const LogRectConfigurationBC& );

      void fillInflowData( const MagGeom&  a_geometry, 
                           const Real      a_time );
      
      void setAllBcType( const BoundaryBoxLayoutPtrVect&  a_bdry_layout );

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

      inline GridFunction& radialInflowFunc( const Side::LoHiSide& a_side);
      inline GridFunction& poloidalInflowFunc( const Side::LoHiSide& a_side);
#if CFG_DIM==3
      inline GridFunction& toroidalInflowFunc( const Side::LoHiSide& a_side);
#endif
      inline GridFunction& inflowFunc( const int& a_dir,
                                        const Side::LoHiSide& a_side );
      inline std::string radialBcType( const Side::LoHiSide& a_side);
      inline std::string poloidalBcType( const Side::LoHiSide& a_side);
#if CFG_DIM==3
      inline std::string toroidalBcType( const Side::LoHiSide& a_side);
#endif
      inline std::string getBcType( const int& a_dir,
                                    const Side::LoHiSide& a_side );
   
      const std::string m_species_name;
      const std::string m_variable_name;
      int m_verbosity;

  enum {INVALID=-1,
	RADIAL_LOWER, RADIAL_UPPER,
	POLOIDAL_LOWER, POLOIDAL_UPPER,
#if CFG_DIM==3
	TOROIDAL_LOWER, TOROIDAL_UPPER,
#endif
	NUM_INFLOW};

      Vector<RefCountedPtr<GridFunction> > m_inflow_function;
      Vector<std::string> m_bdry_name;
      Vector<std::string> m_bc_type;

      Vector<std::string> m_all_bc_type;
      BoundaryBoxLayoutPtrVect m_all_bdry_layouts;
      BoundaryBoxLayoutPtrVect m_all_bdry_layouts_flux;
      BoundaryBoxLayoutPtrVect m_all_bdry_layouts_edge;
      Vector< RefCountedPtr<LevelData<FArrayBox>> > m_all_bdry_data;
      bool m_all_bdry_defined=false;
      bool m_all_bdry_layouts_flux_defined=false;
      bool m_all_bdry_layouts_edge_defined=false;

};

#include "NamespaceFooter.H"

#endif

