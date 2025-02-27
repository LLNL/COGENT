#ifndef _FLUID_VAR_BC_H_
#define _FLUID_VAR_BC_H_

#include "FluidSpecies.H"
#include "EdgeDataBox.H"
#include "GridFunction.H"
#include "GridFunctionLibrary.H"
#include "Directions.H"
#include "InsulatorConductorBC.H"

#include "NamespaceHeader.H"

class FluidVarBC
{
   public:
   
      FluidVarBC(const std::string&  species_name,
                 const std::string&  variable_name,
                 const int&          verbosity,
                 const int           num_boundaries);

      virtual ~FluidVarBC();
   
      /// Set the boundary conditions.
      /**
       * Sets the boundary condition on the species for the given geometry
       * at the given time.  Assumes distribution function is in computational
       * coordinates.
       *
       * @param[out] species_phys Fluid species to set.
       * @param[in] time          Time value at which to evaluate the function.
       */
      
      virtual void apply( FluidSpecies&  species_phys,
                          const Real&    time );

      virtual void applyCellBC( const FluidSpecies&  species_phys,
                                LevelData<FArrayBox>&  dst,
                                const Real&    time );
      
      virtual void applyPhysBC( const FluidSpecies&  species_phys,
                                LevelData<FArrayBox>&  dst );
       
      virtual void applyFluxBC( const FluidSpecies&  species_phys,
                                LevelData<FluxBox>&  dst,
                          const Real                 time );
      
      virtual void applyEdgeBC( const FluidSpecies&      species_phys,
                                LevelData<EdgeDataBox>&  dst,
                          const Real                     time );
      
      virtual void applyNodeBC( const FluidSpecies&        species_phys,
                                LevelData<NodeFArrayBox>&  dst,
                          const Real                       time );      
            
      virtual void setFluxBC( const FluidSpecies&  species_phys,
                              LevelData<FluxBox>&  dst,
                        const LevelData<FluxBox>&  src,
                        const Real                 time );
      
      virtual void setEdgeBC( const FluidSpecies&     species_phys,
                              LevelData<EdgeDataBox>& dst,
                        const LevelData<EdgeDataBox>& src,
                        const Real                    time );
      
      virtual void setNodeBC( const FluidSpecies&     species_phys,
                              LevelData<NodeFArrayBox>& dst,
                        const LevelData<NodeFArrayBox>& src,
                        const Real                    time );
      
      virtual void setOnAxisCurlBC( LevelData<EdgeDataBox>&  a_curl_covar,
                              const LevelData<FArrayBox>&    a_By_phys,
                              const FluidSpecies&            a_species_comp );
      
      virtual void setInsulatorBC( const FluidSpecies&    a_species_comp,
                                   LevelData<FArrayBox>&  a_dst,
                             const LevelData<FArrayBox>&  a_src,
                             const Real                   a_time );
   
      virtual void applyRecyclingBC(LevelData<FluxBox>&         outflux,
                                    const LevelData<FluxBox>&   influx,
                                    const Real                  time );
   
      virtual void applyRecyclingBC(EllipticOpBC&               elliptic_bc,
                                    const LevelData<FluxBox>&   face_val,
                                    const Real                  time );
   
      virtual void applyRecyclingBC(LevelData<FArrayBox>&       dst,
                                    const LevelData<FArrayBox>& src,
                                    const Real                  time );
      
      virtual void defineBoundaryLDT( const FluidSpecies&  a_species_phys );
   
      virtual void fillInflowData(const MagGeom&  geometry,
                                  const Real      time );
   
      virtual void setAllBcType(const MagGeom&                    geometry,
                                const BoundaryBoxLayoutPtrVect&   bdry_layouts );
      
      bool isRecyclingBC() const {return m_recycling_bc;}
      
      bool isInsulatorConductorBC() const {return m_insulator_conductor_bc;}

      /// Check association with variable
      /**
       * Returns true if this boundary condition corresponds to the variable
       * identified by the variable_name parameter.  This assumes that all derived
       * types store the name of the variable to which they apply.
       *
       * @param[in] variable_name String identifier of the variable name.
       */
      bool isForVariable( const std::string& a_variable_name ) const
      {
         return (m_variable_name==a_variable_name);
      }
   
   protected:
   
      virtual std::string getBCType(const int block_number,
                                    const int dir,
                                    const int side) = 0;
   
      virtual RefCountedPtr<GridFunction> getBCFunction(const int block_number,
                                                        const int dir,
                                                        const int side ) = 0;
   
      virtual void setNames() = 0;

   
      const std::string m_species_name;
      const std::string m_variable_name;
   
      Vector<RefCountedPtr<GridFunction> > m_bc_function;
      Vector<std::string> m_bdry_name;
      Vector<std::string> m_bc_type;
      Vector<std::string> m_block_bc_type;
   
      BoundaryBoxLayoutPtrVect m_block_bdry_layouts;
   
      Vector< RefCountedPtr<LevelData<FArrayBox>> > m_block_bdry_data;
      
      bool m_block_bdry_defined;
      bool m_inflow_data_defined;
      bool m_recycling_bc;
      
      RefCountedPtr<InsulatorConductorBC> m_InsulatorConductorBC;
      bool m_insulator_conductor_bc;
   
      Real m_recycling_coefficient;
   
   private:
   
      int m_verbosity;
};


#include "NamespaceFooter.H"

#endif
