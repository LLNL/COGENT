#include "FluidVarBC.H"
#include "CodimBC.H"
#include "FluidBCUtils.H"
#include "EllipticOpBC.H"
#include "RecyclingBC.H"

#include "NamespaceHeader.H"


FluidVarBC::FluidVarBC(const std::string&  a_species_name,
                       const std::string&  a_variable_name,
                       const int&          a_verbosity,
                       const int           a_num_boundaries)
   : m_species_name(a_species_name),
     m_variable_name(a_variable_name),
     m_recycling_bc(false),
     m_insulator_conductor_bc(false),
     m_recycling_coefficient(1.0),
     m_verbosity(a_verbosity)
{
   m_bc_function.resize(a_num_boundaries);
   m_bc_type.resize(a_num_boundaries);
   m_bdry_name.resize(a_num_boundaries);
   
   m_block_bdry_defined = false;
   m_inflow_data_defined = false;
   
}


FluidVarBC::~FluidVarBC()
{
   m_bc_function.resize(0);
   m_bdry_name.resize(0);
   m_bc_type.resize(0);
}


void FluidVarBC::apply(FluidSpecies&  a_species_phys,
                       const Real&    a_time )
{
   CH_TIME("FluidVarBC::apply()");
   const MagGeom& geometry( a_species_phys.configurationSpaceGeometry() );
   const MagCoordSys& coord_sys( dynamic_cast<const MagCoordSys&>( *geometry.getCoordSys()) );

   // define boundary Layout, Data, and bc Types
   defineBoundaryLDT( a_species_phys );

   // apply BCs to cell variable
   LevelData<FArrayBox>& cell_var( a_species_phys.cell_var( m_variable_name ) );
   const LevelData<FluxBox>& velocity = a_species_phys.velocity();
   FluidBCUtils::setCellBC( cell_var,
                            velocity,
                            m_block_bdry_layouts,
                            m_block_bdry_data,
                            m_block_bc_type,
                            geometry );
   
   if(m_insulator_conductor_bc) { // apply InsulatorConductorBC if defined
      m_InsulatorConductorBC->applyBC( cell_var, m_block_bdry_layouts, m_block_bc_type, geometry, a_time );
   }
   
//#define NEW_CODIM_CORNER_BC_METHOD
#ifdef NEW_CODIM_CORNER_BC_METHOD
   // interpolate to fill physical corner boundaries
   // JRA, something wrong here for multi-block? (see two_field_neutrals)
   // JRA, Also, code is actually slower in 3D. While time spent here goes
   //      down, for some reason time spent doing exchange in convertToPhysical
   //      goes up by a lot. Don't understand this.
   CodimBC::setCodimCornerValues( cell_var, coord_sys );
#else
   // fills both physical and internal corners with extrapolation
   CodimBC::setCodimBoundaryValues( cell_var, coord_sys );
   
   // copy internal corners with non-corner values via exchange
   geometry.fillCorners(cell_var, cell_var.ghostVect(), SpaceDim); 
#endif
  
}

void FluidVarBC::applyFluxBC(const FluidSpecies&  a_species_phys,
                             LevelData<FluxBox>&  a_dst,
                             const Real           a_time )
{
   CH_TIME("FluidVarBC::applyFluxBC()");
   const MagGeom& geometry( a_species_phys.configurationSpaceGeometry() );
  
   // define boundary Layout, Data, and bc Types
   defineBoundaryLDT( a_species_phys );

   // apply BCs to a_dst
   FluidBCUtils::setFluxBC( a_dst,
                            m_block_bdry_layouts,
                            m_block_bc_type,
                            geometry );

}

void FluidVarBC::applyEdgeBC(const FluidSpecies&      a_species_phys,
                             LevelData<EdgeDataBox>&  a_dst,
                             const Real               a_time )
{
   CH_TIME("FluidVarBC::applyEdgeBC()");
   
   // define boundary Layout, Data, and bc Types
   defineBoundaryLDT( a_species_phys );
   
   const IntVect& dst_ghost_vect( a_dst.ghostVect() );
   if(dst_ghost_vect > IntVect::Zero) {
      FluidBCUtils::setEdgeBC( a_dst,
                               m_block_bdry_layouts,
                               m_block_bc_type );
   }
   
   if(m_insulator_conductor_bc) {
      const MagGeom& geometry( a_species_phys.configurationSpaceGeometry() );
      m_InsulatorConductorBC->applyEdgeBC( a_dst, m_block_bdry_layouts, m_block_bc_type, geometry );
   }

}

void FluidVarBC::setFluxBC(const FluidSpecies&        a_species_phys,
                           LevelData<FluxBox>&        a_dst,
                           const LevelData<FluxBox>&  a_src,
                           const Real                 a_time )
{
   CH_TIME("FluidVarBC::setFluxBC()");
      
   // check that src and dst ghost vects are what they need to be
   const IntVect& dst_ghost_vect( a_dst.ghostVect() );
   const IntVect& src_ghost_vect( a_src.ghostVect() );
   CH_assert(dst_ghost_vect==src_ghost_vect || src_ghost_vect==IntVect::Zero);   

   // define boundary Layout, Data, and bc Types
   defineBoundaryLDT( a_species_phys );
   
   // copy a_src to a_dst in boundary region
   FluidBCUtils::setFluxBC( a_dst,
                            m_block_bdry_layouts,
                            a_src ) ;
   
}

void FluidVarBC::setEdgeBC(const FluidSpecies&            a_species_phys,
                           LevelData<EdgeDataBox>&        a_dst,
                           const LevelData<EdgeDataBox>&  a_src,
                           const Real                     a_time )
{
   CH_TIME("FluidVarBC::setEdgeBC()");
   
   // check that src and dst ghost vects are what they need to be
   const IntVect& dst_ghost_vect( a_dst.ghostVect() );
   const IntVect& src_ghost_vect( a_src.ghostVect() );
   CH_assert(dst_ghost_vect==src_ghost_vect || src_ghost_vect==IntVect::Zero);   
   
   // define boundary Layout, Data, and bc Types
   defineBoundaryLDT( a_species_phys );
   
   // copy a_src to a_dst in boundary region
   FluidBCUtils::setEdgeBC( a_dst,
                            m_block_bdry_layouts,
                            a_src );
   
}

void FluidVarBC::setOnAxisCurlBC( LevelData<EdgeDataBox>&  a_curl_covar,
                            const LevelData<FArrayBox>&    a_By_phys,
                            const FluidSpecies&            a_species_phys )
{
   CH_TIME("FluidVarBC::setOnAxisCurlBC()");
   
   // Compute the appropriate BC on axis for the curl of a function
   // that is odd wrt the axis for axisymmetric systems where the Jacobian is zero.
   //
   // Assumed a_curl is covariant curl(a_By_phys) on cell edges   
   
   const MagGeom& geometry( a_species_phys.configurationSpaceGeometry() );
   
   // define boundary Layout, Data, and bc Types
   defineBoundaryLDT( a_species_phys );

   FluidBCUtils::setOnAxisCurlBC( a_curl_covar,
                                  a_By_phys,
                                  geometry,
                                  m_block_bdry_layouts,
                                  m_block_bc_type );

}

void FluidVarBC::setInsulatorBC( const FluidSpecies&    a_species_comp,
                                 LevelData<FArrayBox>&  a_dst,
                           const LevelData<FArrayBox>&  a_src,
                           const Real                   a_time )
{
   CH_TIME("FluidVarBC::setInsulatorBC()");
  
   // set a_dst equal to a_src in boundary region where an insulator exists
   
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   
   if(m_insulator_conductor_bc) {
      m_InsulatorConductorBC->setInsulatorBC( a_dst, a_src, m_block_bdry_layouts, m_block_bc_type, geometry, a_time );
   }
   
}
      
void FluidVarBC::defineBoundaryLDT( const FluidSpecies&  a_species_phys )
{

   const MagGeom& geometry( a_species_phys.configurationSpaceGeometry() );
   if(m_block_bdry_defined==false) {
   
      const MagCoordSys& coord_sys( dynamic_cast<const MagCoordSys&>( *geometry.getCoordSys()) );
 
      const LevelData<FArrayBox>& u( a_species_phys.cell_var(0) );
      const DisjointBoxLayout& grids( u.disjointBoxLayout() );
      const IntVect& ghost_vect( u.ghostVect() );
      
      FluidBCUtils::defineBoundaryBoxLayouts( m_block_bdry_layouts,
                                              grids,
                                              coord_sys,
                                              ghost_vect );
      
      setAllBcType( geometry, m_block_bdry_layouts );
      m_block_bdry_defined = true;
   }
 
   const bool is_cell_var = a_species_phys.is_cell_var(m_variable_name);
   if(m_inflow_data_defined==false && is_cell_var) {
      FluidBCUtils::defineInflowDataStorage( m_block_bdry_data,
                                             m_block_bdry_layouts,
                                             m_variable_name,
                                             a_species_phys );
      fillInflowData( geometry, 0.0 );
      m_inflow_data_defined = true;
   }

}

void FluidVarBC::fillInflowData( const MagGeom&  a_geometry,
                                 const Real      a_time )
{
   for (int i(0); i<m_block_bdry_layouts.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(m_block_bdry_layouts[i]) );
      const int& block( bdry_layout.block() );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      RefCountedPtr<GridFunction> inflow_func( getBCFunction( block, dir, side ) );
      if(!inflow_func.isNull()) {
         LevelData<FArrayBox>& bdry_data( *(m_block_bdry_data[i]) );
         inflow_func->assign( bdry_data, a_geometry, bdry_layout, a_time );
         bdry_data.exchange();
      }
      else {
         m_block_bdry_data[i] = RefCountedPtr<LevelData<FArrayBox>>(NULL);
      }
   }
}

void FluidVarBC::setAllBcType(const MagGeom&                   a_geometry,
                              const BoundaryBoxLayoutPtrVect&  a_bdry_layouts )
{
   m_block_bc_type.resize(a_bdry_layouts.size());
   for (int i(0); i<a_bdry_layouts.size(); i++) {
      const BoundaryBoxLayout& bdry_layout( *(a_bdry_layouts[i]) );
      const int& block( bdry_layout.block() );
      const int& dir( bdry_layout.dir() );
      const Side::LoHiSide& side( bdry_layout.side() );
      m_block_bc_type[i] = getBCType(block, dir, side);
   }
}

void FluidVarBC::applyRecyclingBC(LevelData<FluxBox>&         a_outflux,
                                  const LevelData<FluxBox>&   a_influx,
                                  const Real                  a_time )
{
   CH_TIME("FluidVarBC::applyRecyclingBC()");
   
   if(m_block_bdry_defined == false) {
      MayDay::Error( "FluidVarBC::applyRecyclingBC: boundary layouts should be defined!" );
   }
   
   for (int i(0); i<m_block_bdry_layouts.size(); i++) {

      const std::string this_bc_type (m_block_bc_type[i]);

      if (this_bc_type == "recycling") {
         
         string pp_prefix = "BC." + m_species_name + "." + m_variable_name;
         ParmParse pp(pp_prefix.c_str());
         
         RecyclingBC recyclingBC(m_block_bdry_layouts[i], pp);
         recyclingBC.applyBC(a_outflux, a_influx, m_recycling_coefficient, a_time);
      }
   }
}

void FluidVarBC::applyRecyclingBC(EllipticOpBC&             a_elliptic_bc,
                                  const LevelData<FluxBox>& a_face_dist,
                                  const Real                a_time )
{
   CH_TIME("FluidVarBC::applyRecyclingBC()");
   
   if(m_block_bdry_defined == false) {
      MayDay::Error( "FluidVarBC::applyRecyclingBC: boundary layouts should be defined!" );
   }
   
   for (int i(0); i<m_block_bdry_layouts.size(); i++) {

      const std::string this_bc_type (m_block_bc_type[i]);

      if (this_bc_type == "recycling") {

         string pp_prefix = "BC." + m_species_name + "." + m_variable_name;
         ParmParse pp(pp_prefix.c_str());

         RecyclingBC recyclingBC(m_block_bdry_layouts[i], pp);
         recyclingBC.applyBC(a_elliptic_bc, a_face_dist, m_recycling_coefficient, a_time);
      }
   }
}

void FluidVarBC::applyRecyclingBC(LevelData<FArrayBox>&       a_dst,
                                  const LevelData<FArrayBox>& a_src,
                                  const Real                  a_time )
{
   CH_TIME("FluidVarBC::applyRecyclingBC()");
   
   if(m_block_bdry_defined == false) {
      MayDay::Error( "FluidVarBC::applyRecyclingBC: boundary layouts should be defined!" );
   }
   
   for (int i(0); i<m_block_bdry_layouts.size(); i++) {

      const std::string this_bc_type (m_block_bc_type[i]);

      if (this_bc_type == "recycling_boundary") {

         string pp_prefix = "BC." + m_species_name + "." + m_variable_name;
         ParmParse pp(pp_prefix.c_str());

         RecyclingBC recyclingBC(m_block_bdry_layouts[i], pp);
         recyclingBC.applyBC(a_dst, a_src, m_recycling_coefficient, a_time);
      }
   }
}

#include "NamespaceFooter.H"
