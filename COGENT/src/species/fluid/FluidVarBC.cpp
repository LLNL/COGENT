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
     m_verbosity(a_verbosity),
     m_recycling_coefficient(1.0),
     m_recycling_bc(false)
{
   m_bc_function.resize(a_num_boundaries);
   m_bc_type.resize(a_num_boundaries);
   m_bdry_name.resize(a_num_boundaries);
   
   m_block_bdry_defined=false;
   m_block_bdry_layouts_flux_defined=false;
   m_block_bdry_layouts_edge_defined=false;
   
}


FluidVarBC::~FluidVarBC()
{
   m_bc_function.resize(0);
   m_bdry_name.resize(0);
   m_bc_type.resize(0);
}


void FluidVarBC::apply(FluidSpecies&  a_species_comp,
                       const Real&    a_time )
{
   CH_TIME("FluidVarBC::apply()");
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const MagCoordSys& coord_sys( dynamic_cast<const MagCoordSys&>( *geometry.getCoordSys()) );

   LevelData<FArrayBox>& u( a_species_comp.cell_var( m_variable_name ) );
   const DisjointBoxLayout& grids( u.disjointBoxLayout() );
   const IntVect& ghost_vect( u.ghostVect() );

   if(m_block_bdry_defined==false) {
      FluidBCUtils::defineBoundaryBoxLayouts( m_block_bdry_layouts,
                                              grids,
                                              coord_sys,
                                              ghost_vect );
   
      FluidBCUtils::defineInflowDataStorage( m_block_bdry_data,
                                             m_block_bdry_layouts,
                                             m_variable_name,
                                             a_species_comp );
   
      fillInflowData( geometry, a_time );
      setAllBcType( geometry, m_block_bdry_layouts );
      m_block_bdry_defined = true;
   }

   const LevelData<FluxBox>& velocity= a_species_comp.velocity();
   FluidBCUtils::setInflowOutflowBC( u,
                                     m_block_bdry_layouts,
                                     m_block_bdry_data,
                                     m_block_bc_type,
                                     geometry,
                                     velocity );
   
   // interpolate all other codim boundaries
   CodimBC::setCodimBoundaryValues( u, coord_sys );
  
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
      LevelData<FArrayBox>& bdry_data( *(m_block_bdry_data[i]) );
      inflow_func->assign( bdry_data, a_geometry, bdry_layout, a_time );
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

void FluidVarBC::applyFluxBC(const FluidSpecies&  a_species_comp,
                             LevelData<FluxBox>&  a_dst,
                             const Real           a_time )
{
   CH_TIME("FluidVarBC::applyFluxBC()");
      
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const MagCoordSys& coord_sys( dynamic_cast<const MagCoordSys&>( *geometry.getCoordSys()) );
  
   const DisjointBoxLayout& grids( a_dst.disjointBoxLayout() );
   const IntVect& ghost_vect( a_dst.ghostVect() );
   
   if(m_block_bdry_layouts_flux_defined==false) {
      FluidBCUtils::defineBoundaryBoxLayouts( m_block_bdry_layouts_flux,
                                              grids,
                                              coord_sys,
                                              ghost_vect );
      m_block_bdry_layouts_flux_defined = true;
   }
   setAllBcType( geometry, m_block_bdry_layouts_flux );

   FluidBCUtils::setFluxBC( a_dst,
                            m_block_bdry_layouts_flux,
                            m_block_bc_type,
                            geometry ) ;
   
}

void FluidVarBC::applyEdgeBC(const FluidSpecies&      a_species_comp,
                             LevelData<EdgeDataBox>&  a_dst,
                             const Real               a_time )
{
   CH_TIME("FluidVarBC::applyEdgeBC()");
   
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const MagCoordSys& coord_sys( dynamic_cast<const MagCoordSys&>( *geometry.getCoordSys()) );
   
   const DisjointBoxLayout& grids( a_dst.disjointBoxLayout() );
   const IntVect& dst_ghost_vect( a_dst.ghostVect() );
   if(dst_ghost_vect > IntVect::Zero) {

      if(m_block_bdry_layouts_edge_defined==false) {
         FluidBCUtils::defineBoundaryBoxLayouts( m_block_bdry_layouts_edge,
                                                 grids,
                                                 coord_sys,
                                                 dst_ghost_vect );
         
         m_block_bdry_layouts_edge_defined = true;
      }
      setAllBcType( geometry, m_block_bdry_layouts_edge );

      FluidBCUtils::setEdgeBC( a_dst,
                               m_block_bdry_layouts_edge,
                               m_block_bc_type );

   }

}

void FluidVarBC::setFluxBC(const FluidSpecies&        a_species_comp,
                           LevelData<FluxBox>&        a_dst,
                           const LevelData<FluxBox>&  a_src,
                           const Real                 a_time )
{
   CH_TIME("FluidVarBC::setFluxBC()");
      
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const MagCoordSys& coord_sys( dynamic_cast<const MagCoordSys&>( *geometry.getCoordSys()) );
   
   const DisjointBoxLayout& grids( a_dst.disjointBoxLayout() );
   const IntVect& ghost_vect( a_dst.ghostVect() );
   
   if(m_block_bdry_layouts_flux_defined==false) {
      FluidBCUtils::defineBoundaryBoxLayouts( m_block_bdry_layouts_flux,
                                              grids,
                                              coord_sys,
                                              ghost_vect );
      m_block_bdry_layouts_flux_defined = true;
   }
   setAllBcType( geometry, m_block_bdry_layouts_flux );

   FluidBCUtils::setFluxBC( a_dst,
                            m_block_bdry_layouts_flux,
                            a_src ) ;
   
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

void FluidVarBC::setEdgeBC(const FluidSpecies&            a_species_comp,
                           LevelData<EdgeDataBox>&        a_dst,
                           const LevelData<EdgeDataBox>&  a_src,
                           const Real                     a_time )
{
   CH_TIME("FluidVarBC::setEdgeBC()");
   
   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const MagCoordSys& coord_sys( dynamic_cast<const MagCoordSys&>( *geometry.getCoordSys()) );
   
   const DisjointBoxLayout& grids( a_dst.disjointBoxLayout() );
   const IntVect& ghost_vect( a_dst.ghostVect() );
   
   if(m_block_bdry_layouts_edge_defined==false) {
      FluidBCUtils::defineBoundaryBoxLayouts( m_block_bdry_layouts_edge,
                                              grids,
                                              coord_sys,
                                              ghost_vect );
   
      m_block_bdry_layouts_edge_defined = true;
   }
   setAllBcType( geometry, m_block_bdry_layouts_edge );

   FluidBCUtils::setEdgeBC( a_dst,
                            m_block_bdry_layouts_edge,
                            a_src );
   
}

void FluidVarBC::applyBC(const FluidSpecies&    a_species_comp,
                         LevelData<FArrayBox>&  a_dst,
                         const Real             a_time )
{
   CH_TIME("FluidVarBC::applyBC()");

   const MagGeom& geometry( a_species_comp.configurationSpaceGeometry() );
   const MagCoordSys& coord_sys( dynamic_cast<const MagCoordSys&>( *geometry.getCoordSys()) );
   
   const DisjointBoxLayout& grids( a_dst.disjointBoxLayout() );
   const IntVect& ghost_vect( a_dst.ghostVect() );

   if(m_block_bdry_defined==false) {
      FluidBCUtils::defineBoundaryBoxLayouts( m_block_bdry_layouts,
                                              grids,
                                              coord_sys,
                                              ghost_vect );
   
      setAllBcType( geometry, m_block_bdry_layouts );
      m_block_bdry_defined = true;
   }
   setAllBcType( geometry, m_block_bdry_layouts );
   
   FluidBCUtils::setBC( a_dst,
                        m_block_bdry_layouts,
                        m_block_bc_type,
                        coord_sys );

}


#include "NamespaceFooter.H"
