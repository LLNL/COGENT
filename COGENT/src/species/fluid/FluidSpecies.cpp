#include "FluidSpecies.H"

#include "NamespaceHeader.H"


FluidSpecies::FluidSpecies( const string&     a_pp_prefix,
                            const string&     a_name,
                            const MagGeom&    a_geometry,
                            const IntVect&    a_ghost_vect )
   : CFGVars(a_pp_prefix, a_name, a_geometry)
{
   ParmParse pp(a_pp_prefix.c_str());

   if ( pp.contains("mass") ) {
      pp.get("mass", m_mass);
   }
   else {
      MayDay::Error("FluidSpecies::FluidSpecies(): mass not specified");
   }

   if ( pp.contains("charge")) {
      pp.get("charge", m_charge);
   }
   else {
      MayDay::Error("FluidSpecies::FluidSpecies(): charge not specified");
   }

   if ( pp.contains("operator_type") ) {
      pp.get("operator_type", m_op_type);
   }
   else {
      MayDay::Error("FluidSpecies::FluidSpecies(): operator_type not specified");
   }
      
   addCellVar("density", 1, a_ghost_vect);

   m_velocity.define(a_geometry.gridsFull(), SpaceDim, a_ghost_vect);
   for (DataIterator dit(m_velocity.dataIterator()); dit.ok(); ++dit) {
      m_velocity[dit].setVal(0.);
   }
}


FluidSpecies::FluidSpecies( const FluidSpecies& a_foo )
   : CFGVars(a_foo.pp_prefix(), a_foo.name(), a_foo.configurationSpaceGeometry()),
     m_op_type(a_foo.m_op_type),
     m_mass(a_foo.m_mass),
     m_charge(a_foo.m_charge)
{
   for (int n=0; n<num_cell_vars(); ++n) {
      const LevelData<FArrayBox>& data = a_foo.cell_var(n);
      const string& name = a_foo.cell_var_name(n);
      addCellVar(name, data.nComp(), data.ghostVect());
   }

   for (int n=0; n<num_face_vars(); ++n) {
      const LevelData<FluxBox>& data = a_foo.face_var(n);
      const string& name = a_foo.face_var_name(n);
      addFaceVar(name, data.nComp(), data.ghostVect());
   }
}


void FluidSpecies::numberDensity( LevelData<FArrayBox>& a_rho ) const
{
   massDensity(a_rho);
   for (DataIterator dit(a_rho.dataIterator()); dit.ok(); ++dit ) {
      a_rho[dit].divide( m_mass );
   }
}


void FluidSpecies::massDensity( LevelData<FArrayBox>& a_rho ) const
{
   DisjointBoxLayout dbl( cell_var(0).disjointBoxLayout() );
   CH_assert( dbl.compatible( a_rho.disjointBoxLayout() ) );
   Interval src_interval( MASS_DENSITY, MASS_DENSITY );
   Interval dst_interval( 0 , 0 );
   cell_var(0).copyTo( src_interval, a_rho, dst_interval );
}


void FluidSpecies::chargeDensity( LevelData<FArrayBox>& a_rho ) const
{
   massDensity( a_rho );
   Real ratio( m_charge / m_mass );
   for (DataIterator dit( a_rho.dataIterator()); dit.ok(); ++dit ) {
      a_rho[dit].mult( ratio );
   }
}


void FluidSpecies::pressure( LevelData<FArrayBox>& a_pressure ) const
{
}


void FluidSpecies::temperature( LevelData<FArrayBox>& a_temperature ) const
{
}


void FluidSpecies::velocity( LevelData<FArrayBox>& a_velocity ) const
{
}


bool FluidSpecies::isSpecies( const string& a_name ) const
{
   if (name() == a_name) return true;
   return false;
}


const FluidSpecies& FluidSpecies::operator=( const FluidSpecies& a_rhs )
{
   if (&a_rhs != this)
   {
      CH_assert( (cell_var(0).disjointBoxLayout()).compatible( a_rhs.cell_var(0).disjointBoxLayout() ) );
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;
      m_op_type = a_rhs.m_op_type;
      cell_var(0).define( a_rhs.cell_var(0) );
   }
   return *this;
}


void FluidSpecies::copy( const FluidSpecies& a_rhs )
{
   if (&a_rhs != this) {
      CFGVars::copy(a_rhs);

      m_op_type = a_rhs.m_op_type;
      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;
   }
}


Real FluidSpecies::maxValue() const
{
   const LevelData<FArrayBox>& density = cell_var(0);

   const DisjointBoxLayout& grids(density.getBoxes());
   Real local_maximum( -CH_BADVAL );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      Box box( grids[dit] );
      Real box_max( density[dit].max( box ) );
      local_maximum = Max( local_maximum, box_max );
   }

   Real maximum( local_maximum );
#ifdef CH_MPI
   MPI_Allreduce( &local_maximum, &maximum, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif

   return maximum;
}


Real FluidSpecies::minValue() const
{
   const LevelData<FArrayBox>& density = cell_var(0);

   const DisjointBoxLayout& grids(density.getBoxes());
   Real local_minimum( CH_BADVAL );
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      Box box( grids[dit] );
      Real box_min(density[dit].min( box ) );
      if (box_min<0.0) {
         IntVect box_min_loc = density[dit].minIndex( box );
         pout() << box_min_loc << ":  " << box_min << endl;
      }
      local_minimum = Min( local_minimum, box_min );
   }
   Real minimum(local_minimum);
#ifdef CH_MPI
   MPI_Allreduce( &local_minimum, &minimum, 1, MPI_CH_REAL, MPI_MIN, MPI_COMM_WORLD );
#endif
   return minimum;
}


RefCountedPtr<CFGVars>
FluidSpecies::convertToPhysical( const IntVect& a_ghost_vect ) const
{
   RefCountedPtr<FluidSpecies> physical_ptr = clone(a_ghost_vect, true);

   LevelData<FArrayBox>& cell_data = physical_ptr->cell_var(0);

   if ( cell_data.isDefined() ) {
      configurationSpaceGeometry().divideJonValid(cell_data);
      cell_data.exchange();
   }

   return physical_ptr;
}


void
FluidSpecies::convertFromPhysical()
{
   if ( cell_var(0).isDefined() ) {
      configurationSpaceGeometry().multJonValid(cell_var(0));
      cell_var(0).exchange();
   }
}
   
#if 0
int
FluidSpecies::cell_component(const string& a_variable_name) const
{
   int component = -1;

   for (int n=0; n<num_cell_vars(); ++n) {
      if ( a_variable_name == cell_var_name(n) ) {
         component = n;
         break;
      }
   }

   if ( component == -1 ) {
      MayDay::Error("FluidSpecies::cell_component(): variable_name not found");
   }

   return component;
}
#endif


RefCountedPtr<CFGVars>
FluidSpecies::clone( const IntVect& a_ghost_vect,
                     const bool     a_copy_soln_data ) const
{
   RefCountedPtr<FluidSpecies> result
      = RefCountedPtr<FluidSpecies>(
            new FluidSpecies( pp_prefix(), name(), configurationSpaceGeometry(), a_ghost_vect ) );

   if (a_copy_soln_data) {
      const LevelData<FArrayBox>& this_density = cell_var(0);
      LevelData<FArrayBox>& result_density(result->cell_var(0));

      for (DataIterator dit(result_density.dataIterator() ); dit.ok(); ++dit) {
         result_density[dit].copy(this_density[dit] );
      }
   }

   return result;
}


#include "NamespaceFooter.H"
