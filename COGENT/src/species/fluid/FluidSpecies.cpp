#include "FluidSpecies.H"
#include "SpaceUtils.H"

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
 
   if ( pp.contains("momentumDensity") ) {
      pp.get("momentumDensity", m_evolve_momentumDensity);
   }
   if ( pp.contains("energyDensity") ) {
      pp.get("energyDensity", m_evolve_energyDensity);
   }
   if ( pp.contains("energyDensity_ele") ) {
      pp.get("energyDensity_ele", m_evolve_energyDensity_ele);
   }
   if ( pp.contains("momentumDensity_virtual") ) {
      pp.get("momentumDensity_virtual", m_evolve_momentumDensity_virtual);
   }
   if ( pp.contains("magneticField") ) {
      pp.get("magneticField", m_evolve_magneticField);
   }
   if ( pp.contains("magneticField_virtual") ) {
      pp.get("magneticField_virtual", m_evolve_magneticField_virtual);
   }
   if ( pp.contains("electricField") ) {
      pp.get("electricField", m_evolve_electricField);
   }
   //if ( pp.contains("electricField_virtual") ) {
   //   pp.get("electricField_virtual", m_evolve_electricField_virtual);
   //}
   if ( pp.contains("currentDensity") ) {
      pp.get("currentDensity", m_evolve_currentDensity);
   }
   //if ( pp.contains("currentDensity_virtual") ) {
   //   pp.get("currentDensity_virtual", m_evolve_currentDensity_virtual);
   //}
   if ( pp.contains("plotMemberVars") ) {
      pp.get("plotMemberVars", m_plotMemberVars);
      if(m_plotMemberVars) {
         int n(0);
         pp.query( "plotMemberVars.number", n );
         m_plotMemberVarNames.resize( n );
         pp.queryarr( "plotMemberVars.list", m_plotMemberVarNames, 0, n );
      }
   }
   if ( pp.contains("adiabatic_coefficient") ) {
      pp.get("adiabatic_coefficient", m_gamma);
   }
   
   addCellVar("density", 1, a_ghost_vect);
   if(m_evolve_momentumDensity)   addCellVar("momentumDensity", SpaceDim, a_ghost_vect);
   if(m_evolve_energyDensity)     addCellVar("energyDensity",     1, a_ghost_vect);
   if(m_evolve_energyDensity_ele) addCellVar("energyDensity_ele", 1, a_ghost_vect);
   if(m_evolve_momentumDensity_virtual)  addCellVar("momentumDensity_virtual", 1, a_ghost_vect);
   if(m_evolve_magneticField)    addFaceVar("magneticField", 1, a_ghost_vect);
   if(m_evolve_magneticField_virtual)    addCellVar("magneticField_virtual", 1, a_ghost_vect);
   if(m_evolve_electricField)    addEdgeVar("electricField", 1, a_ghost_vect);
   //if(m_evolve_electricField_virtual)    addNodeVar("electricField_virtual", 1, a_ghost_vect);
   if(m_evolve_currentDensity)   addEdgeVar("currentDensity", 1, a_ghost_vect);
   //if(m_evolve_currentDensity_virtual)   addNodeVar("currentDensity_virtual", 1, a_ghost_vect);


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
   
   for (int n=0; n<num_edge_vars(); ++n) {
      const LevelData<EdgeDataBox>& data = a_foo.edge_var(n);
      const string& name = a_foo.edge_var_name(n);
      addEdgeVar(name, data.nComp(), data.ghostVect());
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
   //Interval src_interval( MASS_DENSITY, MASS_DENSITY );
   //Interval dst_interval( 0 , 0 );
   //cell_var(0).copyTo( src_interval, a_rho, dst_interval );
   //cell_var(0).copyTo( a_rho ); 
   // JRA, copyTo doesnt fill ghost unless passing Copier
   // the following does ghost cells. 
   for (DataIterator dit( a_rho.dataIterator()); dit.ok(); ++dit ) {
      //a_rho[dit].copy( cell_var(0)[dit] );
      a_rho[dit].copy( cell_var("density")[dit] );
   }
}


void FluidSpecies::chargeDensity( LevelData<FArrayBox>& a_rho ) const
{
   massDensity( a_rho );
   Real ratio( m_charge / m_mass );
   for (DataIterator dit( a_rho.dataIterator()); dit.ok(); ++dit ) {
      a_rho[dit].mult( ratio );
   }
}

void FluidSpecies::velocity_virtual( LevelData<FArrayBox>& a_velocity ) const
{
   CH_TIME("FluidSpecies::velocity_virtual()");
   
   if(m_evolve_momentumDensity_virtual==1 ) {
      DisjointBoxLayout dbl( cell_var("momentumDensity_virtual").disjointBoxLayout() );
      CH_assert( dbl.compatible( a_velocity.disjointBoxLayout() ) );
      LevelData<FArrayBox> a_massDensity(dbl, cell_var("density").nComp(), a_velocity.ghostVect());
      //LevelData<FArrayBox> a_massDensity;
      //a_massDensity.define( cell_var("density") );
      massDensity( a_massDensity );
      for (DataIterator dit( a_velocity.dataIterator()); dit.ok(); ++dit ) {
         a_velocity[dit].copy( cell_var("momentumDensity_virtual")[dit] );
         a_velocity[dit].divide( a_massDensity[dit] );
      }
   }
   else {
      for (DataIterator dit( a_velocity.dataIterator()); dit.ok(); ++dit ) {
         a_velocity[dit].setVal( 0.0 );
      }
   }
}

void FluidSpecies::meanEnergy( LevelData<FArrayBox>& a_meanEnergy ) const
{
   CH_TIME("FluidSpecies::meanEnergy()");
   
   DisjointBoxLayout dbl( cell_var(0).disjointBoxLayout() );
   CH_assert( dbl.compatible( a_meanEnergy.disjointBoxLayout() ) );
   for (DataIterator dit( a_meanEnergy.dataIterator()); dit.ok(); ++dit ) {
      a_meanEnergy[dit].setVal( 0.0 );
   }

   
   if(m_evolve_momentumDensity) {
      LevelData<FArrayBox> a_vel(dbl, cell_var("momentumDensity").nComp(), a_meanEnergy.ghostVect());
      //LevelData<FArrayBox> a_vel;
      //a_vel.define( cell_var("momentumDensity") );
      velocity( a_vel );
      for (DataIterator dit( a_meanEnergy.dataIterator()); dit.ok(); ++dit ) {
         a_vel[dit].mult( cell_var("momentumDensity")[dit],0,0,SpaceDim );
         a_vel[dit].mult( 0.5 );
         for (int n=0; n<SpaceDim; ++n) {
            a_meanEnergy[dit].plus( a_vel[dit],n,0,1 );
         }
      } 
   } 
   if(m_evolve_momentumDensity_virtual) {
      LevelData<FArrayBox> a_vel_virtual(dbl, cell_var("momentumDensity_virtual").nComp(), a_meanEnergy.ghostVect());
      //LevelData<FArrayBox> a_vel_virtual;
      //a_vel_virtual.define( cell_var("momentumDensity_virtual") );
      velocity_virtual( a_vel_virtual );
      for (DataIterator dit( a_meanEnergy.dataIterator()); dit.ok(); ++dit ) {
         a_vel_virtual[dit].mult( cell_var("momentumDensity_virtual")[dit] );
         a_vel_virtual[dit].mult( 0.5 );
         a_meanEnergy[dit].plus( a_vel_virtual[dit] );
      } 
   } 

}

void FluidSpecies::pressure( LevelData<FArrayBox>& a_pressure ) const
{
   CH_TIME("FluidSpecies::pressure()");

   if( m_evolve_energyDensity==1 ) {
      DisjointBoxLayout dbl( cell_var("energyDensity").disjointBoxLayout() );
      CH_assert( dbl.compatible( a_pressure.disjointBoxLayout() ) );
   
      LevelData<FArrayBox> a_meanEnergy(dbl, a_pressure.nComp(), a_pressure.ghostVect());
      //LevelData<FArrayBox> a_meanEnergy;
      //a_meanEnergy.define( a_pressure );
      meanEnergy( a_meanEnergy );

      for (DataIterator dit( a_pressure.dataIterator()); dit.ok(); ++dit ) {
         a_pressure[dit].copy( cell_var("energyDensity")[dit] );
         a_pressure[dit].minus( a_meanEnergy[dit] );
         a_pressure[dit].mult( m_gamma-1.0 );
      }
   }
   else {
      for (DataIterator dit( a_pressure.dataIterator()); dit.ok(); ++dit ) {
         a_pressure[dit].setVal( 0.0 );
      }
   }
}

void FluidSpecies::pressure_ele( LevelData<FArrayBox>& a_pressure_ele ) const
{
   CH_TIME("FluidSpecies::pressure_ele()");

   if( m_evolve_energyDensity_ele==1 ) {
      DisjointBoxLayout dbl( cell_var("energyDensity_ele").disjointBoxLayout() );
      CH_assert( dbl.compatible( a_pressure_ele.disjointBoxLayout() ) );
   
      for (DataIterator dit( a_pressure_ele.dataIterator()); dit.ok(); ++dit ) {
         a_pressure_ele[dit].copy( cell_var("energyDensity_ele")[dit] );
         a_pressure_ele[dit].mult( m_gamma-1.0 );
      }
   }
   else {
      for (DataIterator dit( a_pressure_ele.dataIterator()); dit.ok(); ++dit ) {
         a_pressure_ele[dit].setVal( 0.0 );
      }
   }
}

void FluidSpecies::Bpressure( LevelData<FArrayBox>& a_Bpressure ) const
{
   CH_TIME("FluidSpecies::Bpressure()");
   
   if( m_evolve_magneticField_virtual ) {
      DisjointBoxLayout dbl( cell_var("magneticField_virtual").disjointBoxLayout() );
      CH_assert( dbl.compatible( a_Bpressure.disjointBoxLayout() ) );
   
      for (DataIterator dit( a_Bpressure.dataIterator()); dit.ok(); ++dit ) {
         a_Bpressure[dit].copy( cell_var("magneticField_virtual")[dit] );
         a_Bpressure[dit] *= cell_var("magneticField_virtual")[dit];
         a_Bpressure[dit].mult(0.5);
      }
   }
   else {
      for (DataIterator dit( a_Bpressure.dataIterator()); dit.ok(); ++dit ) {
         a_Bpressure[dit].setVal( 0.0 );
      }
   }
}


void FluidSpecies::temperature( LevelData<FArrayBox>& a_temperature ) const
{
}


void FluidSpecies::velocity( LevelData<FArrayBox>& a_velocity ) const
{
   CH_TIME("FluidSpecies::velocity()");
   
   if(m_evolve_momentumDensity) {
      DisjointBoxLayout dbl( cell_var("momentumDensity").disjointBoxLayout() );
      CH_assert( dbl.compatible( a_velocity.disjointBoxLayout() ) );
      LevelData<FArrayBox> a_massDensity(dbl, 1, a_velocity.ghostVect());
      //LevelData<FArrayBox> a_massDensity;
      //a_massDensity.define( cell_var("density") );
      massDensity( a_massDensity );
      for (DataIterator dit( a_velocity.dataIterator()); dit.ok(); ++dit ) {
         a_velocity[dit].copy( cell_var("momentumDensity")[dit] );
         for (int n=0; n<SpaceDim; ++n) {
            a_velocity[dit].divide( a_massDensity[dit],0,n,1 );
         }
      }
   }
   else {
      for (DataIterator dit( a_velocity.dataIterator()); dit.ok(); ++dit ) {
         a_velocity[dit].setVal( 0.0 );
      }
   }
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
      m_mass    = a_rhs.m_mass;
      m_charge  = a_rhs.m_charge;
      m_gamma   = a_rhs.m_gamma;
      m_op_type = a_rhs.m_op_type;
      m_evolve_momentumDensity = a_rhs.m_evolve_momentumDensity;
      m_evolve_energyDensity = a_rhs.m_evolve_energyDensity;
      m_evolve_energyDensity_ele = a_rhs.m_evolve_energyDensity_ele;
      m_evolve_momentumDensity_virtual = a_rhs.m_evolve_momentumDensity_virtual;
      m_evolve_magneticField = a_rhs.m_evolve_magneticField;
      m_evolve_magneticField_virtual = a_rhs.m_evolve_magneticField_virtual;
      m_evolve_electricField = a_rhs.m_evolve_electricField;
      //m_evolve_electricField_virtual = a_rhs.m_evolve_electricField_virtual;
      m_evolve_currentDensity = a_rhs.m_evolve_currentDensity;
      //m_evolve_currentDensity_virtual = a_rhs.m_evolve_currentDensity_virtual;
      for (int n=0; n<a_rhs.num_cell_vars(); ++n) {
         cell_var(n).define( a_rhs.cell_var(n) );
      }
      for (int n=0; n<a_rhs.num_face_vars(); ++n) {
         face_var(n).define( a_rhs.face_var(n) );
      }
      for (int n=0; n<a_rhs.num_edge_vars(); ++n) {
         edge_var(n).define( a_rhs.edge_var(n) );
      }
   }
   return *this;
}

void FluidSpecies::copy( const FluidSpecies& a_rhs )
{
   if (&a_rhs != this) {
      CFGVars::copy(a_rhs);

      m_mass = a_rhs.m_mass;
      m_charge = a_rhs.m_charge;
      m_gamma   = a_rhs.m_gamma;
      m_op_type = a_rhs.m_op_type;
      m_evolve_momentumDensity = a_rhs.m_evolve_momentumDensity;
      m_evolve_energyDensity = a_rhs.m_evolve_energyDensity;
      m_evolve_energyDensity_ele = a_rhs.m_evolve_energyDensity_ele;
      m_evolve_momentumDensity_virtual = a_rhs.m_evolve_momentumDensity_virtual;
      m_evolve_magneticField = a_rhs.m_evolve_magneticField;
      m_evolve_magneticField_virtual = a_rhs.m_evolve_magneticField_virtual;

      m_evolve_electricField = a_rhs.m_evolve_electricField;
      //m_evolve_electricField_virtual = a_rhs.m_evolve_electricField_virtual;
      m_evolve_currentDensity = a_rhs.m_evolve_currentDensity;
      //m_evolve_currentDensity_virtual = a_rhs.m_evolve_currentDensity_virtual;

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

void
FluidSpecies::convertToPhysical()
{

   for (int n=0; n<num_cell_vars(); ++n) {
      if ( cell_var(n).isDefined() ) {
         configurationSpaceGeometry().divideJonValid(cell_var(n));
         cell_var(n).exchange();
      }
   }
   
   for (int n=0; n<num_face_vars(); ++n) {
      if ( face_var(n).isDefined() ) {
         configurationSpaceGeometry().divideJonFaces(face_var(n));
         face_var(n).exchange();
      }
   }
   
   for (int n=0; n<num_edge_vars(); ++n) {
      if ( edge_var(n).isDefined() ) {
         //configurationSpaceGeometry().divideJonEdges(edge_var(n));
         //edge_var(n).exchange();
         SpaceUtils::exchangeEdgeDataBox(edge_var(n)); 
      }
   }


}


void
FluidSpecies::interpFaceVarToCell( LevelData<FArrayBox>&  a_cell_var,
                             const string&                a_face_var_name ) const
{

   CH_assert( face_var(a_face_var_name).isDefined() );
   SpaceUtils::interpFacesToCell(a_cell_var,face_var(a_face_var_name),"c2");
   
}


void
FluidSpecies::interpEdgeVarToCell( LevelData<FArrayBox>&  a_cell_var,
                             const string&                a_edge_var_name ) const
{

   CH_assert( edge_var(a_edge_var_name).isDefined() );
   SpaceUtils::interpEdgesToCell(a_cell_var,edge_var(a_edge_var_name),"c2");
   
}
   
   
void
FluidSpecies::convertFromPhysical()
{

   for (int n=0; n<num_cell_vars(); ++n) {
      if ( cell_var(n).isDefined() ) {
         configurationSpaceGeometry().multJonValid(cell_var(n));
         cell_var(n).exchange();
      }
   }
   for (int n=0; n<num_face_vars(); ++n) {
      if ( face_var(n).isDefined() ) {
         configurationSpaceGeometry().multJonFaces(face_var(n));
         face_var(n).exchange();
      }
   }
   for (int n=0; n<num_edge_vars(); ++n) {
      if ( edge_var(n).isDefined() ) {
         //configurationSpaceGeometry().multJonEdges(edge_var(n));
         //edge_var(n).exchange();
      }
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
   CH_TIME("FluidSpecies::clone");
   RefCountedPtr<FluidSpecies> result
      = RefCountedPtr<FluidSpecies>(
            new FluidSpecies( pp_prefix(), name(), configurationSpaceGeometry(), a_ghost_vect ) );

   if (a_copy_soln_data) {

      const LevelData<FArrayBox>& this_density = cell_var(0);
      LevelData<FArrayBox>& result_density(result->cell_var(0));
      for (DataIterator dit(result_density.dataIterator() ); dit.ok(); ++dit) {
         result_density[dit].copy(this_density[dit] );
      }
      if(m_evolve_energyDensity) {
         const LevelData<FArrayBox>& this_energyDensity = cell_var("energyDensity");
         LevelData<FArrayBox>& result_energyDensity(result->cell_var("energyDensity"));
         for (DataIterator dit(result_energyDensity.dataIterator() ); dit.ok(); ++dit) {
            result_energyDensity[dit].copy(this_energyDensity[dit] );
         }
      }
      if(m_evolve_energyDensity_ele) {
         const LevelData<FArrayBox>& this_energyDensity_ele = cell_var("energyDensity_ele");
         LevelData<FArrayBox>& result_energyDensity_ele(result->cell_var("energyDensity_ele"));
         for (DataIterator dit(result_energyDensity_ele.dataIterator() ); dit.ok(); ++dit) {
            result_energyDensity_ele[dit].copy(this_energyDensity_ele[dit] );
         }
      }
      if(m_evolve_momentumDensity) {
         const LevelData<FArrayBox>& this_momentumDensity = cell_var("momentumDensity");
         LevelData<FArrayBox>& result_momentumDensity(result->cell_var("momentumDensity"));
         for (DataIterator dit(result_momentumDensity.dataIterator() ); dit.ok(); ++dit) {
            result_momentumDensity[dit].copy(this_momentumDensity[dit] );
         }
      }
      if(m_evolve_momentumDensity_virtual) {
         const LevelData<FArrayBox>& this_momentumDensity_v = cell_var("momentumDensity_virtual");
         LevelData<FArrayBox>& result_momentumDensity_v(result->cell_var("momentumDensity_virtual"));
         for (DataIterator dit(result_momentumDensity_v.dataIterator() ); dit.ok(); ++dit) {
            result_momentumDensity_v[dit].copy(this_momentumDensity_v[dit] );
         }
      }
      if(m_evolve_magneticField) {
         const LevelData<FluxBox>& this_magneticField = face_var("magneticField");
         LevelData<FluxBox>& result_magneticField(result->face_var("magneticField"));
         for (DataIterator dit(result_magneticField.dataIterator() ); dit.ok(); ++dit) {
            result_magneticField[dit].copy(this_magneticField[dit] );
         }
      }
      if(m_evolve_magneticField_virtual) {
         const LevelData<FArrayBox>& this_magneticField_v = cell_var("magneticField_virtual");
         LevelData<FArrayBox>& result_magneticField_v(result->cell_var("magneticField_virtual"));
         for (DataIterator dit(result_magneticField_v.dataIterator() ); dit.ok(); ++dit) {
            result_magneticField_v[dit].copy(this_magneticField_v[dit] );
         }
      }
      int this_edge_nComp; 
      if(m_evolve_electricField) {
         const LevelData<EdgeDataBox>& this_electricField = edge_var("electricField");
         LevelData<EdgeDataBox>& result_electricField(result->edge_var("electricField"));
         this_edge_nComp = this_electricField.nComp(); 
         for (DataIterator dit(result_electricField.dataIterator() ); dit.ok(); ++dit) {
            result_electricField[dit].copy(this_electricField[dit],0,0,this_edge_nComp );
         }
      }
      if(m_evolve_currentDensity) {
         const LevelData<EdgeDataBox>& this_currentDensity = edge_var("currentDensity");
         LevelData<EdgeDataBox>& result_currentDensity(result->edge_var("currentDensity"));
         this_edge_nComp = this_currentDensity.nComp(); 
         for (DataIterator dit(result_currentDensity.dataIterator() ); dit.ok(); ++dit) {
            result_currentDensity[dit].copy(this_currentDensity[dit],0,0,this_edge_nComp );
         }
      }

   }

   return result;
}


#include "NamespaceFooter.H"
