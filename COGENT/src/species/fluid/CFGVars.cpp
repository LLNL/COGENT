#include "CFGVars.H"
#include "SpaceUtils.H.multidim"

#include "NamespaceHeader.H"


CFGVars::CFGVars( const string&       a_pp_prefix,
                  const std::string&  a_name,
                  const MagGeom&      a_geometry )
   : m_pp_prefix(a_pp_prefix),
     m_name(a_name),
     m_geometry(a_geometry)
{
   ParmParse pp(a_pp_prefix.c_str());

   pp.query( "plot_mapped_face_vars" , m_plot_mapped_face_vars );
   pp.query( "plot_mapped_edge_vars" , m_plot_mapped_edge_vars );

   if ( pp.contains("plotMemberVars") ) {
      pp.get("plotMemberVars", m_plotMemberVars);
      if(m_plotMemberVars) {
         int n(0);
         pp.query( "plotMemberVars.number", n );
         m_plotMemberVarNames.resize( n );
         pp.queryarr( "plotMemberVars.list", m_plotMemberVarNames, 0, n );
      }
   }
}


const CFGVars& CFGVars::operator=( const CFGVars& a_rhs )
{
   if (&a_rhs != this)
   {
      m_pp_prefix = a_rhs.m_pp_prefix;
      m_name = a_rhs.m_name;

      m_plotMemberVarNames = a_rhs.m_plotMemberVarNames;

      m_cell_data.resize(0);
      for (int i=0; i<a_rhs.num_cell_vars(); ++i) {
         const LevelData<FArrayBox>& data = a_rhs.cell_var(i);
         addCellVar(cell_var_name(i), data.nComp(), data.ghostVect());
         m_cell_data[i]->define(data);
      }

      m_face_data.resize(0);
      for (int i=0; i<a_rhs.num_face_vars(); ++i) {
         const LevelData<FluxBox>& data = a_rhs.face_var(i);
         addFaceVar(face_var_name(i), data.nComp(), data.ghostVect());
         m_face_data[i]->define(data);
      }
      
      m_edge_data.resize(0);
      for (int i=0; i<a_rhs.num_edge_vars(); ++i) {
         const LevelData<EdgeDataBox>& data = a_rhs.edge_var(i);
         addEdgeVar(edge_var_name(i), data.nComp(), data.ghostVect());
         m_edge_data[i]->define(data);
      }

      m_node_data.resize(0);
      for (int i=0; i<a_rhs.num_node_vars(); ++i) {
         const LevelData<NodeFArrayBox>& data = a_rhs.node_var(i);
         addNodeVar(node_var_name(i), data.nComp(), data.ghostVect());
         m_node_data[i]->define(data);
      }

   }
   return *this;
}


void CFGVars::copy( const CFGVars& a_rhs )
{
   if ( &a_rhs != this ) {
      m_pp_prefix = a_rhs.m_pp_prefix;
      m_name = a_rhs.m_name;

      int num_cell_data = a_rhs.num_cell_vars();
      CH_assert(num_cell_data == num_cell_vars());
      for (int i=0; i<num_cell_data; ++i) {
         LevelData<FArrayBox>& this_cell_data = cell_var(i);
         const LevelData<FArrayBox>& that_cell_data = a_rhs.cell_var(i);
         CH_assert( (this_cell_data.disjointBoxLayout()).compatible( that_cell_data.disjointBoxLayout() ) );
         for (DataIterator dit( this_cell_data.dataIterator() ); dit.ok(); ++dit) {
            this_cell_data[dit].copy(that_cell_data[dit] );
         }
      }

      int num_face_data = a_rhs.num_face_vars();
      CH_assert(num_face_data == num_face_vars());
      for (int i=0; i<num_face_data; ++i) {
         LevelData<FluxBox>& this_face_data = face_var(i);
         const LevelData<FluxBox>& that_face_data = a_rhs.face_var(i);
         CH_assert( (this_face_data.disjointBoxLayout()).compatible( that_face_data.disjointBoxLayout() ) );
         for (DataIterator dit( this_face_data.dataIterator() ); dit.ok(); ++dit) {
            this_face_data[dit].copy(that_face_data[dit] );
         }
      }
      
      int num_edge_data = a_rhs.num_edge_vars();
      CH_assert(num_edge_data == num_edge_vars());
      for (int i=0; i<num_edge_data; ++i) {
         LevelData<EdgeDataBox>& this_edge_data = edge_var(i);
         const LevelData<EdgeDataBox>& that_edge_data = a_rhs.edge_var(i);
         CH_assert( (this_edge_data.disjointBoxLayout()).compatible( that_edge_data.disjointBoxLayout() ) );
         for (DataIterator dit( this_edge_data.dataIterator() ); dit.ok(); ++dit) {
            //this_edge_data[dit].copy(that_edge_data[dit],0,0,this_edge_data.nComp() );
            SpaceUtils::copyEdgeDataBox(this_edge_data[dit],that_edge_data[dit]); 
         }
      }

      int num_node_data = a_rhs.num_node_vars();
      CH_assert(num_node_data == num_node_vars());
      for (int i=0; i<num_node_data; ++i) {
         LevelData<NodeFArrayBox>& this_node_data = node_var(i);
         const LevelData<NodeFArrayBox>& that_node_data = a_rhs.node_var(i);
         CH_assert( (this_node_data.disjointBoxLayout()).compatible( that_node_data.disjointBoxLayout() ) );
         for (DataIterator dit( this_node_data.dataIterator() ); dit.ok(); ++dit) {
            this_node_data[dit].getFab().copy(that_node_data[dit].getFab());
         }
      }

   }
}


void CFGVars::zeroData()
{
   for (int i=0; i<num_cell_vars(); ++i) {
      LevelData<FArrayBox>& this_cell_data = cell_var(i);
      for (DataIterator dit( this_cell_data.dataIterator() ); dit.ok(); ++dit) {
         this_cell_data[dit].setVal(0.);
      }
   }

   for (int i=0; i<num_face_vars(); ++i) {
      LevelData<FluxBox>& this_face_data = face_var(i);
      for (DataIterator dit( this_face_data.dataIterator() ); dit.ok(); ++dit) {
         this_face_data[dit].setVal(0.);
      }
   }
   
   for (int i=0; i<num_edge_vars(); ++i) {
      LevelData<EdgeDataBox>& this_edge_data = edge_var(i);
      for (DataIterator dit( this_edge_data.dataIterator() ); dit.ok(); ++dit) {
         this_edge_data[dit].setVal(0.);
      }
   }

   for (int i=0; i<num_node_vars(); ++i) {
      LevelData<NodeFArrayBox>& this_node_data = node_var(i);
      for (DataIterator dit( this_node_data.dataIterator() ); dit.ok(); ++dit) {
         this_node_data[dit].getFab().setVal(0.);
      }
   }

}


void CFGVars::addData( const CFGVars&  a_rhs,
                       const Real      a_factor )
{
   int num_cell_data = a_rhs.num_cell_vars();
   CH_assert(num_cell_data == num_cell_vars());
   for (int i=0; i<num_cell_data; ++i) {
      LevelData<FArrayBox>& this_cell_data = cell_var(i);
      const LevelData<FArrayBox>& that_cell_data = a_rhs.cell_var(i);
      CH_assert( (this_cell_data.disjointBoxLayout()).compatible( that_cell_data.disjointBoxLayout() ) );
      for (DataIterator dit( this_cell_data.dataIterator() ); dit.ok(); ++dit) {
         this_cell_data[dit].plus(that_cell_data[dit], a_factor);
      }
   }

   int num_face_data = a_rhs.num_face_vars();
   CH_assert(num_face_data == num_face_vars());
   for (int i=0; i<num_face_data; ++i) {
      LevelData<FluxBox>& this_face_data = face_var(i);
      const LevelData<FluxBox>& that_face_data = a_rhs.face_var(i);
      CH_assert( (this_face_data.disjointBoxLayout()).compatible( that_face_data.disjointBoxLayout() ) );
      for (DataIterator dit( this_face_data.dataIterator() ); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            this_face_data[dit][dir].plus(that_face_data[dit][dir], a_factor);
         }
      }
   }
   
   int num_edge_data = a_rhs.num_edge_vars();
   CH_assert(num_edge_data == num_edge_vars());
   for (int i=0; i<num_edge_data; ++i) {
      LevelData<EdgeDataBox>& this_edge_data = edge_var(i);
      const LevelData<EdgeDataBox>& that_edge_data = a_rhs.edge_var(i);
      CH_assert( (this_edge_data.disjointBoxLayout()).compatible( that_edge_data.disjointBoxLayout() ) );
      for (DataIterator dit( this_edge_data.dataIterator() ); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            this_edge_data[dit][dir].plus(that_edge_data[dit][dir], a_factor);
         }
      }
   }

   int num_node_data = a_rhs.num_node_vars();
   CH_assert(num_node_data == num_node_vars());
   for (int i=0; i<num_node_data; ++i) {
      LevelData<NodeFArrayBox>& this_node_data = node_var(i);
      const LevelData<NodeFArrayBox>& that_node_data = a_rhs.node_var(i);
      CH_assert( (this_node_data.disjointBoxLayout()).compatible( that_node_data.disjointBoxLayout() ) );
      for (DataIterator dit( this_node_data.dataIterator() ); dit.ok(); ++dit) {
         this_node_data[dit].getFab().plus(that_node_data[dit].getFab(), a_factor);
      }
   }

}


int CFGVars::size( bool a_count_ghosts )
{
   int size(0);

   for (int i=0; i<num_cell_vars(); ++i) {
      LevelData<FArrayBox>& this_cell_data = cell_var(i);
      const DisjointBoxLayout& dbl( this_cell_data.disjointBoxLayout() );
      for (DataIterator dit( this_cell_data.dataIterator() ); dit.ok(); ++dit) {
         const Box& box( a_count_ghosts ? this_cell_data[dit].box() : dbl[dit] );
         size += box.numPts() * this_cell_data.nComp();
      }
   }

   for (int i=0; i<num_face_vars(); ++i) {
      LevelData<FluxBox>& this_face_data = face_var(i);
      const DisjointBoxLayout& dbl( this_face_data.disjointBoxLayout() );
      for (DataIterator dit( this_face_data.dataIterator() ); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            const Box& box( a_count_ghosts ? this_face_data[dit][dir].box() : surroundingNodes(dbl[dit],dir) );
            size += box.numPts() * this_face_data.nComp();
         }
      }
   }
   
   for (int i=0; i<num_edge_vars(); ++i) {
      LevelData<EdgeDataBox>& this_edge_data = edge_var(i);
      const DisjointBoxLayout& dbl( this_edge_data.disjointBoxLayout() );
      for (DataIterator dit( this_edge_data.dataIterator() ); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; ++dir) {
            const Box& box( a_count_ghosts ? this_edge_data[dit][dir].box() : enclosedCells(surroundingNodes(dbl[dit]),dir) );
            size += box.numPts() * this_edge_data.nComp();
         }
      }
   }

   for (int i=0; i<num_node_vars(); ++i) {
      LevelData<NodeFArrayBox>& this_node_data = node_var(i);
      const DisjointBoxLayout& dbl( this_node_data.disjointBoxLayout() );
      for (DataIterator dit( this_node_data.dataIterator() ); dit.ok(); ++dit) {
         const Box& box( a_count_ghosts ? this_node_data[dit].getFab().box() : surroundingNodes(dbl[dit]) );
         size += box.numPts() * this_node_data.nComp();
      }
   }

   return size;
}


bool CFGVars::conformsTo( const CFGVars&  a_rhs,
                          const bool      a_include_ghost_cells ) const
{
   bool status( true );

   for (int i=0; i<num_cell_vars(); ++i) {
      const LevelData<FArrayBox>& this_cell_data = cell_var(i);
      const LevelData<FArrayBox>& that_cell_data = a_rhs.cell_var(i);

      const DisjointBoxLayout& thisBoxes( this_cell_data.disjointBoxLayout() );
      const DisjointBoxLayout& rhsBoxes( that_cell_data.disjointBoxLayout() );
   
      status &= thisBoxes.compatible( rhsBoxes );
      status &= ( this_cell_data.nComp() == that_cell_data.nComp() );
   
      if ( a_include_ghost_cells ) {
         status &= ( this_cell_data.ghostVect() == that_cell_data.ghostVect() );
      }

   }

   for (int i=0; i<num_face_vars(); ++i) {
      const LevelData<FluxBox>& this_face_data = face_var(i);
      const LevelData<FluxBox>& that_face_data = a_rhs.face_var(i);

      const DisjointBoxLayout& thisBoxes( this_face_data.disjointBoxLayout() );
      const DisjointBoxLayout& rhsBoxes( that_face_data.disjointBoxLayout() );
   
      status &= thisBoxes.compatible( rhsBoxes );
      status &= ( this_face_data.nComp() == that_face_data.nComp() );
   
      if ( a_include_ghost_cells ) {
         status &= ( this_face_data.ghostVect() == that_face_data.ghostVect() );
      }
   }
   
   for (int i=0; i<num_edge_vars(); ++i) {
      const LevelData<EdgeDataBox>& this_edge_data = edge_var(i);
      const LevelData<EdgeDataBox>& that_edge_data = a_rhs.edge_var(i);

      const DisjointBoxLayout& thisBoxes( this_edge_data.disjointBoxLayout() );
      const DisjointBoxLayout& rhsBoxes( that_edge_data.disjointBoxLayout() );
   
      status &= thisBoxes.compatible( rhsBoxes );
      status &= ( this_edge_data.nComp() == that_edge_data.nComp() );
   
      if ( a_include_ghost_cells ) {
         status &= ( this_edge_data.ghostVect() == that_edge_data.ghostVect() );
      }
   }

   for (int i=0; i<num_node_vars(); ++i) {
      const LevelData<NodeFArrayBox>& this_node_data = node_var(i);
      const LevelData<NodeFArrayBox>& that_node_data = a_rhs.node_var(i);

      const DisjointBoxLayout& thisBoxes( this_node_data.disjointBoxLayout() );
      const DisjointBoxLayout& rhsBoxes( that_node_data.disjointBoxLayout() );

      status &= thisBoxes.compatible( rhsBoxes );
      status &= ( this_node_data.nComp() == that_node_data.nComp() );

      if ( a_include_ghost_cells ) {
         status &= ( this_node_data.ghostVect() == that_node_data.ghostVect() );
      }
   }

   return status;
}


void
CFGVars::interpFaceVarToCell( LevelData<FArrayBox>&  a_cell_var,
                        const string&                a_face_var_name ) const
{
   CH_assert( face_var(a_face_var_name).isDefined() );
   SpaceUtils::interpFaceVectorToCell(a_cell_var,face_var(a_face_var_name),"c2");
}


void
CFGVars::interpEdgeVarToCell( LevelData<FArrayBox>&  a_cell_var,
                        const string&                a_edge_var_name ) const
{
   CH_assert( edge_var(a_edge_var_name).isDefined() );
   SpaceUtils::interpEdgesToCell(a_cell_var,edge_var(a_edge_var_name),"c2");
}

void
CFGVars::interpNodeVarToCell( LevelData<FArrayBox>&  a_cell_var,
                        const string&                a_node_var_name ) const
{
   CH_assert( node_var(a_node_var_name).isDefined() );
   SpaceUtils::interpNodesToCells(a_cell_var,node_var(a_node_var_name),"c2");
}

#include "NamespaceFooter.H"

