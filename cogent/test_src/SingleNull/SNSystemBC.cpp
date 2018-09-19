#include "SNSystemBC.H"

#include "inspect.H"

#include "BlockBoundary.H"
#include "Box.H"
#include "DataIterator.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "IntVect.H"
#include "LoHiSide.H"
#include "MultiBlockCoordSys.H"
#include "MagCoordSys.H"
#include "MagBlockCoordSys.H"
#include "Tuple.H"
#include "Vector.H"

#include "SingleNullBlockCoordSys.H"

#include "SingleNullBC.H"

#include "NamespaceHeader.H"


inline
std::string determineCoordSysType( const MagGeom& a_geometry )
{
   const MagBlockCoordSys& block_coord_sys( a_geometry.getBlockCoordSys(1) );
   return block_coord_sys.geometryType();
}


SNSystemBC::SNSystemBC( ParmParse& a_pp,
                        const MagGeom& a_geometry,
                        const LevelData<FArrayBox>& a_solution )
   : m_verbosity(0),
     m_geometry( a_geometry ),
     m_bcs(NULL),
     m_mblexPtr(NULL)

{
   a_pp.query( "gksystem.verbosity", m_verbosity );

   const std::string coord_sys_type( determineCoordSysType( a_geometry ) );

   parseInput( a_pp, coord_sys_type, a_solution );

   m_mblexPtr = new MultiBlockLevelExchangeAverage();

   int numGhost(4);
   int spaceOrder(4);
#ifdef TIME_MULTIBLOCK
   MPI_Barrier( MPI_COMM_WORLD );
   double define_time( MPI_Wtime() / 60.0 );
#endif
   if (procID()==0) cout << "Beginning exchange define" << endl;
   m_mblexPtr->define( &m_geometry, numGhost, spaceOrder );
#ifdef TIME_MULTIBLOCK
   define_time = MPI_Wtime() / 60.0 - define_time;
   cout << "Exchange define time on proc " << procID()
        << " = " << define_time << endl;
#endif
 }


SNSystemBC::~SNSystemBC()
{
   delete m_bcs;
   delete m_mblexPtr;
}


inline
void SNSystemBC::executeInternalExchanges( LevelData<FArrayBox>& a_solution ) const
{
  // Fill ghost cells except for those on block boundaries
   // (aka extra-block ghost cells)
   a_solution.exchange();

   // Fill extra-block ghost cells
#ifdef TIME_MULTIBLOCK
   double interpGhosts_time( MPI_Wtime() / 60.0 );
#endif

   m_mblexPtr->interpGhosts( a_solution );

#ifdef TIME_MULTIBLOCK
   interpGhosts_time = MPI_Wtime() / 60.0 - interpGhosts_time;
   cout << "interpGhosts time on proc " << procID() << " = "
        << interpGhosts_time << endl;
#endif
}

#if 0
void SNSystemBC::fillBlockBoundaryCellsInCompSpace(
   LevelData<FArrayBox>& a_solution_comp,
   const LevelData<FArrayBox>& a_solution_phys ) const
{
   const MultiBlockCoordSys& coord_sys( *(m_geometry.coordSysPtr()) );
   typedef Tuple<BlockBoundary,2*SpaceDim> BlockBoundaryTuple;
   const Vector<BlockBoundaryTuple>& boundaries( coord_sys.boundaries() );

   LevelData<FArrayBox>& uJ( a_solution_comp );
   const DisjointBoxLayout& grid( uJ.disjointBoxLayout() );
   const LevelData<FArrayBox>& u( a_solution_phys );
   for (DataIterator dit( grid ); dit.ok(); ++dit) {

      const int i_block( coord_sys.whichBlock( grid[dit] ) );
      const BlockBoundaryTuple block_bdry( boundaries[i_block] );

      for (int dir(0); dir<CH_SPACEDIM; dir++) {
         for (SideIterator sit; sit.ok(); ++sit) {
            int i( dir + CH_SPACEDIM * sit() );
            if (block_bdry[i].isInterface()) {
               const IntVect& ghost_vect( uJ.ghostVect() );
               Box box( grid[dit] );
               box = adjCellBox( box, dir, sit(), ghost_vect[dir] );
               for (int ddir(0); ddir<CH_SPACEDIM; ddir++) {
                  if (ddir!=dir) {
                     for (SideIterator ssit; ssit.ok(); ++ssit) {
                        int j( ddir + CH_SPACEDIM * ssit() );
                        if (!block_bdry[j].isDomainBoundary()) {
                           box.growDir( ddir, ssit(), ghost_vect[ddir] );
                        }
                     }
                  }
               }
               Box domain_box( box );
               domain_box.growDir( dir, Side::flip( sit() ), 1 );
               ProblemDomain domain( domain_box );
               const FArrayBox& thisU( u[dit] );
               FArrayBox& thisUJ( uJ[dit] );
               const MagBlockCoordSys&
               block_coord_sys( m_geometry.getBlockCoordSys( grid[dit] ) );
//               inspect( thisU );
               block_coord_sys.multJ( thisU, box, domain, thisUJ );
//               inspect( thisUJ );
            }
         }
      }
   }
//   u.copyTo( uJ );
//   m_geometry.multJ( thisU, box, domain, thisUJ );
   uJ.exchange();
}
#endif


#if 0
void SNSystemBC::fillDomainGhostCellsInPhysSpace(
   LevelData<FArrayBox>& a_solution_phys,
   const LevelData<FArrayBox>& a_solution_comp ) const
{
   const MultiBlockCoordSys& coord_sys( *(m_geometry.coordSysPtr()) );
   typedef Tuple<BlockBoundary,2*SpaceDim> BlockBoundaryTuple;
   const Vector<BlockBoundaryTuple>& boundaries( coord_sys.boundaries() );

   LevelData<FArrayBox>& u( a_solution_phys );
   const DisjointBoxLayout& grid( u.disjointBoxLayout() );
   const LevelData<FArrayBox>& uJ( a_solution_comp );

   for (DataIterator dit( grid ); dit.ok(); ++dit) {

      const int i_block( coord_sys.whichBlock( grid[dit] ) );
      const BlockBoundaryTuple block_bdry( boundaries[i_block] );

      for (int dir(0); dir<CH_SPACEDIM; dir++) {
         for (SideIterator sit; sit.ok(); ++sit) {
            int i( dir + CH_SPACEDIM * sit() );
            if (block_bdry[i].isDomainBoundary()) {
               Box box( grid[dit] );
               const IntVect& ghost_vect( u.ghostVect() );
               box.grow( ghost_vect );
               box.growDir( dir, sit(), -ghost_vect[dir] );
               box = adjCellBox( box, dir, sit(), ghost_vect[dir] );
               FArrayBox& thisU( u[dit] );
               const FArrayBox& thisUJ( uJ[dit] );
               const MagBlockCoordSys&
               block_coord_sys( m_geometry.getBlockCoordSys( grid[dit] ) );
//               inspect( thisU );
               block_coord_sys.divideJ( thisUJ, box, thisU );
//               inspect( thisUJ );
            }
         }
      }
   }
   u.exchange();
}
#endif

void SNSystemBC::fillGhostCells(
   LevelData<FArrayBox>& a_u,
   const LevelData<FArrayBox>& a_uJ,
   const LevelData<FluxBox>& a_velocity,
   const Real& a_time ) const
{
#if 0
   LevelData<FArrayBox> solution_comp;
   solution_comp.define( a_solution_comp );
   a_solution_comp.copyTo( solution_comp );
   
   LevelData<FArrayBox>& solution_phys( a_solution_phys );
   executeInternalExchanges( solution_phys );
   
   fillBlockBoundaryCellsInCompSpace( solution_comp, solution_phys );
   
   m_bcs->apply( solution_comp, m_geometry, a_velocity, a_time );
   
   fillDomainGhostCellsInPhysSpace( solution_phys, solution_comp );
#else
   executeInternalExchanges( a_u );
   m_bcs->apply( a_u, m_geometry, a_velocity, a_time );
#endif   
}


void SNSystemBC::parseInput( ParmParse& a_pp,
                             const std::string& a_coord_sys_type,
                             const LevelData<FArrayBox>& a_solution )
{
   ParmParse ppsp( "BC" );
   m_bcs = new SingleNullBC( "solution", ppsp, m_verbosity );
}


#include "NamespaceFooter.H"
