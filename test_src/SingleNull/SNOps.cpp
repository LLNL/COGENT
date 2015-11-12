#include "SNOps.H"
#include "CONSTANTS.H"
#include "FORT_PROTO.H"
#include "inspect.H"

#include <fstream>
#include <sstream>

#include "LocalizedF_F.H"

#if 1  // warning, OS dependencies, will not work on all platforms
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "CH_HDF5.H"
#include "HDF5Portable.H"

#include "BoxIterator.H"
#include "LoadBalance.H"
#include "MagBlockCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "newMappedGridIO.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "inspect.H"

#include "AMRIO.H"

#include "EdgeToCell.H"
#include "FourthOrderUtil.H"
//#include "Directions.H"

#include "NamespaceHeader.H"

SNOps::SNOps( ParmParse& a_pp,
              const MagGeom& a_geometry,
              const SNSystemBC& a_boundary_conditions,
              const SNSystemIC& a_initial_conditions,
              const int a_verbosity  )
   : m_operator(NULL),
     m_initial_conditions( a_initial_conditions ),
     m_boundary_conditions( a_boundary_conditions ),
     m_geometry( a_geometry ),
     m_ghost_vect( 4*IntVect::Unit ),
     m_verbosity( a_verbosity )
{
   parseParameters( a_pp );
   m_operator = new Advect( ParmParse( Advect::pp_name ) );
}


SNOps::~SNOps()
{
   delete m_operator;
}


Real SNOps::stableDt( const LevelData<FArrayBox>& a_soln,
                      const int a_step_number )
{
   return m_operator->computeDt( m_geometry, a_soln );
}


void SNOps::evalRHS( LevelData<FArrayBox>&       a_rhs,
                     const LevelData<FArrayBox>& a_soln_comp,
                     const int                   a_step_number,
                     const Real                  a_time,
                     const int                   a_stage )
{
   LevelData<FArrayBox> a_soln_phys( a_soln_comp.disjointBoxLayout(),
                                     a_soln_comp.nComp(),
                                     a_soln_comp.ghostVect() );
   a_soln_comp.copyTo( a_soln_phys );
   m_geometry.divideJonValid( a_soln_phys );

   LevelData<FluxBox> velocity( a_soln_comp.disjointBoxLayout(),
                                SpaceDim,
                                IntVect::Unit );
   m_operator->updateVelocity( m_geometry, velocity, true );
   m_boundary_conditions.fillGhostCells( a_soln_phys,
                                         a_soln_comp,
                                         velocity,
                                         a_time );

   m_operator->evalRHS( m_geometry, a_rhs, a_soln_phys, a_time );
}


void SNOps::divideJ( const LevelData<FArrayBox>& a_soln_mapped,
                     LevelData<FArrayBox>&       a_soln_physical,
                     bool                        a_restrict_to_valid )
{
   CH_assert( a_soln_mapped.ghostVect() >= a_soln_physical.ghostVect() );

   const DisjointBoxLayout& grids = a_soln_mapped.disjointBoxLayout();

   IntVect grown_vect = IntVect::Unit;

   LevelData<FArrayBox> J(grids, 1, grown_vect);
   m_geometry.getJ( J );

   LevelData<FArrayBox> grown_soln(grids, 1, grown_vect);

   DataIterator dit = a_soln_physical.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     grown_soln[dit].copy(a_soln_mapped[dit]);
   }
   grown_soln.exchange();

   // We now have soln and J with one layer of ghost cells that have been filled at least
   // on interior box sides.  Compute fourth-order quotient using one-sided differences at
   // physical boundaries.

   const MagCoordSys* coords( m_geometry.getCoordSys() );
   for (dit.begin(); dit.ok(); ++dit) {
      int block_number = coords->whichBlock( grids[dit] );
      const SingleNullBlockCoordSys* coord_sys
         = dynamic_cast<const SingleNullBlockCoordSys*>( coords->getCoordSys(block_number) );

      cellFGToCellF(a_soln_physical[dit], grown_soln[dit], J[dit], grids[dit], coord_sys->domain(), true);
   }

   a_soln_physical.exchange();
}


void SNOps::parseParameters( ParmParse& a_ppgksys )
{
}


void SNOps::writeCheckpointFile( HDF5Handle& a_handle ) const
{
}
   
void SNOps::readCheckpointFile( HDF5Handle& a_handle, const int& a_cur_step )
{
}


#include "NamespaceFooter.H"

