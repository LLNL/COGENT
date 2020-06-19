#include "LogicalSheathBC.H"
#include "Directions.H"
#include "FlipGrids.H"
#include "LogicalSheathBCF_F.H"

#undef CH_SPACEDIM
#include "ReductionOps.H.multidim"
#include "ReductionCopier.H.multidim"
#include "SpreadingCopier.H.multidim"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif

#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

const string LogicalSheathBC::pp_name = "logical_sheath";

LogicalSheathBC::LogicalSheathBC( const BoundaryBoxLayoutPtr&  a_bdry_layout,
                                  const ParmParse&             a_pp)

: m_bdry_layout(a_bdry_layout),
  m_compute_potential_BC(false),
  m_sheath_bc_type("first_order")
{
   // Create boundary dbl, which is one-cell-wide in vpar and mu
   const DisjointBoxLayout& grids_full( a_bdry_layout->disjointBoxLayout() );
   DisjointBoxLayout grids_tmp;
   adjCellLo(grids_tmp, grids_full, VPARALLEL_DIR, -1);
   adjCellLo(m_grids_inj, grids_tmp, MU_DIR, -1);
   
   parseParameters( a_pp );
}

void
LogicalSheathBC::computeSheathBC(LevelData<FArrayBox>& a_phi_bc,
                                 const KineticSpecies& a_species) const
{
   // Get coordinate system parameters
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   const MultiBlockCoordSys& coord_sys( *(geometry.coordSysPtr()) );
   const PhaseGrid& phase_grid = geometry.phaseGrid();
   
   // Get dfn and create bdry_data container
   const DisjointBoxLayout& bdry_dbl( m_bdry_layout->disjointBoxLayout() );
   LevelData<FArrayBox> bdry_data(bdry_dbl, 1, IntVect::Zero);
   const LevelData<FArrayBox>& dfn = a_species.distributionFunction();

   // Solve for potential BC
   Real residual = DBL_MAX;
   while (residual > 1.0e-6) {

      // Fill boundary data object with dfn
      fillBoundaryData(bdry_data, dfn);

      for (DataIterator dit(bdry_dbl); dit.ok(); ++dit) {
         const Box& bdry_box = bdry_dbl[dit];
         const FArrayBox& this_data = bdry_data[dit];
         FArrayBox& this_phi = a_phi_bc[dit];
         /*
            Here, do something to the phase-space bdry_data using
            the current iteration of a_phi_bc. Note that this_phi is defined
            with m_grids_inj, which are flattened in vpar and mu.
         */
      }
      
      // Integrate over velocity space
      LevelData<FArrayBox> bdry_data_summed(m_grids_inj, 1, IntVect::Zero);
      sum(bdry_data, bdry_data_summed);

      Real local_residual = DBL_MAX;
      for (DataIterator dit(m_grids_inj); dit.ok(); ++dit) {
         const FArrayBox& this_data = bdry_data_summed[dit];
         FArrayBox& this_phi = a_phi_bc[dit];
         FArrayBox this_residual(m_grids_inj[dit], 1);
         /*
            Here, compute residual and update a_phi_bc
         */
         local_residual = this_residual.norm();
      }
      
#ifdef CH_MPI
      MPI_Allreduce(&local_residual, &residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
      residual = local_residual;
#endif
   }
}

void
LogicalSheathBC::fillBoundaryData(LevelData<FArrayBox>& a_bdry_data,
                                  const LevelData<FArrayBox>& a_dfn) const
{
   /*
    Copy distribution function from the last valid cell
    into the bdry_data object that occupies one-cell-wide boundary layer.
    */
   const DisjointBoxLayout& bdry_grids( a_bdry_data.disjointBoxLayout() );
   for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

      const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
      const FArrayBox& this_dfn( a_dfn[interior_dit] );
      FArrayBox& this_data( a_bdry_data[dit] );
      
      FArrayBox tmp(this_dfn.box(), 1);
      tmp.copy(this_dfn);
      const int& dir( m_bdry_layout->dir() );
      const Side::LoHiSide& side( m_bdry_layout->side() );
      tmp.shift(dir, sign(side));
      this_data.copy(tmp);
   }
}

void
LogicalSheathBC::sum( const LevelData<FArrayBox>& a_src,
                      LevelData<FArrayBox>&       a_dst ) const
{
   /*
    Sum over the velocity space
    */
   
   CH_assert(a_src.nComp() == a_dst.nComp());

   const DisjointBoxLayout& dst_grids = a_dst.getBoxes();
   const DisjointBoxLayout& src_grids = a_src.getBoxes();
   const ProblemDomain& problem_domain = src_grids.physDomain();

   // Initialize the destination, since SumOp does not do that for us.
   DataIterator dit = a_dst.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_dst[dit].setVal(0.);
   }
   
   DisjointBoxLayout grids_tmp;
   adjCellLo(grids_tmp, src_grids, VPARALLEL_DIR, -1);

   // Initialize the destination, since SumOp does not do that for us.
   LevelData<FArrayBox> tmp(grids_tmp, a_src.nComp(), IntVect::Zero);
   for (DataIterator dit(grids_tmp); dit.ok(); ++dit) {
      tmp[dit].setVal(0.);
   }
   
   // Define ReductionCopier to compute intersections (sum in the poloidal direction)
   ReductionCopier reduceCopierVp(src_grids, grids_tmp, problem_domain, VPARALLEL_DIR);
   
   SumOp opVp(VPARALLEL_DIR);
   opVp.scale = 1.0;

   // Do the summing operation -- sums data in src along the vparallel direction
   a_src.copyTo(a_src.interval(), tmp, tmp.interval(), reduceCopierVp, opVp);

   // Define ReductionCopier to compute intersections (sum in the mu  direction)
   ReductionCopier reduceCopierMu(grids_tmp, dst_grids, problem_domain, MU_DIR);
   
   SumOp opMu(MU_DIR);
   opMu.scale = 1.0;

   // Do the summing operation -- sums data in src along the nu direction
   tmp.copyTo(tmp.interval(), a_dst, a_dst.interval(), reduceCopierMu, opMu);
}

void LogicalSheathBC::applyBC( KineticSpecies& a_species,
                               const LevelData<FluxBox>& a_velocity,
                               const CFG::LevelData<CFG::FArrayBox>& a_phi) const
{
   /*
    This BC reflects particles that cannot make the potential barrier
    */
   
   const int& dir( m_bdry_layout->dir() );
   const Side::LoHiSide& side( m_bdry_layout->side() );
   
   LevelData<FArrayBox>& soln( a_species.distributionFunction() );
   const double mass = a_species.mass();
   const double charge = a_species.charge();
            
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
         
   LevelData<FArrayBox> phi_injected;
   if (!m_compute_potential_BC) {
      geometry.injectConfigurationToPhase(a_phi, phi_injected);
   }
   else {
      phi_injected.define(m_grids_inj, 1, IntVect::Zero);
      computeSheathBC(phi_injected, a_species);
   }
         
   // Create valid-cell box that extends ghost_vect number
   // of cells away from the boundary
   const Box& domain_box = (geometry.domain()).domainBox();
   const IntVect& ghost_vect( soln.ghostVect() );
   const DisjointBoxLayout& dbl = soln.getBoxes();
         
   IntVect lo_end(domain_box.smallEnd());
   IntVect hi_end(domain_box.bigEnd());

   if (side == Side::Lo) hi_end[dir] = lo_end[dir] + ghost_vect[dir];
   if (side == Side::Hi) lo_end[dir] = hi_end[dir] - ghost_vect[dir];
         
   Box refl_bnd_box(lo_end, hi_end);
         
   // Create the flipped data object we will need for reflecting particles (nominally, electrons)
   //  below the potential barrier
   int reflectDir = VPARALLEL_DIR;
   int reflectCoord = 0;
   Vector<Tuple<DataIndex, 2> > boxCorrelation;
   
   DisjointBoxLayout flippedGrids;
   getFlippedGrids(flippedGrids, boxCorrelation, dbl, refl_bnd_box, reflectDir, reflectCoord);
   LevelData<FArrayBox> flippedData(flippedGrids, 1);
   soln.copyTo(flippedData);

   // Iterate over patches of flipped data, a small subset of the full data
   DataIterator fdit = flippedData.dataIterator();
   for (fdit.begin(); fdit.ok(); ++fdit) {
            
      // find the iterator value for the UNFLIPPED data corresponding to this flipped data
      DataIndex regDataIndex;
      for (int n=0; n<boxCorrelation.size(); n++)
      {
         if (boxCorrelation[n][1] == fdit() )
         {
            regDataIndex = boxCorrelation[n][0];
         }
      }
            
      FArrayBox& this_dfn = soln[regDataIndex];
      const Box& this_box = dbl[regDataIndex];
      FArrayBox& this_phi = phi_injected[regDataIndex];
      
      // Because m_grids_inj corresponds to ghost cell layer,
      // need to shift phi_bc to the last valid cell
      if (m_compute_potential_BC) {
         this_phi.shift(dir, -sign(side));
      }

      const FluxBox& this_vel = a_velocity[regDataIndex];
                  
      Box boundaryBox = adjCellBox(this_box, dir, side, ghost_vect[dir]);
            
      FArrayBox velocityRealCoords(boundaryBox, VEL_DIM);
      const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(this_box);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);
            
      const int SIDE(side);
      
      if (m_sheath_bc_type == "first_order") {
         FORT_SET_LOGICAL_SHEATH_BC(CHF_FRA1(this_dfn,0),
                                    CHF_BOX(boundaryBox),
                                    CHF_CONST_FRA1(flippedData[fdit],0),
                                    CHF_CONST_FRA(velocityRealCoords),
                                    CHF_CONST_FRA1(this_vel[dir],dir),
                                    CHF_CONST_FRA1(this_phi,0),
                                    CHF_CONST_REAL(mass),
                                    CHF_CONST_REAL(charge),
                                    CHF_CONST_INT(SIDE) );
      }
      
      else if (m_sheath_bc_type == "second_order") {
      
         FORT_SET_LOGICAL_SHEATH_BC_ORDER2(CHF_FRA1(this_dfn,0),
                                           CHF_BOX(boundaryBox),
                                           CHF_CONST_FRA1(flippedData[fdit],0),
                                           CHF_CONST_FRA(velocityRealCoords),
                                           CHF_CONST_FRA1(this_vel[dir],dir),
                                           CHF_CONST_FRA1(this_phi,0),
                                           CHF_CONST_REAL(mass),
                                           CHF_CONST_REAL(charge),
                                           CHF_CONST_INT(SIDE) );
      }
      
      else {
         MayDay::Error( "LogicalSheathBC::applyBC: unknown BC type" );
      }
      
   }
}

void
LogicalSheathBC::parseParameters(const ParmParse& a_pp)
{
   a_pp.query("compute_potential_BC", m_compute_potential_BC);
   a_pp.query("sheath_bc_type", m_sheath_bc_type);
}

#if 0
LogicalSheathBC::computeSheathBC(LevelData<FArrayBox> a_phi_bc,
                                 const KineticSpecies& a_species)
{
   /*
      Experimental option in progress. This one uses low-level MPI stuff.
      May want to continue development if the working option gets expensive.
    */
   
   // Get coordinate system parameters
   const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   const MultiBlockCoordSys& coord_sys( *(geometry.coordSysPtr()) );
   const PhaseGrid& phase_grid = geometry.phaseGrid();
   
   // Get decomposition data
   int num_config_boxes = phase_grid.numConfigBoxes();
   
   // Fill boundary data object with dfn
   const DisjointBoxLayout& bdry_dbl( m_bdry_layout->disjointBoxLayout() );
   LevelData<FArrayBox> bdry_data(bdry_dbl, 1, IntVect::Zero);
   const LevelData<FArrayBox>& dfn = a_species.distributionFunction();
   fillBoundaryData(bdry_data, dfn);

   // Loop over patches of boundary dbl and search for a valid config_box
   // neighbouring the 4D boundary patch, then get info about procIDs (MPI_comm)
   // that contain velocity slices for this cfg patch.
   for (DataIterator dit(bdry_dbl); dit.ok(); ++dit) {
      for (int k=0; k<num_config_boxes; ++k) {
         const CFG::Box& config_box = phase_grid.configBox(k);
         const Box& bdry_box = bdry_dbl[dit];
         const Side::LoHiSide& side( bdry_layout.side() );
         if ((side == Side::Lo && config_box.smallEnd(SHEATH_DIR) == bdry_box.bigEnd(SHEATH_DIR) + 1) ||
             (side == Side::Hi && config_box.bigEnd(SHEATH_DIR) == bdry_box.smallEnd(SHEATH_DIR) - 1)) {
            
            const FArrayBox& this_data = bdry_data[dit];
            FArrayBox& this_phi = a_phi_bc[dit];
            
            FArrayBox residual(m_grids_inj[dit], 1);
            
            residual.setVal(1.0);
            while (residual.norm() > 1.0e-6) {
            
               FArrayBox tmp(this_data);
               /*
                  Here, we do something to the phase-space bdry_data using
                  current iteration of phi.
               */
            
               FArrayBox reduced_data(bdry_dbl[dit], 1);
            
               int size = bdry_dbl[dit].numPts();
               const MPI_Comm& config_box_comm = phase_grid.configBoxComm(k);
               MPI_Allreduce(this_data.dataPtr(), reduced_data.dataPtr(), size, MPI_DOUBLE, MPI_SUM, config_box_comm);
            
               Box
               
               /*
                reduced_data contains sum over different velocity slices; in order to get vpar and mu inegrals we
                now need to sum over vpar and mu indices. Then we should check if the obtained function satisfies
                our criteria...
                */
               
               
            }
         }
      }
   }
   
}
#endif

#include "NamespaceFooter.H"

