#include "LogRectCoordSys.H"
#include "MagBlockCoordSys.H"
#include "MillerBlockCoordSys.H"
#include "SlabBlockCoordSys.H"
#include "ToroidalBlockCoordSys.H"
#include "CylindricalBlockCoordSys.H"
#include "FourthOrderUtil.H"
#include "BlockBoundary.H"
#include "SpreadingCopier.H"
#include "ReductionOps.H"
#include "Directions.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

const std::string LogRectCoordSys::pp_name = "LogicallyRectangular";

LogRectCoordSys::LogRectCoordSys(ParmParse&               a_pp,
                                 const string&            a_mag_geom_type,
                                 const std::vector<int>&  a_numcells,
                                 const std::vector<bool>& a_is_periodic,
                                 const std::vector<int>&  a_decomposition)
  : m_decomposition(a_decomposition),
    m_num_blocks(1)
{
  m_spread_radially = false;
  m_mag_geom_type = a_mag_geom_type;

  //mb_dir multiblock direction (y_dir for 3D, z_dir for 2D)
#if CFG_DIM == 3
    m_mb_dir = TOROIDAL_DIR;
#else
   //this options has not been fully developed yet
    m_mb_dir = POLOIDAL_DIR;
#endif

  a_pp.query( "number_of_blocks", m_num_blocks);
  m_decomposition[m_mb_dir] /= m_num_blocks;

  if (m_num_blocks != 1 && m_decomposition[m_mb_dir] != 1) {
    MayDay::Error("LogRectCoordSys:: current multiblock sheared geometry implementation only works for the same number of blocks and proc in the toroidal dir");  
  }
  
  int num_block_cells = a_numcells[m_mb_dir] / m_num_blocks;
  
  Vector<Box> domain_boxes(m_num_blocks);
  for ( int block_number=0; block_number<m_num_blocks; ++block_number ) {
    IntVect lo = IntVect::Zero;
    IntVect hi;
    for (int dir=0; dir<SpaceDim; ++dir) {
      hi[dir] = a_numcells[dir] - 1;
      if (dir == m_mb_dir ) {
        lo[dir] = (BLOCK_SEPARATION + num_block_cells) * block_number;
        hi[dir] = lo[dir] + num_block_cells - 1;
      }
    }
    domain_boxes[block_number] = Box(lo, hi);
  }

  bool is_periodic[SpaceDim];
  for (int dir=0; dir<SpaceDim; ++dir) {
    is_periodic[dir] = a_is_periodic[dir];
  }

  for ( int block_number = 0; block_number < m_num_blocks; ++block_number ) {

     if (m_mag_geom_type == "miller") {
        MillerBlockCoordSys* block_coords = new MillerBlockCoordSys( a_pp, ProblemDomain(domain_boxes[block_number], is_periodic));
        m_coord_vec.push_back(block_coords);
        m_spread_radially = block_coords->getConstMinorrad() != 0;
     }
     else if (m_mag_geom_type == "slab") {
        SlabBlockCoordSys* block_coords = new SlabBlockCoordSys( a_pp, ProblemDomain(domain_boxes[block_number], is_periodic),
                                                                 block_number, m_num_blocks, BLOCK_SEPARATION);
        m_coord_vec.push_back(block_coords);
     }
#if CFG_DIM==3
     else if (m_mag_geom_type == "toroidal") {
       ToroidalBlockCoordSys* block_coords = new ToroidalBlockCoordSys( a_pp, ProblemDomain(domain_boxes[block_number], is_periodic),
                                                                    block_number, m_num_blocks, BLOCK_SEPARATION);
       m_coord_vec.push_back(block_coords);
     }
#endif
     else if (m_mag_geom_type == "cylindrical") {
        CylindricalBlockCoordSys* block_coords = new CylindricalBlockCoordSys( a_pp, ProblemDomain(domain_boxes[block_number], is_periodic));
        m_coord_vec.push_back(block_coords);
     }
  }

  defineCoordSystemsAndBoundaries(m_coord_vec);

#if 1
  // Define the boundary conditions for divergence cleaning (whether or not they're used)
  // Does not check for periodicity. FIX IF NEEDED.
  for (int block = 0; block < m_num_blocks; block++) {
    Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[block];
    const MagBlockCoordSys& coord_sys( *(MagBlockCoordSys*)getCoordSys(block) );
    const ProblemDomain& domain( coord_sys.domain() );

    for (int dir=0; dir<SpaceDim; ++dir) {
      for (int side=0; side<2; ++side) {
        if ( blockBoundaries[dir + side*SpaceDim].isDomainBoundary() && !domain.isPeriodic(dir)) {
          double bc_value = 0.;
          int bc_type = EllipticOpBC::DIRICHLET;     // Homogeneous Dirichlet
          m_divergence_cleaning_bc.setBCType(block, dir, side, bc_type);
          m_divergence_cleaning_bc.setBCValue(block, dir, side, bc_value);
        }
      }
    }
  }
#endif
}



LogRectCoordSys::~LogRectCoordSys()
{
   for (int i=0; i<m_coord_vec.size(); ++i) {
      delete m_coord_vec[i];
   }
}



void
LogRectCoordSys::defineBoundaries()
{
   CH_assert(gotMappingBlocks());
   m_boundaries.resize(numBlocks());

   for (int iblock(0); iblock<numBlocks(); iblock++) {

      const MagBlockCoordSys& coord_sys( *(MagBlockCoordSys*)getCoordSys(iblock) );
      const ProblemDomain& domain( coord_sys.domain() );
      Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries( m_boundaries[iblock] );
     
      for (int idir(RADIAL_DIR); idir<SpaceDim; idir++) {
         if ( domain.isPeriodic(idir) ) {
	   //We have to put something at the periodic domain boundary. 
	   //The periodicity is used by exchange(), which is informed by the 
	   //isPeriodic member of the ProblemDomain object. Here, we basically  
	   //indicate that the boundary is not the one where physical BCs are 
	   //applied. The value of the "shift" variable in the code below 
	   //does not matter. It is however would be  important for the 
	   //multiblock exchange operations requiring the info about the connectivities.
	   //NB: blockBoundaries[idir] = BlockBoundary::CONFORMAL did not work as
	   // a workaround. M.A.D. 03/20/18
	   IntVect shift( m_mappingBlocks[iblock].size(idir) * BASISV(idir) );
	   IndicesTransformation it;
	   it.defineFromTranslation( shift );
	   blockBoundaries[idir].define( it, iblock );
	   it.defineFromTranslation( -shift );
	   blockBoundaries[idir+SpaceDim].define( it, iblock );
         }
         else if (idir == m_mb_dir && numBlocks()>1) {
	   IntVect shift( m_mappingBlocks[iblock].size(idir) * BASISV(idir) );
	   IndicesTransformation it;
	   it.defineFromTranslation( shift );
	   blockBoundaries[idir].define( it, iblock );
	   it.defineFromTranslation( -shift );
	   blockBoundaries[idir+SpaceDim].define( it, iblock );
         }
	 else {
            blockBoundaries[idir].define(0);
            blockBoundaries[idir+SpaceDim].define(0);
	 }
      }

#if CGF_DIM==3
      // Overwrite toroidal boundary coupling
      // Needed to properly construct stencils for MBSolver 

      if ( m_num_blocks == 1 ) {  // Assuming periodic coupling 
	shift = m_mappingBlocks[iblock].size(TOROIDAL_DIR) * BASISV(TOROIDAL_DIR);
	it.defineFromTranslation( shift );
	blockBoundaries[TOROIDAL_DIR].define( it, iblock );
	it.defineFromTranslation( -shift );
	blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, iblock );
      }

      else {
	// Lower face coupling
	if ( iblock > 0 ) {  // Couple to the block on the lower boundary
	  shift = -BLOCK_SEPARATION * BASISV(TOROIDAL_DIR);
	  it.defineFromTranslation( shift );
	  blockBoundaries[TOROIDAL_DIR].define( it, iblock - 1);
	}

	else {  // First toroidal block: couple it to the last one
	  shift = (m_num_blocks * (BLOCK_SEPARATION + m_mappingBlocks[iblock].size(TOROIDAL_DIR))
		   - BLOCK_SEPARATION) * BASISV(TOROIDAL_DIR);
	  it.defineFromTranslation( shift );
	  blockBoundaries[TOROIDAL_DIR].define( it, iblock + m_num_blocks - 1);
	}
	
	// Upper face coupling
	if ( toroidal_sector < m_num_blocks-1 ) {  // Couple to the block on the upper toroidal boundary
	  shift = BLOCK_SEPARATION * BASISV(TOROIDAL_DIR);
	  it.defineFromTranslation( shift );
	  blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, iblock + 1);
	}

	else {  // Last toroidal block: couple it to the first one
	  shift = -(m_num_blocks * (BLOCK_SEPARATION + m_mappingBlocks[iblock].size(TOROIDAL_DIR))
		    - BLOCK_SEPARATION) * BASISV(TOROIDAL_DIR);
	  it.defineFromTranslation( shift );
	  blockBoundaries[TOROIDAL_DIR+SpaceDim].define( it, iblock - m_num_blocks + 1);
	}
      }
#endif
   }
}



void
LogRectCoordSys::blockRemapping(RealVect&       a_xi_valid,
                                int&            a_n_valid,
                                const RealVect& a_xiSrc,
                                int             a_nSrc) const
{
   //Need to have it here; otherwise have a failing assertion in MultiBlockLevelGeom::define 
   a_n_valid = a_nSrc;
   a_xi_valid = a_xiSrc;
}

void
LogRectCoordSys::postProcessMetricData( LevelData<FArrayBox>& a_data ) const
{
   if (m_spread_radially) {
      spreadRadially(a_data);
   }
}



void
LogRectCoordSys::postProcessMetricData( LevelData<FluxBox>& a_data ) const
{
   if (m_spread_radially) {
      spreadRadially(a_data);
   }
}



void
LogRectCoordSys::spreadRadially( LevelData<FArrayBox>& a_data ) const
{
   int radial_dir = 0;
   const MagBlockCoordSys& coord_sys( *(MagBlockCoordSys*)getCoordSys(0) );
   const ProblemDomain& block_zero_domain( coord_sys.domain() );

   Box domainBox = block_zero_domain.domainBox();
   int rmid = (domainBox.smallEnd(radial_dir) + domainBox.bigEnd(radial_dir) + 1)/2;
   
   LevelData<FArrayBox> grown_data(a_data.disjointBoxLayout(), a_data.nComp(), a_data.ghostVect()+2*IntVect::Unit);
   
   DataIterator dit = a_data.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      grown_data[dit].copy(a_data[dit]);
   }
   grown_data.exchange();
   
   LevelData<FluxBox> face_data(grown_data.disjointBoxLayout(), grown_data.nComp(), a_data.ghostVect());
   fourthOrderCellToFace(face_data, grown_data);
   
   LevelData<FArrayBox> data2(a_data.disjointBoxLayout(), a_data.nComp(), a_data.ghostVect());
   for (dit.begin(); dit.ok(); ++dit) {
      FArrayBox& radial_face = face_data[dit][radial_dir];
      radial_face.shiftHalf(radial_dir,1);
      data2[dit].copy(radial_face);
   }
   
   // Make a poloidal slice at the radial midpoint
   //LevelData<FArrayBox> temp;
   const SliceSpec radial_slice(radial_dir,rmid);
   //a_data.degenerate(temp, radial_slice);
   // do this explicitly here
   //   DisjointBoxLayout sliceGrids;
   Box sliceBox(domainBox);
   sliceBox.setSmall(radial_dir,rmid);
   sliceBox.setBig(radial_dir, rmid);
   Vector<Box> sliceBoxes(1,sliceBox);
   Vector<int> bogusProcAssign(1,0);
   DisjointBoxLayout sliceGrids(sliceBoxes, bogusProcAssign, block_zero_domain);
#if 0
   const DisjointBoxLayout& grids = a_data.getBoxes();
   grids.degenerate(sliceGrids, radial_slice);
#endif
   
   LevelData<FArrayBox> temp(sliceGrids, a_data.nComp(), a_data.ghostVect());
   data2.copyTo(temp);
   
   // Spread the slice radially
   Vector<int> spreadingDirs;
   spreadingDirs.push_back(radial_dir);
   
   SpreadingCopier spreadCopier(temp.getBoxes(),
                                a_data.getBoxes(),
                                block_zero_domain,
                                a_data.ghostVect(),
                                spreadingDirs);
   const SpreadingOp spreadOp(spreadingDirs);
   
   temp.copyTo(temp.interval(),
               a_data, a_data.interval(),
               spreadCopier,
               spreadOp);
}



void
LogRectCoordSys::spreadRadially( LevelData<FluxBox>& a_data ) const
{
   int radial_dir = 0;
   const MagBlockCoordSys& coord_sys( *(MagBlockCoordSys*)getCoordSys(0) );
   const ProblemDomain& block_zero_domain( coord_sys.domain() );

   Box domainBox = block_zero_domain.domainBox();
   int rmid = (domainBox.smallEnd(radial_dir) + domainBox.bigEnd(radial_dir) + 1)/2;
   
   LevelData<FluxBox> grown_data(a_data.disjointBoxLayout(), a_data.nComp(), a_data.ghostVect()+2*IntVect::Unit);
   
   DataIterator dit = a_data.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      grown_data[dit].copy(a_data[dit]);
   }
   grown_data.exchange();
   
   // Temporaries
   LevelData<FArrayBox> tmp_cell(grown_data.disjointBoxLayout(), grown_data.nComp(), grown_data.ghostVect());
   LevelData<FluxBox> tmp_face(a_data.disjointBoxLayout(), a_data.nComp(), a_data.ghostVect());
   
   // Average data on non-radial faces.  We do this by shifting to cell
   // centers, apply the 4th-order averaging, then copy the resulting
   // average on the middle flux surface
   for (int dir=0; dir<SpaceDim; ++dir) {
      if (dir != radial_dir) {
         for (int side=-1; side<2; side+=2) {
            for (dit.begin(); dit.ok(); ++dit) {
               FArrayBox& face(grown_data[dit][dir]);
               face.shiftHalf(dir,-side);
               tmp_cell[dit].copy(face);
               face.shiftHalf(dir,side);  // put it back
            }
            tmp_cell.exchange();
            fourthOrderCellToFace(tmp_face, tmp_cell);
            
            for (dit.begin(); dit.ok(); ++dit) {
               FArrayBox& face(tmp_face[dit][radial_dir]);
               face.shiftHalf(radial_dir,1);
               face.shiftHalf(dir,side);
               a_data[dit][dir].copy(face);
               face.shiftHalf(dir,-side); // put it back
               face.shiftHalf(radial_dir,-1);  // put it back
            }
         }
      }
   }
   a_data.exchange();
   
   // Make a poloidal slice at the radial midpoint
   //LevelData<FluxBox> temp;
   const SliceSpec radial_slice(radial_dir,rmid);
   //a_data.degenerae(temp, radial_slice);
   
   // do this explicitly here
#if 0
   DisjointBoxLayout sliceGrids;
   const DisjointBoxLayout& grids = a_data.getBoxes();
   grids.degenerate(sliceGrids, radial_slice);
#endif
   
   Box sliceBox(domainBox);
   sliceBox.setSmall(radial_dir,rmid);
   sliceBox.setBig(radial_dir, rmid);
   Vector<Box> sliceBoxes(1,sliceBox);
   Vector<int> bogusProcAssign(1,0);
   DisjointBoxLayout sliceGrids(sliceBoxes, bogusProcAssign, block_zero_domain);
   
   LevelData<FluxBox> temp(sliceGrids, a_data.nComp(), a_data.ghostVect());
   a_data.copyTo(temp);
   
   // Spread the slice radially
   Vector<int> spreadingDirs;
   spreadingDirs.push_back(radial_dir);
   
   SpreadingCopier spreadCopier(temp.getBoxes(),
                                a_data.getBoxes(),
                                block_zero_domain,
                                a_data.ghostVect(),
                                spreadingDirs);
   const FaceSpreadingOp spreadOp(spreadingDirs);
   
   temp.copyTo(temp.interval(),
               a_data, a_data.interval(),
               spreadCopier,
               spreadOp);
}

#if CFG_DIM == 3
void
LogRectCoordSys::toroidalBlockRemapping(IntVect& a_ivDst,
                                        Vector<Real>& a_interpStecil,
                                        const RealVect& a_xiSrc,
                                        const int a_nSrc,
                                        const Side::LoHiSide& a_side) const
{

   //Getting the destination block number (nDst)
   int nDst;
   if (a_side == Side::LoHiSide::Lo) nDst = a_nSrc - 1;
   if (a_side == Side::LoHiSide::Hi) nDst = a_nSrc + 1;
   
   //Correction for periodicity
   if (nDst == numBlocks()) {
     nDst = 0;
   }

   if (nDst == -1) {
     nDst = numBlocks() - 1;
   }
   
   const MagBlockCoordSys* src_coord_sys = (MagBlockCoordSys*)getCoordSys(a_nSrc);
   const RealVect& dx_src = src_coord_sys->getMappedCellSize();
   
   const MagBlockCoordSys* dst_coord_sys = (MagBlockCoordSys*)getCoordSys(nDst);
   const RealVect& dx_dst = dst_coord_sys->getMappedCellSize();
   
   CH_assert(dx_src==dx_dst);
   
   //Get mapped cell-centered coord in dst block
   RealVect x_src = src_coord_sys->realCoord(a_xiSrc);
   applyPeriodicity(x_src);
   RealVect xi0_dst = dst_coord_sys->mappedCoord(x_src);
   
   //Get mapped lower-left coord in dst block
   RealVect xiSW_src(a_xiSrc);
   xiSW_src[RADIAL_DIR] += -0.5 * dx_src[RADIAL_DIR];
   xiSW_src[POLOIDAL_DIR] += -0.5 * dx_src[POLOIDAL_DIR];

   RealVect xSW_src = src_coord_sys->realCoord(xiSW_src);
   applyPeriodicity(xSW_src);
   double xiSW_dst = (dst_coord_sys->mappedCoord(xSW_src))[POLOIDAL_DIR];
   
   //Get mapped upper-left coord in dst block
   RealVect xiNW_src(a_xiSrc);
   xiNW_src[RADIAL_DIR] += -0.5 * dx_src[RADIAL_DIR];
   xiNW_src[POLOIDAL_DIR] += 0.5 * dx_src[POLOIDAL_DIR];

   RealVect xNW_src = src_coord_sys->realCoord(xiNW_src);
   applyPeriodicity(xNW_src);
   double xiNW_dst = (dst_coord_sys->mappedCoord(xNW_src))[POLOIDAL_DIR];
   
   //Get mapped upper-right coord in dst block
   RealVect xiNE_src(a_xiSrc);
   xiNE_src[RADIAL_DIR] += 0.5 * dx_src[RADIAL_DIR];
   xiNE_src[POLOIDAL_DIR] += 0.5 * dx_src[POLOIDAL_DIR];

   RealVect xNE_src = src_coord_sys->realCoord(xiNE_src);
   applyPeriodicity(xNE_src);
   double xiNE_dst = (dst_coord_sys->mappedCoord(xNE_src))[POLOIDAL_DIR];
   
   //Get mapped lower-right coord in dst block
   RealVect xiSE_src(a_xiSrc);
   xiSE_src[RADIAL_DIR] += 0.5 * dx_src[RADIAL_DIR];
   xiSE_src[POLOIDAL_DIR] += -0.5 * dx_src[POLOIDAL_DIR];
   
   RealVect xSE_src = src_coord_sys->realCoord(xiSE_src);
   applyPeriodicity(xSE_src);
   double xiSE_dst = (dst_coord_sys->mappedCoord(xSE_src))[POLOIDAL_DIR];
  
  
   if (m_mag_geom_type == "toroidal") {
     enforcePoloidalCut(a_xiSrc, xiNW_dst, xiNE_dst, xiSW_dst, xiSE_dst, xi0_dst[POLOIDAL_DIR]);
   }
  
   //Check the remapping properties
   if (   (xiNE_dst - xiSE_dst) > 3.0 * dx_dst[POLOIDAL_DIR]
       || (xiNW_dst - xiSW_dst) > 3.0 * dx_dst[POLOIDAL_DIR]
       || (xiNE_dst - xiSE_dst) < 0.0
       || (xiNW_dst - xiSW_dst) < 0.0 )
   {
      MayDay::Error("LogRectCoordSys::toroidalBlockRemapping(): remapping failure; try decreasing radial and/or toroidal cell cize");
   }
   
   //Get the global index of the dst cell center
   IntVect iv0_dst;
   for (int dir=0; dir<SpaceDim; ++dir) {
      iv0_dst[dir] = floor(xi0_dst[dir]/dx_dst[dir]);
   }

   for (int dir=0; dir<SpaceDim; ++dir) {
      a_ivDst[dir] = iv0_dst[dir];
   }
   
   //Get interpolation coefficients
   int order = 3;
   getInterpolationCoefficients(a_interpStecil, xi0_dst, iv0_dst, dx_dst, order);

}

void
LogRectCoordSys::applyPeriodicity(RealVect& a_x) const                              
{

  const MagBlockCoordSys* src_coord_0 = (MagBlockCoordSys*)getCoordSys(0);
  const MagBlockCoordSys* src_coord_N = (MagBlockCoordSys*)getCoordSys(numBlocks()-1);

  if (m_mag_geom_type == "slab") {
    RealVect lowerMappedCoord;
    RealVect upperMappedCoord;

    for (int dir=0; dir<SpaceDim; ++dir) {
        lowerMappedCoord[dir] = src_coord_0->lowerMappedCoordinate(dir);
        upperMappedCoord[dir] = src_coord_N->upperMappedCoordinate(dir);
    }
  
    RealVect lowerCoord = src_coord_0->realCoord(lowerMappedCoord);
    RealVect upperCoord = src_coord_N->realCoord(upperMappedCoord);

    if (a_x[TOROIDAL_DIR]<lowerCoord[TOROIDAL_DIR]) {
      double shift = lowerCoord[TOROIDAL_DIR] - a_x[TOROIDAL_DIR];
      a_x[TOROIDAL_DIR] = upperCoord[TOROIDAL_DIR] - shift;
    }

    if (a_x[TOROIDAL_DIR]>upperCoord[TOROIDAL_DIR]) {
      double shift = a_x[TOROIDAL_DIR] - upperCoord[TOROIDAL_DIR];
      a_x[TOROIDAL_DIR] =	lowerCoord[TOROIDAL_DIR] + shift;
    }
  }
  
  else if (m_mag_geom_type == "toroidal") {
    double phi_max = ((const ToroidalBlockCoordSys*)src_coord_0)->getToroidalWedgeFraction();
    phi_max *= 2. * Pi;
    double phi = atan2(a_x[1], a_x[0]);
    double r = sqrt(pow(a_x[0],2) + pow(a_x[1],2));
    
    if (phi < 0.0) {
      double phi_reflected = phi + phi_max;
      a_x[0] = r * cos(phi_reflected);
      a_x[1] = r * sin(phi_reflected);
    }

    if (phi > phi_max) {
      double phi_reflected = phi - phi_max;
      a_x[0] = r * cos(phi_reflected);
      a_x[1] = r * sin(phi_reflected);
    }
  }
  
  else {
    MayDay::Error("LogRectCoordSys::toroidalBlockRemapping(): remapping is not implemented for this geometry");
  }
  
}


void
LogRectCoordSys::enforcePoloidalCut(const RealVect& a_xiSrc,
                                    double& a_xiNW_dst,
                                    double& a_xiNE_dst,
                                    double& a_xiSW_dst,
                                    double& a_xiSE_dst,
                                    double& a_xi0_dst) const
{
   
   Vector<double> poloidalNodes(5,0);
   poloidalNodes[0] = a_xiNW_dst;
   poloidalNodes[1] = a_xiNE_dst;
   poloidalNodes[2] = a_xiSW_dst;
   poloidalNodes[3] = a_xiSE_dst;
   poloidalNodes[4] = a_xi0_dst;

   
   for (int n=0; n<5; ++n) {
      if (fabs(a_xiSrc[POLOIDAL_DIR] - poloidalNodes[n]) > Pi) {
         int sign = (a_xiSrc[POLOIDAL_DIR] - Pi > 0.) ? 1 : -1;
         poloidalNodes[n] += sign * 2.0 * Pi;
      }
   }
   
   a_xiNW_dst = poloidalNodes[0];
   a_xiNE_dst = poloidalNodes[1];
   a_xiSW_dst = poloidalNodes[2];
   a_xiSE_dst = poloidalNodes[3];
   a_xi0_dst = poloidalNodes[4];
}


void
LogRectCoordSys::getInterpolationCoefficients(Vector<Real>&    a_coeff,
                                              const RealVect&  a_xi0_dst,
                                              const IntVect&   a_iv0_dst,
					      const RealVect&  a_dx_dst,
                                              const int        a_order) const

{
   Real h = a_dx_dst[POLOIDAL_DIR];
   
   Real cent = (a_iv0_dst[POLOIDAL_DIR] + 0.5) * a_dx_dst[POLOIDAL_DIR];
   Real lo = cent - h;
   Real hi = cent + h;

   Real coeffLo(0.0);
   Real coeffCent(1.0);
   Real coeffHi(0.0);
   
   Real point = a_xi0_dst[POLOIDAL_DIR];

   if (a_order == 3) {
      if (abs(cent - point) > 1.0e-10) {

         Real d0 = lo - point;
         Real d1 = cent - point;
         Real d2 = hi - point;
   
         coeffLo = d1*d2/(d0-d1)/(d0-d2);
         coeffCent = d2*d0/(d1-d2)/(d1-d0);
         coeffHi = d0*d1/(d2-d0)/(d2-d1);
      }
   }
   
   else if (a_order == 2) {
      if (point > cent) {
         coeffLo = 0.0;
         coeffCent = 1.0 - (point - cent)/h;
         coeffHi = (point - cent)/h;
      }
      else {
         coeffLo = (cent - point)/h;
         coeffCent =  1.0 - (cent - point)/h;
         coeffHi = 0.0;
      }
   }
   
   else {
      MayDay::Error("getInterpolationCoefficients:: only order = 2 or 3 is currently supported");
   }
   
   a_coeff[0] = coeffLo;
   a_coeff[1] = coeffCent;
   a_coeff[2] = coeffHi;

   
}

/*
 This is an experimental implementation (presently not used). It seems to be 
 only first order in the poloidal dir, but includes some radial variations
 in the remapping coefficeints (by taking advantage of the nodal values).
 */
void
LogRectCoordSys::getInterpolationCoefficients(Vector<Real>&   a_coeff,
                                              const Real&     a_xiNW_dst,
                                              const Real&     a_xiNE_dst,
                                              const Real&     a_xiSW_dst,
                                              const Real&     a_xiSE_dst,
                                              const IntVect&  a_iv0_dst,
					      const RealVect& a_dx_dst) const
{
   //High and low poloida_dir ends of the central dst cell
   double hi_end = (a_iv0_dst[POLOIDAL_DIR] + 1) * a_dx_dst[POLOIDAL_DIR];
   double lo_end = (a_iv0_dst[POLOIDAL_DIR]) * a_dx_dst[POLOIDAL_DIR];
   
   //second-order interpolation stencil
   double interpCoeffHi(0.0);
   double interpCoeffCent(a_dx_dst[POLOIDAL_DIR]);
   double interpCoeffLo(0.0);
   
   //Handle hi and cent cells
   if (((a_xiNW_dst - hi_end) >= 0) && ((a_xiNE_dst - hi_end) > 0)) {
      
      double area1 = (a_xiNW_dst - hi_end);
      double area2 = (a_xiNE_dst - hi_end);
      
      interpCoeffHi += (area2 + area1)/2.0 ;
      
   }
   
   else if (((a_xiNW_dst - hi_end) <= 0) && ((a_xiNE_dst - hi_end) < 0)) {
      
      double area1 = (hi_end - a_xiNW_dst);
      double area2 = (hi_end - a_xiNE_dst);
      
      interpCoeffCent -= (area2 + area1)/2.0 ;
   }
   
   else {
      
      double rad_inters_frac = abs(a_xiNE_dst - hi_end)/abs(a_xiNW_dst-hi_end);
      double rad_inters_dist = 1.0 / (1.0 + rad_inters_frac);
      
      double area1 = 0.5 * (a_xiNW_dst - hi_end) * rad_inters_dist;
      if (area1 < 0) {
         interpCoeffCent += area1;
      }
      else {
         interpCoeffHi += area1;
      }
      
      double area2 = 0.5 * (a_xiNE_dst - hi_end) * (1.0 - rad_inters_dist);
      if (area2 < 0) {
         interpCoeffCent += area2;
      }
      else {
         interpCoeffHi += area2;
      }
      
   }
   
   //Handle lo and cent cells
   if (((a_xiSW_dst - lo_end) <= 0) && ((a_xiSE_dst - lo_end) < 0)) {
      
      double area1 = (lo_end - a_xiSW_dst);
      double area2 = (lo_end - a_xiSE_dst);
      
      interpCoeffLo += (area2 + area1)/2.0 ;
      
   }
   
   else if (((a_xiSW_dst - lo_end) >= 0) && ((a_xiSE_dst - lo_end) > 0)) {
      
      double area1 = (a_xiSW_dst - lo_end);
      double area2 = (a_xiSE_dst - lo_end);
      
      interpCoeffCent -= (area2 + area1)/2.0 ;
   }
   
   else {
      
      double rad_inters_frac = abs(a_xiSE_dst - lo_end)/abs(a_xiSW_dst-lo_end);
      double rad_inters_dist = 1.0 / (1.0 + rad_inters_frac);
      
      double area1 = 0.5 * (lo_end - a_xiSW_dst) * rad_inters_dist;
      if (area1 < 0) {
         interpCoeffCent += area1;
      }
      else {
         interpCoeffLo += area1;
      }
      
      double area2 = 0.5 * (lo_end - a_xiSE_dst) * (1.0 - rad_inters_dist);
      if (area2 < 0) {
         interpCoeffCent += area2;
      }
      else {
         interpCoeffLo += area2;
      }
      
   }
   
   //Use this overwrite for the case where the ghost cell is
   //conformal to the destination cell (e.g., By = 0 case)
   if ((abs(a_xiSW_dst - lo_end) < 1.0e-10) && (abs(a_xiSE_dst - lo_end) < 1.0e-10)
       && (abs(a_xiNW_dst - hi_end) < 1.0e-10) && (abs(a_xiNE_dst - hi_end) < 1.0e-10)) {
      
      interpCoeffLo = 0.0;
      interpCoeffHi = 0.0;
      interpCoeffCent = 1.0;
   }

   double interpSum = interpCoeffLo + interpCoeffCent + interpCoeffHi;

   a_coeff[0] = interpCoeffLo/interpSum;
   a_coeff[1] = interpCoeffCent/interpSum;
   a_coeff[2] = interpCoeffHi/interpSum;
}

#endif

#include "NamespaceFooter.H"
