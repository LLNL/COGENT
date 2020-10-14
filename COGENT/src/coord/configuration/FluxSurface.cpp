#include "FluxSurface.H"
#include "FluxSurfaceF_F.H"
#include "MagCoordSys.H"
#include "SingleNullBlockCoordSys.H"
#include "SingleNullCoordSys.H"
#include "SNCoreCoordSys.H"

#ifdef CH_MULTIDIM

#undef CH_SPACEDIM
#include "ReductionOps.H.multidim"
#include "ReductionCopier.H.multidim"
#include "SpreadingCopier.H.multidim"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif
#define CH_SPACEDIM CFG_DIM

#include "Directions.H"

#include "NamespaceHeader.H"



FluxSurface::FluxSurface( const MagGeom& a_geom,
                          const bool     a_shell_average )
   : m_magnetic_geometry(&a_geom),
     m_shell_average(a_shell_average)
{
   m_volume.define(m_magnetic_geometry->grids(), 1, 2*IntVect::Unit);
   if (m_shell_average) {
      m_magnetic_geometry->getCellVolumes(m_volume);
   }
   else {
      LevelData<FArrayBox> face_areas(m_magnetic_geometry->grids(), SpaceDim, 2*IntVect::Unit);
      m_magnetic_geometry->getPointwiseFaceAreas(face_areas);
      for (DataIterator dit(face_areas.dataIterator()); dit.ok(); ++dit) {
        m_volume[dit].copy(face_areas[dit],RADIAL_DIR, 0, 1);
      }
   }
   
#if CFG_DIM == 2
   adjCellLo(m_grids, m_magnetic_geometry->grids(), POLOIDAL_DIR, -1);
#else
   adjCellLo(m_grids2D, m_magnetic_geometry->grids(), POLOIDAL_DIR, -1);
   adjCellLo(m_grids, m_grids2D, TOROIDAL_DIR, -1);
#endif

   m_areas.define(m_grids, 1, IntVect::Zero);
   
   const MagCoordSys& sn_coord_sys( *(m_magnetic_geometry->getCoordSys()) );
   m_single_null = ( typeid(sn_coord_sys) == typeid(SingleNullCoordSys) );

   if (m_single_null) {
      m_areas_core.define(m_grids, 1, IntVect::Zero);
      m_areas_pf.define(m_grids, 1, IntVect::Zero);
   } 

   computeAreas();

}


void
FluxSurface::computeAreas()
{
   LevelData<FArrayBox> ones(m_magnetic_geometry->grids(), 1, IntVect::Zero);
   DataIterator dit = ones.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
        ones[dit].setVal(1.);
   }
   
   integrate(ones, m_areas);

   //get flux surface areas for the core and private flux regions
   if ( m_single_null ) {
    
     LevelData<FArrayBox> ones_tmp(m_magnetic_geometry->grids(), 1, IntVect::Zero);

     zeroBlockData(PF, ones, ones_tmp);
     integrate(ones_tmp, m_areas_core);

     zeroBlockData(CORE, ones, ones_tmp);
     integrate(ones_tmp, m_areas_pf);
   }

}


void
FluxSurface::sum( const LevelData<FArrayBox>& a_src,
                  LevelData<FArrayBox>&       a_dst ) const
{
   CH_assert(a_src.nComp() == a_dst.nComp());

   const DisjointBoxLayout& dst_grids = a_dst.getBoxes();
   const DisjointBoxLayout& src_grids = a_src.getBoxes();
   const ProblemDomain& problem_domain = src_grids.physDomain();

   // Initialize the destination, since SumOp does not do that for us.
   DataIterator dit = a_dst.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_dst[dit].setVal(0.);
   }

#if CFG_DIM == 2

   int poloidal_dir = POLOIDAL_DIR;
   // Define ReductionCopier to compute intersections (sum in the y-direction)
   ReductionCopier reduceCopier(src_grids, dst_grids, problem_domain, poloidal_dir);

   SumOp op(poloidal_dir);
   op.scale = 1.0;

   // Do the summing operation -- sums data in src along the poloidal direction
   // and places the result in dst
   a_src.copyTo(a_src.interval(), a_dst, a_dst.interval(), reduceCopier, op);

#else

   // Initialize the destination, since SumOp does not do that for us.
   LevelData<FArrayBox> tmp_2d(m_grids2D, a_src.nComp(), IntVect::Zero);
   for (DataIterator dit(m_grids2D); dit.ok(); ++dit) {
      tmp_2d[dit].setVal(0.);
   }
   
   // Define ReductionCopier to compute intersections (sum in the poloidal direction)
   ReductionCopier reduceCopier2D(src_grids, m_grids2D, problem_domain, POLOIDAL_DIR);
   
   SumOp op2D(POLOIDAL_DIR);
   op2D.scale = 1.0;

   // Do the summing operation -- sums data in src along the poloidal direction
   a_src.copyTo(a_src.interval(), tmp_2d, tmp_2d.interval(), reduceCopier2D, op2D);

   // Define ReductionCopier to compute intersections (sum in the toroidal  direction)
   ReductionCopier reduceCopier3D(m_grids2D, dst_grids, problem_domain, TOROIDAL_DIR);
   
   SumOp op3D(TOROIDAL_DIR);
   op3D.scale = 1.0;

   // Do the summing operation -- sums data in src along the toroidal direction
   tmp_2d.copyTo(tmp_2d.interval(), a_dst, a_dst.interval(), reduceCopier3D, op3D);
   
#endif
}



void
FluxSurface::spread( const LevelData<FArrayBox>& a_src,
                     LevelData<FArrayBox>&       a_dst ) const
{
   CH_assert(a_src.nComp() == a_dst.nComp());

   const DisjointBoxLayout& dst_grids = a_dst.getBoxes();
   const DisjointBoxLayout& src_grids = a_src.getBoxes();
   const ProblemDomain& problem_domain = dst_grids.physDomain();

#if CFG_DIM == 2

   int poloidal_dir = POLOIDAL_DIR;
   // Define SpreadingCopier to spread in the poloidal direction
   SpreadingCopier spreadCopier(src_grids, dst_grids, problem_domain, poloidal_dir);

   const SpreadingOp spreadOp(poloidal_dir);

   // Do the spreading
   a_src.copyTo(a_src.interval(), a_dst, a_dst.interval(), spreadCopier, spreadOp);

#else
   
   LevelData<FArrayBox> tmp_2d(m_grids2D, a_src.nComp(), IntVect::Zero);
   
   // Define SpreadingCopier to spread in the poloidal direction
   SpreadingCopier spreadCopier2D(src_grids, m_grids2D, problem_domain, TOROIDAL_DIR);
   
   const SpreadingOp spreadOp2D(TOROIDAL_DIR);
   
   // Do the spreading in the poloidal direction
   a_src.copyTo(a_src.interval(), tmp_2d, tmp_2d.interval(), spreadCopier2D, spreadOp2D);
   
   // Define SpreadingCopier to spread in the toroidal  direction
   SpreadingCopier spreadCopier3D(m_grids2D, dst_grids, problem_domain, POLOIDAL_DIR);
   
   const SpreadingOp spreadOp3D(POLOIDAL_DIR);
   
   // Do the remaining spreading in the toroidal direction
   tmp_2d.copyTo(tmp_2d.interval(), a_dst, a_dst.interval(), spreadCopier3D, spreadOp3D);

#endif
}



void
FluxSurface::integrate( const LevelData<FArrayBox>& a_src,
                        LevelData<FArrayBox>&       a_dst ) const
{
   CH_assert(a_src.nComp() == a_dst.nComp());

   int ncomp = a_src.nComp();

   LevelData<FArrayBox> tmp(m_magnetic_geometry->grids(), ncomp, IntVect::Zero);
   DataIterator dit = tmp.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      tmp[dit].copy(a_src[dit]);
      for (int comp=0; comp<ncomp; ++comp) {
         tmp[dit].mult(m_volume[dit],0,comp,1);
      }
   }

   sum(tmp, a_dst);
}



void
FluxSurface::average( const LevelData<FArrayBox>& a_src,
                      LevelData<FArrayBox>&       a_dst ) const
{

   CH_assert(a_src.nComp() == a_dst.nComp());

   integrate(a_src, a_dst);

   int ncomp = a_src.nComp();

   DataIterator dit = a_dst.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      for (int comp=0; comp<ncomp; ++comp) {
         a_dst[dit].divide(m_areas[dit],0,comp,1);
      }
   }
}

void
FluxSurface::averageAndSpread( const LevelData<FArrayBox>& a_src,
                               LevelData<FArrayBox>&       a_dst ) const
{
   LevelData<FArrayBox> flux_aver_tmp(m_grids, a_src.nComp(), a_src.ghostVect());
   average(a_src, flux_aver_tmp);
   spread(flux_aver_tmp, a_dst);  

   if (m_single_null) {

     int ncomp = a_src.nComp();
     CH_assert(a_src.nComp() == a_dst.nComp());

     LevelData<FArrayBox> src_tmp(a_src.disjointBoxLayout(),ncomp, a_src.ghostVect());
     LevelData<FArrayBox> dst_tmp(a_dst.disjointBoxLayout(),ncomp, a_dst.ghostVect());
     LevelData<FArrayBox> flux_aver_tmp(m_grids, ncomp, a_src.ghostVect());

     DataIterator dit = a_dst.dataIterator();

     //Perform the average-and-spreading in the core region
     zeroBlockData(PF, a_src, src_tmp);
     integrate(src_tmp, flux_aver_tmp);
     for (dit.begin(); dit.ok(); ++dit) {
      for (int comp=0; comp<ncomp; ++comp) {
         flux_aver_tmp[dit].divide(m_areas_core[dit],0,comp,1);
      }
     }
     spread(flux_aver_tmp, dst_tmp);  
     zeroBlockData(PF, dst_tmp, src_tmp); 
     //Add flux-surface averages from the core region
     for (dit.begin(); dit.ok(); ++dit) {
         a_dst[dit].copy(src_tmp[dit]);
     }


     //Perform the average-and-spreading in the private flux region
     zeroBlockData(CORE, a_src, src_tmp);
     integrate(src_tmp, flux_aver_tmp);
     for (dit.begin(); dit.ok(); ++dit) {
      for (int comp=0; comp<ncomp; ++comp) {
         flux_aver_tmp[dit].divide(m_areas_pf[dit],0,comp,1);
      }
     }
     spread(flux_aver_tmp, dst_tmp);  
     zeroBlockData(SOL, dst_tmp, src_tmp); 
     zeroBlockData(CORE, src_tmp, dst_tmp); 
     //Add flux-surface averages from the prvate-flux region
     for (dit.begin(); dit.ok(); ++dit) {
         a_dst[dit].plus(dst_tmp[dit]);
     }

   }

}



void
FluxSurface::add( const LevelData<FArrayBox>& a_flux_surface_data,
                  LevelData<FArrayBox>&       a_data              ) const
{
   const DisjointBoxLayout& grids = a_data.getBoxes();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FORT_ADD_FLUX_SURFACE_ARRAY(CHF_BOX(grids[dit]),
                                  CHF_CONST_FRA1(a_flux_surface_data[dit],0),
                                  CHF_FRA(a_data[dit]));
   }
}



void
FluxSurface::subtract( const LevelData<FArrayBox>& a_flux_surface_data,
                       LevelData<FArrayBox>&       a_data              ) const
{
   const DisjointBoxLayout& grids = a_data.getBoxes();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FORT_SUBTRACT_FLUX_SURFACE_ARRAY(CHF_BOX(grids[dit]),
                                       CHF_CONST_FRA1(a_flux_surface_data[dit],0),
                                       CHF_FRA(a_data[dit]));
   }
}



void
FluxSurface::multiply( const LevelData<FArrayBox>& a_flux_surface_data,
                       LevelData<FArrayBox>&       a_data              ) const
{
   const DisjointBoxLayout& grids = a_data.getBoxes();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FORT_MULTIPLY_FLUX_SURFACE_ARRAY(CHF_BOX(grids[dit]),
                                       CHF_CONST_FRA1(a_flux_surface_data[dit],0),
                                       CHF_FRA(a_data[dit]));
   }
}



void
FluxSurface::divide( const LevelData<FArrayBox>& a_flux_surface_data,
                     LevelData<FArrayBox>&       a_data              ) const
{
   const DisjointBoxLayout& grids = a_data.getBoxes();

   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      FORT_DIVIDE_FLUX_SURFACE_ARRAY(CHF_BOX(grids[dit]),
                                     CHF_CONST_FRA1(a_flux_surface_data[dit],0),
                                     CHF_FRA(a_data[dit]));
   }
}


void
FluxSurface::zeroBlockData( const int                   region_type,
                            const LevelData<FArrayBox>& a_src,
                            LevelData<FArrayBox>&       a_dst       ) const
{
   const DisjointBoxLayout& grids = a_src.getBoxes();
   const MagCoordSys& coord_sys = *(m_magnetic_geometry->getCoordSys());

   CH_assert( typeid(coord_sys) == typeid(SingleNullCoordSys) );
   
   DataIterator dit = grids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
     a_dst[dit].copy(a_src[dit]);
     int block_number = coord_sys.whichBlock( grids[dit] );
     switch (region_type) 
      { 
      case CORE:
	if ( ((const SingleNullCoordSys&)coord_sys).isCORE(block_number) ) {
	  a_dst[dit].setVal(0.0);
	}
        break;
      case PF:
	if ( ((const SingleNullCoordSys&)coord_sys).isPF(block_number) ) {
          a_dst[dit].setVal(0.0);
        }
	break;    
      case SOL:
	if ( ((const SingleNullCoordSys&)coord_sys).isSOL(block_number) ) {
          a_dst[dit].setVal(0.0);
        }
        break;    
      default:
         MayDay::Error("FluxSurface::zeroBlockData(): Invalid region_type encountered");
      }
   }
}

#include "NamespaceFooter.H"

#endif // end ifdef Multidim
