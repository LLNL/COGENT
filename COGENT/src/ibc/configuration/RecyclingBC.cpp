#include "RecyclingBC.H"
#include "Directions.H"
#include "RecyclingBCF_F.H"
#include "DataArray.H"
#include "inspect.H"

#include "NamespaceHeader.H"

const string RecyclingBC::pp_name = "recycling_boundary";

RecyclingBC::RecyclingBC(const BoundaryBoxLayoutPtr&  a_bdry_layout,
                         const ParmParse&             a_pp)

: m_bdry_layout(a_bdry_layout)
{
   parseParameters( a_pp );
}

void RecyclingBC::applyBC(LevelData<FluxBox>&         a_outflux,
                          const LevelData<FluxBox>&   a_influx,
                          const Real                  a_R,
                          const Real                  a_time ) const
{
   /*
    This BC reflects the incoming normal flux with a_R weight.
    */
   
   const int& dir( m_bdry_layout->dir() );
   const Side::LoHiSide& side( m_bdry_layout->side() );
   const int iSIDE(side);
   
   CH_assert(a_outflux.nComp() == 1);
   CH_assert(a_influx.nComp() == 1);
   const DisjointBoxLayout& bdry_grids( m_bdry_layout->disjointBoxLayout() );

   for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

      const Box fill_box( bdry_grids[dit] );

      //We only need one layer of ghosts in direction dir to loop over
      Box boundaryBox(bdry_grids[dit]);
      if (boundaryBox.size(dir) > 1) {
         boundaryBox.growDir(dir, side, 1 - boundaryBox.size(dir) );
      }
      
      const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
      const FArrayBox& this_influx(a_influx[interior_dit][dir]);
      
      FArrayBox& this_outflux(a_outflux[interior_dit][dir]);

      FORT_REFLECT_FLUX(CHF_FRA1(this_outflux,0),
                        CHF_BOX(boundaryBox),
                        CHF_CONST_FRA1(this_influx,0),
                        CHF_CONST_REAL(a_R),
                        CHF_CONST_INT(dir),
                        CHF_CONST_INT(iSIDE) );
   }
}

void RecyclingBC::applyBC(EllipticOpBC&             a_elliptic_bc,
                          const LevelData<FluxBox>& a_face_val,
                          const Real                a_R,
                          const Real                a_time ) const
{
   /*
    Set ellipticOp BCs based on the boundary values of a_face_val
    */

   CH_assert(a_face_val.ghostVect() == IntVect::Zero);
  
   const int& dir( m_bdry_layout->dir() );
   const Side::LoHiSide& side( m_bdry_layout->side() );
   const int iSIDE(side);
   
   CH_assert(a_face_val.nComp() == 1);
   const DisjointBoxLayout& bdry_grids( m_bdry_layout->disjointBoxLayout() );
   
   const int& block( m_bdry_layout->block() );
   
   // Create a single-layer cell-centered boundary box
   const ProblemDomain& bdry_domain = bdry_grids.physDomain();
   Box bdry_box = bdry_domain.domainBox();
   if (bdry_box.size(dir) > 1) {
       bdry_box.growDir(dir, side, 1 - bdry_box.size(dir) );
   }
   bdry_box.shiftHalf(dir, -sign(side));
   
   // Create a local (this proc) temporary data_array storage of BCs
   FArrayBox data_array_local(bdry_box, 1);
   data_array_local.setVal(0.);
   
   for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

      //We only need one layer of ghosts in direction dir to loop over
      Box this_boundaryBox = bdry_grids[dit];
      if (this_boundaryBox.size(dir) > 1) {
         this_boundaryBox.growDir(dir, side, 1 - this_boundaryBox.size(dir) );
      }
      this_boundaryBox.shiftHalf(dir, -sign(side));

      
      const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
      const FArrayBox& this_face_val(a_face_val[interior_dit][dir]);

      FArrayBox this_bc_data(this_boundaryBox, 1);

      for (BoxIterator bit(this_boundaryBox); bit.ok(); ++bit) {
	IntVect iv = bit();
	iv[dir] = (this_face_val.box()).sideEnd(side)[dir];
	this_bc_data(bit(),0) = this_face_val(iv,0);
      }      
      
      // Fill a part of the data_array object
      data_array_local.copy(this_bc_data, this_boundaryBox);
   }
   
   // Create a global temporary data_array storage of BCs
   FArrayBox data_array_global(bdry_box, 1);
   data_array_global.setVal(0.);

#ifdef CH_MPI
   //Broadcast data_array object across processors
   int size = bdry_box.numPts();
   MPI_Allreduce(data_array_local.dataPtr(), data_array_global.dataPtr(), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
   data_array_global.copy(data_array_local);
#endif
   
   // Multiply by recycling coeff
   data_array_global.mult(a_R);

   RefCountedPtr<GridFunction> bc_func = a_elliptic_bc.getBCFunction(block, dir, side);
   CH_assert((typeid(*bc_func) == typeid(DataArray)));
   ((RefCountedPtr<DataArray>)bc_func)->setData(data_array_global, false);
}

void RecyclingBC::applyBC(LevelData<FArrayBox>&       a_dst,
                          const LevelData<FArrayBox>& a_src,
                          const Real                  a_R,
                          const Real                  a_time ) const
{
   /*
    This BC reflects extrapolates a_src into the ghost cells of a_dst with a_R weight
    */
   
   const int& dir( m_bdry_layout->dir() );
   const Side::LoHiSide& side( m_bdry_layout->side() );
   const int iSIDE(side);
      
   CH_assert(a_src.nComp() == a_dst.nComp());
   const DisjointBoxLayout& bdry_grids( m_bdry_layout->disjointBoxLayout() );

   for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

      const Box fill_box( bdry_grids[dit] );

      //We only need one layer of ghosts in direction dir to loop over
      Box boundaryBox(bdry_grids[dit]);
      if (boundaryBox.size(dir) > 1) {
         boundaryBox.growDir(dir, side, 1 - boundaryBox.size(dir) );
      }
      
      const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
      const FArrayBox& this_src(a_src[interior_dit]);

      FArrayBox& this_dst(a_dst[interior_dit]);

      FORT_EXTRAP_TO_GHOSTS(CHF_FRA(this_dst),
                            CHF_BOX(boundaryBox),
                            CHF_CONST_FRA(this_src),
                            CHF_CONST_REAL(a_R),
                            CHF_CONST_INT(dir),
                            CHF_CONST_INT(iSIDE) );


   }
}

void
RecyclingBC::parseParameters(const ParmParse& a_pp)
{
   //Parse input here
}

#include "NamespaceFooter.H"

