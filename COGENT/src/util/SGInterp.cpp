#include "SGInterp.H"

#include "NamespaceHeader.H"

void SGInterp::coarsen1DCC( FArrayBox&        a_dst,
                            const Box&        a_dst_bx,
                            const FArrayBox&  a_src,
                            const Box&        a_src_bx,
                            const int         a_dir )
{
  CH_assert(a_dst.nComp() == a_src.nComp());
  CH_assert( (a_dir >= 0) && (a_dir < SpaceDim) );

  static const int npts_interp = 6;

  const IntVect& src_lo(a_src_bx.smallEnd());
  const IntVect& src_hi(a_src_bx.bigEnd());
  const IntVect& dst_lo(a_dst_bx.smallEnd());
  const IntVect& dst_hi(a_dst_bx.bigEnd());

  if ((src_lo == dst_lo) && (src_hi == dst_hi)) {
    a_dst.copy(a_src);
    return;
  }

  for (int d=0; d<SpaceDim; d++) {
    if (d != a_dir) {
      CH_assert(src_lo[d] == dst_lo[d]);
      CH_assert(src_hi[d] == dst_hi[d]);
    }
  }

  int n_src = src_hi[a_dir] - src_lo[a_dir] + 1;
  int n_dst = dst_hi[a_dir] - dst_lo[a_dir] + 1;
  CH_assert(n_src > n_dst);

  double fac = ((double) n_src) / ((double) n_dst);
  int stride = (int) fac;

  /* interp coefficients */
  static std::vector<double> c(npts_interp);
  c[0] = c[5] = 3.0/256.0;
  c[1] = c[4] = -25.0/256.0;
  c[2] = c[3] = 150.0/256.0;

  a_dst.setVal(0.0);
  for (BoxIterator bit(a_dst_bx); bit.ok(); ++bit) {

    IntVect iv_dst(bit());

    std::vector<IntVect> iv_src(npts_interp, bit());
    iv_src[2][a_dir] = src_lo[a_dir]
                       + (iv_dst[a_dir]-dst_lo[a_dir])*stride + (stride/2-1);
    iv_src[1][a_dir] = iv_src[2][a_dir]-1;
    iv_src[0][a_dir] = iv_src[2][a_dir]-2;
    iv_src[3][a_dir] = iv_src[2][a_dir]+1;
    iv_src[4][a_dir] = iv_src[2][a_dir]+2;
    iv_src[5][a_dir] = iv_src[2][a_dir]+3;

    for (int n = 0; n < a_dst.nComp(); n++) {
      for (int i = 0; i < npts_interp; i++) {
        a_dst(iv_dst,n) += c[i] * a_src(iv_src[i],n);
      }
    }

  }

  return;
}

void SGInterp::coarsen1DFC( FArrayBox&        a_dst,
                            const Box&        a_dst_bx,
                            const FArrayBox&  a_src,
                            const Box&        a_src_bx,
                            const int         a_dir )
{
  CH_assert(a_dst.nComp() == a_src.nComp());
  CH_assert( (a_dir >= 0) && (a_dir < SpaceDim) );
  CH_assert(a_src_bx.type(a_dir) == IndexType::NODE);
  CH_assert(a_dst_bx.type(a_dir) == IndexType::NODE);

  Box src_bx_cc = enclosedCells(a_src_bx, a_dir);
  Box dst_bx_cc = enclosedCells(a_dst_bx, a_dir);

  const IntVect& src_lo(src_bx_cc.smallEnd());
  const IntVect& src_hi(src_bx_cc.bigEnd());
  const IntVect& dst_lo(dst_bx_cc.smallEnd());
  const IntVect& dst_hi(dst_bx_cc.bigEnd());

  if ((src_lo == dst_lo) && (src_hi == dst_hi)) {
    a_dst.copy(a_src);
    return;
  }

  for (int d=0; d<SpaceDim; d++) {
    if (d != a_dir) {
      CH_assert(src_lo[d] == dst_lo[d]);
      CH_assert(src_hi[d] == dst_hi[d]);
    }
  }

  int n_src = src_hi[a_dir] - src_lo[a_dir] + 1;
  int n_dst = dst_hi[a_dir] - dst_lo[a_dir] + 1;
  CH_assert(n_src > n_dst);

  double fac = ((double) n_src) / ((double) n_dst);
  int stride = (int) fac;

  a_dst.setVal(0.0);
  for (BoxIterator bit(a_dst_bx); bit.ok(); ++bit) {

    IntVect iv_dst(bit());

    IntVect iv_src(bit());
    iv_src[a_dir] = src_lo[a_dir] + (iv_dst[a_dir]-dst_lo[a_dir])*stride;

    for (int n = 0; n < a_dst.nComp(); n++) {
      a_dst(iv_dst,n) = a_src(iv_src,n);
    }

  }

  return;
}

void SGInterp::refine1DCC(  FArrayBox&        a_dst,
                            const Box&        a_dst_bx,
                            const FArrayBox&  a_src,
                            const Box&        a_src_bx,
                            const int         a_dir )
{
  CH_assert(a_dst.nComp() == a_src.nComp());
  CH_assert( (a_dir >= 0) && (a_dir < SpaceDim) );

  static const int npts_interp = 7;

  const IntVect& src_lo(a_src_bx.smallEnd());
  const IntVect& src_hi(a_src_bx.bigEnd());
  const IntVect& dst_lo(a_dst_bx.smallEnd());
  const IntVect& dst_hi(a_dst_bx.bigEnd());

  for (int d=0; d<SpaceDim; d++) {
    if (d != a_dir) {
      CH_assert(src_lo[d] == dst_lo[d]);
      CH_assert(src_hi[d] == dst_hi[d]);
    }
  }

  int n_src = src_hi[a_dir] - src_lo[a_dir] + 1;
  int n_dst = dst_hi[a_dir] - dst_lo[a_dir] + 1;
  CH_assert(n_src < n_dst);

  double fac = ((double) n_dst) / ((double) n_src);
  int stride = (int) fac;

  a_dst.setVal(0.0);
  for (BoxIterator bit(a_src_bx); bit.ok(); ++bit) {

    std::vector<IntVect> iv_src(npts_interp, bit());
    iv_src[0][a_dir] -= 3;
    iv_src[1][a_dir] -= 2;
    iv_src[2][a_dir] -= 1;
    iv_src[4][a_dir] += 1;
    iv_src[5][a_dir] += 2;
    iv_src[6][a_dir] += 3;

    for (int j = 0; j < stride; j++) {

      double xi = -0.5 + ((double)j + 0.5) * (1.0/fac);

      IntVect iv_dst(bit());
      iv_dst[a_dir] = dst_lo[a_dir]
                      + stride*(iv_src[3][a_dir]-src_lo[a_dir]) + j;
      CH_assert(a_dst_bx.contains(iv_dst));

      std::vector<double> c(npts_interp);
      c[0] = ( 1.0/720.0) * (xi-3)*(xi-2)*(xi-1)*xi*(xi+1)*(xi+2);
      c[1] = (-1.0/120.0) * (xi-3)*(xi-2)*(xi-1)*xi*(xi+1)*(xi+3);
      c[2] = ( 1.0/ 48.0) * (xi-3)*(xi-2)*(xi-1)*xi*(xi+2)*(xi+3);
      c[3] = ( 1.0/ 36.0) * (36 - xi*xi*(xi*xi-7)*(xi*xi-7));
      c[4] = ( 1.0/ 48.0) * (xi-3)*(xi-2)*xi*(xi+1)*(xi+2)*(xi+3);
      c[5] = (-1.0/120.0) * (xi-3)*(xi-1)*xi*(xi+1)*(xi+2)*(xi+3);
      c[6] = ( 1.0/720.0) * (xi-2)*(xi-1)*xi*(xi+1)*(xi+2)*(xi+3);

      for (int n = 0; n < a_dst.nComp(); n++) {
        for (int i = 0; i < npts_interp; i++) {
          a_dst(iv_dst,n) += c[i] * a_src(iv_src[i],n);
        }
      }
    }

  }

  return;
}

void SGInterp::refine1DFC(  FArrayBox&        a_dst,
                            const Box&        a_dst_bx,
                            const FArrayBox&  a_src,
                            const Box&        a_src_bx,
                            const int         a_dir )
{
  CH_assert(a_dst.nComp() == a_src.nComp());
  CH_assert( (a_dir >= 0) && (a_dir < SpaceDim) );
  CH_assert(a_src_bx.type(a_dir) == IndexType::NODE);
  CH_assert(a_dst_bx.type(a_dir) == IndexType::NODE);

  static const int npts_interp = 7;

  Box src_bx_cc = enclosedCells(a_src_bx, a_dir);
  Box dst_bx_cc = enclosedCells(a_dst_bx, a_dir);

  const IntVect& src_lo(src_bx_cc.smallEnd());
  const IntVect& src_hi(src_bx_cc.bigEnd());
  const IntVect& dst_lo(dst_bx_cc.smallEnd());
  const IntVect& dst_hi(dst_bx_cc.bigEnd());

  for (int d=0; d<SpaceDim; d++) {
    if (d != a_dir) {
      CH_assert(src_lo[d] == dst_lo[d]);
      CH_assert(src_hi[d] == dst_hi[d]);
    }
  }

  int n_src = src_hi[a_dir] - src_lo[a_dir] + 1;
  int n_dst = dst_hi[a_dir] - dst_lo[a_dir] + 1;
  CH_assert(n_src < n_dst);

  double fac = ((double) n_dst) / ((double) n_src);
  int stride = (int) fac;

  a_dst.setVal(0.0);
  for (BoxIterator bit(a_src_bx); bit.ok(); ++bit) {

    std::vector<IntVect> iv_src(npts_interp, bit());
    iv_src[0][a_dir] -= 3;
    iv_src[1][a_dir] -= 2;
    iv_src[2][a_dir] -= 1;
    iv_src[4][a_dir] += 1;
    iv_src[5][a_dir] += 2;
    iv_src[6][a_dir] += 3;

    for (int j = -stride/2; j <= stride/2; j++) {

      double xi = -0.5 + ((double)(j+stride/2)) * (1.0/fac);

      IntVect iv_dst(bit());
      iv_dst[a_dir] = dst_lo[a_dir]
                      + stride*(iv_src[3][a_dir]-src_lo[a_dir]) + j;

      if (a_dst_bx.contains(iv_dst)) {

        std::vector<double> c(npts_interp);
        c[0] = ( 1.0/720.0) * (xi-3)*(xi-2)*(xi-1)*xi*(xi+1)*(xi+2);
        c[1] = (-1.0/120.0) * (xi-3)*(xi-2)*(xi-1)*xi*(xi+1)*(xi+3);
        c[2] = ( 1.0/ 48.0) * (xi-3)*(xi-2)*(xi-1)*xi*(xi+2)*(xi+3);
        c[3] = ( 1.0/ 36.0) * (36 - xi*xi*(xi*xi-7)*(xi*xi-7));
        c[4] = ( 1.0/ 48.0) * (xi-3)*(xi-2)*xi*(xi+1)*(xi+2)*(xi+3);
        c[5] = (-1.0/120.0) * (xi-3)*(xi-1)*xi*(xi+1)*(xi+2)*(xi+3);
        c[6] = ( 1.0/720.0) * (xi-2)*(xi-1)*xi*(xi+1)*(xi+2)*(xi+3);

        for (int n = 0; n < a_dst.nComp(); n++) {
          for (int i = 0; i < npts_interp; i++) {
            a_dst(iv_dst,n) += c[i] * a_src(iv_src[i],n);
          }
        }
      }
    }

  }

  return;
}

void SGInterp::interpolate( FArrayBox&        a_dst,
                            const Box&        a_dst_bx,
                            const FArrayBox&  a_src,
                            const Box&        a_src_int_bx )
{
  CH_assert(a_dst.nComp() == a_src.nComp());

  const IntVect n_src = a_src_int_bx.bigEnd() 
                        - a_src_int_bx.smallEnd() 
                        + IntVect::Unit;
  const IntVect n_dst = a_dst_bx.bigEnd() 
                        - a_dst_bx.smallEnd() 
                        + IntVect::Unit;

  FArrayBox* fab_from(nullptr);
  FArrayBox* fab_to(nullptr);
  Box* bx_from(nullptr);
  Box* bx_to(nullptr);
  IntVect n_from(IntVect::Zero), n_to(IntVect::Zero);

  const Box& src_bx_wg = a_src.box();
  CH_assert( src_bx_wg.contains(grow(a_src_int_bx,num_gpt)) );

  bx_to = new Box( src_bx_wg );
  fab_to = new FArrayBox( *bx_to, a_src.nComp() );
  fab_to->copy(a_src);
  n_to = n_src;

  IntVect gvec = num_gpt*IntVect::Unit;

  for (int dir = 0; dir < SpaceDim; dir++) {

    if (n_src[dir] < n_dst[dir]) {
      
      /* refine */
      n_from = n_to;
      bx_from = bx_to;
      fab_from = fab_to;

      n_to[dir] = n_dst[dir];

      IntVect ref_vec(IntVect::Unit);
      ref_vec[dir] = (int) n_to[dir]/n_from[dir];

      Box bx_from_no_ghosts = grow( *bx_from, -1*gvec );
      Box bx_to_no_ghosts = refine( bx_from_no_ghosts, ref_vec );
      gvec[dir] = 0;

      bx_to = new Box( grow( bx_to_no_ghosts, gvec ) );
      fab_to = new FArrayBox ( *bx_to, a_dst.nComp() );

      refine1DCC( *fab_to, 
                  *bx_to, 
                  *fab_from, 
                  grow( bx_from_no_ghosts, gvec ), 
                  dir);

      delete bx_from;
      delete fab_from;

    } else if (n_src[dir] > n_dst[dir]) {

      /* coarsen */
      n_from = n_to;
      bx_from = bx_to;
      fab_from = fab_to;

      n_to[dir] = n_dst[dir];

      IntVect ref_vec(IntVect::Unit);
      ref_vec[dir] = (int) n_from[dir]/n_to[dir];

      Box bx_from_no_ghosts = grow( *bx_from, -1*gvec );
      Box bx_to_no_ghosts = refine( bx_from_no_ghosts, ref_vec );
      gvec[dir] = 0;

      bx_to = new Box( grow( bx_to_no_ghosts, gvec ) );
      fab_to = new FArrayBox ( *bx_to, a_dst.nComp() );

      coarsen1DCC(  *fab_to, 
                    *bx_to, 
                    *fab_from, 
                    grow( bx_from_no_ghosts, gvec ), 
                    dir );

      delete bx_from;
      delete fab_from;

    } else {

      /* copy */
      n_from = n_to;
      bx_from = bx_to;
      fab_from = fab_to;

      n_to[dir] = n_dst[dir];

      IntVect ref_vec(IntVect::Unit);
      ref_vec[dir] = (int) n_to[dir]/n_from[dir];

      Box bx_from_no_ghosts = grow( *bx_from, -1*gvec );
      Box bx_to_no_ghosts = refine( bx_from_no_ghosts, ref_vec );
      gvec[dir] = 0;

      bx_to = new Box( grow( bx_to_no_ghosts, gvec ) );
      fab_to = new FArrayBox ( *bx_to, a_dst.nComp() );

      fab_to->copy( *fab_from, *bx_to );

      delete bx_from;
      delete fab_from;

    }

  }

  CH_assert( fab_to != nullptr );
  CH_assert( bx_to != nullptr );
  CH_assert( (*bx_to) == a_dst_bx );
  CH_assert( n_to == n_dst );

  a_dst.copy( *fab_to );
  delete fab_to;
  delete bx_to;

  return;
}

void SGInterp::interpolate( FluxBox&        a_dst,
                            const Box&      a_dst_bx,
                            const FluxBox&  a_src,
                            const Box&      a_src_int_bx )
{
  CH_assert(a_dst.nComp() == a_src.nComp());

  const IntVect n_src = a_src_int_bx.bigEnd() 
                        - a_src_int_bx.smallEnd() 
                        + IntVect::Unit;
  const IntVect n_dst = a_dst_bx.bigEnd() 
                        - a_dst_bx.smallEnd() 
                        + IntVect::Unit;

  for (int dim = 0; dim < SpaceDim; dim++) {

    FArrayBox* fab_from(nullptr);
    FArrayBox* fab_to(nullptr);
    Box* bx_from(nullptr);
    Box* bx_to(nullptr);
    IntVect n_from(IntVect::Zero), n_to(IntVect::Zero);
  
    const Box& src_bx_wg = a_src.box();
    CH_assert( src_bx_wg.contains(grow(a_src_int_bx,num_gpt)) );
  
    bx_to = new Box( src_bx_wg );
    fab_to = new FArrayBox( surroundingNodes(*bx_to,dim), a_src.nComp() );
    fab_to->copy(a_src[dim]);
    n_to = n_src;
  
    IntVect gvec = num_gpt*IntVect::Unit;
  
    for (int dir = 0; dir < SpaceDim; dir++) {
  
      if (n_src[dir] < n_dst[dir]) {
        
        /* refine */
        n_from = n_to;
        bx_from = bx_to;
        fab_from = fab_to;
  
        n_to[dir] = n_dst[dir];
  
        IntVect ref_vec(IntVect::Unit);
        ref_vec[dir] = (int) n_to[dir]/n_from[dir];
  
        Box bx_from_no_ghosts = grow( *bx_from, -1*gvec );
        Box bx_to_no_ghosts = refine( bx_from_no_ghosts, ref_vec );
        gvec[dir] = 0;
  
        bx_to = new Box( grow( bx_to_no_ghosts, gvec ) );
        fab_to = new FArrayBox( surroundingNodes(*bx_to,dim), a_dst.nComp() );
 
        if (dim == dir ) {
          refine1DFC( *fab_to, 
                      surroundingNodes(*bx_to,dim), 
                      *fab_from, 
                      surroundingNodes(grow(bx_from_no_ghosts,gvec), dim), 
                      dir);
        } else {
          refine1DCC( *fab_to, 
                      surroundingNodes(*bx_to,dim), 
                      *fab_from, 
                      surroundingNodes(grow(bx_from_no_ghosts,gvec), dim), 
                      dir);
        }
  
        delete bx_from;
        delete fab_from;
  
      } else if (n_src[dir] > n_dst[dir]) {
  
        /* coarsen */
        n_from = n_to;
        bx_from = bx_to;
        fab_from = fab_to;
  
        n_to[dir] = n_dst[dir];
  
        IntVect ref_vec(IntVect::Unit);
        ref_vec[dir] = (int) n_from[dir]/n_to[dir];
  
        Box bx_from_no_ghosts = grow( *bx_from, -1*gvec );
        Box bx_to_no_ghosts = refine( bx_from_no_ghosts, ref_vec );
        gvec[dir] = 0;
  
        bx_to = new Box( grow( bx_to_no_ghosts, gvec ) );
        fab_to = new FArrayBox( surroundingNodes(*bx_to,dim), a_dst.nComp() );
  
        if (dim == dir) {
          coarsen1DFC(  *fab_to, 
                        surroundingNodes(*bx_to,dim), 
                        *fab_from, 
                        surroundingNodes(grow(bx_from_no_ghosts,gvec), dim), 
                        dir );
        } else {
          coarsen1DCC(  *fab_to, 
                        surroundingNodes(*bx_to,dim), 
                        *fab_from, 
                        surroundingNodes(grow(bx_from_no_ghosts,gvec), dim), 
                        dir );
        }
  
        delete bx_from;
        delete fab_from;
  
      } else {
  
        /* copy */
        n_from = n_to;
        bx_from = bx_to;
        fab_from = fab_to;
  
        n_to[dir] = n_dst[dir];
  
        IntVect ref_vec(IntVect::Unit);
        ref_vec[dir] = (int) n_to[dir]/n_from[dir];
  
        Box bx_from_no_ghosts = grow( *bx_from, -1*gvec );
        Box bx_to_no_ghosts = refine( bx_from_no_ghosts, ref_vec );
        gvec[dir] = 0;
  
        bx_to = new Box( grow( bx_to_no_ghosts, gvec ) );
        fab_to = new FArrayBox( surroundingNodes(*bx_to,dim), a_dst.nComp() );
  
        fab_to->copy( *fab_from );
  
        delete bx_from;
        delete fab_from;
  
      }
  
    }
  
    CH_assert( fab_to != nullptr );
    CH_assert( bx_to != nullptr );
    CH_assert( (*bx_to) == a_dst_bx );
    CH_assert( n_to == n_dst );
  
    a_dst[dim].copy( *fab_to );
    delete fab_to;
    delete bx_to;

  }

  return;
}

void SGInterp::interpolate( EdgeDataBox&        a_dst,
                            const Box&          a_dst_bx,
                            const EdgeDataBox&  a_src,
                            const Box&          a_src_int_bx )
{
  CH_assert(a_dst.nComp() == a_src.nComp());

  const IntVect n_src = a_src_int_bx.bigEnd() 
                        - a_src_int_bx.smallEnd() 
                        + IntVect::Unit;
  const IntVect n_dst = a_dst_bx.bigEnd() 
                        - a_dst_bx.smallEnd() 
                        + IntVect::Unit;

  for (int dim = 0; dim < SpaceDim; dim++) {

    FArrayBox* fab_from(nullptr);
    FArrayBox* fab_to(nullptr);
    Box* bx_from(nullptr);
    Box* bx_to(nullptr);
    IntVect n_from(IntVect::Zero), n_to(IntVect::Zero);
  
    const Box& src_bx_wg = a_src.box();
    CH_assert( src_bx_wg.contains(grow(a_src_int_bx,num_gpt)) );
  
    bx_to = new Box( src_bx_wg );
    fab_to = new FArrayBox( enclosedCells(surroundingNodes(*bx_to),dim), 
                            a_src.nComp() );
    fab_to->copy(a_src[dim]);
    n_to = n_src;
  
    IntVect gvec = num_gpt*IntVect::Unit;
    IntVect cvec = IntVect::Unit; cvec[dim] = 0;
  
    for (int dir = 0; dir < SpaceDim; dir++) {
  
      if (n_src[dir] < n_dst[dir]) {
        
        /* refine */
        n_from = n_to;
        bx_from = bx_to;
        fab_from = fab_to;
  
        n_to[dir] = n_dst[dir];
  
        IntVect ref_vec(IntVect::Unit);
        ref_vec[dir] = (int) n_to[dir]/n_from[dir];
  
        Box bx_from_no_ghosts = grow( *bx_from, -1*gvec );
        Box bx_to_no_ghosts = refine( bx_from_no_ghosts, ref_vec );
        gvec[dir] = 0;
  
        bx_to = new Box( grow( bx_to_no_ghosts, gvec ) );
        fab_to = new FArrayBox( enclosedCells(surroundingNodes(*bx_to),dim), 
                                a_dst.nComp() );
 
        if (dim == dir ) {
          refine1DCC( *fab_to, 
                      enclosedCells(surroundingNodes(*bx_to),dim), 
                      *fab_from, 
                      enclosedCells(
                        surroundingNodes(grow(bx_from_no_ghosts,gvec)),dim), 
                      dir);
        } else {
          refine1DFC( *fab_to, 
                      enclosedCells(surroundingNodes(*bx_to),dim), 
                      *fab_from, 
                      enclosedCells (
                        surroundingNodes(grow(bx_from_no_ghosts,gvec)),dim), 
                      dir);
        }
  
        delete bx_from;
        delete fab_from;
  
      } else if (n_src[dir] > n_dst[dir]) {
  
        /* coarsen */
        n_from = n_to;
        bx_from = bx_to;
        fab_from = fab_to;
  
        n_to[dir] = n_dst[dir];
  
        IntVect ref_vec(IntVect::Unit);
        ref_vec[dir] = (int) n_from[dir]/n_to[dir];
  
        Box bx_from_no_ghosts = grow( *bx_from, -1*gvec );
        Box bx_to_no_ghosts = refine( bx_from_no_ghosts, ref_vec );
        gvec[dir] = 0;
  
        bx_to = new Box( grow( bx_to_no_ghosts, gvec ) );
        fab_to = new FArrayBox( enclosedCells(surroundingNodes(*bx_to),dim),
                                a_dst.nComp() );
  
        if (dim == dir) {
          coarsen1DCC(  *fab_to, 
                        enclosedCells(surroundingNodes(*bx_to),dim), 
                        *fab_from, 
                        enclosedCells (
                          surroundingNodes(grow(bx_from_no_ghosts,gvec)),dim),
                        dir );
        } else {
          coarsen1DFC(  *fab_to, 
                        enclosedCells(surroundingNodes(*bx_to),dim), 
                        *fab_from, 
                        enclosedCells (
                          surroundingNodes(grow(bx_from_no_ghosts,gvec)),dim),
                        dir );
        }
  
        delete bx_from;
        delete fab_from;
  
      } else {
  
        /* copy */
        n_from = n_to;
        bx_from = bx_to;
        fab_from = fab_to;
  
        n_to[dir] = n_dst[dir];
  
        IntVect ref_vec(IntVect::Unit);
        ref_vec[dir] = (int) n_to[dir]/n_from[dir];
  
        Box bx_from_no_ghosts = grow( *bx_from, -1*gvec );
        Box bx_to_no_ghosts = refine( bx_from_no_ghosts, ref_vec );
        gvec[dir] = 0;
  
        bx_to = new Box( grow( bx_to_no_ghosts, gvec ) );
        fab_to = new FArrayBox( enclosedCells(surroundingNodes(*bx_to),dim), 
                                a_dst.nComp() );
  
        fab_to->copy( *fab_from );
  
        delete bx_from;
        delete fab_from;
  
      }
  
    }
  
    CH_assert( fab_to != nullptr );
    CH_assert( bx_to != nullptr );
    CH_assert( (*bx_to) == a_dst_bx );
    CH_assert( n_to == n_dst );
  
    a_dst[dim].copy( *fab_to );
    delete fab_to;
    delete bx_to;

  }

  return;
}

#include "NamespaceFooter.H"

