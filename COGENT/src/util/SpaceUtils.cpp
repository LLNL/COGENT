#include "SpaceUtils.H"
#include "SpaceUtilsF_F.H"

#include "NamespaceHeader.H"

void
SpaceUtils::faceInterpolate(  const int         a_dir,
                              const Box&        a_fc_box,
                              const int         a_order,
                              const FArrayBox&  a_var,
                              FArrayBox&        a_face_var )
{
  int ncomp = a_var.nComp();
  CH_assert(a_face_var.nComp() == ncomp);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert (a_dir >= 0 && a_dir < SpaceDim);

  for (int n=0; n<ncomp; n++) {
    FORT_FACE_INTERPOLATE(  CHF_CONST_INT(a_dir),
                            CHF_BOX(a_fc_box),
                            CHF_CONST_INT(a_order),
                            CHF_CONST_FRA1(a_var, n),
                            CHF_FRA1(a_face_var, n) );
  }

  return;
}


void
SpaceUtils::faceInterpolate(  const int         a_dir,
                              const Box&        a_fc_box,
                              const Box&        a_cc_box,
                              const int         a_order,
                              const FArrayBox&  a_var,
                              FArrayBox&        a_face_var )
{
  int ncomp = a_var.nComp();
  CH_assert(a_face_var.nComp() == ncomp);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert (a_dir >= 0 && a_dir < SpaceDim);

  int imin = a_cc_box.smallEnd(a_dir),
      imax = a_cc_box.bigEnd(a_dir);

  if (imin == imax) {

    /* box is flat in along this dimension,
     * so copy instead of interpolate */
    BoxIterator bit(a_cc_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      IntVect iw = iv; iw[a_dir]++;

      for (int n=0; n<ncomp; n++) {
        a_face_var(iv,n) = a_var(iv,n);
        a_face_var(iw,n) = a_var(iv,n);
      }
    }

  } else {

    faceInterpolate(a_dir, a_fc_box, a_order, a_var, a_face_var);

  }

  return;
}


void
SpaceUtils::cellCenteredGradientComponent(  const Box&        a_box,
                                            const int         a_dir,
                                            const FArrayBox&  a_var,
                                            const RealVect&   a_dx,
                                            const int         a_order,
                                            FArrayBox&        a_grad_var )
{
  CH_assert(a_var.nComp() == 1);
  CH_assert(a_grad_var.nComp() > a_dir);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  FORT_CELL_CENTERED_GRAD_COMPONENT(  CHF_BOX(a_box),
                                      CHF_CONST_INT(a_dir),
                                      CHF_CONST_FRA1(a_var,0),
                                      CHF_CONST_REALVECT(a_dx),
                                      CHF_CONST_INT(a_order),
                                      CHF_CONST_FRA1(a_grad_var, a_dir) );

  return;
}

void
SpaceUtils::faceCenteredGradientComponent(  const Box&        a_box,
                                            const int         a_dir,
                                            const FArrayBox&  a_var,
                                            const RealVect&   a_dx,
                                            const int         a_order,
                                            FArrayBox&        a_grad_var )
{
  CH_assert(a_var.nComp() == 1);
  CH_assert(a_grad_var.nComp() > a_dir);
  CH_assert(a_order == 2 || a_order == 4);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  FORT_FACE_CENTERED_GRAD_COMPONENT(  CHF_BOX(a_box),
                                      CHF_CONST_INT(a_dir),
                                      CHF_CONST_FRA1(a_var, 0),
                                      CHF_CONST_REALVECT(a_dx),
                                      CHF_CONST_INT(a_order),
                                      CHF_FRA1(a_grad_var, a_dir) );

  return;
}

void
SpaceUtils::extrapBoundaryGhostsForCC( FArrayBox&                              a_data,
                                       const Box&                              a_interiorbox,
                                       const Box&                              a_domain_box,
                                       const int                               a_dir,
                                       const int                               a_order,
                                       const Tuple<BlockBoundary, 2*SpaceDim>& a_block_boundaries)
{
   // If a_order = 2, this function second-order extrapolates to fill two layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   // If a_order = 4, this function fourth-order extrapolates to fill three layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));

   int imin = a_interiorbox.smallEnd(a_dir);
   int imax = a_interiorbox.bigEnd(a_dir);

      int depth = (a_order==2)? 2: 3;
   
   for (int side=-1; side<2; side+=2) {

      if ( a_block_boundaries[a_dir + (side+1)*SpaceDim/2].isDomainBoundary() ) {

         Box dstbox;
         bool isbc = false;
      
         switch(side) 
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(a_dir) == a_domain_box.smallEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellLo(a_interiorbox, a_dir, depth);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(a_dir) == a_domain_box.bigEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellHi(a_interiorbox, a_dir, depth);
               }
               break;
            }

         if (isbc) {

            CH_assert(bx.contains(dstbox));

            if (imin == imax) {

               BoxIterator bit(dstbox);
               for (bit.begin(); bit.ok(); ++bit) {
                 IntVect i1 = bit();
                 IntVect i2 = i1; i2[a_dir] = imin;
                 for (int n=0; n<a_data.nComp(); n++) {
                    a_data(i1,n) = a_data(i2,n);
                 }
               }

            } else {

               CH_assert( (a_order==4 && a_interiorbox.size(a_dir)>=5)
                       || (a_order==2 && a_interiorbox.size(a_dir)>=3));
               FORT_EXTRAP_FOR_CC_OPS(CHF_CONST_INT(a_dir),
                                      CHF_CONST_INT(side),
                                      CHF_CONST_INT(a_order),
                                      CHF_BOX(dstbox),
                                      CHF_BOX(a_interiorbox),
                                      CHF_FRA(a_data));

            }
         }
      }
   }

   return;
} 

void
SpaceUtils::extrapBoundaryGhostsForFC( FArrayBox&                              a_data,
                                       const Box&                              a_interiorbox,
                                       const Box&                              a_domain_box,
                                       const int                               a_dir,
                                       const int                               a_order,
                                       const Tuple<BlockBoundary, 2*SpaceDim>& a_block_boundaries)
{
   // This function fourth-order extrapolates to fill two layers of a_data ghost cells
   // at physical boundaries in the direction a_dir.

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));

   int imin = a_interiorbox.smallEnd(a_dir);
   int imax = a_interiorbox.bigEnd(a_dir);

   for (int side=-1; side<2; side+=2) {

      if ( a_block_boundaries[a_dir + (side+1)*SpaceDim/2].isDomainBoundary() ) {

         Box dstbox;
         bool isbc = false;
      
         switch(side) 
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(a_dir) == a_domain_box.smallEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellLo(a_interiorbox, a_dir, 2);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(a_dir) == a_domain_box.bigEnd(a_dir));
               if (isbc) {
                  dstbox = adjCellHi(a_interiorbox, a_dir, 2);
               }
               break;
            }

         if (isbc) {

            CH_assert(bx.contains(dstbox));

            if (imin == imax) {

               BoxIterator bit(dstbox);
               for (bit.begin(); bit.ok(); ++bit) {
                 IntVect i1 = bit();
                 IntVect i2 = i1; i2[a_dir] = imin;
                 for (int n=0; n<a_data.nComp(); n++) {
                    a_data(i1,n) = a_data(i2,n);
                 }
               }

            } else {

               CH_assert( (a_order==4 && a_interiorbox.size(a_dir)>=5)
                       || (a_order==2 && a_interiorbox.size(a_dir)>=3));
               FORT_EXTRAP_FOR_FC_OPS(CHF_CONST_INT(a_dir),
                                      CHF_CONST_INT(side),
                                      CHF_CONST_INT(a_order),
                                      CHF_BOX(dstbox),
                                      CHF_BOX(a_interiorbox),
                                      CHF_FRA(a_data));

            }
         }
      }
   }

   return;
} 


void
SpaceUtils::secondOrderTransExtrapAtDomainBdry(FArrayBox&           a_data,
                                               const int            a_dir,
                                               const Box&           a_interiorbox,
                                               const ProblemDomain& a_domain,
                                               const int            a_maxdim )
{
   CH_TIME("secondOrderTransExtrapAtDomainBdry");

   const Box& bx= a_data.box();
   CH_assert(bx.contains(a_interiorbox));
   const Box& domBox = a_domain.domainBox();

   for (int tdir=0; tdir<a_maxdim; ++tdir)
   {
      if (tdir != a_dir && !a_domain.isPeriodic(tdir))
      {
         CH_assert(a_interiorbox.size(tdir)>=3);

         for (int side=-1; side<2; side+=2)
         {
            Box dstbox;
            bool isbc = false;

            switch(side)
            {
            case -1:
               // Handle the low side
               isbc = (a_interiorbox.smallEnd(tdir) == domBox.smallEnd(tdir));
               if (isbc)
               {
                  dstbox = adjCellLo(a_interiorbox,tdir,1);
               }
               break;
            case 1:
               // Handle high side
               isbc = (a_interiorbox.bigEnd(tdir) == domBox.bigEnd(tdir));
               if (isbc)
               {
                  dstbox = adjCellHi(a_interiorbox,tdir,1);
               }
               break;
            }

            if (isbc)
            {
               CH_assert(bx.contains(dstbox));
               FORT_SECOND_ORDER_EXTRAPOLATION( CHF_CONST_INT(tdir),
                                                CHF_CONST_INT(side),
                                                CHF_BOX(a_interiorbox),
                                                CHF_BOX(dstbox),
                                                CHF_FRA(a_data) );
            }
         }
      }
   }

   return;
}



#include "NamespaceFooter.H"
