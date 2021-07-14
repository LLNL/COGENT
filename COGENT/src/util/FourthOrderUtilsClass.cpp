
#include "FourthOrderUtilsClass.H"
#include "FourthOrderUtil.H"
#include "BoxIterator.H"

#if CH_SPACEDIM==2
#include "FourthOrderUtilF_F.H.2D"
#endif
#if CH_SPACEDIM==3
#include "FourthOrderUtilF_F.H.3D"
#endif
#if CH_SPACEDIM==3
#include "FourthOrderUtilF_F.H.3D"
#endif
#if CH_SPACEDIM==4
#include "FourthOrderUtilF_F.H.4D"
#endif
#if CH_SPACEDIM==5
#include "FourthOrderUtilF_F.H.5D"
#endif

#include "NamespaceHeader.H"

//********************************************************************************//
// Wrappers for calling the Chombo routines if sparse grids are 
// turned off and calling the new SG-friendly routines if they're
// turned on.
//********************************************************************************//
void FourthOrderUtil::fourthOrderAverageGen( LevelData<FArrayBox>& a_phi, 
                                             int a_sgn)
{
   
   if (!useSG) { fourthOrderAverage(a_phi, a_sgn); }
   else
   {
      FourthOrderUtil::fourthOrderAverageSG(a_phi, a_sgn);
   }

}
//********************************************************************************//

void FourthOrderUtil::fourthOrderAverageGen( LevelData<FluxBox>& a_phi, 
                                          int a_sgn)
{
   
   if (!useSG) { fourthOrderAverage(a_phi, a_sgn); }
   else
   {
      FourthOrderUtil::fourthOrderAverageSG(a_phi, a_sgn);
   }

}
//********************************************************************************//
// This routine is already sparse friendly (the operation only happens in one 
// direction, so just call the old routine regardless of whether SG is turned on)
void FourthOrderUtil::fourthOrderCellToFace( LevelData<FluxBox>& a_faceData, 
                                             const LevelData<FArrayBox>& a_cellData)
{
   
   fourthOrderCellToFace(a_faceData, a_cellData);

}
//********************************************************************************//
// This routine is also already sparse friendly, so just call the original
void FourthOrderUtil::fourthOrderCellToFaceCenters( LevelData<FluxBox>& a_faceData, 
                                                    const LevelData<FArrayBox>& a_cellData)
{
   
   fourthOrderCellToFaceCenters(a_faceData, a_cellData);

}
//********************************************************************************//
void FourthOrderUtil::fourthOrderMult( LevelData<FluxBox>& a_uTimesV, 
                                       const LevelData<FluxBox>& a_u, 
                                       const LevelData<FluxBox>& a_v,
                                       int a_sgn)
{
   if (!useSG) { fourthOrderMult(a_uTimesV, a_u, a_v, a_sgn); }
   else
   {
      FourthOrderUtil::fourthOrderMultSG(a_uTimesV, a_u, a_v, a_sgn);
   }
}
//********************************************************************************//
void FourthOrderUtil::fourthOrderCellProd( LevelData<FArrayBox>& a_uTimesV, 
                                           const LevelData<FArrayBox>& a_u, 
                                           const LevelData<FArrayBox>& a_v)
{
   if (!useSG) { fourthOrderCellProd( a_uTimesV, a_u, a_v ); }
   else
   {
      FourthOrderUtil::fourthOrderCellProdSG( a_uTimesV, a_u, a_v ); 
   }
}
//********************************************************************************//
//********************************************************************************//
// Highest level new SG-friendly routines
//********************************************************************************//
//********************************************************************************//

void FourthOrderUtil::fourthOrderAverageSG( LevelData<FArrayBox>& a_phi, 
                                            int a_sgn)
{
   CH_TIME("fourthOrderAverageSG");

   DataIterator dit = a_phi.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FArrayBox& thisPhi = a_phi[dit];

      fourthOrderAverageCellSG(thisPhi, a_sgn);
   }
}

//********************************************************************************//
void FourthOrderUtil::fourthOrderAverageSG( LevelData<FluxBox>& a_phi, 
                                            int a_sgn)
{
   CH_TIME("fourthOrderAverageSG_Flux");

   DataIterator dit = a_phi.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FluxBox& thisPhi = a_phi[dit];

      fourthOrderAverageFaceSG(thisPhi, a_sgn);
   }
}

//********************************************************************************//
void FourthOrderUtil::fourthOrderMultSG( LevelData<FluxBox>& a_uTimesV, 
                                         const LevelData<FluxBox>& a_u, 
                                         const LevelData<FluxBox>& a_v, 
                                         int a_sgn)
{
   int ncomp = a_u.nComp();
   CH_assert(a_v.nComp() == ncomp);
   CH_assert(a_uTimesV.nComp() == ncomp);

   DataIterator dit = a_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FourthOrderUtil::fourthOrderMultFaceSG(a_uTimesV[dit], a_u[dit], a_v[dit], a_sgn);
   }
}
//********************************************************************************//
void FourthOrderUtil::fourthOrderCellProdSG( LevelData<FArrayBox>& a_uTimesV, 
                                             const LevelData<FArrayBox>& a_u, 
                                             const LevelData<FArrayBox>& a_v)
{
   int ncomp = a_u.nComp();
   CH_assert(a_v.nComp() == ncomp);
   CH_assert(a_uTimesV.nComp() == ncomp);

   DataIterator dit = a_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FourthOrderUtil::fourthOrderCellProdSG(a_uTimesV[dit], a_u[dit], a_v[dit]);
   }
}
//********************************************************************************//
//********************************************************************************//
// One level lower new SG-friendly routines
//********************************************************************************//
//********************************************************************************//


void FourthOrderUtil::fourthOrderAverageCellGen( FArrayBox& a_phi, 
								   				 const ProblemDomain& a_domain, 
								   				 const Box& 			a_bx, 
                                          		 int a_sgn)
{
   
   if (!useSG) { fourthOrderAverageCell(a_phi, a_domain, a_bx, a_sgn); }
   else
   {
      FourthOrderUtil::fourthOrderAverageCellSG(a_phi, a_domain, a_bx, a_sgn);
   }
}


void FourthOrderUtil::fourthOrderAverageCellSG( FArrayBox& a_phi,
                                            int a_sgn )
{
   //derivBox shrunk by 1
   Box derivBox(a_phi.box());
   derivBox.grow(-1);

   FArrayBox tempLap(derivBox, a_phi.nComp());

   Real factor = a_sgn * (1.0/24.0);
   for (int dir=0; dir < SpaceDim; dir++) 
   {
      tempLap.setVal(0.0);
      FORT_INCREMENTLAPLACIAN(CHF_FRA(tempLap), 
                              CHF_CONST_FRA(a_phi), 
                              CHF_BOX(derivBox), 
                              CHF_CONST_INT(dir), 
                              CHF_CONST_REAL(factor));
      a_phi.plus(tempLap, derivBox, 0, 0, tempLap.nComp());
   }
}

//********************************************************************************//
void FourthOrderUtil::fourthOrderAverageCellSG( FArrayBox& a_phi,
								 				const ProblemDomain& a_domain,
								 				const Box& 		  a_bx,
								 				int a_sgn)
{
   CH_TIME("fourthOrderAverageCellSG domain");
   
   CH_assert(a_phi.box().contains(a_bx));
   Box derivBox(a_bx);
   derivBox.grow(1);
   derivBox &= a_domain;
   CH_assert(a_phi.box().contains(derivBox));

   FArrayBox tempLap(a_bx, a_phi.nComp());

   Real factor = a_sgn * (1.0/24.0);
   for (int dir=0; dir <SpaceDim; dir++)
   {
      tempLap.setVal(0.0);
      Box derivDirBox(a_bx);
	  if (!a_domain.isPeriodic(dir) &&
		  (a_domain.domainBox().smallEnd(dir) == a_bx.smallEnd(dir)))
	  {
	     // one-sided derivative, low end
	     Box loBox(a_bx);
		 loBox.setRange(dir, a_bx.smallEnd(dir));
		 FORT_INCREMENTLOSIDELAPLACIAN(CHF_FRA(tempLap),
									   CHF_FRA(a_phi),
									   CHF_BOX(loBox),
									   CHF_INT(dir),
									   CHF_REAL(factor));
		 derivDirBox.growLo(dir,-1);
	  }

	  if (!a_domain.isPeriodic(dir) &&
		  (a_domain.domainBox().bigEnd(dir) == a_bx.bigEnd(dir)))
	  {
	     // one-sided derivative, low end
	     Box hiBox(a_bx);
		 hiBox.setRange(dir, a_bx.bigEnd(dir));
		 FORT_INCREMENTHISIDELAPLACIAN(CHF_FRA(tempLap),
									   CHF_FRA(a_phi),
									   CHF_BOX(hiBox),
									   CHF_INT(dir),
									   CHF_REAL(factor));
		 derivDirBox.growHi(dir,-1);
	  }

	  FORT_INCREMENTHISIDELAPLACIAN(CHF_FRA(tempLap),
									CHF_FRA(a_phi),
									CHF_BOX(derivDirBox),
									CHF_INT(dir),
									CHF_REAL(factor));

	  a_phi.plus(tempLap, a_bx, 0, 0, tempLap.nComp());
	  

   }
}
//********************************************************************************//
void FourthOrderUtil::fourthOrderAverageFaceSG( FluxBox& a_phi, 
                                                int a_sgn )
{
   Real factor = a_sgn * (1.0/24.0);
   for (int dir=0; dir<SpaceDim; dir++)
   {
      FArrayBox& thisPhiDir = a_phi[dir];
      // derivBox shrunk by 1 in tangential directions
      Box derivBox(thisPhiDir.box());
      derivBox.grow(-1);
      derivBox.grow(dir,1);
      
      FArrayBox tempLap(derivBox, thisPhiDir.nComp());

      for (int tanDir=0; tanDir < SpaceDim; tanDir++)
      {
         if (tanDir != dir)
         {
            tempLap.setVal(0.0);
            FORT_INCREMENTLAPLACIAN(CHF_FRA(tempLap), 
                                    CHF_CONST_FRA(thisPhiDir),
                                    CHF_BOX(derivBox),
                                    CHF_CONST_INT(tanDir), 
                                    CHF_CONST_REAL(factor));
            thisPhiDir.plus(tempLap, derivBox, 0, 0, tempLap.nComp());
         }
      }
   }
}

//********************************************************************************//
void FourthOrderUtil::fourthOrderAverageFaceSG( FArrayBox& a_phi, 
												int a_dir, 
                                                int a_sgn )
{
   Real factor = a_sgn * (1.0/24.0);
   Box derivBox(a_phi.box());
   derivBox.grow(-1);
   derivBox.grow(a_dir,1);

   FArrayBox tempLap(derivBox, a_phi.nComp());

   for (int tanDir=0; tanDir < SpaceDim; tanDir++) {
      if (tanDir != a_dir) {
         tempLap.setVal(0.0);
         FORT_INCREMENTLAPLACIAN(CHF_FRA(tempLap),
  								 CHF_CONST_FRA(a_phi),
                				 CHF_BOX(derivBox), 
								 CHF_CONST_INT(tanDir), 
								 CHF_CONST_REAL(factor));
         a_phi.plus(tempLap, derivBox, 0, 0, tempLap.nComp());
      }
   }
}
//********************************************************************************//

//********************************************************************************//
// Already SG friendly, so just call original
void FourthOrderUtil::fourthOrderCellToFace( FluxBox& a_faceData, 
                                             const FArrayBox& a_cellData )
{
   fourthOrderCellToFace(a_faceData, a_cellData);
}

//********************************************************************************//
// Already SG friendly, so just call original
void FourthOrderUtil::fourthOrderCellToFaceCenters( FluxBox& a_faceData, 
                                                    const FArrayBox& a_cellData )
{
   fourthOrderCellToFaceCenters(a_faceData, a_cellData);
}
//********************************************************************************//
void FourthOrderUtil::fourthOrderMultFaceSG( FluxBox& a_uTimesV, 
                                             const FluxBox& a_u, 
                                             const FluxBox& a_v,
                                             int a_sgn)
{
   CH_TIME("computeCompFaceFluxes_SGversion");
   int ncomp = a_u.nComp();
   CH_assert(a_v.nComp() == ncomp);
   CH_assert(a_uTimesV.nComp() == ncomp);

   int deavg_sign = -a_sgn;
   int reavg_sign = a_sgn;
   
   // loop over directions
   for (int faceDir=0; faceDir<SpaceDim; faceDir++)
   {
      // pluck out the components for the present direction
      FArrayBox& thisUVdir = a_uTimesV[faceDir];
      const FArrayBox& thisUdir = a_u[faceDir];
      const FArrayBox& thisVdir = a_v[faceDir];

      // fetch the appropriate box
      Box intersectBox(thisUdir.box());
      intersectBox &= thisVdir.box();
      intersectBox &= thisUVdir.box();

      // temporary storage for point-values of u and v
      FArrayBox u_ptvals(intersectBox, thisUdir.nComp());
      FArrayBox v_ptvals(intersectBox, thisVdir.nComp());
      u_ptvals.copy(thisUdir, intersectBox);
      v_ptvals.copy(thisVdir, intersectBox);

      // compute point values of u and v
      FourthOrderUtil::fourthOrderAverageFaceSG( u_ptvals, faceDir, deavg_sign );
      FourthOrderUtil::fourthOrderAverageFaceSG( v_ptvals, faceDir, deavg_sign );

      // multiply point values
      thisUVdir.copy(u_ptvals, intersectBox);
      thisUVdir.mult(v_ptvals, intersectBox, 0, 0, ncomp);

      // reaverage to get <uv>
      FourthOrderUtil::fourthOrderAverageFaceSG( thisUVdir, faceDir, reavg_sign );

   }

}
//********************************************************************************//
void FourthOrderUtil::fourthOrderCellProdSG( FArrayBox& a_uTimesV, 
                                             const FArrayBox& a_u, 
                                             const FArrayBox& a_v)
{
   int ncomp = a_u.nComp();
   CH_assert(a_v.nComp() == ncomp);
   CH_assert(a_uTimesV.nComp() == ncomp);

   int deavg_sign = -1;
   int reavg_sign = 1;

   Box intersectBox(a_u.box());
   intersectBox &= a_v.box();
   intersectBox &= a_uTimesV.box();

   FArrayBox u_ptvals(intersectBox, ncomp);
   FArrayBox v_ptvals(intersectBox, ncomp);
   u_ptvals.copy(a_u, intersectBox);
   v_ptvals.copy(a_v, intersectBox);

   FourthOrderUtil::fourthOrderAverageCellSG( u_ptvals, deavg_sign );
   FourthOrderUtil::fourthOrderAverageCellSG( v_ptvals, deavg_sign );

   a_uTimesV.copy( u_ptvals, intersectBox );
   a_uTimesV.mult( v_ptvals, intersectBox, 0, 0, ncomp );

   FourthOrderUtil::fourthOrderAverageCellSG( a_uTimesV, reavg_sign );
}


#include "NamespaceFooter.H"
