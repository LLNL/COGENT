#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif



#include "mappedAdvectionFlux.H"
#include "mappedLimiter.H"
#include "FourthOrderUtil.H"
#include "BoxIterator.H"

#include "mappedAdvectionF_F.H"
#include "FourthOrderUtil.H"
#include "FourthOrderUtilsClass.H.multidim"

#include "NamespaceHeader.H"


void
cellUJToCellU(
              LevelData<FArrayBox>& a_cell_u,
              const LevelData<FArrayBox>& a_cell_uJ,
              const FourthOrderCoordSys* a_FOCS)
  {
   int ncomp = a_cell_u.nComp();
   CH_assert(a_cell_uJ.nComp() == ncomp);

   // Obtain J
   const LevelData<FArrayBox>& J = a_FOCS->getJ();

   DataIterator dit = a_cell_uJ.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FArrayBox& thisU = a_cell_u[dit];
      const FArrayBox& thisUJ = a_cell_uJ[dit];
      const FArrayBox& thisJ = J[dit];
      FArrayBox thisJinv(thisJ.box(), 1);
      thisJinv.copy(thisJ);
      thisJinv.invert(1.0);

      Box intersectBox(thisUJ.box());
      intersectBox &= thisJinv.box();
      intersectBox &= thisU.box();

      thisU.copy(thisUJ, intersectBox);
      BoxIterator ibit(intersectBox);
      for (ibit.begin(); ibit.ok(); ++ibit)
      {
         const IntVect& iv = ibit();
         const Real val = thisJinv(iv);
         for (int n=0; n<ncomp; n++)
            thisU(iv,n) *= val;
      }

      // compute box over which we can do this
      Box gradIntersectBox(thisU.box());
      gradIntersectBox.grow(-1);
      gradIntersectBox &= thisJ.box();
      FArrayBox gradProduct(gradIntersectBox, ncomp);
      gradProduct.setVal( 0.0 );

      // since dx's cancel, use 1 here
      Real fakeDx = 1.0;
      Real factor = -1.0/12.0;
      for (int dir=0; dir<SpaceDim; dir++)
        {
          // call fortran wrapper function to handle multidim case
          incrementGradProduct(gradProduct,
                               thisU,
                               thisJ,
                               gradIntersectBox,
                               fakeDx,
                               factor,
                               dir);

        }

      BoxIterator jbit(gradIntersectBox);
      for (jbit.begin(); jbit.ok(); ++jbit)
      {
         const IntVect& iv = jbit();
         const Real val = thisJinv(iv);
         for (int n=0; n<ncomp; n++)
            gradProduct(iv,n) *= val;
      }

      thisU.plus(gradProduct, gradIntersectBox, 0, 0, ncomp);
   } // end loop over boxes
}




void
computeCompFaceFluxes( LevelData<FluxBox>& a_uTimesV,
                       const LevelData<FluxBox>& a_u,
                       const LevelData<FluxBox>& a_v,
                       bool a_useFourthOrder)
{
  CH_TIME("computeCompFaceFluxes");

   // Compute the SpaceDim-by-nComp face-averaged fluxes in computational
   // space, where a_v is the SpaceDim-dimensional velocity vector and
   // a_u is the nComp-dim state vector
   int ncomp = a_u.nComp();
   CH_assert(a_v.nComp() == SpaceDim);
   CH_assert(a_uTimesV.nComp() == SpaceDim * ncomp);

   // loop over boxes
   DataIterator dit = a_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FluxBox& thisUV = a_uTimesV[dit];
      const FluxBox& thisU = a_u[dit];
      const FluxBox& thisV = a_v[dit];

      // loop over faces (index "d" in notes)
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
      {
         FArrayBox& thisUVdir = thisUV[faceDir];
         const FArrayBox& thisUdir = thisU[faceDir];
         const FArrayBox& thisVdir = thisV[faceDir];

         // compute <u_p><v_d> tensor
         Box intersectBox(thisUdir.box());
         intersectBox &= thisVdir.box();
         intersectBox &= thisUVdir.box();

         thisUVdir.setVal(0.0);
         FORT_INCREMENTFACEPROD(CHF_FRA(thisUVdir),
                                CHF_CONST_FRA(thisUdir),
                                CHF_CONST_FRA(thisVdir),
                                CHF_BOX(intersectBox));

         if (a_useFourthOrder)
         {
            // now increment with product of tangential gradients
            // (sum over index "d'" in notes)
            for (int tanDir=0; tanDir<SpaceDim; tanDir++)
            {
               if (tanDir!=faceDir)
               {
                  // since dx's cancel, use 1 here
                  Real fakeDx = 1.0;
                  Real factor = 1.0/12.0;

                  // compute box over which we can do this
                  Box gradIntersectBox(thisUdir.box());
                  gradIntersectBox &= thisVdir.box();
                  gradIntersectBox.grow(tanDir,-1);
                  gradIntersectBox &= thisUVdir.box();

                  FORT_INCREMENTFACEPRODGRAD(CHF_FRA(thisUVdir),
                                             CHF_CONST_FRA(thisUdir),
                                             CHF_CONST_FRA(thisVdir),
                                             CHF_BOX(gradIntersectBox),
                                             CHF_REAL(fakeDx),
                                             CHF_REAL(factor),
                                             CHF_INT(tanDir));

               } // end if tangential direction
            } // end loop over tanDir
         } // end if fourth order
      } // end loop over faceDir
   } // end loop over boxes
}

// The SG-friendly version
void
computeCompFaceFluxesSG( LevelData<FluxBox>& a_uTimesV,
                       const LevelData<FluxBox>& a_u,
                       const LevelData<FluxBox>& a_v,
                       bool a_useFourthOrder)
{
  CH_TIME("computeCompFaceFluxes SG version");

   // Compute the SpaceDim-by-nComp face-averaged fluxes in computational
   // space, where a_v is the SpaceDim-dimensional velocity vector and
   // a_u is the nComp-dim state vector
   int ncomp = a_u.nComp();
   CH_assert(a_v.nComp() == SpaceDim);
   CH_assert(a_uTimesV.nComp() == SpaceDim * ncomp);
   CH_assert(a_uTimesV.ghostVect() == IntVect::Unit);
	
	FourthOrderUtil FourthOrderOperators; //Object that holds various fourth-order operations 
	FourthOrderOperators.setSG(true); //Whether to use the SG versions of fourth order stencils
	int deavg_sign = -1.0;
   int reavg_sign = 1.0;
  
	LevelData<FluxBox> u_ptvals( a_u.disjointBoxLayout(), ncomp, a_u.ghostVect() );
	LevelData<FluxBox> v_ptvals( a_v.disjointBoxLayout(), SpaceDim, a_v.ghostVect() );

	for (auto dit(a_u.dataIterator()); dit.ok(); ++dit)
	{
		for (int j=0; j < SpaceDim; ++j)
		{
			u_ptvals[dit][j].copy(a_u[dit][j]);
		}
	}
	
	for (auto dit(a_v.dataIterator()); dit.ok(); ++dit)
	{
		for (int j=0; j < SpaceDim; ++j)
		{
			v_ptvals[dit][j].copy(a_v[dit][j]);
		}
	}

	FourthOrderOperators.fourthOrderAverageGen(u_ptvals, deavg_sign);
	FourthOrderOperators.fourthOrderAverageGen(v_ptvals, deavg_sign);

	u_ptvals.exchange();
	v_ptvals.exchange();

   // loop over boxes
   DataIterator dit = u_ptvals.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FluxBox& thisUV = a_uTimesV[dit];
      const FluxBox& thisU = u_ptvals[dit];
      const FluxBox& thisV = v_ptvals[dit];

      // loop over faces (index "d" in notes)
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
      {
         FArrayBox& thisUVdir = thisUV[faceDir];
         const FArrayBox& thisUdir = thisU[faceDir];
         const FArrayBox& thisVdir = thisV[faceDir];

         // compute <u_p><v_d> tensor
         Box intersectBox(thisUdir.box());
         intersectBox &= thisVdir.box();
         intersectBox &= thisUVdir.box();

			thisUVdir.setVal(0.0);
         FORT_INCREMENTFACEPROD(CHF_FRA(thisUVdir),
                                CHF_CONST_FRA(thisUdir),
                                CHF_CONST_FRA(thisVdir),
                                CHF_BOX(intersectBox));

		}
	}
	

	FourthOrderOperators.fourthOrderAverageGen(a_uTimesV, reavg_sign);

}

void
computeCompFaceFluxesGen( LevelData<FluxBox>& a_uTimesV,
                       const LevelData<FluxBox>& a_u,
                       const LevelData<FluxBox>& a_v,
                       bool a_useFourthOrder,
							  bool a_useSG)
{
   if (a_useSG) { computeCompFaceFluxesSG(a_uTimesV, a_u, a_v, a_useFourthOrder); }
	else { computeCompFaceFluxes(a_uTimesV, a_u, a_v, a_useFourthOrder); }
}

void
computeCompFaceFluxesNormal( LevelData<FluxBox>& a_uTimesV,
                             const LevelData<FluxBox>& a_u,
                             const LevelData<FluxBox>& a_v,
                             bool a_useFourthOrder)
{
   // Compute the SpaceDim-by-nComp face-averaged fluxes in computational
   // space, where a_v is the SpaceDim-dimensional velocity vector and
   // a_u is the nComp-dim state vector

   // loop over boxes
   DataIterator dit = a_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FluxBox& thisUV = a_uTimesV[dit];
      const FluxBox& thisU = a_u[dit];
      const FluxBox& thisV = a_v[dit];

      // loop over faces (index "d" in notes)
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
      {
         FArrayBox& thisUVdir = thisUV[faceDir];
         const FArrayBox& thisUdir = thisU[faceDir];
         const FArrayBox& thisVdir = thisV[faceDir];

         // compute <u_p><v_d> tensor
         Box intersectBox(thisUdir.box());
         intersectBox &= thisVdir.box();
         intersectBox &= thisUVdir.box();

         thisUVdir.setVal(0.0);
         FORT_INCREMENTFACEPRODNORMAL(CHF_FRA(thisUVdir),
                                CHF_CONST_FRA(thisUdir),
                                CHF_CONST_FRA(thisVdir),
                                CHF_BOX(intersectBox));

         if (a_useFourthOrder)
         {
            // now increment with product of tangential gradients
            // (sum over index "d'" in notes)
            for (int tanDir=0; tanDir<SpaceDim; tanDir++)
            {
               if (tanDir!=faceDir)
               {
                  // since dx's cancel, use 1 here
                  Real fakeDx = 1.0;
                  Real factor = 1.0/12.0;

                  // compute box over which we can do this
                  Box gradIntersectBox(thisUdir.box());
                  gradIntersectBox &= thisVdir.box();
                  gradIntersectBox.grow(tanDir,-1);
                  gradIntersectBox &= thisUVdir.box();

                  FORT_INCREMENTFACEPRODGRADNORMAL(CHF_FRA(thisUVdir),
                                             CHF_CONST_FRA(thisUdir),
                                             CHF_CONST_FRA(thisVdir),
                                             CHF_BOX(gradIntersectBox),
                                             CHF_REAL(fakeDx),
                                             CHF_REAL(factor),
                                             CHF_INT(tanDir));

               } // end if tangential direction
            } // end loop over tanDir
         } // end if fourth order
      } // end loop over faceDir
   } // end loop over boxes
}


#include "NamespaceFooter.H"
