#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include <cstdio>

#include "BoxIterator.H"
#include "ParmParse.H"
#include "DebugDump.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "AMRIO.H"
#include "FluxBox.H"

#include <fstream>
#include <float.h>
#include "FieldData.H"

#include "UsingNamespace.H"

#define LARGEINTVAL -999999999
#define LARGEREALVAL 123456789.0

void oneSidedPureDerivs(Real               & a_firstDeriv,
                        Real               & a_secDeriv  , 
                        const Real         & a_evalPt    ,
                        const Vector<Real> & a_yVal      ,
                        const Vector<Real> & a_xVal      )
{
  // interpolate quadratic (Lagrange formula for evenly spaced points). Differentiate twice. 
  CH_assert (a_yVal.size() == a_xVal.size());
  
  // order of extrapolation
  int order = a_yVal.size();
  
  if (order == 0)
    {
      MayDay::Abort("a_points has size 0");
    }
  else
    {
      // elements of Lagrange formula for evenly spaced points
      Vector<Real>prodTerm(order);
      Vector<Real>sumTerm (order);

      // set to 1.0 for ease of applying a product formula
      for(int j = 0; j < order; ++ j)
        {
          prodTerm[j] = 1.0;
        }
      
      // set to 0.0 for summing
      for(int j = 0; j < order; ++ j)
        {
          sumTerm[j] = 0.0;
        }

      for (int j = 0; j < order; ++j)
        {
          for(int m = 0; m < order; ++ m)
            {
              if (j != m)
                {
                  prodTerm[j] *= (a_xVal[j] - a_xVal[m]); 
                  sumTerm [j] += a_xVal[m];
                }
            }
        }

      // elements of Lagrange formula for evenly spaced points
      Vector<Real>firstDerivSummand(order);
      Vector<Real>secDerivSummand  (order);
      for (int j = 0; j < order; ++j)
        {
          firstDerivSummand[j] = (2.0*a_evalPt - sumTerm[j])/prodTerm[j]; 
          secDerivSummand  [j] = (2.0                      )/prodTerm[j]; 
        }
      
      // return values
      a_firstDeriv = 0.0;
      a_secDeriv   = 0.0;

      for (int j = 0; j < order; ++j)
        {
          a_firstDeriv += a_yVal[j]*firstDerivSummand[j];
          a_secDeriv   += a_yVal[j]*secDerivSummand  [j];
        }
    }
}
void firstDerivCenteredDiff(Real                    & a_xDeriv         ,
                            Real                    & a_yDeriv         ,
                            const IntVect      & a_iv                 ,
                            const int              & a_dir               ,
                            const RealVect    & a_dx               ,
                            const FArrayBox & a_curr            ,
                            const bool            & a_writeInfo    )
{
  // only called for non-cornerCells
  Real retvalX = LARGEREALVAL;
  Real retvalY = LARGEREALVAL;
  
  int xComp = 0;
  int yComp = 1;

  // stencil ivs for centered difference
  IntVect hi = a_iv;
  IntVect lo = a_iv;
  
  hi += BASISV(a_dir);
  lo -= BASISV(a_dir);
 
  // for differences with respect to the a_dir direction
  Real dx = a_dx[a_dir];

  // derivative of X-component in the a_dir direction
  retvalX = (a_curr(hi,xComp) - a_curr(lo,xComp))/(2.0*dx);
   
  //retval
  a_xDeriv = retvalX;
  
  // derivative of Y-component in the a_dir direction
  retvalY = (a_curr(hi,yComp) - a_curr(lo,yComp))/(2.0*dx);
    
  //retval
  a_yDeriv = retvalY;
}

void calculateBdDerivIv(FArrayBox       & a_mapPartials,
                        const IntVect   & a_iv         ,
                        const FArrayBox & a_curr       ,
                        const int       & a_oneSideDir ,
                        const IntVect   & a_hiLoSide   ,
                        const RealVect  & a_dx         ,
                        const bool      & a_cornerIv   ,
                        const bool      & a_BFieldCase )
{
  Real deriv1 = LARGEREALVAL;
  Real deriv2 = LARGEREALVAL;

  Real evalPt = 0.0;

  Vector<Real> xVal(3);
  Vector<Real> yVal(3);

  for (int xyDir = 0; xyDir < SpaceDim; ++ xyDir)
    {
      for(int xiEtaDir = 0; xiEtaDir< SpaceDim; ++ xiEtaDir)
        {
          if((a_oneSideDir == xiEtaDir) || (a_cornerIv))
            {
              // x-values for quadratic interpolation
              xVal[0] = evalPt;
              xVal[1] = evalPt - a_hiLoSide[xiEtaDir]*    a_dx[xiEtaDir];
              xVal[2] = evalPt - a_hiLoSide[xiEtaDir]*2.0*a_dx[xiEtaDir];
              
              // y-values for quadratic interpolation
              yVal[0] = a_curr(a_iv,xyDir);
              yVal[1] = a_curr(a_iv - a_hiLoSide[xiEtaDir]*     BASISV(xiEtaDir),xyDir);
              yVal[2] = a_curr(a_iv - a_hiLoSide[xiEtaDir]* 2 * BASISV(xiEtaDir),xyDir);
              
              // differentiate the quadratic
              oneSidedPureDerivs(deriv1,
                                 deriv2, 
                                 evalPt,
                                 yVal  ,
                                 xVal  );
            }
          else
            {
              // use centered differencing
              deriv1 = (a_curr(a_iv + BASISV(xiEtaDir),xyDir) -                          a_curr(a_iv - BASISV(xiEtaDir),xyDir))/(2.0           *a_dx[xiEtaDir]);
              deriv2 = (a_curr(a_iv + BASISV(xiEtaDir),xyDir) - 2.0*a_curr(a_iv,xyDir) + a_curr(a_iv - BASISV(xiEtaDir),xyDir))/(a_dx[xiEtaDir]*a_dx[xiEtaDir]);
            }

          // assign first derivatives
          a_mapPartials(a_iv,xyDir +        xiEtaDir *SpaceDim) = deriv1;
          
          // this function is reused from mapPartials in the case of setting the perpField.
          if (!a_BFieldCase)
            {
              // assign second derivatives
              a_mapPartials(a_iv,xyDir +(2 + xiEtaDir)*SpaceDim) = deriv2;
            }
        }    
    }
}
void computeSecondDeriv( Real            &a_x2Xi      , 
                         Real            &a_y2Xi      ,   
                         Real            &a_x2Eta     ,  
                         Real            &a_y2Eta     ,
                         Real            &a_xXiEta    ,
                         Real            &a_yXiEta    ,
                         const IntVect   & a_iv       ,
                         const RealVect  & a_dx       ,
                         const FArrayBox & a_curr     ,
                         const bool      & a_writeInfo)      
{
  // indexing macros
  int xDir = 0;
  int yDir = 1;
  
  // stencil ivs for pure derivatives
  IntVect xHi( 1, 0);
  IntVect xLo(-1, 0);
  IntVect yHi( 0, 1);
  IntVect yLo( 0,-1);
  
  xHi += a_iv;
  xLo += a_iv;
  yHi += a_iv;
  yLo += a_iv;

  // stencil ivs for cross derivatives
  IntVect hiHi( 1, 1);
  IntVect loLo(-1,-1);
  IntVect hiLo( 1,-1);
  IntVect loHi(-1, 1);

  hiHi += a_iv;
  loLo += a_iv;
  hiLo += a_iv;
  loHi += a_iv;

  // dx
  Real dx   = a_dx[xDir];
  Real dy   = a_dx[yDir];
  Real dx2  = dx*dx;
  Real dy2  = dy*dy;
  Real dxdy = dx*dy;
  
  // second partials for R and Z
  a_x2Xi  = (a_curr(xHi,xDir) - 2*a_curr(a_iv,xDir) + a_curr(xLo,xDir))/dx2; 
  a_x2Eta = (a_curr(yHi,xDir) - 2*a_curr(a_iv,xDir) + a_curr(yLo,xDir))/dy2;
  
  a_y2Xi  = (a_curr(xHi,yDir) - 2*a_curr(a_iv,yDir) + a_curr(xLo,yDir))/dx2; 
  a_y2Eta = (a_curr(yHi,yDir) - 2*a_curr(a_iv,yDir) + a_curr(yLo,yDir))/dy2;
  
  a_xXiEta = (a_curr(hiHi,xDir) + a_curr(loLo,xDir) - a_curr(hiLo,xDir) - a_curr(loHi,xDir))/(4*dxdy);
  a_yXiEta = (a_curr(hiHi,yDir) + a_curr(loLo,yDir) - a_curr(hiLo,yDir) - a_curr(loHi,yDir))/(4*dxdy);
}

void computeMapPartialsBdIv(FArrayBox       & a_mapPartials, 
                            const FArrayBox & a_curr       ,
                            const Box       & a_validBox   ,
                            const RealVect  & a_dx         ,
                            const bool      & a_BFieldCase ) 
{
  // macros for components of the mapping
  int xDir = 0;
  int yDir = 1;

  Box domainBox = a_validBox;

  // distinguish different stencils for the boundary ivs
  bool cornerIv = true;
  IntVect hiLoSide;
   
  // not used for corner ivs
  int oneSideDir = LARGEINTVAL;
   
  // lower left corner
  IntVect loLo = domainBox.smallEnd();
  hiLoSide[xDir] = -1;
  hiLoSide[yDir] = -1; 
   
  // this variable reused for all boundary cases
  IntVect iv = LARGEREALVAL*IntVect::Unit;
   
  // IntVect
  iv = loLo;
   
  // corner iv
  cornerIv = true;
   
  // fill first derivatives and the pure secondDerivatives 
  calculateBdDerivIv(a_mapPartials,
                     iv           ,
                     a_curr       ,
                     oneSideDir   ,
                     hiLoSide     ,
                     a_dx         ,
                     cornerIv     ,
                     a_BFieldCase );
   
  // upper right corner
  IntVect hiHi = domainBox.bigEnd  ();
  hiLoSide[xDir] = 1;
  hiLoSide[yDir] = 1; 
   
  // IntVect
  iv = hiHi;
   
  // corner iv
  cornerIv = true;
   
  // fill first derivatives and the pure secondDerivatives 
  calculateBdDerivIv(a_mapPartials,
                     iv           ,
                     a_curr       ,
                     oneSideDir   ,
                     hiLoSide     ,
                     a_dx         ,
                     cornerIv     ,
                     a_BFieldCase );
   
   
  // upper left corner
  IntVect loHi(loLo[0],hiHi[1]);
  hiLoSide[xDir] = -1;
  hiLoSide[yDir] =  1; 
   
  // IntVect
  iv = loHi;
   
  // corner iv
  cornerIv = true;
   
  // fill first derivatives and the pure secondDerivatives 
  calculateBdDerivIv(a_mapPartials,
                     iv           ,
                     a_curr       ,
                     oneSideDir   ,
                     hiLoSide     ,
                     a_dx         ,
                     cornerIv     ,
                     a_BFieldCase );
   
  // lower right corner
  IntVect hiLo(hiHi[0],loLo[1]);
  hiLoSide[xDir] =  1;
  hiLoSide[yDir] = -1; 
   
  // IntVect
  iv = hiLo;
   
  // corner iv
  cornerIv = true;
   
  // fill first derivatives and the pure second derivatives 
  calculateBdDerivIv(a_mapPartials,
                     iv           ,
                     a_curr       ,
                     oneSideDir   ,
                     hiLoSide     ,
                     a_dx         ,
                     cornerIv     ,
                     a_BFieldCase );
   
  // set the corner mixed second derivatives arbitrarily if not the BField case
  if(!a_BFieldCase)
    {
      a_mapPartials(loLo,xDir + 4*SpaceDim) = 1.0;
      a_mapPartials(loLo,yDir + 4*SpaceDim) = 1.0;
      
      a_mapPartials(hiHi,xDir + 4*SpaceDim) = 1.0;
      a_mapPartials(hiHi,yDir + 4*SpaceDim) = 1.0;
      
      a_mapPartials(loHi,xDir + 4*SpaceDim) = 1.0;
      a_mapPartials(loHi,yDir + 4*SpaceDim) = 1.0;
      
      a_mapPartials(hiLo,xDir + 4*SpaceDim) = 1.0;
      a_mapPartials(hiLo,yDir + 4*SpaceDim) = 1.0;
    }

  // iteration boxes around the edges
  Box xEdgeLo(loLo,loHi);
  Box xEdgeHi(hiLo,hiHi);
   
  Box yEdgeLo(loLo,hiLo);
  Box yEdgeHi(loHi,hiHi);
   
  // iterate over the interior of the edges
  xEdgeLo.grow(-BASISV(yDir));
  xEdgeHi.grow(-BASISV(yDir));
   
  yEdgeLo.grow(-BASISV(xDir));
  yEdgeHi.grow(-BASISV(xDir));
   
  Vector<Real> xVal(3);
  Vector<Real> yVal(3);
   
  // iterate over x-lo side
  for (BoxIterator bit(xEdgeLo); bit.ok(); ++bit)
    {
      // use for case statement for stencil-type
      int oneSideDir = xDir;
      int centerDir  = (oneSideDir + 1)%2;
       
      // two-components for reuse in case iv is in a corner
      IntVect hiLoSide     = LARGEREALVAL*IntVect::Unit;
      hiLoSide[oneSideDir] = -1; 
       
      // IntVect
      iv = bit();
       
      // this box exludes corners
      cornerIv = false;
       
      // fill first derivatives and the pure second derivatives 
      calculateBdDerivIv(a_mapPartials,
                         iv           ,
                         a_curr       ,
                         oneSideDir   ,
                         hiLoSide     ,
                         a_dx         ,
                         cornerIv     ,
                         a_BFieldCase );
       
      // fill mixed partial using centered differences if not BField
      if (!a_BFieldCase)
        {
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              a_mapPartials(iv,idir + 4*SpaceDim)  = a_mapPartials(iv + BASISV(centerDir),idir + 4*SpaceDim); 
              a_mapPartials(iv,idir + 4*SpaceDim) -= a_mapPartials(iv - BASISV(centerDir),idir + 4*SpaceDim);
              a_mapPartials(iv,idir + 4*SpaceDim) /= (2.0*a_dx[centerDir])                                  ;
            }
        }
    }  

  // iterate over lo side
  for (BoxIterator bit(yEdgeLo); bit.ok(); ++bit)
    {
      // use for case statement for stencil-type
      int oneSideDir = yDir              ;
      int centerDir  = (oneSideDir + 1)%2;
       
      // two-components for reuse in case iv is in a corner
      IntVect hiLoSide     = LARGEREALVAL*IntVect::Unit;
      hiLoSide[oneSideDir] = -1                        ;
       
      // IntVect
      IntVect iv = bit();
       
      // this box exludes corners
      cornerIv = false;
       
      // fill first derivatives and the pure second derivatives 
      calculateBdDerivIv(a_mapPartials,
                         iv           ,
                         a_curr       ,
                         oneSideDir   ,
                         hiLoSide     ,
                         a_dx         ,
                         cornerIv     ,
                         a_BFieldCase );
       
      // fill mixed partial using centered differences if not BField
      if (!a_BFieldCase)
        {
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              a_mapPartials(iv,idir + 4*SpaceDim)  = a_mapPartials(iv + BASISV(centerDir),idir + 4*SpaceDim); 
              a_mapPartials(iv,idir + 4*SpaceDim) -= a_mapPartials(iv - BASISV(centerDir),idir + 4*SpaceDim);
              a_mapPartials(iv,idir + 4*SpaceDim) /= (2.0*a_dx[centerDir])                                  ;
            }
        }
    }

  // iterate over hi side
  for (BoxIterator bit(xEdgeHi); bit.ok(); ++bit)
    { 
      // use for case statement for stencil-type
      int oneSideDir = xDir;
      int centerDir  = (oneSideDir + 1)%2;
       
      // two-components for reuse in case iv is in a corner
      IntVect hiLoSide     = LARGEREALVAL*IntVect::Unit;
      hiLoSide[oneSideDir] = 1                         ;
       
      // IntVect
      iv = bit();
       
      // this box exludes corners
      cornerIv = false;
       
      // fill first derivatives and the pure second derivatives 
      calculateBdDerivIv(a_mapPartials,
                         iv           ,
                         a_curr       ,
                         oneSideDir   ,
                         hiLoSide     ,
                         a_dx         ,
                         cornerIv     ,
                         a_BFieldCase );
       
      // fill mixed partial using centered differences if not BField
      if (!a_BFieldCase)
        {
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              a_mapPartials(iv,idir + 4*SpaceDim)  = a_mapPartials(iv + BASISV(centerDir),idir + 4*SpaceDim); 
              a_mapPartials(iv,idir + 4*SpaceDim) -= a_mapPartials(iv - BASISV(centerDir),idir + 4*SpaceDim);
              a_mapPartials(iv,idir + 4*SpaceDim) /= (2.0*a_dx[centerDir])                                  ;
            }
        }
    }
   
  for (BoxIterator bit(yEdgeHi); bit.ok(); ++bit)
    {
      // use for case statement for stencil-type
      int oneSideDir = yDir              ;
      int centerDir  = (oneSideDir + 1)%2;
       
      // two-components for reuse in case iv is in a corner
      IntVect hiLoSide     = LARGEREALVAL*IntVect::Unit;
      hiLoSide[oneSideDir] = 1                         ;
       
      // IntVect
      iv = bit();
       
      // this box exludes corners
      cornerIv = false;
       
      // fill first derivatives and the pure second derivatives 
      calculateBdDerivIv(a_mapPartials,
                         iv           ,
                         a_curr       ,
                         oneSideDir   ,
                         hiLoSide     ,
                         a_dx         ,
                         cornerIv     ,
                         a_BFieldCase );
       
      // fill mixed partial using centered differences if not BField
      if (!a_BFieldCase)
        {
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              a_mapPartials(iv,idir + 4*SpaceDim)  = a_mapPartials(iv + BASISV(centerDir),idir + 4*SpaceDim); 
              a_mapPartials(iv,idir + 4*SpaceDim) -= a_mapPartials(iv - BASISV(centerDir),idir + 4*SpaceDim);
              a_mapPartials(iv,idir + 4*SpaceDim) /=                 (2.0*a_dx[centerDir]                  );
            }
        }
    }
}

void reComputeMapPartialsInteriorBd(FArrayBox       & a_mapPartials, 
                                    const FArrayBox & a_curr       ,
                                    const Box       & a_validBox   ,
                                    const RealVect  & a_dx         ,
                                    const bool      & a_BFieldCase ) 
{
  // macros for components of the mapping
  int xDir = 0;
  int yDir = 1;

  Box domainBox = a_validBox;

  // distinguish different stencils for the boundary ivs
  bool cornerIv = true;
  IntVect hiLoSide;
  
  // lower left corner
  IntVect loLo = domainBox.smallEnd();
  hiLoSide[xDir] = -1;
  hiLoSide[yDir] = -1; 
   
  // this variable reused for all boundary cases
  IntVect iv = LARGEREALVAL*IntVect::Unit;
   
  // IntVect
  iv = loLo;
   
  // corner iv
  cornerIv = true;
  
  // upper right corner
  IntVect hiHi = domainBox.bigEnd  ();
  hiLoSide[xDir] = 1;
  hiLoSide[yDir] = 1; 
   
  // IntVect
  iv = hiHi;
   
  // upper left corner
  IntVect loHi(loLo[0],hiHi[1]);
  hiLoSide[xDir] = -1;
  hiLoSide[yDir] =  1; 
   
  // IntVect
  iv = loHi;
   
  // lower right corner
  IntVect hiLo(hiHi[0],loLo[1]);
  hiLoSide[xDir] =  1;
  hiLoSide[yDir] = -1; 
   
  // IntVect
  iv = hiLo;
     
  // iteration boxes around the edges
  Box xEdgeLo(loLo,loHi);
  Box xEdgeHi(hiLo,hiHi);
   
  Box yEdgeLo(loLo,hiLo);
  Box yEdgeHi(loHi,hiHi);
   
  // iterate over the interior of the edges
  // xEdgeLo.grow(-BASISV(yDir));
  //xEdgeHi.grow(-BASISV(yDir));
   
  // yEdgeLo.grow(-BASISV(xDir));
  //yEdgeHi.grow(-BASISV(xDir));
   
  Vector<Real> xVal(3);
  Vector<Real> yVal(3);
   
  // iterate over x-lo side
  for (BoxIterator bit(xEdgeLo); bit.ok(); ++bit)
    {
      // use for case statement for stencil-type
      int oneSideDir = xDir;
      int centerDir  = (oneSideDir + 1)%2;
       
      // two-components for reuse in case iv is in a corner
      IntVect hiLoSide     = LARGEREALVAL*IntVect::Unit;
      hiLoSide[oneSideDir] = -1; 
       
      // IntVect
      iv = bit();
       
      // this box excludes corners
      cornerIv = false;
       
      // fill first derivatives and the pure second derivatives 
      calculateBdDerivIv(a_mapPartials,
                         iv           ,
                         a_curr       ,
                         oneSideDir   ,
                         hiLoSide     ,
                         a_dx         ,
                         cornerIv     ,
                         a_BFieldCase );
       
      // fill mixed partial using centered differences if not BField
      if (!a_BFieldCase)
        {
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              a_mapPartials(iv,idir + 4*SpaceDim)  = a_mapPartials(iv + BASISV(centerDir),idir + 4*SpaceDim); 
              a_mapPartials(iv,idir + 4*SpaceDim) -= a_mapPartials(iv - BASISV(centerDir),idir + 4*SpaceDim);
              a_mapPartials(iv,idir + 4*SpaceDim) /= (2.0*a_dx[centerDir])                                  ;
            }
        }
    }  

  // iterate over lo side
  for (BoxIterator bit(yEdgeLo); bit.ok(); ++bit)
    {
      // use for case statement for stencil-type
      int oneSideDir = yDir              ;
      int centerDir  = (oneSideDir + 1)%2;
       
      // two-components for reuse in case iv is in a corner
      IntVect hiLoSide     = LARGEREALVAL*IntVect::Unit;
      hiLoSide[oneSideDir] = -1                        ;
       
      // IntVect
      IntVect iv = bit();
       
      // this box exludes corners
      cornerIv = false;
       
      // fill first derivatives and the pure second derivatives 
      calculateBdDerivIv(a_mapPartials,
                         iv           ,
                         a_curr       ,
                         oneSideDir   ,
                         hiLoSide     ,
                         a_dx         ,
                         cornerIv     ,
                         a_BFieldCase );
       
      // fill mixed partial using centered differences if not BField
      if (!a_BFieldCase)
        {
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              a_mapPartials(iv,idir + 4*SpaceDim)  = a_mapPartials(iv + BASISV(centerDir),idir + 4*SpaceDim); 
              a_mapPartials(iv,idir + 4*SpaceDim) -= a_mapPartials(iv - BASISV(centerDir),idir + 4*SpaceDim);
              a_mapPartials(iv,idir + 4*SpaceDim) /= (2.0*a_dx[centerDir])                                  ;
            }
        }
    }

  // iterate over hi side
  for (BoxIterator bit(xEdgeHi); bit.ok(); ++bit)
    { 
      // use for case statement for stencil-type
      int oneSideDir = xDir;
      int centerDir  = (oneSideDir + 1)%2;
       
      // two-components for reuse in case iv is in a corner
      IntVect hiLoSide     = LARGEREALVAL*IntVect::Unit;
      hiLoSide[oneSideDir] = 1                         ;
       
      // IntVect
      iv = bit();
       
      // this box exludes corners
      cornerIv = false;
       
      // fill first derivatives and the pure second derivatives 
      calculateBdDerivIv(a_mapPartials,
                         iv           ,
                         a_curr       ,
                         oneSideDir   ,
                         hiLoSide     ,
                         a_dx         ,
                         cornerIv     ,
                         a_BFieldCase );
       
      // fill mixed partial using centered differences if not BField
      if (!a_BFieldCase)
        {
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              a_mapPartials(iv,idir + 4*SpaceDim)  = a_mapPartials(iv + BASISV(centerDir),idir + 4*SpaceDim); 
              a_mapPartials(iv,idir + 4*SpaceDim) -= a_mapPartials(iv - BASISV(centerDir),idir + 4*SpaceDim);
              a_mapPartials(iv,idir + 4*SpaceDim) /= (2.0*a_dx[centerDir])                                  ;
            }
        }
    }
   
  for (BoxIterator bit(yEdgeHi); bit.ok(); ++bit)
    {
      // use for case statement for stencil-type
      int oneSideDir = yDir              ;
      int centerDir  = (oneSideDir + 1)%2;
       
      // two-components for reuse in case iv is in a corner
      IntVect hiLoSide     = LARGEREALVAL*IntVect::Unit;
      hiLoSide[oneSideDir] = 1                         ;
       
      // IntVect
      iv = bit();
       
      // this box exludes corners
      cornerIv = false;
       
      // fill first derivatives and the pure second derivatives 
      calculateBdDerivIv(a_mapPartials,
                         iv           ,
                         a_curr       ,
                         oneSideDir   ,
                         hiLoSide     ,
                         a_dx         ,
                         cornerIv     ,
                         a_BFieldCase );
       
      // fill mixed partial using centered differences if not BField
      if (!a_BFieldCase)
        {
          for (int idir = 0; idir < SpaceDim; ++idir)
            {
              a_mapPartials(iv,idir + 4*SpaceDim)  = a_mapPartials(iv + BASISV(centerDir),idir + 4*SpaceDim); 
              a_mapPartials(iv,idir + 4*SpaceDim) -= a_mapPartials(iv - BASISV(centerDir),idir + 4*SpaceDim);
              a_mapPartials(iv,idir + 4*SpaceDim) /=                 (2.0*a_dx[centerDir]                  );
            }
        }
    }
}

void computeMapPartials(FArrayBox       & a_mapPartials,
                        const FArrayBox & a_extendedMap,
                        const Box       & a_validBox   ,
                        const RealVect  & a_dx         ,
                        const bool      & a_writeInfo  )
{
  // macros for components of the mapping
  int xDir = 0;
  int yDir = 1;

  Box domainBox = a_extendedMap.box();
 
  Real xDeriv = LARGEREALVAL;
  Real yDeriv = LARGEREALVAL;

  // shrink the box in the x- and y-directions
  Box intBox = domainBox;
  intBox.grow(-1);

  // iterate
  for (BoxIterator bit(intBox); bit.ok(); ++bit)
    { 
      // IntVect
      const IntVect iv = bit();
      
      for (int xiEtaDir = 0;xiEtaDir < SpaceDim; ++xiEtaDir)
        {
          firstDerivCenteredDiff(xDeriv       ,
                                 yDeriv       ,
                                 iv           ,
                                 xiEtaDir     ,
                                 a_dx         ,
                                 a_extendedMap,
                                 a_writeInfo  );
                      
          // assign first derivative
          a_mapPartials(iv,xDir + xiEtaDir*SpaceDim) = xDeriv;
          a_mapPartials(iv,yDir + xiEtaDir*SpaceDim) = yDeriv;
          
        }
    }

  // iterate over the interior
  for (BoxIterator bit(intBox); bit.ok(); ++bit)
    {
      // IntVect
      const IntVect iv = bit(); 
      
      // derivatives of the first two components: xXi and yXi 
      Real x2Xi  ;
      Real y2Xi  ;
      Real x2Eta ;
      Real y2Eta ;
      Real xXiEta;
      Real yXiEta;
      computeSecondDeriv(x2Xi         ,
                         y2Xi         ,
                         x2Eta        ,
                         y2Eta        ,
                         xXiEta       ,
                         yXiEta       ,
                         iv           ,
                         a_dx         ,
                         a_extendedMap,
                         a_writeInfo  );
      
      // xi deriv
      a_mapPartials(iv,xDir + 2*SpaceDim) = x2Xi ;
      a_mapPartials(iv,yDir + 2*SpaceDim) = y2Xi ;

      //eta deriv
      a_mapPartials(iv,xDir + 3*SpaceDim) = x2Eta;
      a_mapPartials(iv,yDir + 3*SpaceDim) = y2Eta;

      // mixed deriv
      a_mapPartials(iv,xDir + 4*SpaceDim) = xXiEta;
      a_mapPartials(iv,yDir + 4*SpaceDim) = yXiEta;
    }

  // flag for computeMapPartialsBdIv
  bool BFieldCase = false;

  // set the corners, iterate over four edge boxes using sometimes centered and sometimes one-sided stencils (which are based on quadratic extrapolation)
  computeMapPartialsBdIv(a_mapPartials      , 
                         a_extendedMap      ,
                         a_extendedMap.box(),
                         a_dx               ,
                         BFieldCase         );


  Box grownValidBox = a_validBox;
  //grownValidBox.grow(1);

  //recompute derivatives around the boundary of the valid box
  reComputeMapPartialsInteriorBd(a_mapPartials, 
                                 a_extendedMap,
                                 grownValidBox,
                                 a_dx         ,
                                 BFieldCase   );
#if 0
  grownValidBox.grow(1);
  reComputeMapPartialsInteriorBd(a_mapPartials, 
                                 a_extendedMap,
                                 grownValidBox,
                                 a_dx         ,
                                 BFieldCase   );
#endif
}

void SmJnRes(Real            & a_res0       ,
             Real            & a_res1       , 
             const IntVect   & a_iv         ,
             const FArrayBox & a_mapPartials,
             const FArrayBox & a_Jacobian   ,
             const bool      & a_writeInfo  )      
{
  // indexing macros
  //int xiSmooth  = 0;
  //int etaSmooth = 1;
  
  int xDir = 0;
  int yDir = 1;

  int xiDir = 0;
  int etaDir = 1;
  
  int m = LARGEREALVAL;
  int n = LARGEREALVAL;
    
  ParmParse pp;
  pp.get("Spower",m);
  pp.get("Jpower",n);
  
  Real x_xi    = a_mapPartials(a_iv,xDir + 0*SpaceDim);
  Real y_xi    = a_mapPartials(a_iv,yDir + 0*SpaceDim);
  
  Real x_eta   = a_mapPartials(a_iv,xDir + 1*SpaceDim);
  Real y_eta   = a_mapPartials(a_iv,yDir + 1*SpaceDim);
  
  Real x_2xi   = a_mapPartials(a_iv,xDir + 2*SpaceDim);
  Real y_2xi   = a_mapPartials(a_iv,yDir + 2*SpaceDim);
  
  Real x_2eta  = a_mapPartials(a_iv,xDir + 3*SpaceDim);
  Real y_2eta  = a_mapPartials(a_iv,yDir + 3*SpaceDim);
  
  Real x_xi_eta = a_mapPartials(a_iv,xDir + 4*SpaceDim);
  Real y_xi_eta = a_mapPartials(a_iv,yDir + 4*SpaceDim);
  
  // mixed partials are equal
  Real x_eta_xi = x_xi_eta;
  Real y_eta_xi = y_xi_eta;

  // sum of squares
  Real S  = (x_xi*x_xi + y_xi*y_xi) + (x_eta*x_eta + y_eta*y_eta);
  
  Real Sx_xi  = 2.0*x_xi ;
  Real Sy_xi  = 2.0*y_xi ;
  Real Sx_eta = 2.0*x_eta;
  Real Sy_eta = 2.0*y_eta;

  Real Sx_xi_xi   =  2.0*x_2xi ;
  Real Sy_xi_xi   =  2.0*y_2xi ;
  Real Sx_eta_eta =  2.0*x_2eta;
  Real Sy_eta_eta =  2.0*y_2eta;

  Real S_xi  = Sx_xi*x_2xi    + Sy_xi*y_2xi    + Sx_eta*x_eta_xi + Sy_eta*y_eta_xi;
  Real S_eta = Sx_xi*x_xi_eta + Sy_xi*y_xi_eta + Sx_eta*x_2eta   + Sy_eta*y_2eta  ;
   
  // Jacobian
  Real J  = a_Jacobian(a_iv,0);                           //0
   
  // derivatives with respect to derivatives of the mapping
  Real Jx_xi  = a_Jacobian(a_iv,1 + 0*SpaceDim + xDir    );//1
  Real Jy_xi  = a_Jacobian(a_iv,1 + 0*SpaceDim + yDir    );//2
  Real Jx_eta = a_Jacobian(a_iv,1 + 1*SpaceDim + xDir    );//3
  Real Jy_eta = a_Jacobian(a_iv,1 + 1*SpaceDim + yDir    );//4

  // derivatives of the above with respect to xi and eta
  Real Jx_xi_xi   = a_Jacobian(a_iv,1 + 2*SpaceDim + xDir);//5
  Real Jy_xi_xi   = a_Jacobian(a_iv,1 + 2*SpaceDim + yDir);//6
  Real Jx_eta_eta = a_Jacobian(a_iv,1 + 3*SpaceDim + xDir);//7
  Real Jy_eta_eta = a_Jacobian(a_iv,1 + 3*SpaceDim + yDir);//8

  // derivatives with respect to xi and eta
  Real J_xi  = a_Jacobian(a_iv,1 + 4*SpaceDim + xiDir    );//9
  Real J_eta = a_Jacobian(a_iv,1 + 4*SpaceDim + etaDir   );//10

  Real L1 = pow(S,m);
  Real L2 = pow(J,n);
    
  //Lagrangian; L = L1*L2;
  Real L1_xi  = m*pow(S,m-1)*S_xi ;
  Real L1_eta = m*pow(S,m-1)*S_eta;

  Real L2_xi  = n*pow(J,n-1)*J_xi ;
  Real L2_eta = n*pow(J,n-1)*J_eta;

  Real L1x_xi  = m*pow(S,m-1)*Sx_xi;
  Real L1y_xi  = m*pow(S,m-1)*Sy_xi;
  
  Real L1x_eta = m*pow(S,m-1)*Sx_eta;
  Real L1y_eta = m*pow(S,m-1)*Sy_eta;

  Real L1x_xi_xi   = m*(m-1)*pow(S,m-2)*S_xi *Sx_xi  + m*pow(S,m-1)*Sx_xi_xi  ;
  Real L1y_xi_xi   = m*(m-1)*pow(S,m-2)*S_xi *Sy_xi  + m*pow(S,m-1)*Sy_xi_xi  ;
  Real L1x_eta_eta = m*(m-1)*pow(S,m-2)*S_eta*Sx_eta + m*pow(S,m-1)*Sx_eta_eta;
  Real L1y_eta_eta = m*(m-1)*pow(S,m-2)*S_eta*Sy_eta + m*pow(S,m-1)*Sy_eta_eta;

  Real L2x_xi  = n*pow(J,n-1)*Jx_xi;
  Real L2y_xi  = n*pow(J,n-1)*Jy_xi;
  
  Real L2x_eta = n*pow(J,n-1)*Jx_eta;
  Real L2y_eta = n*pow(J,n-1)*Jy_eta;

  Real L2x_xi_xi   = n*(n-1)*pow(J,n-2)*J_xi *Jx_xi  + n*pow(J,n-1)*Jx_xi_xi ;
  Real L2y_xi_xi   = n*(n-1)*pow(J,n-2)*J_xi *Jy_xi  + n*pow(J,n-1)*Jy_xi_xi ;
  Real L2x_eta_eta = n*(n-1)*pow(J,n-2)*J_eta*Jx_eta + n*pow(J,n-1)*Jx_eta_eta;
  Real L2y_eta_eta = n*(n-1)*pow(J,n-2)*J_eta*Jy_eta + n*pow(J,n-1)*Jy_eta_eta;

  //Real Lx_xi  = L1x_xi* L2 + L1*L2x_xi ;
  //Real Ly_xi  = L1y_xi* L2 + L1*L2y_xi ;

  //Real Lx_eta = L1x_eta*L2 + L1*L2x_eta;
  //Real Ly_eta = L1y_eta*L2 + L1*L2y_eta;

  Real Lx_xi_xi   = L1x_xi_xi*L2   + L1x_xi*L2_xi   + L1_xi*L2x_xi   + L1*L2x_xi_xi  ;
  Real Ly_xi_xi   = L1y_xi_xi*L2   + L1y_xi*L2_xi   + L1_xi*L2y_xi   + L1*L2y_xi_xi  ;
  Real Lx_eta_eta = L1x_eta_eta*L2 + L1x_eta*L2_eta + L1_eta*L2x_eta + L1*L2x_eta_eta;
  Real Ly_eta_eta = L1y_eta_eta*L2 + L1y_eta*L2_eta + L1_eta*L2y_eta + L1*L2y_eta_eta;

  // no term in S includes x or includes y (explicitly)
  Real Lx = 0;
  Real Ly = 0;

  // return residual
  a_res0 = -(Lx_xi_xi + Lx_eta_eta) + Lx;
  a_res1 = -(Ly_xi_xi + Ly_eta_eta) + Ly;
}

void computeInc(FArrayBox       & a_inc           ,
                FArrayBox       & a_smoothRes     ,
                const Box       & a_fixedBox      ,
                const FArrayBox & a_curr          ,
                const FArrayBox & a_mapPartials   ,
                const FArrayBox & a_smoothDiag    ,
                const FArrayBox & a_Jacobian      ,
                const bool      & a_writeInfo     ,
                const RealVect  & a_dx            ,
                const int       & a_iter          )
{
  // iteration boxes
  Box box    = a_inc.box();
  Box intBox = box.grow(-1);

  // iterate over the interior to get data for boundary op
  for (BoxIterator bit(intBox); bit.ok(); ++bit)
    {
      // IntVect
      const IntVect& iv = bit();  
      if (!a_fixedBox.contains(iv))
        {
          // indexing macros
          int xDir = 0;
          int yDir = 1;
          
          // smooth res
          Real smoothRes0 = LARGEREALVAL;
          Real smoothRes1 = LARGEREALVAL;
          
          // smooth diag term
          Real smoothDiag0 = a_smoothDiag(iv,xDir);
          Real smoothDiag1 = a_smoothDiag(iv,yDir);
          
          SmJnRes(smoothRes0   ,
                  smoothRes1   , 
                  iv           ,
                  a_mapPartials,
                  a_Jacobian   ,
                  a_writeInfo );
          
          a_smoothRes(iv,xDir) = smoothRes0;
          a_smoothRes(iv,yDir) = smoothRes1;
          
          Real diagTerm0 = smoothDiag0;
          Real diagTerm1 = smoothDiag1;
          
          if (diagTerm0 == 0.0 || diagTerm1 == 0)
            {
              MayDay::Abort("No diagonal term for Euler-Lagrange equation");
            }
          
          Real lambda0 = -0.5/diagTerm0;
          Real lambda1 = -0.5/diagTerm1;
          
          if (a_writeInfo)
            {
              pout()<<"lambda0 = "<< lambda0 << endl;
              pout()<<"lambda1 = "<< lambda1 << endl;
            }
          
          a_inc(iv,xDir) = lambda0*smoothRes0;
          a_inc(iv,yDir) = lambda1*smoothRes1;
          
          a_smoothRes(iv,xDir) = smoothRes0;
          a_smoothRes(iv,yDir) = smoothRes1;
        }
    }
}
  
IntVect findXiEtaXPt(bool      & a_hasXPt,
                     const int & a_block )
{
  a_hasXPt = true;
  int xDir = 0;
  int yDir = 1;

  IntVect retval;
  if (a_block == 0 )
    {
      retval = BASISV(xDir) + BASISV(yDir);
    }
  else if (a_block == 1)
    {
      retval = BASISV(xDir);
    }
  else if (a_block == 2)
    {
      retval = BASISV(yDir);
    }
  else if (a_block == 3)
    {
      retval = IntVect::Zero;
    }
  else if (a_block == 4)
    {
      retval = IntVect::Zero;
    }
  else if (a_block == 5)
    {
      retval = BASISV(yDir);
    }
  else if (a_block == 6)
    {
      retval = BASISV(xDir);
    }
  else if (a_block == 7)
    {
      retval = BASISV(xDir) + BASISV(yDir);
    }
  else
    {
      MayDay::Abort("bad block number");
    }
  return retval;
}

bool edgeContainsXPoint(bool              & a_xPointHiSide,
                        const std::string & a_blockName   ,
                        const std::string & a_edge        )
{
  bool retval;
  
  if (a_blockName == "lcore")
    {
      if (a_edge == "bottom")
        {
          retval = false;
        }
      else if (a_edge == "right")
        { 
          a_xPointHiSide = true;
          retval         = true; 
        }
      else if (a_edge == "top")
        {
          a_xPointHiSide = true;
          retval         = true;
        }
      else if (a_edge == "left")
        {
          retval = false;
        }
      else
        {
          MayDay::Abort("bad string for edge");
        }
    }
  else if (a_blockName == "rcore")  
    {
      if (a_edge == "bottom")
        {
          a_xPointHiSide = true;
          retval         = true;
        }
      else if (a_edge == "right")
        {
          a_xPointHiSide = false;
          retval         = true;
        }
      else if (a_edge == "top")
        {
          retval = false;
        }
      else if (a_edge == "left")
        {
          retval = false;
        }
      else
        {
          MayDay::Abort("bad string for edge");
        }
    }
  else if (a_blockName == "lcsol")
    {
      if (a_edge == "bottom")
        {
          retval = false;
        }
      else if (a_edge == "right")
        {
          retval = false;
        }
      else if (a_edge == "top")
        {
          a_xPointHiSide = false;
          retval         = true;
        }
      else if (a_edge == "left")
        {
          a_xPointHiSide = true;
          retval         = true;
        }
      else
        {
          MayDay::Abort("bad string for edge");
        }
    }
  else  if (a_blockName == "rcsol") 
    {
      if (a_edge == "bottom")
        {
          a_xPointHiSide = false;
          retval         = true;
        }
      else if (a_edge == "right")
        {
          retval = false;
        }
      else if (a_edge == "top")
        {
          retval = false;
        }
      else if (a_edge == "left")
        {
          a_xPointHiSide = false;
          retval         = true;
        }
      else
        {
          MayDay::Abort("bad string for edge");
        }
    }
  else if (a_blockName == "lsol")
    {
      if (a_edge == "bottom")
        {
          a_xPointHiSide = false;
          retval         = true;
        }
      else if (a_edge == "right")
        {
          retval = false;
        }
      else if (a_edge == "top")
        {
          retval = false;
        }
      else if (a_edge == "left")
        {
          a_xPointHiSide = false;
          retval         =  true;
        }
      else
        {
          MayDay::Abort("bad string for edge");
        }
    }
  else  if (a_blockName == "rsol")   
    {
      if (a_edge == "bottom")
        {
          retval = false;
        }
      else if (a_edge == "right")
        {
          retval = false;
        }
      else if (a_edge == "top")
        {
          a_xPointHiSide = false;
          retval         = true;
        }
      else if (a_edge == "left")
        {
          a_xPointHiSide = true;
          retval         = true;
        }
      else
        {
          MayDay::Abort("bad string for edge");
        }
    }
  else if (a_blockName == "lpf") 
    {
      if (a_edge == "bottom")
        {
          a_xPointHiSide = true;
          retval         =  true;
        }
      else if (a_edge == "right")
        {
          a_xPointHiSide = false;
          retval         = true;
        }
      else if (a_edge == "top")
        {
          
          retval = false;
        }
      else if (a_edge == "left")
        {
          retval = false;
        }
      else
        {
          MayDay::Abort("bad string for edge");
        }
    }
  else if (a_blockName == "rpf")  
    {
      if (a_edge == "bottom")
        {
          retval = false;
        }
      else if (a_edge == "right")
        { 
          a_xPointHiSide = true;
          retval         = true;
        }
      else if (a_edge == "top")
        {
          a_xPointHiSide = true;
          retval         = true;
        }
      else if (a_edge == "left")
        {
          retval = false;
        }
      else
        {
          MayDay::Abort("bad string for edge");
        }
    }
  else if (a_blockName == "mcore")  
    {
      if (a_edge == "bottom")
        {
          retval = false;
        }
      else if (a_edge == "right")
        { 
          retval = false;
        }
      else if (a_edge == "top")
        {
          retval = false;
        }
      else if (a_edge == "left")
        {
          retval = false;
        }
      else
        {
          MayDay::Abort("bad string for edge");
        }
    }
  else if (a_blockName == "mcsol")  
    {
      if (a_edge == "bottom")
        {
          retval = false;
        }
      else if (a_edge == "right")
        { 
          retval = false;
        }
      else if (a_edge == "top")
        {
          retval = false;
        }
      else if (a_edge == "left")
        {
          retval = false;
        }
      else
        {
          MayDay::Abort("bad string for edge");
        }
    }
  else
    {
      MayDay::Abort("bad string for block");
    }
  return retval;
}

void ghostCellPlan(int       & a_getGhostFromBlockNumber,
                   int       & a_getGhostFromFace       , 
                   bool      & a_useExtrapolation       ,
                   const int & a_needGhostForBlock      ,
                   const int & a_needGhostForFace       )
{
  // faces in (xi,eta) blocks are numbered 0,1,2,3 (counter-clockwise from bottom): 
  int bottom = 0 ; 
  int right  = 1;
  int top    = 2;
  int left   = 3;

  // face value for the extrapolation case
  int invalid = LARGEINTVAL;

  // block numbering
  int lcore = 0;
  int rcore = 1;
  int lcsol = 2;
  int rcsol = 3;
  int lsol  = 4;
  int rsol  = 5;
  int lpf   = 6;
  int rpf   = 7;
  int mcore = 8;
  int mcsol = 9;
  
  // block 0
  if (a_needGhostForBlock == lcore)
    {
      if (a_needGhostForFace == bottom)
        {
          a_getGhostFromBlockNumber = mcore;
          a_getGhostFromFace        = top  ;
          a_useExtrapolation        = true;
        }
      if (a_needGhostForFace == right)
        {
          a_getGhostFromBlockNumber = lcsol;
          a_getGhostFromFace        = left ;
          a_useExtrapolation        = true ;
        } 
      if (a_needGhostForFace == top)
        {
          a_getGhostFromBlockNumber = rcore ;
          a_getGhostFromFace        = bottom;
          a_useExtrapolation        = true  ;
        }

      if (a_needGhostForFace == left)
        {
          a_getGhostFromBlockNumber = invalid;       
          a_getGhostFromFace        = invalid;      
          a_useExtrapolation        = true   ;      
        } 
    }

  // block 1
  if (a_needGhostForBlock == rcore)
    {
      if (a_needGhostForFace == bottom)
        {
          a_getGhostFromBlockNumber = lcore;
          a_getGhostFromFace        = top  ;
          a_useExtrapolation        = true ;
        }
      if (a_needGhostForFace == right)
        {
          a_getGhostFromBlockNumber = rcsol;
          a_getGhostFromFace        = left ;
          a_useExtrapolation        = true ;
        } 

      if (a_needGhostForFace == top)
        {
          a_getGhostFromBlockNumber = mcore ;
          a_getGhostFromFace        = bottom;
          a_useExtrapolation        = true  ;
        }

      if (a_needGhostForFace == left)
        {
          a_getGhostFromBlockNumber = invalid;       
          a_getGhostFromFace        = invalid;      
          a_useExtrapolation        = true   ;     
        } 
    }

  // block 2
  if (a_needGhostForBlock == lcsol)
    {
      if (a_needGhostForFace == bottom)
        {
          a_getGhostFromBlockNumber = mcsol;
          a_getGhostFromFace        = top  ;
          a_useExtrapolation        = true;
        }
      if (a_needGhostForFace == right)
        {
          a_getGhostFromBlockNumber = invalid;       
          a_getGhostFromFace        = invalid;      
          a_useExtrapolation        = true   ;    
        } 

      if (a_needGhostForFace == top)
        {
          a_getGhostFromBlockNumber = lsol  ;
          a_getGhostFromFace        = bottom;
          a_useExtrapolation        = true  ;
        }

      if (a_needGhostForFace == left)
        {
          a_getGhostFromBlockNumber = lcore;       
          a_getGhostFromFace        = right;      
          a_useExtrapolation        = true ;
        } 
    }
  
  // block 3
  if (a_needGhostForBlock == rcsol)
    {
      if (a_needGhostForFace == bottom)
        {
          a_getGhostFromBlockNumber = rsol ;
          a_getGhostFromFace        = top  ;
          a_useExtrapolation        = true ;
        }
      if (a_needGhostForFace == right)
        {
          a_getGhostFromBlockNumber = invalid;       
          a_getGhostFromFace        = invalid;      
          a_useExtrapolation        = true   ;      
        } 

      if (a_needGhostForFace == top)
        {
          a_getGhostFromBlockNumber = mcsol ;
          a_getGhostFromFace        = bottom;
          a_useExtrapolation        = true  ;
        }

      if (a_needGhostForFace == left)
        {
          a_getGhostFromBlockNumber = rcore;       
          a_getGhostFromFace        = right;      
          a_useExtrapolation        = true ;
        } 
    }

  //block 4
  if (a_needGhostForBlock == lsol)
    {
      if (a_needGhostForFace == bottom)
        {
          a_getGhostFromBlockNumber = lcsol;
          a_getGhostFromFace        = top  ;
          a_useExtrapolation        = true ;
        }
      if (a_needGhostForFace == right)
        {
          a_getGhostFromBlockNumber = invalid;       
          a_getGhostFromFace        = invalid;      
          a_useExtrapolation        = true   ;      
        } 

      if (a_needGhostForFace == top)
        {
          a_getGhostFromBlockNumber = invalid;       
          a_getGhostFromFace         = invalid;      
          a_useExtrapolation         = true   ;      
        }

      if (a_needGhostForFace == left)
        {
          a_getGhostFromBlockNumber = lpf  ;
          a_getGhostFromFace        = right;
          a_useExtrapolation        = true ;
        } 
    }

  //block 5
  if (a_needGhostForBlock == rsol)
    {
      if (a_needGhostForFace == bottom)
        {
          a_getGhostFromBlockNumber = invalid;       
          a_getGhostFromFace        = invalid;      
          a_useExtrapolation        = true   ;      
        }
      if (a_needGhostForFace == right)
        {
          a_getGhostFromBlockNumber = invalid;       
          a_getGhostFromFace        = invalid;      
          a_useExtrapolation        = true   ;      
        } 

      if (a_needGhostForFace == top)
        {
          a_getGhostFromBlockNumber = rcsol ;
          a_getGhostFromFace        = bottom;
          a_useExtrapolation        = true  ;
        }

      if (a_needGhostForFace == left)
        {
          a_getGhostFromBlockNumber = rpf  ;       
          a_getGhostFromFace        = right;      
          a_useExtrapolation        = true ;
        } 
    }

  // block 6
  if (a_needGhostForBlock == lpf)
    {
      if (a_needGhostForFace == bottom)
        {
          a_getGhostFromBlockNumber = rpf ;
          a_getGhostFromFace        = top ;
          a_useExtrapolation        = true;
        }
      if (a_needGhostForFace == right)
        {
          a_getGhostFromBlockNumber = lsol;
          a_getGhostFromFace        = left;
          a_useExtrapolation        = true;
        } 

      if (a_needGhostForFace == top)
        {
          a_getGhostFromBlockNumber = invalid;
          a_getGhostFromFace        = invalid;
          a_useExtrapolation        = true   ;
        }

      if (a_needGhostForFace == left)
        {
          a_getGhostFromBlockNumber = invalid;       
          a_getGhostFromFace        = invalid;      
          a_useExtrapolation        = true   ;      
        } 
    }

  //block 7 
  if (a_needGhostForBlock == rpf)
    {
      if (a_needGhostForFace == bottom)
        {
          a_getGhostFromBlockNumber = invalid;
          a_getGhostFromFace        = invalid;
          a_useExtrapolation        = true   ;
        }
      if (a_needGhostForFace == right)
        {
          a_getGhostFromBlockNumber = rsol;
          a_getGhostFromFace        = left;
          a_useExtrapolation        = true;
        } 

      if (a_needGhostForFace == top)
        {
          a_getGhostFromBlockNumber = lpf   ;
          a_getGhostFromFace        = bottom;
          a_useExtrapolation        = true  ;
        }

      if (a_needGhostForFace == left)
        {
          a_getGhostFromBlockNumber = invalid;       
          a_getGhostFromFace        = invalid;      
          a_useExtrapolation        = true   ;      
        } 
    }

  //block 8
  if (a_needGhostForBlock == mcore)
    {
      if (a_needGhostForFace == bottom)
        {
          a_getGhostFromBlockNumber = rcore;
          a_getGhostFromFace        = top  ;
          a_useExtrapolation        = true ;
        }
      if (a_needGhostForFace == right)
        {
          a_getGhostFromBlockNumber = mcsol;
          a_getGhostFromFace        = left ;
          a_useExtrapolation        = true;
        } 

      if (a_needGhostForFace == top)
        {
          a_getGhostFromBlockNumber = lcore ;
          a_getGhostFromFace        = bottom;
          a_useExtrapolation        = true  ;
        }

      if (a_needGhostForFace == left)
        {
          a_getGhostFromBlockNumber = invalid;       
          a_getGhostFromFace        = invalid;      
          a_useExtrapolation        = true   ;      
        } 
    }

  // block 9
  if (a_needGhostForBlock == mcsol)
    {
      if (a_needGhostForFace == bottom)
        {
          a_getGhostFromBlockNumber = rcsol;
          a_getGhostFromFace        = top  ;
          a_useExtrapolation        = true ;
        }
      if (a_needGhostForFace == right)
        {
          a_getGhostFromBlockNumber = invalid;
          a_getGhostFromFace        = invalid;
          a_useExtrapolation        = true   ;
        } 

      if (a_needGhostForFace == top)
        {
          a_getGhostFromBlockNumber = lcsol ;
          a_getGhostFromFace        = bottom;
          a_useExtrapolation        = true  ;
        }

      if (a_needGhostForFace == left)
        {
          a_getGhostFromBlockNumber = mcore;       
          a_getGhostFromFace        = right;      
          a_useExtrapolation        = true;      
        } 
    }
}

void computeError(FArrayBox       & a_error      ,
                  const FArrayBox & a_exactAnswer,
                  const FArrayBox & a_extendedMap)
{
  for (BoxIterator bit(a_error.box());bit.ok();++bit) 
    {
      IntVect iv = bit();
      for (int iComp = 0; iComp < SpaceDim; ++ iComp)
        {
          a_error(iv,iComp) = a_exactAnswer(iv,iComp) - a_extendedMap(iv,iComp);
        }
    }
}

void testFunction(FArrayBox     & a_fab        ,
                  FArrayBox     & a_exactAnswer,
                  const Box     & a_validBox   ,
                  const RealVect& a_origin     ,
                  const RealVect& a_dx         )
 
{
  RealVect pt;
  
  Vector<Real> a3;
  Vector<Real> a2;
  Vector<Real> a1;
  Vector<Real> a0;

  ParmParse pp;
  pp.getarr("a3",a3,0,SpaceDim);
  pp.getarr("a2",a2,0,SpaceDim);
  pp.getarr("a1",a1,0,SpaceDim);
  pp.getarr("a0",a0,0,SpaceDim);
  
  for (BoxIterator bit(a_validBox);bit.ok();++bit) 
    {
      IntVect iv = bit();
      for (int idir = 0; idir < SpaceDim; ++ idir)
        {
          pt[idir] = a_origin[idir] + a_dx[idir]*iv[idir]; 
        }

      for (int iComp = 0; iComp < SpaceDim;++iComp)
        {
          a_fab(iv,iComp)   = a3[iComp]*pt[iComp]*pt[iComp]*pt[iComp] + a2[iComp]*pt[iComp]*pt[iComp] + a1[iComp]*pt[iComp] + a0[iComp];
        }
    }

  for (BoxIterator bit(a_exactAnswer.box());bit.ok();++bit) 
    {
      IntVect iv = bit();
      for (int idir = 0; idir < SpaceDim; ++ idir)
        {
          pt[idir] = a_origin[idir] + a_dx[idir]*iv[idir]; 
        }

      for (int iComp = 0; iComp < SpaceDim;++iComp)
        {
          a_exactAnswer(iv,iComp) = a3[iComp]*pt[iComp]*pt[iComp]*pt[iComp] + a2[iComp]*pt[iComp]*pt[iComp] + a1[iComp]*pt[iComp] + a0[iComp];
          
        }
    }
}

Real extrapolate(Real               & a_extrapPt,
                 const Vector<Real> & a_yVal    ,
                 const Vector<Real> & a_xPts    ,
                 const Real         & a_dx      )
{
  Real retval = LARGEREALVAL;

  CH_assert (a_yVal.size() == a_xPts.size());
  
  // order of extrapolation
  int order = a_yVal.size();
  
  if (order == 0)
    {
      MayDay::Warning("a_points has size 0");
      retval = -LARGEREALVAL;
    }
  else
    {
      // elements of Lagrange formula for evenly spaced points
      Vector<Real>basisPoly(order);
      
      // set to 1.0 for ease of applying a product formula
      for(int j = 0; j < order; ++ j)
        {
          basisPoly[j] = 1.0;
        }
      
      for (int j = 0; j < order; ++j)
        {
          for(int m = 0; m < order; ++ m)
            {
              if (j != m)
                {
                  basisPoly[j] *= (a_extrapPt - a_xPts[m])/(a_xPts[j] - a_xPts[m]); 
                }
            }
        }
      
      Real LagrangePoly = 0.0;
      for (int j = 0; j < order; ++j)
        {
          LagrangePoly += a_yVal[j]*basisPoly[j];
        }
      retval = LagrangePoly;
    }
  
  return retval;
}


IntVect cornerExtrapDirections(const std::string& a_blockName, 
                               const std::string& a_corner   )
{
  // for indexing
  int xDir = 0;
  int yDir = 1;

  IntVect xOnly = BASISV(xDir);
  IntVect yOnly = BASISV(yDir);
  IntVect xAndY = IntVect::Unit;

  IntVect retval;
  if (a_blockName == "lcore")
    {
      if (a_corner == "lowerLeft")
        {
          retval = xOnly;
        }
      else if (a_corner == "lowerRight")
        {
          retval = xOnly;
        }
      else if (a_corner == "upperRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperLeft")
        {
          retval = xAndY;
        }
      else
        {
          MayDay::Abort("bad value for corner");
        }
    }
  else if (a_blockName == "rcore")
    {
      if (a_corner == "lowerLeft")
        {
          retval = xOnly;
        }
      else if (a_corner == "lowerRight")
        {
          retval = xOnly;
        }
      else if (a_corner == "upperRight")
        {
          retval = yOnly;
        }
      else if (a_corner == "upperLeft")
        {
          retval = yOnly;
        }
      else
        {
          MayDay::Abort("bad value for corner");
        }
    }
  else if (a_blockName == "lcsol")
    {
      if (a_corner == "lowerLeft")
        {
          retval = xAndY;
        }
      else if (a_corner == "lowerRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperLeft")
        { retval = xAndY;
        }
      else
        {
          MayDay::Abort("bad value for corner");
        }
    }
  else if (a_blockName == "rcsol")
    {
      
      if (a_corner == "lowerLeft")
        {
          retval = xAndY;
        }
      else if (a_corner == "lowerRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperLeft")
        {
          retval = xAndY;
        }
      else
        {
          MayDay::Abort("bad value for corner");
        }
    }
  else if (a_blockName == "lsol")
    {
      if (a_corner == "lowerLeft")
        {
          retval = xAndY;
        }
      else if (a_corner == "lowerRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperLeft")
        {
          retval = xAndY;
        }
      else
        {
          MayDay::Abort("bad value for corner");
        }
    }
  else if (a_blockName == "rsol")
    {
      if (a_corner == "lowerLeft")
        {
          retval = xAndY;
        }
      else if (a_corner == "lowerRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperLeft")
        {
          retval = xAndY;
        }
      else
        {
          MayDay::Abort("bad value for corner");
        }
    }
  else if (a_blockName == "lpf")
    {
      if (a_corner == "lowerLeft")
        {
          retval = xAndY;
        }
      else if (a_corner == "lowerRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperLeft")
        {
          retval = xAndY;
        }
      else
        {
          MayDay::Abort("bad value for corner");
        }
    }
  else if (a_blockName == "rpf")
    {
      if (a_corner == "lowerLeft")
        {
          retval = xAndY;
        }
      else if (a_corner == "lowerRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperRight")
        {
          retval = xAndY;
        }
      else if (a_corner == "upperLeft")
        {
          retval = xAndY;
        }
      else
        {
          MayDay::Abort("bad value for corner");
        }
    }
  else 
    {
      MayDay::Abort("bad name for blockname");
    }

  return retval;
}

RealVect cornerExtrap(const IntVect     & a_baseIv   ,
                      const IntVect     & a_offsetDir,
                      const FArrayBox   & a_fab      ,
                      const int         & a_order    ,
                      const RealVect    & a_dx       ,
                      const std::string & a_blockName,
                      const std::string & a_corner   )
{    
  // use one or both of the horizontal and vertical extrapolations
  IntVect directions = cornerExtrapDirections(a_blockName,
                                              a_corner   );
    
  // for indexing
  int xDir = 0;
  int yDir = 1;

  RealVect retval = LARGEREALVAL*RealVect::Unit;

  for (int iComp = 0; iComp < a_fab.nComp(); ++iComp)
    {
      //use stencil based on previously extrapolated vertical or horizontal cells
      RealVect extrapVal;
      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          Vector<Real> funcVal(a_order);
          Vector<Real> argPts (a_order);
            
          for (int iExtrap = 0; iExtrap < a_order; ++iExtrap)
            {
              IntVect extrapOffset = a_offsetDir[idir]*(iExtrap + 1)*BASISV(idir);
              funcVal[iExtrap] = a_fab(a_baseIv + extrapOffset,iComp);
              argPts [iExtrap] = a_dx[idir]*(a_baseIv[idir] + extrapOffset[idir]);
            }
            
          Real extrapPt  =  a_dx[idir]*(a_baseIv[idir]);
          extrapVal[idir] =  extrapolate(extrapPt   ,
                                         funcVal    ,
                                         argPts     ,
                                         a_dx[idir]);
        }
      
        
      if (directions == IntVect::Unit)
        {
          retval[iComp] = 0.5*(extrapVal[0] + extrapVal[1]);
        }
      else if (directions == BASISV(xDir))
        {
          retval[iComp] = extrapVal[0];
        }
      else if(directions == BASISV(yDir))
        {
          retval[iComp] = extrapVal[1];
        }
    }

  return retval;
}

void fillInCorners(FArrayBox         & a_extendedFab ,
                   const Box         & a_validBox    ,
                   const int         & a_order       ,
                   const RealVect    & a_dx          ,
                   const std::string & a_blockName   )
{
  IntVect lo = a_validBox.smallEnd();
  IntVect offsetDirLo(1,1);

  IntVect hi = a_validBox.bigEnd  ();
  IntVect offsetDirHi(-1,-1);

  IntVect loHi(lo[0],hi[1]);
  IntVect offsetDirLoHi(1,-1);

  IntVect hiLo(hi[0],lo[1]);
  IntVect offsetDirHiLo(-1,1);

  RealVect avgVal = LARGEREALVAL*RealVect::Unit;
  
  // for selecting extrapolation direction(s) in the corner cell
  std::string corner;

  // lower left
  IntVect baseIv  = lo;
  corner = "lowerLeft";
  avgVal = cornerExtrap(baseIv       ,
                        offsetDirLo  ,
                        a_extendedFab,
                        a_order      ,
                        a_dx         ,
                        a_blockName  ,
                        corner       );
  
  for (int iComp = 0; iComp < SpaceDim;++iComp)
    {
      a_extendedFab(lo,iComp) = avgVal[iComp];
    }

  // hi (upper right)
  baseIv  = hi;
  corner = "upperRight";
  avgVal = cornerExtrap(baseIv       ,
                        offsetDirHi  ,
                        a_extendedFab,
                        a_order      ,
                        a_dx         ,
                        a_blockName  ,
                        corner       );
  
  for (int iComp = 0; iComp < SpaceDim;++iComp)
    {
      a_extendedFab(hi,iComp) = avgVal[iComp];
    }

  // loHi (lower right)
  corner = "upperLeft";
  baseIv  = loHi;
  avgVal = cornerExtrap(baseIv       ,
                        offsetDirLoHi,
                        a_extendedFab,
                        a_order      ,
                        a_dx         ,
                        a_blockName  ,
                        corner       );
  
  for (int iComp = 0; iComp < SpaceDim;++iComp)
    {
      a_extendedFab(loHi,iComp) = avgVal[iComp];
    }

  // hiLo (upperLeft)
  baseIv  = hiLo;
  corner = "lowerRight"; 
  avgVal = cornerExtrap(baseIv       ,
                        offsetDirHiLo,
                        a_extendedFab,
                        a_order      ,
                        a_dx         ,
                        a_blockName  ,
                        corner       );
  
  for (int iComp = 0; iComp < SpaceDim;++iComp)
    {
      a_extendedFab(hiLo,iComp) = avgVal[iComp];
    }
}
void extendMapOneSide(FArrayBox            & a_extendedFab        ,
                      const FArrayBox      & a_smallerMap         ,  
                      const Side::LoHiSide & a_sit                ,
                      const bool             a_useExtrapolation[2],
                      const int            & a_extrapDir          ,
                      const RealVect       & a_dx                 ,
                      const int            & a_order              ,
                      const Real           & a_weight             )
{
  IntVect hi = a_smallerMap.bigEnd  ();
  IntVect lo = a_smallerMap.smallEnd();
  
  int otherDir = (a_extrapDir + 1) % 2;
  
  IntVect numPts = hi - lo;
      
  int hiLo = 2*a_sit - 1;
  IntVect baseIv = lo + a_sit*numPts[a_extrapDir]*BASISV(a_extrapDir);
  for (int iComp = 0; iComp < a_smallerMap.nComp(); ++iComp)
    {
      for (int iPt = 0; iPt <= numPts[otherDir]; ++ iPt)
        {
          IntVect offsetIv = BASISV(otherDir)*iPt;
          Vector<Real> funcVal(a_order);
          Vector<Real> argPts (a_order);

          for (int iExtrap = 0; iExtrap < a_order; ++iExtrap)
            {
              IntVect extrapOffset = -hiLo*iExtrap*BASISV(a_extrapDir);
              funcVal[iExtrap] = a_smallerMap(baseIv + offsetIv + extrapOffset,iComp);
              argPts [iExtrap] = a_dx[a_extrapDir]*(baseIv[a_extrapDir] + offsetIv[a_extrapDir] + extrapOffset[a_extrapDir]);
    
#if 0
              //debug
              pout()<< "baseIv       = "<<baseIv      <<endl;
              pout()<< "offsetIv     = "<<offsetIv    <<endl;
              pout()<< "extrapOffset = "<<extrapOffset<<endl;
              pout()<< "baseIv + offsetIv + extrapOffset = "<<baseIv + offsetIv + extrapOffset<<endl;
                      
              pout()<< "funcVal[0] = "<<funcVal[0]<<endl;
              pout()<< "funcVal[1] = "<<funcVal[1]<<endl;
                      
              pout()<< "argPts[0] = "<<argPts[0]<<endl;
              pout()<< "argPts[1] = "<<argPts[1]<<endl;
#endif
            }

          Real extrapPt  =  a_dx[a_extrapDir]*(baseIv[a_extrapDir] + hiLo);
          Real extrapVal =  extrapolate(extrapPt         ,
                                        funcVal          ,
                                        argPts           ,
                                        a_dx[a_extrapDir]);
                  
          a_extendedFab(baseIv + offsetIv + hiLo*BASISV(a_extrapDir),iComp) += a_weight*extrapVal;
        }
    }
}

std::string convertEdgeToString(const int& a_side)
{
  std::string retval = "invalid";
  if (a_side == 0)
    {
      retval = "bottom";
    }
  else if(a_side == 1)
    {
      retval = "right";
    }
  else if(a_side == 2)
    {
      
      retval = "top";
    }
  else if(a_side == 3)
    {
      retval = "left";
    }
  else
    {
      MayDay::Abort("bad value for side");
    }
  return retval;
}

std::string convertEdgeToString(const Side::LoHiSide& a_sit,
                                const int           & a_dir)
{
  std::string retval = "invalid";
  if (a_dir == 0)
    {
      if (a_sit == Side::Lo)
        {
          retval = "left";
        }
      else if (a_sit == Side::Hi)
        {
          retval = "right";
        }
      else
        {
          MayDay::Abort("bad value for side");
        }
    }
  else if (a_dir == 1)
    {
      if (a_sit == Side::Lo)
        {
          retval = "bottom";
        }
      else if (a_sit == Side::Hi)
        {
          retval = "top";
        }
      else
        {
          MayDay::Abort("bad value for side");
        }
    }
  else
    {
      MayDay::Abort("bad value for dir");
    }
  return retval;
}

// respace interior ghost nodes
void reSpaceInteriorGhostNodes(FArrayBox & a_extendedMap,
                               const Box & a_iterBox    ,
                               const int & a_tanDir     )

{
  // find top and bottom iv and linearly interpolate
  IntVect topIv = a_iterBox.bigEnd();
  IntVect botIv = a_iterBox.smallEnd();
  
  int numNodes =  topIv[a_tanDir];  
  for (BoxIterator bit(a_iterBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      Real alpha = (iv[a_tanDir]/(Real)numNodes);
      
      for (int iComp = 0; iComp < a_extendedMap.nComp(); ++iComp)
        {
          a_extendedMap(iv,iComp) = alpha*a_extendedMap(topIv,iComp) + (1.0 - alpha)*a_extendedMap(botIv,iComp);
        }
    }
}
  
void reSpaceInteriorGhostNodes(FArrayBox            & a_extendedMap,
                               const Box            & a_extendedBox,
                               const Box            & a_currBox    ,
                               const int            & a_numGhost   ,
                               const std::string    & a_blockName  ,
                               const Side::LoHiSide & a_sit        ,
                               const int            & a_dir        )
  
{
  bool hasXPoint;
  int xPointEdge;
  int iEdge = LARGEINTVAL;
  if (a_sit == 1 && a_dir == 0)
    {
      iEdge = 1;
    }
  else if (a_sit == 1 && a_dir == 1)
    {
      iEdge = 2;
    }
  else if (a_sit == 0 && a_dir == 0)
    {
      iEdge = 3;
    }
  else
    {
      iEdge = 0;
    }
  
  std::string edge = convertEdgeToString(iEdge); 
  bool xPointHiSide;
  hasXPoint = edgeContainsXPoint(xPointHiSide,
                                 a_blockName ,
                                 edge        );

  if(hasXPoint)
    {
      xPointEdge = iEdge;
      
      Box xPointBox = a_currBox;
      if (xPointEdge == 0)
        {
          xPointBox.growLo(1,a_numGhost);
        }
      else if (xPointEdge == 1)
        {
          xPointBox.growHi(0,a_numGhost);
        }
      else if (xPointEdge == 2)
        {
          xPointBox.growHi(1,a_numGhost);
        }
      else if (xPointEdge == 3)
        {
          xPointBox.growLo(0,a_numGhost);
        }
      
      int xDir = 0;
      int yDir = 1;
      
      IntVect numNodes = a_currBox.size();
      
      // grow currBox in xDir and yDir
      Box xBox = a_currBox;
      xBox.grow(a_numGhost*BASISV(xDir));
      
      Box yBox = a_currBox;
      yBox.grow(a_numGhost*BASISV(yDir));
      
      // box iterate over extended and intersect with grown curr (xDir and yDir)
      for (BoxIterator bit(a_extendedBox); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          bool xTrue = xBox.contains(iv);
          bool yTrue = yBox.contains(iv);
          bool xPointSide = xPointBox.contains(iv);
          bool ghostCellTrue = !a_currBox.contains(iv);
          CH_assert(!(xTrue && yTrue && ghostCellTrue));
          
          int extendDir = 1;
          if (xTrue)
            {
              extendDir = 0;
            }
          
          int tanDir = (extendDir + 1)% SpaceDim;
          if ( (xTrue || yTrue) && (ghostCellTrue) &&(xPointSide) && (xPointEdge == 0 || xPointEdge == 2))
            {
              bool hiSide = false;
              if (iv[extendDir] > numNodes[extendDir] - 1)
                {
                  hiSide = true;
                }
              // find top and bottom iv and linearly interpolate
              IntVect topIv ;
              IntVect botIv;
              
              topIv[extendDir] = iv[extendDir];
              topIv[tanDir   ] = numNodes[tanDir] - 1;
              
              botIv[extendDir] = iv[extendDir];
              botIv[tanDir   ] = 0;
              
              Real alpha = LARGEREALVAL;
              
              alpha = (iv[tanDir]/(Real)numNodes[tanDir]);
              
              for (int iComp = 0; iComp < a_extendedMap.nComp(); ++iComp)
                {
                  a_extendedMap(iv,iComp) = alpha*a_extendedMap(topIv,iComp) + (1.0 - alpha)*a_extendedMap(botIv,iComp);
                }
            }
          
        }
    }
}
void extendMapDir(FArrayBox         & a_extendedMap        ,
                  const bool        & a_cornerCase         ,
                  const int         & a_dir                ,     
                  const bool          a_useExtrapolation[2],
                  const Box         & a_stencilBox         ,
                  const int         & a_numGhost           ,
                  const RealVect    & a_dx                 ,  
                  const int         & a_order              ,
                  const Real        & a_weight             )
{
  for (SideIterator sit; sit.ok(); ++sit)
    {
      if (a_useExtrapolation[sit()])
        {  
          Box validBox = a_stencilBox;
                    
          for (int iGhost = 1; iGhost <= a_numGhost; ++iGhost)
            {
              FArrayBox smallerMap(validBox,SpaceDim);
              smallerMap.copy(a_extendedMap);
              
              // lagrange formula for evenly spaced points applied hi and lo
              extendMapOneSide(a_extendedMap     ,
                               smallerMap        ,
                               sit()             ,
                               a_useExtrapolation,
                               a_dir             ,
                               a_dx              ,
                               a_order           ,
                               a_weight          );
              
              // grow the valid region for the next layer of ghost cells
              int oneCell = 1;
              validBox.growDir(a_dir,sit(),oneCell);
            }
        }
    }
}

void extendMap(FArrayBox         & a_extendedMap        ,       
               const bool          a_useExtrapolation[4],
               const Box         & a_currBox            ,
               const int         & a_numGhost           ,
               const RealVect    & a_dx                 ,  
               const int         & a_order              ,
               const std::string & a_blockName          ,
               const Real        & a_weight             )
{
  // extendedMap agrees with curr on currBox
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      // indexing for a_useExtrapolation
      // bottom = 0
      // right  = 1 
      // top    = 2
      // left   = 3
          
      // indexing for useExtrapolation
      //lo    = 0
      //high  = 1 
      bool useExtrapolation[2];
          
      if (idir == 0)
        {
          useExtrapolation[0] = a_useExtrapolation[3];
          useExtrapolation[1] = a_useExtrapolation[1]; 
        }
      else
        {
          useExtrapolation[0] = a_useExtrapolation[0];
          useExtrapolation[1] = a_useExtrapolation[2]; 
        }

      // not doing corner ghost cells in this call
      bool cornerCase = false;

      extendMapDir(a_extendedMap   ,
                   cornerCase      ,
                   idir            ,
                   useExtrapolation,
                   a_currBox       ,
                   a_numGhost      ,
                   a_dx            ,  
                   a_order         ,
                   a_weight        );
    }
}

void getUnitSquareGrid(const IntVect& a_hiCorner,
                       FArrayBox    & a_xi      )
{
  CH_assert(a_hiCorner[0] != 0);
  CH_assert(a_hiCorner[1] != 0);

  Box box(IntVect::Zero, a_hiCorner);
  RealVect dx;
  dx[0] = 1.0/a_hiCorner[0];
  dx[1] = 1.0/a_hiCorner[1];
  
  a_xi.define(box,SpaceDim);
  
  for (BoxIterator bit(box);bit.ok();++bit) 
    {
      IntVect iv = bit();
      for (int idir = 0; idir<SpaceDim; ++idir) 
        {
          a_xi(iv,idir) = iv[idir]*dx[idir];
        }
    }
}

void computeFirstDerivDir(Real            & a_xDeriv   ,
                          Real            & a_yDeriv   ,
                          const IntVect   & a_iv       ,
                          const int       & a_dir      ,
                          const RealVect  & a_dx       ,
                          const FArrayBox & a_curr     ,
                          const bool      & a_writeInfo)
{
  int xDir = 0;
  int yDir = 1;

  // stencil ivs
  IntVect hi = a_iv;
  IntVect lo = a_iv;
  
  hi += BASISV(a_dir);
  lo -= BASISV(a_dir);

  Real dx = a_dx[a_dir];
 
  a_xDeriv = (a_curr(hi,xDir) - a_curr(lo,xDir))/(2.0*dx);
  a_yDeriv = (a_curr(hi,yDir) - a_curr(lo,yDir))/(2.0*dx);
}

void setBField(FArrayBox         & a_BField          ,
               const FArrayBox   & a_newPhysicalCoord,
               const RealVect    & a_dx              ,
               const std::string & a_geometryFileName,
               const std::string & a_fieldFileName   , 
               const std::string & a_blockName       )
{
  ParmParse pp;
    
  bool simpleTestCase = false;
  pp.query("simpleTestCase",simpleTestCase);

  if(simpleTestCase)
    {
      int xDir = 0;
      int yDir = 1;
     
      Box box = a_BField.box();
      
      Real BPerturb = LARGEREALVAL;
      pp.get("BPerturb",BPerturb);
     
      // x-component
      a_BField.setVal(1.0,xDir);
      
      // y-component
      a_BField.setVal(0.0,yDir);
      for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          // IntVect
          const IntVect& iv = bit();
          a_BField(iv,yDir) += BPerturb*a_dx[yDir]*M_PI*cos(M_PI*iv[xDir]*a_dx[xDir])*sin(M_PI*iv[yDir]*a_dx[yDir]);
        }
    }
  else
    {
      // Construct the field data object
      FieldData field_data(a_geometryFileName,
                           a_fieldFileName   ,
                           a_blockName       );
      
      // normalized B field and four first partials  
      field_data.getFieldUnitVectorFromDCT(a_newPhysicalCoord,
                                           a_BField          );
      Box box = a_BField.box();
      for (BoxIterator bit(box);bit.ok();++bit) 
        {
          IntVect iv = bit();
          Real mag = 0.0;
          for (int idir = 0; idir < SpaceDim; ++idir) 
            {
              mag += a_BField(iv,idir)*a_BField(iv,idir);
            }
          if (mag != 0.0)
            {
              mag = sqrt(mag);
              for (int iComp = 0; iComp < SpaceDim; ++iComp)
                {
                  a_BField(iv,iComp) *= 1/mag;
                }
            }
          else
            {
              pout()<<"BField = (0,0) at iv = "<<iv<<endl;
            }
        }
    }
}


void allocArray(const int& a_rows,
                const int& a_cols,
                Real**   & a_A   )
{
  a_A = new Real* [a_rows];

  for (int i = 0; i < a_rows;i++)
    {
      a_A[i] = new Real [a_cols];
      Real* scanA = a_A[i];

      for (int j = 0; j < a_cols; j++)
        {
          // faster than the more obvious A[i][j] = 0.0;
          *(scanA++) = 0.0;
        }
    }
}

void freeArray(const int& a_rows,
               const int& a_cols,
               Real**   & a_A   )
{
  for (int i = 0; i < a_rows; i++)
    {
      delete[] a_A[i];
    }

  delete[] a_A;
}

void outputArray(const int& a_rows,
                 const int& a_cols,
                 Real**   & a_A   )
{
  for (int i = a_rows -1 ; i >= 0; i--)
    {
      for (int j = 0; j < a_cols; j++)
        {
          pout() << a_A[i][j] << " ";
        }
      pout() << endl;
    }
  pout()<<endl;
}





int getBlockNumber(const std::string & a_sourceBlockName)
{
  int retval = LARGEINTVAL;
  if (a_sourceBlockName == "lcore")
    {
      retval = 0;
    }
  else if (a_sourceBlockName == "rcore")
    {
      retval = 1;
    }
  else if (a_sourceBlockName == "lcsol")
    {
      retval = 2;
    }
  else if (a_sourceBlockName == "rcsol")
    {
      retval = 3;
    }
  else if (a_sourceBlockName == "lsol")
    {
      retval = 4;
    }
  else if (a_sourceBlockName == "rsol")
    {
      retval = 5;
    }
  else if (a_sourceBlockName == "lpf")
    {
      retval = 6;
    }
  else if (a_sourceBlockName == "rpf")
    {
      retval = 7;
    }
  else if (a_sourceBlockName == "mcore")
    {
      retval = 8;
    }
  else if (a_sourceBlockName == "mcsol")
    {
      retval = 9;
    }
  else
    {
      MayDay::Abort("bad block name");
    }
  return retval;
}




void setPhysicalCoordinates (FArrayBox         & a_oldPhysicalCoord ,
                             const IntVect     & a_hiCorner         ,
                             const std::string & a_blockName        , 
                             const std::string & a_geometryFileName ,
                             const std::string & a_fieldFileName    )
{
  // Construct the field data object
  FieldData field_data(a_geometryFileName,
                       a_fieldFileName   ,
                       a_blockName       );

  /*
    Fill an FArrayBox with physical coordinates.  In the
    new mesh generator code using the Brackbill-Saltzmann
    approach, the Euler equation will solved for these
    coordinates.  Here, however, we simply create a grid
    in the unit square, which is the assumed mapped
    domain used internally in the FieldData object and
    get the corresponding physical coordinates.
  */
  
  // Get a grid on the unit square.  The argument sets the size.
  FArrayBox xi;
  getUnitSquareGrid(a_hiCorner, 
                    xi        );
 
  bool simpleTestCase = false;
  ParmParse pp;
  pp.query("simpleTestCase",simpleTestCase);

  CH_assert(xi.box() == a_oldPhysicalCoord.box());
  if (simpleTestCase)
    {
      a_oldPhysicalCoord.copy(xi);
    }
  else
    {
      field_data.getPhysicalCoordinates(xi                , 
                                        a_oldPhysicalCoord);

    }
}

void clearExtendedMappingFile(const std::string& a_outputMap,
                              const int        & a_numRow   )
{
  char suffix[80];
  sprintf(suffix,".%03d.2d.txt",a_numRow + 1);
  
  std::string outputMap = a_outputMap + suffix;
  const char* fileName = outputMap.c_str();
  
  FILE* fd;
  fd = fopen (fileName, "w");
  if (fd == NULL) 
    {
      perror ("Error opening file");
    }
  fclose(fd);
}

void outputMapping(const FArrayBox   & a_extendedMap,
                   const FArrayBox   & a_BField     ,
                   const IntVect     & a_numGhost   ,
                   const std::string & a_blockName  ,
                   const std::string & a_outputMap  )
{
  IntVect boxSize = a_extendedMap.box().size();
  int row = boxSize[0];
  int col = boxSize[1];
  int numRad = row - 2*a_numGhost[0];
  int numPol = col - 2*a_numGhost[1];

  pout()<<"Writing extended map for base resolution = (" << numRad<<","<<numPol<<")"<<endl;
  
  char suffix[80];
  sprintf(suffix,".%03d.2d.txt",numRad);
  
  std::string outputMap = a_outputMap + suffix;
  const char* fileName = outputMap.c_str();

  FILE* fd;
  fd = fopen (fileName, "a");
  if (fd == NULL) 
    {
      perror ("Error opening file");
    }

  fprintf(fd," %s %d %d %d %d\n",a_blockName.c_str(),
          numRad             ,
          a_numGhost[0]      ,
          numPol             ,
          a_numGhost[1]      );

  for (BoxIterator bit(a_extendedMap.box());bit.ok();++bit) 
    {
      IntVect iv = bit();
      for (int iComp = 0; iComp < SpaceDim; ++ iComp)
        {
          fprintf(fd,"%20.13e ",a_extendedMap(iv,iComp));
        }
      for (int iComp = 0; iComp < SpaceDim; ++ iComp)
        {
          fprintf(fd,"%20.13e ",a_BField(iv,iComp)*a_extendedMap(iv,0));
        }
      fprintf(fd,"\n");
    }
  
  fclose(fd);
 }

void getBlockName(std::string& a_blockName,
                  const int  & a_iBlock  )
{
  if (a_iBlock == 0)
    {
      a_blockName = "lcore";
    }
  else if (a_iBlock == 1)
    {
      a_blockName = "rcore";
    }
  else if (a_iBlock == 2)
    {
      a_blockName = "lcsol";
    }
  else if (a_iBlock == 3)
    {
      a_blockName = "rcsol";
    }
  else if (a_iBlock == 4)
    {
      a_blockName = "lsol";
    }
  else if (a_iBlock == 5)
    {
      a_blockName = "rsol";
    }
  else if (a_iBlock == 6)
    {
      a_blockName = "lpf";
    }
  else if (a_iBlock == 7)
    {
      a_blockName = "rpf";
    }
  else if (a_iBlock == 8)
    {
      a_blockName = "mcore";
    }
  else if (a_iBlock == 9)
    {
      a_blockName = "mcsol";
    }
}

Box setSourceGhostBox(const int & a_numGhost        ,
                      const int & a_getGhostFromFace,
                      const Box & a_box     )
{

  int bottom = 0;
  int right  = 1;
  int top    = 2;
  int left   = 3;
  
  int xDir   = 0;
  int yDir   = 1;

  int numCol = a_box.bigEnd()[xDir] - a_box.smallEnd()[xDir];
  int numRow = a_box.bigEnd()[yDir] - a_box.smallEnd()[yDir];

  IntVect ghostWidth  = (a_numGhost - 1)*BASISV(xDir);
  IntVect ghostHeight = (a_numGhost - 1)*BASISV(yDir);
  
  Box retval;

  if (a_getGhostFromFace == bottom)
    {
      IntVect smallEnd = a_box.smallEnd() + BASISV(yDir);
      IntVect bigEnd   =  smallEnd + ghostHeight + numCol*BASISV(xDir);
      retval.define(smallEnd,bigEnd);
    }
  else if (a_getGhostFromFace == right)
    {
      IntVect smallEnd = a_box.smallEnd() + (numCol - ghostWidth)*BASISV(xDir) - BASISV(xDir);
      IntVect bigEnd   = smallEnd + ghostWidth + numRow*BASISV(yDir);
      retval.define(smallEnd,bigEnd);
    }
  else if (a_getGhostFromFace == top)
    {
      IntVect smallEnd = a_box.smallEnd() + (numRow - ghostHeight)*BASISV(yDir) - BASISV(yDir);
      IntVect bigEnd   =  smallEnd + ghostHeight + numCol*BASISV(xDir);
      retval.define(smallEnd,bigEnd);
    }
  else if (a_getGhostFromFace == left)
    {
      IntVect smallEnd = a_box.smallEnd() + BASISV(xDir);
      IntVect bigEnd = smallEnd + ghostWidth + numRow*BASISV(yDir);
      retval.define(smallEnd,bigEnd);
    }
  else
    {
      MayDay::Abort("bad value in getGhostFromFace");
    }

  // box corresponding to ghost cell data
  return retval;

}

Box  setDestGhostBox(const int & a_numGhost        ,
                     const int & a_needGhostForFace,
                     const Box & a_box             )
{
  int bottom = 0;
  int right  = 1;
  int top    = 2;
  int left   = 3;
  
  int xDir   = 0;
  int yDir   = 1;

  int numCol = a_box.bigEnd()[xDir] - a_box.smallEnd()[xDir];
  int numRow = a_box.bigEnd()[yDir] - a_box.smallEnd()[yDir];
  
  IntVect ghostWidth  = (a_numGhost - 1)*BASISV(xDir);
  IntVect ghostHeight = (a_numGhost - 1)*BASISV(yDir);

  Box retval;
 
  if (a_needGhostForFace == bottom)
    {
      IntVect smallEnd = a_box.smallEnd() - a_numGhost*BASISV(yDir);
      IntVect bigEnd   = smallEnd + ghostHeight + numCol*BASISV(xDir);
      retval.define(smallEnd,bigEnd);
    }
  else if (a_needGhostForFace == right)
    {
      IntVect smallEnd = a_box.smallEnd() + (numCol + 1)*BASISV(xDir);
      IntVect bigEnd   = smallEnd + ghostWidth + numRow*BASISV(yDir);
      retval.define(smallEnd,bigEnd);
      
    }
  else if (a_needGhostForFace == top)
    {
      IntVect smallEnd = a_box.smallEnd() + (numRow + 1)*BASISV(yDir);
      IntVect bigEnd   = smallEnd + ghostHeight + numCol*BASISV(xDir);
      retval.define(smallEnd,bigEnd);
    }
  else if (a_needGhostForFace == left)
    {
      IntVect smallEnd = a_box.smallEnd() - a_numGhost*BASISV(xDir);
      IntVect bigEnd   = smallEnd + ghostWidth + numRow*BASISV(yDir);
      retval.define(smallEnd,bigEnd);
    }
  else
    {
      MayDay::Abort("bad value in needGhostForFace");
    }

  // box corresponding to ghost cell data
  return retval;
}
void  matchData(FArrayBox       & a_destFab         ,
                const int       & a_blockNumber     ,
                const int       & a_ghostBlockNumber,
                const FArrayBox & a_ghostCellFab    ,
                const int       & a_needGhostForFace,
                const int       & a_getGhostFromFace)
{
  // for indexing and separating cases
  int bottom = 0;
  int right  = 1;
  int top    = 2;
  int left   = 3;
  
  //int xDir   = 0;
  //int yDir   = 1;

  // for a call to copy
  Interval comps(0,1);
  
  if (a_needGhostForFace == bottom)
    {
      if (a_getGhostFromFace == bottom)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else if (a_getGhostFromFace == right)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else if (a_getGhostFromFace == top)
        {
          a_destFab.copy(a_ghostCellFab.box(),
                         comps               ,
                         a_destFab.box()     ,
                         a_ghostCellFab      ,
                         comps               );
        
        }
      else if (a_getGhostFromFace == left)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else
        {
          MayDay::Abort("bad value for needGhostFromFace");
        }
    }
  else if (a_needGhostForFace == right)
    {
      if (a_getGhostFromFace == bottom)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else if (a_getGhostFromFace == right)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else if (a_getGhostFromFace == top)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else if (a_getGhostFromFace == left)
        {
          a_destFab.copy(a_ghostCellFab.box(),
                         comps               ,
                         a_destFab.box()     ,
                         a_ghostCellFab      ,
                         comps               );

        }
      else
        {
          MayDay::Abort("bad value for needGhostFromFace");
        }
    }
  else if (a_needGhostForFace == top)
    {
      if (a_getGhostFromFace == bottom)
        {
          a_destFab.copy(a_ghostCellFab.box(),
                         comps               ,
                         a_destFab.box()     ,
                         a_ghostCellFab      ,
                         comps               );     
        }
      else if (a_getGhostFromFace == right)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else if (a_getGhostFromFace == top)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else if (a_getGhostFromFace == left)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else
        {
          MayDay::Abort("bad value for needGhostFromFace");
        }
    }
  else if (a_needGhostForFace == left)
    {
      if (a_getGhostFromFace == bottom)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else if (a_getGhostFromFace == right)
        {
          a_destFab.copy(a_ghostCellFab.box(),
                         comps               ,
                         a_destFab.box()     ,
                         a_ghostCellFab      ,
                         comps               );     
        }
      else if (a_getGhostFromFace == top)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else if (a_getGhostFromFace == left)
        {
          MayDay::Abort("This case never happens in our geometry");
        }
      else
        {
          MayDay::Abort("bad value for needGhostFromFace");
        }
    }  
  
  else
    {
      MayDay::Abort("bad value for needGhostForFace");
    }
}

void copyGhostCell(FArrayBox              & a_extendedMapFab     ,
                   bool                     a_useExtrapolation[4],
                   const Vector<FArrayBox*> a_extendedMaps       , 
                   const std::string      & a_blockName          ,
                   const int              & a_blockNumber        ,
                   const std::string      & a_geometryFileName   ,
                   const std::string      & a_fieldFileName      ,
                   const Box              & a_currBox            ,
                   const int              & a_numComp            ,
                   const int              & a_numGhost           )
{
  for (int needGhostForFace = 0;needGhostForFace < 4;++needGhostForFace)
    {
      //call ghostCellPlan for this block
      int needGhostForBlock = a_blockNumber;
      int getGhostFromBlock = LARGEINTVAL;
      int getGhostFromFace  = LARGEINTVAL;

      ghostCellPlan(getGhostFromBlock                   ,
                    getGhostFromFace                    , 
                    a_useExtrapolation[needGhostForFace],
                    needGhostForBlock                   ,
                    needGhostForFace                    );
      
      std::string ghostBlockName;
      getBlockName(ghostBlockName   ,
                   getGhostFromBlock);
      
      std::string edge = convertEdgeToString(needGhostForFace); 
      bool xPointHiSide;
      bool xPointEdge = edgeContainsXPoint(xPointHiSide,
                                           a_blockName ,
                                           edge        );

      // copy on the entire xPoint edge, then overwrite through extrapolation within the radius
      if (a_useExtrapolation[needGhostForFace] == false || xPointEdge == true)
        {
          //  overwrite a portion later for xPoint edge
          Box sourceGhostBox = setSourceGhostBox(a_numGhost      ,
                                                 getGhostFromFace,
                                                 a_currBox       );
          
          // overwrite a portion later for xPoint edge
          Box destGhostBox = setDestGhostBox(a_numGhost      ,
                                             needGhostForFace,
                                             a_currBox       );
          
          FArrayBox sourceFab(sourceGhostBox,a_numComp);
          sourceFab.copy(*(a_extendedMaps[getGhostFromBlock]));
          
          FArrayBox destFab(destGhostBox,a_numComp);
          matchData(destFab          ,
                    needGhostForBlock,
                    getGhostFromBlock,
                    sourceFab        ,
                    needGhostForFace ,
                    getGhostFromFace);
          
          a_extendedMapFab.copy(destFab);
        }
    }
}

Box getCornerRegion(const int         & a_numGhost , 
                    const std::string & a_corner   ,
                    const Box         & a_currBox  )
{
  Box retval;

  int xDir = 0;
  int yDir = 1;
  
  int numRow = a_currBox.size()[yDir] - 1;
  int numCol = a_currBox.size()[xDir] - 1;

  // spanning vectors of the valid cell box
  IntVect xVec = numCol*BASISV(xDir);
  IntVect yVec = numRow*BASISV(yDir);

  // corners of the box containing valid cells
  IntVect lowerLeft  = a_currBox.smallEnd();
  IntVect lowerRight = lowerLeft + xVec   ;
  IntVect upperRight = a_currBox.bigEnd() ;
  IntVect upperLeft  = lowerLeft + yVec   ;

  // spanning vectors of the ghost box
  IntVect xGhost = a_numGhost*BASISV(xDir);
  IntVect yGhost = a_numGhost*BASISV(yDir);

  // parameters of the retval box
  IntVect lo;
  IntVect hi;

  if (a_corner == "lowerLeft")
    {
      lo = lowerLeft - xGhost - yGhost;
    }
  else if (a_corner == "lowerRight")
    {
      lo = lowerRight + BASISV(xDir) - yGhost;
    }
  else if (a_corner == "upperRight")
    {
      lo = upperRight + IntVect::Unit;
    }
  else if (a_corner == "upperLeft")
    {
      lo = upperLeft - xGhost + BASISV(yDir);
    }
  else 
    {
      MayDay::Abort("bad string value for corner");
    }

  hi = lo + xGhost + yGhost - IntVect::Unit;
  
  // return box in the destination extended box
  retval.define(lo,hi);
  return retval;
}

bool getSourceRegion(std::string       & a_srcBlockName   ,
                     Box               & a_srcBox         ,
                     const Box         & a_destBox        ,
                     const std::string & a_destBlockName  ,
                     const int         & a_numGhost       ,  
                     const std::string & a_corner         ,
                     const Box         & a_currBox        ,
                     const bool        & a_blending          )
{
  bool extrapolate = false;
  
  int xDir = 0;
  int yDir = 1;
  
  int numRow = a_currBox.size()[yDir];
  int numCol = a_currBox.size()[xDir];
  
  // modified spanning vectors of the valid cell box
  IntVect xVec = (numCol - 1)*BASISV(xDir);
  IntVect yVec = (numRow - 1)*BASISV(yDir);
  
  // corners of the box containing valid cells
  IntVect lowerLeft  = a_currBox.smallEnd();
  IntVect lowerRight = lowerLeft + xVec   ;
  IntVect upperRight = a_currBox.bigEnd() ;
  IntVect upperLeft  = lowerLeft + yVec   ;
  
  // spanning vectors of the ghost box
  IntVect xGhost = a_numGhost*BASISV(xDir);
  IntVect yGhost = a_numGhost*BASISV(yDir);
  
  // parameters of the retval box
  IntVect lo = LARGEINTVAL*IntVect::Unit;
  IntVect hi = LARGEINTVAL*IntVect::Unit;
  
  IntVect loDest = a_destBox.smallEnd();
  IntVect hiDest = a_destBox.bigEnd  ();
  
  IntVect sourceDir = LARGEINTVAL*IntVect::Unit;
  
  // all blocks share of row or column of values, which creates an offset in the calculations below
  if (a_destBlockName == "lcore")
    {
      if (a_corner == "lowerLeft")
        {
          // lo depends on srcBlock and a_corner. In particular, "upperLeft" below depends on srcBlock as well as a_corner.
          lo             = upperLeft - xGhost - yGhost;
          a_srcBlockName = "mcore";
          extrapolate    = true   ;
        }
      else if (a_corner == "lowerRight")
        {
          lo             = upperRight + BASISV(xDir) - yGhost;
          a_srcBlockName = "mcore";
          extrapolate    = true   ;
        }
      else if (a_corner == "upperRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else 
        {
          MayDay::Abort("bad string value for corner");
        }
      
      // true for all cases, but largeIntVal if extrapolate = true
      hi = lo + xGhost + yGhost - IntVect::Unit;
    }
  else  if (a_destBlockName == "rcore")
    {
      if (a_corner == "lowerLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "lowerRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperRight")
        {
          lo              = upperLeft + BASISV(xDir) + BASISV(yDir);                       
          a_srcBlockName  = "rcsol";
          extrapolate     = true;
        }
      else if (a_corner == "upperLeft")
        {
          lo              = lowerLeft - xGhost + BASISV(yDir);
          a_srcBlockName  = "mcore";
          extrapolate     = true;
        }
      else 
        {
          MayDay::Abort("bad string value for corner");
        }
      
      // true for all cases, but largeIntVal if extrapolate = true
      hi = lo + xGhost + yGhost - IntVect::Unit;
    }
  else if (a_destBlockName == "lcsol")
    {
      if (a_corner == "lowerLeft")
        {
          lo              = lowerRight - xGhost - yGhost;
          a_srcBlockName  = "mcore";
          extrapolate     = true;
        }
      else if (a_corner == "lowerRight")
        {
          lo              = upperRight + BASISV(xDir) - yGhost;
          a_srcBlockName  = "mcore";
          extrapolate     = true;
        }
      else if (a_corner == "upperRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else 
        {
          MayDay::Abort("bad string value for corner");
        }
      
      // true for all cases, but largeIntVal if extrapolate = true
      hi = lo + xGhost + yGhost - IntVect::Unit;
    }
  else if (a_destBlockName == "rcsol")
    {
      if (a_corner == "lowerLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "lowerRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperRight")
        {
          lo              = lowerRight + BASISV(xDir) + BASISV(yDir);
          a_srcBlockName  = "mcsol";
          extrapolate     = true;
        }
      else if (a_corner == "upperLeft")
        {
          lo              = lowerLeft - xGhost + BASISV(yDir);
          a_srcBlockName  = "mcsol";
          extrapolate     = true;
        }
      else 
        {
          MayDay::Abort("bad string value for corner");
        }
      
      // true for all cases, but largeIntVal if extrapolate = true
      hi = lo + xGhost + yGhost - IntVect::Unit;
    }
  else if (a_destBlockName == "lsol")
    {
      if (a_corner == "lowerLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "lowerRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else 
        {
          MayDay::Abort("bad string value for corner");
        }
      
      // true for all cases, but largeIntVal if extrapolate = true
      hi = lo + xGhost + yGhost - IntVect::Unit;
    }
  else if (a_destBlockName == "rsol")
    {
      if (a_corner == "lowerLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "lowerRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else 
        {
          MayDay::Abort("bad string value for corner");
        }

      // true for all cases, but largeIntVal if extrapolate = true
      hi = lo + xGhost + yGhost - IntVect::Unit;
    }
  else if (a_destBlockName == "lpf")
    {
      if (a_corner == "lowerLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "lowerRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else 
        {
          MayDay::Abort("bad string value for corner");
        }
      
      // true for all cases, but largeIntVal if extrapolate = true
      hi = lo + xGhost + yGhost - IntVect::Unit;
    }
  else if (a_destBlockName == "rpf")
    {
      if (a_corner == "lowerLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "lowerRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperLeft")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else 
        {
          MayDay::Abort("bad string value for corner");
        }
      
      // true for all cases, but largeIntVal if extrapolate = true
      hi = lo + xGhost + yGhost - IntVect::Unit;
    }
  else if (a_destBlockName == "mcore")
    {
      if (a_corner == "lowerLeft")
        {
          lo             = upperLeft - xGhost - yGhost;
          a_srcBlockName = "rcore";
          extrapolate       = true; 
        }
      else if (a_corner == "lowerRight")
        {
          lo             = upperRight + BASISV(xDir) - yGhost;
          a_srcBlockName = "rcore";
          extrapolate       = true;
        }
      else if (a_corner == "upperRight")
        {
          lo              = upperLeft + BASISV(xDir) + BASISV(yDir);                       
          a_srcBlockName = "mcsol";
          extrapolate   = true;
        }
      else if (a_corner == "upperLeft")
        {
          lo              = lowerLeft - xGhost + BASISV(yDir);
          a_srcBlockName = "lcore";
          extrapolate       = true;
        }
      else 
        {
          MayDay::Abort("bad string value for corner");
        }

      // true for all cases, but largeIntVal if extrapolate = true
      hi = lo + xGhost + yGhost - IntVect::Unit;
    }
  else if (a_destBlockName == "mcsol")
    {
      if (a_corner == "lowerLeft")
        {
          lo              = lowerRight - xGhost - yGhost;
          a_srcBlockName  = "mcore";
          extrapolate    = true;
        }
      else if (a_corner == "lowerRight")
        {
          lo             = upperRight + BASISV(xDir) - yGhost;
          a_srcBlockName = "rcsol";
          extrapolate       = true;
        }
      else if (a_corner == "upperRight")
        {
          lo             = lowerRight + BASISV(xDir) + BASISV(yDir);
          a_srcBlockName = "lcsol";
          extrapolate       = true;
        }
      else if (a_corner == "upperLeft")
        {
          lo             = lowerLeft - xGhost + BASISV(yDir);
          a_srcBlockName = "lcsol";
          extrapolate       = true;
        }
      else 
        {
          MayDay::Abort("bad string value for corner");
        }

      // true for all cases, but largeIntVal if extrapolate = true
      hi = lo + xGhost + yGhost - IntVect::Unit;
    }
  else
    {
      MayDay::Abort("bad string value for blockname");
    }
  
  // return box in the source extended box
  a_srcBox.define(lo,hi);
  
  return extrapolate;
}

void cornerGhostCells(Vector<FArrayBox*>        a_extendedMap   ,
                      const std::string       & a_corner        ,
                      const Vector<FArrayBox*>  a_curr          ,
                      const int               & a_numGhost      ,
                      const std::string       & a_blockName     ,
                      const RealVect          & a_dx            ,
                      const int               & a_order         ,
                      const bool                a_useHiLoXDir[2],
                      const bool                a_useHiLoYDir[2],
                      const int               & a_blockNumber   ,
                      const bool              & a_blending      )
{
  bool cornerCase = true;

  std::string sourceBlockName;
  Box sourceBox;
  Interval comps(0,1);
  Box cornerBox = getCornerRegion(a_numGhost                  ,
                                  a_corner                    ,
                                  a_curr[a_blockNumber]->box());

  //get the source region or learn that this region will be extrapolated using data in the extrapDirection
  bool extrap = getSourceRegion(sourceBlockName             ,
                                sourceBox                   ,
                                cornerBox                   ,
                                a_blockName                 ,
                                a_numGhost                  , 
                                a_corner                    ,
                                a_curr[a_blockNumber]->box(),
                                a_blending                     );

  if (extrap == false)
    {
      // block number of data that will be copied
      int sourceBlock = getBlockNumber(sourceBlockName);

      // fill in the ghost cells 
      a_extendedMap[a_blockNumber]->copy(sourceBox                   ,
                                         comps                       ,
                                         cornerBox                   ,
                                         *(a_extendedMap[sourceBlock]),
                                         comps                       );
    }
  else
    {
    
      int xDir = 0;
      int yDir = 1;
                    
      // extrapolate in each direction and average
      Real weight = 1.0;
          
      // extrapolate in y-direction
      Box smallerBoxXShift;
      if (a_corner == "lowerLeft" || a_corner == "upperLeft")
        {
          smallerBoxXShift = adjCellLo(a_curr[a_blockNumber]->box(),xDir,a_numGhost);
        }
      else if (a_corner == "lowerRight" || a_corner == "upperRight")
        {
          smallerBoxXShift = adjCellHi(a_curr[a_blockNumber]->box(),xDir,a_numGhost);
        }
      else
        {
          MayDay::Abort("bad corner name");
        }
          
      // use yDir appropriate flags (in particular, shift the box in xDir)
      extendMapDir(*(a_extendedMap  [a_blockNumber]),
                   cornerCase                       ,
                   yDir                             ,
                   a_useHiLoYDir                    ,
                   smallerBoxXShift                 ,
                   a_numGhost                       ,
                   a_dx                             ,
                   a_order                          ,
                   weight                           );

      // average the extrapolations in the xDir and yDir directions
      for (BoxIterator bit(cornerBox);bit.ok();++bit) 
        {
          IntVect iv = bit();
          for (int iComp = 0; iComp < a_extendedMap[a_blockNumber]->nComp(); ++iComp)
            {
              (*(a_extendedMap [a_blockNumber]))(iv,iComp) *= 0.5;
            }
        }
          
      // extrapolate in x-direction
      Box smallerBoxYShift;
      if (a_corner == "lowerLeft" || a_corner == "lowerRight")
        {
          smallerBoxYShift = adjCellLo(a_curr[a_blockNumber]->box(),yDir,a_numGhost);
        }
      else if (a_corner == "upperLeft" || a_corner == "upperRight")
        {
          smallerBoxYShift = adjCellHi(a_curr[a_blockNumber]->box(),yDir,a_numGhost);
        }
          
      weight = 0.5;
          
      // use xDir appropriate flags (in particular, shift the box in yDir)
      extendMapDir(*(a_extendedMap [a_blockNumber]),
                   cornerCase                      ,
                   xDir                            ,
                   a_useHiLoXDir                   ,
                   smallerBoxYShift                ,
                   a_numGhost                      ,
                   a_dx                            ,
                   a_order                         ,
                   weight                          );
    }
}

void blend(Vector<RealVect>       & a_blend    ,
           const Vector<RealVect> & a_extrap   ,
           const Vector<RealVect> & a_copy     ,
           const Real             & a_startFrac,
           const Real             & a_endFrac  )
{
  int length  = a_extrap.size();
  int start   = floor(a_startFrac*length);
  int end     = floor(a_endFrac*length);
  
  Real normalize = end - start;
  for (int idir = 0; idir < SpaceDim ;++idir)
    {
      for (int ival = 0; ival < length; ++ival)
        {
          if (ival < start)
            {
              a_blend[ival][idir] = a_extrap[ival][idir];
            }
          else if (ival >= start && ival <= end)
            {
              Real alpha = 1.0 - (ival - start)/normalize;
              a_blend[ival][idir] = alpha*a_extrap[ival][idir] + (1.0 - alpha)*a_copy[ival][idir];
            }
          else if (ival > end)
            {
              a_blend[ival][idir] = a_copy[ival][idir];
            }
        }
    }
}

void blend(FArrayBox         & a_extendedMap,
           const Box         & a_currBox    ,
           const FArrayBox   & a_copyTemp   ,
           const FArrayBox   & a_extrapTemp ,
           const std::string & a_blockName  ,
           const int         & a_iBlock     ,
           const int         & a_numGhost   ,
           const Real        & a_startFrac ,
           const Real        & a_endFrac  )
{
  int nComp = a_extendedMap.nComp();

  for (int iFace = 0; iFace < 4; ++iFace)
    {
      int dir = 0;
      int otherDir = 1;
      Side::LoHiSide side = Side::Lo;
      
      if (iFace == 0 || iFace == 2)
        {
          dir      = 1;
          otherDir = 0;
        }
      
      if (iFace == 1 || iFace == 2)
        {
          side = Side::Hi;
        }
      
      //call ghostCellPlan for this block
      int sourceBlock       = LARGEINTVAL;
      int sourceFace        = LARGEINTVAL;
      int sourceBlockNumber = LARGEINTVAL;
      bool extrapolation;
      ghostCellPlan(sourceBlock  ,
                    sourceFace   , 
                    extrapolation,
                    a_iBlock     ,
                    iFace        );
      
      std::string sourceBlockName;
      getBlockName(sourceBlockName   ,
                   sourceBlockNumber);
      
      std::string edge = convertEdgeToString(iFace); 
      bool xPointHiSide;
      bool xPointEdge = edgeContainsXPoint(xPointHiSide,
                                           a_blockName ,
                                           edge        );
      if (xPointEdge)
        {
          int length = a_currBox.size(otherDir);
          
          for (int iGhost = 0; iGhost < a_numGhost; ++iGhost) 
            {
              //blend
              Vector<RealVect> extrapVec(length);
              Vector<RealVect> copyVec  (length);
              Vector<RealVect> blendVec (length);

              int iStart;
              int iEnd;
              int iInc;

              int jStart;
              int jEnd;
              int jInc;
             
              if (dir == 0)
                {
                  if (xPointHiSide)
                    {
                      iStart  = a_currBox.smallEnd()[dir] + sign(side)*(1 + iGhost);
                      iStart += (1 + sign(side))/2 *( a_currBox.bigEnd()[dir] -  a_currBox.smallEnd()[dir]);
                      iEnd    = iStart;
                      iInc    = 1;
                      
                      jStart = a_currBox.bigEnd  ()[otherDir];
                      jEnd   = a_currBox.smallEnd()[otherDir];
                      jInc   = -1;
                    }
                  else
                    {
                      iStart  = a_currBox.smallEnd()[dir] + sign(side)*(1 + iGhost);
                      iStart += (1 + sign(side))/2 *( a_currBox.bigEnd()[dir] -  a_currBox.smallEnd()[dir]);
                      iEnd    = iStart;
                      iInc    = 1;
                      
                      jStart = a_currBox.smallEnd()[otherDir];
                      jEnd   = a_currBox.bigEnd  ()[otherDir];
                      jInc   = 1;
                    }
                }
              else
                {
                  if (xPointHiSide)
                    {
                      jStart  = a_currBox.smallEnd()[dir] + sign(side)*(1 + iGhost);
                      jStart += (1 + sign(side))/2 *( a_currBox.bigEnd()[dir] -  a_currBox.smallEnd()[dir]);
                      jEnd    = jStart;
                      jInc    = 1;
                      
                      iStart = a_currBox.bigEnd  ()[otherDir];
                      iEnd   = a_currBox.smallEnd()[otherDir];
                      iInc   = -1;
                    }
                  else
                    {
                      jStart  = a_currBox.smallEnd()[dir] + sign(side)*(1 + iGhost);
                      jStart += (1 + sign(side))/2 *( a_currBox.bigEnd()[dir] -  a_currBox.smallEnd()[dir]);
                      jEnd    = jStart;
                      jInc    = 1;
                      
                      iStart = a_currBox.smallEnd()[otherDir];
                      iEnd   = a_currBox.bigEnd  ()[otherDir];
                      iInc   = 1;
                    }
                }

              int iVec = 0; 
              for (int i = iStart; i*iInc <= iEnd*iInc; i += iInc)
                {
                  for (int j = jStart; j*jInc <= jEnd*jInc; j += jInc)
                    {
                      IntVect iv(i,j);
                      for (int iComp = 0; iComp < nComp; ++iComp)
                        {
                          extrapVec[iVec][iComp] = a_extrapTemp(iv,iComp);
                          copyVec  [iVec][iComp] = a_copyTemp  (iv,iComp);
                        }
                      
                      iVec += 1;
                    }
                }
              
              blend(blendVec   ,
                    extrapVec  ,
                    copyVec    ,
                    a_startFrac,
                    a_endFrac  );
         
              // fill extendedMap
              iVec = 0; 
              for (int i = iStart; i*iInc <= iEnd*iInc; i += iInc)
                {
                  for (int j = jStart; j*jInc <= jEnd*jInc; j += jInc)
                    {
                      IntVect iv(i,j);
                      for (int iComp = 0; iComp < nComp; ++iComp)
                        {
                          a_extendedMap(iv,iComp) = blendVec[iVec][iComp];
                        }
                      
                      iVec += 1;
                    }
                }
              
            }
        }
      else 
        {
          // box containing ghost cells
          Box validBox = adjCellBox(a_currBox , 
                                    dir       , 
                                    side      ,
                                    a_numGhost);
          
          // for copy call
          Interval nVar(0,nComp - 1);
          
          if (extrapolation)
            {
              // copy extrapTemp
              a_extendedMap.copy(validBox    ,
                                 nVar        ,
                                 validBox    ,
                                 a_extrapTemp,
                                 nVar        );
            }
          else 
            {
              //copy copyTemp
              a_extendedMap.copy(validBox  ,
                                 nVar      ,
                                 validBox  ,
                                 a_copyTemp,
                                 nVar      );
            }
        }
    }
}
void computeJacobianPartials(FArrayBox       & a_Jacobian   ,
                             const FArrayBox & a_mapPartials,
                             const bool      & a_writeInfo  )
{
  // iindexing macros
  int xDir   = 0;
  int yDir   = 1;
  int xiDir  = 0;
  int etaDir = 1;
  
  Box box = a_Jacobian.box();
  
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      // IntVect
      const IntVect& iv = bit();
    
      Real xXi    = a_mapPartials(iv,xDir + 0*SpaceDim); //0
      Real yXi    = a_mapPartials(iv,yDir + 0*SpaceDim); //1
      
      Real xEta   = a_mapPartials(iv,xDir + 1*SpaceDim); //2
      Real yEta   = a_mapPartials(iv,yDir + 1*SpaceDim); //3
      
      Real x2Xi   = a_mapPartials(iv,xDir + 2*SpaceDim); //4
      Real y2Xi   = a_mapPartials(iv,yDir + 2*SpaceDim); //5
      
      Real x2Eta  = a_mapPartials(iv,xDir + 3*SpaceDim); //6
      Real y2Eta  = a_mapPartials(iv,yDir + 3*SpaceDim); //7
      
      Real xXiEta = a_mapPartials(iv,xDir + 4*SpaceDim); //8
      Real yXiEta = a_mapPartials(iv,yDir + 4*SpaceDim); //9
      
      // Jacobian
      Real J = xXi*yEta - xEta*yXi;
      
      // derivatives with respect to derivatives of the mapping
      Real JxXi  =  yEta;
      Real JyXi  = -xEta;
      Real JxEta = -yXi;
      Real JyEta =  xXi;
      
      Real JxXiXi   =  yXiEta;
      Real JyXiXi   = -xXiEta;
      Real JxEtaEta = -yXiEta;
      Real JyEtaEta =  xXiEta;
      
      // derivatives with respect to xi and eta
      Real Jxi  = x2Xi  *yEta + xXi*yXiEta - xXiEta*yXi - xEta*y2Xi  ;
      Real Jeta = xXiEta*yEta + xXi*y2Eta  - x2Eta *yXi - xEta*yXiEta;
      
      //J
      a_Jacobian(iv,0                      ) = J       ; //0
      
      //JxXi,JyXi,JxEta,JyEta
      a_Jacobian(iv,1 + 0*SpaceDim + xDir  ) = JxXi    ; //1
      a_Jacobian(iv,1 + 0*SpaceDim + yDir  ) = JyXi    ; //2
      a_Jacobian(iv,1 + 1*SpaceDim + xDir  ) = JxEta   ; //3
      a_Jacobian(iv,1 + 1*SpaceDim + yDir  ) = JyEta   ; //4
      
      //JxXiXi,JyXiXi,JxEtaEta,JyEtaEta
      a_Jacobian(iv,1 + 2*SpaceDim + xDir  ) = JxXiXi  ; //5
      a_Jacobian(iv,1 + 2*SpaceDim + yDir  ) = JyXiXi  ; //6
      a_Jacobian(iv,1 + 3*SpaceDim + xDir  ) = JxEtaEta; //7
      a_Jacobian(iv,1 + 3*SpaceDim + yDir  ) = JyEtaEta; //8
      
      //Jxi,Jeta
      a_Jacobian(iv,1 + 4*SpaceDim + xiDir ) = Jxi     ; //9
      a_Jacobian(iv,1 + 4*SpaceDim + etaDir) = Jeta    ; //10
      
      
      // smaller than zero jacobian => we're searching outside the space of 1-1, orientation preserving maps.
      if (J <= 0.0)
        {
          if(a_writeInfo)
            {
              pout()<<"iv = "      <<iv<<endl;
              pout()<<"Jacobian = "<<J <<endl; 
            }
          
          //MayDay::Abort("Jacobian <= 0.0");    
        }
    }
}
void computeSmoothnessPartials(FArrayBox       & a_smoothness ,
                               const FArrayBox & a_mapPartials,
                               const bool      & a_writeInfo  )
{
  // indexing macros
  int xDir = 0;
  int yDir = 1;

  int xiDir  = 0;
  int etaDir = 1;

  int xXiDir  = 0;
  int yXiDir  = 1;
  
  int xEtaDir = 0;
  int yEtaDir = 1;
  
  int xiSmooth  = 0;
  int etaSmooth = 1;

  Box box = a_smoothness.box();
  Box intBox = box;
  intBox.grow(-1);
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      // IntVect
      const IntVect& iv = bit();
      
      // set domain boundary cells to zero for the color map
      if (!intBox.contains(iv))
        {
          // xi smoothing 
          a_smoothness(iv,    7*xiSmooth                       ) = 0.0; //0
          a_smoothness(iv,1 + 7*xiSmooth  +   xXiDir           ) = 0.0; //1
          a_smoothness(iv,1 + 7*xiSmooth  +   yXiDir           ) = 0.0; //2
          a_smoothness(iv,1 + 7*xiSmooth  +   SpaceDim + xXiDir) = 0.0; //3
          a_smoothness(iv,1 + 7*xiSmooth  +   SpaceDim + yXiDir) = 0.0; //4
          a_smoothness(iv,1 + 7*xiSmooth  + 2*SpaceDim + xiDir ) = 0.0; //5
          a_smoothness(iv,1 + 7*xiSmooth  + 2*SpaceDim + etaDir) = 0.0; //6

          // eta smoothing
          a_smoothness(iv,    7*etaSmooth                       ) = 0.0; //7
          a_smoothness(iv,1 + 7*etaSmooth +   xEtaDir           ) = 0.0; //8
          a_smoothness(iv,1 + 7*etaSmooth +   yEtaDir           ) = 0.0; //9
          a_smoothness(iv,1 + 7*etaSmooth +   SpaceDim + xEtaDir) = 0.0; //10 
          a_smoothness(iv,1 + 7*etaSmooth +   SpaceDim + yEtaDir) = 0.0; //11
          a_smoothness(iv,1 + 7*etaSmooth + 2*SpaceDim + xiDir  ) = 0.0; //12 
          a_smoothness(iv,1 + 7*etaSmooth + 2*SpaceDim + etaDir ) = 0.0; //13 
        }
      else
        {
          //   SmoothnessMeasure = xXi *xXi  + yXi*yXi   (first  six entries are S and various derivatives) 
          //or SmoothnessMeasure = xEta*xEta + yEta*yEta (second six entries are S and various derivatives)

          Real xXi    = a_mapPartials(iv,xDir + 0*SpaceDim);
          Real yXi    = a_mapPartials(iv,yDir + 0*SpaceDim);
          
          Real xEta   = a_mapPartials(iv,xDir + 1*SpaceDim);
          Real yEta   = a_mapPartials(iv,yDir + 1*SpaceDim);
          
          Real x2Xi   = a_mapPartials(iv,xDir + 2*SpaceDim);
          Real y2Xi   = a_mapPartials(iv,yDir + 2*SpaceDim);
          
          Real x2Eta  = a_mapPartials(iv,xDir + 3*SpaceDim);
          Real y2Eta  = a_mapPartials(iv,yDir + 3*SpaceDim);
          
          Real xXiEta = a_mapPartials(iv,xDir + 4*SpaceDim);
          Real yXiEta = a_mapPartials(iv,yDir + 4*SpaceDim);

          // xi smoothing 
          Real xiS      = xXi*xXi + yXi*yXi            ;
          Real xiSxXi   = 2.0*xXi                      ;
          Real xiSyXi   = 2.0*yXi                      ;
          Real xiSxXiXi = 2.0*x2Xi                     ;
          Real xiSyXiXi = 2.0*y2Xi                     ;
          Real xiSxi    = 2.0*(xXi*x2Xi   + yXi*y2Xi  );
          Real xiSeta   = 2.0*(xXi*xXiEta + yXi*yXiEta);

          // eta Smoothing
          Real etaS        = xEta*xEta + yEta*yEta          ;
          Real etaSxEta    = 2.0*xEta                       ;
          Real etaSyEta    = 2.0*yEta                       ;
          Real etaSxEtaEta = 2.0*x2Eta                      ;
          Real etaSyEtaEta = 2.0*y2Eta                      ;
          Real etaSxi      = 2.0*(xEta*xXiEta + yEta*yXiEta);
          Real etaSeta     = 2.0*(xEta*x2Eta  + yEta*y2Eta );

          // fill data
          a_smoothness(iv,     7*xiSmooth                       ) = xiS       ;//0
          a_smoothness(iv, 1 + 7*xiSmooth  +   xXiDir           ) = xiSxXi    ;//1
          a_smoothness(iv, 1 + 7*xiSmooth  +   yXiDir           ) = xiSyXi    ;//2
          a_smoothness(iv, 1 + 7*xiSmooth  +   SpaceDim + xXiDir) = xiSxXiXi  ;//3
          a_smoothness(iv, 1 + 7*xiSmooth  +   SpaceDim + yXiDir) = xiSyXiXi  ;//4
          a_smoothness(iv, 1 + 7*xiSmooth  + 2*SpaceDim + xiDir ) = xiSxi     ;//5
          a_smoothness(iv, 1 + 7*xiSmooth  + 2*SpaceDim + etaDir) = xiSeta    ;//6
          
          // eta smoothing
          a_smoothness(iv,     7*etaSmooth                       ) = etaS       ;//7
          a_smoothness(iv, 1 + 7*etaSmooth +   xEtaDir           ) = etaSxEta   ;//8
          a_smoothness(iv, 1 + 7*etaSmooth +   yEtaDir           ) = etaSyEta   ;//9
          a_smoothness(iv, 1 + 7*etaSmooth +   SpaceDim + xEtaDir) = etaSxEtaEta;//10
          a_smoothness(iv, 1 + 7*etaSmooth +   SpaceDim + yEtaDir) = etaSyEtaEta;//11
          a_smoothness(iv, 1 + 7*etaSmooth + 2*SpaceDim + xiDir  ) = etaSxi     ;//12
          a_smoothness(iv, 1 + 7*etaSmooth + 2*SpaceDim + etaDir ) = etaSeta    ;//13
        }
    }
}

void computeSmJnDiagonalTerm(FArrayBox       & a_smoothDiag ,
                             const RealVect  & a_dx         ,
                             const FArrayBox & a_Jacobian   ,
                             const FArrayBox & a_mapPartials)
{
  // indexing macros
  int firstEq = 0;
  int secEq   = 1;

  //int xiEtaComp        = 0;
  //int smoothVolVarComp = 1;

  // indexing macros
  int xDir = 0;
  int yDir = 1;
  
  //int xiDir  = 0;
  //int etaDir = 1;
  
  int m = LARGEREALVAL;
  int n = LARGEREALVAL;
  
  ParmParse pp;
  pp.get("Spower",m);
  pp.get("Jpower",n);
  
  // iterate
  Box box = a_smoothDiag.box();
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      const IntVect iv = bit(); 

      Real alpha = 1.0;
      Real beta  = 1.0;
      
      Real dx = a_dx[xDir];
      Real dy = a_dx[yDir];
      
      Real dx2 = dx*dx;
      Real dy2 = dy*dy;
      
      Real x_xi  = a_mapPartials(iv,xDir + 0*SpaceDim); //0
      Real y_xi  = a_mapPartials(iv,yDir + 0*SpaceDim); //1
      
      Real x_eta = a_mapPartials(iv,xDir + 1*SpaceDim); //2
      Real y_eta = a_mapPartials(iv,yDir + 1*SpaceDim); //3
        
      Real S  = x_xi*x_xi + y_xi*y_xi + x_eta*x_eta + y_eta*y_eta;
      Real Sx_xi  = 2.0*x_xi ;
      Real Sy_xi  = 2.0*y_xi ;
      Real Sx_eta = 2.0*x_eta;
      Real Sy_eta = 2.0*y_eta;
      
      //jacobian
      Real J  = a_Jacobian(iv,0                        );    //0
      
      // derivatives with respect to derivatives of the mapping
      Real Jx_xi  = a_Jacobian(iv,1 + 0*SpaceDim + xDir  );   //1
      Real Jy_xi  = a_Jacobian(iv,1 + 0*SpaceDim + yDir  );   //2
      Real Jx_eta = a_Jacobian(iv,1 + 1*SpaceDim + xDir  );   //3
      Real Jy_eta = a_Jacobian(iv,1 + 1*SpaceDim + yDir  );   //4
      
      // derivatives of the above with respect to xi and eta
      //Real JxXiXi   = a_Jacobian(iv,1 + 2*SpaceDim + xDir  );//5
      //Real JyXiXi   = a_Jacobian(iv,1 + 2*SpaceDim + yDir  );//6
      //Real JxEtaEta = a_Jacobian(iv,1 + 3*SpaceDim + xDir  );//7
      //Real JyEtaEta = a_Jacobian(iv,1 + 3*SpaceDim + yDir  );//8
      
      // derivatives with respect to xi and eta
      //Real Jxi  = a_Jacobian(iv,1 + 4*SpaceDim + xiDir );    //9
      //Real Jeta = a_Jacobian(iv,1 + 4*SpaceDim + etaDir);    //10
      
      // Diagonal contribution of linearized operator
      Real x_2xi_D  = -2.0/dx2;
      Real x_2eta_D = -2.0/dy2;
      
      Real y_2xi_D  = -2.0/dx2;
      Real y_2eta_D = -2.0/dy2;      
      
      //Lagrangian = pow(S,m)*pow(J,n) = L1*L2
      
      Real L1 = pow(S,m);
      Real L2 = pow(J,n);

      //Real L1_xi  = alpha*(m*pow(S,m-1)*S_xi );
      //Real L1_eta = beta *(m*pow(S,m-1)*S_eta); 

      Real Eq1L1_xi_D  = m*pow(S,m-1)*2.0*x_xi *x_2xi_D ;
      Real Eq2L1_xi_D  = m*pow(S,m-1)*2.0*y_xi *y_2xi_D ;
      Real Eq1L1_eta_D = m*pow(S,m-1)*2.0*x_eta*x_2eta_D;
      Real Eq2L1_eta_D = m*pow(S,m-1)*2.0*y_eta*y_2eta_D; 
     
      Real L1x_xi  = m*pow(S,m-1)*Sx_xi ;
      Real L1y_xi  = m*pow(S,m-1)*Sy_xi ;
      Real L1x_eta = m*pow(S,m-1)*Sx_eta; 
      Real L1y_eta = m*pow(S,m-1)*Sy_eta;

      //Real L1x_xi_xi   = alpha*(m*(m-1)*pow(S,m-2)*Sx_xi  + m*pow(S,m-1)*2.0*x_2xi );
      //Real L1y_xi_xi   = alpha*(m*(m-1)*pow(S,m-2)*Sy_xi  + m*pow(S,m-1)*2.0*y_2xi );
      //Real L1x_eta_eta = beta *(m*(m-1)*pow(S,m-2)*Sx_eta + m*pow(S,m-1)*2.0*x_2eta); 
      //Real L1y_eta_eta = beta *(m*(m-1)*pow(S,m-2)*Sy_eta + m*pow(S,m-1)*2.0*y_2eta);

      Real L1x_xi_xi_D   = m*pow(S,m-1)*2.0*x_2xi_D ;
      Real L1y_xi_xi_D   = m*pow(S,m-1)*2.0*y_2xi_D ;
      Real L1x_eta_eta_D = m*pow(S,m-1)*2.0*x_2eta_D; 
      Real L1y_eta_eta_D = m*pow(S,m-1)*2.0*y_2eta_D;

      // Real L2_xi  = alpha*(n*pow(J,n-1)*J_xi );
      // Real L2_eta = beta *(n*pow(J,n-1)*J_eta); 

      Real Eq1L2_xi_D  =  n*pow(J,n-1)*( y_eta*x_2xi_D); 
      Real Eq2L2_xi_D  =  n*pow(J,n-1)*(-x_eta*y_2xi_D); 
      Real Eq1L2_eta_D =  n*pow(J,n-1)*(-y_xi*x_2eta_D); 
      Real Eq2L2_eta_D =  n*pow(J,n-1)*( x_xi*y_2eta_D); 
      
      Real L2x_xi  = n*pow(J,n-1)*Jx_xi ;
      Real L2y_xi  = n*pow(J,n-1)*Jy_xi ;
      Real L2x_eta = n*pow(J,n-1)*Jx_eta; 
      Real L2y_eta = n*pow(J,n-1)*Jy_eta; 

      // L2x_xi_xi   = n*(n-1)*pow(J,n-2)*J_xi *Jx_xi  + n*pow(J,n-1)*( y_eta_xi);
      // L2y_xi_xi   = n*(n-1)*pow(J,n-2)*J_xi *Jy_xi  + n*pow(J,n-1)*(-x_eta_xi);
      // L2x_eta_eta = n*(n-1)*pow(J,n-2)*J_eta*Jx_eta + n*pow(J,n-1)*(-y_xi_eta); 
      // L2y_eta_eta = n*(n-1)*pow(J,n-2)*J_eta*Jy_eta + n*pow(J,n-1)*( x_xi_eta); 
      Real Eq1J_xi_D  =  y_eta*x_2xi_D ;
      Real Eq2J_xi_D  = -x_eta*y_2xi_D ;
      Real Eq1J_eta_D = -y_xi *x_2eta_D;
      Real Eq2J_eta_D =  x_xi *y_2eta_D;

      Real L2x_xi_xi_D   = n*(n-1)*pow(J,n-2)*Eq1J_xi_D *Jx_xi ;
      Real L2y_xi_xi_D   = n*(n-1)*pow(J,n-2)*Eq2J_xi_D *Jy_xi ; 
      Real L2x_eta_eta_D = n*(n-1)*pow(J,n-2)*Eq1J_eta_D*Jx_eta;
      Real L2y_eta_eta_D = n*(n-1)*pow(J,n-2)*Eq2J_eta_D*Jy_eta;
          
      // Lx_xi  = alpha*(L1x_xi *L2 + L1*L2x_xi );     
      // Ly_xi  = alpha*(L1y_xi *L2 + L1*L2y_xi ); 
      // Lx_eta = alpha*(L1x_eta*L2 + L1*L2x_eta);     
      // Ly_eta = alpha*(L1y_eta*L2 + L1*L2y_eta); 

      // Lx_xi_xi   = alpha*(L1x_xi_xi  *L2 + L1x_xi *L2_xi  + L1_xi *L2x_xi  + L1*L2x_xi_xi  );     
      // Ly_xi_xi   = alpha*(L1y_xi_xi  *L2 + L1y_xi *L2_xi  + L1_xi *L2y_xi  + L1*L2y_xi_xi  ); 
      // Lx_eta_eta = beta *(L1x_eta_eta*L2 + L1x_eta*L2_eta + L1_eta*L2x_eta + L1*L2x_eta_eta);
      // Ly_eta_eta = beta *(L1y_eta_eta*L2 + L1y_eta*L2_eta + L1_eta*L2y_eta + L1*L2y_eta_eta);
       
      Real Lx_xi_xi_D   = alpha*(L1x_xi_xi_D  *L2 + L1x_xi *Eq1L2_xi_D  + Eq1L1_xi_D *L2x_xi  + L1*L2x_xi_xi_D);
      Real Ly_xi_xi_D   = alpha*(L1y_xi_xi_D  *L2 + L1y_xi *Eq2L2_xi_D  + Eq2L1_xi_D *L2y_xi  + L1*L2y_xi_xi_D);
      
      Real Lx_eta_eta_D = beta*(L1x_eta_eta_D*L2 + L1x_eta*Eq1L2_eta_D + Eq1L1_eta_D*L2x_eta + L1*L2x_eta_eta_D);
      Real Ly_eta_eta_D = beta*(L1y_eta_eta_D*L2 + L1y_eta*Eq2L2_eta_D + Eq2L1_eta_D*L2y_eta + L1*L2y_eta_eta_D);
      
      // assign data 
      a_smoothDiag(iv,firstEq) = -(Lx_xi_xi_D + Lx_eta_eta_D); 
      a_smoothDiag(iv,secEq  ) = -(Ly_xi_xi_D + Ly_eta_eta_D);
    }
}
void computeBdSmJnDiagonalTerm(FArrayBox                & a_smoothDiag ,
                               const RealVect           & a_dx         ,
                               const FArrayBox          & a_mapPartials,
                               const FArrayBox          & a_Jacobian   ,
                               const Vector<FArrayBox*> & a_tangent    ,
                               const Vector<int>        & a_tanDir     )
{
  int m = LARGEREALVAL;
  int n = LARGEREALVAL;

  ParmParse pp;
  pp.get("Spower",m);
  pp.get("Jpower",n);

  int xDir = 0;
  int yDir = 1;

  //int xiDir = 0;
  //int etaDir = 1;
  
  //int xiEtaComp = 0;
 
  // iterate
  for (int iEdge = 0; iEdge < 4; ++iEdge)
    {
      int hiLo = 1;
      if ((iEdge == 0) || (iEdge == 3))
        {
          hiLo = -1;
        }
      
      int tanDir = a_tanDir[iEdge];

      Box tanBox = a_tangent[iEdge]->box();
      Box tanBoxInt = tanBox; 
      tanBoxInt.grow(tanDir,-1);

      for (BoxIterator bit(tanBoxInt); bit.ok(); ++bit)
        {
          const IntVect iv = bit();
         
          Real xiSWeight  = 1.0;
          Real etaSWeight = 1.0;
      
          Real alpha = xiSWeight;
          Real beta  = etaSWeight;

          Real x_xi    = a_mapPartials(iv,xDir + 0*SpaceDim); //0
          Real y_xi    = a_mapPartials(iv,yDir + 0*SpaceDim); //1
          
          Real x_eta   = a_mapPartials(iv,xDir + 1*SpaceDim); //2
          Real y_eta   = a_mapPartials(iv,yDir + 1*SpaceDim); //3
          
          Real x_xi_D  = hiLo*1.5/a_dx[xDir];
          Real y_xi_D  = hiLo*1.5/a_dx[xDir];
          Real x_eta_D = hiLo*1.5/a_dx[yDir];
          Real y_eta_D = hiLo*1.5/a_dx[yDir];

          Real S = x_xi*x_xi + y_xi*y_xi + x_eta*x_eta + y_eta*y_eta;
          Real Sx_xi  = 2.0*x_xi ;
          Real Sy_xi  = 2.0*y_xi ;
          Real Sx_eta = 2.0*x_eta;
          Real Sy_eta = 2.0*y_eta;
           
          Real L1 = pow(S,m);

          //J
          Real J = a_Jacobian(iv,0                      ); //0
          
          //Jx_xi,Jy_xi,Jx_eta,Jy_eta
          Real Jx_xi  = a_Jacobian(iv,1 + 0*SpaceDim + xDir  ) ; //1
          Real Jy_xi  = a_Jacobian(iv,1 + 0*SpaceDim + yDir  ) ; //2
          Real Jx_eta = a_Jacobian(iv,1 + 1*SpaceDim + xDir  ) ; //3
          Real Jy_eta = a_Jacobian(iv,1 + 1*SpaceDim + yDir  ); //4
          
          //Jx_xi_xi,Jy_xi_xi,Jx_eta_eta,Jy_eta_eta
          //Real Jx_xi_xi   = a_Jacobian(iv,1 + 2*SpaceDim + xDir  ); //5
          //Real Jy_xi_xi   = a_Jacobian(iv,1 + 2*SpaceDim + yDir  ) ; //6
          //Real Jx_eta_eta = a_Jacobian(iv,1 + 3*SpaceDim + xDir  ); //7
          //Real Jy_eta_eta = a_Jacobian(iv,1 + 3*SpaceDim + yDir  ); //8
          
          //Jxi,Jeta
          //Real J_xi  = a_Jacobian(iv,1 + 4*SpaceDim + xiDir ); //9
          //Real J_eta = a_Jacobian(iv,1 + 4*SpaceDim + etaDir); //10
          Real L2 = pow(J,n);

          // iEdge = 1 and 3
          // L1x_xi = alpha*(m*pow(S,m-1)*Sx_xi) 
          // L1y_xi = alpha*(m*pow(S,m-1)*Sx_xi) 
          
          // L2x_xi = alpha*(n*pow(J,n-1)*Jx_xi) 
          // L2y_xi = alpha*(n*pow(J,n-1)*Jy_xi)

          //Real L = L1*L2
          // Lx_xi = L1x_xi*L2 + L1*L2x_xi;
          // Ly_xi = L1y_xi*L2 + L1*L2y_xi;

          // linearize boundary operator to extract diagonal term
          Real L1x_xi_D = x_xi_D*alpha*(m*(m-1)*pow(S,m-2)*Sx_xi*Sx_xi + m*pow(S,m-1)*2.0);
          Real L1y_xi_D = y_xi_D*alpha*(m*(m-1)*pow(S,m-2)*Sy_xi*Sy_xi + m*pow(S,m-1)*2.0);

          Real L2x_xi_D = x_xi_D*alpha*(n*(n-1)*pow(J,n-2)*Jx_xi*Jx_xi);
          Real L2y_xi_D = y_xi_D*alpha*(n*(n-1)*pow(J,n-2)*Jy_xi*Jy_xi);

          Real Lx_xi_D = L1x_xi_D*L2 + L1*L2x_xi_D;
          Real Ly_xi_D = L1y_xi_D*L2 + L1*L2y_xi_D;

          // iEdges = 0 and 2 
          // L1x_eta = beta*(m*pow(S,m-1)*Sx_eta); 
          // L1y_eta = beta*(m*pow(S,m-1)*Sy_eta);
 
          // L2x_eta = alpha*(n*pow(J,n-1)*Jx_eta); 
          // L2y_eta = alpha*(n*pow(J,n-1)*Jy_eta);

          // L = L1*L2
          // Lx_eta = L1x_eta*L2 + L1*L2x_eta;
          // Ly_eta = L1y_eta*L2 + L1*L2y_eta;

          // linearize boundary operator to extract diagonal term
          Real L1x_eta_D = x_eta_D*alpha*(m*(m-1)*pow(S,m-2)*Sx_eta*Sx_eta + m*pow(S,m-1)*2.0);
          Real L1y_eta_D = y_eta_D*alpha*(m*(m-1)*pow(S,m-2)*Sy_eta*Sy_eta + m*pow(S,m-1)*2.0);
          
          Real L2x_eta_D = x_eta_D*alpha*(n*(n-1)*pow(J,n-2)*Jx_eta*Jx_eta);
          Real L2y_eta_D = y_eta_D*alpha*(n*(n-1)*pow(J,n-2)*Jy_eta*Jy_eta);

          Real Lx_eta_D = L1x_eta_D*L2 + L1*L2x_eta_D;
          Real Ly_eta_D = L1y_eta_D*L2 + L1*L2y_eta_D;

          // left or right edges
          if ((iEdge == 1) || (iEdge == 3))
            {
              a_smoothDiag(iv,xDir) = alpha*Lx_xi_D;
              a_smoothDiag(iv,yDir) = alpha*Ly_xi_D;
            }
          // top or bottom edges
          else
            {
              if (beta != 0)
                {
                  a_smoothDiag(iv,xDir) = beta*Lx_eta_D;
                  a_smoothDiag(iv,yDir) = beta*Ly_eta_D;
                }
              else
                {
                  ParmParse pp;
                  Real incScale;
                  pp.get("incScale",incScale);
                  if (incScale == 0.0)
                    {
                      a_smoothDiag(iv,xDir) = alpha*Lx_eta_D;
                      a_smoothDiag(iv,yDir) = alpha*Ly_eta_D;
                    }
                  else
                    {
                      MayDay::Abort("beta = 0.0 doesn't work when moving boundary");
                    }
                }
            }
        }
    }
}
void incrementBdSpacing(Vector<RealVect>       & a_movedVec         ,
                        Vector<Real>           & a_arcLengthVec     ,
                        const Vector<Real>     & a_targetSpacingVec ,
                        const Vector<Real>     & a_incrementVec     ,
                        const Vector<RealVect> & a_hiResVec         ,
                        const Vector<Real>     & a_hiResArcLengthVec,
                        const Box              & a_validBox         ,
                        const bool             & a_writeInfo        ,
                        const int              & a_iter             )
{
  int xDir = 0;
  int yDir = 1;
  
  for (int iVal = 1; iVal < a_arcLengthVec.size(); ++iVal)
    {
      a_arcLengthVec[iVal]  = a_targetSpacingVec[iVal]; 
    }
  
  //int lastShiftedCell = LARGEREALVAL;
  
  // end points of the edge
  // RealVect xPoint      = a_hiResVec[0];
  //RealVect otherCorner = a_hiResVec[a_hiResVec.size() -1];

  IntVect smallEnd = a_validBox.smallEnd();
  IntVect bigEnd   = a_validBox.bigEnd  ();
  
  int numRow = bigEnd[yDir] - smallEnd[yDir] + 1;
  int numCol = bigEnd[xDir] - smallEnd[xDir] + 1;
  RealVect dx(1.0/(numCol - 1),1.0/(numRow - 1));

  int numFaceIv = a_hiResVec.size();
  int currNewCell = 1;
  
  // arcLengthVec[0] is undefined. arcLengthVec[j] is the spacing between nodes (j-1) and j. 
  Real inc = a_arcLengthVec[currNewCell];

  for(int iOldCell  = 1; iOldCell< numFaceIv; iOldCell++)
    {
      while(inc < a_hiResArcLengthVec[iOldCell])
        {
          // relative discrepancy between hi and lo cells
          Real frac = inc/a_hiResArcLengthVec[iOldCell];
          
          // assign new node
          a_movedVec[currNewCell] = (1.0 - frac)*a_hiResVec[iOldCell - 1] + frac*a_hiResVec[iOldCell];
          
          // next new node
          currNewCell += 1;
                
          if (iOldCell < numFaceIv - 1)
            {
              if (currNewCell < a_arcLengthVec.size())
                {
                  inc += a_arcLengthVec[currNewCell];
                }
              else
                {
                  MayDay::Abort("check total arc length");
                }
            }
          else
            {
              break;
            }
        }

      // measure where even spacing multiples land, relative to current arc length increment 
      inc -= a_hiResArcLengthVec[iOldCell];
    }
 
  a_movedVec[0] = a_hiResVec[0];
  a_movedVec[a_movedVec.size() - 1] = a_hiResVec[numFaceIv - 1];

  // debugging output
  Vector<Real> newSpacing(a_movedVec.size());

  for (int iCell = 1; iCell < a_movedVec.size(); ++iCell)
    {
      newSpacing[iCell]  = (a_movedVec[iCell][0] - a_movedVec[iCell - 1][0])*(a_movedVec[iCell][0] - a_movedVec[iCell - 1][0]);
      newSpacing[iCell] += (a_movedVec[iCell][1] - a_movedVec[iCell - 1][1])*(a_movedVec[iCell][1] - a_movedVec[iCell - 1][1]);
      newSpacing[iCell] = sqrt(newSpacing[iCell]);
    }
    
  newSpacing[0] = 0;
#if 1
  if (a_writeInfo)
    {
      for (int iSpacing = 0; iSpacing < newSpacing.size(); ++iSpacing)
        {
          pout()<<"newSpacing["<<iSpacing<<"] = "<<newSpacing[iSpacing]<<endl;
        }
      pout()<<endl;
    }
#endif
}

void incrementBd(FArrayBox                & a_curr          ,
                 FArrayBox                & a_xi_eta        ,
                 Vector<FArrayBox*>       & a_arcLength     ,
                 const Vector<FArrayBox*> & a_targetSpacing , 
                 const Vector<FArrayBox*> & a_increment     ,
                 const FArrayBox          & a_hiRes         ,
                 const int                & a_refRatio      ,
                 const Vector<FArrayBox*> & a_hiResArcLength,
                 const bool               & a_writeInfo     ,
                 const std::string        & a_blockName     ,
                 const int                & a_iter          )
{
  Box box = a_curr.box     ();
  IntVect lo = box.smallEnd();
  IntVect hi = box.bigEnd  ();
  
  IntVect hiResLo = lo*a_refRatio;
  IntVect hiResHi = hi*a_refRatio;
  IntVect numHiResIv = hiResHi - hiResLo + IntVect::Unit;
    
  int xDir = 0;
  int yDir = 1;
  
  int numCol = 1 + hi[xDir] - lo[xDir];
  int numRow = 1 + hi[yDir] - lo[yDir];

  for (int iEdge = 0; iEdge < 4; ++iEdge)
    {
      const FArrayBox& targetSpacing = *(a_targetSpacing [iEdge]);

      FArrayBox& arcLength           = *(a_arcLength     [iEdge]);
      FArrayBox& increment           = *(a_increment     [iEdge]);
      FArrayBox& hiResArcLength      = *(a_hiResArcLength[iEdge]);
    
      int tanDir    = 1     ;
      int norDir    = 0     ;
      int numFaceIv = numRow;
      
      if (iEdge == 0 || iEdge == 2)
        {
          tanDir = 0;
          norDir = 1;
          numFaceIv = numCol;
        }
      
      Side::LoHiSide side = Side::Lo;
      if (iEdge == 1 || iEdge == 2)
        {
          side = Side::Hi;
        }
      
      // 0 or 1 valued
      int hilo = (1 + sign(side))/2;
      
      // x,y
      Vector<RealVect> hiResVec     (numHiResIv[tanDir]);
      Vector<Real> hiResArcLengthVec(numHiResIv[tanDir]);
      
      Vector<RealVect> movedVec     (numFaceIv);
      Vector<Real> incrementVec     (numFaceIv);
      Vector<Real> arcLengthVec     (numFaceIv);
      Vector<Real> targetSpacingVec (numFaceIv);
          
      int iStart;
      int iEnd;
      int iInc;
      
      int jStart;
      int jEnd;
      int jInc;
      
      // hi res version
      if (norDir == 0)
        {
          iStart  = hiResLo[norDir];
          iStart += hilo*(hiResHi[norDir] - hiResLo[norDir]);
          iEnd    = iStart;
          iInc    = 1;
          
          jStart = hiResLo[tanDir];
          jEnd   = hiResHi[tanDir];
          jInc   = 1;
        }
      else
        {
          jStart  = hiResLo[norDir];
          jStart += hilo*(hiResHi[norDir] - hiResLo[norDir]);
          jEnd    = jStart;
          jInc    = 1;
          
          iStart = hiResLo[tanDir];
          iEnd   = hiResHi[tanDir];
          iInc   = 1;
        }
      
      int iVec = 0; 
      for (int i = iStart; i*iInc <= iEnd*iInc; i += iInc)
        {
          for (int j = jStart; j*jInc <= jEnd*jInc; j += jInc)
            {
              IntVect hiResIv(i,j);
                  
              hiResArcLengthVec[iVec] = hiResArcLength(hiResIv,0);
                  
              for (int iComp = 0; iComp < a_curr.nComp(); ++iComp)
                {
                  hiResVec[iVec][iComp] = a_hiRes(hiResIv,iComp);
                }
                  
              iVec += 1;
            }
        }
          
      if (norDir == 0)
        {
          iStart  = lo[norDir];
          iStart += hilo*(hi[norDir] - lo[norDir]);
          iEnd    = iStart;
          iInc    = 1;
          
          jStart = lo[tanDir];
          jEnd   = hi[tanDir];
          jInc   = 1;
        }
      else
        {
          jStart  = lo[norDir];
          jStart += hilo*(hi[norDir] - lo[norDir]);
          jEnd    = jStart;
          jInc    = 1;
          
          iStart = lo[tanDir];
          iEnd   = hi[tanDir];
          iInc   = 1;
        }
      
      iVec = 0; 
      for (int i = iStart; i*iInc <= iEnd*iInc; i += iInc)
        {
          for (int j = jStart; j*jInc <= jEnd*jInc; j += jInc)
            {
              IntVect iv(i,j);
              
              arcLengthVec     [iVec] = arcLength     (iv,0);
              targetSpacingVec [iVec] = targetSpacing (iv,0); 
              incrementVec     [iVec] = increment     (iv,0);
                           
              // diagnostic debugging
              Real incScale = LARGEREALVAL;
              ParmParse pp;
              pp.get ("incScale",incScale);
              incrementVec[iVec] *= incScale;
                               
              iVec += 1;
            }
        }
           
      // change spacing
      incrementBdSpacing(movedVec         ,
                         arcLengthVec     ,
                         targetSpacingVec ,
                         incrementVec     ,
                         hiResVec         ,
                         hiResArcLengthVec,
                         box              ,
                         a_writeInfo      ,
                         a_iter           );
      
      // fill unit square boundary values
      iVec = 0; 
      for (int i = iStart; i*iInc <= iEnd*iInc; i += iInc)
        {
          for (int j = jStart; j*jInc <= jEnd*jInc; j += jInc)
            {
              IntVect iv(i,j);
              arcLength(iv,0) = arcLengthVec[iVec];
              for (int iComp = 0; iComp < SpaceDim; ++iComp)
                {
                  a_curr(iv,iComp) = movedVec[iVec][iComp];
                }

              iVec += 1;
            }
        }
    }
}

void pointJacobi(FArrayBox         & a_curr       ,
                 FArrayBox         & a_mapPartials,  
                 FArrayBox         & a_smoothDiag ,
                 FArrayBox         & a_Jacobian   ,
                 Real              & a_maxNorm    ,
                 const Box         & a_validBox   ,
                 const int         & a_iter       ,
                 const bool        & a_writeInfo  ,
                 const bool        & a_writePlot  ,
                 const RealVect    & a_dx         ,
                 const int         & a_baseRes    ,
                 const std::string & a_blockName  )
                 
{
  int numComp = a_curr.nComp();
  Box box = a_curr.box();
 
  // allocate increments
  FArrayBox inc      (box,numComp);
  FArrayBox smoothRes(box,numComp);
      
  // set the increments to zero
  inc      .setVal(0.0);
  smoothRes.setVal(0.0);
     
  // compute change (up to a factor of lambda) between two iterations
  computeInc(inc           ,
             smoothRes     ,
             a_validBox    ,
             a_curr        , 
             a_mapPartials ,
             a_smoothDiag  ,
             a_Jacobian    ,
             a_writeInfo   ,
             a_dx          ,
             a_iter        );

  // increment interior
  Box intBox = box;
  intBox.grow(-1);

  // comp
  int srcComp  = 0;
  int destComp = 0;
  
  // update curr interior
  a_curr.plus(inc     ,
              intBox  ,
              srcComp ,
              destComp,
              numComp );
          
  // norm of R,Z increments
  a_maxNorm = inc.norm(0,0,2);

  Real maxInc = 0.0;
  maxInc = Max(Max(Abs(inc.min(0)),Abs(inc.max(0))),Max(Abs(inc.min(1)),Abs(inc.max(1))));
          
  Real maxSmoothRes = 0.0;
  maxSmoothRes = Max(Max(Abs(smoothRes.min(0)),Abs(smoothRes.max(0))),Max(Abs(smoothRes.min(1)),Abs(smoothRes.max(1))));
          
  pout() << "max inc            : " << maxInc       << endl;
  pout() << "max smooth residual: " << maxSmoothRes << endl;
          
  if(a_writePlot) 
    {  
      int printFreq = 1;
      ParmParse pp;
      pp.get("printFreq_ghostCellExtrapBS",printFreq);
      if ((a_iter % printFreq) == 0)
        {  
          // 
          string prefix("diagnostic/");
          std::string suffix(".2d.hdf5");
              
          // dummy default argument
          Vector<string> componentNames;

          char name[1000];
          // increment
          sprintf(name,"diagnostic/%03d.inc.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&inc,name,componentNames,a_dx[0]);
              
          // smooth residual
          sprintf(name,"diagnostic/%03d.smoothRes.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&smoothRes,name,componentNames,a_dx[0]);
              
          // Jacobian
          sprintf(name,"diagnostic/%03d.Jac.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_Jacobian,name,componentNames,a_dx[0]);
              
          // current map
          sprintf(name,"diagnostic/%03d.curr.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_curr,name,componentNames,a_dx[0]);
        }
    }
  
  // recompute unit derivative
  computeMapPartials(a_mapPartials,
                     a_curr       ,
                     a_validBox   ,
                     a_dx         ,
                     a_writeInfo  );

  computeJacobianPartials(a_Jacobian   ,
                          a_mapPartials,
                          a_writeInfo  );
            
  computeSmJnDiagonalTerm(a_smoothDiag ,
                          a_dx         ,
                          a_Jacobian   ,
                          a_mapPartials);
          
  if(a_writePlot) 
    { 
      int printFreq = 1;
      ParmParse pp;
      pp.get("printFreq_ghostCellExtrapBS",printFreq);
      if ((a_iter % printFreq) == 0)
        {  
          // 
          string prefix("diagnostic/");
          std::string suffix(".2d.hdf5");
                  
          // dummy default argument
          Vector<string> componentNames;
                  
          char name[1000];

          // current map partials
          sprintf(name,"diagnostic/%03d.mapPartials.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_mapPartials,name,componentNames,a_dx[0]);

          // current smoothDiag
          sprintf(name,"diagnostic/%03d.smoothDiag.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_smoothDiag,name,componentNames,a_dx[0]);
                  
          // current Jacobian
          sprintf(name,"diagnostic/%03d.Jacobian.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_Jacobian,name,componentNames,a_dx[0]);
        }
    }
}

void refineExtendedMap(FArrayBox      & a_hiResFab   ,
                       const FArrayBox& a_extendedFab,
                       const int      & a_refRatio   )

{
  Box hiResBox = a_hiResFab.box();
   
  IntVect hiVal = a_hiResFab.box().bigEnd  ();
  IntVect loVal = a_hiResFab.box().smallEnd();

  for (BoxIterator bit(a_hiResFab.box());bit.ok();++bit) 
    {
      IntVect iv = bit();
      
      // rely on automatic flooring of rational number
      IntVect coarseIvLoLo = iv/a_refRatio;
      
      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          if (iv[idir] == hiVal[idir])
            {
              coarseIvLoLo -= BASISV(idir);
            }
        }
    
      IntVect coarseIvHiLo = coarseIvLoLo + BASISV(0);
      IntVect coarseIvLoHi = coarseIvLoLo + BASISV(1);

      IntVect coarseIvHiHi = coarseIvLoLo + BASISV(0)+ BASISV(1);

      // alpha is a value between 0 and 1
      RealVect alpha(iv[0]/(Real)a_refRatio,iv[1]/(Real)a_refRatio);
      alpha -= coarseIvLoLo;
      
      for (int iComp = 0; iComp < SpaceDim; ++iComp)
        {
          a_hiResFab(iv,iComp) = (1.0 - alpha[1])*((1.0 - alpha[0])*a_extendedFab(coarseIvLoLo,iComp) + alpha[0]*a_extendedFab(coarseIvHiLo,iComp))
            +        alpha[1] *((1.0 - alpha[0])*a_extendedFab(coarseIvLoHi,iComp) + alpha[0]*a_extendedFab(coarseIvHiHi,iComp));
        }
    }
} 

void computeArcLength(FArrayBox       & a_arcLength     ,
                      FArrayBox       & a_hiResArcLength,
                      Real            & a_evenSpacingVal,
                      Real            & a_totalArcLength,
                      const FArrayBox & a_hiRes         ,
                      const int       & a_tanDir        ,
                      const int       & a_refRatio      ,
                      const Box       & a_hiResEdgeBox  )
{
  int norDir = (a_tanDir + 1) % SpaceDim;
  int zeroComp = 0;
  Box intBox = a_hiResEdgeBox;

  //shrink the box
  intBox.growLo(a_tanDir,-1);
  
  // incremental counters
  Real totalHiResLength = 0.0;
  Real totalLoResLength = 0.0;
  Real loResSegLength   = 0.0;
 
  for (BoxIterator bit(intBox);bit.ok();++bit) 
    {
      IntVect iv = bit();
      
      // seg length at high resolution
      Real hiResSegLength = 0.0;
      for (int iComp = 0; iComp < SpaceDim; ++iComp)
        {
          hiResSegLength += (a_hiRes(iv,iComp) - a_hiRes(iv - BASISV(a_tanDir),iComp))*(a_hiRes(iv,iComp) - a_hiRes(iv - BASISV(a_tanDir),iComp));
        }
      
      hiResSegLength   = sqrt(hiResSegLength);
      
      // increment counters
      loResSegLength   += hiResSegLength;
      totalHiResLength += hiResSegLength;

      // fill data holder
      a_hiResArcLength(iv,zeroComp) = hiResSegLength;
      
      // sub-sample for lo res
      if (((iv[a_tanDir]) % a_refRatio) == 0)
        {
          int loResTan = iv[a_tanDir];
          int loResNor = iv[norDir  ];
          loResTan   /= a_refRatio;
          loResNor   /= a_refRatio;
          IntVect loResIv = loResTan*BASISV(a_tanDir) + loResNor*BASISV(norDir);
          a_arcLength(loResIv,zeroComp) = loResSegLength;
          totalLoResLength             += loResSegLength;

          // reset for the refRatio number of terms
          loResSegLength = 0.0;
        }
    }

  // use to calculate the spacing for evenly spaced boundary nodes
  int numSpaces = a_arcLength.box().size()[a_tanDir] - 1;
 
  // spacing
  a_evenSpacingVal = totalHiResLength/numSpaces;
  
  a_totalArcLength = totalHiResLength;
  // hi and lo res values should be equal
  pout()<<"totalHiResLength = "<<totalHiResLength<<endl;
  pout()<<"totalLoResLength = "<<totalLoResLength<<endl;
}

void setTargetSpacing(FArrayBox        & a_targetSpacing  ,  
                      const Real       & a_outerTotalLength,
                      const Real       & a_innerTotalLength,
                      const FArrayBox  & a_innerArcLength  ,
                      const int        & a_tanDir          ,
                      const IntVect    & a_offset          )
{
  // indexing macro
  int compZero = 0;

  Box edgeBoxInt = a_targetSpacing.box();
  
  //shrink the box
  edgeBoxInt.growLo(a_tanDir,-1);
  
  Real stretchingFactor = a_outerTotalLength/a_innerTotalLength;
  Real totalNewLength = 0;
  for (BoxIterator bit(edgeBoxInt);bit.ok();++bit) 
    {
      IntVect iv = bit();
      IntVect innerIv = iv + a_offset;
      a_targetSpacing(iv,compZero) = stretchingFactor*a_innerArcLength(innerIv,compZero);
      totalNewLength +=  a_targetSpacing(iv,compZero);
    }

  for (BoxIterator bit(edgeBoxInt);bit.ok();++bit) 
    {
      IntVect iv = bit();
      a_targetSpacing(iv,compZero) *= a_outerTotalLength/totalNewLength;
    }
  Real totalAdjustedLength = 0;
  for (BoxIterator bit(edgeBoxInt);bit.ok();++bit) 
    {
      IntVect iv = bit();
      totalAdjustedLength += a_targetSpacing(iv,compZero);
    }
  pout()<< totalAdjustedLength - a_outerTotalLength<<endl;
}

void moveExtendedBdNodes(Vector<FArrayBox*>         & a_extendedMap        , 
                         const int                  & a_refRatio           ,
                         const int                  & a_numBlocks          ,
                         const Vector<int>          & a_tanDir             ,
                         const bool                 & a_writeInfo          ,
                         const Vector<Box>          & a_box                ,
                         const Vector<Vector<Box> > & a_outerEdgeBoxes     ,
                         const Vector<Vector<Box> > & a_innerEdgeBoxes     ,
                         const Vector<Box>          & a_hiResBox           ,
                         const Vector<Vector<Box> > & a_outerHiResEdgeBoxes,
                         const Vector<Vector<Box> > & a_innerHiResEdgeBoxes,
                         const int                  & a_beginBlock         ,
                         const int                  & a_endBlock           ,
                         const int                  & a_iter               ,
                         const int                  & a_numGhost           )
{
  int numComp = 2;
  int oneComp = 1;

  // Define hiRes and xi_eta
  Vector<FArrayBox*> hiRes(a_numBlocks);

  Vector<FArrayBox*> xi_eta(a_numBlocks);
  
  // Define arcLength, hiResArcLength, and targetSpacing
  Vector<Vector< FArrayBox*> >  outerArcLength     (a_numBlocks);
  Vector<Vector< FArrayBox*> >  outerHiResArcLength(a_numBlocks);
  Vector<Vector< FArrayBox*> >  outerTargetSpacing (a_numBlocks);

  // Define arcLength, hiResArcLength, and targetSpacing
  Vector<Vector< FArrayBox*> >  innerArcLength     (a_numBlocks);
  Vector<Vector< FArrayBox*> >  innerHiResArcLength(a_numBlocks);

  //allocate memory
  for (int iBlock = a_beginBlock; iBlock <= a_endBlock; ++iBlock)
    {
      // FArray boxes
      hiRes           [iBlock] = new FArrayBox(a_hiResBox[iBlock],numComp);
      xi_eta          [iBlock] = new FArrayBox(a_box     [iBlock],numComp);
     
      // data holders for edges
      outerArcLength        [iBlock].resize(4);
      outerTargetSpacing    [iBlock].resize(4);
      outerHiResArcLength   [iBlock].resize(4);

      innerArcLength        [iBlock].resize(4);
      innerHiResArcLength   [iBlock].resize(4);
      
      for(int iEdge = 0;iEdge< 4; ++iEdge)
        {
          // ordering of edges: bottom,right,top,left
          outerArcLength     [iBlock][iEdge] = new FArrayBox(a_outerEdgeBoxes     [iBlock] [iEdge],oneComp);
          outerTargetSpacing [iBlock][iEdge] = new FArrayBox(a_outerEdgeBoxes     [iBlock] [iEdge],oneComp);
          outerHiResArcLength[iBlock][iEdge] = new FArrayBox(a_outerHiResEdgeBoxes[iBlock] [iEdge],oneComp);

          innerArcLength     [iBlock][iEdge] = new FArrayBox(a_innerEdgeBoxes     [iBlock] [iEdge],oneComp);
          innerHiResArcLength[iBlock][iEdge] = new FArrayBox(a_innerHiResEdgeBoxes[iBlock] [iEdge],oneComp);
        }
    }

  for (int iBlock = 0; iBlock < a_numBlocks; ++iBlock)
    {
      std::string blockName;
      getBlockName(blockName,
                   iBlock  );
      
      // extended data used for making hi res interpolated data
      FArrayBox & extendedFab = *(a_extendedMap[iBlock]);
      FArrayBox & hiResFab    = *(hiRes        [iBlock]);
     
      Vector<Real>evenSpacingVal(4);
      
      //refine extendedMap to make a hi res version
      refineExtendedMap(hiResFab   ,
                        extendedFab,
                        a_refRatio ); 
            
      // discretize tangent to the boundary 
      Vector <Real> outerTotalLength(4);

      for (int iEdge = 0; iEdge < 4;++iEdge)    
        {
          computeArcLength(*(outerArcLength     [iBlock][iEdge]),
                           *(outerHiResArcLength[iBlock][iEdge]),
                           evenSpacingVal               [iEdge] ,
                           outerTotalLength             [iEdge] ,  
                           *(hiRes              [iBlock]       ),
                           a_tanDir                     [iEdge] ,
                           a_refRatio                           ,
                           a_outerHiResEdgeBoxes[iBlock][iEdge]);
        }
    
      // discretize tangent to the boundary 
      Vector <Real> innerTotalLength(4);

      for (int iEdge = 0; iEdge < 4;++iEdge)    
        {
          computeArcLength(*(innerArcLength      [iBlock][iEdge]),
                           *(innerHiResArcLength [iBlock][iEdge]),
                           evenSpacingVal                [iEdge] ,
                           innerTotalLength              [iEdge] ,  
                           *(hiRes               [iBlock]       ),
                           a_tanDir                      [iEdge] ,
                           a_refRatio                            ,
                           a_innerHiResEdgeBoxes [iBlock][iEdge]);
        }

      // data holder for calculating the even spacing map
      Vector<FArrayBox*> arcLengthDummy(4);
      
      // dummy data holder
      Vector<FArrayBox*> incrementDummy(4);

      for (int iEdge = 0; iEdge < 4; ++iEdge)
        {
          if ( (iBlock < 8) || (iEdge == 0) || (iEdge == 2) ||(true))
            {
              arcLengthDummy           [iEdge] = new FArrayBox((outerArcLength[iBlock])[iEdge]->box(),1);
              (*(arcLengthDummy        [iEdge])).setVal(LARGEREALVAL);
              
              incrementDummy           [iEdge] = new FArrayBox((outerArcLength[iBlock])[iEdge]->box(),1);
              (*(incrementDummy        [iEdge])).setVal(LARGEREALVAL);
                       
              //replicate inner spacing on ghost cell boundary, compensating for a different overall boundary length
              int hiLo = -1;
              if (iEdge == 0 || iEdge == 3)
                {
                  hiLo = 1;
                }

              int norDir = 1;
              if (iEdge == 1 || iEdge == 3)
                {
                  norDir = 0;
                }
              
              IntVect offset = hiLo*a_numGhost*BASISV(norDir); 
              setTargetSpacing(*(outerTargetSpacing[iBlock][iEdge]),  
                               outerTotalLength            [iEdge] ,
                               innerTotalLength            [iEdge] ,
                               *(innerArcLength    [iBlock][iEdge]),
                               a_tanDir                    [iEdge] ,
                               offset                              );
            }
          else
            {
              arcLengthDummy           [iEdge] = new FArrayBox((outerArcLength[iBlock])[iEdge]->box(),1);
              (*(arcLengthDummy        [iEdge])).setVal(LARGEREALVAL);

              incrementDummy           [iEdge] = new FArrayBox((outerArcLength[iBlock])[iEdge]->box(),1);
              (*(incrementDummy        [iEdge])).setVal(LARGEREALVAL);

              (*(outerTargetSpacing[iBlock])[iEdge] ).copy(*((outerArcLength[iBlock])[iEdge]));
            }

        }

      // move boundary points
      incrementBd(*(a_extendedMap    [iBlock]),
                  *(xi_eta           [iBlock]),
                  arcLengthDummy              ,
                  outerTargetSpacing [iBlock] ,
                  incrementDummy              ,
                  *(hiRes            [iBlock]),
                  a_refRatio                  ,
                  outerHiResArcLength[iBlock] ,
                  a_writeInfo                 ,
                  blockName                   ,
                  a_iter                      );

    }
}

void createEdgeBoxes(Vector<Box> & a_hiResEdgeBoxes,
                     const Box   & a_domainBox     )
{
  // four edges in each block
  CH_assert(a_hiResEdgeBoxes.size() == 4);
  
  // lower left corner
  IntVect loLo = a_domainBox.smallEnd();
     
  // upper right corner
  IntVect hiHi = a_domainBox.bigEnd  ();
     
  // upper left corner
  IntVect loHi(loLo[0],hiHi[1]);
     
  // lower right corner
  IntVect hiLo(hiHi[0],loLo[1]);
   
  // iteration boxes around the edges xEdge (yEdge) means normal direction is x (y)
  Box xEdgeLo(loLo,loHi);
  Box xEdgeHi(hiLo,hiHi);
   
  Box yEdgeLo(loLo,hiLo);
  Box yEdgeHi(loHi,hiHi);
  
  // the order used is: bottom,right,top,left
  a_hiResEdgeBoxes[0] = yEdgeLo;
  a_hiResEdgeBoxes[1] = xEdgeHi;
  a_hiResEdgeBoxes[2] = yEdgeHi;
  a_hiResEdgeBoxes[3] = xEdgeLo;
}

void createInnerHiResEdgeBoxes(Vector<Box> & a_innerHiResEdgeBoxes,
                               Vector<Box> & a_outerHiResEdgeBoxes,
                               const int   & a_numGhost           ,
                               const int   & a_refRatio           )
{
  int xDir =0;
  int yDir =1;

  // the order used is: bottom,right,top,left
  Box yEdgeLo = a_outerHiResEdgeBoxes[0];
  Box xEdgeHi = a_outerHiResEdgeBoxes[1];
  Box yEdgeHi = a_outerHiResEdgeBoxes[2];
  Box xEdgeLo = a_outerHiResEdgeBoxes[3];

  // iteration boxes around the edges xEdge (yEdge) means normal direction is x (y)
  yEdgeLo += a_numGhost*a_refRatio*BASISV(yDir);
  xEdgeHi -= a_numGhost*a_refRatio*BASISV(xDir);
  yEdgeHi -= a_numGhost*a_refRatio*BASISV(yDir);
  xEdgeLo += a_numGhost*a_refRatio*BASISV(xDir);

  // the order used is: bottom,right,top,left
  a_innerHiResEdgeBoxes[0] = yEdgeLo;
  a_innerHiResEdgeBoxes[1] = xEdgeHi;
  a_innerHiResEdgeBoxes[2] = yEdgeHi;
  a_innerHiResEdgeBoxes[3] = xEdgeLo; 
}
void createInnerEdgeBoxes(Vector<Box> & a_innerEdgeBoxes,
                          Vector<Box> & a_outerEdgeBoxes,
                          const int   & a_numGhost           )
{
  int xDir =0;
  int yDir =1;

  // the order used is: bottom,right,top,left
  Box yEdgeLo = a_outerEdgeBoxes[0];
  Box xEdgeHi = a_outerEdgeBoxes[1];
  Box yEdgeHi = a_outerEdgeBoxes[2];
  Box xEdgeLo = a_outerEdgeBoxes[3];

  // iteration boxes around the edges xEdge (yEdge) means normal direction is x (y)
  yEdgeLo += a_numGhost*BASISV(yDir);
  xEdgeHi -= a_numGhost*BASISV(xDir);
  yEdgeHi -= a_numGhost*BASISV(yDir);
  xEdgeLo += a_numGhost*BASISV(xDir);

  // the order used is: bottom,right,top,left
  a_innerEdgeBoxes[0] = yEdgeLo;
  a_innerEdgeBoxes[1] = xEdgeHi;
  a_innerEdgeBoxes[2] = yEdgeHi;
  a_innerEdgeBoxes[3] = xEdgeLo; 
}

void makeGhostCells(Vector<FArrayBox*>       & a_extendedMap     ,
                    Vector<FArrayBox*>       & a_extendedBField  ,
                    const std::string        & a_geometryFileName,
                    const std::string        & a_fieldFileName   ,
                    const Vector<FArrayBox*> & a_curr            ,
                    const int                & a_numGhost        ,
                    const int                & a_numBlocks       ,
                    const int                & a_numComp         ,
                    const Vector<RealVect>   & a_dx              ,
                    const int                & a_order           ,
                    const Real               & a_startFrac       ,
                    const Real               & a_endFrac         ,
                    const bool               & a_blending        ,
                    const int                & a_num1Deriv       )
{
    
  for (int iGhost = 1; iGhost <= a_numGhost; ++iGhost)
    {
      // initialize ghost cells to 0 and interior to map
      for(int iBlock = 0; iBlock < a_numBlocks; ++iBlock)
        {
          Box extendedBox = a_curr[iBlock]->box();
          extendedBox.grow(a_numGhost);

          // allocate memory outside the function
          a_extendedMap    [iBlock] = new FArrayBox(extendedBox,a_numComp              );
          a_extendedBField [iBlock] = new FArrayBox(extendedBox,a_numComp + a_num1Deriv);
           
          // set to arbitrary value
          a_extendedBField[iBlock]->setVal(0.13579);
            
          std::string blockName;
                        
          getBlockName(blockName,
                       iBlock  );
            
          // set the extended map to zero (including ghost region)
          a_extendedMap[iBlock]->setVal(0.0);

          // copy the valid cells onto the extendedMap
          a_extendedMap[iBlock]->copy(*(a_curr[iBlock]));
        }
        
      for (int iBlock = 0; iBlock < a_numBlocks; ++iBlock)
        {
          std::string blockName;
          getBlockName(blockName,
                       iBlock  );
            
          // bottom = 0, right = 1, top =2, left = 3
          bool useExtrapolation[4];
            
          Box extendedBox = a_extendedMap[iBlock]->box();
            
          // temporary variables for blending 
          FArrayBox copyTemp(extendedBox,a_numComp);
          copyTemp.copy(*(a_extendedMap[iBlock]));
            
          FArrayBox extrapTemp(extendedBox,a_numComp);
          extrapTemp.copy(*(a_extendedMap[iBlock]));
            
          Box stencilBox  = a_curr[iBlock]->box();
          stencilBox.grow(1 - iGhost);

          if (a_numGhost > 0)        
            {
              // copy where the ghost cells are the valid cells of neighboring blocks 
              copyGhostCell(copyTemp            ,
                            useExtrapolation    ,
                            a_extendedMap         ,
                            blockName           ,
                            iBlock              ,
                            a_geometryFileName,
                            a_fieldFileName       ,
                            stencilBox          ,
                            a_numComp             , 
                            a_numGhost            );
                
                
              // here weight is dummy parameter; elsewhere extrapolation comes from two directions and solution is incremented with weight = 0.5 
              Real dummyWeight = 1.0;
              extendMap(extrapTemp           ,
                        useExtrapolation     ,
                        a_curr[iBlock]->box(),
                        a_numGhost           ,
                        a_dx[iBlock]         ,
                        a_order              ,
                        blockName            ,
                        dummyWeight          );
                
              blend(*(a_extendedMap[iBlock]),
                    a_curr[iBlock]->box()   ,
                    copyTemp                ,
                    extrapTemp              ,
                    blockName               ,
                    iBlock                  ,
                    a_numGhost              ,
                    a_startFrac             ,
                    a_endFrac               );
                
            }
        }
        
      for (int iBlock = 0; iBlock < a_numBlocks; ++iBlock)
        {
          if (a_numGhost > 0)
            {
              Interval comps(0,1);
                
              Box sourceBox;
              std::string sourceBlockName;
                
              std::string blockName;
              getBlockName(blockName,
                           iBlock   );
                
              // selective use of the extend function
              bool extrapHiLoXDir[2];
              bool extrapHiLoYDir[2];
                
              // for indexing in useHiLo
              int loSide = 0;
              int hiSide = 1;
                
              // iterate through the corners
              std::string corner;
                
              //ghost cells in lower left square
              corner = "lowerLeft";
                
              //xDir
              extrapHiLoXDir[loSide] = true ;
              extrapHiLoXDir[hiSide] = false;
                
              //yDir
              extrapHiLoYDir[loSide] = true ;
              extrapHiLoYDir[hiSide] = false;
                 
              // fill in the corners
              cornerGhostCells(a_extendedMap ,
                               corner        ,
                               a_curr        ,
                               a_numGhost    ,
                               blockName     ,
                               a_dx[iBlock]  ,
                               a_order       ,
                               extrapHiLoXDir,
                               extrapHiLoYDir,
                               iBlock        ,
                               a_blending    );
                
              // going counter-clockwise, starting lower left
              corner = "lowerRight";
                
              //xDir
              extrapHiLoXDir[loSide] = false;
              extrapHiLoXDir[hiSide] = true ;
                
              //yDir
              extrapHiLoYDir[loSide] = true;
              extrapHiLoYDir[hiSide] = false;
                
              // fill in the corners
              cornerGhostCells(a_extendedMap ,
                               corner        ,
                               a_curr        ,
                               a_numGhost    ,
                               blockName     ,
                               a_dx[iBlock]  ,
                               a_order       ,
                               extrapHiLoXDir,
                               extrapHiLoYDir,
                               iBlock        ,
                               a_blending    );
                
              // upper right
              corner = "upperRight";
                
              //xDir
              extrapHiLoXDir[loSide] = false;
              extrapHiLoXDir[hiSide] = true ;
            
              //yDir
              extrapHiLoYDir[loSide] = false;
              extrapHiLoYDir[hiSide] = true ;
                
              // fill in the corners
              cornerGhostCells(a_extendedMap ,
                               corner        ,
                               a_curr        ,
                               a_numGhost    ,
                               blockName     ,
                               a_dx [iBlock] ,
                               a_order       ,
                               extrapHiLoXDir,
                               extrapHiLoYDir,
                               iBlock        ,
                               a_blending    );
                
              //upper left
              corner = "upperLeft";
                
              //xDir
              extrapHiLoXDir[loSide] = true ;
              extrapHiLoXDir[hiSide] = false;
                
              //yDir
              extrapHiLoYDir[loSide] = false;
              extrapHiLoYDir[hiSide] = true ;
            
              // fill in the corners
              cornerGhostCells(a_extendedMap ,
                               corner        ,
                               a_curr        ,
                               a_numGhost    ,
                               blockName     ,
                               a_dx[iBlock]  ,
                               a_order       ,
                               extrapHiLoXDir,
                               extrapHiLoYDir,
                               iBlock        ,
                               a_blending    );
            }
            
        }
    }
}

#if 0
void makeWeightFunction(FArrayBox & a_blendingWeights,
                        const int & a_iBlock         ,
                        const Box & a_validBox       )
{
  for (BoxIterator bit(extendedBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (!currValidBox.contains(iv))
        {
          for (int iComp = 0; iComp<SpaceDim; ++iComp)
            {
              // at the first of four ghost cells take 
              if (whichGhostMethod == "aligned")
                {
                  blendedFab(iv,iComp) =  alignedFab(iv,iComp);
                }
              else if(whichGhostMethod == "even")
                {
                  blendedFab(iv,iComp) =  evenFab(iv,iComp);
                }
              else if (whichGhostMethod == "weighted")
                {
                  if (a_iBlock == 8)
                    {
                      int xiDir = 0;
                      int etaDir = 1;
                      IntVect lo = validBox.smallEnd();
                      IntVect hi = validBox.bigEnd();
                      IntVect loHi(lo[xiDir],hi[etaDir]);
                      IntVect hiLo(hi[xiDir],lo[etaDir]);
                      Real dx = 1.0/(hi[xiDir ] - lo[xiDir ]);
                      Real dy = 1.0/(hi[etaDir] - lo[etaDir]);
                      
                      IntVect base1(-a_numGhost,-a_numGhost);
                      IntVect base2(-a_numGhost, loHi[etaDir]+ a_numGhost);
                      IntVect base3(-a_numGhost, 36);
                      
                      Real dist1 = ((iv[0] - base1[0])*(iv[0] - base1[0])) +  ((iv[1] - base1[1])*(iv[1] - base1[1]));
                      dist1 = sqrt(dist1);
                      Real sd1 = 10;
                      int expo1 = 2;
                      Real gaussianWeight1 =  exp(-pow(dist1/sd1,expo1));
                      
                      Real dist2 = ((iv[0] - base2[0])*(iv[0] - base2[0])) +  ((iv[1] - base2[1])*(iv[1] - base2[1]));
                      dist1 = sqrt(dist2);
                      Real sd2 = 35;
                      int expo2 = 2;
                      Real gaussianWeight2 =  exp(-pow(dist2/sd2,expo2));
                      
                      Real dist3 = ((iv[0] - base3[0])*(iv[0] - base3[0])) +  ((iv[1] - base3[1])*(iv[1] - base3[1]));
                      dist1 = sqrt(dist3);
                      Real sd3 = 12;
                      int expo3 = 2;
                      Real gaussianWeight3 =  exp(-pow(dist3/sd3,expo3));
                      
                      Real gaussianWeight = gaussianWeight1 + gaussianWeight2 + gaussianWeight3;
                      blendingWeights(iv,0) = gaussianWeight;
                      
                      blendedFab(iv,iComp) = (1.0 - gaussianWeight)*alignedFab(iv,iComp) + gaussianWeight*evenFab(iv,iComp);
                    }
                }
            }
        }
    }
}
#endif
            
void blendGhostCells(Vector<FArrayBox*> &a_blendedExtendedMap,
                     Vector<FArrayBox*> &a_alignedExtendedMap,
                     Vector<FArrayBox*> &a_evenExtendedMap   ,
                     Vector<FArrayBox*> &a_alignedCurr       ,
                     int                &a_numGhost          ,
                     int                &a_numBlocks         ,
                     Vector<int>        &a_whichGhostMethod  )
{
  for (int iBlock = 0; iBlock < a_numBlocks; ++iBlock)
    {
      std::string whichGhostMethod;
      if (a_whichGhostMethod[iBlock] == 0)
        {
          whichGhostMethod = "aligned"; 
        }
      else if(a_whichGhostMethod[iBlock] == 1)
        {
          whichGhostMethod = "even"; 
        }
      else if(a_whichGhostMethod[iBlock] == 2)
        {
          whichGhostMethod = "weighted"; 
        }
      else
        {
          MayDay::Abort("bad value for which ghost method");
        }
      
      FArrayBox& blendedFab = *(a_blendedExtendedMap[iBlock]);
      FArrayBox& alignedFab = *(a_alignedExtendedMap[iBlock]);
      FArrayBox& evenFab    = *(a_evenExtendedMap   [iBlock]);
          
      Box validBox   = a_alignedCurr[iBlock]->box();
      Box ghostBox = validBox;
      ghostBox.grow(a_numGhost);

      FArrayBox gWeight(ghostBox,1);
        
      for(int iGhost = 1; iGhost<= a_numGhost; ++iGhost)
        {
          // blending weight
          Real weight = ((Real)iGhost)/(Real)a_numGhost;

          // identifies which layer of ghost cells we're blending
          Box extendedBox  = validBox;
          Box currValidBox = validBox;
          extendedBox. grow(iGhost    );
          currValidBox.grow(iGhost - 1);
          
          for (BoxIterator bit(extendedBox); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              if (!currValidBox.contains(iv))
                {
                  for (int iComp = 0; iComp<SpaceDim; ++iComp)
                    {
                      // at the first of four ghost cells take 
                      if (whichGhostMethod == "aligned")
                        {
                          blendedFab(iv,iComp) =  alignedFab(iv,iComp);
                        }
                      else if(whichGhostMethod == "even")
                        {
                          blendedFab(iv,iComp) =  evenFab(iv,iComp);
                        }
                      else if (whichGhostMethod == "weighted")
                        {
                          if (iBlock == 8)
                            {
                              int xiDir = 0;
                              int etaDir = 1;
                              IntVect lo = validBox.smallEnd();
                              IntVect hi = validBox.bigEnd();
                              IntVect loHi(lo[xiDir],hi[etaDir]);
                              IntVect hiLo(hi[xiDir],lo[etaDir]);
                              Real dx = 1.0/(hi[xiDir ] - lo[xiDir ]);
                              Real dy = 1.0/(hi[etaDir] - lo[etaDir]);
 
                              IntVect base1(-a_numGhost,-a_numGhost);
                              IntVect base2(-a_numGhost, loHi[etaDir]+ a_numGhost);
                              IntVect base3(-a_numGhost, 36);

                              Real dist1 = ((iv[0] - base1[0])*(iv[0] - base1[0])) +  ((iv[1] - base1[1])*(iv[1] - base1[1]));
                              dist1 = sqrt(dist1);
                              Real sd1 = 10;
                              int expo1 = 2;
                              Real gaussianWeight1 =  exp(-pow(dist1/sd1,expo1));
                            
                              Real dist2 = ((iv[0] - base2[0])*(iv[0] - base2[0])) +  ((iv[1] - base2[1])*(iv[1] - base2[1]));
                              dist1 = sqrt(dist2);
                              Real sd2 = 35;
                              int expo2 = 2;
                              Real gaussianWeight2 =  exp(-pow(dist2/sd2,expo2));
                            
                              Real dist3 = ((iv[0] - base3[0])*(iv[0] - base3[0])) +  ((iv[1] - base3[1])*(iv[1] - base3[1]));
                              dist1 = sqrt(dist3);
                              Real sd3 = 12;
                              int expo3 = 2;
                              Real gaussianWeight3 =  exp(-pow(dist3/sd3,expo3));
                            
                              Real gaussianWeight = gaussianWeight1 + gaussianWeight2 + gaussianWeight3;
                              gWeight(iv,0) = gaussianWeight;
                              
                              blendedFab(iv,iComp) = (1.0 - gaussianWeight)*alignedFab(iv,iComp) + gaussianWeight*evenFab(iv,iComp);

                            }
                          else
                            {
                              blendedFab(iv,iComp) =  (1.0 - weight)*alignedFab(iv,iComp) +  weight*evenFab(iv,iComp);
                            }
                        }
                      else
                        {
                          MayDay::Abort("bad values for whichGhostMethod");
                        }
                    }
                } 
            }
        }
      if (iBlock == 8)
        {
           char fileName[80];
        
               
        // mapping with ghost cells 
        sprintf(fileName,"diagnostic/gWeight.%s.2d.hdf5","mcore");
      
        Real dummy = 1;
        
        // dummy default argument
        Vector<string> componentNames;
          // mapping with ghost cells
        writeFABname(&gWeight      ,
                     fileName      ,
                     componentNames,
                     dummy         );
        }
    } 
}

void chopOffExtraGhost(Vector<FArrayBox*>       & a_extendedMap      ,
                       const Vector<FArrayBox*>& a_blendedExtendedMap,
                       const int               & a_numBlocks         )
{
  for (int iBlock = 0; iBlock< a_numBlocks;++iBlock)
    {
      a_extendedMap[iBlock]->copy(*(a_blendedExtendedMap[iBlock]));
    }
}
    
int main(int    a_argc,
         char** a_argv)
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
#endif
  
  // Begin forever present scoping trick
  {
    const char* inFile = NULL;
    
    if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
    else
      {
        pout() << "Usage: ghostCell...ex <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }
   
    // indexing macros
    int xDir = 0;
    int yDir =1;

    // Parse input file
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);
    
    // prolixity 
    bool writeInfo = false;
    pp.query("writeInfo",writeInfo) ;
    
    // prolixity 
    bool writePlot = false;
    pp.query("writePlot",writePlot) ;
        
    // nCell is a parameter for cell-centered data holders, which simulate node-centered data holders
    Vector<int> numCell(2);
    Real log2MinRes;
    Real log2MaxRes;
    pp.get("log2MinRes_ghostCellExtrapBS",log2MinRes);
    pp.get("log2MaxRes_ghostCellExtrapBS",log2MaxRes);
    
    int numBlocks = -LARGEINTVAL;
    pp.get("numBlocks", numBlocks);

    int refRatioEtaInput = LARGEINTVAL;
    pp.get("refRatioEta",refRatioEtaInput);
    refRatioEtaInput -= 1;

    Vector<int> refRatioXi (numBlocks);
    Vector<int> refRatioEta(numBlocks);
    for (int iBlock = 0; iBlock < numBlocks;++iBlock)
      {
        refRatioXi [iBlock] = 0;
        refRatioEta[iBlock] = 0;
      }

    // set mcore and mcsol
    refRatioEta[8] = refRatioEtaInput;
    refRatioEta[9] = refRatioEtaInput;
   
    int numGhost = LARGEINTVAL;
    int order    = LARGEINTVAL;
    pp.get("numGhost",numGhost);
    pp.get("order"   ,order   );

    std::string outputFileName;
    pp.get("outputFileName",outputFileName);
             
    bool blending;
    pp.get("blending",blending);
    
    Real startFrac = LARGEREALVAL;
    Real endFrac = LARGEREALVAL;
    if (blending)
      {
        pp.get("startFrac",startFrac);
        pp.get("endFrac",endFrac);
      }
    else
      {
        startFrac = 1.0;
        endFrac   = 1.0;
      }

    int numExtraGhost = LARGEINTVAL;
    pp.get("numExtraGhost",numExtraGhost);
    int numTrueGhost  = numGhost;

    numGhost += numExtraGhost;
    
    if (numGhost > 0 && numBlocks < 10)
      {
        MayDay::Abort("must do all ten blocks to get ghost cells");
      }
    
    int refRatio = LARGEREALVAL;
    pp.get("refRatio",refRatio);

    int beginBlock = -LARGEINTVAL;
    pp.get("beginBlock", beginBlock);
    
    int endBlock = -LARGEINTVAL;
    pp.get("endBlock", endBlock);

    Vector<int> numRow                     (numBlocks);
    Vector<int> numCol                     (numBlocks);
    Vector<IntVect> loEnd                  (numBlocks);
    Vector<IntVect> hiEnd                  (numBlocks);
    Vector<Box> box                        (numBlocks);
    Vector<Box> extendedBox                (numBlocks);
    Vector<Box> trueExtendedBox            (numBlocks);
    Vector<IntVect> hiCorner               (numBlocks);
    Vector<RealVect> dx                    (numBlocks);
    Vector<IntVect> exLoEnd                (numBlocks);
    Vector<IntVect> exHiEnd                (numBlocks); 
    Vector<Box> hiResBox                   (numBlocks);
    Vector<IntVect> hiResHiCorner          (numBlocks); 
    Vector<Vector<Box> >outerEdgeBoxes     (numBlocks);    
    Vector<Vector<Box> >innerEdgeBoxes     (numBlocks);    
    Vector<Vector<Box> >outerHiResEdgeBoxes(numBlocks);
    Vector<Vector<Box> >innerHiResEdgeBoxes(numBlocks);
    
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        numRow[iBlock] = pow(2,log2MinRes + refRatioEta[iBlock]);
        numCol[iBlock] = pow(2,log2MinRes + refRatioXi [iBlock]);

        loEnd[iBlock] = IntVect::Zero;
        hiEnd[iBlock][0] = numCol[iBlock];
        hiEnd[iBlock][1] = numRow[iBlock];

        box  [iBlock].define(loEnd[iBlock],hiEnd[iBlock]);
        hiCorner[iBlock]= box[iBlock].bigEnd();
        
        extendedBox[iBlock] = box[iBlock];
        extendedBox[iBlock].grow(numGhost);

        dx[iBlock][0] = 1.0/Real(numCol[iBlock]);
        dx[iBlock][1] = 1.0/Real(numRow[iBlock]);

        trueExtendedBox[iBlock] = box[iBlock];
        trueExtendedBox[iBlock].grow(numTrueGhost);
               
        exLoEnd [iBlock]    = extendedBox[iBlock].smallEnd();
        exHiEnd [iBlock]    = extendedBox[iBlock].bigEnd  ();

        hiResBox     [iBlock].define(refRatio*exLoEnd[iBlock],refRatio*exHiEnd[iBlock]);
        hiResHiCorner[iBlock]= hiResBox[iBlock].bigEnd();
        
        outerEdgeBoxes     [iBlock].resize(4);
        innerEdgeBoxes     [iBlock].resize(4);
        innerHiResEdgeBoxes[iBlock].resize(4);
        outerHiResEdgeBoxes[iBlock].resize(4);

        createEdgeBoxes(outerEdgeBoxes[iBlock],
                        extendedBox   [iBlock]);

        createEdgeBoxes(outerHiResEdgeBoxes[iBlock],
                        hiResBox           [iBlock]);
        
        createInnerHiResEdgeBoxes(innerHiResEdgeBoxes[iBlock],
                                  outerHiResEdgeBoxes[iBlock],
                                  numGhost           ,
                                  refRatio           );
        
        createInnerEdgeBoxes(innerEdgeBoxes[iBlock],
                             outerEdgeBoxes[iBlock],
                             numGhost              );
      }

    // no loop over resolutions yet in ghost cell code.
    int baseRes = pow(2,log2MinRes);

    char suffix[1000];
    sprintf(suffix,".%03d.2d.txt",baseRes + 1);
    
    // global fieldData object: block independent
    std::string alignedGeometryBaseName;
    pp.get("alignedGeometryBaseNameGhost", alignedGeometryBaseName);
    std::string alignedGeometryFileName = alignedGeometryBaseName + suffix;

    // global fieldData object: block independent
    std::string evenGeometryBaseName;
    pp.get("evenGeometryBaseNameGhost", evenGeometryBaseName);
    std::string evenGeometryFileName = evenGeometryBaseName + suffix;
    
    std::string fieldFileName;
    pp.get("fieldFileNameGhost", fieldFileName); 
    // mapped quantities have two components
    int numComp    = 2 ;
 
    // clean start for mapping file output function    
    clearExtendedMappingFile(outputFileName,
                             baseRes      );
      
    //data holders for valid cells on all blocks 
    Vector<FArrayBox*> evenCurr      (numBlocks);
    Vector<FArrayBox*> alignedCurr   (numBlocks);
    
    //data holders with ghost cells on all blocks 
    Vector<FArrayBox*> evenExtendedMap     (numBlocks);
    Vector<FArrayBox*> alignedExtendedMap  (numBlocks);
    Vector<FArrayBox*> blendedExtendedMap  (numBlocks);
    Vector<FArrayBox*> extendedMap         (numBlocks);
    
       // indexing macro for component number
    int num1Deriv = 4;
    
   // boxes used for tangent vector FArrayBoxes
    
    Vector<FArrayBox*> extendedBField    (numBlocks);
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        // data holders for valid cells
        evenCurr         [iBlock] = new FArrayBox(box[iBlock],numComp);
        alignedCurr      [iBlock] = new FArrayBox(box[iBlock],numComp);
        
        // data holders with ghost cells
        evenExtendedMap   [iBlock] = new FArrayBox(extendedBox    [iBlock],numComp);
        alignedExtendedMap[iBlock] = new FArrayBox(extendedBox    [iBlock],numComp);
        blendedExtendedMap[iBlock] = new FArrayBox(extendedBox    [iBlock],numComp);
        extendedMap       [iBlock] = new FArrayBox(trueExtendedBox[iBlock],numComp);

        extendedBField    [iBlock] = new FArrayBox(extendedBox    [iBlock],numComp + num1Deriv);
      }
       
    // fill in curr for each block copy onto currIter
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        std::string blockName;
        getBlockName(blockName,
                     iBlock  );
            
        // Construct the field data object
        FieldData field_data(evenGeometryFileName,
                             fieldFileName       ,
                             blockName           );
            
        setPhysicalCoordinates(*evenCurr[iBlock]    ,
                               hiCorner [iBlock]    ,
                               blockName            ,
                               evenGeometryFileName ,
                               fieldFileName        );
        
        // write even spaced grid before extending
        // output
        char fileName[80];
        
        // dummy default argument
        Vector<string> componentNames;
        
        // mapping with ghost cells 
        sprintf(fileName,"diagnostic/%03d.even.final.%s.2d.hdf5",baseRes,blockName.c_str());
      
        // mapping with ghost cells
        writeFABname(evenCurr[iBlock],
                     fileName        ,
                     componentNames  ,
                     dx[iBlock][0]   );
      }
   
    makeGhostCells(evenExtendedMap     ,
                   extendedBField      ,
                   evenGeometryFileName,
                   fieldFileName       ,
                   evenCurr            ,
                   numGhost            ,
                   numBlocks           ,
                   numComp             ,
                   dx                  ,
                   order               ,
                   startFrac           ,
                   endFrac             ,
                   blending            ,
                   num1Deriv           );

    //brackbill saltzman on ghost region
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        std::string blockName;
        getBlockName(blockName,
                     iBlock  );
        
        // Construct the field data object
        FieldData field_data(alignedGeometryFileName,
                             fieldFileName          ,
                             blockName              );
        
        setPhysicalCoordinates(*alignedCurr[iBlock]   ,
                               hiCorner    [iBlock]   ,
                               blockName              ,
                               alignedGeometryFileName,
                               fieldFileName          );
        
        // output
        char fileName[80];
        
        // dummy default argument
        Vector<string> componentNames;
        
        // mapping with ghost cells 
        sprintf(fileName,"diagnostic/%03d.evenExtended.final%s.2d.hdf5",baseRes,blockName.c_str());
      
        // mapping with ghost cells
        writeFABname(evenExtendedMap[iBlock],
                     fileName               ,
                     componentNames         ,
                     dx[iBlock][0]          );
        
        // copy the aligned valid cells into the evenExtendedMap, in preparation for Brackbill-Salzman on ghost cells
        evenExtendedMap[iBlock]->copy(*(alignedCurr[iBlock]));

        // mapping with ghost cells 
        sprintf(fileName,"diagnostic/%03d.evenAndAlignedBeforeBS.%s.2d.hdf5",baseRes,blockName.c_str());
      
        // mapping with ghost cells
        writeFABname(evenExtendedMap[iBlock],
                     fileName               ,
                     componentNames         ,
                     dx[iBlock][0]          );
      }
   
    // ordering of edges is bottom,right,top,left
    Vector<int> tanDir(4); 
    tanDir[0] = xDir;
    tanDir[1] = yDir;
    tanDir[2] = xDir;
    tanDir[3] = yDir;

    int dummy = LARGEINTVAL;

    // respace the nodes on the extended boundary
    moveExtendedBdNodes(evenExtendedMap    , 
                        refRatio           ,
                        numBlocks          ,
                        tanDir             ,
                        writeInfo          ,
                        box                ,
                        outerEdgeBoxes     ,
                        innerEdgeBoxes     ,
                        hiResBox           ,
                        outerHiResEdgeBoxes,
                        innerHiResEdgeBoxes,
                        beginBlock         ,
                        endBlock           ,
                        dummy              ,
                        numGhost           );

    // define data holders for pre-calculated derivatives
    Vector<FArrayBox*> mapPartials(numBlocks);
    Vector<FArrayBox*> smoothDiag (numBlocks);
    Vector<FArrayBox*> Jacobian   (numBlocks);
    Real maxNorm;
    
    // mapped quantities have two components
    int numJacComp = 11;
    int num2Deriv  = 6 ;
    
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        Box extendedBox = (*(evenExtendedMap[iBlock])).box();

        smoothDiag  [iBlock] = new FArrayBox(extendedBox,numComp              );
        Jacobian    [iBlock] = new FArrayBox(extendedBox,numJacComp           );
        mapPartials [iBlock] = new FArrayBox(extendedBox,num1Deriv + num2Deriv);
      }
    
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        
        computeMapPartials(*(mapPartials    [iBlock])       ,
                           *(evenExtendedMap[iBlock])       ,
                           (*(alignedCurr   [iBlock])).box(),
                           dx               [iBlock]        ,
                           writeInfo                        );
        
        computeJacobianPartials(*(Jacobian   [iBlock]),
                                *(mapPartials[iBlock]),
                                writeInfo             );
        
        computeSmJnDiagonalTerm(*(smoothDiag [iBlock]),
                                dx           [iBlock] ,
                                *(Jacobian   [iBlock]),
                                *(mapPartials[iBlock]));
      }
    
    int maxIter_ghostCellExtrapBS;
    pp.get("maxIter_ghostCellExtrapBS",maxIter_ghostCellExtrapBS);
     
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        for (int iter = 0; iter <maxIter_ghostCellExtrapBS; ++ iter)
          {
            std::string blockName;
            getBlockName(blockName,
                         iBlock   );
                
            Box validBox = alignedCurr[iBlock]->box();
                
            // call pointJacobi on extendedMap with fixed central square
            pointJacobi(*(evenExtendedMap[iBlock]),
                        *(mapPartials[iBlock])    ,  
                        *(smoothDiag [iBlock])    ,
                        *(Jacobian   [iBlock])    ,
                        maxNorm                   ,
                        validBox                  ,
                        iter                      ,
                        writeInfo                 ,
                        writePlot                 ,
                        dx           [iBlock]     ,
                        baseRes                   ,
                        blockName                 );
          }
      }
    
    // extrapolate from field aligned valid points
    // fill in curr for each block copy onto currIter
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        std::string blockName;
        getBlockName(blockName,
                     iBlock  );
            
        // Construct the field data object
        FieldData field_data(alignedGeometryFileName,
                             fieldFileName          ,
                             blockName              );
            
        setPhysicalCoordinates(*(alignedCurr[iBlock]) ,
                               hiCorner     [iBlock]  ,
                               blockName              ,
                               alignedGeometryFileName,
                               fieldFileName          );
        
        blendedExtendedMap[iBlock]->copy(*(alignedCurr[iBlock]));
      }
   
    makeGhostCells(alignedExtendedMap     ,
                   extendedBField         ,
                   alignedGeometryFileName,
                   fieldFileName          ,
                   alignedCurr            ,
                   numGhost               ,
                   numBlocks              ,
                   numComp                ,
                   dx                     ,
                   order                  ,
                   startFrac              ,
                   endFrac                ,
                   blending               ,
                   num1Deriv              );
    
    // combine the ghost cells extrapolated from aligned curr and the Brackbill ghost cells
    Vector<int> whichGhostMethod(numBlocks);
    pp.getarr ("whichGhostMethod",whichGhostMethod,0,numBlocks - 1);
  
    blendGhostCells(blendedExtendedMap,
                    alignedExtendedMap,
                    evenExtendedMap   ,
                    alignedCurr       ,
                    numGhost          ,
                    numBlocks         , 
                    whichGhostMethod  );
    
    chopOffExtraGhost(extendedMap       ,
                      blendedExtendedMap,
                      numBlocks         );

    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        std::string blockName;
        getBlockName(blockName,
                     iBlock   );
        
        // use FieldData clas for magnetic field and first partials
        setBField(*(extendedBField      [iBlock]),
                  *(extendedMap         [iBlock]),
                  dx                    [iBlock] ,
                  alignedGeometryFileName        ,
                  fieldFileName                  ,
                  blockName                      );
        
        // write rawMapping file
        outputMapping(*(extendedMap   [iBlock]),
                      *(extendedBField[iBlock]),
                      numGhost*IntVect::Unit   ,
                      blockName                ,
                      outputFileName           );
        
        // output
        char fileName[80];
        
        // dummy default argument
        Vector<string> componentNames;
        
        // mapping with ghost cells 
        sprintf(fileName,"diagnostic/%03d.extended.final.%s.2d.hdf5",baseRes,blockName.c_str());
        
        // mapping with ghost cells
        writeFABname(extendedMap[iBlock],
                     fileName           ,
                     componentNames     ,
                     dx[iBlock][0]      );
        // mapping with ghost cells 
        sprintf(fileName,"diagnostic/%03d.alignedExtended.final.%s.2d.hdf5",baseRes,blockName.c_str());
        
        // mapping with ghost cells
        writeFABname(alignedExtendedMap[iBlock],
                     fileName                  ,
                     componentNames            ,
                     dx[iBlock][0]             );

        // mapping with ghost cells 
        sprintf(fileName,"diagnostic/%03d.evenAndAligned.final.%s.2d.hdf5",baseRes,blockName.c_str());
        
        // mapping with ghost cells
        writeFABname(evenExtendedMap[iBlock],
                     fileName      ,
                     componentNames,
                     dx[iBlock][0] );

        // mapping with ghost cells 
        sprintf(fileName,"diagnostic/%03d.aligned.final.%s.2d.hdf5",baseRes,blockName.c_str());
      
        // mapping with ghost cells
        writeFABname(alignedCurr[iBlock],
                     fileName           ,
                     componentNames     ,
                     dx[iBlock][0]      );
      
        // mapping with ghost cells 
        sprintf(fileName,"diagnostic/%03d.extendedBField.final.%s.2d.hdf5",baseRes,blockName.c_str());
      
        // mapping with ghost cells
        writeFABname(extendedBField[iBlock]   ,
                     fileName                 ,
                     componentNames           ,
                     dx            [iBlock][0]);
      }
  
    //delete pointers (last in first out)
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        delete evenExtendedMap    [iBlock];
        delete alignedExtendedMap [iBlock];
        delete extendedBField     [iBlock];
        delete evenCurr           [iBlock];
        delete alignedCurr        [iBlock];
        delete blendedExtendedMap [iBlock];
        delete extendedMap        [iBlock];
      }
  
    

  }// End scoping trick
 
#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
 
  return 0;
  
}

