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

void cornerAngles(Vector<Real>   & a_angles       ,
                  const FArrayBox& a_physicalCoord,
                  const int      & a_blockNumber  )
{
  int xDir = 0;
  int yDir = 1;
  Box validBox = a_physicalCoord.box();

  int firstRow = validBox.smallEnd()[xDir];
  int firstCol = validBox.smallEnd()[yDir]; 

  int lastRow = validBox.bigEnd()[xDir];
  int lastCol = validBox.bigEnd()[yDir]; 
  
  IntVect P00(firstRow,firstCol);
  IntVect P000 = P00 + BASISV(xDir);
  IntVect P001 = P00 + BASISV(yDir);

  IntVect P10(lastCol,firstRow);
  IntVect P100 = P10 - BASISV(xDir);
  IntVect P101 = P10 + BASISV(yDir);

  IntVect P11(lastCol,lastRow);
  IntVect P110 = P11 - BASISV(xDir);
  IntVect P111 = P11 - BASISV(yDir);

  IntVect P01(firstCol,lastRow);
  IntVect P010 = P01 + BASISV(xDir);
  IntVect P011 = P01 - BASISV(yDir);
  
  
  RealVect P00Pt (a_physicalCoord(P00 ,xDir),a_physicalCoord(P00 ,yDir));
  RealVect P000Pt(a_physicalCoord(P000,xDir),a_physicalCoord(P000,yDir));
  RealVect P001Pt(a_physicalCoord(P001,xDir),a_physicalCoord(P001,yDir));

  RealVect P10Pt (a_physicalCoord(P10 ,xDir),a_physicalCoord(P10 ,yDir));
  RealVect P100Pt(a_physicalCoord(P100,xDir),a_physicalCoord(P100,yDir));
  RealVect P101Pt(a_physicalCoord(P101,xDir),a_physicalCoord(P101,yDir));

  RealVect P11Pt (a_physicalCoord(P11 ,xDir),a_physicalCoord(P11 ,yDir));
  RealVect P110Pt(a_physicalCoord(P110,xDir),a_physicalCoord(P110,yDir));
  RealVect P111Pt(a_physicalCoord(P111,xDir),a_physicalCoord(P111,yDir));
  
  RealVect P01Pt (a_physicalCoord(P01 ,xDir),a_physicalCoord(P01 ,yDir));
  RealVect P010Pt(a_physicalCoord(P010,xDir),a_physicalCoord(P010,yDir));
  RealVect P011Pt(a_physicalCoord(P011,xDir),a_physicalCoord(P011,yDir));
    
  RealVect V000 =  (P00Pt - P000Pt);
  RealVect V001 =  (P00Pt - P001Pt);

  RealVect V100 =  (P10Pt - P100Pt);
  RealVect V101 =  (P10Pt - P101Pt);
  
  RealVect V110 =  (P11Pt - P110Pt);
  RealVect V111 =  (P11Pt - P111Pt);

  RealVect V010 =  (P01Pt - P010Pt);
  RealVect V011 =  (P01Pt - P011Pt);
 
  Real M000 = sqrt(V000[0]*V000[0] + V000[1]*V000[1]);
  Real M001 = sqrt(V001[0]*V001[0] + V001[1]*V001[1]);

  Real M100 = sqrt(V100[0]*V100[0] + V100[1]*V100[1]);
  Real M101 = sqrt(V101[0]*V101[0] + V101[1]*V101[1]);

  Real M110 = sqrt(V110[0]*V110[0] + V110[1]*V110[1]);
  Real M111 = sqrt(V111[0]*V111[0] + V111[1]*V111[1]);

  Real M010 = sqrt(V010[0]*V010[0] + V010[1]*V010[1]);
  Real M011 = sqrt(V011[0]*V011[0] + V011[1]*V011[1]);

  V000 /= M000;
  V001 /= M001;

  V100 /= M100;
  V101 /= M101;

  V110 /= M110;
  V111 /= M111;

  V010 /= M010;
  V011 /= M011;
  
  Real cosA00 = V000[0]*V001[0] + V000[1]*V001[1];
  Real cosA11 = V110[0]*V111[0] + V110[1]*V111[1];
  Real cosA01 = V010[0]*V011[0] + V010[1]*V011[1];
  Real cosA10 = V100[0]*V101[0] + V100[1]*V101[1];

  // looking for the angle in the third quadrant
  Real A00 = M_PI - acos(cosA00);
  Real A10 =        acos(cosA10 );
  Real A11 = M_PI - acos(cosA11);
  Real A01 =        acos(cosA01 );
 
  // lower left
  a_angles[0] = A00;

  // lower right
  a_angles[1] = A10;

  //  upper right
  a_angles[2] = A11;

  // upper left
  a_angles[3] = A01;
}

RealVect interpolateAngle(const FArrayBox    & a_curr      ,
                          const IntVect      & a_iv        , 
                          const int          & a_firstRow  ,
                          const int          & a_lastRow   ,
                          const bool         & a_goingLower)
{ 
  RealVect retVal;
  
  int xDir = 0;
  int yDir = 1;
  
  int hiLo = 1;
  if (a_goingLower)
    {
      hiLo = -1;
    }
  Real yDirFrac = a_iv[yDir]/((Real)a_lastRow);

  IntVect bdLo0(a_iv[xDir] - hiLo,a_firstRow);
  IntVect bdLo1(a_iv[xDir]       ,a_firstRow);
  IntVect bdHi0(a_iv[xDir] - hiLo,a_lastRow );
  IntVect bdHi1(a_iv[xDir]       ,a_lastRow );

  RealVect ptLo0(a_curr(bdLo0,xDir),a_curr(bdLo0,yDir));
  RealVect ptLo1(a_curr(bdLo1,xDir),a_curr(bdLo1,yDir));
  RealVect ptHi0(a_curr(bdHi0,xDir),a_curr(bdHi0,yDir));
  RealVect ptHi1(a_curr(bdHi1,xDir),a_curr(bdHi1,yDir));
  
  RealVect loDir = ptLo1;
  loDir         -= ptLo0;
 
  RealVect hiDir = ptHi1;
  hiDir         -= ptHi0;
   
  retVal = (1.0 - yDirFrac)*loDir + (yDirFrac)*hiDir;

  Real mag = sqrt(retVal[0]*retVal[0] +retVal[1]*retVal[1]); 

  retVal /= mag;
  retVal *= -1.0;
  return (retVal);
}

RealVect interpolateAngle(const Vector<Real> & a_angles     ,
                          const RealVect     & a_normalDir  ,
                          const IntVect      & a_iv         ,
                          const int          & a_lastRow    ,
                          const int          & a_lastCol    ,
                          const int          & a_blockNumber)
{ 
  RealVect retVal;
  int xDir = 0;
  int yDir = 1;

  Real xDirFrac = a_iv[xDir]/((Real)a_lastCol);
  Real yDirFrac = a_iv[yDir]/((Real)a_lastRow);

  Real lowerLeft  = a_angles[0];
  Real lowerRight = a_angles[1];
  Real upperRight = a_angles[2]; 
  Real upperLeft  = a_angles[3]; 

  Real loAvg = lowerLeft*(1.0 - xDirFrac) + lowerRight*xDirFrac;
  Real hiAvg = upperLeft*(1.0 - xDirFrac) + upperRight*xDirFrac;
  
  Real theta = (1.0 - yDirFrac)*loAvg + (yDirFrac)*hiAvg;

  theta -= 0.5*M_PI;
    
  retVal[xDir] = cos(theta)*a_normalDir[xDir] - sin(theta)*a_normalDir[yDir];
  retVal[yDir] = sin(theta)*a_normalDir[xDir] + cos(theta)*a_normalDir[yDir];
  
  return (retVal);
}
void  fillTargetFlux(FArrayBox       & a_magFluxFab,
                     const FArrayBox & a_curr      ,
                     const FieldData & a_field_data)
{
  int xDir = 0;
  int yDir = 1;

  Box box = a_magFluxFab.box();
  for (BoxIterator bit(box);bit.ok();++bit) 
    {
      IntVect iv = bit();
      RealVect pt(a_curr(iv,xDir),a_curr(iv,yDir));
      Real flux = a_field_data.getMagneticFluxFromDCT(pt);
      a_magFluxFab(iv,0) = flux;
    }
}

Real  getTargetFlux(const IntVect   & a_iv        ,
                    const FArrayBox & a_curr      ,
                    const FieldData & a_field_data)
{
  int xDir = 0;
  int yDir = 1;

  IntVect bigEnd = a_curr.box().bigEnd();
  IntVect smallEnd = a_curr.box().smallEnd();

  IntVect loEndCol(a_iv[xDir],smallEnd[yDir]);
  IntVect hiEndCol(a_iv[xDir],bigEnd  [yDir]);

  Real frac = ((Real)a_iv[yDir])/((Real)bigEnd[yDir]);
  RealVect loPt(a_curr(loEndCol,xDir),a_curr(loEndCol,yDir));
  RealVect hiPt(a_curr(hiEndCol,xDir),a_curr(hiEndCol,yDir));
 
  Real retval =  (1.0 - frac)*a_field_data.getMagneticFluxFromDCT(loPt) + frac*a_field_data.getMagneticFluxFromDCT(hiPt);

  return retval;
}

RealVect findRoot(const RealVect  & a_loPt       , 
                  const RealVect  & a_hiPt       ,
                  const Real      & a_targetFlux ,
                  const FieldData & a_field_data ,
                  const Real      & a_tolerance  ,
                  const int       & a_blockNumber,
                  const IntVect   & a_iv         )
{
  RealVect newPt = LARGEREALVAL*RealVect::Unit;
  
  // lo point
  RealVect loPt = a_loPt;

  // hi point
  RealVect hiPt = a_hiPt;
   
  bool findingRoot = true;
  while(findingRoot)
    {
      // mid point
      RealVect midPt = hiPt;
      midPt         += loPt;
      midPt         *= 0.5;
      
      Real fluxLo   = a_field_data.getMagneticFluxFromDCT(loPt );
      Real fluxHi   = a_field_data.getMagneticFluxFromDCT(hiPt );
      Real fluxMid  = a_field_data.getMagneticFluxFromDCT(midPt);
      
      Real diffLo   = a_targetFlux - fluxLo ;
      Real diffHi   = a_targetFlux - fluxHi ;
      Real diffMid  = a_targetFlux - fluxMid;
      
      if (diffHi * diffLo < 0.0)
        {
          if (diffHi * diffMid >= 0.0)
            {
              hiPt = midPt;
            }
          else if (diffMid * diffLo >= 0)
            {
              loPt = midPt;
            }
          else
            {
              pout()<< "Problem in block "<<a_blockNumber<<", at iv = "<<a_iv<<endl; 
              FILE* fd;
              fd = fopen ("done","w");
              if (fd == NULL) 
                {
                  perror ("Error opening file");
                }
              fclose(fd);
              MayDay::Abort("points must bracket target flux value");
            }
        }
      else
        {
          pout()<< "Problem in block "<<a_blockNumber<<", at iv = "<<a_iv<<endl;
          FILE* fd;
          fd = fopen ("done","w");
          if (fd == NULL) 
            {
              perror ("Error opening file");
            }
          fclose(fd);
          MayDay::Abort("points must bracket target flux value");
        }
      
      pout()<< "diffMid = "<<diffMid<<endl;
      if (Abs(diffMid) < a_tolerance)
        {
          findingRoot = false;  
          newPt = midPt;
        }
    }
  return newPt;
}

void orthoMarch(FArrayBox         & a_curr            ,
                const FArrayBox   & a_BField          ,
                const FArrayBox   & a_perpField       ,
                const FArrayBox   & a_magFluxFab      ,
                const RealVect    & a_dx              ,
                const std::string & a_geometryFileName,
                const std::string & a_fieldFileName   ,
                const std::string & a_blockName       ,
                const bool        & a_writeInfo       ,
                const int         & a_startHiLo       ,
                const int         & a_marchDirection  ,
                const Real        & a_tolerance       ,
                const int         & a_blockNumber     ,
                const int         & a_numSubSteps     )
{
  // Construct the field data object
  FieldData field_data(a_geometryFileName,
                       a_fieldFileName   ,
                       a_blockName       );

  int xDir = 0;
  int yDir = 1;
 
  // iterate
  Box validBox  = a_curr.box       ();
  IntVect loEnd = validBox.smallEnd();
  IntVect hiEnd = validBox.bigEnd  ();
  
  int firstRow = loEnd[ a_marchDirection              ];
  int firstCol = loEnd[(a_marchDirection + 1)%SpaceDim];
  
  int lastCol  = hiEnd[ a_marchDirection              ];
  int lastRow  = hiEnd[(a_marchDirection + 1)%SpaceDim];

  FArrayBox currLo(validBox,SpaceDim);
  FArrayBox currHi(validBox,SpaceDim);
  currLo.setVal(LARGEREALVAL);
  currHi.setVal(LARGEREALVAL);
  
  // initialize boundary values for currLo and currHi
  currLo.copy(a_curr);
  currHi.copy(a_curr);
  Real dist;
  Real prevDist;
  RealVect loPt;
  RealVect direction;
#if 0  
  Vector<Real> angles;
  angles.resize(4);
  cornerAngles(angles       ,
               currLo       ,
               a_blockNumber);
#endif
#if 1
  for (int iRow = firstRow + 1; iRow < lastRow; ++iRow)
    {
      for (int iCol = firstCol + 1; iCol <= lastCol; ++iCol)
        {
          IntVect iv{iCol,iRow};
          
          // (possibly interpolated) flux value for for this column
          Real targetFlux = getTargetFlux(iv        ,
                                          currLo    ,
                                          field_data);

          if (iCol < lastCol)
            {
              // flux on the edge where we begin marching
              IntVect firstColIv(firstCol,iRow);
              RealVect firstColPt(currLo(firstColIv,xDir),currLo(firstColIv,yDir));
              Real firstColFlux = field_data.getMagneticFluxFromDCT(firstColPt);
              
              // find the column where the flux crosses the target flux
              int baseCol;
              for (baseCol = firstCol; baseCol <= lastCol; ++baseCol)
                {
                  // get flux in this column
                  IntVect baseColIv(baseCol,iRow);
                  RealVect baseColPt(currLo(baseColIv,xDir),currLo(baseColIv,yDir));
                  Real baseColFlux = field_data.getMagneticFluxFromDCT(baseColPt);
                  
                  if (((targetFlux - baseColFlux)*(targetFlux - firstColFlux)) < 0.0)
                    {
                      break;
                    }
                }
              if (baseCol == lastCol + 1)
                {
                  pout()<< "can't bracket flux at iv = "<<iv<<" and iBlock = "<<a_blockNumber<<endl;
                  MayDay::Abort("must bracket flux");
                }

              // calculate a sufficiently far distance to look for a crossing point
              IntVect loIv (iCol - 1,iRow);
              IntVect hiIv (baseCol ,iRow);
              
              RealVect loPt  (currLo(loIv,xDir),currLo(loIv,yDir));
              RealVect loHiPt(currLo(hiIv,xDir),currLo(hiIv,yDir));
              
              dist =  (loHiPt[xDir] - loPt[xDir])*(loHiPt[xDir] - loPt[xDir]); 
              dist += (loHiPt[yDir] - loPt[yDir])*(loHiPt[yDir] - loPt[yDir]);
              dist = sqrt(dist);
              
              // stretch in case the path is curved
              dist *= 2.0;
              prevDist = dist;
            }
          else
            {
              dist = prevDist;
            }

          Real subDist = dist/(Real)(a_numSubSteps); 
          
          RealVect fieldDir;
          IntVect baseIv(iCol - 1,iRow);
          RealVect basePt(currLo(baseIv,xDir),currLo(baseIv,yDir));
          RealVect hiPt = basePt;

          bool success;
          for (int iStep = 0; iStep < a_numSubSteps; ++iStep)
            {
              success = false;
              fieldDir = field_data.getFieldUnitVectorFromDCT(basePt);

              RealVect normalDir;
              normalDir[0] = -fieldDir[1];
              normalDir[1] =  fieldDir[0];
                                 
              bool goingLower = false;
              
 
#if 0  
              pout()<< "iv = "<<iv<<", iStep = "<< iStep;
              pout()<< ", normalDir = "<<normalDir;
              pout()<< ", direction = "<<direction<<endl;
              pout()<< "basePt = "<<basePt<<endl;
#endif              
              if (a_blockNumber < 8)
                {
                  direction =  interpolateAngle(currLo,
                                            iv        , 
                                            firstRow  ,
                                            lastRow   ,
                                            goingLower);
                }
              else
                {
                  direction = normalDir;
                }
              hiPt -= subDist*direction;

              // check whether points bracket target flux
              Real loFlux = field_data.getMagneticFluxFromDCT(basePt);
              Real hiFlux = field_data.getMagneticFluxFromDCT(hiPt  );

              Real loDiff = loFlux - targetFlux;
              Real hiDiff = hiFlux - targetFlux;

              if (loDiff*hiDiff < 0.0)
                {
                  success = true;
                  break;
                }
              else
                {
                  basePt = hiPt;
                }
            }
          
          if (!success)
            {
              pout()<<"Points don't bracket flux for currLo, iv = "<<iv<<"in block "<<a_blockNumber<<endl;
              MayDay::Abort("points don't bracket flux");
            } 
                                                
          RealVect newPt = findRoot(basePt       ,
                                    hiPt         ,
                                    targetFlux   ,
                                    field_data   ,
                                    a_tolerance  ,
                                    a_blockNumber,
                                    iv           );
                    
          currLo(iv,xDir) = newPt[xDir];
          currLo(iv,yDir) = newPt[yDir];
        }
    }
  #endif  
#if 1
  for (int iRow = firstRow + 1; iRow < lastRow ; ++iRow)
    {
      for (int iCol = lastCol -1; iCol >= firstCol; --iCol)
        {
          IntVect iv{iCol,iRow};
          
          // (possibly interpolated) flux value for for this column
          Real targetFlux = getTargetFlux(iv        ,
                                          currHi    ,
                                          field_data);
          if (iCol > firstCol)
            {
              // flux on the edge where we begin marching
              IntVect lastColIv(lastCol,iRow);
              RealVect lastColPt(currHi(lastColIv,xDir),currHi(lastColIv,yDir));
              Real lastColFlux = field_data.getMagneticFluxFromDCT(lastColPt);
              
              // find the column where the flux crosses the target flux
              int baseCol;
              for (baseCol = lastCol; baseCol >= firstCol; --baseCol)
                {
                  // get flux in this column
                  IntVect baseColIv(baseCol,iRow);
                  RealVect baseColPt(currHi(baseColIv,xDir),currHi(baseColIv,yDir));
                  Real baseColFlux = field_data.getMagneticFluxFromDCT(baseColPt);
                  
                  if (((targetFlux - baseColFlux)*(targetFlux - lastColFlux)) < 0.0)
                    {
                      break;
                    }
                }
              if (baseCol == firstCol - 1)
                {
                  pout()<< "can't bracket flux at iv = "<<iv<<" and iBlock = "<<a_blockNumber<<endl;
                  FILE* fd;
                  fd = fopen ("done","w");
                  if (fd == NULL) 
                    {
                      perror ("Error opening file");
                    }
                  fclose(fd);
                  MayDay::Abort("must bracket flux");
                }

              // calculate a sufficiently far distance to look for a crossing point
              IntVect loIv (iCol + 1,iRow);
              IntVect hiIv (baseCol ,iRow);
              
              RealVect loPt  (currHi(loIv,xDir),currHi(loIv,yDir));
              RealVect loHiPt(currHi(hiIv,xDir),currHi(hiIv,yDir));
              
              dist = (loHiPt[xDir] - loPt[xDir])*(loHiPt[xDir] - loPt[xDir]);
              dist     += (loHiPt[yDir] - loPt[yDir])*(loHiPt[yDir] - loPt[yDir]);
              dist = sqrt(dist);
              
              // stretch in case the path is curved
              dist *= 2;
              prevDist = dist;
            }
          else
            {
              dist = prevDist;
            }
          
          Real subDist = dist/((Real)a_numSubSteps); 
          
          RealVect fieldDir;
          IntVect baseIv(iCol + 1,iRow);
          RealVect basePt(currHi(baseIv,xDir),currHi(baseIv,yDir));
          RealVect hiPt = basePt;
          
          bool success;
          for (int iStep = 0; iStep < a_numSubSteps; ++iStep)
            {
              success = false;
              fieldDir = field_data.getFieldUnitVectorFromDCT(basePt);

              RealVect normalDir;
              normalDir[0] = -fieldDir[1];
              normalDir[1] =  fieldDir[0];
              normalDir *= -1.0;

              bool goingLower = true;
              if (a_blockNumber < 8)
                {
                  direction =  interpolateAngle(currLo,
                                            iv        , 
                                            firstRow  ,
                                            lastRow   ,
                                            goingLower);
                }
              else
                {
                  direction = normalDir;
                }

#if 0
              pout()<< "iRow = "<<iRow<<", iCol = "<< iCol<<", iStep = "<< iStep;
              pout()<< ", normalDir = "<<normalDir<<endl;
              pout()<< "basePt = "<<basePt<<endl;
#endif
              hiPt -= subDist*direction;
              
              // check whether points bracket target flux
              Real loFlux = field_data.getMagneticFluxFromDCT(basePt);
              Real hiFlux = field_data.getMagneticFluxFromDCT(hiPt  );

              Real loDiff = loFlux - targetFlux;
              Real hiDiff = hiFlux - targetFlux;
              if (loDiff*hiDiff < 0.0)
                {
                  success = true;
                  break;
                }
              else
                {
                  basePt = hiPt;
                }
            }
          if (!success)
            {
              pout()<<"Points don't bracket flux for currHi, iv =  "<<iv<<"in block "<<a_blockNumber<<endl;
              FILE* fd;
              fd = fopen ("done","w");
              if (fd == NULL) 
                {
                  perror ("Error opening file");
                }
              fclose(fd);
              MayDay::Abort("points don't bracket flux");
            }
          RealVect newPt = findRoot(basePt       ,
                                    hiPt         ,
                                    targetFlux   ,
                                    field_data   ,
                                    a_tolerance  ,
                                    a_blockNumber,
                                    iv           );
                    
          currHi(iv,xDir) = newPt[xDir];
          currHi(iv,yDir) = newPt[yDir];
        }
    }
#endif
   for (int iRow = 1; iRow < lastRow; ++iRow)
    {
      for (int iCol = 0; iCol <=lastCol; ++iCol)
        {
          IntVect iv{iCol,iRow};
          Real frac = ((Real)iv[xDir])/((Real)lastCol);
          //frac = 1.0;
          a_curr(iv,xDir) = (1.0 - frac)*currLo(iv,xDir) + frac*currHi(iv,xDir);
          a_curr(iv,yDir) = (1.0 - frac)*currLo(iv,yDir) + frac*currHi(iv,yDir);
        }
    }
}

void sepGrid(FArrayBox         & a_curr            ,
             const FArrayBox   & a_BField          ,
             const FArrayBox   & a_perpField       ,
             const FArrayBox   & a_magFluxFab      ,
             const RealVect    & a_dx              ,
             const std::string & a_geometryFileName,
             const std::string & a_fieldFileName   ,
             const std::string & a_blockName       ,
             const bool        & a_writeInfo       ,
             const int         & a_startHiLo       ,
             const int         & a_marchDirection  ,
             const Real        & a_tolerance       ,
             const int         & a_blockNumber     ,
             const int         & a_numSubSteps     )
{
  // Construct the field data object
  FieldData field_data(a_geometryFileName,
                       a_fieldFileName   ,
                       a_blockName       );

  int xDir = 0;
  int yDir = 1;
 
  // iterate
  Box validBox  = a_curr.box       ();
  IntVect loEnd = validBox.smallEnd();
  IntVect hiEnd = validBox.bigEnd  ();
  
  int firstRow = loEnd[ a_marchDirection              ];
  int firstCol = loEnd[(a_marchDirection + 1)%SpaceDim];
  
  int lastRow  = hiEnd[ a_marchDirection              ];
  int lastCol  = hiEnd[(a_marchDirection + 1)%SpaceDim];

  FArrayBox currLo(validBox,SpaceDim);
  FArrayBox currHi(validBox,SpaceDim);
  currLo.setVal(LARGEREALVAL);
  currHi.setVal(LARGEREALVAL);
  
  // initialize boundary values for currLo and currHi
  currLo.copy(a_curr);
  currHi.copy(a_curr);
  Real dist;
  Real prevDist;
  RealVect loPt;
  RealVect direction;
  Vector<Real> angles;
  angles.resize(4);
  cornerAngles(angles       ,
               currLo       ,
               a_blockNumber);
#if 1
  for (int iRow = firstRow + 1; iRow < lastRow; ++iRow)
    {
      for (int iCol = firstCol + 1; iCol <= lastCol; ++iCol)
        {
          IntVect iv{iCol,iRow};
          
          // (possibly interpolated) flux value for for this column
          Real targetFlux = getTargetFlux(iv        ,
                                          currLo    ,
                                          field_data);

          if (iCol < lastCol)
            {
              // flux on the edge where we begin marching
              IntVect firstColIv(firstCol,iRow);
              RealVect firstColPt(currLo(firstColIv,xDir),currLo(firstColIv,yDir));
              Real firstColFlux = field_data.getMagneticFluxFromDCT(firstColPt);
              
              // find the column where the flux crosses the target flux
              int baseCol;
              for (baseCol = firstCol; baseCol <= lastCol; ++baseCol)
                {
                  // get flux in this column
                  IntVect baseColIv(baseCol,iRow);
                  RealVect baseColPt(currLo(baseColIv,xDir),currLo(baseColIv,yDir));
                  Real baseColFlux = field_data.getMagneticFluxFromDCT(baseColPt);
                  
                  if (((targetFlux - baseColFlux)*(targetFlux - firstColFlux)) < 0.0)
                    {
                      break;
                    }
                }
              if (baseCol == lastCol + 1)
                {
                  pout()<< "can't bracket flux at iv = "<<iv<<" and iBlock = "<<a_blockNumber<<endl;
                  MayDay::Abort("must bracket flux");
                }

              // calculate a sufficiently far distance to look for a crossing point
              IntVect loIv (iCol - 1,iRow);
              IntVect hiIv (baseCol ,iRow);
              
              RealVect loPt  (currLo(loIv,xDir),currLo(loIv,yDir));
              RealVect loHiPt(currLo(hiIv,xDir),currLo(hiIv,yDir));
              
              dist =  (loHiPt[xDir] - loPt[xDir])*(loHiPt[xDir] - loPt[xDir]); 
              dist += (loHiPt[yDir] - loPt[yDir])*(loHiPt[yDir] - loPt[yDir]);
              dist = sqrt(dist);
              
              // stretch in case the path is curved
              dist *= 2.0;
              prevDist = dist;
            }
          else
            {
              dist = prevDist;
            }

          Real subDist = dist/(Real)(a_numSubSteps); 
          
          RealVect fieldDir;
          IntVect baseIv(iCol - 1,iRow);
          RealVect basePt(currLo(baseIv,xDir),currLo(baseIv,yDir));
          RealVect hiPt = basePt;

          bool success;
          for (int iStep = 0; iStep < a_numSubSteps; ++iStep)
            {
              success = false;
              fieldDir = field_data.getFieldUnitVectorFromDCT(basePt);

              RealVect normalDir;
              normalDir[0] = -fieldDir[1];
              normalDir[1] =  fieldDir[0];
                                 
              bool goingLower = false;
              
 
#if 1  
              pout()<< "iv = "<<iv<<", iStep = "<< iStep;
              pout()<< ", normalDir = "<<normalDir;
              pout()<< ", direction = "<<direction<<endl;
              pout()<< "basePt = "<<basePt<<endl;
#endif              
              if (a_blockNumber < 8)
                {
                  direction =  interpolateAngle(currLo,
                                            iv        , 
                                            firstRow  ,
                                            lastRow   ,
                                            goingLower);
                }
              else
                {
                  direction = normalDir;
                }
              hiPt -= subDist*direction;

              // check whether points bracket target flux
              Real loFlux = field_data.getMagneticFluxFromDCT(basePt);
              Real hiFlux = field_data.getMagneticFluxFromDCT(hiPt  );

              Real loDiff = loFlux - targetFlux;
              Real hiDiff = hiFlux - targetFlux;

              if (loDiff*hiDiff < 0.0)
                {
                  success = true;
                  break;
                }
              else
                {
                  basePt = hiPt;
                }
            }
          
          if (!success)
            {
              pout()<<"Points don't bracket flux for currLo, iv = "<<iv<<"in block "<<a_blockNumber<<endl;
              MayDay::Abort("points don't bracket flux");
            } 
                                                
          RealVect newPt = findRoot(basePt       ,
                                    hiPt         ,
                                    targetFlux   ,
                                    field_data   ,
                                    a_tolerance  ,
                                    a_blockNumber,
                                    iv           );
                    
          currLo(iv,xDir) = newPt[xDir];
          currLo(iv,yDir) = newPt[yDir];
        }
    }
  #endif  
#if 1
  for (int iRow = firstRow + 1; iRow < lastRow ; ++iRow)
    {
      for (int iCol = lastCol -1; iCol >= firstCol; --iCol)
        {
          IntVect iv{iCol,iRow};
          
          // (possibly interpolated) flux value for for this column
          Real targetFlux = getTargetFlux(iv        ,
                                          currHi    ,
                                          field_data);
          if (iCol > firstCol)
            {
              // flux on the edge where we begin marching
              IntVect lastColIv(lastCol,iRow);
              RealVect lastColPt(currHi(lastColIv,xDir),currHi(lastColIv,yDir));
              Real lastColFlux = field_data.getMagneticFluxFromDCT(lastColPt);
              
              // find the column where the flux crosses the target flux
              int baseCol;
              for (baseCol = lastCol; baseCol >= firstCol; --baseCol)
                {
                  // get flux in this column
                  IntVect baseColIv(baseCol,iRow);
                  RealVect baseColPt(currHi(baseColIv,xDir),currHi(baseColIv,yDir));
                  Real baseColFlux = field_data.getMagneticFluxFromDCT(baseColPt);
                  
                  if (((targetFlux - baseColFlux)*(targetFlux - lastColFlux)) < 0.0)
                    {
                      break;
                    }
                }
              if (baseCol == firstCol - 1)
                {
                  pout()<< "can't bracket flux at iv = "<<iv<<" and iBlock = "<<a_blockNumber<<endl;
                  FILE* fd;
                  fd = fopen ("done","w");
                  if (fd == NULL) 
                    {
                      perror ("Error opening file");
                    }
                  fclose(fd);
                  MayDay::Abort("must bracket flux");
                }

              // calculate a sufficiently far distance to look for a crossing point
              IntVect loIv (iCol + 1,iRow);
              IntVect hiIv (baseCol ,iRow);
              
              RealVect loPt  (currHi(loIv,xDir),currHi(loIv,yDir));
              RealVect loHiPt(currHi(hiIv,xDir),currHi(hiIv,yDir));
              
              dist = (loHiPt[xDir] - loPt[xDir])*(loHiPt[xDir] - loPt[xDir]);
              dist     += (loHiPt[yDir] - loPt[yDir])*(loHiPt[yDir] - loPt[yDir]);
              dist = sqrt(dist);
              
              // stretch in case the path is curved
              dist *= 2;
              prevDist = dist;
            }
          else
            {
              dist = prevDist;
            }
          
          Real subDist = dist/((Real)a_numSubSteps); 
          
          RealVect fieldDir;
          IntVect baseIv(iCol + 1,iRow);
          RealVect basePt(currHi(baseIv,xDir),currHi(baseIv,yDir));
          RealVect hiPt = basePt;
          
          bool success;
          for (int iStep = 0; iStep < a_numSubSteps; ++iStep)
            {
              success = false;
              fieldDir = field_data.getFieldUnitVectorFromDCT(basePt);

              RealVect normalDir;
              normalDir[0] = -fieldDir[1];
              normalDir[1] =  fieldDir[0];
              normalDir *= -1.0;

              bool goingLower = true;
              if (a_blockNumber < 8)
                {
                  direction =  interpolateAngle(currLo,
                                            iv        , 
                                            firstRow  ,
                                            lastRow   ,
                                            goingLower);
                }
              else
                {
                  direction = normalDir;
                }

#if 0
              pout()<< "iRow = "<<iRow<<", iCol = "<< iCol<<", iStep = "<< iStep;
              pout()<< ", normalDir = "<<normalDir<<endl;
              pout()<< "basePt = "<<basePt<<endl;
#endif
              hiPt -= subDist*direction;
              
              // check whether points bracket target flux
              Real loFlux = field_data.getMagneticFluxFromDCT(basePt);
              Real hiFlux = field_data.getMagneticFluxFromDCT(hiPt  );

              Real loDiff = loFlux - targetFlux;
              Real hiDiff = hiFlux - targetFlux;
              if (loDiff*hiDiff < 0.0)
                {
                  success = true;
                  break;
                }
              else
                {
                  basePt = hiPt;
                }
            }
          if (!success)
            {
              pout()<<"Points don't bracket flux for currHi, iv =  "<<iv<<"in block "<<a_blockNumber<<endl;
              FILE* fd;
              fd = fopen ("done","w");
              if (fd == NULL) 
                {
                  perror ("Error opening file");
                }
              fclose(fd);
              MayDay::Abort("points don't bracket flux");
            }
          RealVect newPt = findRoot(basePt       ,
                                    hiPt         ,
                                    targetFlux   ,
                                    field_data   ,
                                    a_tolerance  ,
                                    a_blockNumber,
                                    iv           );
                    
          currHi(iv,xDir) = newPt[xDir];
          currHi(iv,yDir) = newPt[yDir];
        }
    }
#endif
   for (int iRow = 1; iRow < lastRow; ++iRow)
    {
      for (int iCol = 0; iCol <=lastCol; ++iCol)
        {
          IntVect iv{iCol,iRow};
          Real frac = ((Real)iv[xDir])/((Real)lastCol);
          //frac = 1.0;
          a_curr(iv,xDir) = (1.0 - frac)*currLo(iv,xDir) + frac*currHi(iv,xDir);
          a_curr(iv,yDir) = (1.0 - frac)*currLo(iv,yDir) + frac*currHi(iv,yDir);
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
  else if (a_block == 8)
    {
      retval = LARGEREALVAL*IntVect::Unit;
      a_hasXPt =  false;
    }
  else if (a_block == 9)
    {
      retval = LARGEREALVAL*IntVect::Unit;
      a_hasXPt =  false;
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
      retval = false;
    }
  else if (a_blockName == "mcsol")  
    {
      retval = false;
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
  
  // block 0
  if (a_needGhostForBlock == lcore)
    {
      if (a_needGhostForFace == bottom)
        {
          a_getGhostFromBlockNumber = rcore;
          a_getGhostFromFace        = top  ;
          a_useExtrapolation        = false;
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
          a_getGhostFromBlockNumber = lcore ;
          a_getGhostFromFace        = bottom;
          a_useExtrapolation        = false ;
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
          a_getGhostFromBlockNumber = rcsol;
          a_getGhostFromFace        = top  ;
          a_useExtrapolation        = false;
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
          a_getGhostFromBlockNumber = lcsol ;
          a_getGhostFromFace        = bottom;
          a_useExtrapolation        = false ;
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
      
      for (int m = 0; m < order; ++m)
        {
          for(int j = 0; j < order; ++ j)
            {
              if (j != m)
                {
                  basisPoly[j]  *= (a_extrapPt - a_xPts[m])/(a_xPts[j] - a_xPts[m]); 
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

void extendMapDir(FArrayBox         & a_extendedMap        ,
                  const bool        & a_cornerCase         ,
                  const int         & a_dir                ,     
                  const bool          a_useExtrapolation[2],
                  const Box         & a_stencilBox         ,
                  const int         & a_numGhost           ,
                  const RealVect    & a_dx                 ,  
                  const int         & a_order              ,
                  const std::string & a_blockName          ,
                  const Real        & a_weight             )
{
  for (SideIterator sit; sit.ok(); ++sit)
    {
      if (a_useExtrapolation[sit()])
        {  
          Box validBox = a_stencilBox;
                    
          for (int iGhost = 0; iGhost < a_numGhost; ++iGhost)
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
                   a_blockName     ,
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

void firstDerivCenteredDiff(Real            & a_xDeriv       ,
                            Real            & a_yDeriv       ,
                            const IntVect   & a_iv           ,
                            const int       & a_dir          ,
                            const RealVect  & a_dx           ,
                            const FArrayBox & a_curr         ,
                            const bool      & a_writeInfo    )
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

void setBField(FArrayBox         & a_BField          ,
               const FArrayBox   & a_newPhysicalCoord,
               const RealVect    & a_dx              ,
               const std::string & a_geometryFileName,
               const std::string & a_fieldFileName   , 
               const std::string & a_blockName       )
{
  // magnetic flux is the last comp of BField
  int fluxComp = a_BField.nComp() - 3;
  int fluxDiffHi = fluxComp + 1;
  int fluxDiffLo = fluxComp + 2;
  int xComp    = 0;
  int yComp    = 1;
  int xiDir    = 0;
  int etaDir   = 1;

  // Construct the field data object
  FieldData field_data(a_geometryFileName,
                       a_fieldFileName   ,
                       a_blockName       );
  
  // normalized B field and four first partials  
#ifdef  USE_DCT_FIELD
  field_data.getFieldUnitVectorFromDCT(a_newPhysicalCoord,
                                       a_BField          );
  
  field_data.convertToMappedDerivatives(a_newPhysicalCoord,
                                        a_BField          );
#else
  field_data.getFieldUnitVectorFromMappingFile(a_newPhysicalCoord,
                                               a_BField          );
#endif

  Box box = a_BField.box();
  for (BoxIterator bit(box);bit.ok();++bit) 
    {
      IntVect iv = bit();
      Real mag = 0.0;
      for (int idir = 0; idir < SpaceDim; ++idir) 
        {
          mag += a_BField(iv,idir)*a_BField(iv,idir);
        }
#if 0
      if (mag != 0.0)
        {
          mag = sqrt(mag);
          pout()<<"mag = "<<mag<<endl;
          for (int iComp = 0; iComp < SpaceDim; ++iComp)
            {
              a_BField(iv,iComp) *= 1/mag;
            }
        }

      if (mag != 1.0)
        {
           pout()<<"BField not normalized at iv = "<<iv <<endl;
           pout()<<"magnitude                   = "<<mag<<endl;
        }
      
#endif     
      // set magnetic flux
      Real x = a_newPhysicalCoord(iv,xComp);
      Real y = a_newPhysicalCoord(iv,yComp);
      RealVect pt(x,y);
 
      // iv at high anchor point
      int hiCell = box.bigEnd()  [etaDir];
      IntVect hiIv(iv[xiDir],hiCell);
      
      // physical coord at high anchor point
      Real xHi = a_newPhysicalCoord(hiIv,xComp);
      Real yHi = a_newPhysicalCoord(hiIv,yComp);
      RealVect ptHi(xHi,yHi);
      
      // flux at high anchor point
      Real hiFlux = field_data.getMagneticFluxFromDCT(ptHi);
      
      // iv at low anchor point
      int loCell = box.smallEnd()[etaDir];
      IntVect loIv(iv[xiDir],loCell);
      
      // physical coord at low anchor point
      Real xLo = a_newPhysicalCoord(loIv,xComp);
      Real yLo = a_newPhysicalCoord(loIv,yComp);
      RealVect ptLo(xLo,yLo);

      // flux at low anchor point
      Real loFlux = field_data.getMagneticFluxFromDCT(ptLo);
      
      // assign flux val and flux differences
      a_BField(iv,fluxComp  ) = field_data.getMagneticFluxFromDCT(pt);
      a_BField(iv,fluxDiffHi) = hiFlux - a_BField(iv,fluxComp); 
      a_BField(iv,fluxDiffLo) = loFlux - a_BField(iv,fluxComp); 
    }
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

void computeMapPartialsBdIv(FArrayBox       & a_mapPartials, 
                            const FArrayBox & a_curr       ,
                            const RealVect  & a_dx         ,
                            const bool      & a_BFieldCase ) 
{
  // macros for components of the mapping
  int xDir = 0;
  int yDir = 1;

  Box domainBox = a_curr.box();

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

void  setPerpField(FArrayBox         & a_perpField       ,
                   const FArrayBox   & a_newPhysicalCoord,  
                   const FArrayBox   & a_BField          ,
                   const RealVect    & a_dx              ,
                   const std::string & a_geometryFileName,
                   const std::string & a_fieldFileName   ,
                   const std::string & a_blockName       ,
                   const bool        & a_writeInfo       )
{
  Box box = a_BField.box();
  int xDir = 0;
  int yDir = 1;
   
  // make a normal vector: that is, switch the components, negate one component, and normalize
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      // IntVect
      const IntVect& iv = bit();
      Real b0 = a_BField(iv,xDir);
      Real b1 = a_BField(iv,yDir);
      Real norm = sqrt(b0*b0 + b1*b1);
      
      if (norm > 0.0)
        {
          // perpField = normal vector 
          a_perpField(iv,xDir) = -b1/norm;
          a_perpField(iv,yDir) =  b0/norm;
        }
      else
        {
          a_perpField(iv,xDir) = LARGEREALVAL;
          a_perpField(iv,yDir) = LARGEREALVAL;
        }
    }

  //  iterate over the cells with a neighbor in the x-direction
  Box intBox = box;
  intBox.grow(-1);
  for (BoxIterator bit(intBox); bit.ok(); ++bit)
    {
      // IntVect
      const IntVect& iv = bit();
            
      // partials of the first and second components of the mapping
      Real xDeriv = LARGEREALVAL;    
      Real yDeriv = LARGEREALVAL;    
      
      // four first partials of the vector field perpendicular to the magnetic field lines 
      for (int idir = 0; idir <SpaceDim; ++idir)
        {
          // idir denotes a direction in (xi,eta) space.
          firstDerivCenteredDiff(xDeriv     ,
                                 yDeriv     ,
                                 iv         ,
                                 idir       ,
                                 a_dx       ,
                                 a_perpField,
                                 a_writeInfo);
          
          // component-0            first  component of perp field
          // component-1            second component of perp field
          
          // component-2 partial of first  component with respect to xi
          // component-3 partial of second component with respect to xi
          
          // component-4 partial of first  component with respect to eta
          // component-5 partial of second component with respect to eta
          a_perpField(iv,xDir + (1 + idir)*SpaceDim) = xDeriv;
          a_perpField(iv,yDir + (1 + idir)*SpaceDim) = yDeriv;
        }
    }

  //fix bug in this function. It seems wrongly to set the first two components. 
#if 0  
  // set bBield deriv on boundary the data and derivatives are both in a_perpField. Hence the repeated arguement
  bool BFieldCase = true;
  computeMapPartialsBdIv(a_perpField, 
                         a_perpField,
                         a_dx       ,
                         BFieldCase );
#endif
}

void computeBxBy(Real       & a_Bx  ,
                  Real      & a_By  ,
                  const Real& a_BXi ,
                  const Real& a_BEta,
                  const Real& a_xXi ,
                  const Real& a_yXi ,
                  const Real& a_xEta,
                  const Real& a_yEta)
 {
   // matrix equation for Bx and By
   // |xXi  yXi | |Bx|   |BXi |
   // |         | |  | = |    |
   // |XEta YEta| |By|   |BEta|

   // shorter names
   Real a = a_xXi;
   Real b = a_yXi;
   Real c = a_xEta;
   Real d = a_yEta;

   //discrimant
   Real disc = a*d-b*c;

   // discriminant zero is not possible in a valid calculation (non-zero Jacobian)
   if (disc == 0.0)
     {
       MayDay::Abort("Discriminant zero in computeBxBy");
     }
   
   // invert the matrix
   // |Bx|    1   |d -b| |BXi | 
   // |  | = ---- |    | |    |
   // |By|   disc |-c a| |BEta|

   a_Bx =  d*a_BXi - b*a_BEta;
   a_By = -c*a_BXi + a*a_BEta;

   // complete the inversion
   a_Bx /= disc;
   a_By /= disc;
 } 

void setdXdY(FArrayBox        & a_perpXperpY ,
              const FArrayBox & a_perpField  ,
              const FArrayBox & a_mapPartials,
              const RealVect  & a_dx         )
{

  Box box = a_perpXperpY.box();
  
  int xDir = 0;
  int yDir = 1;
  // iterate over the interior cells
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      
      // Component 0: First component  of unit vector pointing in the magnetic field direction
      // Component 1: Second component of unit vector pointing in the magnetic field direction
      // Component 2: derivative wrt xi[0] of the first component of target
      // Component 3: derivative wrt xi[1] of the first component of target
      // Component 4: derivative wrt xi[0] of the second component of target
      // Component 5: derivative wrt xi[1] of the second component of target
      
      Real B0Xi   = a_perpField(iv,xDir +   SpaceDim);
      Real B1Xi   = a_perpField(iv,yDir +   SpaceDim);
      Real B0Eta  = a_perpField(iv,xDir + 2*SpaceDim);
      Real B1Eta  = a_perpField(iv,yDir + 2*SpaceDim);
      
      // xi partials of the mapping
      Real xXi  = a_mapPartials(iv,xDir             );
      Real yXi  = a_mapPartials(iv,yDir             );
      Real xEta = a_mapPartials(iv,xDir +   SpaceDim);
      Real yEta = a_mapPartials(iv,yDir +   SpaceDim);
      
      // compute the derivative of B0 with repect to X and Y
      Real B0x = LARGEREALVAL;
      Real B0y = LARGEREALVAL;
      computeBxBy(B0x  ,
                  B0y  ,
                  B0Xi ,
                  B0Eta,
                  xXi  ,
                  yXi  ,
                  xEta ,
                  yEta );
      
      // compute the derivative of B1 with repect to X and Y
      Real B1x = LARGEREALVAL;
      Real B1y = LARGEREALVAL;
      computeBxBy(B1x  ,
                  B1y  ,
                  B1Xi ,
                  B1Eta,
                  xXi  ,
                  yXi  ,
                  xEta ,
                  yEta );
      
      a_perpXperpY(iv,xDir)            = B0x;
      a_perpXperpY(iv,yDir)            = B0y;
      a_perpXperpY(iv,xDir + SpaceDim) = B1x;
      a_perpXperpY(iv,yDir + SpaceDim) = B1y;
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

void computeBdSmoothDiagonalTerm(FArrayBox                & a_smoothDiag ,
                                 const RealVect           & a_dx         ,
                                 const FArrayBox          & a_weights    ,
                                 const FArrayBox          & a_Jacobian   ,
                                 const FArrayBox          & a_mapPartials,
                                 const Vector<FArrayBox*> & a_tangent    ,
                                 const Vector<int>        & a_tanDir     )
{
  int xDir = 0;
  int yDir = 1;
  
  int xiEtaComp        = 0;
  //int smoothVolVarComp = 1;
 
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
         
          Real xiSWeight  = a_weights(iv,xiEtaComp);
          Real etaSWeight = 1.0 - xiSWeight;
      
          Real alpha = xiSWeight;
          Real beta  = etaSWeight;

          // Real J = xXi*yEta - xEta*yXi;
                   
          Real J  = a_Jacobian(iv,0);   
          Real J3 = J*J*J           ;

          Real xXi    = a_mapPartials(iv,xDir + 0*SpaceDim); //0
          Real yXi    = a_mapPartials(iv,yDir + 0*SpaceDim); //1
          
          Real xEta   = a_mapPartials(iv,xDir + 1*SpaceDim); //2
          Real yEta   = a_mapPartials(iv,yDir + 1*SpaceDim); //3

          RealVect vXi (xXi ,yXi );
          RealVect vEta(xEta,yEta);

          RealVect vX  (xXi ,xEta);
          RealVect vY  (yXi ,yEta);
          
          Real vXi2  = vXi .dotProduct(vXi );
          Real vEta2 = vEta.dotProduct(vEta);
          
          Real vX2   = vX.dotProduct(vX);
          //Real vXvY  = vX.dotProduct(vY);
          Real vY2   = vY.dotProduct(vY);
          
          Real xXi_D  = hiLo*1.5/a_dx[xDir];
          Real yXi_D  = hiLo*1.5/a_dx[xDir];
          Real xEta_D = hiLo*1.5/a_dx[yDir];
          Real yEta_D = hiLo*1.5/a_dx[yDir];

          // iEdge = 1 and 3
          // SxXi = alpha*(2*xXi*J - JxXi*(xXi^2 + xEta^2 + yXi^2 + yEta^2)))/ J^2 
          // SyXi = alpha*(2*yXi*J - JyXi*(xXi^2 + xEta^2 + yXi^2 + yEta^2)))/ J^2 
          
          // linearize boundary operator to extract diagonal term
          Real SxXi_D = xXi_D*( 2.0*vEta2*vY2 /J3);
          Real SyXi_D = yXi_D*( 2.0*vEta2*vX2 /J3);
         
          // iEdges = 0 and 2 
          // SxEta = beta*(2*xEta*J - JxEta*(xXi^2 + xEta^2 + yXi^2 + yEta^2)))/ J^2 
          // SyEta = beta*(2*yEta*J - JyEta*(xXi^2 + xEta^2 + yXi^2 + yEta^2)))/ J^2 

          // linearize boundary operator to extract diagonal term
          Real SxEta_D = xEta_D*( 2.0*vXi2*vY2 /J3);
          Real SyEta_D = yEta_D*( 2.0*vXi2*vX2 /J3);          
          
          // left or right edges
          if ((iEdge == 1) || (iEdge == 3))
            {
              a_smoothDiag(iv,xDir) = alpha*SxXi_D;
              a_smoothDiag(iv,yDir) = alpha*SyXi_D;
            }
          
          // top or bottom edges
          else
            {
              a_smoothDiag(iv,xDir) = beta*SxEta_D;
              a_smoothDiag(iv,yDir) = beta*SyEta_D;
            }
        }
    } 
}

void computeBdOneOverJDiagonalTerm(FArrayBox                & a_smoothDiag ,
                                   const RealVect           & a_dx         ,
                                   const FArrayBox          & a_weights    ,
                                   const FArrayBox          & a_Jacobian   ,
                                   const FArrayBox          & a_mapPartials,
                                   const Vector<FArrayBox*> & a_tangent    ,
                                   const Vector<int>        & a_tanDir     )
{
  int xDir = 0;
  int yDir = 1;
  
  int xiEtaComp        = 0;
  //int smoothVolVarComp = 1;
 
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
         
          Real xiSWeight  = a_weights(iv,xiEtaComp);
          Real etaSWeight = 1.0 - xiSWeight;
      
          Real alpha = xiSWeight;
          Real beta  = etaSWeight;

          // Real J = xXi*yEta - xEta*yXi;
                   
          Real J  = a_Jacobian(iv,0);   
          Real J3 = J*J*J           ;

          Real xXi    = a_mapPartials(iv,xDir + 0*SpaceDim); //0
          Real yXi    = a_mapPartials(iv,yDir + 0*SpaceDim); //1
          
          Real xEta   = a_mapPartials(iv,xDir + 1*SpaceDim); //2
          Real yEta   = a_mapPartials(iv,yDir + 1*SpaceDim); //3
          
          Real xXi_D  = hiLo*1.5/a_dx[xDir];
          Real yXi_D  = hiLo*1.5/a_dx[xDir];
          Real xEta_D = hiLo*1.5/a_dx[yDir];
          Real yEta_D = hiLo*1.5/a_dx[yDir];

          // iEdge = 1 and 3
          // LxXi = alpha*(-JxXi/J^2) 
          // LyXi = alpha*(-JyXi/J^2)
          
          // linearize boundary operator to extract diagonal term
          Real LxXi_D = xXi_D*4.0*yEta*yEta/J3;
          Real LyXi_D = yXi_D*4.0*xEta*xEta/J3;
         
          // iEdges = 0 and 2 
          // LxEta = beta*(-JxEta/J^2) 
          // LyEta = beta*(-JyEta/J^2) 

          // linearize boundary operator to extract diagonal term
          Real LxEta_D = xEta_D*4.0*yXi*yXi/J3;
          Real LyEta_D = yEta_D*4.0*xXi*xXi/J3;          
          
          // left or right edges
          if ((iEdge == 1) || (iEdge == 3))
            {
              a_smoothDiag(iv,xDir) = alpha*LxXi_D;
              a_smoothDiag(iv,yDir) = alpha*LyXi_D;
            }
          
          // top or bottom edges
          else
            {
              a_smoothDiag(iv,xDir) = beta*LxEta_D;
              a_smoothDiag(iv,yDir) = beta*LyEta_D;
            }
        }
    } 
}

void computeBdOneOverSDiagonalTerm(FArrayBox                & a_smoothDiag ,
                                   const RealVect           & a_dx         ,
                                   const FArrayBox          & a_weights    ,
                                   const FArrayBox          & a_mapPartials,
                                   const Vector<FArrayBox*> & a_tangent    ,
                                   const Vector<int>        & a_tanDir     )
{
  int xDir = 0;
  int yDir = 1;
  
  int xiEtaComp        = 0;
  //int smoothVolVarComp = 1;
 
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
         
          Real xiSWeight  = a_weights(iv,xiEtaComp);
          Real etaSWeight = 1.0 - xiSWeight;
      
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
          Real S2 = S*S;
          Real S4 = S2*S2;

          // iEdge = 1 and 3
          // LxXi = alpha*(-2.0*x_xi/S^2) 
          // LyXi = alpha*(-2.0*y_xi/S^2)
          
          // linearize boundary operator to extract diagonal term
          Real LxXi_D = x_xi_D*(-2.0*S + 8.0*x_xi*x_xi*S)/S4;
          Real LyXi_D = y_xi_D*(-2.0*S + 8.0*y_xi*y_xi*S)/S4;
         
          // iEdges = 0 and 2 
          // LxEta = beta*(-2.0*x_eta/S^2) 
          // LyEta = beta*(-2.0*y_eta/S^2) 

          // linearize boundary operator to extract diagonal term
          Real LxEta_D = x_eta_D*(-2.0*S + 8.0*x_eta*x_eta*S)/S4;
          Real LyEta_D = y_eta_D*(-2.0*S + 8.0*y_eta*y_eta*S)/S4;          
          
          // left or right edges
          if ((iEdge == 1) || (iEdge == 3))
            {
              a_smoothDiag(iv,xDir) = alpha*LxXi_D;
              a_smoothDiag(iv,yDir) = alpha*LyXi_D;
            }
          
          // top or bottom edges
          else
            {
              a_smoothDiag(iv,xDir) = beta*LxEta_D;
              a_smoothDiag(iv,yDir) = beta*LyEta_D;
            }
        }
    } 
}

void computeBdSmJnDiagonalTerm(FArrayBox                & a_smoothDiag ,
                               const RealVect           & a_dx         ,
                               const FArrayBox          & a_weights    ,
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
  
  int xiEtaComp = 0;
 
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
         
          Real xiSWeight  = a_weights(iv,xiEtaComp);
          Real etaSWeight = 1.0 - xiSWeight        ;
      
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

void computeSmoothDiagonalTerm(FArrayBox       & a_smoothDiag ,
                               const RealVect  & a_dx         ,
                               const FArrayBox & a_weights    ,
                               const FArrayBox & a_Jacobian   ,
                               const FArrayBox & a_mapPartials)
{
  // indexing macros
  int firstEq = 0;
  int secEq   = 1;

  int xiEtaComp        = 0;

  // indexing macros
  int xDir = 0;
  int yDir = 1;
  
  //int xiDir  = 0;
  //int etaDir = 1;
  
  // iterate
  Box box = a_smoothDiag.box();
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      const IntVect iv = bit(); 

      Real xiSWeight  = a_weights(iv,xiEtaComp);
      Real etaSWeight = 1.0 - xiSWeight;
      
      Real alpha = xiSWeight;
      Real beta  = etaSWeight;
      
      Real dx = a_dx[xDir];
      Real dy = a_dx[yDir];
      
      Real dx2 = dx*dx;
      Real dy2 = dy*dy;
      
      Real xXi    = a_mapPartials(iv,xDir + 0*SpaceDim); //0
      Real yXi    = a_mapPartials(iv,yDir + 0*SpaceDim); //1
      
      Real xEta   = a_mapPartials(iv,xDir + 1*SpaceDim); //2
      Real yEta   = a_mapPartials(iv,yDir + 1*SpaceDim); //3
      
     
      Real J  = a_Jacobian(iv,0); 
      Real J3 = J*J*J           ;
           
      RealVect vXi (xXi ,yXi );
      RealVect vEta(xEta,yEta);
      
      RealVect vX  (xXi ,xEta);
      RealVect vY  (yXi ,yEta);
      
      Real vXi2  = vXi .dotProduct(vXi );
      Real vEta2 = vEta.dotProduct(vEta);
      
      Real vX2   = vX.dotProduct(vX);
      //Real vXvY  = vX.dotProduct(vY);
      Real vY2   = vY.dotProduct(vY);
      
      //Lagrangian = S/J
      
      // contributions to Euler-Lagrange equations
      //S = alpha*xiS + beta*etaS;
      
      //LxXi  = (SxXi* J - S*JxXi )/J2;
      //LyXi  = (SyXi* J - S*JyXi )/J2;
      //LxEta = (SxEta*J - S*JxEta)/J2;
      //LxEta = (SyEta*J - S*JyEta)/J2;
      
      // terms in Euler-Lagrange equation
      //Real LxXiXi   = ((SxXiXi*J    + SxXi*Jxi      - Sxi*JxXi    - S*JxXiXi   )*J2 - 2*J*Jxi   *(SxXi*J   - S*JxXi ))/J4;
      //Real LyXiXi   = ((SyXiXi  *J  + SyXi*Jxi      - Sxi *JyXi   - S*JyXiXi   )*J2 - 2*J*Jxi   *(SyXi*J   - S*JyXi ))/J4;
      //Real LxEtaEta = ((SxEtaEta*J  + SxEta*Jeta    - Seta*JxEta  - S*JxEtaEta )*J2 - 2*J*Jeta  *(SxEta*J - S*JxEta)) /J4;
      //Real LyEtaEta = ((SyEtaEta*J  + SyEta*Jeta    - Seta*JyEta  - S*JyEtaEta )*J2 - 2*J*Jeta  *(SyEta*J - S*JyEta)) /J4;         
     
      // Diagonal contribution of linearized operator

      Real xXiXi_D   = -2.0/dx2;
      Real xEtaEta_D = -2.0/dy2;
      
      Real yXiXi_D   = -2.0/dx2;
      Real yEtaEta_D = -2.0/dy2;
      
      Real LxXiXi_D   = alpha*xXiXi_D  *2.0*vEta2*vY2/J3;
      Real LxEtaEta_D = beta *xEtaEta_D*2.0*vXi2 *vY2/J3;
      
      Real LyXiXi_D   = alpha*yXiXi_D  *2.0*vEta2*vX2/J3;
      Real LyEtaEta_D = beta *yEtaEta_D*2.0*vXi2 *vX2/J3;

      // assign data 
      a_smoothDiag(iv,firstEq) = -(LxXiXi_D + LxEtaEta_D); 
      a_smoothDiag(iv,secEq  ) = -(LyXiXi_D + LyEtaEta_D);
    }
}

void computeOneOverJDiagonalTerm(FArrayBox       & a_smoothDiag ,
                                 const RealVect  & a_dx         ,
                                 const FArrayBox & a_weights    ,
                                 const FArrayBox & a_Jacobian   ,
                                 const FArrayBox & a_mapPartials)
{
  // indexing macros
  int firstEq = 0;
  int secEq   = 1;

  int xiEtaComp  = 0;

  // indexing macros
  int xDir = 0;
  int yDir = 1;
  
  //int xiDir  = 0;
  //int etaDir = 1;
  
  // iterate
  Box box = a_smoothDiag.box();
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      const IntVect iv = bit(); 

      Real xiSWeight  = a_weights(iv,xiEtaComp);
      Real etaSWeight = 1.0 - xiSWeight        ;
      
      Real alpha = xiSWeight;
      Real beta  = etaSWeight;
      
      Real dx = a_dx[xDir];
      Real dy = a_dx[yDir];
      
      Real dx2 = dx*dx;
      Real dy2 = dy*dy;
      
      Real xXi  = a_mapPartials(iv,xDir + 0*SpaceDim); //0
      Real yXi  = a_mapPartials(iv,yDir + 0*SpaceDim); //1
      
      Real xEta = a_mapPartials(iv,xDir + 1*SpaceDim); //2
      Real yEta = a_mapPartials(iv,yDir + 1*SpaceDim); //3
           
      Real J  = a_Jacobian(iv,0); 
      Real J3 = J*J*J           ;
     
      //Lagrangian = 1.0/J

      // LxXi = alpha*(-JxXi/J^2) 
      // LyXi = alpha*(-JyXi/J^2)
      // LxEta = beta*(-JxEta/J^2) 
      // LyEta = beta*(-JyEta/J^2) 
          
      // LxXiXi   = (-yXiEta*J + 2.0*yEta*(xXiXi *yEta + xXi*yxiEta  - xXiEta *yXi - xEta*yXiXi ))/J3
      // LxEtaEta = ( yXiEta*J - 2.0*yXi *(xXiEta*yEta + xXi*yEtaEta - xEtaEta*yXi - xEta*yXiEta))/J3

      // LyXiXi   = (-xXiEta*J + 2.0*xEta*(xXiXi *yEta + xXi*yxiEta  - xXiEta *yXi - xEta*yXiXi ))/J3
      // LyEtaEta = ( xXiEta*J - 2.0*xXi *(xXiEta*yEta + xXi*yEtaEta - xEtaEta*yXi - xEta*yXiEta))/J3

      // Diagonal contribution of linearized operator
      Real xXiXi_D   = -2.0/dx2;
      Real xEtaEta_D = -2.0/dy2;
      
      Real yXiXi_D   = -2.0/dx2;
      Real yEtaEta_D = -2.0/dy2;
       
      Real LxXiXi_D   = alpha*(xXiXi_D  *2.0*yEta*yEta/J3);
      Real LxEtaEta_D = beta *(xEtaEta_D*2.0*yXi *yXi /J3);
      
      Real LyXiXi_D   = alpha*(yXiXi_D  *2.0*xEta*xEta/J3);
      Real LyEtaEta_D = beta *(yEtaEta_D*2.0*xXi *xXi /J3);
      
      // assign data 
      a_smoothDiag(iv,firstEq) = -(LxXiXi_D + LxEtaEta_D); 
      a_smoothDiag(iv,secEq  ) = -(LyXiXi_D + LyEtaEta_D);
    }
}

void computeOneOverSDiagonalTerm(FArrayBox       & a_smoothDiag ,
                                 const RealVect  & a_dx         ,
                                 const FArrayBox & a_weights    ,
                                 const FArrayBox & a_Jacobian   ,
                                 const FArrayBox & a_mapPartials)
{
  // indexing macros
  int firstEq = 0;
  int secEq   = 1;

  int xiEtaComp        = 0;

  // indexing macros
  int xDir = 0;
  int yDir = 1;
  
  //int xiDir  = 0;
  //int etaDir = 1;
  
  // iterate
  Box box = a_smoothDiag.box();
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      const IntVect iv = bit(); 

      Real xiSWeight  = a_weights(iv,xiEtaComp);
      Real etaSWeight = (1.0 - xiSWeight);
      
      Real alpha = xiSWeight;
      Real beta  = etaSWeight;
      
      Real dx = a_dx[xDir];
      Real dy = a_dx[yDir];
      
      Real dx2 = dx*dx;
      Real dy2 = dy*dy;
      
      Real x_xi  = a_mapPartials(iv,xDir + 0*SpaceDim); //0
      Real y_xi  = a_mapPartials(iv,yDir + 0*SpaceDim); //1
      
      Real x_eta = a_mapPartials(iv,xDir + 1*SpaceDim); //2
      Real y_eta = a_mapPartials(iv,yDir + 1*SpaceDim); //3
        
      Real S  = x_xi*x_xi + y_xi*y_xi + x_eta*x_eta + y_eta*y_eta;
      Real S2 = S*S;
      Real S4 = S2*S2;

      //Lagrangian = 1.0/S

      // Lx_xi  = alpha*(-2.0*x_xi/S^2) 
      // Ly_xi  = alpha*(-2.0*y_xi/S^2)
      // Lx_eta = beta *(-2.0*x_eta/S^2) 
      // Ly_eta = beta *(-2.0*y_eta/S^2) 
          
      // Lx_xi_xi   = alpha*(-2.0*x_2xi*S2  - 4.0*x_xi *S*(x_2xi    + y_2xi    + x_eta_xi  + y_eta_xi)/S4
      // Ly_xi_xi   = alpha*(-2.0*y_2xi*S2  - 4.0*y_xi *S*(x_2xi    + y_2xi    + x_eta_xi  + y_eta_xi)/S4 

      // Lx_eta_eta = beta *(-2.0*x_2eta*S2 - 4.0*x_eta*S*(x_xi_eta + y_xi_eta + x_2eta    + y_2eta)  /S4
      // Ly_eta_eta = beta *(-2.0*y_2eta*S2 - 4.0*y_eta*S*(x_xi_eta + y_xi_eta + x_2eta    + y_2eta)  /S4

      // Diagonal contribution of linearized operator
      Real x_2xi_D  = -2.0/dx2;
      Real x_2eta_D = -2.0/dy2;
      
      Real y_2xi_D  = -2.0/dx2;
      Real y_2eta_D = -2.0/dy2;
       
      Real Lx_xi_xi_D = alpha*(x_2xi_D *(-2.0*S2 - 4.0*x_xi*S))/S4;
      Real Ly_xi_xi_D = alpha*(y_2xi_D *(-2.0*S2 - 4.0*y_xi*S))/S4;
      
      Real Lx_eta_eta_D = beta*(x_2eta_D *(-2.0*S2 - 4.0*x_eta*S))/S4;
      Real Ly_eta_eta_D = beta*(y_2eta_D *(-2.0*S2 - 4.0*y_eta*S))/S4;
      
      // assign data 
      a_smoothDiag(iv,firstEq) = -(Lx_xi_xi_D + Lx_eta_eta_D); 
      a_smoothDiag(iv,secEq  ) = -(Ly_xi_xi_D + Ly_eta_eta_D);
    }
}

void computeSmJnDiagonalTerm(FArrayBox       & a_smoothDiag ,
                             const RealVect  & a_dx         ,
                             const FArrayBox & a_weights    ,
                             const FArrayBox & a_Jacobian   ,
                             const FArrayBox & a_mapPartials)
{
  // indexing macros
  int firstEq = 0;
  int secEq   = 1;

  int xiEtaComp        = 0;
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

      Real xiSWeight  = a_weights(iv,xiEtaComp);
      Real etaSWeight = 1.0 - xiSWeight;
      
      Real alpha = xiSWeight;
      Real beta  = etaSWeight;
      
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

void smoothRes(Real            & a_res0      ,
               Real            & a_res1      , 
               const Real      & a_xiSWeight ,
               const Real      & a_etaSWeight,
               const IntVect   & a_iv        ,
               const FArrayBox & a_smoothness,
               const FArrayBox & a_Jacobian  ,
               const bool      & a_writeInfo )      
{
  //indexing macros
  int xiSmooth  = 0;
  int etaSmooth = 1;

  int xiDir  = 0;
  int etaDir = 1;

  int xXiDir = 0;
  int yXiDir = 1;
  
  int xEtaDir = 0;
  int yEtaDir = 1;
  
   //weights
  Real alpha = a_xiSWeight;
  Real beta  = a_etaSWeight;

  // xi smoothing 
  Real xiS        = a_smoothness(a_iv,     7*xiSmooth                       );//0
  Real xiSxXi     = a_smoothness(a_iv, 1 + 7*xiSmooth  +   xXiDir           );//1
  Real xiSyXi     = a_smoothness(a_iv, 1 + 7*xiSmooth  +   yXiDir           );//2
  Real xiSxXiXi   = a_smoothness(a_iv, 1 + 7*xiSmooth  +   SpaceDim + xXiDir);//3
  Real xiSyXiXi   = a_smoothness(a_iv, 1 + 7*xiSmooth  +   SpaceDim + yXiDir);//4
  Real xiSxi      = a_smoothness(a_iv, 1 + 7*xiSmooth  + 2*SpaceDim + xiDir );//5
  Real xiSeta     = a_smoothness(a_iv, 1 + 7*xiSmooth + 2*SpaceDim + etaDir);//6
  
  // eta smoothing
  Real etaS        = a_smoothness(a_iv,     7*etaSmooth                       );//7
  Real etaSxEta    = a_smoothness(a_iv, 1 + 7*etaSmooth +   xEtaDir           );//8
  Real etaSyEta    = a_smoothness(a_iv, 1 + 7*etaSmooth +   yEtaDir           );//9
  Real etaSxEtaEta = a_smoothness(a_iv, 1 + 7*etaSmooth +   SpaceDim + xEtaDir);//10
  Real etaSyEtaEta = a_smoothness(a_iv, 1 + 7*etaSmooth +   SpaceDim + yEtaDir);//11
  Real etaSxi      = a_smoothness(a_iv, 1 + 7*etaSmooth + 2*SpaceDim + xiDir  );//12
  Real etaSeta     = a_smoothness(a_iv, 1 + 7*etaSmooth + 2*SpaceDim + etaDir );//13

  //indexing macros
  int xDir = 0;
  int yDir = 1;

  //jacobian
  Real J  = a_Jacobian(a_iv,0                        );    //0
  Real J2 = J*J                                      ;
  Real J4 = J2*J2                                    ;
  
  // derivatives with respect to derivatives of the mapping
  Real JxXi  = a_Jacobian(a_iv,1 + 0*SpaceDim + xDir  );   //1
  Real JyXi  = a_Jacobian(a_iv,1 + 0*SpaceDim + yDir  );   //2
  Real JxEta = a_Jacobian(a_iv,1 + 1*SpaceDim + xDir  );   //3
  Real JyEta = a_Jacobian(a_iv,1 + 1*SpaceDim + yDir  );   //4

  // derivatives of the above with respect to xi and eta
  Real JxXiXi   = a_Jacobian(a_iv,1 + 2*SpaceDim + xDir  );//5
  Real JyXiXi   = a_Jacobian(a_iv,1 + 2*SpaceDim + yDir  );//6
  Real JxEtaEta = a_Jacobian(a_iv,1 + 3*SpaceDim + xDir  );//7
  Real JyEtaEta = a_Jacobian(a_iv,1 + 3*SpaceDim + yDir  );//8

  // derivatives with respect to xi and eta
  Real Jxi  = a_Jacobian(a_iv,1 + 4*SpaceDim + xiDir );    //9
  Real Jeta = a_Jacobian(a_iv,1 + 4*SpaceDim + etaDir);    //10
  
  Real S = alpha*xiS + beta*etaS;
  
  Real SxXi = alpha*xiSxXi;
  Real SyXi = alpha*xiSyXi;

  Real SxEta = beta*etaSxEta;
  Real SyEta = beta*etaSyEta;
  
  Real SxXiXi   = alpha*xiSxXiXi   ;
  Real SyXiXi   = alpha*xiSyXiXi   ;
  Real SxEtaEta = beta *etaSxEtaEta;
  Real SyEtaEta = beta *etaSyEtaEta;

  Real Sxi  = alpha*xiSxi  + beta*etaSxi ;
  Real Seta = alpha*xiSeta + beta*etaSeta;
  
  //Lagrangian = S/J
  //LxXi  = (SxXi *J - S*JxXi )/J2;
  //LyXi  = (SyXi *J - S*JyXi )/J2;
  //LxEta = (SxEta*J - S*JxEta)/J2;
  //LyEta = (SyEta*J - S*JyEta)/J2;

  Real LxXiXi   = ((SxXiXi  *J + SxXi *Jxi  - Sxi *JxXi  - S*JxXiXi  )*J2 - 2*J*Jxi *(SxXi *J - S*JxXi ))/J4;
  Real LyXiXi   = ((SyXiXi  *J + SyXi *Jxi  - Sxi *JyXi  - S*JyXiXi  )*J2 - 2*J*Jxi *(SyXi *J - S*JyXi ))/J4;
  Real LxEtaEta = ((SxEtaEta*J + SxEta*Jeta - Seta*JxEta - S*JxEtaEta)*J2 - 2*J*Jeta*(SxEta*J - S*JxEta))/J4;
  Real LyEtaEta = ((SyEtaEta*J + SyEta*Jeta - Seta*JyEta - S*JyEtaEta)*J2 - 2*J*Jeta*(SyEta*J - S*JyEta))/J4;         

 // no term in S or J includes x or includes y (explicitly)
  Real Lx = 0;
  Real Ly = 0;

  // return residual
  a_res0 = -(LxXiXi + LxEtaEta) + Lx;
  a_res1 = -(LyXiXi + LyEtaEta) + Ly;
}

void oneOverJRes(Real            & a_res0      ,
                 Real            & a_res1      , 
                 const Real      & a_xiSWeight ,
                 const Real      & a_etaSWeight,
                 const IntVect   & a_iv        ,
                 const FArrayBox & a_smoothness,
                 const FArrayBox & a_Jacobian  ,
                 const bool      & a_writeInfo )      
{
  //indexing macros
  //int xiSmooth  = 0;
  //int etaSmooth = 1;

  int xiDir  = 0;
  int etaDir = 1;

  //int xXiDir = 0;
  //int yXiDir = 1;
  
  //int xEtaDir = 0;
  //int yEtaDir = 1;
  
   //weights
  Real alpha = a_xiSWeight;
  Real beta  = a_etaSWeight;

  // xi smoothing 
  //Real xiS        = a_smoothness(a_iv,     7*xiSmooth                       );//0
  //Real xiSxXi     = a_smoothness(a_iv, 1 + 7*xiSmooth  +   xXiDir           );//1
  //Real xiSyXi     = a_smoothness(a_iv, 1 + 7*xiSmooth  +   yXiDir           );//2
  //Real xiSxXiXi   = a_smoothness(a_iv, 1 + 7*xiSmooth  +   SpaceDim + xXiDir);//3
  //Real xiSyXiXi   = a_smoothness(a_iv, 1 + 7*xiSmooth  +   SpaceDim + yXiDir);//4
  //Real xiSxi      = a_smoothness(a_iv, 1 + 7*xiSmooth  + 2*SpaceDim + xiDir );//5
  //Real xiSeta     = a_smoothness(a_iv, 1 + 7*xiSmooth + 2*SpaceDim + etaDir);//6
  
  // eta smoothing
  //Real etaS        = a_smoothness(a_iv,     7*etaSmooth                       );//7
  //Real etaSxEta    = a_smoothness(a_iv, 1 + 7*etaSmooth +   xEtaDir           );//8
  //Real etaSyEta    = a_smoothness(a_iv, 1 + 7*etaSmooth +   yEtaDir           );//9
  //Real etaSxEtaEta = a_smoothness(a_iv, 1 + 7*etaSmooth +   SpaceDim + xEtaDir);//10
  //Real etaSyEtaEta = a_smoothness(a_iv, 1 + 7*etaSmooth +   SpaceDim + yEtaDir);//11
  //Real etaSxi      = a_smoothness(a_iv, 1 + 7*etaSmooth + 2*SpaceDim + xiDir  );//12
  //Real etaSeta     = a_smoothness(a_iv, 1 + 7*etaSmooth + 2*SpaceDim + etaDir );//13

  //indexing macros
  int xDir = 0;
  int yDir = 1;

  //jacobian
  Real J  = a_Jacobian(a_iv,0                        );    //0
  Real J2 = J*J                                      ;
  Real J4 = J2*J2                                    ;
  
  // derivatives with respect to derivatives of the mapping
  Real JxXi  = a_Jacobian(a_iv,1 + 0*SpaceDim + xDir  );   //1
  Real JyXi  = a_Jacobian(a_iv,1 + 0*SpaceDim + yDir  );   //2
  Real JxEta = a_Jacobian(a_iv,1 + 1*SpaceDim + xDir  );   //3
  Real JyEta = a_Jacobian(a_iv,1 + 1*SpaceDim + yDir  );   //4

  // derivatives of the above with respect to xi and eta
  Real JxXiXi   = a_Jacobian(a_iv,1 + 2*SpaceDim + xDir  );//5
  Real JyXiXi   = a_Jacobian(a_iv,1 + 2*SpaceDim + yDir  );//6
  Real JxEtaEta = a_Jacobian(a_iv,1 + 3*SpaceDim + xDir  );//7
  Real JyEtaEta = a_Jacobian(a_iv,1 + 3*SpaceDim + yDir  );//8

  // derivatives with respect to xi and eta
  Real Jxi  = a_Jacobian(a_iv,1 + 4*SpaceDim + xiDir );    //9
  Real Jeta = a_Jacobian(a_iv,1 + 4*SpaceDim + etaDir);    //10
  
  //Real S = 1.0/J;
    
  //Lagrangian = 1.0/J
  //LxXi  = -JxXi /J2;
  //LyXi  = -JyXi /J2;
  //LxEta = -JxEta/J2;
  //LyEta = -JyEta/J2;

  Real LxXiXi   = alpha*(-JxXiXi*J2 + JxXi*2*J*Jxi)/J4;
  Real LyXiXi   = alpha*(-JyXiXi*J2 + JyXi*2*J*Jxi)/J4;

  Real LxEtaEta = beta*(-JxEtaEta*J2 + JxEta*2*J*Jeta)/J4;
  Real LyEtaEta = beta*(-JyEtaEta*J2 + JyEta*2*J*Jeta)/J4;         

 // no term in S or J includes x or includes y (explicitly)
  Real Lx = 0;
  Real Ly = 0;

  // return residual
  a_res0 = -(LxXiXi + LxEtaEta) + Lx;
  a_res1 = -(LyXiXi + LyEtaEta) + Ly;
}

void oneOverSRes(Real            & a_res0       ,
                 Real            & a_res1       , 
                 const Real      & a_xiSWeight  ,
                 const Real      & a_etaSWeight ,
                 const IntVect   & a_iv         ,
                 const FArrayBox & a_mapPartials,
                 const bool      & a_writeInfo  )      
{
  // indexing macros
  //int xiSmooth  = 0;
  //int etaSmooth = 1;
  
  int xDir = 0;
  int yDir = 1;
  
  // weights
  Real alpha = a_xiSWeight;
  Real beta  = a_etaSWeight;

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

  Real S = x_xi*x_xi + y_xi*y_xi + x_eta*x_eta + y_eta*y_eta;
  //Real S2 = S*S;
  Real S3 = S*S*S;

  //Real L = 1.0/S;
    
  //Lagrangian = 1.0/S
  //Lx_xi  = -2.0*x_xi /S2
  //Ly_xi  = -2.0*y_xi /S2
  //Lx_eta = -2.0*x_eta/S2
  //Ly_eta = -2.0*y_eta/S2


  Real Lx_xi_xi   = alpha *(8.0*(x_xi*x_2xi    + y_xi*y_2xi    + x_eta*x_eta_xi + y_eta*y_eta_xi)*x_xi  - 2*S*x_2xi )/S3;
  Real Ly_xi_xi   = alpha *(8.0*(x_xi*x_2xi    + y_xi*y_2xi    + x_eta*x_eta_xi + y_eta*y_eta_xi)*y_xi  - 2*S*y_2xi )/S3;
  Real Lx_eta_eta = beta  *(8.0*(x_xi*x_xi_eta + y_xi*y_xi_eta + x_eta*x_2eta   + y_eta*y_2eta  )*x_eta - 2*S*x_2eta)/S3;
  Real Ly_eta_eta = beta  *(8.0*(x_xi*x_xi_eta + y_xi*y_xi_eta + x_eta*x_2eta   + y_eta*y_2eta  )*y_eta - 2*S*y_2eta)/S3;
  
  // no term in S includes x or includes y (explicitly)
  Real Lx = 0;
  Real Ly = 0;

  // return residual
  a_res0 = -(Lx_xi_xi + Lx_eta_eta) + Lx;
  a_res1 = -(Ly_xi_xi + Ly_eta_eta) + Ly;
}

void SmJnRes(Real            & a_res0       ,
             Real            & a_res1       , 
             const Real      & a_xiSWeight  ,
             const Real      & a_etaSWeight ,
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
  
  // weights
  Real alpha = a_xiSWeight;
  Real beta  = a_etaSWeight;

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
  Real S  = alpha*(x_xi*x_xi + y_xi*y_xi) + beta*(x_eta*x_eta + y_eta*y_eta);
  
  Real Sx_xi  = alpha*2.0*x_xi ;
  Real Sy_xi  = alpha*2.0*y_xi ;
  Real Sx_eta = beta *2.0*x_eta;
  Real Sy_eta = beta *2.0*y_eta;

  Real Sx_xi_xi   =  alpha*2.0*x_2xi ;
  Real Sy_xi_xi   =  alpha*2.0*y_2xi ;
  Real Sx_eta_eta =  beta *2.0*x_2eta;
  Real Sy_eta_eta =  beta *2.0*y_2eta;

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

void computeFirstDerivCoeff(Real            & a_A       ,
                            Real            & a_B       ,
                            Real            & a_C       ,   
                            Real            & a_alpha   ,
                            Real            & a_beta    ,
                            Real            & a_gamma   ,
                            Real            & a_Jacobian,
                            const IntVect   & a_iv      ,
                            const RealVect  & a_dx      ,
                            const FArrayBox & a_curr    ,
                            const bool      & a_writeInfo)

{
  // indexing macros
  int xDir = 0;
  int yDir = 1;
  
  // stencil ivs
  IntVect xHi = a_iv;
  IntVect xLo = a_iv;
  IntVect yHi = a_iv;
  IntVect yLo = a_iv;
  
  xHi += BASISV(xDir);
  xLo -= BASISV(xDir);
  yHi += BASISV(yDir);
  yLo -= BASISV(yDir);

  Real dx = a_dx[0];
  Real dy = a_dx[1];

  Real x_xi = (a_curr(xHi,xDir) - a_curr(xLo,xDir))/(2.0*dx);
  Real y_xi = (a_curr(xHi,yDir) - a_curr(xLo,yDir))/(2.0*dx);

  Real x_eta  = (a_curr(yHi,xDir) - a_curr(yLo,xDir))/(2.0*dy);
  Real y_eta  = (a_curr(yHi,yDir) - a_curr(yLo,yDir))/(2.0*dy);

  // Jacobian
  Real J = x_xi*y_eta - x_eta*y_xi;

  // smaller than zero jacobian => we're searching outside the space of 1-1, orientation preserving maps.
  if (J <= 0.0)
    {
      if(a_writeInfo)
        {
          pout()<<"iv = "      <<a_iv<<endl;
          pout()<<"Jacobian = "<<J <<endl; 
        }
      MayDay::Warning("Jacobian <= 0.0");    
    }

  // coefficients for the Euler equations
  a_A = x_xi*y_xi + x_eta*y_eta;
  a_B = y_xi*y_xi + y_eta*y_eta;
  a_C = x_xi*x_xi + x_eta*x_eta;
  
  a_alpha = (x_eta *x_eta  + y_eta *y_eta )/(J*J*J);
  a_beta  = (x_xi*x_eta  + y_xi*y_eta )/(J*J*J);
  a_gamma = (x_xi*x_xi + y_xi*y_xi)/(J*J*J);

  // return J for diagnostics
  a_Jacobian = J;
}

void smoothRes(Real            & a_inc0     ,
               Real            & a_inc1     , 
               Real            & a_diagTerm0,
               Real            & a_diagTerm1, 
               const Real      & a_A        ,
               const Real      & a_B        ,
               const Real      & a_C        ,
               const Real      & a_alpha    , 
               const Real      & a_beta     ,  
               const Real      & a_gamma    ,
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
  Real dx = a_dx[xDir];
  Real dy = a_dx[yDir];
  Real dx2 = dx*dx;
  Real dy2 = dy*dy;
  
  // second partials for R and Z
  Real x_2xi  = (a_curr(xHi,xDir) - 2*a_curr(a_iv,xDir) + a_curr(xLo,xDir))/(dx*dx); 
  Real x_2eta = (a_curr(yHi,xDir) - 2*a_curr(a_iv,xDir) + a_curr(yLo,xDir))/(dy*dy);
  
  Real y_2xi  = (a_curr(xHi,yDir) - 2*a_curr(a_iv,yDir) + a_curr(xLo,yDir))/(dx*dx); 
  Real y_2eta = (a_curr(yHi,yDir) - 2*a_curr(a_iv,yDir) + a_curr(yLo,yDir))/(dy*dy);
  
  Real x_xiEta = (a_curr(hiHi,xDir) + a_curr(loLo,xDir) - a_curr(hiLo,xDir) - a_curr(loHi,xDir))/(4*dx*dy);
  Real y_xiEta = (a_curr(hiHi,yDir) + a_curr(loLo,yDir) - a_curr(hiLo,yDir) - a_curr(loHi,yDir))/(4*dx*dy);

  // increment for R
  a_inc0 = +a_B*(a_alpha*x_2xi - 2*a_beta*x_xiEta + a_gamma*x_2eta) 
           -a_A*(a_alpha*y_2xi - 2*a_beta*y_xiEta + a_gamma*y_2eta);

  // increment for Z
  a_inc1 = -a_A*(a_alpha*x_2xi - 2*a_beta*x_xiEta + a_gamma*x_2eta) 
           +a_C*(a_alpha*y_2xi - 2*a_beta*y_xiEta + a_gamma*y_2eta);
  
  // relaxation parameter = (-1/2)*(reciprocal of the diagonal term of the operator)
  a_diagTerm0 = -2.0*a_B*(a_alpha/dx2 + a_gamma/dy2);
  a_diagTerm1 = -2.0*a_C*(a_alpha/dx2 + a_gamma/dy2);
    
  if (a_writeInfo)
    {
      pout()<< "iv     = "     <<a_iv    <<endl;
      pout()<< "smoothRes0 = " <<a_inc0  <<endl;
      pout()<< "smoothRes1 = " <<a_inc1  <<endl;
  
      pout()<<endl;
    }    
}

void computePartial(Real            & a_xPartial,
                    Real            & a_yPartial,
                    const RealVect  & a_dx      ,
                    const IntVect   & a_iv      ,
                    const int       & a_dir     ,
                    const FArrayBox & a_fab     )

{
  // indexing macros
  int xDir = 0;
  int yDir = 1;
  
  // stencil ivs
  IntVect hi = a_iv;
  IntVect lo = a_iv;
  hi += BASISV(a_dir);
  lo -= BASISV(a_dir);
 
  // dx
  Real dx = a_dx[a_dir];

  a_xPartial = (a_fab(hi,xDir) - a_fab(lo,xDir))/(2.0*dx);
  a_yPartial = (a_fab(hi,yDir) - a_fab(lo,yDir))/(2.0*dx);
}

void computeSecondDeriv( Real            &a_x2Xi    , 
                         Real            &a_y2Xi    ,   
                         Real            &a_x2Eta     ,  
                         Real            &a_y2Eta     ,
                         Real            &a_xXiEta  ,
                         Real            &a_yXiEta  ,
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


void xiVolVarRes(Real            & a_volVarRes0 ,
                 Real            & a_volVarRes1 ,
                 Real            & a_diagTerm0  ,
                 Real            & a_diagTerm1  , 
                 const IntVect   & a_iv         ,
                 const RealVect  & a_dx         ,
                 const FArrayBox & a_curr       ,
                 const FArrayBox & a_mapPartials,
                 const FArrayBox & a_perpXperpY ,
                 const FArrayBox & a_perpField  ,
                 const bool      & a_writeInfo  )
{
  int xDir = 0;
  int yDir = 1;
  IntVect loEnd = a_curr.box().smallEnd();
  IntVect hiEnd = a_curr.box().bigEnd();
  
  Real dx2 = a_dx[xDir]*a_dx[xDir];
  Real dy2 = a_dx[yDir]*a_dx[yDir];

  // Component 0: First component  of unit vector pointing in the magnetic field direction
  // Component 1: Second component of unit vector pointing in the magnetic field direction
  // Component 2: derivative wrt xi[0] of the first component of target
  // Component 3: derivative wrt xi[1] of the first component of target
  // Component 4: derivative wrt xi[0] of the second component of target
  // Component 5: derivative wrt xi[1] of the second component of target

  Real B0     = a_perpField(a_iv,xDir             );
  Real B1     = a_perpField(a_iv,yDir             );
  Real B0Xi   = a_perpField(a_iv,xDir +   SpaceDim);
  Real B1Xi   = a_perpField(a_iv,yDir +   SpaceDim);
      
  // xi partials of the mapping
  Real xXi  = a_mapPartials(a_iv,xDir             ); //0
  Real yXi  = a_mapPartials(a_iv,yDir             ); //1
  Real x2Xi = a_mapPartials(a_iv,xDir + 2*SpaceDim); //4
  Real y2Xi = a_mapPartials(a_iv,yDir + 2*SpaceDim); //5
  
  // compute the derivative of B0 with repect to X and Y
  Real B0x = a_perpXperpY(a_iv,xDir           );
  Real B0y = a_perpXperpY(a_iv,yDir           );
  Real B1x = a_perpXperpY(a_iv,xDir + SpaceDim);
  Real B1y = a_perpXperpY(a_iv,yDir + SpaceDim);
  
  // q = <T_xi,B> = (xXi,yXi).dot.(B0,B1)
  Real q = xXi*B0 + yXi*B1;
  Real q2 = q*q;

  // differentiate q with respect to x and y
  Real qX = xXi*B0x + yXi*B1x; 
  Real qY = xXi*B0y + yXi*B1y;
  
  // differentiate q with respect to xXi and yXi
  //Real qXxi  = B0; 
  //Real qYxi  = B1;

  // differentiate q with respect to Xi 
  Real qXi  = x2Xi*B0 + xXi*B0Xi  + y2Xi*B1 + yXi*B1Xi; 
  
  // length squared of the (xi)tangent vector
  Real p = xXi*xXi + yXi*yXi;
  Real p2 = p*p;
  Real p4 = p2*p2;

  // differentiate p with respect to XXi and with respect to YXi
  //Real pXxi = 2.0*xXi;
  //Real pYxi = 2.0*yXi;
  
  // differentiate p with respect to Xi
  Real pXi = 2.0*(xXi*x2Xi + yXi*y2Xi);
  
  // Lagrangian: L = q*q/p;

  // differentiate L with respect to x and y (Note: p does not depend on x or y)
  Real Lx = 2.0*q*qX/p;
  Real Ly = 2.0*q*qY/p;

  // differentiate L with respect to XXi and with respect to YXi
  //LXxi  = (2.0*q*qXxi*p - pXxi*q*q)/(p*p) ;
  //LYxi  = (2.0*q*qYxi*p - pYxi*q*q)/(p*p) ;
  
  // simplify the above
  //LXxi  = 2.0*(q*B0*p - xXi*q2)/p2;
  //LYxi  = 2.0*(q*B1*p - yXi*q2)/p2 ;
  
  //Real LXEta = 0.0 ; 
  //Real LYEta = 0.0 ;

  Real LXEtaEta = 0.0 ; 
  Real LYEtaEta = 0.0 ;

  // differentiate first with respect to Xxi or YXi then differentiate again with respect to Xi
  Real LXxiXi = ((2.0*p2*((qXi*B0*p + q*B0Xi*p + q*B0*pXi) - (x2Xi*q2 + xXi*2.0*q*qXi))) - (2.0*(q*B0*p - xXi*q2)*2.0*p*pXi))/p4;
  Real LYxiXi = ((2.0*p2*((qXi*B1*p + q*B1Xi*p + q*B1*pXi) - (y2Xi*q2 + yXi*2.0*q*qXi))) - (2.0*(q*B1*p - yXi*q2)*2.0*p*pXi))/p4;
    
  //Real LXEtaEta = 0.0 ; 
  //Real LYEtaEta = 0.0 ;
  
  // Euler-Lagrange system: -(LXxiXi + LXetaEta) + Lx = 0  
  //                        -(LYxiXi + LYetaEta) + Ly = 0  
  a_volVarRes0 = -(LXxiXi + LXEtaEta) + Lx;
  a_volVarRes1 = -(LYxiXi + LYEtaEta) + Ly;

  // contributions to the diagonal term in the stencil with respect to the first  equation
  Real x2XiDiag = -2.0/(4.0*dx2);
  Real y2XiDiag = -2.0/(4.0*dy2);
  
  Real pXiDiag0 = 2.0*xXi*x2XiDiag;
  Real pXiDiag1 = 2.0*yXi*y2XiDiag;

  Real qXiDiag0 = 2.0*B0*x2XiDiag;
  Real qXiDiag1 = 2.0*B1*y2XiDiag;

  Real LXxiXiDiag = ((2.0*p2*((qXiDiag0*B0*p + q*B0*pXiDiag0) - (x2XiDiag*q2 + xXi*2.0*q*qXiDiag0))) - (2.0*(-xXi*q2)*2.0*p*pXiDiag0))/p4;
  Real LYxiXiDiag = ((2.0*p2*((qXiDiag1*B1*p + q*B1*pXiDiag1) - (y2XiDiag*q2 + yXi*2.0*q*qXiDiag1))) - (2.0*(-yXi*q2)*2.0*p*pXiDiag1))/p4;

  Real LXetaEtaDiag = 0.0;
  Real LYetaEtaDiag = 0.0;
  
  a_diagTerm0 = LXxiXiDiag + LXetaEtaDiag;
  a_diagTerm1 = LYxiXiDiag + LYetaEtaDiag;
   
  if(a_writeInfo) 
    {
      pout()<< "iv            = " <<a_iv         <<endl;
      pout()<< "a_volVarRes0  = " <<a_volVarRes0 <<endl;
      pout()<< "a_volVarRes1  = " <<a_volVarRes1 <<endl;
      pout()                                     <<endl;
    }
}

void etaVolVarRes(Real            & a_volVarRes0    ,
                  Real            & a_volVarRes1    ,
                  Real            & a_diagTerm0     ,
                  Real            & a_diagTerm1     , 
                  const IntVect   & a_iv            ,
                  const RealVect  & a_dx            ,
                  const FArrayBox & a_curr          ,
                  const FArrayBox & a_mapPartials   ,
                  const FArrayBox & a_perpXperpY    ,
                  const FArrayBox & a_perpField     ,
                  const bool      & a_writeInfo     )
{
  int xDir = 0;
  int yDir = 1;

  //int xiDir  = 0;
  int etaDir = 1;
  
  IntVect loEnd = a_curr.box().smallEnd();
  IntVect hiEnd = a_curr.box().bigEnd  ();
  
  // transition between one-sided derivatives
  int numCells = hiEnd[etaDir];
  Real delta   = a_iv[etaDir];
  delta       /= (numCells - 1);
  
  //Real dx2 = a_dx[xDir]*a_dx[xDir];
  Real dy  = a_dx[yDir];
  Real dy2 = dy*dy     ;
  
  // component-0            first  component of perp field
  // component-1            second component of perp field
  
  // component-2 partial of first  component with respect to xi
  // component-3 partial of second component with respect to xi
  
  // component-4 partial of first  component with respect to eta
  // component-5 partial of second component with respect to eta
  
  Real B0     = a_perpField(a_iv,xDir             );
  Real B1     = a_perpField(a_iv,yDir             );
  Real B0_eta  = a_perpField(a_iv,xDir + 2*SpaceDim);
  Real B1_eta  = a_perpField(a_iv,yDir + 2*SpaceDim);
      
  // eta partials of the mapping
  Real x_eta  = a_mapPartials(a_iv,xDir + 1*SpaceDim); //0
  Real y_eta  = a_mapPartials(a_iv,yDir + 1*SpaceDim); //1
  Real x_2eta = a_mapPartials(a_iv,xDir + 3*SpaceDim); //4
  Real y_2eta = a_mapPartials(a_iv,yDir + 3*SpaceDim); //5
  bool oneSidedEtaDiff = false;
  ParmParse pp;
  pp.get("oneSidedEtaDiff",oneSidedEtaDiff);

  if (oneSidedEtaDiff)
    {
      Real x_etaTemp = x_eta;
      Real y_etaTemp = y_eta;

      // the other iv is futher away from the x-point
      IntVect hiIv = a_iv + BASISV(etaDir);
      IntVect loIv = a_iv - BASISV(etaDir);
      
      Real x_etaHi =  (a_curr(hiIv, xDir) - a_curr(a_iv, xDir))/a_dx[etaDir];
      //Real y_etaHi =  (a_curr(hiIv, yDir) - a_curr(a_iv, yDir))/a_dx[etaDir];

      Real x_etaLo = -(a_curr(loIv, xDir) - a_curr(a_iv, xDir))/a_dx[etaDir];
      Real y_etaLo = -(a_curr(loIv, yDir) - a_curr(a_iv, yDir))/a_dx[etaDir];
      
      
      x_eta = delta*x_etaHi + (1.0 - delta)*x_etaLo;
      y_eta = delta*x_etaHi + (1.0 - delta)*y_etaLo;
      
      if ((Abs(x_etaTemp - x_eta) > 1.0e-03) ||(Abs(y_etaTemp - y_eta) > 1.0e-03))
      {
        pout()<<"x_eta discrepancy = "<< x_etaTemp - x_eta<<endl;
        pout()<<"y_eta discrepancy = "<< y_etaTemp - y_eta<<endl;
      }
    }
 
  // compute the derivative of B0 with repect to X and Y
  Real B0_x = a_perpXperpY(a_iv,xDir           );
  Real B0_y = a_perpXperpY(a_iv,yDir           );
  Real B1_x = a_perpXperpY(a_iv,xDir + SpaceDim);
  Real B1_y = a_perpXperpY(a_iv,yDir + SpaceDim);
  
  // q = <T_eta,B> = (x_eta,y_eta).dot.(B0,B1)
  Real q = x_eta*B0 + y_eta*B1;
  Real q2 = q*q;

  // differentiate q with respect to x and y
  Real qX = x_eta*B0_x + y_eta*B1_x; 
  Real qY = x_eta*B0_y + y_eta*B1_y;
  
  // differentiate q with respect to xXi and yXi
  //Real qX_eta  = B0; 
  //Real qY_eta  = B1;

  // differentiate q with respect to eta
  Real q_eta  = x_2eta*B0 + x_eta*B0_eta  + y_2eta*B1 + y_eta*B1_eta; 
  
  // length squared of the (eta)tangent vector
  Real p = x_eta*x_eta + y_eta*y_eta;
  Real p2 = p*p ;
  Real p3 = p*p*p;
  Real p4 = p2*p2 ;

  // differentiate p with respect to X_eta and with respect to Y_eta
  //Real pXEta = 2.0*x_eta;
  //Real pYEta = 2.0*y_eta;
  
  // differentiate p with respect to eta
  Real p_eta = 2.0*(x_eta*x_2eta + y_eta*y_2eta);
  
  // Lagrangian: L = q*q/p;

  // differentiate L with respect to x and y (Note: p does not depend on x or y)
  Real Lx = 2.0*q*qX/p;
  Real Ly = 2.0*q*qY/p;

  // differentiate L with respect to X_eta and with respect to Y_eta
  //LX_eta  = (2.0*q*qX_eta*p - pX_eta*q*q)/(p*p) ;
  //LY_eta  = (2.0*q*qY_eta*p - pY_eta*q*q)/(p*p) ;
  
  // simplify the above
  //LX_eta  = 2.0*(q*B0*p - x_eta*q2)/p2;
  //LY_eta  = 2.0*(q*B1*p - y_eta*q2)/p2 ;
  
  //Real LXxi = 0.0 ; 
  //Real LYXi = 0.0 ;

  Real LX_xi_xi = 0.0 ; 
  Real LY_xi_xi = 0.0 ;

  // differentiate first with respect to Xeta or Yeta then differentiate again with respect to eta
  Real LX_eta_eta = ((2.0*p2*((q_eta*B0*p + q*B0_eta*p + q*B0*p_eta) - (x_2eta*q2 + x_eta*2.0*q*q_eta))) - (2.0*(q*B0*p - x_eta*q2)*2.0*p*p_eta))/p4;
  Real LY_eta_eta = ((2.0*p2*((q_eta*B1*p + q*B1_eta*p + q*B1*p_eta) - (y_2eta*q2 + y_eta*2.0*q*q_eta))) - (2.0*(q*B1*p - y_eta*q2)*2.0*p*p_eta))/p4;
    
  // Euler-Lagrange system: -(LX_xi_xi + LX_eta_eta) + Lx = 0  
  //                        -(LY_xi_xi + LY_eta_eta) + Ly = 0  
  a_volVarRes0 = -(LX_xi_xi + LX_eta_eta) + Lx;
  a_volVarRes1 = -(LY_xi_xi + LY_eta_eta) + Ly;

  // contributions to the diagonal term in the stencil with respect to the first  equation
  if (oneSidedEtaDiff)
    {
      Real x_etaDiag = (1.0 - 2.0*delta)/dy;
      Real y_etaDiag = (1.0 - 2.0*delta)/dy;

      Real x_2etaDiag = -2.0/(4.0*dy2);
      Real y_2etaDiag = -2.0/(4.0*dy2);
     
      // linearize first equation around x_eta: take numerator to be three terms separated by the minus signs
      Real I0 = (q_eta*B0*p + q*B0_eta*p + q*B0*p_eta);
      
      Real q_eta_x_eta = B0_eta   ;
      Real p_x_eta     = 2.0*x_eta; 
      Real q_x_eta     = B0       ;
      Real p_eta_x_eta = 2.0*x_2eta;
      
      Real I_x_eta = q_eta_x_eta*B0*p + q_eta*B0*p_x_eta + q_x_eta*B0_eta*p + q*B0_eta*p_x_eta + q_x_eta*B0*p_eta + q*B0*p_eta_x_eta;

      Real II0 = -(x_2eta*q2 + x_eta*2.0*q*q_eta);
      
      //Real q_x_eta     = B0    ;
      Real q2_x_eta    = q*B0  ;
      //Real q_eta_x_eta = B0_eta;
      
      Real II_x_eta = -(x_2eta*q2_x_eta + 2.0*q*q_eta + x_eta*2.0*q_x_eta*q_eta + x_eta*2.0*q*q_eta_x_eta);
      

      Real III0 = -(2.0*(q*B0*p - x_eta*q2)*2.0*p*p_eta);
      //Real q_x_eta     = B0;
      //Real p_x_eta     = 2.0*x_eta; 
      //Real q2_x_eta    = q*B0;
      //Real p_eta_x_eta = 2.0*x_2eta;
      Real III_x_eta   =  -(2.0*(q_x_eta*B0*p + q*B0*p_x_eta - q2 - x_eta*q2_x_eta)*2.0*p*p_eta) - (2.0*(q*B0*p - x_eta*q2)*2.0*(p_x_eta*p_eta + p*2.0*p_eta_x_eta));

      Real LX_eta_eta_firstDerDiag = x_etaDiag*(2.0*(I_x_eta + II_x_eta)*p2 - ((I0 + II0)*2.0*p*p_x_eta)/p4 + (III_x_eta*p4 - III0*4.0*p3*p_x_eta)/(p4*p4));
     
      // linearize second equation around x_eta: take numerator to be three terms separated by the minus signs
      Real I1 = (q_eta*B1*p + q*B1_eta*p + q*B1*p_eta);
      
      Real q_eta_y_eta = B1_eta   ;
      Real p_y_eta     = 2.0*y_eta; 
      Real q_y_eta     = B1       ;
      Real p_eta_y_eta = 2.0*y_2eta;
      
      Real I_y_eta = q_eta_y_eta*B1*p + q_eta*B1*p_y_eta + q_y_eta*B1_eta*p + q*B1_eta*p_y_eta + q_y_eta*B1*p_eta + q*B1*p_eta_y_eta;

      Real II1 = -(y_2eta*q2 + y_eta*2.0*q*q_eta);
      
      //Real q_y_eta     = B1    ;
      Real q2_y_eta    = q*B1  ;
      //Real q_eta_y_eta = B1_eta;
      
      Real II_y_eta = -(y_2eta*q2_y_eta + 2.0*q*q_eta + y_eta*2.0*q_y_eta*q_eta + y_eta*2.0*q*q_eta_y_eta);
      

      Real III1 = -(2.0*(q*B1*p - y_eta*q2)*2.0*p*p_eta);
      //Real q_y_eta     = B1;
      //Real p_y_eta     = 2.0*y_eta; 
      //Real q2_y_eta    = q*B1;
      //Real p_eta_y_eta = 2.0*y_2eta;
      Real III_y_eta   =  -(2.0*(q_y_eta*B1*p + q*B1*p_y_eta - q2 - y_eta*q2_y_eta)*2.0*p*p_eta) - (2.0*(q*B1*p - y_eta*q2)*2.0*(p_y_eta*p_eta + p*2.0*p_eta_y_eta));

      Real LY_eta_eta_firstDerDiag = y_etaDiag*(2.0*(I_y_eta + II_y_eta)*p2 - ((I1 + II1)*2.0*p*p_y_eta)/p4 + (III_y_eta*p4 - III1*4.0*p3*p_y_eta)/(p4*p4));
     
      // linearize the first equation around  x_2eta and the second equation around y_2eta
      Real p_etaDiag0 = 2.0*x_eta*x_2etaDiag;
      Real p_etaDiag1 = 2.0*y_eta*y_2etaDiag;
      
      Real q_etaDiag0 = 2.0*B0*x_2etaDiag;
      Real q_etaDiag1 = 2.0*B1*y_2etaDiag;
      Real LX_eta_etaDiag = ((2.0*p2*((q_etaDiag0*B0*p + q*B1*p_etaDiag0) - (x_2etaDiag*q2 + x_eta*2.0*q*q_etaDiag0))) - (2.0*(-x_eta*q2)*2.0*p*p_etaDiag0))/p4;
      Real LY_eta_etaDiag = ((2.0*p2*((q_etaDiag1*B1*p + q*B1*p_etaDiag1) - (y_2etaDiag*q2 + y_eta*2.0*q*q_etaDiag1))) - (2.0*(-y_eta*q2)*2.0*p*p_etaDiag1))/p4;
      
      Real LX_xi_xiDiag = 0.0;
      Real LY_xi_xiDiag = 0.0;
      a_diagTerm0 = LX_xi_xiDiag + LX_eta_eta_firstDerDiag + LX_eta_etaDiag;
      a_diagTerm1 = LY_xi_xiDiag + LY_eta_eta_firstDerDiag + LY_eta_etaDiag;
      
    }
  else
    {
      Real x_2etaDiag = -2.0/(4.0*dy2);
      Real y_2etaDiag = -2.0/(4.0*dy2);
      
      Real p_etaDiag0 = 2.0*x_eta*x_2etaDiag;
      Real p_etaDiag1 = 2.0*y_eta*y_2etaDiag;
      
      Real q_etaDiag0 = 2.0*B0*x_2etaDiag;
      Real q_etaDiag1 = 2.0*B1*y_2etaDiag;
      
      Real LX_eta_etaDiag = ((2.0*p2*((q_etaDiag0*B0*p + q*B0*p_etaDiag0) - (x_2etaDiag*q2 + x_eta*2.0*q*q_etaDiag0))) - (2.0*(-x_eta*q2)*2.0*p*p_etaDiag0))/p4;
      Real LY_eta_etaDiag = ((2.0*p2*((q_etaDiag1*B1*p + q*B1*p_etaDiag1) - (y_2etaDiag*q2 + y_eta*2.0*q*q_etaDiag1))) - (2.0*(-y_eta*q2)*2.0*p*p_etaDiag1))/p4;
      
      Real LX_xi_xiDiag = 0.0;
      Real LY_xi_xiDiag = 0.0;
      a_diagTerm0 = LX_xi_xiDiag + LX_eta_etaDiag;
      a_diagTerm1 = LY_xi_xiDiag + LY_eta_etaDiag;
    }
  
   
  if(a_writeInfo) 
    {
      pout()<< "iv            = " <<a_iv         <<endl;
      pout()<< "a_volVarRes0  = " <<a_volVarRes0 <<endl;
      pout()<< "a_volVarRes1  = " <<a_volVarRes1 <<endl;
      pout()                                     <<endl;
    }
}

void volVarResOld(Real            & a_volVarRes0   ,
                  Real            & a_volVarRes1    ,
                  Real            & a_diagTerm0     ,
                  Real            & a_diagTerm1     ,   
                  const IntVect   & a_iv            ,
                  const RealVect  & a_dx            ,
                  const FArrayBox & a_curr          ,
                  const FArrayBox & a_mapPartials   ,
                  const FArrayBox & a_perpXperpY    ,
                  const FArrayBox & a_perpField     ,
                  const bool      & a_writeInfo     )
{
  ParmParse pp;
  bool xiTangent  = false;
  bool etaTangent = false;
  pp.get("xiTangent" ,xiTangent);
  pp.get("etaTangent",etaTangent);

  CH_assert((xiTangent == false && etaTangent == true) ||(xiTangent == true && etaTangent == false ));

  if (xiTangent)
    {
      xiVolVarRes(a_volVarRes0 ,
                  a_volVarRes1 ,
                  a_diagTerm0  ,
                  a_diagTerm1  , 
                  a_iv         ,
                  a_dx         ,
                  a_curr       ,
                  a_mapPartials,
                  a_perpXperpY ,
                  a_perpField  ,
                  a_writeInfo  );
    }
  else
    {
      etaVolVarRes(a_volVarRes0    ,
                   a_volVarRes1    ,
                   a_diagTerm0     ,
                   a_diagTerm1     ,   
                   a_iv            ,
                   a_dx            ,
                   a_curr          ,
                   a_mapPartials   ,
                   a_perpXperpY    ,
                   a_perpField     ,
                   a_writeInfo     );
    }

}

void volVarResNew(Real            & a_volVarResXi0  ,
                  Real            & a_volVarResXi1  ,
                  Real            & a_volVarResEta0 ,
                  Real            & a_volVarResEta1 ,
                  Real            & a_diagTermXi0   ,
                  Real            & a_diagTermXi1   ,   
                  Real            & a_diagTermEta0  ,
                  Real            & a_diagTermEta1  ,   
                  const IntVect   & a_iv            ,
                  const RealVect  & a_dx            ,
                  const FArrayBox & a_curr          ,
                  const FArrayBox & a_mapPartials   ,
                  const FArrayBox & a_perpXperpY    ,
                  const FArrayBox & a_perpField     ,
                  const FArrayBox & a_BxBy          ,
                  const FArrayBox & a_BField        ,
                  const bool      & a_writeInfo     )
{
  
  xiVolVarRes(a_volVarResXi0,
              a_volVarResXi1,
              a_diagTermXi0 ,
              a_diagTermXi1 , 
              a_iv          ,
              a_dx          ,
              a_curr        ,
              a_mapPartials ,
              a_BxBy        ,
              a_BField      ,
              a_writeInfo   );

  etaVolVarRes(a_volVarResEta0 ,
               a_volVarResEta1 ,
               a_diagTermEta0  ,
               a_diagTermEta1  ,   
               a_iv            ,
               a_dx            ,
               a_curr          ,
               a_mapPartials   ,
               a_perpXperpY    ,
               a_perpField     ,
               a_writeInfo     );
}

Real setSmoothVolVarWeight(const IntVect   & a_iv   , 
                           const Real      & a_sigma,
                           const Real      & a_expo ,
                           const RealVect  & a_pt   ,
                           const RealVect  & a_xPt  )
{
  Real dist = 0.0;
  for (int idir = 0;idir <SpaceDim;++idir)
    {
      dist += (a_pt[idir] - a_xPt[idir])*(a_pt[idir] - a_xPt[idir]);
    }
  
  dist = sqrt(dist);
  
  Real weight = LARGEREALVAL;
  if (a_sigma <= 0.0)
    {
      if (a_sigma < 0.0)
        {
          MayDay::Abort("sigma must be non-negative");
        }
      else
        {
          weight = 0.999;
        }
    }
  else
    {
      weight = exp(-pow(dist/a_sigma,a_expo));
    }
  // pout()<<"xPt    = "<<a_xPt<<endl;
  // pout()<<"pt     = "<<a_pt<<endl;
  // pout()<<"weight = "<<weight<<endl;
  
  return weight;
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


void setWeights(FArrayBox      & a_weights    ,
                const IntVect  & a_nCell      ,
                const Real     & a_xiEtaWeight,
                const Real     & a_sigma      ,
                const Real     & a_expo       ,
                const RealVect & a_dx         ,
                const bool     & a_hasXPt     ,
                const IntVect  & a_xiEtaXPt   )
 {
   int xiEtaComp        = 0;
   int smoothVolVarComp = 1;
   
   a_weights.setVal(a_xiEtaWeight,xiEtaComp       );
   a_weights.setVal(LARGEREALVAL ,smoothVolVarComp);
   
   // iteration box
   Box box = a_weights.box();
    
   int xComp = 0;
   int yComp = 1;
 
   RealVect xPt = LARGEREALVAL*RealVect::Unit;

   if (a_hasXPt)
     {
       for(int idir = 0;idir < SpaceDim;++idir)
         { 
           xPt[idir]   = a_xiEtaXPt[idir]*a_nCell[idir]*a_dx[idir];
         }

       // iterate over the interior cells
       for (BoxIterator bit(box); bit.ok(); ++bit)
         {
           const IntVect iv = bit();
          
           
           RealVect pt;
           pt[0] = iv[xComp]*a_dx[xComp];
           pt[1] = iv[yComp]*a_dx[yComp];
           
           a_weights(iv,smoothVolVarComp) = setSmoothVolVarWeight(iv     ,
                                                                  a_sigma,  
                                                                  a_expo ,
                                                                  pt     ,
                                                                  xPt    );
         }
     }
   else
     {
       pout()<<"No xpoint!!"<<endl;
       a_weights.setVal(0.9999999999,smoothVolVarComp);
     }
   
 }

void computeIncOld(FArrayBox       & a_inc           ,
                   FArrayBox       & a_smoothRes     ,
                   FArrayBox       & a_volVarRes     ,
                   FArrayBox       & a_totalRes      ,
                   const FArrayBox & a_curr          ,
                   const FArrayBox & a_mapPartials   ,
                   const FArrayBox & a_smoothness    ,
                   const FArrayBox & a_smoothDiag    ,
                   const FArrayBox & a_Jacobian      ,
                   const FArrayBox & a_weights       ,
                   const FArrayBox & a_perpXperpY    ,
                   const FArrayBox & a_perpField     ,
                   const bool      & a_writeInfo     ,
                   const RealVect  & a_dx            ,
                   const int       & a_xPointHiLoSide,
                   const int       & a_iter          ,
                   const int       & a_numInitSmooth )
{
  // iteration boxes
  Box box    = a_inc.box();
  Box intBox = box.grow(-1);

  // iterate over the interior to get data for boundary op
  for (BoxIterator bit(intBox); bit.ok(); ++bit)
    {
      // IntVect
      const IntVect& iv = bit();   
      
      // indexing macros
      int xDir = 0;
      int yDir = 1;

      int xiEtaComp        = 0;
      int smoothVolVarComp = 1;

      Real xiEtaWeight        = a_weights(iv,xiEtaComp       );
      Real smoothVolVarWeight = a_weights(iv,smoothVolVarComp);
      if (a_iter < a_numInitSmooth)
        {
          smoothVolVarWeight = 1.0;
        }

      Real xiSWeight  = xiEtaWeight;
      Real etaSWeight = (1.0 - xiEtaWeight);
      
      // smooth res
      Real smoothRes0 = LARGEREALVAL;
      Real smoothRes1 = LARGEREALVAL;
      
      // smooth diag term
      Real smoothDiag0 = a_smoothDiag(iv,xDir);
      Real smoothDiag1 = a_smoothDiag(iv,yDir);
            
      SmJnRes(smoothRes0   ,
              smoothRes1   , 
              xiSWeight    ,
              etaSWeight   ,
              iv           ,
              a_mapPartials,
              a_Jacobian   ,
              a_writeInfo );

      a_smoothRes(iv,xDir) = smoothRes0;
      a_smoothRes(iv,yDir) = smoothRes1;
     
      Real volVarRes0      = LARGEREALVAL;
      Real volVarRes1      = LARGEREALVAL;
      Real volVarDiagTerm0 = LARGEREALVAL;
      Real volVarDiagTerm1 = LARGEREALVAL;

      volVarResOld(volVarRes0      , 
                   volVarRes1      ,
                   volVarDiagTerm0 ,
                   volVarDiagTerm1 ,
                   iv              ,
                   a_dx            ,
                   a_curr          ,
                   a_mapPartials   ,
                   a_perpXperpY    ,
                   a_perpField     ,
                   a_writeInfo     );

      // todo:: where should this go or come from
      volVarRes0 *= -1.0;
      volVarRes1 *= -1.0;

      a_volVarRes(iv,xDir) = volVarRes0;
      a_volVarRes(iv,yDir) = volVarRes1;

      Real totalRes0 = LARGEREALVAL;
      Real totalRes1 = LARGEREALVAL;
      
      // residual and diagonal term for all the principles combined
      totalRes0 = smoothRes0*smoothVolVarWeight + volVarRes0*(1.0 - smoothVolVarWeight);
      totalRes1 = smoothRes1*smoothVolVarWeight + volVarRes1*(1.0 - smoothVolVarWeight);
          
      Real diagTerm0 = smoothDiag0*smoothVolVarWeight + volVarDiagTerm0*(1.0 - smoothVolVarWeight);
      Real diagTerm1 = smoothDiag1*smoothVolVarWeight + volVarDiagTerm1*(1.0 - smoothVolVarWeight);
          
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
      
      a_inc(iv,xDir) = lambda0*totalRes0;
      a_inc(iv,yDir) = lambda1*totalRes1;

      //a_inc(iv,xDir) = (-0.5*smoothRes0/smoothDiag0) + (-0.000001*volVarRes0/volVarDiagTerm0);
      //a_inc(iv,yDir) = (-0.5*smoothRes1/smoothDiag1) + (-0.000001*volVarRes1/volVarDiagTerm1);
     
      a_totalRes(iv,xDir) = totalRes0;
      a_totalRes(iv,yDir) = totalRes1;
      
      a_smoothRes(iv,xDir) = smoothRes0;
      a_smoothRes(iv,yDir) = smoothRes1;
      
      a_volVarRes(iv,xDir) = volVarRes0;
      a_volVarRes(iv,yDir) = volVarRes1;
    } 
}

void computeIncNew(FArrayBox       & a_inc           ,
                   FArrayBox       & a_smoothRes     ,
                   FArrayBox       & a_volVarRes     ,
                   FArrayBox       & a_totalRes      ,
                   const FArrayBox & a_curr          ,
                   const FArrayBox & a_mapPartials   ,
                   const FArrayBox & a_smoothness    ,
                   const FArrayBox & a_smoothDiag    ,
                   const FArrayBox & a_Jacobian      ,
                   const FArrayBox & a_weights       ,
                   const FArrayBox & a_perpXperpY    ,
                   const FArrayBox & a_perpField     ,
                   const FArrayBox & a_BxBy          ,
                   const FArrayBox & a_BField        ,
                   const bool      & a_writeInfo     ,
                   const RealVect  & a_dx            ,
                   const int       & a_iter          ,
                   const int       & a_numInitSmooth )
{
  // iteration boxes
  Box box    = a_inc.box();
  Box intBox = box.grow(-1);

  // iterate over the interior to get data for boundary op
  for (BoxIterator bit(intBox); bit.ok(); ++bit)
    {
      // IntVect
      const IntVect& iv = bit();   
      
      // indexing macros
      int xDir = 0;
      int yDir = 1;

      int xiEtaComp        = 0;
      int smoothVolVarComp = 1;

      Real xiEtaWeight        = a_weights(iv,xiEtaComp       );
      Real smoothVolVarWeight = a_weights(iv,smoothVolVarComp);
      if (a_iter < a_numInitSmooth)
        {
          smoothVolVarWeight = 1.0;
        }

      Real xiSWeight  = xiEtaWeight;
      Real etaSWeight = (1.0 - xiEtaWeight);
      
      // smooth res
      Real smoothRes0 = LARGEREALVAL;
      Real smoothRes1 = LARGEREALVAL;
      
      // smooth diag term
      Real smoothDiag0 = a_smoothDiag(iv,xDir);
      Real smoothDiag1 = a_smoothDiag(iv,yDir);
            
      SmJnRes(smoothRes0   ,
              smoothRes1   , 
              xiSWeight    ,
              etaSWeight   ,
              iv           ,
              a_mapPartials,
              a_Jacobian   ,
              a_writeInfo );

      a_smoothRes(iv,xDir) = smoothRes0;
      a_smoothRes(iv,yDir) = smoothRes1;
     
      Real volVarResXi0      = LARGEREALVAL;
      Real volVarResXi1      = LARGEREALVAL;

      Real volVarResEta0      = LARGEREALVAL;
      Real volVarResEta1      = LARGEREALVAL;
      
      Real volVarDiagTermXi0 = LARGEREALVAL;
      Real volVarDiagTermXi1 = LARGEREALVAL;
      
      Real volVarDiagTermEta0 = LARGEREALVAL;
      Real volVarDiagTermEta1 = LARGEREALVAL;

      volVarResNew(volVarResXi0      ,   
                   volVarResXi1      ,
                   volVarResEta0     , 
                   volVarResEta1     ,
                   volVarDiagTermXi0 ,
                   volVarDiagTermXi1 ,
                   volVarDiagTermEta0,
                   volVarDiagTermEta1,
                   iv                ,
                   a_dx              ,
                   a_curr            ,
                   a_mapPartials     ,
                   a_perpXperpY      ,
                   a_perpField       ,
                   a_BxBy            ,
                   a_BField          ,
                   a_writeInfo       );

      // todo:: where should this go or come from
      volVarResXi0 *= -1.0;
      volVarResXi1 *= -1.0;

      volVarResEta0 *= -1.0;
      volVarResEta1 *= -1.0;
      
      int xiTan  = 0;
      int etaTan = 1;
      
      a_volVarRes(iv,xDir + xiTan) = volVarResXi0;
      a_volVarRes(iv,yDir + xiTan) = volVarResXi1;
      
      a_volVarRes(iv,xDir + etaTan*SpaceDim) = volVarResEta0;
      a_volVarRes(iv,yDir + etaTan*SpaceDim) = volVarResEta1;
      

      Real totalRes0 = LARGEREALVAL;
      Real totalRes1 = LARGEREALVAL;
      
      // residual and diagonal term for all the principles combined
      totalRes0 = smoothRes0*smoothVolVarWeight + (volVarResXi0 + volVarResEta0)*(1.0 - smoothVolVarWeight);
      totalRes1 = smoothRes1*smoothVolVarWeight + (volVarResXi1 + volVarResEta1)*(1.0 - smoothVolVarWeight);
          
      Real diagTerm0 = smoothDiag0*smoothVolVarWeight + (volVarDiagTermXi0 + volVarDiagTermEta0)*(1.0 - smoothVolVarWeight);
      Real diagTerm1 = smoothDiag1*smoothVolVarWeight + (volVarDiagTermXi1 + volVarDiagTermEta1)*(1.0 - smoothVolVarWeight);
          
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
      
      a_inc(iv,xDir) = lambda0*totalRes0;
      a_inc(iv,yDir) = lambda1*totalRes1;

      a_totalRes(iv,xDir) = totalRes0;
      a_totalRes(iv,yDir) = totalRes1;
      
      a_smoothRes(iv,xDir) = smoothRes0;
      a_smoothRes(iv,yDir) = smoothRes1;
      
      a_volVarRes(iv,xDir) = volVarResXi0;
      a_volVarRes(iv,yDir) = volVarResXi1;
    } 
}

bool stoppingCriterion(const Real& a_maxNorm                  ,
                       const Real& a_maxTol                   ,
                       const int & a_iter                     ,
                       const int & a_maxIter_BrackbillSalzman)
{
  bool keepWorking = true;
  if (a_maxNorm < a_maxTol || a_iter >= a_maxIter_BrackbillSalzman)
    {
      keepWorking = false;
    }

  return keepWorking;
}

void computeMapPartials(FArrayBox       &a_mapPartials,
                        const FArrayBox &a_curr       ,
                        const RealVect  &a_dx         ,
                        const bool      &a_writeInfo  )
{
  // macros for components of the mapping
  int xDir = 0;
  int yDir = 1;

  Box domainBox = a_curr.box();
 
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
          firstDerivCenteredDiff(xDeriv     ,
                                 yDeriv     ,
                                 iv         ,
                                 xiEtaDir   ,
                                 a_dx       ,
                                 a_curr     ,
                                 a_writeInfo);
                      
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
      computeSecondDeriv(x2Xi       ,
                         y2Xi       ,
                         x2Eta      ,
                         y2Eta      ,
                         xXiEta     ,
                         yXiEta     ,
                         iv         ,
                         a_dx       ,
                         a_curr     ,
                         a_writeInfo);
      
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
  
  // set the corners, iterate over four edge boxes using sometimes centered and sometimes one-sided stencils (which are based on quadratic extrapolation)
  bool BFieldCase = false;
  computeMapPartialsBdIv(a_mapPartials, 
                         a_curr       ,
                         a_dx         ,
                         BFieldCase   );

}

void xiEtaTranspose(FArrayBox & a_mappedFab)
{
  int xDir = 0;
  int yDir = 1;
  
  Box box = a_mappedFab.box();
  int nComp = a_mappedFab.nComp();
  FArrayBox transposedMap(box,nComp);

    for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      // IntVect
      const IntVect& iv = bit();
      IntVect vi;
      vi[0] = iv[1];
      vi[1] = iv[0];
      transposedMap(iv,xDir) = a_mappedFab(vi,xDir);
      transposedMap(iv,yDir) = a_mappedFab(vi,yDir);
    }
    a_mappedFab.copy(transposedMap);
}

Real checkAnnulus (const RealVect & a_res       ,
                   const Real     & a_xXi       , 
                   const Real     & a_yXi       ,
                   const Real     & a_xEta      , 
                   const Real     & a_yEta      ,
                   const IntVect  & a_iv        ,
                   const RealVect & a_dx        ,
                   const int      & a_iEdge     ,
                   const Real     & a_xiSWeight ,
                   const Real     & a_etaSWeight,
                   const Real     & a_J         ,
                   const Real     & a_JxXi      ,
                   const Real     & a_JyXi      ,
                   const Real     & a_JxEta     ,
                   const Real     & a_JyEta     ,
                   const RealVect & a_tanVec    )
{
  Real pi = 4.0*atan(1.0);
  
  Real xi  = a_iv[0]*a_dx[0];
  Real eta = a_iv[1]*a_dx[1]; 
 
  Real r   = 2.0*xi + 1;
  Real t   = pi*eta + pi/2.0;

  //Real x = r*cos(t);
  //Real y = r*sin(t);

  Real x_xi  = 2.0*cos(t);
  Real x_eta = -r*pi*sin(t);

  Real y_xi  = 2.0*sin(t);
  Real y_eta = r*pi*cos(t);

  Real Je = x_xi*y_eta - x_eta*y_xi;
  Real Je2 = Je*Je;

  Real Se = x_xi *x_xi  + y_xi *y_xi  + x_eta *x_eta  + y_eta *y_eta ;
  //Real Sc = a_xXi*a_xXi + a_yXi*a_yXi + a_xEta*a_xEta + a_yEta*a_yEta;

  RealVect res;
  
  if (a_iEdge == 1 || a_iEdge == 3)
    {
      res[0] = a_xiSWeight *(2.0*x_xi*Je - y_eta*Se)/Je2;
      res[1] = a_etaSWeight*(2.0*y_xi*Je + x_eta*Se)/Je2;
    }
  else
    {
      res[0] = a_xiSWeight *(2.0*x_eta*Je + y_xi*Se)/Je2;
      res[1] = a_etaSWeight*(2.0*y_eta*Je - x_xi*Se)/Je2;
    }

  pout()<<endl;
  pout()<<"exact res    = "<<   res <<endl;
  pout()<<"computed res = "<< a_res <<endl;
  pout()<<endl;

  Real retval = sqrt(a_tanVec[0]*(res[0] - a_res[0])*a_tanVec[0]*(res[0] - a_res[0]) + a_tanVec[1]*(res[1] - a_res[1])*a_tanVec[1]*(res[1] - a_res[1]));
  return retval;
}


void SmJnBoundaryOperator(Vector<FArrayBox*>       & a_arcLengthInc  ,
                          FArrayBox                & a_smoothRes     ,
                          const FArrayBox          & a_Jacobian      ,
                          const FArrayBox          & a_mapPartials   ,
                          const Vector<int>        & a_tanDir        ,
                          const Vector<FArrayBox*> & a_tangent       ,
                          const FArrayBox          & a_smoothDiag    ,
                          const RealVect           & a_dx            ,
                          const FArrayBox          & a_weights       ,
                          const int                & a_iter          )
{
  int xiEtaComp        = 0;
  //int smoothVolVarComp = 1;

  // J = Xxi*Yeta - Xeta*Yxi
  
  //S = (xiSWeight*(Xxi*Xxi + Yxi*Yxi) +etaSWeight*(Xeta*Xeta + Yeta*Yeta))
  //smooth boundary operator: hilo*(SXxi,SYxi) on left and right edges or hilo*(SXeta,SYeta) on top and bottom edges  

  // indexing macros
  //int xiDir  = 0;
  //int etaDir = 1;
  
  int xDir = 0;
  int yDir = 1;

  int zerothComp = 0;

  // debugging output
  Vector<Real> incTotal(4);

  // iterate over edges: low, right, top, left.
  for (int iEdge = 0; iEdge < 4;++iEdge)
    {
      // debugging output
      incTotal[iEdge] = 0.0;
      
      const FArrayBox& tangent   = *(a_tangent     [iEdge]);
      FArrayBox& arcLengthInc    = *(a_arcLengthInc[iEdge]);
      int tanDir                 = a_tanDir        [iEdge] ;

      // set the increments to zero at the corners
      arcLengthInc(tangent.box().smallEnd(),zerothComp) = 0.0;
      arcLengthInc(tangent.box().bigEnd  (),zerothComp) = 0.0;
      
      // iv interior to this edge
      Box intBox = tangent.box();
      intBox.grow(tanDir,-1);
      
      //debugging output for annulus 
      Real totalDisc  = 0.0;
      Real maxDisc    = 0.0;
      IntVect maxIv   = LARGEREALVAL*IntVect::Unit;
      bool printTotal = false;

      // iterate
      for (BoxIterator bit(intBox); bit.ok(); ++bit)
        { 
          // IntVect
          const IntVect iv = bit();
          
          Real xiEtaWeight        = a_weights(iv,xiEtaComp       );
          //Real smoothVolVarWeight = a_weights(iv,smoothVolVarComp);
          
          Real xiSWeight  = xiEtaWeight      ;
          Real etaSWeight = 1.0 - xiEtaWeight;
          
          //tangent vector
          RealVect tanVec;
          tanVec[0] = tangent(iv,xDir); 
          tanVec[1] = tangent(iv,yDir);
          
          // diagonal term
          Real diag0 = a_smoothDiag(iv,xDir);
          Real diag1 = a_smoothDiag(iv,yDir);
          
          // lambda
          Real lambda0 = -0.5/diag0;
          Real lambda1 = -0.5/diag1;

          Real signFix0 = 1.0;
          if (lambda0 < 0.0)
            {
              signFix0 = -1.0;
            }
          
          Real signFix1 = 1.0;
          if (lambda1 < 0.0)
            {
              signFix1 = -1.0;
            }
          Real minLambda = Min(Abs(lambda0),Abs(lambda1));
          lambda0 = signFix0*minLambda;
          lambda1 = signFix1*minLambda;
          
           // indexing macros
          //int xiSmooth  = 0;
          //int etaSmooth = 1;
          
          int xDir = 0;
          int yDir = 1;
          
          int xiDir = 0;
          int etaDir = 1;
          
          Real m = LARGEREALVAL;
          Real n = LARGEREALVAL;
           
          ParmParse pp;
          pp.get("Spower",m);
          pp.get("Jpower",n);
  
          // weights
          //Real alpha = xiSWeight;
          //Real beta  = etaSWeight;
          
          Real x_xi    = a_mapPartials(iv,xDir + 0*SpaceDim);
          Real y_xi    = a_mapPartials(iv,yDir + 0*SpaceDim);
          
          Real x_eta   = a_mapPartials(iv,xDir + 1*SpaceDim);
          Real y_eta   = a_mapPartials(iv,yDir + 1*SpaceDim);
          
          //Real x_2xi   = a_mapPartials(iv,xDir + 2*SpaceDim);
          //Real y_2xi   = a_mapPartials(iv,yDir + 2*SpaceDim);
          
          //Real x_2eta  = a_mapPartials(iv,xDir + 3*SpaceDim);
          //Real y_2eta  = a_mapPartials(iv,yDir + 3*SpaceDim);
          
          //Real x_xi_eta = a_mapPartials(iv,xDir + 4*SpaceDim);
          //Real y_xi_eta = a_mapPartials(iv,yDir + 4*SpaceDim);
          
          // mixed partials are equal
          //Real x_eta_xi = x_xi_eta;
          //Real y_eta_xi = y_xi_eta;
          
          // sum of squares
          Real S  = x_xi*x_xi + y_xi*y_xi + x_eta*x_eta + y_eta*y_eta;
          
          Real Sx_xi  = 2.0*x_xi ;
          Real Sy_xi  = 2.0*y_xi ;
          Real Sx_eta = 2.0*x_eta;
          Real Sy_eta = 2.0*y_eta;
          
          //Real Sx_xi_xi   = 2.0*x_2xi ;
          //Real Sy_xi_xi   = 2.0*y_2xi ;
          //Real Sx_eta_eta = 2.0*x_2eta;
          //Real Sy_eta_eta = 2.0*y_2eta;
          
          //Real S_xi  = Sx_xi*x_2xi    + Sy_xi*y_2xi    + Sx_eta*x_eta_xi + Sy_eta*y_eta_xi;
          //Real S_eta = Sx_xi*x_xi_eta + Sy_xi*y_xi_eta + Sx_eta*x_2eta   + Sy_eta*y_2eta  ;
          
          // Jacobian
          Real J  = a_Jacobian(iv,0);                           //0
          
          // derivatives with respect to derivatives of the mapping
          Real Jx_xi  = a_Jacobian(iv,1 + 0*SpaceDim + xDir    );//1
          Real Jy_xi  = a_Jacobian(iv,1 + 0*SpaceDim + yDir    );//2
          Real Jx_eta = a_Jacobian(iv,1 + 1*SpaceDim + xDir    );//3
          Real Jy_eta = a_Jacobian(iv,1 + 1*SpaceDim + yDir    );//4
          
          // derivatives of the above with respect to xi and eta
          //Real Jx_xi_xi   = a_Jacobian(iv,1 + 2*SpaceDim + xDir);//5
          //Real Jy_xi_xi   = a_Jacobian(iv,1 + 2*SpaceDim + yDir);//6
          //Real Jx_eta_eta = a_Jacobian(iv,1 + 3*SpaceDim + xDir);//7
          //Real Jy_eta_eta = a_Jacobian(iv,1 + 3*SpaceDim + yDir);//8
          
          // derivatives with respect to xi and eta
          //Real J_xi  = a_Jacobian(iv,1 + 4*SpaceDim + xiDir    );//9
          //Real J_eta = a_Jacobian(iv,1 + 4*SpaceDim + etaDir   );//10
          
          Real L1 = pow(S,m);
          Real L2 = pow(J,n);
          
          //Lagrangian; L = L1*L2;
          //Real L1_xi  = m*pow(S,m-1)*S_xi ;
          //Real L1_eta = m*pow(S,m-1)*S_eta;
          
          //Real L2_xi  = n*pow(J,n-1)*J_xi ;
          //Real L2_eta = n*pow(J,n-1)*J_eta;
          
          Real L1x_xi  = m*pow(S,m-1)*Sx_xi;
          Real L1y_xi  = m*pow(S,m-1)*Sy_xi;
          
          Real L1x_eta = m*pow(S,m-1)*Sx_eta;
          Real L1y_eta = m*pow(S,m-1)*Sy_eta;
  
          //Real L1x_xi_xi   = m*(m-1)*pow(S,m-2)*S_xi *Sx_xi  + m*pow(S,m-1)*Sx_xi *Sx_xi_xi  ;
          //Real L1y_xi_xi   = m*(m-1)*pow(S,m-2)*S_xi *Sy_xi  + m*pow(S,m-1)*Sy_xi *Sy_xi_xi  ;
          //Real L1x_eta_eta = m*(m-1)*pow(S,m-2)*S_eta*Sx_eta + m*pow(S,m-1)*Sx_eta*Sx_eta_eta;
          //Real L1y_eta_eta = m*(m-1)*pow(S,m-2)*S_eta*Sy_eta + m*pow(S,m-1)*Sy_eta*Sy_eta_eta;
          
          Real L2x_xi  = n*pow(J,n-1)*Jx_xi;
          Real L2y_xi  = n*pow(J,n-1)*Jy_xi;
          
          Real L2x_eta = n*pow(J,n-1)*Jx_eta;
          Real L2y_eta = n*pow(J,n-1)*Jy_eta;
          
          //Real L2x_xi_xi   = n*(n-1)*pow(J,n-2)*J_xi *Jx_xi  + n*pow(J,n-1)*Jx_xi *Jx_xi_xi ;
          //Real L2y_xi_xi   = n*(n-1)*pow(J,n-2)*J_xi *Jy_xi  + n*pow(J,n-1)*Jy_xi *Jy_xi_xi ;
          //Real L2x_eta_eta = n*(n-1)*pow(J,n-2)*J_eta*Jx_eta + n*pow(J,n-1)*Jx_eta*Jx_eta_eta;
          //Real L2y_eta_eta = n*(n-1)*pow(J,n-2)*J_eta*Jy_eta + n*pow(J,n-1)*Jy_eta*Jy_eta_eta;
          
          Real Lx_xi  = L1x_xi* L2 + L1*L2x_xi ;
          Real Ly_xi  = L1y_xi* L2 + L1*L2y_xi ;

          Real Lx_eta = L1x_eta*L2 + L1*L2x_eta;
          Real Ly_eta = L1y_eta*L2 + L1*L2y_eta;

          // residual
          RealVect res = RealVect::Unit*LARGEREALVAL;
                  
          // left or right edges
          if ((iEdge == 1) || (iEdge == 3))
            {
              res[xiDir ] = xiSWeight*Lx_xi;
              res[etaDir] = xiSWeight*Ly_xi;
            }
           
          // top or bottom edges
           else
            {
              res[xiDir ] = etaSWeight*Lx_eta;
              res[etaDir] = etaSWeight*Ly_eta;
            }
      
          // project residual onto the boundary
          Real inc  = lambda0*res[xiDir ]*tanVec[xiDir ];
          inc      += lambda1*res[etaDir]*tanVec[etaDir];
       
          // the new spacing between iv and iv - 1 is old spacing + inc
          arcLengthInc(iv,0)= inc;

          // for diagnostics
          a_smoothRes(iv,xiDir ) = res[xiDir ]*tanVec[xiDir ];
          a_smoothRes(iv,etaDir) = res[etaDir]*tanVec[etaDir];
          
          pout()<<"iv = "<<iv << ", inc = "<<inc<<endl;
        }
      pout()<<endl;
      
      if (printTotal)
        {
          pout()<<"Iteration num = "<<a_iter           <<endl;
          pout()<<"Data for edge = "<<iEdge            <<endl;
          pout()<<"totalDisc     = "<<totalDisc*a_dx[0]<<endl;
          pout()<<"maxDisc       = "<<maxDisc          <<endl;
          pout()<<"maxIv         = "<<maxIv            <<endl;
          pout()<<"incTotal      = "<<incTotal[iEdge]  <<endl;
          pout()<<"                "                   <<endl; 
        }
    }
}

Real evalSmoothLagrangian(FArrayBox       & a_smoothLagrangian,
                          const FArrayBox & a_mapPartials     ,
                          const FArrayBox & a_Jacobian        ,
                          const RealVect  & a_dx              )
{
  int xDir = 0;
  int yDir = 1;
  
  int xiDir  = 0;
  int etaDir = 1;
  
  Box box = a_smoothLagrangian.box();
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      // first partials of the map
      Real Xxi  = a_mapPartials(iv,xDir + xiDir *SpaceDim);
      Real Yxi  = a_mapPartials(iv,yDir + xiDir *SpaceDim); 
      Real Xeta = a_mapPartials(iv,xDir + etaDir*SpaceDim);
      Real Yeta = a_mapPartials(iv,yDir + etaDir*SpaceDim);
      
      //Jacobian in the zeroth component
      Real J = a_Jacobian(iv,0);

      a_smoothLagrangian(iv,0) = (Xxi*Xxi + Xeta*Xeta + Yxi*Yxi + Yeta*Yeta)/J;
    }
  
  return a_smoothLagrangian.norm(1)*a_dx[xiDir]*a_dx[etaDir];
}

Real evalOneOverJLagrangian(FArrayBox       & a_smoothLagrangian,
                            const FArrayBox & a_mapPartials     ,
                            const FArrayBox & a_Jacobian        ,
                            const RealVect  & a_dx              )
{
  //int xDir = 0;
  //int yDir = 1;
  
  int xiDir  = 0;
  int etaDir = 1;
  
  Box box = a_smoothLagrangian.box();
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      
      //Jacobian in the zeroth component
      Real J = a_Jacobian(iv,0);

      a_smoothLagrangian(iv,0) = 1.0/J;
    }
  
  return a_smoothLagrangian.norm(1)*a_dx[xiDir]*a_dx[etaDir];
}

Real evalOneOverSLagrangian(FArrayBox       & a_smoothLagrangian,
                            const FArrayBox & a_mapPartials     ,
                            const FArrayBox & a_Jacobian        ,
                            const RealVect  & a_dx              )
{
  int xDir = 0;
  int yDir = 1;
  
  int xiDir  = 0;
  int etaDir = 1;
  
  Box box = a_smoothLagrangian.box();
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      
      // first partials of the map
      Real x_xi  = a_mapPartials(iv,xDir + xiDir *SpaceDim);
      Real x_eta = a_mapPartials(iv,xDir + etaDir*SpaceDim);
      Real y_xi  = a_mapPartials(iv,yDir + xiDir *SpaceDim); 
      Real y_eta = a_mapPartials(iv,yDir + etaDir*SpaceDim);    
      
      //sum of squares of first derivatives
      Real S = x_xi*x_xi + y_xi*y_xi + x_eta*x_eta + y_eta*y_eta;

      a_smoothLagrangian(iv,0) = 1/S;
    }
  
  return a_smoothLagrangian.norm(1)*a_dx[xiDir]*a_dx[etaDir];
}

void  xPointHiLo(int               & a_hiLo     ,
                 const std::string & a_blockName)
{
  int retval = LARGEINTVAL;
  int blockNumber = getBlockNumber(a_blockName);
  if (blockNumber == 0)
    {
      retval = 1.0;
    }
  else if (blockNumber == 1)
    {
      retval =  -1.0;
    }
   else if (blockNumber == 2)
    {
      retval =  1.0;
    }
  else if (blockNumber == 3)
    {
      retval = -1.0;
    }
   else if (blockNumber == 4)
    {
      retval = -1.0;
    }
   else if (blockNumber == 5)
    {
      retval = 1.0;
    }
   else if (blockNumber == 6)
    {
      retval = -1.0;
    }
  else if (blockNumber == 7)
    {
      retval = 1.0;
    }
   else if (blockNumber == 8)
    {
      retval = LARGEINTVAL;
    }
   else if (blockNumber == 9)
    {
      retval = LARGEINTVAL;
    }
   else
     {
       MayDay::Abort("bad blockName");
     }
  a_hiLo = retval;
}

void pointJacobiOld(FArrayBox                & a_curr            ,
                    FArrayBox                & a_BField          ,
                    FArrayBox                & a_perpField       ,
                    FArrayBox                & a_perpXperpY      , 
                    FArrayBox                & a_mapPartials     ,  
                    FArrayBox                & a_smoothness      ,
                    FArrayBox                & a_smoothDiag      ,
                    FArrayBox                & a_Jacobian        ,
                    FArrayBox                & a_weights         ,
                    Vector<FArrayBox*>       & a_arcLength       , 
                    const Vector<FArrayBox*> & a_targetSpacing   ,
                    Real                     & a_maxNorm         ,
                    const FArrayBox          & a_hiRes           , 
                    const Vector<int>        & a_tanDir          ,
                    const Vector<FArrayBox*> & a_tangent         , 
                    const Vector<FArrayBox*> & a_hiResArcLength  , 
                    const int                & a_iter            ,
                    const bool               & a_writeInfo       ,
                    const bool               & a_writePlot       ,
                    const RealVect           & a_dx              ,
                    const int                & a_baseRes         ,
                    const std::string        & a_geometryFileName,
                    const std::string        & a_fieldFileName   ,
                    const std::string        & a_blockName       ,
                    const Real               & a_xiEtaWeight     ,   
                    const Real               & a_sigma           ,   
                    const Real               & a_expo            ,
                    const int                & a_refRatio        ,
                    const int                & a_numInitSmooth   )
{
  //fieldData object, which is required for updating the magnetic field
    
  // Construct the field data object
  FieldData fieldData(a_geometryFileName,
                      a_fieldFileName   ,
                      a_blockName       );

  int numComp = a_curr.nComp();
  Box box = a_curr.box();
 
  // allocate increments
  FArrayBox inc      (box,numComp);
  FArrayBox smoothRes(box,numComp);
  FArrayBox totalRes (box,numComp);
  FArrayBox volVarRes(box,numComp);
  
  // set the increments to zero
  inc      .setVal(0.0);
  smoothRes.setVal(0.0);
  volVarRes.setVal(0.0);
  totalRes .setVal(0.0);
   

  // track whether the xPoint is on the hi side or low side with respect to eta
  int xPointHiLoSide = LARGEINTVAL;
  xPointHiLo(xPointHiLoSide,
             a_blockName   );

  // compute change (up to a factor of lambda) between two iterations
  computeIncOld(inc            ,
                smoothRes      ,
                volVarRes      ,
                totalRes       ,
                a_curr         , 
                a_mapPartials  ,
                a_smoothness   ,
                a_smoothDiag   ,
                a_Jacobian     ,
                a_weights      ,
                a_perpXperpY   ,
                a_perpField    ,
                a_writeInfo    ,
                a_dx           ,
                xPointHiLoSide ,
                a_iter         ,
                a_numInitSmooth);

  // increment interior
  Box intBox = box;
  intBox.grow(-1);

  // comp
  int srcComp  = 0;
  int destComp = 0;

  ParmParse pp;
  bool moveInterior;
  pp.get("moveInterior",moveInterior);

  if (moveInterior)
    {
      // update curr interior
      a_curr.plus(inc     ,
                  intBox  ,
                  srcComp ,
                  destComp,
                  numComp );
    }

      // data holder for boundary increments
      Vector<FArrayBox*> arcLengthInc(4);
      for (int iEdge = 0; iEdge < 4; ++iEdge)
        {
          arcLengthInc [iEdge] = new FArrayBox(a_tangent[iEdge]->box(),1);
        }
    
  // norm of R,Z increments
  a_maxNorm = inc.norm(0,0,2);

  Real maxInc = 0.0;
  maxInc = Max(Max(Abs(inc.min(0)),Abs(inc.max(0))),Max(Abs(inc.min(1)),Abs(inc.max(1))));
          
  Real maxSmoothRes = 0.0;
  maxSmoothRes = Max(Max(Abs(smoothRes.min(0)),Abs(smoothRes.max(0))),Max(Abs(smoothRes.min(1)),Abs(smoothRes.max(1))));
          
  Real maxVolVarRes = 0.0;
  maxVolVarRes = Max(Max(Abs(volVarRes.min(0)),Abs(volVarRes.max(0))),Max(Abs(volVarRes.min(1)),Abs(volVarRes.max(1))));

  Real maxTotalRes = 0.0;
  maxTotalRes = Max(Max(Abs(totalRes.min(0)),Abs(totalRes.max(0))),Max(Abs(totalRes.min(1)),Abs(totalRes.max(1))));
          
  pout() << "max inc            : " << maxInc       << endl;
  pout() << "max smooth residual: " << maxSmoothRes << endl;
  pout() << "max volVar residual: " << maxVolVarRes << endl;
  pout() << "max total  residual: " << maxTotalRes  << endl;
        
  if(a_writePlot) 
    {  
      int printFreq = 1;
      ParmParse pp;
      pp.get("printFreq_BrackbillSalzman",printFreq);
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
              
          // volVar residual
          sprintf(name,"diagnostic/%03d.volVarRes.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&volVarRes,name,componentNames,a_dx[0]);
              
          // residual
          sprintf(name,"diagnostic/%03d.totalRes.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&totalRes ,name,componentNames,a_dx[0]);
              
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
                     a_dx         ,
                     a_writeInfo  );

  computeJacobianPartials(a_Jacobian   ,
                          a_mapPartials,
                          a_writeInfo  );
            
  computeSmJnDiagonalTerm(a_smoothDiag ,
                          a_dx         ,
                          a_weights    ,
                          a_Jacobian   ,
                          a_mapPartials);

  computeBdSmJnDiagonalTerm(a_smoothDiag ,
                            a_dx         ,
                            a_weights    ,
                            a_mapPartials,
                            a_Jacobian   ,
                            a_tangent    ,
                            a_tanDir     );
          
  // recompute magnetic field lines
  setBField(a_BField          ,
            a_curr            ,
            a_dx              ,
            a_geometryFileName,
            a_fieldFileName   ,
            a_blockName       );

  // recompute perpendicular field
  setPerpField(a_perpField       ,
               a_curr            , 
               a_BField          ,
               a_dx              ,
               a_geometryFileName,
               a_fieldFileName   ,
               a_blockName       ,
               a_writeInfo       );

  // recompute x and y derivatives
  setdXdY(a_perpXperpY ,
          a_perpField  ,
          a_mapPartials,
          a_dx         );
  
  if(a_writePlot) 
    { 
      int printFreq = 1;
      ParmParse pp;
      pp.get("printFreq_BrackbillSalzman",printFreq);
      if ((a_iter % printFreq) == 0)
        {  
          // 
          string prefix("diagnostic/");
          std::string suffix(".2d.hdf5");
                  
          // dummy default argument
          Vector<string> componentNames;
                  
          char name[1000];
                  
          // BField
          sprintf(name,"diagnostic/%03d.BField.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_BField,name,componentNames,a_dx[0]);

          // current weights
          sprintf(name,"diagnostic/%03d.weights.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_weights,name,componentNames,a_dx[0]);

          // perpField
          sprintf(name,"diagnostic/%03d.perpField.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_perpField,name,componentNames,a_dx[0]);
                  
          // x and y derivatives of perp field
          sprintf(name,"diagnostic/%03d.perpXperpY.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_perpXperpY,name,componentNames,a_dx[0]);

          // current map partials
          sprintf(name,"diagnostic/%03d.mapPartials.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_mapPartials,name,componentNames,a_dx[0]);

          // current smoothness
          sprintf(name,"diagnostic/%03d.smoothness.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_smoothness,name,componentNames,a_dx[0]);

          // current smoothDiag
          sprintf(name,"diagnostic/%03d.smoothDiag.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_smoothDiag,name,componentNames,a_dx[0]);
                  
          // current Jacobian
          sprintf(name,"diagnostic/%03d.Jacobian.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_Jacobian,name,componentNames,a_dx[0]);
        }
    }
}

void pointJacobiNew(FArrayBox                & a_curr            ,
                    FArrayBox                & a_BField          ,
                    FArrayBox                & a_BxBy            ,
                    FArrayBox                & a_perpField       ,
                    FArrayBox                & a_perpXperpY      , 
                    FArrayBox                & a_mapPartials     ,  
                    FArrayBox                & a_smoothness      ,
                    FArrayBox                & a_smoothDiag      ,
                    FArrayBox                & a_Jacobian        ,
                    FArrayBox                & a_weights         ,
                    Vector<FArrayBox*>       & a_arcLength       , 
                    const Vector<FArrayBox*> & a_targetSpacing   ,
                    Real                     & a_maxNorm         ,
                    const FArrayBox          & a_hiRes           , 
                    const Vector<int>        & a_tanDir          ,
                    const Vector<FArrayBox*> & a_tangent         , 
                    const Vector<FArrayBox*> & a_hiResArcLength  , 
                    const int                & a_iter            ,
                    const bool               & a_writeInfo       ,
                    const bool               & a_writePlot       ,
                    const RealVect           & a_dx              ,
                    const int                & a_baseRes         ,
                    const std::string        & a_geometryFileName,
                    const std::string        & a_fieldFileName   ,
                    const std::string        & a_blockName       ,
                    const Real               & a_xiEtaWeight     ,   
                    const Real               & a_sigma           ,   
                    const Real               & a_expo            ,
                    const int                & a_refRatio        ,
                    const int                & a_numInitSmooth   )
{
  //fieldData object, which is required for updating the magnetic field
    
  // Construct the field data object
  FieldData fieldData(a_geometryFileName,
                      a_fieldFileName   ,
                      a_blockName       );

  int numComp = a_curr.nComp();
  Box box = a_curr.box();
 
  // allocate increments
  FArrayBox inc      (box,numComp);
  FArrayBox smoothRes(box,numComp);
  FArrayBox totalRes (box,numComp);

  // SpaceDim components each for xi tangent and eta tangent
  FArrayBox volVarRes(box,numComp*SpaceDim);
  
  // set the increments to zero
  inc      .setVal(0.0);
  smoothRes.setVal(0.0);
  volVarRes.setVal(0.0);
  totalRes .setVal(0.0);
   

  // track whether the xPoint is on the hi side or low side with respect to eta
  int xPointHiLoSide = LARGEINTVAL;
  xPointHiLo(xPointHiLoSide,
             a_blockName   );

  // compute change (up to a factor of lambda) between two iterations
  computeIncNew(inc            ,
                smoothRes      ,
                volVarRes      ,
                totalRes       ,
                a_curr         , 
                a_mapPartials  ,
                a_smoothness   ,
                a_smoothDiag   ,
                a_Jacobian     ,
                a_weights      ,
                a_perpXperpY   ,
                a_perpField    ,
                a_BxBy         ,
                a_BField       ,
                a_writeInfo    ,
                a_dx           ,
                a_iter         ,
                a_numInitSmooth);

  // increment interior
  Box intBox = box;
  intBox.grow(-1);

  // comp
  int srcComp  = 0;
  int destComp = 0;

  ParmParse pp;
  bool moveInterior;
  pp.get("moveInterior",moveInterior);

  if (moveInterior)
    {
      // update curr interior
      a_curr.plus(inc     ,
                  intBox  ,
                  srcComp ,
                  destComp,
                  numComp );
    }

      // data holder for boundary increments
      Vector<FArrayBox*> arcLengthInc(4);
      for (int iEdge = 0; iEdge < 4; ++iEdge)
        {
          arcLengthInc [iEdge] = new FArrayBox(a_tangent[iEdge]->box(),1);
        }
    
  // norm of R,Z increments
  a_maxNorm = inc.norm(0,0,2);

  Real maxInc = 0.0;
  maxInc = Max(Max(Abs(inc.min(0)),Abs(inc.max(0))),Max(Abs(inc.min(1)),Abs(inc.max(1))));
          
  Real maxSmoothRes = 0.0;
  maxSmoothRes = Max(Max(Abs(smoothRes.min(0)),Abs(smoothRes.max(0))),Max(Abs(smoothRes.min(1)),Abs(smoothRes.max(1))));
          
  Real maxVolVarRes = 0.0;
  maxVolVarRes = Max(Max(Abs(volVarRes.min(0)),Abs(volVarRes.max(0))),Max(Abs(volVarRes.min(1)),Abs(volVarRes.max(1))));

  Real maxTotalRes = 0.0;
  maxTotalRes = Max(Max(Abs(totalRes.min(0)),Abs(totalRes.max(0))),Max(Abs(totalRes.min(1)),Abs(totalRes.max(1))));
          
  pout() << "max inc            : " << maxInc       << endl;
  pout() << "max smooth residual: " << maxSmoothRes << endl;
  pout() << "max volVar residual: " << maxVolVarRes << endl;
  pout() << "max total  residual: " << maxTotalRes  << endl;
        
  if(a_writePlot) 
    {  
      int printFreq = 1;
      ParmParse pp;
      pp.get("printFreq_BrackbillSalzman",printFreq);
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
              
          // volVar residual
          sprintf(name,"diagnostic/%03d.volVarRes.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&volVarRes,name,componentNames,a_dx[0]);
              
          // residual
          sprintf(name,"diagnostic/%03d.totalRes.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&totalRes ,name,componentNames,a_dx[0]);
              
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
                     a_dx         ,
                     a_writeInfo  );

  computeJacobianPartials(a_Jacobian   ,
                          a_mapPartials,
                          a_writeInfo  );
            
  computeSmJnDiagonalTerm(a_smoothDiag ,
                          a_dx         ,
                          a_weights    ,
                          a_Jacobian   ,
                          a_mapPartials);

  computeBdSmJnDiagonalTerm(a_smoothDiag ,
                            a_dx         ,
                            a_weights    ,
                            a_mapPartials,
                            a_Jacobian   ,
                            a_tangent    ,
                            a_tanDir     );
          
  // recompute magnetic field lines
  setBField(a_BField          ,
            a_curr            ,
            a_dx              ,
            a_geometryFileName,
            a_fieldFileName   ,
            a_blockName       );

  // recompute x and y derivatives
  setdXdY(a_BxBy       ,
          a_BField     ,
          a_mapPartials,
          a_dx         );
          
  // recompute perpendicular field
  setPerpField(a_perpField       ,
               a_curr            , 
               a_BField          ,
               a_dx              ,
               a_geometryFileName,
               a_fieldFileName   ,
               a_blockName       ,
               a_writeInfo       );

  // recompute x and y derivatives
  setdXdY(a_perpXperpY ,
           a_perpField,
           a_mapPartials,
           a_dx       );

  
          
  if(a_writePlot) 
    { 
      int printFreq = 1;
      ParmParse pp;
      pp.get("printFreq_BrackbillSalzman",printFreq);
      if ((a_iter % printFreq) == 0)
        {  
          // 
          string prefix("diagnostic/");
          std::string suffix(".2d.hdf5");
                  
          // dummy default argument
          Vector<string> componentNames;
                  
          char name[1000];
                  
          // BField
          sprintf(name,"diagnostic/%03d.BField.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_BField,name,componentNames,a_dx[0]);

          sprintf(name,"diagnostic/%03d.BxBy.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_BxBy,name,componentNames,a_dx[0]);

          // current weights
          sprintf(name,"diagnostic/%03d.weights.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_weights,name,componentNames,a_dx[0]);

          // perpField
          sprintf(name,"diagnostic/%03d.perpField.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_perpField,name,componentNames,a_dx[0]);
                  
          // x and y derivatives of perp field
          sprintf(name,"diagnostic/%03d.perpXperpY.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_perpXperpY,name,componentNames,a_dx[0]);

          // current map partials
          sprintf(name,"diagnostic/%03d.mapPartials.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_mapPartials,name,componentNames,a_dx[0]);

          // current smoothness
          sprintf(name,"diagnostic/%03d.smoothness.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_smoothness,name,componentNames,a_dx[0]);

          // current smoothDiag
          sprintf(name,"diagnostic/%03d.smoothDiag.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_smoothDiag,name,componentNames,a_dx[0]);
                  
          // current Jacobian
          sprintf(name,"diagnostic/%03d.Jacobian.iter.%06d.%s.2d.hdf5",a_baseRes,a_iter,a_blockName.c_str());
          writeFABname(&a_Jacobian,name,componentNames,a_dx[0]);
        }
    }
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


void interpolateInteriorXiEta(FArrayBox & a_initXiEta)
{
  
  // indexing macros
  int comp0 = 0;
  int comp1 = 1;

  int xDir = 0;
  int yDir = 1;

  // iterate over the interior
  Box domain      = a_initXiEta.box();

  IntVect loCorner = domain.smallEnd();
  IntVect hiCorner = domain.bigEnd();
  // fix for non square
  int loVal = loCorner[xDir];
  int hiVal = hiCorner[yDir];
  Box interiorBox = domain;
  interiorBox.grow(-1);

  // iterate over interior 
  for (BoxIterator bit(interiorBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      // xi stencil values
      IntVect loXi(iv[xDir],loVal);
      IntVect hiXi(iv[xDir],hiVal);
      Real xi0  = a_initXiEta(loXi ,comp0);
      Real xi1  = a_initXiEta(hiXi ,comp0);

      // eta stencil values
      IntVect loEta(loVal,iv[yDir]);
      IntVect hiEta(hiVal,iv[yDir]);
      Real eta0 = a_initXiEta(loEta,comp1);
      Real eta1 = a_initXiEta(hiEta,comp1);
      
      // intermediate quantities
      Real deltaXi  = xi1  - xi0 ;
      Real deltaEta = eta1 - eta0;
      Real denom    = 1.0 - deltaEta*deltaXi;

      // interpolated value
      RealVect retVal;
      retVal[0] = ( xi0 + deltaXi *eta0)/denom;
      retVal[1] = (eta0 + deltaEta* xi0)/denom;

      // assign to data holder
      a_initXiEta(iv,comp0) = retVal[0];
      a_initXiEta(iv,comp1) = retVal[1];           
    }
}



void sample(FArrayBox      & a_subSample,
            const FArrayBox& a_fineRes  )
{
  // box from coarser data holder
  Box coarseBox = a_subSample.box();
  
  // box iterate over finer data holder
  for (BoxIterator bit(coarseBox);bit.ok();++bit) 
    {
      IntVect iv = bit();
      for (int idir = 0; idir<SpaceDim; ++idir) 
        {
          // assign from even coord
          a_subSample(iv,idir) = a_fineRes(2*iv,idir);
        }
    }
}

void clearExtendedMappingFile(const char* a_fileName )
{
  
  FILE* fd;
  fd = fopen (a_fileName,"w");
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
                   const char*         a_fileName   )
{
  IntVect boxSize = a_extendedMap.box().size();
  int row = boxSize[0];
  int col = boxSize[1];
  int numRad = row - 2*a_numGhost[0];
  int numPol = col - 2*a_numGhost[1];

  pout()<<"Writing extended map for base resolution = (" << numRad<<","<<numPol<<")"<<endl;

  FILE* fd;
  fd = fopen (a_fileName, "a");
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
  else 
    {
      MayDay::Abort("bad block number");
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
  
  int numRow = a_currBox.bigEnd()[yDir];
  int numCol = a_currBox.bigEnd()[xDir];

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
          a_srcBlockName = "rcore";
        }
      else if (a_corner == "lowerRight")
        {
          lo             = upperRight + BASISV(xDir) - yGhost;
          a_srcBlockName = "rcore";
        }
      else if (a_corner == "upperRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperLeft")
        {
          if(a_blending)
            {
              lo             = lowerLeft - xGhost + BASISV(yDir);
              a_srcBlockName = "rcore";
            }
          else
            {
              extrapolate       = true;
              a_srcBlockName    = "invalid";
            }
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
          if(a_blending)
            {
              lo              = upperLeft - xGhost - yGhost;
              a_srcBlockName = "lcore";
            }
          else
            {
              extrapolate       = true;
              a_srcBlockName    = "invalid";
            }
        }
      else if (a_corner == "lowerRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperRight")
        {
          lo              = upperLeft + BASISV(xDir) + BASISV(yDir);                       
          a_srcBlockName = "rcsol";
        }
      else if (a_corner == "upperLeft")
        {
          lo              = lowerLeft - xGhost + BASISV(yDir);
          a_srcBlockName = "lcore";
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
          a_srcBlockName = "lcore";
        }
      else if (a_corner == "lowerRight")
        {
          lo              = upperRight + BASISV(xDir) - yGhost;
          a_srcBlockName = "rcsol";
        }
      else if (a_corner == "upperRight")
        {
          if(a_blending)
            {
              lo             = lowerRight + BASISV(xDir) + BASISV(yDir);
              a_srcBlockName = "lsol";
            }
          else
            {
              extrapolate       = true;
              a_srcBlockName    = "invalid";
            }
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
          if(a_blending)
            {
              lo             = upperRight + BASISV(xDir) - yGhost;
              a_srcBlockName = "rsol";
            }
          else
            {
              extrapolate       = true;
              a_srcBlockName    = "invalid";
            }
        }
      else if (a_corner == "upperRight")
        {
          lo              = lowerRight + BASISV(xDir) + BASISV(yDir);
          a_srcBlockName = "lcsol";
        }
      else if (a_corner == "upperLeft")
        {
          lo              = lowerLeft - xGhost + BASISV(yDir);
          a_srcBlockName = "lcsol";
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
          if(a_blending)
            {
              lo = upperRight + BASISV(xDir) - yGhost;
              a_srcBlockName = "lcsol";
            }
          else
            {
              extrapolate       = true;
              a_srcBlockName    = "invalid";
            }
        }
      else if (a_corner == "upperRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperLeft")
        {
          //lo = upperRight - xGhost + BASISV(yDir);
          //a_srcBlockName = "lpf";
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
          if(a_blending)
            {
              lo              = lowerRight + BASISV(xDir) + BASISV(yDir);
              a_srcBlockName = "rcsol";
            }
          else
            {
              extrapolate       = true;
              a_srcBlockName    = "invalid";
            }
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
          if (a_blending)
            {
              lo              = upperLeft - xGhost - yGhost;
              a_srcBlockName = "rpf";
            }
          else
            {
              extrapolate       = true;
              a_srcBlockName    = "invalid";
            }
        }
      else if (a_corner == "lowerRight")
        {
          extrapolate       = true;
          a_srcBlockName    = "invalid";
        }
      else if (a_corner == "upperRight")
        {
          //lo              = upperLeft + BASISV(xDir) + BASISV(yDir);
          //a_srcBlockName = "lsol";
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
          //lo             = lowerLeft + BASISV(xDir) - yGhost;
          //a_srcBlockName = "rsol";
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
          if(a_blending)
            {
              lo             = lowerLeft - xGhost + BASISV(yDir);
              a_srcBlockName = "lpf";
            }
          else
            {
              extrapolate       = true;
              a_srcBlockName    = "invalid";
            }
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
                       a_blockName                      ,
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
                       a_blockName                     ,
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

  for (int iEdge = 0; iEdge < 4; ++iEdge)
    {
      int dir = 0;
      int otherDir = 1;
      Side::LoHiSide side = Side::Lo;
      
      if (iEdge == 0 || iEdge == 2)
        {
          dir      = 1;
          otherDir = 0;
        }
      
      if (iEdge == 1 || iEdge == 2)
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
                    iEdge        );
      
      std::string sourceBlockName;
      getBlockName(sourceBlockName   ,
                   sourceBlockNumber);
      
      std::string edge = convertEdgeToString(iEdge); 
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

void computeArcLength(FArrayBox       & a_arcLength     ,
                      FArrayBox       & a_hiResArcLength,
                      Real            & a_evenSpacingVal,
                      Real            & a_totalArcLength,
                      const FArrayBox & a_hiRes         ,
                      const int       & a_tanDir        ,
                      const int       & a_refRatio      )
{
  int norDir = (a_tanDir + 1) % SpaceDim;
  int zeroComp = 0;
  Box edgeBoxInt = a_hiResArcLength.box();
  edgeBoxInt.growLo(a_tanDir,-1);
  
  // incremental counters
  Real totalHiResLength = 0.0;
  Real totalLoResLength = 0.0;
  Real loResSegLength   = 0.0;
 
  for (BoxIterator bit(edgeBoxInt);bit.ok();++bit) 
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
      if ((iv[a_tanDir] % a_refRatio) == 0)
        {
          CH_assert(iv[norDir] % a_refRatio == 0);
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

void computeBoundaryTangentVector(FArrayBox       & a_tangent,
                                  const FArrayBox & a_curr   ,
                                  const int       & a_tanDir ,
                                  const RealVect  & a_dx     )
{
  Box edgeBoxInt = a_tangent.box();
  edgeBoxInt.grow(-BASISV(a_tanDir));
  for (BoxIterator bit(edgeBoxInt);bit.ok();++bit) 
    {
      IntVect iv = bit();
      for (int iComp = 0; iComp < SpaceDim; ++ iComp)
        {
          a_tangent(iv,iComp) = a_curr(iv + BASISV(a_tanDir),iComp) - a_curr(iv - BASISV(a_tanDir),iComp);
          a_tangent(iv,iComp) /= (2.0*a_dx[a_tanDir]);
        }
    }
  
  // paramteters for calculating tangent vector in the corner
  Real hiLoSide = LARGEREALVAL;
  Real evalPt = 0.0;
  Real deriv1 = LARGEREALVAL;
  Real deriv2 = LARGEREALVAL;
  
  // calculate tangent vector in the low corner through quadratic extrapolation
  IntVect loCorner = a_tangent.box().smallEnd();
  hiLoSide = -1;
  for (int iComp = 0;  iComp < SpaceDim; ++iComp)
    {
      Vector<Real> xVal(3);
      Vector<Real> yVal(3);
      
      // x-values for quadratic interpolation
      xVal[0] = evalPt;
      xVal[1] = evalPt - hiLoSide*    a_dx[a_tanDir];
      xVal[2] = evalPt - hiLoSide*2.0*a_dx[a_tanDir];
      
      // y-values for quadratic interpolation
      yVal[0] = a_curr(loCorner,iComp);
      yVal[1] = a_curr(loCorner - hiLoSide*     BASISV(a_tanDir),iComp);
      yVal[2] = a_curr(loCorner - hiLoSide* 2 * BASISV(a_tanDir),iComp);
      
      // differentiate the quadratic
      oneSidedPureDerivs(deriv1,
                         deriv2, 
                         evalPt,
                         yVal  ,
                         xVal  );

      // deriv2 comes for free, but not used here
      a_tangent(loCorner,iComp) = deriv1;
      //pout() << "loCorner deriv = "<<deriv1<<endl;
    }

  // calculate tangent vector in the hi corner through quadratic extrapolation
  IntVect hiCorner = a_tangent.box().bigEnd  ();
  hiLoSide = 1;
  for (int iComp = 0;  iComp < SpaceDim; ++iComp)
    {
      Vector<Real> xVal(3);
      Vector<Real> yVal(3);
      
      // x-values for quadratic interpolation
      xVal[0] = evalPt;
      xVal[1] = evalPt - hiLoSide*    a_dx[a_tanDir];
      xVal[2] = evalPt - hiLoSide*2.0*a_dx[a_tanDir];
      
      // y-values for quadratic interpolation
      yVal[0] = a_curr(hiCorner,iComp);
      yVal[1] = a_curr(hiCorner - hiLoSide*     BASISV(a_tanDir),iComp);
      yVal[2] = a_curr(hiCorner - hiLoSide* 2 * BASISV(a_tanDir),iComp);
      
      // differentiate the quadratic
      oneSidedPureDerivs(deriv1,
                         deriv2, 
                         evalPt,
                         yVal  ,
                         xVal  );

      // deriv2 comes for free, but not used here
      a_tangent(hiCorner,iComp) = deriv1;
      //pout() << "hiCorner deriv = "<<deriv1<<endl;
    }
  for (BoxIterator bit(a_tangent.box());bit.ok();++bit) 
    {
      IntVect iv = bit();
      Real len = 0.0;
      for (int iComp = 0; iComp < SpaceDim; ++ iComp)
        {
          len += a_tangent(iv,iComp)*a_tangent(iv,iComp);
        }
      len = sqrt(len);
      if (len > 0.0)
        {
          for (int iComp = 0; iComp < SpaceDim; ++ iComp)
            {
              a_tangent(iv,iComp) /= len;
            }
        }
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

void checkPoint(const int                & a_beginBlock        ,
                const int                & a_endBlock          ,
                const int                & a_numBlocks         ,
                const std::string        & a_geometryFileName  ,
                const std::string        & a_fieldFileName     ,
                const Vector<FArrayBox*> & a_curr              ,
                const int                & a_numRad            ,
                const Vector<IntVect>    & a_hiCorner          ,
                const Vector<FArrayBox*> & a_BField            ,
                const Vector<RealVect>   & a_dx                ,
                const int                & a_numGhost          ,
                const int                & a_iter              ,
                const std::string        & a_checkpointBaseName,
                const std::string        & a_filenameSuffix    )
  
{
  //int xiDir  = 0;
  //int etaDir = 1;
  for (int iBlock = 0; iBlock < a_beginBlock; ++iBlock)
    {
      
      std::string blockName;
      getBlockName(blockName,
                   iBlock  );
      
      // Construct the field data object
      FieldData field_data(a_geometryFileName,
                           a_fieldFileName   ,
                           blockName          );
      
      
      // mapped coordinates of the unit square
      setPhysicalCoordinates(*(a_curr  [iBlock]) ,
                             a_hiCorner[iBlock]  ,
                             blockName           ,
                             a_geometryFileName  ,
                             a_fieldFileName     );
      
      
      // use FieldData class for magnetic field and first partials
      setBField(*(a_BField[iBlock]),
                *(a_curr  [iBlock]),
                a_dx      [iBlock] ,
                a_geometryFileName ,
                a_fieldFileName    ,
                blockName          );

    }
 
  for (int iBlock = a_endBlock + 1; iBlock < a_numBlocks; ++iBlock)
    {
      std::string blockName;
      getBlockName(blockName,
                   iBlock   );
      
      // Construct the field data object
      FieldData field_data(a_geometryFileName,
                           a_fieldFileName   ,
                           blockName       );
      
      
      // mapped coordinates of the unit square
      setPhysicalCoordinates(*(a_curr  [iBlock]),
                             a_hiCorner[iBlock] ,
                             blockName          ,
                             a_geometryFileName ,
                             a_fieldFileName    );
      
      // use FieldData class for magnetic field and first partials
      setBField(*(a_BField[iBlock]),
                *(a_curr  [iBlock]),
                a_dx      [iBlock] ,
                a_geometryFileName ,
                a_fieldFileName    ,
                blockName          );
    }
  
  // clean start for mapping file output function  
  char fileName[1000];

  // assumes no ghost cells
  sprintf(fileName,"%s.%d.",a_checkpointBaseName.c_str(),a_iter);
  std::string name = fileName;
    
  clearExtendedMappingFile(name.c_str());
  
  for (int iBlock = 0; iBlock < a_numBlocks; ++iBlock)
    {
      std::string blockName;
      getBlockName(blockName,
                   iBlock   );
            
      outputMapping(*(a_curr   [iBlock]),
                    *(a_BField [iBlock]),
                    IntVect::Zero       ,
                    blockName           ,
                    name.c_str()        );
    }

}

void computeFluxMatchedVal(FArrayBox                 & a_fluxMatchedVal  ,
                           const FArrayBox           & a_evenSpaceCurr   ,
                           const Vector<FArrayBox *> & a_hiResArcLength  ,
                           const FArrayBox           & a_hiRes           ,
                           const Vector<Real>        & a_evenSpacingVal  ,
                           const int                 & a_iEdge           ,
                           const int                 & a_iBlock          ,
                           const RealVect            & a_hiResDx         ,
                           const int                 & a_refRatio        ,
                           const std::string         & a_geometryFileName,
                           const std::string         & a_fieldFileName   )
{
  Real totalLength = 0.0;
  int zeroComp = 0;
  int xComp = 0;
  int yComp = 1;
  
  std::string blockName;
  getBlockName(blockName,
               a_iBlock  );
  
  // Construct the field data object
  FieldData field_data(a_geometryFileName,
                       a_fieldFileName   ,
                       blockName         );

    int loHi = 1;
  if ((a_iEdge == 0) || a_iEdge == 1)
    {
      loHi = 0;
    }

  int farEdge = (a_iEdge + 2) % 4;
  int tanDir = 1;
  if ((a_iEdge == 0) || (a_iEdge == 2))
  {
    tanDir = 0;
  }
  int norDir = (tanDir + 1)%2;
  
  //Real evenSpacingValFar = a_evenSpacingVal[farEdge];

  const FArrayBox& hiResArcLengthNear = *(a_hiResArcLength[a_iEdge]);
  const FArrayBox& hiResArcLengthFar  = *(a_hiResArcLength[farEdge]);

  Box nearEdgeBox = hiResArcLengthFar.box();
  Box farEdgeBox  = hiResArcLengthFar.box();
    
  int numCells      = a_fluxMatchedVal  .box().size()[tanDir] - 1;
  int numHiResCells = hiResArcLengthNear.box().size()[tanDir] - 1;
  bool hasXPt;
  IntVect xPtIv = findXiEtaXPt(hasXPt,
                               a_iBlock);
  
  bool xPointEdge = false;

  //Real xPt;
  if (hasXPt)
    {
      std::string edge = convertEdgeToString(a_iEdge); 
      bool xPointHiSide;
      
      std::string blockName;
      getBlockName(blockName,
                   a_iBlock );
      
      xPointEdge = edgeContainsXPoint(xPointHiSide,
                                      blockName   ,
                                      edge        );
    }
  
  int currHiRes = 0;
  Real currArcLength = 0.0;
  Box box = a_fluxMatchedVal.box();
  IntVect smallEnd = box.smallEnd();
  IntVect bigEnd   = box.bigEnd  ();

  if (((a_iEdge == 0) || (a_iEdge == 2)) && xPointEdge)
   {
     for (BoxIterator bit(box);bit.ok();++bit) 
       {
         IntVect iv = bit();
         if (iv != smallEnd)
           {
             if (iv != bigEnd)
               {
                 IntVect farIv = iv + (1 - loHi)*numCells*BASISV(norDir);
                 Real x = a_evenSpaceCurr(farIv,xComp);
                 Real y = a_evenSpaceCurr(farIv,yComp);
                 RealVect farPt(x,y);
                 Real farFlux = field_data.getMagneticFluxFromDCT(farPt);
                 
                 // iterate through fine data holders to find bracketing fluxes
                 for (int hiResCounter = currHiRes; hiResCounter < numHiResCells; ++hiResCounter) 
                   {
                     // condition for exiting the loop
                     bool foundBracketing = false;
                     
                     // lo member of bracketing ivs
                     IntVect hiResIvLo;
                     hiResIvLo[tanDir] = hiResCounter;
                     hiResIvLo[norDir] = loHi*numHiResCells;
                     
                     // flux for lo iv
                     Real xLo = a_hiRes(hiResIvLo,xComp);
                     Real yLo = a_hiRes(hiResIvLo,yComp);
                     RealVect hiResPtLo(xLo,yLo);
                     Real nearFluxLo = field_data.getMagneticFluxFromDCT(hiResPtLo);
                     
                     // hi member of bracketing ivs
                     IntVect hiResIvHi = hiResIvLo + BASISV(tanDir);
                     
                     // flux for hi iv
                     Real xHi = a_hiRes(hiResIvHi,xComp);
                     Real yHi = a_hiRes(hiResIvHi,yComp);
                     RealVect hiResPtHi(xHi,yHi);
                     Real nearFluxHi = field_data.getMagneticFluxFromDCT(hiResPtHi);
                     
                     Real fluxDiffHi = nearFluxHi - farFlux;
                     Real fluxDiffLo = nearFluxLo - farFlux;
                     
                     Real tol = 1.0e-10;

                     // condition to check for bracketing
                     if (fluxDiffHi*fluxDiffLo <= 0.0)
                       {
                         RealVect newHi = hiResPtHi;
                         RealVect newLo = hiResPtLo;
                             
                         RealVect newPt = LARGEREALVAL*RealVect::Unit;
                         int counter = 0;
                         while (((Abs(fluxDiffHi)  > tol) || (Abs(fluxDiffLo) > tol)) && (counter < 1000))
                           {
                             // mid-point
                             newPt = 0.5*(newHi + newLo);
                             
                             // flux at mid-point
                             Real newFlux = field_data.getMagneticFluxFromDCT(newPt);
                             
                             // use the flux at mid-point to further restrict
                             Real fluxDiffNew = newFlux - farFlux;
                             
                             if (fluxDiffNew*fluxDiffLo <= 0.0)
                               {
                                 newLo      = newLo;
                                 newHi      = newPt;
                                 fluxDiffHi = fluxDiffNew;
                               }
                             else if (fluxDiffNew*fluxDiffHi <= 0.0)
                               {
                                 newLo      = newPt;
                                 newHi      = newHi;
                                 fluxDiffLo = fluxDiffNew;
                               }
                             else
                               {
                                 FILE* fd;
                                 fd = fopen ("done","w");
                                 if (fd == NULL) 
                                   {
                                     perror ("Error opening file");
                                   }
                                 fclose(fd);
                                 MayDay::Abort("not bracketed");
                               }

                             // double check in case the process hangs
                             counter += 1;
                             if (counter == 1000)
                               {
                                 pout()<<"bisection failed"<<endl;
                                 FILE* fd;
                                 fd = fopen ("done","w");
                                 if (fd == NULL) 
                                   {
                                     perror ("Error opening file");
                                   }
                                 fclose(fd);
                                 MayDay::Abort("bisection failed");
                               }
                           }
                         
                         Real dist2Hi  = sqrt((newPt[0] - hiResPtHi[0])*(newPt[0] - hiResPtHi[0]) + (newPt[1] - hiResPtHi[1])*(newPt[1] - hiResPtHi[1])); 
                         Real dist2Lo  = sqrt((newPt[0] - hiResPtLo[0])*(newPt[0] - hiResPtLo[0]) + (newPt[1] - hiResPtLo[1])*(newPt[1] - hiResPtLo[1]));
                                                                                                              
                         Real inc      = currArcLength + dist2Lo;
                         a_fluxMatchedVal(iv,zeroComp) = inc;
                         totalLength                  += inc;
                         currArcLength = dist2Hi;

                         pout()<<currArcLength<<endl;
                         currHiRes = hiResCounter + 1;
                         foundBracketing = true;
                       }
                     
                     if (foundBracketing)
                       {
                         break;
                       }
                     
                     // start looking for bracketing ivs at this index
                     currHiRes = hiResCounter + 1;
                     //Real temp = hiResArcLengthNear(hiResIvHi,zeroComp);
                     currArcLength += hiResArcLengthNear(hiResIvHi,zeroComp);
                   }
               }
             else
               {
                  for (int hiResCounter = currHiRes; hiResCounter < numHiResCells; ++hiResCounter) 
                   {
                     IntVect hiResIvHi;
                     hiResIvHi[tanDir] = hiResCounter + 1;
                     hiResIvHi[norDir] = loHi*numHiResCells;
                     currArcLength += hiResArcLengthNear(hiResIvHi,zeroComp);
                   }
                  
                  a_fluxMatchedVal(iv,zeroComp) = currArcLength;
                  totalLength                  += currArcLength;
               }
           }
       }
   }
 else
   {
     a_fluxMatchedVal.setVal(LARGEREALVAL);
   }
  pout()<<"totalLength = "<<totalLength <<endl;
}

void computeTargetSpacingVal(FArrayBox       & a_targetSpacing           ,
                             Real            & a_totalTargetSpacingLength,
                             const FArrayBox & a_arcLength               ,
                             const Real      & a_evenSpacingVal          ,
                             const FArrayBox & a_fluxMatchedVal          ,
                             const int       & a_iEdge                   ,
                             const Real      & a_sigma                   ,
                             const Real      & a_expo                    ,
                             const int       & a_iBlock                  ,
                             const RealVect  & a_dx                      ,
                             const int       & a_tanDir                  ,
                             const Real      & a_totalArcLength          )
{
  // macros 
  int evenComp = 0;
  int fluxComp = 1;
  int twoComp  = 2;
  
  // initialize to 0.0
  a_totalTargetSpacingLength = 0.0;
  
  int numCells = a_targetSpacing.box().size()[a_tanDir] - 1;
  FArrayBox cumulativeInc(a_fluxMatchedVal.box(),twoComp);
  int index = -1;
  Real cumulativeFluxInc = 0.0;
  for (BoxIterator bit(a_targetSpacing.box());bit.ok();++bit) 
     {
       index += 1;
       IntVect iv = bit();
       if (iv[a_tanDir] ==  0)
         {
           cumulativeInc(iv,evenComp) = index*a_evenSpacingVal;
           cumulativeInc(iv,fluxComp) = cumulativeFluxInc     ;
         }
       else 
         {
           cumulativeFluxInc         += a_fluxMatchedVal(iv,0);
           cumulativeInc(iv,evenComp) = index*a_evenSpacingVal;
           cumulativeInc(iv,fluxComp) = cumulativeFluxInc     ;
           
         }
     }
  bool hasXPt;

  // indicates lo corner hi corner info
  IntVect xPtLogicalIv = findXiEtaXPt(hasXPt,
                                      a_iBlock);

  IntVect xPtIv = LARGEINTVAL*IntVect::Unit; 
  Real xPt;

  bool xPointEdge = false;
  
  if (hasXPt)
    {
      std::string edge = convertEdgeToString(a_iEdge); 
      bool xPointHiSide;
      
      std::string blockName;
      getBlockName(blockName,
                   a_iBlock );
      
      xPointEdge = edgeContainsXPoint(xPointHiSide,
                                      blockName   ,
                                      edge        );
      
      if (xPointEdge && !xPointHiSide)
        {
          CH_assert(xPtLogicalIv[a_tanDir] == 0);
        }

      if (xPointEdge)
        {
          xPtIv = xPtLogicalIv*numCells         ;
          xPt   = xPtIv[a_tanDir]*a_dx[a_tanDir];
        }
    }
  
  for (BoxIterator bit(a_targetSpacing.box());bit.ok();++bit) 
     {
       IntVect iv = bit();
       if (iv[a_tanDir] > 0)
         {
           Real pt = (iv[a_tanDir] - 0.5)*a_dx[a_tanDir];
           
           if ((a_iEdge == 0) || (a_iEdge == 2))
             {
               if (!xPointEdge)
                 {
                   a_targetSpacing(iv,0)       = a_evenSpacingVal;
                   a_totalTargetSpacingLength += a_evenSpacingVal;
                 }
               else
                 {
                   IntVect preIv = iv - BASISV(a_tanDir);
                   Real    prePt = (preIv[a_tanDir] - 0.5)*a_dx[a_tanDir];
                   
                   Real dist    = (pt    - xPt)*(pt    - xPt);
                   Real preDist = (prePt - xPt)*(prePt - xPt);
                   
                   dist    = sqrt(dist);
                   preDist = sqrt(preDist);
                   Real weight    = LARGEREALVAL;
                   Real preWeight = LARGEREALVAL;
                   if (a_sigma <= 0.0)
                     {
                       FILE* fd;
                       fd = fopen ("done","w");
                       if (fd == NULL) 
                         {
                           perror ("Error opening file");
                         }
                       fclose(fd);
                       MayDay::Abort("sigma must be positive");
                     }
                   else
                     {
                       weight    = exp(-pow(   dist/a_sigma,a_expo));
                       preWeight = exp(-pow(preDist/a_sigma,a_expo));
                     }
                   
                   // check the weight
                   pout()<<"weight = "<<weight<<endl;

                   a_targetSpacing(iv,0)       = weight   *cumulativeInc(iv   ,evenComp) + (1.0 - weight   )*cumulativeInc(iv   ,fluxComp);
                   a_targetSpacing(iv,0)      -= preWeight*cumulativeInc(preIv,evenComp) + (1.0 - preWeight)*cumulativeInc(preIv,fluxComp);
                   a_totalTargetSpacingLength += a_targetSpacing(iv,0);
                 }
             }
           else
             {
               if (a_iBlock < 8)
                 {
                   a_targetSpacing(iv,0)       = a_evenSpacingVal;
                   a_totalTargetSpacingLength += a_evenSpacingVal;
                 }
               else
                 {
                   a_targetSpacing(iv,0)       = a_arcLength(iv,0);
                   a_totalTargetSpacingLength += a_arcLength(iv,0);
                 }
             }
           
         }
     }
}

void straightenEtaBd(Vector <FArrayBox* > & a_curr            ,
                     Vector <FArrayBox* > & a_BField          ,
                     const std::string    & a_geometryFileName,
                     const std::string    & a_fieldFileName   ,
                     const int            & a_numStraighten   ,
                     const RealVect       & a_dx              ,
                     const int            & a_numBlocks       )
{
  int xiDir  = 0;
  int etaDir = 1;
  
  int xDir = 0;
  int yDir = 1;
  
  // initialize ghost cells to 0 and interior to map
  for(int iBlock = 0; iBlock < a_numBlocks; ++iBlock)
    {
      Box box          = a_curr[iBlock]->box();
      IntVect hiCorner = box.bigEnd();
      int numCol       = hiCorner[etaDir];
      int numRow       = hiCorner[xiDir] ;

      // set to arbitrary value
      a_BField[iBlock]->setVal(0.13579);
      
      std::string blockName;
      getBlockName(blockName,
                   iBlock  );
      
      // Construct the field data object
      FieldData field_data(a_geometryFileName,
                           a_fieldFileName   ,
                           blockName       );
      
      FArrayBox& currFab =  *a_curr[iBlock];
      setPhysicalCoordinates(currFab         ,
                             hiCorner        ,
                             blockName       ,
                             a_geometryFileName,
                             a_fieldFileName   );
      
      if (iBlock < a_numStraighten)
        {
          // corner IntVects
          IntVect lo = box.smallEnd();
          IntVect hi = box.bigEnd  ();
          IntVect loHi(lo[0],hi[1]);
          IntVect hiLo(hi[0],lo[1]);
          
          // lo physical coord
          Real loValx = currFab(lo,xiDir );
          Real loValy = currFab(lo,etaDir);
          
          // hi  physical coord
          Real hiValx = currFab(hi,xiDir );
          Real hiValy = currFab(hi,etaDir);
          
          // loHi physical coord
          Real loHiValx = currFab(loHi,xiDir );
          Real loHiValy = currFab(loHi,etaDir);
          
          // hiLo physical coord 
          Real hiLoValx = currFab(hiLo,xiDir );
          Real hiLoValy = currFab(hiLo,etaDir);
          
          // four corners of the map
          RealVect loVal  (loValx  ,loValy  );
          RealVect hiVal  (hiValx  ,hiValy  );
          RealVect loHiVal(loHiValx,loHiValy);
          RealVect hiLoVal(hiLoValx,hiLoValy);
          
          for (BoxIterator bit(box);bit.ok();++bit) 
            {
              IntVect iv = bit();
              if ( (iv != lo) && (iv != hi) && (iv != loHi) && (iv != hiLo) )
                {
                  Real frac = iv[xiDir];
                  frac     /= numCol    ;
                    
                  // only change the bottom or top edges
                  if(iv[etaDir] == 0)
                    {
                      Real x1 = currFab(iv + BASISV(etaDir),xDir);
                      Real y1 = currFab(iv + BASISV(etaDir),yDir);
                        
                      Real x2 = currFab(iv + 2*BASISV(etaDir),xDir);
                      Real y2 = currFab(iv + 2*BASISV(etaDir),yDir);
                        
                      Real x3 = hiLoVal[xDir];
                      Real y3 = hiLoVal[yDir];
                        
                      Real x4 = loVal[xDir];
                      Real y4 = loVal[yDir];
                        
                      Real a = x3 - x4;
                      Real b = x2 - x1;
                      Real c = y3 - y4;
                      Real d = y2 - y1;
                        
                      Real b1 = x1 - x4;
                      Real b2 = y1 - y4;
                      Real disc = a*d - b*c;
                        
                      Real s = (1.0/disc)*( d*b1 - b*b2);
                      Real t = (1.0/disc)*(-c*b1 + a*b2);
                        
                      Real xVal1 = x4 + s*(x3 - x4);
                      Real yVal1 = y4 + s*(y3 - y4);
                        
                      Real xVal2 = x1 + t*(x1 - x2);
                      Real yVal2 = y1 + t*(y1 - y2);

                      RealVect P1(x1,y1);
                      RealVect P2(x2,y2);
                      RealVect P3(x3,y3);
                      RealVect P4(x4,y4);
                      RealVect V1(xVal1,yVal1);
                      RealVect V2(xVal2,yVal2);

                      RealVect D1 = V1 - P1;
                      RealVect D2 = V1 - P2;
                      RealVect D3 = V1 - P3;
                      RealVect D4 = V1 - P4;

                      RealVect E1 = V2 - P1;
                      RealVect E2 = V2 - P2;
                      RealVect E3 = V2 - P3;
                      RealVect E4 = V2 - P4;

                      Real Dratio1 = D1[1]/D1[0];
                      Real Dratio2 = D2[1]/D2[0];
                      Real Dratio3 = D3[1]/D3[0];
                      Real Dratio4 = D4[1]/D4[0];

                      Real Eratio1 = E1[1]/E1[0];
                      Real Eratio2 = E2[1]/E2[0];
                      Real Eratio3 = E3[1]/E3[0];
                      Real Eratio4 = E4[1]/E4[0];

                      Real Disc1 = Dratio3 - Dratio4;
                      Real Disc2 = Dratio2 - Dratio1;

                      Real Eisc1 = Eratio3 - Eratio4;
                      Real Eisc2 = Eratio2 - Eratio1;
                       
                      if ( (Abs (Disc1) > 1.0e-8) || (Abs (Disc2) > 1.0e-8) ||(Abs (Eisc1) > 1.0e-8) || (Abs (Eisc2) > 1.0e-8) )
                        {
                          pout()<<" Disc1 = "<<Disc1 << endl;
                          pout()<<" Disc2 = "<<Disc2 << endl;
                            
                          pout()<<" Eisc1 = "<<Eisc1 << endl;
                          pout()<<" Eisc2 = "<<Eisc2 << endl;
                            
                          pout() <<"iv = "<< iv << endl;
                            
                          MayDay::Abort("new point not at intersection");
                        }
                      pout()<<" Disc1 = "<<Disc1 << endl;
                      pout()<<" Disc2 = "<<Disc2 << endl;
                        
                      pout()<<" Eisc1 = "<<Eisc1 << endl;
                      pout()<<" Eisc2 = "<<Eisc2 << endl;
                        
                      pout() <<"iv = "<< iv << endl;
                      currFab(iv,xDir) = xVal1;
                      currFab(iv,yDir) = yVal1;
                    }
                  else if(iv[etaDir] == numRow)
                    {
                      Real x1 = currFab(iv - BASISV(etaDir),xDir);
                      Real y1 = currFab(iv - BASISV(etaDir),yDir);
                        
                      Real x2 = currFab(iv - 2*BASISV(etaDir),xDir);
                      Real y2 = currFab(iv - 2*BASISV(etaDir),yDir);
                        
                      Real x3 = hiVal[xDir];
                      Real y3 = hiVal[yDir];
                        
                      Real x4 = loHiVal[xDir];
                      Real y4 = loHiVal[yDir];
                        
                      Real a = x3 - x4;
                      Real b = x2 - x1;
                      Real c = y3 - y4;
                      Real d = y2 - y1;
                        
                      Real b1 = x1 - x4;
                      Real b2 = y1 - y4;
                      Real disc = a*d - b*c;
                        
                      Real s = (1.0/disc)*( d*b1 - b*b2);
                      Real t = (1.0/disc)*(-c*b1 + a*b2);
                        
                      Real xVal1 = x4 + s*(x3 - x4);
                      Real yVal1 = y4 + s*(y3 - y4);
                        
                      Real xVal2 = x1 + t*(x1 - x2);
                      Real yVal2 = y1 + t*(y1 - y2);
                                                
                      RealVect P1(x1,y1);
                      RealVect P2(x2,y2);
                      RealVect P3(x3,y3);
                      RealVect P4(x4,y4);
                      RealVect V1(xVal1,yVal1);
                      RealVect V2(xVal2,yVal2);

                      RealVect D1 = V1 - P1;
                      RealVect D2 = V1 - P2;
                      RealVect D3 = V1 - P3;
                      RealVect D4 = V1 - P4;

                      RealVect E1 = V2 - P1;
                      RealVect E2 = V2 - P2;
                      RealVect E3 = V2 - P3;
                      RealVect E4 = V2 - P4;

                      Real Dratio1 = D1[1]/D1[0];
                      Real Dratio2 = D2[1]/D2[0];
                      Real Dratio3 = D3[1]/D3[0];
                      Real Dratio4 = D4[1]/D4[0];

                      Real Eratio1 = E1[1]/E1[0];
                      Real Eratio2 = E2[1]/E2[0];
                      Real Eratio3 = E3[1]/E3[0];
                      Real Eratio4 = E4[1]/E4[0];

                      Real Disc1 = Dratio3 - Dratio4;
                      Real Disc2 = Dratio2 - Dratio1;

                      Real Eisc1 = Eratio3 - Eratio4;
                      Real Eisc2 = Eratio2 - Eratio1;
                       
                      if ( (Abs (Disc1) > 1.0e-8) || (Abs (Disc2) > 1.0e-8) ||(Abs (Eisc1) > 1.0e-8) || (Abs (Eisc2) > 1.0e-8) )
                        {
                          pout()<<" Disc1 = "<<Disc1 << endl;
                          pout()<<" Disc2 = "<<Disc2 << endl;
                            
                          pout()<<" Eisc1 = "<<Eisc1 << endl;
                          pout()<<" Eisc2 = "<<Eisc2 << endl;

                          pout() <<"iv = "<< iv << endl;
                            
                          MayDay::Abort("new point not at intersection");
                        }
                      pout()<<" Disc1 = "<<Disc1 << endl;
                      pout()<<" Disc2 = "<<Disc2 << endl;
                        
                      pout()<<" Eisc1 = "<<Eisc1 << endl;
                      pout()<<" Eisc2 = "<<Eisc2 << endl;
                        
                      pout() <<"iv = "<< iv << endl;
                      currFab(iv,xDir) = xVal2;
                      currFab(iv,yDir) = yVal2;
                    }
                }
            }
        }
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
        pout() << "Usage:  bs...ex <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }
   
    // Parse input file
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);
    
    // macros for indexing
    int xDir = 0;
    int yDir = 1;
    int oneComp = 1;

    //begin loop over resolution here
    int log2MinRes = LARGEINTVAL;
    int log2MaxRes = LARGEINTVAL;
    pp.get("log2MinRes_BrackbillSalzman",log2MinRes);
    pp.get("log2MaxRes_BrackbillSalzman",log2MaxRes);

    // value for stopping criterion
    Real maxTol = -LARGEREALVAL;
    pp.get("maxTol",maxTol);
    
    // parameter for point Jacobi 
    int maxIter_BrackbillSalzman = LARGEREALVAL;
    pp.get("maxIter_BrackbillSalzman",maxIter_BrackbillSalzman) ;
    
    // parameter for smoothnessRes 
    Real xiEtaWeight = LARGEREALVAL;
    pp.get("xiEtaWeight",xiEtaWeight) ;
    
    // prolixity 
    bool writeInfo = false;
    pp.query("writeInfo",writeInfo) ;
    
    // prolixity 
    bool writePlot = false;
    pp.query("writePlot",writePlot) ;

    // checkpoint interval
    int checkPointInt = LARGEREALVAL;
    pp.get("checkPointInt",checkPointInt) ;
    
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
    
    // box for hi res discretization of the boundary
    int refRatio = LARGEREALVAL;
    pp.get("refRatio",refRatio);

    for (int numberCells = log2MinRes; numberCells <= log2MaxRes; ++ numberCells)
      {
        // resolution before refRatioXi and refRatioEta
        int baseRes = pow(2,numberCells);

        // different box data for each block
        Vector<int> numRow           (numBlocks); 
        Vector<int> numCol           (numBlocks);
        Vector<IntVect> nCell        (numBlocks);
        Vector<IntVect> loEnd        (numBlocks);
        Vector<IntVect> hiEnd        (numBlocks);
        Vector<Box> box              (numBlocks);
        Vector<IntVect> hiCorner     (numBlocks);
        Vector<RealVect> dx          (numBlocks);

        Vector<Box> hiResBox         (numBlocks); 
        Vector<IntVect> hiResHiCorner(numBlocks);
        Vector<int> hiResNumRow      (numBlocks);
        Vector<int> hiResNumCol      (numBlocks);
        Vector<RealVect> hiResDx     (numBlocks);
        
        // boxes coresponding to boundary values on four sides
        Vector<Vector<Box> >edgeBoxes     (numBlocks);
        Vector<Vector<Box> >hiResEdgeBoxes(numBlocks);

        for (int iBlock = 0; iBlock < numBlocks;++ iBlock)
          {
            numRow[iBlock] = pow(2,numberCells + refRatioEta[iBlock]);
            numCol[iBlock] = pow(2,numberCells + refRatioXi [iBlock] );

            nCell[iBlock][0] = numCol[iBlock];
            nCell[iBlock][1] = numRow[iBlock];

            loEnd[iBlock] = IntVect::Zero;
            hiEnd[iBlock][0] = numCol[iBlock];
            hiEnd[iBlock][1] = numRow[iBlock];

            box     [iBlock].define(loEnd[iBlock],hiEnd         [iBlock]);
            hiCorner[iBlock] = box[iBlock].bigEnd();

            dx[iBlock][0] = 1.0/((Real)numCol[iBlock]);
            dx[iBlock][1] = 1.0/((Real)numRow[iBlock]);

            // hi res versions
            hiResBox     [iBlock].define(loEnd     [iBlock],refRatio*hiEnd[iBlock]);
            hiResHiCorner[iBlock] = hiResBox       [iBlock].bigEnd();
            hiResNumRow  [iBlock] = refRatio*numRow[iBlock];
            hiResNumCol  [iBlock] = refRatio*numCol[iBlock];

            hiResDx[iBlock][0] = 1.0/((Real)hiResNumCol[iBlock]);
            hiResDx[iBlock][1] = 1.0/((Real)hiResNumRow[iBlock]);

            edgeBoxes     [iBlock].resize(4);
            hiResEdgeBoxes[iBlock].resize(4);
            
            // boxes used for tangent vector FArrayBoxes
            createEdgeBoxes(edgeBoxes     [iBlock],box     [iBlock]);
            createEdgeBoxes(hiResEdgeBoxes[iBlock],hiResBox[iBlock]);
          }
                        
        int beginBlock = -LARGEINTVAL;
        pp.get("beginBlockBS", beginBlock);

        int endBlock = -LARGEINTVAL;
        pp.get("endBlockBS", endBlock);
        
        // mapped quantities have two components
        int numComp    = 2 ;
        int numJacComp = 11;
        int numWeights = 2 ;      

        Vector<FArrayBox*>curr            (numBlocks);
        Vector<FArrayBox*>magFluxFab      (numBlocks);
        Vector<FArrayBox*>currTemp        (numBlocks);
        Vector<FArrayBox*>hiRes           (numBlocks);
        Vector<FArrayBox*>initDisc        (numBlocks);
        Vector<FArrayBox*>oldPhysicalCoord(numBlocks);
        Vector<FArrayBox*>newPhysicalCoord(numBlocks);
        Vector<FArrayBox*>smoothDiag      (numBlocks);
        Vector<FArrayBox*>Jacobian        (numBlocks);
        Vector<FArrayBox*>weights         (numBlocks);
        Vector<FArrayBox*>smoothLagrangian(numBlocks);
        
        for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
          {
            //discretizes the boundary at high resolution
            hiRes           [iBlock] = new FArrayBox(hiResBox[iBlock],numComp   );
            curr            [iBlock] = new FArrayBox(box     [iBlock],numComp   );
            magFluxFab      [iBlock] = new FArrayBox(box     [iBlock],1         );
            currTemp        [iBlock] = new FArrayBox(box     [iBlock],numComp   );
            initDisc        [iBlock] = new FArrayBox(box     [iBlock],numComp   );
            oldPhysicalCoord[iBlock] = new FArrayBox(box     [iBlock],numComp   );
            newPhysicalCoord[iBlock] = new FArrayBox(box     [iBlock],numComp   );
            smoothDiag      [iBlock] = new FArrayBox(box     [iBlock],numComp   );
            Jacobian        [iBlock] = new FArrayBox(box     [iBlock],numJacComp);
            weights         [iBlock] = new FArrayBox(box     [iBlock],numWeights);
            smoothLagrangian[iBlock] = new FArrayBox(box     [iBlock],oneComp   );
          }
    
        // unit vector field perpendicular to magnetic field lines combined with (all) four first partial derivatives
        int num1Deriv      = 4 ;
        int num2Deriv      = 6 ;
        int numSmoothTerms = 14;
        int fluxComp       = 3 ;
        
        Vector<FArrayBox*>perpField  (numBlocks);
        Vector<FArrayBox*>BField     (numBlocks);
        
        Vector<FArrayBox*>mapPartials(numBlocks);
        Vector<FArrayBox*>BxBy       (numBlocks);
        Vector<FArrayBox*>perpXperpY (numBlocks);
        Vector<FArrayBox*>smoothness (numBlocks);
           
        for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
          {
            perpField  [iBlock] = new FArrayBox(box[iBlock],numComp   + num1Deriv           );
            BField     [iBlock] = new FArrayBox(box[iBlock],numComp   + num1Deriv + fluxComp);
            BxBy       [iBlock] = new FArrayBox(box[iBlock],num1Deriv                       );
            mapPartials[iBlock] = new FArrayBox(box[iBlock],num1Deriv + num2Deriv           );
            perpXperpY [iBlock] = new FArrayBox(box[iBlock],num1Deriv                       );
            smoothness [iBlock] = new FArrayBox(box[iBlock],numSmoothTerms                  );
          }
         
        // paths for filenames
        std::string geometryBaseName;
        std::string checkpointBaseName;
        std::string outputFileBaseName;
        
        // even spaced boundary nodes or more aligned
        bool evenSpacingCase;
        pp.get ("evenSpacingCase",evenSpacingCase);
        bool readFromCheckpoint;
        pp.get ("readFromCheckpoint",readFromCheckpoint);

        if (evenSpacingCase)
          {
            if (readFromCheckpoint)
              {
                pp.get("evenCheckpointBaseName",geometryBaseName  );
                pp.get("evenOutputBaseName"    ,outputFileBaseName);              
              }
            else
              {
                pp.get("evenGeometryBaseName",geometryBaseName  );
                pp.get("evenOutputBaseName"  ,outputFileBaseName);
              }
          }
        else
          {
            if (readFromCheckpoint)
              {
                pp.get("alignedCheckpointBaseName",geometryBaseName  );
                pp.get("alignedOutputBaseName"    ,outputFileBaseName);              
              }
            else
              {
                pp.get("alignedGeometryBaseName",geometryBaseName);
                pp.get("alignedOutputBaseName"  ,outputFileBaseName);
              }
          }

        char fileNameSuffix[1000];
        sprintf(fileNameSuffix,".%03d.2d.txt",baseRes + 1);

        std::string geometryFileName   = geometryBaseName   + fileNameSuffix;
        std::string outputFileName     = outputFileBaseName + fileNameSuffix;
       
        std::string fieldFileName;
        pp.get("fieldFileName", fieldFileName);

        Real sigma;
        pp.get("sigma",sigma);
      
        Real expo;
        pp.get("expo", expo);
        if (expo < 0.0)
          {
            MayDay::Abort("exponent < 0");
          }

        Vector<Vector<FArrayBox*> >tangent;
        tangent.resize(numBlocks);
        
        Vector<Vector<FArrayBox*> >arcLength;
        arcLength.resize(numBlocks);

        Vector<Vector<FArrayBox*> >targetSpacing;
        targetSpacing.resize(numBlocks);

        Vector<Vector<FArrayBox*> >targetSpacingTemp;
        targetSpacingTemp.resize(numBlocks);

        Vector<Vector<FArrayBox*> >hiResArcLength;
        hiResArcLength.resize(numBlocks);

        Vector<Vector<FArrayBox*> >fluxMatchedVal;
        fluxMatchedVal.resize(numBlocks);
                
        // ordering of edges is bottom,right,top,left
        Vector<int> tanDir(4); 
        tanDir[0] = xDir;
        tanDir[1] = yDir;
        tanDir[2] = xDir;
        tanDir[3] = yDir;
        
        // FArrayBoxes for tangent vector to the (hiRes representation of the) boundary
        for (int iBlock = beginBlock; iBlock <= endBlock; ++iBlock)
          {
            tangent          [iBlock].resize(4);
            arcLength        [iBlock].resize(4);
            targetSpacing    [iBlock].resize(4);
            targetSpacingTemp[iBlock].resize(4);
            fluxMatchedVal   [iBlock].resize(4);
            hiResArcLength   [iBlock].resize(4);
            
            for(int iEdge = 0;iEdge< 4; ++iEdge)
              {
                // ordering of edges: bottom,right,top,left
                tangent          [iBlock][iEdge] = new FArrayBox(edgeBoxes     [iBlock][iEdge],SpaceDim);
                arcLength        [iBlock][iEdge] = new FArrayBox(edgeBoxes     [iBlock][iEdge],oneComp );
                targetSpacing    [iBlock][iEdge] = new FArrayBox(edgeBoxes     [iBlock][iEdge],oneComp );
                targetSpacingTemp[iBlock][iEdge] = new FArrayBox(edgeBoxes     [iBlock][iEdge],oneComp );
                fluxMatchedVal   [iBlock][iEdge] = new FArrayBox(edgeBoxes     [iBlock][iEdge],oneComp );
                hiResArcLength   [iBlock][iEdge] = new FArrayBox(hiResEdgeBoxes[iBlock][iEdge],oneComp );
              }
          }
        
        Vector<int> keepWorking;
        keepWorking.resize(numBlocks);
        for (int iBlock = beginBlock; iBlock <= endBlock; ++iBlock)
          {
            keepWorking[iBlock] = 1;
            Vector<Real>evenSpacingVal  (4);
            
            for (int iter = 0; iter < maxIter_BrackbillSalzman && (keepWorking[iBlock] == 1);++iter)
              {
                 if (iter % checkPointInt == 0)
                  {
                    // always write without ghost cells
                    int numGhost = 0;
                    int numRad = 1 + numRow[iBlock] - 2*numGhost;

                    checkPoint(beginBlock        ,
                               endBlock          ,
                               numBlocks         ,
                               geometryFileName  ,
                               fieldFileName     ,
                               curr              ,
                               numRad            ,
                               hiCorner          ,
                               BField            ,
                               dx                ,
                               iter              ,
                               numGhost          ,
                               checkpointBaseName,
                               fileNameSuffix    );
                  }

                std::string blockName;
                getBlockName(blockName,
                             iBlock  );
                    
                // Construct the field data object
                FieldData field_data(geometryFileName,
                                     fieldFileName   ,
                                     blockName       );

                if (iter == 0)
                  { 
                    // hi Res mapped coordinates of the unit square
                    setPhysicalCoordinates(*hiRes       [iBlock],
                                           hiResHiCorner[iBlock],
                                           blockName            ,
                                           geometryFileName     ,
                                           fieldFileName        );
                    
                    // mapped coordinates of the unit square
                    setPhysicalCoordinates(*newPhysicalCoord[iBlock],
                                           hiCorner         [iBlock],
                                           blockName                ,
                                           geometryFileName         ,
                                           fieldFileName            );
                    
                    // initialize with output from FieldData interpolation
                    curr[iBlock]->copy(*newPhysicalCoord[iBlock]);
                    
                    // parameters that balance between variational principles
                    bool hasXPt;
                    IntVect xPt = findXiEtaXPt(hasXPt,
                                               iBlock);
                    
                    setWeights(*(weights[iBlock]),
                               nCell    [iBlock] ,
                               xiEtaWeight       ,
                               sigma             ,
                               expo              ,
                               dx       [iBlock] ,
                               hasXPt            ,
                               xPt               );

                    computeMapPartials(*(mapPartials[iBlock]),
                                       *(curr       [iBlock]),
                                       dx           [iBlock] ,
                                       writeInfo             );
                    
                    computeJacobianPartials(*(Jacobian   [iBlock]),
                                            *(mapPartials[iBlock]),
                                            writeInfo             );
                                                               
                    // use FieldData class for magnetic field and first partials
                    setBField(*(BField[iBlock]),
                              *(curr  [iBlock]),
                              dx      [iBlock] ,
                              geometryFileName ,
                              fieldFileName    ,
                              blockName        );
                                                
                    // recompute x and y derivatives
                    setdXdY(*(BxBy       [iBlock]),
                            *(BField     [iBlock]),
                            *(mapPartials[iBlock]),
                            dx           [iBlock]);

                    // diagonal term for the smooth operator 
                    computeSmJnDiagonalTerm(*(smoothDiag [iBlock]),
                                            dx           [iBlock] ,
                                            *(weights    [iBlock]),
                                            *(Jacobian   [iBlock]),
                                            *(mapPartials[iBlock]));
                    
                    computeBdSmJnDiagonalTerm(*(smoothDiag [iBlock]),
                                              dx           [iBlock] ,
                                              *(weights    [iBlock]),
                                              *(mapPartials[iBlock]),
                                              *(Jacobian   [iBlock]),
                                              tangent      [iBlock] ,
                                              tanDir               );
                    
                    // switch x for y and negate zeroth component
                    setPerpField(*(perpField[iBlock]),
                                 *(curr     [iBlock]),
                                 *(BField   [iBlock]),
                                 dx         [iBlock] ,
                                 geometryFileName    ,
                                 fieldFileName       ,
                                 blockName           ,
                                 writeInfo           );
                        
                    // recompute x and y derivatives
                    setdXdY(*(perpXperpY [iBlock]),
                            *(perpField  [iBlock]),
                            *(mapPartials[iBlock]),
                            dx           [iBlock]);
      
                    if (writePlot)
                      {
                        pout()<<"Write target and initial" <<endl;
                        char filename[1000];
                            
                        // dummy default argument
                        Vector<string> componentNames;
                            
                        // initial mapping 
                        sprintf(filename,"diagnostic/%03d.initial.init.%s.2d.hdf5",baseRes,blockName.c_str());
                        writeFABname(curr[iBlock],filename,componentNames,dx[iBlock][0] );
                            
                        // map derivatives
                        sprintf(filename,"diagnostic/%03d.mapPartials.init.%s.2d.hdf5",baseRes,blockName.c_str());
                        writeFABname(mapPartials[iBlock],filename,componentNames,dx[iBlock][0] );
                            
                        // map derivatives
                        sprintf(filename,"diagnostic/%03d.smoothness.init.%s.2d.hdf5",baseRes,blockName.c_str());
                        writeFABname(smoothness[iBlock],filename,componentNames,dx[iBlock][0] );
                            
                        // original physical coordinates of mapped unit square 
                        sprintf(filename,"diagnostic/%03d.oldPhysicalCoord.curr.%s.2d.hdf5",baseRes,blockName.c_str());
                        writeFABname(oldPhysicalCoord[iBlock],filename,componentNames,dx[iBlock][0] );
                            
                        // modified physical coordinates of mapped unit square 
                        sprintf(filename,"diagnostic/%03d.newPhysicalCoord.init.%s.2d.hdf5",baseRes,blockName.c_str());
                        writeFABname(newPhysicalCoord[iBlock],filename,componentNames,dx[iBlock][0] );
                            
                        // B field
                        sprintf(filename,"diagnostic/%03d.BField.init.%s.2d.hdf5",baseRes,blockName.c_str());
                        writeFABname(BField[iBlock],filename,componentNames,dx[iBlock][0] );

                        // x and y derivatives of B field
                        sprintf(filename,"diagnostic/%03d.BxBy.init.%s.2d.hdf5",baseRes,blockName.c_str());
                        writeFABname(BxBy[iBlock],filename,componentNames,dx[iBlock][0] );
                        
                        // perp field with all first partials in xi and eta
                        sprintf(filename,"diagnostic/%03d.perpField.init.%s.2d.hdf5",baseRes,blockName.c_str());
                        writeFABname(perpField[iBlock],filename,componentNames,dx[iBlock][0] );
                            
                        // x and y derivatives of perp field
                        sprintf(filename,"diagnostic/%03d.perpXperpY.init.%s.2d.hdf5",baseRes,blockName.c_str());
                        writeFABname(perpXperpY[iBlock],filename,componentNames,dx[iBlock][0] );
                            
                        // weights for parameterizing the various variational principles
                        sprintf(filename,"diagnostic/%03d.weights.init.%s.2d.hdf5",baseRes,blockName.c_str());
                        writeFABname(weights[iBlock],filename,componentNames,dx[iBlock][0] );
                      }
                  }

                Real maxNorm = LARGEREALVAL;
                bool minGprime = true;
                bool separatrixGrid = false;
                if (minGprime)
                  {
                    int startHiLo = 0;
                    int marchDir  = 0;
                    
                    Real tolerance = 1.0e-10;
                    
                    fillTargetFlux(*(magFluxFab[iBlock]),
                                   *(curr      [iBlock]),
                                   field_data          );
                    
                    int numSubSteps = 10*(nCell[iBlock][0]);
                   
                    orthoMarch(*(curr      [iBlock]),
                               *(BField    [iBlock]),
                               *(perpField [iBlock]),
                               *(magFluxFab[iBlock]),
                               dx          [iBlock] ,
                               geometryFileName     ,
                               fieldFileName        ,
                               blockName            ,
                               writeInfo            ,
                               startHiLo            ,
                               marchDir             ,
                               tolerance            ,
                               iBlock               ,
                               numSubSteps          );
                  }
                else if(separatrixGrid)
                  {
                    int startHiLo = 0;
                    int marchDir  = 0;
                    
                    Real tolerance = 1.0e-10;
                    
                    fillTargetFlux(*(magFluxFab[iBlock]),
                                   *(curr      [iBlock]),
                                   field_data          );
                    
                    int numSubSteps = 10*(nCell[iBlock][0]);
                    sepGrid(*(curr      [iBlock]),
                            *(BField    [iBlock]),
                            *(perpField [iBlock]),
                            *(magFluxFab[iBlock]),
                            dx          [iBlock] ,
                            geometryFileName     ,
                            fieldFileName        ,
                            blockName            ,
                            writeInfo            ,
                            startHiLo            ,
                            marchDir             ,
                            tolerance            ,
                            iBlock               ,
                            numSubSteps          );
                  }
                else
                  {
                    // variational grid generation from BrackBillSalzman (JCP 1982)
                    int numInitSmooth;
                    pp.get ("numInitSmooth",numInitSmooth);
                    
                    // in this code smoothness is not used; set to zero
                    pointJacobiNew(*(curr        [iBlock]),
                                   *(BField      [iBlock]),
                                   *(BxBy        [iBlock]),
                                   *(perpField   [iBlock]),
                                   *(perpXperpY  [iBlock]),
                                   *(mapPartials [iBlock]),
                                   *(smoothness  [iBlock]),
                                   *(smoothDiag  [iBlock]),
                                   *(Jacobian    [iBlock]),
                                   *(weights     [iBlock]),
                                   arcLength     [iBlock] ,
                                   targetSpacing [iBlock] ,
                                   maxNorm                ,
                                   *(hiRes       [iBlock]),
                                   tanDir                 ,
                                   tangent       [iBlock] ,
                                   hiResArcLength[iBlock] ,
                                   iter                   ,
                                   writeInfo              ,
                                   writePlot              ,
                                   dx            [iBlock] ,
                                   baseRes                ,
                                   geometryFileName       ,
                                   fieldFileName          ,
                                   blockName              ,
                                   xiEtaWeight            ,
                                   sigma                  ,
                                   expo                   ,
                                   refRatio               ,
                                   numInitSmooth          );
                  }
                
                // one round of relaxation for iBlock data
                pout()<<"finished pointJacobi for iter = "<<iter<<endl;
                
                // dummy default argument
                Vector<string> componentNames;
                
                
                bool keepWorkingBool = (stoppingCriterion(maxNorm                  ,
                                                          maxTol                   ,
                                                          iter + 1                 ,
                                                          maxIter_BrackbillSalzman));
                if (!keepWorkingBool)
                  {
                    keepWorking[iBlock] = -1;
                  }
                if (keepWorking[iBlock] != 1)
                  {
                    string prefix("diagnostic/");
                    std::string suffix(".2d.hdf5");
                    
                    char name[1000];
                    
                    // dummy default argument
                    Vector<string> componentNames;
                    
                    // final map
                    //xiEtaTranspose(a_curr);
                    sprintf(name,"diagnostic/%03d.curr.final.%s.2d.hdf5",baseRes,blockName.c_str());
                    writeFABname(curr[iBlock],name,componentNames,dx[iBlock][0]);
                    
                    // final map partials
                    sprintf(name,"diagnostic/%03d.mapPartials.final.%s.2d.hdf5",baseRes,blockName.c_str());
                    writeFABname(mapPartials[iBlock],name,componentNames,dx[iBlock][0]);
                    
                    // final smoothness
                    sprintf(name,"diagnostic/%03d.smoothness.final.%s.2d.hdf5",baseRes,blockName.c_str());
                    writeFABname(smoothness[iBlock],name,componentNames,dx[iBlock][0]);
                    
                    // final smoothDiag
                    sprintf(name,"diagnostic/%03d.smoothDiag.final.%s.2d.hdf5",baseRes,blockName.c_str());
                    writeFABname(smoothDiag[iBlock],name,componentNames,dx[iBlock][0]);
                    
                    // final Jacobian
                    sprintf(name,"diagnostic/%03d.Jacobian.final.%s.2d.hdf5",baseRes,blockName.c_str());
                    writeFABname(Jacobian[iBlock],name,componentNames,dx[iBlock][0]);
                    
                    // final smoothLagrangian
                    sprintf(name,"diagnostic/%03d.smoothLagrangian.final.%s.2d.hdf5",baseRes,blockName.c_str());
                    writeFABname(smoothLagrangian[iBlock],name,componentNames,dx[iBlock][0]);
                    
                    // BField
                    sprintf(name,"diagnostic/%03d.BField.final.%s.2d.hdf5",baseRes,blockName.c_str());
                    writeFABname(BField[iBlock],name,componentNames,dx[iBlock][0]); 
                    
                        // perpField
                    sprintf(name,"diagnostic/%03d.perpField.final.%s.2d.hdf5",baseRes,blockName.c_str());
                    writeFABname(perpField[iBlock],name,componentNames,dx[iBlock][0]); 
                    
                    // perpXperpY
                    sprintf(name,"diagnostic/%03d.perpXperpY.final.%s.2d.hdf5",baseRes,blockName.c_str());
                    writeFABname(perpXperpY[iBlock],name,componentNames,dx[iBlock][0]); 
                  }
              }
          }
        
        for (int iBlock = 0; iBlock < beginBlock; ++iBlock)
          {
            
            std::string blockName;
            getBlockName(blockName,
                         iBlock  );
            
            // Construct the field data object
            FieldData field_data(geometryFileName,
                                 fieldFileName   ,
                                 blockName       );
            
            
            // mapped coordinates of the unit square
            setPhysicalCoordinates(*(curr  [iBlock]),
                                   hiCorner[iBlock] ,
                                   blockName        ,
                                   geometryFileName ,
                                   fieldFileName    );
            
            
            // use FieldData class for magnetic field and first partials
            setBField(*(BField[iBlock]),
                      *(curr  [iBlock]),
                      dx       [iBlock],
                      geometryFileName ,
                      fieldFileName    ,
                      blockName        );
          }
        
        for (int iBlock = endBlock + 1; iBlock < numBlocks; ++iBlock)
          {
            std::string blockName;
            getBlockName(blockName,
                         iBlock   );
            
            // Construct the field data object
            FieldData field_data(geometryFileName,
                                 fieldFileName   ,
                                 blockName       );
            
                
            // mapped coordinates of the unit square
            setPhysicalCoordinates(*(curr  [iBlock]),
                                   hiCorner[iBlock] ,
                                   blockName        ,
                                   geometryFileName ,
                                   fieldFileName    );
            
            // use FieldData class for magnetic field and first partials
            setBField(*(BField[iBlock]),
                      *(curr  [iBlock]),
                      dx      [iBlock] ,
                      geometryFileName ,
                      fieldFileName    ,
                      blockName        );
          }
        
        // clean start for mapping file output function    
        clearExtendedMappingFile(outputFileName.c_str());

        for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
          {
            std::string blockName;
            getBlockName(blockName,
                         iBlock   );
           
            // write unghosted mapping file
            outputMapping(*(curr   [iBlock])    ,
                          *(BField[iBlock])     ,
                          IntVect::Zero         ,
                          blockName             ,
                          outputFileName.c_str());
          }
        
        // write coarser versions of curr and BField
        
        for (int iBlock = beginBlock; iBlock <= endBlock; ++iBlock)
          {
            Box validBox = curr[iBlock]->box();
            Box coarseBox = validBox; 
            int coarseningFactor = 2;
            for (int numberCells = log2MaxRes - 1; numberCells > 2; -- numberCells)
              {
                // make coarse curr
                IntVect coarseBigEnd = validBox.bigEnd();
                coarseBigEnd = coarseBigEnd/2;
                coarseBox.define(IntVect::Zero,coarseBigEnd);
                FArrayBox coarseCurr(coarseBox,SpaceDim);
    
                std::string blockName;
                getBlockName(blockName,
                             iBlock   );

                int currRes = pow(2,numberCells);
                Real dx = 1.0/((Real)currRes);
               
                for (BoxIterator bit(coarseBox);bit.ok();++bit) 
                  {
                    IntVect coarseIv = bit();
                    IntVect fineIv   = coarseIv*coarseningFactor;
                    coarseCurr(coarseIv,xDir) = (*curr[iBlock])(fineIv,xDir);
                    coarseCurr(coarseIv,yDir) = (*curr[iBlock])(fineIv,yDir);               
                  }
                string prefix("diagnostic/");
                std::string suffix(".2d.hdf5");
                
                char name[1000];
                
                // dummy default argument
                Vector<string> componentNames;
                int coarseRes = pow(2,numberCells);
                
                sprintf(name,"diagnostic/%03d.curr.final.%s.2d.hdf5",coarseRes,blockName.c_str());
                writeFABname(&coarseCurr,name,componentNames,dx);
                
                // coarsen valid box
                validBox = coarseBox;
                coarseningFactor *= 2;
              }
          }
            //delete pointers (last in first out)
        for (int iBlock = beginBlock; iBlock <= endBlock; ++iBlock)
          {
            delete perpField  [iBlock];
            delete BField     [iBlock];
            delete BxBy       [iBlock];
            delete mapPartials[iBlock];
            delete perpXperpY [iBlock];
            delete smoothness [iBlock];
            
            delete curr            [iBlock];
            delete initDisc        [iBlock];
            delete oldPhysicalCoord[iBlock];
            delete newPhysicalCoord[iBlock];
            delete smoothDiag      [iBlock];
            delete Jacobian        [iBlock];
            delete weights         [iBlock];
            delete smoothLagrangian[iBlock];
          }
       
        // write an empty file to signal python script that bs.cpp is finished.
        FILE* fd;
        fd = fopen ("done","w");
        if (fd == NULL) 
          {
            perror ("Error opening file");
          }
        fclose(fd);
      }
  }// End scoping trick 
  
#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
  
  return 0;
}


