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



void setBField(FArrayBox         & a_BField              ,
               const FArrayBox   & a_newPhysicalCoord    ,
               const RealVect    & a_dx                  ,
               const std::string & a_geometryFileNameNewBd,
               const std::string & a_fieldFileNameNewBd   , 
               const std::string & a_blockName            )
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
  FieldData field_data(a_geometryFileNameNewBd,
                       a_fieldFileNameNewBd   ,
                       a_blockName            );
  
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



void incrementBdSpacing(Vector<RealVect>       & a_movedVec         ,
                        Vector<Real>           & a_moved_xi_eta     ,
                        Vector<Real>           & a_arcLengthVec     ,
                        const Vector<Real>     & a_targetSpacingVec ,
                        const Vector<Real>     & a_incrementVec     ,
                        const Vector<RealVect> & a_hiResVec         ,
                        const Vector<Real>     & a_hiResArcLengthVec,
                        const Box              & a_validBox         ,
                        const bool             & a_writeInfo        ,
                        const int              & a_refRatio         )
{
  int xDir = 0;
  int yDir = 1;
  
  for (int iVal = 1; iVal < a_arcLengthVec.size(); ++iVal)
    {
      //a_arcLengthVec[iVal] += (a_incrementVec[iVal] - a_incrementVec[iVal -1]);
      //Real gamma = 0.001;
      //gamma *= .01;
      Real gamma = 1.0;  
      a_arcLengthVec[iVal]  = gamma*a_targetSpacingVec[iVal] + (1.0 - gamma)*a_arcLengthVec[iVal]; 
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
          a_movedVec    [currNewCell] = (1.0 - frac)*a_hiResVec[iOldCell - 1] + frac*a_hiResVec[iOldCell];
          a_moved_xi_eta[currNewCell] = (1.0 - frac)*          (iOldCell - 1) + frac*           iOldCell;

          
          CH_assert(a_refRatio > 0);
          a_moved_xi_eta[currNewCell]/= numFaceIv;
          
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

  a_moved_xi_eta[0] = 0;
  a_moved_xi_eta[a_movedVec.size() - 1] = 1.0;
  
    // debugging output
  Vector<Real> newSpacing(a_movedVec.size());

    for (int iCell = 1; iCell < a_movedVec.size(); ++iCell)
      {
        newSpacing[iCell]  = (a_movedVec[iCell][0] - a_movedVec[iCell - 1][0])*(a_movedVec[iCell][0] - a_movedVec[iCell - 1][0]);
        newSpacing[iCell] += (a_movedVec[iCell][1] - a_movedVec[iCell - 1][1])*(a_movedVec[iCell][1] - a_movedVec[iCell - 1][1]);
        newSpacing[iCell] = sqrt(newSpacing[iCell]);
      }
    
    newSpacing[0] = 0;
#if 0
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
                 const std::string        & a_blockName     )
{
  Box box = a_curr.box     ();
  IntVect lo = box.smallEnd();
  IntVect hi = box.bigEnd  ();
  
  IntVect hiResLo = lo;
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

      //      pout()<<"iEdge = "<<iEdge<<endl;
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
      Vector<Real>  moved_xi_eta    (numFaceIv);
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
                         moved_xi_eta     ,
                         arcLengthVec     ,
                         targetSpacingVec ,
                         incrementVec     ,
                         hiResVec         ,
                         hiResArcLengthVec,
                         box              ,
                         a_writeInfo      ,
                         a_refRatio       );
      
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
              
              // new values for (xi,eta) on this edge
              a_xi_eta(iv,tanDir) = moved_xi_eta[iVec];
              a_xi_eta(iv,norDir) = hilo;

              iVec += 1;
            }
        }
    }
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

void setPhysicalCoordinates (FArrayBox         & a_oldPhysicalCoord      ,
                             const IntVect     & a_hiCorner              ,
                             const std::string & a_blockName             , 
                             const std::string & a_geometryFileNameNewBd ,
                             const std::string & a_fieldFileNameNewBd    )
{
  // Construct the field data object
  FieldData field_data(a_geometryFileNameNewBd,
                       a_fieldFileNameNewBd   ,
                       a_blockName            );

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
  int xiLoVal = loCorner[yDir];
  int xiHiVal = hiCorner[xDir];
  
  int etaLoVal = loCorner[xDir];
  int etaHiVal = hiCorner[yDir];

  Box interiorBox = domain;
  interiorBox.grow(-1);

  // iterate over interior 
  for (BoxIterator bit(interiorBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      // xi stencil values
      IntVect loXi(iv[xDir],etaLoVal);
      IntVect hiXi(iv[xDir],etaHiVal);
      Real xi0  = a_initXiEta(loXi ,comp0);
      Real xi1  = a_initXiEta(hiXi ,comp0);

      // eta stencil values
      IntVect loEta(xiLoVal,iv[yDir]);
      IntVect hiEta(xiHiVal,iv[yDir]);
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

void moveBdPoints(Vector<RealVect>       & a_movedVec  ,
                  Vector<Real>           & a_xiEtaBCVal,
                  const Vector<RealVect> & a_givenVec  ,
                  const Vector<Real>     & a_givenXiEta,
                  const Real             & a_shiftDist ,
                  const Box              & a_validBox  ,
                  const bool             & a_moveAll   ,
                  const bool             & a_writeInfo )
{
  int xDir = 0;
  int yDir = 1;

  // end points of the edge
  RealVect xPoint      = a_givenVec[0];
  RealVect otherCorner = a_givenVec[a_givenVec.size() -1];

  IntVect smallEnd = a_validBox.smallEnd();
  IntVect bigEnd   = a_validBox.bigEnd  ();
  
  int numRow = bigEnd[yDir] - smallEnd[yDir] + 1;
  int numCol = bigEnd[xDir] - smallEnd[xDir] + 1;
  RealVect dx(1.0/(numCol - 1),1.0/(numRow - 1));

  int numFaceIv = a_givenVec.size();
  Vector<Real> increments  (numFaceIv);
  Vector<Real> distToX     (numFaceIv);
  Vector<Real> arcLengthToX(numFaceIv);
  Vector<Real> xiEtaBcVal  (numFaceIv);

  increments  [0] = 0.0;
  distToX     [0] = 0.0;
  arcLengthToX[0] = 0.0;

  int lastShiftedCell = 0;
  for (int iCell = 1; iCell < numFaceIv; ++iCell)
    {
      RealVect curr = a_givenVec[iCell    ];
      RealVect prev = a_givenVec[iCell - 1];
      
      Real segLength = 0.0;
      segLength += (prev[xDir] - curr[xDir])*(prev[xDir] - curr[xDir]);
      segLength += (prev[yDir] - curr[yDir])*(prev[yDir] - curr[yDir]);
      segLength  = sqrt(segLength);
      
      Real dist = 0;
      dist  += (a_givenVec[0][xDir] - curr[xDir])*(a_givenVec[0][xDir] - curr[xDir]);
      dist  += (a_givenVec[0][yDir] - curr[yDir])*(a_givenVec[0][yDir] - curr[yDir]);
      dist   = sqrt(dist);

      increments  [iCell] = segLength;
      distToX     [iCell] = dist     ;
      arcLengthToX[iCell] = arcLengthToX[iCell - 1] + segLength;

      if ((distToX[iCell] < a_shiftDist) || a_moveAll)
        {
          lastShiftedCell = iCell;
        }
    }
 
  // can only move up to the last cell
  if (lastShiftedCell == numFaceIv - 1)
    {
      lastShiftedCell -= 1;
    }
  
  Real evenSpacingValue = arcLengthToX[1 + lastShiftedCell]/(1 + lastShiftedCell);
  Real disc             = evenSpacingValue;
 
  int currNewCell = 1;
  
  for(int iOldCell  = 1; iOldCell< numFaceIv; iOldCell++)
    {
      if ((iOldCell <= lastShiftedCell))
        {
          while(disc < increments[iOldCell])
            {
              // relative discrepancy between hi and lo cells
              Real frac = disc/increments[iOldCell];
              
              // assign new node
              a_movedVec[currNewCell] = (1.0 - frac)*a_givenVec[iOldCell - 1] + frac*a_givenVec[iOldCell];
              
              // invert the new point
              a_xiEtaBCVal[currNewCell] = (1.0 - frac)*a_givenXiEta[iOldCell - 1] + frac*a_givenXiEta[iOldCell];

              // next new node
              currNewCell += 1;
              disc += evenSpacingValue;
            }

          // measure where even spacing multiples land, relative to current arc length increment 
          disc -= increments[iOldCell];
        }
      else
        {
          // iOldCell and currNewCell agree when ioldCell is beyond the lastShiftedCell
          a_movedVec  [iOldCell] = a_givenVec  [iOldCell];
          a_xiEtaBCVal[iOldCell] = a_givenXiEta[iOldCell];
        }
    }
  
    // the x-Point location never changes.
    a_movedVec  [0] = a_givenVec  [0];
    a_xiEtaBCVal[0] = a_givenXiEta[0];

    // last node never moves
    a_movedVec  [numFaceIv - 1] = a_givenVec  [numFaceIv - 1]; 
    a_xiEtaBCVal[numFaceIv - 1] = a_givenXiEta[numFaceIv - 1];

    // debugging output
    Vector<Real> newSpacing(numFaceIv);
    Vector<Real> oldSpacing(numFaceIv);
    for (int iCell = 1; iCell < numFaceIv; ++iCell)
      {
        newSpacing[iCell]  = (a_movedVec[iCell][0] - a_movedVec[iCell - 1][0])*(a_movedVec[iCell][0] - a_movedVec[iCell - 1][0]);
        newSpacing[iCell] += (a_movedVec[iCell][1] - a_movedVec[iCell - 1][1])*(a_movedVec[iCell][1] - a_movedVec[iCell - 1][1]);
        newSpacing[iCell] = sqrt(newSpacing[iCell]);

        oldSpacing[iCell]  = (a_givenVec[iCell][0] - a_givenVec[iCell - 1][0])*(a_givenVec[iCell][0] - a_givenVec[iCell - 1][0]);
        oldSpacing[iCell] += (a_givenVec[iCell][1] - a_givenVec[iCell - 1][1])*(a_givenVec[iCell][1] - a_givenVec[iCell - 1][1]);
        oldSpacing[iCell] = sqrt(oldSpacing[iCell]);
      }
    
    newSpacing[0] = 0;
    oldSpacing[0] = 0;
    if (a_writeInfo)
      {
        pout()<<"evenSpacingValue = "<<evenSpacingValue<<endl;
        pout()<<endl;
        pout()<<"lastShiftedCell = "<<lastShiftedCell<<endl;
        pout()<<endl;
        pout()<<"newSpacing[1 ] = "<<newSpacing[1]<<endl;
        pout()<<"newSpacing[2 ] = "<<newSpacing[2]<<endl;
        pout()<<"newSpacing[3 ] = "<<newSpacing[3]<<endl;
        pout()<<"newSpacing[4 ] = "<<newSpacing[4]<<endl;
        pout()<<"newSpacing[5 ] = "<<newSpacing[5]<<endl;
        pout()<<"newSpacing[6 ] = "<<newSpacing[6]<<endl;
        pout()<<"newSpacing[7 ] = "<<newSpacing[7]<<endl;
        pout()<<"newSpacing[8 ] = "<<newSpacing[8]<<endl;
        pout()<<"newSpacing[9 ] = "<<newSpacing[6]<<endl;
        pout()<<"newSpacing[10] = "<<newSpacing[7]<<endl;
        pout()<<"newSpacing[11] = "<<newSpacing[8]<<endl;
        pout()<<endl;
        pout()<<"oldSpacing[1 ] = "<<oldSpacing[1]<<endl;
        pout()<<"oldSpacing[2 ] = "<<oldSpacing[2]<<endl;
        pout()<<"oldSpacing[3 ] = "<<oldSpacing[3]<<endl;
        pout()<<"oldSpacing[4 ] = "<<oldSpacing[4]<<endl;
        pout()<<"oldSpacing[5 ] = "<<oldSpacing[5]<<endl;
        pout()<<"oldSpacing[6 ] = "<<oldSpacing[6]<<endl;
        pout()<<"oldSpacing[7 ] = "<<oldSpacing[7]<<endl;
        pout()<<"oldSpacing[8 ] = "<<oldSpacing[8]<<endl;
        pout()<<"oldSpacing[9 ] = "<<oldSpacing[6]<<endl;
        pout()<<"oldSpacing[10] = "<<oldSpacing[7]<<endl;
        pout()<<"oldSpacing[11] = "<<oldSpacing[8]<<endl;
        pout()<<endl;
        pout()<<"shift 0  = "<<a_givenVec[0 ] - a_movedVec[0 ]<<endl;
        pout()<<"shift 1  = "<<a_givenVec[1 ] - a_movedVec[1 ]<<endl;
        pout()<<"shift 2  = "<<a_givenVec[2 ] - a_movedVec[2 ]<<endl;
        pout()<<"shift 3  = "<<a_givenVec[3 ] - a_movedVec[3 ]<<endl;
        pout()<<"shift 4  = "<<a_givenVec[4 ] - a_movedVec[4 ]<<endl;
        pout()<<"shift 5  = "<<a_givenVec[5 ] - a_movedVec[5 ]<<endl;
        pout()<<"shift 6  = "<<a_givenVec[6 ] - a_movedVec[6 ]<<endl;
        pout()<<"shift 7  = "<<a_givenVec[7 ] - a_movedVec[7 ]<<endl;
        pout()<<"shift 8  = "<<a_givenVec[8 ] - a_movedVec[8 ]<<endl;
        pout()<<"shift 9  = "<<a_givenVec[9 ] - a_movedVec[9 ]<<endl;
        pout()<<"shift 10 = "<<a_givenVec[10] - a_movedVec[10]<<endl;
        pout()<<"shift 11 = "<<a_givenVec[11] - a_movedVec[11]<<endl;
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



void computeFluxMatchedVal(FArrayBox                   & a_fluxMatchedVal       ,
                           const FArrayBox             & a_evenSpaceCurr        ,
                           const Vector<FArrayBox *>   & a_hiResArcLength       ,
                           const FArrayBox             & a_hiRes                ,
                           const Vector<Real>          & a_evenSpacingVal       ,
                           const int                   & a_iEdge                ,
                           const int                   & a_iBlock               ,
                           const RealVect              & a_hiResDx              ,
                           const int                   & a_refRatio             ,
                           const std::string           & a_geometryFileNameNewBd,
                           const std::string           & a_fieldFileNameNewBd   ,
                           const Vector<Vector<Real> > & a_fluxValues           )
{
  Real totalLength = 0.0;
  int zeroComp = 0;
  int xComp = 0;
  int yComp = 1;
  
  std::string blockName;
  getBlockName(blockName,
               a_iBlock );

  // Construct the field data object
  FieldData field_data(a_geometryFileNameNewBd,
                       a_fieldFileNameNewBd   ,
                       blockName              );
  
  int loHi = 1;
  if ((a_iEdge == 0) || a_iEdge == 1)
    {
      loHi = 0;
    }
  
  int tanDir = 1;
  if ((a_iEdge == 0) || (a_iEdge == 2))
    {
      tanDir = 0;
    }
  int norDir = (tanDir + 1)%2;
  
  const FArrayBox& hiResArcLengthNear = *(a_hiResArcLength[a_iEdge]);
  Box nearEdgeBox = hiResArcLengthNear.box();
    
  // hi res number of cells
  int numHiResCells = hiResArcLengthNear.box().size()[tanDir] - 1;
  
  //index in normal direction
  int normalIndex = hiResArcLengthNear.box().smallEnd()[norDir];

  // int normalIndexCheck = hiResArcLengthNear.box().bigEnd()[norDir];
    
  int currHiRes = 0;
  Real currArcLength = 0.0;
  Box box = a_fluxMatchedVal.box();
  IntVect smallEnd = box.smallEnd();
  IntVect bigEnd   = box.bigEnd  ();
  
  if (a_iEdge == 0 || a_iEdge == 2) 
    {
      for (BoxIterator bit(box);bit.ok();++bit) 
        {
          IntVect iv = bit();
          if (iv != smallEnd)
            {
              if (iv != bigEnd)
                {
                  Real farFlux = (a_fluxValues[a_iBlock])[iv[tanDir]];
                 
                  // iterate through fine data holders to find bracketing fluxes
                  for (int hiResCounter = currHiRes; hiResCounter < numHiResCells; ++hiResCounter) 
                    {
                      // condition for exiting the loop
                      bool foundBracketing = false;
                     
                      // lo member of bracketing ivs
                      IntVect hiResIvLo;
                      hiResIvLo[tanDir] = hiResCounter;
                      hiResIvLo[norDir] = loHi*normalIndex;
                     
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
                                  MayDay::Abort("not bracketed");
                                }

                              // double check in case the process hangs
                              counter += 1;
                              if (counter == 1000)
                                {
                                  pout()<<"bisection failed"<<endl;
                                  MayDay::Abort("bisection failed");
                                }
                            }
                         
                          Real dist2Hi  = sqrt((newPt[0] - hiResPtHi[0])*(newPt[0] - hiResPtHi[0]) + (newPt[1] - hiResPtHi[1])*(newPt[1] - hiResPtHi[1])); 
                          Real dist2Lo  = sqrt((newPt[0] - hiResPtLo[0])*(newPt[0] - hiResPtLo[0]) + (newPt[1] - hiResPtLo[1])*(newPt[1] - hiResPtLo[1]));
                                                                                                              
                          Real inc      = currArcLength + dist2Lo;
                          a_fluxMatchedVal(iv,zeroComp) = inc;
                          totalLength                  += inc;
                          currArcLength = dist2Hi;

                          //pout()<<currArcLength<<endl;
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
                      hiResIvHi[norDir] = loHi*normalIndex;
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
    
  //pout()<<"totalLength = "<<totalLength <<endl;
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
                             const Real      & a_totalArcLength          ,
                             const bool      & a_moveMcoreMcsol          )
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
                   a_targetSpacing(iv,0)       = a_fluxMatchedVal(iv,0);
                   a_totalTargetSpacingLength += a_fluxMatchedVal(iv,0);
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
                       MayDay::Abort("sigma must be positive");
                     }
                   else
                     {
                       weight    = exp(-pow(   dist/a_sigma,a_expo));
                       preWeight = exp(-pow(preDist/a_sigma,a_expo));
                     }
                   
                   // check the weight
                   //pout()<<"weight = "<<weight<<endl;

                   a_targetSpacing(iv,0)       = weight   *cumulativeInc(iv   ,evenComp) + (1.0 - weight   )*cumulativeInc(iv   ,fluxComp);
                   a_targetSpacing(iv,0)      -= preWeight*cumulativeInc(preIv,evenComp) + (1.0 - preWeight)*cumulativeInc(preIv,fluxComp);
                   a_totalTargetSpacingLength += a_targetSpacing(iv,0);
                 }
             }
           else
             {
               if ((a_iBlock < 8) || a_moveMcoreMcsol)
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

void checkCornerFlux(FArrayBox   &a_curr                 ,
                     std::string &a_blockName            ,
                     std::string &a_geometryFileNameNewBd,
                     std::string &a_fieldFileNameNewBd   )
{
  int xiDir = 0;
  Box box = a_curr.box();
  int numCell = box.size()[xiDir] - 1;
  IntVect ll = box.smallEnd();
  IntVect ur = box.bigEnd()  ; 
  IntVect lr(ll[0],ur[1])    ;
  IntVect ul(ur[0],ll[1])    ;
  
  // Construct the field data object
  FieldData field_data(a_geometryFileNameNewBd,
                       a_fieldFileNameNewBd   ,
                       a_blockName            );
  
  RealVect llr(a_curr(ll,0),a_curr(ll,1));
  RealVect urr(a_curr(ur,0),a_curr(ur,1));
  RealVect lrr(a_curr(lr,0),a_curr(lr,1));
  RealVect ulr(a_curr(ul,0),a_curr(ul,1));

  Real llFlux = field_data.getMagneticFluxFromDCT(llr);
  Real urFlux = field_data.getMagneticFluxFromDCT(urr);
  Real lrFlux = field_data.getMagneticFluxFromDCT(lrr);
  Real ulFlux = field_data.getMagneticFluxFromDCT(ulr);

  pout()<<"llFlux = "<<llFlux<<endl;
  pout()<<"lrFlux = "<<lrFlux<<endl;
  pout()<<"ulFlux = "<<ulFlux<<endl;
  pout()<<"urFlux = "<<urFlux<<endl;

  pout()<<"lowerFluxDiff = "<<llFlux - lrFlux<<endl;
  pout()<<"upperFluxDiff = "<<ulFlux - urFlux<<endl;

   for (int iCell = 0; iCell < numCell; ++iCell)
    {
      
      // matching iv
      IntVect hiIv = lr + iCell*BASISV(xiDir);
      IntVect loIv = ll + iCell*BASISV(xiDir);
      
      // mapped values that should have same flux (away from x point)
      RealVect hiSide(a_curr(hiIv,0),a_curr(hiIv,1));
      RealVect loSide(a_curr(loIv,0),a_curr(loIv,1));
      
      // magetic flux
      Real hiFlux =  field_data.getMagneticFluxFromDCT(hiSide);
      Real loFlux =  field_data.getMagneticFluxFromDCT(loSide);

      
      pout()<<"hiIv     = "<<hiIv<<endl;
      pout()<<"loIv     = "<<loIv<<endl;
      pout()<<endl;
      pout()<<"hiFlux   = "<<hiFlux<<endl;
      pout()<<"loFlux   = "<<loFlux<<endl;
      pout()<<endl;
      pout()<<"FluxDiff = "<<hiFlux - loFlux<<endl;
      
    }
  
}

void straightenEtaBd(Vector <FArrayBox* > & a_curr                 ,
                     Vector <FArrayBox* > & a_BField               ,
                     const std::string    & a_geometryFileNameNewBd,
                     const std::string    & a_fieldFileNameNewBd   ,
                     const int            & a_numStraighten        ,
                     const RealVect       & a_dx                   ,
                     const int            & a_beginBlockInit       ,
                     const int            & a_endBlockInit         )
{
  int xiDir  = 0;
  int etaDir = 1;
  
  int xDir = 0;
  int yDir = 1;
  
  // initialize ghost cells to 0 and interior to map
  for(int iBlock = a_beginBlockInit; iBlock <= a_endBlockInit; ++iBlock)
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
      FieldData field_data(a_geometryFileNameNewBd,
                           a_fieldFileNameNewBd   ,
                           blockName              );
      
      FArrayBox& currFab =  *a_curr[iBlock];
      setPhysicalCoordinates(currFab                ,
                             hiCorner               ,
                             blockName              ,
                             a_geometryFileNameNewBd,
                             a_fieldFileNameNewBd   );
      
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
                       
                      if ( (Abs (Disc1) > 1.0e-07) || (Abs (Disc2) > 1.0e-07) ||(Abs (Eisc1) > 1.0e-07) || (Abs (Eisc2) > 1.0e-07) )
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
                       
                      if ( (Abs (Disc1) > 1.0e-07) || (Abs (Disc2) > 1.0e-07) ||(Abs (Eisc1) > 1.0e-07) || (Abs (Eisc2) > 1.0e-07) )
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

void getFluxValueVector(Vector<Vector<Real> >     & a_fluxValues           ,
                        const Vector<FArrayBox*>  & a_evenSpacedCurr       ,
                        const Vector<Vector<Box> >& a_edgeBoxes            ,
                        const Vector <int>        & a_numCol               ,
                        const std::string         & a_geometryFileNameNewBd,
                        const std::string         & a_fieldFileNameNewBd   ,
                        const int                 & a_beginBlockInit       ,
                        const int                 & a_endBlockInit         )
{ 
  // macros for indexing
  int tanDir = 0;
  int xComp  = 0;
  int yComp  = 1;

  // number of blocks determines the size of vectors
  int numBlocks = a_evenSpacedCurr.size();

  // need to ask for magnetic flux from these blocks
  FieldData field_data9(a_geometryFileNameNewBd,
                        a_fieldFileNameNewBd   ,
                       "mcsol"                 );
  
  FieldData field_data8(a_geometryFileNameNewBd,
                       a_fieldFileNameNewBd    ,
                       "mcore"                 );
  
  FieldData field_data7(a_geometryFileNameNewBd,
                       a_fieldFileNameNewBd    ,
                       "rpf"                   );

  // resize the flux value vector
  a_fluxValues.resize(numBlocks);
  for (int iBlock = a_beginBlockInit; iBlock <= a_endBlockInit; ++iBlock)
    {
      a_fluxValues[iBlock].resize(a_numCol[iBlock] + 1);
    }
  
  // iteration boxes
  Box topEdge9    = a_edgeBoxes[9][2];
  Box topEdge8    = a_edgeBoxes[8][2];
  Box topEdge7    = a_edgeBoxes[7][2];

  Box bottomEdge9 = a_edgeBoxes[9][0];
  Box bottomEdge8 = a_edgeBoxes[8][0];
  Box bottomEdge7 = a_edgeBoxes[7][0];

  // block 9 has values corresponding to evenly spaced points on one (arbitrary, say the top) eta edge
  FArrayBox& mcsolFab  = *(a_evenSpacedCurr[9]);
  
  // block 8 has values corresponding to evenly spaced points on one (let's match block 9, so the top) eta edge
  FArrayBox& mcoreFab  = *(a_evenSpacedCurr[8]);
  
  // block 7 has values corresponding to evenly spaced points on the bottom   
  FArrayBox& rpfFab    = *(a_evenSpacedCurr[7]);
  
  for (BoxIterator bit(topEdge9);bit.ok();++bit) 
    {
      IntVect iv = bit();
      
      Real x9 = mcsolFab(iv,xComp);
      Real y9 = mcsolFab(iv,yComp);
      RealVect pt9(x9,y9);
             
      // get the flux at this node
      Real fluxValue9 =  field_data9.getMagneticFluxFromDCT(pt9);

      // block 9 has values corresponding to evenly spaced points on one (arbitrary, say the top) eta edge
      (a_fluxValues[9])[iv[tanDir]] = fluxValue9;
      
      // block 5 has values that match block 9
      // block 4 has values that match block 9
      // block 3 has values that match block 9
      // block 2 has values that match block 9
   
      (a_fluxValues[5])[iv[tanDir]] = fluxValue9;
      (a_fluxValues[4])[iv[tanDir]] = fluxValue9;
      (a_fluxValues[3])[iv[tanDir]] = fluxValue9;
      (a_fluxValues[2])[iv[tanDir]] = fluxValue9; 

      Real x8 = mcoreFab(iv,xComp);
      Real y8 = mcoreFab(iv,yComp);
      RealVect pt8(x8,y8);

       // get the flux at this node
      Real fluxValue8 =  field_data8.getMagneticFluxFromDCT(pt8);

      // block 8 has values corresponding to evenly spaced points on one (let's match block 9, so the top) eta edge
      (a_fluxValues[8])[iv[tanDir]] = fluxValue8;
      
      // block 1 has values that match block 8
      // block 0 has values that match block 8
      (a_fluxValues[1])[iv[tanDir]] = fluxValue8; 
      (a_fluxValues[0])[iv[tanDir]] = fluxValue8; 
    }
  
  for (BoxIterator bit(bottomEdge7);bit.ok();++bit) 
    {
      IntVect iv = bit();
      
      Real x7 = rpfFab(iv,xComp);
      Real y7 = rpfFab(iv,yComp);
      RealVect pt7(x7,y7);
      
      // get the flux at this node
       Real fluxValue7 =  field_data7.getMagneticFluxFromDCT(pt7);

       // block 7 has values corresponding to evenly spaced points on the bottom   
       (a_fluxValues[7])[iv[tanDir]] = fluxValue7;
     
       // block 6 has valuess that match block 7
       (a_fluxValues[6])[iv[tanDir]] = fluxValue7;
    }
}

void copyBoundary(FArrayBox& a_newCurr,
                  FArrayBox& a_oldCurr)
{
  Box interiorBox = a_newCurr.box();
  interiorBox.grow(-1);
  for (BoxIterator bit(a_newCurr.box());bit.ok();++bit) 
    {
      IntVect iv = bit();
      
      if(!interiorBox.contains(iv))
        {
          for (int iComp = 0; iComp < SpaceDim; ++iComp)
            {
              a_newCurr(iv,iComp) = a_oldCurr(iv,iComp);
            }
          
        }
    }
  
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
      
      //pout()<< "diffMid = "<<diffMid<<endl;
      if (Abs(diffMid) < a_tolerance)
        {
          findingRoot = false;  
          newPt = midPt;
        }
    }
  return newPt;
}

void alignBdPt(FArrayBox       & a_hiRes      ,
               const IntVect   & a_iv         ,
               const FieldData & a_field_data ,
               const int       & a_numSteps   ,
               const Real      & a_dist       , 
               const Real      & a_edgeFlux   ,
               const int       & a_blockNumber,
               const Real      & a_tolerance  )
{
  int xDir = 0;
  int yDir = 1;

  RealVect fieldDir  = LARGEREALVAL*RealVect::Unit;
  RealVect normalDir = LARGEREALVAL*RealVect::Unit;

  RealVect pt(a_hiRes(a_iv,xDir),a_hiRes(a_iv,yDir));
  RealVect basePt = pt;
  RealVect hiPt   = basePt;
  RealVect loPt   = basePt;
  RealVect bracketPt = LARGEREALVAL*RealVect::Unit;

  bool success = false;
  for (int iStep = 0; iStep < a_numSteps; ++iStep)
    {
      fieldDir = a_field_data.getFieldUnitVectorFromDCT(basePt);
      normalDir[0] = -fieldDir[1];
      normalDir[1] =  fieldDir[0];
      
      hiPt = basePt + normalDir*a_dist;
      loPt = basePt - normalDir*a_dist;
      
      Real hiFlux = a_field_data.getMagneticFluxFromDCT(hiPt  );
      Real loFlux = a_field_data.getMagneticFluxFromDCT(loPt  );
      
      if ((hiFlux - a_edgeFlux)*(loFlux - a_edgeFlux) < 0)
         {
           success = true;
           break;
         }

      Real hiDiff   = a_edgeFlux - hiFlux;
      Real loDiff   = a_edgeFlux - loFlux;

      if (Abs(hiDiff) < Abs(loDiff))
        {
          basePt = hiPt;
        }
      else
        {
          basePt = loPt;
        }
    }

  if (success)
    {
      RealVect newPt = findRoot(loPt         ,
                                hiPt         ,
                                a_edgeFlux   ,
                                a_field_data ,
                                a_tolerance  ,
                                a_blockNumber,
                                a_iv         );
      
      a_hiRes(a_iv,xDir) = newPt[xDir];
      a_hiRes(a_iv,yDir) = newPt[yDir];
    }
  else
    {
      pout()<<"blockNumber = "<<a_blockNumber<<endl;
      pout()<<"iv          = "<< a_iv<<endl;
      pout()<<"hiPt        = "<<hiPt<<endl;
      pout()<<"loPt        = "<<loPt<<endl;
      pout()<<"pt          = "<<pt<<endl;
      pout()<<"edgeFlux    = "<<a_edgeFlux<<endl;
      pout()<<"normalDir   = "<<normalDir<<endl;
      
      MayDay::Abort("couldn't bracket edgeFlux");
    }
}

void checkFluxAlignment(FArrayBox         & a_mapping  ,
                        const Vector<Box> & a_edgeBoxes,    
                        const std::string & a_blockName,
                        const FieldData   & a_field_data)
{
  int xComp = 0;
  int yComp = 1;
  Real maxFlux = -LARGEREALVAL;
  Real minFlux = LARGEREALVAL;
  IntVect maxFluxIv;
  IntVect minFluxIv;

  // iterate over edge
  Box leftEdge  = a_edgeBoxes[3];
  Box rightEdge = a_edgeBoxes[1];

  IntVect lowerLeft  = a_edgeBoxes[3].smallEnd();
  IntVect lowerRight = a_edgeBoxes[1].smallEnd();

  IntVect upperRight = a_edgeBoxes[1].bigEnd();
  IntVect upperLeft  = a_edgeBoxes[3].bigEnd();

  RealVect  lowerLeftPt (a_mapping(lowerLeft ,xComp),a_mapping(lowerLeft ,yComp));
  RealVect  lowerRightPt(a_mapping(lowerRight,xComp),a_mapping(lowerRight,yComp));
  
  RealVect  upperRightPt(a_mapping(upperRight,xComp),a_mapping(upperRight,yComp));
  RealVect  upperLeftPt (a_mapping(upperLeft ,xComp),a_mapping(upperLeft ,yComp)); 

  for (BoxIterator bit(leftEdge);bit.ok();++bit) 
    {
      // for each point move normal to the boundary towards alignment
      IntVect iv = bit();
      Real x = a_mapping(iv,xComp);
      Real y = a_mapping(iv,yComp);
      RealVect pt(x,y);
      Real thisFlux = a_field_data.getMagneticFluxFromDCT(pt);
      pout()<<"iv = "<< iv<<"(x,y) = "<<pt<<" flux = "<<thisFlux<<endl;
      if (thisFlux > maxFlux)
        {
          maxFlux = thisFlux;
          maxFluxIv = iv;
        }

       if (thisFlux < minFlux)
        {
          minFlux = thisFlux;
          minFluxIv = iv;
        }
    }
   
  pout()<<"left side for block = "<<a_blockName<<endl;
  pout()<<"maxFlux minus minFlux = "<<maxFlux - minFlux<<endl;
  pout()<<"maxFlux = "<<maxFlux<<" at iv = "<< maxFluxIv <<endl;
  pout()<<"minFlux = "<<minFlux<<" at iv = "<< minFluxIv <<endl;

  
  maxFlux = -LARGEREALVAL;
  minFlux = LARGEREALVAL;
  for (BoxIterator bit(rightEdge);bit.ok();++bit) 
    {
      // for each point move normal to the boundary towards alignment
      IntVect iv = bit();
      Real x = a_mapping(iv,xComp);
      Real y = a_mapping(iv,yComp);
      RealVect pt(x,y);
      Real thisFlux = a_field_data.getMagneticFluxFromDCT(pt);
      if (thisFlux > maxFlux)
        {
          maxFlux = thisFlux;
          maxFluxIv = iv;
        }

       if (thisFlux < minFlux)
        {
          minFlux = thisFlux;
          minFluxIv = iv;
        }
    }
  
  pout()<<"right side for block = "<<a_blockName<<endl;
  pout()<<"maxFlux minus minFlux = "<<maxFlux - minFlux<<endl;
  pout()<<"maxFlux = "<<maxFlux<<" at iv = "<< maxFluxIv <<endl;
  pout()<<"minFlux = "<<minFlux<<" at iv = "<< minFluxIv <<endl;
}

void alignBoundary(FArrayBox         & a_hiRes                ,
                   const Vector<Box> & a_hiResEdgeBoxes       ,    
                   const std::string & a_geometryFileNameNewBd,
                   const std::string & a_fieldFileNameNewBd   , 
                   const std::string & a_blockName            )
{
  // indexing macros
  int xDir = 0;
  int yDir = 1;

  //int xiDir   = 0;
  int etaDir = 1;

  Real tolerance = 1.0e-12;
  int numSteps = 50000;

  int blockNumber = getBlockNumber(a_blockName);

  Box topBox =  a_hiResEdgeBoxes[2];
  IntVect smallEnd = topBox.smallEnd();
  IntVect bigEnd   = topBox.bigEnd  ();

  RealVect smallEndPt (a_hiRes(smallEnd,xDir),a_hiRes(smallEnd,yDir));
  RealVect bigEndPt   (a_hiRes(bigEnd  ,xDir),a_hiRes(bigEnd  ,yDir));

  Real dist = (smallEndPt[xDir] - bigEndPt[xDir])*(smallEndPt[xDir] - bigEndPt[xDir]);
  dist     += (smallEndPt[yDir] - bigEndPt[yDir])*(smallEndPt[yDir] - bigEndPt[yDir]);
  dist = sqrt(dist);
  
  int numCells = smallEnd[0] - bigEnd[0];
  dist /= (numCells*1000);

  RealVect fieldDir  = LARGEREALVAL*RealVect::Unit;
  RealVect normalDir = LARGEREALVAL*RealVect::Unit;

 // Construct the field data object
  FieldData field_data(a_geometryFileNameNewBd,
                       a_fieldFileNameNewBd   ,
                       a_blockName            );

  // field aligned edges
  Box leftEdge  = a_hiResEdgeBoxes[3];
  Box rightEdge = a_hiResEdgeBoxes[1];

  // corners are set correctly
  leftEdge.grow (etaDir,-1);
  rightEdge.grow(etaDir,-1);

  Real leftEdgeFlux  = field_data.getMagneticFluxFromDCT(smallEndPt);
  Real rightEdgeFlux = field_data.getMagneticFluxFromDCT(bigEndPt  );

  // iterate through the const xi boundaries
  
  for (BoxIterator bit(leftEdge);bit.ok();++bit) 
    {
      // for each point move normal to the boundary towards alignment
      IntVect iv = bit();
      
      alignBdPt(a_hiRes      ,
                iv           ,
                field_data   ,
                numSteps     ,
                dist         , 
                leftEdgeFlux ,
                blockNumber  ,
                tolerance    );
    }     
  for (BoxIterator bit(rightEdge);bit.ok();++bit) 
    {
      // for each point move normal to the boundary towards alignment
      IntVect iv = bit();
      alignBdPt(a_hiRes      ,
                iv           ,
                field_data   ,
                numSteps     ,
                dist         , 
                rightEdgeFlux,
                blockNumber  ,
                tolerance    );
    }
}

void writeFluxFab(const FArrayBox   & a_mapping   ,
                  const std::string & a_fileName  ,
                  const int         & a_baseRes   ,
                  const std::string & a_blockName ,
                  const RealVect    & a_dx        ,
                  const FieldData   & a_field_data)
{
  int numComp = 1;
  int xComp = 0;
  int yComp = 1;
  int fluxComp = 0;

  FArrayBox fluxFab(a_mapping.box(),numComp);
  for (BoxIterator bit(a_mapping.box());bit.ok();++bit) 
    {
      IntVect iv = bit();
      Real x = a_mapping(iv,xComp);
      Real y = a_mapping(iv,yComp);
      RealVect pt(x,y);
      fluxFab(iv,fluxComp) = a_field_data.getMagneticFluxFromDCT(pt);
    }

    // dummy default argument
  Vector<string> componentNames;
  char fileName[1000];

  sprintf(fileName,"diagnostic/%03d.%s.%s.2d.hdf5",a_baseRes,a_fileName.c_str(),a_blockName.c_str());
  writeFABname(&fluxFab,fileName,componentNames,a_dx[0]);
}


void computeEvenSpacingCurr(Vector<FArrayBox*>         & a_evenSpacingCurr      , 
                            Vector<FArrayBox*>         & a_evenSpacingFlux      , 
                            const int                  & a_refRatio             ,
                            const Vector<IntVect>      & a_hiCorner             ,
                            const int                  & a_numBlocks            ,
                            const Vector<IntVect>      & a_hiResHiCorner        ,
                            const Vector<int>          & a_tanDir               ,
                            const Vector<RealVect>     & a_dx                   ,
                            const std::string          & a_geometryFileNameNewBd,
                            const std::string          & a_fieldFileNameNewBd   ,
                            const bool                 & a_writeInfo            ,
                            const Vector<Box>          & a_box                  ,
                            const Vector<Vector<Box> > & a_edgeBoxes            ,
                            const Vector <Box>         & a_hiResBox             ,
                            const Vector<Vector<Box> > & a_hiResEdgeBoxes       ,
                            const int                  & a_baseRes              ,
                            const bool                 & a_moveMcoreMcsol       ,
                            const int                  & a_beginBlockInit       ,
                            const int                  & a_endBlockInit         )
{
  int numComp = 2;
  int oneComp = 1;

  // Define hiRes and xi_eta
  Vector<FArrayBox*> hiRes(a_numBlocks);
  Vector<FArrayBox*> xi_eta(a_numBlocks);
  
  // Define arcLength, hiResArcLength, and targetSpacing
  Vector<Vector< FArrayBox*> > arcLength     (a_numBlocks);
  Vector<Vector< FArrayBox*> > hiResArcLength(a_numBlocks);
  Vector<Vector< FArrayBox*> > targetSpacing (a_numBlocks);

  //allocate memory
  for (int iBlock = a_beginBlockInit; iBlock < a_endBlockInit; ++iBlock)
    {
      // FArray boxes
      hiRes           [iBlock] = new FArrayBox(a_hiResBox[iBlock],numComp);
      xi_eta          [iBlock] = new FArrayBox(a_box     [iBlock],numComp);
      
      // data holders for edges
      arcLength        [iBlock].resize(4);
      targetSpacing    [iBlock].resize(4);
      hiResArcLength   [iBlock].resize(4);
      
      for(int iEdge = 0;iEdge< 4; ++iEdge)
        {
          // ordering of edges: bottom,right,top,left
          arcLength        [iBlock][iEdge] = new FArrayBox(a_edgeBoxes     [iBlock][iEdge],oneComp);
          targetSpacing    [iBlock][iEdge] = new FArrayBox(a_edgeBoxes     [iBlock][iEdge],oneComp);
          hiResArcLength   [iBlock][iEdge] = new FArrayBox(a_hiResEdgeBoxes[iBlock][iEdge],oneComp);
        }
    }

  for (int iBlock = a_beginBlockInit; iBlock < a_endBlockInit; ++iBlock)
    { 
      

      
      Vector<Real>evenSpacingVal(4);
      
      std::string blockName;
      getBlockName(blockName,
                   iBlock  );
      
      // Construct the field data object
      FieldData field_data(a_geometryFileNameNewBd,
                           a_fieldFileNameNewBd   ,
                           blockName              );
      
      // hi Res mapped coordinates of the unit square
      setPhysicalCoordinates(*hiRes         [iBlock],
                             a_hiResHiCorner[iBlock],
                             blockName              ,
                             a_geometryFileNameNewBd,
                             a_fieldFileNameNewBd   );

      writeFluxFab(*(hiRes   [iBlock])   ,
                   "hiResFluxBeforeAlign",
                   a_hiCorner[iBlock][1] ,
                   blockName             ,  
                   a_dx[iBlock]          ,
                   field_data            );   

      // set the constant-xi block boundary edges to have the same 
      alignBoundary(*(hiRes         [iBlock]),
                    a_hiResEdgeBoxes[iBlock] , 
                    a_geometryFileNameNewBd  ,
                    a_fieldFileNameNewBd     ,
                    blockName                );

      //      writeFluxFab(*(hiRes   [iBlock])  ,
      //            "hiResFluxAfterAlign"    ,
      //             a_hiCorner[iBlock][1],
      //             blockName            ,  
      //             a_dx[iBlock]         ,
      //             field_data           );   
      
      // mapped coordinates of the unit square
      setPhysicalCoordinates(*a_evenSpacingCurr[iBlock],
                             a_hiCorner        [iBlock],
                             blockName                 ,
                             a_geometryFileNameNewBd   ,
                             a_fieldFileNameNewBd      );
                
      
      // discretize tangent to the boundary 
      Vector <Real> totalArcLength(4);
      for (int iEdge = 0; iEdge < 4;++iEdge)    
        {
          computeArcLength(*(arcLength     [iBlock][iEdge]),
                           *(hiResArcLength[iBlock][iEdge]),
                           evenSpacingVal          [iEdge] ,
                           totalArcLength          [iEdge] ,  
                           *(hiRes         [iBlock]       ),
                           a_tanDir                [iEdge] ,
                           a_refRatio                      );
        }
      
      // data holder for calculating the even spacing map
      Vector<FArrayBox*> arcLengthInc(4);
      
      // dummy data holder
      Vector<FArrayBox*> incrementDummy(4);

      for (int iEdge = 0; iEdge < 4; ++iEdge)
        {
          if ( (iBlock < 8) || (iEdge == 0) || (iEdge == 2) || (a_moveMcoreMcsol) ||true)
            {
              // always do this.Clean up the code after checking this change is good.
              arcLengthInc             [iEdge] = new FArrayBox((arcLength[iBlock])[iEdge]->box(),1);
              (*(arcLengthInc          [iEdge])).setVal(evenSpacingVal[iEdge]);
              
              incrementDummy           [iEdge] = new FArrayBox((arcLength[iBlock])[iEdge]->box(),1);
              (*(incrementDummy        [iEdge])).setVal(LARGEREALVAL);
              
              (*(targetSpacing[iBlock])[iEdge] ).setVal( evenSpacingVal[iEdge]);
            }
          else
            {                                                           
              // Never do this else. Clean up the code after checking this change is good.
              arcLengthInc             [iEdge] = new FArrayBox((arcLength[iBlock])[iEdge]->box(),1);
              (*(arcLengthInc          [iEdge])).copy(*((arcLength[iBlock])[iEdge]));

              incrementDummy           [iEdge] = new FArrayBox((arcLength[iBlock])[iEdge]->box(),1);
              (*(incrementDummy        [iEdge])).setVal(LARGEREALVAL);

              (*(targetSpacing[iBlock])[iEdge] ).copy(*((arcLength[iBlock])[iEdge]));
            }

        }

      // move boundary points
      incrementBd(*(a_evenSpacingCurr[iBlock]),
                  *(xi_eta           [iBlock]),
                  arcLengthInc                ,
                  targetSpacing      [iBlock] ,
                  incrementDummy              ,
                  *(hiRes            [iBlock]),
                  a_refRatio                  ,
                  hiResArcLength     [iBlock] ,
                  a_writeInfo                 ,
                  blockName                   );
      
      char fileName[1000];
      Vector<string> componentNames;
      sprintf(fileName,"diagnostic/%03d.evenBeforeInterpolaton.%s.2d.hdf5",a_baseRes,blockName.c_str());
      writeFABname(a_evenSpacingCurr[iBlock],fileName,componentNames,a_dx[iBlock][0]);
      
      interpolateInteriorXiEta(*(xi_eta[iBlock]));
      
      // use interior xi_eta to upgrade initial guess on interior
      FArrayBox tempEvenFab((a_evenSpacingCurr[iBlock])->box(),(a_evenSpacingCurr[iBlock])->nComp());
      field_data.getPhysicalCoordinates(*(xi_eta[iBlock]), 
                                        tempEvenFab      );
      
      // restore values on the boundary to the values established in the the call to incrementBd
      copyBoundary(tempEvenFab                ,
                   *(a_evenSpacingCurr[iBlock]));
      
      // rename newCurr curr
      a_evenSpacingCurr[iBlock]->copy(tempEvenFab);
      
      //fill in magnetic flux
      FArrayBox& evenSpacingFab = *(a_evenSpacingCurr[iBlock]);
      
      // writeFluxFab(evenSpacingFab          ,
      //             "evenSpacingFluxBeforeAlign",
      //             a_hiCorner[iBlock][1]   ,
      //             blockName               ,  
      //             a_dx[iBlock]            ,
      //             field_data              );

      // set the constant-xi block boundary edges to have the same 
      alignBoundary(evenSpacingFab          ,
                    a_edgeBoxes[iBlock]     , 
                    a_geometryFileNameNewBd ,
                    a_fieldFileNameNewBd    ,
                    blockName               );

      writeFluxFab(evenSpacingFab          ,
                   "evenSpacingFluxAfterAlign",
                   a_hiCorner[iBlock][1]   ,
                   blockName               ,  
                   a_dx[iBlock]            ,
                   field_data              );

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
    pp.get("log2MinRes_newBdNodesIC",log2MinRes);
    pp.get("log2MaxRes_newBdNodesIC",log2MaxRes);

    // prolixity 
    bool writeInfo = false;
    pp.query("writeInfo",writeInfo) ;
    
    // prolixity 
    bool writePlot = false;
    pp.query("writePlot",writePlot) ;

    int numBlocks = -LARGEINTVAL;
    pp.get("numBlocks", numBlocks);

    int refRatioEtaInput = LARGEINTVAL;
    pp.get("refRatioEta",refRatioEtaInput);
    refRatioEtaInput -= 1;

    Vector<int> refRatioXi (numBlocks);
    Vector<int> refRatioEta(numBlocks);
    int beginBlockInit;           
    int endBlockInit  ;
    pp.get("beginBlockInit",beginBlockInit);
    pp.get("endBlockInit",endBlockInit);
    
    for (int iBlock = beginBlockInit; iBlock <= endBlockInit;++iBlock)
      {
        refRatioXi [iBlock] = 0;
        refRatioEta[iBlock] = 0;
      }

    // set mcore and mcsol
    refRatioEta[8] = refRatioEtaInput;
    refRatioEta[9] = refRatioEtaInput;
    
    for (int numberCells = log2MinRes; numberCells <= log2MaxRes; ++ numberCells)
      {
        // different numRow and numCol for each block
        Vector<int> numRow(numBlocks);
        Vector<int> numCol(numBlocks);

        // baseRes used for making file names
        int baseRes = pow(2,numberCells);
        
        // different hi corner for each block
        Vector<IntVect> hiCorner(numBlocks);
        
         // dx for each block
        Vector<RealVect> dx(numBlocks);

        // different loEnd and hiEnd for each block 
        Vector<IntVect> loEnd(numBlocks);
        Vector<IntVect> hiEnd(numBlocks);

        // different box for each block
        Vector < Box > box(numBlocks);
        
        // refRatio for hi res discretization of the boundary
        int refRatio = LARGEREALVAL;
        pp.get("refRatio",refRatio);

        Vector<RealVect> hiResDx(numBlocks);

        // differnt hiResBox for each block
        Vector<Box> hiResBox(numBlocks);

        // different hiRes hi corner for each block
        Vector<IntVect> hiResHiCorner(numBlocks);

        // different hi Res numRow numCol for each block
        Vector<int> hiResNumRow(numBlocks);
        Vector<int> hiResNumCol(numBlocks);

        for (int iBlock = beginBlockInit; iBlock <= endBlockInit; ++ iBlock)
          {
            numRow[iBlock] = pow(2,numberCells + refRatioEta[iBlock]);
            numCol[iBlock] = pow(2,numberCells + refRatioXi [iBlock]);

            dx[iBlock][0] = 1.0/Real(numCol[iBlock]);
            dx[iBlock][1] = 1.0/Real(numRow[iBlock]);

            loEnd[iBlock] = IntVect::Zero;
            hiEnd[iBlock][0] = numCol[iBlock];
            hiEnd[iBlock][1] = numRow[iBlock];
            
            box     [iBlock].define(loEnd[iBlock],hiEnd[iBlock]);
            hiCorner[iBlock]= box[iBlock].bigEnd();

            hiResNumRow  [iBlock]    = numRow[iBlock]*refRatio;
            hiResNumCol  [iBlock]    = numCol[iBlock]*refRatio;

            hiResDx      [iBlock][0] = 1.0/Real(hiResNumCol[iBlock]);
            hiResDx      [iBlock][1] = 1.0/Real(hiResNumRow[iBlock]);

            hiResBox     [iBlock].define(loEnd[iBlock],refRatio*hiEnd[iBlock]);
            hiResHiCorner[iBlock]= hiResBox[iBlock].bigEnd();
          }          
    
        // mapped quantities have two components
        int numComp    = 2 ;
        
        Vector<FArrayBox*>curr            (numBlocks);
        Vector<FArrayBox*>xi_eta          (numBlocks);
        Vector<FArrayBox*>evenSpacingCurr (numBlocks);
        Vector<FArrayBox*>evenSpacingFlux (numBlocks);
        Vector<FArrayBox*>hiRes           (numBlocks);
        Vector<FArrayBox*>initDisc        (numBlocks);
        Vector<FArrayBox*>oldPhysicalCoord(numBlocks);
        Vector<FArrayBox*>newPhysicalCoord(numBlocks);

        // Currently: COGENT_mapping_3-22-16
        std::string geometryFileNameNewBd;
        pp.get("geometryFileNameNewBd", geometryFileNameNewBd);

        //Currently: DCT_coefficients_3-22-16
        std::string fieldFileNameNewBd;
        pp.get("fieldFileNameNewBd", fieldFileNameNewBd);
        
        for (int iBlock = beginBlockInit; iBlock <= endBlockInit; ++iBlock)
          {
            //discretizes the boundary at high resolution
            hiRes           [iBlock] = new FArrayBox(hiResBox[iBlock],numComp     );
            curr            [iBlock] = new FArrayBox(box     [iBlock],numComp     );
            xi_eta          [iBlock] = new FArrayBox(box     [iBlock],numComp     );
            evenSpacingCurr [iBlock] = new FArrayBox(box     [iBlock],numComp     );
            evenSpacingFlux [iBlock] = new FArrayBox(box     [iBlock],oneComp     );
            initDisc        [iBlock] = new FArrayBox(box     [iBlock],numComp     );
            oldPhysicalCoord[iBlock] = new FArrayBox(box     [iBlock],numComp     );
            newPhysicalCoord[iBlock] = new FArrayBox(box     [iBlock],numComp     );
          }

        // write fabs corresponding to the input map, called "inputMap"
        for (int iBlock = beginBlockInit; iBlock <= endBlockInit; ++iBlock)
          {
            std::string blockName;
            getBlockName(blockName,
                         iBlock  );
                    
            // Construct the field data object
            FieldData field_data(geometryFileNameNewBd,
                                 fieldFileNameNewBd   ,
                                 blockName            );
                
            // mapped coordinates of the unit square
            setPhysicalCoordinates(*(oldPhysicalCoord[iBlock]),
                                   hiEnd             [iBlock] ,
                                   blockName                  ,
                                   geometryFileNameNewBd      ,
                                   fieldFileNameNewBd         );
            
            char filename[1000];
            // dummy default argument
            Vector<string> componentNames;
            sprintf(filename,"diagnostic/%03d.inputMap.%s.2d.hdf5",baseRes,blockName.c_str());
            writeFABname(oldPhysicalCoord[iBlock],filename,componentNames,dx[iBlock][0]);
          }
    
        // unit vector field perpendicular to magnetic field lines combined with (all) four first partial derivatives
        int num1Deriv      = 4 ;
        int fluxComp       = 3 ;
        
        Vector<FArrayBox*>BField   (numBlocks);
                   
        for (int iBlock = beginBlockInit; iBlock <= endBlockInit; ++iBlock)
          {
            BField     [iBlock] = new FArrayBox(box[iBlock],numComp   + num1Deriv + fluxComp);
          }
        Real sigma_core;
        pp.get("sigma_core",sigma_core);
        if (sigma_core < 0.0)
          {
            MayDay::Abort("sigma < 0");
          }

        Real sigma_scrape_off;
        pp.get("sigma_scrape_off",sigma_scrape_off);
        if (sigma_scrape_off < 0.0)
          {
            MayDay::Abort("sigma < 0");
          }

        Real sigma_private_flux;
        pp.get("sigma_private_flux",sigma_private_flux);
        if (sigma_private_flux < 0.0)
          {
            MayDay::Abort("sigma < 0");
          }

        Real expo_core;
        pp.get("expo_core", expo_core);
        if (expo_core < 0.0)
          {
            MayDay::Abort("exponent < 0");
          }

        Real expo_scrape_off;
        pp.get("expo_scrape_off", expo_scrape_off);
        if (expo_scrape_off < 0.0)
          {
            MayDay::Abort("exponent < 0");
          }
        
        Real expo_private_flux;
        pp.get("expo_private_flux", expo_private_flux);
        if (expo_private_flux < 0.0)
          {
            MayDay::Abort("exponent < 0");
          }
        
        Vector<Real> sigmaBlock(numBlocks);
        Vector<Real> expoBlock(numBlocks);
        for (int iBlock = beginBlockInit; iBlock <= endBlockInit; ++iBlock)
          {
            if ((iBlock == 0) || (iBlock == 1))
              {
                sigmaBlock[iBlock] = sigma_core;
                expoBlock[iBlock]  = expo_core;
              }
            else if (iBlock== 2 || iBlock == 3 || iBlock == 4 || iBlock == 5)
              {
                sigmaBlock[iBlock] = sigma_scrape_off;
                expoBlock[iBlock]  = expo_scrape_off;
              }
            else if (iBlock == 6 || iBlock == 7)
              {
                sigmaBlock[iBlock] = sigma_private_flux;
                expoBlock[iBlock]  = expo_private_flux;
              }
          }
      
        // evenly space boundary nodes along the long edges of mcore and mcsol
        bool moveMcoreMcsol = false;
        pp.get("moveMcoreMcsol",moveMcoreMcsol);
        
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
        
        // boxes coresponding to boundary values on four sides
        Vector<Vector<Box> > edgeBoxes     (numBlocks);
        Vector<Vector<Box> > hiResEdgeBoxes(numBlocks);
        for (int iBlock = beginBlockInit; iBlock <= endBlockInit; ++iBlock)
          {
            edgeBoxes[iBlock]     .resize(4);
            hiResEdgeBoxes[iBlock].resize(4);

            // boxes used for tangent vector FArrayBoxes
            createEdgeBoxes(edgeBoxes     [iBlock],box     [iBlock]);
            createEdgeBoxes(hiResEdgeBoxes[iBlock],hiResBox[iBlock]);
          }
        
        // ordering of edges is bottom,right,top,left
        Vector<int> tanDir(4); 
        tanDir[0] = xDir;
        tanDir[1] = yDir;
        tanDir[2] = xDir;
        tanDir[3] = yDir;
        
        // FArrayBoxes for tangent vector to the (hiRes representation of the) boundary
        for (int iBlock = beginBlockInit; iBlock <= endBlockInit; ++iBlock)
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

        computeEvenSpacingCurr(evenSpacingCurr      ,
                               evenSpacingFlux      ,
                               refRatio             ,
                               hiCorner             ,
                               numBlocks            ,
                               hiResHiCorner        ,
                               tanDir               ,
                               dx                   ,
                               geometryFileNameNewBd,
                               fieldFileNameNewBd   ,
                               writeInfo            ,
                               box                  ,
                               edgeBoxes            ,
                               hiResBox             ,
                               hiResEdgeBoxes       ,
                               baseRes              ,
                               moveMcoreMcsol       ,
                               beginBlockInit       ,
                               endBlockInit         );
        
        // clean start for mapping file output function    
        char fileName0[1000];
        sprintf(fileName0,"outputMapping/even.curr.%03d.2d.txt",numRow[0] + 1);
        clearExtendedMappingFile(fileName0);
       
        // write fabs corresponding to the map based on even spaced bd nodes, called evenCurr
        for (int iBlock = beginBlockInit; iBlock <= endBlockInit; ++iBlock)
          {
            std::string blockName;
            getBlockName(blockName,
                         iBlock   );
           
            // write unghosted mapping file
            outputMapping(*(evenSpacingCurr[iBlock]),
                          *(BField         [iBlock]),
                          IntVect::Zero             ,
                          blockName                 ,
                          fileName0                 );
        
            // initial mapping
            char fileName[1000];
            // dummy default argument
            Vector<string> componentNames;
            sprintf(fileName,"diagnostic/%03d.evenAfterInterpolation.%s.2d.hdf5",baseRes,blockName.c_str());
            writeFABname(evenSpacingCurr[iBlock],fileName,componentNames,dx[iBlock][0]);

            sprintf(fileName,"diagnostic/%03d.evenFlux.%s.2d.hdf5",baseRes,blockName.c_str());
            writeFABname(evenSpacingFlux[iBlock],fileName,componentNames,dx[iBlock][0]);
          }
        
        for (int iBlock = beginBlockInit; iBlock <= endBlockInit; ++iBlock)
          {
            pout()<< "starting block number = "<<iBlock<<endl;

            Vector<Real>evenSpacingVal(4);

            std::string blockName;
            getBlockName(blockName,
                         iBlock  );
                    
            // Construct the field data object
            FieldData field_data(geometryFileNameNewBd,
                                 fieldFileNameNewBd   ,
                                 blockName            );
                
              
            // hi Res mapped coordinates of the unit square
            setPhysicalCoordinates(*hiRes       [iBlock],
                                   hiResHiCorner[iBlock],
                                   blockName            ,
                                   geometryFileNameNewBd,
                                   fieldFileNameNewBd   );
                    
            // mapped coordinates of the unit square
            setPhysicalCoordinates(*newPhysicalCoord[iBlock],
                                   hiCorner         [iBlock],
                                   blockName                ,
                                   geometryFileNameNewBd    ,
                                   fieldFileNameNewBd       );
                    
            // initialize with output from FieldData interpolation
            curr[iBlock]->copy(*newPhysicalCoord[iBlock]);
                    
            checkCornerFlux(*(curr[iBlock])      ,
                            blockName            ,
                            geometryFileNameNewBd,
                            fieldFileNameNewBd   );
    
            // discretize tangent to the boundary 
            for (int iEdge = 0; iEdge < 4;++iEdge)    
              {
                computeBoundaryTangentVector(*(tangent[iBlock][iEdge]),
                                             *(curr   [iBlock]       ),
                                             tanDir           [iEdge] ,
                                             dx       [iBlock]        );
              }
                    
            // arc lengths on all four edges 
            Vector <Real> totalArcLength(4);
            for (int iEdge = 0; iEdge < 4;++iEdge)    
              {
                computeArcLength(*(arcLength     [iBlock][iEdge])  ,
                                 *(hiResArcLength[iBlock][iEdge])  ,
                                 evenSpacingVal          [iEdge]   ,
                                 totalArcLength          [iEdge]   ,  
                                 *(hiRes         [iBlock]       )  ,
                                 tanDir                  [iEdge]   ,
                                 refRatio                          );
              }
            
            // record the magentic flux at the boundaries with evenly spaced nodes on the top edges of blocks 8 and 9 and the bottom of block 7
            Vector<Vector<Real> > fluxValues;
            getFluxValueVector(fluxValues           ,
                               evenSpacingCurr      ,
                               edgeBoxes            ,
                               numCol               ,
                               geometryFileNameNewBd,
                               fieldFileNameNewBd   ,
                               beginBlockInit       ,
                               endBlockInit         );

            // compute increments for a weighted combination of flux matching and evenly spaced nodes
            for (int iEdge = 0; iEdge < 4;++iEdge)    
              {
                computeFluxMatchedVal(*(fluxMatchedVal [iBlock][iEdge]),
                                      *(evenSpacingCurr[iBlock]       ),
                                      hiResArcLength   [iBlock]        ,
                                      *(hiRes          [iBlock]       ),
                                      evenSpacingVal                   ,
                                      iEdge                            ,
                                      iBlock                           ,
                                      hiResDx          [iBlock]        ,
                                      refRatio                         ,
                                      geometryFileNameNewBd            ,
                                      fieldFileNameNewBd               ,
                                      fluxValues                       );
                        
                Real totalTargetArcLength = LARGEREALVAL;
                computeTargetSpacingVal(*(targetSpacing   [iBlock][iEdge]),
                                        totalTargetArcLength              ,
                                        *(arcLength       [iBlock][iEdge]),
                                        evenSpacingVal            [iEdge] ,
                                        *(fluxMatchedVal  [iBlock][iEdge]),
                                        iEdge                             ,
                                        sigmaBlock        [iBlock]        ,
                                        expoBlock         [iBlock]        ,
                                        iBlock                            ,
                                        dx                [iBlock]        ,
                                        tanDir                    [iEdge] ,
                                        totalArcLength            [iEdge] ,
                                        moveMcoreMcsol                    );
                          
                Real ratio = totalArcLength[iEdge]/totalTargetArcLength;
                if (Abs(ratio - 1.0) > 1.0e-13)
                  {
                    pout()<<"ratio - 1.0 = "<<ratio - 1.0 <<endl;
                    MayDay::Abort("check arc length");
                  }

                //*(targetSpacing [iBlock][iEdge])*= ratio;
              }
                             
            // data holder for boundary increments
            Vector<FArrayBox*> arcLengthInc(4);            for (int iEdge = 0; iEdge < 4; ++iEdge)
              {
                arcLengthInc [iEdge] = new FArrayBox(tangent[iBlock][iEdge]->box(),1);
              }
                
            incrementBd(*(curr        [iBlock]),
                        *(xi_eta      [iBlock]), 
                        arcLength     [iBlock] ,
                        targetSpacing [iBlock] ,
                        arcLengthInc           ,
                        *(hiRes       [iBlock]),
                        refRatio               ,
                        hiResArcLength[iBlock] ,
                        writeInfo              ,
                        blockName              );
            
             char fileName[1000];
             Vector<string> componentNames;
            sprintf(fileName,"diagnostic/%03d.alignedBeforeInterpolaton.%s.2d.hdf5",baseRes,blockName.c_str());
            writeFABname(curr[iBlock],fileName,componentNames,dx[iBlock][0]);

            interpolateInteriorXiEta(*(xi_eta[iBlock]));
            
            // use interior xi_eta to upgrade initial guess on interior
            FArrayBox tempCurrFab((curr[iBlock])->box(),(curr[iBlock])->nComp());
            field_data.getPhysicalCoordinates(*(xi_eta[iBlock]), 
                                              tempCurrFab      );

            // restore values on the boundary to the values established in the the call to incrementBd
            copyBoundary(tempCurrFab    ,
                         *(curr[iBlock]));

            // rename newCurr curr
            curr[iBlock]->copy(tempCurrFab);

            // set the constant-xi block boundary edges to have the same flux
            alignBoundary(*(curr   [iBlock])   ,
                          edgeBoxes[iBlock]    , 
                          geometryFileNameNewBd,
                          fieldFileNameNewBd   ,
                          blockName            );

            writeFluxFab(*(curr   [iBlock]) ,
                         "finalFluxValues"  ,
                         hiCorner[iBlock][1],
                         blockName          ,  
                         dx[iBlock]         ,
                         field_data         );   


            //            tempCurrFab -= *(curr[iBlock]);
            
            //checkFluxAlignment(*(curr   [iBlock]),
            //                   edgeBoxes[iBlock] , 
            //                   blockName         ,
            //                   field_data        );
          }
        
       

        // clean start for mapping file output function    
        char fileName[1000];
        sprintf(fileName,"outputMapping/aligned.curr.%03d.2d.txt",numRow[0] + 1);
        clearExtendedMappingFile(fileName);

        for (int iBlock = beginBlockInit; iBlock <= endBlockInit; ++iBlock)
          {
            std::string blockName;
            getBlockName(blockName,
                         iBlock   );
            
            // write unghosted mapping file
            outputMapping(*(curr   [iBlock]),
                          *(BField[iBlock]) ,
                          IntVect::Zero     ,
                          blockName         ,
                          fileName          );
            
            // aligned mapping
            char fileName[1000];
            // dummy default argument
            Vector<string> componentNames;
          
            sprintf(fileName,"diagnostic/%03d.alignedAfterInterpolaton.%s.2d.hdf5",baseRes,blockName.c_str());
            writeFABname(curr[iBlock],fileName,componentNames,dx[iBlock][0]);
          }

        // end loop over resolution          
        //delete pointers (last in first out)
        for (int iBlock = beginBlockInit; iBlock <= endBlockInit; ++iBlock)
          {
            delete BField          [iBlock];
            delete curr            [iBlock];
            delete initDisc        [iBlock];
            delete oldPhysicalCoord[iBlock];
            delete newPhysicalCoord[iBlock];
          }
      }
  }// End scoping trick 
  
#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
  
  return 0;
}


