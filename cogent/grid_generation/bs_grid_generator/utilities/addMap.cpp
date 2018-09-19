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

void setBField(FArrayBox         & a_BField          ,
               const FArrayBox   & a_newPhysicalCoord,
               const RealVect    & a_dx              ,
               const std::string & a_geometryFileName,
               const std::string & a_fieldFileName   , 
               const std::string & a_blockName       )
{
  // Construct the field data object
  FieldData field_data(a_geometryFileName,
                       a_fieldFileName   ,
                       a_blockName       );
  
  // normalized B field and four first partials  
  field_data.getFieldUnitVectorFromMappingFile(a_newPhysicalCoord,
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



void clearExtendedMappingFile(const FArrayBox & a_extendedMap    ,
                              const bool      & a_evenSpacingCase)
{
  IntVect boxSize = a_extendedMap.box().size();
  int row = boxSize[0];
    
  char filename[80];

  if (!a_evenSpacingCase)
    {
      sprintf(filename,"noGhostMap.curr.%03d.2d.txt",row);
      FILE* fd;
      fd = fopen (filename, "w");
      if (fd == NULL) 
        {
          perror ("Error opening file");
        }
      fclose(fd);
    }
  else
    {
      sprintf(filename,"even.curr.%03d.2d.txt",row);
      FILE* fd;
      fd = fopen (filename, "w");
      if (fd == NULL) 
        {
          perror ("Error opening file");
        }
      fclose(fd);
    }
}
void outputMapping(const FArrayBox   & a_extendedMap    ,
                   const FArrayBox   & a_BField         ,
                   const IntVect     & a_numGhost       ,
                   const std::string & a_blockName      ,
                   const int         & a_iter           ,
                   const bool        & a_evenSpacingCase)
{
  
  IntVect boxSize = a_extendedMap.box().size();
  int row = boxSize[0];
  int col = boxSize[1];
  int numRad = row - 2*a_numGhost[0];
  int numPol = col - 2*a_numGhost[1];

 char filename[80];
  if (!a_evenSpacingCase)
    {
      pout()<<"Writing extended map for base resololution = (" << numRad<<","<<numPol<<")"<<endl;
     
      
      if (a_iter == LARGEINTVAL)
        {
          sprintf(filename,"outputMapping/aligned.curr.%03d.2d.txt",numRad);
        }
      else
        {
          sprintf(filename,"outputMapping/aligned.iter.%d.curr.%03d.2d.txt",a_iter,numRad);
        }
    }
  else
    {
      pout()<<"Writing even space map for base resolution (" << numRad<<","<<numPol<<")"<<endl;
      
      if (a_iter == LARGEINTVAL)
        {
          sprintf(filename,"outputMapping/even.curr.%03d.2d.txt",numRad);
        }
      else
        {
          sprintf(filename,"outputMapping/even.iter.%d.curr.%03d.2d.txt",a_iter,numRad);
        }
    }

      FILE* fd;
      fd = fopen (filename, "a");
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
    
    int log2MinCells;
    int log2MaxCells;
    pp.get("log2MinRes_BrackbillSalzman",log2MinCells);
    pp.get("log2MaxRes_BrackbillSalzman",log2MaxCells);

    int baseRes = pow(2,log2MinCells);

    int beginBlock = -LARGEINTVAL;
    pp.get("beginBlock", beginBlock);
    
    int endBlock = -LARGEINTVAL;
    pp.get("endBlock", endBlock);
    
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

    // prolixity 
    bool writeInfo = false;
    pp.query("writeInfo",writeInfo) ;
     
    Vector<int> numRow      (numBlocks);
    Vector<int> numCol      (numBlocks);
    Vector<IntVect> loEnd   (numBlocks);
    Vector<IntVect> hiEnd   (numBlocks);
    Vector<Box> box         (numBlocks);
    Vector<IntVect> hiCorner(numBlocks);
    Vector<RealVect> dx     (numBlocks);

    for(int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        numRow[iBlock] = pow(2,log2MinCells + refRatioEta[iBlock]);
        numCol[iBlock] = pow(2,log2MinCells + refRatioXi [iBlock]);

        loEnd[iBlock] = IntVect::Zero;
        hiEnd[iBlock][0] = numCol[iBlock]; 
        hiEnd[iBlock][1] = numRow[iBlock]; 
        
        box     [iBlock].define(loEnd[iBlock],hiEnd[iBlock]);
        hiCorner[iBlock]= box[iBlock].bigEnd();
               
        dx[iBlock][0] = 1.0/Real(numCol[iBlock]);
        dx[iBlock][1] = 1.0/Real(numRow[iBlock]);
      }

    // mapped quantities have two components
    int numComp    = 2 ;
      
    Vector<FArrayBox*>curr(numBlocks);
               
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        curr[iBlock] = new FArrayBox(box[iBlock],numComp);
      }
    
    // unit vector field perpendicular to magnetic field lines combined with (all) four first partial derivatives
    int num1Deriv      = 4 ;
    Vector<FArrayBox*>BField(numBlocks);
        
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        BField[iBlock] = new FArrayBox(box[iBlock],numComp + num1Deriv);
        BField[iBlock]->setVal(0.0);
      }
    
    std::string geometryFileName;
    pp.get("geometryFileNameAddMap", geometryFileName);
    
    std::string tempName;
    std::string fieldFileName;
    pp.get("fieldFileNameAddMap", fieldFileName);

    bool evenSpacingCase;
    pp.get("evenSpacingCase",evenSpacingCase);
   
    for (int iBlock = beginBlock; iBlock <= endBlock; ++iBlock)
      { 
        char iGeometryFileName[1000];
        if (evenSpacingCase)
          {
            sprintf(iGeometryFileName,"../block%d/outputMapping/even.curr.%03d.2d.txt",iBlock,baseRes + 1);
          }
        else
          {
            sprintf(iGeometryFileName,"../block%d/outputMapping/aligned.curr.%03d.2d.txt",iBlock,baseRes + 1);
          }
        std::string blockName;
        getBlockName(blockName,
                     iBlock  );
        
        // Construct the field data objects
        FieldData field_data(iGeometryFileName,
                             fieldFileName    ,
                             blockName        );
        
        // mapped coordinates of the unit square
        setPhysicalCoordinates(*(curr   [iBlock]),
                               hiCorner [iBlock] ,
                               blockName         ,
                               iGeometryFileName ,
                               fieldFileName     );
      }

    for (int iBlock = endBlock + 1; iBlock < numBlocks; ++iBlock)
      { 
        std::string blockName;
        getBlockName(blockName,
                     iBlock  );
        
        // Construct the field data objects
        FieldData field_data(geometryFileName,
                             fieldFileName   ,
                             blockName      );
        
        // mapped coordinates of the unit square
        setPhysicalCoordinates(*(curr  [iBlock]),
                               hiCorner[iBlock] ,
                               blockName        ,
                               geometryFileName ,
                               fieldFileName    );
      }


    clearExtendedMappingFile(*(curr[0])     ,
                             evenSpacingCase);
    
    for (int iBlock = 0; iBlock < numBlocks; ++iBlock)
      {
        std::string blockName;
        getBlockName(blockName,
                     iBlock   );
        
     int dummyVal = LARGEINTVAL;
        // write unghosted mapping file
        outputMapping(*(curr  [iBlock]),
                      *(BField[iBlock]) ,
                      IntVect::Zero     ,
                      blockName         ,
                      dummyVal          ,
                      evenSpacingCase   );
      }
    // end loop over resolution          
    //delete pointers (last in first out)
    for (int iBlock = beginBlock; iBlock < endBlock; ++iBlock)
      {
        delete BField[iBlock];
        delete curr  [iBlock];
      }
    
  }// End scoping trick 
  
#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
  
  return 0;
}


