#include <cmath>
#include <cstdio>


#include "BoxIterator.H"
#include "ParmParse.H"
#include "DebugDump.H"
#include "FABView.H"
#include "LoadBalance.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "AMRIO.H"

#include "fftw3.h"

#include <fstream>
#include <float.h>

#include "UsingNamespace.H"

#define LARGEINTVAL -999999999
#define LARGEREALVAL 123456789.0

void myWriteFab(const FArrayBox & a_fab     , 
                const Real      & a_sVal    ,
                const bool      & a_inputFab)
{
   // dummy default arguments
  Real dx = 1.0;
  Vector<string> componentNames;
  
  // for sprintf
  char name[1000];
  
  if (a_inputFab)
    {
      sprintf(name,"noisyFabIn.2d.hdf5");
      writeFABname(&a_fab,name,componentNames,dx);
    }
  else
    {
      int iDir = 0;
      int jDir = 1;
           
      int ires = a_fab.box().size()[iDir];
      int jres = a_fab.box().size()[jDir];
      
     
      sprintf(name,"smoothFabsOut/%03dx%03d.smoothFab.sVal.%g.2d.hdf5",ires,jres,a_sVal);
      writeFABname(&a_fab,name,componentNames,dx);
    }
}

void IDCT(FArrayBox & a_yOut,
          FArrayBox & a_yIn )
{
  IntVect size = a_yIn.box().size();
  int nRow = size[1];
  int nCol  = size[0];

  // ALWAYS USE THIS BEFORE INITIALIZING  the 'in' array
  fftw_plan plan = fftw_plan_r2r_2d(nRow            , 
                                    nCol            , 
                                    a_yIn.dataPtr() ,  
                                    a_yOut.dataPtr(), 
                                    FFTW_REDFT00    , 
                                    FFTW_REDFT00    ,
                                    FFTW_MEASURE    );
  // execute
  fftw_execute(plan);
  a_yOut.divide(4*(nRow - 1)*(nCol - 1));
}

void DCT(FArrayBox & a_yOut,
         FArrayBox & a_yIn )
{
  IntVect size = a_yIn.box().size();
  int nRow = size[1];
  int nCol = size[0];

  // ALWAYS USE THIS BEFORE INITIALIZING  the 'in' array
  fftw_plan plan = fftw_plan_r2r_2d(nRow            , 
                                    nCol            , 
                                    a_yIn.dataPtr() , 
                                    a_yOut.dataPtr(), 
                                    FFTW_REDFT00    , 
                                    FFTW_REDFT00    ,
                                    FFTW_MEASURE    );
  
  // execute
  fftw_execute(plan);
}

// z,y and s are the variables in Garcia's code
void smoothn(FArrayBox & a_z,
             FArrayBox & a_y,
             const Real& a_s)
{
  // use this box for new data holders
  Box box = a_y.box();
  
  //one component data
  int nComp = 1;
  int comp  = 0;
    
  //DCT of y
  FArrayBox yDCT(box,nComp);
  DCT(yDCT,a_y);

  // make Lambda, Gamma, and yDCT times Gamma
  FArrayBox Lambda   (box,nComp);
  FArrayBox Gamma    (box,nComp);
  FArrayBox yDCTGamma(box,nComp);
  
  // parameters used in loop
  IntVect boxSize = box.size();
  Real iMax = boxSize[0];
  Real jMax = boxSize[1];
  
  // loop over all iv
  for (BoxIterator bit(box);bit.ok();++bit)
    {
      IntVect iv = bit();
      Real i = iv[0];
      Real j = iv[1];

      // set Lambda
      Lambda(iv,comp) = -2.0*(2.0 - cos(M_PI*i/iMax) - cos(M_PI*j/jMax)); 
      
      // set Gamma
      Gamma (iv,comp) = 1.0/ (1.0 + a_s*Lambda(iv,0)*Lambda(iv,0)); 
        
      // multiply yDCT by Gamma
      yDCTGamma(iv,comp) = Gamma(iv,0)*yDCT(iv,0); 
    }

  // z is the smoothed result
  IDCT(a_z,yDCTGamma);
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
        pout() << "Usage:  garcia...ex <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }

    // Parse input file
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);
   
    // high corner in (i,j)-space
    Vector<int> nCell (SpaceDim);
    cin >> nCell[0];
    cin >> nCell[1];

    IntVect bigEnd(nCell[0] - 1,nCell[1] - 1);
    IntVect smallEnd = IntVect::Zero;

    // box for Fabs
    Box box(smallEnd,bigEnd);
    
    // one component data
    int nComp    = 1;
    int comp = 0;

    // input noisy array 
    FArrayBox noisyFab(box,nComp);

    for (int i = smallEnd[1]; i <= bigEnd[1]; ++i)
      {
        for (int j = smallEnd[0]; j <= bigEnd[0]; ++j)
          {
            Real yScalar;
            cin >> yScalar;
            //pout() << "{i,j] = "<<"{"<<i<<","<<j<<"} = "<< yScalar<<endl;
            noisyFab(IntVect(i,j),comp) = yScalar;
          }
      }

    // write out noisyInputFab as a FArrayBox
    bool inputFabBool = true;
    Real dummy;
    myWriteFab(noisyFab    ,
               dummy       ,
               inputFabBool);

    //set to false for later writing of smooth fabs
    inputFabBool = false;

     // smooth arrays are parameterized by a smoothing parameter, smoothingParam. 
    int numSVal;
    Real sInit;
    Real deltaS;

    // read in parameters
    pp.get("numSVal",numSVal);
    pp.get("sInit"  ,sInit);
    pp.get("deltaS",deltaS);

    // make numSVal increasingly smoothed versions
     for (int iSVal = 0 ;iSVal < numSVal;++iSVal)
      {
        // allocate data holder
        FArrayBox smoothFab(box,nComp);
        
        // set smoothing parameter
        Real smoothingParam = sInit + iSVal*deltaS;

        // smooth yn and put the result in zVec
        smoothn(smoothFab     ,
                noisyFab      ,
                smoothingParam);

        // output smoothFab
        myWriteFab(smoothFab     ,
                   smoothingParam,
                   inputFabBool  );
      }

  }// End scoping trick 
  
#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
  
  return 0;
}
