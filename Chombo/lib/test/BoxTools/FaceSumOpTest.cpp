
#include "LevelData.H"
#include "ReductionCopier.H"
#include "ReductionOps.H"
#include "FABView.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "UsingNamespace.H"
int main( int a_argc, char* a_argv[] )
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init( &a_argc, &a_argv );
  setChomboMPIErrorHandler();
#endif

  if (SpaceDim > 1)
    {
      //domain box parameters in up to five dimensions
      IntVect smallEnd = IntVect::Zero;
      IntVect bigEnd   = IntVect::Unit;
      bigEnd *= 7;
  
      // make an 8x8x8x... domainBox
      pout()<<"bigEnd = "<<bigEnd<<endl;
  
      Box domainBox(smallEnd, bigEnd);
      ProblemDomain probDom(domainBox);

      pout()<<"domainBox = "<<domainBox<<endl;
      pout()<<endl;

      // to do: make this a default argument of the copier constructor
      bool exchange = false;
  
      // reduce in the 0 and 1 directions 
      IntVect reduceDir = IntVect::Zero;
      reduceDir[0] = 1;
      reduceDir[1] = 1;
  
      int numReduceDir = 2;
  
      // Make two domain filling adjacent boxes
      int n = 4;

      IntVect big1 = 7*IntVect::Unit;
      big1[SpaceDim - 1] = n-1;

      IntVect small2 = IntVect::Zero;
      small2[SpaceDim -1] = n;
  

      Box oldBox1(IntVect::Zero,big1);
      Box oldBox2(small2,domainBox.bigEnd());

      // print the boxes
      pout() << "oldBox 1 = " << oldBox1 << ", oldBox 2 = " << oldBox2 << endl;
      pout()<<endl;
  
      Vector<Box> oldBoxes;
      oldBoxes.push_back(oldBox1);
      oldBoxes.push_back(oldBox2);

      // set the processor ids
      Vector<int> oldPids;
      oldPids.push_back(0);
      oldPids.push_back(1);
      DisjointBoxLayout oldDbl(oldBoxes, oldPids, domainBox);
  
      // define the source level data
      LevelData<FluxBox> ld(oldDbl       , 
                            1            , 
                            IntVect::Zero);
  
      for (DataIterator dit(oldDbl); dit.ok(); ++dit) 
        {
          pout()<<"oldDbl[dit] = "<<oldDbl[dit] <<endl;
          pout()<<endl;
            
          ld[dit].setVal(1.0);
        }

      //use this parameter to call the new constructor
      bool usingFaceSumOp = true; 

      // Make an LD<FluxBox> that is (in 2D) degenerate in the 0 direction or (in 3D) degenerate in the 0 and 1 directions;
      Vector<DisjointBoxLayout> tempDestGrids(numReduceDir + 1);

      tempDestGrids[0] = ld.getBoxes();

      int index = -1;
      for (int idir = 0; idir <SpaceDim; ++idir)
        {
          if (reduceDir[idir] == 1)
            {
              index += 1;
              adjCellLo(tempDestGrids[index + 1], 
                        tempDestGrids[index    ], 
                        idir                   ,
                        -1                     );
          
         
            }
        }
  
      //  DisjointBoxLayout newDestGrids(tempDestGrids[numReduceDir]);
      DisjointBoxLayout newDestGrids(tempDestGrids[index + 1]);
  
      pout()<<"newDestGrids.physDomain() =" <<newDestGrids.physDomain() << endl;
      pout()<<endl;
  
      // Make a level data to test the new method
      IntVect newDestGhosts = IntVect::Unit;

      // no ghost cell in the reduction directions
      IntVect reductionVect = IntVect::Zero;
      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          reductionVect += BASISV(idir)*reduceDir[idir];
        }
   
      // vector of all ones minus the vector which is one in the reduction directions
      newDestGhosts -= reductionVect;

      //temporary
      newDestGhosts = IntVect::Zero;
      LevelData<FluxBox> newDestLD(newDestGrids , 
                                   ld.nComp()   ,  
                                   newDestGhosts);

      // Initialize the reduction destination for the new method
      for (DataIterator dit(newDestGrids); dit.ok(); ++dit) 
        {
          newDestLD[dit].setVal(0.0);
          pout()<<"newDestGrids[dit] = "<<newDestGrids[dit] <<endl;
          pout()<<endl;
        }

      pout()<<"calling constructor for new Reduction Copier";
      const ReductionCopier newReduceCopier(usingFaceSumOp      ,  
                                            ld.getBoxes()       ,
                                            newDestLD.getBoxes(),
                                            oldDbl.physDomain   (),
                                            newDestGhosts        ,
                                            reduceDir           ,
                                            exchange            );

      pout() <<"finished calling constructor for ReductionCopier"<<endl;
 
      FaceSumOp op(reduceDir);

      pout()<<"Starting the copy"<<endl;
      ld.copyTo( ld.interval()       ,
                 newDestLD           ,
                 newDestLD.interval(),
                 newReduceCopier     ,
                 op );
      pout()<<"Ending the copy"<<endl;
   
      for (DataIterator dit(newDestGrids); dit.ok(); ++dit) 
        {
          for (int dir=0; dir<SpaceDim; ++dir) 
            {
              pout()<<endl;
              pout()<<"newDestLD[dit][dir].box() = "<<newDestLD[dit][dir].box()<<endl;
              pout() << "direction = " << dir <<endl;
              pout()<<endl;
              for (BoxIterator bit(newDestLD[dit][dir].box()); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  pout()<<"newDestLD[dit][dir]("<<iv<<",0) = "<< newDestLD[dit][dir](iv,0) << endl;
                }
            }
        }
    }
#ifdef CH_MPI
      MPI_Finalize();
#endif
    
  return 0;
}

