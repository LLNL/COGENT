#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// #include <cstdio>

#include "FourthOrderTensorInterp.H"
#include "FourthOrderFineInterp.H"
#include "BoxIterator.H"


#include "NamespaceHeader.H"
//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
FourthOrderTensorInterp::FourthOrderTensorInterp()
{
  m_defined = false;
}




void FourthOrderTensorInterp::define(LevelData<FArrayBox>      & a_fineData        ,
                                     const int                 & a_refRatio        ,
                                     const Box                 & a_fineDomain      ,
                                     const LevelData<FArrayBox>& a_coarseData      ,
                                     const Box                 & a_coarseDomain    ,
                                     const int                 & a_coarseGhostsFill,
                                     const Interval            & a_fixedDims       )
{
  // make domain and interior levelData
  Box smallDomain = a_fineDomain;
  smallDomain.grow(-2*a_refRatio);
  
  //layout of boxes that are within 2*refRatio fine cells of the domain boundary
  DisjointBoxLayout   domainDbl    ;
  Vector<Box>         domainBoxes  ;
  Vector<int>         domainProcIDs;
  LevelData<FArrayBox>domainData   ;

  ProblemDomain probFine(a_fineDomain);
  
  ////layout of boxes that are farther than 2*refRatio fine cells of the domain boundary
  DisjointBoxLayout   interiorDbl    ;
  Vector<Box>         interiorBoxes  ;
  Vector<int>         interiorProcIDs;
  LevelData<FArrayBox>interiorData   ;
  
  //fine dbl
  DisjointBoxLayout fineDbl = a_fineData.disjointBoxLayout();
  Vector<Box> fineBoxes     = fineDbl.boxArray();
  Vector<int> fineProcIDs   = fineDbl.procIDs();
  
  for (int iBox = 0; iBox < fineBoxes.size(); ++ iBox)
    {
      Box thisBox = fineBoxes[iBox];
      thisBox     &= smallDomain;
      
      if (thisBox == fineBoxes[iBox])
        {
          interiorBoxes  .push_back(fineBoxes  [iBox]);
          interiorProcIDs.push_back(fineProcIDs[iBox]);
        }
      else
        {
          domainBoxes  .push_back(fineBoxes  [iBox]);
          domainProcIDs.push_back(fineProcIDs[iBox]);
        }
    }

  domainDbl  .define(domainBoxes  ,domainProcIDs  ,a_fineDomain );
  interiorDbl.define(interiorBoxes,interiorProcIDs,a_fineDomain );

  domainData  .define(domainDbl  ,a_fineData.nComp(),a_fineData.ghostVect());
  interiorData.define(interiorDbl,a_fineData.nComp(),a_fineData.ghostVect());

  // fine data parameters
  Interval compInterval(0,a_fineData.nComp() - 1);
  int nComp =  a_fineData.nComp();

   // Interior boxes: use the newer Tensor algorithm designed by Phil C.
   Vector < DisjointBoxLayout>dblVec;
   makeDBLVec(dblVec    ,
             interiorDbl,
             a_refRatio );

  Vector < FourthOrderOneDInterp*> fineInterpVec;
  IntVect coarseGhostVect = 2*IntVect::Unit;
  makeSemiCoarsenVec(fineInterpVec  ,
                     nComp          ,
                     a_refRatio     ,
                     coarseGhostVect,
                     a_coarseDomain );
  
  Vector <LevelData <FArrayBox>* > levelDataVec;
  defineLevelDataVec(levelDataVec  ,
                     dblVec        ,
                     nComp         ,
                     coarseGhostVect);

  levelDataVec[0] = &(interiorData);
  //the directional update is a += operation, which requires a 0 initial value
  for (int iComp = 0; iComp < nComp; ++iComp)
    {
      for (int iVec = 0; iVec< levelDataVec.size(); ++iVec)
        {
          for (DataIterator dit = (*(levelDataVec[iVec])).dataIterator(); dit.ok(); ++dit)
            {
              FArrayBox& fab = (*(levelDataVec[iVec]))[dit];
              fab.setVal(0.0);
            }
        }
    }

  LevelData<FArrayBox>* coarsenedFine = (levelDataVec[SpaceDim]);
  a_coarseData.copyTo(compInterval  ,
                      *coarsenedFine,
                      compInterval  );

  semiCoarsenInterpLoop(levelDataVec ,
                        fineInterpVec);

  // use the older algorithm implemented by Peter M. for boxes near the domain.
   FourthOrderFineInterp interpolator;
   interpolator.define(domainDbl         ,
                       nComp             ,
                       a_refRatio        ,
                       probFine          ,
                       a_coarseGhostsFill,
                       a_fixedDims       );

   interpolator.interpToFine(domainData  ,
                             a_coarseData);
   
  // copy the answers, which live on disjoint layouts, onto a_fineData
  interiorData.copyTo(compInterval, a_fineData, compInterval);
  domainData  .copyTo(compInterval, a_fineData, compInterval);
}

void FourthOrderTensorInterp::makeDBLVec(Vector < DisjointBoxLayout>& a_dblVec  ,
                                         const DisjointBoxLayout    & a_fineDbl ,
                                         const int                  & a_refRatio)
{

  a_dblVec.resize(SpaceDim + 1);
  a_dblVec[0] = a_fineDbl;
  
  //coarsen the fine dbl in the x-direction, and call it a_dblVec[1], then increment the coarsening direction and increment the vector index 
  for (int idir = 0; idir< SpaceDim; ++idir)
    {
      DisjointBoxLayout semiCoarsenDbl;
      
      //coarsening vector
      IntVect changeVecDx = IntVect::Unit;
      changeVecDx[idir] = a_refRatio;

      coarsen(semiCoarsenDbl,
              a_dblVec[idir],
              changeVecDx   );

      //a_dblVec[0] = fine; a_dblVec[1] = x-coarsen; a_dblVec[2] = x,y-coarsen; a_dblVec[2] = x,y,z-coarsen; 
      a_dblVec[idir + 1] = semiCoarsenDbl;
      
    }
  
}
void FourthOrderTensorInterp::makeSemiCoarsenVec(Vector < FourthOrderOneDInterp*>& a_fineInterpVec  ,
                                                 const int                       & a_nComp          ,
                                                 const int                       & a_refRatio       ,
                                                 const IntVect                   & a_coarseGhostVect,
                                                 const Box                       & a_coarseDomain   )
                        
{

  //a_fineInterpVec.resize(SpaceDim);
  // Vector < FourthOrderOneDInterp*> fineInterpVec;
  a_fineInterpVec.resize(SpaceDim,NULL);
 
   
  for (int idir = 0; idir< SpaceDim; ++idir)
    {
       //coarsening/refining vector
      IntVect changeVecDx = IntVect::Unit;
      
      //data for turning a 1x1 coarse-cell box into a semi-refined box
      int interpDir     = idir;
      changeVecDx[idir] = a_refRatio;

      a_fineInterpVec[idir] = new FourthOrderOneDInterp;
      a_fineInterpVec[idir]->define(a_nComp          ,
                                    changeVecDx      ,
                                    interpDir        ,
                                    a_coarseGhostVect,
                                    a_refRatio       ,
                                    a_coarseDomain   );   
    }
  
}
void FourthOrderTensorInterp::defineLevelDataVec(Vector <LevelData <FArrayBox>*> & a_levelDataVec,
                                                 const Vector< DisjointBoxLayout>& a_dblVec      ,
                                                 const int                       & a_nComp       ,
                                                 const IntVect                   & a_ghostVect   )

{
  a_levelDataVec.resize(SpaceDim + 1,NULL);

  for (int idir = 0; idir<=SpaceDim; ++idir)
    {
      a_levelDataVec[idir] = new LevelData <FArrayBox>;
      a_levelDataVec[idir]->define(a_dblVec[idir],
                                   a_nComp       ,
                                   a_ghostVect   );   
    }
}

void FourthOrderTensorInterp::semiCoarsenInterpLoop(Vector <LevelData <FArrayBox>* >      & a_levelDataVec ,
                                                    const Vector < FourthOrderOneDInterp*>& a_fineInterpVec)
{
  // idir = 0,  growDir = (2,0) or (2,2,0),
  // idir = 1, growDir  = (0,0) or (2,0,0).
  //If SpaceDim   == 3, then growDir = (0,0,0)  
   IntVect growDir = 2*IntVect::Unit;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      int finerIndex   = SpaceDim - (idir + 1);
      int coarserIndex = SpaceDim - idir      ;
      growDir -= 2*BASISV(SpaceDim - 1 - idir);

      a_fineInterpVec[SpaceDim - 1 - idir]->interpToFine(*(a_levelDataVec[finerIndex  ]),
                                                         *(a_levelDataVec[coarserIndex]),
                                                         growDir);
      
    }
  
}

 

#include "NamespaceFooter.H"
