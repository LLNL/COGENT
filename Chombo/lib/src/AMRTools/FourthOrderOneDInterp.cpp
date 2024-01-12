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

#include "FourthOrderOneDInterp.H"
#include "NamespaceHeader.H"

// Constructor 
FourthOrderOneDInterp::FourthOrderOneDInterp()
{
  m_defined = false;
  
}

// Destructor
FourthOrderOneDInterp::~FourthOrderOneDInterp()
{
}

//define
void FourthOrderOneDInterp::define(const int              & a_nComp          ,
                                   const IntVect          & a_changeVecDx    ,
                                   const int              & a_interpDir      ,
                                   const IntVect          & a_coarseGhostVect,
                                   const int              & a_refRatio       ,
                                   const Box              & a_coarseDomain   )
{
  m_nComp            = a_nComp          ;
  m_changeVecDx      = a_changeVecDx    ;
  m_interpDir        = a_interpDir      ;
  m_coarseGhostVect  = a_coarseGhostVect;
  m_refRatio         = a_refRatio       ;
  m_coarseDomain     = a_coarseDomain   ; 

  int order = 4;
  m_stencil.define(m_refRatio,
                   order     );
  
  m_defined = true;
}



FourthOrderPatchInterp& FourthOrderOneDInterp::patchInterp()
{
  CH_assert(m_defined);
  return m_patchInterp;
}

void FourthOrderOneDInterp::interpToFine(LevelData<FArrayBox>      & a_fine      ,
                                         const LevelData<FArrayBox>& a_coarseData,
                                         const IntVect             & a_growDir)
{
  CH_assert(m_defined);

   DataIterator dit = a_coarseData.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box coarseValidBox = a_coarseData.disjointBoxLayout()[dit];
      const Box fineValidBox   = a_fine.      disjointBoxLayout()[dit];

      //pout()<<"coarseValidBox  = "<< coarseValidBox <<endl << endl;
      
      // pout()<<"fineValidBox    = "<<fineValidBox    <<endl << endl;        

      const FArrayBox& coarseFab = a_coarseData[dit];
      FArrayBox      & fineFab   = a_fine      [dit];

      const Box coarseBox = coarseFab.box();
      const Box fineBox   = fineFab  .box();

      Box boxIter = coarseValidBox;
      boxIter.grow(a_growDir);

      //pout()<<"coarseBox        = "<< coarseBox<<endl;
      //pout()<<"fineBox          = "<< fineBox   <<endl;             
      //pout()<<"boxIter          = "<< boxIter   <<endl;
      
      // interpolate one direction at a time
      patchInterp(fineFab       ,
                  coarseFab     ,
                  boxIter  );

     
      
      //else
      //  {
      //   cout<<"Whoops"<<endl;
          //call patchInterp from FourthOrderFineInterp (ie, the non-tensor version)
    }
    
}

void FourthOrderOneDInterp::patchInterp(FArrayBox            & a_fineFab  ,
                                              const FArrayBox& a_coarseFab,
                                              const Box      & a_coarseBox)
{
    for (BoxIterator bit(a_coarseBox); bit.ok(); ++bit)
      {
        IntVect coarseIv = bit();
        m_stencil.apply(a_fineFab     ,
                        a_coarseFab   ,
                        coarseIv      ,
                        m_interpDir   ,
                        m_changeVecDx);  
      }
}

#include "NamespaceFooter.H"
