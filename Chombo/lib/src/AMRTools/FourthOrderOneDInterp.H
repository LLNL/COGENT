#ifdef CH_LANG_CC
/*
 *      _______              __1
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERONEDINTERP_H_
#define _FOURTHORDERONEDINTERP_H_

#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "FourthOrderOneDStencil.H"
#include "FourthOrderPatchInterp.H"

#include "NamespaceHeader.H"

/// Fourth-order interpolation in space

/**
 */
class FourthOrderOneDInterp
{
public:
  FourthOrderOneDInterp();
  ~FourthOrderOneDInterp();

  void define(const int              & a_nComp          ,
              const IntVect          & a_changeVecDx    ,
              const int              & a_interpDir      ,
              const IntVect          & a_coarseGhostVect,
              const int              & a_refRatio       ,
              const Box              & a_coarseDomain   );

  
  LevelData<FArrayBox>  & coarsenedFineData();
  FourthOrderPatchInterp& patchInterp      ();
  
  void interpToFine(LevelData<FArrayBox>      & a_fine  ,
                    const LevelData<FArrayBox>& a_coarse,
                    const IntVect             & a_growDir);
  
  
  
  void patchInterp(FArrayBox      & a_fineFab   ,
                   const FArrayBox& a_coarseFab ,
                   const Box      & a_coarseBox );
  //protected:
  
  //level data on the coarse level
    IntVect              m_coarseGhostVect;

  // intermediate classes
  FourthOrderPatchInterp m_patchInterp;
  FourthOrderOneDStencil m_stencil    ;

  // This vector coarsens a_fineLayout to m_layoutCoarsened
  int     m_nComp       ;
  IntVect m_changeVecDx ;
  int     m_interpDir   ;
  int     m_refRatio    ;
  bool    m_defined     ;
  Box     m_coarseDomain;

  

private:

  // Disallowed for all the usual reasons
  void operator=(const FourthOrderOneDInterp&);
  FourthOrderOneDInterp(const FourthOrderOneDInterp&);
};

#include "NamespaceFooter.H"
#endif
