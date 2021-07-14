#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCOARSENSTUFF_H_
#define _EBCOARSENSTUFF_H_

#include "LevelData.H"
#include "EBCellFAB.H"
#include "NamespaceHeader.H"

void
coarsenStuff(LevelData<EBCellFAB>               & a_acoefCoar,
             LevelData<EBFluxFAB>               & a_etaCoar,
             LevelData<EBFluxFAB>               & a_lambdaCoar,
             LevelData<BaseIVFAB<Real> >        & a_etaCoarIrreg,
             LevelData<BaseIVFAB<Real> >        & a_lambdaCoarIrreg,
             const EBLevelGrid                  & a_eblgFine,
             const EBLevelGrid                  & a_eblgCoar,
             const LevelData<EBCellFAB>         & a_acoefFine,
             const LevelData<EBFluxFAB>         & a_etaFine,
             const LevelData<EBFluxFAB>         & a_lambdaFine,
             const LevelData<BaseIVFAB<Real> >  & a_etaFineIrreg,
             const LevelData<BaseIVFAB<Real> >  & a_lambdaFineIrreg,
             const int                          & a_refToDepth);

#include "NamespaceFooter.H"
#endif
