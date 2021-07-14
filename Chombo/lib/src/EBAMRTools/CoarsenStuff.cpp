#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LevelData.H"
#include "EBCellFAB.H"

#include "EBLevelGrid.H"
#include "EBISLayout.H"
#include "EBCoarseAverage.H"
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
             const int                          & a_refToDepth)
{
  CH_assert(a_refToDepth > 0);

  Interval interv(0, 0);
  if (a_refToDepth == 1)
    {
      a_acoefFine. copyTo(interv,  a_acoefCoar, interv);
      a_etaFine.   copyTo(interv,    a_etaCoar, interv);
      a_lambdaFine.copyTo(interv, a_lambdaCoar, interv);

      a_etaFineIrreg.   copyTo(interv,    a_etaCoarIrreg, interv);
      a_lambdaFineIrreg.copyTo(interv, a_lambdaCoarIrreg, interv);
    }
  else
    {
      EBCoarseAverage averageOp(a_eblgFine.getDBL(),    a_eblgCoar.getDBL(),
                                a_eblgFine.getEBISL(),  a_eblgCoar.getEBISL(),
                                a_eblgCoar.getDomain(), a_refToDepth, 1,
                                a_eblgCoar.getEBIS());

      //MayDay::Warning("might want to figure out what harmonic averaging is in this context");
      averageOp.average( a_acoefCoar     ,  a_acoefFine     , interv);
      averageOp.average(   a_etaCoar     ,    a_etaFine     , interv);
      averageOp.average(   a_etaCoarIrreg,    a_etaFineIrreg, interv);
      averageOp.average(a_lambdaCoar     , a_lambdaFine     , interv);
      averageOp.average(a_lambdaCoarIrreg, a_lambdaFineIrreg, interv);

    }
  a_acoefCoar.exchange(interv);
  a_etaCoar.exchange(interv);
  a_lambdaCoar.exchange(interv);
  a_etaCoarIrreg.exchange(interv);
  a_lambdaCoarIrreg.exchange(interv);
}

#include "NamespaceFooter.H"
