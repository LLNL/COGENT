#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LoadBalance.H"
#include "EBArith.H"

#include "EBViscousTensorOp.H"
#include "EBArith.H"

#include "CH_Timer.H"
#include "EBViscousTensorOpFactory.H"
#include "EBCoarseAverage.H"
#include "CoarsenStuff.H"
#include "NamespaceHeader.H"

/****/
EBViscousTensorOpFactory::~EBViscousTensorOpFactory()
{
}
/****/
int
EBViscousTensorOpFactory::
refToFiner(const ProblemDomain& a_domain) const
{
  int retval = -1;
  bool found = false;
  for (int ilev = 0; ilev < m_eblgs.size(); ilev++)
    {
      if (m_eblgs[ilev].getDomain() == a_domain)
        {
          retval = m_refRatio[ilev];
          found = true;
        }
    }
  if (!found)
    {
      MayDay::Error("Domain not found in AMR hierarchy");
    }
  return retval;
}

//-----------------------------------------------------------------------
EBViscousTensorOpFactory::
EBViscousTensorOpFactory(const Vector<EBLevelGrid>&                                  a_eblgs,
                         const Real&                                                 a_alpha,
                         const Real&                                                 a_beta,
                         const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
                         const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_eta,
                         const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_lambda,
                         const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_etaIrreg,
                         const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_lambdaIrreg,
                         const Real&                                                 a_dxCoarse,
                         const Vector<int>&                                          a_refRatio,
                         const RefCountedPtr<BaseDomainBCFactory>&                   a_domainBCFactory,
                         const RefCountedPtr<BaseEBBCFactory>    &                   a_ebBCFactory,
                         const IntVect&                                              a_ghostCellsPhi,
                         const IntVect&                                              a_ghostCellsRhs,
                         int a_numLevels,
                         bool a_noMG)
{
  CH_assert(a_eblgs.size() <= a_refRatio.size());
  m_noMG = a_noMG;
  m_ghostCellsRhs = a_ghostCellsRhs;
  m_ghostCellsPhi = a_ghostCellsPhi;
  m_lambda        = a_lambda;
  m_acoef         = a_acoef;
  m_eta           = a_eta;
  m_lambdaIrreg   = a_lambdaIrreg;
  m_etaIrreg      = a_etaIrreg;
  m_alpha         = a_alpha;
  m_beta          = a_beta;
  if (a_numLevels > 0)
  {
    m_numLevels = a_numLevels;
  }
  else
  {
    m_numLevels = a_eblgs.size();
  }

  m_domainBCFactory = a_domainBCFactory;
  m_ebBCFactory     = a_ebBCFactory;

  m_eblgs           = a_eblgs;
  m_refRatio        = a_refRatio;
  m_dx.resize(m_numLevels);

  m_dx[0] = a_dxCoarse;
  for (int ilev = 1; ilev < m_numLevels; ilev++)
  {
    m_dx[ilev] = m_dx[ilev-1];
    m_dx[ilev] /= m_refRatio[ilev-1];
  }
  m_alpha = a_alpha;
  m_beta = a_beta;

  m_eblgsMG      .resize(m_numLevels);
  m_acoefMG      .resize(m_numLevels);
  m_etaMG        .resize(m_numLevels);
  m_lambdaMG     .resize(m_numLevels);
  m_etaIrregMG   .resize(m_numLevels);
  m_lambdaIrregMG.resize(m_numLevels);
  m_hasMGObjects .resize(m_numLevels);
  for (int ilev = 0; ilev < m_numLevels; ilev++)
  {
    if ((ilev==0) || (m_refRatio[ilev] > 2))
    {
      m_hasMGObjects[ilev] = true;

      int mgRef = 2;
      m_eblgsMG      [ilev].resize(0);
      m_acoefMG      [ilev].resize(0);
      m_etaMG        [ilev].resize(0);
      m_lambdaMG     [ilev].resize(0);
      m_etaIrregMG   [ilev].resize(0);
      m_lambdaIrregMG[ilev].resize(0);

      m_eblgsMG      [ilev].push_back(m_eblgs      [ilev]);
      m_acoefMG      [ilev].push_back(m_acoef      [ilev]);
      m_etaMG        [ilev].push_back(m_eta        [ilev]);
      m_lambdaMG     [ilev].push_back(m_lambda     [ilev]);
      m_etaIrregMG   [ilev].push_back(m_etaIrreg   [ilev]);
      m_lambdaIrregMG[ilev].push_back(m_lambdaIrreg[ilev]);

      bool hasCoarser = true;
      while (hasCoarser)
      {
        int imgsize = m_eblgsMG[ilev].size();
        const EBLevelGrid& eblgFine=  m_eblgsMG[ilev][imgsize-1];
        DisjointBoxLayout dblCoarMG;
        ProblemDomain  domainCoarMG;
        int maxBoxSize = 32;
        bool dumbool;
        int testRef = 2;
        hasCoarser = EBAMRPoissonOp::getCoarserLayouts(dblCoarMG,
                                                       domainCoarMG,
                                                       eblgFine.getDBL(),
                                                       eblgFine.getEBISL(),
                                                       eblgFine.getDomain(),
                                                       mgRef,
                                                       eblgFine.getEBIS(),
                                                       maxBoxSize, dumbool, testRef);

        if (hasCoarser)
        {
          m_eblgsMG[ilev].push_back(EBLevelGrid(dblCoarMG, domainCoarMG, 4, eblgFine.getEBIS()));
          int img = m_eblgsMG[ilev].size() - 1;
          const EBLevelGrid& eblgCoar = m_eblgsMG[ilev][img  ];
          const EBLevelGrid& eblgFine = m_eblgsMG[ilev][img-1];

          LayoutData<IntVectSet> irregSets(eblgCoar.getDBL());
          int nghost = 1;
          for (DataIterator dit = eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
          {
            Box grownBox = grow(eblgCoar.getDBL().get(dit()), nghost);
            grownBox &= domainCoarMG;
            irregSets[dit()] = eblgCoar.getEBISL()[dit()].getIrregIVS(grownBox);
          }
          EBFluxFactory       ebfluxfact(eblgCoar.getEBISL());
          EBCellFactory       ebcellfact(eblgCoar.getEBISL());
          BaseIVFactory<Real> baseivfact(eblgCoar.getEBISL(), irregSets);

          RefCountedPtr<LevelData<BaseIVFAB<Real> > >    etaIrregCoar( new LevelData<BaseIVFAB<Real> >(eblgCoar.getDBL(), 1, nghost*IntVect::Unit, baseivfact) );
          RefCountedPtr<LevelData<BaseIVFAB<Real> > > lambdaIrregCoar( new LevelData<BaseIVFAB<Real> >(eblgCoar.getDBL(), 1, nghost*IntVect::Unit, baseivfact) );
          RefCountedPtr<LevelData<EBCellFAB> >              acoefCoar( new LevelData<EBCellFAB>       (eblgCoar.getDBL(), 1, nghost*IntVect::Unit, ebcellfact) );
          RefCountedPtr<LevelData<EBFluxFAB> >                etaCoar( new LevelData<EBFluxFAB>       (eblgCoar.getDBL(), 1, nghost*IntVect::Unit, ebfluxfact) );
          RefCountedPtr<LevelData<EBFluxFAB> >             lambdaCoar( new LevelData<EBFluxFAB>       (eblgCoar.getDBL(), 1, nghost*IntVect::Unit, ebfluxfact) );

          coarsenStuff(*acoefCoar, *etaCoar, *lambdaCoar, *etaIrregCoar, *lambdaIrregCoar, eblgFine, eblgCoar,
              *m_acoefMG[ilev][img-1], *m_etaMG[ilev][img-1], *m_lambdaMG[ilev][img-1],
              *m_etaIrregMG[ilev][img-1], *m_lambdaIrregMG[ilev][img-1], mgRef);

          m_acoefMG      [ilev].push_back(      acoefCoar);
          m_etaMG        [ilev].push_back(        etaCoar);
          m_lambdaMG     [ilev].push_back(     lambdaCoar);
          m_etaIrregMG   [ilev].push_back(   etaIrregCoar);
          m_lambdaIrregMG[ilev].push_back(lambdaIrregCoar);
        }
      }

    }
    else
    {
      m_hasMGObjects[ilev] = false;
    }
  }
}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
EBViscousTensorOp*
EBViscousTensorOpFactory::
MGnewOp(const ProblemDomain& a_domainFine,
        int                  a_depth,
        bool                 a_homoOnly)
{
  //debug turn off multigrid
  //if (a_depth > 2)
  //  return NULL;
  //end debug
  //find out if there is a real starting point here.
  int ref=-1;
  bool found = false;
  //to turn off multigrid all together
  if (m_noMG && (a_depth > 0))
    {
      return NULL;
    }

  RefCountedPtr<LevelData<BaseIVFAB<Real> > >     etaIrreg;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  lambdaIrreg;
  RefCountedPtr<LevelData<EBCellFAB> >               acoef;
  RefCountedPtr<LevelData<EBFluxFAB> >                 eta;
  RefCountedPtr<LevelData<EBFluxFAB> >              lambda;

  for (int ilev = 0; ilev < m_numLevels && !found; ilev++)
  {
    if (a_domainFine == m_eblgs[ilev].getDomain())
    {
      found = true;
      ref = ilev ;
    }
  }
  if (!found)
  {
    MayDay::Error("No corresponding AMRLevel to starting point of MGnewOp");
  }

  //multigrid operator.  coarn by two from depth  no
  //coarse or finer to worry about.
  Real          dxMGLevel;
  EBLevelGrid eblgMGLevel;
  EBLevelGrid eblgCoarMG;
  Real dxCoar = 1.0;
  dxCoar *= -1.0;
  int refToDepth = 1;
  if (ref > 0)
  {
    dxCoar = m_dx[ref-1];
  }
  bool hasCoarMGObjects = false;
  if (a_depth == 0)
  {
    eblgMGLevel =       m_eblgs[ref];
    dxMGLevel   =          m_dx[ref];
    etaIrreg    =    m_etaIrreg[ref];
    lambdaIrreg = m_lambdaIrreg[ref];
    acoef       =       m_acoef[ref];
    eta         =         m_eta[ref];
    lambda      =      m_lambda[ref];

    hasCoarMGObjects = m_hasMGObjects[ref];
    if (hasCoarMGObjects)
    {
      eblgCoarMG = m_eblgsMG[ref][1];
    }
  }
  else
  {
    int icoar = 1;
    for (int idep = 0; idep < a_depth; idep++)
    {
      icoar *= 2;
    }
    refToDepth = icoar;
    const ProblemDomain domainFine = m_eblgs[ref].getDomain();
    ProblemDomain domainBoxMGLevel = coarsen(domainFine, icoar);
    bool foundMGLevel = false;
    int numMGLevels = m_eblgsMG[ref].size();
    for (int img = 0; img < numMGLevels; img++)
    {
      if (m_eblgsMG[ref][img].getDomain() == domainBoxMGLevel)
      {
        eblgMGLevel =       m_eblgsMG[ref][img];
        etaIrreg    =    m_etaIrregMG[ref][img];
        lambdaIrreg = m_lambdaIrregMG[ref][img];
        acoef     =       m_acoefMG[ref][img];
        eta       =         m_etaMG[ref][img];
        lambda    =      m_lambdaMG[ref][img];

        foundMGLevel = true;

        hasCoarMGObjects = ((img+1) < (numMGLevels));
        if (hasCoarMGObjects)
        {
          eblgCoarMG = m_eblgsMG[ref][img+1];
        }
        break;
      }
    }
    bool coarsenable = foundMGLevel;

    dxMGLevel = m_dx[ref];
    dxMGLevel *= Real(icoar);

    if (!coarsenable)
    {
      //not coarsenable.
      //return null
      return NULL;
    }
  }

  ViscousBaseEBBC*     viscEBBC  = (ViscousBaseEBBC*)         m_ebBCFactory->create(eblgMGLevel.getDomain(), eblgMGLevel.getEBISL(),
      dxMGLevel*RealVect::Unit,
      &m_ghostCellsPhi, &m_ghostCellsRhs);
  ViscousBaseDomainBC* viscDomBC = (ViscousBaseDomainBC*) m_domainBCFactory->create(eblgMGLevel.getDomain(), eblgMGLevel.getEBISL(),
      dxMGLevel*RealVect::Unit);
  RefCountedPtr<ViscousBaseEBBC>      ebbc(viscEBBC);
  RefCountedPtr<ViscousBaseDomainBC> dombc(viscDomBC);

  //no need for finer or coarser amr levels here
  int bogRef = 2;
  EBViscousTensorOp* newOp = NULL;
  // Time-independent a coefficient.
  newOp = new EBViscousTensorOp(EBLevelGrid(), eblgMGLevel, EBLevelGrid(), eblgCoarMG,
                                m_alpha, m_beta, acoef, eta, lambda, etaIrreg, lambdaIrreg,
                                dxMGLevel, dxCoar,  bogRef, bogRef,  dombc, ebbc,
                                hasCoarMGObjects, m_ghostCellsPhi, m_ghostCellsRhs);

  // We're finished.
  return newOp;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
EBViscousTensorOp*
EBViscousTensorOpFactory::
AMRnewOp(const ProblemDomain& a_domainFine)
{
  //figure out which level we are at.
  int ref=-1;
  bool found = false;
  EBLevelGrid eblgFine, eblgCoar;

  for (int ilev = 0; ilev < m_numLevels && !found; ilev++)
    {
      if (a_domainFine == m_eblgs[ilev].getDomain())
        {
          found = true;
          ref = ilev ;
        }
    }
  if (!found)
    {
      MayDay::Error("No corresponding AMRLevel to starting point of AMRnewOp");
    }

  int refToFiner   = 2;
  int refToCoarser = 2;
  Real dxCoar = -1.0;
  if (ref > 0)
    {
      dxCoar = m_dx[ref-1];
      eblgCoar = m_eblgs[ref-1];
      refToCoarser= m_refRatio[ref-1];
    }
  if (ref < m_numLevels-1)
    {
      eblgFine = m_eblgs[ref+1];
      refToFiner = m_refRatio[ref];
    }
  //creates coarse and finer info and bcs and all that
  EBLevelGrid      eblgMGLevel = m_eblgs[ref];
  Real               dxMGLevel =    m_dx[ref];

  bool hasCoarMGObjects = m_hasMGObjects[ref];
  EBLevelGrid eblgCoarMG;
  if (hasCoarMGObjects)
     {
      eblgCoarMG = m_eblgsMG[ref][1];
    }
  ViscousBaseEBBC*     viscEBBC  = (ViscousBaseEBBC*)     m_ebBCFactory->create(    m_eblgs[ref].getDomain(),
                                                                                    m_eblgs[ref].getEBISL(), dxMGLevel*RealVect::Unit,
                                                                                    &m_ghostCellsPhi, &m_ghostCellsRhs);
  ViscousBaseDomainBC* viscDomBC = (ViscousBaseDomainBC*) m_domainBCFactory->create(m_eblgs[ref].getDomain(),
                                                                     m_eblgs[ref].getEBISL(), dxMGLevel*RealVect::Unit);
  RefCountedPtr<ViscousBaseEBBC>      ebbc(viscEBBC);
  RefCountedPtr<ViscousBaseDomainBC> dombc(viscDomBC);

  EBViscousTensorOp* newOp = NULL;
  newOp = new EBViscousTensorOp(eblgFine, eblgMGLevel, eblgCoar, eblgCoarMG,
                                m_alpha, m_beta, m_acoef[ref], m_eta[ref], m_lambda[ref],
                                m_etaIrreg[ref], m_lambdaIrreg[ref], dxMGLevel, dxCoar,
                                refToFiner, refToCoarser,  dombc, ebbc, hasCoarMGObjects,
                                m_ghostCellsPhi, m_ghostCellsRhs);

  return newOp;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBViscousTensorOpFactory::
reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim)
{
  delete a_reclaim;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBViscousTensorOpFactory::
AMRreclaim(EBViscousTensorOp* a_reclaim)
{
  delete a_reclaim;
}
//-----------------------------------------------------------------------
#include "NamespaceFooter.H"
