#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DataIterator.H"
#include "IntVect.H"
#include "Copier.H"
#include "MayDay.H"
#include "LayoutIterator.H"
#include "NeighborIterator.H"
#include "BoxIterator.H"
#include "SPMD.H"
#include "CH_Timer.H"
#include "memtrack.H"
#include "parstream.H"
#include <chrono>

#include <vector>
#include "NamespaceHeader.H"

using std::ostream;

Pool Copier::s_motionItemPool(sizeof(MotionItem), "Copier::MotionItem", 500);

CopierBuffer::~CopierBuffer()
{
  clear();
}

void CopierBuffer::clear()
{
  if (m_sendbuffer != NULL) freeMT(m_sendbuffer);
  if (m_recbuffer  != NULL) freeMT(m_recbuffer);
  m_sendbuffer = NULL;
  m_recbuffer  = NULL;
  m_sendcapacity = 0;
  m_reccapacity = 0;
  m_ncomps = 0;
}

Copier::Copier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               bool a_exchange,
               IntVect a_shift)
{
  const ProblemDomain& domain = a_level.physDomain();
  define(a_level, a_dest, domain, IntVect::Zero, a_exchange, a_shift);
}

Copier::Copier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               const ProblemDomain& a_domain,
               bool a_exchange,
               IntVect a_shift)
{
  define(a_level, a_dest, a_domain, IntVect::Zero, a_exchange, a_shift);
}

Copier::Copier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               const IntVect& a_destGhost,
               bool a_exchange,
               IntVect a_shift)
{
  const ProblemDomain domain = a_level.physDomain();
  define(a_level, a_dest, domain, a_destGhost, a_exchange, a_shift);
}

Copier::Copier(const DisjointBoxLayout& a_level,
               const BoxLayout& a_dest,
               const ProblemDomain& a_domain,
               const IntVect& a_destGhost,
               bool a_exchange,
               IntVect a_shift)
{
  define(a_level, a_dest, a_domain, a_destGhost, a_exchange, a_shift);
}

void Copier::clear()
{
  for (unsigned int i = 0; i < m_localMotionPlan.size(); ++i)
    {
      s_motionItemPool.returnPtr( m_localMotionPlan[i] );
    }
  for (unsigned int i = 0; i < m_fromMotionPlan.size(); ++i)
    {
      s_motionItemPool.returnPtr( m_fromMotionPlan[i] );
    }
  for (unsigned int i = 0; i < m_toMotionPlan.size(); ++i)
    {
      s_motionItemPool.returnPtr( m_toMotionPlan[i] );
    }

  m_localMotionPlan.resize(0);
  m_fromMotionPlan.resize(0);
  m_toMotionPlan.resize(0);
  m_isDefined = false;
  m_buffers.clear();
}

Copier::Copier(const Copier& a_rhs)
{
  *this = a_rhs;
}

Copier& Copier::operator= (const Copier& b)
{
  clear();

  // only do any of this if the source is defined
  if (b.isDefined() )
    {
      m_localMotionPlan.resize(b.m_localMotionPlan.size());
      for (int i = 0; i < m_localMotionPlan.size(); ++i)
        {
          m_localMotionPlan[i] = new (s_motionItemPool.getPtr()) MotionItem(*(b.m_localMotionPlan[i]));
        }
      
      m_fromMotionPlan.resize(b.m_fromMotionPlan.size());
      for (int i = 0; i < m_fromMotionPlan.size(); ++i)
        {
          m_fromMotionPlan[i] = new (s_motionItemPool.getPtr()) MotionItem(*(b.m_fromMotionPlan[i]));
        }
      
      m_toMotionPlan.resize(b.m_toMotionPlan.size());
      for (int i = 0; i < m_toMotionPlan.size(); ++i)
        {
          m_toMotionPlan[i] = new (s_motionItemPool.getPtr()) MotionItem(*(b.m_toMotionPlan[i]));
        }
      
      m_isDefined = true;
    }
  return *this;
}

void Copier::trimMotion(const DisjointBoxLayout& a_exchangedLayout, const IntVect& a_ghost,
                        const Vector<MotionItem*>& a_oldItems, Vector<MotionItem*>& a_newItems)
{
  for (int i = 0; i < a_oldItems.size(); ++i)
    {
      const MotionItem& item = *(a_oldItems[i]);
      const Box& b = a_exchangedLayout[item.toIndex];
      const Box& c = item.toRegion;
      for (int d = 0; d < CH_SPACEDIM; ++d)
        {
          Box g(b);
          g.grow(d, a_ghost[d]);
          if (g.intersectsNotEmpty(c))
            {
              MotionItem* it = new (s_motionItemPool.getPtr()) MotionItem(item);
              a_newItems.push_back(it);
              break;
            }
        }
    }
}

///
void Copier::trimEdges(const DisjointBoxLayout& a_exchangedLayout, const IntVect& a_ghost)
{
  Copier oldCopier = *this;
  clear();

  m_isDefined = oldCopier.m_isDefined;
  
  trimMotion(a_exchangedLayout, a_ghost, oldCopier.m_localMotionPlan, m_localMotionPlan);
  //   pout() << "old Copy operations:" << oldCopier.m_localMotionPlan.size() << "  "
  //         << "new Copy operations:" << m_localMotionPlan.size() << "\n";
  trimMotion(a_exchangedLayout, a_ghost, oldCopier.m_fromMotionPlan, m_fromMotionPlan);
  trimMotion(a_exchangedLayout, a_ghost, oldCopier.m_toMotionPlan, m_toMotionPlan);
}

void Copier::reverse()
{
  for (int i = 0; i < m_localMotionPlan.size(); ++i)
    {
      m_localMotionPlan[i]->reverse();
    }
  for (int i = 0; i < m_fromMotionPlan.size(); ++i)
    {
      m_fromMotionPlan[i]->reverse();
    }
  for (int i = 0; i < m_toMotionPlan.size(); ++i)
    {
      m_toMotionPlan[i]->reverse();
    }
  m_fromMotionPlan.swap(m_toMotionPlan);

}

bool Copier::operator==(const Copier& rhs) const
{
  if(m_localMotionPlan.size() != rhs.m_localMotionPlan.size()) return false;
  if(m_fromMotionPlan.size() != rhs.m_fromMotionPlan.size()) return false;
  if(m_toMotionPlan.size() != rhs.m_toMotionPlan.size()) return false;
  for(int i=0; i<m_localMotionPlan.size(); ++i)
    {
      if(!(*m_localMotionPlan[i] == *rhs.m_localMotionPlan[i]))
        {
          return false;
        }
    }
  return true;

}

void Copier::coarsen(int a_refRatio)
{
  for (int i = 0; i < m_localMotionPlan.size(); ++i)
    {
      m_localMotionPlan[i]->fromRegion.coarsen(a_refRatio);
      m_localMotionPlan[i]->toRegion.coarsen(a_refRatio);
    }
  for (int i = 0; i < m_fromMotionPlan.size(); ++i)
    {
      m_fromMotionPlan[i]->fromRegion.coarsen(a_refRatio);
      m_fromMotionPlan[i]->toRegion.coarsen(a_refRatio);
    }
  for (int i = 0; i < m_toMotionPlan.size(); ++i)
    {
      m_toMotionPlan[i]->fromRegion.coarsen(a_refRatio);
      m_toMotionPlan[i]->toRegion.coarsen(a_refRatio);
    }

}
void Copier::define(const DisjointBoxLayout& a_level,
                    const BoxLayout& a_dest,
                    bool a_exchange,
                    IntVect a_shift)
{
  const ProblemDomain domain = a_level.physDomain();
  define(a_level, a_dest, domain, IntVect::Zero, a_exchange, a_shift);
}

void Copier::define(const DisjointBoxLayout& a_level,
                    const BoxLayout& a_dest,
                    const ProblemDomain& a_domain,
                    bool a_exchange,
                    IntVect a_shift)
{
  define(a_level, a_dest, a_domain, IntVect::Zero, a_exchange, a_shift);
}

void Copier::define(const DisjointBoxLayout& a_level,
                    const BoxLayout& a_dest,
                    const IntVect& a_ghost,
                    bool a_exchange,
                    IntVect a_shift)
{
  const ProblemDomain domain = a_level.physDomain();
  define(a_level, a_dest, domain, a_ghost, a_exchange, a_shift);
}

void Copier::defineFixedBoxSize(const DisjointBoxLayout& a_src,
                                const LMap&  a_lmap,
                                const IntVect&  a_ghost,
                                const ProblemDomain& a_domain,
                                bool  a_includeSelf,
                                bool  a_reverse)
{
  CH_TIME("Copier::defineFixedBoxSize same layout");
  CH_assert(a_src.isClosed());
  auto t0 = ch_ticks();
  auto th = ch_ticks();
  auto ti = ch_ticks();
  auto ts = std::chrono::system_clock::now();
  th=0; ti=0;
  clear();
  m_isDefined=true;
  buffersAllocated = false;

  DataIterator dit = a_src.dataIterator();
  if (dit.size() == 0)
    {
      // just go ahead and return, this processor owns no boxes, and thus will not participate in this Copier operation
      return;
    }

  Box domainBox(a_domain.domainBox());
  bool isPeriodic = false;
  if (!domainBox.isEmpty())
    {
      isPeriodic = a_domain.isPeriodic();
    }

  IntVect origin = domainBox.smallEnd();
  dit.begin();
  int myprocID = procID();
  IntVect fsize = a_src[dit].size();
  IntVect hghost =(a_ghost + fsize - IntVect::Unit)/fsize; // rounding up integer division
  for (; dit.ok(); ++dit)
    {
      const Box& b = a_src[dit];
      CH_assert(b.size() == fsize); // verify every box is the same size for this function
      Box bghost(b);
      bghost.grow(a_ghost);
      IntVect loEnd = b.smallEnd();
      uint64_t hash = loEnd.hash(origin, fsize);
      loEnd/=fsize;
      CH_assert(loEnd*fsize == b.smallEnd());// verify coarsen-refine is grid aligned
 
      // verify the hash retrieves this box from the layout
      auto th0 = ch_ticks();
      auto m = a_lmap.find(hash);
      th += ch_ticks()-th0;
      auto e = a_lmap.end();
      CH_assert(m != e && dit() == m->second);
      BoxIterator bit(Box(loEnd-hghost, loEnd+hghost));
      for(bit.begin(); bit.ok(); ++bit)
      {
        IntVect i = bit()*fsize;
        if(a_domain.image(i)) // code still won't do periodic correctly yet.....(bvs)
          {
            uint64_t h = i.hash(origin, fsize);
            if(h == hash)  { if(a_includeSelf)
              {
                MotionItem* item = new (s_motionItemPool.getPtr())
                  MotionItem(dit(), dit(), bghost);
                m_localMotionPlan.push_back(item);
              }}
            else
              {
                th0=ch_ticks();
                auto index = a_lmap.find(h);
                th+=ch_ticks()-th0;
                if(index != e) // got a hit
                  {
                    Box d = a_src[index->second];
                    Box dghost(d);
                    dghost.grow(a_ghost);
                    if(b.intersectsNotEmpty(dghost))
                      {
                        int proc = a_src.procID(index->second);
                        if(myprocID == proc) // local copy operation
                          {
                            auto ti0=ch_ticks();
                            void* spot = s_motionItemPool.getPtr();
                            ti+=ch_ticks()-ti0;
                            MotionItem* item = new (spot)
                              MotionItem(dit(), DataIndex(index->second), b & dghost);
        
                            m_localMotionPlan.push_back(item);
     
                          }
                        else
                          {
                            auto ti0 = ch_ticks();
                            auto spot = s_motionItemPool.getPtr();
                            ti+= ch_ticks() - ti0;
                            MotionItem* item1 = new (spot)
                              MotionItem(dit(), DataIndex(index->second), b & dghost);
                            item1->procID=proc;
                            m_fromMotionPlan.push_back(item1);
                            ti0 = ch_ticks();
                            spot = s_motionItemPool.getPtr();
                            ti+= ch_ticks() - ti0;
                            MotionItem* item2 = new (spot)
                              MotionItem(DataIndex(index->second), dit(), d & bghost);
                            item2->procID=proc;
                            m_toMotionPlan.push_back(item2);
                          }
                      }
                  }
              }
          }
      }
    }
  if (a_reverse)
    {
      reverse();
    }
  auto t1 = ch_ticks();
  sort();
  auto t2 = ch_ticks();
  std::chrono::duration<double> diff = std::chrono::system_clock::now()-ts;
  double rate = diff.count()/(t2-t1);
  pout()<<"total:"<<(t2-t0)*rate<<" sort:"<<(t2-t1)*rate<<" map.find:"<<th*rate
        <<" pool acquire:"<<ti*rate<<"\n";
}

void Copier::defineFixedBoxSize(const DisjointBoxLayout& a_src,
                                const LMap&  a_lmapSrc,
                                const IntVect&  a_refRatioSrc,
                                const DisjointBoxLayout& a_dst,
                                const LMap&  a_lmapDst,
                                const IntVect&  a_refRatioDst,
                                const IntVect&  a_ghost,
                                const ProblemDomain& a_domain,
                                bool  a_reverse)
{
  CH_TIME("Copier::defineFixedBoxSize diff layouts");
  CH_assert(a_src.isClosed());
  CH_assert(a_dst.isClosed());
  int myprocID = procID();
  clear();
  m_isDefined=true;
  buffersAllocated = false;
  Box domainBox(a_domain.domainBox());
  bool isPeriodic = false;
  if (!domainBox.isEmpty())
    {
      isPeriodic = a_domain.isPeriodic();
    }

  DataIterator ditSrc = a_src.dataIterator();
  if (ditSrc.size() == 0)
    {
      // just go ahead and return, this processor owns no boxes, and thus will not participate in this Copier operation
      return;
    }

  // fsizeSrc and fsizeDst are the sizes of the boxes in the layouts
  // BEFORE any coarsening of either.

  ditSrc.begin();
  IntVect fsizeSrc = a_src[ditSrc].size();
  // IntVect fsizeSrcFull = fsizeSrc * a_refRatioSrc;
  IntVect fsizeSrcReal = fsizeSrc / a_refRatioSrc;
  ProblemDomain domainSrc(a_domain);
  domainSrc.refine(a_refRatioSrc);

  DataIterator ditDst = a_dst.dataIterator();
  ditDst.begin();
  IntVect fsizeDst = a_dst[ditDst].size();
  // IntVect fsizeDstFull = fsizeDst * a_refRatioDst;
  IntVect fsizeDstReal = fsizeDst / a_refRatioDst;
  ProblemDomain domainDst(a_domain);
  domainDst.refine(a_refRatioDst);

  // IntVect origin = domainBox.smallEnd();

  IntVect originSrc = domainSrc.domainBox().smallEnd();
  IntVect originDst = domainDst.domainBox().smallEnd();

  auto eSrc = a_lmapSrc.end();
  auto eDst = a_lmapDst.end();

  for (; ditSrc.ok(); ++ditSrc)
    {
      const Box& b = a_src[ditSrc];
      CH_assert(b.size() == fsizeSrc); // verify every source box is the same size for this function
      IntVect loSrc = b.smallEnd();
      { // Verify coarsen-refine is grid-aligned.
        IntVect boxIndSrc = loSrc / fsizeSrc; // index of source box
        CH_assert(boxIndSrc * fsizeSrc == loSrc);
      }
      { // verify the hash retrieves this box from the layout
        uint64_t hashSrc = loSrc.hash(originSrc, fsizeSrc);
        auto mSrc = a_lmapSrc.find(hashSrc);
        CH_assert(mSrc != eSrc);
        CH_assert(ditSrc() == mSrc->second);
      }

      Box bReal(b);
      bReal.coarsen(a_refRatioSrc);

      // the grown source box: cells of source box + ghosts,
      // coarsened to give
      // indices of dest boxes that overlap the grown source box
      Box boxIndsDstCover(bReal);
      boxIndsDstCover.grow(a_ghost);
      boxIndsDstCover.coarsen(fsizeDstReal);
      BoxIterator bit(boxIndsDstCover);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect loDst = bit()*fsizeDst; // in uncoarsened index space
          if (domainDst.image(loDst)) // code still won't do periodic correctly yet.....(bvs)
            {
              uint64_t hashDst = loDst.hash(originDst, fsizeDst);
              auto indexNbr = a_lmapDst.find(hashDst);
              if (indexNbr != eDst) // got a hit
                {
                  auto nbr = indexNbr->second;
                  Box dReal = a_dst[nbr];
                  dReal.coarsen(a_refRatioDst);
                  Box dRealghost = grow(dReal, a_ghost);
                  if (bReal.intersectsNotEmpty(dRealghost))
                    {
                      Box intersection = dRealghost & bReal;
                      unsigned int nbrProc = a_dst.procID(nbr);
                      MotionItem* item = new (s_motionItemPool.getPtr())
                        MotionItem(ditSrc(), // source on this proc
                                   DataIndex(nbr), // dest
                                   intersection); // region
                      if (myprocID == nbrProc) // local copy operation
                        {
                          m_localMotionPlan.push_back(item);
                        }
                      else
                        { // TO nbrProc FROM myprocID
                          item->procID = nbrProc;
                          // Add to list of things that this proc will send
                          m_fromMotionPlan.push_back(item);
                          // Need to add to m_toMotionPlan from nbrProc,
                          // which we'll get to in the ditDst loop.
                        }
                    }
                }
            }
        }
    }
  for (; ditDst.ok(); ++ditDst)
    {
      const Box& d = a_dst[ditDst];
      CH_assert(d.size() == fsizeDst); // verify every dest box is the same size for this function
      IntVect loDst = d.smallEnd();
      { // Verify coarsen-refine is grid-aligned.
        IntVect boxIndDst = loDst / fsizeDst; // index of source box
        CH_assert(boxIndDst * fsizeDst == loDst);
      }
      { // verify the hash retrieves this box from the layout
        uint64_t hashDst = loDst.hash(originDst, fsizeDst);
        auto mDst = a_lmapDst.find(hashDst);
        CH_assert(mDst != eDst);
        CH_assert(ditDst() == mDst->second);
      }

      Box dReal(d);
      dReal.coarsen(a_refRatioDst);
      Box dRealghost = grow(dReal, a_ghost);
      IntVect loEndDst = dReal.smallEnd();
      IntVect hiEndDst = dReal.bigEnd();

      // the grown dest box: cells of dest box + ghosts,
      // coarsened to give
      // indices of source boxes that overlap the grown dest box
      Box boxIndsSrcCover(loEndDst - a_ghost,
                          hiEndDst + a_ghost,
                          IndexType::TheCellType());
      boxIndsSrcCover.coarsen(fsizeSrcReal);
      BoxIterator bit(boxIndsSrcCover);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect loSrc = bit()*fsizeSrc; // in uncoarsened index space
          if (domainSrc.image(loSrc)) // code still won't do periodic correctly yet.....(bvs)
            {
              uint64_t hashSrc = loSrc.hash(originSrc, fsizeSrc);
              auto indexNbr = a_lmapSrc.find(hashSrc);
              if (indexNbr != eSrc) // got a hit
                {
                  auto nbr = indexNbr->second;
                  Box bReal = a_src[nbr];
                  bReal.coarsen(a_refRatioSrc);
                  // Box bRealghost = grow(bReal, a_ghost);
                  if (bReal.intersectsNotEmpty(dRealghost))
                    {
                      Box intersection = dRealghost & bReal;
                      unsigned int nbrProc = a_src.procID(nbr);
                      // Case (myprocID == nbrProc) already handled above.
                      if (myprocID != nbrProc)
                        {
                          MotionItem* item = new (s_motionItemPool.getPtr())
                            MotionItem(DataIndex(nbr), // source
                                       ditDst(), // dest on this proc
                                       intersection); // region
                          // FROM nbrProc TO myprocID
                          item->procID = nbrProc;
                          // Add to list of things that this proc will receive
                          m_toMotionPlan.push_back(item);
                        }
                    }
                }
            }
        }
    }
  if (a_reverse)
    {
      reverse();
    }
  sort();
}

void Copier::defineFixedSizeNodesCollect(const DisjointBoxLayout& a_layout,
                                         const LMap&  a_lmap,
                                         const IntVect& a_ghostSrc,
                                         const IntVect& a_ghostDst,
                                         const ProblemDomain& a_domain)
{
  CH_TIME("Copier::defineFixedSizeNodesCollect one layout");
  CH_assert(a_layout.isClosed());

  // This does not yet work with periodic domain.
  CH_assert(!a_domain.isPeriodic());
  // Maybe a flag later if you want a Copier that doesn't include self.
  bool includeSelf = true;

  clear();
  m_isDefined = true;
  buffersAllocated = false; // mutable member data of Copier

  DataIterator dit = a_layout.dataIterator();
  if (dit.size() == 0)
    {
      // just go ahead and return: this processor owns no boxes,
      // and thus will not participate in this Copier operation
      return;
    }

  Box domainBox(a_domain.domainBox());
  IntVect origin = domainBox.smallEnd();
  dit.begin();
  unsigned int thisProc = procID();
  // fixedBoxSize = dimensions of first box; should be same for all boxes
  IntVect fixedBoxSize = a_layout[dit].size();
  IntVect ghostMax = max(a_ghostSrc, a_ghostDst);
  IntVect ghostMin = min(a_ghostSrc, a_ghostDst);
  // length of ghost vector PLUS 1 -- so that you get all NODEs --
  // in units of fixedBoxSize, rounding up
  IntVect hghost = (ghostMax + fixedBoxSize) / fixedBoxSize;
  auto indexEnd = a_lmap.end();
  for (; dit.ok(); ++dit)
    {
      const Box& baseBox = a_layout[dit];
      // Verify every box is the same size.
      CH_assert(baseBox.size() == fixedBoxSize);
      IntVect baseBoxLo = baseBox.smallEnd();
      // Result of hash function on low corner of box
      // in the layout with the given origin and blocking factor.
      auto hashBase = baseBoxLo.hash(origin, fixedBoxSize);

      baseBoxLo /= fixedBoxSize;
      // From here on, baseBoxLo is in units of fixedBoxSize.
      // Verify coarsen-refine is grid-aligned.
      CH_assert(baseBoxLo*fixedBoxSize == baseBox.smallEnd());

      {
        // Verify the hash retrieves this box from the layout.
        auto indexBase = a_lmap.find(hashBase);
        CH_assert(indexBase != indexEnd);
        CH_assert(dit() == indexBase->second);
      }

      Box thisBaseNodes = surroundingNodes(baseBox);
      Box thisGhostSrcNodes = grow(thisBaseNodes, a_ghostSrc);
      Box thisGhostDstNodes = grow(thisBaseNodes, a_ghostDst);
      Box thisGhostMinNodes = grow(thisBaseNodes, ghostMin);

      BoxIterator bit(Box(baseBoxLo - hghost, baseBoxLo + hghost));
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect nbrLo = bit()*fixedBoxSize;
          if (a_domain.image(nbrLo)) // nbrLo set to periodic image in a_domain
            {
              auto hashNbr = nbrLo.hash(origin, fixedBoxSize);
              if (hashNbr == hashBase)
                {
                  if (includeSelf)
                    {
                      MotionItem* item = new (s_motionItemPool.getPtr())
                        MotionItem(dit(), // source
                                   dit(), // dest
                                   thisGhostMinNodes); // region
                      m_localMotionPlan.push_back(item);
                    }
                }
              else
                { // Not the same patch.
                  auto indexNbr = a_lmap.find(hashNbr);
                  if (indexNbr != indexEnd) // got a hit
                    {
                      auto nbr = indexNbr->second;
                      Box nbrBox = a_layout[nbr];
                      Box nbrBaseNodes = surroundingNodes(nbrBox);
                      Box nbrGhostSrcNodes = grow(nbrBaseNodes, a_ghostSrc);
                      if (nbrGhostSrcNodes.intersectsNotEmpty(thisGhostDstNodes))
                        { // Communicate from nbr patch to this patch.
                          Box intersectNodes = nbrGhostSrcNodes & thisGhostDstNodes;
                          unsigned int nbrProc = a_layout.procID(nbr);
                          // Copy intersectNodes from nbr patch to this patch.
                          MotionItem* item = new (s_motionItemPool.getPtr())
                            MotionItem(DataIndex(nbr), // source
                                       dit(), // dest
                                       intersectNodes); // region
                          if (thisProc == nbrProc) // local copy operation
                            {
                              m_localMotionPlan.push_back(item);
                            }
                          else
                            { // thisProc receives from nbrProc 
                              item->procID = nbrProc;
                              m_toMotionPlan.push_back(item);
                            }
                        }
                      Box nbrGhostDstNodes = grow(nbrBaseNodes, a_ghostDst);
                      if (thisGhostSrcNodes.intersectsNotEmpty(nbrGhostDstNodes))
                        { // Communicate from this patch to nbr patch.
                          Box intersectNodes = thisGhostSrcNodes & nbrGhostDstNodes;
                          unsigned int nbrProc = a_layout.procID(nbr);
                          if (thisProc != nbrProc)
                            { // thisProc sends to nbrProc 
                              MotionItem* item = new (s_motionItemPool.getPtr())
                                MotionItem(dit(), // source
                                           DataIndex(nbr), // dest
                                           intersectNodes); // region
                              item->procID = nbrProc;
                              m_fromMotionPlan.push_back(item);
                            }
                        }
                    }
                }
            }
        }
    }
  sort();
}


void Copier::defineFixedSizeNodesCollect(const DisjointBoxLayout& a_src,
                                         const LMap&  a_lmapSrc,
                                         const IntVect&  a_refRatioSrc,
                                         const IntVect&  a_ghostSrc,
                                         const DisjointBoxLayout& a_dst,
                                         const LMap&  a_lmapDst,
                                         const IntVect&  a_refRatioDst,
                                         const IntVect&  a_ghostDst,
                                         const ProblemDomain& a_domain)
{
  CH_TIME("Copier::defineFixedSizeNodesCollect diff layouts");
  CH_assert(a_src.isClosed());
  CH_assert(a_dst.isClosed());

  // This does not yet work with periodic domain.
  CH_assert(!a_domain.isPeriodic());

  int myprocID = procID();
  clear();
  m_isDefined = true;
  buffersAllocated = false;
  
  DataIterator ditSrc = a_src.dataIterator();
  if (ditSrc.size() == 0)
    {
      // just go ahead and return, this processor owns no boxes, and thus will not participate in this Copier operation
      return;
    }

  Box domainBox(a_domain.domainBox());

  // fsizeSrc and fsizeDst are the sizes of the boxes in the layouts
  // BEFORE any coarsening of either.

  ditSrc.begin();
  IntVect fsizeSrc = a_src[ditSrc].size();
  // IntVect fsizeSrcFull = fsizeSrc * a_refRatioSrc;
  IntVect fsizeSrcReal = fsizeSrc / a_refRatioSrc;
  ProblemDomain domainSrc(a_domain);
  domainSrc.refine(a_refRatioSrc);

  DataIterator ditDst = a_dst.dataIterator();
  ditDst.begin();
  IntVect fsizeDst = a_dst[ditDst].size();
  // IntVect fsizeDstFull = fsizeDst * a_refRatioDst;
  IntVect fsizeDstReal = fsizeDst / a_refRatioDst;
  ProblemDomain domainDst(a_domain);
  domainDst.refine(a_refRatioDst);

  // IntVect origin = domainBox.smallEnd();

  IntVect originSrc = domainSrc.domainBox().smallEnd();
  IntVect originDst = domainDst.domainBox().smallEnd();

  auto eSrc = a_lmapSrc.end();
  auto eDst = a_lmapDst.end();

  IntVect ghostMax = max(a_ghostSrc, a_ghostDst);
  // PLUS 1 so that you get all NODEs
  IntVect ghost1 = ghostMax + IntVect::Unit;

  for (; ditSrc.ok(); ++ditSrc)
    {
      const Box& b = a_src[ditSrc];
      CH_assert(b.size() == fsizeSrc); // verify every source box is the same size for this function
      IntVect loSrc = b.smallEnd();
      { // Verify coarsen-refine is grid-aligned.
        IntVect boxIndSrc = loSrc / fsizeSrc; // index of source box
        CH_assert(boxIndSrc * fsizeSrc == loSrc);
      }
      { // verify the hash retrieves this box from the layout
        uint64_t hashSrc = loSrc.hash(originSrc, fsizeSrc);
        auto mSrc = a_lmapSrc.find(hashSrc);
        CH_assert(mSrc != eSrc);
        CH_assert(ditSrc() == mSrc->second);
      }

      Box bRealCells(b);
      bRealCells.coarsen(a_refRatioSrc);
      Box bRealNodes = surroundingNodes(bRealCells);
      Box bRealNodesGhost = grow(bRealNodes, a_ghostSrc);
      
      // the grown source box: cells of source box + ghosts,
      // coarsened to give
      // indices of dest boxes that overlap the grown source box.
      Box boxIndsDstCover(bRealCells);
      boxIndsDstCover.grow(ghost1);
      boxIndsDstCover.coarsen(fsizeDstReal);
      BoxIterator bit(boxIndsDstCover);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect loDst = bit()*fsizeDst; // in uncoarsened index space
          uint64_t hashDst = loDst.hash(originDst, fsizeDst);
          auto indexNbr = a_lmapDst.find(hashDst);
          if (indexNbr != eDst) // got a hit
            {
              auto nbr = indexNbr->second;
              Box dRealCells = a_dst[nbr];
              dRealCells.coarsen(a_refRatioDst);
              Box dRealNodes = surroundingNodes(dRealCells);
              Box dRealNodesGhost = grow(dRealNodes, a_ghostDst);
              if (bRealNodesGhost.intersectsNotEmpty(dRealNodesGhost))
                {
                  Box intersection = dRealNodesGhost & bRealNodesGhost;
                  unsigned int nbrProc = a_dst.procID(nbr);
                  MotionItem* item = new (s_motionItemPool.getPtr())
                    MotionItem(ditSrc(), // source on this proc
                               DataIndex(nbr), // dest
                               intersection); // region
                  if (myprocID == nbrProc) // local copy operation
                    {
                      m_localMotionPlan.push_back(item);
                    }
                  else
                    { // TO nbrProc FROM myprocID
                      item->procID = nbrProc;
                      // Add to list of things that this proc will send
                      m_fromMotionPlan.push_back(item);
                      // Need to add to m_toMotionPlan from nbrProc,
                      // which we'll get to in the ditDst loop.
                    }
                }
            }
        }
    }
  for (; ditDst.ok(); ++ditDst)
    {
      const Box& d = a_dst[ditDst];
      CH_assert(d.size() == fsizeDst); // verify every dest box is the same size for this function
      IntVect loDst = d.smallEnd();
      { // Verify coarsen-refine is grid-aligned.
        IntVect boxIndDst = loDst / fsizeDst; // index of source box
        CH_assert(boxIndDst * fsizeDst == loDst);
      }
      { // verify the hash retrieves this box from the layout
        uint64_t hashDst = loDst.hash(originDst, fsizeDst);
        auto mDst = a_lmapDst.find(hashDst);
        CH_assert(mDst != eDst);
        CH_assert(ditDst() == mDst->second);
      }

      Box dRealCells(d);
      dRealCells.coarsen(a_refRatioDst);
      Box dRealNodes = surroundingNodes(dRealCells);
      Box dRealNodesGhost = grow(dRealNodes, a_ghostDst);

      // the grown dest box: cells of dest box + ghosts,
      // coarsened to give
      // indices of source boxes that overlap the grown dest box
      Box boxIndsSrcCover(dRealCells);
      boxIndsSrcCover.grow(ghost1);
      boxIndsSrcCover.coarsen(fsizeSrcReal);
      BoxIterator bit(boxIndsSrcCover);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect loSrc = bit()*fsizeSrc; // in uncoarsened index space
          uint64_t hashSrc = loSrc.hash(originSrc, fsizeSrc);
          auto indexNbr = a_lmapSrc.find(hashSrc);
          if (indexNbr != eSrc) // got a hit
            {
              auto nbr = indexNbr->second;
              Box bRealCells = a_src[nbr];
              bRealCells.coarsen(a_refRatioSrc);
              Box bRealNodes = surroundingNodes(bRealCells);
              Box bRealNodesGhost = grow(bRealNodes, a_ghostSrc);
              if (bRealNodesGhost.intersectsNotEmpty(dRealNodesGhost))
                {
                  Box intersection = dRealNodesGhost & bRealNodesGhost;
                  unsigned int nbrProc = a_src.procID(nbr);
                  // Case (myprocID == nbrProc) already handled above.
                  if (myprocID != nbrProc)
                    {
                      MotionItem* item = new (s_motionItemPool.getPtr())
                        MotionItem(DataIndex(nbr), // source
                                   ditDst(), // dest on this proc
                                   intersection); // region
                      // FROM nbrProc TO myprocID
                      item->procID = nbrProc;
                      // Add to list of things that this proc will receive
                      m_toMotionPlan.push_back(item);
                    }
                }
            }
        }
    }
  sort();
}


void Copier::define(const BoxLayout& a_level,
                    const BoxLayout& a_dest,
                    const ProblemDomain& a_domain,
                    const IntVect& a_ghost,
                    bool a_exchange,
                    IntVect a_shift)
{
  CH_TIME("Copier::define");
  CH_assert(a_level.isClosed());
  CH_assert(a_dest.isClosed());
  //  CH_assert(a_level.checkPeriodic(a_domain));

  clear();
  m_isDefined = true;
  buffersAllocated = false;
  //bool self = a_dest == a_level;
  const BoxLayout& level = a_level;
  const BoxLayout& dest = a_dest;

  // set up vector of dataIndexes to keep track of which
  // "to" boxes are not completely contained within the primary
  // domain.  these boxes are then candidates for filling by
  // periodic images of the "from" data.
  Vector<DataIndex> periodicallyFilledToVect;

  // in order to cull which "from" data may be needed to
  // fill the "to" data, keep track of the radius around the
  // primary domain in which all these cells lie.
  // do this by incrementally growing the domain box and
  // keeping track of what this radius is.
  // just to make things simpler, start off with a radius of one
  Box grownDomainCheckBox = a_domain.domainBox();
  grownDomainCheckBox.grow(1);
  int periodicCheckRadius = 1;
  // keep track of whether any boxes are in the negative domain, since they will require an extra shift
  // due to integer math
  CH_assert(SpaceDim <= 6);
  bool negativeStretch[6] = {false};

  // since valid regions of the "from" DBL may also be outside
  // the primary domain, need to keep track of whether any of these
  // need to be checked separately.
  Vector<DataIndex> periodicFromVect;
  // use same domain trick here as well
  Box grownFromDomainCheckBox = a_domain.domainBox();
  int periodicFromCheckRadius = 1;

  Box domainBox(a_domain.domainBox());
  bool isPeriodic = false;
  if (!domainBox.isEmpty())
    isPeriodic = a_domain.isPeriodic();

  // (dfm -- 9/13/05) as currently written, the Copier won't correctly
  // handle periodic cases where the number of ghost cells is greater
  // than the width of the domain.  We _should_ do multiple wraparounds,
  // but we don't. So, put in this assertion. We can revisit this if it
  // becomes an issue
  bool multipleWraps = false;
  if (isPeriodic)
    {
      for (int dir = 0; dir < SpaceDim; dir++)
        {
          if (a_domain.isPeriodic(dir))
            {
              // first, check to see if ghost cells require multiple wraps
              if (a_ghost[dir] > domainBox.size(dir))
                {
                  MayDay::Warning("nGhost > domainBox size in periodic direction - requires multiple wraps");
                  multipleWraps = true;
                }
              // then, check to see if valid regions require multiple shifts
              // don't bother if we've already decided to do multiple wraps
              if (!multipleWraps)
                {
                  int singleWrapHi = domainBox.bigEnd(dir) +domainBox.size(dir);
                  int singleWrapLo = domainBox.smallEnd(dir) -domainBox.size(dir);

                  // for now, do this for all boxes (no need for communication
                  // for really big layouts, may want to distribute this
                  LayoutIterator lit = dest.layoutIterator();             
                  for (lit.begin(); lit.ok(); ++lit)
                    {
                      // check to see if this box requires multiple wraps
                      Box destBox = dest[lit];
                      if ((destBox.bigEnd(dir) > singleWrapHi) ||
                          (destBox.smallEnd(dir) < singleWrapLo))
                        {
                          MayDay::Warning("Copier::define -- multiple wraps required for valid-region copying");
                          multipleWraps = true;
                        } // end if box doesn't fit in single-wrapped domain
                    } // end loop over dest boxes
                } // end if we don't already have multiple wraps

              // if we switch to distrubuted testing for the last loop, 
              // would need a reduction here...
            }
        }
    }

  unsigned int myprocID = procID();

  // The following 4 for loops are the result of a performance optimization.
  // When increasing the size of the problem, we found that the code was
  // looping over every destination box for every source box which was N1*N2
  // loop iterations (essentially an N-squared approach).
  // The following code attempts to simply reduce N1 and N2 by first separating
  // the boxes (or LayoutIndexes to boxes) that reside on the current processor.
  // Then the loop to determine which boxes of the first list intersect with
  // which boxes of the second list can be done in N1' * N2' iterations,
  // where N1' is the reduced N1 and N2' is the reduced N2.
  // We have to break up the assigning of MotionItems into two separate
  // loops and be careful about the local copies.  These 4 loops are
  // significantly faster than the original for loop -- _especially_
  // for large problems.  (ndk)

#ifdef CH_MPI  // don't need to do this in serial
  // make a vector of boxes (or LayoutIndexes to boxes) from destination layout
  // that are known to reside on this processor.
  vector<DataIndex> vectorDestDI;
  vector<DataIndex> vectorDestOnProcDI;
  for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
  {
    vectorDestDI.push_back(DataIndex(to()));
    if (myprocID == dest.procID(to()))
    {
      vectorDestOnProcDI.push_back(DataIndex(to()));
    }
  }

  // make a vector of boxes (or LayoutIndexes to boxes) from "level"/src layout
  // that are known to reside on this processor.
  vector<DataIndex> vectorLevelDI;
  vector<DataIndex> vectorLevelOnProcDI;
  for (LayoutIterator from(a_level.layoutIterator()); from.ok(); ++from)
  {
    vectorLevelDI.push_back(DataIndex(from()));
    if (myprocID == level.procID(from()))
    {
      vectorLevelOnProcDI.push_back(DataIndex(from()));
    }
  }
#else
  // in serial, it's not very interesting as it's all of them.
  vector<DataIndex> vectorDestOnProcDI;
  vector<DataIndex> vectorLevelDI;
  for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
  {
    vectorDestOnProcDI.push_back(DataIndex(to()));
  }
  for (LayoutIterator from(a_level.layoutIterator()); from.ok(); ++from)
  {
    vectorLevelDI.push_back(DataIndex(from()));
  }
#endif

  bool isSorted = (a_level.isSorted() && a_dest.isSorted());

  // loop over all dest/to DI's on my processor
  for (vector<DataIndex>::iterator vdi = vectorDestOnProcDI.begin();
      vdi != vectorDestOnProcDI.end(); ++vdi)
  {
    // at this point, i know myprocID == toProcID
    const DataIndex todi(*vdi);

    Box ghost(dest[todi]);
    ghost -= a_shift;

    ghost.grow(a_ghost);

    //bool isSorted = (a_level.isSorted() && a_dest.isSorted());
    // then for each level/from DI, see if they intersect
    for (vector<DataIndex>::iterator vli = vectorLevelDI.begin();
        vli != vectorLevelDI.end(); ++vli)
    {
      const DataIndex fromdi(*vli);
      const unsigned int fromProcID = level.procID(fromdi);
      const Box& fromBox = level[fromdi];
      if ((fromBox.bigEnd(0) < ghost.smallEnd(0)) && isSorted )
      {
        //can skip rest cuz we haven't gotten to something interesting
        continue;
      }

      if (ghost.intersectsNotEmpty(fromBox))
      {
        Box srcBox(ghost); // ??
        srcBox &= fromBox; // ??

        Box destBox = srcBox + a_shift;

        MotionItem* item = new (s_motionItemPool.getPtr())
          MotionItem(fromdi, todi, srcBox, destBox);
        if (item == NULL)
        {
          MayDay::Error("Out of Memory in copier::define");
        }
        if (fromProcID == myprocID)
        { // local move
          if (a_exchange && fromdi == todi)
            s_motionItemPool.returnPtr(item);
          else
            m_localMotionPlan.push_back(item);
        }
        else
        {
          item->procID = fromProcID;
          m_toMotionPlan.push_back(item);
        }
      }
      if ((fromBox.smallEnd(0) > ghost.bigEnd(0)) && isSorted)
      {
        //can break out of loop, since we know that the smallEnd
        // of all the remaining boxes are lexigraphically beyond this ghosted box.
        break;
      }

    }
  }

  // Don't need to worry about this in serial as we already
  // took care of the local copy motion items just above.  skip this.
#ifdef CH_MPI
  // loop over all dest/to DI's
  for (vector<DataIndex>::iterator vdi = vectorDestDI.begin();
      vdi != vectorDestDI.end(); ++vdi)
  {

    const DataIndex todi(*vdi);

    Box ghost(dest[todi]);
    ghost -= a_shift;

    ghost.grow(a_ghost);

    const unsigned int toProcID = dest.procID(todi);

    // then for each level/from DI on this processor, see if they intersect
    for (vector<DataIndex>::iterator vli = vectorLevelOnProcDI.begin();
        vli != vectorLevelOnProcDI.end(); ++vli)
    {

      // at this point, i know myprocID == fromProcID

      const DataIndex fromdi(*vli);
      const Box& fromBox = level[fromdi];

      if ((fromBox.bigEnd(0) < ghost.smallEnd(0)) && isSorted)
      {
        //can skip rest cuz we haven't gotten to something interesting
        continue;
      }

      if (ghost.intersectsNotEmpty(fromBox))
      {
        Box srcBox(ghost); // ??
        srcBox &= fromBox; // ??

        Box destBox = srcBox + a_shift;

        if (toProcID == myprocID)
        { // local move
          // don't push back here!  or you will get two.
          //     we already did it above...
          //m_localMotionPlan.push_back(item);
        }
        else
        {
          MotionItem* item = new (s_motionItemPool.getPtr())
            MotionItem(fromdi, todi, srcBox, destBox);
          if (item == NULL)
          {
            MayDay::Error("Out of Memory in copier::define");
          }

          item->procID = toProcID;
          m_fromMotionPlan.push_back(item);
        }
      }
      if ((fromBox.smallEnd(0) > ghost.bigEnd(0)) && isSorted)
      {
        //can break out of loop, since we know that the smallEnd
        // of all the remaining boxes are lexigraphically beyond this ghosted box.
        break;
      }
    }
  }
#endif

  if (isPeriodic && a_shift != IntVect::Zero)
  {
    MayDay::Error("Copier::define - domain periodic and a non-zero shift in the copy is not implemented");
  }

  // put periodic intersection checking in here for "to" boxes
  if (isPeriodic)
  {
    for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
    {
      Box ghost(dest[to()]);
      ghost.grow(a_ghost);
      //unsigned int toProcID = dest.procID(to());  // unused variable

      // only do this if ghost box hangs over domain edge
      if (!domainBox.contains(ghost))
      {
        // add the dataIndex for this box to the list
        // of boxes which we need to come back to
        periodicallyFilledToVect.push_back(DataIndex(to()));
        // now check to see if we need to grow the
        // periodic check radius
        if (!grownDomainCheckBox.contains(ghost))
        {
          // grow the domainCheckBox until it contains ghost
          while (!grownDomainCheckBox.contains(ghost))
          {
            grownDomainCheckBox.grow(1);
            periodicCheckRadius++;
            for (int dir=0; dir<SpaceDim; dir++)
              {
                if (ghost.smallEnd()[dir] < 0) negativeStretch[dir] = true;
              }
          }
        } // end if we need to grow radius around domain

      } //end if ghost box is not contained in domain
    } // end if periodic
  }

  // Here ends the so-called N-squared optimizations.  the rest is unchanged. (ndk)

  // now do periodic checking, if necessary
  if (isPeriodic)
    {

      // the only "from" boxes we will need to check
      // will be those within periodicCheckRadius of the
      // domain boundary. so, create a box to screen out
      // those which we will need to check.
      Box shrunkDomainBox = a_domain.domainBox();
      shrunkDomainBox.grow(-periodicCheckRadius);

      IntVect numWraps = IntVect::Unit;
      IntVect shiftMult(domainBox.size());

      
      if (multipleWraps)
        {
          // make this be a function of direction in case 
          // our domains are very lopsided
          for (int dir=0; dir<SpaceDim; dir++)
            {
              numWraps[dir] = Max(numWraps[dir], periodicCheckRadius/shiftMult[dir]);
              // if we have boxes in the negative plane, need to add an extra wrap
              if (negativeStretch[dir]) numWraps[dir] += 1;
            }
        }

      
      ShiftIterator shiftIt = a_domain.shiftIterator();
      // if multipleWraps, make a custom ShiftIterator here
      if (multipleWraps)
        {
          const bool* isPeriodic = a_domain.isPeriodicVect();
          shiftIt.computeShifts(isPeriodic, numWraps);
        }

      // now loop over "from" boxes
      for (LayoutIterator from(a_level.layoutIterator()); from.ok(); ++from)
        {
          // first check to see whether we need to look at this box
          const Box& fromBox = level[from()];

          if (!shrunkDomainBox.contains(fromBox))
            {
              unsigned int fromProcID = level.procID(from());

              // check to see if fromBox is contained in domain,
              // if not, add it to the list of fromBoxes we need to
              // go back and check separately to see if it will
              // fill one of the "to" boxes
              if (!domainBox.contains(fromBox))
                {
                  periodicFromVect.push_back(DataIndex(from()));

                  if (!grownFromDomainCheckBox.contains(fromBox))
                    {
                      while (!grownFromDomainCheckBox.contains(fromBox))
                        {
                          grownFromDomainCheckBox.grow(1);
                          periodicFromCheckRadius++;
                        }
                    } // end if we need to grow domain check box
                } // end if fromBox is outside domain

              // now loop over those "to" boxes which were not contained
              // in the domain
              for (int toRef = 0; toRef < periodicallyFilledToVect.size(); toRef++)
                {
                  DataIndex toIndex = periodicallyFilledToVect[toRef];
                  unsigned int toProcID = dest.procID(toIndex);

                  // don't worry about anything that doesn't involve this proc
                  if (toProcID != myprocID && fromProcID != myprocID)
                    {
                      // do nothing
                    }
                  else
                    {
                      Box ghost(dest[toIndex]);
                      ghost.grow(a_ghost);
                      // now need to loop over shift vectors and look at images
                      for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                        {
                          IntVect shiftVect(shiftIt()*shiftMult);
                          ghost.shift(shiftVect);
                          if (ghost.intersectsNotEmpty(fromBox)) // rarely happens
                            {
                              Box intersectBox(ghost);
                              intersectBox &= fromBox;
                              Box toBox(intersectBox);
                              toBox.shift(-shiftVect);
                              MotionItem* item = new (s_motionItemPool.getPtr())
                                MotionItem(DataIndex(from()), DataIndex(toIndex),
                                           intersectBox, toBox);
                              if (item == NULL)
                                {
                                  MayDay::Error("Out of Memory in copier::define");
                                }
                              if (toProcID == fromProcID) // local move
                                m_localMotionPlan.push_back(item);
                              else if (fromProcID == myprocID)
                                {
                                  item->procID = toProcID;
                                  m_fromMotionPlan.push_back(item);
                                }
                              else
                                {
                                  item->procID = fromProcID;
                                  m_toMotionPlan.push_back(item);
                                }

                            } // end if shifted box intersects

                          ghost.shift(-shiftVect);
                        } // end loop over shift vectors
                    } // end if either from box or to box are on this proc
                } // end loop over destination boxes
            } // end if source box is close to domain boundary
        } // end loop over destination boxes

      // now go back through the "from" boxes which were outside
      // the domain and see if they intersect any toBoxes
      if (periodicFromVect.size() != 0)
        {
          // the only "to" boxes we will need to check
          // will be those within periodicCheckRadius of the
          // domain boundary. so, create a box to screen out
          // those which we will need to check.
          shrunkDomainBox = a_domain.domainBox();
          shrunkDomainBox.grow(-periodicFromCheckRadius);

          // now loop over the "to" boxes
          for (LayoutIterator to(a_dest.layoutIterator()); to.ok(); ++to)
            {
              // first check to see whether we need to look at this box
              Box ghost(dest[to()]);
              ghost.grow(a_ghost);

              if (!shrunkDomainBox.contains(ghost))
                {
                  unsigned int toProcID = a_dest.procID(to());

                  // now loop over those "from" boxes which are not
                  // contained by the domain
                  for (int fromRef = 0; fromRef < periodicFromVect.size(); fromRef++)
                    {
                      DataIndex fromIndex = periodicFromVect[fromRef];
                      const Box& fromBox = level[fromIndex];
                      unsigned int fromProcID = level.procID(fromIndex);

                      // don't worry about anything which doesn't involve
                      // this proc
                      if (toProcID != myprocID && fromProcID != myprocID)
                        {
                          // do nothing
                        }
                      else
                        {
                          // now need to loop over shift vectors and look at images
                          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                            {
                              IntVect shiftVect(shiftIt()*shiftMult);
                              ghost.shift(shiftVect);
                              if (ghost.intersectsNotEmpty(fromBox))
                                {
                                  Box intersectBox(ghost);
                                  intersectBox &= fromBox;
                                  Box toBox(intersectBox);
                                  toBox.shift(-shiftVect);
                                  MotionItem* item = new (s_motionItemPool.getPtr())
                                    MotionItem(DataIndex(fromIndex), DataIndex(to()),
                                               intersectBox, toBox);
                                  if (item == NULL)
                                    {
                                      MayDay::Error("Out of Memory in copier::define");
                                    }
                                  if (toProcID == fromProcID) // local move
                                    m_localMotionPlan.push_back(item);
                                  else if (fromProcID == myprocID)
                                    {
                                      item->procID = toProcID;
                                      m_fromMotionPlan.push_back(item);
                                    }
                                  else
                                    {
                                      item->procID = fromProcID;
                                      m_toMotionPlan.push_back(item);
                                    }

                                } // end if shifted box intersects

                              ghost.shift(-shiftVect);
                            } // end loop over shift vectors
                        } // end if either from box or to box are on this proc
                    } // end loop over "from" boxes
                } // end if destination box is close to domain boundary
            } // end loop over destination boxes
        } // end if any of the "From" boxes were outside the domain

    } // end if we need to do anything for periodicity
  sort();
}

void Copier::ghostDefine(const DisjointBoxLayout& a_src,
                         const DisjointBoxLayout& a_dest,
                         const ProblemDomain& a_domain,
                         const IntVect& a_srcGhost)
{
  //first, define a regular copier operation
  define(a_dest, a_src,  a_domain, a_srcGhost);

  //now, reverse the direction of the operation.
  reverse();

  sort();
}

void Copier::ghostDefine(const DisjointBoxLayout& a_src,
                         const DisjointBoxLayout& a_dest,
                         const ProblemDomain& a_domain,
                         const IntVect& a_srcGhost,
                         const IntVect& a_destGhost)
{
  CH_TIME("Copier::ghostDefine");
  clear();
  BoxLayout src ; src.deepCopy(a_src);
  src.grow(a_srcGhost);
  src.close();
  define(src, a_dest, a_domain, a_destGhost);
}

void Copier::exchangeDefine(const DisjointBoxLayout& a_grids,
                            const IntVect& a_ghost, bool a_includeSelf)
{
  CH_TIME("Copier::exchangeDefine");
  clear();
  m_isDefined=true;  
  DataIterator dit = a_grids.dataIterator();
  NeighborIterator nit(a_grids);
  int myprocID = procID();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& b = a_grids[dit];
      Box bghost(b);
      bghost.grow(a_ghost);
      if(a_includeSelf)
      {
        MotionItem* item = new (s_motionItemPool.getPtr())
                MotionItem(dit(), dit(), bghost);
        m_localMotionPlan.push_back(item);
      } 
      for (nit.begin(dit()); nit.ok(); ++nit)
        {
          Box neighbor = nit.box();
          int fromProcID = a_grids.procID(nit());
          if (neighbor.intersectsNotEmpty(bghost))
            {
              Box box(neighbor & bghost);

              MotionItem* item = new (s_motionItemPool.getPtr())
                MotionItem(DataIndex(nit()), dit(), nit.unshift(box), box);
              if (fromProcID == myprocID)
              { // local move
                m_localMotionPlan.push_back(item);
              }
              else
              {
                item->procID = fromProcID;
                m_toMotionPlan.push_back(item);
              }
            }
          neighbor.grow(a_ghost);
          if (neighbor.intersectsNotEmpty(b) && fromProcID != myprocID)
            {
              Box box(neighbor & b);
              MotionItem* item = new (s_motionItemPool.getPtr())
                MotionItem(dit(), DataIndex(nit()), box, nit.unshift(box) );
              item->procID = fromProcID;
              m_fromMotionPlan.push_back(item);
            }
        }

    }
  sort();
}

class MotionItemSorter
{
public:
  inline bool operator()(MotionItem* const & lhs, MotionItem* const& rhs) const;
};

inline bool  MotionItemSorter::operator()(MotionItem* const & lhs, MotionItem* const& rhs) const
{
  if (lhs->procID == rhs->procID)
    {
      const Box& left = lhs->toRegion;
      const Box& right = rhs->toRegion;
      if (left.smallEnd() == right.smallEnd())
        {
          return left.bigEnd().lexLT(right.bigEnd());
        }
      else
        {
          return left < right;
        }
    }
  return lhs->procID < rhs->procID;
}

void Copier::sort()
{
  //return;
  /*
  std::vector<MotionItem*>& vlocal = m_localMotionPlan.stdVector();
  std::sort(vlocal.begin(), vlocal.end(), [](MotionItem* const & lhs, MotionItem* const& rhs)->bool{ return lhs->toIndex.intCode()<rhs->toIndex.intCode();});
  m_range.resize(0);
  int items = vlocal.size();
  if(items > 0)
    {
      m_range.push_back(IndexTM<int,2>::Zero);
       for (int n=0; n<items-1; n++)
         {
           const MotionItem& item = *(vlocal[n]);
           //  pout()<<item.toIndex.intCode()<<" ";
           const MotionItem& itemNext = *(vlocal[n+1]);
           if(item.toIndex != itemNext.toIndex)
             {
               m_range.back()[1]=n+1;
               m_range.push_back(m_range.back());
               m_range.back()[0]=n+1;
             }
         }
       m_range.back()[1]=items;
    }
  */
  std::vector<MotionItem*>& vfrom  = m_fromMotionPlan.stdVector();
  std::sort(vfrom.begin(), vfrom.end(), MotionItemSorter());
  std::vector<MotionItem*>& vto = m_toMotionPlan.stdVector();
  std::sort(vto.begin(), vto.end(), MotionItemSorter());
}

int Copier::print() const
{
  pout()  << *this;
  return 0;
}

int Copier::numLocalCellsToCopy() const
{
  int sum = 0;
  for (unsigned int i = 0; i < m_localMotionPlan.size(); ++i)
    {
      sum += m_localMotionPlan[i]->fromRegion.numPts();
    }
  return sum;
}

int Copier::numFromCellsToCopy() const
{
  int sum = 0;
  for (unsigned int i = 0; i < m_fromMotionPlan.size(); ++i)
    {
      sum += m_fromMotionPlan[i]->fromRegion.numPts();
    }
  return sum;
}

int Copier::numToCellsToCopy() const
{
  int sum = 0;
  for (unsigned int i = 0; i < m_toMotionPlan.size(); ++i)
    {
      sum += m_toMotionPlan[i]->fromRegion.numPts();
    }
  return sum;
}

ostream& operator<< (ostream& os, const Copier& copier)
{
  os << "local(" << procID() << "): ";
  for (CopyIterator it(copier, CopyIterator::LOCAL); it.ok(); ++it)
    {
      os << " from " << it().fromIndex << it().fromRegion << " to " << it().toIndex << it().toRegion << '\n'
         << "          ";
    }
  os << "\nfrom(" << procID() << "): ";
  for (CopyIterator it(copier, CopyIterator::FROM); it.ok(); ++it)
    {
      os << " from " << it().fromIndex << it().fromRegion << " to " << it().toIndex << it().toRegion << "[" << it().procID << "]" << '\n'
         << "         ";
    }
  os << "\nto(" << procID() << "): ";
  for (CopyIterator it(copier, CopyIterator::TO); it.ok(); ++it)
    {
      os << " from " << it().fromIndex << it().fromRegion << " to " << it().toIndex << it().toRegion << "[" <<it().procID << "]" << '\n'
         << "       ";
    }
  os << "\n";
  return os;
}

Copier::~Copier()
{
  CH_TIME("~Copier");
  clear();
}

void Copier::setBufferAllocated(bool arg) const
{
  buffersAllocated  = arg;
}

bool Copier::bufferAllocated() const
{
  return buffersAllocated;
}

const ProblemDomain&
Copier::getPhysDomain(const DisjointBoxLayout& a_level) const
{
  return a_level.physDomain();
}

#include "NamespaceFooter.H"
