#include "CFGBlock.H"

#include "NamespaceHeader.H"

CFGBlock::CFGBlock(const IntVect& lo_mapped_index,
                   const IntVect& hi_mapped_index,
                   const IntVect& is_periodic,
                   const IntVect& decomp)

   : m_lo_mapped_index(lo_mapped_index),
     m_hi_mapped_index(hi_mapped_index),
     m_is_periodic(is_periodic),
     m_decomp(decomp)
{
}


ProblemDomain CFGBlock::getConfigurationDomain() const
{
   bool isPeriodic[SpaceDim];

   for (int dir=0; dir<SpaceDim; ++dir) {
      isPeriodic[dir] = m_is_periodic[dir];
   }

   return ProblemDomain(Box(m_lo_mapped_index,m_hi_mapped_index), isPeriodic);
}

#include "NamespaceFooter.H"
