#ifndef _CFGBLOCK_H_
#define _CFGBLOCK_H_

#include "IntVect.H"
#include "ProblemDomain.H"

#include "Box.H"
#include "ProblemDomain.H"

#include "NamespaceHeader.H"

class CFGBlock
{
public:

   CFGBlock(const IntVect& lo_mapped_index,
            const IntVect& hi_mapped_index,
            const IntVect& is_periodic,
            const IntVect& decomp);

   ~CFGBlock() {;}

   ProblemDomain getConfigurationDomain() const;

   IntVect cfgDecomp() const {return m_decomp;}

   IntVect loMappedIndex() const {return m_lo_mapped_index;}

   IntVect hiMappedIndex() const {return m_hi_mapped_index;}

   IntVect isPeriodic() const {return m_is_periodic;}

private:

   IntVect m_lo_mapped_index;
   IntVect m_hi_mapped_index;
   IntVect m_is_periodic;
   IntVect m_decomp;
};

#include "NamespaceFooter.H"

#endif
