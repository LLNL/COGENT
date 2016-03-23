#ifndef _Mapping_H_
#define _Mapping_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"

#include "NamespaceHeader.H"

class Mapping
{
  public:

    Mapping() { m_isDefined = false; }
    ~Mapping() { }

    void define(const KineticSpeciesPtrVect&);

    inline bool isDefined() { return m_isDefined; }
    inline const LevelData<FArrayBox>& getPointMapping() const 
      { 
        CH_assert(m_isDefined);
        return m_pointMapping; 
      }

  private:
    bool m_isDefined;
    LevelData<FArrayBox> m_pointMapping;
};

/* Fortran ordering */
#define _BoxIndexMapping_(p,i,N,ndims) \
  { \
    p = i[ndims-1]; \
    for (int c=ndims-2; c>=0; c--) p = (p*N[c]) + i[c];\
  }

#include "NamespaceFooter.H"

#endif
