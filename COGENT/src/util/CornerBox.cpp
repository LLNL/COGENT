#include "CornerBox.H"
#include "NamespaceHeader.H"

// first do simple access functions
// ---------------------------------------------------------
int
CornerBox::nComp() const
{
  return m_nvar;
}

// ---------------------------------------------------------
const Box&
CornerBox::box() const
{
  return m_bx;
}

// ---------------------------------------------------------
FArrayBox&
CornerBox::data()
{
  CH_assert(m_nvar >0);
  CH_assert(m_data != NULL);

  return *m_data;
}

// ---------------------------------------------------------
const FArrayBox&
CornerBox::data() const
{
  CH_assert(m_nvar >0);
  CH_assert(m_data != NULL);

  return *m_data;
}

// ---------------------------------------------------------
FArrayBox*
CornerBox::dataPtr()
{
  CH_assert(m_nvar > 0);
  CH_assert(m_data != NULL);

  return m_data;
}

// ---------------------------------------------------------
// constructors and destructors
// ---------------------------------------------------------
CornerBox::CornerBox() : m_data(NULL)
{
  m_nvar = -1;
}

// ---------------------------------------------------------
CornerBox::CornerBox(const Box& a_bx, int a_nComp)  : m_data(NULL)
{
  define(a_bx, a_nComp);
}

// ---------------------------------------------------------
CornerBox::CornerBox( const Box&   a_bx,
                      const int    a_nComp,
                      Real *const  a_alias)
  :
  m_bx(a_bx),
  m_nvar(a_nComp),
  m_data(NULL)
{
  CH_assert(a_nComp > 0);
  CH_assert(a_alias != NULL);

  Box cornerBox(surroundingNodes(a_bx));
  m_data = new FArrayBox(cornerBox, a_nComp, a_alias);
}

// ---------------------------------------------------------
CornerBox::~CornerBox()
{
  clear();
}

// ---------------------------------------------------------
void
CornerBox::clear()
{
  // first delete storage
  if (m_data != NULL)
  {
    delete m_data;
    m_data = NULL;
  }

  // now reset all other variables
  m_nvar = -1;

  // set the box to the empty box...
  m_bx = Box();
}

// ---------------------------------------------------------
// define function
void
CornerBox::define(const Box& a_bx, int a_nComp)
{
  CH_assert(a_nComp > 0);

  m_bx = a_bx;
  m_nvar = a_nComp;

  Box cornerBox(surroundingNodes(m_bx));
  if (m_data != NULL) {
    delete m_data;
  }
  m_data = new FArrayBox(cornerBox,m_nvar);
}

// ---------------------------------------------------------
// should resize fluxes in space (could be faster than re-allocating
// storage)
void
CornerBox::resize(const Box& a_bx, int a_nComp)
{
  // if this object has not already been defined, call define fn.
  if (m_nvar < 0)
    {
      define(a_bx, a_nComp);
    }
  else
    {
      CH_assert(a_nComp > 0);
      m_bx = a_bx;
      m_nvar = a_nComp;
      
      Box cornerBox(surroundingNodes(m_bx));
      if (m_data != NULL)
      {
        m_data->resize(cornerBox,m_nvar);
      }
      else
      {
        FArrayBox* newFabPtr = new FArrayBox(cornerBox, m_nvar);
        m_data = newFabPtr;
      }
    }
}

// ---------------------------------------------------------
void
CornerBox::setVal(const Real val)
{
  CH_assert(m_nvar > 0);
  CH_assert(m_data != NULL);
  m_data->setVal(val);
}

// ---------------------------------------------------------
void
CornerBox::setVal(const Real val, const int startComp, const int nComp)
{
  CH_assert(startComp >-1);
  CH_assert(startComp + nComp <= m_nvar);
  CH_assert(m_data != NULL);

  for (int comp=startComp; comp < startComp+nComp; comp++)
    {
      m_data->setVal(val,comp);
    }

}

// ---------------------------------------------------------
void
CornerBox::setVal(const Real val, const Box& bx)
{
  CH_assert(m_bx.contains(bx));
  CH_assert(m_data != NULL);
  // move cell-centered box to appropriate edge
  Box cornerBox(surroundingNodes(bx));
  m_data->setVal(val,cornerBox,0,m_nvar);
}

// ---------------------------------------------------------
void
CornerBox::setVal(const Real val, const Box& bx, const int startComp, const int nComp)
{
  CH_assert(m_bx.contains(bx));
  CH_assert(m_data != NULL);
  CH_assert(startComp > -1);
  CH_assert(startComp +nComp <=m_nvar);

  Box cornerBox(surroundingNodes(bx));
  m_data->setVal(val, cornerBox, startComp, nComp);
}

// ---------------------------------------------------------
// copy on intersection
void
CornerBox::copy(const CornerBox& a_src)
{
  CH_assert(a_src.nComp() == m_nvar);
  CH_assert(m_data != NULL);

  // copy on intersection of cell-centered boxes
  Box overlap(m_bx);
  overlap &= a_src.m_bx;
  
  Box cornerBox(overlap);
  cornerBox.surroundingNodes();
  m_data->copy(a_src.data(), cornerBox);
}

// ---------------------------------------------------------
void
CornerBox::copy(const CornerBox& a_src, const int srcComp,
                const int destComp, const int numComp)
{
  // to ensure that neither comp is negative
  CH_assert(srcComp*destComp > -1);
  CH_assert(srcComp+numComp <= a_src.nComp());
  CH_assert(destComp+numComp <= m_nvar);
  CH_assert(m_data != 0);
  
  const FArrayBox& srcFab = a_src.data();
  m_data->copy(srcFab,srcComp, destComp, numComp);
}

// ---------------------------------------------------------
void
CornerBox::copy(const CornerBox& a_src,
                const Box&     a_destbox)
{
  CH_assert(a_src.nComp() == m_nvar);
  CH_assert(m_bx.contains(a_destbox));
  CH_assert(a_src.box().contains(a_destbox));

  Interval comps(0, m_nvar-1);
  copy(a_destbox, comps, a_src, comps);
}

// ---------------------------------------------------------
void
CornerBox::copy(const Box& R, const Interval& Cdest, const CornerBox& a_src,
                const Interval& Csrc)
{
  CH_assert(m_data != NULL);
  Box Rcorner(R);
  Rcorner.surroundingNodes();
  const FArrayBox& srcFab = a_src.data();
  // all this intersecting is necessary due to the node-centered nature of things
  Rcorner &= srcFab.box();
  Rcorner &= m_data->box();
  if (!Rcorner.isEmpty())
  {
    //this is probably wrong in periodic---dtg
    m_data->copy(Rcorner,Cdest,Rcorner,srcFab,Csrc);
  }
}

// ---------------------------------------------------------
void
CornerBox::copy(const Box& srcbox,
                const Interval& destcomps,
                const Box& destbox,
                const CornerBox& a_src,
                const Interval& srccomps)
{
  // boxes do need to be the same size
  CH_assert(srcbox.sameSize(destbox));
  CH_assert (m_data != NULL);
  const FArrayBox& srcFab = a_src.data();
  
  Box srcCornerBox;
  Box destCornerBox(destbox);
  destCornerBox.surroundingNodes();
  // this is somewhat different if a_src and dest boxes don't coincide...
  if (srcbox == destbox)
  {
    // safety check -- due to edge-centered nature of things,
    // destbox may not be contained in m_dataes[dir]
    // (DFM 11-4-09) however, both a_src and dest boxes need to
    // be the same size assume, however, that we only need to
    // intersect with the dest, and not the a_src
    destCornerBox &= (m_data->box());
    srcCornerBox = destCornerBox;
  }
  else
  {
    // (DFM 11-4-09) if a_src and dest boxes are not the same
    // (probably due to a periodic wrapping), then this is a bit
    // more complicated.
    
    // first compute the shift
    IntVect shiftVect(srcbox.smallEnd());
    shiftVect -= destbox.smallEnd();
    // intersect destBox with dest
    destCornerBox &= (m_data->box());
    srcCornerBox = destCornerBox;
    srcCornerBox.shift(shiftVect);
  }
  
  // don't even bother to do this for an empty box
  if (!destCornerBox.isEmpty())
  {
    m_data->copy(srcCornerBox, destcomps,destCornerBox, srcFab, srccomps);
  }

}

void
CornerBox::copyTo(FArrayBox& a_fab) const
{
  const FArrayBox& my_fab = data();
  const Box& cc_grid = a_fab.box();

  BoxIterator bit(cc_grid);
  for (bit.begin(); bit.ok(); ++bit) {
    IntVect ic = bit();

    IntVect i1(ic), i2(ic), i3(ic), i4(ic);
    i2[2]++;
    i3[3]++;
    i4[2]++; i4[3]++;

    Real val1, val2, val3, val4;
    my_fab.getVal(&val1, i1);
    my_fab.getVal(&val2, i2);
    my_fab.getVal(&val3, i3);
    my_fab.getVal(&val4, i4);

    Real val = 0.25 * (val1 + val2 + val3 + val4);
    a_fab.set(ic, 0, val);
  }

}

// ---------------------------------------------------------
CornerBox&
CornerBox::negate()
{
  CH_assert(m_data != NULL);
  m_data->negate();
  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::negate(int        comp,
                  int        numcomp)
{
  CH_assert(m_data != NULL);
  m_data->negate(comp,numcomp);
  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::negate(const Box& subbox,
                  int        comp,
                  int        numcomp)
{
  CH_assert(m_data != NULL);
  m_data->negate(subbox,comp,numcomp);
  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::plus(const CornerBox&  a_src,
                const Box&        a_subbox,
                int               a_srccomp,
                int               a_destcomp,
                int               a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));

  Box cornerBox(surroundingNodes(a_subbox));
  // now call corresponding FArrayBox function
  m_data->plus(a_src.data(), cornerBox, a_srccomp, a_destcomp, a_numcomp);

  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::minus( const CornerBox&  a_src,
                  const Box&        a_subbox,
                  int               a_srccomp,
                  int               a_destcomp,
                  int               a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));
  
  Box cornerBox(surroundingNodes(a_subbox));
  // now call corresponding FArrayBox function
  m_data->minus(a_src.data(), cornerBox, a_srccomp, a_destcomp, a_numcomp);

  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::mult(const CornerBox&  a_src,
                const Box&        a_subbox,
                int               a_srccomp,
                int               a_destcomp,
                int               a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));
  
  Box cornerBox(surroundingNodes(a_subbox));
  // now call corresponding FArrayBox function
  m_data->mult(a_src.data(), cornerBox, a_srccomp, a_destcomp, a_numcomp);

  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::divide(const CornerBox&  a_src,
                  const Box&        a_subbox,
                  int               a_srccomp,
                  int               a_destcomp,
                  int               a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));
  
  Box cornerBox(surroundingNodes(a_subbox));
  // now call corresponding FArrayBox function
  m_data->divide(a_src.data(), cornerBox, a_srccomp, a_destcomp, a_numcomp);

  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::operator+= (Real r)
{
  CH_assert(m_data != NULL);
  m_data->plus(r);
  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::operator+= (const CornerBox& f)
{
  CH_assert(m_data != NULL);
  m_data->plus(f.data());
  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::operator-= (Real r)
{
  CH_assert(m_data != NULL);
  m_data->plus(-r);
  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::operator-= (const CornerBox& f)
{
  CH_assert(m_data != NULL);
  m_data->minus(f.data());
  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::operator*= (Real r)
{
  CH_assert(m_data != NULL);
  m_data->mult(r);
  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::operator*= (const CornerBox& f)
{
  CH_assert(m_data != NULL);
  m_data->mult(f.data());
  return *this;
}

// ---------------------------------------------------------
CornerBox&
CornerBox::shift(const IntVect& iv)
{
  m_bx.shift(iv);
  m_data->shift(iv);

  return *this;
}

// ---------------------------------------------------------
int
CornerBox::size(const Box& bx, const Interval& comps) const
{
  FArrayBox tempFab;
  Box cornerBox(surroundingNodes(bx));
  return tempFab.size(cornerBox,comps);
}

// ---------------------------------------------------------
void
CornerBox::linearOut(void* buf, const Box& R, const Interval& comps) const
{
  linearOut2(buf, R, comps);
}

void*
CornerBox::linearOut2(void* buf, const Box& R, const Interval& comps) const
{
  Real* buffer = (Real*) buf;
  
  CH_assert(m_data != NULL);
  Box cornerBox(surroundingNodes(R));
  //int dirSize = m_data->size(cornerBox, comps);
  void* newBuf = m_data->linearOut2(buffer, cornerBox, comps);
  buffer = (Real*) newBuf;
  //buffer += dirSize/sizeof(Real);

  return (void*) buffer;
}

// ---------------------------------------------------------
void
CornerBox::linearIn(void* buf, const Box& R, const Interval& comps)
{
  linearIn2(buf, R, comps);
}

void*
CornerBox::linearIn2(void* buf, const Box& R, const Interval& comps)
{
  Real* buffer = (Real*) buf;
  
  CH_assert(m_data != NULL);
  Box cornerBox(surroundingNodes(R));
  //int dirSize = m_data->size(cornerBox, comps);
  void* newBuf = m_data->linearIn2(buffer,cornerBox,comps);
  buffer = (Real*) newBuf;
  //buffer += dirSize/sizeof(Real);
  
  return (void*) buffer;
}

#include "NamespaceFooter.H"
