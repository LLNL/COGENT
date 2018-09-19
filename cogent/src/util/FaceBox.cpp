#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FaceBox.H"
#include "NamespaceHeader.H"

// first do simple access functions
// ---------------------------------------------------------
int
FaceBox::nComp() const
{
  return m_nvar;
}

int
FaceBox::dir() const
{
  return m_dir;
}

// ---------------------------------------------------------
const Box&
FaceBox::box() const
{
  return m_bx;
}

// ---------------------------------------------------------
FArrayBox&
FaceBox::data()
{
  CH_assert(m_nvar >0);
  CH_assert(m_data != NULL);

  return *m_data;
}

// ---------------------------------------------------------
const FArrayBox&
FaceBox::data() const
{
  CH_assert(m_nvar >0);
  CH_assert(m_data != NULL);

  return *m_data;
}

// ---------------------------------------------------------
FArrayBox*
FaceBox::dataPtr()
{
  CH_assert(m_nvar >0);
  CH_assert(m_data != NULL);

  return m_data;
}

// ---------------------------------------------------------
// constructors and destructors
// ---------------------------------------------------------
FaceBox::FaceBox() : m_data(NULL)
{
  m_nvar = -1;
}

// ---------------------------------------------------------
FaceBox::FaceBox(const Box& a_bx, int a_dir, int a_nComp)  : m_data(NULL)
{
  define(a_bx, a_dir, a_nComp);
}

FaceBox::FaceBox(const Box& a_bx, int a_nComp)  : m_data(NULL)
{
  define(a_bx, a_nComp);
}

// ---------------------------------------------------------
FaceBox::FaceBox(const Box&   a_bx,
                 const int    a_dir,
                 const int    a_nComp,
                 Real *const  a_alias)
  :
  m_bx(a_bx),
  m_nvar(a_nComp),
  m_dir(a_dir),
  m_data(NULL)
{
  CH_assert(a_nComp > 0);
  CH_assert((a_dir >=0) && (a_dir < SpaceDim));
  CH_assert(a_alias != NULL);
  Box edgeBox(surroundingNodes(a_bx, a_dir));
  m_data = new FArrayBox(edgeBox, a_nComp, a_alias);
}

FaceBox::FaceBox(const Box&   a_bx,
                 const int    a_nComp,
                 Real *const  a_alias)
  :
  m_bx(a_bx),
  m_nvar(a_nComp),
  m_dir(0),
  m_data(NULL)
{
  CH_assert(a_nComp > 0);
  CH_assert(a_alias != NULL);
  Box edgeBox(surroundingNodes(a_bx, m_dir));
  m_data = new FArrayBox(edgeBox, a_nComp, a_alias);
}

// ---------------------------------------------------------
FaceBox::~FaceBox()
{
  clear();
}

// ---------------------------------------------------------
void
FaceBox::clear()
{
  // first delete storage
  if (m_data != NULL)
  {
    delete m_data;
    m_data = NULL;
  }

  // now reset all other variables
  m_nvar = -1;
  m_dir = -1;

  // set the box to the empty box...
  m_bx = Box();
}

// ---------------------------------------------------------
// define function
void
FaceBox::define(const Box& a_bx, int a_dir, int a_nComp)
{
  CH_assert(a_nComp > 0);
  CH_assert((a_dir >= 0) && (a_dir < SpaceDim));

  m_bx = a_bx;
  m_nvar = a_nComp;
  m_dir = a_dir;

  Box edgeBox(surroundingNodes(m_bx,a_dir));
  if (m_data != NULL) {
    delete m_data;
  }
  m_data = new FArrayBox(edgeBox,m_nvar);
}

void
FaceBox::define(const Box& a_bx, int a_nComp)
{
  CH_assert(a_nComp > 0);

  m_bx = a_bx;
  m_nvar = a_nComp;
  m_dir = 0;

  Box edgeBox(surroundingNodes(m_bx,m_dir));
  if (m_data != NULL) {
    delete m_data;
  }
  m_data = new FArrayBox(edgeBox,m_nvar);
}

// ---------------------------------------------------------
// should resize fluxes in space (could be faster than re-allocating
// storage)
void
FaceBox::resize(const Box& a_bx, int a_dir, int a_nComp)
{
  // if this object has not already been defined, call define fn.
  if (m_nvar < 0)
    {
      define(a_bx, a_dir, a_nComp);
    }
  else
    {
      CH_assert(a_nComp > 0);
      CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
      m_bx = a_bx;
      m_nvar = a_nComp;
      m_dir = a_dir;
      
      Box edgeBox(surroundingNodes(m_bx, a_dir));
      if (m_data != NULL)
      {
        m_data->resize(edgeBox,m_nvar);
      }
      else
      {
        FArrayBox* newFabPtr = new FArrayBox(edgeBox, m_nvar);
        m_data = newFabPtr;
      }
    }
}

void
FaceBox::resize(const Box& a_bx, int a_nComp)
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
      
      Box edgeBox(surroundingNodes(m_bx, m_dir));
      if (m_data != NULL)
      {
        m_data->resize(edgeBox,m_nvar);
      }
      else
      {
        FArrayBox* newFabPtr = new FArrayBox(edgeBox, m_nvar);
        m_data = newFabPtr;
      }
    }
}

void 
FaceBox::reDir(int a_dir)
{
  CH_assert(m_nvar != -1);
  m_dir = a_dir;

  if (m_data != NULL)
  {
    delete m_data;
    m_data = NULL;
  }

  CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
  Box edgeBox(surroundingNodes(m_bx, a_dir));
  m_data = new FArrayBox(edgeBox, m_nvar);
}

// ---------------------------------------------------------
void
FaceBox::setVal(const Real val)
{
  CH_assert(m_nvar > 0);
  CH_assert(m_data != NULL);
  m_data->setVal(val);
}

// ---------------------------------------------------------
void
FaceBox::setVal(const Real val, const int startComp, const int nComp)
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
FaceBox::setVal(const Real val, const Box& bx)
{
  CH_assert(m_bx.contains(bx));
  CH_assert(m_data != NULL);
  // move cell-centered box to appropriate edge
  Box edgeBox(surroundingNodes(bx, m_dir));
  m_data->setVal(val,edgeBox,0,m_nvar);

}

// ---------------------------------------------------------
void
FaceBox::setVal(const Real val, const Box& bx, const int startComp, const int nComp)
{
  CH_assert(m_bx.contains(bx));
  CH_assert(m_data != NULL);
  CH_assert(startComp > -1);
  CH_assert(startComp +nComp <=m_nvar);

  Box edgeBox(surroundingNodes(bx, m_dir));
  m_data->setVal(val, edgeBox, startComp, nComp);
}

// ---------------------------------------------------------
// copy on intersection
void
FaceBox::copy(const FaceBox& a_src)
{
  CH_assert(a_src.nComp() == m_nvar);
  CH_assert(a_src.dir() == m_dir);
  CH_assert(m_data != NULL);

  // copy on intersection of cell-centered boxes
  Box overlap(m_bx);
  overlap &= a_src.m_bx;
  
  Box faceBox(overlap);
  faceBox.surroundingNodes(m_dir);
  m_data->copy(a_src.data(), faceBox);
}

// ---------------------------------------------------------
void
FaceBox::copy(const FaceBox& a_src, const int srcComp,
              const int destComp, const int numComp)
{
  // to ensure that neither comp is negative
  CH_assert(srcComp*destComp > -1);
  CH_assert(a_src.dir() == m_dir);
  CH_assert(srcComp+numComp <= a_src.nComp());
  CH_assert(destComp+numComp <= m_nvar);
  CH_assert(m_data != 0);
  
  const FArrayBox& srcFab = a_src.data();
  m_data->copy(srcFab,srcComp, destComp, numComp);
}

// ---------------------------------------------------------
void
FaceBox::copy(const FaceBox& a_src,
              const Box&     a_destbox)
{
  CH_assert(a_src.nComp() == m_nvar);
  CH_assert(a_src.dir() == m_dir);
  CH_assert(m_bx.contains(a_destbox));
  CH_assert(a_src.box().contains(a_destbox));

  Interval comps(0, m_nvar-1);
  copy(a_destbox, comps, a_src, comps);
}

// ---------------------------------------------------------
void
FaceBox::copy(const Box& R, const Interval& Cdest, const FaceBox& a_src,
              const Interval& Csrc)
{
  CH_assert(m_data != NULL);
  CH_assert(a_src.dir() == m_dir);
  Box Redge(R);
  Redge.surroundingNodes(m_dir);
  const FArrayBox& srcFab = a_src.data();
  // all this intersecting is necessary due to the face-centered nature of things
  Redge &= srcFab.box();
  Redge &= m_data->box();
  if (!Redge.isEmpty())
  {
    //this is probably wrong in periodic---dtg
    m_data->copy(Redge,Cdest,Redge,srcFab,Csrc);
  }
}

// ---------------------------------------------------------
void
FaceBox::copy(const Box& srcbox,
              const Interval& destcomps,
              const Box& destbox,
              const FaceBox& a_src,
              const Interval& srccomps)
{
  // boxes do need to be the same size
  CH_assert(srcbox.sameSize(destbox));
  CH_assert(a_src.dir() == m_dir);
  CH_assert (m_data != NULL);
  const FArrayBox& srcFab = a_src.data();
  
  Box srcEdgeBox;
  Box destEdgeBox(destbox);
  destEdgeBox.surroundingNodes(m_dir);
  // this is somewhat different if a_src and dest boxes don't coincide...
  if (srcbox == destbox)
  {
    // safety check -- due to edge-centered nature of things,
    // destbox may not be contained in m_dataes[dir]
    // (DFM 11-4-09) however, both a_src and dest boxes need to
    // be the same size assume, however, that we only need to
    // intersect with the dest, and not the a_src
    destEdgeBox &= (m_data->box());
    srcEdgeBox = destEdgeBox;
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
    destEdgeBox &= (m_data->box());
    srcEdgeBox = destEdgeBox;
    srcEdgeBox.shift(shiftVect);
  }
  
  // don't even bother to do this for an empty box
  if (!destEdgeBox.isEmpty())
  {
    m_data->copy(srcEdgeBox, destcomps,destEdgeBox, srcFab, srccomps);
  }

}

void
FaceBox::copyTo(FArrayBox& a_fab) const
{
  const FArrayBox& my_fab = data();
  const Box& cc_grid = a_fab.box();

  BoxIterator bit(cc_grid);
  for (bit.begin(); bit.ok(); ++bit) {
    IntVect ic = bit();
    IntVect ip(ic); ip[m_dir]++;

    Real val1, val2;
    my_fab.getVal(&val1, ic);
    my_fab.getVal(&val2, ip);

    Real val = 0.5 * (val1 + val2);
    a_fab.set(ic, 0, val);
  }

}

// ---------------------------------------------------------
FaceBox&
FaceBox::negate()
{
  CH_assert(m_data != NULL);
  m_data->negate();
  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::negate(int        comp,
                int        numcomp)
{
  CH_assert(m_data != NULL);
  m_data->negate(comp,numcomp);
  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::negate(const Box& subbox,
                int        comp,
                int        numcomp)
{
  CH_assert(m_data != NULL);
  m_data->negate(subbox,comp,numcomp);
  return *this;
}

// ---------------------------------------------------------
void
FaceBox::plus(const FaceBox& a_src)
{
  CH_assert(a_src.dir() == m_dir);
  m_data->plus(a_src.data());
}

FaceBox&
FaceBox::plus(const FaceBox& a_src,
              const Box&     a_subbox,
              int            a_srccomp,
              int            a_destcomp,
              int            a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));
  CH_assert(a_src.dir() == m_dir);

  Box faceBox(surroundingNodes(a_subbox,m_dir));
  // now call corresponding FArrayBox function
  m_data->plus(a_src.data(), faceBox, a_srccomp, a_destcomp, a_numcomp);

  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::minus(const FaceBox& a_src,
               const Box&     a_subbox,
               int            a_srccomp,
               int            a_destcomp,
               int            a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));
  CH_assert(a_src.dir() == m_dir);
  
  Box faceBox(surroundingNodes(a_subbox,m_dir));
  // now call corresponding FArrayBox function
  m_data->minus(a_src.data(), faceBox, a_srccomp, a_destcomp, a_numcomp);

  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::mult(const FaceBox& a_src,
              const Box&     a_subbox,
              int            a_srccomp,
              int            a_destcomp,
              int            a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));
  CH_assert(a_src.dir() == m_dir);
  
  Box faceBox(surroundingNodes(a_subbox,m_dir));
  // now call corresponding FArrayBox function
  m_data->mult(a_src.data(), faceBox, a_srccomp, a_destcomp, a_numcomp);

  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::divide(const FaceBox& a_src,
                const Box&     a_subbox,
                int            a_srccomp,
                int            a_destcomp,
                int            a_numcomp)
{
  CH_assert(m_bx.contains(a_subbox));
  CH_assert(a_src.box().contains(a_subbox));
  CH_assert(a_src.dir() == m_dir);
  
  Box faceBox(surroundingNodes(a_subbox,m_dir));
  // now call corresponding FArrayBox function
  m_data->divide(a_src.data(), faceBox, a_srccomp, a_destcomp, a_numcomp);

  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::operator+= (Real r)
{
  CH_assert(m_data != NULL);
  m_data->plus(r);
  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::operator+= (const FaceBox& f)
{
  CH_assert(f.dir() == m_dir);
  CH_assert(m_data != NULL);
  m_data->plus(f.data());
  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::operator-= (Real r)
{
  CH_assert(m_data != NULL);
  m_data->plus(-r);
  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::operator-= (const FaceBox& f)
{
  CH_assert(f.dir() == m_dir);
  CH_assert(m_data != NULL);
  m_data->minus(f.data());
  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::operator*= (Real r)
{
  CH_assert(m_data != NULL);
  m_data->mult(r);
  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::operator*= (const FaceBox& f)
{
  CH_assert(f.dir() == m_dir);
  CH_assert(m_data != NULL);
  m_data->mult(f.data());
  return *this;
}

// ---------------------------------------------------------
FaceBox&
FaceBox::shift(const IntVect& iv)
{
  m_bx.shift(iv);
  m_data->shift(iv);

  return *this;
}

// ---------------------------------------------------------
int
FaceBox::size(const Box& bx, const Interval& comps) const
{
  FArrayBox tempFab;
  Box edgeBox(surroundingNodes(bx,m_dir));
  return tempFab.size(edgeBox,comps);
}

// ---------------------------------------------------------
void
FaceBox::linearOut(void* buf, const Box& R, const Interval& comps) const
{
  linearOut2(buf, R, comps);
}

void*
FaceBox::linearOut2(void* buf, const Box& R, const Interval& comps) const
{
  Real* buffer = (Real*) buf;
  
  CH_assert(m_data != NULL);
  Box dirBox(surroundingNodes(R,m_dir));
  //int dirSize = m_data->size(dirBox, comps);
  void* newBuf = m_data->linearOut2(buffer, dirBox, comps);
  buffer = (Real*) newBuf;
  //buffer += dirSize/sizeof(Real);

  return (void*) buffer;
}

// ---------------------------------------------------------
void
FaceBox::linearIn(void* buf, const Box& R, const Interval& comps)
{
  linearIn2(buf, R, comps);
}

void*
FaceBox::linearIn2(void* buf, const Box& R, const Interval& comps)
{
  Real* buffer = (Real*) buf;
  
  CH_assert(m_data != NULL);
  Box dirBox(surroundingNodes(R,m_dir));
  //int dirSize = m_data->size(dirBox, comps);
  void* newBuf = m_data->linearIn2(buffer,dirBox,comps);
  buffer = (Real*) newBuf;
  //buffer += dirSize/sizeof(Real);
  
  return (void*) buffer;
}

#include "NamespaceFooter.H"
