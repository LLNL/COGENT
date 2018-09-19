#ifndef _CORNERBOX_H_
#define _CORNERBOX_H_

#include "Box.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "NamespaceHeader.H"

class CornerBox
{

public:
  /// Default constructor
  CornerBox();

  /// Constructs CornerBox on cell-centered box with n components
  CornerBox(const Box& bx, int n=1);

  /// Constructs CornerBox aliasing a memory pointer for the FArrayBox
  CornerBox(const Box& a_bx, const int a_nComp, Real *const a_alias);

  /// Destructor
  ~CornerBox();

  /// Resize CornerBox similar to BaseFab::resize()
  void resize(const Box& bx, int n=1);

  /// Define function
  void define(const Box& bx, int n=1);
  
  void define(const CornerBox& a_fb)
  {
    define(a_fb.box(), a_fb.nComp());
  }

  /// Returns the CornerBox to the undefined state
  void clear();

  //  access functions

  /// Number of components
  int nComp() const;

  /// Returns cell-centered box which defines cornerBox
  const Box& box() const;

  /// Returns node-centered data
  FArrayBox& data();

  /// Returns const reference to node-centered data
  const FArrayBox& data() const;

  /// Returns a pointer
  FArrayBox* dataPtr();

  //  data modification functions

  /// Set all data to val
  void setVal(const Real val);

  /// More specific setVal
  void setVal(const Real val, const int startComp, const int nComp);

  /// Sets data on faces surrounding cell-centered box bx
  void setVal(const Real val, const Box& bx);

  /// Most specific setVal
  /**
     Sets data on faces surrounding cell-centered box bx
  */
  void setVal(const Real val, const Box& bx, const int startComp, const int nComp);

  /// Copy from src to this CornerBox -- sizes must be identical
  void copy(const CornerBox& src);

  /// Copy on overlap, for all directions
  void copy(const CornerBox& src, const int srcComp,
            const int destComp, const int numComp);

    /**
       Modifies this CornerBox by copying the contents of the argument src
       into it.  A copy within the intersecting region of the domains of the
       two CornerBoxes and the specified Box a_destbox is performed.  All
       components are copied.
    */
  void copy(const CornerBox& a_src,
            const Box&     a_destbox);

  /// Copies from a subsection of one box into another.
  /**
     Assumes the boxes are both in the same index space, and that
     box R is completely contained in both the src and destination
     boxes.  This function required by BoxLayoutData
  */
  void copy(const Box& R, const Interval& Cdest,  const CornerBox& src,
            const Interval& Csrc);

  /// Modifies this CornerBox by copying the contents of \em src into it.
  /**
     This, the most general form
     of copy, specifies the contents of any sub-box srcbox
     in `CornerBox' src may be copied into a (possibly
     different) destbox in the destination `CornerBox'.  Note
     that although the srcbox and the destbox may be
     disjoint, they must be the same size and shape.  If the
     sizes differ, the copy is undefined and a runtime error
     results.  This copy function is the only one of the
     copy functions to allow a copy between differing
     boxes. The user also specifies how many components are
     copied, starting at component srccomp in src and stored
     starting at component destcomp.  The results are
     UNDEFINED if the src and dest CornerBoxes are the same and
     the srcbox and destbox overlap.
  */
  void copy (const Box&        srcbox,
             const Interval&   destcomps,
             const Box&        destbox,
             const CornerBox&    src,
             const Interval&   srccomps);

  /* Copy to a cell-centered FArrayBox */
  void copyTo(FArrayBox&) const;

  /// Modifies this CornerBox to its additive inverse.
  /**
     Modifies this CornerBox by replacing each value with its additive
     inverse, for the given range of components and within the given subbox.
     Returns *this.
  */
  CornerBox& negate ( const Box& subbox,
                      int   comp=0,
                      int   numcomp=1);

  /// Modifies this CornerBox to its additive inverse.
  /**
     Modifies this CornerBox by replacing each value with its additive
     inverse, for the given range of components over the whole domain of
     the CornerBox.  Returns *this.
  */
  CornerBox& negate ( int comp,
                      int numcomp=1);

  /// Modifies this CornerBox to its additive inverse.
  /**
     Modifies this CornerBox by replacing each value with its additive
     inverse for all components over the whole domain of the CornerBox.
     Returns *this.
  */
  CornerBox& negate ();

  /// Modifies this CornerBox by adding src in the CELL-CENTERED sub-box
  /**
     Modifies this CornerBox by pointwise addition of values in the argument
     FArrayBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
     to this CornerBox's components (a_destcomp : a_destcomp+numcomp-1)
     where the domain of this CornerBox intersects the a_subbox.  The actual
     directionally-dependent subbox is a_subbox.surroundingNodes(dir);
     NOTE:
     a_subbox must be contained in the cell-centered Box of this CornerBox.
     Returns *this
   */
  CornerBox& plus(const CornerBox&  a_src,
                  const Box&        a_subbox,
                  int               a_srccomp,
                  int               a_destcomp,
                  int               a_numcomp = 1);

  /// Modifies this CornerBox by subtracting src in the CELL-CENTERED sub-box
  /**
     Modifies this CornerBox by pointwise addition of values in the argument
     FArrayBox.  Subtracts src's components (a_srccomp : a_srccomp+a_numcomp-1)
     from this CornerBox's components (a_destcomp : a_destcomp+numcomp-1)
     where the domain of this CornerBox intersects the a_subbox.  The actual
     directionally-dependent subbox is a_subbox.surroundingNodes(dir);
     NOTE:
     a_subbox must be contained in the cell-centered Box of this CornerBox.
     Returns *this
   */
  CornerBox& minus( const CornerBox&  a_src,
                    const Box&        a_subbox,
                    int               a_srccomp,
                    int               a_destcomp,
                    int               a_numcomp = 1);

  /// Modifies this CornerBox by multiplying src in the CELL-CENTERED sub-box
  /**
     Modifies this CornerBox by pointwise multiplication of values in the
     argument FArrayBox.  Multiplies src's components
     (a_srccomp : a_srccomp+a_numcomp-1) with this CornerBox's components
     (a_destcomp : a_destcomp+numcomp-1) where the domain of this CornerBox
     intersects the a_subbox.  The actual directionally-dependent subbox
     is a_subbox.surroundingNodes(dir);
     NOTE:
     a_subbox must be contained in the cell-centered Box of this CornerBox.
     Returns *this
   */
  CornerBox& mult(const CornerBox&  a_src,
                  const Box&        a_subbox,
                  int               a_srccomp,
                  int               a_destcomp,
                  int               a_numcomp = 1);

  /// Modifies this CornerBox by dividing src in the CELL-CENTERED sub-box
  /**
     Modifies this CornerBox by pointwise division of values in the argument
     FArrayBox.  Divides src's components (a_srccomp : a_srccomp+a_numcomp-1)
     into this CornerBox's components (a_destcomp : a_destcomp+numcomp-1)
     where the domain of this CornerBox intersects the a_subbox.  The actual
     directionally-dependent subbox is a_subbox.surroundingNodes(dir);
     NOTE:
     a_subbox must be contained in the cell-centered Box of this CornerBox.
     Returns *this
   */
  CornerBox& divide(const CornerBox&  a_src,
                    const Box&        a_subbox,
                    int               a_srccomp,
                    int               a_destcomp,
                    int               a_numcomp = 1);

  /// Modifies this CornerBox by adding the scalar Real r to all values.
  CornerBox& operator+= (Real r);

  /// Modifies this CornerBox by incrementing with the argument CornerBox.
  /**
     Modifies this CornerBox by pointwise addition of the values of the
     argument CornerBox.  You might come to grief if the domains of the
     FArrayBoxes don't match, just as in FArrayBox::plus().
  */
  CornerBox& operator+= (const CornerBox& f);

  /// Modifies this CornerBox by subtracting the scalar Real r to all values.
  CornerBox& operator-= (Real r);

  /// Modifies this CornerBox by decrementing with the argument CornerBox.
  /**
     Modifies this CornerBox by pointwise subtraction of the values of the
     argument CornerBox.  You might come to grief if the domains of the
     CornerBoxes don't match, just as in FArrayBox::minus().
  */
  CornerBox& operator-= (const CornerBox& f);

  /// Modifies this CornerBox by multiplying all values by the scalar Real r.
  CornerBox& operator*= (Real r);

  /// Modifies this CornerBox by multiplying by the argument CornerBox.
  /**
     Modifies this CornerBox by pointwise multiplication of the values by the
     argument CornerBox.  You might come to grief if the domains of the
     CornerBoxes don't match, just as in FArrayBox::mult().
  */
  CornerBox& operator*= (const CornerBox& f);

  /// Modifies this CornerBox by shifting its domain box.
  CornerBox& shift (const IntVect& v);

  /// Returns size of linearized data over \em bx.
  /**
     Returns size, in number of bytes, of a flat linear
     representation of data in components comps in faces around
     cell-centered box R
  */
  int size(const Box& bx, const Interval& comps) const;

  /// Writes a linear representation of this CornerBox
  /**
     Write a linear representaion of the internal data for
     the faces surrounding cell-centered box R.  Assumes that sufficient
     memory for the buffer has already been allocated by the caller
  */
  void linearOut(void*buf, const Box& R, const Interval& comps) const;

  /// like linearOut, but returns current position in the buffer
  void* linearOut2(void*buf, const Box& R, const Interval& comps) const;

  /// Read a linear representation of the data over the Box \rm R
  /**
      Reads in the output of linearOut
  */
  void linearIn(void*buf, const Box& R, const Interval& comps);

  /// like linearInt, but returns current position in the buffer
  void* linearIn2(void*buf, const Box& R, const Interval& comps);

  /// Helper function for linearization
  static int preAllocatable()
  {
    return 0;
  }

protected:
  /**
     Box over which this CornerBox is defined
  */
  Box m_bx;

  /**
     Number of variables on each face
  */
  int m_nvar;

  /**
     CH_SPACEDIM FArrayBox which holds data
  */
  FArrayBox* m_data;


private:
  // these are disallowed

  /// Disallowed
  CornerBox (const CornerBox&);
  /// Disallowed
  CornerBox& operator = (const CornerBox&);

};

#include "NamespaceFooter.H"
#endif
