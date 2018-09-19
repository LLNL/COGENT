#ifndef _FACEBOX_H_
#define _FACEBOX_H_

#include "Box.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "NamespaceHeader.H"

/// Chombo's FluxBox but only along one dimension
/** 
 * This is a copy of Chombo's FluxBox, but defined along
 * only one dimension.
 */
class FaceBox
{

public:
  /// Default constructor
  FaceBox();

  /// Constructs FaceBox on cell-centered box with n components
  FaceBox(const Box& bx, int a_dir, int n);
  FaceBox(const Box& bx, int n=1);

  /// Constructs FaceBox aliasing a memory pointer for the FArrayBox
  FaceBox(const Box& a_bx, const int a_dir, const int a_nComp, Real *const a_alias);
  FaceBox(const Box& a_bx, const int a_nComp, Real *const a_alias);

  /// Destructor
  ~FaceBox();

  /// Resize FaceBox similar to BaseFab::resize()
  void resize(const Box& bx, int dir, int n);
  void resize(const Box& bx, int n=1);

  /// Define function
  void define(const Box& bx, int d, int n);
  void define(const Box& bx, int n=1);
  
  void define(const FaceBox& a_fb)
  {
    define(a_fb.box(), a_fb.dir(), a_fb.nComp());
  }

  /// Modify the dimension along which this FaceBox represents faces
  void reDir(int a_dir);

  /// Returns the FaceBox to the undefined state
  void clear();

  //  access functions

  /// Number of components
  int nComp() const;

  /// Dimension
  int dir() const;

  /// Returns cell-centered box which defines faceBox
  const Box& box() const;

  /// Returns face-centered flux in direction \em dir
  FArrayBox& data();

  /// Returns const reference to face-centered flux in direction \em dir
  const FArrayBox& data() const;

  /// Returns a pointer
  FArrayBox* dataPtr();

  //  data modification functions

  /// Set all fluxes to val
  void setVal(const Real val);

  /// More specific setVal
  void setVal(const Real val, const int startComp, const int nComp);

  /// Sets fluxes on faces surrounding cell-centered box bx
  void setVal(const Real val, const Box& bx);

  /// Most specific setVal
  /**
     Sets fluxes on faces surrounding cell-centered box bx
  */
  void setVal(const Real val, const Box& bx, const int startComp, const int nComp);

  /// Copy from src to this FaceBox -- sizes must be identical
  void copy(const FaceBox& src);

  /// Copy on overlap, for all directions
  void copy(const FaceBox& src, const int srcComp,
            const int destComp, const int numComp);

    /**
       Modifies this FaceBox by copying the contents of the argument src
       into it.  A copy within the intersecting region of the domains of the
       two FaceBoxes and the specified Box a_destbox is performed.  All
       components are copied.
    */
  void copy(const FaceBox& a_src,
            const Box&     a_destbox);

  /// Copies from a subsection of one box into another.
  /**
     Assumes the boxes are both in the same index space, and that
     box R is completely contained in both the src and destination
     boxes.  This function required by BoxLayoutData
  */
  void copy(const Box& R, const Interval& Cdest,  const FaceBox& src,
            const Interval& Csrc);

  /* Copy to a cell-centered FArrayBox */
  void copyTo(FArrayBox&) const;

  /// Modifies this FaceBox by copying the contents of \em src into it.
  /**
     This, the most general form
     of copy, specifies the contents of any sub-box srcbox
     in `FaceBox' src may be copied into a (possibly
     different) destbox in the destination `FaceBox'.  Note
     that although the srcbox and the destbox may be
     disjoint, they must be the same size and shape.  If the
     sizes differ, the copy is undefined and a runtime error
     results.  This copy function is the only one of the
     copy functions to allow a copy between differing
     boxes. The user also specifies how many components are
     copied, starting at component srccomp in src and stored
     starting at component destcomp.  The results are
     UNDEFINED if the src and dest FaceBoxes are the same and
     the srcbox and destbox overlap.
  */
  void copy (const Box&        srcbox,
             const Interval&   destcomps,
             const Box&        destbox,
             const FaceBox&    src,
             const Interval&   srccomps);

  /// Modifies this FaceBox to its additive inverse.
  /**
     Modifies this FaceBox by replacing each value with its additive
     inverse, for the given range of components and within the given subbox.
     Returns *this.
  */
  FaceBox& negate (const Box& subbox,
                   int        comp=0,
                   int        numcomp=1);

  /// Modifies this FaceBox to its additive inverse.
  /**
     Modifies this FaceBox by replacing each value with its additive
     inverse, for the given range of components over the whole domain of
     the FaceBox.  Returns *this.
  */
  FaceBox& negate (int comp,
                   int numcomp=1);

  /// Modifies this FaceBox to its additive inverse.
  /**
     Modifies this FaceBox by replacing each value with its additive
     inverse for all components over the whole domain of the FaceBox.
     Returns *this.
  */
  FaceBox& negate ();

  void plus(const FaceBox& a_src);

  /// Modifies this FaceBox by adding src in the CELL-CENTERED sub-box
  /**
     Modifies this FaceBox by pointwise addition of values in the argument
     FArrayBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
     to this FaceBox's components (a_destcomp : a_destcomp+numcomp-1)
     where the domain of this FaceBox intersects the a_subbox.  The actual
     directionally-dependent subbox is a_subbox.surroundingNodes(dir);
     NOTE:
     a_subbox must be contained in the cell-centered Box of this FaceBox.
     Returns *this
   */
  FaceBox& plus(const FaceBox& a_src,
                const Box&     a_subbox,
                int            a_srccomp,
                int            a_destcomp,
                int            a_numcomp = 1);

  /// Modifies this FaceBox by subtracting src in the CELL-CENTERED sub-box
  /**
     Modifies this FaceBox by pointwise addition of values in the argument
     FArrayBox.  Subtracts src's components (a_srccomp : a_srccomp+a_numcomp-1)
     from this FaceBox's components (a_destcomp : a_destcomp+numcomp-1)
     where the domain of this FaceBox intersects the a_subbox.  The actual
     directionally-dependent subbox is a_subbox.surroundingNodes(dir);
     NOTE:
     a_subbox must be contained in the cell-centered Box of this FaceBox.
     Returns *this
   */
  FaceBox& minus(const FaceBox& a_src,
                 const Box&     a_subbox,
                 int            a_srccomp,
                 int            a_destcomp,
                 int            a_numcomp = 1);

  /// Modifies this FaceBox by multiplying src in the CELL-CENTERED sub-box
  /**
     Modifies this FaceBox by pointwise multiplication of values in the
     argument FArrayBox.  Multiplies src's components
     (a_srccomp : a_srccomp+a_numcomp-1) with this FaceBox's components
     (a_destcomp : a_destcomp+numcomp-1) where the domain of this FaceBox
     intersects the a_subbox.  The actual directionally-dependent subbox
     is a_subbox.surroundingNodes(dir);
     NOTE:
     a_subbox must be contained in the cell-centered Box of this FaceBox.
     Returns *this
   */
  FaceBox& mult(const FaceBox& a_src,
                const Box&     a_subbox,
                int            a_srccomp,
                int            a_destcomp,
                int            a_numcomp = 1);

  /// Modifies this FaceBox by dividing src in the CELL-CENTERED sub-box
  /**
     Modifies this FaceBox by pointwise division of values in the argument
     FArrayBox.  Divides src's components (a_srccomp : a_srccomp+a_numcomp-1)
     into this FaceBox's components (a_destcomp : a_destcomp+numcomp-1)
     where the domain of this FaceBox intersects the a_subbox.  The actual
     directionally-dependent subbox is a_subbox.surroundingNodes(dir);
     NOTE:
     a_subbox must be contained in the cell-centered Box of this FaceBox.
     Returns *this
   */
  FaceBox& divide(const FaceBox& a_src,
                  const Box&     a_subbox,
                  int            a_srccomp,
                  int            a_destcomp,
                  int            a_numcomp = 1);

  /// Modifies this FaceBox by adding the scalar Real r to all values.
  FaceBox& operator+= (Real r);

  /// Modifies this FaceBox by incrementing with the argument FaceBox.
  /**
     Modifies this FaceBox by pointwise addition of the values of the
     argument FaceBox.  You might come to grief if the domains of the
     FArrayBoxes don't match, just as in FArrayBox::plus().
  */
  FaceBox& operator+= (const FaceBox& f);

  /// Modifies this FaceBox by subtracting the scalar Real r to all values.
  FaceBox& operator-= (Real r);

  /// Modifies this FaceBox by decrementing with the argument FaceBox.
  /**
     Modifies this FaceBox by pointwise subtraction of the values of the
     argument FaceBox.  You might come to grief if the domains of the
     FaceBoxes don't match, just as in FArrayBox::minus().
  */
  FaceBox& operator-= (const FaceBox& f);

  /// Modifies this FaceBox by multiplying all values by the scalar Real r.
  FaceBox& operator*= (Real r);

  /// Modifies this FaceBox by multiplying by the argument FaceBox.
  /**
     Modifies this FaceBox by pointwise multiplication of the values by the
     argument FaceBox.  You might come to grief if the domains of the
     FaceBoxes don't match, just as in FArrayBox::mult().
  */
  FaceBox& operator*= (const FaceBox& f);

  /// Modifies this FaceBox by shifting its domain box.
  FaceBox& shift (const IntVect& v);

  /// Returns size of linearized data over \em bx.
  /**
     Returns size, in number of bytes, of a flat linear
     representation of data in components comps in faces around
     cell-centered box R
  */
  int size(const Box& bx, const Interval& comps) const;

  /// Writes a linear representation of this FaceBox
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
     Cell-centered Box over which this FaceBox is defined
  */
  Box m_bx;

  /**
     Number of variables on each face
  */
  int m_nvar;

  /**
    Dimension of this facebox (i.e., along which
    dimension it represents faces).
  */
  int m_dir;

  /**
     CH_SPACEDIM FArrayBoxes which hold fluxes
  */
  FArrayBox* m_data;


private:
  // these are disallowed

  /// Disallowed
  FaceBox (const FaceBox&);
  /// Disallowed
  FaceBox& operator = (const FaceBox&);

};

#include "NamespaceFooter.H"
#endif
