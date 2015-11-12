#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SPARSE_COUPLING_H_
#define _SPARSE_COUPLING_H_

#include <iostream>

#include "Vector.H"
#include "IntVect.H"
#include "REAL.H"
#include "SPMD.H"
#include "BaseFab.H"
#include "NamespaceHeader.H"

class SparseCoupling;

class SparseCouplingIterator
{
public:
  ///
  /**
     Default constructor.  This constructs an invalid iterator.
     The user must call define before using.
  */
  SparseCouplingIterator();

  ///
  /**
     Constructs a BoxIterator and associates it with a Box.
     Arguments:
     a_bx (not modified) the Box to iterate over.
  */
  SparseCouplingIterator( const SparseCoupling& stencil );

  ///
  /**
     Associates a Box with this BoxIterator.
     Arguments:
     a_bx (not modified) the Box to iterate over.
  */
  void define( const SparseCoupling& stencil );

  ///
  /**
     Copy constructor.
     Arguments:
     a_iterIn (not modified) the BoxIterator to copy.
  */
  SparseCouplingIterator( const SparseCouplingIterator& iterator );

  ///
  ~SparseCouplingIterator ()
  {
  }

  ///
  /**
     Sets this BoxIterator to the first IntVect in its Box.  The
     definition of the "first" IntVect is
     implementation-dependent.
  */
  void begin();

  ///
  /**
     Sets this BoxIterator to the first IntVect in its Box.  The
     definition of the "first" IntVect is
     implementation-dependent.
  */
  void reset();

  ///
  /**
     Modifies this BoxIterator to set it to the next location in its
     Box.  The definition of the "next location" of a Box is
     implementation-dependent.
  */
  void operator++();

  void next();

  ///
  /**
     Returns the value of the InVect for the current location of this
     BoxIterator.
  */
  const int& operator()() const;

  ///
  /**
     Returns true if this BoxIterator's location is within its Box.
  */
  bool ok();

protected:
  int m_current;
  int m_lo;
  int m_hi;
};


class SparseCoupling
{
public:
   ///
   SparseCoupling();
   
   ///
   virtual ~SparseCoupling();
   
   ///
   /**
      add a cell index dependency and its weight
      If the cell index is already in the stencil, accumulate the weights.
   */
   void add( const IntVect& index, const Real& weight );
   
   ///
   /**
      number of dependencies in the stencil
   */
   int size() const;
   
   ///
   /** access an index
    */
   const IntVect& operator[]( const int i ) const;
   
   ///
   /**
      access a weight
   */
   /**
    */
   const Real& weight( const IntVect& index ) const;

   ///
   /**
      access a weight
   */
   /**
    */
   const Real& weight( const int& index ) const;

   ///
   /**
    */
   SparseCouplingIterator iterator() const;
   
   ///
   /**
    */
   SparseCoupling( const SparseCoupling& stencil );
   
   ///
   /**
    */
   const SparseCoupling& operator=( const SparseCoupling& stencil );

   
   ///
   /**
    */
   bool operator!=( const SparseCoupling& stencil ) const;
   
   ///
   /**
    */
   bool operator==( const SparseCoupling& stencil ) const;

   ///
   /**
    */
   void print( std::ostream& out, const std::string& prefix = "" ) const;

   ///
   /**
    */
   int linearSize() const;

   ///
   /**
    */
   void linearIn( const void* const inBuf );

   ///
   /**
    */
   void linearOut( void* const outBuf ) const;

   ///
   /**
    */
   friend std::ostream& operator<<( std::ostream& out, const SparseCoupling& spc );

private:
   
   ///
   /**
    */
   int find( const IntVect& a_index ) const;
   
   Vector<IntVect> m_index;
   Vector<Real> m_weight;
};


template < >
int BaseFab<SparseCoupling>::preAllocatable();

template < >
int BaseFab<SparseCoupling>::size(const Box&      a_box,
                                  const Interval& a_comps) const;

template < >
void BaseFab<SparseCoupling>::linearOut(void*           a_buf,
                                        const Box&      a_R,
                                        const Interval& a_comps) const;
template < >
void BaseFab<SparseCoupling>::linearIn(void*           a_buf,
                                       const Box&      a_R,
                                       const Interval& a_comps);

#include "SparseCouplingI.H"

#include "NamespaceFooter.H"

#include "BaseNamespaceHeader.H"

#include "NamespaceVar.H"

template < > int linearSize( const CH_XDIR::SparseCoupling& a_input );

template < > void linearIn( CH_XDIR::SparseCoupling& a_outputT, const void* const a_inBuf );
template < > void linearOut( void* const a_outBuf, const CH_XDIR::SparseCoupling& a_inputT );

#include "BaseNamespaceFooter.H"


#endif


