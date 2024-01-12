#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERONEDSTENCIL_H_
#define _FOURTHORDERONEDSTENCIL_H_

#include "FArrayBox.H"
#include "IntVect.H"
#include "Vector.H"
#include "LAPACKMatrix.H"
#include "NamespaceHeader.H"

/// Class to manage coarse-to-fine spatial interpolation to all fine cells within a single particular coarse cell

/**
 */
class FourthOrderOneDStencil
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
  */
  FourthOrderOneDStencil();

    

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
  */
  ~FourthOrderOneDStencil()
  {
    if (m_weights != NULL)
      {
        for (int iWeight = 0; iWeight <=3; ++ iWeight)
          {
            delete m_weights[iWeight];
          }
        delete m_weights;
      }
  }

  /// Actual constructor.
  /**
     Set up object.
  */
  void define(const int& a_refRatio,
              const int& a_order   );

  void makeMomentMatrix(LAPACKMatrix & a_matrix,
                        const Real   & a_dx    ,
                        const Real   & a_vol   );
  
  void calculateStencilWeights(LAPACKMatrix       & a_weights  ,
                               const LAPACKMatrix & a_matrix   ,
                               const LAPACKMatrix & a_coarseData);

  void calculateMoments(Vector <Real> & a_oneDMomentMap,
                        const Real    & a_dx           ,
                        const Real    & a_vol          );

  void apply(FArrayBox      & a_fineFab         ,
             const FArrayBox& a_coarseFab       ,
             const IntVect  & a_coarseIv        ,
             const int      & a_interpDir       ,
             const IntVect  & a_coarseRefineVect);

protected:
 
public:
 int m_refRatio;
  int m_order  ;

  Real ** m_weights;
protected:

  /// whether define() has been called
  bool m_defined;
  
  
private:

  // Disallowed for all the usual reasons
  void operator=(const FourthOrderOneDStencil& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  FourthOrderOneDStencil(const FourthOrderOneDStencil& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
