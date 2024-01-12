#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERTENSORINTERP_H_
#define _FOURTHORDERTENSORINTERP_H_


#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
using std::string;
#include  <iostream>
#include "parstream.H"
#include "CONSTANTS.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "ParmParse.H"
#include "LoadBalance.H"
#include "FourthOrderOneDInterp.H"
#include "FArrayBox.H"
#include "IntVect.H"
#include "Vector.H"
#include "LAPACKMatrix.H"

#include "NamespaceHeader.H"

/// Class to manage coarse-to-fine spatial interpolation to all fine cells within a single particular coarse cell

/**
 */
class FourthOrderTensorInterp
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
  */
  FourthOrderTensorInterp();

    

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
  */
  ~FourthOrderTensorInterp()
  {
    
  }

  /// Actual constructor.
  /**
     Set up object.
  */
  void define(LevelData<FArrayBox>      & a_fineData    ,
              const int                 & a_refRatio    ,
              const Box                 & a_fineDomain  ,
              const LevelData<FArrayBox>& a_coarseData  ,
              const Box                 & a_coarseDomain);

  void define(LevelData<FArrayBox>      & a_fineData        ,
              const int                 & a_refRatio        ,
              const Box                 & a_fineDomain      ,
              const LevelData<FArrayBox>& a_coarseData      ,
              const Box                 & a_coarseDomain    ,
              const int                 & a_coarseGhostsFill,
              const Interval            & a_fixedDims       );
  

  void makeDBLVec(Vector < DisjointBoxLayout>& a_dblVec  ,
                  const DisjointBoxLayout    & a_fineDbl ,
                  const int                  & a_refRatio);
    
  void makeSemiCoarsenVec(Vector < FourthOrderOneDInterp*>& a_fineInterpVec  ,
                          const int                       & a_nComp          ,
                          const int                       & a_refRatio       ,
                          const IntVect                   & a_coarseGhostVect,
                          const Box                       & a_coarseDomain  );

  void defineLevelDataVec(Vector <LevelData <FArrayBox>* >& a_levelDataVec,
                          const Vector< DisjointBoxLayout>& a_dblVec      ,
                          const int                       & a_nComp       ,
                          const IntVect                   & a_ghostVect   );

  void semiCoarsenInterpLoop(Vector <LevelData <FArrayBox>* >      & a_levelDataVec ,
                             const Vector < FourthOrderOneDInterp*>& a_fineInterpVec);


protected:
 
public:
 
protected:

  /// whether define() has been called
  bool m_defined;
  
  
private:

  // Disallowed for all the usual reasons
  void operator=(const FourthOrderTensorInterp& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  FourthOrderTensorInterp(const FourthOrderTensorInterp& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
