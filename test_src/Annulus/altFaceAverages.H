#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _ALTFACEAVERAGES_H_
#define _ALTFACEAVERAGES_H_

#include "MultiBlockCoordSys.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

void
altFaceAverages( LevelData<FluxBox>& facePhi,
                 const LevelData<FArrayBox>& cellPhi,
                 const LevelData<FluxBox>& faceVel,
                 const MultiBlockCoordSys& coords);

#include "NamespaceFooter.H"

#endif
