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

#include "MagGeom.H"

#include "NamespaceHeader.H"

void
bwenoFaceAverages( LevelData<FluxBox>& facePhi,
                   const LevelData<FArrayBox>& cellPhi,
                   const LevelData<FluxBox>& faceVel,
                   const MagGeom& coord );

void
uw5FaceAverages( LevelData<FluxBox>&           a_face_phi,
                 const LevelData<FArrayBox>&   a_cell_phi,
                 const LevelData<FluxBox>&     a_face_vel,
                 const MagGeom&                a_geom );

void
weno5FaceAverages( LevelData<FluxBox>&           a_face_phi,
                   const LevelData<FArrayBox>&   a_cell_phi,
                   const LevelData<FluxBox>&     a_face_vel,
                   const MagGeom&                a_geom );

#include "NamespaceFooter.H"

#endif
