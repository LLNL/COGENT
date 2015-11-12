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

#include "PhaseGeom.H"

#include "NamespaceHeader.H"

void
uw1FaceAverages( LevelData<FluxBox>&         a_face_phi,
                 const LevelData<FArrayBox>& a_cell_phi,
                 const LevelData<FluxBox>&   a_face_vel,
                 const PhaseGeom&            a_geom );

void
uw3FaceAverages( LevelData<FluxBox>&         a_face_phi,
                 const LevelData<FArrayBox>& a_cell_phi,
                 const LevelData<FluxBox>&   a_face_vel,
                 const PhaseGeom&            a_geom );

void
uw5FaceAverages( LevelData<FluxBox>&           a_face_phi,
                 const LevelData<FArrayBox>&   a_cell_phi,
                 const LevelData<FluxBox>&     a_face_vel,
                 const PhaseGeom&              a_geom );

void
weno5FaceAverages( LevelData<FluxBox>&         a_face_phi,
                   const LevelData<FArrayBox>& a_cell_phi,
                   const LevelData<FluxBox>&   a_face_vel,
                   const PhaseGeom&            a_geom );

void
bwenoFaceAverages( LevelData<FluxBox>&         a_face_phi,
                   const LevelData<FArrayBox>& a_cell_phi,
                   const LevelData<FluxBox>&   a_face_vel,
                   const PhaseGeom&            a_geom);

#include "NamespaceFooter.H"

#endif
