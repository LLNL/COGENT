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
uw1FaceAverages( LevelData<FluxBox>&         face_phi,
                 const LevelData<FArrayBox>& cell_phi,
                 const LevelData<FluxBox>&   face_vel,
                 const PhaseGeom&            geom,
                 const bool                  mapped_velocity );

void
uw3FaceAverages( LevelData<FluxBox>&         face_phi,
                 const LevelData<FArrayBox>& cell_phi,
                 const LevelData<FluxBox>&   face_vel,
                 const PhaseGeom&            geom,
                 const bool                  mapped_velocity );

void
uw5FaceAverages( LevelData<FluxBox>&           face_phi,
                 const LevelData<FArrayBox>&   cell_phi,
                 const LevelData<FluxBox>&     face_vel,
                 const PhaseGeom&              geom,
                 const bool                    mapped_velocity );

void
weno5FaceAverages( LevelData<FluxBox>&         face_phi,
                   const LevelData<FArrayBox>& cell_phi,
                   const LevelData<FluxBox>&   face_vel,
                   const PhaseGeom&            geom,
                   const bool                  mapped_velocity );

void
bwenoFaceAverages( LevelData<FluxBox>&         face_phi,
                   const LevelData<FArrayBox>& cell_phi,
                   const LevelData<FluxBox>&   face_vel,
                   const PhaseGeom&            geom,
                   const bool                  mapped_velocity );

#include "NamespaceFooter.H"

#endif
