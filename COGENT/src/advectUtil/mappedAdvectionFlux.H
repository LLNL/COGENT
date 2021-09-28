#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _MAPPEDADVECTIONFLUX_H_
#define _MAPPEDADVECTIONFLUX_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "FourthOrderCoordSys.H"
#include "mappedLimiter.H"

#include "NamespaceHeader.H"

/// compute high-order flux for advection equation on mapped grids.
/** Assumes all physical boundary conditions have been set.
    -- also returns face-averaged face advection velocity in faceVel
 */


void
cellUJToCellU(
              LevelData<FArrayBox>& a_cell_u,
              const LevelData<FArrayBox>& a_cell_uJ,
              const FourthOrderCoordSys* a_FOCS);


// original, SG-unfriendly version - implements the 
// "average of product is product of averages plus gradient correction"
// Is equivalent to the SG-friendly version if a_useFourthOrder = false
void
computeCompFaceFluxes( LevelData<FluxBox>& a_uTimesV,
                       const LevelData<FluxBox>& a_u,
                       const LevelData<FluxBox>& a_v,
                       bool a_useFourthOrder);

// New, SG-friendly version - implements the 
// "de-average to face centers, multiply, then re-average" method
// Is equivalent to the SG-unfriendly version if a_useFourthOrder = false
void
computeCompFaceFluxesSG( LevelData<FluxBox>& a_uTimesV,
                       const LevelData<FluxBox>& a_u,
                       const LevelData<FluxBox>& a_v,
                       bool a_useFourthOrder);

void
computeCompFaceFluxesGen( LevelData<FluxBox>& a_uTimesV,
                       const LevelData<FluxBox>& a_u,
                       const LevelData<FluxBox>& a_v,
                       bool a_useFourthOrder,
							  bool a_useSG);


void
computeCompFaceFluxesNormal( LevelData<FluxBox>& a_uTimesV,
                             const LevelData<FluxBox>& a_u,
                             const LevelData<FluxBox>& a_v,
                             bool a_useFourthOrder);

#include "NamespaceFooter.H"

#endif
