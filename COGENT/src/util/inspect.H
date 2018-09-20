#include <LevelData.H>
#include <FArrayBox.H>
#include <FluxBox.H>

#include "NamespaceHeader.H"

/*
  Functions to facilitate (except for the fact they need to be compiled
  in a priori) the inspection of FArrayBox, FluxBox and LevelData objects
  in TotalView.
*/

void inspect(const FArrayBox& data);
void inspect(const FluxBox& data);
void inspect(const LevelData<FArrayBox>& data);
void inspect(const LevelData<FluxBox>& data);

#include "NamespaceFooter.H"
