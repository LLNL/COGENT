#ifndef RECTANGULARTORUSDIRECTIONS_H
#define RECTANGULARTORUSDIRECTIONS_H
#include "Directions.H"

#if PDIM==4
#define R_DIR       RADIAL_DIR
#define Z_DIR       POLOIDAL_DIR
#elif PDIM==5
#define R_DIR       RADIAL_DIR
#define PHI_DIR     POLOIDAL_DIR
#define Z_DIR       TOROIDAL_DIR
#endif

#endif


