#ifndef _DIRECTIONS_H_
#define _DIRECTIONS_H_

#include "Dimensions.H"
#include "SPACE.H"

#if PDIM==2
enum directions {RADIAL_DIR, POLOIDAL_DIR};
#endif
#if PDIM==4
enum directions {RADIAL_DIR, POLOIDAL_DIR, VPARALLEL_DIR, MU_DIR};
#endif
#if PDIM==5
enum directions {RADIAL_DIR, POLOIDAL_DIR, TOROIDAL_DIR, VPARALLEL_DIR, MU_DIR};
#endif

#endif
