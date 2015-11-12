#ifndef _DIMENSIONS_H_
#define _DIMENSIONS_H_

#if CFG_DIM==2
#define CFG_NAMESPACE D2
#define CP1_NAMESPACE D3
#define PS_NAMESPACE D4
#define PDIM 4
#endif
#if CFG_DIM==3
#define CFG_NAMESPACE D3
#define CP1_NAMESPACE D4
#define PS_NAMESPACE D5
#define PDIM 5
#endif
#define VEL_DIM 2
#define VEL_NAMESPACE D2

#endif
