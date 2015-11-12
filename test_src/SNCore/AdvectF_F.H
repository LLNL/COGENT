#ifndef _ADVECTF_F_H_
#define _ADVECTF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDVEL_TEST 
#define GUARDVEL_TEST 
// Prototype for Fortran procedure vel_test ...
//
void FORTRAN_NAME( VEL_TEST ,vel_test )(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_REALVECT(vvec)
      ,CHFp_FRA(velocity) );

#define FORT_VEL_TEST FORTRAN_NAME( inlineVEL_TEST, inlineVEL_TEST)
#define FORTNT_VEL_TEST FORTRAN_NAME( VEL_TEST, vel_test)

inline void FORTRAN_NAME(inlineVEL_TEST, inlineVEL_TEST)(
      CHFp_CONST_INT(dir)
      ,CHFp_BOX(gridbox)
      ,CHFp_REALVECT(vvec)
      ,CHFp_FRA(velocity) )
{
 CH_TIMELEAF("FORT_VEL_TEST");
 FORTRAN_NAME( VEL_TEST ,vel_test )(
      CHFt_CONST_INT(dir)
      ,CHFt_BOX(gridbox)
      ,CHFt_REALVECT(vvec)
      ,CHFt_FRA(velocity) );
}
#endif  // GUARDVEL_TEST 

#ifndef GUARDMULT_CFG_NT 
#define GUARDMULT_CFG_NT 
// Prototype for Fortran procedure mult_cfg_nt ...
//
void FORTRAN_NAME( MULT_CFG_NT ,mult_cfg_nt )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(N)
      ,CHFp_FRA(velocity) );

#define FORT_MULT_CFG_NT FORTRAN_NAME( inlineMULT_CFG_NT, inlineMULT_CFG_NT)
#define FORTNT_MULT_CFG_NT FORTRAN_NAME( MULT_CFG_NT, mult_cfg_nt)

inline void FORTRAN_NAME(inlineMULT_CFG_NT, inlineMULT_CFG_NT)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_FRA(N)
      ,CHFp_FRA(velocity) )
{
 CH_TIMELEAF("FORT_MULT_CFG_NT");
 FORTRAN_NAME( MULT_CFG_NT ,mult_cfg_nt )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_FRA(N)
      ,CHFt_FRA(velocity) );
}
#endif  // GUARDMULT_CFG_NT 

#ifndef GUARDCFGDIVERGENCE 
#define GUARDCFGDIVERGENCE 
// Prototype for Fortran procedure cfgdivergence ...
//
void FORTRAN_NAME( CFGDIVERGENCE ,cfgdivergence )(
      CHFp_CONST_FRA(uEdge)
      ,CHFp_FRA(div)
      ,CHFp_BOX(gridInt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(dir) );

#define FORT_CFGDIVERGENCE FORTRAN_NAME( inlineCFGDIVERGENCE, inlineCFGDIVERGENCE)
#define FORTNT_CFGDIVERGENCE FORTRAN_NAME( CFGDIVERGENCE, cfgdivergence)

inline void FORTRAN_NAME(inlineCFGDIVERGENCE, inlineCFGDIVERGENCE)(
      CHFp_CONST_FRA(uEdge)
      ,CHFp_FRA(div)
      ,CHFp_BOX(gridInt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(dir) )
{
 CH_TIMELEAF("FORT_CFGDIVERGENCE");
 FORTRAN_NAME( CFGDIVERGENCE ,cfgdivergence )(
      CHFt_CONST_FRA(uEdge)
      ,CHFt_FRA(div)
      ,CHFt_BOX(gridInt)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_INT(dir) );
}
#endif  // GUARDCFGDIVERGENCE 

}

#endif
