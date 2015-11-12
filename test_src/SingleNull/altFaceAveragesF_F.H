#ifndef _ALTFACEAVERAGESF_F_H_
#define _ALTFACEAVERAGESF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDBWENOFACEVALUES 
#define GUARDBWENOFACEVALUES 
// Prototype for Fortran procedure BWENOFACEVALUES ...
//
void FORTRAN_NAME( BWENOFACEVALUES ,bwenofacevalues )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_BWENOFACEVALUES FORTRAN_NAME( inlineBWENOFACEVALUES, inlineBWENOFACEVALUES)
#define FORTNT_BWENOFACEVALUES FORTRAN_NAME( BWENOFACEVALUES, bwenofacevalues)

inline void FORTRAN_NAME(inlineBWENOFACEVALUES, inlineBWENOFACEVALUES)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_BWENOFACEVALUES");
 FORTRAN_NAME( BWENOFACEVALUES ,bwenofacevalues )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDBWENOFACEVALUES 

#ifndef GUARDWENO5FACEVALUES 
#define GUARDWENO5FACEVALUES 
// Prototype for Fortran procedure WENO5FACEVALUES ...
//
void FORTRAN_NAME( WENO5FACEVALUES ,weno5facevalues )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_WENO5FACEVALUES FORTRAN_NAME( inlineWENO5FACEVALUES, inlineWENO5FACEVALUES)
#define FORTNT_WENO5FACEVALUES FORTRAN_NAME( WENO5FACEVALUES, weno5facevalues)

inline void FORTRAN_NAME(inlineWENO5FACEVALUES, inlineWENO5FACEVALUES)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_WENO5FACEVALUES");
 FORTRAN_NAME( WENO5FACEVALUES ,weno5facevalues )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDWENO5FACEVALUES 

#ifndef GUARDUW5FACEVALUES 
#define GUARDUW5FACEVALUES 
// Prototype for Fortran procedure UW5FACEVALUES ...
//
void FORTRAN_NAME( UW5FACEVALUES ,uw5facevalues )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) );

#define FORT_UW5FACEVALUES FORTRAN_NAME( inlineUW5FACEVALUES, inlineUW5FACEVALUES)
#define FORTNT_UW5FACEVALUES FORTRAN_NAME( UW5FACEVALUES, uw5facevalues)

inline void FORTRAN_NAME(inlineUW5FACEVALUES, inlineUW5FACEVALUES)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_UW5FACEVALUES");
 FORTRAN_NAME( UW5FACEVALUES ,uw5facevalues )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDUW5FACEVALUES 

}

#endif
