#ifndef _ALTFACEAVERAGESF_F_H_
#define _ALTFACEAVERAGESF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDALTFACEVALUES 
#define GUARDALTFACEVALUES 
// Prototype for Fortran procedure ALTFACEVALUES ...
//
void FORTRAN_NAME( ALTFACEVALUES ,altfacevalues )(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_REAL(dx) );

#define FORT_ALTFACEVALUES FORTRAN_NAME( inlineALTFACEVALUES, inlineALTFACEVALUES)
#define FORTNT_ALTFACEVALUES FORTRAN_NAME( ALTFACEVALUES, altfacevalues)

inline void FORTRAN_NAME(inlineALTFACEVALUES, inlineALTFACEVALUES)(
      CHFp_FRA(facePhi)
      ,CHFp_CONST_FRA(cellPhi)
      ,CHFp_CONST_FRA1(faceVel)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_ALTFACEVALUES");
 FORTRAN_NAME( ALTFACEVALUES ,altfacevalues )(
      CHFt_FRA(facePhi)
      ,CHFt_CONST_FRA(cellPhi)
      ,CHFt_CONST_FRA1(faceVel)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDALTFACEVALUES 

}

#endif
