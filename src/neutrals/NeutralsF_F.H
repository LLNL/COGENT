#ifndef _NEUTRALSF_F_H_
#define _NEUTRALSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_SC_NTR_FREQ 
#define GUARDCOMPUTE_SC_NTR_FREQ 
// Prototype for Fortran procedure compute_sc_ntr_freq ...
//
void FORTRAN_NAME( COMPUTE_SC_NTR_FREQ ,compute_sc_ntr_freq )(
      CHFp_BOX(box)
      ,CHFp_FRA1(ntr_freq)
      ,CHFp_CONST_FRA1(n)
      ,CHFp_CONST_FRA1(T) );

#define FORT_COMPUTE_SC_NTR_FREQ FORTRAN_NAME( inlineCOMPUTE_SC_NTR_FREQ, inlineCOMPUTE_SC_NTR_FREQ)
#define FORTNT_COMPUTE_SC_NTR_FREQ FORTRAN_NAME( COMPUTE_SC_NTR_FREQ, compute_sc_ntr_freq)

inline void FORTRAN_NAME(inlineCOMPUTE_SC_NTR_FREQ, inlineCOMPUTE_SC_NTR_FREQ)(
      CHFp_BOX(box)
      ,CHFp_FRA1(ntr_freq)
      ,CHFp_CONST_FRA1(n)
      ,CHFp_CONST_FRA1(T) )
{
 CH_TIMELEAF("FORT_COMPUTE_SC_NTR_FREQ");
 FORTRAN_NAME( COMPUTE_SC_NTR_FREQ ,compute_sc_ntr_freq )(
      CHFt_BOX(box)
      ,CHFt_FRA1(ntr_freq)
      ,CHFt_CONST_FRA1(n)
      ,CHFt_CONST_FRA1(T) );
}
#endif  // GUARDCOMPUTE_SC_NTR_FREQ 

}

#endif
