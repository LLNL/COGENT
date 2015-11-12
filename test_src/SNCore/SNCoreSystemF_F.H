#ifndef _SNCORESYSTEMF_F_H_
#define _SNCORESYSTEMF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDBLOB_INITIAL_CONDITION 
#define GUARDBLOB_INITIAL_CONDITION 
// Prototype for Fortran procedure blob_initial_condition ...
//
void FORTRAN_NAME( BLOB_INITIAL_CONDITION ,blob_initial_condition )(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_FRA1(data) );

#define FORT_BLOB_INITIAL_CONDITION FORTRAN_NAME( inlineBLOB_INITIAL_CONDITION, inlineBLOB_INITIAL_CONDITION)
#define FORTNT_BLOB_INITIAL_CONDITION FORTRAN_NAME( BLOB_INITIAL_CONDITION, blob_initial_condition)

inline void FORTRAN_NAME(inlineBLOB_INITIAL_CONDITION, inlineBLOB_INITIAL_CONDITION)(
      CHFp_BOX(gridbox)
      ,CHFp_CONST_REALVECT(h)
      ,CHFp_FRA1(data) )
{
 CH_TIMELEAF("FORT_BLOB_INITIAL_CONDITION");
 FORTRAN_NAME( BLOB_INITIAL_CONDITION ,blob_initial_condition )(
      CHFt_BOX(gridbox)
      ,CHFt_CONST_REALVECT(h)
      ,CHFt_FRA1(data) );
}
#endif  // GUARDBLOB_INITIAL_CONDITION 

}

#endif
