#ifndef _CHKPT_H_
#define _CHKPT_H_

#include "NamespaceHeader.H"

typedef enum {
  chkpt_pre_time_step,  /* pre-timestep */
  chkpt_ti_advance_1,   /* in time step, before stage calculation */
  chkpt_ti_advance_2,   /* in time step, after stage calculation */
  chkpt_ti_advance_3,   /* in time step, before step completion */
  chkpt_stage_func_0,   /* stage function evaluation for first Newton iteration */
  chkpt_stage_func_n,   /* stage function evaluation for Newton iteration > 1 */
  chkpt_stage_jac       /* stage Jacobian evaluation */
} Checkpoint;


#include "NamespaceFooter.H"

#endif
