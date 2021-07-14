#ifndef _TIDEFINITIONS_H_
#define _TIDEFINITIONS_H_

#include "NamespaceHeader.H"

typedef enum {ti_explicit, ti_imex, ti_implicit} TimeIntegratorType;

typedef enum {_explicit_op_, _implicit_op_} tiOpType;

typedef enum {
  chkpt_naive,          /* no checkpoint-based optimizations */
  chkpt_pre_time_step,  /* pre-timestep */
  chkpt_ti_init,        /* initializeing the time integrators */
  chkpt_ti_advance_1,   /* in time step, before stage calculation */
  chkpt_ti_advance_2,   /* in time step, after stage calculation */
  chkpt_ti_advance_3,   /* in time step, before step completion */
  chkpt_stage_func_0,   /* stage function evaluation for first Newton iteration */
  chkpt_stage_func_n,   /* stage function evaluation for Newton iteration > 1 */
  chkpt_stage_jac       /* stage Jacobian evaluation */
} Checkpoint;

#include "NamespaceFooter.H"
#endif
