/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision: 1.5 $
 ***********************************************************************EHEADER*/




/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the New stuff                                         */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __MLMATRIX__
#define __MLMATRIX__

#include "ml_common.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


extern int oldML_Mdfy_Prolongator_DirBdry(ML *, int , double *, double *);
extern int ML_Compute_Coarse_Bdry(ML *ml_handle, int level, int size, 
           int fine_size);
extern int ML_Mdfy_Prolongator_DirBdry(ML *ml_handle, int level2, int size,
     int fine_size );

extern ML_Operator *ML_Operator_ExplicitlyScale(ML_Operator *matrix,
					 double scalar);

extern int ML_Operator_ChangeToSinglePrecision(ML_Operator *matrix);
extern int ML_Operator_ChangeToChar(ML_Operator *matrix);
extern int ML_Operator_ImplicitTranspose(ML_Operator *Rmat, 
					 ML_Operator *Pmat,
					 int PostCommAlreadySet);
extern int ML_Gen_Restrictor_TransP(ML *, int, int);
extern int ML_Gen_Prolongator_Getrow(ML *, int , int , int , int ,
            int (*)(void* , int , int *, int , int *, double *, int *),
            int (*)(double *, void*), void *data, int);
  /*
extern int ML_Operator_Transpose(ML_Operator *Amat, ML_Operator *Amat_trans );
  */
extern int ML_Operator_ColPartition2RowPartition(ML_Operator *A, ML_Operator *Atrans);
extern int eye_getrows(ML_Operator *data, int N_requested_rows, int requested_rows[],
		       int allocated_space, int columns[], double values[],
		       int row_lengths[]);
extern	int eye_matvec(ML_Operator *Amat_in, int ilen, double p[], int olen, double ap[]);
extern int ML_Operator_Transpose_byrow(ML_Operator *A, ML_Operator *Atrans);
extern int ML_Operator_Getrow_Diag(ML_Operator *Amat, double **diagonal);
extern int ML_build_overlapped_pre_comm(ML_Operator *tempA, ML_CommInfoOP
					*old_comm, int max_per_proc,
					int *hash_list, int hash_length, 
					int *hash_used, int old_Nrows, 
					int *Nexternal, int *external[],
					int *Nexternal_allocated);
extern int ML_Operator_HashGlobalRcvList(ML_CommInfoOP *pre_comm, int Nrows, 
					 int hash_list[], int hash_length, 
					 int *hash_used, int Gid_assigned_to_proc[], 
					 ML_Comm *comm, 
					 int *Nexternal, int **external,
					 int *Nexternal_allocated);
extern int ML_overlap(ML_Operator *oldA, ML_Operator *newA, int overlap,
		      ML_CommInfoOP **nonOverlapped_2_Overlapped);
extern void ML_Operator_ReportStatistics(ML_Operator *mat, char *appendlabel,
                                         int perfAndCommStats);
extern void ML_Operator_Profile(ML_Operator *A, char *appendlabel, int numits);
extern void ML_Operator_Profile_SetIterations(int numits);
extern int ML_Operator_Profile_GetIterations();
extern int ML_Operator_Get_Nnz(ML_Operator *A);

#ifndef ML_CPP
#ifdef __cplusplus
  }
#endif
#endif

#endif
