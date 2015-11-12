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




/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

#ifndef ML_IFPACK_H
#define ML_IFPACK_H

#include "ml_include.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_TEUCHOS)

int ML_Smoother_Ifpack(ML_Smoother *sm,int inlen,double x[],int outlen,
			      double rhs[]);

void ML_Smoother_Clean_Ifpack(void * Ifpack_Handle);

/** Solves using IFPACK */
int ML_Ifpack_Solve(void * Ifpack_Handle, double * x, double * rhs);

/** Destroy all data associated to the IFPACK smoother. */
void ML_Ifpack_Destroy(void * Ifpack_Handle);

#endif

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef ML_IFPACK_H */
