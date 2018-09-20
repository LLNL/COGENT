/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision: 1.1 $
 ***********************************************************************EHEADER*/

#ifndef __MLI_SOLVER_KACZMARZ_H__
#define __MLI_SOLVER_KACZMARZ_H__

#include <stdio.h>
#include "_hypre_parcsr_mv.h"
#include "mli_matrix.h"
#include "mli_vector.h"
#include "mli_solver.h"

/******************************************************************************
 * data structure for the Kaczmarz relaxation scheme
 *---------------------------------------------------------------------------*/

class MLI_Solver_Kaczmarz : public MLI_Solver
{
   MLI_Matrix  *Amat_;
   int         nSweeps_;
   double      *AsqDiag_;
   int         zeroInitialGuess_;

public :

   MLI_Solver_Kaczmarz(char *name);
   ~MLI_Solver_Kaczmarz();
   int setup(MLI_Matrix *Amat);
   int solve(MLI_Vector *f, MLI_Vector *u);

   int setParams(char *paramString, int argc, char **argv);
};

#endif

