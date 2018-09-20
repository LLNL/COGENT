/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision: 2.10 $
 ***********************************************************************EHEADER*/

#include "_hypre_struct_ls.h"

/*==========================================================================*/

HYPRE_Int
HYPRE_StructGMRESCreate( MPI_Comm comm, HYPRE_StructSolver *solver )
{
   hypre_GMRESFunctions * gmres_functions =
      hypre_GMRESFunctionsCreate(
         hypre_CAlloc, hypre_StructKrylovFree,
         hypre_StructKrylovCommInfo,
         hypre_StructKrylovCreateVector,
         hypre_StructKrylovCreateVectorArray,
         hypre_StructKrylovDestroyVector, hypre_StructKrylovMatvecCreate,
         hypre_StructKrylovMatvec, hypre_StructKrylovMatvecDestroy,
         hypre_StructKrylovInnerProd, hypre_StructKrylovCopyVector,
         hypre_StructKrylovClearVector,
         hypre_StructKrylovScaleVector, hypre_StructKrylovAxpy,
         hypre_StructKrylovIdentitySetup, hypre_StructKrylovIdentity );

   *solver = ( (HYPRE_StructSolver) hypre_GMRESCreate( gmres_functions ) );

   return hypre_error_flag;
}

/*==========================================================================*/

HYPRE_Int 
HYPRE_StructGMRESDestroy( HYPRE_StructSolver solver )
{
   return( hypre_GMRESDestroy( (void *) solver ) );
}

/*==========================================================================*/

HYPRE_Int 
HYPRE_StructGMRESSetup( HYPRE_StructSolver solver,
                        HYPRE_StructMatrix A,
                        HYPRE_StructVector b,
                        HYPRE_StructVector x      )
{
   return( HYPRE_GMRESSetup( (HYPRE_Solver) solver,
                             (HYPRE_Matrix) A,
                             (HYPRE_Vector) b,
                             (HYPRE_Vector) x ) );
}

/*==========================================================================*/

HYPRE_Int 
HYPRE_StructGMRESSolve( HYPRE_StructSolver solver,
                        HYPRE_StructMatrix A,
                        HYPRE_StructVector b,
                        HYPRE_StructVector x      )
{
   return( HYPRE_GMRESSolve( (HYPRE_Solver) solver,
                             (HYPRE_Matrix) A,
                             (HYPRE_Vector) b,
                             (HYPRE_Vector) x ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructGMRESSetTol( HYPRE_StructSolver solver,
                         double             tol    )
{
   return( HYPRE_GMRESSetTol( (HYPRE_Solver) solver, tol ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructGMRESSetAbsoluteTol( HYPRE_StructSolver solver,
                                 double             atol    )
{
   return( HYPRE_GMRESSetAbsoluteTol( (HYPRE_Solver) solver, atol ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructGMRESSetMaxIter( HYPRE_StructSolver solver,
                             HYPRE_Int          max_iter )
{
   return( HYPRE_GMRESSetMaxIter( (HYPRE_Solver) solver, max_iter ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructGMRESSetKDim( HYPRE_StructSolver solver,
                          HYPRE_Int          k_dim )
{
   return( HYPRE_GMRESSetKDim( (HYPRE_Solver) solver, k_dim ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructGMRESSetPrecond( HYPRE_StructSolver         solver,
                             HYPRE_PtrToStructSolverFcn precond,
                             HYPRE_PtrToStructSolverFcn precond_setup,
                             HYPRE_StructSolver         precond_solver )
{
   return( HYPRE_GMRESSetPrecond( (HYPRE_Solver) solver,
                                  (HYPRE_PtrToSolverFcn) precond,
                                  (HYPRE_PtrToSolverFcn) precond_setup,
                                  (HYPRE_Solver) precond_solver ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructGMRESSetLogging( HYPRE_StructSolver solver,
                             HYPRE_Int          logging )
{
   return( HYPRE_GMRESSetLogging( (HYPRE_Solver) solver, logging ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructGMRESSetPrintLevel( HYPRE_StructSolver solver,
                                HYPRE_Int          print_level )
{
   return( HYPRE_GMRESSetPrintLevel( (HYPRE_Solver) solver, print_level ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructGMRESGetNumIterations( HYPRE_StructSolver  solver,
                                   HYPRE_Int          *num_iterations )
{
   return( HYPRE_GMRESGetNumIterations( (HYPRE_Solver) solver, num_iterations ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructGMRESGetFinalRelativeResidualNorm( HYPRE_StructSolver  solver,
                                               double             *norm   )
{
   return( HYPRE_GMRESGetFinalRelativeResidualNorm( (HYPRE_Solver) solver, norm ) );
}


