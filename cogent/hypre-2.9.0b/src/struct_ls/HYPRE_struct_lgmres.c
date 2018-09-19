/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision: 2.6 $
 ***********************************************************************EHEADER*/

#include "_hypre_struct_ls.h"

/*==========================================================================*/

HYPRE_Int
HYPRE_StructLGMRESCreate( MPI_Comm comm, HYPRE_StructSolver *solver )
{
   hypre_LGMRESFunctions * lgmres_functions =
      hypre_LGMRESFunctionsCreate(
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

   *solver = ( (HYPRE_StructSolver) hypre_LGMRESCreate( lgmres_functions ) );

   return hypre_error_flag;
}

/*==========================================================================*/

HYPRE_Int 
HYPRE_StructLGMRESDestroy( HYPRE_StructSolver solver )
{
   return( hypre_LGMRESDestroy( (void *) solver ) );
}

/*==========================================================================*/

HYPRE_Int 
HYPRE_StructLGMRESSetup( HYPRE_StructSolver solver,
                         HYPRE_StructMatrix A,
                         HYPRE_StructVector b,
                         HYPRE_StructVector x      )
{
   return( HYPRE_LGMRESSetup( (HYPRE_Solver) solver,
                              (HYPRE_Matrix) A,
                              (HYPRE_Vector) b,
                              (HYPRE_Vector) x ) );
}

/*==========================================================================*/

HYPRE_Int 
HYPRE_StructLGMRESSolve( HYPRE_StructSolver solver,
                         HYPRE_StructMatrix A,
                         HYPRE_StructVector b,
                         HYPRE_StructVector x      )
{
   return( HYPRE_LGMRESSolve( (HYPRE_Solver) solver,
                              (HYPRE_Matrix) A,
                              (HYPRE_Vector) b,
                              (HYPRE_Vector) x ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructLGMRESSetTol( HYPRE_StructSolver solver,
                          double             tol    )
{
   return( HYPRE_LGMRESSetTol( (HYPRE_Solver) solver, tol ) );
}
/*==========================================================================*/

HYPRE_Int
HYPRE_StructLGMRESSetAbsoluteTol( HYPRE_StructSolver solver,
                                  double             tol    )
{
   return( HYPRE_LGMRESSetAbsoluteTol( (HYPRE_Solver) solver, tol ) );
}
/*==========================================================================*/

HYPRE_Int
HYPRE_StructLGMRESSetMaxIter( HYPRE_StructSolver solver,
                              HYPRE_Int          max_iter )
{
   return( HYPRE_LGMRESSetMaxIter( (HYPRE_Solver) solver, max_iter ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructLGMRESSetKDim( HYPRE_StructSolver solver,
                           HYPRE_Int          k_dim )
{
   return( HYPRE_LGMRESSetKDim( (HYPRE_Solver) solver, k_dim ) );
}



/*==========================================================================*/

HYPRE_Int
HYPRE_StructLGMRESSetAugDim( HYPRE_StructSolver solver,
                             HYPRE_Int          aug_dim )
{
   return( HYPRE_LGMRESSetAugDim( (HYPRE_Solver) solver, aug_dim ) );
}


/*==========================================================================*/

HYPRE_Int
HYPRE_StructLGMRESSetPrecond( HYPRE_StructSolver         solver,
                              HYPRE_PtrToStructSolverFcn precond,
                              HYPRE_PtrToStructSolverFcn precond_setup,
                              HYPRE_StructSolver         precond_solver )
{
   return( HYPRE_LGMRESSetPrecond( (HYPRE_Solver) solver,
                                   (HYPRE_PtrToSolverFcn) precond,
                                   (HYPRE_PtrToSolverFcn) precond_setup,
                                   (HYPRE_Solver) precond_solver ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructLGMRESSetLogging( HYPRE_StructSolver solver,
                              HYPRE_Int          logging )
{
   return( HYPRE_LGMRESSetLogging( (HYPRE_Solver) solver, logging ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructLGMRESSetPrintLevel( HYPRE_StructSolver solver,
                                 HYPRE_Int          print_level )
{
   return( HYPRE_LGMRESSetPrintLevel( (HYPRE_Solver) solver, print_level ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructLGMRESGetNumIterations( HYPRE_StructSolver  solver,
                                    HYPRE_Int          *num_iterations )
{
   return( HYPRE_LGMRESGetNumIterations( (HYPRE_Solver) solver, num_iterations ) );
}

/*==========================================================================*/

HYPRE_Int
HYPRE_StructLGMRESGetFinalRelativeResidualNorm( HYPRE_StructSolver  solver,
                                                double             *norm   )
{
   return( HYPRE_LGMRESGetFinalRelativeResidualNorm( (HYPRE_Solver) solver,
                                                     norm ) );
}


