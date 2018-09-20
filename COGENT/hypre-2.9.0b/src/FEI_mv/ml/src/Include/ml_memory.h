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
/* ML memory management functions                                       */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : April, 1998                                          */
/* ******************************************************************** */

#ifndef __MLMEM__
#define __MLMEM__

#include <stdio.h>
#include <stdlib.h>
#include "ml_common.h"

#define MAX_MALLOC_LOG 1000

#ifdef size_t
#define ml_size_t size_t
#else
#define ml_size_t int
#endif
#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif
#ifndef ML_MEM_CHECK
extern void ML_free(void *);
extern void *ML_allocate(ml_size_t size);
#endif
#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#ifndef ML_MEM_CHECK
#define ML_allocate(i)    malloc((i + sizeof(double) ))
#define ML_realloc(i,j)   realloc(i,j)
extern void *ml_void_mem_ptr;
#define ML_free(i)        { ml_void_mem_ptr = (void *) i;  if (ml_void_mem_ptr != NULL) {free( (void*) i); i = NULL;} }
#else
#define ML_free(i)        { ML_myfree(i); i = NULL; }
#define ML_realloc(i,j)    ML_myrealloc(i,j)
#endif
#define ML_allocate_check(ptr_to_check) \
                         {if ((ptr_to_check) == NULL) {\
                            printf("In file %s (line %d): memory allocation failed for pointer #%lu\n", __FILE__, __LINE__, (long unsigned int) ptr_to_check);\
                            exit(1);\
                            }\
                         }

#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif
  extern char *ML_memory_check(char *fmt, ...);


extern  int  ML_memory_alloc( void **, unsigned int, char * );
extern  int  ML_memory_free( void ** );
extern  int  ML_memory_check_var(void *);
extern  int  ML_memory_inquire(void);
extern  int  ML_memory_inquire_short( int );
extern  int  ML_memory_clean( char *, int );
#ifdef ML_MEM_CHECK
extern void  ML_print_it();
extern char *ML_allocate(unsigned int isize);
extern void  ML_myfree(void *vptr);
extern char *ML_myrealloc(void *vptr, unsigned int new_size);
extern void ML_spit_it_out();
#endif
extern int ML_MaxAllocatableSize();
extern int ML_MaxMemorySize();

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif

