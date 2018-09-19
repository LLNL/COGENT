/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision: 2.19 $
 ***********************************************************************EHEADER*/

/******************************************************************************
 *
 * Header info for the Box structures
 *
 *****************************************************************************/

#ifndef hypre_BOX_HEADER
#define hypre_BOX_HEADER

#ifndef HYPRE_MAXDIM
#define HYPRE_MAXDIM 3
#endif

/*--------------------------------------------------------------------------
 * hypre_Index:
 *   This is used to define indices in index space, or dimension
 *   sizes of boxes.
 *
 *   The spatial dimensions x, y, and z may be specified by the
 *   integers 0, 1, and 2, respectively (see the hypre_IndexD macro below).
 *   This simplifies the code in the hypre_Box class by reducing code
 *   replication.
 *--------------------------------------------------------------------------*/

typedef HYPRE_Int  hypre_Index[3];
typedef HYPRE_Int *hypre_IndexRef;

/*--------------------------------------------------------------------------
 * hypre_Box:
 *--------------------------------------------------------------------------*/

typedef struct hypre_Box_struct
{
   hypre_Index imin;           /* min bounding indices */
   hypre_Index imax;           /* max bounding indices */

} hypre_Box;

/*--------------------------------------------------------------------------
 * hypre_BoxArray:
 *   An array of boxes.
 *--------------------------------------------------------------------------*/

typedef struct hypre_BoxArray_struct
{
   hypre_Box  *boxes;         /* Array of boxes */
   HYPRE_Int   size;          /* Size of box array */
   HYPRE_Int   alloc_size;    /* Size of currently alloced space */

} hypre_BoxArray;

#define hypre_BoxArrayExcess 10

/*--------------------------------------------------------------------------
 * hypre_BoxArrayArray:
 *   An array of box arrays.
 *--------------------------------------------------------------------------*/

typedef struct hypre_BoxArrayArray_struct
{
   hypre_BoxArray  **box_arrays;    /* Array of pointers to box arrays */
   HYPRE_Int         size;          /* Size of box array array */

} hypre_BoxArrayArray;


/*--------------------------------------------------------------------------
 * Accessor macros: hypre_Index
 *--------------------------------------------------------------------------*/

#define hypre_IndexD(index, d)  (index[d])

#define hypre_IndexX(index)     hypre_IndexD(index, 0)
#define hypre_IndexY(index)     hypre_IndexD(index, 1)
#define hypre_IndexZ(index)     hypre_IndexD(index, 2)

/*--------------------------------------------------------------------------
 * Member functions: hypre_Index
 *--------------------------------------------------------------------------*/

#define hypre_SetIndex(index, ix, iy, iz) \
( hypre_IndexX(index) = ix,\
  hypre_IndexY(index) = iy,\
  hypre_IndexZ(index) = iz )

#define hypre_ClearIndex(index)  hypre_SetIndex(index, 0, 0, 0)

#define hypre_IndexZero(index)\
   (hypre_IndexX(index) == 0 &&  hypre_IndexY(index) == 0 \
    && hypre_IndexZ(index) == 0)

#define hypre_IndexGTESize(index, size) \
(hypre_IndexX(index) >= size &&  hypre_IndexY(index) >= size \
    && hypre_IndexZ(index) >= size)

#define hypre_CopyIndex(index1, index2) \
( hypre_IndexX(index2) = hypre_IndexX(index1),\
  hypre_IndexY(index2) = hypre_IndexY(index1),\
  hypre_IndexZ(index2) = hypre_IndexZ(index1) )

#define hypre_CopyToCleanIndex(in_index, ndim, out_index) \
{\
   HYPRE_Int d;\
   for (d = 0; d < ndim; d++)\
   {\
      hypre_IndexD(out_index, d) = hypre_IndexD(in_index, d);\
   }\
   for (d = ndim; d < 3; d++)\
   {\
      hypre_IndexD(out_index, d) = 0;\
   }\
}

/*--------------------------------------------------------------------------
 * Accessor macros: hypre_Box
 *--------------------------------------------------------------------------*/

#define hypre_BoxIMin(box)     ((box) -> imin)
#define hypre_BoxIMax(box)     ((box) -> imax)

#define hypre_AddIndex(index1, index2, index3) \
( hypre_IndexX(index3) = hypre_IndexX(index2) + hypre_IndexX(index1),\
  hypre_IndexY(index3) = hypre_IndexY(index2) + hypre_IndexY(index1),\
  hypre_IndexZ(index3) = hypre_IndexZ(index2) + hypre_IndexZ(index1) )

#define hypre_SubtractIndex(index1, index2, index3) \
( hypre_IndexX(index3) = hypre_IndexX(index1) - hypre_IndexX(index2),\
  hypre_IndexY(index3) = hypre_IndexY(index1) - hypre_IndexY(index2),\
  hypre_IndexZ(index3) = hypre_IndexZ(index1) - hypre_IndexZ(index2) )

#define hypre_BoxIMinD(box, d) (hypre_IndexD(hypre_BoxIMin(box), d))
#define hypre_BoxIMaxD(box, d) (hypre_IndexD(hypre_BoxIMax(box), d))
#define hypre_BoxSizeD(box, d) \
hypre_max(0, (hypre_BoxIMaxD(box, d) - hypre_BoxIMinD(box, d) + 1))

#define hypre_BoxIMinX(box)    hypre_BoxIMinD(box, 0)
#define hypre_BoxIMinY(box)    hypre_BoxIMinD(box, 1)
#define hypre_BoxIMinZ(box)    hypre_BoxIMinD(box, 2)

#define hypre_BoxIMaxX(box)    hypre_BoxIMaxD(box, 0)
#define hypre_BoxIMaxY(box)    hypre_BoxIMaxD(box, 1)
#define hypre_BoxIMaxZ(box)    hypre_BoxIMaxD(box, 2)

#define hypre_BoxSizeX(box)    hypre_BoxSizeD(box, 0)
#define hypre_BoxSizeY(box)    hypre_BoxSizeD(box, 1)
#define hypre_BoxSizeZ(box)    hypre_BoxSizeD(box, 2)

#define hypre_BoxEqualP( box1, box2 ) (\
 hypre_BoxIMinX(box1)==hypre_BoxIMinX(box2) &&\
 hypre_BoxIMaxX(box1)==hypre_BoxIMaxX(box2) &&\
 hypre_BoxIMinY(box1)==hypre_BoxIMinY(box2) &&\
 hypre_BoxIMaxY(box1)==hypre_BoxIMaxY(box2) &&\
 hypre_BoxIMinZ(box1)==hypre_BoxIMinZ(box2) &&\
 hypre_BoxIMaxZ(box1)==hypre_BoxIMaxZ(box2) )

#define hypre_IndexInBoxP( index, box ) (\
 hypre_IndexX(index)>=hypre_BoxIMinX(box) &&\
 hypre_IndexX(index)<=hypre_BoxIMaxX(box) &&\
 hypre_IndexY(index)>=hypre_BoxIMinY(box) &&\
 hypre_IndexY(index)<=hypre_BoxIMaxY(box) &&\
 hypre_IndexZ(index)>=hypre_BoxIMinZ(box) &&\
 hypre_IndexZ(index)<=hypre_BoxIMaxZ(box) )


#define hypre_IndexDInBoxP( index, d, box ) (\
 hypre_IndexD(index, d)>=hypre_BoxIMinD(box, d) &&\
 hypre_IndexD(index, d)<=hypre_BoxIMaxD(box, d) )

#define hypre_CopyBox(box1, box2) \
( hypre_CopyIndex(hypre_BoxIMin(box1), hypre_BoxIMin(box2)),\
  hypre_CopyIndex(hypre_BoxIMax(box1), hypre_BoxIMax(box2)) )

#define hypre_BoxVolume(box) \
(hypre_BoxSizeX(box) * hypre_BoxSizeY(box) * hypre_BoxSizeZ(box))

#define hypre_BoxShiftPos(box, shift) \
{\
   hypre_BoxIMinX(box) += hypre_IndexX(shift);\
   hypre_BoxIMinY(box) += hypre_IndexY(shift);\
   hypre_BoxIMinZ(box) += hypre_IndexZ(shift);\
   hypre_BoxIMaxX(box) += hypre_IndexX(shift);\
   hypre_BoxIMaxY(box) += hypre_IndexY(shift);\
   hypre_BoxIMaxZ(box) += hypre_IndexZ(shift);\
}

#define hypre_BoxShiftNeg(box, shift) \
{\
   hypre_BoxIMinX(box) -= hypre_IndexX(shift);\
   hypre_BoxIMinY(box) -= hypre_IndexY(shift);\
   hypre_BoxIMinZ(box) -= hypre_IndexZ(shift);\
   hypre_BoxIMaxX(box) -= hypre_IndexX(shift);\
   hypre_BoxIMaxY(box) -= hypre_IndexY(shift);\
   hypre_BoxIMaxZ(box) -= hypre_IndexZ(shift);\
}

#define hypre_BoxIndexRank(box, index) \
((hypre_IndexX(index) - hypre_BoxIMinX(box)) + \
 ((hypre_IndexY(index) - hypre_BoxIMinY(box)) + \
   ((hypre_IndexZ(index) - hypre_BoxIMinZ(box)) * \
    hypre_BoxSizeY(box))) * \
  hypre_BoxSizeX(box))

/* The first hypre_CCBoxIndexRank is better style because it keeps
   its similarity to the variable coefficient hypre_BoxIndexRank.
   The second one sometimes avoids compiler warnings...*/
#define hypre_CCBoxIndexRank(box, index) 0
#define hypre_CCBoxIndexRank_noargs() 0

#define hypre_BoxOffsetDistance(box, index) \
(hypre_IndexX(index) + \
 (hypre_IndexY(index) + \
  (hypre_IndexZ(index) * \
   hypre_BoxSizeY(box))) * \
 hypre_BoxSizeX(box))

#define hypre_CCBoxOffsetDistance(box, index) 0
  
/*--------------------------------------------------------------------------
 * Accessor macros: hypre_BoxArray
 *--------------------------------------------------------------------------*/

#define hypre_BoxArrayBoxes(box_array)     ((box_array) -> boxes)
#define hypre_BoxArrayBox(box_array, i)    &((box_array) -> boxes[(i)])
#define hypre_BoxArraySize(box_array)      ((box_array) -> size)
#define hypre_BoxArrayAllocSize(box_array) ((box_array) -> alloc_size)

/*--------------------------------------------------------------------------
 * Accessor macros: hypre_BoxArrayArray
 *--------------------------------------------------------------------------*/

#define hypre_BoxArrayArrayBoxArrays(box_array_array) \
((box_array_array) -> box_arrays)
#define hypre_BoxArrayArrayBoxArray(box_array_array, i) \
((box_array_array) -> box_arrays[(i)])
#define hypre_BoxArrayArraySize(box_array_array) \
((box_array_array) -> size)

/*--------------------------------------------------------------------------
 * Looping macros:
 *--------------------------------------------------------------------------*/

#define hypre_ForBoxI(i, box_array) \
for (i = 0; i < hypre_BoxArraySize(box_array); i++)

#define hypre_ForBoxArrayI(i, box_array_array) \
for (i = 0; i < hypre_BoxArrayArraySize(box_array_array); i++)

/*--------------------------------------------------------------------------
 * BoxLoop macros:
 *
 * NOTE: PThreads version of BoxLoop looping macros are in `box_pthreads.h'.
 *--------------------------------------------------------------------------*/

#ifndef HYPRE_USE_PTHREADS

#define HYPRE_BOX_PRIVATE hypre__nx,hypre__ny,hypre__nz,hypre__i,hypre__j,hypre__k

#define hypre_BoxLoopDeclareS(dbox, stride, sx, sy, sz) \
HYPRE_Int  sx = (hypre_IndexX(stride));\
HYPRE_Int  sy = (hypre_IndexY(stride)*hypre_BoxSizeX(dbox));\
HYPRE_Int  sz = (hypre_IndexZ(stride)*\
           hypre_BoxSizeX(dbox)*hypre_BoxSizeY(dbox))

#define hypre_BoxLoopDeclareN(loop_size) \
HYPRE_Int  hypre__i, hypre__j, hypre__k;\
HYPRE_Int  hypre__nx = hypre_IndexX(loop_size);\
HYPRE_Int  hypre__ny = hypre_IndexY(loop_size);\
HYPRE_Int  hypre__nz = hypre_IndexZ(loop_size);\
HYPRE_Int  hypre__mx = hypre__nx;\
HYPRE_Int  hypre__my = hypre__ny;\
HYPRE_Int  hypre__mz = hypre__nz;\
HYPRE_Int  hypre__dir, hypre__max;\
HYPRE_Int  hypre__div, hypre__mod;\
HYPRE_Int  hypre__block, hypre__num_blocks;\
hypre__dir = 0;\
hypre__max = hypre__nx;\
if (hypre__ny > hypre__max)\
{\
   hypre__dir = 1;\
   hypre__max = hypre__ny;\
}\
if (hypre__nz > hypre__max)\
{\
   hypre__dir = 2;\
   hypre__max = hypre__nz;\
}\
hypre__num_blocks = hypre_NumThreads();\
if (hypre__max < hypre__num_blocks)\
{\
   hypre__num_blocks = hypre__max;\
}\
if (hypre__num_blocks > 0)\
{\
   hypre__div = hypre__max / hypre__num_blocks;\
   hypre__mod = hypre__max % hypre__num_blocks;\
}

#define hypre_BoxLoopSet(i, j, k) \
i = 0;\
j = 0;\
k = 0;\
hypre__nx = hypre__mx;\
hypre__ny = hypre__my;\
hypre__nz = hypre__mz;\
if (hypre__num_blocks > 1)\
{\
   if (hypre__dir == 0)\
   {\
      i = hypre__block * hypre__div + hypre_min(hypre__mod, hypre__block);\
      hypre__nx = hypre__div + ((hypre__mod > hypre__block) ? 1 : 0);\
   }\
   else if (hypre__dir == 1)\
   {\
      j = hypre__block * hypre__div + hypre_min(hypre__mod, hypre__block);\
      hypre__ny = hypre__div + ((hypre__mod > hypre__block) ? 1 : 0);\
   }\
   else if (hypre__dir == 2)\
   {\
      k = hypre__block * hypre__div + hypre_min(hypre__mod, hypre__block);\
      hypre__nz = hypre__div + ((hypre__mod > hypre__block) ? 1 : 0);\
   }\
}

#define hypre_BoxLoopGetIndex(index) \
index[0] = hypre__i; index[1] = hypre__j; index[2] = hypre__k

/* Use this before the For macros below to force only one block */
#define hypre_BoxLoopSetOneBlock() hypre__num_blocks = 1

/* Use this to get the block iteration inside a BoxLoop */
#define hypre_BoxLoopBlock() hypre__block

/*-----------------------------------*/

#define hypre_BoxLoop0Begin(ndim, loop_size)\
{\
   hypre_BoxLoopDeclareN(loop_size);

#define hypre_BoxLoop0For()\
   hypre__BoxLoop0For(hypre__i, hypre__j, hypre__k)
#define hypre__BoxLoop0For(i, j, k)\
   for (hypre__block = 0; hypre__block < hypre__num_blocks; hypre__block++)\
   {\
   hypre_BoxLoopSet(i, j, k);\
   for (k = 0; k < hypre__nz; k++)\
   {\
      for (j = 0; j < hypre__ny; j++)\
      {\
         for (i = 0; i < hypre__nx; i++)\
         {

#define hypre_BoxLoop0End()\
         }\
      }\
   }\
   }\
}
  
/*-----------------------------------*/

#define hypre_BoxLoop1Begin(ndim, loop_size,\
                            dbox1, start1, stride1, i1)\
{\
   HYPRE_Int  hypre__i1start = hypre_BoxIndexRank(dbox1, start1);\
   hypre_BoxLoopDeclareS(dbox1, stride1, hypre__sx1, hypre__sy1, hypre__sz1);\
   hypre_BoxLoopDeclareN(loop_size);

#define hypre_BoxLoop1For(i1)\
   hypre__BoxLoop1For(hypre__i, hypre__j, hypre__k, i1)
#define hypre__BoxLoop1For(i, j, k, i1)\
   for (hypre__block = 0; hypre__block < hypre__num_blocks; hypre__block++)\
   {\
   hypre_BoxLoopSet(i, j, k);\
   i1 = hypre__i1start + i*hypre__sx1 + j*hypre__sy1 + k*hypre__sz1;\
   for (k = 0; k < hypre__nz; k++)\
   {\
      for (j = 0; j < hypre__ny; j++)\
      {\
         for (i = 0; i < hypre__nx; i++)\
         {

#define hypre_BoxLoop1End(i1)\
            i1 += hypre__sx1;\
         }\
         i1 += hypre__sy1 - hypre__nx*hypre__sx1;\
      }\
      i1 += hypre__sz1 - hypre__ny*hypre__sy1;\
   }\
   }\
}
  
/*-----------------------------------*/

#define hypre_BoxLoop2Begin(ndim,loop_size,\
                            dbox1, start1, stride1, i1,\
                            dbox2, start2, stride2, i2)\
{\
   HYPRE_Int  hypre__i1start = hypre_BoxIndexRank(dbox1, start1);\
   HYPRE_Int  hypre__i2start = hypre_BoxIndexRank(dbox2, start2);\
   hypre_BoxLoopDeclareS(dbox1, stride1, hypre__sx1, hypre__sy1, hypre__sz1);\
   hypre_BoxLoopDeclareS(dbox2, stride2, hypre__sx2, hypre__sy2, hypre__sz2);\
   hypre_BoxLoopDeclareN(loop_size);

#define hypre_BoxLoop2For(i1, i2)\
   hypre__BoxLoop2For(hypre__i, hypre__j, hypre__k, i1, i2)
#define hypre__BoxLoop2For(i, j, k, i1, i2)\
   for (hypre__block = 0; hypre__block < hypre__num_blocks; hypre__block++)\
   {\
   hypre_BoxLoopSet(i, j, k);\
   i1 = hypre__i1start + i*hypre__sx1 + j*hypre__sy1 + k*hypre__sz1;\
   i2 = hypre__i2start + i*hypre__sx2 + j*hypre__sy2 + k*hypre__sz2;\
   for (k = 0; k < hypre__nz; k++)\
   {\
      for (j = 0; j < hypre__ny; j++)\
      {\
         for (i = 0; i < hypre__nx; i++)\
         {

#define hypre_BoxLoop2End(i1, i2)\
            i1 += hypre__sx1;\
            i2 += hypre__sx2;\
         }\
         i1 += hypre__sy1 - hypre__nx*hypre__sx1;\
         i2 += hypre__sy2 - hypre__nx*hypre__sx2;\
      }\
      i1 += hypre__sz1 - hypre__ny*hypre__sy1;\
      i2 += hypre__sz2 - hypre__ny*hypre__sy2;\
   }\
   }\
}

/*-----------------------------------*/

#define hypre_BoxLoop3Begin(ndim, loop_size,\
                            dbox1, start1, stride1, i1,\
                            dbox2, start2, stride2, i2,\
                            dbox3, start3, stride3, i3)\
{\
   HYPRE_Int  hypre__i1start = hypre_BoxIndexRank(dbox1, start1);\
   HYPRE_Int  hypre__i2start = hypre_BoxIndexRank(dbox2, start2);\
   HYPRE_Int  hypre__i3start = hypre_BoxIndexRank(dbox3, start3);\
   hypre_BoxLoopDeclareS(dbox1, stride1, hypre__sx1, hypre__sy1, hypre__sz1);\
   hypre_BoxLoopDeclareS(dbox2, stride2, hypre__sx2, hypre__sy2, hypre__sz2);\
   hypre_BoxLoopDeclareS(dbox3, stride3, hypre__sx3, hypre__sy3, hypre__sz3);\
   hypre_BoxLoopDeclareN(loop_size);

#define hypre_BoxLoop3For(i1, i2, i3)\
   hypre__BoxLoop3For(hypre__i, hypre__j, hypre__k, i1, i2, i3)
#define hypre__BoxLoop3For(i, j, k, i1, i2, i3)\
   for (hypre__block = 0; hypre__block < hypre__num_blocks; hypre__block++)\
   {\
   hypre_BoxLoopSet(i, j, k);\
   i1 = hypre__i1start + i*hypre__sx1 + j*hypre__sy1 + k*hypre__sz1;\
   i2 = hypre__i2start + i*hypre__sx2 + j*hypre__sy2 + k*hypre__sz2;\
   i3 = hypre__i3start + i*hypre__sx3 + j*hypre__sy3 + k*hypre__sz3;\
   for (k = 0; k < hypre__nz; k++)\
   {\
      for (j = 0; j < hypre__ny; j++)\
      {\
         for (i = 0; i < hypre__nx; i++)\
         {

#define hypre_BoxLoop3End(i1, i2, i3)\
            i1 += hypre__sx1;\
            i2 += hypre__sx2;\
            i3 += hypre__sx3;\
         }\
         i1 += hypre__sy1 - hypre__nx*hypre__sx1;\
         i2 += hypre__sy2 - hypre__nx*hypre__sx2;\
         i3 += hypre__sy3 - hypre__nx*hypre__sx3;\
      }\
      i1 += hypre__sz1 - hypre__ny*hypre__sy1;\
      i2 += hypre__sz2 - hypre__ny*hypre__sy2;\
      i3 += hypre__sz3 - hypre__ny*hypre__sy3;\
   }\
   }\
}

/*-----------------------------------*/

#define hypre_BoxLoop4Begin(ndim, loop_size,\
                            dbox1, start1, stride1, i1,\
                            dbox2, start2, stride2, i2,\
                            dbox3, start3, stride3, i3,\
                            dbox4, start4, stride4, i4)\
{\
   HYPRE_Int  hypre__i1start = hypre_BoxIndexRank(dbox1, start1);\
   HYPRE_Int  hypre__i2start = hypre_BoxIndexRank(dbox2, start2);\
   HYPRE_Int  hypre__i3start = hypre_BoxIndexRank(dbox3, start3);\
   HYPRE_Int  hypre__i4start = hypre_BoxIndexRank(dbox4, start4);\
   hypre_BoxLoopDeclareS(dbox1, stride1, hypre__sx1, hypre__sy1, hypre__sz1);\
   hypre_BoxLoopDeclareS(dbox2, stride2, hypre__sx2, hypre__sy2, hypre__sz2);\
   hypre_BoxLoopDeclareS(dbox3, stride3, hypre__sx3, hypre__sy3, hypre__sz3);\
   hypre_BoxLoopDeclareS(dbox4, stride4, hypre__sx4, hypre__sy4, hypre__sz4);\
   hypre_BoxLoopDeclareN(loop_size);

#define hypre_BoxLoop4For(i1, i2, i3, i4)\
   hypre__BoxLoop4For(hypre__i, hypre__j, hypre__k, i1, i2, i3, i4)
#define hypre__BoxLoop4For(i, j, k, i1, i2, i3, i4)\
   for (hypre__block = 0; hypre__block < hypre__num_blocks; hypre__block++)\
   {\
   hypre_BoxLoopSet(i, j, k);\
   i1 = hypre__i1start + i*hypre__sx1 + j*hypre__sy1 + k*hypre__sz1;\
   i2 = hypre__i2start + i*hypre__sx2 + j*hypre__sy2 + k*hypre__sz2;\
   i3 = hypre__i3start + i*hypre__sx3 + j*hypre__sy3 + k*hypre__sz3;\
   i4 = hypre__i4start + i*hypre__sx4 + j*hypre__sy4 + k*hypre__sz4;\
   for (k = 0; k < hypre__nz; k++)\
   {\
      for (j = 0; j < hypre__ny; j++)\
      {\
         for (i = 0; i < hypre__nx; i++)\
         {

#define hypre_BoxLoop4End(i1, i2, i3, i4)\
            i1 += hypre__sx1;\
            i2 += hypre__sx2;\
            i3 += hypre__sx3;\
            i4 += hypre__sx4;\
         }\
         i1 += hypre__sy1 - hypre__nx*hypre__sx1;\
         i2 += hypre__sy2 - hypre__nx*hypre__sx2;\
         i3 += hypre__sy3 - hypre__nx*hypre__sx3;\
         i4 += hypre__sy4 - hypre__nx*hypre__sx4;\
      }\
      i1 += hypre__sz1 - hypre__ny*hypre__sy1;\
      i2 += hypre__sz2 - hypre__ny*hypre__sy2;\
      i3 += hypre__sz3 - hypre__ny*hypre__sy3;\
      i4 += hypre__sz4 - hypre__ny*hypre__sy4;\
   }\
   }\
}

/*-----------------------------------*/

#endif  /* ifndef HYPRE_USE_PTHREADS */

#endif

/******************************************************************************
 *
 * NEW BoxLoop STUFF
 *
 *****************************************************************************/

#ifndef hypre_ZBOX_HEADER
#define hypre_ZBOX_HEADER

/*--------------------------------------------------------------------------
 * BoxLoop macros:
 *--------------------------------------------------------------------------*/

#define zypre_BoxLoopDeclare() \
HYPRE_Int  hypre__i[HYPRE_MAXDIM+1];\
HYPRE_Int  hypre__n[HYPRE_MAXDIM+1];\
HYPRE_Int  hypre__m[HYPRE_MAXDIM];\
HYPRE_Int  hypre__d, hypre__ndim;\
HYPRE_Int  hypre__dir, hypre__max;\
HYPRE_Int  hypre__div, hypre__mod;\
HYPRE_Int  hypre__block, hypre__num_blocks;\
HYPRE_Int  hypre__I, hypre__J, hypre__IN, hypre__JN

#define zypre_BoxLoopDeclareK(k) \
HYPRE_Int  hypre__sk[k][HYPRE_MAXDIM];\
HYPRE_Int  hypre__ikstart[k], hypre__ikinc[k][HYPRE_MAXDIM+1], hypre__i0inc[k]

#define zypre_BoxLoopInit(ndim, loop_size) \
hypre__ndim = ndim;\
for (hypre__d = 0; hypre__d < hypre__ndim; hypre__d++)\
{\
   hypre__n[hypre__d] = hypre__m[hypre__d] = loop_size[hypre__d];\
}\
hypre__dir = 0;\
hypre__max = hypre__n[0];\
for (hypre__d = 1; hypre__d < hypre__ndim; hypre__d++)\
{\
   if (hypre__n[hypre__d] > hypre__max)\
   {\
      hypre__dir = hypre__d;\
      hypre__max = hypre__n[hypre__d];\
   }\
}\
hypre__num_blocks = hypre_NumThreads();\
if (hypre__max < hypre__num_blocks)\
{\
   hypre__num_blocks = hypre__max;\
}\
if (hypre__num_blocks > 0)\
{\
   hypre__div = hypre__max / hypre__num_blocks;\
   hypre__mod = hypre__max % hypre__num_blocks;\
}

#define zypre_BoxLoopInitK(k, dboxk, startk, stridek, ik) \
hypre__sk[k][0] = 1;\
for (hypre__d = 1; hypre__d < hypre__ndim; hypre__d++)\
{\
   hypre__sk[k][hypre__d] =\
      hypre__sk[k][hypre__d-1]*hypre_BoxSizeD(dboxk, hypre__d-1);\
}\
for (hypre__d = 0; hypre__d < hypre__ndim; hypre__d++)\
{\
   hypre__sk[k][hypre__d] = hypre__sk[k][hypre__d]*stridek[hypre__d];\
}\
hypre__ikstart[k] = hypre_BoxIndexRank(dboxk, startk)

#define zypre_BoxLoopSet() \
for (hypre__d = 0; hypre__d < hypre__ndim; hypre__d++)\
{\
   hypre__n[hypre__d] = hypre__m[hypre__d];\
}\
if (hypre__num_blocks > 1)\
{\
   hypre__i[hypre__dir] =\
      hypre__block * hypre__div + hypre_min(hypre__mod, hypre__block);\
   hypre__n[hypre__dir] =\
      hypre__div + ((hypre__mod > hypre__block) ? 1 : 0);\
}

#define zypre_BoxLoopSetK(k, ik) \
ik = hypre__ikstart[k];\
if (hypre__num_blocks > 1)\
{\
   ik += hypre__i[hypre__dir]*hypre__sk[k][hypre__dir];\
}\
hypre__ikinc[k][0] = hypre__sk[k][0];\
hypre__i0inc[k] = hypre__ikinc[k][0];\
if (hypre__ndim > 1)\
{\
   hypre__ikinc[k][1] = hypre__sk[k][1] - hypre__n[0]*hypre__sk[k][0];\
}\
for (hypre__d = 2; hypre__d < hypre__ndim; hypre__d++)\
{\
   hypre__ikinc[k][hypre__d] = hypre__ikinc[k][hypre__d-1] +\
      hypre__sk[k][hypre__d] - hypre__n[hypre__d-1]*hypre__sk[k][hypre__d-1];\
}\
hypre__ikinc[k][hypre__ndim] = 0

#define zypre_BoxLoopSetLoop() \
hypre__IN = hypre__n[0];\
hypre__JN = 1;\
for (hypre__d = 1; hypre__d < hypre__ndim; hypre__d++)\
{\
   hypre__JN *= hypre__n[hypre__d];\
   hypre__i[hypre__d] = 0;\
   hypre__n[hypre__d] -= 2; /* this produces a simpler comparison below */\
}\
hypre__i[hypre__ndim] = 0;\
hypre__n[hypre__ndim] = 0

#define zypre_BoxLoopIncLoop() \
for (hypre__d = 1; hypre__i[hypre__d] > hypre__n[hypre__d]; hypre__d++)\
{\
   hypre__i[hypre__d] = 0;\
}\
hypre__i[hypre__d]++

#define zypre_BoxLoopGetIndex(index) \
index[0] = hypre__I;\
for (hypre__d = 1; hypre__d < hypre__ndim; hypre__d++)\
{\
   index[hypre__d] = hypre__i[hypre__d];\
}

/*-----------------------------------*/

#define zypre_BoxLoop0Begin(ndim, loop_size)\
{\
   zypre_BoxLoopDeclare();\
   zypre_BoxLoopInit(ndim, loop_size);

#define zypre_BoxLoop0For()\
   for (hypre__block = 0; hypre__block < hypre__num_blocks; hypre__block++)\
   {\
      zypre_BoxLoopSet();\
      zypre_BoxLoopSetLoop();\
      for (hypre__J = 0; hypre__J < hypre__JN; hypre__J++)\
      {\
         for (hypre__I = 0; hypre__I < hypre__IN; hypre__I++)\
         {

#define zypre_BoxLoop0End()\
         }\
         zypre_BoxLoopIncLoop();\
      }\
   }\
}

/*-----------------------------------*/

#define zypre_BoxLoop1Begin(ndim, loop_size,\
                            dbox1, start1, stride1, i1)\
{\
   zypre_BoxLoopDeclare();\
   zypre_BoxLoopDeclareK(1);\
   zypre_BoxLoopInit(ndim, loop_size);\
   zypre_BoxLoopInitK(0, dbox1, start1, stride1, i1);

#define zypre_BoxLoop1For(i1)\
   for (hypre__block = 0; hypre__block < hypre__num_blocks; hypre__block++)\
   {\
      zypre_BoxLoopSet();\
      zypre_BoxLoopSetK(0, i1);\
      zypre_BoxLoopSetLoop();\
      for (hypre__J = 0; hypre__J < hypre__JN; hypre__J++)\
      {\
         for (hypre__I = 0; hypre__I < hypre__IN; hypre__I++)\
         {

#define zypre_BoxLoop1End(i1)\
            i1 += hypre__i0inc[0];\
         }\
         zypre_BoxLoopIncLoop();\
         i1 += hypre__ikinc[0][hypre__d];\
      }\
   }\
}

/*-----------------------------------*/

#define zypre_BoxLoop2Begin(ndim, loop_size,\
                            dbox1, start1, stride1, i1,\
                            dbox2, start2, stride2, i2)\
{\
   zypre_BoxLoopDeclare();\
   zypre_BoxLoopDeclareK(2);\
   zypre_BoxLoopInit(ndim, loop_size);\
   zypre_BoxLoopInitK(0, dbox1, start1, stride1, i1);\
   zypre_BoxLoopInitK(1, dbox2, start2, stride2, i2);

#define zypre_BoxLoop2For(i1, i2)\
   for (hypre__block = 0; hypre__block < hypre__num_blocks; hypre__block++)\
   {\
      zypre_BoxLoopSet();\
      zypre_BoxLoopSetK(0, i1);\
      zypre_BoxLoopSetK(1, i2);\
      zypre_BoxLoopSetLoop();\
      for (hypre__J = 0; hypre__J < hypre__JN; hypre__J++)\
      {\
         for (hypre__I = 0; hypre__I < hypre__IN; hypre__I++)\
         {

#define zypre_BoxLoop2End(i1, i2)\
            i1 += hypre__i0inc[0];\
            i2 += hypre__i0inc[1];\
         }\
         zypre_BoxLoopIncLoop();\
         i1 += hypre__ikinc[0][hypre__d];\
         i2 += hypre__ikinc[1][hypre__d];\
      }\
   }\
}

/*-----------------------------------*/

#define zypre_BoxLoop3Begin(ndim, loop_size,\
                            dbox1, start1, stride1, i1,\
                            dbox2, start2, stride2, i2,\
                            dbox3, start3, stride3, i3)\
{\
   zypre_BoxLoopDeclare();\
   zypre_BoxLoopDeclareK(3);\
   zypre_BoxLoopInit(ndim, loop_size);\
   zypre_BoxLoopInitK(0, dbox1, start1, stride1, i1);\
   zypre_BoxLoopInitK(1, dbox2, start2, stride2, i2);\
   zypre_BoxLoopInitK(2, dbox3, start3, stride3, i3);

#define zypre_BoxLoop3For(i1, i2, i3)\
   for (hypre__block = 0; hypre__block < hypre__num_blocks; hypre__block++)\
   {\
      zypre_BoxLoopSet();\
      zypre_BoxLoopSetK(0, i1);\
      zypre_BoxLoopSetK(1, i2);\
      zypre_BoxLoopSetK(2, i3);\
      zypre_BoxLoopSetLoop();\
      for (hypre__J = 0; hypre__J < hypre__JN; hypre__J++)\
      {\
         for (hypre__I = 0; hypre__I < hypre__IN; hypre__I++)\
         {

#define zypre_BoxLoop3End(i1, i2, i3)\
            i1 += hypre__i0inc[0];\
            i2 += hypre__i0inc[1];\
            i3 += hypre__i0inc[2];\
/*            i1 ++;*/\
/*            i2 ++;*/\
/*            i3 ++;*/\
         }\
         zypre_BoxLoopIncLoop();\
         i1 += hypre__ikinc[0][hypre__d];\
         i2 += hypre__ikinc[1][hypre__d];\
         i3 += hypre__ikinc[2][hypre__d];\
      }\
   }\
}

/*-----------------------------------*/

#define zypre_BoxLoop4Begin(ndim, loop_size,\
                            dbox1, start1, stride1, i1,\
                            dbox2, start2, stride2, i2,\
                            dbox3, start3, stride3, i3,\
                            dbox4, start4, stride4, i4)\
{\
   zypre_BoxLoopDeclare();\
   zypre_BoxLoopDeclareK(4);\
   zypre_BoxLoopInit(ndim, loop_size);\
   zypre_BoxLoopInitK(0, dbox1, start1, stride1, i1);\
   zypre_BoxLoopInitK(1, dbox2, start2, stride2, i2);\
   zypre_BoxLoopInitK(2, dbox3, start3, stride3, i3);\
   zypre_BoxLoopInitK(3, dbox4, start4, stride4, i4);

#define zypre_BoxLoop4For(i1, i2, i3, i4)\
   for (hypre__block = 0; hypre__block < hypre__num_blocks; hypre__block++)\
   {\
      zypre_BoxLoopSet();\
      zypre_BoxLoopSetK(0, i1);\
      zypre_BoxLoopSetK(1, i2);\
      zypre_BoxLoopSetK(2, i3);\
      zypre_BoxLoopSetK(3, i4);\
      zypre_BoxLoopSetLoop();\
      for (hypre__J = 0; hypre__J < hypre__JN; hypre__J++)\
      {\
         for (hypre__I = 0; hypre__I < hypre__IN; hypre__I++)\
         {

#define zypre_BoxLoop4End(i1, i2, i3, i4)\
            i1 += hypre__i0inc[0];\
            i2 += hypre__i0inc[1];\
            i3 += hypre__i0inc[2];\
            i4 += hypre__i0inc[3];\
         }\
         zypre_BoxLoopIncLoop();\
         i1 += hypre__ikinc[0][hypre__d];\
         i2 += hypre__ikinc[1][hypre__d];\
         i3 += hypre__ikinc[2][hypre__d];\
         i4 += hypre__ikinc[3][hypre__d];\
      }\
   }\
}

/*-----------------------------------*/

#endif




#if 0

/*- NOTES ------------------------------*/

#define hypre_BoxLoop2Begin(loop_size,
                            dbox1, start1, stride1, i1,
                            dbox2, start2, stride2, i2)
{
   /* init hypre__i1start */
   HYPRE_Int  hypre__i1start = hypre_BoxIndexRank(dbox1, start1);
   HYPRE_Int  hypre__i2start = hypre_BoxIndexRank(dbox2, start2);
   /* declare and set hypre__s1 */
   hypre_BoxLoopDeclareS(dbox1, stride1, hypre__sx1, hypre__sy1, hypre__sz1);
   hypre_BoxLoopDeclareS(dbox2, stride2, hypre__sx2, hypre__sy2, hypre__sz2);
   /* declare and set hypre__n, hypre__m, hypre__dir, hypre__max,
    *                 hypre__div, hypre__mod, hypre__block, hypre__num_blocks */
   hypre_BoxLoopDeclareN(loop_size);

#define hypre_BoxLoop2For(i, j, k, i1, i2)
   for (hypre__block = 0; hypre__block < hypre__num_blocks; hypre__block++)
   {
   /* set i and hypre__n */
   hypre_BoxLoopSet(i, j, k);
   /* set i1 */
   i1 = hypre__i1start + i*hypre__sx1 + j*hypre__sy1 + k*hypre__sz1;
   i2 = hypre__i2start + i*hypre__sx2 + j*hypre__sy2 + k*hypre__sz2;
   for (k = 0; k < hypre__nz; k++)
   {
      for (j = 0; j < hypre__ny; j++)
      {
         for (i = 0; i < hypre__nx; i++)
         {

#define hypre_BoxLoop2End(i1, i2)
            i1 += hypre__sx1;
            i2 += hypre__sx2;
         }
         i1 += hypre__sy1 - hypre__nx*hypre__sx1;
         i2 += hypre__sy2 - hypre__nx*hypre__sx2;
      }
      i1 += hypre__sz1 - hypre__ny*hypre__sy1;
      i2 += hypre__sz2 - hypre__ny*hypre__sy2;
   }
   }
}

/*-----------------------------------*/

/* Idea 2 */

N = 1;
for (d = 0; d < ndim; d++)
{
   N *= n[d];
   i[d] = 0;
   n[d] -= 2; /* this produces a simpler comparison below */
}
i[ndim] = 0;
n[ndim] = 0;
for (I = 0; I < N; I++)
{
   /* loop body */

   for (d = 0; i[d] > n[d]; d++)
   {
      i[d] = 0;
   }
   i[d]++;
   i1 += s1[d]; /* NOTE: These are different from hypre__sx1, etc. above */
   i2 += s2[d]; /* The lengths of i, n, and s must be (ndim+1) */
}

/*-----------------------------------*/

/* Idea 3 */

N = 1;
for (d = 1; d < ndim; d++)
{
   N *= n[d];
   i[d] = 0;
   n[d] -= 2; /* this produces a simpler comparison below */
}
i[ndim] = 0;
n[ndim] = 0;
for (J = 0; J < N; J++)
{
   for (I = 0; I < n[0]; I++)
   {
      /* loop body */

      i1 += s1[0];
      i2 += s2[0];
   }
   for (d = 1; i[d] > n[d]; d++)
   {
      i[d] = 0;
   }
   i[d]++;
   i1 += s1[d]; /* NOTE: These are different from hypre__sx1, etc. above */
   i2 += s2[d]; /* The lengths of i, n, and s must be (ndim+1) */
}

#endif
