/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision: 2.5 $
 ***********************************************************************EHEADER*/




 
#include "headers.h"
 
/*--------------------------------------------------------------------------
 * hypre_GenerateLaplacian27pt
 *--------------------------------------------------------------------------*/

hypre_CSRMatrix *
hypre_GenerateLaplacian27pt(HYPRE_Int      nx,
                            HYPRE_Int      ny,
                            HYPRE_Int      nz,
                            HYPRE_Int      P,
                            HYPRE_Int      Q,
                            HYPRE_Int      R,
                            double  *value )
{
   hypre_CSRMatrix *A;

   HYPRE_Int    *A_i;
   HYPRE_Int    *A_j;
   double *A_data;

   HYPRE_Int *global_part;
   HYPRE_Int ix, iy, iz;
   HYPRE_Int cnt;
   HYPRE_Int row_index;
   HYPRE_Int p, q, r;

   HYPRE_Int nx_local, ny_local;
   HYPRE_Int nx_size, ny_size, nz_size;
   HYPRE_Int nxy;
   HYPRE_Int grid_size;

   HYPRE_Int *nx_part;
   HYPRE_Int *ny_part;
   HYPRE_Int *nz_part;

   grid_size = nx*ny*nz;

   hypre_GeneratePartitioning(nx,P,&nx_part);
   hypre_GeneratePartitioning(ny,Q,&ny_part);
   hypre_GeneratePartitioning(nz,R,&nz_part);

   global_part = hypre_CTAlloc(HYPRE_Int,P*Q*R+1);

   global_part[0] = 0;
   cnt = 1;
   for (iz = 0; iz < R; iz++)
   {
      nz_size = nz_part[iz+1]-nz_part[iz];
      for (iy = 0; iy < Q; iy++)
      {
         ny_size = (ny_part[iy+1]-ny_part[iy])*nz_size;
         for (ix = 0; ix < P; ix++)
         {
            nx_size = nx_part[ix+1] - nx_part[ix];
            global_part[cnt] = global_part[cnt-1];
            global_part[cnt++] += nx_size*ny_size;
         }
      }
   }
   
   A_i = hypre_CTAlloc(HYPRE_Int, grid_size+1);

   cnt = 0;
   A_i[0] = 0;
   for (r = 0; r < R; r++)
   {
   for (q = 0; q < Q; q++)
   {
   for (p = 0; p < P; p++)
   {
   for (iz = nz_part[r];  iz < nz_part[r+1]; iz++)
   {
      for (iy = ny_part[q];  iy < ny_part[q+1]; iy++)
      {
         for (ix = nx_part[p]; ix < nx_part[p+1]; ix++)
         {
            cnt++;
            A_i[cnt] = A_i[cnt-1];
            A_i[cnt]++;
            if (iz > nz_part[r]) 
            {
               A_i[cnt]++;
               if (iy > ny_part[q]) 
               {
                  A_i[cnt]++;
      	          if (ix > nx_part[p])
      	          {
      	             A_i[cnt]++;
      	          }
      	          else
      	          {
      	             if (ix) 
      		        A_i[cnt]++;
      	          }
      	          if (ix < nx_part[p+1]-1)
      	          {
      	             A_i[cnt]++;
      	          }
      	          else
      	          {
      	             if (ix+1 < nx) 
      		        A_i[cnt]++;
      	          }
               }
               else
               {
                  if (iy) 
                  {
                     A_i[cnt]++;
      	             if (ix > nx_part[p])
      	             {
      	                A_i[cnt]++;
      	             }
      	             else if (ix)
      	             {
      	                A_i[cnt]++;
      	             }
      	             if (ix < nx_part[p+1]-1)
      	             {
      	                A_i[cnt]++;
      	             }
      	             else if (ix < nx-1)
      	             {
      	                A_i[cnt]++;
      	             }
                  }
               }
               if (ix > nx_part[p]) 
                  A_i[cnt]++;
               else
               {
                  if (ix) 
                  {
                     A_i[cnt]++; 
                  }
               }
               if (ix+1 < nx_part[p+1]) 
                  A_i[cnt]++;
               else
               {
                  if (ix+1 < nx) 
                  {
                     A_i[cnt]++; 
                  }
               }
               if (iy+1 < ny_part[q+1]) 
               {
                  A_i[cnt]++;
      	          if (ix > nx_part[p])
      	          {
      	             A_i[cnt]++;
      	          }
      	          else
      	          {
      	             if (ix) 
      		        A_i[cnt]++;
      	          }
      	          if (ix < nx_part[p+1]-1)
      	          {
      	             A_i[cnt]++;
      	          }
      	          else
      	          {
      	             if (ix+1 < nx) 
      		        A_i[cnt]++;
      	          }
               }
               else
               {
                  if (iy+1 < ny) 
                  {
                     A_i[cnt]++;
      	             if (ix > nx_part[p])
      	             {
      	                A_i[cnt]++;
      	             }
      	             else if (ix)
      	             {
      	                A_i[cnt]++;
      	             }
      	             if (ix < nx_part[p+1]-1)
      	             {
      	                A_i[cnt]++;
      	             }
      	             else if (ix < nx-1)
      	             {
      	                A_i[cnt]++;
      	             }
                  }
               }
            }
            else
            {
               if (iz)
	       {
		  A_i[cnt]++;
                  if (iy > ny_part[q]) 
                  {
                     A_i[cnt]++;
      	             if (ix > nx_part[p])
      	             {
      	                A_i[cnt]++;
      	             }
      	             else
      	             {
      	                if (ix) 
      	   	        A_i[cnt]++;
      	             }
      	             if (ix < nx_part[p+1]-1)
      	             {
      	                A_i[cnt]++;
      	             }
      	             else
      	             {
      	                if (ix+1 < nx) 
      	   	        A_i[cnt]++;
      	             }
                  }
                  else
                  {
                     if (iy) 
                     {
                        A_i[cnt]++;
      	                if (ix > nx_part[p])
      	                {
      	                   A_i[cnt]++;
      	                }
      	                else if (ix)
      	                {
      	                   A_i[cnt]++;
      	                }
      	                if (ix < nx_part[p+1]-1)
      	                {
      	                   A_i[cnt]++;
      	                }
      	                else if (ix < nx-1)
      	                {
      	                   A_i[cnt]++;
      	                }
                     }
                  }
                  if (ix > nx_part[p]) 
                     A_i[cnt]++;
                  else
                  {
                     if (ix) 
                     {
                        A_i[cnt]++; 
                     }
                  }
                  if (ix+1 < nx_part[p+1]) 
                     A_i[cnt]++;
                  else
                  {
                     if (ix+1 < nx) 
                     {
                        A_i[cnt]++; 
                     }
                  }
                  if (iy+1 < ny_part[q+1]) 
                  {
                     A_i[cnt]++;
      	             if (ix > nx_part[p])
      	             {
      	                A_i[cnt]++;
      	             }
      	             else
      	             {
      	                if (ix) 
      	   	           A_i[cnt]++;
      	             }
      	             if (ix < nx_part[p+1]-1)
      	             {
      	                A_i[cnt]++;
      	             }
      	             else
      	             {
      	                if (ix+1 < nx) 
      	   	           A_i[cnt]++;
      	             }
                  }
                  else
                  {
                     if (iy+1 < ny) 
                     {
                        A_i[cnt]++;
      	                if (ix > nx_part[p])
      	                {
      	                   A_i[cnt]++;
      	                }
      	                else if (ix)
      	                {
      	                   A_i[cnt]++;
      	                }
      	                if (ix < nx_part[p+1]-1)
      	                {
      	                   A_i[cnt]++;
      	                }
      	                else if (ix < nx-1)
      	                {
      	                   A_i[cnt]++;
      	                }
                     }
                  }
               }
            }
            if (iy > ny_part[q]) 
            {
               A_i[cnt]++;
   	       if (ix > nx_part[p])
   	       {
   	          A_i[cnt]++;
   	       }
   	       else
   	       {
   	          if (ix) 
   		     A_i[cnt]++;
   	       }
   	       if (ix < nx_part[p+1]-1)
   	       {
   	          A_i[cnt]++;
   	       }
   	       else
   	       {
   	          if (ix+1 < nx) 
   		     A_i[cnt]++;
   	       }
            }
            else
            {
               if (iy) 
               {
                  A_i[cnt]++;
   	          if (ix > nx_part[p])
   	          {
   	             A_i[cnt]++;
   	          }
   	          else if (ix)
   	          {
   	             A_i[cnt]++;
   	          }
   	          if (ix < nx_part[p+1]-1)
   	          {
   	             A_i[cnt]++;
   	          }
   	          else if (ix < nx-1)
   	          {
   	             A_i[cnt]++;
   	          }
               }
            }
            if (ix > nx_part[p]) 
               A_i[cnt]++;
            else
            {
               if (ix) 
               {
                  A_i[cnt]++; 
               }
            }
            if (ix+1 < nx_part[p+1]) 
               A_i[cnt]++;
            else
            {
               if (ix+1 < nx) 
               {
                  A_i[cnt]++; 
               }
            }
            if (iy+1 < ny_part[q+1]) 
            {
               A_i[cnt]++;
   	       if (ix > nx_part[p])
   	       {
   	          A_i[cnt]++;
   	       }
   	       else
   	       {
   	          if (ix) 
   		     A_i[cnt]++;
   	       }
   	       if (ix < nx_part[p+1]-1)
   	       {
   	          A_i[cnt]++;
   	       }
   	       else
   	       {
   	          if (ix+1 < nx) 
   		     A_i[cnt]++;
   	       }
            }
            else
            {
               if (iy+1 < ny) 
               {
                  A_i[cnt]++;
   	          if (ix > nx_part[p])
   	          {
   	             A_i[cnt]++;
   	          }
   	          else if (ix)
   	          {
   	             A_i[cnt]++;
   	          }
   	          if (ix < nx_part[p+1]-1)
   	          {
   	             A_i[cnt]++;
   	          }
   	          else if (ix < nx-1)
   	          {
   	             A_i[cnt]++;
   	          }
               }
            }
            if (iz+1 < nz_part[r+1]) 
            {
               A_i[cnt]++;
               if (iy > ny_part[q]) 
               {
                  A_i[cnt]++;
      	          if (ix > nx_part[p])
      	          {
      	             A_i[cnt]++;
      	          }
      	          else
      	          {
      	             if (ix) 
      		        A_i[cnt]++;
      	          }
      	          if (ix < nx_part[p+1]-1)
      	          {
      	             A_i[cnt]++;
      	          }
      	          else
      	          {
      	             if (ix+1 < nx) 
      		        A_i[cnt]++;
      	          }
               }
               else
               {
                  if (iy) 
                  {
                     A_i[cnt]++;
      	             if (ix > nx_part[p])
      	             {
      	                A_i[cnt]++;
      	             }
      	             else if (ix)
      	             {
      	                A_i[cnt]++;
      	             }
      	             if (ix < nx_part[p+1]-1)
      	             {
      	                A_i[cnt]++;
      	             }
      	             else if (ix < nx-1)
      	             {
      	                A_i[cnt]++;
      	             }
                  }
               }
               if (ix > nx_part[p]) 
                  A_i[cnt]++;
               else
               {
                  if (ix) 
                  {
                     A_i[cnt]++; 
                  }
               }
               if (ix+1 < nx_part[p+1]) 
                  A_i[cnt]++;
               else
               {
                  if (ix+1 < nx) 
                  {
                     A_i[cnt]++; 
                  }
               }
               if (iy+1 < ny_part[q+1]) 
               {
                  A_i[cnt]++;
      	          if (ix > nx_part[p])
      	          {
      	             A_i[cnt]++;
      	          }
      	          else
      	          {
      	             if (ix) 
      		        A_i[cnt]++;
      	          }
      	          if (ix < nx_part[p+1]-1)
      	          {
      	             A_i[cnt]++;
      	          }
      	          else
      	          {
      	             if (ix+1 < nx) 
      		        A_i[cnt]++;
      	          }
               }
               else
               {
                  if (iy+1 < ny) 
                  {
                     A_i[cnt]++;
      	             if (ix > nx_part[p])
      	             {
      	                A_i[cnt]++;
      	             }
      	             else if (ix)
      	             {
      	                A_i[cnt]++;
      	             }
      	             if (ix < nx_part[p+1]-1)
      	             {
      	                A_i[cnt]++;
      	             }
      	             else if (ix < nx-1)
      	             {
      	                A_i[cnt]++;
      	             }
                  }
               }
            }
            else
            {
               if (iz+1 < nz)
	       {
		  A_i[cnt]++;
                  if (iy > ny_part[q]) 
                  {
                     A_i[cnt]++;
      	             if (ix > nx_part[p])
      	             {
      	                A_i[cnt]++;
      	             }
      	             else
      	             {
      	                if (ix) 
      	   	        A_i[cnt]++;
      	             }
      	             if (ix < nx_part[p+1]-1)
      	             {
      	                A_i[cnt]++;
      	             }
      	             else
      	             {
      	                if (ix+1 < nx) 
      	   	        A_i[cnt]++;
      	             }
                  }
                  else
                  {
                     if (iy) 
                     {
                        A_i[cnt]++;
      	                if (ix > nx_part[p])
      	                {
      	                   A_i[cnt]++;
      	                }
      	                else if (ix)
      	                {
      	                   A_i[cnt]++;
      	                }
      	                if (ix < nx_part[p+1]-1)
      	                {
      	                   A_i[cnt]++;
      	                }
      	                else if (ix < nx-1)
      	                {
      	                   A_i[cnt]++;
      	                }
                     }
                  }
                  if (ix > nx_part[p]) 
                     A_i[cnt]++;
                  else
                  {
                     if (ix) 
                     {
                        A_i[cnt]++; 
                     }
                  }
                  if (ix+1 < nx_part[p+1]) 
                     A_i[cnt]++;
                  else
                  {
                     if (ix+1 < nx) 
                     {
                        A_i[cnt]++; 
                     }
                  }
                  if (iy+1 < ny_part[q+1]) 
                  {
                     A_i[cnt]++;
      	             if (ix > nx_part[p])
      	             {
      	                A_i[cnt]++;
      	             }
      	             else
      	             {
      	                if (ix) 
      	   	           A_i[cnt]++;
      	             }
      	             if (ix < nx_part[p+1]-1)
      	             {
      	                A_i[cnt]++;
      	             }
      	             else
      	             {
      	                if (ix+1 < nx) 
      	   	           A_i[cnt]++;
      	             }
                  }
                  else
                  {
                     if (iy+1 < ny) 
                     {
                        A_i[cnt]++;
      	                if (ix > nx_part[p])
      	                {
      	                   A_i[cnt]++;
      	                }
      	                else if (ix)
      	                {
      	                   A_i[cnt]++;
      	                }
      	                if (ix < nx_part[p+1]-1)
      	                {
      	                   A_i[cnt]++;
      	                }
      	                else if (ix < nx-1)
      	                {
      	                   A_i[cnt]++;
      	                }
                     }
                  }
               }
            }
         }
      }
   }
   }
   }
   }

   A_j = hypre_CTAlloc(HYPRE_Int, A_i[grid_size]);
   A_data = hypre_CTAlloc(double, A_i[grid_size]);

   row_index = 0;
   cnt = 0;
   for (r=0; r < R; r++)
   {
   for (q=0; q < Q; q++)
   {
   for (p=0; p < P; p++)
   {
   for (iz = nz_part[r];  iz < nz_part[r+1]; iz++)
   {
      for (iy = ny_part[q];  iy < ny_part[q+1]; iy++)
      {
   	 ny_local = ny_part[q+1] - ny_part[q];
         for (ix = nx_part[p]; ix < nx_part[p+1]; ix++)
         {
   	    nx_local = nx_part[p+1] - nx_part[p];
   	    nxy = nx_local*ny_local;
            A_j[cnt] = row_index;
            A_data[cnt++] = value[0];
            if (iz > nz_part[r]) 
            {
               if (iy > ny_part[q]) 
               {
      	          if (ix > nx_part[p])
      	          {
      	             A_j[cnt] = row_index-nxy-nx_local-1;
      	             A_data[cnt++] = value[1];
      	          }
      	          else
      	          {
      	             if (ix) 
      	             {
      		        A_j[cnt] = map3(ix-1,iy-1,iz-1,p-1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	          }
      	          A_j[cnt] = row_index-nxy-nx_local;
      	          A_data[cnt++] = value[1];
      	          if (ix < nx_part[p+1]-1)
      	          {
      	             A_j[cnt] = row_index-nxy-nx_local+1;
      	             A_data[cnt++] = value[1];
      	          }
      	          else
      	          {
      	             if (ix+1 < nx) 
      	             {
      		        A_j[cnt] = map3(ix+1,iy-1,iz-1,p+1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	          }
               }
               else
               {
                  if (iy) 
                  {
      	             if (ix > nx_part[p])
      	             {
      		        A_j[cnt] = map3(ix-1,iy-1,iz-1,p,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else if (ix)
      	             {
      		        A_j[cnt] = map3(ix-1,iy-1,iz-1,p-1,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      		     A_j[cnt] = map3(ix,iy-1,iz-1,p,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
      	             if (ix < nx_part[p+1]-1)
      	             {
      		        A_j[cnt] = map3(ix+1,iy-1,iz-1,p,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else if (ix < nx-1)
      	             {
      		        A_j[cnt] = map3(ix+1,iy-1,iz-1,p+1,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
                  }
               }
               if (ix > nx_part[p]) 
      	       {   
      	          A_j[cnt] = row_index-nxy-1;
      	          A_data[cnt++] = value[1];
      	       }   
               else
               {
                  if (ix) 
                  {
      		     A_j[cnt] = map3(ix-1,iy,iz-1,p-1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
                  }
               }
      	       A_j[cnt] = row_index-nxy;
      	       A_data[cnt++] = value[1];
               if (ix+1 < nx_part[p+1]) 
      	       {   
      	          A_j[cnt] = row_index-nxy+1;
      	          A_data[cnt++] = value[1];
      	       }   
               else
               {
                  if (ix+1 < nx) 
                  {
      		     A_j[cnt] = map3(ix+1,iy,iz-1,p+1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
                  }
               }
               if (iy+1 < ny_part[q+1]) 
               {
      	          if (ix > nx_part[p])
      	          {
      	             A_j[cnt] = row_index-nxy+nx_local-1;
      	             A_data[cnt++] = value[1];
      	          }
      	          else
      	          {
      	             if (ix) 
                     {
      		        A_j[cnt] = map3(ix-1,iy+1,iz-1,p-1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
                     }
      	          }
      	          A_j[cnt] = row_index-nxy+nx_local;
      	          A_data[cnt++] = value[1];
      	          if (ix < nx_part[p+1]-1)
      	          {
      	             A_j[cnt] = row_index-nxy+nx_local+1;
      	             A_data[cnt++] = value[1];
      	          }
      	          else
      	          {
      	             if (ix+1 < nx) 
                     {
      		        A_j[cnt] = map3(ix+1,iy+1,iz-1,p+1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
                     }
      	          }
               }
               else
               {
                  if (iy+1 < ny) 
                  {
      	             if (ix > nx_part[p])
      	             {
      		        A_j[cnt] = map3(ix-1,iy+1,iz-1,p,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else if (ix)
      	             {
      		        A_j[cnt] = map3(ix-1,iy+1,iz-1,p-1,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      		     A_j[cnt] = map3(ix,iy+1,iz-1,p,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
      	             if (ix < nx_part[p+1]-1)
      	             {
      		        A_j[cnt] = map3(ix+1,iy+1,iz-1,p,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else if (ix < nx-1)
      	             {
      		        A_j[cnt] = map3(ix+1,iy+1,iz-1,p+1,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
                  }
               }
            }
            else
            {
               if (iz)
	       {
                  if (iy > ny_part[q]) 
                  {
      	             if (ix > nx_part[p])
      	             {
      		        A_j[cnt] = map3(ix-1,iy-1,iz-1,p,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else
      	             {
      	                if (ix) 
      	                {
      		           A_j[cnt] = map3(ix-1,iy-1,iz-1,p-1,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	             }
      		     A_j[cnt] = map3(ix,iy-1,iz-1,p,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
      	             if (ix < nx_part[p+1]-1)
      	             {
      		        A_j[cnt] = map3(ix+1,iy-1,iz-1,p,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else
      	             {
      	                if (ix+1 < nx) 
      	                {
      		           A_j[cnt] = map3(ix+1,iy-1,iz-1,p+1,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	             }
                  }
                  else
                  {
                     if (iy) 
                     {
      	                if (ix > nx_part[p])
      	                {
      		           A_j[cnt] = map3(ix-1,iy-1,iz-1,p,q-1,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	                else if (ix)
      	                {
      		           A_j[cnt] =map3(ix-1,iy-1,iz-1,p-1,q-1,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      		        A_j[cnt] = map3(ix,iy-1,iz-1,p,q-1,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	                if (ix < nx_part[p+1]-1)
      	                {
      		           A_j[cnt] = map3(ix+1,iy-1,iz-1,p,q-1,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	                else if (ix < nx-1)
      	                {
      		           A_j[cnt] =map3(ix+1,iy-1,iz-1,p+1,q-1,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
                     }
                  }
                  if (ix > nx_part[p]) 
                  {
      		     A_j[cnt] =map3(ix-1,iy,iz-1,p,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
                  }
                  else
                  {
                     if (ix) 
                     {
      		        A_j[cnt] =map3(ix-1,iy,iz-1,p-1,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
                     }
                  }
      		  A_j[cnt] =map3(ix,iy,iz-1,p,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		  A_data[cnt++] = value[1];
                  if (ix+1 < nx_part[p+1]) 
                  {
      		     A_j[cnt] =map3(ix+1,iy,iz-1,p,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
                  }
                  else
                  {
                     if (ix+1 < nx) 
                     {
      		        A_j[cnt] =map3(ix+1,iy,iz-1,p+1,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
                     }
                  }
                  if (iy+1 < ny_part[q+1]) 
                  {
      	             if (ix > nx_part[p])
      	             {
      		        A_j[cnt] =map3(ix-1,iy+1,iz-1,p,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else
      	             {
      	                if (ix) 
      	                {
      		           A_j[cnt] =map3(ix-1,iy+1,iz-1,p-1,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	             }
      		     A_j[cnt] =map3(ix,iy+1,iz-1,p,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
      	             if (ix < nx_part[p+1]-1)
      	             {
      		        A_j[cnt] =map3(ix+1,iy+1,iz-1,p,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else
      	             {
      	                if (ix+1 < nx) 
      	                {
      		           A_j[cnt] =map3(ix+1,iy+1,iz-1,p+1,q,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	             }
                  }
                  else
                  {
                     if (iy+1 < ny) 
                     {
      	                if (ix > nx_part[p])
      	                {
      		           A_j[cnt] =map3(ix-1,iy+1,iz-1,p,q+1,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	                else if (ix)
      	                {
      		           A_j[cnt] =map3(ix-1,iy+1,iz-1,p-1,q+1,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      		        A_j[cnt] =map3(ix,iy+1,iz-1,p,q+1,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	                if (ix < nx_part[p+1]-1)
      	                {
      		           A_j[cnt] =map3(ix+1,iy+1,iz-1,p,q+1,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	                else if (ix < nx-1)
      	                {
      		           A_j[cnt] =map3(ix+1,iy+1,iz-1,p+1,q+1,r-1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
                     }
                  }
               }
            }
            if (iy > ny_part[q]) 
            {
   	       if (ix > nx_part[p])
   	       {
   	          A_j[cnt] = row_index-nx_local-1;
   	          A_data[cnt++] = value[1];
   	       }
   	       else
   	       {
   	          if (ix) 
   	          {
      		     A_j[cnt] =map3(ix-1,iy-1,iz,p-1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
   	          }
   	       }
   	       A_j[cnt] = row_index-nx_local;
   	       A_data[cnt++] = value[1];
   	       if (ix < nx_part[p+1]-1)
   	       {
   	          A_j[cnt] = row_index-nx_local+1;
   	          A_data[cnt++] = value[1];
   	       }
   	       else
   	       {
   	          if (ix+1 < nx) 
   	          {
      		     A_j[cnt] =map3(ix+1,iy-1,iz,p+1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
   	          }
   	       }
            }
            else
            {
               if (iy) 
               {
   	          if (ix > nx_part[p])
   	          {
      		     A_j[cnt] =map3(ix-1,iy-1,iz,p,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
   	          }
   	          else if (ix)
   	          {
      		     A_j[cnt] =map3(ix-1,iy-1,iz,p-1,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
   	          }
      		  A_j[cnt] =map3(ix,iy-1,iz,p,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		  A_data[cnt++] = value[1];
   	          if (ix < nx_part[p+1]-1)
   	          {
      		     A_j[cnt] =map3(ix+1,iy-1,iz,p,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
   	          }
   	          else if (ix < nx-1)
   	          {
      		     A_j[cnt] =map3(ix+1,iy-1,iz,p+1,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
   	          }
               }
            }
            if (ix > nx_part[p]) 
            {
               A_j[cnt] = row_index-1;
               A_data[cnt++] = value[1];
            }
            else
            {
               if (ix) 
               {
      		  A_j[cnt] =map3(ix-1,iy,iz,p-1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		  A_data[cnt++] = value[1];
               }
            }
            if (ix+1 < nx_part[p+1]) 
            {
               A_j[cnt] = row_index+1;
               A_data[cnt++] = value[1];
            }
            else
            {
               if (ix+1 < nx) 
               {
      		  A_j[cnt] =map3(ix+1,iy,iz,p+1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		  A_data[cnt++] = value[1];
               }
            }
            if (iy+1 < ny_part[q+1]) 
            {
   	       if (ix > nx_part[p])
   	       {
                  A_j[cnt] = row_index+nx_local-1;
                  A_data[cnt++] = value[1];
   	       }
   	       else
   	       {
   	          if (ix) 
                  {
      		     A_j[cnt] =map3(ix-1,iy+1,iz,p-1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
                  }
   	       }
               A_j[cnt] = row_index+nx_local;
               A_data[cnt++] = value[1];
   	       if (ix < nx_part[p+1]-1)
   	       {
                  A_j[cnt] = row_index+nx_local+1;
                  A_data[cnt++] = value[1];
   	       }
   	       else
   	       {
   	          if (ix+1 < nx) 
                  {
      		     A_j[cnt] =map3(ix+1,iy+1,iz,p+1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
                  }
   	       }
            }
            else
            {
               if (iy+1 < ny) 
               {
   	          if (ix > nx_part[p])
   	          {
      		     A_j[cnt] =map3(ix-1,iy+1,iz,p,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
   	          }
   	          else if (ix)
   	          {
      		     A_j[cnt] =map3(ix-1,iy+1,iz,p-1,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
   	          }
      		  A_j[cnt] =map3(ix,iy+1,iz,p,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		  A_data[cnt++] = value[1];
   	          if (ix < nx_part[p+1]-1)
   	          {
      		     A_j[cnt] =map3(ix+1,iy+1,iz,p,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
   	          }
   	          else if (ix < nx-1)
   	          {
      		     A_j[cnt] =map3(ix+1,iy+1,iz,p+1,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
   	          }
               }
            }
            if (iz+1 < nz_part[r+1]) 
            {
               if (iy > ny_part[q]) 
               {
      	          if (ix > nx_part[p])
      	          {
      	             A_j[cnt] = row_index+nxy-nx_local-1;
      	             A_data[cnt++] = value[1];
      	          }
      	          else
      	          {
      	             if (ix) 
   	             {
      		        A_j[cnt] =map3(ix-1,iy-1,iz+1,p-1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
   	             }
      	          }
      	          A_j[cnt] = row_index+nxy-nx_local;
      	          A_data[cnt++] = value[1];
      	          if (ix < nx_part[p+1]-1)
      	          {
      	             A_j[cnt] = row_index+nxy-nx_local+1;
      	             A_data[cnt++] = value[1];
      	          }
      	          else
      	          {
      	             if (ix+1 < nx) 
   	             {
      		        A_j[cnt] =map3(ix+1,iy-1,iz+1,p+1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
   	             }
      	          }
               }
               else
               {
                  if (iy) 
                  {
      	             if (ix > nx_part[p])
      	             {
      		        A_j[cnt] =map3(ix-1,iy-1,iz+1,p,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else if (ix)
      	             {
      		        A_j[cnt] =map3(ix-1,iy-1,iz+1,p-1,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      		     A_j[cnt] =map3(ix,iy-1,iz+1,p,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
      	             if (ix < nx_part[p+1]-1)
      	             {
      		        A_j[cnt] =map3(ix+1,iy-1,iz+1,p,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else if (ix < nx-1)
      	             {
      		        A_j[cnt] =map3(ix+1,iy-1,iz+1,p+1,q-1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
                  }
               }
               if (ix > nx_part[p]) 
               {
                  A_j[cnt] = row_index+nxy-1;
                  A_data[cnt++] = value[1];
               }
               else
               {
                  if (ix) 
                  {
      		     A_j[cnt] =map3(ix-1,iy,iz+1,p-1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
                  }
               }
               A_j[cnt] = row_index+nxy;
               A_data[cnt++] = value[1];
               if (ix+1 < nx_part[p+1]) 
               {
                  A_j[cnt] = row_index+nxy+1;
                  A_data[cnt++] = value[1];
               }
               else
               {
                  if (ix+1 < nx) 
                  {
      		     A_j[cnt] =map3(ix+1,iy,iz+1,p+1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
                  }
               }
               if (iy+1 < ny_part[q+1]) 
               {
      	          if (ix > nx_part[p])
      	          {
                     A_j[cnt] = row_index+nxy+nx_local-1;
                     A_data[cnt++] = value[1];
      	          }
      	          else
      	          {
      	             if (ix) 
                     {
      		        A_j[cnt] =map3(ix-1,iy+1,iz+1,p-1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
                     }
      	          }
                  A_j[cnt] = row_index+nxy+nx_local;
                  A_data[cnt++] = value[1];
      	          if (ix < nx_part[p+1]-1)
      	          {
                     A_j[cnt] = row_index+nxy+nx_local+1;
                     A_data[cnt++] = value[1];
      	          }
      	          else
      	          {
      	             if (ix+1 < nx) 
                     {
      		        A_j[cnt] =map3(ix+1,iy+1,iz+1,p+1,q,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
                     }
      	          }
               }
               else
               {
                  if (iy+1 < ny) 
                  {
      	             if (ix > nx_part[p])
      	             {
      		        A_j[cnt] =map3(ix-1,iy+1,iz+1,p,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else if (ix)
      	             {
      		        A_j[cnt] =map3(ix-1,iy+1,iz+1,p-1,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      		     A_j[cnt] =map3(ix,iy+1,iz+1,p,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
      	             if (ix < nx_part[p+1]-1)
      	             {
      		        A_j[cnt] =map3(ix+1,iy+1,iz+1,p,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else if (ix < nx-1)
      	             {
      		        A_j[cnt] =map3(ix+1,iy+1,iz+1,p+1,q+1,r,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
                  }
               }
            }
            else
            {
               if (iz+1 < nz)
	       {
                  if (iy > ny_part[q]) 
                  {
      	             if (ix > nx_part[p])
      	             {
      		        A_j[cnt] =map3(ix-1,iy-1,iz+1,p,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else
      	             {
      	                if (ix) 
      	                {
      		           A_j[cnt] =map3(ix-1,iy-1,iz+1,p-1,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	             }
      		     A_j[cnt] =map3(ix,iy-1,iz+1,p,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
      	             if (ix < nx_part[p+1]-1)
      	             {
      		        A_j[cnt] =map3(ix+1,iy-1,iz+1,p,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else
      	             {
      	                if (ix+1 < nx) 
      	                {
      		           A_j[cnt] =map3(ix+1,iy-1,iz+1,p+1,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	             }
                  }
                  else
                  {
                     if (iy) 
                     {
      	                if (ix > nx_part[p])
      	                {
      		           A_j[cnt] =map3(ix-1,iy-1,iz+1,p,q-1,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	                else if (ix)
      	                {
      		           A_j[cnt] =map3(ix-1,iy-1,iz+1,p-1,q-1,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      		        A_j[cnt] =map3(ix,iy-1,iz+1,p,q-1,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	                if (ix < nx_part[p+1]-1)
      	                {
      		           A_j[cnt] =map3(ix+1,iy-1,iz+1,p,q-1,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	                else if (ix < nx-1)
      	                {
      		           A_j[cnt] =map3(ix+1,iy-1,iz+1,p+1,q-1,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
                     }
                  }
                  if (ix > nx_part[p]) 
                  {
      		     A_j[cnt] =map3(ix-1,iy,iz+1,p,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
                  }
                  else
                  {
                     if (ix) 
                     {
      		        A_j[cnt] =map3(ix-1,iy,iz+1,p-1,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
                     }
                  }
      		  A_j[cnt] =map3(ix,iy,iz+1,p,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		  A_data[cnt++] = value[1];
                  if (ix+1 < nx_part[p+1]) 
                  {
      		     A_j[cnt] =map3(ix+1,iy,iz+1,p,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
                  }
                  else
                  {
                     if (ix+1 < nx) 
                     {
      		        A_j[cnt] =map3(ix+1,iy,iz+1,p+1,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
                     }
                  }
                  if (iy+1 < ny_part[q+1]) 
                  {
      	             if (ix > nx_part[p])
      	             {
      		        A_j[cnt] =map3(ix-1,iy+1,iz+1,p,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else
      	             {
      	                if (ix) 
      	                {
      		           A_j[cnt] =map3(ix-1,iy+1,iz+1,p-1,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	             }
      		     A_j[cnt] =map3(ix,iy+1,iz+1,p,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		     A_data[cnt++] = value[1];
      	             if (ix < nx_part[p+1]-1)
      	             {
      		        A_j[cnt] =map3(ix+1,iy+1,iz+1,p,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	             }
      	             else
      	             {
      	                if (ix+1 < nx) 
      	                {
      		           A_j[cnt] =map3(ix+1,iy+1,iz+1,p+1,q,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	             }
                  }
                  else
                  {
                     if (iy+1 < ny) 
                     {
      	                if (ix > nx_part[p])
      	                {
      		           A_j[cnt] =map3(ix-1,iy+1,iz+1,p,q+1,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	                else if (ix)
      	                {
      		           A_j[cnt] =map3(ix-1,iy+1,iz+1,p-1,q+1,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      		        A_j[cnt] =map3(ix,iy+1,iz+1,p,q+1,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		        A_data[cnt++] = value[1];
      	                if (ix < nx_part[p+1]-1)
      	                {
      		           A_j[cnt] =map3(ix+1,iy+1,iz+1,p,q+1,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
      	                else if (ix < nx-1)
      	                {
      		           A_j[cnt] =map3(ix+1,iy+1,iz+1,p+1,q+1,r+1,P,Q,R,
					nx_part,ny_part,nz_part,global_part);
      		           A_data[cnt++] = value[1];
      	                }
                     }
                  }
               }
            }
            row_index++;
         }
      }
   }
   }
   }
   }

   A = hypre_CSRMatrixCreate(grid_size, grid_size, A_i[grid_size]);

   hypre_CSRMatrixI(A) = A_i;
   hypre_CSRMatrixJ(A) = A_j;
   hypre_CSRMatrixData(A) = A_data;

   hypre_TFree(nx_part);
   hypre_TFree(ny_part);
   hypre_TFree(nz_part);
   hypre_TFree(global_part);

   return A;
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

HYPRE_Int
map3( HYPRE_Int  ix,
      HYPRE_Int  iy,
      HYPRE_Int  iz,
      HYPRE_Int  p,
      HYPRE_Int  q,
      HYPRE_Int  r,
      HYPRE_Int  P,
      HYPRE_Int  Q,
      HYPRE_Int  R,
      HYPRE_Int *nx_part,
      HYPRE_Int *ny_part,
      HYPRE_Int *nz_part,
      HYPRE_Int *global_part )
{
   HYPRE_Int nx_local;
   HYPRE_Int ix_local;
   HYPRE_Int iy_local;
   HYPRE_Int iz_local;
   HYPRE_Int nxy;
   HYPRE_Int global_index;
   HYPRE_Int proc_num;
 
   proc_num = r*P*Q + q*P + p;
   nx_local = nx_part[p+1] - nx_part[p];
   nxy = nx_local*(ny_part[q+1] - ny_part[q]);
   ix_local = ix - nx_part[p];
   iy_local = iy - ny_part[q];
   iz_local = iz - nz_part[r];
   global_index = global_part[proc_num] 
      + iz_local*nxy + iy_local*nx_local + ix_local;

   return global_index;
}
