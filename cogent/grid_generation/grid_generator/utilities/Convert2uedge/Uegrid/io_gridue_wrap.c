/*
C wrappers for for fortran routines READG, WRITEG etc
Fortran interface to IDL is platform specific, it
is easier to link via C wrappers
*/

#include <stdio.h>  
#include <stdlib.h>
//#include "export.h" /* IDL export definitions */
#include "idl_export.h" /* IDL export definitions */


#if UNDERSCORE
 #define WRITE_GRIDUE write_gridue_
 #define READ_GRIDUE read_gridue_
 #define READ_GRIDUE_DIMS read_gridue_dims_
#else
 #define WRITE_GRIDUE write_gridue
 #define READ_GRIDUE read_gridue
 #define READ_GRIDUE_DIMS read_gridue_dims
#endif




int write_gridue(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void WRITE_GRIDUE(int nxm, int nym, int ixpt1, int ixpt2, int iysptrx1,
			double* rm, double* zm, double* psi, 
			double* br, double* bz, double* bpol, 
			double* bphi, double* b);    
  //int ix,iy;
  //double (*rm)[13][3][5]=argv[5];

  if (argc != 13){
    fprintf(stderr,"write_gridue: Incorrect number of arguments\n");
    return(0);
  }

  //printf("In C: argv[0]=%d\n", *((int *) argv[0]));
  //printf("In C: argv[1]=%d\n", *((int *) argv[1]));
  //printf("In C: argv[2]=%d\n", *((int *) argv[2]));
  //printf("In C: argv[3]=%d\n", *((int *) argv[3]));
  //printf("In C: argv[4]=%d\n", *((int *) argv[4]));

  //printf("In C: rm[0]=%f\n", ((double*) argv[5])[0]);



  WRITE_GRIDUE(
	    (int) argv[0], (int) argv[1], (int) argv[2], (int) argv[3], (int) argv[4],
	    (double*) argv[5], (double*) argv[6], (double*) argv[7], 
	    (double*) argv[8], (double*) argv[9], (double*) argv[10], 
	    (double*) argv[11], (double*) argv[12]
	    );

  return (1);
}  



int read_gridue(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void READ_GRIDUE(int nxm, int nym, int ixpt1, int ixpt2, int iysptrx1,
			double* rm, double* zm, double* psi, 
			double* br, double* bz, double* bpol, 
			double* bphi, double* b);    

  if (argc != 13){
    fprintf(stderr,"read_gridue: Incorrect number of arguments\n");
    return(0);
  }

  READ_GRIDUE(
	    (int) argv[0], (int) argv[1], (int) argv[2], (int) argv[3], (int) argv[4],
	    (double*) argv[5], (double*) argv[6], (double*) argv[7], 
	    (double*) argv[8], (double*) argv[9], (double*) argv[10], 
	    (double*) argv[11], (double*) argv[12]
	    );

  //printf("...back to C...\n");

  return (1);
}  



int read_gridue_dims(int argc, void *argv[])  
{  

  /* Fortran routine, may need _ at the end */  
  extern void READ_GRIDUE_DIMS(int nxm, int nym);    

  //printf("In read_gridue_dims C wrapper\n");

  if (argc != 2){
    printf("read_gridue_dims: Incorrect number of arguments\n");
    fprintf(stderr,"read_gridue_dims: Incorrect number of arguments\n");
    return(0);
  }

  //printf("Before READ_GRIDUE_DIMS...\n");
  READ_GRIDUE_DIMS( (int) argv[0], (int) argv[1]);

  return (1);
}  
