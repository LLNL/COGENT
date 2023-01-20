/*
  This code reads in ODEVectors from two files
  (ASCII or binary) and computes the norms of
  the difference between the two.
  
  The vectors must be of the same size.
  
  The code will read an input file "diff.inp"
  with the following inputs:
  - name of first ODEVector file 
  - name of second ODEVector file 
  - format (ascii/binary)

  The norms are reported in the following order:
  L1, L2, Linf
  
  This is a serial code.
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>

int readODEVector(  const std::string&,
                    std::vector<double>&,
                    const bool );

void computeNorms(  const std::vector<double>&,
                    const std::vector<double>&,
                    double* const );

int main() 
{
  std::string fname1, fname2, format;

  FILE* inp;
  inp = fopen("diff.inp", "r");
  if (!inp) {
    fprintf(stderr, "Error: diff.inp not found.\n");
    return 0;
  }
  {
    char c_str[100];
    fscanf(inp, "%s", c_str);
    fname1 = std::string(c_str);
    fscanf(inp, "%s", c_str);
    fname2 = std::string(c_str);
    fscanf(inp, "%s", c_str);
    format = std::string(c_str);
  }
  fclose(inp);

  printf("ODEVector filenames are:\n");
  printf("  %s\n", fname1.c_str() );
  printf("  %s\n", fname2.c_str() );
  printf(" Format: %s\n", format.c_str() );

  bool is_ascii( format == "ascii");

  int ierr;
  std::vector<double> vec1(0), vec2(0);
  ierr = readODEVector( fname1, vec1, is_ascii );
  if (ierr) return ierr;
  ierr = readODEVector( fname2, vec2, is_ascii );
  if (ierr) return ierr;

  if ( vec1.size() != vec2.size() ) {
    fprintf(stderr, "Error: vector sizes differ (%d, %d)\n",
            vec1.size(), vec2.size() );
  }
  printf("ODEVector size: %d\n", vec1.size());

  double norms[3];
  computeNorms(vec1, vec2, norms);
  printf( "diff: %1.6e, %1.6e, %1.6e\n",
          norms[0], norms[1], norms[2] );

  return 0;
}

int readODEVector(  const std::string& a_fname,
                    std::vector<double>& a_vec,
                    const bool a_is_ascii )
{
  FILE* in;

  if (a_is_ascii) {
    in = fopen( a_fname.c_str(), "r" );
  } else {
    in = fopen( a_fname.c_str(), "rb" );
  }

  if (!in) {
    fprintf(stderr, "Error: unable to open file %s.\n", a_fname.c_str() );
    return 1;
  }
  printf("Reading file %s.\n", a_fname.c_str() );

  int vecsize;
  if (a_is_ascii) {
    fscanf( in, "%d", &vecsize ); 
  } else {
    fread( &vecsize, sizeof(int), 1, in );
  }

  a_vec.resize(vecsize);

  if (a_is_ascii) {
    for (int i=0; i<vecsize; i++) {
      fscanf(in, "%lf", &a_vec.data()[i] );
    }
  } else {
    fread(a_vec.data(), sizeof(double), vecsize, in );
  }

  fclose(in);

  return 0;
}

void computeNorms(  const std::vector<double>& a_v1,
                    const std::vector<double>& a_v2,
                    double* const a_norms )
{
  if (a_v1.size() != a_v2.size()) {
    a_norms[0] = a_norms[1] = a_norms[2] = -1.0;
  }

  double array_sum_abs = 0;
  double array_sum_square = 0;
  double array_max = 0;
  for (int i=0; i<a_v1.size(); i++) {
    double delta(a_v1[i]-a_v2[i]);
    array_sum_abs += (delta < 0 ? -delta : delta );
    array_sum_square += delta*delta;
    if (delta > array_max) array_max = delta;
  }
  a_norms[0] = array_sum_abs / ((double)a_v1.size());
  a_norms[1] = sqrt(array_sum_square/ ((double)a_v1.size()));
  a_norms[2] = array_max;

  return;
}
