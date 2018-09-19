#include <iostream>
#include <stdio.h>
#include <cmath>

using namespace std;
 
// function declaration:
double getPsi(double r, double z);

int main ()
{

  /* 
   Compute DCT coefficients for the "model" geometry
  */

  double PI = 3.141592;

  char file_name[80];
  sprintf(file_name, "DCT_coefficients.txt");
  FILE* fd = fopen(file_name, "w");

  int N = 128;
  double Rmin = 0.5; 
  double Rmax = 2.6; 
  double Zmin = -1.0;
  double Zmax = 3.5;

  double m_a =  1.2;
  double m_b =  0.9;
  double m_c =  0.7;
  double m_R0 = 1.6;
  double m_Zx = -acos(m_c/m_b);
  double m_Zc  = acos(m_c/m_b);
  double m_Z0  = 0.4-m_Zx;

  double m_Zaxis = m_Zc + m_Z0;
  double m_Zxpt = m_Zx + m_Z0;

  cout<<getPsi(m_R0,m_Zxpt)<<" "<<getPsi(Rmin,Zmin)<<" "<<getPsi(m_R0,m_Zaxis)<<endl;

  fprintf(fd, "%20.7e %20.7e %20.7e %20.7e %d %d %20.7e %20.7e %20.7e %20.7e \n", Rmin, Rmax, Zmin, Zmax, N, N, m_R0, m_Zaxis, m_R0, m_Zxpt);

  double dct_coeff;

  for (int v = 0; v < N; v++) {
    for (int u = 0; u < N; u++) {

      dct_coeff = 0.0;
      double Lambda_u(1.0);
      double Lambda_v(1.0);

      if (u == 0) Lambda_u = 1.0/sqrt(2.0);
      if (v == 0) Lambda_v = 1.0/sqrt(2.0);

      for (int i = 0; i < N; i++) {
	for (int j = 0; j < N; j++) {

	  double r = Rmin + (Rmax-Rmin)/double(N-1)*double(i);
	  double z = Zmin + (Zmax-Zmin)/double(N-1)*double(j);

	  double fac = getPsi(r,z);
	  fac *= cos(PI*double(u)/(2.0*double(N))*(2.0*double(i)+1.0));
	  fac *= cos(PI*double(v)/(2.0*double(N))*(2.0*double(j)+1.0));

	  dct_coeff += 2.0/double(N)*Lambda_u*Lambda_v*fac;

	}
      }
      fprintf(fd, "%20.12e ", dct_coeff);
    }
    fprintf(fd, "\n ");
  }

  fclose(fd);

  return 0;
}

double getPsi(double RR, double ZZ)
{

  double scaling_fac = 0.501984;
  double m_a =  1.2;
  double m_b =  0.9;
  double m_c =  0.7;
  double m_R0 = 1.6;
  double m_Zx = -acos(m_c/m_b);
  double m_Zc  = acos(m_c/m_b);
  double m_Z0  = 0.4-m_Zx;

  double result;
  double x = RR - m_R0;
  double z = ZZ - m_Z0;
   
  if (z<m_Zx) z = 2.0 * m_Zx - z;
   
  double psiMin = m_b * sin(m_Zc) - m_c * m_Zc;
  double psi = cos(m_a * x) + m_b * sin(z) - m_c * z;
   
  //Bacause psi is periodic (though periodicity is seen only far outside the domain)
  //we flatten it out here, baising to a constant after the "first period boundary".
  //This is the extra cautionary measure for various Newton solvers.
  double psiTruncated = (psi < psiMin) ? psiMin : psi;

  psi *= scaling_fac;
  return psi; //psiTruncated;

}
