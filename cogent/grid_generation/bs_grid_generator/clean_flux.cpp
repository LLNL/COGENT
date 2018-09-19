#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NUM_BLOCKS 10
#define MAX_NEWTON_ITER  20
#define NEWTON_TOL       1.e-13

#define Pi (3.14159265358979323846264338327950288e0)

using namespace std;


class Flux {
 public:

   ~Flux()
      {
         if (data) delete [] data;
      }

   void initialize(char *flux_file)
   {
      FILE* fd_flux;

      if ( (fd_flux = fopen( flux_file, "r" )) == NULL) {
         std::cout << "Cannot open flux file" << std::endl;
         exit(1);
      }

      double magAxisR, magAxisZ, XpointR, XpointZ;

      fscanf(fd_flux, "%lf %lf %lf %lf %d %d", &Rmin, &Rmax, &Zmin, &Zmax, &NR, &NZ);
      fscanf(fd_flux, "%lf %lf %lf %lf", &magAxisR, &magAxisZ, &XpointR, &XpointZ);

      coef = new double*[NR];
      for (int i=0; i<NR; ++i) {
         coef[i] = new double[NZ];
      }

      for (int j=0; j<NZ; ++j) {
         for (int i=0; i<NR; ++i) {
            fscanf(fd_flux, "%lf", &coef[i][j]);
         }
      }

      fclose(fd_flux);
      
      int llen = NR>NZ? NR: NZ;
      double * data = new double[7*llen];
      lambda = data;
      facR = lambda + llen;
      facZ = facR + llen;
      sinfacR = facZ + llen;
      cosfacR = sinfacR + llen;
      sinfacZ = cosfacR + llen;
      cosfacZ = sinfacZ + llen;
      
      lambda[0] = 1. / sqrt(2.);
      for (int l=1; l<llen; ++l) {
         lambda[l] = 1.;
      }
      
      for (int i=0; i<NR; ++i) {
         facR[i] = i * Pi / NR;
      }

      for (int j=0; j<NZ; ++j) {
         facZ[j] = j * Pi / NZ;
      }

      Rscale = (NR-1)/(Rmax - Rmin);
      Zscale = (NZ-1)/(Zmax - Zmin);
   }

   double evaluate( const double&  R,
                    const double&  Z,
                    const int      derivR,
                    const int      derivZ)
   {
      double sR = (R - Rmin) * Rscale;
      double sZ = (Z - Zmin) * Zscale;

      for (int i=0; i<NR; ++i) {
         double t = facR[i] * (sR + 0.5);
         sinfacR[i] = sin(t);
         cosfacR[i] = cos(t);
      }

      for (int j=0; j<NZ; ++j) {
         double t = facZ[j] * (sZ + 0.5);
         sinfacZ[j] = sin(t);
         cosfacZ[j] = cos(t);
      }

      double value = 0.;

      if (derivR == 0 && derivZ == 0) {
         for (int j=0; j<NZ; ++j) {
            for (int i=0; i<NR; ++i) {
               value += lambda[i] * lambda[j] * coef[i][j] * cosfacR[i] * cosfacZ[j];
            }
         }
      }
      else if (derivR == 1 && derivZ == 0) {
         for (int j=0; j<NZ; ++j) {
            for (int i=0; i<NR; ++i) {
               value -= lambda[i] * lambda[j] * coef[i][j] * sinfacR[i] * cosfacZ[j] * facR[i];
            }
         }
         value *= Rscale;
      }
      else if (derivR == 0 && derivZ == 1) {
         for (int j=0; j<NZ; ++j) {
            for (int i=0; i<NR; ++i) {
               value -= lambda[i] * lambda[j] * coef[i][j] * cosfacR[i] * sinfacZ[j] * facZ[j];
            }
         }
         value *= Zscale;
      }
      else {
         std::cerr << "Illegal derivative flag in evalflux()" << std::endl;
         exit(1);
      }
   
      value *= 2. / sqrt(NR*NZ);

      return value;
   };

   bool pointValid( const double& R,
                    const double& Z ) const
   {
      return R >= Rmin && R <= Rmax && Z >= Zmin && Z <= Zmax;
   }

 private:

   double*  data;
   double*  facR;
   double*  facZ;
   double*  lambda;
   double** coef;
   double*  sinfacR;
   double*  cosfacR;
   double*  sinfacZ;
   double*  cosfacZ;

   int      NR;
   int      NZ;
   double   Rmin;
   double   Rmax;
   double   Zmin;
   double   Zmax;
   double   Rscale;
   double   Zscale;
} flux_data;


bool isCoreBlock(const char* block_name)
{
   return strcmp(block_name,"lcore") == 0 ||
          strcmp(block_name,"rcore") == 0 ||
          strcmp(block_name,"mcore") == 0;
}


bool isSOLBlock(const char* block_name)
{
   return strcmp(block_name,"lsol")  == 0 ||
          strcmp(block_name,"lcsol") == 0 ||
          strcmp(block_name,"mcsol") == 0 ||
          strcmp(block_name,"rcsol") == 0 ||
          strcmp(block_name,"rsol")  == 0;
}


bool isPfBlock(const char* block_name)
{
   return strcmp(block_name,"lpf") == 0 ||
          strcmp(block_name,"rpf") == 0;
}


void adjustFlux(double&  R,
                double&  Z,
                const    double& target_flux,
                int&     niter,
                double&  change)
{
   double Rorig = R;
   double Zorig = Z;

   double Rnew = Rorig;
   double Znew = Zorig;
   double s = 0.;

   niter = 0;
   change = 0.;

   // Compute flux surface normal
   double normalR = flux_data.evaluate(Rnew, Znew, 1, 0);
   double normalZ = flux_data.evaluate(Rnew, Znew, 0, 1);
   double len = sqrt(normalR*normalR + normalZ*normalZ);
   normalR /= len;
   normalZ /= len;

   while(niter < MAX_NEWTON_ITER) {

      double psi = flux_data.evaluate(Rnew, Znew, 0, 0);

      double diff = psi - target_flux;

      change = sqrt(pow(Rnew - Rorig,2) + pow(Znew - Zorig,2));

      if (fabs(diff) < NEWTON_TOL * fabs(target_flux)) {
         break;
      }
      else {
         niter++;

         double dpsidR = flux_data.evaluate(Rnew, Znew, 1, 0);
         double dpsidZ = flux_data.evaluate(Rnew, Znew, 0, 1);
                     
         s -= diff / (dpsidR*normalR + dpsidZ*normalZ);

         Rnew = Rorig + s*normalR;
         Znew = Zorig + s*normalZ;

         if ( !flux_data.pointValid(Rnew,Znew) ) {
            cerr << "A Newton iterate has fallen outside of the function domain" << endl;
            exit(1);
         }
      }
   }

   R = Rnew;
   Z = Znew;
}


int main(int argc, char** argv)
{
   if (argc != 3) {
      printf("Usage: fixflux COGENT_mapping flux_coefficients\n");
      exit(1);
   }

   // Read and initialize the flux data
   flux_data.initialize(argv[2]);

   FILE* fd_mapping;

   if ( (fd_mapping = fopen( argv[1], "r" )) == NULL) {
      cerr << "Cannot open mapping file" << endl;
      exit(1);
   }

   char block[NUM_BLOCKS][80];

   double*** R; 
   double*** Z; 

   R = new double**[NUM_BLOCKS];
   Z = new double**[NUM_BLOCKS];

   int nradial[NUM_BLOCKS], nradial_extend[NUM_BLOCKS], npoloidal[NUM_BLOCKS], npoloidal_extend[NUM_BLOCKS];

   for (int block_number=0; block_number<NUM_BLOCKS; ++block_number) {

      fscanf(fd_mapping, "%s %d %d %d %d\n", block[block_number], &nradial[block_number], &nradial_extend[block_number],
             &npoloidal[block_number], &npoloidal_extend[block_number]);

      //      cout << block[block_number] << " " << nradial[block_number] << " " << nradial_extend[block_number] <<
      //         " " << npoloidal[block_number] << " " << npoloidal_extend[block_number] << endl;

      R[block_number] = new double*[nradial[block_number]];
      Z[block_number] = new double*[nradial[block_number]];
      for (int i=0; i<nradial[block_number]; ++i) {
         R[block_number][i] = new double[npoloidal[block_number]];
         Z[block_number][i] = new double[npoloidal[block_number]];
      }

      for (int j=0; j<npoloidal[block_number]; ++j) {
         for (int i=0; i<nradial[block_number]; ++i) {
            double dummy1, dummy2;
            fscanf(fd_mapping, "%lf %lf %lf %lf", &R[block_number][i][j], &Z[block_number][i][j], &dummy1, &dummy2);
            //         cout << R[block_number][i][j] << " " << Z[block_number][i][j] << " " << dummy1 << " " << dummy2 << endl;
         }
      }
   }

   fclose(fd_mapping);

   // Average the core fluxes

   int core_nradial = -1;

   for (int block_number=0; block_number<NUM_BLOCKS; ++block_number) {
      if ( isCoreBlock(block[block_number]) ) {
         if (core_nradial == -1) {
            core_nradial = nradial[block_number];
         } else if (core_nradial != nradial[block_number]) {
            cerr << "Inconsistent radial dimension" << endl;
         }
      }
   }

   double* core_flux = new double[core_nradial];

   for (int i=0; i<core_nradial; ++i) {
      core_flux[i] = 0.;
      int poloidal_vals = 0;

      for (int block_number=0; block_number<NUM_BLOCKS; ++block_number) {
         if ( isCoreBlock(block[block_number]) ) {
            for (int j=0; j<npoloidal[block_number]; ++j) {
               core_flux[i] += flux_data.evaluate(R[block_number][i][j], Z[block_number][i][j], 0, 0);
               poloidal_vals++;
            }
         }
      }

      core_flux[i] /= (double)poloidal_vals;
   }

   // Average the scrape-off fluxes

   int sol_nradial = -1;

   for (int block_number=0; block_number<NUM_BLOCKS; ++block_number) {
      if ( isSOLBlock(block[block_number]) ) {
         if (sol_nradial == -1) {
            sol_nradial = nradial[block_number];
         } else if (sol_nradial != nradial[block_number]) {
            cerr << "Inconsistent radial dimension" << endl;
         }
      }
   }

   double* sol_flux = new double[sol_nradial];

   for (int i=0; i<sol_nradial; ++i) {
      sol_flux[i] = 0.;
      int poloidal_vals = 0;

      for (int block_number=0; block_number<NUM_BLOCKS; ++block_number) {
         if ( isSOLBlock(block[block_number]) ) {
            for (int j=0; j<npoloidal[block_number]; ++j) {
               sol_flux[i] += flux_data.evaluate(R[block_number][i][j], Z[block_number][i][j], 0, 0);
               poloidal_vals++;
            }
         }
      }

      sol_flux[i] /= (double)poloidal_vals;
   }

   // Average the private flux fluxes

   int pf_nradial = -1;

   for (int block_number=0; block_number<NUM_BLOCKS; ++block_number) {
      if ( isPfBlock(block[block_number]) ) {
         if (pf_nradial == -1) {
            pf_nradial = nradial[block_number];
         } else if (pf_nradial != nradial[block_number]) {
            cerr << "Inconsistent radial dimension" << endl;
         }
      }
   }

   double* pf_flux = new double[pf_nradial];

   for (int i=0; i<pf_nradial; ++i) {
      pf_flux[i] = 0.;
      int poloidal_vals = 0;

      for (int block_number=0; block_number<NUM_BLOCKS; ++block_number) {
         if ( isPfBlock(block[block_number]) ) {
            for (int j=0; j<npoloidal[block_number]; ++j) {
               pf_flux[i] += flux_data.evaluate(R[block_number][i][j], Z[block_number][i][j], 0, 0);
               poloidal_vals++;
            }
         }
      }

      pf_flux[i] /= (double)poloidal_vals;
   }

   // Ensure that core and pf separatrix flux agrees with the scape-off
   core_flux[core_nradial-1] = sol_flux[0];
   pf_flux[pf_nradial-1]     = sol_flux[0];

   int max_newton_iter = 0;
   double max_change = 0.;

   cout << "Adjusting core points" << endl;

   for (int i=0; i<core_nradial; ++i) {
      for (int block_number=0; block_number<NUM_BLOCKS; ++block_number) {
         if ( isCoreBlock(block[block_number]) ) {
            for (int j=0; j<npoloidal[block_number]; ++j) {

               int niter;
               double change;
               adjustFlux(R[block_number][i][j], Z[block_number][i][j], core_flux[i], niter, change);

               if (niter > max_newton_iter) max_newton_iter = niter;
               if (change > max_change) max_change = change;
            }
         }
      }
   }

   cout << "Adjusting scrape-off points" << endl;

   for (int i=0; i<sol_nradial; ++i) {
      for (int block_number=0; block_number<NUM_BLOCKS; ++block_number) {
         if ( isSOLBlock(block[block_number]) ) {
            for (int j=0; j<npoloidal[block_number]; ++j) {

               int niter;
               double change;
               adjustFlux(R[block_number][i][j], Z[block_number][i][j], sol_flux[i], niter, change);

               if (niter > max_newton_iter) max_newton_iter = niter;
               if (change > max_change) max_change = change;
            }
         }
      }
   }

   cout << "Adjusting private flux points" << endl;

   for (int i=0; i<pf_nradial; ++i) {
      for (int block_number=0; block_number<NUM_BLOCKS; ++block_number) {
         if ( isPfBlock(block[block_number]) ) {
            for (int j=0; j<npoloidal[block_number]; ++j) {

               int niter;
               double change;
               adjustFlux(R[block_number][i][j], Z[block_number][i][j], pf_flux[i], niter, change);

               if (niter > max_newton_iter) max_newton_iter = niter;
               if (change > max_change) max_change = change;
            }
         }
      }
   }

   cout << "Maximum number of Newton iterations performed = " << max_newton_iter << endl;
   cout << "Maximum point change = " << max_change << endl;

   FILE* fd_output;

   char output_file[80];
   strcpy(output_file, argv[1]);
   strcat(output_file, "_cf");

   if ( (fd_output = fopen( output_file, "w" )) == NULL) {
      cerr << "Cannot open output mapping file" << endl;
      exit(1);
   }

   for (int block_number=0; block_number<NUM_BLOCKS; ++block_number) {
      fprintf(fd_output, "%s  %d  %d  %d  %d\n", block[block_number],
              nradial[block_number], nradial_extend[block_number],
              npoloidal[block_number], npoloidal_extend[block_number]);
      for (int j=0; j<npoloidal[block_number]; ++j) {
         for (int i=0; i<nradial[block_number]; ++i) {
            fprintf(fd_output, "%20.13e  %20.13e  0.  0.\n", R[block_number][i][j], Z[block_number][i][j]);
         }
      }
   }

   fclose(fd_output);
}
