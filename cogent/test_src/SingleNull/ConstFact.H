#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <cmath>

#include "NamespaceHeader.H"

/**
   A namespace providing definitions of commomly used physical constants
   and conversion factors.
 */
namespace Constants {

   /**
     Fundamental Physical Constants
   */
   const double SPEED_OF_LIGHT      = 2.99792458e+08; // m/s
   const double MASS_OF_ELECTRON    = 9.1093826e-31;  // kg
   const double MASS_OF_NEUTRON     = 1.67492728e-27; // kg
   const double MASS_OF_PROTON      = 1.67262171e-27; // kg
   const double ELEMENTARY_CHARGE   = 1.60217653e-19; // C
   const double VACUUM_PERMITTIVITY = 8.85418782e-12; // C^2/N-m^2
   const double BOLTZMANN_CONSTANT  = 1.3806505e-23;  // J/K
   const double PLANCK_CONSTANT     = 6.6260693e-34;  // J s

   /**
     Mathematical constants
   */
   const double PI     = M_PI;
   const double TWO_PI = 2 * M_PI;
   const double E      = exp( 1.0 );

   /**
     Conversion Factors
   */
   const double JOULE_PER_KEV     = 1.60217646e-16;
   const double KEV_PER_JOULE     = 1.0 / JOULE_PER_KEV;
   const double CM3_PER_M3        = 1.0e+06;
   const double M3_PER_CM3        = 1.0e-06;
   const double SECOND_PER_HOUR   = 3.6e+03;
   const double SECOND_PER_MINUTE = 6.0e+01;

}

#include "NamespaceFooter.H"

#endif
