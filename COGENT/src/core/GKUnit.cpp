#include "GKUnit.H"
#include "ConstFact.H"
#include <assert.h>

#include "NamespaceHeader.H"

inline void getPosDefUnit(
   Real& a_val, const string &a_name, ParmParse& a_pp )
{
   Real val = 0.0;
   a_pp.get( a_name.c_str(), val ); // a_name must exist or abort
   assert( val > 0.0 );
   a_val = val;
}

GKUnits::GKUnits()
{
   // Fundamental Characteristic Scales
   ParmParse ppunits( "units" );
   getPosDefUnit( m_scale[NUMBER_DENSITY], "number_density", ppunits );
   getPosDefUnit( m_scale[TEMPERATURE],    "temperature",    ppunits );
   getPosDefUnit( m_scale[LENGTH],         "length",         ppunits );
   getPosDefUnit( m_scale[MASS],           "mass",           ppunits );
   getPosDefUnit( m_scale[MAGNETIC_FIELD], "magnetic_field", ppunits );

   // Universal Constants
   Real pi = Constants::PI;
   m_scale[CHARGE]       = Constants::ELEMENTARY_CHARGE;
   m_scale[BOLTZMANN]    = Constants::BOLTZMANN_CONSTANT;
   m_scale[PERMITTIVITY] = Constants::VACUUM_PERMITTIVITY;
   m_scale[PERMEABILITY] = pi * 4.0e-7;

   // Derived Characteristic Scales
   Real masskg = m_scale[MASS]*Constants::MASS_OF_PROTON;
   Real tempJoules = m_scale[CHARGE] * m_scale[TEMPERATURE];
   m_scale[THERMAL_SPEED] = sqrt( tempJoules / masskg );
   m_scale[TRANSIT_TIME] = m_scale[LENGTH] / m_scale[THERMAL_SPEED];
   m_scale[MAGNETIC_MOMENT] = tempJoules / ( 2.0 * m_scale[MAGNETIC_FIELD] );
   m_scale[POTENTIAL] = tempJoules / m_scale[CHARGE];

   Real denom = pi * pow(m_scale[THERMAL_SPEED],3);
   m_scale[DISTRIBUTION_FUNCTION] = m_scale[NUMBER_DENSITY] / denom;

   m_scale[GYROFREQUENCY] = m_scale[CHARGE] * m_scale[MAGNETIC_FIELD] / masskg;
   m_scale[GYRORADIUS] = m_scale[THERMAL_SPEED] / m_scale[GYROFREQUENCY];

   denom = m_scale[NUMBER_DENSITY] * m_scale[CHARGE] * m_scale[CHARGE];
   m_scale[DEBYE_LENGTH] = sqrt( m_scale[PERMITTIVITY] * tempJoules / denom );

   // Nondmensional Parameters
   m_larmor_number = m_scale[GYRORADIUS] / m_scale[LENGTH];
   m_debye_number = m_scale[DEBYE_LENGTH] / m_scale[LENGTH];

}


inline void multLDFArrayBox( LevelData<FArrayBox>& a_data,
                             const Interval&       a_components,
                             const Real&           a_factor )
{
   DataIterator dit = a_data.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      a_data[dit].mult( a_factor, a_components.begin(), a_components.size() );
   }
}

void GKUnits::makeDimensionless( LevelData<FArrayBox>& a_data,
                                 const Interval&       a_components,
                                 const CHAR_SCALE&     a_type ) const
{
   checkType( a_type );
   multLDFArrayBox( a_data, a_components, 1.0 / m_scale[a_type] );
}

void GKUnits::makeDimensional( LevelData<FArrayBox>& a_data,
                               const Interval&       a_components,
                               const CHAR_SCALE&     a_type ) const
{
   checkType( a_type );
   multLDFArrayBox( a_data, a_components, m_scale[a_type] );
}


inline void multLDFluxBox( LevelData<FluxBox>& a_data,
                           const Interval      a_components,
                           const Real&         a_factor )
{
   DataIterator dit = a_data.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      for (int dir=0; dir<CH_SPACEDIM; dir++)
      {
         a_data[dit][dir].mult( a_factor, a_components.begin(), a_components.size() );
      }
   }
}

void GKUnits::makeDimensionless( LevelData<FluxBox>& a_data,
                                 const Interval&     a_components,
                                 const CHAR_SCALE&   a_type ) const
{
   checkType( a_type );
   multLDFluxBox( a_data, a_components, 1.0 / m_scale[a_type] );
}

void GKUnits::makeDimensional( LevelData<FluxBox>& a_data,
                               const Interval&     a_components,
                               const CHAR_SCALE&   a_type ) const
{
   checkType( a_type );
   multLDFluxBox( a_data, a_components, m_scale[a_type] );
}

void GKUnits::print( std::ostream& a_out ) const
{
   a_out << "** GKUnits ***********************************************************" << endl;
   a_out << endl;
   a_out << "==== Fundamental Units ===============================================" << endl;
   a_out << "  NUMBER DENSITY    [1/m^3]: " << m_scale[NUMBER_DENSITY] << endl;
   a_out << "  TEMPERATURE           [eV]: " << m_scale[TEMPERATURE] << endl;
   a_out << "  LENGTH                [m]: " << m_scale[LENGTH] << endl;
   a_out << "  MASS                 [AMU]: " << m_scale[MASS] << endl;
   a_out << "  MAGNETIC FIELD        [T]: " << m_scale[MAGNETIC_FIELD] << endl;
   a_out << "  CHARGE                [C]: " << m_scale[CHARGE] << endl;
   a_out << "  BOLTZMANN CONSTANT  [J/K]: " << m_scale[BOLTZMANN] << endl;
   a_out << "  ELECTRIC CONSTANT   [F/m]: " << m_scale[PERMITTIVITY] << endl;
   a_out << "  MAGNETIC CONSTANT   [H/m]: " << m_scale[PERMEABILITY] << endl;
   a_out << endl;
   a_out << "==== Derived Units ===================================================" << endl;
   a_out << endl;
   a_out << "  THERMAL SPEED       [m/s]: " << m_scale[THERMAL_SPEED] << endl;
   a_out << "  TRANSIT TIME          [s]: " << m_scale[TRANSIT_TIME] << endl;
   a_out << "  MAGNETIC MOMENT     [J/T]: " << m_scale[MAGNETIC_MOMENT] << endl;
   a_out << "  DIST. FUNCTION [s/(Am^6)]: " << m_scale[DISTRIBUTION_FUNCTION] << endl;
   a_out << "  POTENTIAL             [V]: " << m_scale[POTENTIAL] << endl;
   a_out << "  GYROFREQUENCY       [1/s]: " << m_scale[GYROFREQUENCY] << endl;
   a_out << "  GYRORADIUS            [m]: " << m_scale[GYRORADIUS] << endl;
   a_out << "  DEBYE LENGTH          [m]: " << m_scale[DEBYE_LENGTH] << endl;
   a_out << endl;
   a_out << "==== Dimensionless Parameters =========================================" << endl;
   a_out << endl;
   a_out << "  LARMOR NUMBER            : " << m_larmor_number << endl;
   a_out << "  DEBYE NUMBER             : " << m_debye_number << endl;
   a_out << endl;
   a_out << "***********************************************************************" << endl;
}

#include "NamespaceFooter.H"
