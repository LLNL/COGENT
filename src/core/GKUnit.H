#ifndef  _GKUNITS_H_
#define  _GKUNITS_H_

#include <iostream>

#include "REAL.H"
#include "ParmParse.H"
#include "Interval.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

/**
 * Gyrokinetic units class.
 *
 * Responsible for converting between dimensionless and dimensional
 * representations.
 *
 * Required input keys and data types:
 *
 *    -\b number_density
 *      real value specifying the characteristic number density in [\f$m^{-3}\f$]
 *
 *    -\b temperature
 *      real value specifying the characteristic temperature in [\f$eV\f$]
 *
 *    -\b length
 *      real value specifying the characteristic length in [\f$m\f$]
 *
 *    -\b mass
 *      real value specifying the characteristic mass in [\f$AMU\f$]
 *
 *    -\b magnetic_field
 *      real value specifying the characteristic magnetic field in [\f$T\f$]
 *
 * The following represents a sample input entry:
 *
 * \verbatim
 *  units.number_density = 1.0
 *  units.temperature = 1.0
 *  units.length = 1.0
 *  units.mass = 1.0
 *  units.magnetic_field = 1.0
 * \endverbatim
*/
class GKUnits
{
   public:

      /// Enumerated characteristic scales.
      /**
       * The associate mks units are shown.
       */
      enum CHAR_SCALE {
         NUMBER_DENSITY = 0,    /**< [\f$m^{-3}\f$] */
         TEMPERATURE,           /**< [\f$eV\f$] */
         LENGTH,                /**< [\f$m\f$] */
         MASS,                  /**< [\f$AMU\f$] */
         MAGNETIC_FIELD,        /**< [\f$T\f$] = [\f$Wb m^{-2}\f$] = [\f$kg s^{-2} A^{-1}\f$] */
         CHARGE,                /**< [\f$C\f$] = [\f$s A\f$] */
         BOLTZMANN,             /**< [\f$J K^{-1}\f$] = [\f$m^2 kg s^{-2} K^{-1}\f$] */
         PERMITTIVITY,          /**< [\f$F m^{-1}\f$] = [\f$C (V m)^{-1}\f$] = [\f$s^4 A^2 m^{-2} kg^{-1}\f$]*/
         PERMEABILITY,          /**< [\f$H m^{-1}\f$] = [\f$Wb (A m)^{-1}\f$] = [\f$m^2 kg (A s)^{-2}\f$]*/
         THERMAL_SPEED,         /**< [\f$m s^{-1}\f$] */
         TRANSIT_TIME,          /**< [\f$s\f$] */
         MAGNETIC_MOMENT,       /**< [\f$J T^{-1}\f$] = [\f$J m^4 Wb^{-1}\f$] = [\f$A m^2\f$]*/
         DISTRIBUTION_FUNCTION, /**< [\f$s A^{-1} m^{-6}\f$] */
         POTENTIAL,             /**< [\f$V\f$] = [\f$J(As)^{-1}\f$] = [\f$kg m^2 A^{-1} s^{-2}\f$]*/
         GYROFREQUENCY,         /**< [\f$s^{-1}\f$] */
         GYRORADIUS,            /**< [\f$m\f$] */
         DEBYE_LENGTH,          /**< [\f$m\f$] */
         COUNT};

      /// Constructor.
      /**
       * Initializes the units.
       *
       *  @param[in] parm_parse the ParmParse database.
       */
      GKUnits( ParmParse& parm_parse );

      /// Destructor.
      /**
       */
      ~GKUnits() {;}

      /// Make data dimensionless.
      /**
       * Scales out the units from the components of data.
       *
       *  @param[in,out] data the data to be rescaled.
       *  @param[in] components the range of components to be rescales.
       *  @param[in] type the type of characteristic scale by which to normalize.
       */
      void makeDimensionless( LevelData<FArrayBox>& data,
                              const Interval&       components,
                              const CHAR_SCALE&     type ) const;

      /// Make data dimensionless.
      /**
       * Scales out the units from the components of data.
       *
       *  @param[in,out] data the data to be rescaled.
       *  @param[in] components the range of components to be rescales.
       *  @param[in] type the type of characteristic scale by which to normalize.
       */
      void makeDimensionless( LevelData<FluxBox>& data,
                              const Interval&     components,
                              const CHAR_SCALE&   type ) const;

      /// Make data dimensionless.
      /**
       * Scales out the units from the components of data.
       *
       *  @param[in,out] data the data to be rescaled.
       *  @param[in] type the type of characteristic scale by which to normalize.
       */
      void makeDimensionless( Real& data, const CHAR_SCALE& type ) const
      { checkType( type ); data /= m_scale[type]; }

      /// Make data dimensional.
      /**
       * Rescales the components of data to have mks units.
       *
       *  @param[in,out] data the data to be rescaled.
       *  @param[in] components the range of components to be rescales.
       *  @param[in] type the type of characteristic scale from which to provide units.
       */
      void makeDimensional( LevelData<FArrayBox>& data,
                            const Interval&       components,
                            const CHAR_SCALE&     type ) const;

      /// Make data dimensional.
      /**
       * Rescales the components of data to have mks units.
       *
       *  @param[in,out] data the data to be rescaled.
       *  @param[in] components the range of components to be rescales.
       *  @param[in] type the type of characteristic scale from which to provide units.
       */
      void makeDimensional( LevelData<FluxBox>& data,
                            const Interval&     components,
                            const CHAR_SCALE&   type ) const;

      /// Make data dimensional.
      /**
       * Rescales the components of data to have mks units.
       *
       *  @param[in,out] data the data to be rescaled.
       *  @param[in] type the type of characteristic scale from which to provide units.
       */
      void makeDimensional( Real& data, const CHAR_SCALE& type ) const
      { checkType( type ); data *= m_scale[type]; }

      /// Return the characterisitc Larmor number.
      /**
       * Returns the ratio of the characteristic Larmor radius to the
       * characteristic length.
       */
      Real larmorNumber() const { return m_larmor_number; }

      /// Return the charactersitic Debye number.
      /**
       * Returns the ratio of the characteristic Debye length to the
       * characteristic length.
       */
      Real debyeNumber() const { return m_debye_number; }

      /// Print charactersitic scales to output stream.
      /**
       *  @param[in] out the output stream.
       */
      void print( std::ostream& out ) const;

   private:

      Real m_scale[COUNT];
      Real m_larmor_number;
      Real m_debye_number;

      /// Check that the input type is valid.
      /**
       */
      void checkType( const CHAR_SCALE& type ) const
      { CH_assert( type>=0 && type<COUNT ) };
};

#include "NamespaceFooter.H"

#endif
