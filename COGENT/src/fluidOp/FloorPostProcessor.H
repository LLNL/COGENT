#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _FLOORPOSTPROCESSOR_H_
#define _FLOORPOSTPROCESSOR_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluidSpecies.H"

#include "NamespaceHeader.H"

  /// Enforces floor values for fluid quantities
  /**
  */

class FloorPostProcessor
{
   public:

      /// Default Constructor.
      /**
       */
      FloorPostProcessor() :
         m_value(0),
         m_absolute(true),
         m_verbose(false)
      {;}

      /// Constructor.
      /**
       * @param[in] a_value  Floor value
       * @param[in] a_absolute If true/false use the floor value in the absolute/relative sence.
       */
      FloorPostProcessor(const Real& a_value,
                         const bool& a_absolute,
                         const bool& a_verbose = false);

      /// Destructor.
      /**
       */
      ~FloorPostProcessor()
      {
      }

      /// Define the enforcer.
      /**
       * @param[in] a_value  Floor value
       * @param[in] a_absolute If true/false use the floor value in the absolute/relative sence.
       */
      void define(const Real& a_value,
                  const bool& a_absolute,
                  const bool& a_verbose = false);

      /// Enforces the floor value through local redistribution of deficit.
      /**
       * @param[in,out] a_phi Scalar or vector cell-averaged field on which
       *   to enforce the ref value.
       * @param[in] a_ref_val Reference value used to estimate zero, i.e.
       *   the value that is to be considered O(1).
       */
      void enforce(FluidSpeciesPtrVect& a_data) const;

   private:

      void getMinMax(const LevelData<FArrayBox>& data,
                     double& min,
                     double& max) const;
   
      FloorPostProcessor( const FloorPostProcessor& foo );

      FloorPostProcessor& operator=( const FloorPostProcessor& foo );
   
      Real m_value;
      bool m_absolute;
      bool m_verbose;
};

#include "NamespaceFooter.H"

#endif
