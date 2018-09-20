#ifndef _MOMENT_OP_H_
#define _MOMENT_OP_H_


#include "Dimensions.H"
#include "Kernels.H"

#undef CH_SPACEDIM
#include "FArrayBox.H.multidim"
#include "LevelData.H.multidim"
#define CH_SPACEDIM PDIM

#include "SliceSpec.H"

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;
namespace CP1 = CP1_NAMESPACE;


class KineticSpecies;

/// Moment Operator class.
/**
 * This is a singleton that provides the basic partial intergration
 * over phase space to produce moments of KineticSpecies distribution
 * functions.
 *
 * NB: This implementation is not thread safe.
 */
class MomentOp
{
   public:

      /// Aquire a reference.
      /**
       * Aquire a reference to the Singleton object.
       */
      static MomentOp& instance();

      /// Computes the moment.
      /**
       * @param[out] result Configuration-space LevelData into which
       * cell-averaged result is placed.  The component layout of the
       * result is determined by the kernel used to take the moment
       * @param[in] kinetic_species Kinetic species of which to take moment
       * @param[in] kernel Kernel object that evaluates the kernel
       */
      void compute( CFG::LevelData<CFG::FArrayBox>& result,
                    const KineticSpecies& kinetic_species,
                    const Kernel& kernel ) const;

      // modified compute to take moments of arbitrary LevelData<FArrayBox>
      void compute( CFG::LevelData<CFG::FArrayBox>& result,
                    const KineticSpecies& kinetic_species,
                    const LevelData<FArrayBox>& a_function,
                    const Kernel& kernel ) const;
   private:

      /// Default Constructor.
      /**
       */
      MomentOp() {;}

      /// Copy Constructor.
      /**
       * @param[in] arg MomentOp object to copy
       */
      MomentOp( const MomentOp& arg ) {;}

      /// Assignment operator.
      /**
       * @param[in] arg MomentOp object to copy
       */
      MomentOp& operator=( const MomentOp& rhs ) { return instance(); }

      /// Compute integrand.
      /**
       * Computes the integrand using the distribution function of the
       * kinetic species and the kernel function object specified at creation.
       *
       * @param[out] integrand LevelData into which cell-averaged result is
       * placed
       * @param[in] kinetic_species Kinetic species to use
       * @param[in] kernel Kernel object that evaluates the kernel
       */
      void computeIntegrand( LevelData<FArrayBox>&  integrand,
                             const KineticSpecies&  kinetic_species,
                             const Kernel&          kernel ) const;

      // modifed computeIntegrand to work for arbitrary LevelData<FArrayBox>
      void computeIntegrand( LevelData<FArrayBox>&       integrand,
                             const KineticSpecies&       kinetic_species,
                             const LevelData<FArrayBox>& function,
                             const Kernel&               kernel ) const;


      /// Partial integration over \f$\mu\f$ direction.
      /**
       * Partial integration over \f$\mu\f$ direction.
       *
       * @param[out] result Result of partial integration
       * @param[in] integrand Full dimension LevelData
       * @param[in] domain Phase space ProblemDomain
       */
      void partialIntegralMu( CP1::LevelData<CP1::FArrayBox>& result,
                              const LevelData<FArrayBox>& integrand,
                              const ProblemDomain& domain,
                              const SliceSpec& slice_mu) const;

      /// Partial integration over \f$v_{\parallel}\f$ direction.
      /**
       * Partial integration over \f$v_{\parallel}\f$ direction.
       *
       * @param[out] result Result of partial integration
       * @param[in] integrand Full dimension LevelData
       * @param[in] domain Phase space ProblemDomain
       */
      void partialIntegralVp( CFG::LevelData<CFG::FArrayBox>& result,
                              const CP1::LevelData<CP1::FArrayBox>& integrand,
                              const CP1::ProblemDomain& domain,
                              const CP1::SliceSpec& slice_vp) const;
};

#include "NamespaceFooter.H"

#endif
